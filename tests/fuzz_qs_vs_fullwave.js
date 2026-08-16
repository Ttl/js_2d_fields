// Fuzzer: quasi-static (rectilinear FDM) vs full-wave (triangular FEM).
//
// Generates random transmission-line geometries across ALL UI line types and advanced
// features, solves each with BOTH backends, and flags:
//   • large discrepancies in the field-driven quantities (Z0, eps_eff, C) the two
//     methods should agree on,
//   • geometries whose triangular mesh comes out structurally poor (degenerate/NaN,
//     constraint crossings, missing constraint edges, extreme area ratio, >5% bad
//     triangles, or an extreme Q>2000 sliver),
//   • geometries that make ONE backend throw while the other succeeds (a backend bug).
//
// Coverage:
//   line types — microstrip, diff_microstrip, stripline, diff_stripline,
//                gcpw, diff_gcpw, broadside_stripline
//   features   — solder mask (incl. coplanar mask on gcpw; er 2–10, thicknesses
//                10–30 µm, coverage style: full conformal / trace top exposed /
//                no side coverage / trace-only), top dielectric (material from
//                the substrate range), ground cutout, enclosure (closed gnd box),
//                surface plating (σ 1e6–6e7 independent of the bulk — both
//                better- and worse-conducting; thickness 1–8 µm; random face
//                combo top/sides/bottom; plating rq up to 2 µm), surface
//                roughness (rq up to 3 µm), differential gaps, causal
//                (Djordjevic-Sarkar) dielectric dispersion
//   materials  — bulk conductor σ 1e6–6e7 S/m; substrate er 2.2–10 (broadside
//                layers 2.2–6), tanδ 0.001–0.02
//   modes      — single-ended lines compare the one quasi-TEM mode; differential pairs
//                compare BOTH the odd and the even mode (the worst of the two is flagged)
//
// Conductor loss (R per mode) IS compared (since 2026-08-16, after the multi-drive-MQS
// broadside reference and the skin-transition/plating/GCPW loss fixes): gate R_THRESH
// (default 25%), RELAXED to R_RELAX (default 60%) when a solve carries a loss-accuracy
// note — QS reason 'skin-transition' (δ deep into the conductor thickness) or
// 'broadside-proximity' (strongly coupled pair, QS reads up to ~50% high), or a tri
// mqs-* fallback warning (reference degraded to perturbation). Certificate warnings do
// NOT relax the gate (they describe C, not R). Geometries that BOTH backends reject
// (invalid random combos) are skipped, not flagged.
//
// Usage:   node tests/fuzz_qs_vs_fullwave.js [N] [seed] [threshold%] [R-threshold%]
//   e.g.   node tests/fuzz_qs_vs_fullwave.js 40 1 15 25
// Reproduce a flagged case: rerun with the same N and seed (generation is deterministic;
// NOTE the 2026-08-16 generator extension shifted every seed's stream — older seed/case
// references, incl. the ones in docs/backend_agreement_handover.md, no longer reproduce).

import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';

const N      = parseInt(process.argv[2]) || 40;
const SEED   = process.argv[3] !== undefined ? parseInt(process.argv[3]) : 1;
const THRESH = (parseFloat(process.argv[4]) || 15) / 100;
// Bad-mesh gate: structural errors (degenerate/NaN/constraint crossings/missing
// constraint edges/extreme area ratio), systemic quality collapse (>5% of
// triangles with Q>5 — checkMeshQuality's own warning threshold), or a single
// extreme sliver. The singleton cap was 100 and is now 2000: coarse-budget
// meshes on wide-domain thin-layer geometries (wide trace + solder mask)
// legitimately stretch far-field cap triangles to Q ~120-300 inside the 10-25 µm
// mask/trace bands, where the field is ~0 — verified benign (C/Z0 within 0.06%,
// R within 0.5% of a 60k sliver-free mesh, all structural checks clean) and
// self-healing with budget. The structural signals above are the ones a genuine
// mesh breakdown trips.
const BAD_Q  = 2000;
// Conductor-loss (R) gates: base, and relaxed for solves carrying a
// loss-accuracy note (see header). Calibrated 2026-08-16 on seeds 1-13
// (360 specs, app settings): unwarned median ~6% / max 17.2% (a low-coupling
// broadside; the sweep found the broadside bias crossing 26% at w/gap 1.70,
// which is why the proximity-warning threshold sits at 1.5); warned max 47.8%
// (deep-skin-transition plated stripline — QS ~2× high — and strongly-coupled
// broadside). New unwarned rows above 25% are worth investigating, not just
// gate-bumping: they'd extend the known bias envelope (issue-9 pour faces,
// VACUUM_LOSS_CAL tilt, broadside proximity below the warning threshold).
const R_THRESH = (parseFloat(process.argv[5]) || 25) / 100;
const R_RELAX  = 0.60;

// --- Seeded PRNG (deterministic, reproducible) ---
export function createRng(seed) {
    let s = (seed | 0) || 1;
    const next = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; };
    return {
        f: (lo, hi) => lo + next() * (hi - lo),
        logf: (lo, hi) => lo * Math.pow(hi / lo, next()),
        pick: arr => arr[Math.floor(next() * arr.length)],
        bool: p => next() < p,
    };
}

const TL_TYPES = ['microstrip', 'diff_microstrip', 'stripline', 'diff_stripline',
                  'gcpw', 'diff_gcpw', 'broadside_stripline'];

// --- Random spec generator (geometry + which features are enabled) ---
export function randomSpec(rng) {
    const tl = rng.pick(TL_TYPES);
    // rq up to 3 µm: at copper GHz skin depths that is rq/δ > 2, where the
    // roughness models diverge the most — exactly the regime worth comparing.
    const spec = { tl, freq: rng.logf(0.5e9, 6e9), rq: rng.bool(0.3) ? rng.logf(0.1e-6, 3e-6) : 0 };
    // Plating σ is drawn INDEPENDENTLY of the bulk σ over the same 1e6–6e7 range,
    // so plating lands both better- and worse-conducting than the bulk (silver
    // over copper vs nickel over copper); thickness spans thinner and thicker
    // than typical GHz skin depths, exercising the layered-impedance mix.
    const rollPlating = () => {
        if (!rng.bool(0.2)) return false;
        spec.plating_sigma = rng.logf(1e6, 6e7);
        spec.plating_t = rng.logf(1e-6, 8e-6);
        // Face combo (top/sides/bottom independent, at least one guaranteed) and
        // plated-surface roughness — both were fixed (top+sides, rq 0) before
        // 2026-08-16; the app exposes all of them.
        spec.plating_top = rng.bool(0.8);
        spec.plating_sides = rng.bool(0.7);
        spec.plating_bottom = rng.bool(0.25);
        if (!spec.plating_top && !spec.plating_sides && !spec.plating_bottom) spec.plating_top = true;
        spec.plating_rq = rng.bool(0.4) ? rng.logf(0.1e-6, 2e-6) : 0;
        return true;
    };
    // Causal (Djordjevic-Sarkar) dielectric dispersion — applies to every line type and to
    // both backends, so it stays inside the cross-backend comparison (it shifts er(f), hence
    // C and eps_eff, which the fuzzer DOES compare). See solveOn for how it is evaluated.
    // Kept a modest fraction (not 0.3+): a causal full-wave case is ~20× slower than a plain
    // one because the triangular backend bypasses its eigensolve cache on every refinement
    // pass when causal is on, so too many causal cases make a 100-run impractically long.
    spec.use_causal = rng.bool(0.2);

    if (tl === 'broadside_stripline') {
        const bs_w = rng.logf(0.05e-3, 1e-3), bs_t = rng.logf(10e-6, 50e-6);
        Object.assign(spec, {
            bs_w, bs_t,
            bs_x_offset: rng.bool(0.5) ? rng.f(-0.15e-3, 0.15e-3) : 0,
            bs_h_bottom: rng.logf(0.05e-3, 0.5e-3), bs_er_bottom: rng.f(2.2, 6), bs_tand_bottom: rng.f(0.001, 0.02),
            // Floor of 2.5t: the middle layer must clear 2× trace_thickness or the
            // broadside pair collides and validation rejects the case on both backends.
            bs_h_middle: rng.logf(Math.max(0.05e-3, 2.5 * bs_t), 0.5e-3), bs_er_middle: rng.f(2.2, 6), bs_tand_middle: rng.f(0.001, 0.02),
            bs_h_top: rng.logf(0.05e-3, 0.5e-3), bs_er_top: rng.f(2.2, 6), bs_tand_top: rng.f(0.001, 0.02),
            bs_sigma: rng.logf(1e6, 6e7),
        });
        spec.use_enclosure = rng.bool(0.3);
        // An enclosure always closes its side walls (see the enclosure block below).
        spec.use_side_gnd = spec.use_enclosure;
        if (spec.use_enclosure) spec.enclosure_width = spec.bs_w * rng.f(3, 8);
        spec.use_plating = rollPlating();
        return spec;
    }

    // Bulk conductor σ spans 1e6–6e7 S/m (poor alloys up to copper): at fixed
    // frequency this sweeps δ ~8×, walking the skin-transition/thickness ratio.
    const w = rng.logf(0.05e-3, 3e-3), h = rng.logf(0.05e-3, 1.5e-3), t = rng.logf(10e-6, 70e-6);
    Object.assign(spec, {
        w, h, t, er: rng.f(2.2, 10), tand: rng.f(0.001, 0.02), sigma: rng.logf(1e6, 6e7),
        gnd_thickness: rng.logf(10e-6, 35e-6),
    });
    const isDiff = tl.includes('diff');
    if (isDiff) spec.trace_spacing = rng.logf(0.02e-3, 0.6e-3);
    // Top height floor of 2t: a cover thinner than the trace always fails
    // parameter validation (trace_thickness > enclosure_height) on both backends.
    if (tl.includes('stripline')) { spec.er_top = spec.er; spec.tand_top = spec.tand; spec.stripline_top_h = rng.logf(Math.max(0.05e-3, 2 * t), 1.5e-3); }
    if (tl.includes('gcpw')) { spec.gap = rng.logf(0.05e-3, 0.5e-3); spec.via_gap = rng.logf(0.05e-3, 0.5e-3); }

    // Advanced features. Solder mask applies to the microstrip AND gcpw families
    // (gcpw takes the coplanar solder mask path in MicrostripSolver); its material
    // and thicknesses are randomized — er well past the usual 3.5 to make the
    // C/eps shift it causes clearly visible in the comparison. Top dielectric /
    // ground cutout stay microstrip-family-only, as in the app.
    if (tl === 'microstrip' || tl === 'diff_microstrip' || tl === 'gcpw' || tl === 'diff_gcpw') {
        spec.use_sm = rng.bool(0.3);
        if (spec.use_sm) {
            spec.sm_er = rng.f(2, 10); spec.sm_tand = rng.f(0.001, 0.03);
            spec.sm_t_sub = rng.f(10e-6, 30e-6); spec.sm_t_trace = rng.f(10e-6, 30e-6);
            spec.sm_t_side = rng.f(10e-6, 30e-6);
            // Coverage style (was always full conformal): half the time zero one
            // component — trace top exposed, no side coverage, or trace-only
            // (bare substrate). Both backends verified to handle the degenerate
            // (zero-thickness) mask rects.
            spec.sm_style = rng.pick(['full', 'full', 'full', 'no-trace-top', 'no-sides', 'trace-only']);
            if (spec.sm_style === 'no-trace-top') spec.sm_t_trace = 0;
            else if (spec.sm_style === 'no-sides') spec.sm_t_side = 0;
            else if (spec.sm_style === 'trace-only') spec.sm_t_sub = 0;
        }
    }
    if (tl === 'microstrip' || tl === 'diff_microstrip') {
        // Top dielectric material drawn from the same range as the substrate.
        if (rng.bool(0.3)) {
            spec.use_top_diel = true; spec.top_diel_h = rng.logf(0.05e-3, 0.5e-3);
            spec.top_diel_er = rng.f(2.2, 10); spec.top_diel_tand = rng.f(0.001, 0.02);
        }
        if (rng.bool(0.2)) { spec.use_gnd_cut = true; spec.gnd_cut_w = w * rng.f(0.2, 1.5); spec.gnd_cut_h = h * rng.f(0.1, 0.7); }
    }
    // Enclosure walls apply to any of these. An enclosure is always a CLOSED box
    // (gnd sides + gnd top): it shrinks the domain to a few trace widths, and an
    // OPEN boundary that close to the trace is a domain truncation the two
    // backends legitimately approximate differently (the tri backend's natural-BC
    // walls act as magnetic mirrors on a tight box and inflate C ~1/width, while
    // the QS open treatment is near-transparent) — a modeling artifact, not a
    // backend bug, so those combos are excluded by construction instead of
    // surfacing as spurious discrepancies. Open boundaries stay covered by the
    // non-enclosure cases, whose auto-sized domains keep the walls far away.
    if (rng.bool(0.3)) {
        spec.use_enclosure = true;
        let span = isDiff ? (2 * w + spec.trace_spacing) : w;
        // GCPW: the enclosure must clear the full coplanar active width (trace +
        // gaps + via fences, matching MicrostripSolver's active_width), not just
        // the trace span — otherwise the combo always fails validation.
        if (tl.includes('gcpw')) span += (isDiff ? 0 : w) + 2 * (spec.gap + spec.via_gap);
        spec.enclosure_width = span * rng.f(2.5, 6) + 2 * spec.gnd_thickness;
        spec.enclosure_height = (h + t) * rng.f(1.5, 4);
        spec.use_side_gnd = true;
        spec.use_top_gnd = true;
    }
    spec.use_plating = rollPlating();
    return spec;
}

// --- Build solver options for a spec on a backend (mirrors app_solver.updateGeometry) ---
function addCommon(o, spec) {
    if (spec.use_sm) {
        o.use_sm = true;
        o.sm_t_sub = spec.sm_t_sub; o.sm_t_trace = spec.sm_t_trace; o.sm_t_side = spec.sm_t_side;
        o.sm_er = spec.sm_er; o.sm_tand = spec.sm_tand;
    }
    if (spec.use_top_diel) { o.top_diel_h = spec.top_diel_h; o.top_diel_er = spec.top_diel_er; o.top_diel_tand = spec.top_diel_tand; }
    if (spec.use_gnd_cut) { o.gnd_cut_width = spec.gnd_cut_w; o.gnd_cut_sub_h = spec.gnd_cut_h; }
    if (spec.use_enclosure) {
        o.enclosure_width = spec.enclosure_width;
        if (o.enclosure_height === undefined) o.enclosure_height = spec.enclosure_height;
        const lr = spec.use_side_gnd ? 'gnd' : 'open';
        const top = spec.use_top_gnd ? 'gnd' : (o.boundaries ? o.boundaries[2] : 'open');
        const bot = o.boundaries ? o.boundaries[3] : 'gnd';
        o.boundaries = [lr, lr, top, bot];
    }
    if (spec.use_plating) o.plating = { sigma: spec.plating_sigma, thickness: spec.plating_t, rq: spec.plating_rq ?? 0,
        top: spec.plating_top ?? true, sides: spec.plating_sides ?? true, bottom: spec.plating_bottom ?? false, thick_corners: true };
}

export function buildSolver(spec, backend) {
    if (spec.tl === 'broadside_stripline') {
        const o = {
            trace_width: spec.bs_w, trace_thickness: spec.bs_t, x_offset: spec.bs_x_offset, sigma_cond: spec.bs_sigma,
            h_bottom: spec.bs_h_bottom, er_bottom: spec.bs_er_bottom, tand_bottom: spec.bs_tand_bottom,
            h_middle: spec.bs_h_middle, er_middle: spec.bs_er_middle, tand_middle: spec.bs_tand_middle,
            h_top: spec.bs_h_top, er_top: spec.bs_er_top, tand_top: spec.bs_tand_top,
            freq: spec.freq, rq: spec.rq, boundaries: ['open', 'open', 'gnd', 'gnd'], mesh_backend: backend,
            nx: 30, ny: 30,   // coarse initial grid (matches the app) so adaptive refinement has headroom
        };
        if (spec.use_enclosure) { o.enclosure_width = spec.enclosure_width; if (spec.use_side_gnd) o.boundaries = ['gnd', 'gnd', 'gnd', 'gnd']; }
        if (spec.use_plating) o.plating = { sigma: spec.plating_sigma, thickness: spec.plating_t, rq: spec.plating_rq ?? 0,
        top: spec.plating_top ?? true, sides: spec.plating_sides ?? true, bottom: spec.plating_bottom ?? false, thick_corners: true };
        const s = new BroadsideStriplineSolver(o);
        s.use_causal_materials = !!spec.use_causal;
        return s;
    }
    const base = {
        trace_width: spec.w, substrate_height: spec.h, trace_thickness: spec.t, gnd_thickness: spec.gnd_thickness,
        epsilon_r: spec.er, tan_delta: spec.tand, sigma_cond: spec.sigma, freq: spec.freq, rq: spec.rq,
        mesh_backend: backend,
        nx: 30, ny: 30,   // coarse initial grid (matches the app) so adaptive refinement has headroom
    };
    const o = { ...base };
    if (spec.trace_spacing) o.trace_spacing = spec.trace_spacing;
    if (spec.tl.includes('gcpw')) { o.boundaries = ['open', 'open', 'open', 'gnd']; o.use_coplanar_gnd = true; o.gap = spec.gap; o.via_gap = spec.via_gap; o.use_vias = true; }
    else if (spec.tl.includes('stripline')) { o.epsilon_r_top = spec.er_top; o.tan_delta_top = spec.tand_top; o.enclosure_height = spec.stripline_top_h; o.boundaries = ['open', 'open', 'gnd', 'gnd']; }
    else o.boundaries = ['open', 'open', 'open', 'gnd'];
    addCommon(o, spec);
    const s = new MicrostripSolver(o);
    s.use_causal_materials = !!spec.use_causal;
    return s;
}

async function solveOn(spec, backend) {
    const log = console.log;
    console.log = () => {};
    try {
        const s = buildSolver(spec, backend);
        // Mirror the app's adaptive controls (src/app_solver.js) so the fuzzer measures
        // real app behaviour, not an unrefined coarse grid. Without these the rectilinear
        // backend inherited nx=300 (a >90k-node initial grid that exceeds the 20k budget →
        // zero refinement) and min_converged_passes=1 (stops at the first stably-wrong pass),
        // which manufactured spurious QS-vs-FW discrepancies on tight differential gaps.
        let r = await s.solve_adaptive({
            max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2,
        });
        // Causal materials: solve_adaptive refines the mesh with NON-causal er, so — exactly as
        // the app does (app_solver.runSimulation re-evaluates the solve point through
        // computeAtFrequency when use_causal_materials is set) — re-evaluate at the solve
        // frequency to apply the Djordjevic-Sarkar model. This applies the SAME er(f) shift to
        // both backends (rectilinear re-solves Laplace with the causal er; triangular re-runs
        // solveAt → _applyCausal), keeping the cross-backend comparison fair.
        if (spec.use_causal) r = await s.computeAtFrequency(s.freq, r);
        // Compare every mode the backend returns: single-ended → [single]; differential →
        // [odd, even], in that fixed order from the shared FieldSolver2D path on both backends.
        return {
            modes: r.modes.map(m => ({ Z0: m.Z0, eps_eff: m.eps_eff, C: m.RLGC.C, R: m.RLGC.R, mode: m.mode })),
            // Warning tags for the R-gate relaxation: QS loss-accuracy reasons
            // ('skin-transition' / 'broadside-proximity') and tri fallback types.
            warns: (r.warnings || []).map(w => ({ type: w.type, reason: w.reason || null })),
            maxQ: s.meshQuality ? s.meshQuality.maxQ : null,
            degenerate: s.meshQuality ? s.meshQuality.degenerateCount : 0,
            nan: s.meshQuality ? s.meshQuality.nanNodes : 0,
            // Structural quality signals for the bad-mesh gate (see BAD_Q note).
            crossings: s.meshQuality ? s.meshQuality.crossings : 0,
            missingEdges: s.meshQuality ? s.meshQuality.missingEdges : 0,
            areaRatio: s.meshQuality ? s.meshQuality.areaRatio : 0,
            badFraction: s.meshQuality ? s.meshQuality.badFraction : 0,
        };
    } catch (e) {
        return { error: e.message || String(e) };
    } finally {
        console.log = log;
    }
}

function fmtSpec(spec) {
    const u = (x, s = 1e3, d = 3) => (x * s).toFixed(d);
    let s = spec.tl;
    if (spec.tl === 'broadside_stripline') s += ` w=${u(spec.bs_w)}mm t=${u(spec.bs_t, 1e6, 1)}µm off=${u(spec.bs_x_offset)}mm hM=${u(spec.bs_h_middle)}mm`;
    else s += ` w=${u(spec.w)}mm h=${u(spec.h)}mm t=${u(spec.t, 1e6, 1)}µm er=${spec.er.toFixed(2)}`;
    s += ` sig=${((spec.sigma ?? spec.bs_sigma) / 1e6).toFixed(1)}e6 f=${(spec.freq / 1e9).toFixed(2)}GHz`;
    if (spec.trace_spacing) s += ` gap=${u(spec.trace_spacing)}mm`;
    const feats = ['use_sm', 'use_top_diel', 'use_gnd_cut', 'use_enclosure', 'use_side_gnd', 'use_top_gnd', 'use_plating', 'use_causal']
        .filter(k => spec[k]).map(k => {
            if (k === 'use_sm' && spec.sm_style && spec.sm_style !== 'full') return `sm:${spec.sm_style}`;
            if (k === 'use_plating') {
                const faces = [spec.plating_top && 't', spec.plating_sides && 's', spec.plating_bottom && 'b']
                    .filter(Boolean).join('');
                return `plating:${faces || 'ts'}${spec.plating_rq ? '+rq' : ''}`;
            }
            return k.replace('use_', '');
        });
    if (feats.length) s += ` [${feats.join(',')}]`;
    return s;
}

const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);

async function main() {
    console.log(`\n### Fuzz QS vs full-wave — N=${N} seed=${SEED} flag>${(THRESH * 100).toFixed(0)}% ` +
        `R>${(R_THRESH * 100).toFixed(0)}% (warned ${(R_RELAX * 100).toFixed(0)}%) Qmax>${BAD_Q} ###`);
    console.log(`covers: ${TL_TYPES.join(', ')} + solder-mask/top-diel/gnd-cut/enclosure/plating/roughness/causal; diff pairs check odd+even\n`);
    const rng = createRng(SEED);
    const flagged = { discrepancy: [], loss: [], badMesh: [], crash: [] };
    let ok = 0, skipped = 0, rejected = 0;
    const byType = {};
    for (let i = 0; i < N; i++) {
        const spec = randomSpec(rng);
        byType[spec.tl] = (byType[spec.tl] || 0) + 1;
        const qs = await solveOn(spec, 'rectilinear');
        const fw = await solveOn(spec, 'triangular');

        // Both reject → invalid random geometry (validation is backend-independent); skip.
        if (qs.error && fw.error) { skipped++; continue; }
        // Only one backend fails → a real backend bug — UNLESS it's the pre-mesh
        // meshability guard cleanly refusing a geometry one backend's cost model
        // can't afford (working as intended; the guard has its own dedicated test).
        if (qs.error || fw.error) {
            const err = String(qs.error || fw.error);
            const who = qs.error ? `quasi-static threw (full-wave OK): ${qs.error}` : `full-wave threw (quasi-static OK): ${fw.error}`;
            if (/cannot be meshed/i.test(err)) {
                rejected++;
                console.log(`[${i}] ○ MESHABILITY REJECTED (${qs.error ? 'quasi-static' : 'full-wave'} guard)\n      ${fmtSpec(spec)}`);
            } else {
                flagged.crash.push({ spec, who });
                console.log(`[${i}] ✗ ONE-SIDED FAILURE — ${who}\n      ${fmtSpec(spec)}`);
            }
            continue;
        }

        // areaRatio is deliberately NOT gated: adaptive skin-band refinement
        // routinely produces max/min element-area ratios of 1e6-1e7 (µm-scale
        // skin elements vs mm-scale far field) on perfectly healthy meshes —
        // measured on cases with Qmax 3-6 and badFraction ~0 — and the direct
        // factorization handles it; the ratio separates nothing here.
        const badMesh = (fw.maxQ != null && fw.maxQ > BAD_Q) || fw.degenerate > 0 || fw.nan > 0
            || fw.crossings > 0 || fw.missingEdges > 0 || fw.badFraction > 0.05;
        if (badMesh) {
            flagged.badMesh.push({ spec, maxQ: fw.maxQ });
            console.log(`[${i}] ✗ BAD MESH (Qmax=${fw.maxQ?.toFixed(0)} degen=${fw.degenerate} nan=${fw.nan} ` +
                `cross=${fw.crossings} missEdge=${fw.missingEdges} areaRatio=${fw.areaRatio?.toExponential(1)} ` +
                `badFrac=${(fw.badFraction * 100).toFixed(1)}%)\n      ${fmtSpec(spec)}`);
        }

        // Flag on the field-driven quantities Z0 and C only. eps_eff is reported with
        // DIFFERENT conventions by the two backends — the full-wave includes the
        // conductor internal inductance (eps_eff = c²·L·C), the quasi-static reports the
        // pure C/C0 — so it differs by the (loss-dependent) L_internal even when the
        // electrostatic field agrees. That's the loss-modeling difference the fuzzer is
        // meant to exclude; eps is still printed for context.
        // Compare each matched mode pair and flag the worst. Single-ended lines have one
        // quasi-TEM mode; differential pairs have two (odd then even), so the even mode — which
        // the old fuzzer never checked — is now validated cross-backend alongside the odd.
        const nModes = Math.min(qs.modes.length, fw.modes.length);
        let worst = null;
        for (let mi = 0; mi < nModes; mi++) {
            const q = qs.modes[mi], w = fw.modes[mi];
            const dZ = relDiff(q.Z0, w.Z0), dE = relDiff(q.eps_eff, w.eps_eff), dC = relDiff(q.C, w.C);
            const dMax = Math.max(dZ, dC);
            if (!worst || dMax > worst.dMax) worst = { mode: q.mode, dZ, dE, dC, dMax, q, w };
        }
        if (worst && worst.dMax > THRESH) {
            flagged.discrepancy.push({ i, spec, mode: worst.mode, dZ: worst.dZ, dE: worst.dE, dC: worst.dC });
            const tag = nModes > 1 ? ` [${worst.mode}]` : '';
            console.log(`[${i}] ✗ DISCREPANCY ${(worst.dMax * 100).toFixed(0)}%${tag} (Z0 ${(worst.dZ * 100).toFixed(0)}% eps ${(worst.dE * 100).toFixed(0)}% C ${(worst.dC * 100).toFixed(0)}%)\n` +
                `      ${fmtSpec(spec)}\n` +
                `      qs: Z0=${worst.q.Z0.toFixed(2)} eps=${worst.q.eps_eff.toFixed(3)} C=${(worst.q.C * 1e12).toFixed(1)}pF | ` +
                `fw: Z0=${worst.w.Z0.toFixed(2)} eps=${worst.w.eps_eff.toFixed(3)} C=${(worst.w.C * 1e12).toFixed(1)}pF`);
        }

        // Conductor-loss gate (R per mode). Relaxed when either backend flagged
        // reduced loss accuracy for this solve: QS reasons 'skin-transition' /
        // 'broadside-proximity', or the tri backend downgrading MQS to
        // perturbation (mqs-* fallback warnings). Certificate warnings describe
        // C, not R, and do NOT relax the gate.
        const lossRelaxed =
            qs.warns.some(w => w.reason === 'skin-transition' || w.reason === 'broadside-proximity')
            || fw.warns.some(w => /^mqs-/.test(w.type));
        const rGate = lossRelaxed ? R_RELAX : R_THRESH;
        let worstR = null;
        for (let mi = 0; mi < nModes; mi++) {
            const q = qs.modes[mi], w = fw.modes[mi];
            if (!(q.R > 0) || !(w.R > 0)) continue;
            const dR = relDiff(q.R, w.R);
            if (!worstR || dR > worstR.dR) worstR = { mode: q.mode, dR, q, w };
        }
        if (worstR && worstR.dR > rGate) {
            flagged.loss.push({ i, spec, mode: worstR.mode, dR: worstR.dR, relaxed: lossRelaxed });
            const tag = nModes > 1 ? ` [${worstR.mode}]` : '';
            console.log(`[${i}] ✗ LOSS DISCREPANCY ${(worstR.dR * 100).toFixed(0)}%${tag}` +
                `${lossRelaxed ? ' (relaxed gate)' : ''}\n` +
                `      ${fmtSpec(spec)}\n` +
                `      R qs=${worstR.q.R.toFixed(2)} vs fw=${worstR.w.R.toFixed(2)} Ω/m`);
        }
        if (!badMesh && (!worst || worst.dMax <= THRESH) && (!worstR || worstR.dR <= rGate)) ok++;
    }

    console.log(`\n${'='.repeat(72)}`);
    console.log(`coverage: ${Object.entries(byType).map(([k, v]) => `${k}:${v}`).join('  ')}`);
    // Split the loss failures by gate: an UNWARNED one is the interesting kind — it
    // extends the known bias envelope rather than re-confirming a regime both backends
    // already flag as approximate (see the R_THRESH note above).
    const lossUnwarned = flagged.loss.filter(f => !f.relaxed).length;
    const lossWarned = flagged.loss.length - lossUnwarned;
    console.log(`SUMMARY: ${ok}/${N} clean | ${flagged.discrepancy.length} discrepancy | ` +
        `${flagged.loss.length} loss-discrepancy (${lossUnwarned} unwarned, ${lossWarned} warned) | ` +
        `${flagged.badMesh.length} bad-mesh | ` +
        `${flagged.crash.length} one-sided-fail | ${rejected} meshability-rejected | ${skipped} skipped(invalid)`);
    if (lossUnwarned) {
        console.log(`unwarned loss rows (investigate, do not just raise the gate):`);
        for (const f of flagged.loss.filter(x => !x.relaxed))
            console.log(`  [${f.i}] ${(100 * f.dR).toFixed(0)}%${f.mode ? ` [${f.mode}]` : ''}  ${fmtSpec(f.spec)}`);
    }
    console.log(`${'='.repeat(72)}`);
    if (flagged.discrepancy.length || flagged.loss.length || flagged.badMesh.length || flagged.crash.length) process.exitCode = 1;
}

import { pathToFileURL } from 'url';
if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) main();
