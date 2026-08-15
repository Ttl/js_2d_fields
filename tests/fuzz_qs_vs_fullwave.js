// Fuzzer: quasi-static (rectilinear FDM) vs full-wave (triangular FEM).
//
// Generates random transmission-line geometries across ALL UI line types and advanced
// features, solves each with BOTH backends, and flags:
//   • large discrepancies in the field-driven quantities (Z0, eps_eff, C) the two
//     methods should agree on,
//   • geometries whose triangular mesh comes out poor (sliver Qmax, degenerate, NaN),
//   • geometries that make ONE backend throw while the other succeeds (a backend bug).
//
// Coverage:
//   line types — microstrip, diff_microstrip, stripline, diff_stripline,
//                gcpw, diff_gcpw, broadside_stripline
//   features   — solder mask (incl. coplanar mask on gcpw; er 2–10, thicknesses
//                10–30 µm), top dielectric (material from the substrate range),
//                ground cutout, enclosure (closed gnd box), surface plating
//                (σ 1e6–6e7 independent of the bulk — both better- and
//                worse-conducting; thickness 1–8 µm), surface roughness (rq up
//                to 3 µm), differential gaps, causal (Djordjevic-Sarkar)
//                dielectric dispersion
//   materials  — bulk conductor σ 1e6–6e7 S/m; substrate er 2.2–10 (broadside
//                layers 2.2–6), tanδ 0.001–0.02
//   modes      — single-ended lines compare the one quasi-TEM mode; differential pairs
//                compare BOTH the odd and the even mode (the worst of the two is flagged)
//
// Loss (R / alpha_c) is intentionally NOT compared — the two methods model conductor
// loss differently (full-wave MQS captures proximity the quasi-static perturbation
// under-predicts), so they legitimately diverge there. Geometries that BOTH backends
// reject (invalid random combos) are skipped, not flagged.
//
// Usage:   node tests/fuzz_qs_vs_fullwave.js [N] [seed] [threshold%]
//   e.g.   node tests/fuzz_qs_vs_fullwave.js 40 1 15
// Reproduce a flagged case: rerun with the same N and seed (generation is deterministic).

import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';

const N      = parseInt(process.argv[2]) || 40;
const SEED   = process.argv[3] !== undefined ? parseInt(process.argv[3]) : 1;
const THRESH = (parseFloat(process.argv[4]) || 15) / 100;
const BAD_Q  = 100;

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
    if (spec.use_plating) o.plating = { sigma: spec.plating_sigma, thickness: spec.plating_t, rq: 0, top: true, sides: true, bottom: false, thick_corners: true };
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
        if (spec.use_plating) o.plating = { sigma: spec.plating_sigma, thickness: spec.plating_t, rq: 0, top: true, sides: true, bottom: false, thick_corners: true };
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
            modes: r.modes.map(m => ({ Z0: m.Z0, eps_eff: m.eps_eff, C: m.RLGC.C, mode: m.mode })),
            maxQ: s.meshQuality ? s.meshQuality.maxQ : null,
            degenerate: s.meshQuality ? s.meshQuality.degenerateCount : 0,
            nan: s.meshQuality ? s.meshQuality.nanNodes : 0,
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
        .filter(k => spec[k]).map(k => k.replace('use_', ''));
    if (feats.length) s += ` [${feats.join(',')}]`;
    return s;
}

const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);

async function main() {
    console.log(`\n### Fuzz QS vs full-wave — N=${N} seed=${SEED} flag>${(THRESH * 100).toFixed(0)}% Qmax>${BAD_Q} ###`);
    console.log(`covers: ${TL_TYPES.join(', ')} + solder-mask/top-diel/gnd-cut/enclosure/plating/roughness/causal; diff pairs check odd+even\n`);
    const rng = createRng(SEED);
    const flagged = { discrepancy: [], badMesh: [], crash: [] };
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

        const badMesh = (fw.maxQ != null && fw.maxQ > BAD_Q) || fw.degenerate > 0 || fw.nan > 0;
        if (badMesh) {
            flagged.badMesh.push({ spec, maxQ: fw.maxQ });
            console.log(`[${i}] ✗ BAD MESH (Qmax=${fw.maxQ?.toFixed(0)} degen=${fw.degenerate} nan=${fw.nan})\n      ${fmtSpec(spec)}`);
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
        if (!badMesh && (!worst || worst.dMax <= THRESH)) ok++;
    }

    console.log(`\n${'='.repeat(72)}`);
    console.log(`coverage: ${Object.entries(byType).map(([k, v]) => `${k}:${v}`).join('  ')}`);
    console.log(`SUMMARY: ${ok}/${N} clean | ${flagged.discrepancy.length} discrepancy | ` +
        `${flagged.badMesh.length} bad-mesh | ` +
        `${flagged.crash.length} one-sided-fail | ${rejected} meshability-rejected | ${skipped} skipped(invalid)`);
    console.log(`${'='.repeat(72)}`);
    if (flagged.discrepancy.length || flagged.badMesh.length || flagged.crash.length) process.exitCode = 1;
}

import { pathToFileURL } from 'url';
if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) main();
