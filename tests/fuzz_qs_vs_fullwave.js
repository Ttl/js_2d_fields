// Fuzzer: quasi-static (rectilinear FDM) vs full-wave (triangular FEM).
//
// Generates random transmission-line geometries across ALL UI line types and advanced
// features, solves each with BOTH backends, and flags:
//   • large discrepancies in the field-driven quantities (Z0, eps_eff, C) the two
//     methods should agree on,
//   • geometries whose triangular mesh comes out structurally poor (degenerate/NaN,
//     >5% bad triangles, or an extreme Q>2000 sliver),
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
//   materials  — bulk conductor σ 5e5–6.3e7 S/m; substrate er 2.2–10 (broadside
//                layers 2.2–6) with ~10% low-eps draws (1–2.2), tanδ 0.001–0.02
//                with ~10% lossy draws (0.02–0.2)
//   conductors — thickness 10–70 µm with ~15% thin films (0.5–10 µm); embedded
//                (negative-thickness) traces on every family incl. broadside
//   frequency  — 0.5–6 GHz, and ~30% of draws at 1–500 MHz (skin depth at and
//                beyond the conductor thickness)
//   modes      — single-ended lines compare the one quasi-TEM mode; differential pairs
//                compare BOTH the odd and the even mode (the worst of the two is flagged)
//
// Conductor loss (R per mode) is compared: gate R_THRESH
// (default 25%), RELAXED to R_RELAX (default 60%) when a solve carries a loss-accuracy
// note. QS reason 'skin-transition' (δ deep into the conductor thickness) or
// 'narrow-gap', or any tri
// mqs-* warning: either a fallback (reference degraded to perturbation) or
// 'mqs-band-capped' (MQS ran, but on a skin band that stopped short of its size
// target). Certificate warnings do NOT relax the gate (they describe C, not R).
// Geometries that BOTH backends reject (invalid random combos) are skipped, not flagged.
//
// Usage:   node tests/fuzz_qs_vs_fullwave.js [N] [seed] [threshold%] [R-threshold%]
//   e.g.   node tests/fuzz_qs_vs_fullwave.js 40 1 15 25
// Parallel: node tests/fuzz_parallel.mjs [-j jobs] [N] [seed] [threshold%] [R-threshold%]
//   (same cases and output, the cases spread over worker threads)
// Reproduce a flagged case: rerun with the same N and seed (generation is deterministic;

import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';
import { CustomGeometrySolver } from '../src/custom_geometry.js';

// Run settings from the command line [N] [seed] [threshold%] [R-threshold%], and
// ONLY=4,11, which solves just those case indices (the others are drawn and skipped).
export function parseArgs(args) {
    return {
        N: parseInt(args[0]) || 40,
        seed: args[1] !== undefined ? parseInt(args[1]) : 1,
        thresh: (parseFloat(args[2]) || 15) / 100,
        rThresh: (parseFloat(args[3]) || 25) / 100,
        only: process.env.ONLY ? new Set(process.env.ONLY.split(',').map(Number)) : null,
    };
}
// Bad-mesh gate: structural errors (degenerate/NaN), systemic quality collapse (>5% of
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
// (360 specs, app settings): unwarned median ~6% / max 17.2%; warned max 47.8%
// (deep-skin-transition plated stripline, QS ~2x high). New unwarned rows above
// 25% are worth investigating, not just gate-bumping: they'd extend the known
// bias envelope (issue-9 pour faces, strongly coupled broadside pairs). The base
// gate is a command-line setting (rThresh).
const R_RELAX  = 0.60;
// Shunt-conductance and inductance gates (see the G/L block in main).
const G_THRESH = 0.15;
const L_THRESH = 0.10;
// Calibrated 2026-09-19 on seeds 201-206 (240 specs, 30% below 500 MHz, 15% thin
// films): every L difference above 10% carried a loss-accuracy note, max 28%.
const L_RELAX  = 0.35;

// --- Seeded PRNG (deterministic, reproducible) ---
export function createRng(seed) {
    // mulberry32. The old LCG multiplied in doubles, which dropped the low bits of
    // s * 1103515245 and collapsed every seed onto one 10466-draw cycle (~363 specs).
    let s = seed >>> 0;
    const next = () => {
        s = (s + 0x6D2B79F5) >>> 0;
        let z = Math.imul(s ^ (s >>> 15), s | 1);
        z ^= z + Math.imul(z ^ (z >>> 7), z | 61);
        return ((z ^ (z >>> 14)) >>> 0) / 4294967296;
    };
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
    // Frequency: the GHz band, and in ~30% of draws 1-500 MHz, where the skin depth
    // reaches and passes the conductor thickness for every sigma drawn (35 um copper
    // crosses at ~3.5 MHz): the regime the internal-inductance and R_dc models own.
    const freq = rng.bool(0.3) ? rng.logf(1e6, 0.5e9) : rng.logf(0.5e9, 6e9);
    const spec = { tl, freq, rq: rng.bool(0.3) ? rng.logf(0.1e-6, 3e-6) : 0 };
    // Shared material draws. Substrates: the usual laminate range, ~10% low-eps
    // (foam / air, eps_r 1-2.2). Loss tangent: ~10% lossy (0.02-0.2), so G is compared
    // where it is not a small number. Conductors: resistive alloys up to silver.
    const drawEr = (hi) => rng.bool(0.1) ? rng.f(1, 2.2) : rng.f(2.2, hi);
    const drawTand = () => rng.bool(0.1) ? rng.logf(0.02, 0.2) : rng.f(0.001, 0.02);
    const drawSigma = () => rng.logf(5e5, 6.3e7);
    // Conductor thickness: foil weights, and in ~15% of draws a thin film (0.5-10 um)
    // that sits below the skin depth over most of the band.
    const drawT = (hi) => rng.bool(0.15) ? rng.logf(0.5e-6, 10e-6) : rng.logf(10e-6, hi);
    // Plating σ is drawn INDEPENDENTLY of the bulk σ over the same 1e6–6e7 range,
    // so plating lands both better- and worse-conducting than the bulk (silver
    // over copper vs nickel over copper); thickness spans thinner and thicker
    // than typical GHz skin depths, exercising the layered-impedance mix.
    const rollPlating = () => {
        if (!rng.bool(0.2)) return false;
        spec.plating_sigma = rng.logf(1e6, 6e7);
        // One draw in five puts the plating at or beyond the trace thickness: solid
        // plating metal on both backends (see S4 in docs/untested_edge_cases.md).
        const tAbs = Math.abs(spec.t ?? spec.bs_t);
        spec.plating_t = rng.bool(0.2) ? tAbs * rng.f(1, 1.5) : rng.logf(1e-6, 8e-6);
        spec.plating_thick_corners = rng.bool(0.7);
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
        const bs_w = rng.logf(0.05e-3, 1e-3), bs_t = drawT(50e-6);
        Object.assign(spec, {
            bs_w, bs_t,
            bs_x_offset: rng.bool(0.5) ? rng.f(-0.15e-3, 0.15e-3) : 0,
            bs_h_bottom: rng.logf(0.05e-3, 0.5e-3), bs_er_bottom: drawEr(6), bs_tand_bottom: drawTand(),
            // Floor of 2.5t: the middle layer must clear 2× trace_thickness or the
            // broadside pair collides and validation rejects the case on both backends.
            bs_h_middle: rng.logf(Math.max(0.05e-3, 2.5 * bs_t), 0.5e-3), bs_er_middle: drawEr(6), bs_tand_middle: drawTand(),
            bs_h_top: rng.logf(0.05e-3, 0.5e-3), bs_er_top: drawEr(6), bs_tand_top: drawTand(),
            bs_sigma: drawSigma(),
        });
        // Embedded (negative-thickness) traces in ~15% of draws: they grow into the
        // outer layers, which validation requires to be thicker than the trace.
        if (rng.bool(0.15)) spec.bs_t = -Math.min(bs_t, 0.8 * Math.min(spec.bs_h_bottom, spec.bs_h_top));
        spec.use_enclosure = rng.bool(0.3);
        // An enclosure always closes its side walls (see the enclosure block below).
        spec.use_side_gnd = spec.use_enclosure;
        // Sized from the offset pair's span: the offset upper trace must stay inside
        // the walls or validation rejects the case.
        if (spec.use_enclosure) spec.enclosure_width = (spec.bs_w + 2 * Math.abs(spec.bs_x_offset)) * rng.f(3, 8);
        spec.use_plating = rollPlating();
        return spec;
    }

    // Bulk conductor σ spans 1e6–6e7 S/m (poor alloys up to copper): at fixed
    // frequency this sweeps δ ~8×, walking the skin-transition/thickness ratio.
    const w = rng.logf(0.05e-3, 3e-3), h = rng.logf(0.05e-3, 1.5e-3);
    // Embedded (negative-thickness) traces in ~15% of draws, kept above the
    // substrate floor (validation requires t > -h).
    let t = drawT(70e-6);
    if (rng.bool(0.15)) t = -Math.min(t, 0.8 * h);
    Object.assign(spec, {
        w, h, t, er: drawEr(10), tand: drawTand(), sigma: drawSigma(),
        gnd_thickness: rng.logf(10e-6, 35e-6),
    });
    const isDiff = tl.includes('diff');
    if (isDiff) spec.trace_spacing = rng.logf(0.02e-3, 0.6e-3);
    // Top height floor of 2t: a cover thinner than the trace always fails
    // parameter validation (trace_thickness > enclosure_height) on both backends.
    // The cover material matches the substrate half the time (the classic
    // homogeneous stripline); otherwise it is air (eps_r 1) or an independent
    // dielectric, which floats the two grounds against each other in the full-wave
    // eigenproblem (see TriBackend._eigenAbc).
    if (tl.includes('stripline')) {
        if (rng.bool(0.5)) { spec.er_top = spec.er; spec.tand_top = spec.tand; }
        else if (rng.bool(0.4)) { spec.er_top = 1; spec.tand_top = 0; }
        else { spec.er_top = drawEr(10); spec.tand_top = drawTand(); }
        spec.stripline_top_h = rng.logf(Math.max(0.05e-3, 2 * Math.abs(t)), 1.5e-3);
    }
    if (tl.includes('gcpw')) {
        spec.gap = rng.logf(0.05e-3, 0.5e-3); spec.via_gap = rng.logf(0.05e-3, 0.5e-3);
        // Finite coplanar grounds in ~40% of draws, via fence at least 50 um wide
        if (rng.bool(0.4)) spec.gnd_width = 2 * spec.via_gap + rng.logf(0.05e-3, 3e-3);
    }

    // Advanced features. Solder mask applies to the microstrip AND gcpw families
    // (gcpw takes the coplanar solder mask path in MicrostripSolver); its material
    // and thicknesses are randomized — er well past the usual 3.5 to make the
    // C/eps shift it causes clearly visible in the comparison. Top dielectric and
    // ground cutout apply to the microstrip AND gcpw families, as in the app.
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
    if (tl === 'microstrip' || tl === 'diff_microstrip' || tl === 'gcpw' || tl === 'diff_gcpw') {
        // Top dielectric material drawn from the same range as the substrate.
        if (rng.bool(0.3)) {
            spec.use_top_diel = true; spec.top_diel_h = rng.logf(0.05e-3, 0.5e-3);
            spec.top_diel_er = rng.f(2.2, 10); spec.top_diel_tand = rng.f(0.001, 0.02);
        }
        if (rng.bool(0.2)) { spec.use_gnd_cut = true; spec.gnd_cut_w = w * rng.f(0.2, 1.5); spec.gnd_cut_h = h * rng.f(0.1, 0.7); }
    }
    if (tl === 'stripline' || tl === 'diff_stripline') {
        // Top dielectric layer carved out of the stripline cover (as in the app):
        // the ground-to-ground spacing is unchanged, the layer takes the lowest part.
        // Half the draws are air (eps_r = 1), the stripline-with-air-gap use case.
        if (rng.bool(0.3)) {
            spec.use_top_diel = true; spec.top_diel_h = spec.stripline_top_h * rng.f(0.2, 0.8);
            spec.top_diel_er = rng.bool(0.5) ? 1 : rng.f(2.2, 10);
            spec.top_diel_tand = spec.top_diel_er === 1 ? 0 : rng.f(0.001, 0.02);
        }
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
    // Besides the closed box, the app's other enclosure shapes: a lid alone over
    // the auto-sized (far) open sides, side walls alone under the auto-sized (far)
    // open top, and side walls on a stripline (its lid is intrinsic). A lid without
    // side walls leaves the grounds floating in the full-wave eigenproblem.
    if (rng.bool(0.35)) {
        spec.use_enclosure = true;
        let span = isDiff ? (2 * w + spec.trace_spacing) : w;
        // GCPW: the enclosure must clear the full coplanar active width (trace +
        // gaps + via fences, matching MicrostripSolver's active_width), not just
        // the trace span — otherwise the combo always fails validation.
        if (tl.includes('gcpw')) span += (isDiff ? 0 : w) + 2 * (spec.gap + (spec.gnd_width ?? spec.via_gap));
        const style = tl.includes('stripline') ? 'sides' : rng.pick(['box', 'box', 'lid', 'sides']);
        spec.enclosure_style = style;
        spec.use_side_gnd = style !== 'lid';
        spec.use_top_gnd = style !== 'sides';
        spec.enclosure_width = style === 'lid' ? 'auto' : span * rng.f(2.5, 6) + 2 * spec.gnd_thickness;
        spec.enclosure_height = style === 'sides' ? 'auto' : (h + t) * rng.f(1.5, 4);
    }
    spec.use_plating = rollPlating();
    return spec;
}

// --- Random custom geometry (geometry text, rectangles only: both backends solve it) ---
//
// A stack of full-width dielectric layers with the signal (or a pair) on one of its
// interfaces, embedded or on the surface, then:
//   grounds    - bottom: a gnd wall, a full-width plane over an open bottom, a finite
//                plane or none; coplanar grounds beside the signals (to the domain edge
//                or finite, one or both sides); a ground buried at another interface;
//                a top ground over the stack (wall, plane or finite). Finite grounds are
//                meshed metal in the full-wave loss solve, the ones reaching the domain
//                edges on three sides are wall grounds. Separate ground bodies float
//                against each other in the full-wave eigenproblem.
//   pairs      - equal or unequal widths and thicknesses, the two traces on the same
//                interface or on different ones (laterally offset)
//   dielectric - random rectangles painted over the stack (later ones win), some of
//                them over the trace; in ~12% of cases one dielectric conducts
//   metal      - stacked touching blocks for a signal (T shapes), per-conductor
//                conductivity, roughness and plating material on any conductor
//   bounds     - open sides, or ground side walls with an explicit domain width
// Sizes are drawn in um, the text is written in um.
export function randomCustomSpec(rng) {
    const freq = rng.bool(0.3) ? rng.logf(1e6, 0.5e9) : rng.logf(0.5e9, 6e9);
    const spec = { tl: 'custom', freq, sigma: rng.logf(5e5, 6.3e7), rq: rng.bool(0.3) ? rng.logf(0.1e-6, 3e-6) : 0,
        thick_plating: rng.bool(0.7), use_causal: rng.bool(0.2), tags: [] };
    const tags = spec.tags;
    const n = v => +v.toPrecision(5);
    const drawEr = () => (rng.bool(0.1) ? rng.f(1, 2.2) : rng.f(2.2, 10));
    const drawTand = () => (rng.bool(0.1) ? rng.logf(0.02, 0.2) : rng.f(0.001, 0.02));
    const drawT = () => (rng.bool(0.15) ? rng.logf(0.5, 10) : rng.logf(10, 70));

    // Dielectric layers from y = 0 up, interfaces Y[0..nL].
    const nL = 1 + Math.floor(rng.f(0, 4));
    const H = Array.from({ length: nL }, () => rng.logf(40, 600));
    const Y = [0];
    for (const h of H) Y.push(Y[Y.length - 1] + h);
    const diels = H.map((h, k) => ({ x: '-inf', w: 'inf', y: Y[k], h }));
    // Signal interface: on top of layer ks (ks = nL is the surface).
    const ks = 1 + Math.floor(rng.f(0, nL));
    // A trace inside the stack must fit in the layer it grows into.
    const fit = (k, t) => (k < nL ? Math.min(t, 0.8 * H[k]) : t);
    // A cover layer over a surface trace.
    if (ks === nL && rng.bool(0.3)) {
        const h = rng.logf(30, 500);
        diels.push({ x: '-inf', w: 'inf', y: Y[nL], h });
        H.push(h); Y.push(Y[nL] + h);
        tags.push('cover');
    }
    const nTop = Y.length - 1;

    // Signals: [{ kind, x0, w, y, t }].
    const sigs = [];
    const w1 = rng.logf(50, 2000), t1 = fit(ks, drawT());
    const pair = rng.bool(0.5);
    if (!pair) sigs.push({ kind: 'sig+', x0: -w1 / 2, w: w1, y: Y[ks], t: t1, k: ks });
    else {
        const s = rng.logf(20, 600);
        const w2 = rng.bool(0.6) ? w1 : w1 * rng.f(0.5, 2);
        let k2 = ks;
        if (nTop >= 2 && rng.bool(0.2)) {
            // The other trace on another interface, laterally offset.
            const ks2 = 1 + Math.floor(rng.f(0, nTop));
            if (ks2 !== ks) k2 = ks2;
        }
        const t2 = fit(k2, rng.bool(0.8) ? t1 : drawT());
        const [plus, minus] = rng.bool(0.5) ? ['sig+', 'sig-'] : ['sig-', 'sig+'];
        if (k2 === ks) {
            sigs.push({ kind: minus, x0: -s / 2 - w1, w: w1, y: Y[ks], t: t1, k: ks });
            sigs.push({ kind: plus, x0: s / 2, w: w2, y: Y[ks], t: t2, k: ks });
        } else {
            const lo = Math.min(ks, k2), hi = Math.max(ks, k2);
            // The lower trace must clear the upper interface.
            const tLo = lo === ks ? t1 : t2;
            if (Y[lo] + tLo >= Y[hi] - 5) return randomCustomSpec(rng);
            const off = rng.f(-1, 1) * w1;
            sigs.push({ kind: minus, x0: -w1 / 2, w: w1, y: Y[ks], t: t1, k: ks });
            sigs.push({ kind: plus, x0: -w2 / 2 + off, w: w2, y: Y[k2], t: t2, k: k2 });
            tags.push('stacked-pair');
        }
        if (w2 !== w1 || t2 !== t1 || k2 !== ks) tags.push('asym-pair');
        else tags.push('pair');
    }
    const sx0 = Math.min(...sigs.map(g => g.x0)), sx1 = Math.max(...sigs.map(g => g.x0 + g.w));
    const span = sx1 - sx0;
    const sigK = new Set(sigs.map(g => g.k));

    const conds = [];   // { kind, x, w, y, h, noPlating }
    // Signal metal, a T-shaped stack of two touching blocks in some cases.
    for (const g of sigs) {
        if (g.k === ks && !tags.includes('stacked-pair') && rng.bool(0.15)) {
            const t0 = g.t * rng.f(0.3, 0.7), wTop = g.w * rng.f(0.3, 0.9);
            // The upper block starts at the same expression as the lower one ends, so the
            // two touch exactly after the text rounds each number.
            conds.push({ kind: g.kind, x: g.x0, w: g.w, y: g.y, h: t0, noPlating: true });
            conds.push({ kind: g.kind, x: g.x0 + (g.w - wTop) / 2, w: wTop, y: `${+g.y.toPrecision(5)}+${+t0.toPrecision(5)}`,
                h: g.t - t0, noPlating: true });
            if (!tags.includes('blocks')) tags.push('blocks');
        } else conds.push({ kind: g.kind, x: g.x0, w: g.w, y: g.y, h: g.t });
    }

    // Grounds.
    const bounds = ['open', 'open', 'open', 'open'];
    let domain = null;
    const tg = () => rng.logf(5, 50);
    const gnd = (x, w, y, h, tag) => { conds.push({ kind: 'gnd', x, w, y, h }); tags.push(tag); };
    const bottom = rng.pick(['wall', 'wall', 'wall', 'plane', 'plane', 'finite', 'finite', 'none', 'none']);
    if (bottom === 'wall') { bounds[3] = 'gnd'; tags.push('gnd:wall'); }
    else if (bottom === 'plane') { const t = tg(); gnd('-inf', 'inf', -t, t, 'gnd:plane'); }
    else if (bottom === 'finite') { const t = tg(), W = span * rng.f(1.2, 8); gnd(-W / 2 + (sx0 + sx1) / 2, W, -t, t, 'gnd:finite'); }
    // Coplanar grounds at the signal interface.
    if (rng.bool(0.3)) {
        const t = sigs[0].t, y = Y[ks];
        const gL = rng.logf(20, 500), gR = rng.bool(0.8) ? gL : rng.logf(20, 500);
        const sides = rng.pick(['both', 'both', 'both', 'left', 'right']);
        const toEdge = rng.bool(0.5);
        const wL = rng.logf(50, 2000), wR = rng.logf(50, 2000);
        const lo = Math.min(...sigs.filter(g => g.k === ks).map(g => g.x0));
        const hi = Math.max(...sigs.filter(g => g.k === ks).map(g => g.x0 + g.w));
        if (sides !== 'right') gnd(lo - gL, toEdge ? '-inf' : -wL, y, t, 'gnd:coplanar');
        if (sides !== 'left') gnd(hi + gR, toEdge ? 'inf' : wR, y, t, 'gnd:coplanar');
    }
    // A ground buried at an interface without a signal (a ground plane inside the stack).
    const free = [];
    for (let k = 1; k < nTop; k++) if (!sigK.has(k)) free.push(k);
    if (free.length && rng.bool(0.15)) {
        const kb = rng.pick(free), t = tg();
        const below = [...sigK].some(k => k > kb);
        // Below a signal it hangs under the interface, above one it sits on it.
        const y = below ? Y[kb] - Math.min(t, 0.5 * H[kb - 1]) : Y[kb];
        const h = below ? Y[kb] - y : Math.min(t, 0.5 * H[kb]);
        const clear = sigs.every(g => g.y + g.t < y - 1 || g.y > y + h + 1);
        if (clear) {
            if (rng.bool(0.6)) { const W = span * rng.f(1, 6); gnd(-W / 2 + (sx0 + sx1) / 2, W, y, h, 'gnd:buried'); }
            else gnd('-inf', 'inf', y, h, 'gnd:buried');
        }
    }
    // A ground over the stack when something covers the signals.
    const topSig = Math.max(...sigs.map(g => g.y + g.t));
    if (Y[nTop] > topSig + 1 && rng.bool(0.4)) {
        const top = rng.pick(['wall', 'plane', 'finite']);
        if (top === 'wall') { bounds[2] = 'gnd'; tags.push('gnd:top-wall'); }
        else if (top === 'plane') { const t = tg(); gnd('-inf', 'inf', Y[nTop], t, 'gnd:top-plane'); }
        else { const t = tg(), W = span * rng.f(1.2, 8); gnd(-W / 2 + (sx0 + sx1) / 2, W, Y[nTop], t, 'gnd:top-finite'); }
    }
    if (!conds.some(c => c.kind === 'gnd') && !bounds.includes('gnd')) { bounds[3] = 'gnd'; tags.push('gnd:wall'); }
    // Ground side walls at an explicit domain width.
    if (rng.bool(0.15)) {
        bounds[0] = bounds[1] = 'gnd';
        // The walls clear every finite conductor.
        const cx = (sx0 + sx1) / 2;
        let reach = span / 2;
        for (const c of conds) {
            if (typeof c.x !== 'string') reach = Math.max(reach, Math.abs(c.x - cx));
            if (typeof c.x !== 'string' && typeof c.w !== 'string') reach = Math.max(reach, Math.abs(c.x + c.w - cx));
        }
        const half = reach + span * rng.f(0.5, 3) + 100;
        domain = [n(cx - half), n(cx + half)];
        tags.push('side-walls');
    }

    // Random dielectric rectangles over the stack, some of them over the trace.
    const nR = rng.bool(0.4) ? 0 : 1 + Math.floor(rng.f(0, 3));
    for (let r = 0; r < nR; r++) {
        if (rng.bool(0.35)) {
            const c = rng.logf(5, 0.5 * span + 10);
            diels.push({ x: sx0 - c, w: span + 2 * c, y: Y[ks], h: topSig - Y[ks] + rng.logf(5, 200), extra: true });
            tags.push('diel:over');
        } else {
            const w = span * rng.logf(0.2, 3), x = (sx0 + sx1) / 2 + rng.f(-1, 1) * span - w / 2;
            const y = Y[ks] + rng.f(-1, 1) * Math.max(50, 0.5 * Y[nTop]);
            diels.push({ x, w, y, h: rng.logf(10, 400), extra: true });
            tags.push('diel:rect');
        }
    }
    diels.forEach(d => { d.er = drawEr(); d.tand = drawTand(); });
    // One conducting dielectric in some cases (a doped silicon substrate or layer).
    if (rng.bool(0.12)) {
        rng.pick(diels).sigma = rng.logf(0.01, 300);
        tags.push('cond-diel');
    }

    // Per-conductor metal and finish. Plating only on conductors that touch no other of
    // their kind: single-block signals, and finite grounds clear of the side walls.
    let own = false;
    const sideWalls = bounds[0] === 'gnd';
    for (const c of conds) {
        if (c.kind === 'gnd' && (sideWalls || typeof c.w === 'string')) c.noPlating = true;
        const keys = [];
        if (rng.bool(0.25)) keys.push(`sigma=${n(rng.logf(5e5, 6.3e7))}`);
        if (rng.bool(0.25)) keys.push(`rq=${n(rng.logf(0.1, 3))}`);
        if (!c.noPlating && rng.bool(0.25)) {
            const faces = ['top', 'sides', 'bottom'].filter(() => rng.bool(0.6));
            keys.push(`plating=${faces.length ? faces.join(',') : 'top'}`,
                `plating_sigma=${n(rng.logf(1e6, 6e7))}`, `plating_t=${n(rng.logf(1, 8))}`);
            if (rng.bool(0.3)) keys.push(`plating_rq=${n(rng.logf(0.1, 2))}`);
            if (rng.bool(0.15)) keys.push(`plating_rq_iface=${n(rng.logf(0.1, 2))}`);
            if (!tags.includes('plating')) tags.push('plating');
        }
        if (keys.length) own = true;
        c.keys = keys;
    }
    if (own) tags.push('own-metal');

    const v = x => (typeof x === 'string' ? x : String(n(x)));
    const lines = ['units um', `bounds ${bounds.join(' ')}`];
    if (domain) lines.push(`domain ${domain[0]} ${domain[1]} auto auto`);
    for (const d of diels) {
        lines.push(`diel x=${v(d.x)} w=${v(d.w)} y=${v(d.y)} h=${v(d.h)} er=${n(d.er)} tand=${n(d.tand)}` +
            (d.sigma ? ` sigma=${n(d.sigma)}` : ''));
    }
    for (const c of conds) lines.push(`${c.kind} x=${v(c.x)} w=${v(c.w)} y=${v(c.y)} h=${v(c.h)}${c.keys.length ? ' ' + c.keys.join(' ') : ''}`);
    spec.text = lines.join('\n');
    spec.tags = [...new Set(tags)];
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
        top: spec.plating_top ?? true, sides: spec.plating_sides ?? true, bottom: spec.plating_bottom ?? false,
        thick_corners: spec.plating_thick_corners ?? true };
}

export function buildSolver(spec, backend) {
    if (spec.tl === 'custom') {
        // The app's custom-geometry options (solver_factory.buildSolverFromParams).
        const s = new CustomGeometrySolver({ text: spec.text, sigma_cond: spec.sigma, freq: spec.freq, rq: spec.rq,
            thick_plating: spec.thick_plating, mesh_backend: backend, nx: 30, ny: 30 });
        s.use_causal_materials = !!spec.use_causal;
        return s;
    }
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
        top: spec.plating_top ?? true, sides: spec.plating_sides ?? true, bottom: spec.plating_bottom ?? false,
        thick_corners: spec.plating_thick_corners ?? true };
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
    if (spec.tl.includes('gcpw')) {
        o.boundaries = ['open', 'open', 'open', 'gnd']; o.use_coplanar_gnd = true; o.gap = spec.gap; o.via_gap = spec.via_gap; o.use_vias = true;
        if (spec.gnd_width) o.coplanar_gnd_width = spec.gnd_width;
    }
    else if (spec.tl.includes('stripline')) {
        o.epsilon_r_top = spec.er_top; o.tan_delta_top = spec.tand_top; o.boundaries = ['open', 'open', 'gnd', 'gnd'];
        o.enclosure_height = spec.stripline_top_h - (spec.use_top_diel ? spec.top_diel_h : 0);
    }
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
        // solveAt → _applyMaterials), keeping the cross-backend comparison fair.
        if (spec.use_causal) r = await s.computeAtFrequency(s.freq, r);
        // Compare every mode the backend returns: single-ended → [single]; differential →
        // [odd, even], in that fixed order from the shared FieldSolver2D path on both backends.
        return {
            modes: r.modes.map(m => ({ Z0: m.Z0, eps_eff: m.eps_eff, C: m.RLGC.C, R: m.RLGC.R, G: m.RLGC.G, L: m.RLGC.L, mode: m.mode })),
            // Warning tags for the R-gate relaxation: QS loss-accuracy reasons
            // ('skin-transition' / 'narrow-gap') and tri fallback types.
            warns: (r.warnings || []).map(w => ({ type: w.type, reason: w.reason || null, mode: w.mode ?? null })),
            maxQ: s.meshQuality ? s.meshQuality.maxQ : null,
            degenerate: s.meshQuality ? s.meshQuality.degenerateCount : 0,
            nan: s.meshQuality ? s.meshQuality.nanNodes : 0,
            // Structural quality signals for the bad-mesh gate (see BAD_Q note).
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
    if (spec.tl === 'custom') {
        // The geometry text on one line (';' separates statements), pasteable into the app.
        const f = spec.freq < 0.1e9 ? (spec.freq / 1e6).toFixed(2) + 'MHz' : (spec.freq / 1e9).toFixed(2) + 'GHz';
        const opts = [spec.rq ? `rq=${u(spec.rq, 1e6, 2)}µm` : '', spec.thick_plating ? '' : 'thin-plating',
            spec.use_causal ? 'causal' : ''].filter(Boolean);
        return `custom sig=${(spec.sigma / 1e6).toFixed(1)}e6 f=${f}${opts.length ? ' ' + opts.join(' ') : ''} ` +
            `[${spec.tags.join(',')}]\n      ${spec.text.split('\n').join('; ')}`;
    }
    let s = spec.tl;
    if (spec.tl === 'broadside_stripline') s += ` w=${u(spec.bs_w)}mm t=${u(spec.bs_t, 1e6, 1)}µm off=${u(spec.bs_x_offset)}mm hM=${u(spec.bs_h_middle)}mm`;
    else s += ` w=${u(spec.w)}mm h=${u(spec.h)}mm t=${u(spec.t, 1e6, 1)}µm er=${spec.er.toFixed(2)}`;
    s += ` sig=${((spec.sigma ?? spec.bs_sigma) / 1e6).toFixed(1)}e6 f=${spec.freq < 0.1e9 ? (spec.freq / 1e6).toFixed(2) + 'MHz' : (spec.freq / 1e9).toFixed(2) + 'GHz'}`;
    if (spec.trace_spacing) s += ` gap=${u(spec.trace_spacing)}mm`;
    if (spec.gnd_width) s += ` gnd_w=${u(spec.gnd_width)}mm`;
    if (spec.tl.includes('stripline') && spec.tl !== 'broadside_stripline' && spec.er_top !== spec.er) s += ` er_top=${spec.er_top.toFixed(2)}`;
    const feats = ['use_sm', 'use_top_diel', 'use_gnd_cut', 'use_enclosure', 'use_side_gnd', 'use_top_gnd', 'use_plating', 'use_causal']
        .filter(k => spec[k]).map(k => {
            if (k === 'use_sm' && spec.sm_style && spec.sm_style !== 'full') return `sm:${spec.sm_style}`;
            if (k === 'use_enclosure' && spec.enclosure_style) return `enclosure:${spec.enclosure_style}`;
            if (k === 'use_plating') {
                const faces = [spec.plating_top && 't', spec.plating_sides && 's', spec.plating_bottom && 'b']
                    .filter(Boolean).join('');
                const tAbs = Math.abs(spec.t ?? spec.bs_t);
                return `plating:${faces || 'ts'}${spec.plating_rq ? '+rq' : ''}${spec.plating_thick_corners === false ? '+flat' : ''}${spec.plating_t >= tAbs ? '+solid' : ''}`;
            }
            return k.replace('use_', '');
        });
    if (feats.length) s += ` [${feats.join(',')}]`;
    return s;
}

const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);

export function printHeader(cfg) {
    console.log(`\n### Fuzz QS vs full-wave — N=${cfg.N} seed=${cfg.seed} flag>${(cfg.thresh * 100).toFixed(0)}% ` +
        `R>${(cfg.rThresh * 100).toFixed(0)}% (warned ${(R_RELAX * 100).toFixed(0)}%) Qmax>${BAD_Q} ###`);
    console.log(`covers: ${TL_TYPES.join(', ')}, custom + solder-mask/top-diel/gnd-cut/enclosure/plating/roughness/causal; diff pairs check odd+even\n`);
}

// The specs of a run, in case order: every case is drawn (so the stream does not
// depend on ONLY), and the ones to solve are returned with their index. A second
// stream replaces about one case in five by a custom geometry, so the native cases
// keep the specs of the native stream at their index. FAMILY=custom or native draws
// one family only.
const CUSTOM_SHARE = 0.2;
export function drawCases(cfg) {
    const rng = createRng(cfg.seed);
    const crng = createRng((cfg.seed ^ 0x5bd1e995) >>> 0);
    const family = process.env.FAMILY || 'all';
    const cases = [];
    for (let i = 0; i < cfg.N; i++) {
        let spec = randomSpec(rng);
        const custom = crng.bool(CUSTOM_SHARE);
        const cspec = randomCustomSpec(crng);
        if (family === 'custom' || (family === 'all' && custom)) spec = cspec;
        if (!cfg.only || cfg.only.has(i)) cases.push({ i, spec });
    }
    return cases;
}

// Solve one case on both backends and compare. Returns a plain (cloneable) record:
// the report lines and what was flagged, for printSummary.
export async function evaluateCase(i, spec, cfg) {
    const res = { i, tl: spec.tl, desc: fmtSpec(spec), lines: [], status: 'compared', ok: false,
        badMesh: false, discrepancy: null, loss: null, rlgc: null, floating: null, crash: null };
    const log = line => res.lines.push(line);
    const THRESH = cfg.thresh, R_THRESH = cfg.rThresh;
    const qs = await solveOn(spec, 'rectilinear');
    const fw = await solveOn(spec, 'triangular');

    // Both reject → invalid random geometry (validation is backend-independent); skip.
    if (qs.error && fw.error) { res.status = 'skipped'; return res; }
    // Only one backend fails → a real backend bug — UNLESS it's the pre-mesh
    // meshability guard cleanly refusing a geometry one backend's cost model
    // can't afford (working as intended; the guard has its own dedicated test).
    if (qs.error || fw.error) {
        const err = String(qs.error || fw.error);
        const who = qs.error ? `quasi-static threw (full-wave OK): ${qs.error}` : `full-wave threw (quasi-static OK): ${fw.error}`;
        if (/cannot be meshed/i.test(err)) {
            res.status = 'rejected';
            log(`[${i}] ○ MESHABILITY REJECTED (${qs.error ? 'quasi-static' : 'full-wave'} guard)\n      ${fmtSpec(spec)}`);
        } else {
            res.status = 'crash'; res.crash = who;
            log(`[${i}] ✗ ONE-SIDED FAILURE — ${who}\n      ${fmtSpec(spec)}`);
        }
        return res;
    }

    // areaRatio is deliberately NOT gated: adaptive skin-band refinement
    // routinely produces max/min element-area ratios of 1e6-1e7 (µm-scale
    // skin elements vs mm-scale far field) on perfectly healthy meshes —
    // measured on cases with Qmax 3-6 and badFraction ~0 — and the direct
    // factorization handles it; the ratio separates nothing here.
    const badMesh = (fw.maxQ != null && fw.maxQ > BAD_Q) || fw.degenerate > 0 || fw.nan > 0
        || fw.badFraction > 0.05;
    if (badMesh) {
        res.badMesh = true;
        log(`[${i}] ✗ BAD MESH (Qmax=${fw.maxQ?.toFixed(0)} degen=${fw.degenerate} nan=${fw.nan} ` +
            `areaRatio=${fw.areaRatio?.toExponential(1)} ` +
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
    // A mode the full-wave flags 'floating-grounds' (grounds not connected, fill
    // differing between them) couples to the parallel-plate wave: a different
    // structure from the stitched-ground quasi-static line, with an eps_eff that
    // depends on the domain width (or none at all where the mode leaks). So does an
    // 'eigen-ambiguous' pick (the mode split across near-degenerate eigenpairs, e.g.
    // a box mode near its cutoff), and 'wall-tied-grounds' (grounds connected only
    // through enclosure side walls: the same parallel-plate wave, as box modes). Their
    // Z0/C are not compared, only counted.
    const nModes = Math.min(qs.modes.length, fw.modes.length);
    const SKIP_TYPES = ['floating-grounds', 'eigen-ambiguous', 'wall-tied-grounds'];
    const floatingMode = m => fw.warns.some(w => SKIP_TYPES.includes(w.type) && (w.mode === m || w.mode === 'all'));
    let worst = null, floating = null;
    for (let mi = 0; mi < nModes; mi++) {
        const q = qs.modes[mi], w = fw.modes[mi];
        const dZ = relDiff(q.Z0, w.Z0), dE = relDiff(q.eps_eff, w.eps_eff), dC = relDiff(q.C, w.C);
        const dMax = Math.max(dZ, dC);
        if (floatingMode(q.mode)) {
            if (!floating || dMax > floating.dMax) floating = { mode: q.mode, dMax };
            continue;
        }
        if (!worst || dMax > worst.dMax) worst = { mode: q.mode, dZ, dE, dC, dMax, q, w };
    }
    res.floating = floating;
    // Conductor-loss gate (R per mode). Relaxed when either backend flagged
    // reduced loss accuracy for this solve: QS reasons 'skin-transition' /
    // 'narrow-gap', or any tri mqs-* warning — the backend downgrading
    // MQS to perturbation, or 'mqs-band-capped' (MQS on an under-resolved skin
    // band). Certificate warnings describe C, not R, and do NOT relax the gate.
    const lossRelaxed =
        qs.warns.some(w => w.reason === 'skin-transition' || w.reason === 'narrow-gap'
            || w.reason === 'plating-transition')
        || fw.warns.some(w => /^mqs-/.test(w.type) || w.reason === 'plating-transition');
    // Shunt conductance and total inductance, each at its own gate. G follows the
    // same static loss weighting on both backends, so it tracks C up to the mesh
    // (and the lossless-air-void class Q1 that C cannot see). L carries the
    // backends' different internal-inductance models on top of the shared vacuum
    // C0, a few percent of L at these frequencies. Under a loss-accuracy note the
    // L gate relaxes like the R gate: with the skin depth at or beyond the conductor
    // thickness the quasi-static surface-field L_internal has no lateral current
    // spreading and reads 10-30% of L away from the MQS volume solve, so only gross
    // errors are detectable there.
    const lGate = lossRelaxed ? L_RELAX : L_THRESH;
    let worstGL = null;
    for (let mi = 0; mi < nModes; mi++) {
        const q = qs.modes[mi], w = fw.modes[mi];
        const dG = (q.G > 0 && w.G > 0) ? relDiff(q.G, w.G) : 0;
        const dL = relDiff(q.L, w.L);
        const over = Math.max(dG / G_THRESH, dL / lGate);
        if (!worstGL || over > worstGL.over) worstGL = { mode: q.mode, dG, dL, over, q, w };
    }
    if (worstGL && worstGL.over > 1) {
        res.rlgc = { mode: worstGL.mode, dG: worstGL.dG, dL: worstGL.dL };
        const tag = nModes > 1 ? ` [${worstGL.mode}]` : '';
        log(`[${i}] ✗ G/L DISCREPANCY${tag} (G ${(worstGL.dG * 100).toFixed(0)}% L ${(worstGL.dL * 100).toFixed(0)}%)` +
            `${lossRelaxed ? ' (relaxed L gate)' : ''}\n` +
            `      ${fmtSpec(spec)}\n` +
            `      qs: G=${worstGL.q.G.toExponential(3)} L=${(worstGL.q.L * 1e9).toFixed(2)}nH | ` +
            `fw: G=${worstGL.w.G.toExponential(3)} L=${(worstGL.w.L * 1e9).toFixed(2)}nH`);
    }
    if (worst && worst.dMax > THRESH) {
        res.discrepancy = { mode: worst.mode, dZ: worst.dZ, dE: worst.dE, dC: worst.dC };
        const tag = nModes > 1 ? ` [${worst.mode}]` : '';
        log(`[${i}] ✗ DISCREPANCY ${(worst.dMax * 100).toFixed(0)}%${tag} (Z0 ${(worst.dZ * 100).toFixed(0)}% eps ${(worst.dE * 100).toFixed(0)}% C ${(worst.dC * 100).toFixed(0)}%)\n` +
            `      ${fmtSpec(spec)}\n` +
            `      qs: Z0=${worst.q.Z0.toFixed(2)} eps=${worst.q.eps_eff.toFixed(3)} C=${(worst.q.C * 1e12).toFixed(1)}pF | ` +
            `fw: Z0=${worst.w.Z0.toFixed(2)} eps=${worst.w.eps_eff.toFixed(3)} C=${(worst.w.C * 1e12).toFixed(1)}pF`);
    }

    const rGate = lossRelaxed ? R_RELAX : R_THRESH;
    let worstR = null;
    for (let mi = 0; mi < nModes; mi++) {
        const q = qs.modes[mi], w = fw.modes[mi];
        if (!(q.R > 0) || !(w.R > 0)) continue;
        const dR = relDiff(q.R, w.R);
        if (!worstR || dR > worstR.dR) worstR = { mode: q.mode, dR, q, w };
    }
    if (worstR && worstR.dR > rGate) {
        res.loss = { mode: worstR.mode, dR: worstR.dR, relaxed: lossRelaxed };
        const tag = nModes > 1 ? ` [${worstR.mode}]` : '';
        log(`[${i}] ✗ LOSS DISCREPANCY ${(worstR.dR * 100).toFixed(0)}%${tag}` +
            `${lossRelaxed ? ' (relaxed gate)' : ''}\n` +
            `      ${fmtSpec(spec)}\n` +
            `      R qs=${worstR.q.R.toFixed(2)} vs fw=${worstR.w.R.toFixed(2)} Ω/m`);
    }
    res.ok = !badMesh && (!worst || worst.dMax <= THRESH) && (!worstR || worstR.dR <= rGate)
        && (!worstGL || worstGL.over <= 1);
    return res;
}

// Summary of a run from the evaluateCase records (any order). Sets the exit code
// when anything was flagged.
export function printSummary(results, cfg) {
    const sorted = [...results].sort((a, b) => a.i - b.i);
    const byType = {};
    for (const r of sorted) byType[r.tl] = (byType[r.tl] || 0) + 1;
    const pick = k => sorted.filter(r => r[k]);
    const count = st => sorted.filter(r => r.status === st).length;
    const ok = sorted.filter(r => r.ok).length;
    const loss = pick('loss'), floating = pick('floating');
    console.log(`\n${'='.repeat(72)}`);
    console.log(`coverage: ${Object.entries(byType).map(([k, v]) => `${k}:${v}`).join('  ')}`);
    // Split the loss failures by gate: an UNWARNED one is the interesting kind — it
    // extends the known bias envelope rather than re-confirming a regime both backends
    // already flag as approximate (see the R gate note above).
    const unwarned = loss.filter(r => !r.loss.relaxed);
    console.log(`SUMMARY: ${ok}/${cfg.N} clean | ${pick('discrepancy').length} discrepancy | ` +
        `${loss.length} loss-discrepancy (${unwarned.length} unwarned, ${loss.length - unwarned.length} warned) | ` +
        `${pick('rlgc').length} G/L-discrepancy | ${pick('badMesh').length} bad-mesh | ` +
        `${count('crash')} one-sided-fail | ${count('rejected')} meshability-rejected | ${count('skipped')} skipped(invalid) | ` +
        `${floating.length} floating/wall-tied-grounds/ambiguous (Z0/C not compared, max ` +
        `${(100 * Math.max(0, ...floating.map(r => r.floating.dMax))).toFixed(0)}%)`);
    if (unwarned.length) {
        console.log(`unwarned loss rows (investigate, do not just raise the gate):`);
        for (const r of unwarned)
            console.log(`  [${r.i}] ${(100 * r.loss.dR).toFixed(0)}%${r.loss.mode ? ` [${r.loss.mode}]` : ''}  ${r.desc}`);
    }
    console.log(`${'='.repeat(72)}`);
    if (sorted.some(r => r.discrepancy || r.loss || r.rlgc || r.badMesh || r.crash)) process.exitCode = 1;
}

async function main() {
    const cfg = parseArgs(process.argv.slice(2));
    printHeader(cfg);
    const results = [];
    for (const { i, spec } of drawCases(cfg)) {
        const r = await evaluateCase(i, spec, cfg);
        for (const line of r.lines) console.log(line);
        results.push(r);
    }
    printSummary(results, cfg);
}

import { pathToFileURL } from 'url';
if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) main();
