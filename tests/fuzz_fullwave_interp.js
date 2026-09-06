// Fuzzer: full-wave TYPICAL USE (all accuracy-tradeoff optimizations on, results
// read through an InterpolatingSweep — exactly the app's default path) vs a
// full-wave REFERENCE with every accuracy-tradeoff optimization disabled and
// every compared frequency solved exactly.
//
// The optimizations under test (all default-on in the triangular backend):
//   • ε_eff dispersion anchor cache      (dispTol,       off via 0)
//   • MQS R(f) anchor cache              (mqsInterpTol,  off via 0)
//   • 2·t thickness term + floor in hFine (reference forces the old 1.5·t sizing)
//   • skin-band economy sizing 0.6δ/1.5δ (reference uses the pre-tuning 0.5δ/2δ)
//   • InterpolatingSweep spline itself   (reference calls computeAtFrequency exactly)
// Both sides share the same adaptive-refinement settings (that loop is identical
// code, not an optimization difference) and the same random geometries as the
// QS-vs-fullwave fuzzer (imported generator — same seed ⇒ same geometry list).
//
// Compares the per-mode RLGC — the quantities the sweep actually splines and
// everything downstream (Zc, alpha, eps_eff, S-params) derives from. The
// DERIVED quantities can't be compared across the two paths directly: the
// sweep's buildResults computes Z0/eps_eff from RLGC via the RF formulas
// (gamma and Zc including R and G — significant at low f where R/omega·L is
// not small), while an exact solve reports the static Z0 and the phase
// convention eps_eff = c²·L·C. Identical RLGC ⇒ identical derived values
// within each convention.
//
// Flags a case when any RLGC component at any compare frequency exceeds its
// gate: L/C > 2%, R/G > 4%. The gates are the measured stack-up of the
// individually-validated tolerances (spline 0.5%, anchor caches 0.1–0.2%,
// sizing ≤0.8%, skin band ≤1% on R) with headroom; a genuine regression in
// any layer blows well past them.
//
// A second comparison isolates the anchor caches on the frequency-list path
// (solve_sweep): the adaptive solve at f_max seeds one anchor at the top and
// the list then comes ascending, an order the caches' leave-one-out gate must
// handle without the sweep's own midpoint checks behind it. Typical use is
// compared with the same solver and the caches off (every point exact) on the
// identical mesh, so its gates are the caches' own tolerances with a little
// headroom (see test_anchor_caches.js for the deterministic version).
//
// Usage:   node tests/fuzz_fullwave_interp.js [N] [seed]
//   e.g.   node tests/fuzz_fullwave_interp.js 12 1
// Reproduce a flagged case: rerun with the same N and seed; ONLY=4,11 solves
// just those case indices (the others are drawn and skipped).

import { createRng, randomSpec, buildSolver } from './fuzz_qs_vs_fullwave.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';

const N    = parseInt(process.argv[2]) || 12;
const SEED = process.argv[3] !== undefined ? parseInt(process.argv[3]) : 1;
const ONLY = process.env.ONLY ? new Set(process.env.ONLY.split(',').map(Number)) : null;

const GATES = { R: 0.04, L: 0.02, G: 0.04, C: 0.02 };
// Frequency-list path, caches on vs off on the same mesh: R and L_internal anchors
// are gated at 0.2%, eps_d at 0.1%.
const LIST_GATES = { R: 0.003, L: 0.003, G: 0.0015, C: 0.0015 };
const N_LIST = 9;   // interpolation needs five anchors, so the last four points exercise the gate
const ADAPTIVE = { max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };
const F_START = 0.1e9;

const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
const pickModes = r => r.modes.map(m => ({
    mode: m.mode, R: m.RLGC.R, L: m.RLGC.L, G: m.RLGC.G, C: m.RLGC.C,
}));

// Worst gate ratio over frequencies, modes and RLGC components.
function worstOf(typ, ref, freqs, gates) {
    let worst = null;
    for (let fi = 0; fi < freqs.length; fi++) {
        const nModes = Math.min(typ[fi].length, ref[fi].length);
        for (let mi = 0; mi < nModes; mi++)
            for (const q of Object.keys(gates)) {
                const a = ref[fi][mi][q], b = typ[fi][mi][q];
                if (!isFinite(a) || !isFinite(b)) continue;
                const d = relDiff(a, b);
                const over = d / gates[q];
                if (!worst || over > worst.over)
                    worst = { over, q, d, gate: gates[q], f: freqs[fi], mode: ref[fi][mi].mode, a, b };
            }
    }
    return worst;
}
const fmtWorst = w => `${w.q} ${(w.d * 100).toFixed(2)}%`;
const fmtMismatch = w => `${w.q} ${(w.d * 100).toFixed(2)}% (gate ${(w.gate * 100).toFixed(2)}%) @ ${(w.f / 1e9).toFixed(2)} GHz [${w.mode}] ref=${w.a.toExponential(4)} typ=${w.b.toExponential(4)}`;

// Typical use: adaptive solve + InterpolatingSweep, read results off the splines.
async function solveTypical(spec, fStop, compareFs) {
    const s = buildSolver(spec, 'triangular');
    s.tri_opts = { lossMethod: 'auto' };
    await s.solve_adaptive({ ...ADAPTIVE });
    const sweep = new InterpolatingSweep(s, new Map(), { tolerance: 0.005 });
    await sweep.run(F_START, fStop, {});
    return sweep.buildResults(compareFs).map(r => pickModes(r.result));
}

// Reference: anchor caches off, pre-optimization mesh sizings, exact solve per point.
async function solveReference(spec, fStop, compareFs) {
    const s = buildSolver(spec, 'triangular');
    const tAbs = Math.max(Math.abs(s.t ?? 35e-6), 1e-9);
    const wRef = Math.max(s.w ?? s.domain_width / 10, 1e-9);
    s.tri_opts = {
        lossMethod: 'auto',
        dispTol: 0, mqsInterpTol: 0,                  // anchor caches off → every point exact
        hFine: Math.min(tAbs, wRef / 4) * 1.5,        // pre-2026-07-25 thickness sizing
        mqsBandDelta: 0.5, mqsBand: 2,                // pre-tuning skin band (finer, wider)
    };
    await s.solve_adaptive({ ...ADAPTIVE });
    s._sweepFmax = fStop;    // size the skin band for the whole range, like a sweep does
    const out = [];
    for (const f of compareFs) out.push(pickModes(await s.computeAtFrequency(f, null)));
    return out;
}

// Frequency-list sweep through solve_sweep (f_max solved first, then ascending);
// cachesOff makes every point an exact solve on the same mesh.
async function solveList(spec, listFs, cachesOff) {
    const s = buildSolver(spec, 'triangular');
    s.tri_opts = cachesOff ? { lossMethod: 'auto', dispTol: 0, mqsInterpTol: 0 } : { lossMethod: 'auto' };
    const r = await s.solve_sweep({ frequencies: listFs, energy_tol: ADAPTIVE.energy_tol, max_nodes: ADAPTIVE.max_nodes });
    return listFs.map((_, i) => r.modes.map(m => ({
        mode: m.mode, R: m.RLGC.R[i], L: m.RLGC.L[i], G: m.RLGC.G[i], C: m.RLGC.C[i],
    })));
}

function fmtSpec(spec) {
    const u = (x, sc = 1e3, d = 3) => (x * sc).toFixed(d);
    let s = spec.tl;
    if (spec.tl === 'broadside_stripline') s += ` w=${u(spec.bs_w)}mm t=${u(spec.bs_t, 1e6, 1)}µm`;
    else s += ` w=${u(spec.w)}mm h=${u(spec.h)}mm t=${u(spec.t, 1e6, 1)}µm er=${spec.er.toFixed(2)}`;
    if (spec.trace_spacing) s += ` gap=${u(spec.trace_spacing)}mm`;
    const feats = ['use_sm', 'use_top_diel', 'use_gnd_cut', 'use_enclosure', 'use_side_gnd', 'use_top_gnd', 'use_plating', 'use_causal']
        .filter(k => spec[k]).map(k => k.replace('use_', ''));
    if (feats.length) s += ` [${feats.join(',')}]`;
    return s;
}

async function main() {
    console.log(`\n### Fuzz full-wave typical-use (interp sweep + caches) vs exact reference — N=${N} seed=${SEED} ###`);
    console.log(`gates: ${Object.entries(GATES).map(([q, g]) => `${q} ${(g * 100).toFixed(0)}%`).join('  ')}\n`);
    const rng = createRng(SEED);
    let ok = 0, skipped = 0;
    const flagged = [];
    for (let i = 0; i < N; i++) {
        const spec = randomSpec(rng);
        const fStop = rng.logf(5e9, 20e9);
        const compareFs = [];
        for (let k = 0; k < 5; k++) compareFs.push(F_START * Math.pow(fStop / F_START, k / 4));
        const listFs = [];
        for (let k = 0; k < N_LIST; k++) listFs.push(F_START * Math.pow(fStop / F_START, k / (N_LIST - 1)));
        if (ONLY && !ONLY.has(i)) { skipped++; continue; }
        const quiet = console.log; console.log = () => {};
        let typ = null, ref = null, typErr = null, refErr = null;
        let typL = null, refL = null, listErr = null;
        try { typ = await solveTypical(spec, fStop, compareFs); } catch (e) { typErr = e.message || String(e); }
        try { ref = await solveReference(spec, fStop, compareFs); } catch (e) { refErr = e.message || String(e); }
        if (!typErr && !refErr) {
            try {
                typL = await solveList(spec, listFs, false);
                refL = await solveList(spec, listFs, true);
            } catch (e) { listErr = e.message || String(e); }
        }
        console.log = quiet;
        if (typErr && refErr) { skipped++; continue; }   // invalid random geometry
        if (typErr || refErr) {
            const who = typErr ? `typical-use threw (reference OK): ${typErr}` : `reference threw (typical-use OK): ${refErr}`;
            flagged.push({ i, spec, why: 'one-sided failure: ' + who });
            console.log(`[${i}] ✗ ONE-SIDED FAILURE — ${who}\n      ${fmtSpec(spec)}`);
            continue;
        }
        if (listErr) {
            flagged.push({ i, spec, why: 'frequency-list sweep threw: ' + listErr });
            console.log(`[${i}] ✗ LIST SWEEP FAILURE — ${listErr}\n      ${fmtSpec(spec)}`);
            continue;
        }
        const worst = worstOf(typ, ref, compareFs, GATES);
        const worstL = worstOf(typL, refL, listFs, LIST_GATES);
        const bad = [];
        if (worst && worst.over > 1) bad.push(`MISMATCH ${fmtMismatch(worst)}`);
        if (worstL && worstL.over > 1) bad.push(`LIST MISMATCH ${fmtMismatch(worstL)}`);
        if (bad.length) {
            flagged.push({ i, spec, why: bad.join('; ') });
            console.log(`[${i}] ✗ ${bad.join('\n      ')}\n      ${fmtSpec(spec)}`);
        } else {
            ok++;
            const w = worst ? `worst ${fmtWorst(worst)}` : 'no comparable quantities';
            const wl = worstL ? `list ${fmtWorst(worstL)}` : 'list n/a';
            console.log(`[${i}] ✓ ${fmtSpec(spec)}  (${w} | ${wl})`);
        }
    }
    console.log(`\n${'='.repeat(72)}`);
    console.log(`SUMMARY: ${ok}/${N - skipped} clean | ${flagged.length} flagged | ${skipped} skipped(invalid)`);
    console.log(`${'='.repeat(72)}`);
    if (flagged.length) process.exitCode = 1;
}

import { pathToFileURL } from 'url';
if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) main();
