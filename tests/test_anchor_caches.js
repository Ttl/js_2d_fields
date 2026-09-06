// Full-wave anchor caches vs exact solves on the same mesh.
//
// On a sweep the triangular backend interpolates three per-mode scalars between
// exact anchors: eps_d (dispersion cache, dispTol) and the MQS pair R and
// X_int = omega*L_internal (mqsInterpTol). A leave-one-out check on the anchors
// next to the query decides whether to interpolate, and nothing downstream
// re-verifies the result: the InterpolatingSweep's midpoint check compares its
// spline against computeAtFrequency, which returns the interpolated value once
// the gate has passed. So the gate is the only error control for these
// quantities, and it has to hold for every order the sweep APIs feed
// frequencies in:
//   * solve_sweep with a frequency list: the adaptive solve at f_max seeds one
//     anchor at the top and the points then come ascending, so after four
//     anchors at the bottom every later point is interior to a gap spanning the
//     whole sweep. A gate that only tests the neighbouring anchor is fooled here
//     (a cubic is always accurate next to its own end points) and a thin-film
//     R came out up to 8% off.
//   * InterpolatingSweep: log-uniform initial points, then interval midpoints,
//     so the anchors are always central to their gaps.
// Both orders are checked against the same solver with the caches off
// (dispTol = mqsInterpTol = 0, every point an exact solve) on the identical
// mesh, so the gates are the caches' own tolerances (0.2% on R and X_int, 0.1%
// on eps_d) with a little headroom, not the several-percent stack-up the
// fullwave fuzzer needs. TRI_STATS counts the eigen and MQS solves so the test
// also fails if a change quietly turns a cache off on the bisection order.
//
// Usage: node tests/test_anchor_caches.js

process.env.TRI_STATS = '1';   // before the solver modules load
const { MicrostripSolver } = await import('../src/microstrip.js');
const { InterpolatingSweep } = await import('../src/interpolating_sweep.js');

const F0 = 1e7, F1 = 2e10, N_LIST = 25;
const LIST = Array.from({ length: N_LIST }, (_, i) => F0 * Math.pow(F1 / F0, i / (N_LIST - 1)));
const ADAPT = { energy_tol: 0.01, max_nodes: 18000 };
// R and X_int anchors are gated at 0.2%, eps_d at 0.1%; X_int only reaches L
// through L_internal/L < 1 and eps_d reaches C and G one to one.
const GATES = { R: 0.0025, L: 0.0025, G: 0.0012, C: 0.0012 };

const BASE = { substrate_height: 0.254e-3, gnd_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02,
               sigma_cond: 5.8e7, freq: 1e9, boundaries: ['open', 'open', 'open', 'gnd'],
               nx: 30, ny: 30, mesh_backend: 'triangular' };
const CASES = [
    // The reported case: R changes only 1.6x over the sweep, so a cubic across
    // the whole span passes an edge-anchor check while being 8% off inside.
    // R and L_internal both have a knee where the film starts to shield, so
    // at the sweep's anchor spacing the gate rightly accepts few points.
    { name: 'microstrip 100 nm', o: { trace_width: 0.3e-3, trace_thickness: 100e-9 } },
    // Skin transition (delta = t near 4 GHz) inside the sweep.
    { name: 'microstrip 1 um',   o: { trace_width: 0.3e-3, trace_thickness: 1e-6 } },
    // Two modes on the half domain, smooth sqrt(f) loss: here the caches must
    // actually skip solves.
    { name: 'diff microstrip 35 um', o: { trace_width: 0.2e-3, trace_spacing: 0.15e-3, trace_thickness: 35e-6 }, expectHits: true },
];

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
const quiet = async fn => { const log = console.log; console.log = () => {}; try { return await fn(); } finally { console.log = log; } };
const mk = (o, cachesOff) => {
    const s = new MicrostripSolver({ ...BASE, ...o });
    s.tri_opts = cachesOff ? { dispTol: 0, mqsInterpTol: 0 } : {};
    return s;
};
const pick = r => Object.fromEntries(r.modes.map(m => [m.mode, { R: m.RLGC.R, L: m.RLGC.L, G: m.RLGC.G, C: m.RLGC.C }]));
// The stats object appears when the backend module loads (lazily, at the first solve).
const counts = () => { const st = globalThis.__TRI_STATS__ || { eig: [], lin: [] };
    return { eig: st.eig.length, mqs: st.lin.filter(x => x.complex).length }; };

// Worst relative difference per quantity over all frequencies and modes.
function compare(label, freqs, typ, ref) {
    const worst = {};
    for (const q of Object.keys(GATES)) worst[q] = { d: 0, f: 0, mode: '' };
    for (let i = 0; i < freqs.length; i++)
        for (const mode of Object.keys(ref[i]))
            for (const q of Object.keys(GATES)) {
                const a = ref[i][mode][q], b = typ[i][mode][q];
                const d = Math.abs(a - b) / Math.max(Math.abs(a), 1e-30);
                if (d > worst[q].d) worst[q] = { d, f: freqs[i], mode };
            }
    for (const q of Object.keys(GATES)) {
        const w = worst[q];
        check(`${label}: ${q} within ${(GATES[q] * 100).toFixed(2)}%`, w.d <= GATES[q],
              `worst ${(w.d * 100).toFixed(3)}% at ${(w.f / 1e9).toFixed(3)} GHz [${w.mode}]`);
    }
}

for (const c of CASES) {
    console.log(`\n=== ${c.name} ===`);
    const nModes = c.o.trace_spacing ? 2 : 1;

    // Frequency list through solve_sweep: exact reference and the cached run on
    // the same mesh.
    const ref = mk(c.o, true), typ = mk(c.o, false);
    const listOf = r => LIST.map((_, i) => Object.fromEntries(r.modes.map(m =>
        [m.mode, { R: m.RLGC.R[i], L: m.RLGC.L[i], G: m.RLGC.G[i], C: m.RLGC.C[i] }])));
    const refList = listOf(await quiet(() => ref.solve_sweep({ frequencies: LIST, ...ADAPT })));
    const c0 = counts();
    const typList = listOf(await quiet(() => typ.solve_sweep({ frequencies: LIST, ...ADAPT })));
    const c1 = counts();
    console.log(`  list order: ${N_LIST} points, ${c1.mqs - c0.mqs} MQS solves, ${c1.eig - c0.eig} eigensolves`);
    compare('list order', LIST, typList, refList);

    // InterpolatingSweep order: the samples it solved, each against an exact
    // solve on the reference solver (skin band sized for the sweep like the
    // sweep does).
    const ref2 = mk(c.o, true), typ2 = mk(c.o, false);
    const refInit = await quiet(() => ref2.solve_adaptive({ ...ADAPT }));
    ref2._sweepFmax = F1;
    const typInit = await quiet(() => typ2.solve_adaptive({ ...ADAPT }));
    const s0 = counts();
    const sweep = new InterpolatingSweep(typ2, typInit, { tolerance: 0.005 });
    await quiet(() => sweep.run(F0, F1, {}));
    const s1 = counts();
    const ts = [...sweep.samplePoints.keys()].sort((a, b) => a - b);
    const freqs = ts.map(t => Math.pow(10, t));
    const typS = ts.map(t => Object.fromEntries(sweep.samplePoints.get(t).map(m => [m.mode, m])));
    const refS = [];
    for (const f of freqs) refS.push(pick(await quiet(() => ref2.computeAtFrequency(f, refInit))));
    const nExact = ts.filter(t => Math.pow(10, t) >= 1e8).length;   // eigensolves only above F_STATIC_MAX
    const mqs = s1.mqs - s0.mqs, eig = s1.eig - s0.eig;
    console.log(`  bisection order: ${ts.length} samples, ${mqs} MQS solves, ${eig} eigensolves`);
    compare('bisection order', freqs, typS, refS);
    // The caches must be live on this order: well under one solve per sample
    // and mode (the eigen count includes the static-anchor walk).
    if (c.expectHits) {
        check('bisection order: MQS cache hits', mqs <= 0.75 * nModes * ts.length, `${mqs} solves for ${ts.length} samples`);
        check('bisection order: eps cache hits', eig <= 0.75 * nModes * nExact + 4, `${eig} solves for ${nExact} full-wave samples`);
    }
}

console.log(`\n${failures === 0 ? 'ALL PASSED' : `${failures} FAILED`}`);
process.exit(failures === 0 ? 0 : 1);
