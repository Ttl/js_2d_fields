// Full-wave (triangular FEM) edge cases from docs/untested_edge_cases.md (F1..F10):
// regimes outside the fuzzer's draw ranges where the backend used to return a
// plausible number silently.
//
//   F1 thick wall-absorbed ground at high frequency: the wall slab impedance must not
//      overflow (d/delta > ~355), which used to reject the MQS result and drop to the
//      perturbation estimate with a 2x step in R and no warning.
//   F2 plated conductor in the skin transition: the layered plating-over-bulk surface
//      impedance needs a bulk thick against its skin depth. Both backends must flag a
//      plated film whose bulk under the plating is thinner than two skin depths.
//   F3 perturbation-path internal inductance: bounded by the finite-thickness slab
//      reactance, not by an arbitrary 0.5 L_ext cap. Coax at 1 kHz must sit near the
//      DC internal inductance of the inner wire plus the shield slab, and a
//      forced-perturbation microstrip must stay in the neighbourhood of the MQS value.
//   F8 plating face classification tolerance must scale with the conductor, not the
//      domain: a 20 nm trace in a 50 mm enclosure must keep its top and bottom faces
//      distinct.
//   F10 the eigen-anchor fallback must surface its warning and keep eps_eff sane.
//
// Run: node tests/test_fullwave_edge_cases.js
import { MicrostripSolver } from '../src/microstrip.js';
import { CoaxSolver } from '../src/coax.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    console.log(`  ${cond ? 'PASS' : 'FAIL'}: ${name}${detail ? ` (${detail})` : ''}`);
    cond ? pass++ : fail++;
}
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};
const APP = { max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
const MU0 = 4 * Math.PI * 1e-7;

const MS = {
    trace_width: 0.3e-3, substrate_height: 0.254e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
    epsilon_r: 4, tan_delta: 0.01, sigma_cond: 5.8e7, freq: 1e9, rq: 0, nx: 30, ny: 30,
    boundaries: ['open', 'open', 'open', 'gnd'],
};
function ms(opts, backend = 'triangular', tri = { lossMethod: 'auto' }) {
    const s = new MicrostripSolver({ ...MS, ...opts, mesh_backend: backend });
    s.use_causal_materials = false;
    if (backend === 'triangular') s.tri_opts = tri;
    return s;
}
async function solved(s) { const r = await quiet(() => s.solve_adaptive({ ...APP })); return { s, r }; }
const at = (s, r, f) => quiet(() => s.computeAtFrequency(f, r));
const warnsOf = (s, r) => [...(r.warnings || []), ...(s.modeWarnings || [])];

// F1: 105 um ground at 49 and 50 GHz.
{
    console.log('F1 thick ground plane across the wall-slab overflow point');
    const { s, r } = await solved(ms({ gnd_thickness: 105e-6, freq: 60e9 }));
    const a = await at(s, r, 49e9), b = await at(s, r, 50e9);
    check('MQS path at 49 GHz', a.modes[0].lossVia === 'mqs', a.modes[0].lossVia);
    check('MQS path at 50 GHz', b.modes[0].lossVia === 'mqs', b.modes[0].lossVia);
    check('R continuous across 49 -> 50 GHz (< 4%)', rel(a.modes[0].RLGC.R, b.modes[0].RLGC.R) < 0.04,
        `${a.modes[0].RLGC.R.toFixed(2)} vs ${b.modes[0].RLGC.R.toFixed(2)} ohm/m`);
    check('no MQS rejection warning', !warnsOf(s, b).some(w => w.type === 'mqs-rejected'));
}

// F2: plating-transition warning on both backends.
for (const backend of ['rectilinear', 'triangular']) {
    console.log(`F2 plating in the skin transition [${backend}]`);
    const tin = (thickness) => ({ sigma: 8.7e6, thickness, rq: 0, top: true, sides: true, bottom: false, thick_corners: true });
    const hasWarn = (s, r) => warnsOf(s, r).some(w => w.reason === 'plating-transition');
    const thin = await solved(ms({ trace_thickness: 1e-6, plating: tin(0.5e-6) }, backend));
    check('1 um Cu + 0.5 um Sn at 1 GHz warns', hasWarn(thin.s, thin.r));
    const thick = await solved(ms({ plating: tin(4e-6) }, backend));
    check('35 um Cu + 4 um Sn at 1 GHz does not warn', !hasWarn(thick.s, thick.r));
    const thickLow = await at(thick.s, thick.r, 1e6);
    check('35 um Cu + 4 um Sn at 1 MHz warns (bulk under the plating thinner than 2 delta)',
        (thickLow.warnings || []).some(w => w.reason === 'plating-transition') || hasWarn(thick.s, thickLow));
}

// F3: perturbation-path internal inductance near DC.
{
    console.log('F3 coax internal inductance near DC (perturbation path)');
    const a = 0.46e-3, b = 1.475e-3;
    const cx = new CoaxSolver({ inner_diameter: 2 * a, dielectric_diameter: 2 * b, epsilon_r: 2.1, tan_delta: 2e-4, sigma_cond: 5.8e7, freq: 1e9 });
    cx.use_causal_materials = false;
    const r = await quiet(() => cx.solve_adaptive({ max_nodes: 20000 }));
    const freqs = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10];
    const lint = [];
    for (const f of freqs) lint.push((await at(cx, r, f)).modes[0].L_internal);
    const Lext = (await at(cx, r, 1e10)).modes[0].L_external;
    // Wire: mu0/(8 pi). Shield: a slab of the modelled shield thickness on the
    // circumference 2 pi b, mu0 t/(3 * 2 pi b).
    const expect = MU0 / (8 * Math.PI) + MU0 * cx.shield_thickness / (3 * 2 * Math.PI * b);
    check('L_int(1 kHz) within 10% of wire + shield slab DC value', rel(lint[0], expect) < 0.10,
        `${(lint[0] * 1e9).toFixed(2)} vs ${(expect * 1e9).toFixed(2)} nH/m`);
    check('L_int(1 kHz) is not the 0.5 L_ext cap', lint[0] < 0.4 * Lext,
        `${(lint[0] * 1e9).toFixed(2)} vs cap ${(0.5 * Lext * 1e9).toFixed(2)} nH/m`);
    let monotone = true;
    for (let k = 1; k < lint.length; k++) if (lint[k] > lint[k - 1] * 1.001) monotone = false;
    check('L_int non-increasing 1 kHz .. 10 GHz', monotone, lint.map(v => (v * 1e9).toFixed(2)).join(' '));

    console.log('F3 microstrip forced perturbation vs MQS at 1 MHz');
    const pert = await solved(ms({}, 'triangular', { lossMethod: 'perturbation' }));
    const mqs = await solved(ms({}, 'triangular', { lossMethod: 'mqs' }));
    const lp = (await at(pert.s, pert.r, 1e6)).modes[0], lm = (await at(mqs.s, mqs.r, 1e6)).modes[0];
    check('perturbation L_int within a factor of two of MQS', lp.L_internal / lm.L_internal > 0.5 && lp.L_internal / lm.L_internal < 2,
        `${(lp.L_internal * 1e9).toFixed(2)} vs ${(lm.L_internal * 1e9).toFixed(2)} nH/m`);
    check('perturbation L_int below the old cap', lp.L_internal < 0.4 * lp.L_external);
}

// F8: plating face classification on a thin trace in a wide domain.
{
    console.log('F8 plating faces on a 20 nm trace in a 50 mm enclosure');
    const pl = (top, bottom) => ({ sigma: 1e7, thickness: 10e-9, rq: 1e-6, top, sides: false, bottom, thick_corners: false });
    const geo = { trace_thickness: 20e-9, enclosure_width: 50e-3, boundaries: ['gnd', 'gnd', 'open', 'gnd'], freq: 5e9 };
    const topOnly = await solved(ms({ ...geo, plating: pl(true, false) }));
    const botOnly = await solved(ms({ ...geo, plating: pl(false, true) }));
    const bare = await solved(ms({ ...geo }));
    const Rt = topOnly.r.modes[0].RLGC.R, Rb = botOnly.r.modes[0].RLGC.R, R0 = bare.r.modes[0].RLGC.R;
    check('rough poor plating on the bottom face costs more than on the top face', Rb > Rt * 1.02,
        `bottom-only ${Rb.toFixed(1)}, top-only ${Rt.toFixed(1)}, bare ${R0.toFixed(1)} ohm/m`);
    check('both plated cases exceed bare', Rt > R0 && Rb > R0);
}

// F10: eigen-anchor fallback warning.
{
    console.log('F10 eigen-anchor fallback');
    const s = ms({ freq: 5e9 });
    const r = await quiet(() => s.solve_adaptive({ ...APP }));
    const tri = await s._ensureTriBackend();
    const orig = tri._eigenPick.bind(tri);
    tri._eigenBiasCache = null;
    tri._eigenPick = (st, fa, phiEps, epsStatic) =>
        fa <= 1e9 ? { fw: null, fwErr: new Error('injected anchor failure') } : orig(st, fa, phiEps, epsStatic);
    const rf = await at(s, r, 5e9);
    tri._eigenPick = orig;
    const w = warnsOf(s, rf).find(w => w.type === 'eigen-anchor');
    check('eigen-anchor warning surfaced', !!w, w ? w.message.slice(0, 80) : 'none');
    check('eps_eff stays within 5% of the static value', rel(rf.modes[0].eps_eff, r.modes[0].eps_eff) < 0.05,
        `${rf.modes[0].eps_eff.toFixed(4)} vs ${r.modes[0].eps_eff.toFixed(4)}`);
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
