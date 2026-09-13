// Robustness at unrealistic inputs (docs/untested_edge_cases.md, X items): the
// solver may be inaccurate there, but it must not throw, return NaN, or return an
// answer that is wrong by orders of magnitude.
//
//   X1 interpolating sweep below 1 Hz: sample keys are the requested log-frequency,
//      so the refinement loop finds its own midpoints (log10(10^t) != t in the last
//      bit for negative t).
//   X2 S-parameters of a line attenuated beyond e^-300: finite, S21 ~ 0 and S11 at the
//      semi-infinite-line reflection, on the single-ended and the 4-port MTL path.
//   X3 causal materials with a loss tangent the wideband Debye fit cannot represent:
//      finite results, the nominal material kept, and a causal-model warning on both
//      backends (the fit used to drive epsilon_r negative and the Laplace solve threw).
//   X4 conductor gaps far below the mesh's merge tolerance keep their own cells on
//      the FDM grid: a 1 nm gap gives a capacitance, not a short and not a snap.
//   X5 the remaining unbounded exponentials: complex tanh saturates instead of
//      dividing two overflowed hyperbolics (layered plating impedance of a metre of
//      plating at 100 GHz), and the causal model's complex log takes no squared
//      frequency (1e200 Hz).
//
// Run: node tests/test_robustness_extremes.js
import { MicrostripSolver } from '../src/microstrip.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';
import { computeSParamsSingleEnded, computeSParamsDifferentialMTL } from '../src/sparameters.js';
import { Complex } from '../src/complex.js';
import { calculate_Zrough, calculate_Zrough_layered } from '../src/surface_roughness.js';
import { djordjevic_sarkar } from '../src/djordjevic_sarkar.js';

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
const fin = z => z && Number.isFinite(z.re) && Number.isFinite(z.im);

const MS = {
    trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0, nx: 30, ny: 30,
    boundaries: ['open', 'open', 'open', 'gnd'],
};
function ms(opts, backend = 'rectilinear', causal = false) {
    const s = new MicrostripSolver({ ...MS, ...opts, mesh_backend: backend });
    s.use_causal_materials = causal;
    if (backend === 'triangular') s.tri_opts = { lossMethod: 'auto' };
    return s;
}
async function solved(s) { const r = await quiet(() => s.solve_adaptive({ ...APP })); return { s, r }; }
const at = (s, r, f) => quiet(() => s.computeAtFrequency(f, r));
const finiteModes = r => (r.modes || []).every(m => [m.Z0, m.eps_eff, m.RLGC.R, m.RLGC.L, m.RLGC.G, m.RLGC.C].every(Number.isFinite));

// X1
{
    console.log('X1 interpolating sweep from 0.1 Hz to 1 Hz');
    const { s, r } = await solved(ms({ freq: 1 }));
    let n = -1, err = null;
    try { const sw = new InterpolatingSweep(s, r, { tolerance: 0.005 }); n = await quiet(() => sw.run(0.1, 1, {})); const b = sw.buildResults([0.3])[0].result; check('interpolated result finite', finiteModes(b)); }
    catch (e) { err = e.message; }
    check('sweep completes', err === null && n > 0, err || `${n} samples`);
}

// X2
{
    console.log('X2 S-parameters beyond e^-300 attenuation');
    const rlgc = { R: 1.5e12, L: 3e-7, G: 0, C: 1.2e-10 };
    const sp = computeSParamsSingleEnded(1e9, rlgc, 0.01, 50);
    check('single-ended entries finite', fin(sp.S11) && fin(sp.S21));
    check('S21 ~ 0', fin(sp.S21) && Math.hypot(sp.S21.re, sp.S21.im) < 1e-100, fin(sp.S21) ? Math.hypot(sp.S21.re, sp.S21.im).toExponential(2) : 'NaN');
    // Semi-infinite line: S11 = (Zc - Zr)/(Zc + Zr).
    const omega = 2 * Math.PI * 1e9;
    const zc = Math.sqrt(rlgc.R / (omega * rlgc.C));   // |Zc| for R >> omega L, with 45 deg phase
    check('S11 magnitude at the semi-infinite-line value', fin(sp.S11) && Math.abs(Math.hypot(sp.S11.re, sp.S11.im) - 1) < 0.05,
        fin(sp.S11) ? `|S11| ${Math.hypot(sp.S11.re, sp.S11.im).toFixed(4)}, |Zc| ${zc.toExponential(2)}` : 'NaN');
    const R2 = [[1.5e12, 1e11], [1e11, 1.6e12]], L2 = [[3e-7, 5e-8], [5e-8, 3.1e-7]], G2 = [[0, 0], [0, 0]], C2 = [[1.2e-10, -2e-11], [-2e-11, 1.3e-10]];
    const sm = computeSParamsDifferentialMTL(1e9, R2, L2, G2, C2, 0.01, 50);
    check('MTL entries finite', sm.S.flat().every(fin));
    check('MTL SDD21 ~ 0', fin(sm.SDD21) && Math.hypot(sm.SDD21.re, sm.SDD21.im) < 1e-100);
    const sLong = computeSParamsSingleEnded(1e9, { R: 20, L: 3e-7, G: 1e-3, C: 1.2e-10 }, 1e6, 50);
    check('1000 km line finite', fin(sLong.S11) && fin(sLong.S21));
}

// X3
for (const backend of ['rectilinear', 'triangular']) {
    console.log(`X3 causal materials with tan_delta = 1 [${backend}]`);
    const { s, r } = await solved(ms({ tan_delta: 1 }, backend, true));
    let rf = null, err = null;
    try { rf = await at(s, r, 10e9); } catch (e) { err = e.message; }
    check('10 GHz point does not throw', err === null, err || '');
    if (rf) {
        check('10 GHz result finite', finiteModes(rf));
        check('eps_eff stays physical', rf.modes[0].eps_eff > 0.9 && rf.modes[0].RLGC.C > 0, `eps ${rf.modes[0].eps_eff.toFixed(3)}`);
        const warns = [...(rf.warnings || []), ...(s.modeWarnings || [])];
        check('causal-model warning present', warns.some(w => w.reason === 'causal-model'), warns.map(w => w.reason || w.type).join(',') || 'none');
    }
}

// X4
{
    console.log('X4 nanometre conductor gaps on the FDM grid');
    // Trace bottom 1 nm above the ground plane: essentially a parallel-plate
    // capacitor eps0 er w / gap; compare against the full-wave result.
    const gap = 1e-9;
    const q = await solved(ms({ trace_thickness: -(MS.substrate_height - gap) }, 'rectilinear'));
    const t = await solved(ms({ trace_thickness: -(MS.substrate_height - gap) }, 'triangular'));
    const Cpp = 8.854e-12 * MS.epsilon_r * MS.trace_width / gap;
    check('FDM capacitance within 2x of the parallel-plate estimate', q.r.modes[0].RLGC.C > 0.5 * Cpp && q.r.modes[0].RLGC.C < 2 * Cpp,
        `${(q.r.modes[0].RLGC.C * 1e6).toFixed(2)} vs ${(Cpp * 1e6).toFixed(2)} uF/m`);
    check('FDM Z0 within 2x of full-wave', rel(q.r.modes[0].Z0, t.r.modes[0].Z0) < 0.5,
        `${q.r.modes[0].Z0.toExponential(3)} vs ${t.r.modes[0].Z0.toExponential(3)}`);
    const gq = await solved(ms({ use_coplanar_gnd: true, use_vias: true, gap, via_gap: 0.1e-3 }, 'rectilinear'));
    const gt = await solved(ms({ use_coplanar_gnd: true, use_vias: true, gap, via_gap: 0.1e-3 }, 'triangular'));
    check('GCPW 1 nm slot: FDM Z0 within 2x of full-wave', rel(gq.r.modes[0].Z0, gt.r.modes[0].Z0) < 0.5,
        `${gq.r.modes[0].Z0.toExponential(3)} vs ${gt.r.modes[0].Z0.toExponential(3)}`);
}

// X5
{
    console.log('X5 unbounded exponentials');
    const th = new Complex(1000, 1).tanh(), tn = new Complex(-1000, 0.5).tanh();
    check('complex tanh saturates at +1 / -1', th.re === 1 && th.im === 0 && tn.re === -1, `${th} ${tn}`);
    const zl = calculate_Zrough_layered(1e11, 5.8e7, 0, 1e7, 1), zs = calculate_Zrough(1e11, 1e7, 0);
    check('1 m plating at 100 GHz: layered Zs finite and equal to the solid plating metal',
        fin(zl) && rel(zl.re, zs.re) < 1e-6 && rel(zl.im, zs.im) < 1e-6, `${zl} vs ${zs}`);
    const d1 = djordjevic_sarkar(1e200, 4.4, 0.02), d2 = djordjevic_sarkar(1e12, 4.4, 0.02);
    check('causal model finite at 1e200 Hz and at its high-frequency limit', Number.isFinite(d1.eps_real) && Number.isFinite(d1.tand_actual)
        && d1.eps_real > 1 && Math.abs(d1.eps_real - d2.eps_real) < 0.1, `${JSON.stringify(d1)} vs ${JSON.stringify(d2)}`);
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
