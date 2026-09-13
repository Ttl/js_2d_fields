// Shared-code edge cases from docs/untested_edge_cases.md (S1..S5): paths both
// backends and the app run through, outside the fuzzer's draw ranges.
//
//   S1 asymmetric pair at DC: the MTL 4-port S-parameters must be finite at f = 0,
//      reciprocal, and continuous with the f -> 0 limit (a series R network).
//   S2 DC point of a sweep: L at exactly f = 0 must carry the DC internal
//      inductance, continuous with the near-DC plateau, on both backends.
//   S3 interpolating sweep vs discrete solve: at a sweep sample frequency the splines
//      are exact, so the derived Z0 and alpha_d must match the discrete solve's
//      conventions (static Z0, alpha_d = G Z0 / 2) at every frequency, including
//      well below the skin regime where Re(Zc) and Z0 differ by tens of percent.
//   S4 plating thicker than the trace: the conductor is entirely plating metal, so
//      R can never fall below the DC resistance of a solid plating-metal trace.
//      Both backends, thin film in the skin transition.
//   S5 rectangular waveguide with b > a: the second cutoff is the lower of TE(1,0)
//      and TE(0,2), and the over-moded warning starts there, not at the fundamental.
//
// Run: node tests/test_shared_edge_cases.js
import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';
import { CoaxSolver } from '../src/coax.js';
import { RectWaveguideSolver } from '../src/rect_waveguide.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';
import { computeSParamsDiffAuto } from '../src/sparameters.js';

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
const C0 = 299792458;

const MS = {
    trace_width: 0.2e-3, substrate_height: 0.1e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
    epsilon_r: 4, tan_delta: 0.01, sigma_cond: 5.8e7, freq: 1e9, rq: 0, nx: 30, ny: 30,
    boundaries: ['open', 'open', 'open', 'gnd'],
};
function ms(opts, backend) {
    const s = new MicrostripSolver({ ...MS, ...opts, mesh_backend: backend });
    if (backend === 'triangular') s.tri_opts = { lossMethod: 'auto' };
    return s;
}
async function solved(s) { const r = await quiet(() => s.solve_adaptive({ ...APP })); return { s, r }; }
const at = (s, r, f) => quiet(() => s.computeAtFrequency(f, r));

// S1: asymmetric broadside pair, MTL 4-port at DC.
{
    console.log('S1 asymmetric pair 4-port S-parameters at DC');
    const bs = new BroadsideStriplineSolver({
        trace_width: 0.2e-3, trace_thickness: 35e-6, x_offset: 50e-6, sigma_cond: 5.8e7,
        h_bottom: 0.2e-3, er_bottom: 4.4, tand_bottom: 0.02, h_middle: 0.15e-3, er_middle: 3, tand_middle: 0.01,
        h_top: 0.3e-3, er_top: 4.4, tand_top: 0.02, freq: 1e9, rq: 0, nx: 30, ny: 30,
        boundaries: ['open', 'open', 'gnd', 'gnd'], mesh_backend: 'rectilinear',
    });
    const { r } = await solved(bs);
    const sp = async (f) => {
        const rf = await at(bs, r, f);
        const odd = rf.modes.find(m => m.mode === 'odd'), even = rf.modes.find(m => m.mode === 'even');
        check(`physMatrix present at ${f} Hz (asymmetric pair routes to the MTL path)`, !!rf.physMatrix);
        return computeSParamsDiffAuto(f, odd.RLGC, even.RLGC, rf.physMatrix, 0.01, 50);
    };
    const s0 = await sp(0), s1 = await sp(1);
    const flat = S => S.flat();
    check('all 16 entries finite at DC', flat(s0.S).every(z => Number.isFinite(z.re) && Number.isFinite(z.im)));
    let maxAsym = 0, maxDiff = 0;
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
        maxAsym = Math.max(maxAsym, Math.hypot(s0.S[i][j].re - s0.S[j][i].re, s0.S[i][j].im - s0.S[j][i].im));
        maxDiff = Math.max(maxDiff, Math.hypot(s0.S[i][j].re - s1.S[i][j].re, s0.S[i][j].im - s1.S[i][j].im));
    }
    check('reciprocal at DC', maxAsym < 1e-12, `max |S_ij - S_ji| = ${maxAsym.toExponential(2)}`);
    check('continuous with f = 1 Hz', maxDiff < 1e-6, `max |S(0) - S(1 Hz)| = ${maxDiff.toExponential(2)}`);
    check('SDD21 ~ 1 for a 10 mm line at DC', Math.abs(s0.SDD21.re - 1) < 0.01 && Math.abs(s0.SDD21.im) < 1e-9,
        `${s0.SDD21.re.toFixed(5)} ${s0.SDD21.im.toExponential(2)}j`);
}

// S2: L at f = 0 continuous with the near-DC plateau.
for (const backend of ['rectilinear', 'triangular']) {
    console.log(`S2 DC inductance [${backend}]`);
    const { s, r } = await solved(ms({}, backend));
    const L0 = (await at(s, r, 0)).modes[0], L1 = (await at(s, r, 1e3)).modes[0];
    check('L_internal(0) > 0', L0.L_internal > 0, `${(L0.L_internal * 1e9).toFixed(2)} nH/m`);
    check('L(0) within 0.5% of L(1 kHz)', rel(L0.RLGC.L, L1.RLGC.L) < 0.005,
        `${(L0.RLGC.L * 1e9).toFixed(2)} vs ${(L1.RLGC.L * 1e9).toFixed(2)} nH/m`);
    check('R(0) is the DC resistance', rel(L0.RLGC.R, L1.RLGC.R) < 0.02,
        `${L0.RLGC.R.toFixed(3)} vs ${L1.RLGC.R.toFixed(3)} ohm/m`);
}

// S3: sweep-built derived quantities vs the discrete solve at sample frequencies.
for (const backend of ['rectilinear', 'triangular']) {
    console.log(`S3 interpolating sweep conventions [${backend}]`);
    const fMin = 1e6, fMax = 10e9;
    const { s, r } = await solved(ms({ freq: fMax }, backend));
    const sweep = new InterpolatingSweep(s, r, { tolerance: 0.005, initialPoints: 6 });
    await quiet(() => sweep.run(fMin, fMax, {}));
    // The endpoints are always sample points, so the splines are exact there and any
    // difference is a definition difference, not interpolation error.
    for (const f of [fMin, fMax]) {
        const sw = sweep.buildResults([f])[0].result.modes[0];
        const ex = (await at(s, r, f)).modes[0];
        check(`Z0 matches the discrete static Z0 at ${(f / 1e9).toFixed(3)} GHz`, rel(sw.Z0, ex.Z0) < 1e-4,
            `sweep ${sw.Z0.toFixed(4)} vs discrete ${ex.Z0.toFixed(4)} (Re Zc ${ex.Zc.re.toFixed(4)})`);
        check(`alpha_d matches at ${(f / 1e9).toFixed(3)} GHz`, rel(sw.alpha_d, ex.alpha_d) < 1e-4,
            `sweep ${sw.alpha_d.toExponential(4)} vs discrete ${ex.alpha_d.toExponential(4)}`);
        check(`L_internal carried at ${(f / 1e9).toFixed(3)} GHz`, rel(sw.L_internal, ex.L_internal) < 1e-3,
            `sweep ${(sw.L_internal * 1e9).toFixed(3)} vs discrete ${(ex.L_internal * 1e9).toFixed(3)} nH/m`);
    }
}

// S4: plating thicker than the trace on a thin film.
{
    console.log('S4 plating thicker than the trace (100 nm Cu + 1 um sigma 1e7 plating, 5 GHz)');
    const t = 100e-9, w = 0.2e-3, sigP = 1e7;
    const opts = { trace_thickness: t, trace_width: w, freq: 5e9,
        plating: { sigma: sigP, thickness: 1e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: true } };
    const Rdc_solid = 1 / (sigP * w * t);
    const R = {};
    for (const backend of ['rectilinear', 'triangular']) {
        const { r } = await solved(ms(opts, backend));
        R[backend] = r.modes[0].RLGC.R;
        check(`[${backend}] R >= DC resistance of the solid plating-metal trace`, R[backend] >= 0.99 * Rdc_solid,
            `${R[backend].toFixed(0)} vs ${Rdc_solid.toFixed(0)} ohm/m`);
    }
    check('backends agree within 25%', rel(R.rectilinear, R.triangular) < 0.25,
        `${R.rectilinear.toFixed(0)} vs ${R.triangular.toFixed(0)} ohm/m`);
    // Plating thinner than the trace keeps the layered model: the same solver at 35 um
    // must not change, the thick-plating tests pin that regime.
    const thick = await solved(ms({ freq: 5e9, plating: { ...opts.plating } }, 'rectilinear'));
    check('35 um trace with 1 um plating still uses the layered model (R well below solid plating)',
        thick.r.modes[0].RLGC.R < 0.5 * (1 / (sigP * w * 35e-6)) + 1e3, `${thick.r.modes[0].RLGC.R.toFixed(1)} ohm/m`);
}

// S5: waveguide second cutoff with b > a.
{
    console.log('S5 rectangular waveguide second cutoff for b > a');
    const a = 10.16e-3, b = 22.86e-3;
    const g = new RectWaveguideSolver({ width: a, height: b, sigma_cond: 5.8e7, freq: 10e9 });
    // Fundamental TE(0,1): kc = pi/b. Next: TE(1,0) at pi/a or TE(0,2) at 2 pi/b.
    const fc2 = C0 * Math.min(Math.PI / a, 2 * Math.PI / b) / (2 * Math.PI);
    check('fc2 is the lower of TE(1,0) and TE(0,2)', rel(g.fc2, fc2) < 1e-9,
        `${(g.fc2 / 1e9).toFixed(3)} vs ${(fc2 / 1e9).toFixed(3)} GHz`);
    const over = f => g.waveguideWarnings(f).some(m => /over-moded/.test(m));
    check('no over-moded warning between the cutoffs', !over(1.2 * g.fc) && !over(0.95 * fc2));
    check('over-moded warning above the second cutoff', over(1.05 * fc2));
    const wide = new RectWaveguideSolver({ width: b, height: a, sigma_cond: 5.8e7, freq: 10e9 });
    check('rotating the guide leaves both cutoffs unchanged', rel(wide.fc, g.fc) < 1e-12 && rel(wide.fc2, g.fc2) < 1e-12);
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
