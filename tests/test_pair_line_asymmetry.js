// Mode conversion of a pair whose traces differ in metal or finish, on both backends.
//
// The geometry is mirror symmetric, so C, G and the external L are too, but R11 != R22
// and L11 != L22 (internal inductance). The loss evaluated with current in one trace
// supplies the two differences (full-wave: the eddy-current solve driven per trace,
// quasi-static: the surface integral of half the sum / difference of the even and odd
// vacuum fields), and the S-parameters go through the 4-port multiconductor routine.
// The quasi-static R reads a few percent low at this mesh, and so do its differences.
//
// Reference: Ansys 2D Extractor, 2 GHz, causal substrate, left trace 5.8e7 S/m, right
// trace 3.8e7 S/m (the case of solve_differential_microstrip_unequal_sigma):
//   R [27.429, 1.4922, 1.4922, 32.306] ohm/m, L [291.58, 19.723, 19.723, 291.96] nH/m.
// and the same with a right trace of 2.3e6 S/m (copper / 25, an asymmetry of the size a
// rough trace causes):
//   R [27.847, 1.4757, 1.4757, 113.11] ohm/m, L [291.61, 19.725, 19.725, 298.16] nH/m.
// G [25.039, -0.10192, -0.10192, 25.038] mS/m and C [116.07, -2.6636, -2.6636, 116.07]
// pF/m in both.
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { computeSParamsDiffAuto, computeSParamsDifferential, computeSParamsDifferentialMTL } from '../src/sparameters.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';
import { Complex } from '../src/complex.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
async function quiet(fn) {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
}
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b));

const pairText = (left, right) => `units mm
bounds gnd gnd gnd gnd
domain -1.5 1.5 0 1.2604
diel x=-inf w=inf y=0 h=0.2104 er=4.4 tand=0.02
sig+ x=-0.6 w=0.35 y=0.2104 h=0.05 ${left}
sig- x=0.25 w=0.35 y=0.2104 h=0.05 ${right}
`;
let BACKEND = 'triangular';
function build(left, right) {
    const s = new CustomGeometrySolver({ text: pairText(left, right), sigma_cond: 5.8e7, freq: 2e9, nx: 10, ny: 10,
        mesh_backend: BACKEND });
    if (BACKEND === 'triangular') s.tri_opts = { lossMethod: 'auto' };
    s.use_causal_materials = true;
    return s;
}
const SOLVE = { energy_tol: 0.002 };
const modesOf = r => [r.modes.find(m => m.mode === 'odd'), r.modes.find(m => m.mode === 'even')];
// Differential in, common out over a line: Scd21 = (S31 - S32 + S41 - S42) / 2.
function scd21(sp) {
    const S = sp.S;
    return S[2][0].sub(S[2][1]).add(S[3][0]).sub(S[3][1]).mul(new Complex(0.5, 0));
}
const dB = z => 20 * Math.log10(Math.max(z.abs(), 1e-15));
const LEN = 0.1, F = 2e9;

const G_REF = [[25.039e-3, -0.10192e-3], [-0.10192e-3, 25.038e-3]], C_REF = [[116.07e-12, -2.6636e-12], [-2.6636e-12, 116.07e-12]];
const REFS = [
    { name: 'right trace 3.8e7 S/m', right: 'sigma=3.8e7',
      R: [[27.429, 1.4922], [1.4922, 32.306]], L: [[291.58e-9, 19.723e-9], [19.723e-9, 291.96e-9]] },
    { name: 'right trace 2.3e6 S/m', right: 'sigma=2.3e6',
      R: [[27.847, 1.4757], [1.4757, 113.11]], L: [[291.61e-9, 19.725e-9], [19.725e-9, 298.16e-9]] },
];
// tolR: R11, R22 and R11 - R22. R12 is a small difference of the mode values.
const TOL = { triangular: { tolR: 0.01, tolR12: 0.03, tolL: 2e-3, tolDL: 0.1, scd: 1 },
              rectilinear: { tolR: 0.06, tolR12: 0.2, tolL: 5e-3, tolDL: 0.1, scd: 1 } };
const rough = {};

for (const [label, backend] of [['full-wave', 'triangular'], ['QS', 'rectilinear']]) {
    BACKEND = backend;
    const T = TOL[backend];
    let first = null;
    for (const ref of REFS) {
        const solver = build('', ref.right);
        const res = await quiet(() => solver.solve_adaptive(SOLVE));
        const [odd, even] = modesOf(res);
        const M = res.RLGC_matrix;
        first = first || { solver, res, odd };
        check(`${label}, ${ref.name}: R11 and R22 match the reference`,
            rel(M.R[0][0], ref.R[0][0]) < T.tolR && rel(M.R[1][1], ref.R[1][1]) < T.tolR
            && rel(M.R[0][0] - M.R[1][1], ref.R[0][0] - ref.R[1][1]) < T.tolR,
            `${M.R[0][0].toFixed(3)} / ${M.R[1][1].toFixed(3)} vs ${ref.R[0][0]} / ${ref.R[1][1]} ohm/m`);
        check(`${label}, ${ref.name}: R12 matches the reference`, rel(M.R[0][1], ref.R[0][1]) < T.tolR12,
            `${M.R[0][1].toFixed(4)} vs ${ref.R[0][1]}`);
        const dLref = ref.L[0][0] - ref.L[1][1];
        check(`${label}, ${ref.name}: L11, L22 and L11 - L22 match the reference`,
            rel(M.L[0][0], ref.L[0][0]) < T.tolL && rel(M.L[1][1], ref.L[1][1]) < T.tolL
            && rel(M.L[0][0] - M.L[1][1], dLref) < T.tolDL,
            `dL ${((M.L[0][0] - M.L[1][1]) * 1e9).toFixed(3)} vs ${(dLref * 1e9).toFixed(3)} nH/m`);
        check(`${label}, ${ref.name}: the diagonal mean is still the mode mean`,
            rel(M.R[0][0] + M.R[1][1], odd.RLGC.R + even.RLGC.R) < 1e-12);
        const refS = computeSParamsDifferentialMTL(F, ref.R, ref.L, G_REF, C_REF, LEN, 50);
        const ourS = computeSParamsDiffAuto(F, odd.RLGC, even.RLGC, res.physMatrix, LEN, 50);
        check(`${label}, ${ref.name}: mode conversion over 100 mm matches the reference matrices`,
            Math.abs(dB(scd21(ourS)) - dB(scd21(refS))) < T.scd,
            `Scd21 ${dB(scd21(ourS)).toFixed(2)} vs ${dB(scd21(refS)).toFixed(2)} dB`);
    }

    // --- Interpolating sweep carries the asymmetry ---
    {
        const { solver, res } = first;
        const sweep = new InterpolatingSweep(solver, res, { tolerance: 0.005 });
        await quiet(() => sweep.run(0.5e9, 10e9, {}));
        const fs = [0.7e9, 1.3e9, 3.1e9, 7.7e9];
        const interp = sweep.buildResults(fs);
        let worst = 0;
        for (let i = 0; i < fs.length; i++) {
            const exact = await quiet(() => solver.computeAtFrequency(fs[i], res));
            const a = modesOf(exact)[0].RLGC, b = modesOf(interp[i].result)[0].RLGC;
            worst = Math.max(worst, rel(a.dR, b.dR), rel(a.dL, b.dL));
        }
        check(`${label}: interpolated sweep points carry R11 - R22 and L11 - L22`, worst < 0.03, `worst rel ${worst.toExponential(2)}`);
        const mat = interp[2].result.RLGC_matrix;
        check(`${label}: and their matrix has unequal diagonals`, mat.R[1][1] > mat.R[0][0] * 1.1);
    }

    // --- Equal traces: nothing changes ---
    {
        const s = build('sigma=3.8e7', 'sigma=3.8e7');
        const r = await quiet(() => s.solve_adaptive(SOLVE));
        const [o, e] = modesOf(r);
        const a = computeSParamsDiffAuto(F, o.RLGC, e.RLGC, r.physMatrix, LEN, 50);
        const b = computeSParamsDifferential(F, o.RLGC, e.RLGC, LEN, 50);
        check(`${label}: equal traces carry no asymmetry and keep the odd/even S-parameters`,
            o.RLGC.dR === undefined && a.S[2][0].re === b.S[2][0].re && a.S[3][0].im === b.S[3][0].im);
    }

    // --- Swapping the traces mirrors the result ---
    {
        const s = build('sigma=3.8e7', '');
        const r = await quiet(() => s.solve_adaptive(SOLVE));
        const [o] = modesOf(r);
        check(`${label}: swapping the two metals flips the sign of the asymmetry`,
            rel(o.RLGC.dR, -first.odd.RLGC.dR) < 0.01 && rel(o.RLGC.dL, -first.odd.RLGC.dL) < 0.02,
            `dR ${o.RLGC.dR.toFixed(3)} vs ${(-first.odd.RLGC.dR).toFixed(3)}`);
    }

    // --- A rough trace converts far more, through the internal inductance ---
    {
        const s = build('', 'rq=0.001');
        const r = await quiet(() => s.solve_adaptive(SOLVE));
        const [o, e] = modesOf(r);
        const sp = computeSParamsDiffAuto(F, o.RLGC, e.RLGC, r.physMatrix, LEN, 50);
        rough[backend] = { dL: o.RLGC.dL, dR: o.RLGC.dR, scd: dB(scd21(sp)) };
        check(`${label}: one rough trace has L11 - L22 of several nH/m and a large conversion`,
            -o.RLGC.dL > 2e-9 && dB(scd21(sp)) > -30,
            `dL ${(o.RLGC.dL * 1e9).toFixed(2)} nH/m, dR ${o.RLGC.dR.toFixed(2)} ohm/m, Scd21 ${dB(scd21(sp)).toFixed(1)} dB`);
    }
}
{
    const q = rough.rectilinear, f = rough.triangular;
    check('rough trace: QS agrees with full-wave', rel(q.dL, f.dL) < 0.08 && rel(q.dR, f.dR) < 0.08 && Math.abs(q.scd - f.scd) < 1,
        `dL ${(q.dL * 1e9).toFixed(2)} vs ${(f.dL * 1e9).toFixed(2)} nH/m, Scd21 ${q.scd.toFixed(1)} vs ${f.scd.toFixed(1)} dB`);
}

console.log(failures === 0 ? '\nALL PAIR LINE ASYMMETRY TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
