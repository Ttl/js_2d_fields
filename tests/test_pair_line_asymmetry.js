// Per-line R and L matrices and mode conversion of a pair whose traces differ in metal,
// finish or width, on both backends.
//
// Both backends evaluate the conductor loss for the trace currents [1, 0], [0, 1] and
// [1, 1] (loss and internal inductance are quadratic forms in the currents) and hand the
// matrices to the 4-port S-parameter routine. A mirror-symmetric pair of two metals keeps
// C, G and the external L symmetric but has R11 != R22 and L11 != L22 (internal
// inductance); it carries the two differences dR, dL on the odd mode instead.
// The quasi-static R reads a few percent low at this mesh, and so do its differences.
//
// References: Ansys 2D Extractor, 2 GHz, causal substrate, 5.8e7 S/m unless noted.
//   right trace 3.8e7 S/m (the case of solve_differential_microstrip_unequal_sigma):
//     R [27.429, 1.4922, 1.4922, 32.306] ohm/m, L [291.58, 19.723, 19.723, 291.96] nH/m
//   right trace 2.3e6 S/m (copper / 25, an asymmetry of the size a rough trace causes):
//     R [27.847, 1.4757, 1.4757, 113.11] ohm/m, L [291.61, 19.725, 19.725, 298.16] nH/m
//   both: G [25.039, -0.10192, -0.10192, 25.038] mS/m, C [116.07, -2.6636, -2.6636, 116.07] pF/m
//   0.5 mm right trace:
//     R [27.39, 1.2808, 1.2808, 21.507] ohm/m, L [291.41, 16.834, 16.834, 240.81] nH/m,
//     G [25.039, -0.10156, -0.10156, 32.001] mS/m, C [116.07, -2.7478, -2.7478, 145.43] pF/m
//     G is the dielectric loss form of the per-trace solves (the eigenvector transform
//     of the modal G read G12 = -0.077).
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { computeSParamsDiffAuto, computeSParamsDifferential, computeSParamsDifferentialMTL, buildPhysicalRLGC } from '../src/sparameters.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';
import { Complex } from '../src/complex.js';
import { check, quiet, rel, done } from './helpers.js';

const pairText = (left, right, wRight) => `units mm
bounds gnd gnd gnd gnd
domain -1.5 1.5 0 1.2604
diel x=-inf w=inf y=0 h=0.2104 er=4.4 tand=0.02
sig+ x=-0.6 w=0.35 y=0.2104 h=0.05 ${left}
sig- x=0.25 w=${wRight} y=0.2104 h=0.05 ${right}
`;
function build(backend, left, right, wRight = 0.35) {
    const s = new CustomGeometrySolver({ text: pairText(left, right, wRight), sigma_cond: 5.8e7, freq: 2e9, nx: 10, ny: 10,
        mesh_backend: backend });
    s.use_causal_materials = true;
    return s;
}
const SOLVE = { energy_tol: 0.002 };
const solve = s => quiet(() => s.solve_adaptive(SOLVE));
const modesOf = r => [r.modes.find(m => m.mode === 'odd'), r.modes.find(m => m.mode === 'even')];
// Differential in, common out over a line: Scd21 = (S31 - S32 + S41 - S42) / 2.
function scd21(sp) {
    const S = sp.S;
    return S[2][0].sub(S[2][1]).add(S[3][0]).sub(S[3][1]).mul(new Complex(0.5, 0));
}
const dB = z => 20 * Math.log10(Math.max(z.abs(), 1e-15));
const LEN = 0.1, F = 2e9;
// Worst relative difference of pick(result) between interpolated and exact sweep points.
async function sweepWorst(solver, res, fs, pick) {
    const sweep = new InterpolatingSweep(solver, res, { tolerance: 0.005 });
    await quiet(() => sweep.run(0.5e9, 10e9, {}));
    const interp = sweep.buildResults(fs);
    let worst = 0;
    for (let i = 0; i < fs.length; i++) {
        const a = pick(await quiet(() => solver.computeAtFrequency(fs[i], res))), b = pick(interp[i].result);
        for (let k = 0; k < a.length; k++) worst = Math.max(worst, rel(a[k], b[k]));
    }
    return { worst, interp };
}

const BACKENDS = [['full-wave', 'triangular'], ['QS', 'rectilinear']];
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

// --- Two metals: the difference path, and the matrix path forced onto it ---
for (const [label, backend] of BACKENDS) {
    const T = TOL[backend];
    let first = null;
    for (const ref of REFS) {
        const solver = build(backend, '', ref.right);
        const res = await solve(solver);
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
        check(`${label}, ${ref.name}: no warning about missing mode conversion`,
            !(res.warnings || []).some(w => w.type === 'line-asymmetry'));
        const refS = computeSParamsDifferentialMTL(F, ref.R, ref.L, G_REF, C_REF, LEN, 50);
        const ourS = computeSParamsDiffAuto(F, odd.RLGC, even.RLGC, res.physMatrix, LEN, 50);
        check(`${label}, ${ref.name}: mode conversion over 100 mm matches the reference matrices`,
            Math.abs(dB(scd21(ourS)) - dB(scd21(refS))) < T.scd,
            `Scd21 ${dB(scd21(ourS)).toFixed(2)} vs ${dB(scd21(refS)).toFixed(2)} dB`);
    }

    // The matrix path, forced onto the same pair, reproduces the reference.
    {
        const ref = REFS[0];
        globalThis.__MODAL_FORCE__ = 'on';
        const r = await solve(build(backend, '', ref.right));
        globalThis.__MODAL_FORCE__ = undefined;
        const M = r.RLGC_matrix;
        check(`${label}: forced matrix path gives the reference R11, R22 and R12`, !!modesOf(r)[0].RLGC.Rm
            && rel(M.R[0][0], ref.R[0][0]) < T.tolR && rel(M.R[1][1], ref.R[1][1]) < T.tolR && rel(M.R[0][1], ref.R[0][1]) < T.tolR12,
            `${M.R[0][0].toFixed(3)} / ${M.R[1][1].toFixed(3)} / ${M.R[0][1].toFixed(4)}`);
        check(`${label}: and L with the internal inductance`, rel(M.L[0][0], ref.L[0][0]) < 5e-3 && rel(M.L[1][1], ref.L[1][1]) < 5e-3
            && rel(M.L[0][0] - M.L[1][1], ref.L[0][0] - ref.L[1][1]) < 0.1,
            `${(M.L[0][0] * 1e9).toFixed(2)} / ${(M.L[1][1] * 1e9).toFixed(2)} nH/m`);
    }

    // The perturbation loss has no per-line data and says so.
    if (backend === 'triangular') {
        const s = build(backend, '', 'sigma=3.8e7');
        s.tri_opts = { lossMethod: 'perturbation' };
        const r = await solve(s);
        check(`${label}: the perturbation loss method warns that mode conversion is missing`,
            (r.warnings || []).some(w => w.type === 'line-asymmetry') && modesOf(r)[0].RLGC.dR === undefined);
    }

    // Interpolating sweep carries the asymmetry.
    {
        const { worst, interp } = await sweepWorst(first.solver, first.res, [0.7e9, 1.3e9, 3.1e9, 7.7e9],
            r => { const m = modesOf(r)[0].RLGC; return [m.dR, m.dL]; });
        check(`${label}: interpolated sweep points carry R11 - R22 and L11 - L22`, worst < 0.03, `worst rel ${worst.toExponential(2)}`);
        const mat = interp[2].result.RLGC_matrix;
        check(`${label}: and their matrix has unequal diagonals`, mat.R[1][1] > mat.R[0][0] * 1.1);
    }

    // Equal traces: nothing changes.
    {
        const r = await solve(build(backend, 'sigma=3.8e7', 'sigma=3.8e7'));
        const [o, e] = modesOf(r);
        const a = computeSParamsDiffAuto(F, o.RLGC, e.RLGC, r.physMatrix, LEN, 50);
        const b = computeSParamsDifferential(F, o.RLGC, e.RLGC, LEN, 50);
        check(`${label}: equal traces carry no asymmetry and keep the odd/even S-parameters`,
            o.RLGC.dR === undefined && a.S[2][0].re === b.S[2][0].re && a.S[3][0].im === b.S[3][0].im);
    }

    // Swapping the traces mirrors the result.
    {
        const [o] = modesOf(await solve(build(backend, 'sigma=3.8e7', '')));
        check(`${label}: swapping the two metals flips the sign of the asymmetry`,
            rel(o.RLGC.dR, -first.odd.RLGC.dR) < 0.01 && rel(o.RLGC.dL, -first.odd.RLGC.dL) < 0.02,
            `dR ${o.RLGC.dR.toFixed(3)} vs ${(-first.odd.RLGC.dR).toFixed(3)}`);
    }

    // A rough trace converts far more, through the internal inductance.
    {
        const r = await solve(build(backend, '', 'rq=0.001'));
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

// --- Unequal widths: the matrices ---
const mats = {};
for (const [label, backend] of BACKENDS) {
    const s = build(backend, '', '', 0.2);
    const r = await solve(s);
    const [odd, even] = modesOf(r);
    const M = r.RLGC_matrix;
    mats[backend] = { M, s, r };
    check(`${label}: an asymmetric pair carries the matrices`, !!r.physMatrix && !!odd.RLGC.Rm && M.R[1][1] > M.R[0][0] * 1.2,
        `R ${M.R.flat().map(v => v.toFixed(2)).join(' ')}`);
    check(`${label}: L is the external matrix plus a positive internal part`,
        M.L[0][0] > r.physMatrix.L[0][0] * 1.002 && M.L[1][1] > r.physMatrix.L[1][1] * 1.002 && M.L[0][0] < r.physMatrix.L[0][0] * 1.03);
    if (backend === 'triangular') {
        // Mode currents i = C v, normalised to |i|^2 = 2, as the mode loss solve uses them.
        for (const m of [odd, even]) {
            const v = s._triBackend._static[m.mode].modalVec, C = r.physMatrix.C;
            let i = [C[0][0] * v[0] + C[0][1] * v[1], C[1][0] * v[0] + C[1][1] * v[1]];
            const n = Math.sqrt(2 / (i[0] ** 2 + i[1] ** 2)); i = i.map(x => x * n);
            const q = (i[0] * i[0] * M.R[0][0] + 2 * i[0] * i[1] * M.R[0][1] + i[1] * i[1] * M.R[1][1]) / 2;
            check(`${label}: ${m.mode}-mode R is the quadratic form of the matrix in the modal currents`, rel(q, m.RLGC.R) < 2e-3,
                `${m.RLGC.R.toFixed(3)} vs ${q.toFixed(3)} ohm/m`);
        }
        // The eigenvector transform it replaces is visibly off.
        const bare = { ...odd.RLGC }; delete bare.Rm; delete bare.Lim;
        const old = buildPhysicalRLGC(bare, even.RLGC, r.physMatrix);
        check(`${label}: the eigenvector transform of the modal R differs from it`, rel(old.R[0][1], M.R[0][1]) > 0.1,
            `R12 ${old.R[0][1].toFixed(3)} transformed vs ${M.R[0][1].toFixed(3)}`);
    }
}
{
    const q = mats.rectilinear.M, f = mats.triangular.M;
    check('QS agrees with full-wave on R11, R22', rel(q.R[0][0], f.R[0][0]) < 0.07 && rel(q.R[1][1], f.R[1][1]) < 0.07,
        `${q.R[0][0].toFixed(2)} / ${q.R[1][1].toFixed(2)} vs ${f.R[0][0].toFixed(2)} / ${f.R[1][1].toFixed(2)}`);
    check('QS agrees with full-wave on R12 and L', rel(q.R[0][1], f.R[0][1]) < 0.2 && rel(q.L[0][0], f.L[0][0]) < 5e-3
        && rel(q.L[1][1], f.L[1][1]) < 5e-3 && rel(q.L[0][1], f.L[0][1]) < 0.01,
        `R12 ${q.R[0][1].toFixed(3)} vs ${f.R[0][1].toFixed(3)}`);
    const { s, r } = mats.rectilinear;
    const { worst } = await sweepWorst(s, r, [0.9e9, 4.3e9],
        res => ['R', 'L'].flatMap(k => [[0, 0], [0, 1], [1, 1]].map(([p, q]) => res.RLGC_matrix[k][p][q])));
    check('interpolated sweep points carry the matrices', worst < 0.01, `worst rel ${worst.toExponential(2)}`);
}
{
    const REF = { R: [27.39, 1.2808, 21.507], L: [291.41e-9, 16.834e-9, 240.81e-9],
        G: [25.039e-3, -0.10156e-3, 32.001e-3], C: [116.07e-12, -2.7478e-12, 145.43e-12] };
    // [diagonal, off-diagonal] tolerance per matrix.
    const TOLS = { triangular: { R: [0.01, 0.02], L: [2e-3, 5e-3], G: [0.01, 0.05], C: [5e-3, 0.02] },
                   rectilinear: { R: [0.05, 0.15], L: [5e-3, 0.01], G: [0.01, 0.06], C: [6e-3, 0.02] } };
    for (const [label, backend] of BACKENDS) {
        const M = (await solve(build(backend, '', '', 0.5))).RLGC_matrix;
        for (const k of ['R', 'L', 'G', 'C']) {
            const got = [M[k][0][0], M[k][0][1], M[k][1][1]], [td, to] = TOLS[backend][k];
            const err = got.map((v, i) => rel(v, REF[k][i]));
            check(`${label}: 0.5 mm right trace, ${k} matrix matches the reference`, err[0] < td && err[2] < td && err[1] < to,
                `${got.map(v => v.toPrecision(5)).join(' ')} vs ${REF[k].map(v => v.toPrecision(5)).join(' ')}`);
        }
    }
}

done();
