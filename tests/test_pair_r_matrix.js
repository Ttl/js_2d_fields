// Per-line R and internal-L matrices of a geometrically asymmetric pair.
//
// Both backends evaluate the conductor loss for the trace currents [1, 0], [0, 1] and
// [1, 1] (loss and internal inductance are quadratic forms in the currents) and hand the
// matrices to the 4-port S-parameter routine, in place of the modal R transformed with
// the voltage eigenvectors and an L matrix without the internal inductance.
//
//   1. the matrix path, forced onto a mirror-symmetric pair of two metals, reproduces
//      the Ansys 2D Extractor matrices of test_pair_line_asymmetry
//   2. unequal trace widths: the full-wave mode R equals the quadratic form of the
//      matrix in the modal currents, QS agrees with full-wave, L carries the internal part
//   3. the interpolating sweep carries the matrices
//   4. Ansys 2D Extractor reference with a 0.5 mm right trace, both 5.8e7 S/m, 2 GHz,
//      causal substrate:
//        R [27.39, 1.2808, 1.2808, 21.507] ohm/m, L [291.41, 16.834, 16.834, 240.81] nH/m,
//        G [25.039, -0.10156, -0.10156, 32.001] mS/m, C [116.07, -2.7478, -2.7478, 145.43] pF/m
//      G is the dielectric loss form of the per-trace solves (the eigenvector transform
//      of the modal G read G12 = -0.077).
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { buildPhysicalRLGC } from '../src/sparameters.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';

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

const pairText = (wRight, right) => `units mm
bounds gnd gnd gnd gnd
domain -1.5 1.5 0 1.2604
diel x=-inf w=inf y=0 h=0.2104 er=4.4 tand=0.02
sig+ x=-0.6 w=0.35 y=0.2104 h=0.05
sig- x=0.25 w=${wRight} y=0.2104 h=0.05 ${right}
`;
function build(backend, wRight, right = '') {
    const s = new CustomGeometrySolver({ text: pairText(wRight, right), sigma_cond: 5.8e7, freq: 2e9, nx: 10, ny: 10,
        mesh_backend: backend });
    if (backend === 'triangular') s.tri_opts = { lossMethod: 'auto' };
    s.use_causal_materials = true;
    return s;
}
const SOLVE = { energy_tol: 0.002 };
const modesOf = r => [r.modes.find(m => m.mode === 'odd'), r.modes.find(m => m.mode === 'even')];
const BACKENDS = [['full-wave', 'triangular', 0.01], ['QS', 'rectilinear', 0.06]];

// --- 1. Forced matrix path against the reference ---
globalThis.__MODAL_FORCE__ = 'on';
for (const [label, backend, tol] of BACKENDS) {
    const r = await quiet(() => build(backend, 0.35, 'sigma=3.8e7').solve_adaptive(SOLVE));
    const M = r.RLGC_matrix;
    check(`${label}: forced matrix path gives the reference R11, R22`, !!modesOf(r)[0].RLGC.Rm
        && rel(M.R[0][0], 27.429) < tol && rel(M.R[1][1], 32.306) < tol,
        `${M.R[0][0].toFixed(3)} / ${M.R[1][1].toFixed(3)} vs 27.429 / 32.306 ohm/m`);
    check(`${label}: and the reference R12`, rel(M.R[0][1], 1.4922) < (backend === 'triangular' ? 0.03 : 0.2), `${M.R[0][1].toFixed(4)} vs 1.4922`);
    check(`${label}: L includes the internal inductance`, rel(M.L[0][0], 291.58e-9) < 5e-3 && rel(M.L[1][1], 291.96e-9) < 5e-3
        && rel(M.L[0][0] - M.L[1][1], -0.38e-9) < 0.1,
        `${(M.L[0][0] * 1e9).toFixed(2)} / ${(M.L[1][1] * 1e9).toFixed(2)} vs 291.58 / 291.96 nH/m`);
}
globalThis.__MODAL_FORCE__ = undefined;

// --- 2. Unequal widths ---
const mats = {};
for (const [label, backend] of BACKENDS) {
    const s = build(backend, 0.2);
    const r = await quiet(() => s.solve_adaptive(SOLVE));
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
}

// --- 3. Interpolating sweep ---
{
    const { s, r } = mats.rectilinear;
    const sweep = new InterpolatingSweep(s, r, { tolerance: 0.005 });
    await quiet(() => sweep.run(0.5e9, 10e9, {}));
    const fs = [0.9e9, 4.3e9];
    const interp = sweep.buildResults(fs);
    let worst = 0;
    for (let i = 0; i < fs.length; i++) {
        const exact = await quiet(() => s.computeAtFrequency(fs[i], r));
        const a = exact.RLGC_matrix, b = interp[i].result.RLGC_matrix;
        for (const k of ['R', 'L']) for (const [p, q] of [[0, 0], [0, 1], [1, 1]]) worst = Math.max(worst, rel(a[k][p][q], b[k][p][q]));
    }
    check('interpolated sweep points carry the matrices', worst < 0.01, `worst rel ${worst.toExponential(2)}`);
}

// --- 4. Unequal widths against the reference ---
{
    const REF = { R: [27.39, 1.2808, 21.507], L: [291.41e-9, 16.834e-9, 240.81e-9],
        G: [25.039e-3, -0.10156e-3, 32.001e-3], C: [116.07e-12, -2.7478e-12, 145.43e-12] };
    // [diagonal, off-diagonal] tolerance per matrix.
    const TOLS = { triangular: { R: [0.01, 0.02], L: [2e-3, 5e-3], G: [0.01, 0.05], C: [5e-3, 0.02] },
                   rectilinear: { R: [0.05, 0.15], L: [5e-3, 0.01], G: [0.01, 0.06], C: [6e-3, 0.02] } };
    for (const [label, backend] of BACKENDS) {
        const r = await quiet(() => build(backend, 0.5).solve_adaptive(SOLVE));
        const M = r.RLGC_matrix;
        for (const k of ['R', 'L', 'G', 'C']) {
            const got = [M[k][0][0], M[k][0][1], M[k][1][1]], [td, to] = TOLS[backend][k];
            const err = got.map((v, i) => rel(v, REF[k][i]));
            check(`${label}: ${k} matrix matches the reference`, err[0] < td && err[2] < td && err[1] < to,
                `${got.map(v => v.toPrecision(5)).join(' ')} vs ${REF[k].map(v => v.toPrecision(5)).join(' ')}`);
        }
    }
}

console.log(failures === 0 ? '\nALL PAIR R MATRIX TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
