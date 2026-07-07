// Eigensolver pencil unit test — no FEM, no mesh: drives solveGeneralized (the WASM
// shift-invert Arnoldi) on a sparse generalized pencil with a KNOWN spectrum and
// verifies every returned pair independently in JS.
//
// Construction: B is real tridiagonal SPD (2 / -0.5); A = B·Λ (column j of B scaled
// by λ_j), so A e_j = λ_j B e_j exactly — eigenvalues are diag(Λ), eigenvectors the
// standard basis. The spectrum mirrors the mode-viewer failure that motivated this
// test (regression 5a34928): a "quasi-TEM" eigenvalue right at the shift, a wanted
// cluster FAR from the shift (cavity modes ~3e6 away at γ²-scale magnitudes), and
// clouds of "spurious" filler beyond them. The solver must return ALL nev nearest
// pairs, not just the near-shift ones a single 20-vector Krylov pass converges —
// with the strict residual gate, that requires the ncvMax subspace growth.
//
// Run: node src/tri_solver/tests/eigen_pencil_test.mjs
import createModule from '../../wasm_solver/eigen_solver.js';
import { createWasmHelpers } from '../fem_core.js';

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}

const NPROB = 400;
const SIGMA = [-3.4e6, 0];          // like -k²·eps_static in the mode viewer
const NEV = 6, NCV = 20, NCV_MAX = 320;

// ---- pencil with known spectrum ----
function buildPencil(evals) {
    const N = evals.length;
    const rowPtr = new Int32Array(N + 1);
    const colIdx = [], bRe = [], aRe = [], aIm = [];
    for (let i = 0; i < N; i++) {
        for (const j of [i - 1, i, i + 1]) {
            if (j < 0 || j >= N) continue;
            const b = j === i ? 2 : -0.5;
            colIdx.push(j); bRe.push(b);
            aRe.push(b * evals[j].re); aIm.push(b * evals[j].im);
        }
        rowPtr[i + 1] = colIdx.length;
    }
    const nnz = colIdx.length;
    return {
        csrA: { rowPtr, colIdx: Int32Array.from(colIdx), valRe: Float64Array.from(aRe), valIm: Float64Array.from(aIm) },
        csrB: { rowPtr, colIdx: Int32Array.from(colIdx), valRe: Float64Array.from(bRe), valIm: new Float64Array(nnz) },
    };
}

// lossTan 0 → fully real pencil (exercises solve_real); >0 → complex (solve_complex),
// like the lossyEpsMap eigensolve the Modes tab always runs.
function makeSpectrum(lossTan) {
    const evals = new Array(NPROB);
    for (let i = 0; i < NPROB; i++) {
        const t = i / (NPROB - 1);
        // "spurious" filler clouds, all further from the shift than the wanted cluster
        const re = (i % 2 === 0) ? -9e6 - 1.1e7 * t : 2.5e6 * Math.pow(40, t);
        evals[i] = { re, im: -Math.abs(re) * lossTan };
    }
    // wanted cluster at the mode-viewer γ² scales: quasi-TEM near the shift, cavity
    // modes ~3e6 away, "evanescent" on the positive side — scattered indices
    const wantedRe = [-3.548e6, -5.447e5, -5.079e5, 4.425e4, 1.256e5, 1.110e6];
    const at = [37, 101, 150, 222, 289, 333];
    wantedRe.forEach((re, k) => { evals[at[k]] = { re, im: -Math.abs(re) * lossTan }; });
    return evals;
}

function nearestToShift(evals, k) {
    const d = (e) => Math.hypot(e.re - SIGMA[0], e.im - SIGMA[1]);
    return evals.map((e, i) => ({ ...e, i })).sort((p, q) => d(p) - d(q)).slice(0, k);
}

// ---- independent verification (complex CSR arithmetic in JS) ----
function csrMatvec(N, csr, xRe, xIm) {
    const yRe = new Float64Array(N), yIm = new Float64Array(N);
    for (let i = 0; i < N; i++) {
        let sr = 0, si = 0;
        for (let p = csr.rowPtr[i]; p < csr.rowPtr[i + 1]; p++) {
            const j = csr.colIdx[p], mr = csr.valRe[p], mi = csr.valIm[p];
            sr += mr * xRe[j] - mi * xIm[j];
            si += mr * xIm[j] + mi * xRe[j];
        }
        yRe[i] = sr; yIm[i] = si;
    }
    return { re: yRe, im: yIm };
}
const norm2 = (re, im) => { let s = 0; for (let i = 0; i < re.length; i++) s += re[i] ** 2 + im[i] ** 2; return Math.sqrt(s); };

// Same normalization as the residual gate in eigen_solver.cpp — recomputed here so
// the solver's own gate is not trusted to verify itself.
function relResidual(N, csrA, csrB, lamRe, lamIm, xRe, xIm) {
    const Ax = csrMatvec(N, csrA, xRe, xIm);
    const Bx = csrMatvec(N, csrB, xRe, xIm);
    const rRe = new Float64Array(N), rIm = new Float64Array(N);
    for (let i = 0; i < N; i++) {
        rRe[i] = Ax.re[i] - (lamRe * Bx.re[i] - lamIm * Bx.im[i]);
        rIm[i] = Ax.im[i] - (lamRe * Bx.im[i] + lamIm * Bx.re[i]);
    }
    return norm2(rRe, rIm) /
        (norm2(xRe, xIm) * (norm2(Ax.re, Ax.im) + Math.hypot(lamRe, lamIm) * norm2(Bx.re, Bx.im)));
}

function runCase(helpers, name, lossTan, useSeed, ncvMax) {
    const evals = makeSpectrum(lossTan);
    const { csrA, csrB } = buildPencil(evals);
    const wanted = nearestToShift(evals, NEV);
    let seed = null;
    if (useSeed) {
        // exact eigenvector of the eigenvalue nearest the shift — mirrors the
        // static-field seed the mode viewer passes (≈ the quasi-TEM eigenvector)
        seed = new Float64Array(NPROB);
        seed[wanted[0].i] = 1;
    }
    const res = helpers.solveGeneralized(NPROB, csrA, csrB, SIGMA, NEV, NCV, seed, ncvMax);

    check(`${name}: returns all ${NEV} requested pairs`, res.nconv === NEV, `nconv=${res.nconv}`);

    // every returned pair: eigenvalue matches a distinct known one; residual re-verified in JS
    const unmatched = wanted.slice();
    let maxRes = 0, maxEvalErr = 0, ok = true;
    for (let k = 0; k < res.nconv; k++) {
        const lr = res.evalsRe[k], li = res.evalsIm[k];
        const xRe = res.evecsRe.subarray(k * NPROB, (k + 1) * NPROB);
        const xIm = res.evecsIm.subarray(k * NPROB, (k + 1) * NPROB);
        const rr = relResidual(NPROB, csrA, csrB, lr, li, xRe, xIm);
        maxRes = Math.max(maxRes, rr);
        let best = -1, bestErr = Infinity;
        unmatched.forEach((w, u) => {
            const err = Math.hypot(lr - w.re, li - w.im) / Math.hypot(w.re, w.im);
            if (err < bestErr) { bestErr = err; best = u; }
        });
        maxEvalErr = Math.max(maxEvalErr, bestErr);
        if (bestErr < 1e-5) unmatched.splice(best, 1); else ok = false;
    }
    check(`${name}: every pair matches a distinct known eigenvalue`, ok && unmatched.length === NEV - res.nconv,
        `max rel eval err=${maxEvalErr.toExponential(1)}`);
    check(`${name}: JS-verified residual < 1e-4 for every pair`, maxRes < 1e-4,
        `max rel res=${maxRes.toExponential(1)}`);
    return res;
}

const M = await createModule();
const helpers = createWasmHelpers(M);

// real pencil → solve_real path; complex (lossy) pencil → solve_complex path,
// which is what the mode viewer always hits (lossyEpsMap floors tanδ at 0.002)
runCase(helpers, 'real, no seed', 0, false, NCV_MAX);
runCase(helpers, 'real, eigenvector seed', 0, true, NCV_MAX);
runCase(helpers, 'complex, no seed', 0.02, false, NCV_MAX);
runCase(helpers, 'complex, eigenvector seed', 0.02, true, NCV_MAX);

// Old fixed-subspace behavior (ncvMax=0): whatever comes back must STILL pass the
// independent residual gate — the count is informational (this is the configuration
// that regressed in 5a34928: far-from-shift pairs don't converge in one 20-vector pass).
{
    const evals = makeSpectrum(0.02);
    const { csrA, csrB } = buildPencil(evals);
    const res = helpers.solveGeneralized(NPROB, csrA, csrB, SIGMA, NEV, NCV, null, 0);
    let maxRes = 0;
    for (let k = 0; k < res.nconv; k++) {
        maxRes = Math.max(maxRes, relResidual(NPROB, csrA, csrB, res.evalsRe[k], res.evalsIm[k],
            res.evecsRe.subarray(k * NPROB, (k + 1) * NPROB), res.evecsIm.subarray(k * NPROB, (k + 1) * NPROB)));
    }
    console.log(`(info)  fixed 20-vector pass converges ${res.nconv}/${NEV} pairs — growth is what finds the rest`);
    check('fixed-subspace pass: returned pairs are genuinely converged', res.nconv === 0 || maxRes < 1e-4,
        `nconv=${res.nconv}, max rel res=${maxRes.toExponential(1)}`);
}

console.log(failures === 0 ? '\nEIGEN PENCIL OK' : `\nEIGEN PENCIL: ${failures} FAILURE(S)`);
process.exit(failures === 0 ? 0 : 1);
