// fem_core.js — Shared FEM utilities for the triangular backend:
// WASM eigensolver/direct-solver helpers, COO→CSR, Jacobi-preconditioned CG,
// complex sqrt, shared Gauss quadrature constants, triangle quality metric.

// ==================== COO to CSR ====================

// Counting sort into per-row buckets + dense-scratch dedup + per-row insertion
// sort on the (small) unique column sets. O(nnz) overall vs the previous global
// comparator sort over an index array, which dominated whole-sweep CPU profiles
// (~44%) at FEM assembly sizes (millions of raw triplets).
export function tripletsToCSR(rows, cols, valsRe, N, valsIm) {
    const hasIm = valsIm != null;
    const nRaw = rows.length;
    // Pass 1: raw entries per row → bucket offsets.
    const rowStart = new Int32Array(N + 1);
    for (let i = 0; i < nRaw; i++) rowStart[rows[i] + 1]++;
    for (let r = 0; r < N; r++) rowStart[r + 1] += rowStart[r];
    // Pass 2: scatter triplets into their row buckets.
    const bCol = new Int32Array(nRaw);
    const bRe = new Float64Array(nRaw);
    const bIm = hasIm ? new Float64Array(nRaw) : null;
    const fill = rowStart.slice(0, N);
    for (let i = 0; i < nRaw; i++) {
        const p = fill[rows[i]]++;
        bCol[p] = cols[i]; bRe[p] = valsRe[i];
        if (hasIm) bIm[p] = valsIm[i];
    }
    // Pass 3: per row, sum duplicates via a dense col→slot scratch map (reset
    // per row by walking only the row's own entries), then insertion-sort the
    // unique columns (rows have tens of entries — insertion sort beats any
    // general sort here and keeps colIdx ascending for the WASM consumers).
    const outCol = new Int32Array(nRaw);
    const outRe = new Float64Array(nRaw);
    const outIm = new Float64Array(nRaw);   // stays zero when !hasIm
    const rowPtr = new Int32Array(N + 1);
    const colPos = new Int32Array(N).fill(-1);
    let nnz = 0;
    for (let r = 0; r < N; r++) {
        const rowOut = nnz;
        for (let p = rowStart[r]; p < rowStart[r + 1]; p++) {
            const c = bCol[p];
            let q = colPos[c];
            if (q >= 0) {
                outRe[q] += bRe[p];
                if (hasIm) outIm[q] += bIm[p];
            } else {
                q = nnz++;
                colPos[c] = q;
                outCol[q] = c; outRe[q] = bRe[p];
                if (hasIm) outIm[q] = bIm[p];
            }
        }
        for (let q = rowOut; q < nnz; q++) colPos[outCol[q]] = -1;
        for (let q = rowOut + 1; q < nnz; q++) {
            const c = outCol[q], vr = outRe[q], vi = outIm[q];
            let j = q - 1;
            while (j >= rowOut && outCol[j] > c) {
                outCol[j + 1] = outCol[j]; outRe[j + 1] = outRe[j]; outIm[j + 1] = outIm[j];
                j--;
            }
            outCol[j + 1] = c; outRe[j + 1] = vr; outIm[j + 1] = vi;
        }
        rowPtr[r + 1] = nnz;
    }
    return {
        rowPtr,
        colIdx: outCol.slice(0, nnz),
        valRe: outRe.slice(0, nnz),
        valIm: outIm.slice(0, nnz),
    };
}

// ==================== Linear Algebra ====================

function csrSpMV(csr, x, N) {
    const y = new Float64Array(N);
    for (let i = 0; i < N; i++) {
        let s = 0;
        for (let k = csr.rowPtr[i]; k < csr.rowPtr[i + 1]; k++)
            s += csr.valRe[k] * x[csr.colIdx[k]];
        y[i] = s;
    }
    return y;
}

export function dot(a, b) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }

export function solveCG(csr, rhs, N, tol = 1e-10, maxIter = 5000) {
    // Convergence is RELATIVE to ‖b‖ (an absolute floor is scale-fragile: the
    // Dirichlet-lift rhs magnitude depends on mesh scale and ε).
    const normB = Math.sqrt(dot(rhs, rhs));
    const stopTol = tol * Math.max(normB, 1e-300);
    // Jacobi (diagonal) preconditioner: M^{-1} where M = diag(A)
    const invDiag = new Float64Array(N);
    for (let i = 0; i < N; i++) {
        for (let k = csr.rowPtr[i]; k < csr.rowPtr[i + 1]; k++) {
            if (csr.colIdx[k] === i) { invDiag[i] = csr.valRe[k] > 1e-30 ? 1 / csr.valRe[k] : 1; break; }
        }
        if (invDiag[i] === 0) invDiag[i] = 1;
    }
    // Preconditioned CG: solve M^{-1}Ax = M^{-1}b
    const x = new Float64Array(N);
    const r = rhs.slice();
    const z = new Float64Array(N);
    for (let i = 0; i < N; i++) z[i] = invDiag[i] * r[i];
    const p = z.slice();
    let rz = dot(r, z);
    for (let iter = 0; iter < maxIter; iter++) {
        const Ap = csrSpMV(csr, p, N);
        const alpha = rz / dot(p, Ap);
        for (let i = 0; i < N; i++) x[i] += alpha * p[i];
        for (let i = 0; i < N; i++) r[i] -= alpha * Ap[i];
        const rnorm = Math.sqrt(dot(r, r));
        if (rnorm < stopTol) return { x, iters: iter + 1, residual: rnorm, converged: true };
        for (let i = 0; i < N; i++) z[i] = invDiag[i] * r[i];
        const rz_new = dot(r, z);
        const beta = rz_new / rz;
        for (let i = 0; i < N; i++) p[i] = z[i] + beta * p[i];
        rz = rz_new;
    }
    return { x, iters: maxIter, residual: Math.sqrt(dot(r, r)), converged: false };
}

export function csqrt(re, im) {
    const r = Math.sqrt(re * re + im * im);
    if (r === 0) return { re: 0, im: 0 };
    const theta = Math.atan2(im, re);
    const sr = Math.sqrt(r);
    return { re: sr * Math.cos(theta / 2), im: sr * Math.sin(theta / 2) };
}

// ==================== WASM Helpers ====================

// Debug hook (node only): EIGEN_WASM_DUMP=<dir> writes each WASM solver call's
// inputs to <dir>/last_<fn>.bin before the call and removes the file on success.
// After an eigen_assert abort() (which kills the WASM instance mid-call), the
// file left behind is the exact failing matrix, feed it to a native harness
// with assertions on to get a real stack trace. Zero cost when the env is unset.
const _dumpDir = (typeof process !== 'undefined' && process.env && process.env.EIGEN_WASM_DUMP) || null;
let _fs = null;
if (_dumpDir) import('fs').then(m => { _fs = m; });   // resolved long before any solve runs
function _dumpCall(fn, meta, arrays) {
    if (!_dumpDir || !_fs) return null;
    arrays = arrays.filter(([, a]) => a != null);
    const header = Buffer.from(JSON.stringify({ ...meta, arrays: arrays.map(([name, a]) => [name, a.constructor.name, a.length]) }));
    const parts = [Buffer.from(new Uint32Array([header.length]).buffer), header];
    for (const [, a] of arrays) parts.push(Buffer.from(a.buffer, a.byteOffset, a.byteLength));
    const path = `${_dumpDir}/last_${fn}.bin`;
    _fs.writeFileSync(path, Buffer.concat(parts));
    return path;
}
function _dumpDone(path) { if (path) try { _fs.unlinkSync(path); } catch { /* keep */ } }

export function createWasmHelpers(M) {
    function allocInt32(arr) {
        const p = M._malloc(4 * arr.length);
        new Int32Array(M.HEAP32.buffer, p, arr.length).set(arr);
        return p;
    }
    function allocFloat64(arr) {
        const p = M._malloc(8 * arr.length);
        new Float64Array(M.HEAPF64.buffer, p, arr.length).set(arr);
        return p;
    }
    function readFloat64(ptr, len) {
        return new Float64Array(M.HEAPF64.buffer, ptr, len).slice();
    }
    // Free every allocation even when the WASM call throws (the frees themselves
    // are guarded: after a genuine abort the module is unusable and _free itself
    // can throw, which must not mask the original error).
    function freeAll(ptrs) {
        for (const p of ptrs) { try { M._free(p); } catch { /* module aborted */ } }
    }
    // ncvMax > ncv lets the WASM solver GROW the Krylov subspace (doubling, reusing
    // the factorization) until nev pairs pass its residual gate or ncvMax columns
    // are reached — needed when wanted modes sit far from the shift (mode viewer).
    // The default (0 → clamped to ncv in C++) keeps the single fixed-size pass.
    function solveGeneralized(N, csrA, csrB, sigma, nev, ncv, initVec, ncvMax = 0) {
        const ptrs = [];
        function ai(a) { const p = allocInt32(a); ptrs.push(p); return p; }
        function af(a) { const p = allocFloat64(a); ptrs.push(p); return p; }
        const dumpPath = _dumpCall('solveGeneralized',
            { N, sigma, nev, ncv, ncvMax, hasInit: !!initVec },
            [['aRowPtr', csrA.rowPtr], ['aColIdx', csrA.colIdx], ['aValRe', csrA.valRe], ['aValIm', csrA.valIm],
             ['bRowPtr', csrB.rowPtr], ['bColIdx', csrB.colIdx], ['bValRe', csrB.valRe], ['bValIm', csrB.valIm],
             ...(initVec ? [['init', initVec]] : [])]);
        try {
            const pAr = ai(csrA.rowPtr), pAc = ai(csrA.colIdx), pAre = af(csrA.valRe), pAim = af(csrA.valIm);
            const pBr = ai(csrB.rowPtr), pBc = ai(csrB.colIdx), pBre = af(csrB.valRe), pBim = af(csrB.valIm);
            const pEvRe = M._malloc(8 * nev); ptrs.push(pEvRe);
            const pEvIm = M._malloc(8 * nev); ptrs.push(pEvIm);
            const pVRe = M._malloc(8 * nev * N); ptrs.push(pVRe);
            const pVIm = M._malloc(8 * nev * N); ptrs.push(pVIm);
            let nc;
            if (initVec) {
                const pInitRe = af(initVec);
                const pInitIm = af(new Float64Array(N));
                nc = M._solve_generalized_eigen_with_init(
                    N, csrA.colIdx.length, pAr, pAc, pAre, pAim,
                    csrB.colIdx.length, pBr, pBc, pBre, pBim,
                    sigma[0], sigma[1], nev, ncv, ncvMax, pEvRe, pEvIm, pVRe, pVIm,
                    pInitRe, pInitIm
                );
            } else {
                nc = M._solve_generalized_eigen(
                    N, csrA.colIdx.length, pAr, pAc, pAre, pAim,
                    csrB.colIdx.length, pBr, pBc, pBre, pBim,
                    sigma[0], sigma[1], nev, ncv, ncvMax, pEvRe, pEvIm, pVRe, pVIm
                );
            }
            // Negative nc = solver-reported failure (factorization failed,
            // exception caught in C++, …). Throw so callers can't mistake it for
            // a truthy "converged" count.
            if (nc < 0) throw new Error(`Eigensolver failed (code ${nc})`);
            _dumpDone(dumpPath);
            return {
                nconv: nc,
                evalsRe: readFloat64(pEvRe, nc > 0 ? nc : 0),
                evalsIm: readFloat64(pEvIm, nc > 0 ? nc : 0),
                evecsRe: readFloat64(pVRe, nc > 0 ? nc * N : 0),
                evecsIm: readFloat64(pVIm, nc > 0 ? nc * N : 0),
            };
        } finally {
            freeAll(ptrs);
        }
    }
    // Direct sparse solver: factorize once, solve for multiple RHS.
    // csr: {rowPtr, colIdx, valRe}, rhsArrays: array of Float64Array(N)
    // Returns array of Float64Array(N) solutions.
    function solveSparseMulti(N, csr, rhsArrays) {
        const ptrs = [];
        function ai(a) { const p = allocInt32(a); ptrs.push(p); return p; }
        function af(a) { const p = allocFloat64(a); ptrs.push(p); return p; }
        const dumpPath = _dumpCall('solveSparseMulti',
            { N, nRhs: rhsArrays.length },
            [['rowPtr', csr.rowPtr], ['colIdx', csr.colIdx], ['valRe', csr.valRe],
             ...rhsArrays.map((r, i) => [`rhs${i}`, r])]);
        try {
            const pR = ai(csr.rowPtr), pC = ai(csr.colIdx), pV = af(csr.valRe);
            const nRhs = rhsArrays.length;
            // Pack RHS into contiguous array
            const rhsPacked = new Float64Array(nRhs * N);
            for (let r = 0; r < nRhs; r++) rhsPacked.set(rhsArrays[r], r * N);
            const pRhs = af(rhsPacked);
            const pX = M._malloc(8 * nRhs * N); ptrs.push(pX);
            const rc = M._solve_sparse_multi(N, csr.colIdx.length, pR, pC, pV, nRhs, pRhs, pX);
            if (rc !== 0) throw new Error(`solve_sparse_multi failed: ${rc}`);
            _dumpDone(dumpPath);
            const results = [];
            for (let r = 0; r < nRhs; r++)
                results.push(readFloat64(pX + r * N * 8, N));
            return results;
        } finally {
            freeAll(ptrs);
        }
    }
    return { allocInt32, allocFloat64, readFloat64, solveGeneralized, solveSparseMulti };
}

// ==================== Shared geometric / quadrature constants ====================

// 3-point Gauss-Legendre rule on [0,1] (for edge / line-contour integrals).
// Points (1 ∓ √(3/5))/2, weights 5/18, 8/18, 5/18 (full double precision —
// 5-digit constants cost ~1e-5 relative error on every contour integral).
const GL3_A = Math.sqrt(3 / 5) / 2;
export const GL3p = [0.5 - GL3_A, 0.5, 0.5 + GL3_A];
export const GL3w = [5 / 18, 8 / 18, 5 / 18];

// Triangle quality metric (∝ circumradius / inradius; ≈1 for an equilateral triangle,
// larger is worse). Returns 1e10 for a degenerate (near-zero-area) triangle. Inputs are
// the three vertex coordinates.
export function triQuality(ax, ay, bx, by, cx, cy) {
    const al = Math.sqrt((bx-cx)**2+(by-cy)**2);
    const bl = Math.sqrt((ax-cx)**2+(ay-cy)**2);
    const cl = Math.sqrt((ax-bx)**2+(ay-by)**2);
    const s = (al+bl+cl)/2;
    const area = Math.abs((bx-ax)*(cy-ay)-(cx-ax)*(by-ay))/2;
    if (area < 1e-30) return 1e10;
    return al*bl*cl/(8*area*(area/s));
}
