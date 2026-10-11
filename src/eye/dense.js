// Small dense linear algebra for the circuit model: matrices are row-major
// Float64Arrays of n*n entries.

// In-place Cholesky factor of a symmetric positive definite matrix: returns L (lower
// triangular, row-major) with A = L L^T, or null when A is not positive definite.
function cholesky(A, n) {
    const L = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
        let d = A[j * n + j];
        for (let k = 0; k < j; k++) d -= L[j * n + k] * L[j * n + k];
        if (!(d > 0)) return null;
        const ljj = Math.sqrt(d);
        L[j * n + j] = ljj;
        for (let i = j + 1; i < n; i++) {
            let s = A[i * n + j];
            for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
            L[i * n + j] = s / ljj;
        }
    }
    return L;
}

// Solve L y = b (forward) for lower triangular L, b of length n. Returns y.
function forwardSub(L, n, b) {
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        let s = b[i];
        for (let k = 0; k < i; k++) s -= L[i * n + k] * y[k];
        y[i] = s / L[i * n + i];
    }
    return y;
}

// Solve L^T x = y (backward) for lower triangular L. Returns x.
function backwardSubT(L, n, y) {
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
        let s = y[i];
        for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k];
        x[i] = s / L[i * n + i];
    }
    return x;
}

// Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations.
// Returns { values, vectors } with vectors column-major in a row-major n*n array
// (column j, vectors[i * n + j], belongs to values[j]).
function jacobiEigen(S, n) {
    const a = Float64Array.from(S);
    const v = new Float64Array(n * n);
    for (let i = 0; i < n; i++) v[i * n + i] = 1;
    for (let sweep = 0; sweep < 100; sweep++) {
        let off = 0, diag = 0;
        for (let i = 0; i < n; i++) {
            diag += a[i * n + i] * a[i * n + i];
            for (let j = i + 1; j < n; j++) off += a[i * n + j] * a[i * n + j];
        }
        if (off <= 1e-32 * diag || off === 0) break;
        for (let p = 0; p < n; p++) {
            for (let q = p + 1; q < n; q++) {
                const apq = a[p * n + q];
                if (apq === 0) continue;
                const app = a[p * n + p], aqq = a[q * n + q];
                const theta = (aqq - app) / (2 * apq);
                const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
                const c = 1 / Math.sqrt(t * t + 1), s = t * c;
                for (let k = 0; k < n; k++) {
                    const akp = a[k * n + p], akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for (let k = 0; k < n; k++) {
                    const apk = a[p * n + k], aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for (let k = 0; k < n; k++) {
                    const vkp = v[k * n + p], vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    const values = new Float64Array(n);
    for (let i = 0; i < n; i++) values[i] = a[i * n + i];
    return { values, vectors: v };
}

// Solve the real system A x = b by Gaussian elimination with partial pivoting.
function solveReal(A, n, b) {
    const m = Float64Array.from(A), x = Float64Array.from(b);
    for (let c = 0; c < n; c++) {
        let p = c, best = Math.abs(m[c * n + c]);
        for (let r = c + 1; r < n; r++) {
            const v = Math.abs(m[r * n + c]);
            if (v > best) { best = v; p = r; }
        }
        if (best === 0) throw new Error('singular circuit matrix');
        if (p !== c) {
            for (let k = 0; k < n; k++) { const t = m[c * n + k]; m[c * n + k] = m[p * n + k]; m[p * n + k] = t; }
            const t = x[c]; x[c] = x[p]; x[p] = t;
        }
        const inv = 1 / m[c * n + c];
        for (let r = c + 1; r < n; r++) {
            const f = m[r * n + c] * inv;
            if (f === 0) continue;
            for (let k = c; k < n; k++) m[r * n + k] -= f * m[c * n + k];
            x[r] -= f * x[c];
        }
    }
    for (let r = n - 1; r >= 0; r--) {
        let s = x[r];
        for (let k = r + 1; k < n; k++) s -= m[r * n + k] * x[k];
        x[r] = s / m[r * n + r];
    }
    return x;
}

// Solve the complex system (Ar + i Ai) X = (Br + i Bi) for nrhs right-hand sides,
// in place. B is n x nrhs row-major; the solution overwrites Br, Bi.
function solveComplexInPlace(Ar, Ai, n, Br, Bi, nrhs) {
    for (let c = 0; c < n; c++) {
        let p = c, best = Ar[c * n + c] ** 2 + Ai[c * n + c] ** 2;
        for (let r = c + 1; r < n; r++) {
            const v = Ar[r * n + c] ** 2 + Ai[r * n + c] ** 2;
            if (v > best) { best = v; p = r; }
        }
        if (best === 0) throw new Error('singular circuit matrix');
        if (p !== c) {
            for (let k = 0; k < n; k++) {
                let t = Ar[c * n + k]; Ar[c * n + k] = Ar[p * n + k]; Ar[p * n + k] = t;
                t = Ai[c * n + k]; Ai[c * n + k] = Ai[p * n + k]; Ai[p * n + k] = t;
            }
            for (let k = 0; k < nrhs; k++) {
                let t = Br[c * nrhs + k]; Br[c * nrhs + k] = Br[p * nrhs + k]; Br[p * nrhs + k] = t;
                t = Bi[c * nrhs + k]; Bi[c * nrhs + k] = Bi[p * nrhs + k]; Bi[p * nrhs + k] = t;
            }
        }
        const pr = Ar[c * n + c], pi = Ai[c * n + c], pd = pr * pr + pi * pi;
        const ir = pr / pd, ii = -pi / pd;
        for (let r = c + 1; r < n; r++) {
            const xr = Ar[r * n + c], xi = Ai[r * n + c];
            if (xr === 0 && xi === 0) continue;
            const fr = xr * ir - xi * ii, fi = xr * ii + xi * ir;
            for (let k = c; k < n; k++) {
                const ur = Ar[c * n + k], ui = Ai[c * n + k];
                Ar[r * n + k] -= fr * ur - fi * ui;
                Ai[r * n + k] -= fr * ui + fi * ur;
            }
            for (let k = 0; k < nrhs; k++) {
                const ur = Br[c * nrhs + k], ui = Bi[c * nrhs + k];
                Br[r * nrhs + k] -= fr * ur - fi * ui;
                Bi[r * nrhs + k] -= fr * ui + fi * ur;
            }
        }
    }
    for (let r = n - 1; r >= 0; r--) {
        const pr = Ar[r * n + r], pi = Ai[r * n + r], pd = pr * pr + pi * pi;
        const ir = pr / pd, ii = -pi / pd;
        for (let k = 0; k < nrhs; k++) {
            let sr = Br[r * nrhs + k], si = Bi[r * nrhs + k];
            for (let c = r + 1; c < n; c++) {
                const ar = Ar[r * n + c], ai = Ai[r * n + c];
                const xr = Br[c * nrhs + k], xi = Bi[c * nrhs + k];
                sr -= ar * xr - ai * xi;
                si -= ar * xi + ai * xr;
            }
            Br[r * nrhs + k] = sr * ir - si * ii;
            Bi[r * nrhs + k] = sr * ii + si * ir;
        }
    }
}

export { cholesky, forwardSub, backwardSubT, jacobiEigen, solveReal, solveComplexInPlace };
