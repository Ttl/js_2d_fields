// fem_core.js — Shared FEM module for p=1 and p=2 Nedelec edge elements
// Provides: WASM helpers, COO→CSR, CG solver, element matrices, Gauss quadrature,
// P2 Nedelec basis construction, Q2 nodal basis.

// ==================== Gauss-Legendre Quadrature ====================

export function gaussLegendre(n) {
    let pts, wts;
    if (n === 2) {
        const s = 1 / Math.sqrt(3);
        pts = [-s, s]; wts = [1, 1];
    } else if (n === 3) {
        const s = Math.sqrt(3 / 5);
        pts = [-s, 0, s]; wts = [5 / 9, 8 / 9, 5 / 9];
    } else if (n === 4) {
        const a = Math.sqrt(3 / 7 - 2 / 7 * Math.sqrt(6 / 5));
        const b = Math.sqrt(3 / 7 + 2 / 7 * Math.sqrt(6 / 5));
        const wa = (18 + Math.sqrt(30)) / 36;
        const wb = (18 - Math.sqrt(30)) / 36;
        pts = [-b, -a, a, b]; wts = [wb, wa, wa, wb];
    } else {
        throw new Error(`Gauss-Legendre not implemented for n=${n}`);
    }
    // Map [-1,1] → [0,1]
    return {
        points: pts.map(s => (s + 1) / 2),
        weights: wts.map(w => w / 2)
    };
}

// ==================== Matrix Utilities ====================

function invertMatrix(M, n) {
    // Gauss-Jordan elimination on n×n matrix (flat array)
    const A = new Array(n);
    for (let i = 0; i < n; i++) {
        A[i] = new Float64Array(2 * n);
        for (let j = 0; j < n; j++) A[i][j] = M[i * n + j];
        A[i][n + i] = 1;
    }
    for (let col = 0; col < n; col++) {
        let maxVal = Math.abs(A[col][col]), maxRow = col;
        for (let row = col + 1; row < n; row++) {
            if (Math.abs(A[row][col]) > maxVal) { maxVal = Math.abs(A[row][col]); maxRow = row; }
        }
        [A[col], A[maxRow]] = [A[maxRow], A[col]];
        const pivot = A[col][col];
        if (Math.abs(pivot) < 1e-15) throw new Error('Singular matrix');
        for (let j = 0; j < 2 * n; j++) A[col][j] /= pivot;
        for (let row = 0; row < n; row++) {
            if (row === col) continue;
            const f = A[row][col];
            for (let j = 0; j < 2 * n; j++) A[row][j] -= f * A[col][j];
        }
    }
    const inv = new Float64Array(n * n);
    for (let i = 0; i < n; i++)
        for (let j = 0; j < n; j++)
            inv[i * n + j] = A[i][n + j];
    return inv;
}

// ==================== P2 Nedelec Basis Construction ====================
//
// Second-order Nedelec first kind on reference element [0,1]²
// E_x ∈ Q_{1,2}: {1, ξ, η, ξη, η², ξη²} — degree ≤1 in ξ, ≤2 in η
// E_y ∈ Q_{2,1}: {1, ξ, ξ², η, ξη, ξ²η} — degree ≤2 in ξ, ≤1 in η
//
// 12 DOFs per element:
//   0,1: bottom h-edge P₀, P₁
//   2,3: top h-edge P₀, P₁
//   4,5: left v-edge P₀, P₁
//   6,7: right v-edge P₀, P₁
//   8,9: interior Ex₀, Ex₁
//   10,11: interior Ey₀, Ey₁
//
// DOFs 0-3,8,9 couple to E_x component; DOFs 4-7,10,11 couple to E_y component.

function buildP2Basis() {
    // Build Vandermonde matrices and invert to get basis coefficients.
    // E_x Vandermonde: V[dof][monomial], 6 DOFs × 6 monomials
    // E_x monomials on ref element: m0=1, m1=ξ, m2=η, m3=ξη, m4=η², m5=ξη²
    // DOF ordering: bottom P₀, bottom P₁, top P₀, top P₁, int Ex₀, int Ex₁
    const Vx = new Float64Array([
        // DOF 0 (bottom P₀): ∫₀¹ E_x(ξ,0) dξ
        1, 1 / 2, 0, 0, 0, 0,
        // DOF 1 (bottom P₁): ∫₀¹ E_x(ξ,0)(2ξ-1) dξ
        0, 1 / 6, 0, 0, 0, 0,
        // DOF 2 (top P₀): ∫₀¹ E_x(ξ,1) dξ
        1, 1 / 2, 1, 1 / 2, 1, 1 / 2,
        // DOF 3 (top P₁): ∫₀¹ E_x(ξ,1)(2ξ-1) dξ
        0, 1 / 6, 0, 1 / 6, 0, 1 / 6,
        // DOF 4 (int Ex₀): ∫∫ E_x·η(1-η) dξdη
        1 / 6, 1 / 12, 1 / 12, 1 / 24, 1 / 20, 1 / 40,
        // DOF 5 (int Ex₁): ∫∫ E_x·ξη(1-η) dξdη
        1 / 12, 1 / 18, 1 / 24, 1 / 36, 1 / 40, 1 / 60,
    ]);

    // E_y Vandermonde: 6 DOFs × 6 monomials
    // E_y monomials: m0=1, m1=ξ, m2=ξ², m3=η, m4=ξη, m5=ξ²η
    // DOF ordering: left P₀, left P₁, right P₀, right P₁, int Ey₀, int Ey₁
    const Vy = new Float64Array([
        // DOF 0 (left P₀): ∫₀¹ E_y(0,η) dη
        1, 0, 0, 1 / 2, 0, 0,
        // DOF 1 (left P₁): ∫₀¹ E_y(0,η)(2η-1) dη
        0, 0, 0, 1 / 6, 0, 0,
        // DOF 2 (right P₀): ∫₀¹ E_y(1,η) dη
        1, 1, 1, 1 / 2, 1 / 2, 1 / 2,
        // DOF 3 (right P₁): ∫₀¹ E_y(1,η)(2η-1) dη
        0, 0, 0, 1 / 6, 1 / 6, 1 / 6,
        // DOF 4 (int Ey₀): ∫∫ E_y·ξ(1-ξ) dξdη
        1 / 6, 1 / 12, 1 / 20, 1 / 12, 1 / 24, 1 / 40,
        // DOF 5 (int Ey₁): ∫∫ E_y·ηξ(1-ξ) dξdη
        1 / 12, 1 / 24, 1 / 40, 1 / 18, 1 / 36, 1 / 60,
    ]);

    // Basis coefficients C = (V^T)^{-1}: C[k][j] = coeff of monomial j in basis function k
    // σ_i(W_k) = Σ_j C[k][j] V[i][j] = (C V^T)[k][i] = δ_{ki}
    // So C V^T = I → C = (V^T)^{-1}
    const rawInvVx = invertMatrix(Vx, 6);
    const rawInvVy = invertMatrix(Vy, 6);
    // Transpose: (V^{-1})^T = (V^T)^{-1}
    const invVx = new Float64Array(36);
    const invVy = new Float64Array(36);
    for (let i = 0; i < 6; i++)
        for (let j = 0; j < 6; j++) {
            invVx[i * 6 + j] = rawInvVx[j * 6 + i];
            invVy[i * 6 + j] = rawInvVy[j * 6 + i];
        }
    return { invVx, invVy };
}

const { invVx, invVy } = buildP2Basis();

// Local DOF → (component, Vandermonde row) mapping
// DOFs 0,1,2,3,8,9 → E_x (Vandermonde rows 0-5)
// DOFs 4,5,6,7,10,11 → E_y (Vandermonde rows 0-5)
const p2DofComp = [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1]; // 0=x, 1=y
const p2DofVrow = [0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 4, 5];

// ==================== Precomputed Quadrature Data ====================

const GL4 = gaussLegendre(4);
const NQ = 16; // 4×4 quadrature points

// P2 basis at area quadrature points: Wx[q*12+i], Wy[q*12+i], curl[q*12+i]
export const p2Wx = new Float64Array(NQ * 12);
export const p2Wy = new Float64Array(NQ * 12);
export const p2Curl = new Float64Array(NQ * 12);
export const qw = new Float64Array(NQ);

// Q2 nodal basis at area quadrature points: N[q*9+k], dN/dxi[q*9+k], dN/deta[q*9+k]
export const q2N = new Float64Array(NQ * 9);
export const q2Dxi = new Float64Array(NQ * 9);
export const q2Deta = new Float64Array(NQ * 9);

// Q2 Lagrange 1D basis on {0, 1/2, 1}
function L0(t) { return 2 * t * t - 3 * t + 1; }
function L1(t) { return -4 * t * t + 4 * t; }
function L2(t) { return 2 * t * t - t; }
function dL0(t) { return 4 * t - 3; }
function dL1(t) { return -8 * t + 4; }
function dL2(t) { return 4 * t - 1; }

// Initialize precomputed data
(function initQuadratureData() {
    for (let qj = 0; qj < 4; qj++) {
        for (let qi = 0; qi < 4; qi++) {
            const q = qj * 4 + qi;
            const xi = GL4.points[qi], eta = GL4.points[qj];
            qw[q] = GL4.weights[qi] * GL4.weights[qj];

            // E_x monomials (Q_{1,2}): 1, ξ, η, ξη, η², ξη²
            const mx = [1, xi, eta, xi * eta, eta * eta, xi * eta * eta];
            const dmx_deta = [0, 0, 1, xi, 2 * eta, 2 * xi * eta];

            // E_y monomials (Q_{2,1}): 1, ξ, ξ², η, ξη, ξ²η
            const my = [1, xi, xi * xi, eta, xi * eta, xi * xi * eta];
            const dmy_dxi = [0, 1, 2 * xi, 0, eta, 2 * xi * eta];

            for (let i = 0; i < 12; i++) {
                const comp = p2DofComp[i], row = p2DofVrow[i];
                if (comp === 0) {
                    // E_x basis function
                    let wx = 0, c = 0;
                    for (let k = 0; k < 6; k++) {
                        wx += invVx[row * 6 + k] * mx[k];
                        c -= invVx[row * 6 + k] * dmx_deta[k];
                    }
                    p2Wx[q * 12 + i] = wx;
                    p2Wy[q * 12 + i] = 0;
                    p2Curl[q * 12 + i] = c; // curl = -∂Wx/∂η
                } else {
                    // E_y basis function
                    let wy = 0, c = 0;
                    for (let k = 0; k < 6; k++) {
                        wy += invVy[row * 6 + k] * my[k];
                        c += invVy[row * 6 + k] * dmy_dxi[k];
                    }
                    p2Wx[q * 12 + i] = 0;
                    p2Wy[q * 12 + i] = wy;
                    p2Curl[q * 12 + i] = c; // curl = ∂Wy/∂ξ
                }
            }

            // Q2 basis: N_k(ξ,η) = L_{ix}(ξ)·L_{jy}(η), k = jy*3 + ix
            const Lx = [L0(xi), L1(xi), L2(xi)];
            const Ly = [L0(eta), L1(eta), L2(eta)];
            const dLx = [dL0(xi), dL1(xi), dL2(xi)];
            const dLy = [dL0(eta), dL1(eta), dL2(eta)];
            for (let jy = 0; jy < 3; jy++) {
                for (let ix = 0; ix < 3; ix++) {
                    const k = jy * 3 + ix;
                    q2N[q * 9 + k] = Lx[ix] * Ly[jy];
                    q2Dxi[q * 9 + k] = dLx[ix] * Ly[jy];
                    q2Deta[q * 9 + k] = Lx[ix] * dLy[jy];
                }
            }
        }
    }
})();

// ==================== P2 Basis Evaluation (for edge integrals) ====================

export function evalP2Basis(xi, eta) {
    const Wx = new Float64Array(12), Wy = new Float64Array(12), curl = new Float64Array(12);
    const mx = [1, xi, eta, xi * eta, eta * eta, xi * eta * eta];
    const dmx_deta = [0, 0, 1, xi, 2 * eta, 2 * xi * eta];
    const my = [1, xi, xi * xi, eta, xi * eta, xi * xi * eta];
    const dmy_dxi = [0, 1, 2 * xi, 0, eta, 2 * xi * eta];
    for (let i = 0; i < 12; i++) {
        const comp = p2DofComp[i], row = p2DofVrow[i];
        if (comp === 0) {
            let wx = 0, c = 0;
            for (let k = 0; k < 6; k++) {
                wx += invVx[row * 6 + k] * mx[k];
                c -= invVx[row * 6 + k] * dmx_deta[k];
            }
            Wx[i] = wx; curl[i] = c;
        } else {
            let wy = 0, c = 0;
            for (let k = 0; k < 6; k++) {
                wy += invVy[row * 6 + k] * my[k];
                c += invVy[row * 6 + k] * dmy_dxi[k];
            }
            Wy[i] = wy; curl[i] = c;
        }
    }
    return { Wx, Wy, curl };
}

export function evalQ2Basis(xi, eta) {
    const N = new Float64Array(9), dxi = new Float64Array(9), deta = new Float64Array(9);
    const Lx = [L0(xi), L1(xi), L2(xi)];
    const Ly = [L0(eta), L1(eta), L2(eta)];
    const dLx = [dL0(xi), dL1(xi), dL2(xi)];
    const dLy = [dL0(eta), dL1(eta), dL2(eta)];
    for (let jy = 0; jy < 3; jy++) {
        for (let ix = 0; ix < 3; ix++) {
            const k = jy * 3 + ix;
            N[k] = Lx[ix] * Ly[jy];
            dxi[k] = dLx[ix] * Ly[jy];
            deta[k] = Lx[ix] * dLy[jy];
        }
    }
    return { N, dxi, deta };
}

// ==================== P2 Element Matrices ====================

export function computeElementMatricesP2(hx, hy) {
    const St = new Float64Array(144);  // 12×12 curl-curl
    const Mt = new Float64Array(144);  // 12×12 edge mass
    const Sz = new Float64Array(81);   // 9×9 nodal stiffness
    const Mz = new Float64Array(81);   // 9×9 nodal mass
    const G = new Float64Array(108);   // 12×9 gradient coupling

    const hyox = hy / hx, hxoy = hx / hy, hxhy = hx * hy;
    const invhxhy = 1 / hxhy;

    for (let q = 0; q < NQ; q++) {
        const w = qw[q];
        for (let i = 0; i < 12; i++) {
            const wxi = p2Wx[q * 12 + i], wyi = p2Wy[q * 12 + i], ci = p2Curl[q * 12 + i];
            for (let j = i; j < 12; j++) {
                const wxj = p2Wx[q * 12 + j], wyj = p2Wy[q * 12 + j], cj = p2Curl[q * 12 + j];
                const st = w * ci * cj * invhxhy;
                const mt = w * (hyox * wxi * wxj + hxoy * wyi * wyj);
                St[i * 12 + j] += st;
                Mt[i * 12 + j] += mt;
                if (j > i) { St[j * 12 + i] += st; Mt[j * 12 + i] += mt; }
            }
        }
        for (let i = 0; i < 9; i++) {
            const dxi_i = q2Dxi[q * 9 + i], deta_i = q2Deta[q * 9 + i], ni = q2N[q * 9 + i];
            for (let j = i; j < 9; j++) {
                const dxi_j = q2Dxi[q * 9 + j], deta_j = q2Deta[q * 9 + j], nj = q2N[q * 9 + j];
                const sz = w * (hyox * dxi_i * dxi_j + hxoy * deta_i * deta_j);
                const mz = w * hxhy * ni * nj;
                Sz[i * 9 + j] += sz;
                Mz[i * 9 + j] += mz;
                if (j > i) { Sz[j * 9 + i] += sz; Mz[j * 9 + i] += mz; }
            }
        }
        for (let i = 0; i < 12; i++) {
            const wxi = p2Wx[q * 12 + i], wyi = p2Wy[q * 12 + i];
            for (let j = 0; j < 9; j++) {
                G[i * 9 + j] += w * (hyox * wxi * q2Dxi[q * 9 + j] + hxoy * wyi * q2Deta[q * 9 + j]);
            }
        }
    }
    return { St, Mt, Sz, Mz, G };
}

// ==================== P1 Element Matrices (analytical) ====================

export function computeElementMatricesP1(hx, hy) {
    const signs = [1, -1, -1, 1];
    const St = new Array(16);
    for (let i = 0; i < 4; i++)
        for (let j = 0; j < 4; j++)
            St[i * 4 + j] = signs[i] * signs[j] / (hx * hy);

    const Mt = new Array(16).fill(0);
    Mt[0] = hy / (3 * hx); Mt[1] = hy / (6 * hx);
    Mt[4] = hy / (6 * hx); Mt[5] = hy / (3 * hx);
    Mt[10] = hx / (3 * hy); Mt[11] = hx / (6 * hy);
    Mt[14] = hx / (6 * hy); Mt[15] = hx / (3 * hy);

    const alpha = hy / (6 * hx), beta = hx / (6 * hy);
    const Sz = [
        2 * alpha + 2 * beta, -2 * alpha + beta, -alpha - beta, alpha - 2 * beta,
        -2 * alpha + beta, 2 * alpha + 2 * beta, alpha - 2 * beta, -alpha - beta,
        -alpha - beta, alpha - 2 * beta, 2 * alpha + 2 * beta, -2 * alpha + beta,
        alpha - 2 * beta, -alpha - beta, -2 * alpha + beta, 2 * alpha + 2 * beta
    ];

    const m = hx * hy / 36;
    const Mz = [4 * m, 2 * m, m, 2 * m, 2 * m, 4 * m, 2 * m, m, m, 2 * m, 4 * m, 2 * m, 2 * m, m, 2 * m, 4 * m];

    const ga = hy / hx, gb = hx / hy;
    const G = [
        -ga / 3, ga / 3, ga / 6, -ga / 6,
        -ga / 6, ga / 6, ga / 3, -ga / 3,
        -gb / 3, -gb / 6, gb / 6, gb / 3,
        -gb / 6, -gb / 3, gb / 3, gb / 6
    ];

    return { St, Mt, Sz, Mz, G };
}

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

export function csrSpMV(csr, x, N) {
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

export function csrSpMVcomplex(csr, x, N) {
    const yr = new Float64Array(N), yi = new Float64Array(N);
    for (let i = 0; i < N; i++) {
        let sr = 0, si = 0;
        for (let k = csr.rowPtr[i]; k < csr.rowPtr[i + 1]; k++) {
            sr += csr.valRe[k] * x[csr.colIdx[k]];
            si += csr.valIm[k] * x[csr.colIdx[k]];
        }
        yr[i] = sr; yi[i] = si;
    }
    return { re: yr, im: yi };
}

export function rayleighQuotient(csrA, csrB, u, N) {
    const Au = csrSpMVcomplex(csrA, u, N), Bu = csrSpMVcomplex(csrB, u, N);
    let numRe = 0, numIm = 0, denRe = 0, denIm = 0;
    for (let i = 0; i < N; i++) {
        numRe += u[i] * Au.re[i]; numIm += u[i] * Au.im[i];
        denRe += u[i] * Bu.re[i]; denIm += u[i] * Bu.im[i];
    }
    const d2 = denRe * denRe + denIm * denIm;
    return { re: (numRe * denRe + numIm * denIm) / d2, im: (numIm * denRe - numRe * denIm) / d2 };
}

export function csqrt(re, im) {
    const r = Math.sqrt(re * re + im * im);
    if (r === 0) return { re: 0, im: 0 };
    const theta = Math.atan2(im, re);
    const sr = Math.sqrt(r);
    return { re: sr * Math.cos(theta / 2), im: sr * Math.sin(theta / 2) };
}

// ==================== WASM Helpers ====================

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
    function solveGeneralized(N, csrA, csrB, sigma, nev, ncv, initVec) {
        const ptrs = [];
        function ai(a) { const p = allocInt32(a); ptrs.push(p); return p; }
        function af(a) { const p = allocFloat64(a); ptrs.push(p); return p; }
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
                    sigma[0], sigma[1], nev, ncv, pEvRe, pEvIm, pVRe, pVIm,
                    pInitRe, pInitIm
                );
            } else {
                nc = M._solve_generalized_eigen(
                    N, csrA.colIdx.length, pAr, pAc, pAre, pAim,
                    csrB.colIdx.length, pBr, pBc, pBre, pBim,
                    sigma[0], sigma[1], nev, ncv, pEvRe, pEvIm, pVRe, pVIm
                );
            }
            // Negative nc = solver-reported failure (factorization failed,
            // exception caught in C++, …). Throw so callers can't mistake it for
            // a truthy "converged" count.
            if (nc < 0) throw new Error(`Eigensolver failed (code ${nc})`);
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
            const results = [];
            for (let r = 0; r < nRhs; r++)
                results.push(readFloat64(pX + r * N * 8, N));
            return results;
        } finally {
            freeAll(ptrs);
        }
    }
    // Rayleigh quotient iteration: refine a known eigenvalue of complex Ax=λBx.
    // csrA, csrB: complex CSR matrices. sigma: [re, im] initial eigenvalue guess.
    // initVec: Float64Array(N) real initial eigenvector.
    // maxIter: number of RQI iterations (1-3 usually enough).
    // Returns { evalRe, evalIm, evecRe, evecIm } or throws on failure.
    function solveRQI(N, csrA, csrB, sigma, initVec, maxIter = 2) {
        const ptrs = [];
        function ai(a) { const p = allocInt32(a); ptrs.push(p); return p; }
        function af(a) { const p = allocFloat64(a); ptrs.push(p); return p; }
        try {
            const pAr = ai(csrA.rowPtr), pAc = ai(csrA.colIdx), pAre = af(csrA.valRe), pAim = af(csrA.valIm);
            const pBr = ai(csrB.rowPtr), pBc = ai(csrB.colIdx), pBre = af(csrB.valRe), pBim = af(csrB.valIm);
            const pInitRe = af(initVec);
            const pEvRe = M._malloc(8); ptrs.push(pEvRe);
            const pEvIm = M._malloc(8); ptrs.push(pEvIm);
            const pVRe = M._malloc(8 * N); ptrs.push(pVRe);
            const pVIm = M._malloc(8 * N); ptrs.push(pVIm);
            const rc = M._solve_rqi(
                N, csrA.colIdx.length, pAr, pAc, pAre, pAim,
                csrB.colIdx.length, pBr, pBc, pBre, pBim,
                sigma[0], sigma[1], maxIter, pInitRe, null,
                pEvRe, pEvIm, pVRe, pVIm
            );
            if (rc < 0) throw new Error(`solve_rqi failed: ${rc}`);
            return {
                evalRe: readFloat64(pEvRe, 1)[0],
                evalIm: readFloat64(pEvIm, 1)[0],
                evecRe: readFloat64(pVRe, N),
                evecIm: readFloat64(pVIm, N),
            };
        } finally {
            freeAll(ptrs);
        }
    }
    // QEP Rayleigh quotient iteration: refine eigenvalue of (A0+λA1+λ²A2)x=0.
    // csr0,csr1,csr2: complex CSR matrices for A0,A1,A2. sigma: [re,im] initial guess.
    // Returns { evalRe, evalIm, evecRe, evecIm } or throws.
    function solveQepRQI(N, csr0, csr1, csr2, sigma, initVec, initVecIm, maxIter = 2) {
        const ptrs = [];
        function ai(a) { const p = allocInt32(a); ptrs.push(p); return p; }
        function af(a) { const p = allocFloat64(a); ptrs.push(p); return p; }
        try {
            const p0r=ai(csr0.rowPtr),p0c=ai(csr0.colIdx),p0re=af(csr0.valRe),p0im=af(csr0.valIm||new Float64Array(csr0.valRe.length));
            const p1r=ai(csr1.rowPtr),p1c=ai(csr1.colIdx),p1re=af(csr1.valRe),p1im=af(csr1.valIm||new Float64Array(csr1.valRe.length));
            const p2r=ai(csr2.rowPtr),p2c=ai(csr2.colIdx),p2re=af(csr2.valRe),p2im=af(csr2.valIm||new Float64Array(csr2.valRe.length));
            const pIRe = af(initVec);
            const pIIm = initVecIm ? af(initVecIm) : null;
            const pEvRe = M._malloc(8); ptrs.push(pEvRe);
            const pEvIm = M._malloc(8); ptrs.push(pEvIm);
            const pVRe = M._malloc(8*N); ptrs.push(pVRe);
            const pVIm = M._malloc(8*N); ptrs.push(pVIm);
            const rc = M._solve_qep_rqi(
                N, csr0.colIdx.length, p0r, p0c, p0re, p0im,
                csr1.colIdx.length, p1r, p1c, p1re, p1im,
                csr2.colIdx.length, p2r, p2c, p2re, p2im,
                sigma[0], sigma[1], maxIter, pIRe, pIIm,
                pEvRe, pEvIm, pVRe, pVIm
            );
            if (rc < 0) throw new Error(`solve_qep_rqi failed: ${rc}`);
            return { evalRe: readFloat64(pEvRe,1)[0], evalIm: readFloat64(pEvIm,1)[0],
                     evecRe: readFloat64(pVRe,N), evecIm: readFloat64(pVIm,N) };
        } finally {
            freeAll(ptrs);
        }
    }
    return { allocInt32, allocFloat64, readFloat64, solveGeneralized, solveSparseMulti, solveRQI, solveQepRQI };
}

// ==================== P2 Basis Verification ====================

export function verifyP2Basis() {
    // Check σ_j(W_i) = δ_{ij} by evaluating DOF functionals on each basis function
    // Uses 4-point Gauss quadrature for edge integrals and area integrals
    const gl = gaussLegendre(4);
    const tol = 1e-10;
    let ok = true;

    for (let i = 0; i < 12; i++) {
        for (let j = 0; j < 12; j++) {
            let val = 0;
            if (j < 4) {
                // Edge DOF: h-edge P₀ or P₁
                const edgeEta = (j < 2) ? 0 : 1;
                const isP1 = (j % 2 === 1);
                for (let q = 0; q < 4; q++) {
                    const xi = gl.points[q], w = gl.weights[q];
                    const b = evalP2Basis(xi, edgeEta);
                    const test = isP1 ? (2 * xi - 1) : 1;
                    val += w * b.Wx[i] * test;
                }
            } else if (j < 8) {
                // Edge DOF: v-edge P₀ or P₁
                const edgeXi = (j < 6) ? 0 : 1;
                const isP1 = (j % 2 === 1);
                for (let q = 0; q < 4; q++) {
                    const eta = gl.points[q], w = gl.weights[q];
                    const b = evalP2Basis(edgeXi, eta);
                    const test = isP1 ? (2 * eta - 1) : 1;
                    val += w * b.Wy[i] * test;
                }
            } else {
                // Interior DOF
                const intIdx = j - 8; // 0,1,2,3
                for (let qj = 0; qj < 4; qj++) {
                    for (let qi = 0; qi < 4; qi++) {
                        const xi = gl.points[qi], eta = gl.points[qj];
                        const w = gl.weights[qi] * gl.weights[qj];
                        const b = evalP2Basis(xi, eta);
                        let test, field;
                        if (intIdx === 0) { test = eta * (1 - eta); field = b.Wx[i]; }
                        else if (intIdx === 1) { test = xi * eta * (1 - eta); field = b.Wx[i]; }
                        else if (intIdx === 2) { test = xi * (1 - xi); field = b.Wy[i]; }
                        else { test = eta * xi * (1 - xi); field = b.Wy[i]; }
                        val += w * field * test;
                    }
                }
            }
            const expected = (i === j) ? 1 : 0;
            if (Math.abs(val - expected) > tol) {
                console.error(`P2 basis check FAIL: σ_${j}(W_${i}) = ${val}, expected ${expected}`);
                ok = false;
            }
        }
    }

    // Check curl-of-gradient = 0: St * grad(φ) = 0 for Q2 nodal φ
    const hx = 1, hy = 1;
    const { St } = computeElementMatricesP2(hx, hy);
    // Gradient of Q2 basis function k → 12 edge DOFs
    for (let kn = 0; kn < 9; kn++) {
        // Compute grad DOFs for Q2 node kn
        const gradDofs = new Float64Array(12);
        for (let q = 0; q < NQ; q++) {
            const w = qw[q];
            const xi = GL4.points[q % 4], eta = GL4.points[Math.floor(q / 4)];
            const dphi_dxi = q2Dxi[q * 9 + kn];
            const dphi_deta = q2Deta[q * 9 + kn];
            // Interior DOFs via quadrature
            gradDofs[8] += w * (-dphi_dxi) * eta * (1 - eta);
            gradDofs[9] += w * (-dphi_dxi) * xi * eta * (1 - eta);
            gradDofs[10] += w * (-dphi_deta) * xi * (1 - xi);
            gradDofs[11] += w * (-dphi_deta) * eta * xi * (1 - xi);
        }
        // Edge DOFs from boundary values
        // Bottom (η=0): P₀ = -(N_k(1,0) - N_k(0,0)), P₁ = (2/3)(-N_k(0,0) + 2N_k(0.5,0) - N_k(1,0))
        const bN = evalQ2Basis(0, 0).N, bM = evalQ2Basis(0.5, 0).N, bR = evalQ2Basis(1, 0).N;
        gradDofs[0] = -(bR[kn] - bN[kn]);
        gradDofs[1] = (2 / 3) * (-bN[kn] + 2 * bM[kn] - bR[kn]);
        // Top (η=1)
        const tN = evalQ2Basis(0, 1).N, tM = evalQ2Basis(0.5, 1).N, tR = evalQ2Basis(1, 1).N;
        gradDofs[2] = -(tR[kn] - tN[kn]);
        gradDofs[3] = (2 / 3) * (-tN[kn] + 2 * tM[kn] - tR[kn]);
        // Left (ξ=0)
        const lB = evalQ2Basis(0, 0).N, lM = evalQ2Basis(0, 0.5).N, lT = evalQ2Basis(0, 1).N;
        gradDofs[4] = -(lT[kn] - lB[kn]);
        gradDofs[5] = (2 / 3) * (-lB[kn] + 2 * lM[kn] - lT[kn]);
        // Right (ξ=1)
        const rB = evalQ2Basis(1, 0).N, rM = evalQ2Basis(1, 0.5).N, rT = evalQ2Basis(1, 1).N;
        gradDofs[6] = -(rT[kn] - rB[kn]);
        gradDofs[7] = (2 / 3) * (-rB[kn] + 2 * rM[kn] - rT[kn]);

        // St * gradDofs should be zero
        let maxErr = 0;
        for (let i = 0; i < 12; i++) {
            let s = 0;
            for (let j = 0; j < 12; j++) s += St[i * 12 + j] * gradDofs[j];
            maxErr = Math.max(maxErr, Math.abs(s));
        }
        if (maxErr > 1e-10) {
            console.error(`P2 curl-grad check FAIL for node ${kn}: max |St·∇N| = ${maxErr.toExponential(3)}`);
            ok = false;
        }
    }

    // Check symmetry
    for (let i = 0; i < 12; i++)
        for (let j = i + 1; j < 12; j++) {
            if (Math.abs(St[i * 12 + j] - St[j * 12 + i]) > 1e-12)
                { console.error(`St not symmetric at (${i},${j})`); ok = false; }
            if (Math.abs(computeElementMatricesP2(1, 1).Mt[i * 12 + j] - computeElementMatricesP2(1, 1).Mt[j * 12 + i]) > 1e-12)
                { console.error(`Mt not symmetric at (${i},${j})`); ok = false; }
        }

    return ok;
}

// ==================== Q2 Robin ABC Mass Matrix (1D, on [0,1]) ====================
// 3×3 mass matrix for L₀, L₁, L₂ on [0,1]
export const q2Mass1D = new Float64Array([
    2 / 15, 1 / 15, -1 / 30,
    1 / 15, 8 / 15, 1 / 15,
    -1 / 30, 1 / 15, 2 / 15
]);

// --- Mode overlap: |<staticVec, eigenvec>| / (|staticVec| * |eigenvec|) ---
// Only over edge DOFs (the transverse E-field part)
export function modeOverlap(staticEdge, evecRe, evecIm, nFreeEdge) {
    let dotRe = 0, dotIm = 0, normS = 0, normE = 0;
    for (let i = 0; i < nFreeEdge; i++) {
        dotRe += staticEdge[i] * evecRe[i];
        dotIm += staticEdge[i] * evecIm[i];
        normS += staticEdge[i] * staticEdge[i];
        normE += evecRe[i] * evecRe[i] + evecIm[i] * evecIm[i];
    }
    return Math.sqrt(dotRe * dotRe + dotIm * dotIm) / (Math.sqrt(normS) * Math.sqrt(normE) + 1e-30);
}

// P2 edge Robin: ∫₀¹ |E_t|² dξ in DOF basis = [[1,0],[0,3]]
// (For E_t = d₀ + d₁(2ξ-1)/... → bilinear form is d₀² + 3d₁²)
// Scale by 1/h_perp for physical integral.

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
