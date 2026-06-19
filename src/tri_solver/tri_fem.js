// Triangular P2 Nedelec FEM — generic mesh, basis, assembly, and solver functions

import { tripletsToCSR, solveCG } from './fem_core.js';

// --- 6-point Gauss quadrature on triangle (barycentric) ---

export const QW = [0.22338159, 0.22338159, 0.22338159, 0.10995174, 0.10995174, 0.10995174];
export const QL1 = [0.10810302, 0.44594849, 0.44594849, 0.81684757, 0.09157621, 0.09157621];
export const QL2 = [0.44594849, 0.44594849, 0.10810302, 0.09157621, 0.09157621, 0.81684757];
export const QL3 = [0.44594849, 0.10810302, 0.44594849, 0.09157621, 0.81684757, 0.09157621];
export const NQ = 6;


// --- Barycentric coefficients ---
// Returns normalized coeff: coeff[i] = [a_i/(2A), b_i/(2A), c_i/(2A)]
// so that λ_i(x,y) = coeff[i][0] + coeff[i][1]*x + coeff[i][2]*y

export function triCoefficients(nodes, v0, v1, v2) {
    const x0 = nodes[2*v0], y0 = nodes[2*v0+1];
    const x1 = nodes[2*v1], y1 = nodes[2*v1+1];
    const x2 = nodes[2*v2], y2 = nodes[2*v2+1];

    let a1 = x1*y2 - y1*x2;
    let a2 = x2*y0 - y2*x0;
    let a3 = x0*y1 - y0*x1;
    let b1 = y1 - y2;
    let b2 = y2 - y0;
    let b3 = y0 - y1;
    let c1 = x2 - x1;
    let c2 = x0 - x2;
    let c3 = x1 - x0;

    const sA = 0.5 * ((x0 - x2) * (y1 - y0) - (x0 - x1) * (y2 - y0));
    const sign = sA >= 0 ? 1 : -1;
    const Area = Math.abs(sA);
    const twoA = 2 * Area;

    // Apply sign correction (same as EMerge tri_coefficients)
    a1 *= sign; a2 *= sign; a3 *= sign;
    b1 *= sign; b2 *= sign; b3 *= sign;
    c1 *= sign; c2 *= sign; c3 *= sign;

    // Normalized coefficients: coeff[vertex][component]
    // coeff[i] = [a_i/(2A), b_i/(2A), c_i/(2A)]
    const coeff = [
        [a1/twoA, b1/twoA, c1/twoA],
        [a2/twoA, b2/twoA, c2/twoA],
        [a3/twoA, b3/twoA, c3/twoA],
    ];

    return { coeff, Area, twoA };
}

// --- P2 Nedelec basis functions (using normalized coeff) ---
// All functions evaluate at a single point (x, y) and return [Wx, Wy] or scalar.
// The coeff is indexed as coeff[vertex_index][0=a, 1=b, 2=c] (normalized).

// _ne1: Whitney edge basis for edge (p, q)
// W_x = -b_p * λ_q + b_q * λ_p
// W_y = -c_p * λ_q + c_q * λ_p
// where λ_i = a_i + b_i*x + c_i*y (all normalized)
export function ne1(coeff, p, q, x, y) {
    const [a1, b1, c1] = coeff[p];
    const [a2, b2, c2] = coeff[q];
    const lam_q = a2 + b2*x + c2*y;
    const lam_p = a1 + b1*x + c1*y;
    return [
        -b1 * lam_q + b2 * lam_p,
        -c1 * lam_q + c2 * lam_p
    ];
}

// _ne2: quadratic edge basis for edge (p, q) — ne1 * (λ_p - λ_q)
export function ne2(coeff, p, q, x, y) {
    const [W0, W1] = ne1(coeff, p, q, x, y);
    const [a1, b1, c1] = coeff[p];
    const [a2, b2, c2] = coeff[q];
    const diff = (a1 - a2) + (b1 - b2)*x + (c1 - c2)*y;
    return [W0 * diff, W1 * diff];
}

// _nf1: face interior basis 1 (uses all 3 vertices)
export function nf1(coeff, x, y) {
    const [a1, b1, c1] = coeff[0];
    const [a2, b2, c2] = coeff[1];
    const [a3, b3, c3] = coeff[2];
    const lam1 = a1 + b1*x + c1*y;
    const lam2 = a2 + b2*x + c2*y;
    const lam3 = a3 + b3*x + c3*y;
    return [
        -(b1*lam3 - b3*lam1)*lam2 - (b2*lam3 - b3*lam2)*lam1,
        -(c1*lam3 - c3*lam1)*lam2 + (-c2*lam3 + c3*lam2)*lam1
    ];
}

// _nf2: face interior basis 2 (uses all 3 vertices)
export function nf2(coeff, x, y) {
    const [a1, b1, c1] = coeff[0];
    const [a2, b2, c2] = coeff[1];
    const [a3, b3, c3] = coeff[2];
    const lam1 = a1 + b1*x + c1*y;
    const lam2 = a2 + b2*x + c2*y;
    const lam3 = a3 + b3*x + c3*y;
    return [
        (b1*lam2 - b2*lam1)*lam3 + (b1*lam3 - b3*lam1)*lam2,
        -(-c1*lam2 + c2*lam1)*lam3 + (c1*lam3 - c3*lam1)*lam2
    ];
}

// _lv: quadratic vertex Lagrange basis for vertex i: -λ_i + 2*λ_i²
export function lv(coeff, i, x, y) {
    const [a1, b1, c1] = coeff[i];
    const lam = a1 + b1*x + c1*y;
    return -lam + 2*lam*lam;
}

// _le: quadratic edge Lagrange basis for edge (i,j): 4*λ_i*λ_j
export function le(coeff, i, j, x, y) {
    const [a1, b1, c1] = coeff[i];
    const [a2, b2, c2] = coeff[j];
    const lam_i = a1 + b1*x + c1*y;
    const lam_j = a2 + b2*x + c2*y;
    return 4*lam_i*lam_j;
}

// _lv_grad: gradient of quadratic vertex basis for vertex i
export function lvGrad(coeff, i, x, y) {
    const [a1, b1, c1] = coeff[i];
    const lam = a1 + b1*x + c1*y;
    return [b1*(4*lam - 1), c1*(4*lam - 1)];
}

// _le_grad: gradient of quadratic edge basis for edge (i,j)
export function leGrad(coeff, i, j, x, y) {
    const [a1, b1, c1] = coeff[i];
    const [a2, b2, c2] = coeff[j];
    const lam_i = a1 + b1*x + c1*y;
    const lam_j = a2 + b2*x + c2*y;
    return [4*b1*lam_j + 4*b2*lam_i, 4*c1*lam_j + 4*c2*lam_i];
}

// _ne1_curl: curl of Whitney edge basis (constant per element)
export function ne1Curl(coeff, p, q) {
    const [, b1, c1] = coeff[p];
    const [, b2, c2] = coeff[q];
    return 2*(b1*c2 - b2*c1);
}

// _ne2_curl: curl of quadratic edge basis
export function ne2Curl(coeff, p, q, x, y) {
    const [a1, b1, c1] = coeff[p];
    const [a2, b2, c2] = coeff[q];
    const lam_p = a1 + b1*x + c1*y;
    const lam_q = a2 + b2*x + c2*y;
    return -(b1 - b2)*(c1*lam_q - c2*lam_p) +
            (c1 - c2)*(b1*lam_q - b2*lam_p) +
            2*(b1*c2 - b2*c1)*(lam_p - lam_q);
}

// _nf1_curl: curl of face basis 1
export function nf1Curl(coeff, x, y) {
    const [a1, b1, c1] = coeff[0];
    const [a2, b2, c2] = coeff[1];
    const [a3, b3, c3] = coeff[2];
    return 3*a1*b2*c3 - 3*a1*b3*c2 + 3*a2*b1*c3 - 3*a2*b3*c1 +
           6*b1*b2*c3*x - 3*b1*b3*c2*x + 3*b1*c2*c3*y -
           3*b2*b3*c1*x + 3*b2*c1*c3*y - 6*b3*c1*c2*y;
}

// _nf2_curl: curl of face basis 2
export function nf2Curl(coeff, x, y) {
    const [a1, b1, c1] = coeff[0];
    const [a2, b2, c2] = coeff[1];
    const [a3, b3, c3] = coeff[2];
    return -3*a2*b1*c3 + 3*a2*b3*c1 - 3*a3*b1*c2 + 3*a3*b2*c1 -
           3*b1*b2*c3*x - 3*b1*b3*c2*x - 6*b1*c2*c3*y +
           6*b2*b3*c1*x + 3*b2*c1*c3*y + 3*b3*c1*c2*y;
}

// --- 12-point Dunavant quadrature on triangle (exact for degree 6) ---
// Used for P3 element matrices where integrands have degree up to 6.
export const QW12 = [
    0.05084490637, 0.05084490637, 0.05084490637,
    0.11678627572, 0.11678627572, 0.11678627572,
    0.08285107561, 0.08285107561, 0.08285107561,
    0.08285107561, 0.08285107561, 0.08285107561
];
// Orbit 1 (3 pts): (1-2a, a, a) permuted, a=0.063089
// Orbit 2 (3 pts): (1-2a, a, a) permuted, a=0.249287
// Orbit 3 (6 pts): all permutations of (a, b, c), a=0.310352, b=0.053145, c=0.636503
export const QL1_12 = [
    0.87382197101, 0.06308901449, 0.06308901449,
    0.50142650966, 0.24928674517, 0.24928674517,
    0.63650249912, 0.31035245103, 0.05314504984,
    0.63650249912, 0.05314504984, 0.31035245103
];
export const QL2_12 = [
    0.06308901449, 0.87382197101, 0.06308901449,
    0.24928674517, 0.50142650966, 0.24928674517,
    0.31035245103, 0.63650249912, 0.63650249912,
    0.05314504984, 0.31035245103, 0.05314504984
];
export const QL3_12 = [
    0.06308901449, 0.06308901449, 0.87382197101,
    0.24928674517, 0.24928674517, 0.50142650966,
    0.05314504984, 0.05314504984, 0.31035245103,
    0.31035245103, 0.63650249912, 0.63650249912
];
export const NQ12 = 12;

// --- P3 Nedelec basis functions (H(curl), hierarchical extension of P2) ---
// P3 adds 7 transverse DOFs: ne3 (3 edge) + nf3..nf6 (4 face)
// Total P3 transverse: 15 DOFs = 9 edge + 6 face
//
// Local DOF ordering for P3 transverse:
// [ne1_0, ne1_1, ne1_2, nf1, ne2_0, ne2_1, ne2_2, nf2,  <- P2 (8)
//  ne3_0, ne3_1, ne3_2, nf3, nf4, nf5, nf6]              <- P3 extra (7)

// _ne3: cubic edge basis for edge (p, q) — ne1 * (λ_p - λ_q)²
// Antisymmetric: ne3(q,p) = -ne3(p,q) — needs sign flip like ne1
export function ne3(coeff, p, q, x, y) {
    const [W0, W1] = ne1(coeff, p, q, x, y);
    const [a1, b1, c1] = coeff[p];
    const [a2, b2, c2] = coeff[q];
    const diff = (a1 - a2) + (b1 - b2)*x + (c1 - c2)*y;
    const d2 = diff * diff;
    return [W0 * d2, W1 * d2];
}

// _ne3_curl: curl of cubic edge basis
// curl(ne1·f²) = 2f·∇f × ne1 + f²·curl(ne1) where f = λ_p - λ_q
// Simplifies to: (λ_p - λ_q)·[(b_p-b_q)ne2_y - (c_p-c_q)ne2_x] + (λ_p-λ_q)·ne2Curl
export function ne3Curl(coeff, p, q, x, y) {
    const [a1, b1, c1] = coeff[p];
    const [a2, b2, c2] = coeff[q];
    const diff = (a1 - a2) + (b1 - b2)*x + (c1 - c2)*y;
    const [ne2x, ne2y] = ne2(coeff, p, q, x, y);
    const cross = (b1 - b2)*ne2y - (c1 - c2)*ne2x;
    return cross + diff * ne2Curl(coeff, p, q, x, y);
}

// _nf3..nf6: P3 face interior basis functions (4 new, hierarchical)
// Constructed as λ_i · nf_j, which vanishes tangentially on all edges.
// nf3 = λ_0 · nf1,  nf4 = λ_1 · nf1,  nf5 = λ_0 · nf2,  nf6 = λ_1 · nf2
export function nf3(coeff, x, y) {
    const lam0 = coeff[0][0] + coeff[0][1]*x + coeff[0][2]*y;
    const [fx, fy] = nf1(coeff, x, y);
    return [lam0 * fx, lam0 * fy];
}
export function nf4(coeff, x, y) {
    const lam1 = coeff[1][0] + coeff[1][1]*x + coeff[1][2]*y;
    const [fx, fy] = nf1(coeff, x, y);
    return [lam1 * fx, lam1 * fy];
}
export function nf5(coeff, x, y) {
    const lam0 = coeff[0][0] + coeff[0][1]*x + coeff[0][2]*y;
    const [fx, fy] = nf2(coeff, x, y);
    return [lam0 * fx, lam0 * fy];
}
export function nf6(coeff, x, y) {
    const lam1 = coeff[1][0] + coeff[1][1]*x + coeff[1][2]*y;
    const [fx, fy] = nf2(coeff, x, y);
    return [lam1 * fx, lam1 * fy];
}

// Curls of P3 face functions: curl(λ_i · F) = b_i·F_y - c_i·F_x + λ_i·curl(F)
export function nf3Curl(coeff, x, y) {
    const [, b0, c0] = coeff[0];
    const lam0 = coeff[0][0] + b0*x + c0*y;
    const [fx, fy] = nf1(coeff, x, y);
    return b0*fy - c0*fx + lam0*nf1Curl(coeff, x, y);
}
export function nf4Curl(coeff, x, y) {
    const [, b1, c1] = coeff[1];
    const lam1 = coeff[1][0] + b1*x + c1*y;
    const [fx, fy] = nf1(coeff, x, y);
    return b1*fy - c1*fx + lam1*nf1Curl(coeff, x, y);
}
export function nf5Curl(coeff, x, y) {
    const [, b0, c0] = coeff[0];
    const lam0 = coeff[0][0] + b0*x + c0*y;
    const [fx, fy] = nf2(coeff, x, y);
    return b0*fy - c0*fx + lam0*nf2Curl(coeff, x, y);
}
export function nf6Curl(coeff, x, y) {
    const [, b1, c1] = coeff[1];
    const lam1 = coeff[1][0] + b1*x + c1*y;
    const [fx, fy] = nf2(coeff, x, y);
    return b1*fy - c1*fx + lam1*nf2Curl(coeff, x, y);
}

// --- P3 Lagrange basis functions (H1, hierarchical extension of P2) ---
// P3 adds 4 longitudinal DOFs: le3 (3 edge extra) + lf3 (1 face bubble)
// Total P3 longitudinal: 10 DOFs = 3 vertex + 3 edge-mid + 3 edge-extra + 1 face
//
// Local DOF ordering for P3 longitudinal:
// [lv_0, lv_1, lv_2, le_0, le_1, le_2,   <- P2 (6)
//  le3_0, le3_1, le3_2, lf3]               <- P3 extra (4)
//
// Hierarchical: P2 vertex (lv) and edge (le) functions are reused unchanged.
// le3(i,j) = λ_i·λ_j·(λ_i - λ_j) — zero at vertices and edge midpoints.
// lf3 = 27·λ_0·λ_1·λ_2 — zero on all edges, = 1 at centroid.

export function le3(coeff, i, j, x, y) {
    const lam_i = coeff[i][0] + coeff[i][1]*x + coeff[i][2]*y;
    const lam_j = coeff[j][0] + coeff[j][1]*x + coeff[j][2]*y;
    return lam_i * lam_j * (lam_i - lam_j);
}

export function le3Grad(coeff, i, j, x, y) {
    const [ai, bi, ci] = coeff[i];
    const [aj, bj, cj] = coeff[j];
    const lam_i = ai + bi*x + ci*y;
    const lam_j = aj + bj*x + cj*y;
    // ∇(λ_i·λ_j·(λ_i-λ_j)) = ∇λ_i·λ_j(2λ_i-λ_j) + ∇λ_j·λ_i(λ_i-2λ_j)
    const fi = lam_j * (2*lam_i - lam_j);
    const fj = lam_i * (lam_i - 2*lam_j);
    return [bi*fi + bj*fj, ci*fi + cj*fj];
}

export function lf3(coeff, x, y) {
    const lam0 = coeff[0][0] + coeff[0][1]*x + coeff[0][2]*y;
    const lam1 = coeff[1][0] + coeff[1][1]*x + coeff[1][2]*y;
    const lam2 = coeff[2][0] + coeff[2][1]*x + coeff[2][2]*y;
    return 27 * lam0 * lam1 * lam2;
}

export function lf3Grad(coeff, x, y) {
    const [a0, b0, c0] = coeff[0];
    const [a1, b1, c1] = coeff[1];
    const [a2, b2, c2] = coeff[2];
    const lam0 = a0 + b0*x + c0*y;
    const lam1 = a1 + b1*x + c1*y;
    const lam2 = a2 + b2*x + c2*y;
    return [
        27*(b0*lam1*lam2 + lam0*b1*lam2 + lam0*lam1*b2),
        27*(c0*lam1*lam2 + lam0*c1*lam2 + lam0*lam1*c2)
    ];
}

// --- P2 element matrices via Gauss quadrature ---
// Returns 14x14 A and B matrices for the generalized eigenvalue problem
// DOF ordering: [ne1_0, ne1_1, ne1_2, nf1, ne2_0, ne2_1, ne2_2, nf2, lv_0, lv_1, lv_2, le_0, le_1, le_2]

export function computeTriP2Matrices(nodes, v0, v1, v2, epsRe, epsIm, k0) {
    const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);

    // Local edge map: edge k → [local vertex p, local vertex q]
    // Edge 0: (v0,v1) → [0,1], Edge 1: (v1,v2) → [1,2], Edge 2: (v2,v0) → [2,0]
    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    const k02 = k0 * k0;

    // Sub-matrices (real and imaginary parts)
    const Att = new Float64Array(64);   // 8x8
    const BttRe = new Float64Array(64); // 8x8
    const BttIm = new Float64Array(64);
    const Dtt = new Float64Array(64);   // 8x8
    const Dzt = new Float64Array(48);   // 6x8
    const Dzz1 = new Float64Array(36);  // 6x6
    const Dzz2Re = new Float64Array(36); // 6x6
    const Dzz2Im = new Float64Array(36);

    // Physical coordinates of triangle vertices
    const txs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
    const tys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];

    // Gauss quadrature loop
    for (let q = 0; q < NQ; q++) {
        const w = QW[q];
        const xq = txs[0]*QL1[q] + txs[1]*QL2[q] + txs[2]*QL3[q];
        const yq = tys[0]*QL1[q] + tys[1]*QL2[q] + tys[2]*QL3[q];

        // Evaluate all 8 transverse basis functions and their curls at (xq, yq)
        // Indices 0-2: ne1 for edges 0,1,2
        // Index 3: nf1
        // Indices 4-6: ne2 for edges 0,1,2
        // Index 7: nf2
        const Wx = new Float64Array(8);
        const Wy = new Float64Array(8);
        const curlW = new Float64Array(8);

        for (let k = 0; k < 3; k++) {
            const [p, qq] = edgeVerts[k];
            const [wx1, wy1] = ne1(coeff, p, qq, xq, yq);
            Wx[k] = wx1; Wy[k] = wy1;
            curlW[k] = ne1Curl(coeff, p, qq);

            const [wx2, wy2] = ne2(coeff, p, qq, xq, yq);
            Wx[k+4] = wx2; Wy[k+4] = wy2;
            curlW[k+4] = ne2Curl(coeff, p, qq, xq, yq);
        }
        {
            const [fx1, fy1] = nf1(coeff, xq, yq);
            Wx[3] = fx1; Wy[3] = fy1;
            curlW[3] = nf1Curl(coeff, xq, yq);

            const [fx2, fy2] = nf2(coeff, xq, yq);
            Wx[7] = fx2; Wy[7] = fy2;
            curlW[7] = nf2Curl(coeff, xq, yq);
        }

        // Evaluate all 6 longitudinal basis functions and their gradients
        // Indices 0-2: lv for vertices 0,1,2
        // Indices 3-5: le for edges 0,1,2
        const Lz = new Float64Array(6);
        const Gx = new Float64Array(6);
        const Gy = new Float64Array(6);

        for (let k = 0; k < 3; k++) {
            Lz[k] = lv(coeff, k, xq, yq);
            const [gx, gy] = lvGrad(coeff, k, xq, yq);
            Gx[k] = gx; Gy[k] = gy;
        }
        for (let k = 0; k < 3; k++) {
            const [p, qq] = edgeVerts[k];
            Lz[k+3] = le(coeff, p, qq, xq, yq);
            const [gx, gy] = leGrad(coeff, p, qq, xq, yq);
            Gx[k+3] = gx; Gy[k+3] = gy;
        }

        // Accumulate sub-matrices with quadrature weight
        // For isotropic materials: Ms = I (ur_inv = 1), Mm = eps * I
        // Msz = 1, Mmz = eps

        // Att[i,j] = ∫ curlW_i * Msz * curlW_j dA = ∫ curlW_i * curlW_j dA
        // Btt[i,j] = ∫ W_i · (Mm * W_j) dA = eps * ∫ W_i · W_j dA
        // Dtt[i,j] = ∫ W_i · (Ms * W_j) dA = ∫ W_i · W_j dA
        for (let i = 0; i < 8; i++) {
            for (let j = 0; j < 8; j++) {
                const idx = i * 8 + j;
                Att[idx] += w * curlW[i] * curlW[j];
                const dot = Wx[i]*Wx[j] + Wy[i]*Wy[j];
                BttRe[idx] += w * epsRe * dot;
                BttIm[idx] += w * epsIm * dot;
                Dtt[idx] += w * dot;
            }
        }

        // Dzt[i,j] = ∫ ∇Lz_i · (Ms * W_j) dA = ∫ ∇Lz_i · W_j dA
        for (let i = 0; i < 6; i++) {
            for (let j = 0; j < 8; j++) {
                Dzt[i * 8 + j] += w * (Gx[i]*Wx[j] + Gy[i]*Wy[j]);
            }
        }

        // Dzz1[i,j] = ∫ ∇Lz_i · (Ms * ∇Lz_j) dA = ∫ ∇Lz_i · ∇Lz_j dA
        // Dzz2[i,j] = ∫ Lz_i * Mmz * Lz_j dA = eps * ∫ Lz_i * Lz_j dA
        for (let i = 0; i < 6; i++) {
            for (let j = 0; j < 6; j++) {
                const idx = i * 6 + j;
                Dzz1[idx] += w * (Gx[i]*Gx[j] + Gy[i]*Gy[j]);
                const prod = Lz[i]*Lz[j];
                Dzz2Re[idx] += w * epsRe * prod;
                Dzz2Im[idx] += w * epsIm * prod;
            }
        }
    }

    // Zero out face-face entries in Dtt (EMerge convention: Dtt[3,3], Dtt[3,7], Dtt[7,3], Dtt[7,7] = 0)
    // This suppresses face DOF self-coupling in the B matrix
    for (const fi of [3, 7]) for (const fj of [3, 7]) Dtt[fi*8+fj] = 0;

    // Multiply by |Area| (Jacobian of the triangle mapping)
    for (let i = 0; i < 64; i++) { Att[i] *= Area; BttRe[i] *= Area; BttIm[i] *= Area; Dtt[i] *= Area; }
    for (let i = 0; i < 48; i++) { Dzt[i] *= Area; }
    for (let i = 0; i < 36; i++) { Dzz1[i] *= Area; Dzz2Re[i] *= Area; Dzz2Im[i] *= Area; }


    // Assemble 14x14 A and B:
    // A[0:8,0:8] = Att - k0² * Btt
    // B[0:8,0:8] = Dtt
    // B[8:14,0:8] = Dzt
    // B[0:8,8:14] = Dzt^T
    // B[8:14,8:14] = Dzz1 - k0² * Dzz2
    // (Note: EMerge uses B[8:,8:] = Dzz1 - k0²*Dzz2 with minus sign)

    const ARe = new Float64Array(196); // 14x14
    const AIm = new Float64Array(196);
    const BRe = new Float64Array(196);
    const BIm = new Float64Array(196);

    // A[0:8,0:8] = Att - k0² * Btt
    for (let i = 0; i < 8; i++) {
        for (let j = 0; j < 8; j++) {
            ARe[i*14+j] = Att[i*8+j] - k02 * BttRe[i*8+j];
            AIm[i*14+j] = -k02 * BttIm[i*8+j];
        }
    }

    // B[0:8,0:8] = Dtt
    for (let i = 0; i < 8; i++)
        for (let j = 0; j < 8; j++)
            BRe[i*14+j] = Dtt[i*8+j];

    // B[8:14,0:8] = Dzt
    for (let i = 0; i < 6; i++)
        for (let j = 0; j < 8; j++)
            BRe[(i+8)*14+j] = Dzt[i*8+j];

    // B[0:8,8:14] = Dzt^T
    for (let i = 0; i < 8; i++)
        for (let j = 0; j < 6; j++)
            BRe[i*14+(j+8)] = Dzt[j*8+i];

    // B[8:14,8:14] = Dzz1 - k0² * Dzz2 (matching EMerge sign convention)
    for (let i = 0; i < 6; i++) {
        for (let j = 0; j < 6; j++) {
            BRe[(i+8)*14+(j+8)] = Dzz1[i*6+j] - k02 * Dzz2Re[i*6+j];
            BIm[(i+8)*14+(j+8)] = -k02 * Dzz2Im[i*6+j];
        }
    }

    // Note: EMerge applies Ls (edge length) scaling to rows and columns.
    // For now, skip Ls scaling and use raw DOF values. The eigenvalues are
    // independent of DOF scaling; only conditioning is affected.

    return { ARe, AIm, BRe, BIm, Area, coeff, Dzz1, Dzz2Re, Dzz2Im, Att, Dtt, Dzt };
}

// --- P2 static element matrices (6x6 nodal only) ---

export function computeTriP2StaticMatrices(nodes, v0, v1, v2) {
    const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    const txs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
    const tys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];

    const Sz = new Float64Array(36); // 6x6 stiffness
    const Mz = new Float64Array(36); // 6x6 mass

    for (let q = 0; q < NQ; q++) {
        const w = QW[q];
        const xq = txs[0]*QL1[q] + txs[1]*QL2[q] + txs[2]*QL3[q];
        const yq = tys[0]*QL1[q] + tys[1]*QL2[q] + tys[2]*QL3[q];

        const Lz = new Float64Array(6);
        const Gx = new Float64Array(6);
        const Gy = new Float64Array(6);

        for (let k = 0; k < 3; k++) {
            Lz[k] = lv(coeff, k, xq, yq);
            const [gx, gy] = lvGrad(coeff, k, xq, yq);
            Gx[k] = gx; Gy[k] = gy;
        }
        for (let k = 0; k < 3; k++) {
            const [p, qq] = edgeVerts[k];
            Lz[k+3] = le(coeff, p, qq, xq, yq);
            const [gx, gy] = leGrad(coeff, p, qq, xq, yq);
            Gx[k+3] = gx; Gy[k+3] = gy;
        }

        for (let i = 0; i < 6; i++) {
            for (let j = 0; j < 6; j++) {
                const idx = i * 6 + j;
                Sz[idx] += w * (Gx[i]*Gx[j] + Gy[i]*Gy[j]);
                Mz[idx] += w * Lz[i]*Lz[j];
            }
        }
    }

    for (let i = 0; i < 36; i++) { Sz[i] *= Area; Mz[i] *= Area; }

    return { Sz, Mz, Area };
}

// --- P3 element matrices via Gauss quadrature ---
// Returns 25x25 A and B matrices for the generalized eigenvalue problem.
// DOF ordering: [ne1_0..2, nf1, ne2_0..2, nf2, ne3_0..2, nf3..nf6, lv_0..2, le_0..2, le3_0..2, lf3]
// First 14 DOFs (8 transverse + 6 longitudinal) match the P2 ordering.
const NT3 = 15; // P3 transverse DOFs
const NZ3 = 10; // P3 longitudinal DOFs
const NDOF3 = NT3 + NZ3; // 25 total

export function computeTriP3Matrices(nodes, v0, v1, v2, epsRe, epsIm, k0) {
    const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
    const k02 = k0 * k0;
    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    const Att = new Float64Array(NT3 * NT3);
    const BttRe = new Float64Array(NT3 * NT3);
    const BttIm = new Float64Array(NT3 * NT3);
    const Dtt = new Float64Array(NT3 * NT3);
    const Dzt = new Float64Array(NZ3 * NT3);
    const Dzz1 = new Float64Array(NZ3 * NZ3);
    const Dzz2Re = new Float64Array(NZ3 * NZ3);
    const Dzz2Im = new Float64Array(NZ3 * NZ3);

    const txs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
    const tys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];

    for (let q = 0; q < NQ12; q++) {
        const w = QW12[q];
        const xq = txs[0]*QL1_12[q] + txs[1]*QL2_12[q] + txs[2]*QL3_12[q];
        const yq = tys[0]*QL1_12[q] + tys[1]*QL2_12[q] + tys[2]*QL3_12[q];

        // Evaluate all 15 transverse basis functions and their curls
        // [0-2: ne1, 3: nf1, 4-6: ne2, 7: nf2, 8-10: ne3, 11-14: nf3..nf6]
        const Wx = new Float64Array(NT3), Wy = new Float64Array(NT3), curlW = new Float64Array(NT3);

        for (let k = 0; k < 3; k++) {
            const [p, qq] = edgeVerts[k];
            [Wx[k], Wy[k]] = ne1(coeff, p, qq, xq, yq);
            curlW[k] = ne1Curl(coeff, p, qq);
            [Wx[k+4], Wy[k+4]] = ne2(coeff, p, qq, xq, yq);
            curlW[k+4] = ne2Curl(coeff, p, qq, xq, yq);
            [Wx[k+8], Wy[k+8]] = ne3(coeff, p, qq, xq, yq);
            curlW[k+8] = ne3Curl(coeff, p, qq, xq, yq);
        }
        [Wx[3], Wy[3]] = nf1(coeff, xq, yq); curlW[3] = nf1Curl(coeff, xq, yq);
        [Wx[7], Wy[7]] = nf2(coeff, xq, yq); curlW[7] = nf2Curl(coeff, xq, yq);
        [Wx[11], Wy[11]] = nf3(coeff, xq, yq); curlW[11] = nf3Curl(coeff, xq, yq);
        [Wx[12], Wy[12]] = nf4(coeff, xq, yq); curlW[12] = nf4Curl(coeff, xq, yq);
        [Wx[13], Wy[13]] = nf5(coeff, xq, yq); curlW[13] = nf5Curl(coeff, xq, yq);
        [Wx[14], Wy[14]] = nf6(coeff, xq, yq); curlW[14] = nf6Curl(coeff, xq, yq);

        // Evaluate all 10 longitudinal basis functions and their gradients
        // [0-2: lv, 3-5: le, 6-8: le3, 9: lf3]
        const Lz = new Float64Array(NZ3), Gx = new Float64Array(NZ3), Gy = new Float64Array(NZ3);

        for (let k = 0; k < 3; k++) {
            Lz[k] = lv(coeff, k, xq, yq);
            [Gx[k], Gy[k]] = lvGrad(coeff, k, xq, yq);
        }
        for (let k = 0; k < 3; k++) {
            const [p, qq] = edgeVerts[k];
            Lz[k+3] = le(coeff, p, qq, xq, yq);
            [Gx[k+3], Gy[k+3]] = leGrad(coeff, p, qq, xq, yq);
            Lz[k+6] = le3(coeff, p, qq, xq, yq);
            [Gx[k+6], Gy[k+6]] = le3Grad(coeff, p, qq, xq, yq);
        }
        Lz[9] = lf3(coeff, xq, yq);
        [Gx[9], Gy[9]] = lf3Grad(coeff, xq, yq);

        // Accumulate sub-matrices
        for (let i = 0; i < NT3; i++) {
            for (let j = 0; j < NT3; j++) {
                const idx = i * NT3 + j;
                Att[idx] += w * curlW[i] * curlW[j];
                const dot = Wx[i]*Wx[j] + Wy[i]*Wy[j];
                BttRe[idx] += w * epsRe * dot;
                BttIm[idx] += w * epsIm * dot;
                Dtt[idx] += w * dot;
            }
        }

        for (let i = 0; i < NZ3; i++)
            for (let j = 0; j < NT3; j++)
                Dzt[i * NT3 + j] += w * (Gx[i]*Wx[j] + Gy[i]*Wy[j]);

        for (let i = 0; i < NZ3; i++) {
            for (let j = 0; j < NZ3; j++) {
                const idx = i * NZ3 + j;
                Dzz1[idx] += w * (Gx[i]*Gx[j] + Gy[i]*Gy[j]);
                const prod = Lz[i]*Lz[j];
                Dzz2Re[idx] += w * epsRe * prod;
                Dzz2Im[idx] += w * epsIm * prod;
            }
        }
    }

    // Zero out P2 face-face entries in Dtt (EMerge convention)
    for (const fi of [3, 7]) for (const fj of [3, 7]) Dtt[fi*NT3+fj] = 0;

    // Multiply by Area
    for (let i = 0; i < NT3*NT3; i++) { Att[i] *= Area; BttRe[i] *= Area; BttIm[i] *= Area; Dtt[i] *= Area; }
    for (let i = 0; i < NZ3*NT3; i++) { Dzt[i] *= Area; }
    for (let i = 0; i < NZ3*NZ3; i++) { Dzz1[i] *= Area; Dzz2Re[i] *= Area; Dzz2Im[i] *= Area; }

    // Assemble 25x25 A and B
    const ARe = new Float64Array(NDOF3 * NDOF3);
    const AIm = new Float64Array(NDOF3 * NDOF3);
    const BRe = new Float64Array(NDOF3 * NDOF3);
    const BIm = new Float64Array(NDOF3 * NDOF3);

    for (let i = 0; i < NT3; i++) {
        for (let j = 0; j < NT3; j++) {
            ARe[i*NDOF3+j] = Att[i*NT3+j] - k02 * BttRe[i*NT3+j];
            AIm[i*NDOF3+j] = -k02 * BttIm[i*NT3+j];
            BRe[i*NDOF3+j] = Dtt[i*NT3+j];
        }
    }

    for (let i = 0; i < NZ3; i++)
        for (let j = 0; j < NT3; j++)
            BRe[(i+NT3)*NDOF3+j] = Dzt[i*NT3+j];

    for (let i = 0; i < NT3; i++)
        for (let j = 0; j < NZ3; j++)
            BRe[i*NDOF3+(j+NT3)] = Dzt[j*NT3+i];

    for (let i = 0; i < NZ3; i++) {
        for (let j = 0; j < NZ3; j++) {
            BRe[(i+NT3)*NDOF3+(j+NT3)] = Dzz1[i*NZ3+j] - k02 * Dzz2Re[i*NZ3+j];
            BIm[(i+NT3)*NDOF3+(j+NT3)] = -k02 * Dzz2Im[i*NZ3+j];
        }
    }

    return { ARe, AIm, BRe, BIm, Area, coeff, Dzz1, Dzz2Re, Dzz2Im, Att, Dtt, Dzt, NT: NT3, NZ: NZ3, NDOF: NDOF3 };
}

// --- P3 static element matrices (10x10 nodal only) ---

export function computeTriP3StaticMatrices(nodes, v0, v1, v2) {
    const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    const txs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
    const tys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];

    const Sz = new Float64Array(100); // 10x10 stiffness
    const Mz = new Float64Array(100); // 10x10 mass

    for (let q = 0; q < NQ12; q++) {
        const w = QW12[q];
        const xq = txs[0]*QL1_12[q] + txs[1]*QL2_12[q] + txs[2]*QL3_12[q];
        const yq = tys[0]*QL1_12[q] + tys[1]*QL2_12[q] + tys[2]*QL3_12[q];

        const Lz = new Float64Array(NZ3), Gx = new Float64Array(NZ3), Gy = new Float64Array(NZ3);

        for (let k = 0; k < 3; k++) {
            Lz[k] = lv(coeff, k, xq, yq);
            [Gx[k], Gy[k]] = lvGrad(coeff, k, xq, yq);
        }
        for (let k = 0; k < 3; k++) {
            const [p, qq] = edgeVerts[k];
            Lz[k+3] = le(coeff, p, qq, xq, yq);
            [Gx[k+3], Gy[k+3]] = leGrad(coeff, p, qq, xq, yq);
            Lz[k+6] = le3(coeff, p, qq, xq, yq);
            [Gx[k+6], Gy[k+6]] = le3Grad(coeff, p, qq, xq, yq);
        }
        Lz[9] = lf3(coeff, xq, yq);
        [Gx[9], Gy[9]] = lf3Grad(coeff, xq, yq);

        for (let i = 0; i < NZ3; i++) {
            for (let j = 0; j < NZ3; j++) {
                Sz[i*NZ3+j] += w * (Gx[i]*Gx[j] + Gy[i]*Gy[j]);
                Mz[i*NZ3+j] += w * Lz[i]*Lz[j];
            }
        }
    }

    for (let i = 0; i < 100; i++) { Sz[i] *= Area; Mz[i] *= Area; }

    return { Sz, Mz, Area };
}

// --- Longitudinal DOF offsets ---
// Returns the global index offsets for each longitudinal DOF group.
export function getLzOffsets(fm) {
    const lzOff = fm.nFreeTransverse;
    const lzEdgeMidOff = lzOff + fm.nFreeVertexDof;
    const lzEdge3Off = lzEdgeMidOff + fm.nFreeEdgeNodeDof;
    const lzFaceNodeOff = lzEdge3Off + (fm.nFreeEdgeNode3Dof || 0);
    return { lzOff, lzEdgeMidOff, lzEdge3Off, lzFaceNodeOff };
}

// --- Freedom map ---
// Global DOF layout (P2-only, backward compatible):
// Transverse: [0, 2*nEdges) edge DOFs + [2*nEdges, 2*nEdges+2*nTris) face DOFs
// Longitudinal: [n_xy, n_xy+nNodes) vertex DOFs + [n_xy+nNodes, n_xy+nNodes+nEdges) edge DOFs
//
// Variable-order (when elemOrder is provided):
// Transverse: [edge Nedelec DOFs (2-3/edge)] [face DOFs (2-6/tri)]
// Longitudinal: [vertex DOFs] [edge midpoint DOFs (1/edge)] [edge extra DOFs (0-1/edge)] [face node DOFs (0-1/tri)]

export function buildTriFreedomMap(mesh, condRect, abc, elemOrder = null) {
    const { nodes, edges, nNodes, nEdges, nTris } = mesh;
    const TOL = 1e-12;
    const xmin = condRect.xmin_domain, xmax = condRect.xmax_domain;
    const ymax = condRect.ymax_domain;
    // Bottom wall: generalized from the legacy "y=0 is PEC ground" assumption to
    // an explicit domain floor. ymin_domain defaults to 0 for backward compat.
    // abc.bottom (like abc.left/right/top) makes the floor a natural (non-PEC)
    // boundary instead of PEC ground.
    const ymin = condRect.ymin_domain ?? 0;
    const condRects = condRect.rects || [condRect];

    // Identify conductor nodes and edges
    const isCondNode = new Uint8Array(nNodes);
    const condNodeGroup = new Int8Array(nNodes);
    for (let n = 0; n < nNodes; n++) {
        const x = nodes[2 * n], y = nodes[2 * n + 1];
        for (let ci = 0; ci < condRects.length; ci++) {
            const cr = condRects[ci];
            if (x >= cr.xmin - TOL && x <= cr.xmax + TOL &&
                y >= cr.ymin - TOL && y <= cr.ymax + TOL) {
                isCondNode[n] = 1;
                condNodeGroup[n] = ci + 1;
                break;
            }
        }
    }

    const isCondEdge = new Uint8Array(nEdges);
    for (let e = 0; e < nEdges; e++) {
        if (isCondNode[edges[2*e]] && isCondNode[edges[2*e+1]]) isCondEdge[e] = 1;
    }

    // Compute edge order: min of adjacent element orders (minimum rule)
    // For edges on a single triangle (boundary), use that triangle's order.
    let edgeOrder = null;
    if (elemOrder) {
        edgeOrder = new Uint8Array(nEdges).fill(255);
        for (let t = 0; t < nTris; t++) {
            const order = elemOrder[t];
            for (let k = 0; k < 3; k++) {
                const e = mesh.triEdges[3*t+k];
                if (order < edgeOrder[e]) edgeOrder[e] = order;
            }
        }
        // Edges not touched by any triangle stay at 255 → default to 2
        for (let e = 0; e < nEdges; e++) if (edgeOrder[e] === 255) edgeOrder[e] = 2;
    }

    // Helper: check if edge is on a PEC boundary (ground, conductor, walls)
    function isEdgePEC(e) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const y0 = nodes[2*n0+1], y1 = nodes[2*n1+1];
        const x0 = nodes[2*n0], x1 = nodes[2*n1];
        if (!abc.bottom && Math.abs(y0-ymin) < TOL && Math.abs(y1-ymin) < TOL) return true;
        if (isCondEdge[e]) return true;
        if (!abc.left && Math.abs(x0-xmin)<TOL && Math.abs(x1-xmin)<TOL) return true;
        if (!abc.right && Math.abs(x0-xmax)<TOL && Math.abs(x1-xmax)<TOL) return true;
        if (!abc.top && Math.abs(y0-ymax)<TOL && Math.abs(y1-ymax)<TOL) return true;
        return false;
    }

    // --- Transverse edge DOFs ---
    // P2: 2 per edge (ne1, ne2). P3: 3 per edge (ne1, ne2, ne3).
    // edgeF[2*e], edgeF[2*e+1]: P2 DOFs (always populated)
    // edgeF3[e]: P3 ne3 DOF index (or -1)
    const edgeF = new Int32Array(2 * nEdges).fill(-1);
    const edgeF3 = elemOrder ? new Int32Array(nEdges).fill(-1) : null;
    let nFreeEdgeDof = 0;
    for (let e = 0; e < nEdges; e++) {
        if (isEdgePEC(e)) continue;
        edgeF[2*e] = nFreeEdgeDof++;     // ne1
        edgeF[2*e+1] = nFreeEdgeDof++;   // ne2
        if (edgeOrder && edgeOrder[e] >= 3) {
            edgeF3[e] = nFreeEdgeDof++;   // ne3
        }
    }

    // --- Transverse face DOFs ---
    // P2: 2 per tri (nf1, nf2). P3: 6 per tri (nf1, nf2, nf3..nf6).
    const faceF = new Int32Array(2 * nTris).fill(-1);
    const faceF3 = elemOrder ? new Int32Array(4 * nTris).fill(-1) : null; // nf3..nf6
    let nFreeFaceDof = 0;
    for (let t = 0; t < nTris; t++) {
        const v0 = mesh.tris[3*t], v1 = mesh.tris[3*t+1], v2 = mesh.tris[3*t+2];
        const yc = (nodes[2*v0+1]+nodes[2*v1+1]+nodes[2*v2+1])/3;
        const xc = (nodes[2*v0]+nodes[2*v1]+nodes[2*v2])/3;
        let inCond = false;
        for (const cr of condRects) {
            if (xc >= cr.xmin-TOL && xc <= cr.xmax+TOL &&
                yc >= cr.ymin-TOL && yc <= cr.ymax+TOL) { inCond = true; break; }
        }
        if (inCond) continue;
        faceF[2*t] = nFreeEdgeDof + nFreeFaceDof++;   // nf1
        faceF[2*t+1] = nFreeEdgeDof + nFreeFaceDof++; // nf2
        if (elemOrder && elemOrder[t] >= 3) {
            for (let k = 0; k < 4; k++) {
                faceF3[4*t+k] = nFreeEdgeDof + nFreeFaceDof++; // nf3..nf6
            }
        }
    }

    const nFreeTransverse = nFreeEdgeDof + nFreeFaceDof;

    // --- Longitudinal vertex DOFs (lv, same for P2/P3) ---
    const nodeF = new Int32Array(nNodes).fill(-1);
    let nFreeVertexDof = 0;
    for (let n = 0; n < nNodes; n++) {
        const x = nodes[2*n], y = nodes[2*n+1];
        if (!abc.bottom && Math.abs(y-ymin) < TOL) continue;
        if (isCondNode[n]) continue;
        if (!abc.left && Math.abs(x-xmin)<TOL) continue;
        if (!abc.right && Math.abs(x-xmax)<TOL) continue;
        if (!abc.top && Math.abs(y-ymax)<TOL) continue;
        nodeF[n] = nFreeVertexDof++;
    }

    // --- Longitudinal edge midpoint DOFs (le, same for P2/P3) ---
    const edgeNodeF = new Int32Array(nEdges).fill(-1);
    let nFreeEdgeNodeDof = 0;
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const ym = (nodes[2*n0+1]+nodes[2*n1+1])/2;
        const xm = (nodes[2*n0]+nodes[2*n1])/2;
        if (!abc.bottom && Math.abs(ym-ymin) < TOL) continue;
        if (isCondNode[n0] && isCondNode[n1]) continue;
        if (!abc.left && Math.abs(xm-xmin)<TOL) continue;
        if (!abc.right && Math.abs(xm-xmax)<TOL) continue;
        if (!abc.top && Math.abs(ym-ymax)<TOL) continue;
        edgeNodeF[e] = nFreeEdgeNodeDof++;
    }

    // --- Longitudinal edge extra DOFs (le3, P3 only) ---
    const edgeNodeF3 = elemOrder ? new Int32Array(nEdges).fill(-1) : null;
    let nFreeEdgeNode3Dof = 0;
    if (elemOrder) {
        for (let e = 0; e < nEdges; e++) {
            if (edgeOrder[e] < 3) continue;
            // Same PEC/conductor conditions as edge midpoint
            if (edgeNodeF[e] < 0) continue; // if midpoint was eliminated, extra is too
            edgeNodeF3[e] = nFreeEdgeNode3Dof++;
        }
    }

    // --- Longitudinal face bubble DOFs (lf3, P3 only) ---
    const faceNodeF = elemOrder ? new Int32Array(nTris).fill(-1) : null;
    let nFreeFaceNodeDof = 0;
    if (elemOrder) {
        for (let t = 0; t < nTris; t++) {
            if (elemOrder[t] < 3) continue;
            if (faceF[2*t] < 0) continue; // conductor interior
            faceNodeF[t] = nFreeFaceNodeDof++;
        }
    }

    const nFreeLongitudinal = nFreeVertexDof + nFreeEdgeNodeDof + nFreeEdgeNode3Dof + nFreeFaceNodeDof;

    return {
        edgeF, faceF, nodeF, edgeNodeF,
        edgeF3, faceF3, edgeNodeF3, faceNodeF,
        elemOrder, edgeOrder,
        nFreeEdgeDof, nFreeFaceDof, nFreeTransverse,
        nFreeVertexDof, nFreeEdgeNodeDof, nFreeEdgeNode3Dof: nFreeEdgeNode3Dof || 0,
        nFreeFaceNodeDof: nFreeFaceNodeDof || 0, nFreeLongitudinal,
        isCondNode, isCondEdge, condNodeGroup,
    };
}

// --- Assembly ---

export function assembleTriFEM(mesh, fm, k2, epsMap, abc, condRect, enrichment) {
    const { tris, edges, triEdges, triSigns, nTris, nEdges } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF, nFreeTransverse, nFreeVertexDof,
            edgeF3, faceF3, edgeNodeF3, faceNodeF, elemOrder,
            nFreeEdgeNodeDof, nFreeEdgeNode3Dof, nFreeFaceNodeDof } = fm;
    const nCorners = enrichment ? enrichment.corners.length : 0;
    const nCornerDofs = nCorners * 2;
    const N = fm.nFreeTransverse + fm.nFreeLongitudinal + nCornerDofs;
    const cornerDofOff = fm.nFreeTransverse + fm.nFreeLongitudinal;

    const nodes = mesh.nodes;
    const k0 = Math.sqrt(k2);

    const { lzOff, lzEdgeMidOff, lzEdge3Off, lzFaceNodeOff } = getLzOffsets(fm);

    // Pre-allocate COO arrays: max 625 nonzeros per P3 element (25×25)
    const maxNnz = nTris * (elemOrder ? 625 : 196);
    const Ar = new Int32Array(maxNnz), Ac = new Int32Array(maxNnz);
    const AvRe = new Float64Array(maxNnz), AvIm = new Float64Array(maxNnz);
    const Br = new Int32Array(maxNnz), Bc = new Int32Array(maxNnz);
    const BvRe = new Float64Array(maxNnz), BvIm = new Float64Array(maxNnz);
    let aNnz = 0, bNnz = 0;
    const RBr = [], RBc = [], RBvRe = [], RBvIm = [];
    const RAr = [], RAc = [], RAvRe = [], RAvIm = [];

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const eps = epsMap[t];
        const order = elemOrder ? elemOrder[t] : 2;

        let ARe, AIm, BRe, BIm, Att, Dtt, nLocal;
        if (order >= 3) {
            const m = computeTriP3Matrices(nodes, v0, v1, v2, eps.re, eps.im, k0);
            ARe = m.ARe; AIm = m.AIm; BRe = m.BRe; BIm = m.BIm; Att = m.Att; Dtt = m.Dtt;
            nLocal = NDOF3; // 25
        } else {
            const m = computeTriP2Matrices(nodes, v0, v1, v2, eps.re, eps.im, k0);
            ARe = m.ARe; AIm = m.AIm; BRe = m.BRe; BIm = m.BIm; Att = m.Att; Dtt = m.Dtt;
            nLocal = 14;
        }

        if (t === 0 && globalThis.__TRI_DEBUG__) {
            const ntrans = order >= 3 ? NT3 : 8;
            let attMax = 0, dttMax = 0;
            for (let i = 0; i < ntrans*ntrans; i++) { attMax = Math.max(attMax, Math.abs(Att[i])); dttMax = Math.max(dttMax, Math.abs(Dtt[i])); }
            let areMax = 0, breMax = 0;
            for (let i = 0; i < nLocal*nLocal; i++) { areMax = Math.max(areMax, Math.abs(ARe[i])); breMax = Math.max(breMax, Math.abs(BRe[i])); }
            globalThis.__TRI_DEBUG__ && console.log(`  Elem 0: Area=${triCoefficients(nodes,v0,v1,v2).Area.toExponential(3)}, |Att|=${attMax.toExponential(3)}, |Dtt|=${dttMax.toExponential(3)}, |A|=${areMax.toExponential(3)}, |B|=${breMax.toExponential(3)}`);
        }

        // Build local-to-global DOF mapping
        const globalDof = new Int32Array(nLocal);
        const signs = new Float64Array(nLocal).fill(1);

        for (let le = 0; le < 3; le++) {
            const eIdx = triEdges[3*t + le];
            const s = triSigns[3*t + le];

            // ne1 DOF (antisymmetric)
            globalDof[le] = edgeF[2*eIdx];
            signs[le] = s;

            // ne2 DOF (symmetric — no sign flip)
            globalDof[le+4] = edgeF[2*eIdx+1];

            // le DOF (edge midpoint nodal)
            const enf = edgeNodeF[eIdx];
            if (order >= 3) {
                // P3 local: [ne1_0..2, nf1, ne2_0..2, nf2, ne3_0..2, nf3..6, lv_0..2, le_0..2, le3_0..2, lf3]
                globalDof[le+18] = enf >= 0 ? lzEdgeMidOff + enf : -1; // le at offset 18
                // ne3 DOF (antisymmetric like ne1)
                globalDof[le+8] = edgeF3 ? edgeF3[eIdx] : -1;
                signs[le+8] = s;
                // le3 DOF (edge extra, antisymmetric — needs sign flip like ne1)
                const enf3 = edgeNodeF3 ? edgeNodeF3[eIdx] : -1;
                globalDof[le+21] = enf3 >= 0 ? lzEdge3Off + enf3 : -1;
                signs[le+21] = s; // le3(i,j) = -le3(j,i)
            } else {
                // P2 local: [ne1_0..2, nf1, ne2_0..2, nf2, lv_0..2, le_0..2]
                globalDof[le+11] = enf >= 0 ? lzEdgeMidOff + enf : -1;
            }
        }

        // Face DOFs
        globalDof[3] = faceF[2*t];   // nf1
        globalDof[7] = faceF[2*t+1]; // nf2
        if (order >= 3 && faceF3) {
            for (let k = 0; k < 4; k++)
                globalDof[11+k] = faceF3[4*t+k]; // nf3..nf6
        }

        // Vertex DOFs (lv)
        const verts = [v0, v1, v2];
        const lvOff = order >= 3 ? 15 : 8;
        for (let k = 0; k < 3; k++) {
            const nf = nodeF[verts[k]];
            globalDof[lvOff + k] = nf >= 0 ? lzOff + nf : -1;
        }
        // le DOFs already set above

        // Face node DOF (lf3, P3 only)
        if (order >= 3 && faceNodeF) {
            const fnf = faceNodeF[t];
            globalDof[24] = fnf >= 0 ? lzFaceNodeOff + fnf : -1;
        }

        // Assemble element matrices into global COO
        for (let li = 0; li < nLocal; li++) {
            const gi = globalDof[li]; if (gi < 0) continue;
            for (let lj = 0; lj < nLocal; lj++) {
                const gj = globalDof[lj]; if (gj < 0) continue;
                const s = signs[li] * signs[lj];
                const idx = li * nLocal + lj;
                const are = s * ARe[idx], aim = s * AIm[idx];
                if (are !== 0 || aim !== 0) {
                    Ar[aNnz] = gi; Ac[aNnz] = gj; AvRe[aNnz] = are; AvIm[aNnz] = aim; aNnz++;
                }
                const bre = s * BRe[idx], bim = s * BIm[idx];
                if (bre !== 0 || bim !== 0) {
                    Br[bNnz] = gi; Bc[bNnz] = gj; BvRe[bNnz] = bre; BvIm[bNnz] = bim; bNnz++;
                }
            }
        }
    }

    // Robin ABC: boundary edges with P2 DOFs (strict equality: 'pmc' skips Robin terms)
    if (abc.top === true || abc.left === true || abc.right === true) {
        const ymax = condRect.ymax_domain;
        const xmin = condRect.xmin_domain, xmax = condRect.xmax_domain;
        const TOL = 1e-12;
        for (let e = 0; e < nEdges; e++) {
            const n0 = edges[2*e], n1 = edges[2*e+1];
            const x0 = nodes[2*n0], y0 = nodes[2*n0+1];
            const x1 = nodes[2*n1], y1 = nodes[2*n1+1];
            const L = Math.sqrt((x1-x0)**2 + (y1-y0)**2);
            let isBoundary = false;
            if (abc.top === true && Math.abs(y0 - ymax) < TOL && Math.abs(y1 - ymax) < TOL) isBoundary = true;
            if (abc.left === true && Math.abs(x0 - xmin) < TOL && Math.abs(x1 - xmin) < TOL) isBoundary = true;
            if (abc.right === true && Math.abs(x0 - xmax) < TOL && Math.abs(x1 - xmax) < TOL) isBoundary = true;
            if (!isBoundary) continue;

            // Robin ABC for transverse edge DOFs: A += jk₀ ∫ W_tang · W_tang dl
            // On edge (p,q) parameterized by s∈[0,1]:
            //   ne1_tang(s) = 1/L (constant), ne2_tang(s) = (1/L)(1-2s)
            // Self-integrals: ∫ne1²dl = 1/L, ∫ne2²dl = 1/(3L), ∫ne1·ne2 dl = 0
            const ef1 = edgeF[2*e];
            if (ef1 >= 0) { RAr.push(ef1); RAc.push(ef1); RAvRe.push(0); RAvIm.push(k0 / L); }

            const ef2 = edgeF[2*e+1];
            if (ef2 >= 0) { RAr.push(ef2); RAc.push(ef2); RAvRe.push(0); RAvIm.push(k0 / (3 * L)); }

            // P2 nodal Robin on boundary: DISABLED — the Ez ABC distorts the
            // eigenvector by forcing ∂Ez/∂n + jk₀Ez = 0, which is a poor approximation
            // for quasi-TEM modes. The transverse ABC alone suppresses cavity modes.
        }
    }

    // Corner singularity enrichment (de Rham compatible).
    // Each corner adds 2 DOFs: 1 transverse (S_t = ∇φ_s) + 1 longitudinal (φ_s).
    // S_t = ∇φ_s ensures grad(Ez_enriched) ⊂ Et_enriched (de Rham).
    // curl(S_t) = 0 (gradient), so St (curl-curl) coupling is zero.
    //
    // DOF ordering: [...standard..., et_corner0, ..., et_cornerN, ez_corner0, ..., ez_cornerN]
    // gt = cornerDofOff + ci (transverse), gz = cornerDofOff + nCorners + ci (longitudinal)
    //
    // Matrix conventions (from computeTriP2Matrices):
    //   A_tt = Att - k₀²ε·Mt       (Att = ∫curl·curl, no ε; ε·Mt in subtracted term)
    //   B_tt = Dtt = ∫W·W           (no ε)
    //   B_zt = Dzt = ∫∇N·W          (no ε)
    //   B_zz = Dzz1 - k₀²·Dzz2     (Dzz1 = ∫∇N·∇N no ε; Dzz2 = ε·∫N·N)
    if (enrichment) {
        const { corners, supports, evalFn, duffyQuadPoints: duffyQP } = enrichment;
        const edgeVerts = [[0, 1], [1, 2], [2, 0]];

        for (let ci = 0; ci < corners.length; ci++) {
            const gt = cornerDofOff + ci;               // transverse enrichment DOF
            const gz = cornerDofOff + nCorners + ci;     // longitudinal enrichment DOF
            const corner = corners[ci];
            const { tris: supportTris, cornerVerts } = supports[ci];

            for (let si = 0; si < supportTris.length; si++) {
                const t = supportTris[si];
                const cornerLocalIdx = cornerVerts[si];
                const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
                if (faceF[2*t] < 0) continue; // skip conductor interior

                const eps = epsMap[t];
                const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
                const txs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
                const tys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];
                const verts = [v0, v1, v2];

                // Quadrature: Duffy for corner-touching, standard otherwise
                const useDuffy = cornerLocalIdx >= 0 && duffyQP;
                const qpts = useDuffy ? duffyQP(cornerLocalIdx) : null;
                const nqp = useDuffy ? qpts.length : NQ;

                for (let q = 0; q < nqp; q++) {
                    let w, xq, yq;
                    if (useDuffy) {
                        const qp = qpts[q];
                        w = qp.w;
                        xq = txs[0]*qp.l1 + txs[1]*qp.l2 + txs[2]*qp.l3;
                        yq = tys[0]*qp.l1 + tys[1]*qp.l2 + tys[2]*qp.l3;
                    } else {
                        w = QW[q];
                        xq = txs[0]*QL1[q] + txs[1]*QL2[q] + txs[2]*QL3[q];
                        yq = tys[0]*QL1[q] + tys[1]*QL2[q] + tys[2]*QL3[q];
                    }

                    // Enrichment: φ_s (scalar) and S_t = ∇φ_s (vector)
                    const { phi, dphidx, dphidy } = evalFn(xq, yq, corner);
                    if (phi === 0 && dphidx === 0 && dphidy === 0) continue;
                    // S_t = (dphidx, dphidy) — the vector enrichment

                    // Standard Nedelec basis at quadrature point
                    const Wx = new Float64Array(8), Wy = new Float64Array(8);
                    for (let k = 0; k < 3; k++) {
                        const [p, qq] = edgeVerts[k];
                        [Wx[k], Wy[k]] = ne1(coeff, p, qq, xq, yq);
                        [Wx[k+4], Wy[k+4]] = ne2(coeff, p, qq, xq, yq);
                    }
                    [Wx[3], Wy[3]] = nf1(coeff, xq, yq);
                    [Wx[7], Wy[7]] = nf2(coeff, xq, yq);

                    const wA = w * Area;
                    const St_dot_St = dphidx*dphidx + dphidy*dphidy; // |S_t|² = |∇φ_s|²

                    // === A matrix entries ===
                    // A[gt, j_t] = -k₀²ε · ∫ S_t·W_j dA  (curl(S_t)=0, so Att=0)
                    // Complex ε: ε = eps.re + j·eps.im (eps.im < 0 for lossy materials)
                    function addA_enr(gi, gj, dot) {
                        const vre = -k2 * wA * (eps.re * dot);
                        const vim = -k2 * wA * (eps.im * dot);
                        if (vre !== 0 || vim !== 0) {
                            RAr.push(gi); RAc.push(gj); RAvRe.push(vre); RAvIm.push(vim);
                            RAr.push(gj); RAc.push(gi); RAvRe.push(vre); RAvIm.push(vim);
                        }
                    }
                    for (let k = 0; k < 3; k++) {
                        const eIdx = triEdges[3*t+k], s = triSigns[3*t+k];
                        const ef1 = edgeF[2*eIdx];
                        if (ef1 >= 0) addA_enr(gt, ef1, dphidx*s*Wx[k] + dphidy*s*Wy[k]);
                        const ef2 = edgeF[2*eIdx+1];
                        if (ef2 >= 0) addA_enr(gt, ef2, dphidx*Wx[k+4] + dphidy*Wy[k+4]);
                    }
                    { const f1=faceF[2*t]; if(f1>=0) addA_enr(gt, f1, dphidx*Wx[3]+dphidy*Wy[3]); }
                    { const f2=faceF[2*t+1]; if(f2>=0) addA_enr(gt, f2, dphidx*Wx[7]+dphidy*Wy[7]); }
                    // A[gt, gt] self
                    { const vre = -k2*wA*eps.re*St_dot_St, vim = -k2*wA*eps.im*St_dot_St;
                      if(vre!==0||vim!==0){ RAr.push(gt);RAc.push(gt);RAvRe.push(vre);RAvIm.push(vim); }}

                    // === B matrix: transverse block (Dtt-like, NO ε) ===
                    // B[gt, j_t] = ∫ S_t·W_j dA
                    for (let k = 0; k < 3; k++) {
                        const eIdx = triEdges[3*t+k], s = triSigns[3*t+k];
                        const ef1 = edgeF[2*eIdx];
                        if (ef1 >= 0) {
                            const v = wA * (dphidx*s*Wx[k] + dphidy*s*Wy[k]);
                            if (v !== 0) { RBr.push(gt);RBc.push(ef1);RBvRe.push(v);RBvIm.push(0);
                                           RBr.push(ef1);RBc.push(gt);RBvRe.push(v);RBvIm.push(0); }
                        }
                        const ef2 = edgeF[2*eIdx+1];
                        if (ef2 >= 0) {
                            const v = wA * (dphidx*Wx[k+4] + dphidy*Wy[k+4]);
                            if (v !== 0) { RBr.push(gt);RBc.push(ef2);RBvRe.push(v);RBvIm.push(0);
                                           RBr.push(ef2);RBc.push(gt);RBvRe.push(v);RBvIm.push(0); }
                        }
                    }
                    { const f1=faceF[2*t]; if(f1>=0){ const v=wA*(dphidx*Wx[3]+dphidy*Wy[3]);
                      if(v!==0){ RBr.push(gt);RBc.push(f1);RBvRe.push(v);RBvIm.push(0); RBr.push(f1);RBc.push(gt);RBvRe.push(v);RBvIm.push(0); }}}
                    { const f2=faceF[2*t+1]; if(f2>=0){ const v=wA*(dphidx*Wx[7]+dphidy*Wy[7]);
                      if(v!==0){ RBr.push(gt);RBc.push(f2);RBvRe.push(v);RBvIm.push(0); RBr.push(f2);RBc.push(gt);RBvRe.push(v);RBvIm.push(0); }}}
                    // B[gt, gt] self = ∫|S_t|² dA (no ε)
                    { const v = wA*St_dot_St;
                      if(v!==0){ RBr.push(gt);RBc.push(gt);RBvRe.push(v);RBvIm.push(0); }}

                    // === B matrix: gradient coupling (Dzt-like, NO ε) ===
                    // B[gz, j_t] = ∫ ∇φ_s·W_j dA (same as B[gt, j_t] since S_t = ∇φ_s)
                    for (let k = 0; k < 3; k++) {
                        const eIdx = triEdges[3*t+k], s = triSigns[3*t+k];
                        const ef1 = edgeF[2*eIdx];
                        if (ef1 >= 0) {
                            const v = wA * (dphidx*s*Wx[k] + dphidy*s*Wy[k]);
                            if (v !== 0) { RBr.push(gz);RBc.push(ef1);RBvRe.push(v);RBvIm.push(0);
                                           RBr.push(ef1);RBc.push(gz);RBvRe.push(v);RBvIm.push(0); }
                        }
                        const ef2 = edgeF[2*eIdx+1];
                        if (ef2 >= 0) {
                            const v = wA * (dphidx*Wx[k+4] + dphidy*Wy[k+4]);
                            if (v !== 0) { RBr.push(gz);RBc.push(ef2);RBvRe.push(v);RBvIm.push(0);
                                           RBr.push(ef2);RBc.push(gz);RBvRe.push(v);RBvIm.push(0); }
                        }
                    }
                    { const f1=faceF[2*t]; if(f1>=0){ const v=wA*(dphidx*Wx[3]+dphidy*Wy[3]);
                      if(v!==0){ RBr.push(gz);RBc.push(f1);RBvRe.push(v);RBvIm.push(0); RBr.push(f1);RBc.push(gz);RBvRe.push(v);RBvIm.push(0); }}}
                    { const f2=faceF[2*t+1]; if(f2>=0){ const v=wA*(dphidx*Wx[7]+dphidy*Wy[7]);
                      if(v!==0){ RBr.push(gz);RBc.push(f2);RBvRe.push(v);RBvIm.push(0); RBr.push(f2);RBc.push(gz);RBvRe.push(v);RBvIm.push(0); }}}

                    // === B matrix: longitudinal block (Dzz-like) ===
                    // B[gz, j_z] = ∫∇φ_s·∇N_j dA - k₀²ε·∫φ_s·N_j dA (Dzz1 - k₀²Dzz2)
                    // Dzz1 has no ε; Dzz2 has complex ε
                    for (let k = 0; k < 3; k++) {
                        const nf = nodeF[verts[k]]; if (nf < 0) continue;
                        const gj = nFreeTransverse + nf;
                        const [gx,gy] = lvGrad(coeff, k, xq, yq);
                        const Nv = lv(coeff, k, xq, yq);
                        const vre = wA * ((dphidx*gx + dphidy*gy) - k2*eps.re*phi*Nv);
                        const vim = wA * (-k2*eps.im*phi*Nv);
                        if (vre !== 0 || vim !== 0) { RBr.push(gz);RBc.push(gj);RBvRe.push(vre);RBvIm.push(vim);
                                       RBr.push(gj);RBc.push(gz);RBvRe.push(vre);RBvIm.push(vim); }
                    }
                    for (let k = 0; k < 3; k++) {
                        const enf = edgeNodeF[triEdges[3*t+k]]; if (enf < 0) continue;
                        const gj = nFreeTransverse + nFreeVertexDof + enf;
                        const [p,qq] = edgeVerts[k];
                        const [gx,gy] = leGrad(coeff, p, qq, xq, yq);
                        const Nv = le(coeff, p, qq, xq, yq);
                        const vre = wA * ((dphidx*gx + dphidy*gy) - k2*eps.re*phi*Nv);
                        const vim = wA * (-k2*eps.im*phi*Nv);
                        if (vre !== 0 || vim !== 0) { RBr.push(gz);RBc.push(gj);RBvRe.push(vre);RBvIm.push(vim);
                                       RBr.push(gj);RBc.push(gz);RBvRe.push(vre);RBvIm.push(vim); }
                    }

                    // === B matrix: self-couplings ===
                    // B[gt, gz] = ∫∇φ_s·∇φ_s dA = ∫|S_t|² dA (Dzt self, no ε)
                    { const v = wA*St_dot_St;
                      if(v!==0){ RBr.push(gt);RBc.push(gz);RBvRe.push(v);RBvIm.push(0);
                                 RBr.push(gz);RBc.push(gt);RBvRe.push(v);RBvIm.push(0); }}
                    // B[gz, gz] = ∫|∇φ_s|² dA - k₀²ε·∫φ_s² dA (Dzz self)
                    { const vre = wA*(St_dot_St - k2*eps.re*phi*phi);
                      const vim = wA*(-k2*eps.im*phi*phi);
                      if(vre!==0||vim!==0){ RBr.push(gz);RBc.push(gz);RBvRe.push(vre);RBvIm.push(vim); }}

                    // === Cross-coupling between corners (if overlapping supports) ===
                    for (let cj = ci + 1; cj < nCorners; cj++) {
                        const gtj = cornerDofOff + cj, gzj = cornerDofOff + nCorners + cj;
                        const { phi:pj, dphidx:dxj, dphidy:dyj } = evalFn(xq, yq, corners[cj]);
                        if (pj === 0 && dxj === 0 && dyj === 0) continue;
                        const dot_ij = dphidx*dxj + dphidy*dyj;
                        // A[gt, gtj] — complex ε
                        { const vre = -k2*wA*eps.re*dot_ij, vim = -k2*wA*eps.im*dot_ij;
                          if(vre!==0||vim!==0){ RAr.push(gt);RAc.push(gtj);RAvRe.push(vre);RAvIm.push(vim);
                                     RAr.push(gtj);RAc.push(gt);RAvRe.push(vre);RAvIm.push(vim); }}
                        // B[gt, gtj] (Dtt cross — no ε)
                        { const v = wA*dot_ij;
                          if(v!==0){ RBr.push(gt);RBc.push(gtj);RBvRe.push(v);RBvIm.push(0);
                                     RBr.push(gtj);RBc.push(gt);RBvRe.push(v);RBvIm.push(0); }}
                        // B[gz, gzj] (Dzz cross — complex ε in mass term)
                        { const vre = wA*(dot_ij - k2*eps.re*phi*pj);
                          const vim = wA*(-k2*eps.im*phi*pj);
                          if(vre!==0||vim!==0){ RBr.push(gz);RBc.push(gzj);RBvRe.push(vre);RBvIm.push(vim);
                                     RBr.push(gzj);RBc.push(gz);RBvRe.push(vre);RBvIm.push(vim); }}
                        // B[gt, gzj] and B[gz, gtj] (Dzt cross — no ε)
                        { const v = wA*dot_ij;
                          if(v!==0){ RBr.push(gt);RBc.push(gzj);RBvRe.push(v);RBvIm.push(0);
                                     RBr.push(gzj);RBc.push(gt);RBvRe.push(v);RBvIm.push(0);
                                     RBr.push(gz);RBc.push(gtj);RBvRe.push(v);RBvIm.push(0);
                                     RBr.push(gtj);RBc.push(gz);RBvRe.push(v);RBvIm.push(0); }}
                    }
                }
            }
        }
    }

    // Merge pre-allocated element arrays with Robin BC dynamic arrays
    const totalA = aNnz + RAr.length;
    const fAr = new Int32Array(totalA), fAc = new Int32Array(totalA);
    const fAvRe = new Float64Array(totalA), fAvIm = new Float64Array(totalA);
    for (let i = 0; i < aNnz; i++) { fAr[i] = Ar[i]; fAc[i] = Ac[i]; fAvRe[i] = AvRe[i]; fAvIm[i] = AvIm[i]; }
    for (let i = 0; i < RAr.length; i++) { fAr[aNnz+i] = RAr[i]; fAc[aNnz+i] = RAc[i]; fAvRe[aNnz+i] = RAvRe[i]; fAvIm[aNnz+i] = RAvIm[i]; }

    const totalB = bNnz + RBr.length;
    const fBr = new Int32Array(totalB), fBc = new Int32Array(totalB);
    const fBvRe = new Float64Array(totalB), fBvIm = new Float64Array(totalB);
    for (let i = 0; i < bNnz; i++) { fBr[i] = Br[i]; fBc[i] = Bc[i]; fBvRe[i] = BvRe[i]; fBvIm[i] = BvIm[i]; }
    for (let i = 0; i < RBr.length; i++) { fBr[bNnz+i] = RBr[i]; fBc[bNnz+i] = RBc[i]; fBvRe[bNnz+i] = RBvRe[i]; fBvIm[bNnz+i] = RBvIm[i]; }

    return {
        csrA: tripletsToCSR(fAr, fAc, fAvRe, N, fAvIm),
        csrB: tripletsToCSR(fBr, fBc, fBvRe, N, fBvIm),
        N
    };
}

// --- P2 Static solver ---
// Uses 6 nodal DOFs per triangle: 3 vertex (lv) + 3 edge midpoint (le)

// enrichment: optional { corners, supports, evalFn, duffyQuadPoints } for corner singularities.
// When provided, adds scalar enrichment DOFs (one per corner) to capture r^ν singular behavior.
// Returns { phiVertex, phiEdge, enrichCoeffs }.
export function solveTriStatic(mesh, fm, epsMap, condPotentials = null, enrichment = null) {
    // condPotentials: array of potentials per conductor group [V1, V2, ...].
    // If null, all conductors get V=1.0.
    const { nodes, tris, nTris } = mesh;
    const { nodeF, edgeNodeF, nFreeVertexDof, nFreeEdgeNodeDof,
            isCondNode, isCondEdge, condNodeGroup,
            elemOrder, edgeOrder, edgeNodeF3, faceNodeF,
            nFreeEdgeNode3Dof, nFreeFaceNodeDof } = fm;
    const nCorners = enrichment ? enrichment.corners.length : 0;
    const nFreeFEM = nFreeVertexDof + nFreeEdgeNodeDof + (nFreeEdgeNode3Dof || 0) + (nFreeFaceNodeDof || 0);
    const nFreeDof = nFreeFEM + nCorners;

    // Longitudinal DOF offsets (within the static system, lzOff=0)
    const edgeMidOff = nFreeVertexDof;
    const edge3Off = edgeMidOff + nFreeEdgeNodeDof;
    const faceNodeOff = edge3Off + (nFreeEdgeNode3Dof || 0);
    // Note: these are offsets within the longitudinal block (no lzOff prefix)
    // because the static solver only has longitudinal DOFs.

    const Rows = [], Cols = [], Vals = [];
    const rhs = new Float64Array(nFreeDof);

    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const order = elemOrder ? elemOrder[t] : 2;
        const nLocal = order >= 3 ? 10 : 6;
        const Sz = order >= 3
            ? computeTriP3StaticMatrices(nodes, v0, v1, v2).Sz
            : computeTriP2StaticMatrices(nodes, v0, v1, v2).Sz;
        const eps = epsMap ? epsMap[t].re : 1.0;
        const verts = [v0, v1, v2];

        const globalDof = new Int32Array(nLocal);
        const isDirichlet = new Uint8Array(nLocal);
        const dirichletVal = new Float64Array(nLocal);
        const signs = new Float64Array(nLocal).fill(1);

        // Vertex DOFs (lv): indices 0-2
        for (let k = 0; k < 3; k++) {
            const nf = nodeF[verts[k]];
            if (nf >= 0) {
                globalDof[k] = nf;
            } else if (isCondNode[verts[k]]) {
                globalDof[k] = -1;
                isDirichlet[k] = 1;
                const cg = condNodeGroup ? condNodeGroup[verts[k]] : 1;
                dirichletVal[k] = condPotentials ? condPotentials[cg - 1] : 1.0;
            } else {
                globalDof[k] = -1;
            }
        }

        // Edge midpoint DOFs (le): indices 3-5
        for (let k = 0; k < 3; k++) {
            const eIdx = mesh.triEdges[3*t + k];
            const enf = edgeNodeF[eIdx];
            if (enf >= 0) {
                globalDof[k+3] = edgeMidOff + enf;
            } else if (isCondEdge[eIdx]) {
                globalDof[k+3] = -1;
                isDirichlet[k+3] = 1;
                const eg = condNodeGroup ? condNodeGroup[mesh.edges[2*eIdx]] : 1;
                dirichletVal[k+3] = condPotentials ? condPotentials[eg - 1] : 1.0;
            } else {
                globalDof[k+3] = -1;
            }
        }

        // P3 extra DOFs
        if (order >= 3) {
            // Edge extra DOFs (le3): indices 6-8 (antisymmetric)
            for (let k = 0; k < 3; k++) {
                const eIdx = mesh.triEdges[3*t + k];
                const enf3 = edgeNodeF3 ? edgeNodeF3[eIdx] : -1;
                if (enf3 >= 0) {
                    globalDof[k+6] = edge3Off + enf3;
                    signs[k+6] = mesh.triSigns[3*t + k]; // antisymmetric
                } else if (isCondEdge[eIdx]) {
                    globalDof[k+6] = -1;
                    isDirichlet[k+6] = 1;
                    // le3 on conductor is zero (not the conductor potential)
                    dirichletVal[k+6] = 0;
                } else {
                    globalDof[k+6] = -1;
                }
            }
            // Face node DOF (lf3): index 9
            const fnf = faceNodeF ? faceNodeF[t] : -1;
            globalDof[9] = fnf >= 0 ? faceNodeOff + fnf : -1;
        }

        // Stiffness assembly
        for (let li = 0; li < nLocal; li++) {
            const gi = globalDof[li]; if (gi < 0 && !isDirichlet[li]) continue;
            for (let lj = 0; lj < nLocal; lj++) {
                const v = eps * signs[li] * signs[lj] * Sz[li * nLocal + lj];
                if (v === 0) continue;
                const gj = globalDof[lj];
                if (gi >= 0 && gj >= 0) {
                    Rows.push(gi); Cols.push(gj); Vals.push(v);
                } else if (gi >= 0 && isDirichlet[lj]) {
                    rhs[gi] -= v * dirichletVal[lj];
                }
            }
        }
    }

    // Enrichment assembly
    if (enrichment) {
        const { corners, supports, evalFn, duffyQuadPoints: duffyQP } = enrichment;
        for (let ci = 0; ci < nCorners; ci++) {
            const gc = nFreeFEM + ci;
            const corner = corners[ci];
            const { tris: supportTris, cornerVerts } = supports[ci];

            for (let si = 0; si < supportTris.length; si++) {
                const t = supportTris[si];
                const cornerLocalIdx = cornerVerts[si];
                const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
                const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
                const txs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
                const tys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];
                const verts = [v0, v1, v2];
                const eps = epsMap ? epsMap[t].re : 1.0;
                const useDuffy = cornerLocalIdx >= 0 && duffyQP;
                const qpts = useDuffy ? duffyQP(cornerLocalIdx) : null;
                const nqp = useDuffy ? qpts.length : NQ;

                // Element DOF map
                const gDof = new Int32Array(6);
                const isDir = new Uint8Array(6);
                const dirVal = new Float64Array(6);
                for (let k = 0; k < 3; k++) {
                    const nf = nodeF[verts[k]];
                    if (nf >= 0) gDof[k] = nf;
                    else if (isCondNode[verts[k]]) { gDof[k]=-1; isDir[k]=1; const cg=condNodeGroup?condNodeGroup[verts[k]]:1; dirVal[k]=condPotentials?condPotentials[cg-1]:1.0; }
                    else gDof[k] = -1;
                }
                for (let k = 0; k < 3; k++) {
                    const eIdx = mesh.triEdges[3*t+k];
                    const enf = edgeNodeF[eIdx];
                    if (enf >= 0) gDof[k+3] = nFreeVertexDof + enf;
                    else if (isCondEdge[eIdx]) { gDof[k+3]=-1; isDir[k+3]=1; const eg=condNodeGroup?condNodeGroup[mesh.edges[2*eIdx]]:1; dirVal[k+3]=condPotentials?condPotentials[eg-1]:1.0; }
                    else gDof[k+3] = -1;
                }

                for (let q = 0; q < nqp; q++) {
                    let w, xq, yq;
                    if (useDuffy) { const qp=qpts[q]; w=qp.w; xq=txs[0]*qp.l1+txs[1]*qp.l2+txs[2]*qp.l3; yq=tys[0]*qp.l1+tys[1]*qp.l2+tys[2]*qp.l3; }
                    else { w=QW[q]; xq=txs[0]*QL1[q]+txs[1]*QL2[q]+txs[2]*QL3[q]; yq=tys[0]*QL1[q]+tys[1]*QL2[q]+tys[2]*QL3[q]; }
                    const { phi: psi, dphidx: dpx, dphidy: dpy } = evalFn(xq, yq, corner);
                    if (psi === 0 && dpx === 0 && dpy === 0) continue;
                    const wA = w * Area;
                    // P2 Lagrange gradients
                    const Gx = new Float64Array(6), Gy = new Float64Array(6);
                    for (let k=0;k<3;k++) { const [gx,gy]=lvGrad(coeff,k,xq,yq); Gx[k]=gx; Gy[k]=gy; }
                    for (let k=0;k<3;k++) { const [p,qq]=edgeVerts[k]; const [gx,gy]=leGrad(coeff,p,qq,xq,yq); Gx[k+3]=gx; Gy[k+3]=gy; }
                    // K[gi, gc] = ε·∫∇N_i·∇ψ dA
                    for (let li=0;li<6;li++) {
                        const v = eps*wA*(Gx[li]*dpx + Gy[li]*dpy);
                        if (v === 0) continue;
                        const gi = gDof[li];
                        if (gi >= 0) { Rows.push(gi);Cols.push(gc);Vals.push(v); Rows.push(gc);Cols.push(gi);Vals.push(v); }
                        else if (isDir[li]) { rhs[gc] -= v * dirVal[li]; }
                    }
                    // K[gc, gc] = ε·∫|∇ψ|² dA
                    { const v=eps*wA*(dpx*dpx+dpy*dpy); if(v!==0){Rows.push(gc);Cols.push(gc);Vals.push(v);} }
                    // Cross-coupling between corners
                    for (let cj=ci+1;cj<nCorners;cj++) {
                        const gcj=nFreeFEM+cj;
                        const {dphidx:dx2,dphidy:dy2}=evalFn(xq,yq,corners[cj]);
                        if(dx2===0&&dy2===0) continue;
                        const v=eps*wA*(dpx*dx2+dpy*dy2);
                        if(v!==0){Rows.push(gc);Cols.push(gcj);Vals.push(v); Rows.push(gcj);Cols.push(gc);Vals.push(v);}
                    }
                }
            }
        }
    }

    const csrS = tripletsToCSR(Rows, Cols, Vals, nFreeDof);
    const { x: phiFree, iters, residual } = solveCG(csrS, rhs, nFreeDof);
    //globalThis.__TRI_DEBUG__ && console.log(`  Static CG: ${iters} iters, residual=${residual.toExponential(2)}`);

    // Reconstruct full potential vector
    const phiVertex = new Float64Array(mesh.nNodes);
    for (let n = 0; n < mesh.nNodes; n++) {
        if (isCondNode[n]) {
            const cg = condNodeGroup ? condNodeGroup[n] : 1;
            phiVertex[n] = condPotentials ? condPotentials[cg - 1] : 1.0;
        } else if (nodeF[n] >= 0) phiVertex[n] = phiFree[nodeF[n]];
    }

    const phiEdge = new Float64Array(mesh.nEdges);
    for (let e = 0; e < mesh.nEdges; e++) {
        if (isCondEdge[e]) {
            const eg = condNodeGroup ? condNodeGroup[mesh.edges[2*e]] : 1;
            phiEdge[e] = condPotentials ? condPotentials[eg - 1] : 1.0;
        } else if (edgeNodeF[e] >= 0) phiEdge[e] = phiFree[edgeMidOff + edgeNodeF[e]];
    }

    // P3 extra DOFs
    const phiEdge3 = elemOrder ? new Float64Array(mesh.nEdges) : null;
    const phiFaceNode = elemOrder ? new Float64Array(nTris) : null;
    if (elemOrder) {
        for (let e = 0; e < mesh.nEdges; e++) {
            const enf3 = edgeNodeF3 ? edgeNodeF3[e] : -1;
            if (enf3 >= 0) phiEdge3[e] = phiFree[edge3Off + enf3];
        }
        for (let t = 0; t < nTris; t++) {
            const fnf = faceNodeF ? faceNodeF[t] : -1;
            if (fnf >= 0) phiFaceNode[t] = phiFree[faceNodeOff + fnf];
        }
    }

    const enrichCoeffs = nCorners > 0 ? new Float64Array(phiFree.buffer, nFreeFEM * 8, nCorners) : null;
    return { phiVertex, phiEdge, phiEdge3, phiFaceNode, enrichCoeffs };
}

// Compute energy W = ½∫ε|∇φ|² dA including enrichment cross terms.
// enrichment: optional { corners, evalFn, coeffs } for including singular contribution.
export function computeTriEnergy(phi, mesh, epsMap, enrichment = null, elemOrder = null) {
    const { nodes, tris, nTris } = mesh;
    const { phiVertex, phiEdge, phiEdge3, phiFaceNode } = phi;
    let W = 0;

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const order = elemOrder ? elemOrder[t] : 2;
        const nLocal = order >= 3 ? 10 : 6;
        const Sz = order >= 3
            ? computeTriP3StaticMatrices(nodes, v0, v1, v2).Sz
            : computeTriP2StaticMatrices(nodes, v0, v1, v2).Sz;
        const eps = epsMap ? epsMap[t].re : 1.0;

        const localPhi = new Float64Array(nLocal);
        const verts = [v0, v1, v2];
        for (let k = 0; k < 3; k++) localPhi[k] = phiVertex[verts[k]];
        for (let k = 0; k < 3; k++) localPhi[k+3] = phiEdge[mesh.triEdges[3*t+k]];
        if (order >= 3) {
            const signs = mesh.triSigns;
            for (let k = 0; k < 3; k++) {
                localPhi[k+6] = phiEdge3 ? signs[3*t+k] * phiEdge3[mesh.triEdges[3*t+k]] : 0;
            }
            localPhi[9] = phiFaceNode ? phiFaceNode[t] : 0;
        }

        for (let li = 0; li < nLocal; li++)
            for (let lj = 0; lj < nLocal; lj++)
                W += eps * Sz[li * nLocal + lj] * localPhi[li] * localPhi[lj];
    }

    // Enrichment cross and self terms: ∫ε(2·∇φ_FEM·c∇ψ + |c∇ψ|²) dA
    if (enrichment && enrichment.coeffs) {
        const { corners, evalFn, coeffs, supports, duffyQuadPoints: duffyQP } = enrichment;
        const { triEdges } = mesh;
        const edgeVerts2 = [[0,1],[1,2],[2,0]];

        for (let ci = 0; ci < corners.length; ci++) {
            if (Math.abs(coeffs[ci]) < 1e-30) continue;
            const corner = corners[ci];
            const supportTris = supports[ci].tris;
            const cornerVerts = supports[ci].cornerVerts;

            for (let si = 0; si < supportTris.length; si++) {
                const t = supportTris[si];
                const cornerLocalIdx = cornerVerts[si];
                const v0=tris[3*t],v1=tris[3*t+1],v2=tris[3*t+2];
                const {coeff,Area}=triCoefficients(nodes,v0,v1,v2);
                const txs=[nodes[2*v0],nodes[2*v1],nodes[2*v2]];
                const tys=[nodes[2*v0+1],nodes[2*v1+1],nodes[2*v2+1]];
                const verts2=[v0,v1,v2];
                const eps = epsMap ? epsMap[t].re : 1.0;
                const lp = new Float64Array(6);
                for(let k=0;k<3;k++) lp[k]=phiVertex[verts2[k]];
                for(let k=0;k<3;k++) lp[k+3]=phiEdge[triEdges[3*t+k]];

                const useDuffy = cornerLocalIdx >= 0 && duffyQP;
                const qpts = useDuffy ? duffyQP(cornerLocalIdx) : null;
                const nqp = useDuffy ? qpts.length : NQ;

                for (let q = 0; q < nqp; q++) {
                    let w, xq, yq;
                    if (useDuffy) { const qp=qpts[q]; w=qp.w; xq=txs[0]*qp.l1+txs[1]*qp.l2+txs[2]*qp.l3; yq=tys[0]*qp.l1+tys[1]*qp.l2+tys[2]*qp.l3; }
                    else { w=QW[q]; xq=txs[0]*QL1[q]+txs[1]*QL2[q]+txs[2]*QL3[q]; yq=tys[0]*QL1[q]+tys[1]*QL2[q]+tys[2]*QL3[q]; }
                    const {dphidx:dpx,dphidy:dpy}=evalFn(xq,yq,corner);
                    if(dpx===0&&dpy===0) continue;
                    // ∇φ_FEM at quadrature point
                    let gx=0,gy=0;
                    for(let k=0;k<3;k++){const [dx,dy]=lvGrad(coeff,k,xq,yq);gx+=dx*lp[k];gy+=dy*lp[k];}
                    for(let k=0;k<3;k++){const [p,qq]=edgeVerts2[k];const [dx,dy]=leGrad(coeff,p,qq,xq,yq);gx+=dx*lp[k+3];gy+=dy*lp[k+3];}
                    const wA = w * Area;
                    // 2·c_i·∫ε∇φ_FEM·∇ψ_i dA
                    W += eps * 2 * coeffs[ci] * wA * (gx*dpx + gy*dpy);
                    // c_i²·∫ε|∇ψ_i|² dA
                    W += eps * coeffs[ci]*coeffs[ci] * wA * (dpx*dpx + dpy*dpy);
                    // Cross: c_i·c_j·∫ε∇ψ_i·∇ψ_j dA
                    for (let cj=ci+1;cj<corners.length;cj++) {
                        if(Math.abs(coeffs[cj])<1e-30) continue;
                        const {dphidx:dx2,dphidy:dy2}=evalFn(xq,yq,corners[cj]);
                        if(dx2===0&&dy2===0) continue;
                        W += eps * 2*coeffs[ci]*coeffs[cj] * wA * (dpx*dx2 + dpy*dy2);
                    }
                }
            }
        }
    }

    return 0.5 * W;
}

// --- Static → all DOFs ---
// Compute initial guess from the static potential

export function staticToEdgeDofs(phi, mesh, fm) {
    const { nodes, edges, nEdges, tris, triEdges, triSigns, nTris } = mesh;
    const { edgeF, faceF, nFreeTransverse, nFreeVertexDof, nodeF, edgeNodeF,
            nFreeLongitudinal, edgeF3, edgeNodeF3, faceNodeF,
            nFreeEdgeNodeDof, nFreeEdgeNode3Dof, nFreeFaceNodeDof } = fm;
    const { phiVertex, phiEdge, phiEdge3, phiFaceNode } = phi;
    const N = nFreeTransverse + nFreeLongitudinal;

    const initVec = new Float64Array(N);

    // ne1 DOFs: Whitney DOF = ∫ E·dl = -(φ_q - φ_p) for edge (p→q)
    for (let e = 0; e < nEdges; e++) {
        const ef = edgeF[2*e]; if (ef < 0) continue;
        const n0 = edges[2*e], n1 = edges[2*e+1];
        initVec[ef] = -(phiVertex[n1] - phiVertex[n0]);
    }

    // ne2 DOFs: ∫ E·t*(λ_p-λ_q) dl = (2/3)(φ_p+φ_q) - (4/3)*φ_pq
    for (let e = 0; e < nEdges; e++) {
        const ef = edgeF[2*e+1]; if (ef < 0) continue;
        const n0 = edges[2*e], n1 = edges[2*e+1];
        initVec[ef] = (2/3)*(phiVertex[n0] + phiVertex[n1]) - (4/3)*phiEdge[e];
    }

    // ne3 DOFs: ∫ E·ne3 dl along each edge.
    // For quasi-TEM static approximation: ne3_DOF ≈ (φp - φq) / 3.
    // (Same convention as ne1_DOF = -(φq - φp) = φp - φq for linear φ.)
    if (edgeF3) {
        for (let e = 0; e < nEdges; e++) {
            const ef3 = edgeF3[e]; if (ef3 < 0) continue;
            const n0 = edges[2*e], n1 = edges[2*e+1];
            initVec[ef3] = (phiVertex[n0] - phiVertex[n1]) / 3;
        }
    }
    // Face DOFs (P2 and P3): set to 0 (interior bubbles, small for quasi-TEM)

    const { lzOff, lzEdgeMidOff, lzEdge3Off, lzFaceNodeOff } = getLzOffsets(fm);

    for (let n = 0; n < mesh.nNodes; n++) {
        const nf = nodeF[n];
        if (nf >= 0) initVec[lzOff + nf] = phiVertex[n];
    }
    for (let e = 0; e < nEdges; e++) {
        const enf = edgeNodeF[e];
        if (enf >= 0) initVec[lzEdgeMidOff + enf] = phiEdge[e];
    }
    // P3 extra longitudinal DOFs
    if (edgeNodeF3 && phiEdge3) {
        for (let e = 0; e < nEdges; e++) {
            const enf3 = edgeNodeF3[e];
            if (enf3 >= 0) initVec[lzEdge3Off + enf3] = phiEdge3[e];
        }
    }
    if (faceNodeF && phiFaceNode) {
        for (let t = 0; t < nTris; t++) {
            const fnf = faceNodeF[t];
            if (fnf >= 0) initVec[lzFaceNodeOff + fnf] = phiFaceNode[t];
        }
    }

    return initVec;
}

// --- Triangular adaptive refinement ---

// Select top refineFrac fraction of triangles by metric value for refinement.
export function markTrianglesForRefinement(metric, refineFrac) {
    const nTris = metric.length;
    const sorted = metric.slice().sort((a, b) => b - a);
    const threshIdx = Math.max(0, Math.round(nTris * refineFrac) - 1);
    const threshold = sorted[threshIdx] || 0;
    if (threshold <= 0) return new Uint8Array(nTris); // nothing to refine
    const marked = new Uint8Array(nTris);
    for (let t = 0; t < nTris; t++) {
        if (metric[t] >= threshold) marked[t] = 1;
    }
    return marked;
}

// Longest-edge bisection with green closure (Rivara-style).
// Marked triangles get their longest edge bisected. Propagation ensures
// conformity: if a triangle has a marked non-longest edge, its longest
// edge also gets marked. This preserves the minimum angle guarantee.
export function refineTriMesh(mesh, marked) {
    const { nodes, tris, edges, triEdges, triSigns, nNodes, nTris, nEdges } = mesh;

    // Compute edge lengths
    const edgeLen2 = new Float64Array(nEdges);
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const dx = nodes[2*n1]-nodes[2*n0], dy = nodes[2*n1+1]-nodes[2*n0+1];
        edgeLen2[e] = dx*dx + dy*dy;
    }

    // Find longest edge of each triangle
    const triLongest = new Int32Array(nTris); // local edge index (0,1,2)
    for (let t = 0; t < nTris; t++) {
        let best = 0;
        if (edgeLen2[triEdges[3*t+1]] > edgeLen2[triEdges[3*t+best]]) best = 1;
        if (edgeLen2[triEdges[3*t+2]] > edgeLen2[triEdges[3*t+best]]) best = 2;
        triLongest[t] = best;
    }

    // Minimum edge length guard
    let domainSize = 0;
    for (let n = 0; n < nNodes; n++) {
        domainSize = Math.max(domainSize, Math.abs(nodes[2*n]), Math.abs(nodes[2*n+1]));
    }
    const minEdgeLen2 = (domainSize * 1e-6) ** 2;

    // Step 1: Mark longest edge of each marked triangle
    const edgeMarked = new Uint8Array(nEdges);
    for (let t = 0; t < nTris; t++) {
        if (!marked[t]) continue;
        const e = triEdges[3*t + triLongest[t]];
        if (edgeLen2[e] > minEdgeLen2) edgeMarked[e] = 1;
    }

    // Longest-edge propagation: if a triangle has a marked non-longest edge,
    // also mark its longest edge. This ensures every green-closure triangle
    // has its longest edge bisected, preventing degenerate sub-triangles.
    let propChanged = true;
    while (propChanged) {
        propChanged = false;
        for (let t = 0; t < nTris; t++) {
            const e0m = edgeMarked[triEdges[3*t]];
            const e1m = edgeMarked[triEdges[3*t+1]];
            const e2m = edgeMarked[triEdges[3*t+2]];
            if (!e0m && !e1m && !e2m) continue; // no marked edges
            const longestE = triEdges[3*t + triLongest[t]];
            if (edgeMarked[longestE]) continue; // longest already marked
            if (edgeLen2[longestE] > minEdgeLen2) {
                edgeMarked[longestE] = 1;
                propChanged = true;
            }
        }
    }

    // Create midpoint nodes for marked edges
    const edgeMidNode = new Int32Array(nEdges).fill(-1);
    let newNodeCount = nNodes;
    const newNodeCoords = [];
    for (let e = 0; e < nEdges; e++) {
        if (!edgeMarked[e]) continue;
        const n0 = edges[2*e], n1 = edges[2*e+1];
        edgeMidNode[e] = newNodeCount++;
        newNodeCoords.push((nodes[2*n0]+nodes[2*n1])/2, (nodes[2*n0+1]+nodes[2*n1+1])/2);
    }

    // Build new nodes array
    const newNodes = new Float64Array(2*newNodeCount);
    newNodes.set(nodes);
    for (let i = 0; i < newNodeCoords.length; i++) newNodes[2*nNodes+i] = newNodeCoords[i];

    // Subdivide triangles
    const newTris = [];
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const e0 = triEdges[3*t], e1 = triEdges[3*t+1], e2 = triEdges[3*t+2];
        const m0 = edgeMidNode[e0], m1 = edgeMidNode[e1], m2 = edgeMidNode[e2];

        const nBisected = (m0>=0?1:0) + (m1>=0?1:0) + (m2>=0?1:0);
        if (nBisected === 3) {
            // Red: 4 sub-triangles
            newTris.push(v0, m0, m2);
            newTris.push(m0, v1, m1);
            newTris.push(m2, m1, v2);
            newTris.push(m0, m1, m2);
        } else if (nBisected === 2) {
            // Green closure: 3 sub-triangles
            if (m0 < 0) {
                newTris.push(v0, v1, m1);
                newTris.push(v0, m1, m2);
                newTris.push(m1, v2, m2);
            } else if (m1 < 0) {
                newTris.push(v0, m0, m2);
                newTris.push(m0, v1, v2);
                newTris.push(m0, v2, m2);
            } else {
                newTris.push(v0, m0, v2);
                newTris.push(m0, v1, m1);
                newTris.push(m0, m1, v2);
            }
        } else if (nBisected === 1) {
            // Bisect into 2 sub-triangles
            if (m0 >= 0) {
                newTris.push(v0, m0, v2);
                newTris.push(m0, v1, v2);
            } else if (m1 >= 0) {
                newTris.push(v0, v1, m1);
                newTris.push(v0, m1, v2);
            } else {
                newTris.push(v0, v1, m2);
                newTris.push(m2, v1, v2);
            }
        } else {
            newTris.push(v0, v1, v2);
        }
    }

    const nNewTris = newTris.length / 3;
    const newTriArr = new Int32Array(newTris);

    // Ensure CCW winding
    for (let t = 0; t < nNewTris; t++) {
        const va = newTriArr[3*t], vb = newTriArr[3*t+1], vc = newTriArr[3*t+2];
        const xa = newNodes[2*va], ya = newNodes[2*va+1];
        const xb = newNodes[2*vb], yb = newNodes[2*vb+1];
        const xc = newNodes[2*vc], yc = newNodes[2*vc+1];
        if ((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya) < 0) {
            newTriArr[3*t+1] = vc; newTriArr[3*t+2] = vb;
        }
    }

    // --- Laplacian smoothing ---
    // Phase 1: smooth new midpoints (as before).
    // Phase 2: quality-targeted smoothing of old vertices adjacent to slivers.
    const constraintYs = Object.keys(mesh.constraintYRanges || {}).map(Number);
    const constraintXs = mesh.constraintXs || [];
    // Domain bounds
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    for (let n = 0; n < newNodeCount; n++) {
        xMin = Math.min(xMin, newNodes[2*n]); xMax = Math.max(xMax, newNodes[2*n]);
        yMin = Math.min(yMin, newNodes[2*n+1]); yMax = Math.max(yMax, newNodes[2*n+1]);
    }
    // Classify ALL vertices: 0=pinned, 1=free, 2=along-Y, 3=along-X
    const canMove = new Uint8Array(newNodeCount);
    for (let n = 0; n < newNodeCount; n++) {
        const mx = newNodes[2*n], my = newNodes[2*n+1];
        let onX = 0, onY = 0;
        if (Math.abs(mx - xMin) < 1e-10 || Math.abs(mx - xMax) < 1e-10) onX = 1;
        if (Math.abs(my - yMin) < 1e-10 || Math.abs(my - yMax) < 1e-10) onY = 1;
        for (const cy of constraintYs) { if (Math.abs(my - cy) < 1e-10) onY = 1; }
        for (const cx of constraintXs) { if (Math.abs(mx - cx) < 1e-10) onX = 1; }
        if (onX && onY) canMove[n] = 0;
        else if (onY) canMove[n] = 3;
        else if (onX) canMove[n] = 2;
        else canMove[n] = 1;
    }
    // Build vertex-to-triangle adjacency
    const vtAdj = new Array(newNodeCount);
    for (let n = 0; n < newNodeCount; n++) vtAdj[n] = [];
    for (let t = 0; t < nNewTris; t++) {
        vtAdj[newTriArr[3*t]].push(t);
        vtAdj[newTriArr[3*t+1]].push(t);
        vtAdj[newTriArr[3*t+2]].push(t);
    }

    function triQuality(t) {
        const va = newTriArr[3*t], vb = newTriArr[3*t+1], vc = newTriArr[3*t+2];
        const ax = newNodes[2*va], ay = newNodes[2*va+1];
        const bx = newNodes[2*vb], by = newNodes[2*vb+1];
        const cx = newNodes[2*vc], cy = newNodes[2*vc+1];
        const al = Math.sqrt((bx-cx)**2+(by-cy)**2);
        const bl = Math.sqrt((ax-cx)**2+(ay-cy)**2);
        const cl = Math.sqrt((ax-bx)**2+(ay-by)**2);
        const s = (al+bl+cl)/2;
        const area = Math.abs((bx-ax)*(cy-ay)-(cx-ax)*(by-ay))/2;
        if (area < 1e-30) return 1e10;
        return al*bl*cl/(8*area*(area/s));
    }

    function smoothVertexSet(eligible, maxPasses, clampLo, clampHi) {
        for (let smoothPass = 0; smoothPass < maxPasses; smoothPass++) {
            const sumX = new Float64Array(newNodeCount), sumY = new Float64Array(newNodeCount);
            const cnt = new Int32Array(newNodeCount);
            for (let t = 0; t < nNewTris; t++) {
                const va = newTriArr[3*t], vb = newTriArr[3*t+1], vc = newTriArr[3*t+2];
                for (const [p, q] of [[va,vb],[va,vc],[vb,vc]]) {
                    sumX[p] += newNodes[2*q]; sumY[p] += newNodes[2*q+1]; cnt[p]++;
                    sumX[q] += newNodes[2*p]; sumY[q] += newNodes[2*p+1]; cnt[q]++;
                }
            }
            let nMoved = 0;
            for (let n = 0; n < newNodeCount; n++) {
                if (!eligible[n] || !canMove[n] || cnt[n] === 0) continue;
                const avgX = sumX[n] / cnt[n], avgY = sumY[n] / cnt[n];
                let nx = canMove[n] === 2 ? newNodes[2*n] : avgX;
                let ny = canMove[n] === 3 ? newNodes[2*n+1] : avgY;
                // Don't let free vertices land on constraint lines (breaks constraint edge check)
                if (canMove[n] === 1) {
                    let snap = false;
                    for (const cy of constraintYs) { if (Math.abs(ny - cy) < 1e-10) { snap = true; break; } }
                    if (!snap) for (const cx of constraintXs) { if (Math.abs(nx - cx) < 1e-10) { snap = true; break; } }
                    if (!snap) { if (Math.abs(nx - xMin) < 1e-10 || Math.abs(nx - xMax) < 1e-10) snap = true; }
                    if (!snap) { if (Math.abs(ny - yMin) < 1e-10 || Math.abs(ny - yMax) < 1e-10) snap = true; }
                    if (snap) continue;
                }
                // Clamp constrained vertices to stay between neighbors on constraint line
                if (clampLo && canMove[n] >= 2) {
                    if (canMove[n] === 3) {
                        // Moving along X: clamp X to [lo, hi]
                        nx = Math.max(clampLo[n] + 1e-12, Math.min(clampHi[n] - 1e-12, nx));
                    } else {
                        // Moving along Y: clamp Y to [lo, hi]
                        ny = Math.max(clampLo[n] + 1e-12, Math.min(clampHi[n] - 1e-12, ny));
                    }
                }
                // Check all adjacent triangles remain positive area
                let valid = true;
                const adj = vtAdj[n];
                for (let i = 0; i < adj.length && valid; i++) {
                    const t = adj[i];
                    const va = newTriArr[3*t], vb = newTriArr[3*t+1], vc = newTriArr[3*t+2];
                    const ox = newNodes[2*n], oy = newNodes[2*n+1];
                    newNodes[2*n] = nx; newNodes[2*n+1] = ny;
                    const area = (newNodes[2*vb]-newNodes[2*va])*(newNodes[2*vc+1]-newNodes[2*va+1]) -
                                 (newNodes[2*vc]-newNodes[2*va])*(newNodes[2*vb+1]-newNodes[2*va+1]);
                    newNodes[2*n] = ox; newNodes[2*n+1] = oy;
                    if (area <= 0) valid = false;
                }
                if (valid) { newNodes[2*n] = nx; newNodes[2*n+1] = ny; nMoved++; }
            }
            if (nMoved === 0) break;
        }
    }

    // Phase 1: smooth new midpoints
    const newMidpoints = new Uint8Array(newNodeCount);
    for (let n = nNodes; n < newNodeCount; n++) newMidpoints[n] = 1;
    smoothVertexSet(newMidpoints, 3);

    // --- Iterative quality improvement: alternating edge swaps and targeted smoothing ---
    function isOnConstraint(na, nb) {
        const ax = newNodes[2*na], ay = newNodes[2*na+1];
        const bx = newNodes[2*nb], by = newNodes[2*nb+1];
        for (const cy of constraintYs) {
            if (Math.abs(ay - cy) < 1e-10 && Math.abs(by - cy) < 1e-10) return true;
        }
        for (const cx of constraintXs) {
            if (Math.abs(ax - cx) < 1e-10 && Math.abs(bx - cx) < 1e-10) return true;
        }
        if (Math.abs(ax - bx) < 1e-10 && (Math.abs(ax - xMin) < 1e-10 || Math.abs(ax - xMax) < 1e-10)) return true;
        if (Math.abs(ay - by) < 1e-10 && (Math.abs(ay - yMin) < 1e-10 || Math.abs(ay - yMax) < 1e-10)) return true;
        return false;
    }

    function edgeSwapPass() {
        const tmpEdgeMap = new Map();
        for (let t = 0; t < nNewTris; t++) {
            const localE = [[0,1],[1,2],[2,0]];
            for (let le = 0; le < 3; le++) {
                const na = newTriArr[3*t+localE[le][0]], nb = newTriArr[3*t+localE[le][1]];
                const n0 = Math.min(na, nb), n1 = Math.max(na, nb);
                const key = n0 * newNodeCount + n1;
                const entry = tmpEdgeMap.get(key);
                if (entry) entry.push(t);
                else tmpEdgeMap.set(key, [t]);
            }
        }
        let nSwaps = 0;
        for (const [key, tris2] of tmpEdgeMap) {
            if (tris2.length !== 2) continue;
            const [t0, t1] = tris2;
            const n0 = Math.floor(key / newNodeCount), n1 = key % newNodeCount;
            if (isOnConstraint(n0, n1)) continue;
            let opp0 = -1, opp1 = -1;
            for (let k = 0; k < 3; k++) {
                if (newTriArr[3*t0+k] !== n0 && newTriArr[3*t0+k] !== n1) opp0 = newTriArr[3*t0+k];
                if (newTriArr[3*t1+k] !== n0 && newTriArr[3*t1+k] !== n1) opp1 = newTriArr[3*t1+k];
            }
            if (opp0 < 0 || opp1 < 0 || opp0 === opp1) continue;
            const qBefore = Math.max(triQuality(t0), triQuality(t1));
            if (qBefore < 10) continue; // only target actual slivers
            // Convexity check: n0 and n1 must be on opposite sides of opp0-opp1
            const cross0 = (newNodes[2*opp1]-newNodes[2*opp0])*(newNodes[2*n0+1]-newNodes[2*opp0+1])
                         - (newNodes[2*n0]-newNodes[2*opp0])*(newNodes[2*opp1+1]-newNodes[2*opp0+1]);
            const cross1 = (newNodes[2*opp1]-newNodes[2*opp0])*(newNodes[2*n1+1]-newNodes[2*opp0+1])
                         - (newNodes[2*n1]-newNodes[2*opp0])*(newNodes[2*opp1+1]-newNodes[2*opp0+1]);
            if (cross0 * cross1 >= 0) continue; // not convex
            // Assign CCW winding based on which side each vertex is on
            const s0 = [newTriArr[3*t0], newTriArr[3*t0+1], newTriArr[3*t0+2]];
            const s1 = [newTriArr[3*t1], newTriArr[3*t1+1], newTriArr[3*t1+2]];
            if (cross0 > 0) {
                // n0 left of opp0→opp1: [opp0, opp1, n0] is CCW
                newTriArr[3*t0] = opp0; newTriArr[3*t0+1] = opp1; newTriArr[3*t0+2] = n0;
                // n1 right: [opp0, n1, opp1] is CCW
                newTriArr[3*t1] = opp0; newTriArr[3*t1+1] = n1; newTriArr[3*t1+2] = opp1;
            } else {
                // n0 right: [opp0, n0, opp1] is CCW
                newTriArr[3*t0] = opp0; newTriArr[3*t0+1] = n0; newTriArr[3*t0+2] = opp1;
                // n1 left: [opp0, opp1, n1] is CCW
                newTriArr[3*t1] = opp0; newTriArr[3*t1+1] = opp1; newTriArr[3*t1+2] = n1;
            }
            const qAfter = Math.max(triQuality(t0), triQuality(t1));
            if (qAfter >= qBefore) {
                // Revert
                newTriArr[3*t0] = s0[0]; newTriArr[3*t0+1] = s0[1]; newTriArr[3*t0+2] = s0[2];
                newTriArr[3*t1] = s1[0]; newTriArr[3*t1+1] = s1[1]; newTriArr[3*t1+2] = s1[2];
            } else {
                nSwaps++;
            }
        }
        return nSwaps;
    }

    // Phase 2: alternating edge swaps + quality-targeted smoothing
    for (let outerPass = 0; outerPass < 3; outerPass++) {
        // Edge swaps
        for (let sp = 0; sp < 5; sp++) { if (edgeSwapPass() === 0) break; }
        // Rebuild vtAdj after swaps
        for (let n = 0; n < newNodeCount; n++) vtAdj[n] = [];
        for (let t = 0; t < nNewTris; t++) {
            vtAdj[newTriArr[3*t]].push(t);
            vtAdj[newTriArr[3*t+1]].push(t);
            vtAdj[newTriArr[3*t+2]].push(t);
        }
        // Quality-targeted smoothing of vertices near slivers
        const Q_THRESH = 10;
        const sliverAdj = new Uint8Array(newNodeCount);
        for (let t = 0; t < nNewTris; t++) {
            if (triQuality(t) > Q_THRESH) {
                for (let k = 0; k < 3; k++) {
                    const v = newTriArr[3*t+k];
                    if (canMove[v] >= 1) sliverAdj[v] = 1;
                }
            }
        }
        // Build constraint-line neighbor bounds for constrained vertices
        // canMove=2: on X boundary, moves along Y; canMove=3: on Y boundary, moves along X
        const clampLo = new Float64Array(newNodeCount).fill(-Infinity);
        const clampHi = new Float64Array(newNodeCount).fill(Infinity);
        for (let n = 0; n < newNodeCount; n++) {
            if (!sliverAdj[n] || canMove[n] < 2) continue;
            const mx = newNodes[2*n], my = newNodes[2*n+1];
            // Find neighbors on same constraint line via mesh edges
            for (const t of vtAdj[n]) {
                for (let k = 0; k < 3; k++) {
                    const nb = newTriArr[3*t+k];
                    if (nb === n) continue;
                    const bx = newNodes[2*nb], by = newNodes[2*nb+1];
                    if (canMove[n] === 3) {
                        // On Y constraint, moving along X: find neighbors at same Y
                        if (Math.abs(by - my) < 1e-10) {
                            if (bx < mx) clampLo[n] = Math.max(clampLo[n], bx);
                            if (bx > mx) clampHi[n] = Math.min(clampHi[n], bx);
                        }
                    } else {
                        // On X constraint, moving along Y: find neighbors at same X
                        if (Math.abs(bx - mx) < 1e-10) {
                            if (by < my) clampLo[n] = Math.max(clampLo[n], by);
                            if (by > my) clampHi[n] = Math.min(clampHi[n], by);
                        }
                    }
                }
            }
        }
        smoothVertexSet(sliverAdj, 5, clampLo, clampHi);
    }
    // Final edge swap pass after last smoothing
    for (let sp = 0; sp < 5; sp++) { if (edgeSwapPass() === 0) break; }

    // --- Rebuild edge list and tri-edge mapping ---
    const edgeMap2 = new Map();
    const edgeList2 = [];
    const newTriEdges = new Int32Array(3*nNewTris);
    const newTriSigns = new Int8Array(3*nNewTris);

    for (let t = 0; t < nNewTris; t++) {
        const localEdges = [[0,1],[1,2],[2,0]];
        for (let le = 0; le < 3; le++) {
            const na = newTriArr[3*t+localEdges[le][0]];
            const nb = newTriArr[3*t+localEdges[le][1]];
            const n0 = Math.min(na,nb), n1 = Math.max(na,nb);
            const key = n0 * newNodeCount + n1;
            let eIdx;
            if (edgeMap2.has(key)) {
                eIdx = edgeMap2.get(key);
            } else {
                eIdx = edgeList2.length;
                edgeList2.push([n0,n1]);
                edgeMap2.set(key, eIdx);
            }
            newTriEdges[3*t+le] = eIdx;
            newTriSigns[3*t+le] = (na === n0) ? 1 : -1;
        }
    }

    const nNewEdges = edgeList2.length;
    const newEdgeArr = new Int32Array(2*nNewEdges);
    for (let e = 0; e < nNewEdges; e++) {
        newEdgeArr[2*e] = edgeList2[e][0]; newEdgeArr[2*e+1] = edgeList2[e][1];
    }

    return {
        nodes: newNodes, tris: newTriArr, edges: newEdgeArr,
        triEdges: newTriEdges, triSigns: newTriSigns,
        nNodes: newNodeCount, nTris: nNewTris, nEdges: nNewEdges,
        Nx: 0, Ny: 0,
        constraintYRanges: mesh.constraintYRanges || {},
        constraintXs: mesh.constraintXs || []
    };
}
