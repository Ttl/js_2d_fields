// Triangular P2 Nedelec FEM — generic mesh, basis, assembly, and solver functions

import { tripletsToCSR, solveCG, triQuality as triQualityXY } from './fem_core.js';

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

    return { ARe, AIm, BRe, BIm, Area, coeff, Dzz1, Dzz2Re, Dzz2Im, Att, Dtt, Dzt, BttRe, BttIm };
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

// --- Longitudinal DOF offsets ---
// Returns the global index offsets for each longitudinal DOF group.
export function getLzOffsets(fm) {
    const lzOff = fm.nFreeTransverse;
    const lzEdgeMidOff = lzOff + fm.nFreeVertexDof;
    return { lzOff, lzEdgeMidOff };
}

// --- Freedom map ---
// Global DOF layout (P2):
// Transverse: [0, 2*nEdges) edge DOFs + [2*nEdges, 2*nEdges+2*nTris) face DOFs
// Longitudinal: [n_xy, n_xy+nNodes) vertex DOFs + [n_xy+nNodes, n_xy+nNodes+nEdges) edge DOFs

export function buildTriFreedomMap(mesh, condRect, abc) {
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

    // An edge is a conductor (PEC) edge only if BOTH endpoints are on metal AND
    // its midpoint lies on/in a conductor rect. Endpoints alone are not enough:
    // an air edge bridging TWO conductors across a one-element gap would be
    // marked PEC and short the gap. The midpoint test also gives the edge its
    // own conductor group (an endpoint shared between touching rects can carry
    // the wrong group).
    const isCondEdge = new Uint8Array(nEdges);
    const condEdgeGroup = new Int8Array(nEdges);
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        if (!(isCondNode[n0] && isCondNode[n1])) continue;
        const xm = (nodes[2*n0] + nodes[2*n1]) / 2;
        const ym = (nodes[2*n0+1] + nodes[2*n1+1]) / 2;
        for (let ci = 0; ci < condRects.length; ci++) {
            const cr = condRects[ci];
            if (xm >= cr.xmin - TOL && xm <= cr.xmax + TOL &&
                ym >= cr.ymin - TOL && ym <= cr.ymax + TOL) {
                isCondEdge[e] = 1;
                condEdgeGroup[e] = ci + 1;
                break;
            }
        }
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

    // --- Transverse edge DOFs (2 per edge: ne1, ne2) ---
    const edgeF = new Int32Array(2 * nEdges).fill(-1);
    let nFreeEdgeDof = 0;
    for (let e = 0; e < nEdges; e++) {
        if (isEdgePEC(e)) continue;
        edgeF[2*e] = nFreeEdgeDof++;     // ne1
        edgeF[2*e+1] = nFreeEdgeDof++;   // ne2
    }

    // --- Transverse face DOFs (2 per tri: nf1, nf2) ---
    const faceF = new Int32Array(2 * nTris).fill(-1);
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
        if (isCondEdge[e]) continue;
        if (!abc.left && Math.abs(xm-xmin)<TOL) continue;
        if (!abc.right && Math.abs(xm-xmax)<TOL) continue;
        if (!abc.top && Math.abs(ym-ymax)<TOL) continue;
        edgeNodeF[e] = nFreeEdgeNodeDof++;
    }

    const nFreeLongitudinal = nFreeVertexDof + nFreeEdgeNodeDof;

    return {
        edgeF, faceF, nodeF, edgeNodeF,
        nFreeEdgeDof, nFreeFaceDof, nFreeTransverse,
        nFreeVertexDof, nFreeEdgeNodeDof, nFreeLongitudinal,
        isCondNode, isCondEdge, condNodeGroup, condEdgeGroup,
    };
}

// --- Assembly ---

// Decomposed assembly. The Lee–Jin system is AFFINE in k² with a Robin (ABC)
// boundary term LINEAR in k₀ = √k²:
//   A(k) = A0 + k²·A1 + j·k0·Ar        (A0 = Att block, A1 = −Btt, Ar = ABC template)
//   B(k) = B0 + k²·B1                  (B0 = Dtt/Dzt/Dzz1 blocks, B1 = −Dzz2)
// This builds the five value templates ONCE on fixed CSR patterns (one pattern
// for A including the ABC entries, one for B); femFromDecomposition() then
// produces the matrices for any frequency in O(nnz). Per-frequency quadrature
// reassembly used to be a dominant sweep cost.
export function assembleTriFEMDecomposed(mesh, fm, epsMap, abc, condRect) {
    const { tris, edges, triEdges, triSigns, nTris, nEdges } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF } = fm;
    const N = fm.nFreeTransverse + fm.nFreeLongitudinal;
    const nodes = mesh.nodes;
    const { lzOff, lzEdgeMidOff } = getLzOffsets(fm);
    const nLocal = 14;

    // Pre-allocate COO arrays: max 196 nonzeros per P2 element (14×14)
    const maxNnz = nTris * 196;
    const eAr = new Int32Array(maxNnz), eAc = new Int32Array(maxNnz);
    const eA0 = new Float64Array(maxNnz);
    const eA1Re = new Float64Array(maxNnz), eA1Im = new Float64Array(maxNnz);
    const eBr = new Int32Array(maxNnz), eBc = new Int32Array(maxNnz);
    const eB0 = new Float64Array(maxNnz);
    const eB1Re = new Float64Array(maxNnz), eB1Im = new Float64Array(maxNnz);
    let aNnz = 0, bNnz = 0;

    // Element templates in the 14×14 local DOF layout
    const A0el = new Float64Array(196);
    const A1elRe = new Float64Array(196), A1elIm = new Float64Array(196);
    const B0el = new Float64Array(196);
    const B1elRe = new Float64Array(196), B1elIm = new Float64Array(196);
    const globalDof = new Int32Array(nLocal);
    const signs = new Float64Array(nLocal);

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const eps = epsMap[t];
        const m = computeTriP2Matrices(nodes, v0, v1, v2, eps.re, eps.im, 0);

        A0el.fill(0); A1elRe.fill(0); A1elIm.fill(0);
        B0el.fill(0); B1elRe.fill(0); B1elIm.fill(0);
        for (let i = 0; i < 8; i++) {
            for (let j = 0; j < 8; j++) {
                A0el[i*14+j] = m.Att[i*8+j];
                A1elRe[i*14+j] = -m.BttRe[i*8+j];
                A1elIm[i*14+j] = -m.BttIm[i*8+j];
                B0el[i*14+j] = m.Dtt[i*8+j];
            }
        }
        for (let i = 0; i < 6; i++) {
            for (let j = 0; j < 8; j++) {
                B0el[(i+8)*14+j] = m.Dzt[i*8+j];     // B[8:14,0:8] = Dzt
                B0el[j*14+(i+8)] = m.Dzt[i*8+j];     // B[0:8,8:14] = Dzt^T
            }
            for (let j = 0; j < 6; j++) {
                B0el[(i+8)*14+(j+8)] = m.Dzz1[i*6+j];
                B1elRe[(i+8)*14+(j+8)] = -m.Dzz2Re[i*6+j];
                B1elIm[(i+8)*14+(j+8)] = -m.Dzz2Im[i*6+j];
            }
        }

        // Local-to-global DOF mapping
        signs.fill(1);
        for (let le = 0; le < 3; le++) {
            const eIdx = triEdges[3*t + le];
            globalDof[le] = edgeF[2*eIdx];
            signs[le] = triSigns[3*t + le];
            globalDof[le+4] = edgeF[2*eIdx+1];
            const enf = edgeNodeF[eIdx];
            globalDof[le+11] = enf >= 0 ? lzEdgeMidOff + enf : -1;
        }
        globalDof[3] = faceF[2*t];
        globalDof[7] = faceF[2*t+1];
        const verts = [v0, v1, v2];
        for (let k = 0; k < 3; k++) {
            const nf = nodeF[verts[k]];
            globalDof[8 + k] = nf >= 0 ? lzOff + nf : -1;
        }

        for (let li = 0; li < nLocal; li++) {
            const gi = globalDof[li]; if (gi < 0) continue;
            for (let lj = 0; lj < nLocal; lj++) {
                const gj = globalDof[lj]; if (gj < 0) continue;
                const s = signs[li] * signs[lj];
                const idx = li * nLocal + lj;
                const a0 = s * A0el[idx], a1re = s * A1elRe[idx], a1im = s * A1elIm[idx];
                if (a0 !== 0 || a1re !== 0 || a1im !== 0) {
                    eAr[aNnz] = gi; eAc[aNnz] = gj;
                    eA0[aNnz] = a0; eA1Re[aNnz] = a1re; eA1Im[aNnz] = a1im; aNnz++;
                }
                const b0 = s * B0el[idx], b1re = s * B1elRe[idx], b1im = s * B1elIm[idx];
                if (b0 !== 0 || b1re !== 0 || b1im !== 0) {
                    eBr[bNnz] = gi; eBc[bNnz] = gj;
                    eB0[bNnz] = b0; eB1Re[bNnz] = b1re; eB1Im[bNnz] = b1im; bNnz++;
                }
            }
        }
    }

    // Robin ABC entries: A += j·k₀·(1/L or 1/(3L)) on boundary edge DOFs — a pure
    // per-unit-k₀ imaginary template. On edge (p,q) parameterized by s∈[0,1] the
    // transverse tangential basis is ne1_tang(s) = 1/L (constant), ne2_tang(s) =
    // (1/L)(1−2s), giving self-integrals ∫ne1²dl = 1/L, ∫ne2²dl = 1/(3L),
    // ∫ne1·ne2 dl = 0. (Strict === true: 'pmc' walls skip the Robin terms. The P2
    // NODAL Robin on Ez is deliberately absent — forcing ∂Ez/∂n + jk₀Ez = 0 distorts
    // quasi-TEM eigenvectors; the transverse ABC alone suppresses cavity modes.)
    const rR = [], rC = [], rV = [];
    if (abc.top === true || abc.left === true || abc.right === true || abc.bottom === true) {
        const ymax = condRect.ymax_domain, ymin = condRect.ymin_domain;
        const xmin = condRect.xmin_domain, xmax = condRect.xmax_domain;
        const TOL = 1e-12;
        for (let e = 0; e < nEdges; e++) {
            const n0 = edges[2*e], n1 = edges[2*e+1];
            const x0 = nodes[2*n0], y0 = nodes[2*n0+1];
            const x1 = nodes[2*n1], y1 = nodes[2*n1+1];
            let isBoundary = false;
            if (abc.top === true && Math.abs(y0 - ymax) < TOL && Math.abs(y1 - ymax) < TOL) isBoundary = true;
            if (abc.bottom === true && Math.abs(y0 - ymin) < TOL && Math.abs(y1 - ymin) < TOL) isBoundary = true;
            if (abc.left === true && Math.abs(x0 - xmin) < TOL && Math.abs(x1 - xmin) < TOL) isBoundary = true;
            if (abc.right === true && Math.abs(x0 - xmax) < TOL && Math.abs(x1 - xmax) < TOL) isBoundary = true;
            if (!isBoundary) continue;
            const L = Math.sqrt((x1-x0)**2 + (y1-y0)**2);
            const ef1 = edgeF[2*e];
            if (ef1 >= 0) { rR.push(ef1); rC.push(ef1); rV.push(1 / L); }
            const ef2 = edgeF[2*e+1];
            if (ef2 >= 0) { rR.push(ef2); rC.push(ef2); rV.push(1 / (3 * L)); }
        }
    }

    // Combined A pattern = element entries + ABC entries. The SAME (row, col)
    // sequence is converted for each value template, so the dedup produces an
    // identical pattern and the templates line up entry-for-entry.
    const nAe = aNnz, nR = rR.length, nA = nAe + nR;
    const rowsA = new Int32Array(nA), colsA = new Int32Array(nA);
    rowsA.set(eAr.subarray(0, nAe)); colsA.set(eAc.subarray(0, nAe));
    const v0A = new Float64Array(nA), v1ReA = new Float64Array(nA), v1ImA = new Float64Array(nA);
    const vrA = new Float64Array(nA);
    v0A.set(eA0.subarray(0, nAe)); v1ReA.set(eA1Re.subarray(0, nAe)); v1ImA.set(eA1Im.subarray(0, nAe));
    for (let i = 0; i < nR; i++) {
        rowsA[nAe + i] = rR[i]; colsA[nAe + i] = rC[i];
        vrA[nAe + i] = rV[i];
    }
    const csrA0 = tripletsToCSR(rowsA, colsA, v0A, N);
    const csrA1 = tripletsToCSR(rowsA, colsA, v1ReA, N, v1ImA);
    const csrAr = tripletsToCSR(rowsA, colsA, new Float64Array(nA), N, vrA);
    const csrB0 = tripletsToCSR(eBr.subarray(0, bNnz), eBc.subarray(0, bNnz), eB0.subarray(0, bNnz), N);
    const csrB1 = tripletsToCSR(eBr.subarray(0, bNnz), eBc.subarray(0, bNnz), eB1Re.subarray(0, bNnz), N, eB1Im.subarray(0, bNnz));

    return {
        N,
        A: { rowPtr: csrA0.rowPtr, colIdx: csrA0.colIdx,
             v0Re: csrA0.valRe, v1Re: csrA1.valRe, v1Im: csrA1.valIm, vrIm: csrAr.valIm },
        B: { rowPtr: csrB0.rowPtr, colIdx: csrB0.colIdx,
             v0Re: csrB0.valRe, v1Re: csrB1.valRe, v1Im: csrB1.valIm },
    };
}

// Combine a decomposition into the concrete A(k), B(k) CSR pair for one k².
export function femFromDecomposition(dec, k2) {
    const k0 = Math.sqrt(k2);
    const A = dec.A, nA = A.v0Re.length;
    const aRe = new Float64Array(nA), aIm = new Float64Array(nA);
    for (let i = 0; i < nA; i++) {
        aRe[i] = A.v0Re[i] + k2 * A.v1Re[i];
        aIm[i] = k2 * A.v1Im[i] + k0 * A.vrIm[i];
    }
    const B = dec.B, nB = B.v0Re.length;
    const bRe = new Float64Array(nB), bIm = new Float64Array(nB);
    for (let i = 0; i < nB; i++) {
        bRe[i] = B.v0Re[i] + k2 * B.v1Re[i];
        bIm[i] = k2 * B.v1Im[i];
    }
    return {
        csrA: { rowPtr: A.rowPtr, colIdx: A.colIdx, valRe: aRe, valIm: aIm },
        csrB: { rowPtr: B.rowPtr, colIdx: B.colIdx, valRe: bRe, valIm: bIm },
        N: dec.N,
    };
}

// Concrete A(k), B(k) CSR pair at a single k² — a thin wrapper over the decomposed
// assembly (verified equivalent to machine precision), for one-shot callers (the
// mode viewer, tests) that don't reuse the templates across frequencies.
export function assembleTriFEM(mesh, fm, k2, epsMap, abc, condRect) {
    return femFromDecomposition(assembleTriFEMDecomposed(mesh, fm, epsMap, abc, condRect), k2);
}

// --- P2 Static solver ---
// Uses 6 nodal DOFs per triangle: 3 vertex (lv) + 3 edge midpoint (le)

// Returns { phiVertex, phiEdge }.
export function solveTriStatic(mesh, fm, epsMap, condPotentials = null) {
    // condPotentials: array of potentials per conductor group [V1, V2, ...].
    // If null, all conductors get V=1.0.
    const { nodes, tris, nTris } = mesh;
    const { nodeF, edgeNodeF, nFreeVertexDof, nFreeEdgeNodeDof,
            isCondNode, isCondEdge, condNodeGroup, condEdgeGroup } = fm;
    const nFreeDof = nFreeVertexDof + nFreeEdgeNodeDof;

    // Longitudinal DOF offsets (within the static system, lzOff=0)
    const edgeMidOff = nFreeVertexDof;
    // Note: these are offsets within the longitudinal block (no lzOff prefix)
    // because the static solver only has longitudinal DOFs.

    const Rows = [], Cols = [], Vals = [];
    const rhs = new Float64Array(nFreeDof);

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const nLocal = 6;
        const Sz = computeTriP2StaticMatrices(nodes, v0, v1, v2).Sz;
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
                const eg = condEdgeGroup ? condEdgeGroup[eIdx]
                         : (condNodeGroup ? condNodeGroup[mesh.edges[2*eIdx]] : 1);
                dirichletVal[k+3] = condPotentials ? condPotentials[eg - 1] : 1.0;
            } else {
                globalDof[k+3] = -1;
            }
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


    const csrS = tripletsToCSR(Rows, Cols, Vals, nFreeDof);
    const { x: phiFree, iters, residual, converged } = solveCG(csrS, rhs, nFreeDof);
    if (!converged) {
        console.warn(`Static CG did not converge in ${iters} iterations ` +
            `(residual ${residual.toExponential(2)}) — static C/Z0 may be inaccurate.`);
    }

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
            const eg = condEdgeGroup ? condEdgeGroup[e]
                     : (condNodeGroup ? condNodeGroup[mesh.edges[2*e]] : 1);
            phiEdge[e] = condPotentials ? condPotentials[eg - 1] : 1.0;
        } else if (edgeNodeF[e] >= 0) phiEdge[e] = phiFree[edgeMidOff + edgeNodeF[e]];
    }

    return { phiVertex, phiEdge };
}

// Compute energy W = ½∫ε|∇φ|² dA.
export function computeTriEnergy(phi, mesh, epsMap) {
    const { nodes, tris, nTris } = mesh;
    const { phiVertex, phiEdge } = phi;
    let W = 0;

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const nLocal = 6;
        const Sz = computeTriP2StaticMatrices(nodes, v0, v1, v2).Sz;
        const eps = epsMap ? epsMap[t].re : 1.0;

        const localPhi = new Float64Array(nLocal);
        const verts = [v0, v1, v2];
        for (let k = 0; k < 3; k++) localPhi[k] = phiVertex[verts[k]];
        for (let k = 0; k < 3; k++) localPhi[k+3] = phiEdge[mesh.triEdges[3*t+k]];

        for (let li = 0; li < nLocal; li++)
            for (let lj = 0; lj < nLocal; lj++)
                W += eps * Sz[li * nLocal + lj] * localPhi[li] * localPhi[lj];
    }

    return 0.5 * W;
}

// --- Static → all DOFs ---
// Compute initial guess from the static potential

export function staticToEdgeDofs(phi, mesh, fm) {
    const { edges, nEdges } = mesh;
    const { edgeF, nFreeTransverse, nodeF, edgeNodeF,
            nFreeLongitudinal } = fm;
    const { phiVertex, phiEdge } = phi;
    const N = nFreeTransverse + nFreeLongitudinal;

    const initVec = new Float64Array(N);

    // ne1 DOFs: Whitney DOF = ∫ E·dl = -(φ_q - φ_p) for edge (p→q)
    for (let e = 0; e < nEdges; e++) {
        const ef = edgeF[2*e]; if (ef < 0) continue;
        const n0 = edges[2*e], n1 = edges[2*e+1];
        initVec[ef] = -(phiVertex[n1] - phiVertex[n0]);
    }

    // ne2 DOFs: the ne2 basis is ne1·(λ_p−λ_q), whose tangential trace gives
    // ∫ ne2·t (λ_p−λ_q) dl = 1/3 — so the DOF dual to it is
    // 3·∫ E·t (λ_p−λ_q) dl = 3·[(2/3)(φ_p+φ_q) − (4/3)φ_pq].
    for (let e = 0; e < nEdges; e++) {
        const ef = edgeF[2*e+1]; if (ef < 0) continue;
        const n0 = edges[2*e], n1 = edges[2*e+1];
        initVec[ef] = 2*(phiVertex[n0] + phiVertex[n1]) - 4*phiEdge[e];
    }
    // Face DOFs: set to 0 (interior bubbles, small for quasi-TEM)

    const { lzOff, lzEdgeMidOff } = getLzOffsets(fm);

    for (let n = 0; n < mesh.nNodes; n++) {
        const nf = nodeF[n];
        if (nf >= 0) initVec[lzOff + nf] = phiVertex[n];
    }
    for (let e = 0; e < nEdges; e++) {
        const enf = edgeNodeF[e];
        if (enf >= 0) initVec[lzEdgeMidOff + enf] = phiEdge[e];
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
    const cyRanges = mesh.constraintYRanges || {};
    const cxRanges = mesh.constraintXRanges || null;
    const constraintYs = Object.keys(cyRanges).map(Number);
    const constraintXs = cxRanges ? Object.keys(cxRanges).map(Number) : (mesh.constraintXs || []);
    // A point is on a constraint line only within the line's extent (a conductor
    // side wall constrains only over the conductor's height, not the whole
    // domain column at that x). Ranges absent → whole line (legacy meshes).
    const onYLine = (x, y) => {
        for (const cy of constraintYs) {
            if (Math.abs(y - cy) >= 1e-10) continue;
            const r = cyRanges[cy];
            if (!r || (x >= r[0] - 1e-10 && x <= r[1] + 1e-10)) return true;
        }
        return false;
    };
    const onXLine = (x, y) => {
        for (const cx of constraintXs) {
            if (Math.abs(x - cx) >= 1e-10) continue;
            const r = cxRanges ? cxRanges[cx] : null;
            if (!r || (y >= r[0] - 1e-10 && y <= r[1] + 1e-10)) return true;
        }
        return false;
    };
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
        if (onYLine(mx, my)) onY = 1;
        if (onXLine(mx, my)) onX = 1;
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
        return triQualityXY(
            newNodes[2*va], newNodes[2*va+1],
            newNodes[2*vb], newNodes[2*vb+1],
            newNodes[2*vc], newNodes[2*vc+1]);
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
                    let snap = onYLine(nx, ny) || onXLine(nx, ny);
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
        if (onYLine(ax, ay) && onYLine(bx, by) && Math.abs(ay - by) < 1e-10) return true;
        if (onXLine(ax, ay) && onXLine(bx, by) && Math.abs(ax - bx) < 1e-10) return true;
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
        // A swap changes BOTH triangles' edge sets, so map entries built at pass
        // start go stale: a later entry referencing a swapped triangle would
        // rewrite triangles that no longer share the edge (mesh corruption).
        // Skip entries touching swapped triangles; the next pass rebuilds the map.
        const dirty = new Uint8Array(nNewTris);
        for (const [key, tris2] of tmpEdgeMap) {
            if (tris2.length !== 2) continue;
            const [t0, t1] = tris2;
            if (dirty[t0] || dirty[t1]) continue;
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
                dirty[t0] = dirty[t1] = 1;
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
        constraintXRanges: mesh.constraintXRanges || null,
        constraintXs: mesh.constraintXs || []
    };
}
