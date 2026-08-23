// Triangular P2 Nedelec FEM — generic mesh, basis, assembly, and solver functions

import { tripletsToCSR, tripletsToCSRMulti, solveCG, triQuality as triQualityXY } from './fem_core.js';
import { shapeContains, REL_SHAPE_TOL } from '../shapes.js';

// 6-point Gauss quadrature on triangle (barycentric)

// Degree-4 exact. Full double precision on purpose. At the 8-digit truncation
// the barycentric triples summed to 0.99999999, so the quadrature point was
// slightly off the triangle resulting in ~1e-5 relative on the element matrices
// of a small element far from the origin (exactly the mesh a PCB cross-section
// produces).
export const QW = [0.223381589678011, 0.223381589678011, 0.223381589678011,
                   0.109951743655322, 0.109951743655322, 0.109951743655322];
export const QL1 = [0.108103018168070, 0.445948490915965, 0.445948490915965,
                    0.816847572980459, 0.091576213509771, 0.091576213509771];
export const QL2 = [0.445948490915965, 0.445948490915965, 0.108103018168070,
                    0.091576213509771, 0.091576213509771, 0.816847572980459];
export const QL3 = [0.445948490915965, 0.108103018168070, 0.445948490915965,
                    0.091576213509771, 0.816847572980459, 0.091576213509771];
export const NQ = 6;


// Barycentric coefficients

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

// P2 Nedelec basis functions (using normalized coeff)

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

// P2 element matrices for the full-wave pencil

// Sub-blocks in the 8 transverse (Nedelec) + 6 longitudinal (nodal) local layout:
//   Att  [8x8]  ∫curlW_i·curlW_j          (ε-free)
//   Btt  [8x8]  ε·∫W_i·W_j                 = ε · Wtt
//   Dtt  [8x8]  ∫W_i·W_j with the face-face entries zeroed
//   Dzt  [6x8]  ∫∇Lz_i·W_j                 (ε-free)
//   Dzz1 [6x6]  ∫∇Lz_i·∇Lz_j               = the P2 stiffness
//   Dzz2 [6x6]  ε·∫Lz_i·Lz_j               = ε · Area · P2_MASS
// Only one integral (Wtt) carries ε, and it does so as a scalar factor, so Btt/Dtt
// share it and Dzz1/Dzz2 come from the closed forms above without quadrature.
//
// The quadrature runs on barycentric coordinates directly. On a straight-sided
// triangle λ_i at quadrature point q is the rule's barycentric weight, so no
// physical point or affine λ evaluation (which amplifies the a_i/(2A) term for a
// small element far from the origin), and no basis-function call. Every basis and
// curl below is written out in terms of (λ, b, c).
//
// The returned object and every array in it are resued across calls (this runs
// once per triangle per assembly): consume the values before the next call.

const _Wx = new Float64Array(8), _Wy = new Float64Array(8), _cW = new Float64Array(8);
const _Gx = new Float64Array(6), _Gy = new Float64Array(6);
const _Att = new Float64Array(64), _Wtt = new Float64Array(64), _Dtt = new Float64Array(64);
const _BttRe = new Float64Array(64), _BttIm = new Float64Array(64);
const _Dzt = new Float64Array(48);
const _Dzz1 = new Float64Array(36), _Dzz2Re = new Float64Array(36), _Dzz2Im = new Float64Array(36);
const _p2out = { Area: 0, Att: _Att, Dtt: _Dtt, Dzt: _Dzt, Dzz1: _Dzz1,
                 Dzz2Re: _Dzz2Re, Dzz2Im: _Dzz2Im, BttRe: _BttRe, BttIm: _BttIm };

export function computeTriP2Matrices(nodes, v0, v1, v2, epsRe, epsIm) {
    const x0 = nodes[2*v0], y0 = nodes[2*v0+1];
    const x1 = nodes[2*v1], y1 = nodes[2*v1+1];
    const x2 = nodes[2*v2], y2 = nodes[2*v2+1];
    const sA = 0.5 * ((x0 - x2) * (y1 - y0) - (x0 - x1) * (y2 - y0));
    const Area = Math.abs(sA);
    const inv = (sA >= 0 ? 1 : -1) / (2 * Area);
    const b0 = (y1 - y2) * inv, b1 = (y2 - y0) * inv, b2 = (y0 - y1) * inv;
    const c0 = (x2 - x1) * inv, c1 = (x0 - x2) * inv, c2 = (x1 - x0) * inv;

    // Whitney curls are constant over the element.
    const cw0 = 2 * (b0*c1 - b1*c0), cw1 = 2 * (b1*c2 - b2*c1), cw2 = 2 * (b2*c0 - b0*c2);

    _Att.fill(0); _Wtt.fill(0); _Dzt.fill(0);

    for (let q = 0; q < NQ; q++) {
        const w = QW[q], l0 = QL1[q], l1 = QL2[q], l2 = QL3[q];

        // Whitney (lowest-order Nedelec) on edges (0,1), (1,2), (2,0):
        // W_pq = λ_p∇λ_q − λ_q∇λ_p.
        const W0x = l0*b1 - l1*b0, W0y = l0*c1 - l1*c0;
        const W1x = l1*b2 - l2*b1, W1y = l1*c2 - l2*c1;
        const W2x = l2*b0 - l0*b2, W2y = l2*c0 - l0*c2;
        _Wx[0] = W0x; _Wy[0] = W0y; _cW[0] = cw0;
        _Wx[1] = W1x; _Wy[1] = W1y; _cW[1] = cw1;
        _Wx[2] = W2x; _Wy[2] = W2y; _cW[2] = cw2;

        // Quadratic edge bases ne2 = W_pq·(λ_p − λ_q); curl via curl(fV) = f·curlV + ∇f×V.
        const d0 = l0 - l1, d1 = l1 - l2, d2 = l2 - l0;
        _Wx[4] = W0x*d0; _Wy[4] = W0y*d0;
        _Wx[5] = W1x*d1; _Wy[5] = W1y*d1;
        _Wx[6] = W2x*d2; _Wy[6] = W2y*d2;
        _cW[4] = cw0*d0 + (b0 - b1)*W0y - (c0 - c1)*W0x;
        _cW[5] = cw1*d1 + (b1 - b2)*W1y - (c1 - c2)*W1x;
        _cW[6] = cw2*d2 + (b2 - b0)*W2y - (c2 - c0)*W2x;

        // Face bubbles: nf1 = λ_0·W_12 − λ_1·W_20,  nf2 = λ_1·W_20 − λ_2·W_01.
        _Wx[3] = l0*W1x - l1*W2x; _Wy[3] = l0*W1y - l1*W2y;
        _cW[3] = l0*cw1 + (b0*W1y - c0*W1x) - l1*cw2 - (b1*W2y - c1*W2x);
        _Wx[7] = l1*W2x - l2*W0x; _Wy[7] = l1*W2y - l2*W0y;
        _cW[7] = l1*cw2 + (b1*W2y - c1*W2x) - l2*cw0 - (b2*W0y - c2*W0x);

        // Nodal P2 gradients: ∇(2λ²−λ) = (4λ−1)∇λ,  ∇(4λ_pλ_q) = 4(λ_q∇λ_p + λ_p∇λ_q).
        const f0 = 4*l0 - 1, f1 = 4*l1 - 1, f2 = 4*l2 - 1;
        _Gx[0] = b0*f0; _Gy[0] = c0*f0;
        _Gx[1] = b1*f1; _Gy[1] = c1*f1;
        _Gx[2] = b2*f2; _Gy[2] = c2*f2;
        _Gx[3] = 4*(b0*l1 + b1*l0); _Gy[3] = 4*(c0*l1 + c1*l0);
        _Gx[4] = 4*(b1*l2 + b2*l1); _Gy[4] = 4*(c1*l2 + c2*l1);
        _Gx[5] = 4*(b2*l0 + b0*l2); _Gy[5] = 4*(c2*l0 + c0*l2);

        for (let i = 0; i < 8; i++) {
            const ci = w * _cW[i], wxi = w * _Wx[i], wyi = w * _Wy[i], row = i * 8;
            for (let j = 0; j < 8; j++) {
                _Att[row + j] += ci * _cW[j];
                _Wtt[row + j] += wxi * _Wx[j] + wyi * _Wy[j];
            }
        }
        for (let i = 0; i < 6; i++) {
            const gxi = w * _Gx[i], gyi = w * _Gy[i], row = i * 8;
            for (let j = 0; j < 8; j++) _Dzt[row + j] += gxi * _Wx[j] + gyi * _Wy[j];
        }
    }

    for (let i = 0; i < 64; i++) {
        const wtt = _Wtt[i] * Area;
        _Att[i] *= Area;
        _Wtt[i] = wtt; _Dtt[i] = wtt;
        _BttRe[i] = epsRe * wtt; _BttIm[i] = epsIm * wtt;
    }
    // No face-DOF self-coupling in the B matrix.
    _Dtt[3*8+3] = _Dtt[3*8+7] = _Dtt[7*8+3] = _Dtt[7*8+7] = 0;
    for (let i = 0; i < 48; i++) _Dzt[i] *= Area;

    triP2Stiffness(nodes, v0, v1, v2, _Dzz1);
    for (let i = 0; i < 36; i++) {
        const m = Area * P2_MASS[i];
        _Dzz2Re[i] = epsRe * m; _Dzz2Im[i] = epsIm * m;
    }

    _p2out.Area = Area;
    return _p2out;
}

// P2 static element matrices (6x6 nodal only)
//
// Closed form instead of a per-triangle quadrature loop. On a straight-sided
// (affine) triangle the barycentric coordinate at quadrature point q is the
// quadrature's barycentric weight, so every P2 shape function and gradient
// there is a constant combination of the element's ∇λ:
//
//   ∇N_i = Σ_m D_im(q)·∇λ_m,  D constant per (i, m, q)   [each i touches ≤ 2 m]
//   ⇒ Sz_ij = Area·Σ_{m,n} (Σ_q w_q D_im D_jn)·(b_m b_n + c_m c_n)
//
// The bracketed tensor is a property of the reference element, built once below
// as a flat (dst, gIdx, coeff) term list, so the element stiffness is 81 fused
// multiply-adds on the six geometric products g_mn, with no quadrature loop, no
// basis-function calls and no per-triangle allocation. Same integral as the
// 6-point rule it replaces (which is exact for the degree-2 integrand), to
// rounding.
//
// The companion mass matrix M_ij = ∫N_i N_j and load vector F_i = ∫N_i are
// likewise Area x a constant reference matrix (P2_MASS / P2_LOAD below), since
// the 6-point rule is exact through degree 4. mqs_loss.js consumes both.

// Atom a = (basis i, gradient index m, λ-profile over the quadrature points).
// Vertex k: ∇N = (4λ_k − 1)∇λ_k. Edge k=(p,q): ∇N = 4λ_q∇λ_p + 4λ_p∇λ_q.
const P2_ATOM_BASIS = [0, 1, 2, 3, 3, 4, 4, 5, 5];
const P2_ATOM_GRAD  = [0, 1, 2, 0, 1, 1, 2, 2, 0];
const { P2_TERM_DST, P2_TERM_G, P2_TERM_C, P2_MASS, P2_LOAD } = (() => {
    const QL = [QL1, QL2, QL3];
    // λ-profile of each atom at every quadrature point.
    const prof = [
        q => 4 * QL1[q] - 1, q => 4 * QL2[q] - 1, q => 4 * QL3[q] - 1,
        q => 4 * QL2[q], q => 4 * QL1[q],        // edge 0 = (v0,v1)
        q => 4 * QL3[q], q => 4 * QL2[q],        // edge 1 = (v1,v2)
        q => 4 * QL1[q], q => 4 * QL3[q],        // edge 2 = (v2,v0)
    ];
    // Merge the 81 atom pairs by (destination entry, g index): duplicates differ
    // only in their coefficient, and summing them up front shrinks the per-triangle
    // loop.
    const merged = new Map();
    for (let a = 0; a < 9; a++) for (let b = 0; b < 9; b++) {
        let c = 0;
        for (let q = 0; q < NQ; q++) c += QW[q] * prof[a](q) * prof[b](q);
        if (c === 0) continue;
        const dst = P2_ATOM_BASIS[a] * 6 + P2_ATOM_BASIS[b];
        const g = P2_ATOM_GRAD[a] * 3 + P2_ATOM_GRAD[b];
        const key = dst * 9 + g;
        merged.set(key, (merged.get(key) || 0) + c);
    }
    const keys = [...merged.keys()].sort((x, y) => x - y);
    const dst = new Int32Array(keys.length), gi = new Int32Array(keys.length);
    const co = new Float64Array(keys.length);
    keys.forEach((k, i) => { dst[i] = (k / 9) | 0; gi[i] = k % 9; co[i] = merged.get(k); });
    // Reference mass matrix and load vector (per unit area).
    const N = (i, q) => (i < 3 ? (l => 2 * l * l - l)(QL[i][q])
                               : 4 * QL[[0, 1, 2][i - 3]][q] * QL[[1, 2, 0][i - 3]][q]);
    const M = new Float64Array(36), L = new Float64Array(6);
    for (let q = 0; q < NQ; q++) for (let i = 0; i < 6; i++) {
        L[i] += QW[q] * N(i, q);
        for (let j = 0; j < 6; j++) M[i * 6 + j] += QW[q] * N(i, q) * N(j, q);
    }
    return { P2_TERM_DST: dst, P2_TERM_G: gi, P2_TERM_C: co, P2_MASS: M, P2_LOAD: L };
})();
export { P2_MASS, P2_LOAD };

const _g9 = new Float64Array(9);

// Element P2 stiffness ∫∇N_i·∇N_j into `Sz` (36 entries, caller-owned so the
// per-triangle loops stay allocation-free). Returns the triangle area.
export function triP2Stiffness(nodes, v0, v1, v2, Sz) {
    const x0 = nodes[2*v0], y0 = nodes[2*v0+1];
    const x1 = nodes[2*v1], y1 = nodes[2*v1+1];
    const x2 = nodes[2*v2], y2 = nodes[2*v2+1];
    const sA = 0.5 * ((x0 - x2) * (y1 - y0) - (x0 - x1) * (y2 - y0));
    const Area = Math.abs(sA);
    // ∇λ_m = sign·(b_m, c_m)/(2A); the sign cancels in every g_mn product, so it
    // is dropped here (triCoefficients keeps it because λ itself needs it).
    const inv = 1 / (2 * Area);
    const b0 = (y1 - y2) * inv, b1 = (y2 - y0) * inv, b2 = (y0 - y1) * inv;
    const c0 = (x2 - x1) * inv, c1 = (x0 - x2) * inv, c2 = (x1 - x0) * inv;
    const g = _g9;
    g[0] = b0*b0 + c0*c0; g[4] = b1*b1 + c1*c1; g[8] = b2*b2 + c2*c2;
    g[1] = g[3] = b0*b1 + c0*c1;
    g[2] = g[6] = b0*b2 + c0*c2;
    g[5] = g[7] = b1*b2 + c1*c2;
    Sz.fill(0);
    for (let k = 0; k < P2_TERM_DST.length; k++)
        Sz[P2_TERM_DST[k]] += P2_TERM_C[k] * g[P2_TERM_G[k]];
    for (let i = 0; i < 36; i++) Sz[i] *= Area;
    return Area;
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
    // Containment test. Shapeless rects take the literal legacy bbox expression at the
    // original TOL, so every rectangular medium is bit-identical. Shaped conductors
    // (coax) test the polygon, at the mesher's geometric tolerance — a boundary node
    // must classify as metal or the PEC surface develops holes.
    const STOL = Math.max(TOL, condRect.geomTol || 0);
    const inCond = (cr, x, y) => (cr.shape
        ? shapeContains(cr, x, y, STOL)
        : (x >= cr.xmin - TOL && x <= cr.xmax + TOL && y >= cr.ymin - TOL && y <= cr.ymax + TOL));

    // Identify conductor nodes and edges
    const isCondNode = new Uint8Array(nNodes);
    const condNodeGroup = new Int8Array(nNodes);
    for (let n = 0; n < nNodes; n++) {
        const x = nodes[2 * n], y = nodes[2 * n + 1];
        for (let ci = 0; ci < condRects.length; ci++) {
            if (inCond(condRects[ci], x, y)) {
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
            if (inCond(condRects[ci], xm, ym)) {
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
        let interior = false;
        for (const cr of condRects) {
            if (inCond(cr, xc, yc)) { interior = true; break; }
        }
        if (interior) continue;
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
        const m = computeTriP2Matrices(nodes, v0, v1, v2, eps.re, eps.im);

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

    // Robin ABC entries: A += j·k₀√ε·(1/L or 1/(3L)) on boundary edge DOFs — a pure
    // per-unit-k₀ imaginary template. On edge (p,q) parameterized by s∈[0,1] the
    // transverse tangential basis is ne1_tang(s) = 1/L (constant), ne2_tang(s) =
    // (1/L)(1−2s), giving self-integrals ∫ne1²dl = 1/L, ∫ne2²dl = 1/(3L),
    // ∫ne1·ne2 dl = 0. (Strict === true: 'pmc' walls skip the Robin terms. The P2
    // NODAL Robin on Ez is deliberately absent — forcing ∂Ez/∂n + jk₀Ez = 0 distorts
    // quasi-TEM eigenvectors; the transverse ABC alone suppresses cavity modes.)
    //
    // The matched coefficient is the local plane-wave wavenumber k*sqrt(er), not vacuum k.
    // A wall edge in dielectric uses one triangle adjacent to the boundary to
    // sample the permittivity.
    const rR = [], rC = [], rV = [];
    if (abc.top === true || abc.left === true || abc.right === true || abc.bottom === true) {
        const ymax = condRect.ymax_domain, ymin = condRect.ymin_domain;
        const xmin = condRect.xmin_domain, xmax = condRect.xmax_domain;
        const TOL = 1e-12;
        // Sample er from adjacent triangle. Built only when the ABC is live.
        const edgeTri = new Int32Array(nEdges).fill(-1);
        for (let t = 0; t < nTris; t++) {
            for (let le = 0; le < 3; le++) {
                const eIdx = triEdges[3*t + le];
                if (edgeTri[eIdx] < 0) edgeTri[eIdx] = t;
            }
        }
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
            const tAdj = edgeTri[e];
            const erAdj = (tAdj >= 0 && epsMap && epsMap[tAdj]) ? epsMap[tAdj].re : 1;
            const nWall = erAdj > 0 ? Math.sqrt(erAdj) : 1;   // local k = k₀·n_wall
            const ef1 = edgeF[2*e];
            if (ef1 >= 0) { rR.push(ef1); rC.push(ef1); rV.push(nWall / L); }
            const ef2 = edgeF[2*e+1];
            if (ef2 >= 0) { rR.push(ef2); rC.push(ef2); rV.push(nWall / (3 * L)); }
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
    // One conversion per pattern (A and B), each carrying all of its aligned value
    // templates, not one conversion per template.
    const csrA = tripletsToCSRMulti(rowsA, colsA, N, [v0A, v1ReA, v1ImA, vrA]);
    const csrB = tripletsToCSRMulti(eBr.subarray(0, bNnz), eBc.subarray(0, bNnz), N,
        [eB0.subarray(0, bNnz), eB1Re.subarray(0, bNnz), eB1Im.subarray(0, bNnz)]);

    return {
        N,
        A: { rowPtr: csrA.rowPtr, colIdx: csrA.colIdx,
             v0Re: csrA.vals[0], v1Re: csrA.vals[1], v1Im: csrA.vals[2], vrIm: csrA.vals[3] },
        B: { rowPtr: csrB.rowPtr, colIdx: csrB.colIdx,
             v0Re: csrB.vals[0], v1Re: csrB.vals[1], v1Im: csrB.vals[2] },
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
export function solveTriStatic(mesh, fm, epsMap, condPotentials = null, directSolver = null) {
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

    // COO in preallocated typed arrays (≤ 36 entries per triangle), and one set of
    // per-element scratch buffers for the whole loop: at mesh sizes here the JS-array
    // pushes and the per-triangle allocations cost more than the arithmetic.
    const Rows = new Int32Array(36 * nTris), Cols = new Int32Array(36 * nTris);
    const Vals = new Float64Array(36 * nTris);
    let nCoo = 0;
    const rhs = new Float64Array(nFreeDof);
    const nLocal = 6;
    const Sz = new Float64Array(36);
    const globalDof = new Int32Array(nLocal);
    const isDirichlet = new Uint8Array(nLocal);
    const dirichletVal = new Float64Array(nLocal);
    const verts = new Int32Array(3);

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        triP2Stiffness(nodes, v0, v1, v2, Sz);
        const eps = epsMap ? epsMap[t].re : 1.0;
        verts[0] = v0; verts[1] = v1; verts[2] = v2;
        isDirichlet.fill(0);

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
            const gi = globalDof[li]; if (gi < 0) continue;
            for (let lj = 0; lj < nLocal; lj++) {
                const v = eps * Sz[li * nLocal + lj];
                if (v === 0) continue;
                const gj = globalDof[lj];
                if (gj >= 0) {
                    Rows[nCoo] = gi; Cols[nCoo] = gj; Vals[nCoo] = v; nCoo++;
                } else if (isDirichlet[lj]) {
                    rhs[gi] -= v * dirichletVal[lj];
                }
            }
        }
    }


    const csrS = tripletsToCSR(Rows.subarray(0, nCoo), Cols.subarray(0, nCoo),
                               Vals.subarray(0, nCoo), nFreeDof);
    // Direct WASM solve (SPD → LDLT, pattern-cached across the eps/air pair on the
    // same freedom map) when a solver is supplied; JS CG otherwise / on failure.
    let phiFree = null;
    if (directSolver && nFreeDof > 0) {
        try {
            const sol = directSolver(nFreeDof, csrS, [rhs])[0];
            if (sol && sol.every(Number.isFinite)) phiFree = sol;
        } catch { phiFree = null; }
    }
    if (!phiFree) {
        const { x, iters, residual, converged } = solveCG(csrS, rhs, nFreeDof);
        if (!converged) {
            console.warn(`Static CG did not converge in ${iters} iterations ` +
                `(residual ${residual.toExponential(2)}) — static C/Z0 may be inaccurate.`);
        }
        phiFree = x;
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
    const { nodes, tris, triEdges, nTris } = mesh;
    const { phiVertex, phiEdge } = phi;
    let W = 0;
    const Sz = new Float64Array(36);
    const localPhi = new Float64Array(6);

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        triP2Stiffness(nodes, v0, v1, v2, Sz);
        const eps = epsMap ? epsMap[t].re : 1.0;

        localPhi[0] = phiVertex[v0]; localPhi[1] = phiVertex[v1]; localPhi[2] = phiVertex[v2];
        localPhi[3] = phiEdge[triEdges[3*t]];
        localPhi[4] = phiEdge[triEdges[3*t+1]];
        localPhi[5] = phiEdge[triEdges[3*t+2]];

        let e = 0;
        for (let li = 0; li < 6; li++) {
            const pi = localPhi[li], row = li * 6;
            for (let lj = 0; lj < 6; lj++) e += Sz[row + lj] * pi * localPhi[lj];
        }
        W += eps * e;
    }

    return 0.5 * W;
}

// Analytic field, Arnoldi start vector
// Project a known transverse E field onto the eigen DOF vector, for geometries that have
// no static solve to seed from. A hollow waveguide is the case: with no second conductor
// drivePotentials() is empty, so there is no potential for staticToEdgeDofs to project.
//
// Only the lowest-order (Whitney / ne1) transverse DOFs are filled, that DOF is the
// tangential line integral over the edge, taken here with the midpoint rule. The
// quadratic ne2 DOFs and the interior face bubbles are left at 0: they are a second-order
// correction to what is only a start vector, and the shift-invert target does the real
// work of selecting the mode. The whole longitudinal (Ez) block is left at 0 too, which
// is exact for a TE mode.
//
// Efun(x, y) -> [Ex, Ey].
export function analyticSeedDofs(mesh, fm, Efun) {
    const { nodes, edges, nEdges } = mesh;
    const { edgeF, nFreeTransverse, nFreeLongitudinal } = fm;
    const initVec = new Float64Array(nFreeTransverse + nFreeLongitudinal);
    for (let e = 0; e < nEdges; e++) {
        const ef = edgeF[2 * e];
        if (ef < 0) continue;
        const n0 = edges[2 * e], n1 = edges[2 * e + 1];
        const x0 = nodes[2 * n0], y0 = nodes[2 * n0 + 1];
        const x1 = nodes[2 * n1], y1 = nodes[2 * n1 + 1];
        const E = Efun(0.5 * (x0 + x1), 0.5 * (y0 + y1));
        initVec[ef] = E[0] * (x1 - x0) + E[1] * (y1 - y0);
    }
    return initVec;
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

    // --- Constraint SEGMENTS (shaped conductors: polygon sides) ---
    // A shaped conductor's boundary is not axis-aligned, so the onXLine/onYLine
    // machinery above cannot see it — and the bbox-wall test cannot either, since a
    // circle only touches its bounding box at four points. Without this, every
    // polygon-boundary node classifies as FREE and the Laplacian smoother pulls it
    // toward the mesh interior on every pass: the conductor shrinks until the
    // containment test stops matching it and the PEC surface develops holes.
    //
    // Such nodes are PINNED rather than allowed to slide. Bisection midpoints of a
    // straight segment land exactly on that segment, so pinning costs nothing
    // geometrically, whereas sliding would need a parametric version of the
    // clampLo/clampHi neighbour bounds and tends to bunch nodes near polygon vertices —
    // degrading exactly the surface resolution the conductor-loss line integral needs.
    //
    // Node classification is done ONCE here and reused by both the smoother's eligibility
    // test and the edge-swap guard, which would otherwise repeat it every pass. The
    // smoother still queries segmentAt for the candidate POSITIONS it proposes (those are
    // not nodes, so there is nothing to precompute for them), hence the union-bbox reject
    // in segmentAt, which is what keeps that per-pass query off the O(segments) path for
    // the majority of vertices, the ones nowhere near a shaped boundary.
    const segs = mesh.constraintSegments || [];
    // Per-segment geometry, precomputed once: padded bbox (the cheap reject), the edge
    // vector, and len·tol so the perpendicular-distance test needs no division. Segments
    // of zero length are dropped here rather than re-checked on every query.
    const segTol = Math.max(xMax - xMin, yMax - yMin) * REL_SHAPE_TOL;
    const segData = [];
    let hullX0 = Infinity, hullX1 = -Infinity, hullY0 = Infinity, hullY1 = -Infinity;
    for (const s of segs) {
        const ex = s.x1 - s.x0, ey = s.y1 - s.y0;
        const len = Math.hypot(ex, ey);
        if (!(len > 0)) continue;
        const d = {
            x0: s.x0, y0: s.y0, ex, ey, crossTol: len * segTol,
            xMin: Math.min(s.x0, s.x1) - segTol, xMax: Math.max(s.x0, s.x1) + segTol,
            yMin: Math.min(s.y0, s.y1) - segTol, yMax: Math.max(s.y0, s.y1) + segTol,
        };
        segData.push(d);
        hullX0 = Math.min(hullX0, d.xMin); hullX1 = Math.max(hullX1, d.xMax);
        hullY0 = Math.min(hullY0, d.yMin); hullY1 = Math.max(hullY1, d.yMax);
    }
    // Index of a constraint segment (px, py) lies on, or -1. `skip` excludes one already
    // found, so a vertex shared by two polygon sides reports both. The union-bbox test
    // first: most queries come from nodes nowhere near a shaped boundary and reject in O(1).
    const segmentAt = (px, py, skip = -1) => {
        if (px < hullX0 || px > hullX1 || py < hullY0 || py > hullY1) return -1;
        for (let si = 0; si < segData.length; si++) {
            if (si === skip) continue;
            const d = segData[si];
            if (px < d.xMin || px > d.xMax || py < d.yMin || py > d.yMax) continue;
            // Perpendicular distance to the (infinite) line; the bbox test above already
            // restricted us to the segment's extent.
            if (Math.abs((px - d.x0) * d.ey - (py - d.y0) * d.ex) > d.crossTol) continue;
            return si;
        }
        return -1;
    };
    // Which segment(s) each EXISTING node sits on. Computed once and reused by the
    // smoother's eligibility test and the edge-swap guard; only the smoother's candidate
    // positions (which are not nodes yet) have to go back through segmentAt.
    const nodeSegA = new Int32Array(newNodeCount).fill(-1);
    const nodeSegB = new Int32Array(newNodeCount).fill(-1);
    if (segData.length) {
        for (let n = 0; n < newNodeCount; n++) {
            const px = newNodes[2*n], py = newNodes[2*n+1];
            const a = segmentAt(px, py);
            if (a < 0) continue;
            nodeSegA[n] = a;
            nodeSegB[n] = segmentAt(px, py, a);
        }
    }
    const onSegment = (n) => nodeSegA[n] >= 0;
    const pointOnAnySegment = (px, py) => segmentAt(px, py) >= 0;

    // Classify ALL vertices: 0=pinned, 1=free, 2=along-Y, 3=along-X
    const canMove = new Uint8Array(newNodeCount);
    for (let n = 0; n < newNodeCount; n++) {
        const mx = newNodes[2*n], my = newNodes[2*n+1];
        if (onSegment(n)) { canMove[n] = 0; continue; }
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
        // Only triangles with an eligible, movable vertex contribute to any average
        // that gets used. The rest of the mesh is skipped on every pass. (The
        // neighbour sums are read only at eligible vertices.)
        const touch = new Uint8Array(nNewTris);
        let anyTouch = false;
        for (let t = 0; t < nNewTris; t++) {
            const va = newTriArr[3*t], vb = newTriArr[3*t+1], vc = newTriArr[3*t+2];
            if ((eligible[va] && canMove[va]) || (eligible[vb] && canMove[vb]) ||
                (eligible[vc] && canMove[vc])) { touch[t] = 1; anyTouch = true; }
        }
        if (!anyTouch) return;
        const sumX = new Float64Array(newNodeCount), sumY = new Float64Array(newNodeCount);
        const cnt = new Int32Array(newNodeCount);
        for (let smoothPass = 0; smoothPass < maxPasses; smoothPass++) {
            sumX.fill(0); sumY.fill(0); cnt.fill(0);
            for (let t = 0; t < nNewTris; t++) {
                if (!touch[t]) continue;
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
                    // Same rule for shaped boundaries: a free vertex that lands on a
                    // polygon side would be misread as a surface node next pass.
                    if (!snap && segData.length) snap = pointOnAnySegment(nx, ny);
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
        // Both endpoints on the SAME polygon side ⇒ this edge lies on a shaped
        // conductor's boundary. Swapping it would tear the surface apart. The domain
        // outline is already safe (one adjacent triangle, so tris2.length !== 2), but
        // an inner conductor's boundary has two — its interior is meshed.
        const a0 = nodeSegA[na], a1 = nodeSegB[na];
        if (a0 >= 0 && (a0 === nodeSegA[nb] || a0 === nodeSegB[nb])) return true;
        if (a1 >= 0 && (a1 === nodeSegA[nb] || a1 === nodeSegB[nb])) return true;
        if (onYLine(ax, ay) && onYLine(bx, by) && Math.abs(ay - by) < 1e-10) return true;
        if (onXLine(ax, ay) && onXLine(bx, by) && Math.abs(ax - bx) < 1e-10) return true;
        if (Math.abs(ax - bx) < 1e-10 && (Math.abs(ax - xMin) < 1e-10 || Math.abs(ax - xMax) < 1e-10)) return true;
        if (Math.abs(ay - by) < 1e-10 && (Math.abs(ay - yMin) < 1e-10 || Math.abs(ay - yMax) < 1e-10)) return true;
        return false;
    }

    // Triangles worth swapping around: only an edge with a sliver (q >= Q_SWAP)
    // on at least one side can be improved by a flip. Computed first so a mesh
    // with no slivers (common case after a band refinement) never pays for the
    // edge map below (3 map ops per triangle, up to 20 passes per refinement).
    const Q_SWAP = 10;
    function anySliverTri() {
        for (let t = 0; t < nNewTris; t++) if (triQuality(t) >= Q_SWAP) return true;
        return false;
    }

    function edgeSwapPass() {
        if (!anySliverTri()) return 0;
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
            if (qBefore < Q_SWAP) continue; // only target actual slivers
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
        // Quality-targeted smoothing of vertices near slivers. Determined
        // before the vtAdj rebuild so a mesh with no slivers left leaves the
        // loop without paying for it (vtAdj is not read after this loop).
        const Q_THRESH = 10;
        const sliverAdj = new Uint8Array(newNodeCount);
        let anySliver = false;
        for (let t = 0; t < nNewTris; t++) {
            if (triQuality(t) > Q_THRESH) {
                for (let k = 0; k < 3; k++) {
                    const v = newTriArr[3*t+k];
                    if (canMove[v] >= 1) { sliverAdj[v] = 1; anySliver = true; }
                }
            }
        }
        if (!anySliver) break;
        // Rebuild vtAdj after swaps
        for (let n = 0; n < newNodeCount; n++) vtAdj[n] = [];
        for (let t = 0; t < nNewTris; t++) {
            vtAdj[newTriArr[3*t]].push(t);
            vtAdj[newTriArr[3*t+1]].push(t);
            vtAdj[newTriArr[3*t+2]].push(t);
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
        constraintXs: mesh.constraintXs || [],
        // Carried forward so every refinement pass re-pins the shaped boundaries.
        constraintSegments: mesh.constraintSegments || [],
        // Bisection never adds elements where there were none, so the refined mesh has a
        // conductor interior exactly when its parent did. Carried forward so the MQS
        // applicability guard still sees it after any number of passes.
        condInteriorMeshed: mesh.condInteriorMeshed
    };
}
