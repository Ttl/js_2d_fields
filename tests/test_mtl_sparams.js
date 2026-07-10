// test_mtl_sparams.js — computeSParamsDifferentialMTL against an independent ground truth.
//
// The MTL 4-port conversion works with 2×2 matrix blocks that do NOT commute for an
// asymmetric pair, so operand order matters everywhere (B = γ⁻¹sinh(γℓ)·[Z], not
// [Z]·γ⁻¹sinh(γℓ); S22 = den⁻¹·(−A+BY−CZ+D), not (…)·den⁻¹). This test checks the full
// 4-port against a ground truth built with none of those formulas:
//
//   1. Chain matrix by LUMPED CASCADE: N = 2^16 symmetric L-sections (series Z·dz/2,
//      shunt Y·dz, series Z·dz/2), chained by plain block multiplication via repeated
//      squaring. Discretization error ~ (γℓ)³/N² ≈ 1e-6 here.
//   2. ABCD → S by solving the port-wave BOUNDARY CONDITIONS as a complex 4×4 linear
//      system per excitation (Gaussian elimination) — no closed-form S expressions.
//
// Properties checked (asymmetric RLGC, several frequencies):
//   1. GROUND TRUTH — full 4×4 S matches the cascade ground truth.
//   2. RECIPROCITY  — S = Sᵀ (the pre-fix code violated this with |S34−S43| ≈ 0.2).
//   3. FLIP SYMMETRY — a uniform line looks the same from both ends: near-near block
//      equals far-far block.
//   4. SYMMETRIC LIMIT — for a symmetric pair the MTL result reduces exactly to the
//      odd/even combination (computeSParamsDifferential), with SDC/SCD = 0.

import { Complex } from '../src/complex.js';
import { Matrix2x2 } from '../src/matrix.js';
import { computeSParamsDifferentialMTL, computeSParamsDifferential } from '../src/sparameters.js';

// --- Asymmetric physical 2×2 p.u.l. matrices (SI), symmetric as matrices ----
const R2 = [[8.0, 1.5], [1.5, 12.0]];             // Ω/m
const L2 = [[3.2e-7, 6.0e-8], [6.0e-8, 2.4e-7]];  // H/m
const G2 = [[2.0e-4, -4.0e-5], [-4.0e-5, 3.5e-4]]; // S/m
const C2 = [[1.3e-10, -2.5e-11], [-2.5e-11, 9.5e-11]]; // F/m
// Homogeneous-dielectric variant: L = με·C⁻¹ makes [L][C] = με·I, i.e. exactly
// velocity-DEGENERATE modes with still-asymmetric matrices (a stripline-like pair).
// This regime is a 0/0 trap for the γ=√(ZY) matrix square root.
const MU_EPS = 4e-7 * Math.PI * 8.854e-12 * 4.0;  // εr = 4
const detC = C2[0][0] * C2[1][1] - C2[0][1] * C2[1][0];
const L2deg = [[MU_EPS * C2[1][1] / detC, -MU_EPS * C2[0][1] / detC],
               [-MU_EPS * C2[1][0] / detC, MU_EPS * C2[0][0] / detC]];
const LENGTH = 0.05, Z_REF = 50;
const FREQS = [1e6, 1e9, 1e10];
const N_LOG2 = 16;                 // cascade sections = 2^16

const TOL_GT = 1e-4;               // vs. cascade ground truth (cascade error ~1e-6)
const TOL_EXACT = 1e-8;            // analytic identities (reciprocity, flip, symmetric limit)

// --- Block-2×2 chain matrices (elements are Matrix2x2) ----------------------
const I2M = new Matrix2x2(new Complex(1), new Complex(0), new Complex(0), new Complex(1));
const Z2M = new Matrix2x2(new Complex(0), new Complex(0), new Complex(0), new Complex(0));
const mC = (re, im) => new Matrix2x2(
    new Complex(re[0][0], im[0][0]), new Complex(re[0][1], im[0][1]),
    new Complex(re[1][0], im[1][0]), new Complex(re[1][1], im[1][1]));

const blockMul = (P, Q) => ({
    A: P.A.mul(Q.A).add(P.B.mul(Q.C)), B: P.A.mul(Q.B).add(P.B.mul(Q.D)),
    C: P.C.mul(Q.A).add(P.D.mul(Q.C)), D: P.C.mul(Q.B).add(P.D.mul(Q.D)),
});

// Cascade chain matrix for the whole line: one symmetric section, squared N_LOG2 times.
function cascadeABCD(freq, Lmat) {
    const omega = 2 * Math.PI * freq;
    const Z = mC(R2, Lmat.map(r => r.map(v => v * omega)));
    const Y = mC(G2, C2.map(r => r.map(v => v * omega)));
    const dz = LENGTH / Math.pow(2, N_LOG2);
    const serHalf = { A: I2M, B: Z.mul(dz / 2), C: Z2M, D: I2M };   // series Z·dz/2
    const shunt = { A: I2M, B: Z2M, C: Y.mul(dz), D: I2M };         // shunt Y·dz
    let sec = blockMul(blockMul(serHalf, shunt), serHalf);
    for (let k = 0; k < N_LOG2; k++) sec = blockMul(sec, sec);
    return sec;
}

// --- Complex 4×4 Gaussian elimination with partial pivoting -----------------
function solve4(M, b) {
    const A = M.map((row, i) => [...row, b[i]]);
    for (let col = 0; col < 4; col++) {
        let piv = col;
        for (let r = col + 1; r < 4; r++) if (A[r][col].abs() > A[piv][col].abs()) piv = r;
        [A[col], A[piv]] = [A[piv], A[col]];
        for (let r = col + 1; r < 4; r++) {
            const f = A[r][col].div(A[col][col]);
            for (let c = col; c < 5; c++) A[r][c] = A[r][c].sub(f.mul(A[col][c]));
        }
    }
    const x = new Array(4);
    for (let r = 3; r >= 0; r--) {
        let s = A[r][4];
        for (let c = r + 1; c < 4; c++) s = s.sub(A[r][c].mul(x[c]));
        x[r] = s.div(A[r][r]);
    }
    return x;
}

// ABCD → S by boundary conditions. Chain convention: [V1; I1] = Φ·[V2; I2] with I2 the
// current flowing OUT of the far port. For incident waves a (ports: 1,2 = near, 3,4 = far):
//   near rows: (Φ.A + Z0·Φ.C)·V2 + (Φ.B + Z0·Φ.D)·I2 = 2√Z0·a_near
//   far rows:  V2 − Z0·I2 = 2√Z0·a_far
// then b_near = ((Φ.A − Z0Φ.C)V2 + (Φ.B − Z0Φ.D)I2)/(2√Z0), b_far = (V2 + Z0·I2)/(2√Z0).
function abcdToS(Phi, Z0) {
    const e = (m, i, j) => [[m.a, m.b], [m.c, m.d]][i][j];
    const S = Array.from({ length: 4 }, () => new Array(4));
    const P = Phi.A.add(Phi.C.mul(Z0)), Q = Phi.B.add(Phi.D.mul(Z0));   // near-row blocks
    const Pm = Phi.A.sub(Phi.C.mul(Z0)), Qm = Phi.B.sub(Phi.D.mul(Z0)); // reflected-wave blocks
    const s = 2 * Math.sqrt(Z0);
    for (let j = 0; j < 4; j++) {
        const a = [0, 0, 0, 0]; a[j] = 1;
        const M = [
            [e(P, 0, 0), e(P, 0, 1), e(Q, 0, 0), e(Q, 0, 1)],
            [e(P, 1, 0), e(P, 1, 1), e(Q, 1, 0), e(Q, 1, 1)],
            [new Complex(1), new Complex(0), new Complex(-Z0), new Complex(0)],
            [new Complex(0), new Complex(1), new Complex(0), new Complex(-Z0)],
        ];
        const rhs = a.map(v => new Complex(s * v));
        const [v1, v2, i1, i2] = solve4(M, rhs);
        const V2 = [v1, v2], I2 = [i1, i2];
        for (let r = 0; r < 2; r++) {
            let bn = new Complex(0), bf = V2[r].add(I2[r].mul(Z0));
            for (let c = 0; c < 2; c++)
                bn = bn.add(e(Pm, r, c).mul(V2[c])).add(e(Qm, r, c).mul(I2[c]));
            S[r][j] = bn.mul(1 / s);
            S[r + 2][j] = bf.mul(1 / s);
        }
    }
    return S;
}

const cdiff = (a, b) => a.sub(b).abs();
const maxDiff = (S, T) => {
    let m = 0;
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) m = Math.max(m, cdiff(S[i][j], T[i][j]));
    return m;
};

// --- Test -------------------------------------------------------------------
let failed = 0;
const check = (name, cond, detail = '') => {
    console.log(`  ${cond ? '✓' : '✗'} ${name}${detail && !cond ? ` — ${detail}` : ''}`);
    if (!cond) failed++;
};

console.log('\n### MTL 4-port S-parameters vs. lumped-cascade ground truth ###');

for (const [Lmat, tag] of [[L2, ''], [L2deg, ', velocity-degenerate']]) for (const freq of FREQS) {
    const label = (freq >= 1e9 ? `${freq / 1e9} GHz` : `${freq / 1e6} MHz`) + tag;
    const sp = computeSParamsDifferentialMTL(freq, R2, Lmat, G2, C2, LENGTH, Z_REF);
    const Sgt = abcdToS(cascadeABCD(freq, Lmat), Z_REF);

    // (1) Ground truth
    const dGT = maxDiff(sp.S, Sgt);
    check(`${label}: matches cascade ground truth (|ΔS| < ${TOL_GT})`, dGT < TOL_GT, `max |ΔS| ${dGT.toExponential(2)}`);

    // (2) Reciprocity: S = Sᵀ
    let dRec = 0;
    for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) dRec = Math.max(dRec, cdiff(sp.S[i][j], sp.S[j][i]));
    check(`${label}: reciprocity S = Sᵀ`, dRec < TOL_EXACT, `max |S_ij − S_ji| ${dRec.toExponential(2)}`);

    // (3) Longitudinal flip symmetry: near-near block == far-far block
    let dFlip = 0;
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) dFlip = Math.max(dFlip, cdiff(sp.S[i][j], sp.S[i + 2][j + 2]));
    check(`${label}: uniform-line flip symmetry (S11 blk = S22 blk)`, dFlip < TOL_EXACT, `max |ΔS| ${dFlip.toExponential(2)}`);
}

// (4) Symmetric limit: equal diagonals → must reduce to the odd/even combination.
console.log('\n### Symmetric-pair limit reduces to odd/even combination ###');
const Rs = [[10, 2], [2, 10]], Ls = [[3.0e-7, 6.0e-8], [6.0e-8, 3.0e-7]];
const Gs = [[2.5e-4, -5.0e-5], [-5.0e-5, 2.5e-4]], Cs = [[1.2e-10, -2.4e-11], [-2.4e-11, 1.2e-10]];
for (const freq of FREQS) {
    const label = freq >= 1e9 ? `${freq / 1e9} GHz` : `${freq / 1e6} MHz`;
    const sp = computeSParamsDifferentialMTL(freq, Rs, Ls, Gs, Cs, LENGTH, Z_REF);
    // Odd mode: opposite drive → X_odd = X11 − X12; even: X_even = X11 + X12.
    const rlgc = sign => ({ R: Rs[0][0] + sign * Rs[0][1], L: Ls[0][0] + sign * Ls[0][1],
                            G: Gs[0][0] + sign * Gs[0][1], C: Cs[0][0] + sign * Cs[0][1] });
    const ref = computeSParamsDifferential(freq, rlgc(-1), rlgc(+1), LENGTH, Z_REF);
    const dMM = Math.max(
        cdiff(sp.SDD11, ref.SDD11), cdiff(sp.SDD21, ref.SDD21),
        cdiff(sp.SCC11, ref.SCC11), cdiff(sp.SCC21, ref.SCC21));
    const dConv = Math.max(sp.SDC11.abs(), sp.SDC21.abs(), sp.SCD11.abs(), sp.SCD21.abs());
    check(`${label}: SDD/SCC match odd/even combination`, dMM < TOL_EXACT, `max |Δ| ${dMM.toExponential(2)}`);
    check(`${label}: SDC/SCD vanish for symmetric pair`, dConv < TOL_EXACT, `max |S| ${dConv.toExponential(2)}`);
}

console.log(failed ? `\nSUMMARY: ${failed} check(s) failed` : '\nSUMMARY: all checks passed');
if (failed) process.exitCode = 1;
