// sparameters.js - S-Parameter Calculations for Transmission Lines
import { Complex } from './complex.js';
import { Matrix2x2 } from './matrix.js';

/**
 * Compute 2-port S-parameters for a single-ended transmission line
 *
 * @param {number} freq - Frequency in Hz
 * @param {object} rlgc - RLGC parameters {R, L, G, C} in SI units
 * @param {number} length - Line length in meters
 * @param {number} Z_ref - Reference impedance in Ohms (typically 50)
 * @returns {object} - {S11, S21, S12, S22} as Complex numbers
 */
function computeSParamsSingleEnded(freq, rlgc, length, Z_ref) {
    const omega = 2 * Math.PI * freq;
    const { R, L, G, C } = rlgc;

    // Handle DC case (frequency = 0)
    // At DC, the transmission line behaves as a simple series resistance R*length
    if (freq === 0 || omega === 0) {
        const R_total = R * length;
        const Zr = Z_ref;

        // ABCD matrix for series resistance:
        // A = 1, B = R_total, C = 0, D = 1
        const A = new Complex(1, 0);
        const B = new Complex(R_total, 0);
        const C_abcd = new Complex(0, 0);
        const D = new Complex(1, 0);

        // Convert to S-parameters
        // den = A + B/Zr + C*Zr + D = 1 + R_total/Zr + 0 + 1 = 2 + R_total/Zr
        const den = new Complex(2 + R_total / Zr, 0);

        // S11 = (A + B/Zr - C*Zr - D) / den = (1 + R_total/Zr - 0 - 1) / den = (R_total/Zr) / den
        const S11 = new Complex(R_total / Zr, 0).div(den);

        // S21 = 2 / den
        const S21 = new Complex(2, 0).div(den);

        const S12 = S21;  // Reciprocal

        // S22 = (-A + B/Zr - C*Zr + D) / den = (-1 + R_total/Zr - 0 + 1) / den = (R_total/Zr) / den
        const S22 = new Complex(R_total / Zr, 0).div(den);

        return { S11, S21, S12, S22 };
    }

    // Series impedance per unit length: Z = R + jwL
    const Z_per_length = new Complex(R, omega * L);

    // Shunt admittance per unit length: Y = G + jwC
    const Y_per_length = new Complex(G, omega * C);

    // Propagation constant: gamma = sqrt(Z * Y)
    const gamma = Z_per_length.mul(Y_per_length).sqrt();

    // Characteristic impedance: Z0 = sqrt(Z / Y)
    const Z0 = Z_per_length.div(Y_per_length).sqrt();

    // gamma * length
    const gl = gamma.mul(length);

    // ABCD matrix elements
    // A = cosh(gamma * l)
    // B = Z0 * sinh(gamma * l)
    // C = sinh(gamma * l) / Z0
    // D = cosh(gamma * l)
    const A = gl.cosh();
    const B = Z0.mul(gl.sinh());
    const C_abcd = gl.sinh().div(Z0);
    const D = gl.cosh();

    // Convert ABCD to S-parameters
    // Reference impedance
    const Zr = new Complex(Z_ref, 0);

    // Common denominator: A + B/Zr + C*Zr + D
    const den = A.add(B.div(Zr)).add(C_abcd.mul(Zr)).add(D);

    // S11 = (A + B/Zr - C*Zr - D) / den
    const S11 = A.add(B.div(Zr)).sub(C_abcd.mul(Zr)).sub(D).div(den);

    // S21 = 2 / den (for reciprocal network)
    const S21 = new Complex(2, 0).div(den);

    // S12 = 2 * (AD - BC) / den = 2 / den (for reciprocal lossless or lossy uniform line, AD-BC=1)
    // For a uniform transmission line, det(ABCD) = 1, so S12 = S21
    const S12 = S21;

    // S22 = (-A + B/Zr - C*Zr + D) / den
    const S22 = A.neg().add(B.div(Zr)).sub(C_abcd.mul(Zr)).add(D).div(den);

    return { S11, S21, S12, S22 };
}

/**
 * Compute modal characteristic impedance from RLGC parameters
 * @param {number} freq - Frequency in Hz
 * @param {object} rlgc - RLGC parameters {R, L, G, C}
 * @returns {number} - Modal Z0 magnitude in Ohms
 */
function computeZ0(freq, rlgc) {
    const omega = 2 * Math.PI * freq;
    const { R, L, G, C } = rlgc;

    // Handle DC case (frequency = 0)
    // At DC, Z0 = sqrt(R/G) with G=0 gives infinity
    if (freq === 0 || omega === 0) {
        return 1e12;  // Return very large impedance for DC
    }

    // Series impedance per unit length: Z = R + jwL
    const Z_per_length = new Complex(R, omega * L);

    // Shunt admittance per unit length: Y = G + jwC
    const Y_per_length = new Complex(G, omega * C);

    // Characteristic impedance: Z0 = sqrt(Z / Y)
    const Z0 = Z_per_length.div(Y_per_length).sqrt();

    return Z0.abs();
}

/**
 * Compute 4-port S-parameters for a differential transmission line
 * Ports 1,3 are at one end (positive and negative), Ports 2,4 at the other end
 *
 * This function uses modal decomposition with odd and even mode parameters.
 *
 * @param {number} freq - Frequency in Hz
 * @param {object} rlgc_odd - RLGC parameters for odd mode {R, L, G, C}
 * @param {object} rlgc_even - RLGC parameters for even mode {R, L, G, C}
 * @param {number} length - Line length in meters
 * @param {number} Z_ref - Reference impedance in Ohms (typically 50)
 * @returns {object} - {S: 4x4 array of Complex, SDD11, SDD21, SCC11, SCC21, SDC11, SDC21, SCD11, SCD21}
 */
function computeSParamsDifferential(freq, rlgc_odd, rlgc_even, length, Z_ref) {
    // Compute single-ended S-params for each mode using system reference impedance
    const S_odd = computeSParamsSingleEnded(freq, rlgc_odd, length, Z_ref);
    const S_even = computeSParamsSingleEnded(freq, rlgc_even, length, Z_ref);

    // For ideal symmetric differential pairs with no coupling between modes,
    // the 4-port S-matrix can be constructed from odd and even mode responses.
    //
    // Port assignment:
    //   Port 1 = near end, trace + (in+)
    //   Port 2 = near end, trace - (in-)
    //   Port 3 = far end, trace + (out+)
    //   Port 4 = far end, trace - (out-)
    //
    // Modal decomposition for symmetric coupled lines:
    //   When exciting port 1 (a1=1, others=0):
    //   - Odd mode excitation: a_odd = (a1 - a2)/√2 = 1/√2
    //   - Even mode excitation: a_even = (a1 + a2)/√2 = 1/√2
    //
    // S-parameter formulas (for symmetric coupled lines):
    //   S11 = S22 = (S_odd_11 + S_even_11) / 2  (reflection at near end)
    //   S33 = S44 = (S_odd_22 + S_even_22) / 2  (reflection at far end)
    //   S21 = S12 = (S_even_11 - S_odd_11) / 2  (near-end coupling, NEXT)
    //   S43 = S34 = (S_even_22 - S_odd_22) / 2  (far-end reflection coupling)
    //   S31 = S13 = S42 = S24 = (S_odd_21 + S_even_21) / 2  (through transmission)
    //   S41 = S14 = S32 = S23 = (S_even_21 - S_odd_21) / 2  (far-end coupling, FEXT)

    const half = new Complex(0.5, 0);

    // Calculate 4-port S-parameters
    const Snn = S_odd.S11.add(S_even.S11).mul(half);   // Near-end reflection (S11, S22)
    const Sff = S_odd.S22.add(S_even.S22).mul(half);   // Far-end reflection (S33, S44)
    const Snext = S_even.S11.sub(S_odd.S11).mul(half); // Near-end coupling (S21, S12)
    const Sfref = S_even.S22.sub(S_odd.S22).mul(half); // Far-end reflection coupling (S43, S34)
    const Sthru = S_odd.S21.add(S_even.S21).mul(half); // Through transmission (S31, S13, S42, S24)
    const Sfext = S_even.S21.sub(S_odd.S21).mul(half); // Far-end coupling (S41, S14, S32, S23)

    // Build 4x4 matrix (symmetric coupled line)
    const S = [
        [Snn,   Snext, Sthru, Sfext],  // Row 1: S11, S12, S13, S14
        [Snext, Snn,   Sfext, Sthru],  // Row 2: S21, S22, S23, S24
        [Sthru, Sfext, Sff,   Sfref],  // Row 3: S31, S32, S33, S34
        [Sfext, Sthru, Sfref, Sff  ]   // Row 4: S41, S42, S43, S44
    ];

    // Mixed-mode S-parameters (for plotting)
    // SDD11 = S_odd_11 (differential reflection)
    // SDD21 = S_odd_21 (differential transmission)
    // SCC11 = S_even_11 (common-mode reflection)
    // SCC21 = S_even_21 (common-mode transmission)
    // SDC, SCD = 0 for symmetric line

    return {
        S,
        SDD11: S_odd.S11,
        SDD21: S_odd.S21,
        SCC11: S_even.S11,
        SCC21: S_even.S21,
        SDC11: new Complex(0, 0),
        SDC21: new Complex(0, 0),
        SCD11: new Complex(0, 0),
        SCD21: new Complex(0, 0)
    };
}

/**
 * Differential S-parameters, automatically using the full MTL 4-port when the physical 2×2
 * [C]/[L] matrices are available (asymmetric-correct) and falling back to the odd/even
 * combination otherwise. Loss [R]/[G] use the symmetric modal reconstruction (the reactive
 * asymmetry that drives S11≠S22 / mode conversion is exact; loss asymmetry is neglected).
 * @param {object} physMatrix {C:[[..]], L:[[..]]} true physical p.u.l. matrices, or null
 */
function computeSParamsDiffAuto(freq, rlgc_odd, rlgc_even, physMatrix, length, Z_ref) {
    if (!physMatrix || !physMatrix.C || !physMatrix.L) {
        return computeSParamsDifferential(freq, rlgc_odd, rlgc_even, length, Z_ref);
    }
    const sym = (o, e) => [[(o + e) / 2, (e - o) / 2], [(e - o) / 2, (o + e) / 2]];
    const R2 = sym(rlgc_odd.R, rlgc_even.R);
    const G2 = sym(rlgc_odd.G, rlgc_even.G);
    return computeSParamsDifferentialMTL(freq, R2, physMatrix.L, G2, physMatrix.C, length, Z_ref);
}

// 2×2 matrix transpose (Matrix2x2 stores a,b,c,d = [[a,b],[c,d]]).
function mTranspose(M) { return new Matrix2x2(M.a, M.c, M.b, M.d); }
// Build a complex 2×2 from real part [[..]] and imaginary part [[..]].
function mComplex(re, im) {
    return new Matrix2x2(
        new Complex(re[0][0], im[0][0]), new Complex(re[0][1], im[0][1]),
        new Complex(re[1][0], im[1][0]), new Complex(re[1][1], im[1][1]));
}

/**
 * Full multiconductor (MTL) 4-port S-parameters for a coupled pair from the physical 2×2
 * per-unit-length matrices [R],[L],[G],[C]. Unlike computeSParamsDifferential (which combines
 * odd/even single-ended responses and is exact only for a SYMMETRIC pair), this handles an
 * asymmetric pair: it yields S11≠S22 and non-zero mode-conversion SDC/SCD.
 *
 * Method: [Z]=[R]+jω[L], [Y]=[G]+jω[C]; [γ]=√([Z][Y]); chain blocks A=cosh(γℓ),
 * B=[Z]·γ⁻¹sinh(γℓ), C=[Y]·γ⁻¹sinh(γℓ), D=Aᵀ (reciprocal); block ABCD→4-port S; mixed-mode
 * transform → SDD/SCC/SDC/SCD. Port order: {1,2}=near ends of trace1/2, {3,4}=far ends.
 *
 * @param {number} freq Hz; @param {number[][]} R2,L2,G2,C2 physical 2×2 p.u.l. matrices (SI)
 * @param {number} length m; @param {number} Z_ref Ω
 */
function computeSParamsDifferentialMTL(freq, R2, L2, G2, C2, length, Z_ref) {
    const omega = 2 * Math.PI * freq;
    const Z = mComplex(R2, L2.map(row => row.map(v => v * omega)));   // [R] + jω[L]
    const Y = mComplex(G2, C2.map(row => row.map(v => v * omega)));   // [G] + jω[C]
    const gamma = Z.mul(Y).sqrt();                                     // [γ] = √([Z][Y])
    const gl = gamma.mul(length);
    const sinch = gamma.inv().mul(gl.sinh());                          // γ⁻¹·sinh(γℓ)
    const A = gl.cosh();
    const B = Z.mul(sinch);
    const Cb = Y.mul(sinch);
    const D = mTranspose(A);                                           // reciprocal: D = Aᵀ
    // Block ABCD → S (post-multiply by den⁻¹), Z_ref·I on every port.
    const Y0 = 1 / Z_ref;
    const BY = B.mul(Y0), CZ = Cb.mul(Z_ref);
    const den = A.add(BY).add(CZ).add(D);
    const denI = den.inv();
    const S11 = A.add(BY).sub(CZ).sub(D).mul(denI);                   // near-near
    const S22 = A.mul(-1).add(BY).sub(CZ).add(D).mul(denI);          // far-far
    const S21 = denI.mul(2);                                           // far-near (reciprocal, det≈1)
    const S12 = S21;
    // Assemble single-ended 4×4 (ports 1,2 near; 3,4 far), then mixed-mode transform.
    const blk = (m) => [[m.a, m.b], [m.c, m.d]];
    const b11 = blk(S11), b12 = blk(S12), b21 = blk(S21), b22 = blk(S22);
    const S = [
        [b11[0][0], b11[0][1], b12[0][0], b12[0][1]],
        [b11[1][0], b11[1][1], b12[1][0], b12[1][1]],
        [b21[0][0], b21[0][1], b22[0][0], b22[0][1]],
        [b21[1][0], b21[1][1], b22[1][0], b22[1][1]],
    ];
    // M: [V1,V2,V3,V4] → [d1,d2,c1,c2], rows = (1/√2)·{[1,-1,0,0],[0,0,1,-1],[1,1,0,0],[0,0,1,1]}.
    const r = 1 / Math.SQRT2;
    const M = [[r, -r, 0, 0], [0, 0, r, -r], [r, r, 0, 0], [0, 0, r, r]];
    const mm = Array.from({ length: 4 }, () => Array.from({ length: 4 }, () => new Complex(0, 0)));
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
        let s = new Complex(0, 0);
        for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++)
            s = s.add(S[a][b].mul(M[i][a] * M[j][b]));                // (M·S·Mᵀ)_ij
        mm[i][j] = s;
    }
    // mm order [d1,d2,c1,c2].
    return {
        S,
        SDD11: mm[0][0], SDD21: mm[1][0],
        SCC11: mm[2][2], SCC21: mm[3][2],
        SDC11: mm[0][2], SDC21: mm[1][2],
        SCD11: mm[2][0], SCD21: mm[3][0],
    };
}

/**
 * Convert Complex S-parameter to dB magnitude
 * @param {Complex} s - S-parameter
 * @returns {number} - Magnitude in dB
 */
function sParamTodB(s) {
    const mag = s.abs();
    if (mag < 1e-15) return -300;  // Prevent log(0)
    return 20 * Math.log10(mag);
}

/**
 * Convert Complex S-parameter to phase in degrees
 * @param {Complex} s - S-parameter
 * @returns {number} - Phase in degrees
 */
function sParamToPhase(s) {
    return s.arg() * 180 / Math.PI;
}

export {
    computeSParamsSingleEnded,
    computeSParamsDifferential,
    computeSParamsDifferentialMTL,
    computeSParamsDiffAuto,
    computeZ0,
    sParamTodB,
    sParamToPhase
};
