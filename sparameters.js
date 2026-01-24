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
    // Compute modal characteristic impedances
    const Z0_odd = computeZ0(freq, rlgc_odd);
    const Z0_even = computeZ0(freq, rlgc_even);

    // Compute single-ended S-params for each mode using modal impedances as reference
    const S_odd = computeSParamsSingleEnded(freq, rlgc_odd, length, Z0_odd);
    const S_even = computeSParamsSingleEnded(freq, rlgc_even, length, Z0_even);

    // For ideal symmetric differential pairs with no coupling between modes,
    // the 4-port S-matrix can be constructed from odd and even mode responses.
    //
    // Mixed-mode S-parameters:
    // SDD = (S_odd) - differential mode
    // SCC = (S_even) - common mode
    // SDC = 0 (no mode conversion for symmetric line)
    // SCD = 0

    // The full 4-port S-matrix (standard single-ended ports):
    // Port assignment: 1=in+, 2=out+, 3=in-, 4=out-
    //
    // For symmetric coupled lines:
    // S11 = S33 = (S_odd_11 + S_even_11) / 2
    // S22 = S44 = (S_odd_22 + S_even_22) / 2
    // S21 = S43 = (S_odd_21 + S_even_21) / 2
    // S12 = S34 = (S_odd_12 + S_even_12) / 2
    // S13 = S31 = (S_even_11 - S_odd_11) / 2
    // S24 = S42 = (S_even_22 - S_odd_22) / 2
    // S23 = S41 = (S_even_21 - S_odd_21) / 2
    // S14 = S32 = (S_even_12 - S_odd_12) / 2

    const half = new Complex(0.5, 0);

    // Calculate 4-port S-parameters
    const S11 = S_odd.S11.add(S_even.S11).mul(half);
    const S22 = S_odd.S22.add(S_even.S22).mul(half);
    const S21 = S_odd.S21.add(S_even.S21).mul(half);
    const S12 = S_odd.S12.add(S_even.S12).mul(half);
    const S13 = S_even.S11.sub(S_odd.S11).mul(half);
    const S24 = S_even.S22.sub(S_odd.S22).mul(half);
    const S23 = S_even.S21.sub(S_odd.S21).mul(half);
    const S14 = S_even.S12.sub(S_odd.S12).mul(half);

    // Build 4x4 matrix (symmetric coupled line)
    // Port assignment: 1=in+, 2=out+, 3=in-, 4=out-
    const S = [
        [S11, S12, S13, S14],  // Row 1: S11, S12, S13, S14
        [S21, S22, S23, S24],  // Row 2: S21, S22, S23, S24
        [S13, S14, S11, S12],  // Row 3: S31, S32, S33, S34
        [S23, S24, S21, S22]   // Row 4: S41, S42, S43, S44
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
    computeZ0,
    sParamTodB,
    sParamToPhase
};
