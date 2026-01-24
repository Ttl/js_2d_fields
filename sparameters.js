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
