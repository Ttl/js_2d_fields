/**
 * Djordjevic-Sarkar (Wideband Debye) causal material model
 *
 * Calculates frequency-dependent complex permittivity using the Djordjevic-Sarkar model
 * to enforce causality in the material parameters. This model adjusts both the real part
 * of permittivity and the loss tangent based on frequency.
 *
 * Reference: Djordjevic, A. R., Bazdar, M. B., Harrington, R. F., & Sarkar, T. K. (1987).
 * "LINPAR for Windows: Matrix Parameters for Multiconductor Transmission Lines,
 * Software and User's Manual".
 */

/**
 * Calculates the frequency-dependent complex permittivity using the
 * Djordjevic-Sarkar (Wideband Debye) model.
 *
 * @param {number} freq - Simulation frequency (Hz)
 * @param {number} eps_ref - Nominal relative permittivity (at f_ref)
 * @param {number} tand_ref - Nominal loss tangent (at f_ref)
 * @param {number} f_ref - Reference frequency (default 1 GHz)
 * @returns {{eps_real: number, tand_actual: number}} - Frequency-dependent permittivity and loss tangent
 */
export function djordjevic_sarkar(freq, eps_ref, tand_ref, f_ref = 1e9) {
    // 1. Define Hidden Constants (Ansys/HFSS Defaults)
    const f_low = 1e3;    // 1 kHz
    const f_high = 1e12;  // 1 THz

    // Prevent potential divide-by-zero at DC
    if (freq < f_low) {
        freq = f_low;
    }

    // 2. Calculate Slope Constant (K)
    const K = (eps_ref * tand_ref * 2) / Math.PI;

    // 3. Calculate Optical Permittivity (eps_inf)
    // We use the reference frequency to anchor the curve to eps_ref
    const term_num_ref = Math.sqrt(f_high**2 + f_ref**2);
    const term_den_ref = Math.sqrt(f_low**2 + f_ref**2);
    const eps_inf = eps_ref - K * Math.log(term_num_ref / term_den_ref);

    // 4. Calculate Complex Log Term at simulation frequency
    // We need the natural log of the complex ratio: ln( (fh + jf) / (fl + jf) )
    const numerator_real = f_high;
    const numerator_imag = freq;
    const denominator_real = f_low;
    const denominator_imag = freq;

    // Complex division: (a + jb) / (c + jd) = ((ac + bd) + j(bc - ad)) / (c^2 + d^2)
    const denom_mag_sq = denominator_real**2 + denominator_imag**2;
    const ratio_real = (numerator_real * denominator_real + numerator_imag * denominator_imag) / denom_mag_sq;
    const ratio_imag = (numerator_imag * denominator_real - numerator_real * denominator_imag) / denom_mag_sq;

    // Complex logarithm: ln(a + jb) = ln(|z|) + j*arg(z)
    const ratio_mag = Math.sqrt(ratio_real**2 + ratio_imag**2);
    const ratio_arg = Math.atan2(ratio_imag, ratio_real);

    const log_real = Math.log(ratio_mag);
    const log_imag = ratio_arg;

    // 5. Final Complex Permittivity: eps_complex = eps_inf + K * log_term
    const eps_complex_real = eps_inf + K * log_real;
    const eps_complex_imag = K * log_imag;

    // Extract results
    const eps_real = eps_complex_real;

    // Loss tangent is ratio of imaginary to real.
    // Note: In most physics conventions, lossy part is negative imaginary.
    // We return positive tand for readability.
    const tand_actual = -eps_complex_imag / eps_complex_real;

    return { eps_real, tand_actual };
}

/**
 * Apply Djordjevic-Sarkar model to all dielectric materials in the solver.
 * Modifies epsilon_r and tand arrays in place based on frequency.
 *
 * @param {FieldSolver2D} solver - The field solver instance
 * @param {number} f_ref - Reference frequency for material parameters (default 1 GHz)
 */
export function applyDjordjevicSarkar(solver, f_ref = 1e9) {
    if (!solver.epsilon_r || !solver.tand) {
        throw new Error("Solver must have epsilon_r and tand arrays initialized");
    }

    const freq = solver.freq;
    const ny = solver.y.length;
    const nx = solver.x.length;

    // Store original values as reference parameters
    // Re-initialize if dimensions have changed (due to mesh refinement)
    const needsInit = !solver._original_epsilon_r ||
                      solver._original_epsilon_r.length !== ny ||
                      (ny > 0 && (!solver._original_epsilon_r[0] || solver._original_epsilon_r[0].length !== nx));

    if (needsInit) {
        solver._original_epsilon_r = solver.epsilon_r.map(row => new Float64Array(row));
        solver._original_tand = solver.tand.map(row => new Float64Array(row));
    }

    // Apply model to each mesh point
    for (let i = 0; i < ny; i++) {
        for (let j = 0; j < nx; j++) {
            const eps_ref = solver._original_epsilon_r[i][j];
            const tand_ref = solver._original_tand[i][j];

            // Skip air/vacuum regions (eps_r = 1)
            if (Math.abs(eps_ref - 1.0) < 1e-6) {
                continue;
            }

            // Skip if loss tangent is zero (lossless material)
            if (Math.abs(tand_ref) < 1e-10) {
                continue;
            }

            // Apply Djordjevic-Sarkar model
            const { eps_real, tand_actual } = djordjevic_sarkar(freq, eps_ref, tand_ref, f_ref);

            solver.epsilon_r[i][j] = eps_real;
            solver.tand[i][j] = tand_actual;
        }
    }
}
