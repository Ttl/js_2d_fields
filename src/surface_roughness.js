import { Complex } from "./complex.js";

const EP0 = 8.854187818814e-12;
const MU0 = 4 * Math.PI * 1e-7;
const C0 = 299792458.0;

/**
 * Computes complex surface impedance using the Gradient Model (Rational Approximation)
 * Reference:
 * D. N. Grujić, "Simple and Accurate Approximation of Rough Conductor Surface
 * Impedance," IEEE Trans. Microwave Theory Tech., vol. 70, no. 4, pp.
 * 2053-2059, April 2022.
 * Implementation is based on:
 * https://github.com/simonp0420/MetalSurfaceImpedance.jl.
 * @param {number} f - Frequency in Hz
 * @param {number} sigma - Bulk conductivity (S/m)
 * @param {number} Rq - RMS Surface roughness (m)
 * @returns {Complex} Complex surface impedance (Re + jIm)
 */
function calculate_Zrough(f, sigma, Rq) {
    // 1. Smooth Case
    const omega = 2 * Math.PI * f;
    const delta = Math.sqrt(2.0 / (omega * MU0 * sigma));
    const R_smooth = 1.0 / (sigma * delta);
    
    // If effectively smooth, return (1+j)*R_smooth
    if (Rq <= 1e-9) { 
        return new Complex(R_smooth, R_smooth); 
    }

    // 2. Gradient Model Constants (Normal Distribution - "Oxide" side default)
    const fz = [8.655e7, 2.3039e9, 4.6915e13, 2.7795e14];
    const fp = [1.7702e9, 7.1614e13, 1.6413e16, 4.9260e12];
    const r_const = [0.50074, 0.45270, 0.43005, 0.29384];

    // 3. Scaling factors
    const Rq_ref = 1e-6;
    const sigma_ref = 58e6;
    const lambda_scale = (Rq * Rq * sigma) / (Rq_ref * Rq_ref * sigma_ref);
    
    const f_ref = lambda_scale * f;
    const omega_ref = 2 * Math.PI * f_ref;
    const delta_ref = Math.sqrt(2.0 / (omega_ref * MU0 * sigma_ref));
    const R_smooth_ref = 1.0 / (sigma_ref * delta_ref);
    
    // Z_smooth_ref = R_smooth_ref + j*R_smooth_ref
    const Z_smooth_ref = new Complex(R_smooth_ref, R_smooth_ref);

    // 4. Compute Psi (Correction Factor)
    // Psi = Product [ (1 + (j*f_ref/fzn)^rn) / (1 + (j*f_ref/fpn)^rn) ]
    let Psi = new Complex(1.0, 0.0);

    for (let k = 0; k < 4; k++) {
        // Term: (j * f_ref / freq)^r
        // This is (f_ref/freq)^r * (j)^r = (ratio)^r * exp(j * pi/2 * r)
        
        // Zero term (numerator)
        const ratio_z = f_ref / fz[k];
        const mag_z = Math.pow(ratio_z, r_const[k]);
        const ang_z = (Math.PI / 2.0) * r_const[k];
        const term_z = new Complex(
            1.0 + mag_z * Math.cos(ang_z), 
            mag_z * Math.sin(ang_z)
        );

        // Pole term (denominator)
        const ratio_p = f_ref / fp[k];
        const mag_p = Math.pow(ratio_p, r_const[k]);
        const ang_p = (Math.PI / 2.0) * r_const[k];
        const term_p = new Complex(
            1.0 + mag_p * Math.cos(ang_p), 
            mag_p * Math.sin(ang_p)
        );

        Psi = Psi.mul(term_z.div(term_p));
    }

    // 5. Final Z_rough = (Rq / (Rq_ref * lambda)) * Psi * Z_smooth_ref
    const scale = Rq / (Rq_ref * lambda_scale);
    return Psi.mul(Z_smooth_ref).mul(scale);
}

// Abramowitz & Stegun approximation 7.1.28, max error ~1.5e-7
function _erf(x) {
    const sign = x >= 0 ? 1 : -1;
    const a = Math.abs(x);
    const t = 1.0 / (1.0 + 0.3275911 * a);
    const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    return sign * (1.0 - poly * Math.exp(-a * a));
}

function _gaussianCDF(x, mean, sigma) {
    return 0.5 * (1.0 + _erf((x - mean) / (Math.SQRT2 * sigma)));
}

/**
 * Layered gradient model using transmission line taper approach.
 * Based on the method described in [1] and generalized for multiple layers in [2].
 * This is a faster and more accurate alternative to the ODE solver method.
 *
 * References:
 * [1] B. Tegowski, T. Jaschke, A. Sieganschin and A. F. Jacob,
 * "A Transmission Line Approach for Rough Conductor Surface Impedance Analysis,"
 * IEEE Trans. Microwave Theory Tech., vol. 71, no. 2, pp. 471-479, Feb. 2023.
 *
 * [2] G. Gold and K. Helmreich, "Modeling of transmission lines with multiple coated conductors,"
 * 2016 46th European Microwave Conference (EuMC), London, UK, 2016, pp. 635-638.
 *
 * @param {number} f - Frequency in Hz
 * @param {number} sigma_bulk - Bulk conductor conductivity (S/m)
 * @param {number} rq - RMS roughness at all interfaces (m)
 * @param {number} sigma_plating - Plating layer conductivity (S/m)
 * @param {number} thickness_plating - Plating layer thickness (m)
 * @param {number} N - Number of points for recursion (default 2048)
 * @returns {Complex} Complex surface impedance
 */
function calculate_Zrough_layered(f, sigma_bulk, rq, sigma_plating, thickness_plating, N = 2048) {
    // Fallback to single-layer if not layered
    if (thickness_plating <= 0 || sigma_plating <= 0) {
        return calculate_Zrough(f, sigma_bulk, rq);
    }

    const omega = 2 * Math.PI * f;

    // Skin depth of the less-conductive material for domain sizing
    const min_sigma = Math.min(sigma_bulk, sigma_plating);
    const skin_depth = Math.sqrt(2.0 / (omega * MU0 * min_sigma));

    // Recursion span: from -5*rq (into air) to well past skin depth
    const epsilon = 1e-9; // Handle rq=0 case
    const recursion_min = -5 * rq - epsilon;
    const recursion_max = Math.max(thickness_plating + 10 * skin_depth, 5e-6) + epsilon;

    // Uniform grid spacing
    const dx = (recursion_max - recursion_min) / (N - 1);

    // Build conductivity profile using CDF approach
    // Two interfaces: air/plating at x=0 and plating/bulk at x=thickness_plating
    const sigma_profile = new Float64Array(N);

    for (let k = 0; k < N; k++) {
        const x_k = recursion_min + k * dx;
        let cdf0, cdf1;

        if (rq <= 1e-12) {
            // Step profile: no roughness smoothing
            cdf0 = x_k >= 0 ? 1.0 : 0.0;
            cdf1 = x_k >= thickness_plating ? 1.0 : 0.0;
        } else {
            // Gaussian CDF for smooth transitions
            cdf0 = _gaussianCDF(x_k, 0, rq);
            cdf1 = _gaussianCDF(x_k, thickness_plating, rq);
        }

        // Build profile from CDFs
        // Region probabilities: air (sigma=0) -> plating -> bulk
        const p_plating = Math.max(0, cdf0 - cdf1);
        const p_bulk = cdf1;

        sigma_profile[k] = sigma_plating * p_plating + sigma_bulk * p_bulk;
    }

    // Compute transmission line properties
    const gamma = new Array(N);
    const Z = new Array(N);

    for (let k = 0; k < N; k++) {
        // Permittivity from conductivity
        const ep = new Complex(EP0, -sigma_profile[k] / omega);

        // Propagation constant
        let g = ep.mul(-omega * omega * MU0).sqrt();
        if (g.re < 0) g = g.neg();  // Ensure positive real part
        gamma[k] = g;

        // Characteristic impedance
        let z = new Complex(MU0, 0).div(ep).sqrt();
        if (z.re < 0) z = z.neg();  // Ensure positive real part
        Z[k] = z;
    }

    // Transmission line recursion (from last to first)
    // Zsi_new = z * (Zsi + z*tanh(g*dx)) / (z + Zsi*tanh(g*dx))
    let Zsi = Z[N - 1];

    for (let k = N - 1; k >= 0; k--) {
        const g = gamma[k];
        const z = Z[k];
        const tanh_gdx = g.mul(dx).tanh();
        const z_tanh = z.mul(tanh_gdx);

        Zsi = z.mul(Zsi.add(z_tanh)).div(z.add(Zsi.mul(tanh_gdx)));
    }

    return Zsi;
}

export { calculate_Zrough, calculate_Zrough_layered };
