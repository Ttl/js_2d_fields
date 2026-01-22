import { Complex } from "./complex.js";

const MU0 = 1.25663706212e-6;
const C0 = 299792458.0;

/**
 * Computes complex surface impedance using the Gradient Model (Rational Approximation)
 * based on Grujić et al. (IEEE T-MTT 2022).
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

export { calculate_Zrough };
