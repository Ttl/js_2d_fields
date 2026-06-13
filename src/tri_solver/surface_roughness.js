// Gradient-model rough-conductor surface impedance (single layer).
//
// Rational approximation of the Gold–Helmreich gradient model:
//   D. N. Grujić, "Simple and Accurate Approximation of Rough Conductor Surface
//   Impedance," IEEE Trans. Microwave Theory Tech., vol. 70, no. 4, 2022.
// Ported from js_2d_fields/src/surface_roughness.js (calculate_Zrough), which
// follows https://github.com/simonp0420/MetalSurfaceImpedance.jl.
//
// The gradient model represents roughness as a depth-graded conductivity
// σ(d) = σ_bulk·Φ(d/Rq) (Gaussian CDF of the surface height distribution) and
// this fit returns the resulting effective surface impedance.
//
// Returns {re, im} in ohms/square. For Rq → 0 this reduces to the smooth
// result Zs = (1+j)·Rs with Rs = 1/(σδ).

const MU0 = 4 * Math.PI * 1e-7;

// Complex helpers on {re, im}
const cmul = (a, b) => ({ re: a.re*b.re - a.im*b.im, im: a.re*b.im + a.im*b.re });
function cdiv(a, b) {
    const d = b.re*b.re + b.im*b.im;
    return { re: (a.re*b.re + a.im*b.im) / d, im: (a.im*b.re - a.re*b.im) / d };
}

export function calculateZrough(f, sigma, Rq) {
    const omega = 2 * Math.PI * f;
    const delta = Math.sqrt(2 / (omega * MU0 * sigma));
    const R_smooth = 1 / (sigma * delta);

    if (!Rq || Rq <= 1e-9) return { re: R_smooth, im: R_smooth };

    // Rational fit constants (normal height distribution, "oxide" side)
    const fz = [8.655e7, 2.3039e9, 4.6915e13, 2.7795e14];
    const fp = [1.7702e9, 7.1614e13, 1.6413e16, 4.9260e12];
    const rc = [0.50074, 0.45270, 0.43005, 0.29384];

    // Scale to the reference case (Rq_ref = 1 µm, σ_ref = 58 MS/m)
    const Rq_ref = 1e-6, sigma_ref = 58e6;
    const lambda = (Rq * Rq * sigma) / (Rq_ref * Rq_ref * sigma_ref);
    const f_ref = lambda * f;
    const delta_ref = Math.sqrt(2 / (2 * Math.PI * f_ref * MU0 * sigma_ref));
    const R_ref = 1 / (sigma_ref * delta_ref);

    // Ψ = Π (1 + (j·f_ref/fz_k)^r_k) / (1 + (j·f_ref/fp_k)^r_k)
    let Psi = { re: 1, im: 0 };
    for (let k = 0; k < 4; k++) {
        const ang = (Math.PI / 2) * rc[k];
        const mz = Math.pow(f_ref / fz[k], rc[k]);
        const mp = Math.pow(f_ref / fp[k], rc[k]);
        const num = { re: 1 + mz * Math.cos(ang), im: mz * Math.sin(ang) };
        const den = { re: 1 + mp * Math.cos(ang), im: mp * Math.sin(ang) };
        Psi = cmul(Psi, cdiv(num, den));
    }

    const scale = Rq / (Rq_ref * lambda);
    return cmul(Psi, { re: R_ref * scale, im: R_ref * scale });
}
