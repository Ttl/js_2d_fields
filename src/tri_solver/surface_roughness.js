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
const EP0 = 8.854187818814e-12;

// Complex helpers on {re, im}
const cmul = (a, b) => ({ re: a.re*b.re - a.im*b.im, im: a.re*b.im + a.im*b.re });
function cdiv(a, b) {
    const d = b.re*b.re + b.im*b.im;
    return { re: (a.re*b.re + a.im*b.im) / d, im: (a.im*b.re - a.re*b.im) / d };
}
const cadd = (a, b) => ({ re: a.re + b.re, im: a.im + b.im });
const cscale = (a, s) => ({ re: a.re * s, im: a.im * s });
const cneg = (a) => ({ re: -a.re, im: -a.im });
function csqrt(a) {
    const r = Math.hypot(a.re, a.im);
    let re = Math.sqrt(Math.max(0, (r + a.re) / 2));
    let im = Math.sqrt(Math.max(0, (r - a.re) / 2));
    if (a.im < 0) im = -im;
    return { re, im };
}
function ctanh(a) {
    // tanh(x+iy) = (sinh 2x + i sin 2y) / (cosh 2x + cos 2y)
    const x2 = 2 * a.re, y2 = 2 * a.im;
    const d = Math.cosh(x2) + Math.cos(y2);
    return { re: Math.sinh(x2) / d, im: Math.sin(y2) / d };
}
function _erf(x) {
    const sign = x >= 0 ? 1 : -1;
    const ax = Math.abs(x);
    const t = 1 / (1 + 0.3275911 * ax);
    const poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    return sign * (1 - poly * Math.exp(-ax * ax));
}
const _gaussianCDF = (x, mean, sigma) => 0.5 * (1 + _erf((x - mean) / (Math.SQRT2 * sigma)));

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

// Layered gradient model: a plating layer (sigma_plating, thickness) over the bulk
// metal, both sharing the roughness Rq (the plating conforms to the bulk surface).
// A transmission-line taper recursion through the depth-graded conductivity
// profile (air → plating → bulk), exactly mirroring the quasi-static backend's
// calculate_Zrough_layered so per-face plating matches between solvers. Returns
// {re, im} ohms/square. Falls back to the single-layer model when there is no
// real plating layer (thickness ≤ 0 or sigma_plating ≤ 0).
export function calculateZroughLayered(f, sigma_bulk, rq, sigma_plating, thickness_plating, N = 2048) {
    if (thickness_plating <= 0 || sigma_plating <= 0) return calculateZrough(f, sigma_bulk, rq);

    const omega = 2 * Math.PI * f;
    const min_sigma = Math.min(sigma_bulk, sigma_plating);
    const skin_depth = Math.sqrt(2 / (omega * MU0 * min_sigma));

    const epsilon = 1e-9;
    const recursion_min = -5 * rq - epsilon;
    const recursion_max = Math.max(thickness_plating + 10 * skin_depth, 5e-6) + epsilon;
    const dx = (recursion_max - recursion_min) / (N - 1);

    // Per-depth propagation constant γ and characteristic impedance Z from the
    // graded conductivity (air σ=0 → plating → bulk via two Gaussian-CDF interfaces).
    const gamma = new Array(N), Z = new Array(N);
    for (let k = 0; k < N; k++) {
        const x_k = recursion_min + k * dx;
        let cdf0, cdf1;
        if (rq <= 1e-12) {
            cdf0 = x_k >= 0 ? 1 : 0;
            cdf1 = x_k >= thickness_plating ? 1 : 0;
        } else {
            cdf0 = _gaussianCDF(x_k, 0, rq);
            cdf1 = _gaussianCDF(x_k, thickness_plating, rq);
        }
        const p_plating = Math.max(0, cdf0 - cdf1), p_bulk = cdf1;
        const sigma_k = sigma_plating * p_plating + sigma_bulk * p_bulk;

        const ep = { re: EP0, im: -sigma_k / omega };
        let g = csqrt(cscale(ep, -omega * omega * MU0));
        if (g.re < 0) g = cneg(g);
        gamma[k] = g;
        let z = csqrt(cdiv({ re: MU0, im: 0 }, ep));
        if (z.re < 0) z = cneg(z);
        Z[k] = z;
    }

    // TL recursion from the bulk back to the surface.
    let Zsi = Z[N - 1];
    for (let k = N - 1; k >= 0; k--) {
        const z = Z[k];
        const tgd = ctanh(cscale(gamma[k], dx));
        const zt = cmul(z, tgd);
        Zsi = cdiv(cmul(z, cadd(Zsi, zt)), cadd(z, cmul(Zsi, tgd)));
    }
    return Zsi;
}
