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
    const Z = Psi.mul(Z_smooth_ref).mul(scale);
    // The excess inductance (Im(Z) - Rs) / omega of the fitted model peaks at
    // f_ref = ROUGH_LX_PEAK and decays below it, negative below about 1 Hz. The
    // roughness is then a thin surface layer far inside the skin depth, whose excess
    // inductance is a constant: it is held at the peak, which keeps L(f) monotone.
    if (f_ref < ROUGH_LX_PEAK && omega > 0) {
        const fp = ROUGH_LX_PEAK / lambda_scale, wp = 2 * Math.PI * fp;
        const Zp = calculate_Zrough(fp, sigma, Rq);
        const Lx = (Zp.im - Math.sqrt(wp * MU0 / (2 * sigma))) / wp;
        return new Complex(Z.re, R_smooth + omega * Lx);
    }
    return Z;
}
const ROUGH_LX_PEAK = 1.2e6;

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

// Antiderivative of the Gaussian CDF:
//   G(x) = ∫ Φ((t−m)/s) dt = (x−m)·Φ((x−m)/s) + s·φ((x−m)/s)
// with φ the standard normal density. At s → 0 the CDF becomes a step and G degenerates
// to the ramp max(x−m, 0), so one expression covers both the roughened interface and the
// perfectly sharp one.
function _cdfIntegral(x, mean, sigma) {
    const d = x - mean;
    if (sigma <= 0) return Math.max(d, 0);
    const z = d / sigma;
    return d * _gaussianCDF(x, mean, sigma)
         + sigma * Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
}

// Mean of that CDF over [a, b], the fraction of the segment lying past the boundary.
// This is what makes the profile below a partial-volume average rather than a point
// sample, so a boundary falling mid-segment is placed exactly instead of snapping to the
// grid. Segments wholly clear of the transition short-circuit to the constant 0 or 1:
// there the antiderivatives are large and nearly equal, and differencing them would both
// cancel and inherit the erf fit's ~1e-7 absolute error.
function _cdfMean(a, b, mean, sigma) {
    if (sigma > 0) {
        if (b <= mean - 6 * sigma) return 0;
        if (a >= mean + 6 * sigma) return 1;
    } else {
        if (b <= mean) return 0;
        if (a >= mean) return 1;
    }
    return (_cdfIntegral(b, mean, sigma) - _cdfIntegral(a, mean, sigma)) / (b - a);
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

    // Recursion span: from 5*rq above the mean surface (out in the air, where the graded
    // profile has died away) to well past the skin depth. Zs is referenced to the top of
    // this span, so that top must sit at the surface for a sharp interface.
    const recursion_min = -5 * rq;
    const recursion_max = Math.max(thickness_plating + 10 * skin_depth, 5e-6);

    // Uniform grid spacing
    const dx = (recursion_max - recursion_min) / (N - 1);

    // Build the conductivity profile from the two interface CDFs: air/plating at x = 0
    // and plating/bulk at x = thickness_plating.
    //
    // Each entry is the segment's partial-volume average over [x_k, x_k + dx], not the
    // profile sampled at x_k, and the recursion below consumes it as a uniform slab of
    // exactly that span. Point sampling made both interfaces snap to the nearest grid
    // node, which cost two distinct errors: the layer thickness was quantized by ±dx
    // (up to 0.48% on Re(Zs), hence on loss, worst for a thin smooth plating), and the
    // topmost segment was read as all-vacuum, standing the whole stack off by up to
    // ~80 nm and inflating Im(Zs), hence the internal inductance. Averaging places both
    // boundaries exactly wherever they fall, so neither survives grid alignment.
    //
    // Averaging the CDFs is enough to average sigma: with thickness_plating > 0 the
    // deeper CDF never exceeds the shallower one, so sigma is a fixed linear combination
    // of the two and the mean passes straight through it.
    const sigma_profile = new Float64Array(N);
    const s_rough = rq <= 1e-12 ? 0 : rq;   // 0 selects the exact-step branch of _cdfMean

    for (let k = 0; k < N; k++) {
        const xa = recursion_min + k * dx, xb = xa + dx;
        const cdf0 = _cdfMean(xa, xb, 0, s_rough);
        const cdf1 = _cdfMean(xa, xb, thickness_plating, s_rough);

        // Region fractions: air (sigma = 0) -> plating -> bulk
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

// Lateral spreading of the return current in a thin wall. A line current over an
// infinitely wide resistive sheet has the sheet current K(k) = -I e^{-|k|h} / (1 - j k Lambda)
// in the transverse wavenumber k, Lambda = delta^2 / d the lateral diffusion length of a
// sheet of thickness d << delta. The dissipation relative to the PEC distribution
// (Lambda -> 0) is g(u) = 2 * integral_0^inf e^{-2s} / (1 + u^2 s^2) ds with u = Lambda / h,
// which is 1 - u^2/2 for small u and pi/u for large u (R ~ 1/(sigma delta^2), vanishing at
// DC, where the PEC distribution would leave the sheet resistance over the PEC current
// width). Written for a general PEC distribution through its effective width
// W_K = I^2 / integral(|K|^2 dl) = 2 pi h for the line source, so u = 2 pi Lambda / W_K.
// Evaluated with t = tan(theta) on theta in [0, pi/2), where the integrand is smooth
// for every u.
function wallSpreadFactor(u) {
    if (!(u > 0.05)) return 1 - u * u / 2;
    const a = 2 / u, n = 2000, h = (Math.PI / 2) / n;
    let acc = 0;
    for (let i = 0; i <= n; i++) {
        const th = i * h;
        const v = i < n ? Math.exp(-a * Math.tan(th)) : 0;
        acc += (i === 0 || i === n) ? v : (i % 2 ? 4 * v : 2 * v);
    }
    return (2 / u) * acc * h / 3;
}

// Blend range of the spreading parameter u = delta^2 / (d W_K) of a ground: below
// SPREAD_U_MIN the confined surface model holds, a decade above it the spread-current
// model (the thin-sheet ground solve, the ideal-ground MQS solve) holds alone, and in
// between the two are blended on log10(u) with this weight of the spreading model.
const SPREAD_U_MIN = 0.02;
function spreadBlendWeight(u) {
    return Math.min(1, Math.max(0, Math.log10(u / SPREAD_U_MIN)));
}

// Lower bound on the return width W_K of the grounds: a quarter of the narrowest trace.
function spreadWidthFloor(traceWidths) {
    return 0.25 * Math.min(...traceWidths);
}

// (1+j) coth((1+j) x) as {re, im}: the surface impedance over Rs of a slab x skin
// depths thick with a field-free back. It is 1+j for a thick slab and tends to
// 1/x + j 2x/3 for a thin one. Taken as 1+j past x = 20, where that holds to double
// precision (cosh(2x) overflows past x ~ 355), and for x <= 0.
function slabCoth(x) {
    if (!(x > 0) || x > 20) return { re: 1, im: 1 };
    const den = Math.cosh(2 * x) - Math.cos(2 * x);
    return { re: (Math.sinh(2 * x) + Math.sin(2 * x)) / den, im: (Math.sinh(2 * x) - Math.sin(2 * x)) / den };
}

export { calculate_Zrough, calculate_Zrough_layered, wallSpreadFactor, slabCoth, SPREAD_U_MIN, spreadBlendWeight,
    spreadWidthFloor };
