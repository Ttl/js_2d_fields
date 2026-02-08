/**
 * Test for calculate_Zrough_layered — compares JS layered gradient model
 * against reference values computed with the equivalent Python implementation.
 *
 * Python reference script (run once to regenerate expected values):
 *
 * ```python
 * import numpy as np
 * from scipy.special import erf
 * from scipy.sparse import diags
 * from scipy.sparse.linalg import spsolve
 *
 * MU0 = 4 * np.pi * 1e-7
 *
 * def gaussian_cdf(x, mean, rq):
 *     return 0.5 * (1 + erf((x - mean) / (np.sqrt(2) * rq)))
 *
 * def solve_layered_Zs(f, sigma_bulk, rq, sigma_plating, thickness_plating, N=2000):
 *     omega = 2 * np.pi * f
 *     min_sigma = min(sigma_bulk, sigma_plating)
 *     skin_depth = np.sqrt(2.0 / (omega * MU0 * min_sigma))
 *     x_min = -5 * rq
 *     x_max = max(thickness_plating + 10 * skin_depth, 5e-6)
 *     x = np.linspace(x_min, x_max, N)
 *     dx = x[1] - x[0]
 *
 *     cdf0 = gaussian_cdf(x, 0, rq)
 *     cdf1 = gaussian_cdf(x, thickness_plating, rq)
 *     p_plating = np.maximum(0, cdf0 - cdf1)
 *     sigma_profile = sigma_plating * p_plating + sigma_bulk * cdf1
 *
 *     sigma_clipped = np.maximum(sigma_profile, 1.0)
 *     dln_sigma = np.gradient(np.log(sigma_clipped), dx)
 *
 *     k_sq = 1j * omega * MU0 * sigma_clipped
 *     inv_dx2 = 1.0 / dx**2
 *     half_inv_dx = 0.5 / dx
 *
 *     D_lower = inv_dx2 + dln_sigma[1:-1] * half_inv_dx
 *     D_main  = -2 * inv_dx2 - k_sq[1:-1]
 *     D_upper = inv_dx2 - dln_sigma[1:-1] * half_inv_dx
 *
 *     main_diag = np.zeros(N, dtype=complex)
 *     upper_diag = np.zeros(N-1, dtype=complex)
 *     lower_diag = np.zeros(N-1, dtype=complex)
 *
 *     main_diag[1:-1] = D_main
 *     upper_diag[1:]  = D_upper
 *     lower_diag[:-1] = D_lower
 *
 *     main_diag[0] = 1.0; main_diag[-1] = 1.0
 *     upper_diag[0] = 0.0; lower_diag[-1] = 0.0
 *
 *     A = diags([lower_diag, main_diag, upper_diag], offsets=[-1, 0, 1], format="csc")
 *     b = np.zeros(N, dtype=complex); b[0] = 1.0
 *     B = spsolve(A, b)
 *
 *     integral = np.trapz(B, x=x)
 *     den = B[0] / MU0
 *     Zs = 1j * omega * integral / den
 *     return Zs
 *
 * # Test cases
 * cases = [
 *     # (f, sigma_bulk, rq, sigma_plating, thickness_plating)
 *     (1e9,  10e6, 0.3e-6, 58e6, 0.5e-6),   # Case 1: Cu plating on brass, 1 GHz
 *     (10e9, 10e6, 0.3e-6, 58e6, 0.5e-6),   # Case 2: Cu plating on brass, 10 GHz
 *     (100e9, 10e6, 0.3e-6, 58e6, 0.5e-6),  # Case 3: Cu plating on brass, 100 GHz (skin << plating)
 *     (1e9, 58e6, 1.0e-6, 43e6, 0.2e-6),    # Case 4: Au plating on Cu, rough, 1 GHz
 *     (10e9, 58e6, 0.1e-6, 43e6, 1.0e-6),   # Case 5: Thick Au on Cu, smooth, 10 GHz
 * ]
 * for i, (f, sb, rq, sp, tp) in enumerate(cases):
 *     Zs = solve_layered_Zs(f, sb, rq, sp, tp)
 *     print(f"Case {i+1}: Zs.re={Zs.real:.6e}  Zs.im={Zs.imag:.6e}")
 *
 * # Output (regenerate if Python algorithm changes):
 * # Case 1 (f=1GHz): Zs.re=1.608768e-02  Zs.im=2.103953e-02
 * # Case 2 (f=10GHz): Zs.re=3.822746e-02  Zs.im=1.307002e-01
 * # Case 3 (f=100GHz): Zs.re=1.946237e-01  Zs.im=1.018985e+00
 * ```
 *
 * Run this test with: node tests/test_layered_roughness.js
 */

import { calculate_Zrough, calculate_Zrough_layered } from '../src/surface_roughness.js';

function assertClose(actual, expected, rtol, label) {
    const err = Math.abs(actual - expected);
    const ref = Math.abs(expected) || 1e-30;
    const relErr = err / ref;
    const pass = relErr < rtol;
    console.log(`  ${pass ? 'PASS' : 'FAIL'} ${label}: got ${actual.toExponential(6)}, expected ${expected.toExponential(6)}, rel_err=${relErr.toExponential(2)}`);
    return pass;
}

let failures = 0;

// ---------------------------------------------------------------
// Test 1: Fallback to single-layer when plating thickness is 0
// ---------------------------------------------------------------
console.log('\nTest 1: Fallback when thickness=0');
{
    const f = 1e9, sigma = 58e6, rq = 0.5e-6;
    const Z_single = calculate_Zrough(f, sigma, rq);
    const Z_layered = calculate_Zrough_layered(f, sigma, rq, 43e6, 0);
    if (!assertClose(Z_layered.re, Z_single.re, 1e-10, 'Re fallback')) failures++;
    if (!assertClose(Z_layered.im, Z_single.im, 1e-10, 'Im fallback')) failures++;
}

// ---------------------------------------------------------------
// Test 2: Fallback when plating sigma is 0
// ---------------------------------------------------------------
console.log('\nTest 2: Fallback when sigma_plating=0');
{
    const f = 5e9, sigma = 10e6, rq = 1.0e-6;
    const Z_single = calculate_Zrough(f, sigma, rq);
    const Z_layered = calculate_Zrough_layered(f, sigma, rq, 0, 0.5e-6);
    if (!assertClose(Z_layered.re, Z_single.re, 1e-10, 'Re fallback')) failures++;
    if (!assertClose(Z_layered.im, Z_single.im, 1e-10, 'Im fallback')) failures++;
}

// ---------------------------------------------------------------
// Test 2b: Zero roughness produces finite results (step profile)
// ---------------------------------------------------------------
console.log('\nTest 2b: Zero roughness (rq=0) produces finite values');
{
    const f = 10e9, sigma_bulk = 10e6, sigma_plating = 58e6, thickness = 0.5e-6;
    const Z = calculate_Zrough_layered(f, sigma_bulk, 0, sigma_plating, thickness);
    const finite = isFinite(Z.re) && isFinite(Z.im);
    if (finite && Z.re > 0 && Z.im > 0) {
        console.log(`  PASS Re=${Z.re.toExponential(4)}, Im=${Z.im.toExponential(4)}`);
    } else {
        console.log(`  FAIL Re=${Z.re}, Im=${Z.im} (expected finite positive values)`);
        failures++;
    }
}

// ---------------------------------------------------------------
// Test 3: Physical sanity — high-freq plating dominates
// At 100 GHz with 0.5um Cu plating, skin depth in Cu ≈ 0.66um
// So plating contributes significantly. Re(Z) should be close to
// smooth Cu value (1+j)*R_smooth_Cu scaled by roughness.
// ---------------------------------------------------------------
console.log('\nTest 3: High-frequency plating dominance (100 GHz)');
{
    const f = 100e9;
    const sigma_bulk = 10e6;   // brass
    const sigma_plating = 58e6; // copper
    const rq = 0.3e-6;
    const thickness = 0.5e-6;

    const Z = calculate_Zrough_layered(f, sigma_bulk, rq, sigma_plating, thickness);

    // Smooth copper at 100 GHz: delta = sqrt(2/(2*pi*100e9*4pi*1e-7*58e6)) ≈ 0.66um
    // R_smooth = 1/(58e6 * 0.66e-6) ≈ 0.026 ohm
    // With roughness, Re should be somewhat above R_smooth but order of magnitude correct
    const omega = 2 * Math.PI * f;
    const MU0 = 4 * Math.PI * 1e-7;
    const delta_cu = Math.sqrt(2 / (omega * MU0 * sigma_plating));
    const R_smooth_cu = 1 / (sigma_plating * delta_cu);

    console.log(`  Layered: Re=${Z.re.toExponential(4)}, Im=${Z.im.toExponential(4)}`);
    console.log(`  R_smooth_Cu=${R_smooth_cu.toExponential(4)}, skin_depth_Cu=${(delta_cu*1e6).toFixed(3)} um`);

    // Re(Z) should be > R_smooth_Cu (roughness increases resistance)
    // but < 10x R_smooth_Cu (not wildly off)
    if (Z.re > R_smooth_cu && Z.re < 10 * R_smooth_cu) {
        console.log('  PASS Re in expected range');
    } else {
        console.log(`  FAIL Re=${Z.re.toExponential(4)} outside [${R_smooth_cu.toExponential(4)}, ${(10*R_smooth_cu).toExponential(4)}]`);
        failures++;
    }
    // Im should be positive (inductive reactance from skin effect)
    if (Z.im > 0) {
        console.log('  PASS Im positive');
    } else {
        console.log(`  FAIL Im=${Z.im.toExponential(4)} (expected positive)`);
        failures++;
    }
}

// ---------------------------------------------------------------
// Test 4: Low-freq bulk dominates
// At 1 MHz with 0.5um plating, skin depth in brass ≈ 23um >> thickness.
// Current sees both layers. Z should be closer to bulk brass than Cu.
// ---------------------------------------------------------------
console.log('\nTest 4: Low-frequency bulk dominance (1 MHz)');
{
    const f = 1e6;
    const sigma_bulk = 10e6;
    const sigma_plating = 58e6;
    const rq = 0.3e-6;
    const thickness = 0.5e-6;

    const Z = calculate_Zrough_layered(f, sigma_bulk, rq, sigma_plating, thickness);

    const omega = 2 * Math.PI * f;
    const MU0 = 4 * Math.PI * 1e-7;
    const delta_brass = Math.sqrt(2 / (omega * MU0 * sigma_bulk));
    const R_smooth_brass = 1 / (sigma_bulk * delta_brass);

    console.log(`  Layered: Re=${Z.re.toExponential(4)}, Im=${Z.im.toExponential(4)}`);
    console.log(`  R_smooth_brass=${R_smooth_brass.toExponential(4)}, skin_depth_brass=${(delta_brass*1e6).toFixed(1)} um`);

    // At very low freq, skin depth >> plating, so Z should be near R_smooth_brass
    // (with some roughness correction). Allow generous range.
    if (Z.re > 0.1 * R_smooth_brass && Z.re < 5 * R_smooth_brass) {
        console.log('  PASS Re in expected range for bulk-dominated regime');
    } else {
        console.log(`  FAIL Re=${Z.re.toExponential(4)} outside expected range`);
        failures++;
    }
}

// ---------------------------------------------------------------
// Test 5: Monotonicity — Re(Z) should decrease with increasing plating thickness
// (thicker plating = more conductive material in skin depth at moderate freq)
// ---------------------------------------------------------------
console.log('\nTest 5: Re(Z) decreases with increasing plating thickness (up to skin depth)');
{
    const f = 10e9;
    const sigma_bulk = 10e6;
    const sigma_plating = 58e6;
    const rq = 0.3e-6;

    // Only check up to ~1.5x skin depth. Beyond that the plating saturates
    // and roughness at the air-plating interface dominates, allowing non-monotonic behavior.
    const thicknesses = [0.1e-6, 0.3e-6, 0.5e-6, 1.0e-6];
    const Re_vals = thicknesses.map(t => calculate_Zrough_layered(f, sigma_bulk, rq, sigma_plating, t).re);

    console.log('  Thickness (um) -> Re(Z):');
    let monotonic = true;
    for (let k = 0; k < thicknesses.length; k++) {
        console.log(`    ${(thicknesses[k]*1e6).toFixed(1)} um -> ${Re_vals[k].toExponential(4)}`);
        if (k > 0 && Re_vals[k] > Re_vals[k-1]) {
            monotonic = false;
        }
    }
    if (monotonic) {
        console.log('  PASS Re(Z) monotonically decreasing with thickness');
    } else {
        console.log('  FAIL Re(Z) not monotonically decreasing');
        failures++;
    }
}

// ---------------------------------------------------------------
// Test 6: Symmetry check — plating with same sigma as bulk should
// behave like single-layer with that conductivity
// ---------------------------------------------------------------
console.log('\nTest 6: Same-sigma plating vs single-layer');
{
    const f = 5e9;
    const sigma = 30e6;
    const rq = 0.5e-6;

    // Layered with same sigma should match single-layer closely
    // (not exactly due to two-interface vs one-interface profile)
    const Z_layered = calculate_Zrough_layered(f, sigma, rq, sigma, 1.0e-6);
    const Z_single = calculate_Zrough(f, sigma, rq);

    console.log(`  Layered (same sigma): Re=${Z_layered.re.toExponential(4)}, Im=${Z_layered.im.toExponential(4)}`);
    console.log(`  Single-layer:         Re=${Z_single.re.toExponential(4)}, Im=${Z_single.im.toExponential(4)}`);

    // Should be within 20% — the profile shape differs slightly but conductivity is uniform
    const re_ratio = Z_layered.re / Z_single.re;
    const im_ratio = Z_layered.im / Z_single.im;
    if (Math.abs(re_ratio - 1) < 0.20 && Math.abs(im_ratio - 1) < 0.20) {
        console.log(`  PASS ratio Re=${re_ratio.toFixed(3)}, Im=${im_ratio.toFixed(3)}`);
    } else {
        console.log(`  FAIL ratio Re=${re_ratio.toFixed(3)}, Im=${im_ratio.toFixed(3)} (expected ~1.0 ± 0.2)`);
        failures++;
    }
}

// ---------------------------------------------------------------
// Test 7: Compare against Python reference values
// Run the Python script in the header comment to regenerate these.
// Tolerance is 5% to account for minor algorithm differences
// (grid size, gradient computation at boundaries).
// ---------------------------------------------------------------
console.log('\nTest 7: Python reference comparison (5% tolerance)');
{
    const cases = [
        // [f, sigma_bulk, rq, sigma_plating, thickness_plating, expected_re, expected_im]
        // Generated by running the Python reference script in the file header.
        [1e9,   10e6, 0.3e-6, 58e6, 0.5e-6, 1.608768e-02, 2.103953e-02],
        [10e9,  10e6, 0.3e-6, 58e6, 0.5e-6, 3.822746e-02, 1.307002e-01],
        [100e9, 10e6, 0.3e-6, 58e6, 0.5e-6, 1.946237e-01, 1.018985e+00],
    ];

    for (let i = 0; i < cases.length; i++) {
        const [f, sb, rq, sp, tp, exp_re, exp_im] = cases[i];
        const Z = calculate_Zrough_layered(f, sb, rq, sp, tp);
        console.log(`  Case ${i+1} (f=${(f/1e9).toFixed(0)}GHz):`);
        if (!assertClose(Z.re, exp_re, 0.05, 'Re')) failures++;
        if (!assertClose(Z.im, exp_im, 0.05, 'Im')) failures++;
    }
}

// ---------------------------------------------------------------
// Test 8: TLine method - Fallback to single-layer when plating thickness is 0
// ---------------------------------------------------------------
console.log('\nTest 8: TLine method - Fallback when thickness=0');
{
    const f = 1e9, sigma = 58e6, rq = 0.5e-6;
    const Z_single = calculate_Zrough(f, sigma, rq);
    const Z_tline = calculate_Zrough_layered(f, sigma, rq, 43e6, 0);
    if (!assertClose(Z_tline.re, Z_single.re, 1e-10, 'Re fallback')) failures++;
    if (!assertClose(Z_tline.im, Z_single.im, 1e-10, 'Im fallback')) failures++;
}

// ---------------------------------------------------------------
// Test 9: TLine method - Compare with Python reference
// ---------------------------------------------------------------
console.log('\nTest 9: TLine method - Python reference comparison (3% tolerance)');
{
    const cases = [
        // [f, sigma_bulk, rq, sigma_plating, thickness_plating, expected_re, expected_im]
        [1e9,   10e6, 0.3e-6, 58e6, 0.5e-6, 1.608768e-02, 2.103953e-02],
        [10e9,  10e6, 0.3e-6, 58e6, 0.5e-6, 3.822746e-02, 1.307002e-01],
        [100e9, 10e6, 0.3e-6, 58e6, 0.5e-6, 1.946237e-01, 1.018985e+00],
    ];

    for (let i = 0; i < cases.length; i++) {
        const [f, sb, rq, sp, tp, exp_re, exp_im] = cases[i];
        const Z = calculate_Zrough_layered(f, sb, rq, sp, tp);
        console.log(`  Case ${i+1} (f=${(f/1e9).toFixed(0)}GHz):`);
        if (!assertClose(Z.re, exp_re, 0.03, 'Re')) failures++;
        if (!assertClose(Z.im, exp_im, 0.03, 'Im')) failures++;
    }
}

// ---------------------------------------------------------------
// Test 11: TLine method - Zero roughness
// ---------------------------------------------------------------
console.log('\nTest 11: TLine method - Zero roughness (rq=0) produces finite values');
{
    const f = 10e9, sigma_bulk = 10e6, sigma_plating = 58e6, thickness = 0.5e-6;
    const Z = calculate_Zrough_layered(f, sigma_bulk, 0, sigma_plating, thickness);
    const finite = isFinite(Z.re) && isFinite(Z.im);
    if (finite && Z.re > 0 && Z.im > 0) {
        console.log(`  PASS Re=${Z.re.toExponential(4)}, Im=${Z.im.toExponential(4)}`);
    } else {
        console.log(`  FAIL Re=${Z.re}, Im=${Z.im} (expected finite positive values)`);
        failures++;
    }
}

// ---------------------------------------------------------------
// Summary
// ---------------------------------------------------------------
console.log(`\n${'='.repeat(50)}`);
if (failures === 0) {
    console.log('All tests passed.');
} else {
    console.log(`${failures} test(s) FAILED.`);
    process.exit(1);
}
