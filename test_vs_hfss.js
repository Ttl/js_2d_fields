import { MicrostripSolver } from './microstrip.js';
import { computeSParamsSingleEnded, computeSParamsDifferential } from './sparameters.js';
import { Complex } from './complex.js';
import { readFileSync } from 'fs';

/**
 * Test microstrip solver results against reference values.
 *
 * @param {Object} solver_results - Object with keys 'Z0', 'eps_eff', 'alpha_diel', 'alpha_cond', 'C', 'R', 'L', 'G'
 * @param {Object} reference - Object with reference values
 * @param {string} test_name - Name of the test for printing
 * @returns {Object} Test results with pass/fail status
 */
function test_microstrip_solution(solver_results, reference, test_name = "Microstrip") {
    // Global error thresholds (relative error in %)
    const MAX_Z0_ERROR = 5.0;
    const MAX_DIEL_LOSS_ERROR = 10.0;
    const MAX_COND_LOSS_ERROR = 50.0;
    const MAX_C_ERROR = 10.0;
    const MAX_R_ERROR = 15.0;
    const MAX_L_ERROR = 10.0;
    const MAX_G_ERROR = 15.0;
    const MAX_EPS_EFF_ERROR = 5.0;

    // Error thresholds mapping
    const error_thresholds = {
        'Z0': MAX_Z0_ERROR,
        'diel_loss': MAX_DIEL_LOSS_ERROR,
        'cond_loss': MAX_COND_LOSS_ERROR,
        'C': MAX_C_ERROR,
        'R': MAX_R_ERROR,
        'L': MAX_L_ERROR,
        'G': MAX_G_ERROR,
        'eps_eff': MAX_EPS_EFF_ERROR
    };

    console.log(`\n${'='.repeat(80)}`);
    console.log(`${test_name.toUpperCase()} VALIDATION TEST`);
    console.log(`${'='.repeat(80)}`);
    console.log(`${'Parameter'.padEnd(15)} ${'Solved'.padEnd(15)} ${'Reference'.padEnd(15)} ${'Error (%)'.padEnd(12)} ${'Status'.padEnd(10)}`);
    console.log(`${'-'.repeat(80)}`);

    let all_passed = true;
    const test_results = {};

    for (const [param, ref_value] of Object.entries(reference)) {
        if (!(param in solver_results)) {
            continue;
        }

        const solved_value = solver_results[param];

        // Calculate relative error
        let rel_error;
        if (ref_value !== 0) {
            rel_error = Math.abs((solved_value - ref_value) / ref_value) * 100;
        } else {
            rel_error = Math.abs(solved_value) * 100;
        }

        // Check against threshold
        const threshold = error_thresholds[param] || 10.0; // default 10%
        const passed = rel_error <= threshold;
        all_passed = all_passed && passed;

        const status = passed ? "✓ PASS" : "✗ FAIL";

        // Format values based on magnitude
        let solved_str, ref_str;
        if (param === 'C') {
            solved_str = `${(solved_value * 1e12).toFixed(2)} pF`;
            ref_str = `${(ref_value * 1e12).toFixed(2)} pF`;
        } else if (param === 'L') {
            solved_str = `${(solved_value * 1e9).toFixed(2)} nH`;
            ref_str = `${(ref_value * 1e9).toFixed(2)} nH`;
        } else if (param === 'G') {
            solved_str = `${(solved_value * 1e3).toFixed(2)} mS`;
            ref_str = `${(ref_value * 1e3).toFixed(2)} mS`;
        } else if (param === 'R') {
            solved_str = `${solved_value.toFixed(2)} Ω`;
            ref_str = `${ref_value.toFixed(2)} Ω`;
        } else if (param.includes('loss')) {
            solved_str = `${solved_value.toFixed(4)} dB/m`;
            ref_str = `${ref_value.toFixed(4)} dB/m`;
        } else {
            solved_str = `${solved_value.toFixed(3)}`;
            ref_str = `${ref_value.toFixed(3)}`;
        }

        console.log(`${param.padEnd(15)} ${solved_str.padEnd(15)} ${ref_str.padEnd(15)} ${rel_error.toFixed(2).padEnd(12)} ${status.padEnd(10)}`);

        test_results[param] = {
            'solved': solved_value,
            'reference': ref_value,
            'error': rel_error,
            'passed': passed
        };
    }

    console.log(`${'-'.repeat(80)}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error(`${test_name} validation failed - see errors above`);
    }

    return test_results;
}

function test_differential_solution(solver_results, reference, test_name = "Differential Microstrip") {
    // Global error thresholds (relative error in %)
    const MAX_Z_DIFF_ERROR = 6.0;
    const MAX_Z_COMMON_ERROR = 6.0;
    const MAX_Z_ODD_ERROR = 6.0;
    const MAX_Z_EVEN_ERROR = 6.0;
    const MAX_EPS_EFF_ERROR = 5.0;
    const MAX_LOSS_ERROR = 50.0;
    const MAX_C_ERROR = 10.0;

    // Error thresholds mapping
    const error_thresholds = {
        'Z_diff': MAX_Z_DIFF_ERROR,
        'Z_common': MAX_Z_COMMON_ERROR,
        'Z_odd': MAX_Z_ODD_ERROR,
        'Z_even': MAX_Z_EVEN_ERROR,
        'eps_eff_odd': MAX_EPS_EFF_ERROR,
        'eps_eff_even': MAX_EPS_EFF_ERROR,
        'C_odd': MAX_C_ERROR,
        'C_even': MAX_C_ERROR,
        'alpha_c_odd': MAX_LOSS_ERROR,
        'alpha_c_even': MAX_LOSS_ERROR,
        'alpha_d_odd': MAX_LOSS_ERROR,
        'alpha_d_even': MAX_LOSS_ERROR,
        'alpha_total_odd': MAX_LOSS_ERROR,
        'alpha_total_even': MAX_LOSS_ERROR,
    };

    console.log(`\n${'='.repeat(80)}`);
    console.log(`${test_name.toUpperCase()} VALIDATION TEST`);
    console.log(`${'='.repeat(80)}`);
    console.log(`${'Parameter'.padEnd(20)} ${'Solved'.padEnd(15)} ${'Reference'.padEnd(15)} ${'Error (%)'.padEnd(12)} ${'Status'.padEnd(10)}`);
    console.log(`${'-'.repeat(80)}`);

    let all_passed = true;
    const test_results = {};

    for (const [param, ref_value] of Object.entries(reference)) {
        if (!(param in solver_results)) {
            continue;
        }

        const solved_value = solver_results[param];

        // Calculate relative error
        let rel_error;
        if (ref_value !== 0) {
            rel_error = Math.abs((solved_value - ref_value) / ref_value) * 100;
        } else {
            rel_error = Math.abs(solved_value) * 100;
        }

        // Check against threshold
        const threshold = error_thresholds[param] || 10.0;
        const passed = rel_error <= threshold;
        all_passed = all_passed && passed;

        const status = passed ? "✓ PASS" : "✗ FAIL";

        // Format values based on parameter type
        let solved_str, ref_str;
        if (param.startsWith('C_')) {
            solved_str = `${(solved_value * 1e12).toFixed(2)} pF`;
            ref_str = `${(ref_value * 1e12).toFixed(2)} pF`;
        } else if (param.startsWith('alpha_')) {
            solved_str = `${solved_value.toFixed(4)} dB/m`;
            ref_str = `${ref_value.toFixed(4)} dB/m`;
        } else if (param.startsWith('Z_')) {
            solved_str = `${solved_value.toFixed(2)} Ω`;
            ref_str = `${ref_value.toFixed(2)} Ω`;
        } else {
            solved_str = `${solved_value.toFixed(3)}`;
            ref_str = `${ref_value.toFixed(3)}`;
        }

        console.log(`${param.padEnd(20)} ${solved_str.padEnd(15)} ${ref_str.padEnd(15)} ${rel_error.toFixed(2).padEnd(12)} ${status.padEnd(10)}`);

        test_results[param] = {
            'solved': solved_value,
            'reference': ref_value,
            'error': rel_error,
            'passed': passed
        };
    }

    console.log(`${'-'.repeat(80)}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error(`${test_name} validation failed - see errors above`);
    }

    return test_results;
}

async function solve_microstrip() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 10,
        ny: 10,
        use_sm: false,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const mode = results.modes[0];

    // Prepare results dictionary matching Python format
    const solver_results = {
        'Z0': mode.Z0,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    // Reference values from HFSS
    const reference = {
        "Z0": 49.8,
        "diel_loss": 2.99,
        "cond_loss": 0.285,
        "C": 123e-12,
        "R": 3.26,
        "G": 13.86e-3,
        "L": 307e-9
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Microstrip 50Ω");

    return solver_results;
}

async function solve_microstrip_1khz() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0,
        sigma_cond: 1e7,
        enclosure_width: 50e-3,
        freq: 1e3,
        nx: 10,
        ny: 10,
        use_sm: false,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const mode = results.modes[0];

    // Prepare results dictionary matching Python format
    const solver_results = {
        'RZc': mode.Zc.re,
        'IZc': mode.Zc.im,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    // Reference values from HFSS
    const reference = {
        "RZc": 866,
        "IZc": -864,
        "cond_loss": 0.0058,
        "C": 123.25e-12,
        "R": 1.01,
        "G": 0,
        // L can't be solved correctly without solving for magnetic field due
        // to current spreading in the ground plane.
        // "L": 523e-9
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Microstrip 1 kHz");

    return solver_results;
}

async function solve_microstrip_embed() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 10,
        ny: 10,
        use_sm: false,
        top_diel_h: 0.2e-3,
        top_diel_er: 4.5,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const mode = results.modes[0];

    const solver_results = {
        'Z0': mode.Z0,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'loss': mode.alpha_total,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    // Reference values from HFSS
    const reference = {
        "Z0": 48.15,
        "eps_eff": 3.621,
        "loss": 3.48,
        "C": 131.8e-12
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Embedded microstrip");

    return solver_results;
}

async function solve_microstrip_cut() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: 1e9,
        nx: 10,
        ny: 10,
        use_sm: false,
        gnd_cut_width: 3e-3,
        gnd_cut_sub_h: 1e-3,
        boundaries: ["open", "open", "open", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const mode = results.modes[0];

    const solver_results = {
        'Z0': mode.Z0,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'loss': mode.alpha_total,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    // Reference values from HFSS
    const reference = {
        "Z0": 55.84,
        "eps_eff": 3.28,
        "loss": 3.19,
        "C": 108.25e-12
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Microstrip with ground cut");

    return solver_results;
}

async function solve_differential_microstrip() {
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.5,
        freq: 1e9,
        nx: 10,
        ny: 10,
        trace_spacing: 1e-3  // This enables differential mode
    });

    const results = await solver.solve_adaptive({energy_tol: 0.01});
    const odd = results.modes.find(m => m.mode === 'odd');
    const even = results.modes.find(m => m.mode === 'even');

    // Map to flat structure for test function
    const solver_results = {
        'Z_odd': odd.Z0,
        'Z_even': even.Z0,
        'Z_diff': results.Z_diff,
        'Z_common': results.Z_common,
        'eps_eff_odd': odd.eps_eff,
        'eps_eff_even': even.eps_eff,
        'alpha_c_odd': odd.alpha_c,
        'alpha_c_even': even.alpha_c,
        'alpha_d_odd': odd.alpha_d,
        'alpha_d_even': even.alpha_d,
        'alpha_total_odd': odd.alpha_total,
        'alpha_total_even': even.alpha_total
    };

    const reference = {
        'Z_odd': 40.23,
        'Z_even': 57.65,
        'eps_eff_even': 3.65,
        'eps_eff_odd': 2.98,
        'alpha_c_odd': 0.363,
        'alpha_c_even': 0.269,
        'alpha_d_odd': 2.67,
        'alpha_d_even': 3.21
    };

    // Test against reference
    test_differential_solution(solver_results, reference, "Differential Microstrip");

    return results;
}

async function solve_stripline() {
    const solver = new MicrostripSolver({
        substrate_height: 0.2e-3,
        trace_width: 0.15e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.1,
        epsilon_r_top: 4.1,
        tan_delta: 0.02,
        tan_delta_top: 0.02,
        enclosure_height: 0.2e-3 + 35e-6,
        freq: 1e9,
        nx: 10,
        ny: 10,
        boundaries: ["open", "open", "gnd", "gnd"]
    });

    const results = await solver.solve_adaptive({energy_tol: 0.01});
    const mode = results.modes[0];

    const solver_results = {
        'Z0': mode.Z0,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'loss': mode.alpha_total,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    // Reference values from HFSS
    const reference = {
        "Z0": 50.2,
        "diel_loss": 3.685,
        "cond_loss": 3.1
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Stripline");

    return results;
}

async function solve_rough_stripline() {
    const solver = new MicrostripSolver({
        substrate_height: 177e-6,
        trace_width: 160e-6,
        trace_thickness: 15e-6,
        gnd_thickness: 15e-6,
        epsilon_r: 3.54,
        epsilon_r_top: 3.48,
        tan_delta: 0.004,
        tan_delta_top: 0.004,
        enclosure_height: 162e-6 + 15e-6,
        freq: 40e9,
        nx: 30,
        ny: 30,
        rq: 0.6e-6,
        boundaries: ["open", "open", "gnd", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const mode = results.modes[0];

    const solver_results = {
        'Z0': mode.Z0,
        'eps_eff': mode.eps_eff,
        'diel_loss': mode.alpha_d,
        'cond_loss': mode.alpha_c,
        'loss': mode.alpha_total,
        'C': mode.RLGC.C,
        'R': mode.RLGC.R,
        'L': mode.RLGC.L,
        'G': mode.RLGC.G
    };

    const phase_delay = 6.3e-9;
    const eps_eff = Math.pow(3e8 * phase_delay, 2);

    // Reference values from Gradient model paper
    // G. Gold and K. Helmreich, "A Physical Surface Roughness Model and Its
    // Applications," in IEEE Transactions on Microwave Theory and Techniques,
    // vol. 65, no. 10, pp. 3720-3732, Oct. 2017.
    const reference = {
        "eps_eff": eps_eff,
        "loss": 80,
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Rough Stripline");

    return results;
}

async function solve_differential_stripline() {
    const solver = new MicrostripSolver({
        substrate_height: 0.2e-3,
        trace_width: 0.15e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.1,
        epsilon_r_top: 4.1,
        tan_delta: 0.02,
        tan_delta_top: 0.02,
        enclosure_height: 0.2e-3 + 35e-6,
        freq: 1e9,
        nx: 10,
        ny: 10,
        trace_spacing: 0.1e-3,  // This enables differential mode
        boundaries: ["open", "open", "gnd", "gnd"]
    });

    const results = await solver.solve_adaptive();
    const odd = results.modes.find(m => m.mode === 'odd');
    const even = results.modes.find(m => m.mode === 'even');

    // Map to flat structure for test function
    const solver_results = {
        'Z_odd': odd.Z0,
        'Z_even': even.Z0,
        'Z_diff': results.Z_diff,
        'Z_common': results.Z_common,
        'eps_eff_odd': odd.eps_eff,
        'eps_eff_even': even.eps_eff,
        'alpha_c_odd': odd.alpha_c,
        'alpha_c_even': even.alpha_c,
        'alpha_d_odd': odd.alpha_d,
        'alpha_d_even': even.alpha_d,
        'alpha_total_odd': odd.alpha_total,
        'alpha_total_even': even.alpha_total
    };

    const reference = {
        'Z_odd': 37.6,
        'Z_even': 61.36,
        'eps_eff_even': 4.162,
        'eps_eff_odd': 4.195,
        'alpha_total_odd': 7.93,
        'alpha_total_even': 6.47,
    };

    // Test against reference
    test_differential_solution(solver_results, reference, "Differential Stripline");

    return results;
}

/**
 * Convert dB and angle (degrees) to Complex number
 * @param {number} dB - Magnitude in dB
 * @param {number} angDeg - Angle in degrees
 * @returns {Complex}
 */
function dbAngToComplex(dB, angDeg) {
    const mag = Math.pow(10, dB / 20);
    const angRad = angDeg * Math.PI / 180;
    return new Complex(mag * Math.cos(angRad), mag * Math.sin(angRad));
}

/**
 * Convert magnitude and angle (degrees) to Complex number
 * @param {number} mag - Linear magnitude
 * @param {number} angDeg - Angle in degrees
 * @returns {Complex}
 */
function magAngToComplex(mag, angDeg) {
    const angRad = angDeg * Math.PI / 180;
    return new Complex(mag * Math.cos(angRad), mag * Math.sin(angRad));
}

/**
 * Calculate absolute difference between two Complex numbers
 * @param {Complex} a
 * @param {Complex} b
 * @returns {number}
 */
function complexAbsDiff(a, b) {
    return a.sub(b).abs();
}

/**
 * Parse Touchstone S2P file in DB format (dB + angle)
 * @param {string} filename - Path to the S2P file
 * @returns {Array<{freq: number, S11: Complex, S21: Complex, S12: Complex, S22: Complex}>}
 */
function parseS2P_DB(filename) {
    const content = readFileSync(filename, 'utf-8');
    const lines = content.split('\n');
    const data = [];

    let freqUnit = 1e9; // Default GHz

    for (const line of lines) {
        const trimmed = line.trim();
        if (!trimmed || trimmed.startsWith('!')) continue;

        // Parse option line
        if (trimmed.startsWith('#')) {
            const parts = trimmed.toUpperCase().split(/\s+/);
            if (parts.includes('HZ')) freqUnit = 1;
            else if (parts.includes('KHZ')) freqUnit = 1e3;
            else if (parts.includes('MHZ')) freqUnit = 1e6;
            else if (parts.includes('GHZ')) freqUnit = 1e9;
            continue;
        }

        // Parse data line
        const parts = trimmed.split(/\s+/).map(parseFloat);
        if (parts.length >= 9 && !isNaN(parts[0])) {
            // Format: freq S11_dB S11_ang S21_dB S21_ang S12_dB S12_ang S22_dB S22_ang
            data.push({
                freq: parts[0] * freqUnit,
                S11: dbAngToComplex(parts[1], parts[2]),
                S21: dbAngToComplex(parts[3], parts[4]),
                S12: dbAngToComplex(parts[5], parts[6]),
                S22: dbAngToComplex(parts[7], parts[8])
            });
        }
    }
    return data;
}

/**
 * Parse Touchstone S4P file in MA format (magnitude + angle)
 * @param {string} filename - Path to the S4P file
 * @returns {Array<{freq: number, S: Array<Array<Complex>>}>}
 */
function parseS4P_MA(filename) {
    const content = readFileSync(filename, 'utf-8');
    const lines = content.split('\n').filter(l => {
        const t = l.trim();
        return t && !t.startsWith('!') && !t.startsWith('#');
    });

    const data = [];
    let i = 0;

    while (i < lines.length) {
        const parts1 = lines[i].trim().split(/\s+/).map(parseFloat);
        if (parts1.length < 9 || isNaN(parts1[0])) {
            i++;
            continue;
        }

        const freq = parts1[0] * 1e9; // GHz to Hz
        const S = [[null, null, null, null],
                   [null, null, null, null],
                   [null, null, null, null],
                   [null, null, null, null]];

        // Row 1: S11 S12 S13 S14
        for (let col = 0; col < 4; col++) {
            S[0][col] = magAngToComplex(parts1[1 + col*2], parts1[2 + col*2]);
        }

        // Rows 2-4
        for (let row = 1; row < 4; row++) {
            i++;
            if (i >= lines.length) break;
            const parts = lines[i].trim().split(/\s+/).map(parseFloat);
            for (let col = 0; col < 4; col++) {
                S[row][col] = magAngToComplex(parts[col*2], parts[col*2 + 1]);
            }
        }

        data.push({ freq, S });
        i++;
    }
    return data;
}

/**
 * Convert magnitude (linear) to dB
 */
function magTodB(mag) {
    return 20 * Math.log10(Math.max(mag, 1e-15));
}

async function test_s2p_generation() {
    console.log(`\n${'='.repeat(80)}`);
    console.log('S2P GENERATION TEST (vs HFSS reference) - Absolute Difference');
    console.log(`${'='.repeat(80)}`);

    // Parse reference file
    const reference = parseS2P_DB('./ms_2d_fr4_sweep.s2p');
    console.log(`Loaded ${reference.length} frequency points from reference`);

    // Setup solver with same parameters as reference
    // t=35um, w=3mm, z_sub=1.6mm
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: 5.8e7,
        freq: reference[0].freq,
        nx: 10,
        ny: 10,
        boundaries: ["open", "open", "open", "gnd"]
    });

    // Build mesh once
    solver.ensure_mesh();

    // Solve for each frequency and compute S-parameters
    const length = 1.0; // 1 m line length (same as HFSS reference)
    const Z_ref = 50;

    let all_passed = true;
    const maxErrors = { S11: 0, S21: 0, S12: 0, S22: 0 };

    // Absolute tolerance for S-parameter comparison
    // S11 (return loss) is sensitive to small Z0 differences
    // S21/S12 (insertion loss) differences accumulate over the 1m line length,
    // so small per-unit-length differences lead to larger S-param differences
    const tolerances = { S11: 0.1, S21: 0.25, S12: 0.25, S22: 0.1 };
    const init_result = await solver.solve_adaptive({energy_tol: 0.01});

    console.log(`\n${'Freq'.padEnd(8)} ${'|ΔS11|'.padEnd(10)} ${'|ΔS21|'.padEnd(10)} ${'|ΔS12|'.padEnd(10)} ${'|ΔS22|'.padEnd(10)} Status`);
    console.log(`${'-'.repeat(70)}`);

    for (const refPoint of reference) {
        // Update frequency and solve
        const freq = refPoint.freq;
        solver.omega = 2 * Math.PI * solver.freq;
        solver.delta_s = Math.sqrt(2 / (solver.omega * 4 * Math.PI * 1e-7 * solver.sigma_cond));

        const result = solver.computeAtFrequency(freq, init_result);
        const mode = result.modes[0];

        // Compute S-parameters
        const sp = computeSParamsSingleEnded(refPoint.freq, mode.RLGC, length, Z_ref);

        // Compute absolute differences for all S-parameters
        const errors = {
            S11: complexAbsDiff(sp.S11, refPoint.S11),
            S21: complexAbsDiff(sp.S21, refPoint.S21),
            S12: complexAbsDiff(sp.S12, refPoint.S12),
            S22: complexAbsDiff(sp.S22, refPoint.S22)
        };

        // Track max errors
        for (const key of Object.keys(errors)) {
            maxErrors[key] = Math.max(maxErrors[key], errors[key]);
        }

        // Check against tolerances
        const passed = errors.S11 < tolerances.S11 &&
                       errors.S21 < tolerances.S21 &&
                       errors.S12 < tolerances.S12 &&
                       errors.S22 < tolerances.S22;
        all_passed = all_passed && passed;

        const freqGHz = (refPoint.freq / 1e9).toFixed(2);
        const status = passed ? '✓' : '✗';
        console.log(`${freqGHz.padEnd(8)} ${errors.S11.toFixed(4).padEnd(10)} ${errors.S21.toFixed(4).padEnd(10)} ${errors.S12.toFixed(4).padEnd(10)} ${errors.S22.toFixed(4).padEnd(10)} ${status}`);
    }

    console.log(`${'-'.repeat(70)}`);
    console.log(`Max absolute errors: S11=${maxErrors.S11.toFixed(4)}, S21=${maxErrors.S21.toFixed(4)}, S12=${maxErrors.S12.toFixed(4)}, S22=${maxErrors.S22.toFixed(4)}`);
    console.log(`Tolerances:          S11=${tolerances.S11}, S21=${tolerances.S21}, S12=${tolerances.S12}, S22=${tolerances.S22}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error('S2P generation test failed - see errors above');
    }
}

/**
 * Test S4P generation against reference file using absolute difference
 */
async function test_s4p_generation() {
    console.log(`\n${'='.repeat(80)}`);
    console.log('S4P GENERATION TEST (vs HFSS reference) - Absolute Difference');
    console.log(`${'='.repeat(80)}`);

    // Parse reference file
    const reference = parseS4P_MA('./stripline_2d_diff_sweep.s4p');
    console.log(`Loaded ${reference.length} frequency points from reference`);

    // Setup solver with same parameters as reference
    // s=100um, t=35um, w=0.15mm, z_sub=200um
    const solver = new MicrostripSolver({
        substrate_height: 0.2e-3,
        trace_width: 0.15e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.1,
        epsilon_r_top: 4.1,
        tan_delta: 0.02,
        tan_delta_top: 0.02,
        enclosure_height: 0.2e-3 + 35e-6,
        freq: reference[0].freq,
        nx: 10,
        ny: 10,
        trace_spacing: 0.1e-3,
        boundaries: ["open", "open", "gnd", "gnd"]
    });

    // Build mesh once
    solver.ensure_mesh();

    // Solve for each frequency and compute S-parameters
    const length = 1.0; // 1 m line length (same as HFSS reference)
    const Z_ref = 50;

    let all_passed = true;

    // Track max errors for all 16 S-parameters
    const maxErrors = Array(4).fill(null).map(() => Array(4).fill(0));

    // Tolerance for absolute difference
    // Similar to S2P, differences accumulate over the 1m line length
    const tolerance = 0.25;
    const init_result = await solver.solve_adaptive({energy_tol: 0.01});

    console.log(`\n${'Freq'.padEnd(8)} ${'Max|ΔSij|'.padEnd(12)} ${'Worst Sij'.padEnd(10)} Status`);
    console.log(`${'-'.repeat(50)}`);

    for (const refPoint of reference) {
        // Update frequency and solve
        const freq = refPoint.freq;
        solver.omega = 2 * Math.PI * solver.freq;
        solver.delta_s = Math.sqrt(2 / (solver.omega * 4 * Math.PI * 1e-7 * solver.sigma_cond));

        const result = solver.computeAtFrequency(freq, init_result);
        const oddMode = result.modes.find(m => m.mode === 'odd');
        const evenMode = result.modes.find(m => m.mode === 'even');

        // Compute 4-port S-parameters
        const sp = computeSParamsDifferential(refPoint.freq, oddMode.RLGC, evenMode.RLGC, length, Z_ref);

        // Compare all 16 S-parameters
        let maxError = 0;
        let worstParam = 'S11';

        for (let row = 0; row < 4; row++) {
            for (let col = 0; col < 4; col++) {
                const error = complexAbsDiff(sp.S[row][col], refPoint.S[row][col]);
                maxErrors[row][col] = Math.max(maxErrors[row][col], error);

                if (error > maxError) {
                    maxError = error;
                    worstParam = `S${row + 1}${col + 1}`;
                }
            }
        }

        const passed = maxError < tolerance;
        all_passed = all_passed && passed;

        const freqGHz = (refPoint.freq / 1e9).toFixed(2);
        const status = passed ? '✓' : '✗';
        console.log(`${freqGHz.padEnd(8)} ${maxError.toFixed(4).padEnd(12)} ${worstParam.padEnd(10)} ${status}`);
    }

    console.log(`${'-'.repeat(50)}`);
    console.log(`\nMax absolute errors per S-parameter:`);
    console.log(`       Port1     Port2     Port3     Port4`);
    for (let row = 0; row < 4; row++) {
        const rowStr = maxErrors[row].map(e => e.toFixed(4).padStart(9)).join(' ');
        console.log(`Port${row + 1} ${rowStr}`);
    }
    console.log(`\nTolerance: ${tolerance}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error('S4P generation test failed - see errors above');
    }
}

// Run tests
async function runTests() {
    await solve_microstrip();
    await solve_microstrip_1khz();
    await solve_microstrip_embed();
    await solve_microstrip_cut();
    await solve_stripline();
    await solve_rough_stripline();
    await solve_differential_stripline();
    await solve_differential_microstrip();
    await test_s2p_generation();
    await test_s4p_generation();
}

runTests();
