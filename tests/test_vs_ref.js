import { MicrostripSolver } from '../src/microstrip.js';
import { computeSParamsSingleEnded, computeSParamsDifferential } from '../src/sparameters.js';
import { Complex } from '../src/complex.js';
import { readFileSync } from 'fs';
import { test_s4p_generation_lossless, test_s4p_generation, test_s2p_generation, test_s2p_generation2 } from './test_snp_export.js';

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
        // RLGC matrix thresholds
        // Note: R has higher tolerance due to differences in conductor loss modeling
        // (surface roughness, skin effect) between solvers
        'R': 20.0,
        'L': 10.0,
        'G': 15.0,
        'C': 10.0,
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

        // Handle array values (for RLGC matrices)
        if (Array.isArray(ref_value)) {
            // Test each element of the array
            for (let i = 0; i < ref_value.length; i++) {
                const elem_ref = ref_value[i];
                const elem_solved = solved_value[i];

                // Calculate relative error
                let rel_error;
                if (elem_ref !== 0) {
                    rel_error = Math.abs((elem_solved - elem_ref) / elem_ref) * 100;
                } else {
                    rel_error = Math.abs(elem_solved) * 100;
                }

                // Check against threshold
                const threshold = error_thresholds[param] || 10.0;
                const passed = rel_error <= threshold;
                all_passed = all_passed && passed;

                const status = passed ? "✓ PASS" : "✗ FAIL";

                // Format values based on parameter type
                let solved_str, ref_str, param_label;
                const indices = ['11', '12', '21', '22'];
                param_label = `${param}[${indices[i]}]`;

                if (param === 'C') {
                    solved_str = `${(elem_solved * 1e12).toFixed(2)} pF/m`;
                    ref_str = `${(elem_ref * 1e12).toFixed(2)} pF/m`;
                } else if (param === 'L') {
                    solved_str = `${(elem_solved * 1e9).toFixed(2)} nH/m`;
                    ref_str = `${(elem_ref * 1e9).toFixed(2)} nH/m`;
                } else if (param === 'G') {
                    solved_str = `${(elem_solved * 1e3).toFixed(2)} mS/m`;
                    ref_str = `${(elem_ref * 1e3).toFixed(2)} mS/m`;
                } else if (param === 'R') {
                    solved_str = `${elem_solved.toFixed(2)} Ω/m`;
                    ref_str = `${elem_ref.toFixed(2)} Ω/m`;
                } else {
                    solved_str = `${elem_solved.toFixed(3)}`;
                    ref_str = `${elem_ref.toFixed(3)}`;
                }

                console.log(`${param_label.padEnd(20)} ${solved_str.padEnd(15)} ${ref_str.padEnd(15)} ${rel_error.toFixed(2).padEnd(12)} ${status.padEnd(10)}`);

                test_results[param_label] = {
                    'solved': elem_solved,
                    'reference': elem_ref,
                    'error': rel_error,
                    'passed': passed
                };
            }
            continue;
        }

        // Handle scalar values
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

async function solve_microstrip_20ghz() {
    const solver = new MicrostripSolver({
        substrate_height: 0.508e-3,
        trace_width: 1.1e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 3.48,
        tan_delta: 0.0037,
        sigma_cond: 5.8e7,
        freq: 20e9,
        rq: 0,
        nx: 10,
        ny: 10,
        use_sm: false,
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

    const reference = {
        "Z0": 50.5,
        "eps_eff": 2.7078,
        "loss": 12.856
    };

    // Test against reference
    test_microstrip_solution(solver_results, reference, "Microstrip 50Ω 20 GHz");

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
        trace_thickness: -15e-6,
        gnd_thickness: 15e-6,
        epsilon_r: 3.48,
        epsilon_r_top: 3.3,
        tan_delta: 0.004,
        tan_delta_top: 0.004,
        enclosure_height: 177e-6,
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

async function solve_differential_stripline_rlgc() {
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
        enclosure_width: 3e-3,
        sigma_cond: 5.8e7,
        freq: 67e9,
        nx: 10,
        ny: 10,
        trace_spacing: 0.1e-3,  // This enables differential mode
        boundaries: ["gnd", "gnd", "gnd", "gnd"]
    });
    solver.use_causal_materials = false;

    const results = await solver.solve_adaptive({energy_tol: 0.001});
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

    // Add 2x2 RLGC matrix values (flattened for comparison)
    if (results.RLGC_matrix) {
        const m = results.RLGC_matrix;
        solver_results.R = [m.R[0][0], m.R[0][1], m.R[1][0], m.R[1][1]];
        solver_results.L = [m.L[0][0], m.L[0][1], m.L[1][0], m.L[1][1]];
        solver_results.G = [m.G[0][0], m.G[0][1], m.G[1][0], m.G[1][1]];
        solver_results.C = [m.C[0][0], m.C[0][1], m.C[1][0], m.C[1][1]];
    }

    const reference = {
        'Z_odd': 37.6,
        'Z_even': 61.36,
        'eps_eff_even': 4.162,
        'eps_eff_odd': 4.195,
        'alpha_total_odd': 278.5,
        'alpha_total_even': 266.2,
        // 2x2 RLGC matrix (left and right traces)
        'R': [297.62, 13.56, 13.56, 297.62],
        'L': [3.32e-7, 8e-8, 8e-8, 3.2e-7],
        'C': [1.46e-10, -3.53e-11, -3.53e-11, 1.46e-10],
        'G': [1.23, -0.297, -0.297, 1.23]
    };

    // Test against reference
    test_differential_solution(solver_results, reference, "Differential Stripline");

    return results;
}

// Run tests
async function runTests() {
    await solve_microstrip();
    await solve_microstrip_1khz();
    await solve_microstrip_embed();
    await solve_microstrip_cut();
    await solve_microstrip_20ghz();
    await solve_stripline();
    await solve_rough_stripline();
    await solve_differential_stripline();
    await solve_differential_stripline_rlgc();
    await solve_differential_microstrip();
    await test_s2p_generation2();
    await test_s2p_generation();
    await test_s4p_generation_lossless();
    await test_s4p_generation();
}

runTests();
