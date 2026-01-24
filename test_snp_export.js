import { MicrostripSolver } from './microstrip.js';
import { computeSParamsSingleEnded, computeSParamsDifferential } from './sparameters.js';
import { Complex } from './complex.js';
import { readFileSync } from 'fs';

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

/**
 * Test S2P generation against reference file using absolute difference
 */
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
        enclosure_width: 3e-3,
        freq: reference[0].freq,
        nx: 10,
        ny: 10,
        trace_spacing: 0.1e-3,
        boundaries: ["gnd", "gnd", "gnd", "gnd"]
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


async function test_s4p_generation_lossless() {
    console.log(`\n${'='.repeat(80)}`);
    console.log('S4P GENERATION TEST (vs HFSS reference) - Absolute Difference');
    console.log(`${'='.repeat(80)}`);

    // Parse reference file
    const reference = parseS4P_MA('./stripline_2d_diff_lossless_fsweep.s4p');
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
        tan_delta: 0,
        tan_delta_top: 0,
        enclosure_height: 0.2e-3 + 35e-6,
        enclosure_width: 3e-3,
        freq: reference[0].freq,
        nx: 10,
        ny: 10,
        trace_spacing: 0.1e-3,
        boundaries: ["gnd", "gnd", "gnd", "gnd"]
    });

    // Build mesh once
    solver.ensure_mesh();

    // Solve for each frequency and compute S-parameters
    const length = 50e-3; // 50 mm line length (same as HFSS reference)
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
    await test_s4p_generation_lossless();
    await test_s4p_generation();
    await test_s2p_generation();
}

runTests();

export { test_s4p_generation_lossless, test_s4p_generation, test_s2p_generation };
