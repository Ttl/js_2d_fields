import { MicrostripSolver } from '../src/microstrip.js';
import { computeSParamsSingleEnded, computeSParamsDifferential } from '../src/sparameters.js';
import { Complex } from '../src/complex.js';
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
    let formatVerified = false;

    for (const line of lines) {
        const trimmed = line.trim();
        if (!trimmed || trimmed.startsWith('!')) continue;

        // Parse option line
        if (trimmed.startsWith('#')) {
            const parts = trimmed.toUpperCase().split(/\s+/);

            // Verify format is DB
            if (!parts.includes('DB')) {
                throw new Error(`Expected DB format but found: ${trimmed}`);
            }
            formatVerified = true;

            if (parts.includes('HZ')) freqUnit = 1;
            else if (parts.includes('KHZ')) freqUnit = 1e3;
            else if (parts.includes('MHZ')) freqUnit = 1e6;
            else if (parts.includes('GHZ')) freqUnit = 1e9;
            continue;
        }

        // Parse data line
        const parts = trimmed.split(/\s+/).map(parseFloat);
        if (parts.length >= 9 && !isNaN(parts[0])) {
            if (!formatVerified) {
                throw new Error('Missing format specification line (# ...)');
            }
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
 * Parse Touchstone S2P file in MA format (magnitude + angle)
 * @param {string} filename - Path to the S2P file
 * @returns {Array<{freq: number, S11: Complex, S21: Complex, S12: Complex, S22: Complex}>}
 */
function parseS2P_MA(filename) {
    const content = readFileSync(filename, 'utf-8');
    const lines = content.split('\n');
    const data = [];

    let freqUnit = 1e9; // Default GHz
    let formatVerified = false;

    for (const line of lines) {
        const trimmed = line.trim();
        if (!trimmed || trimmed.startsWith('!')) continue;

        // Parse option line
        if (trimmed.startsWith('#')) {
            const parts = trimmed.toUpperCase().split(/\s+/);

            // Verify format is MA
            if (!parts.includes('MA')) {
                throw new Error(`Expected MA format but found: ${trimmed}`);
            }
            formatVerified = true;

            if (parts.includes('HZ')) freqUnit = 1;
            else if (parts.includes('KHZ')) freqUnit = 1e3;
            else if (parts.includes('MHZ')) freqUnit = 1e6;
            else if (parts.includes('GHZ')) freqUnit = 1e9;
            continue;
        }

        // Parse data line
        const parts = trimmed.split(/\s+/).map(parseFloat);
        if (parts.length >= 9 && !isNaN(parts[0])) {
            if (!formatVerified) {
                throw new Error('Missing format specification line (# ...)');
            }
            // Format: freq S11_mag S11_ang S21_mag S21_ang S12_mag S12_ang S22_mag S22_ang
            data.push({
                freq: parts[0] * freqUnit,
                S11: magAngToComplex(parts[1], parts[2]),
                S21: magAngToComplex(parts[3], parts[4]),
                S12: magAngToComplex(parts[5], parts[6]),
                S22: magAngToComplex(parts[7], parts[8])
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
    const allLines = content.split('\n');

    let freqUnit = 1e9; // Default GHz
    let formatVerified = false;

    // Check format in header
    for (const line of allLines) {
        const trimmed = line.trim();
        if (trimmed.startsWith('#')) {
            const parts = trimmed.toUpperCase().split(/\s+/);

            // Verify format is MA
            if (!parts.includes('MA')) {
                throw new Error(`Expected MA format but found: ${trimmed}`);
            }
            formatVerified = true;

            if (parts.includes('HZ')) freqUnit = 1;
            else if (parts.includes('KHZ')) freqUnit = 1e3;
            else if (parts.includes('MHZ')) freqUnit = 1e6;
            else if (parts.includes('GHZ')) freqUnit = 1e9;
            break;
        }
    }

    if (!formatVerified) {
        throw new Error('Missing format specification line (# ...)');
    }

    // Filter out comments and option lines for data parsing
    const lines = allLines.filter(l => {
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

        const freq = parts1[0] * freqUnit;
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


async function test_s2p_generation() {
    console.log(`\n${'='.repeat(80)}`);
    console.log('S2P GENERATION TEST');
    console.log(`${'='.repeat(80)}`);

    const reference = parseS2P_MA('./data/ms_2d_fr4_causal_5mm.s2p');
    console.log(`Loaded ${reference.length} frequency points from reference`);

    // Setup solver with same parameters as reference
    const solver = new MicrostripSolver({
        substrate_height: 1.6e-3,
        trace_width: 1e-3,
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
    solver.use_causal_materials = true;

    const sweepResults = await solver.solve_sweep({
        frequencies: reference.map(r => r.freq),
        energy_tol: 0.01
    });

    console.log(`Mesh: ${sweepResults.mesh.nx}x${sweepResults.mesh.ny}`);
    console.log(`Frequencies: ${sweepResults.frequencies.length} points`);

    const length = 5e-3;
    const Z_ref = 50;

    let all_passed = true;
    const maxErrors = { S11: 0, S21: 0, S12: 0, S22: 0 };
    const tolerances = { S11: 0.1, S21: 0.25, S12: 0.25, S22: 0.1 };

    console.log(`\n${'Freq'.padEnd(8)} ${'|ΔS11|'.padEnd(10)} ${'|ΔS21|'.padEnd(10)} ${'|ΔS12|'.padEnd(10)} ${'|ΔS22|'.padEnd(10)} Status`);
    console.log(`${'-'.repeat(70)}`);

    for (let i = 0; i < reference.length; i++) {
        const refPoint = reference[i];
        const freq = sweepResults.frequencies[i];
        const mode = sweepResults.modes[0];

        // Build RLGC object from arrays
        const RLGC = {
            R: mode.RLGC.R[i],
            L: mode.RLGC.L[i],
            G: mode.RLGC.G[i],
            C: mode.RLGC.C[i]
        };

        // Compute S-parameters
        const sp = computeSParamsSingleEnded(freq, RLGC, length, Z_ref);

        // Compute absolute differences
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

        const passed = errors.S11 < tolerances.S11 &&
                       errors.S21 < tolerances.S21 &&
                       errors.S12 < tolerances.S12 &&
                       errors.S22 < tolerances.S22;
        all_passed = all_passed && passed;

        const freqGHz = (freq / 1e9).toFixed(2);
        const status = passed ? '✓' : '✗';
        console.log(`${freqGHz.padEnd(8)} ${errors.S11.toFixed(4).padEnd(10)} ${errors.S21.toFixed(4).padEnd(10)} ${errors.S12.toFixed(4).padEnd(10)} ${errors.S22.toFixed(4).padEnd(10)} ${status}`);
    }

    console.log(`${'-'.repeat(70)}`);
    console.log(`Max absolute errors: S11=${maxErrors.S11.toFixed(4)}, S21=${maxErrors.S21.toFixed(4)}, S12=${maxErrors.S12.toFixed(4)}, S22=${maxErrors.S22.toFixed(4)}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error('S2P solve_sweep test failed');
    }
}

async function test_s2p_generation2() {
    console.log(`\n${'='.repeat(80)}`);
    console.log('S2P GENERATION TEST 2');
    console.log(`${'='.repeat(80)}`);

    const reference = parseS2P_MA('./data/stripline_2d_causal.s2p');
    console.log(`Loaded ${reference.length} frequency points from reference`);

    // Setup solver with same parameters as reference
    const solver = new MicrostripSolver({
        substrate_height: 0.2e-3,
        trace_width: 0.3e-3,
        trace_thickness: 35e-6,
        gnd_thickness: 16e-6,
        epsilon_r: 4.1,
        epsilon_r_top: 4.1,
        tan_delta: 0.02,
        tan_delta_top: 0.02,
        enclosure_height: 0.2e-3 + 35e-6,
        enclosure_width: 3e-3,
        freq: 67e9,
        nx: 10,
        ny: 10,
        boundaries: ["gnd", "gnd", "gnd", "gnd"]
    });

    solver.use_causal_materials = true;

    const sweepResults = await solver.solve_sweep({
        frequencies: reference.map(r => r.freq),
        energy_tol: 0.01
    });

    console.log(`Mesh: ${sweepResults.mesh.nx}x${sweepResults.mesh.ny}`);
    console.log(`Frequencies: ${sweepResults.frequencies.length} points`);

    const length = 50e-3;
    const Z_ref = 50;

    let all_passed = true;
    const maxErrors = { S11: 0, S21: 0, S12: 0, S22: 0 };
    const tolerances = { S11: 0.1, S21: 0.25, S12: 0.25, S22: 0.1 };

    console.log(`\n${'Freq'.padEnd(8)} ${'|ΔS11|'.padEnd(10)} ${'|ΔS21|'.padEnd(10)} ${'|ΔS12|'.padEnd(10)} ${'|ΔS22|'.padEnd(10)} Status`);
    console.log(`${'-'.repeat(70)}`);

    for (let i = 0; i < reference.length; i++) {
        const refPoint = reference[i];
        const freq = sweepResults.frequencies[i];
        const mode = sweepResults.modes[0];

        // Build RLGC object from arrays
        const RLGC = {
            R: mode.RLGC.R[i],
            L: mode.RLGC.L[i],
            G: mode.RLGC.G[i],
            C: mode.RLGC.C[i]
        };

        // Compute S-parameters
        const sp = computeSParamsSingleEnded(freq, RLGC, length, Z_ref);

        // Compute absolute differences
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

        const passed = errors.S11 < tolerances.S11 &&
                       errors.S21 < tolerances.S21 &&
                       errors.S12 < tolerances.S12 &&
                       errors.S22 < tolerances.S22;
        all_passed = all_passed && passed;

        const freqGHz = (freq / 1e9).toFixed(2);
        const status = passed ? '✓' : '✗';
        console.log(`${freqGHz.padEnd(8)} ${errors.S11.toFixed(4).padEnd(10)} ${errors.S21.toFixed(4).padEnd(10)} ${errors.S12.toFixed(4).padEnd(10)} ${errors.S22.toFixed(4).padEnd(10)} ${status}`);
    }

    console.log(`${'-'.repeat(70)}`);
    console.log(`Max absolute errors: S11=${maxErrors.S11.toFixed(4)}, S21=${maxErrors.S21.toFixed(4)}, S12=${maxErrors.S12.toFixed(4)}, S22=${maxErrors.S22.toFixed(4)}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error('S2P solve_sweep test failed');
    }
}

async function test_s4p_generation() {
    console.log(`\n${'='.repeat(80)}`);
    console.log('S4P GENERATION TEST');
    console.log(`${'='.repeat(80)}`);

    const reference = parseS4P_MA('./data/stripline_2d_diff_sweep.s4p');
    console.log(`Loaded ${reference.length} frequency points from reference`);

    // Setup solver
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

    // Use the new solve_sweep() API
    const sweepResults = await solver.solve_sweep({
        frequencies: reference.map(r => r.freq),
        energy_tol: 0.01
    });

    console.log(`Mesh: ${sweepResults.mesh.nx}x${sweepResults.mesh.ny}`);
    console.log(`Frequencies: ${sweepResults.frequencies.length} points`);
    console.log(`Modes: ${sweepResults.modes.map(m => m.mode).join(', ')}`);

    // Verify differential-specific arrays are present
    if (!sweepResults.Z_diff || !sweepResults.Z_common) {
        throw new Error('Z_diff and Z_common should be present for differential mode');
    }

    const length = 1.0;
    const Z_ref = 50;
    const tolerance = 0.25;

    let all_passed = true;
    const maxErrors = Array(4).fill(null).map(() => Array(4).fill(0));

    console.log(`\n${'Freq'.padEnd(8)} ${'Max|ΔSij|'.padEnd(12)} ${'Worst Sij'.padEnd(10)} Status`);
    console.log(`${'-'.repeat(50)}`);

    for (let i = 0; i < reference.length; i++) {
        const refPoint = reference[i];
        const freq = sweepResults.frequencies[i];

        const oddMode = sweepResults.modes.find(m => m.mode === 'odd');
        const evenMode = sweepResults.modes.find(m => m.mode === 'even');

        // Build RLGC objects from arrays
        const oddRLGC = {
            R: oddMode.RLGC.R[i],
            L: oddMode.RLGC.L[i],
            G: oddMode.RLGC.G[i],
            C: oddMode.RLGC.C[i]
        };
        const evenRLGC = {
            R: evenMode.RLGC.R[i],
            L: evenMode.RLGC.L[i],
            G: evenMode.RLGC.G[i],
            C: evenMode.RLGC.C[i]
        };

        // Compute 4-port S-parameters
        const sp = computeSParamsDifferential(freq, oddRLGC, evenRLGC, length, Z_ref);

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

        const freqGHz = (freq / 1e9).toFixed(2);
        const status = passed ? '✓' : '✗';
        console.log(`${freqGHz.padEnd(8)} ${maxError.toFixed(4).padEnd(12)} ${worstParam.padEnd(10)} ${status}`);
    }

    console.log(`${'-'.repeat(50)}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error('S4P solve_sweep test failed');
    }
}

async function test_s4p_generation_lossless() {
    console.log(`\n${'='.repeat(80)}`);
    console.log('S4P GENERATION TEST Lossless');
    console.log(`${'='.repeat(80)}`);

    const reference = parseS4P_MA('./data/stripline_2d_diff_lossless_fsweep.s4p');
    console.log(`Loaded ${reference.length} frequency points from reference`);

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

    const sweepResults = await solver.solve_sweep({
        frequencies: reference.map(r => r.freq),
        energy_tol: 0.01
    });

    console.log(`Mesh: ${sweepResults.mesh.nx}x${sweepResults.mesh.ny}`);
    console.log(`Frequencies: ${sweepResults.frequencies.length} points`);
    console.log(`Modes: ${sweepResults.modes.map(m => m.mode).join(', ')}`);

    // Verify differential-specific arrays are present
    if (!sweepResults.Z_diff || !sweepResults.Z_common) {
        throw new Error('Z_diff and Z_common should be present for differential mode');
    }

    const length = 50e-3;
    const Z_ref = 50;
    const tolerance = 0.25;

    let all_passed = true;
    const maxErrors = Array(4).fill(null).map(() => Array(4).fill(0));

    console.log(`\n${'Freq'.padEnd(8)} ${'Max|ΔSij|'.padEnd(12)} ${'Worst Sij'.padEnd(10)} Status`);
    console.log(`${'-'.repeat(50)}`);

    for (let i = 0; i < reference.length; i++) {
        const refPoint = reference[i];
        const freq = sweepResults.frequencies[i];

        const oddMode = sweepResults.modes.find(m => m.mode === 'odd');
        const evenMode = sweepResults.modes.find(m => m.mode === 'even');

        // Build RLGC objects from arrays
        const oddRLGC = {
            R: oddMode.RLGC.R[i],
            L: oddMode.RLGC.L[i],
            G: oddMode.RLGC.G[i],
            C: oddMode.RLGC.C[i]
        };
        const evenRLGC = {
            R: evenMode.RLGC.R[i],
            L: evenMode.RLGC.L[i],
            G: evenMode.RLGC.G[i],
            C: evenMode.RLGC.C[i]
        };

        // Compute 4-port S-parameters
        const sp = computeSParamsDifferential(freq, oddRLGC, evenRLGC, length, Z_ref);

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

        const freqGHz = (freq / 1e9).toFixed(2);
        const status = passed ? '✓' : '✗';
        console.log(`${freqGHz.padEnd(8)} ${maxError.toFixed(4).padEnd(12)} ${worstParam.padEnd(10)} ${status}`);
    }

    console.log(`${'-'.repeat(50)}`);
    console.log(`Overall Result: ${all_passed ? '✓ ALL TESTS PASSED' : '✗ SOME TESTS FAILED'}`);
    console.log(`${'='.repeat(80)}\n`);

    if (!all_passed) {
        throw new Error('S4P solve_sweep test failed');
    }
}

export { test_s4p_generation_lossless, test_s4p_generation, test_s2p_generation, test_s2p_generation2 };
