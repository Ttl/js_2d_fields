import { MicrostripSolver } from '../src/microstrip.js';

/**
 * Test case from user requirement:
 * "If plating thickness is higher than trace width, and plating rq=0 and bulk rq=0,
 *  and sides=top=enabled but bottom=disabled, then it should equal the case where
 *  plating conductivity was applied to bulk without plating with rq=0."
 */

async function test_plating_thicker_than_trace() {
    console.log('Testing: Plating thickness > trace width');
    console.log('Expected: Should equal solid conductor of plating material\n');

    const freq = 1e9;
    const gold_sigma = 4.1e7;
    const copper_sigma = 5.96e7;
    const trace_width = 100e-6;  // 100 μm

    // Plating thickness GREATER than trace width
    const plating_thickness = 150e-6;  // 150 μm > 100 μm trace width

    console.log(`Trace width: ${(trace_width * 1e6).toFixed(1)} μm`);
    console.log(`Plating thickness: ${(plating_thickness * 1e6).toFixed(1)} μm`);
    console.log(`Ratio: ${(plating_thickness / trace_width).toFixed(2)}×\n`);

    // Case 1: Thick plating on top + sides, bottom disabled
    const thick_plating = {
        sigma: gold_sigma,
        thickness: plating_thickness,
        rq: 0,  // No roughness
        top: true,
        sides: true,
        bottom: false,  // Bottom disabled as per requirement
        thick_corners: true  // Model thick plating covering bottom via sides
    };

    const solver1 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: trace_width,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: copper_sigma,  // Copper bulk
        freq: freq,
        nx: 50,
        ny: 25,
        plating: thick_plating,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0  // No bulk roughness
    });

    // Case 2: Solid gold conductor (no plating)
    const solver2 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: trace_width,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: gold_sigma,  // Gold bulk
        freq: freq,
        nx: 50,
        ny: 25,
        plating: null,  // No plating
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0  // No roughness
    });

    const result1 = await solver1.solve_adaptive({ max_refinements: 0 });
    const mode1 = result1.modes[0];

    const result2 = await solver2.solve_adaptive({ max_refinements: 0 });
    const mode2 = result2.modes[0];

    console.log('Case 1: Plating > trace width, top + sides (bottom disabled)');
    console.log(`  Z0 = ${mode1.Z0.toFixed(4)} Ω`);
    console.log(`  Conductor loss = ${mode1.alpha_c.toFixed(6)} dB/m`);

    console.log('\nCase 2: Solid gold conductor (reference)');
    console.log(`  Z0 = ${mode2.Z0.toFixed(4)} Ω`);
    console.log(`  Conductor loss = ${mode2.alpha_c.toFixed(6)} dB/m`);

    const z0_diff = Math.abs(mode1.Z0 - mode2.Z0);
    const z0_rel_diff = z0_diff / mode2.Z0;
    const alpha_diff = Math.abs(mode1.alpha_c - mode2.alpha_c);
    const alpha_rel_diff = alpha_diff / mode2.alpha_c;

    console.log(`\nZ0 difference: ${(z0_rel_diff * 100).toFixed(3)}%`);
    console.log(`Conductor loss difference: ${(alpha_rel_diff * 100).toFixed(3)}%`);

    const tolerance = 0.01;  // 1%

    console.log('\n' + '='.repeat(60));
    if (alpha_rel_diff <= tolerance) {
        console.log('✓ PASS: Thick plating equals solid conductor');
        console.log(`  (within ${(tolerance * 100)}% tolerance)`);
        return true;
    } else {
        console.log('✗ FAIL: Thick plating does NOT equal solid conductor');
        console.log(`  Expected: < ${(tolerance * 100)}% difference`);
        console.log(`  Actual: ${(alpha_rel_diff * 100).toFixed(3)}% difference`);
        console.log('\nPHYSICAL EXPLANATION:');
        console.log('  When bottom surface is NOT plated, it exposes the copper bulk.');
        console.log('  Current flow at the bottom surface experiences copper conductivity,');
        console.log('  not gold conductivity, leading to higher loss.');
        console.log('\nSOLUTION:');
        console.log('  To achieve equivalence, enable bottom plating:');
        console.log('    plating.bottom = true');
        console.log('  OR model the geometric coverage of thick side plating extending');
        console.log('  across the bottom surface when thickness > trace_width.');
        return false;
    }
}

test_plating_thicker_than_trace()
    .then(passed => {
        console.log('='.repeat(60));
        process.exit(passed ? 0 : 1);
    })
    .catch(err => {
        console.error('\nTest error:', err);
        process.exit(1);
    });
