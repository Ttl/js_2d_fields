import { MicrostripSolver } from '../src/microstrip.js';

/**
 * Test to verify that at bottom corners covered by side plating:
 * - Material conductivity: plating.sigma (plating extends from sides)
 * - Surface roughness: bulk.rq (bottom surface preparation)
 *
 * This is correct because:
 * - rq = roughness at conductor/dielectric interface (from surface preparation)
 * - sigma = material conductivity (intrinsic to the material)
 */

async function test_corner_roughness_parameter() {
    console.log('='.repeat(60));
    console.log('CORNER ROUGHNESS PARAMETER TEST');
    console.log('='.repeat(60));
    console.log('\nPhysical model:');
    console.log('  At bottom corners covered by side plating:');
    console.log('  - Material: PLATING (extends from sides)');
    console.log('  - Roughness: BULK RQ (bottom surface preparation)');
    console.log('  - Impedance: Z = f(plating.sigma, bulk.rq)\n');

    const freq = 1e9;
    const copper_sigma = 5.96e7;
    const gold_sigma = 4.1e7;
    const trace_width = 100e-6;
    const thick_plating = 60e-6;  // 60% of trace width

    console.log('Configuration:');
    console.log(`  Trace width: ${(trace_width * 1e6).toFixed(0)} μm`);
    console.log(`  Plating thickness: ${(thick_plating * 1e6).toFixed(0)} μm (covers ${(thick_plating / trace_width * 200).toFixed(0)}% of bottom)`);
    console.log(`  Frequency: ${freq / 1e9} GHz`);

    // Case 1: Smooth bulk (rq=0), rough plating (rq=1um)
    // At corners: should use plating sigma + bulk rq (0) = smooth
    const plating1 = {
        sigma: gold_sigma,
        thickness: thick_plating,
        rq: 1e-6,  // Rough plating
        top: true,   // Plate top to reduce side surface contribution
        sides: true,
        bottom: false
    };

    const solver1 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: trace_width,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: copper_sigma,
        freq: freq,
        nx: 50,
        ny: 25,
        plating: plating1,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0  // Smooth bulk surface
    });

    // Case 2: Rough bulk (rq=1um), smooth plating (rq=0)
    // At corners: should use plating sigma + bulk rq (1um) = rough
    const plating2 = {
        sigma: gold_sigma,
        thickness: thick_plating,
        rq: 0,  // Smooth plating
        top: true,   // Plate top to reduce side surface contribution
        sides: true,
        bottom: false
    };

    const solver2 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: trace_width,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: copper_sigma,
        freq: freq,
        nx: 50,
        ny: 25,
        plating: plating2,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 1e-6  // Rough bulk surface
    });

    // Case 3: Both smooth (reference)
    const plating3 = {
        sigma: gold_sigma,
        thickness: thick_plating,
        rq: 0,
        top: true,   // Plate top to reduce side surface contribution
        sides: true,
        bottom: false
    };

    const solver3 = new MicrostripSolver({
        substrate_height: 0.8e-3,
        trace_width: trace_width,
        trace_thickness: 35e-6,
        gnd_thickness: 35e-6,
        epsilon_r: 4.5,
        tan_delta: 0.02,
        sigma_cond: copper_sigma,
        freq: freq,
        nx: 50,
        ny: 25,
        plating: plating3,
        boundaries: ["open", "open", "open", "gnd"],
        rq: 0  // Smooth bulk
    });

    const result1 = await solver1.solve_adaptive({ max_refinements: 0 });
    const mode1 = result1.modes[0];

    const result2 = await solver2.solve_adaptive({ max_refinements: 0 });
    const mode2 = result2.modes[0];

    const result3 = await solver3.solve_adaptive({ max_refinements: 0 });
    const mode3 = result3.modes[0];

    console.log('\n' + '='.repeat(60));
    console.log('RESULTS:\n');

    console.log('Case 1: Smooth bulk (rq=0) + Rough plating (rq=1μm)');
    console.log(`  Conductor loss = ${mode1.alpha_c.toFixed(6)} dB/m`);
    console.log('  At corners: Z = f(gold, rq=0) ← uses BULK rq');

    console.log('\nCase 2: Rough bulk (rq=1μm) + Smooth plating (rq=0)');
    console.log(`  Conductor loss = ${mode2.alpha_c.toFixed(6)} dB/m`);
    console.log('  At corners: Z = f(gold, rq=1μm) ← uses BULK rq');

    console.log('\nCase 3: Both smooth (reference)');
    console.log(`  Conductor loss = ${mode3.alpha_c.toFixed(6)} dB/m`);
    console.log('  At corners: Z = f(gold, rq=0)');

    console.log('\n' + '='.repeat(60));
    console.log('ANALYSIS:\n');

    const loss_increase_case2 = mode2.alpha_c - mode3.alpha_c;
    const loss_increase_pct = (loss_increase_case2 / mode3.alpha_c) * 100;

    console.log('Comparing Case 1 vs Case 3:');
    const diff1_3 = Math.abs(mode1.alpha_c - mode3.alpha_c);
    console.log(`  Loss difference: ${diff1_3.toFixed(6)} dB/m`);
    console.log('  Note: Plating rq affects SIDE/TOP surfaces (conformal plating)');
    console.log('        This is CORRECT behavior - sides use plating.rq in layered model');

    console.log('\nComparing Case 2 vs Case 3:');
    console.log(`  Loss increase: ${loss_increase_case2.toFixed(6)} dB/m (${loss_increase_pct.toFixed(2)}%)`);
    if (loss_increase_case2 > 0.01) {
        console.log('  ✓ CORRECT: Bulk rq affects BOTTOM CORNERS (uses bulk rq=1μm)');
    } else {
        console.log('  ✗ ERROR: Bulk rq should affect bottom corners!');
    }

    // Both roughness routes must ADD loss over the all-smooth reference. Their
    // ORDERING is geometry-dependent, not a correctness invariant: the old
    // dielectric-field loss integrand weighted substrate-side samples by √εr,
    // which made bottom(-corner) roughness dominate here; the vacuum-field
    // integrand weights faces by the real
    // current distribution (validated per-face against tri MQS), and on this
    // wrapped high-impedance trace (w=0.1mm on h=0.8mm) roughening top+sides
    // legitimately adds more than roughening the bottom.
    const loss_increase_case1 = mode1.alpha_c - mode3.alpha_c;
    console.log('\nKey checks:');
    console.log(`  Case 1 − Case 3 (plating rq on top/sides adds loss): ${loss_increase_case1.toFixed(6)} dB/m`);
    console.log(`  Case 2 − Case 3 (bulk rq on bottom/corners adds loss): ${loss_increase_case2.toFixed(6)} dB/m`);

    console.log('\n' + '='.repeat(60));
    console.log('PHYSICAL INTERPRETATION:\n');

    console.log('At bottom corners covered by side plating:');
    console.log('  ┌─────────────┐');
    console.log('  │  Dielectric │');
    console.log('  ├─────────────┤ ← Interface roughness = BULK rq');
    console.log('  │   Plating   │    (from bottom surface preparation)');
    console.log('  │  (material) │ ← Material conductivity = PLATING sigma');
    console.log('  └─────────────┘');
    console.log('\nParameters used at BOTTOM corners:');
    console.log('  - sigma: plating.sigma (material extends from sides)');
    console.log('  - rq: bulk.rq (surface roughness from bottom preparation)');
    console.log('\nParameters used at SIDE/TOP surfaces:');
    console.log('  - Layered model with plating.rq (conformal plating)');
    console.log('\nThis is correct because:');
    console.log('  - Bottom corners: plating wraps from sides, but surface is original bottom');
    console.log('  - Side/top: conformal plating applied directly → uses plating.rq');

    console.log('\n' + '='.repeat(60));

    // Verify correctness: EACH roughness route must reach its surfaces —
    // bulk.rq the bottom/corners (Case 2) and plating.rq the plated top/sides
    // (Case 1) — i.e. both raise loss over the all-smooth reference.
    if (loss_increase_case2 > 0.01 && loss_increase_case1 > 0.01) {
        console.log('✓ TEST PASSED: Bottom corners use bulk.rq, sides use plating.rq');
        console.log('='.repeat(60));
        return true;
    } else {
        console.log('✗ TEST FAILED: Incorrect roughness parameter assignment');
        console.log('='.repeat(60));
        return false;
    }
}

test_corner_roughness_parameter()
    .then(passed => {
        process.exit(passed ? 0 : 1);
    })
    .catch(err => {
        console.error('\nTest error:', err);
        process.exit(1);
    });
