/**
 * Test to verify that causal materials correctly affect eps_eff and capacitance
 */

import { MicrostripSolver } from '../src/microstrip.js';
import { applyDjordjevicSarkar } from '../src/djordjevic_sarkar.js';

// Backend selection: `MESH_BACKEND=triangular node tests/test_causal_effect.js`
// runs the triangular full-wave backend; default is rectilinear.
const MESH_BACKEND = process.env.MESH_BACKEND || 'rectilinear';

// Create a simple stripline test case
const options = {
    trace_width: 0.2e-3,
    substrate_height: 0.2e-3,
    trace_thickness: 35e-6,
    epsilon_r: 4.1,
    tan_delta: 0.02,
    sigma_cond: 5.8e7,
    freq: 1e9,
    enclosure_height: 0.2e-3 + 35e-6,
    nx: 30,
    ny: 30,
    boundaries: ["open", "open", "gnd", "gnd"],  // Stripline
    epsilon_r_top: 4.1,
    tan_delta_top: 0.02,
    mesh_backend: MESH_BACKEND,
};

const solver = new MicrostripSolver(options);

(async () => {
console.log('Testing Causal Material Effect on eps_eff and Capacitance');
console.log('==========================================================\n');

// Generate mesh
solver.ensure_mesh();
console.log(solver.x && solver.y
    ? `Mesh generated: ${solver.x.length}x${solver.y.length}`
    : `Mesh generated (backend=${MESH_BACKEND})`);

// Solve at 1 GHz (reference frequency)
console.log('\n1. Solving at 1 GHz (reference frequency)...');
const result_1ghz = await solver.solve_adaptive({
    max_iters: 3,
    energy_tol: 0.05,
    param_tol: 0.05,
    max_nodes: 100000
});

const mode_1ghz = result_1ghz.modes[0];
console.log(`   eps_eff = ${mode_1ghz.eps_eff.toFixed(6)}`);
console.log(`   C = ${(mode_1ghz.C * 1e12).toFixed(6)} pF/m`);
console.log(`   Z0 = ${mode_1ghz.Z0.toFixed(6)} Ω`);

// Now test at 67 GHz without causal materials
console.log('\n2. Computing at 67 GHz WITHOUT causal materials...');
solver.freq = 67e9;
solver.use_causal_materials = false;
const result_67ghz_noncausal = await solver.computeAtFrequency(67e9, result_1ghz);
const mode_67ghz_nc = result_67ghz_noncausal.modes[0];
console.log(`   eps_eff = ${mode_67ghz_nc.eps_eff.toFixed(6)} (should be ~4.1)`);
console.log(`   C = ${(mode_67ghz_nc.C * 1e12).toFixed(6)} pF/m`);
console.log(`   Z0 = ${mode_67ghz_nc.Z0.toFixed(6)} Ω`);

// Now test at 67 GHz WITH causal materials
console.log('\n3. Computing at 67 GHz WITH causal materials...');
solver.freq = 67e9;
solver.use_causal_materials = true;
const result_67ghz_causal = await solver.computeAtFrequency(67e9, result_1ghz);
const mode_67ghz_c = result_67ghz_causal.modes[0];
console.log(`\n   Results with causal materials:`);
console.log(`   eps_eff = ${mode_67ghz_c.eps_eff.toFixed(6)} (should be ~3.89)`);
console.log(`   C = ${(mode_67ghz_c.C * 1e12).toFixed(6)} pF/m`);
console.log(`   Z0 = ${mode_67ghz_c.Z0.toFixed(6)} Ω`);

// Compare
console.log('\n4. Comparison:');
console.log(`   Δeps_eff = ${(mode_67ghz_nc.eps_eff - mode_67ghz_c.eps_eff).toFixed(6)}`);
console.log(`   ΔC = ${((mode_67ghz_nc.C - mode_67ghz_c.C) * 1e12).toFixed(6)} pF/m`);
console.log(`   Causal effect: ${((1 - mode_67ghz_c.eps_eff / mode_67ghz_nc.eps_eff) * 100).toFixed(2)}% reduction in eps_eff`);

// --- Self-checks --------------------------------------------------------
// Causal materials must LOWER eps_eff (and hence C) at 67 GHz relative to the
// non-causal case, matching the ~4.1 → ~3.89 shift the Djordjevic-Sarkar model
// predicts for FR4 (εr=4.1, tanδ=0.02 @ 1 GHz).
let failed = 0;
function check(name, cond) {
    console.log(`  ${cond ? '✓' : '✗'} ${name}`);
    if (!cond) failed++;
}

console.log('\nAssertions:');
check('non-causal eps_eff ≈ 4.1', Math.abs(mode_67ghz_nc.eps_eff - 4.1) < 4.1 * 0.03);
check('causal eps_eff ≈ 3.89', Math.abs(mode_67ghz_c.eps_eff - 3.89) < 3.89 * 0.03);
check('causal eps_eff < non-causal eps_eff', mode_67ghz_c.eps_eff < mode_67ghz_nc.eps_eff);
check('causal C < non-causal C', mode_67ghz_c.C < mode_67ghz_nc.C);
const reduction = 1 - mode_67ghz_c.eps_eff / mode_67ghz_nc.eps_eff;
check('eps_eff reduction in 3%–8% band', reduction > 0.03 && reduction < 0.08);

console.log(failed ? `\n✗ ${failed} check(s) failed` : '\n✓ All checks passed');
process.exit(failed ? 1 : 0);
})();
