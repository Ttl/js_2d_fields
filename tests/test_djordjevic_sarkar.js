/**
 * Test script for Djordjevic-Sarkar causal material model
 * Compare with Python reference implementation
 */

import { djordjevic_sarkar } from '../src/djordjevic_sarkar.js';

// Test cases matching the Python implementation
const testCases = [
    { freq: 1e9, eps_ref: 4.5, tand_ref: 0.02, f_ref: 1e9 },
    { freq: 10e9, eps_ref: 4.5, tand_ref: 0.02, f_ref: 1e9 },
    { freq: 100e6, eps_ref: 4.5, tand_ref: 0.02, f_ref: 1e9 },
    { freq: 1e9, eps_ref: 3.0, tand_ref: 0.01, f_ref: 1e9 },
];

console.log('Testing Djordjevic-Sarkar Model');
console.log('================================\n');

for (const testCase of testCases) {
    const { freq, eps_ref, tand_ref, f_ref } = testCase;
    const result = djordjevic_sarkar(freq, eps_ref, tand_ref, f_ref);

    console.log(`Input: freq=${(freq/1e9).toFixed(3)} GHz, eps_ref=${eps_ref}, tand_ref=${tand_ref}`);
    console.log(`Output: eps_real=${result.eps_real.toFixed(6)}, tand_actual=${result.tand_actual.toFixed(6)}`);
    console.log('');
}

// Test at reference frequency - should match input values
console.log('Verification at reference frequency (1 GHz):');
const refTest = djordjevic_sarkar(1e9, 4.5, 0.02, 1e9);
console.log(`eps_ref=4.5, result=${refTest.eps_real.toFixed(6)} (should be ~4.5)`);
console.log(`tand_ref=0.02, result=${refTest.tand_actual.toFixed(6)} (should be ~0.02)`);
console.log('');

// Test frequency dependence
console.log('Frequency sweep for FR4 (eps_r=4.5, tand=0.02 @ 1GHz):');
console.log('Freq (GHz) | eps_r    | tan(d)');
console.log('-----------|----------|----------');
const frequencies = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0];
const sweep = [];
for (const f_ghz of frequencies) {
    const freq = f_ghz * 1e9;
    const result = djordjevic_sarkar(freq, 4.5, 0.02, 1e9);
    sweep.push(result);
    console.log(`${f_ghz.toString().padStart(10)} | ${result.eps_real.toFixed(6)} | ${result.tand_actual.toFixed(6)}`);
}

// --- Self-checks ------------------------------------------------------------
let failed = 0;
function check(name, cond) {
    console.log(`  ${cond ? '✓' : '✗'} ${name}`);
    if (!cond) failed++;
}

console.log('\nAssertions:');
// At the reference frequency the model must reproduce the input values.
check('eps_real(f_ref) == eps_ref', Math.abs(refTest.eps_real - 4.5) < 4.5 * 1e-3);
check('tand_actual(f_ref) ≈ tand_ref', Math.abs(refTest.tand_actual - 0.02) < 0.02 * 0.02);
// Causal (normal) dispersion: permittivity decreases monotonically with frequency.
let monotonic = true;
for (let i = 1; i < sweep.length; i++) if (sweep[i].eps_real >= sweep[i - 1].eps_real) monotonic = false;
check('eps_real strictly decreasing with frequency', monotonic);
// Everything finite and physical across the band.
check('all eps_real finite and > 1', sweep.every(r => Number.isFinite(r.eps_real) && r.eps_real > 1));
check('all tand_actual finite and > 0', sweep.every(r => Number.isFinite(r.tand_actual) && r.tand_actual > 0));

console.log(failed ? `\n✗ ${failed} check(s) failed` : '\n✓ All checks passed');
process.exit(failed ? 1 : 0);
