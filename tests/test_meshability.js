// Validates the pre-mesh meshability guard (_check_meshability):
//  - normal geometries pass on both backends,
//  - geometries that need an impossibly dense mesh (large domain + high frequency,
//    or an absurd feature/domain ratio) are rejected with a clear error.
import { MicrostripSolver } from '../src/microstrip.js';

let pass = 0, fail = 0;
function check(name, fn) {
    try { fn(); console.log(`  PASS: ${name}`); pass++; }
    catch (e) { console.log(`  FAIL: ${name}\n        ${e.message}`); fail++; }
}

// A normal microstrip should be meshable on both backends and at high frequency.
function normal(backend, freq) {
    return new MicrostripSolver({
        substrate_height: 0.5e-3, trace_width: 0.5e-3, trace_thickness: 35e-6,
        gnd_thickness: 35e-6, epsilon_r: 4.3, tan_delta: 0.02,
        freq, mesh_backend: backend, boundaries: ["open", "open", "open", "gnd"],
    });
}

// A 1 m × (auto) structure at 100 GHz — the classic "too dense to solve" case.
function huge(backend) {
    return new MicrostripSolver({
        substrate_height: 1.0, trace_width: 1.0, trace_thickness: 35e-6,
        gnd_thickness: 35e-6, epsilon_r: 4.3, freq: 100e9,
        enclosure_width: 1.0, enclosure_height: 1.0,
        mesh_backend: backend, boundaries: ["open", "open", "open", "gnd"],
    });
}

const expectsOk = (s) => s._check_meshability(20000);
const expectsReject = (s, label) => {
    try { s._check_meshability(20000); throw new Error(`${label}: expected rejection, got none`); }
    catch (e) {
        if (/expected rejection/.test(e.message)) throw e;
        if (!/cannot be meshed/.test(e.message)) throw new Error(`${label}: wrong error: ${e.message}`);
        console.log(`        rejected as expected: ${e.message.slice(0, 120)}…`);
    }
};

console.log('Normal geometries should pass:');
check('rectilinear microstrip @ 1 GHz', () => expectsOk(normal('rectilinear', 1e9)));
check('rectilinear microstrip @ 100 GHz', () => expectsOk(normal('rectilinear', 100e9)));
check('triangular microstrip @ 1 GHz', () => expectsOk(normal('triangular', 1e9)));
check('triangular microstrip @ 100 GHz', () => expectsOk(normal('triangular', 100e9)));

console.log('\nPathological geometries should be rejected:');
check('rectilinear 1 m @ 100 GHz', () => expectsReject(huge('rectilinear'), 'rect-huge'));
check('triangular 1 m @ 100 GHz', () => expectsReject(huge('triangular'), 'tri-huge'));

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
