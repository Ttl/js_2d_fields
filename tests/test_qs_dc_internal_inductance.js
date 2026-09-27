// Quasi-static internal inductance of the signal traces below and through the skin
// transition, against the full-wave MQS solve. The grounds are made nearly lossless
// (sigma 1e10) so only the traces' own internal inductance is compared, which is what
// the quasi-static DC solve (_dc_signal_inductance) and its blend with the surface
// value provide. Before them the quasi-static value read half the full-wave one at
// low frequency on a microstrip (the surface field of the skin regime, not uniform
// current). Also: the DC limit of a stripline trace, whose faces carry equal fields,
// and the differential modes on the half domain.
import { CustomGeometrySolver } from '../src/custom_geometry.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
async function quiet(fn) {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
}
const rel = (a, b) => Math.abs(a - b) / Math.abs(b);

const SOLVE = { max_iters: 12, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };
// Internal inductance per mode (nH/m) at each frequency.
async function lint(text, freqs, tri) {
    const s = new CustomGeometrySolver({ text, sigma_cond: 1e6, freq: Math.max(...freqs), nx: 30, ny: 30,
        ...(tri ? { mesh_backend: 'triangular' } : {}) });
    if (tri) s.tri_opts = { lossMethod: 'auto' };
    return quiet(async () => {
        const r0 = await s.solve_adaptive(SOLVE);
        const out = [];
        for (const f of freqs) {
            const r = await s.computeAtFrequency(f, r0);
            out.push(Object.fromEntries(r.modes.map(m => [m.mode, m.L_internal * 1e9])));
        }
        return out;
    });
}
async function compare(name, text, freqs, tols) {
    const [q, t] = await Promise.all([lint(text, freqs, false), lint(text, freqs, true)]);
    freqs.forEach((f, k) => {
        for (const mode of Object.keys(t[k])) {
            check(`${name} ${mode} @ ${f / 1e6} MHz: QS internal L within ${100 * tols[k]}% of full-wave`,
                rel(q[k][mode], t[k][mode]) < tols[k], `${q[k][mode].toFixed(2)} vs ${t[k][mode].toFixed(2)} nH/m`);
        }
    });
    return q;
}

const GND = 'gnd x=-4 w=8 y=-0.035 h=0.035 sigma=1e10';
// Microstrip, delta = 159 / 50 / 14.5 um on a 35 um trace.
const ms = await compare('microstrip',
    `units mm\nbounds open open open open\ndiel x=-inf w=inf y=0 h=0.21 er=4.4\n${GND}\nsig+ x=-0.175 w=0.35 y=0.21 h=0.035`,
    [1e7, 1e8, 1.2e9], [0.04, 0.08, 0.04]);
check('microstrip: the DC value is about twice the old surface-model plateau (~16.7 nH/m)', ms[0].single > 30,
    `${ms[0].single.toFixed(2)} nH/m`);
// Differential pair on the half domain: the odd mode has a PEC plane, the even mode a PMC plane.
await compare('differential microstrip',
    `units mm\nbounds open open open open\ndiel x=-inf w=inf y=0 h=0.2 er=4.4\n${GND}\nsig+ x=0.1 w=0.3 y=0.2 h=0.035 mirror=1`,
    [1e7, 1.2e9], [0.05, 0.05]);
// Stripline: equal fields on both faces. The quasi-static surface value itself reads
// ~8% low here in the skin regime (10 GHz), which the 1.2 GHz point inherits.
await compare('stripline',
    `units mm\nbounds open open open open\ndiel x=-4 w=8 y=0 h=0.5 er=4\n${GND}\ngnd x=-4 w=8 y=0.5 h=0.035 sigma=1e10\n` +
    'sig+ x=-0.1 w=0.2 y=0.2325 h=0.035',
    [1e7, 1.2e9], [0.06, 0.10]);

if (failures) { console.log(`\n${failures} CHECK(S) FAILED`); process.exit(1); }
console.log('\nALL CHECKS PASSED');
