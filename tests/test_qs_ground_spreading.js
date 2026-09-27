// Return current spreading sideways in a ground at low frequency, on the quasi-static
// backend (the thin-sheet ground solve, FieldSolver2D._ground_sheet_setup /
// _ground_sheet_impedance) against the full-wave MQS solve, which meshes the ground:
//   1. microstrip on a 50 mm ground in open space (the 1 kHz Ansys case of test_vs_ref),
//      R and L through the transition from full spreading to the surface model
//   2. the blend between the two models leaves L monotone in frequency
//   3. a differential pair on the same ground: half domain = full domain, and full-wave
//   4. a ground plane on the domain edge is a wall on both backends: no spreading
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { MicrostripSolver } from '../src/microstrip.js';

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

// R (ohm/m) and L (nH/m) per mode at each frequency.
async function sweep(solver, freqs, tri) {
    if (tri) solver.tri_opts = { lossMethod: 'auto' };
    return quiet(async () => {
        const r0 = await solver.solve_adaptive();
        const out = [];
        for (const f of freqs) {
            const r = await solver.computeAtFrequency(f, r0);
            out.push(Object.fromEntries(r.modes.map(m => [m.mode, { R: m.RLGC.R, L: m.RLGC.L * 1e9 }])));
        }
        return out;
    });
}
const custom = (text, extra = {}) => new CustomGeometrySolver({ text, sigma_cond: 1e7, freq: 1e6, nx: 10, ny: 10, ...extra });
const TRI = { mesh_backend: 'triangular' };

const MS = `units mm
bounds open open open open
diel x=-25 w=50 y=0 h=1.6 er=4.5
gnd x=-25 w=50 y=-0.035 h=0.035
sig+ x=-1.5 w=3 y=1.6 h=0.035`;

// 1. Full spreading (1 kHz) through the blend (spread 0.2 .. 0.02 near 0.1 .. 1 MHz).
{
    const F = [1e3, 3e4, 1e5, 3e5, 1e6];
    const [q, t] = await Promise.all([sweep(custom(MS), F, false), sweep(custom(MS, TRI), F, true)]);
    F.forEach((f, k) => {
        const a = q[k].single, b = t[k].single;
        check(`microstrip, 50 mm ground @ ${f / 1e3} kHz: L within 2% of full-wave`, rel(a.L, b.L) < 0.02,
            `${a.L.toFixed(1)} vs ${b.L.toFixed(1)} nH/m`);
        check(`microstrip, 50 mm ground @ ${f / 1e3} kHz: R within 3% of full-wave`, rel(a.R, b.R) < 0.03,
            `${a.R.toFixed(4)} vs ${b.R.toFixed(4)} ohm/m`);
    });
    check('microstrip: spreading raises L at 1 kHz well above its 1 MHz value', q[0].single.L > 1.4 * q[4].single.L,
        `${q[0].single.L.toFixed(1)} vs ${q[4].single.L.toFixed(1)} nH/m`);

    // 2. Monotone through the blend, on a fine frequency grid.
    const G = Array.from({ length: 25 }, (_, i) => 1e4 * Math.pow(10, i / 8));
    const g = await sweep(custom(MS), G, false);
    const L = g.map(p => p.single.L);
    check('microstrip: L non-increasing from 10 kHz to 10 MHz across the model blend',
        L.every((v, i) => i === 0 || v <= L[i - 1] * (1 + 1e-6)), L.map(v => v.toFixed(1)).join(' '));
}

// 3. Differential pair: both symmetry-plane conditions of the half domain against the
//    full domain, and full-wave.
{
    const PAIR = MS.replace('sig+ x=-1.5 w=3 y=1.6 h=0.035', 'sig+ x=1 w=3 y=1.6 h=0.035 mirror=1');
    const F = [1e3, 1e5];
    const [half, full, t] = await Promise.all([sweep(custom(PAIR), F, false),
        sweep(custom(PAIR, { symmetry: false }), F, false), sweep(custom(PAIR, TRI), F, true)]);
    F.forEach((f, k) => {
        for (const mode of ['odd', 'even']) {
            const a = half[k][mode], b = full[k][mode], c = t[k][mode];
            check(`pair ${mode} @ ${f / 1e3} kHz: half domain = full domain`, rel(a.L, b.L) < 2e-3 && rel(a.R, b.R) < 2e-3,
                `L ${a.L.toFixed(2)} / ${b.L.toFixed(2)} nH/m, R ${a.R.toFixed(4)} / ${b.R.toFixed(4)}`);
            check(`pair ${mode} @ ${f / 1e3} kHz: L within 2% of full-wave`, rel(a.L, c.L) < 0.02,
                `${a.L.toFixed(1)} vs ${c.L.toFixed(1)} nH/m`);
        }
    });
}

// 4. A ground plane on the domain edge: the triangular backend makes it the wall, and the
//    quasi-static backend keeps its surface model, so neither spreads.
{
    const ms = extra => new MicrostripSolver({ substrate_height: 1.6e-3, trace_width: 3e-3, trace_thickness: 35e-6,
        gnd_thickness: 35e-6, epsilon_r: 4.5, tan_delta: 0, sigma_cond: 1e7, enclosure_width: 50e-3, freq: 1e6,
        nx: 10, ny: 10, boundaries: ['open', 'open', 'open', 'gnd'], ...extra });
    const [q, t] = await Promise.all([sweep(ms({}), [1e3], false), sweep(ms(TRI), [1e3], true)]);
    check('ground plane on the domain edge @ 1 kHz: no spreading, backends agree',
        rel(q[0].single.L, t[0].single.L) < 0.02 && q[0].single.L < 350,
        `${q[0].single.L.toFixed(1)} vs ${t[0].single.L.toFixed(1)} nH/m`);
}

if (failures) { console.log(`\n${failures} CHECK(S) FAILED`); process.exit(1); }
console.log('\nALL CHECKS PASSED');
