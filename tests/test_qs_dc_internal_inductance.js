// Quasi-static internal inductance of the signal traces below and through the skin
// transition, against the full-wave MQS solve. The grounds are made nearly lossless
// (sigma 1e10) so only the traces' own internal inductance is compared, which is what
// the quasi-static DC solve (_dc_signal_inductance) and its blend with the surface
// value provide. Before them the quasi-static value read half the full-wave one at
// low frequency on a microstrip (the surface field of the skin regime, not uniform
// current). Also: the DC limit of a stripline trace, whose faces carry equal fields,
// and the differential modes on the half domain.
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { check, quiet, relErr as rel, APP, done } from './helpers.js';


const SOLVE = { ...APP, max_iters: 12 };
// Internal inductance per mode (nH/m) at each frequency.
async function lint(text, freqs, tri) {
    const s = new CustomGeometrySolver({ text, sigma_cond: 1e6, freq: Math.max(...freqs), nx: 30, ny: 30,
        ...(tri ? { mesh_backend: 'triangular' } : {}) });
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

// Plated and multi-metal traces: faces modelled as plating metal alone (thick corners,
// top plating down the sides) and a thin block of another metal must level off below
// the skin transition like the bare trace, at the full-wave offset over it.
async function lintMs(opts, freqs, tri) {
    const { MicrostripSolver } = await import('../src/microstrip.js');
    const s = new MicrostripSolver({ trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6,
        gnd_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
        boundaries: ['open', 'open', 'open', 'gnd'], nx: 30, ny: 30, ...opts, ...(tri ? { mesh_backend: 'triangular' } : {}) });
    return quiet(async () => {
        const r0 = await s.solve_adaptive(SOLVE);
        const out = [];
        for (const f of freqs) out.push((await s.computeAtFrequency(f, r0)).modes[0].L_internal * 1e9);
        return out;
    });
}
const PLATING = { sigma: 1e7, thickness: 4e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: false };
for (const [name, plating] of [['thick corners', { ...PLATING, thick_corners: true }], ['top-only plating', { ...PLATING, sides: false }]]) {
    const [q, t] = await Promise.all([lintMs({ plating }, [0, 1e3, 1e6], false), lintMs({ plating }, [1e6], true)]);
    check(`${name}: QS internal L flat below the skin transition`, q[0] < 1.1 * q[2] && Math.abs(q[1] - q[0]) < 0.01 * q[0],
        `${q.map(v => v.toFixed(2)).join(' / ')} nH/m at DC / 1 kHz / 1 MHz`);
    // 1 MHz is the start of the skin transition (delta = 1.9 t), where the quasi-static
    // DC/skin blend reads about 5% low on the bare trace too.
    check(`${name} @ 1 MHz: QS internal L within 6% of full-wave`, rel(q[2], t[0]) < 0.06, `${q[2].toFixed(2)} vs ${t[0].toFixed(2)} nH/m`);
}
// Roughness is a surface layer far inside the skin depth at low frequency, a constant
// excess inductance: L_int falls monotonically from DC on both backends.
{
    const FR = [0, 1e2, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9];
    const [q, t] = await Promise.all([lintMs({ rq: 1e-6 }, FR, false), lintMs({ rq: 1e-6 }, FR, true)]);
    for (const [name, v] of [['QS', q], ['full-wave', t]]) {
        // 0.5% slack: the full-wave value from 100 Hz carries its small low-frequency bias over the exact DC one.
        const ok = v.every((x, i) => i === 0 || x <= v[i - 1] * 1.005);
        check(`rough 1 um, ${name}: internal L non-increasing from DC to 1 GHz`, ok, v.map(x => x.toFixed(2)).join(' / '));
    }
}
{
    const text = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.21 er=4.4\n' +
        'sig+ x=-0.175 w=0.35 y=0.21 h=0.031 sigma=5.8e7\nsig+ x=-0.175 w=0.35 y=0.241 h=0.004 sigma=1e7';
    const run = tri => { const s = new CustomGeometrySolver({ text, sigma_cond: 5.8e7, freq: 1e9, nx: 30, ny: 30,
        ...(tri ? { mesh_backend: 'triangular' } : {}) }); return s; };
    const at = async (s, freqs) => quiet(async () => {
        const r0 = await s.solve_adaptive(SOLVE);
        const out = [];
        for (const f of freqs) out.push((await s.computeAtFrequency(f, r0)).modes[0].L_internal * 1e9);
        return out;
    });
    const [q, t] = await Promise.all([at(run(false), [0, 1e6]), at(run(true), [1e6])]);
    check('nickel block on copper: QS internal L flat below the skin transition', q[0] < 1.1 * q[1],
        `${q.map(v => v.toFixed(2)).join(' / ')} nH/m at DC / 1 MHz`);
    check('nickel block on copper @ 1 MHz: QS internal L within 4% of full-wave', rel(q[1], t[0]) < 0.04,
        `${q[1].toFixed(2)} vs ${t[0].toFixed(2)} nH/m`);
}

done();

console.log('\nALL CHECKS PASSED');
