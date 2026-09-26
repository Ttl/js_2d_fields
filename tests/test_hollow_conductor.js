// A hollow stripline trace against the solid one on the full-wave solver.
//
// Inside a closed conductor there is no field, so each wall of a hollow trace is a
// slab driven from one side. Its surface resistance relative to thick metal is
//   F(x) = Re[(1+j) coth((1+j) x)] = (sinh 2x + sin 2x) / (cosh 2x - cos 2x),  x = wall / delta,
// which dips to tanh(pi/2) = 0.917 at x = pi/2: a wall of about 1.6 skin depths loses
// less than solid metal (the reflected wave from the inner surface cancels part of the
// current), and a thinner one loses more. The solid trace of a centred stripline is a
// slab fed from both faces, F(t / 2 delta). Theory: R_hollow / R_solid = F(wall/delta) / F(t/2delta).
//
// The trace is wide (1.2 mm) so the 1-D slab picture holds over most of its perimeter;
// the edges, where the field is two-dimensional, still move the ratio by about 1 %. The
// ground planes are made nearly lossless (solver sigma 1e4 x copper, the trace has its own
// copper sigma) so R is the trace's. A copper shell filled with a poor conductor (two
// touching blocks, meshed metal) and the same trace written as a poor conductor with
// thick copper plating (plating=all with Model Thick Plating, which the full-wave solver
// meshes as a layer of metal)
// must both match the hollow trace. So must the drawing that overlaps a poor core onto
// a copper trace: where conductors of one kind overlap, the later line's metal wins.
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

const SIGMA = 5.8e7, WALL = 5e-6, T = 35e-6;
const F = x => (Math.sinh(2 * x) + Math.sin(2 * x)) / (Math.cosh(2 * x) - Math.cos(2 * x));
const delta = f => Math.sqrt(1 / (Math.PI * f * 4e-7 * Math.PI * SIGMA));
// The dip of F: a wall of pi/2 skin depths.
const fDip = 1 / (Math.PI * 4e-7 * Math.PI * SIGMA * (2 * WALL / Math.PI) ** 2);
const FREQS = [50e6, 200e6, fDip, 1e9, 5e9];

const BASE = `units um\nw = 1200; t = 35; b = 400; tw = 5\nbounds open open gnd gnd\ndomain -2000 2000 0 b\n` +
    'diel x=-inf w=inf y=0 h=b er=3.5\n';
const TRACE = 'x=-w/2 w=w y=b/2-t/2 h=t';
async function rSweep(trace, freqs = FREQS, extra = {}) {
    const s = new CustomGeometrySolver({ text: BASE + trace, freq: 1e9, mesh_backend: 'triangular', sigma_cond: SIGMA * 1e4, ...extra });
    s.tri_opts = { lossMethod: 'auto' };
    await quiet(() => s.solve_adaptive({ max_iters: 8, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 }));
    const out = [];
    for (const f of freqs) {
        const m = (await quiet(() => s.computeAtFrequency(f))).modes[0];
        out.push({ R: m.RLGC.R, via: m.lossVia });
    }
    return out;
}

const solid = await rSweep(`sig+ ${TRACE} sigma=${SIGMA}`);
const hollow = await rSweep(`sig+ ${TRACE} wall=tw sigma=${SIGMA}`);
const shell = await rSweep(`sig+ ${TRACE} wall=tw sigma=${SIGMA}\nsig+ x=-w/2+tw w=w-2*tw y=b/2-t/2+tw h=t-2*tw sigma=${SIGMA / 1000}`);
// Thick plating: meshed. The layered surface model of thin plating needs a thick, good
// bulk and warns here.
const plated = await rSweep(`plating sigma=${SIGMA} t=tw\nsig+ ${TRACE} sigma=${SIGMA / 1000} plating=all`, FREQS, { thick_plating: true });
check('all solves use the MQS loss', [...solid, ...hollow, ...shell, ...plated].every(r => r.via === 'mqs'));

FREQS.forEach((f, i) => {
    const ratio = hollow[i].R / solid[i].R;
    const theory = F(WALL / delta(f)) / F(T / 2 / delta(f));
    const tag = `${(f / 1e6).toFixed(0)} MHz, wall ${(WALL / delta(f)).toFixed(2)} skin depths`;
    const detail = `R_hollow/R_solid ${ratio.toFixed(4)}, theory ${theory.toFixed(4)}`;
    if (i === 0) {
        // Skin depth comparable to the whole trace: the 1-D slab no longer describes the
        // solid trace, only the sign is checked.
        check(`${tag}: the thin wall loses more`, ratio > 1.5, detail);
    } else {
        check(`${tag}: ratio = slab theory`, Math.abs(ratio / theory - 1) < 0.02, detail);
    }
    check(`${tag}: copper shell on a poor core = hollow`, Math.abs(shell[i].R / hollow[i].R - 1) < 0.01,
        `${shell[i].R.toFixed(4)} vs ${hollow[i].R.toFixed(4)} ohm/m`);
    check(`${tag}: copper plating on a poor conductor = hollow`, Math.abs(plated[i].R / hollow[i].R - 1) < 0.01,
        `${plated[i].R.toFixed(4)} vs ${hollow[i].R.toFixed(4)} ohm/m`);
});
// Overlapping drawings: the core over the copper is the plated trace, the copper over
// the core solid copper. At DC the overlap counts once.
const CORE = `sig+ x=-w/2+tw w=w-2*tw y=b/2-t/2+tw h=t-2*tw sigma=${SIGMA / 1000}`;
const over = await rSweep(`sig+ ${TRACE} sigma=${SIGMA}\n${CORE}`, [0, fDip]);
const under = await rSweep(`${CORE}\nsig+ ${TRACE} sigma=${SIGMA}`, [fDip]);
check('poor core drawn over a copper trace = hollow trace', Math.abs(over[1].R / hollow[2].R - 1) < 0.01,
    `${over[1].R.toFixed(4)} vs ${hollow[2].R.toFixed(4)} ohm/m`);
check('copper trace drawn over the core = solid trace', Math.abs(under[0].R / solid[2].R - 1) < 0.01,
    `${under[0].R.toFixed(4)} vs ${solid[2].R.toFixed(4)} ohm/m`);
{
    const w = 1200e-6, aCore = (w - 2 * WALL) * (T - 2 * WALL);
    const rDc = 1 / (SIGMA * (w * T - aCore) + SIGMA / 1000 * aCore);
    check('overlap at DC: the shared area counts once', Math.abs(over[0].R / rDc - 1) < 1e-3,
        `${over[0].R.toFixed(4)} vs ${rDc.toFixed(4)} ohm/m`);
}

// Thin plating over the same poor bulk is a layered surface impedance, which needs a
// thick bulk: it has to say so and point to thick plating.
{
    const s = new CustomGeometrySolver({ text: BASE + `plating sigma=${SIGMA} t=tw\nsig+ ${TRACE} sigma=${SIGMA / 1000} plating=all`,
        freq: 1e9, mesh_backend: 'triangular', sigma_cond: SIGMA * 1e4 });
    s.tri_opts = { lossMethod: 'auto' };
    await quiet(() => s.solve_adaptive({ max_iters: 8, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 }));
    const r = await quiet(() => s.computeAtFrequency(fDip));
    const w = [...(r.warnings || []), ...(s.modeWarnings || [])].find(x => x.reason === 'plating-transition');
    check('thin plating on a poor bulk warns and names thick plating', !!w && /Model Thick Plating/.test(w.message),
        w ? w.message.slice(0, 90) : 'no warning');
}

const dip = hollow[2].R / solid[2].R;
check('a wall of pi/2 skin depths loses less than solid metal', dip < 0.95, `ratio ${dip.toFixed(4)}, theory ${F(Math.PI / 2).toFixed(4)}`);

console.log(failures ? `\n${failures} check(s) failed` : '\nAll checks passed');
process.exit(failures ? 1 : 0);
