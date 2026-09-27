// Line types the fixed templates cannot describe, built as custom geometry and checked
// against a reference each, on both backends:
//   1. microstrip on stacked dielectrics, single-ended and differential
//   2. microstrip on a finite-width ground with air and open boundaries all around
//   3. CPW on a finite substrate over air, against the conformal-mapping closed form
//   4. slotline, full-wave effective permittivity against Janaswamy-Schaubert
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { check, quiet, rel, APP, done } from './helpers.js';

const pct = v => `${(100 * v).toFixed(2)}%`;

const SOLVE = { ...APP, max_iters: 8 };
const TRI = { mesh_backend: 'triangular' };
async function solve(text, extra = {}) {
    const s = new CustomGeometrySolver({ text, nx: 30, ny: 30, freq: 1e9, ...extra });
    const r = await quiet(() => s.solve_adaptive(SOLVE));
    return r.modes.map(m => ({ Z0: m.Z0?.re ?? m.Z0, eps: m.eps_eff, R: m.RLGC.R, G: m.RLGC.G }));
}
// Z0 and eps_eff within `tol`, R within `rTol`, mode by mode.
function agree(name, a, b, tol, rTol) {
    a.forEach((m, i) => {
        const n = a.length > 1 ? `${name} mode ${i}` : name;
        check(`${n}: Z0 and eps_eff`, rel(m.Z0, b[i].Z0) < tol && rel(m.eps, b[i].eps) < tol,
            `Z0 ${m.Z0.toFixed(2)} / ${b[i].Z0.toFixed(2)}, eps ${m.eps.toFixed(4)} / ${b[i].eps.toFixed(4)}`);
        if (rTol) check(`${n}: R`, rel(m.R, b[i].R) < rTol, `${m.R.toFixed(2)} / ${b[i].R.toFixed(2)} ohm/m`);
    });
}

const MS = { trace_width: 0.3e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.3,
    tan_delta: 0.02, sigma_cond: 5.8e7, nx: 30, ny: 30, freq: 1e9 };
const native = await quiet(async () => {
    const r = await new MicrostripSolver(MS).solve_adaptive(SOLVE);
    return r.modes.map(m => ({ Z0: m.Z0?.re ?? m.Z0, eps: m.eps_eff, R: m.RLGC.R, G: m.RLGC.G }));
});

// --- 1. Stacked dielectrics ---
const stack = (l2, trace) => `units mm\nbounds open open open gnd\n` +
    `diel x=-inf w=inf y=0 h=0.1 er=4.3 tand=0.02\ndiel x=-inf w=inf y=0.1 h=0.4 ${l2}\n${trace}`;
const SINGLE = 'sig+ x=-0.15 w=0.3 y=0.5 h=0.035\n';
const PAIR = 'sig- x=-0.4 w=0.3 y=0.5 h=0.035\nsig+ x=0.1 w=0.3 y=0.5 h=0.035\n';
agree('two layers of one material = native microstrip', await solve(stack('er=4.3 tand=0.02', SINGLE)), native, 5e-3, 0.03);
{
    const q = await solve(stack('er=2.2 tand=0.001', SINGLE)), t = await solve(stack('er=2.2 tand=0.001', SINGLE), TRI);
    agree('stacked microstrip, QS vs full-wave', q, t, 0.01, 0.10);
    check('stacked microstrip: low-eps layer lowers eps_eff', q[0].eps < native[0].eps - 0.5, `${q[0].eps.toFixed(3)} vs ${native[0].eps.toFixed(3)}`);
    const qd = await solve(stack('er=2.2 tand=0.001', PAIR)), td = await solve(stack('er=2.2 tand=0.001', PAIR), TRI);
    agree('stacked differential microstrip, QS vs full-wave', qd, td, 0.015, 0.12);
}

// --- 2. Finite ground, open on all sides ---
const finiteGnd = gw => `units mm\nbounds open open open open\ndiel x=-inf w=inf y=0 h=0.5 er=4.3 tand=0.02\n` +
    `gnd x=${-gw / 2} w=${gw} y=-0.035 h=0.035\n${SINGLE}`;
{
    const wide = await solve(finiteGnd(8)), narrow = await solve(finiteGnd(0.6));
    agree('wide finite ground = native microstrip', wide, native, 5e-3, 0.03);
    check('narrow ground raises Z0', narrow[0].Z0 > 1.1 * wide[0].Z0, `${narrow[0].Z0.toFixed(2)} vs ${wide[0].Z0.toFixed(2)}`);
    agree('finite ground 2 mm, QS vs full-wave', await solve(finiteGnd(2)), await solve(finiteGnd(2), TRI), 0.01, 0.10);
    agree('finite ground 0.6 mm, QS vs full-wave', narrow, await solve(finiteGnd(0.6), TRI), 0.01, 0.20);
}

// --- 3. CPW on a finite substrate over air ---
{
    const K = k => { let a = 1, b = Math.sqrt(1 - k * k); for (let i = 0; i < 40; i++) [a, b] = [(a + b) / 2, Math.sqrt(a * b)]; return Math.PI / (2 * a); };
    const Kp = k => K(Math.sqrt(1 - k * k));
    const w = 0.2, g = 0.1, h = 0.635, er = 9.8, a = w / 2, b = w / 2 + g;
    const k0 = a / b, k1 = Math.sinh(Math.PI * a / (2 * h)) / Math.sinh(Math.PI * b / (2 * h));
    const eps = 1 + (er - 1) / 2 * (K(k1) / Kp(k1)) * (Kp(k0) / K(k0));
    const ref = [{ Z0: 30 * Math.PI / Math.sqrt(eps) * Kp(k0) / K(k0), eps }];
    const cpw = (t, gl, gr) => `units mm\nbounds open open open open\ndiel x=-inf w=inf y=${-h} h=${h} er=${er} tand=0.001\n` +
        `gnd ${gl} y=0 h=${t}\ngnd ${gr} y=0 h=${t}\nsig+ x=${-a} w=${w} y=0 h=${t}\n`;
    const full = [`x=${-b} w=-inf`, `x=${b} w=inf`], finite = [`x=${-b - 1} w=1`, `x=${b} w=1`];
    // Thin metal against the zero-thickness closed form: the metal thickness lowers Z0 a little.
    agree('CPW over air, QS vs conformal mapping', await solve(cpw(0.002, ...full)), ref, 0.03);
    agree('CPW over air, full-wave vs conformal mapping', await solve(cpw(0.002, ...full), TRI), ref, 0.03);
    // Finite, floating grounds: the full-wave mode pick has to stay on the CPW mode.
    const o = { freq: 5e9 };
    const qf = await solve(cpw(0.017, ...finite), o), tf = await solve(cpw(0.017, ...finite), { ...TRI, ...o });
    agree('CPW with finite grounds, QS vs full-wave', qf, tf, 0.015, 0.12);
    agree('CPW finite grounds ~ grounds to the wall', qf, await solve(cpw(0.017, ...full), o), 0.015);
}

// --- 4. Slotline ---
{
    const W = 0.1, h = 0.635, er = 9.8, f = 10e9;
    const slot = dom => `units mm\nbounds open open open open\ndomain ${dom}\n` +
        `diel x=-inf w=inf y=${-h} h=${h} er=${er} tand=0.001\nsig+ x=${-W / 2} w=-inf y=0 h=0.005\ngnd x=${W / 2} w=inf y=0 h=0.005\n`;
    // Janaswamy-Schaubert, 9.7 <= er <= 20, 0.02 <= W/h < 0.2 (about 2% accurate).
    const hl = h / (299.792458 / (f / 1e9)), wh = W / h, lg = Math.log10;
    const ratio = 0.923 - 0.448 * lg(er) + 0.2 * wh - (0.29 * wh + 0.047) * lg(hl * 100);
    const epsRef = 1 / (ratio * ratio);
    const s = new CustomGeometrySolver({ text: slot('-5 5 -5 5'), nx: 30, ny: 30, freq: f });
    check('slotline: warns that the line has no quasi-static limit', s.openBoundaryWarnings().some(m => m.includes('slotline')));
    const t5 = await solve(slot('-5 5 -5 5'), { ...TRI, freq: f }), t10 = await solve(slot('-10 10 -10 10'), { ...TRI, freq: f });
    check('slotline: full-wave eps_eff vs Janaswamy-Schaubert', rel(t5[0].eps, epsRef) < 0.06,
        `${t5[0].eps.toFixed(3)} vs ${epsRef.toFixed(3)}, ${pct(rel(t5[0].eps, epsRef))}`);
    check('slotline: full-wave eps_eff does not depend on the domain', rel(t5[0].eps, t10[0].eps) < 0.01,
        `${t5[0].eps.toFixed(4)} vs ${t10[0].eps.toFixed(4)}`);
}

done();
