// Conductive dielectrics (diel sigma=, conductive_dielectric.js): the complex static
// solve with eps* = er (1 - j tand) - j sigma / (omega eps0) on both backends.
//
//   1. text: sigma on a dielectric parses and survives the conversion to text
//   2. parallel plate, oxide over conductive silicon between a full-width strip and
//      a ground wall: the field is one dimensional and C* = eps0 W / sum(h_i / eps_i*)
//      exactly, through the relaxation frequency where C and G change most
//   3. two-layer coax: C* = 2 pi eps0 / sum(ln(r_i+1 / r_i) / eps_i*)
//   4. sigma is the loss tangent sigma / (omega eps0 er): the same material written
//      either way (the complex path forced by a negligible sigma) gives the same line
//   5. far above the relaxation frequency sigma reduces to that loss tangent in the
//      perturbative path, on a microstrip and on symmetric and asymmetric pairs
//   6. backends agree on the on-chip microstrip, half and full domain agree, the QS
//      sweep path matches a direct solve, sigma=0 changes nothing
//   7. a substrate whose skin depth is not large next to its size gets a warning
//   8. the interpolating sweep of an asymmetric pair follows the frequency-dependent
//      [C] of the complex solve between its samples
//   9. DC: a parallel plate on two conductive layers in series has G = W / sum(h_i /
//      sigma_i) and the Maxwell-Wagner C, and an asymmetric pair on them has a finite
//      [G] at f = 0 equal to its low-frequency value
//  10. field plot: the plotted field is the complex static one, the silicon under a
//      coplanar line screens it below f_r and not above; a lossless fill plots no
//      imaginary part
//
// Run: node tests/test_conductive_dielectric.js
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { parseAndEvaluate, solverToGeometryText } from '../src/custom_geometry_text.js';
import { check, quiet, rel, APP, done } from './helpers.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';

const EPS0 = 8.854187817e-12;
const BACKENDS = [['QS', {}], ['full-wave', { mesh_backend: 'triangular' }]];
const fmt = v => v.toExponential(2);

async function solve(text, f, extra = {}, opts = APP) {
    const s = new CustomGeometrySolver({ text, nx: 30, ny: 30, freq: f, ...extra });
    // Nominal materials: the closed forms below have no causal dispersion.
    s.use_causal_materials = false;
    const r = await quiet(() => s.solve_adaptive(opts));
    return { s, r };
}

// Complex arithmetic for the closed forms.
const cdiv = (a, b) => {
    const d = b.re * b.re + b.im * b.im;
    return { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d };
};
const eps = (er, tand, sigma, w) => ({ re: er, im: -(er * tand + sigma / (w * EPS0)) });
// C' and G = omega C'' of C* = num / sum(len_i / eps_i*).
function seriesY(num, layers, w) {
    let s = { re: 0, im: 0 };
    for (const [len, e] of layers) { const q = cdiv({ re: len, im: 0 }, e); s = { re: s.re + q.re, im: s.im + q.im }; }
    const c = cdiv({ re: num, im: 0 }, s);
    return { C: c.re, G: -w * c.im };
}

// --- 1. Text ---
{
    const g = parseAndEvaluate('units um\nr = 0.5\ndiel x=0 y=0 w=1 h=1 er=11.9 sigma=1/(r*1e-2)');
    check('sigma on a dielectric evaluates as a plain number', g.errors.length === 0 && g.rects[0].sigma === 200);
    check('negative sigma on a dielectric is an error',
        parseAndEvaluate('diel x=0 y=0 w=1 h=1 er=2 sigma=-1').errors.length === 1);
    check('rq on a dielectric is still an error', parseAndEvaluate('diel x=0 y=0 w=1 h=1 er=2 rq=1').errors.length === 1);
    const text = 'units um\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=20 er=11.9 sigma=2\n'
        + 'diel x=-inf w=inf y=20 h=5 er=4.1\nsig+ x=-5 w=10 y=25 h=2\n';
    const s = new CustomGeometrySolver({ text });
    const back = new CustomGeometrySolver({ text: solverToGeometryText(s) });
    check('dielectric sigma survives the conversion to text',
        back.dielectrics.map(d => d.sigma).join() === s.dielectrics.map(d => d.sigma).join() && s.dielectrics[0].sigma === 2);
}

// --- 2. Parallel plate ---
// Silicon (20 um, 2 S/m, f_r = 3.0 GHz) under 5 um of oxide with a loss tangent, a strip
// spanning the domain width, ground wall below. Open side walls are natural boundaries,
// so the field is exactly one dimensional. The strip ends on the side walls, where the
// QS charge contour takes the half cells inside the domain.
{
    const W = 100e-6, hsi = 20e-6, hox = 5e-6, sig = 2;
    const text = `units um
bounds open open open gnd
domain -50 50 0 40
diel x=-inf w=inf y=0 h=20 er=11.9 sigma=${sig}
diel x=-inf w=inf y=20 h=5 er=4.1 tand=0.01
sig+ x=-50 w=100 y=25 h=2
`;
    for (const [name, extra] of BACKENDS) {
        for (const f of [0.3e9, 3e9, 30e9]) {
            const w = 2 * Math.PI * f;
            const ex = seriesY(EPS0 * W, [[hox, eps(4.1, 0.01, 0, w)], [hsi, eps(11.9, 0, sig, w)]], w);
            const m = (await solve(text, f, extra)).r.modes[0];
            const C = m.RLGC.C, G = m.RLGC.G;
            check(`${name}: parallel plate at ${f / 1e9} GHz, C and G match the closed form`,
                rel(C, ex.C) < 2e-4 && rel(G, ex.G) < 2e-4,
                `C ${fmt(C)} vs ${fmt(ex.C)}, G ${fmt(G)} vs ${fmt(ex.G)}`);
        }
    }
}

// --- 3. Two-layer coax (full-wave: n-gons) ---
{
    const a = 0.2e-3, rm = 0.5e-3, b = 1.5e-3, sig = 0.05;
    const text = `units mm
bounds open open open open
domain -1.7 1.7 -1.7 1.7
diel ngon x=0 y=0 r=${b * 1e3} n=256 er=11.9 sigma=${sig}
diel ngon x=0 y=0 r=${rm * 1e3} n=256 er=4.1
sig+ ngon x=0 y=0 r=${a * 1e3} n=256
gnd ngon x=0 y=0 r=1.6 r_in=${b * 1e3} n=256
`;
    for (const f of [0.03e9, 0.3e9]) {
        const w = 2 * Math.PI * f;
        const ex = seriesY(2 * Math.PI * EPS0, [[Math.log(rm / a), eps(4.1, 0, 0, w)], [Math.log(b / rm), eps(11.9, 0, sig, w)]], w);
        const m = (await solve(text, f, { mesh_backend: 'triangular' })).r.modes[0];
        check(`full-wave: two-layer coax at ${f / 1e9} GHz, C and G match the closed form`,
            rel(m.RLGC.C, ex.C) < 2e-3 && rel(m.RLGC.G, ex.G) < 2e-3,
            `C ${fmt(m.RLGC.C)} vs ${fmt(ex.C)}, G ${fmt(m.RLGC.G)} vs ${fmt(ex.G)}`);
    }
}

// On-chip microstrip: copper over SiO2, aluminium ground, conductive silicon below.
const chip = (si, extra = '') => `units um
bounds open open open open
diel x=-inf w=inf y=-2.6 h=-inf er=11.9 ${si}
diel x=-inf w=inf y=-2.6 h=12.6 er=4.1 tand=0.001
gnd x=-23 w=46 y=0 h=0.49 sigma=3.5e7
sig+ x=-5.75 w=11.5 y=7.66 h=3
${extra}`;
// Pairs over the same stack, symmetric or with a wider negative trace.
const pair = (si, asym) => `units um
bounds open open open open
diel x=-inf w=inf y=-30 h=30 er=11.9 ${si}
diel x=-inf w=inf y=0 h=8 er=4.1
gnd x=-60 w=120 y=-0.5 h=0.5
sig- x=${asym ? -20 : -14} w=${asym ? 12 : 10} y=8 h=2
sig+ x=4 w=10 y=8 h=2
`;
const SIG = 2, FR = SIG / (2 * Math.PI * EPS0 * 11.9);
const tandOf = f => SIG / (2 * Math.PI * f * EPS0 * 11.9);
const lines = [['microstrip', chip], ['symmetric pair', si => pair(si, false)], ['asymmetric pair', si => pair(si, true)]];
const modesOf = r => r.modes.map(m => m.RLGC);

// --- 4. sigma is a loss tangent in the complex solve ---
for (const [name, extra] of BACKENDS) {
    for (const [lname, geo] of lines) {
        const f = FR;
        const a = modesOf((await solve(geo(`sigma=${SIG}`), f, extra)).r);
        const b = modesOf((await solve(geo(`tand=${tandOf(f)} sigma=1e-30`), f, extra)).r);
        const d = Math.max(...a.flatMap((m, i) => [rel(m.C, b[i].C), rel(m.G, b[i].G)]));
        check(`${name}: ${lname}, sigma and its loss tangent give the same C and G`, d < 1e-6, `max rel ${fmt(d)}`);
    }
}

// --- 5. Far above the relaxation frequency: the perturbative loss tangent ---
for (const [name, extra] of BACKENDS) {
    for (const [lname, geo] of lines) {
        const f = 30 * FR;
        const { r: ra } = await solve(geo(`sigma=${SIG}`), f, extra);
        const { r: rb } = await solve(geo(`tand=${tandOf(f)}`), f, extra);
        const a = modesOf(ra), b = modesOf(rb);
        // The two G estimators of the QS backend (charge flux of the complex solve,
        // cell-wise field integral of the perturbative one) differ by its
        // discretization of G, a few percent on these grids.
        const tolG = name === 'QS' ? 0.05 : 2e-3;
        const dC = Math.max(...a.map((m, i) => rel(m.C, b[i].C)));
        const dG = Math.max(...a.map((m, i) => rel(m.G, b[i].G)));
        check(`${name}: ${lname} at 30 f_r, sigma matches the perturbative loss tangent`,
            dC < 1e-4 && dG < tolG, `C ${fmt(dC)}, G ${fmt(dG)}`);
        if (ra.RLGC_matrix) {
            const G = ra.RLGC_matrix.G, Gb = rb.RLGC_matrix.G;
            const dM = Math.max(rel(G[0][0], Gb[0][0]), rel(G[0][1], Gb[0][1]), rel(G[1][1], Gb[1][1]));
            check(`${name}: ${lname} at 30 f_r, the G matrix matches too`, dM < tolG, `max rel ${fmt(dM)}`);
            if (lname === 'asymmetric pair') check(`${name}: asymmetric pair keeps its physical matrices`, !!ra.physMatrix);
        }
    }
}

// --- 6. Consistency ---
{
    const text = chip(`sigma=${SIG}`);
    const qs = (await solve(text, FR)).r.modes[0].RLGC, fw = (await solve(text, FR, BACKENDS[1][1])).r.modes[0].RLGC;
    const qs0 = (await solve(chip(''), FR)).r.modes[0].RLGC, fw0 = (await solve(chip(''), FR, BACKENDS[1][1])).r.modes[0].RLGC;
    check('on-chip microstrip at f_r, the silicon conductance of the two backends agrees',
        rel(qs.G - qs0.G, fw.G - fw0.G) < 0.03, `QS ${fmt(qs.G - qs0.G)} vs full-wave ${fmt(fw.G - fw0.G)}`);
    check('on-chip microstrip at f_r, C of the two backends agrees', rel(qs.C, fw.C) < 5e-3, `${fmt(qs.C)} vs ${fmt(fw.C)}`);
    // The silicon screens like a floating conductor below f_r and stops conducting above it.
    const lo = (await solve(text, FR / 30)).r.modes[0].RLGC, hi = (await solve(text, FR * 30)).r.modes[0].RLGC;
    const lo0 = (await solve(chip(''), FR / 30)).r.modes[0].RLGC, hi0 = (await solve(chip(''), FR * 30)).r.modes[0].RLGC;
    check('the conductive silicon raises C below f_r and leaves it above', lo.C > lo0.C * (1 + 5e-4) && rel(hi.C, hi0.C) < 5e-5,
        `below ${fmt(lo.C / lo0.C - 1)}, above ${fmt(hi.C / hi0.C - 1)}`);

    for (const [name, extra] of BACKENDS) {
        const half = modesOf((await solve(pair(`sigma=${SIG}`, false), FR, extra)).r);
        const full = modesOf((await solve(pair(`sigma=${SIG}`, false), FR, { ...extra, symmetry: false })).r);
        // The two adaptive solves end on different grids, the lossless pair differs by
        // about 2e-3 the same way.
        const tol = 0.01;
        const d = Math.max(...half.flatMap((m, i) => [rel(m.C, full[i].C), rel(m.G, full[i].G)]));
        check(`${name}: symmetric pair, half and full domain agree`, d < tol, `max rel ${fmt(d)}`);
        const plain = modesOf((await solve(pair('', false), FR, extra)).r);
        const zero = modesOf((await solve(pair('sigma=0', false), FR, extra)).r);
        check(`${name}: sigma=0 changes nothing`, plain.every((m, i) => m.C === zero[i].C && m.G === zero[i].G));
    }

    // QS sweep path (computeAtFrequency re-solves per frequency) against a direct solve.
    for (const [lname, geo] of lines) {
        const { s, r } = await solve(geo(`sigma=${SIG}`), 10e9);
        const sw = modesOf(await quiet(() => s.computeAtFrequency(1e9, r)));
        s.freq = 1e9;
        const direct = modesOf(await quiet(() => s.solve_adaptive({ ...APP, skip_mesh: true })));
        const d = Math.max(...sw.flatMap((m, i) => [rel(m.C, direct[i].C), rel(m.G, direct[i].G)]));
        check(`QS: ${lname}, the sweep path matches a direct solve on the same grid`, d < 1e-9, `max rel ${fmt(d)}`);
    }
}

// --- 7. Warning ---
for (const [name, extra] of BACKENDS) {
    const reasons = r => (r.warnings || []).map(w => w.reason);
    const ok = (await solve(chip(`sigma=${SIG}`), 10e9, extra)).r;
    // 1e5 S/m: skin depth 16 um at 10 GHz, the substrate in the domain is thicker.
    const bad = (await solve(chip('sigma=1e5'), 10e9, extra)).r;
    check(`${name}: a 50 ohm*cm substrate gets no conductive-dielectric warning`, !reasons(ok).includes('conductive-dielectric'));
    check(`${name}: a 0.001 ohm*cm substrate is warned about`, reasons(bad).includes('conductive-dielectric'),
        reasons(bad).join(', '));
}

// --- 8. interpolating sweep carries the per-frequency [C] ---
{
    // No ground between the traces and the silicon: [C] moves through the relaxation.
    const text = 'units um\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=300 er=11.9 sigma=10\n' +
        'diel x=-inf w=inf y=300 h=10 er=3.9\nsig- x=-60 w=40 y=310 h=2\nsig+ x=10 w=15 y=310 h=2\n';
    const s = new CustomGeometrySolver({ text, nx: 30, ny: 30, freq: 1e10 });
    const fs = [2.2e8, 1.3e9, 6e9];
    const { ir, ex } = await quiet(async () => {
        const base = await s.solve_adaptive({ max_iters: 6, max_nodes: 15000, param_tol: 0.03 });
        const sw = new InterpolatingSweep(s, base, { tolerance: 0.01 });
        await sw.run(1e8, 1e10);
        const ex = [];
        for (const f of fs) ex.push(await s.computeAtFrequency(f, base));
        return { ir: sw.buildResults(fs), ex };
    });
    fs.forEach((f, i) => {
        const a = ir[i].result.RLGC_matrix && ir[i].result.RLGC_matrix.C, b = ex[i].RLGC_matrix.C;
        check(`interpolated asymmetric pair [C] at ${fmt(f)} Hz matches a direct solve (0.5%)`,
            !!a && rel(a[0][0], b[0][0]) < 0.005 && rel(a[1][1], b[1][1]) < 0.005
                && Math.abs(a[0][1] - b[0][1]) < 0.005 * b[0][0],
            a ? `C11 ${fmt(a[0][0])} / ${fmt(b[0][0])}, C12 ${fmt(a[0][1])} / ${fmt(b[0][1])}` : 'no RLGC_matrix');
    });
}

// --- 9. DC conductance ---
// Two conductive layers in series under a strip spanning the domain (one dimensional,
// as in 2). At DC the current divides by conductivity: G = W / sum(h_i / sigma_i), and
// the charge on the layer interface leaves C = eps0 W sum(h_i er_i / sigma_i^2) /
// sum(h_i / sigma_i)^2 (Maxwell-Wagner).
{
    const W = 100e-6, L = [[20e-6, 11.9, 2], [5e-6, 4.1, 0.5]];
    const R = L.reduce((a, [h, , sg]) => a + h / sg, 0);
    const Gdc = W / R, Cdc = EPS0 * W * L.reduce((a, [h, er, sg]) => a + h * er / (sg * sg), 0) / (R * R);
    const stack = `units um
bounds open open open gnd
domain -50 50 0 40
diel x=-inf w=inf y=0 h=20 er=11.9 sigma=2
diel x=-inf w=inf y=20 h=5 er=4.1 sigma=0.5
`;
    for (const [name, extra] of BACKENDS) {
        const m = (await solve(stack + 'sig+ x=-50 w=100 y=25 h=2\n', 0, extra)).r.modes[0];
        check(`${name}: parallel plate at DC, G = W / sum(h / sigma) and the Maxwell-Wagner C`,
            rel(m.RLGC.G, Gdc) < 2e-4 && rel(m.RLGC.C, Cdc) < 2e-4,
            `G ${fmt(m.RLGC.G)} vs ${fmt(Gdc)}, C ${fmt(m.RLGC.C)} vs ${fmt(Cdc)}`);
        // Asymmetric pair: the strip split by a 5 um gap into 40 and 55 um traces.
        const pair = stack + 'sig- x=-50 w=40 y=25 h=2\nsig+ x=-5 w=55 y=25 h=2\n';
        const dc = (await solve(pair, 0, extra)).r, lf = (await solve(pair, 1e3, extra)).r;
        const G0 = dc.RLGC_matrix.G, G1 = lf.RLGC_matrix.G;
        check(`${name}: asymmetric pair at DC has the low-frequency [G]`,
            !!dc.physMatrix && [[0, 0], [0, 1], [1, 1]].every(([i, j]) => Math.abs(G0[i][j] - G1[i][j]) < 1e-6 * G1[0][0]),
            `G11 ${fmt(G0[0][0])} / ${fmt(G1[0][0])}, G12 ${fmt(G0[0][1])} / ${fmt(G1[0][1])}`);
        // Both traces at 1 V: the plate less the gap, plus the gap's fringing.
        const both = G0[0][0] + G0[1][1] + 2 * G0[0][1];
        check(`${name}: asymmetric pair at DC, both traces driven conduct between the gapped and the whole plate`,
            both > 0.95 * Gdc && both < Gdc, `${fmt(both)} vs plate ${fmt(Gdc)}`);
    }
}

// --- 10. Field plot ---
{
    // Coplanar line on oxide over conductive silicon, nothing below the silicon.
    const cpw = si => `units um
bounds open open open open
diel x=-inf w=inf y=-50 h=50 er=11.9 ${si}
diel x=-inf w=inf y=0 h=6 er=4.1
gnd x=-30 w=22 y=3 h=1
gnd x=8 w=22 y=3 h=1
sig+ x=-5 w=10 y=3 h=1
`;
    // Peak |E|^2 in the silicon over that in the oxide below the line, and the largest
    // |Im E| over the largest |E|, from the plotted grid of mode 0.
    const plotStats = s => {
        const P = s.getPlotFields(), { x, y } = P;
        let si = 0, ox = 0, im = 0, all = 0;
        for (let i = 0; i < y.length; i++) for (let j = 0; j < x.length; j++) {
            const xi = P.ExIm ? P.ExIm[0][i][j] : 0, yi = P.EyIm ? P.EyIm[0][i][j] : 0;
            const e2 = P.Ex[0][i][j] ** 2 + P.Ey[0][i][j] ** 2 + xi * xi + yi * yi;
            im = Math.max(im, Math.hypot(xi, yi)); all = Math.max(all, Math.sqrt(e2));
            if (Math.abs(x[j]) > 15e-6) continue;
            if (y[i] < -1e-6 && y[i] > -20e-6) si = Math.max(si, e2);
            else if (y[i] > 0.5e-6 && y[i] < 2.5e-6) ox = Math.max(ox, e2);
        }
        return { r: si / ox, im: im / all, mirrored: !P.ExIm || P.ExIm[0][0].length === x.length };
    };
    for (const [name, extra] of BACKENDS) {
        const { s, r } = await solve(cpw(`sigma=${SIG}`), FR, extra);
        const at = async f => { await quiet(() => s.plotFieldsAt(f, r)); return plotStats(s); };
        const lo = await at(FR / 100), mid = await at(FR), hi = await at(FR * 100);
        check(`${name}: the silicon screens the plotted field below f_r and not above`, hi.r > 100 * lo.r,
            `Si / oxide |E|^2 ${fmt(lo.r)} at f_r/100, ${fmt(hi.r)} at 100 f_r`);
        check(`${name}: the plotted field has an imaginary part at f_r`, mid.im > 10 * Math.max(lo.im, hi.im),
            `|Im E| / |E| ${fmt(lo.im)}, ${fmt(mid.im)}, ${fmt(hi.im)}`);
        check(`${name}: the imaginary field is mirrored with the rest`, mid.mirrored);
        const plain = await solve(cpw(''), FR, extra);
        await quiet(() => plain.s.plotFieldsAt(FR, plain.r));
        check(`${name}: a lossless fill plots no imaginary part`, plotStats(plain.s).im < 1e-12,
            fmt(plotStats(plain.s).im));
    }
}

done();
