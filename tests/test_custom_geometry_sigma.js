// Per-conductor conductivity in custom geometry (sigma= on a conductor line), on both
// backends.
//
//   1. a value equal to the solver-wide one changes nothing
//   2. every conductor carrying its own sigma = the solver-wide sigma set to that value
//   3. a poor conductor as the trace only lands between the two uniform lines
//   4. low frequency: R tends to the series DC resistance of the two metals, and
//      plating several skin depths thick equals a trace made of the plating metal
//   5. a pair of two different metals leaves the half domain and reads the mean of the
//      two uniform pairs, whichever trace carries which metal
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { parseAndEvaluate, solverToGeometryText } from '../src/custom_geometry_text.js';

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
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b));

const SOLVE = { max_iters: 8, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };
function build(text, extra = {}) {
    const s = new CustomGeometrySolver({ text, nx: 30, ny: 30, freq: 5e9, ...extra });
    if (s.mesh_backend === 'triangular') s.tri_opts = { lossMethod: 'auto' };
    return s;
}
async function solveR(text, extra = {}) {
    const s = build(text, extra);
    const r = await quiet(() => s.solve_adaptive(SOLVE));
    return r.modes.map(m => ({ R: m.RLGC.R, L: m.RLGC.L, via: m.lossVia }));
}
const TRI = { mesh_backend: 'triangular' };
const CU = 5.8e7, LOW = 1e7;

// --- Text key ---
{
    const g = parseAndEvaluate('units um\ns = 1e7\nsig+ x=0 y=0 w=100 h=35 sigma=s*2');
    check('sigma evaluates as a plain number', g.errors.length === 0 && g.rects[0].sigma === 2e7);
    check('sigma on a dielectric is an error', parseAndEvaluate('diel x=0 y=0 w=1 h=1 er=2 sigma=1e7').errors.length === 1);
    check('sigma must be positive', parseAndEvaluate('sig+ x=0 y=0 w=1 h=1 sigma=0').errors.length === 1);
    const s = build('units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4\nsig+ x=-0.1 w=0.2 y=0.2 h=0.035 sigma=1e7\n');
    const back = build(solverToGeometryText(s));
    check('sigma survives the conversion to text', back.conductors.find(c => c.is_signal).sigma === 1e7);
}

// Microstrip on a finite ground, open on every side: both conductors are rectangles.
const HEAD = 'units mm\nbounds open open open open\ndiel x=-1.5 w=3 y=0 h=0.2 er=4.3 tand=0.02\n';
const line = (sig, gnd) => HEAD + `gnd x=-1 w=2 y=-0.035 h=0.035 ${gnd}\nsig+ x=-0.15 w=0.3 y=0.2 h=0.035 ${sig}\n`;
const PHEAD = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4.3 tand=0.02\n';
const pair = (l, r) => PHEAD + `sig- x=-0.4 w=0.3 y=0.2 h=0.035 ${l}\nsig+ x=0.1 w=0.3 y=0.2 h=0.035 ${r}\n`;

for (const [name, extra, tolSame] of [['QS', {}, 1e-9], ['full-wave', TRI, 0.01]]) {
    // --- 1 ---
    const cu = (await solveR(line('', ''), extra))[0];
    const own = (await solveR(line(`sigma=${CU}`, `sigma=${CU}`), extra))[0];
    check(`${name}: sigma equal to the solver-wide one is a no-op`, cu.R === own.R && cu.L === own.L, `R ${cu.R.toFixed(3)}`);
    if (extra.mesh_backend) check(`${name}: loss comes from the eddy-current solve`, cu.via === 'mqs', cu.via);

    // --- 2 ---
    const low = (await solveR(line('', ''), { ...extra, sigma_cond: LOW }))[0];
    const lowOwn = (await solveR(line(`sigma=${LOW}`, `sigma=${LOW}`), extra))[0];
    check(`${name}: sigma on every conductor = the solver-wide sigma`, rel(low.R, lowOwn.R) < tolSame && rel(low.L, lowOwn.L) < tolSame,
        `R ${lowOwn.R.toFixed(3)} vs ${low.R.toFixed(3)}, rel ${rel(low.R, lowOwn.R).toExponential(1)}`);

    // --- 3 ---
    const sigOnly = (await solveR(line(`sigma=${LOW}`, ''), extra))[0];
    const gndOnly = (await solveR(line('', `sigma=${LOW}`), extra))[0];
    check(`${name}: a poor trace or a poor ground lands between the uniform lines`,
        sigOnly.R > cu.R * 1.1 && sigOnly.R < low.R * 0.99 && gndOnly.R > cu.R * 1.02 && gndOnly.R < sigOnly.R,
        `${cu.R.toFixed(2)} < gnd ${gndOnly.R.toFixed(2)} < trace ${sigOnly.R.toFixed(2)} < ${low.R.toFixed(2)} ohm/m`);
    // Trace and ground losses add, so the two partial changes make up the full one.
    check(`${name}: trace and ground increments add up`, rel(sigOnly.R + gndOnly.R - cu.R, low.R) < 0.03,
        `${(sigOnly.R + gndOnly.R - cu.R).toFixed(2)} vs ${low.R.toFixed(2)}`);

    // --- 4 ---
    const dc = 1 / (LOW * 0.3e-3 * 35e-6) + 1 / (CU * 2e-3 * 35e-6);
    const lf = (await solveR(line(`sigma=${LOW}`, ''), { ...extra, freq: 1e4 }))[0];
    check(`${name}: low-frequency R is the series DC resistance`, rel(lf.R, dc) < 0.03, `${lf.R.toFixed(3)} vs ${dc.toFixed(3)} ohm/m`);

    // --- 4b ---
    // Plating several skin depths thick hides the bulk: the same line as a trace made
    // of the plating metal (delta = 1.2 um at 5 GHz and 3.8e7 S/m).
    const asMetal = (await solveR(line('sigma=3.8e7', ''), extra))[0];
    const asPlating = (await solveR(line('plating=top,sides,bottom plating_sigma=3.8e7 plating_t=0.01', ''), extra))[0];
    check(`${name}: thick plating reads as a trace of the plating metal`,
        rel(asMetal.R, asPlating.R) < 2e-3 && rel(asMetal.L, asPlating.L) < 1e-4 && asMetal.R > cu.R * 1.05,
        `R ${asPlating.R.toFixed(3)} vs ${asMetal.R.toFixed(3)}, bare copper ${cu.R.toFixed(3)}`);

    // --- 5 ---
    const A = `sigma=${LOW}`, B = '';
    const rA = await solveR(pair(A, A), extra), rB = await solveR(pair(B, B), extra);
    const rM = await solveR(pair(A, B), extra), rS = await solveR(pair(B, A), extra);
    for (const [i, mode] of [[0, 'odd'], [1, 'even']]) {
        const mean = 0.5 * (rA[i].R + rB[i].R);
        check(`${name}: mixed pair ${mode}-mode R is the mean of the two uniform pairs`,
            rel(rM[i].R, mean) < 0.03 && Math.abs(rA[i].R - rB[i].R) / mean > 0.2,
            `${rM[i].R.toFixed(2)} vs mean ${mean.toFixed(2)} of ${rA[i].R.toFixed(2)} / ${rB[i].R.toFixed(2)}`);
        check(`${name}: swapping the two metals gives the same ${mode}-mode R`, rel(rM[i].R, rS[i].R) < 0.01,
            `${rM[i].R.toFixed(3)} vs ${rS[i].R.toFixed(3)}`);
    }
}
{
    // The perturbation loss (the fallback when the eddy-current solve cannot run)
    // evaluates each conductor's integral at its own sigma.
    const pert = async (text, extra = {}) => {
        const s = build(text, { ...TRI, ...extra });
        s.tri_opts = { lossMethod: 'perturbation' };
        return (await quiet(() => s.solve_adaptive(SOLVE))).modes[0].RLGC.R;
    };
    const low = await pert(line('', ''), { sigma_cond: LOW });
    const lowOwn = await pert(line(`sigma=${LOW}`, `sigma=${LOW}`));
    check('full-wave perturbation: sigma on every conductor = the solver-wide sigma', rel(low, lowOwn) < 1e-9,
        `${lowOwn.toFixed(3)} vs ${low.toFixed(3)}`);
}
{
    const mixed = build(pair(`sigma=${LOW}`, ''));
    check('traces of different metals solve on the full domain', mixed.sym_half === false
        && build(pair(`sigma=${LOW}`, `sigma=${LOW}`)).sym_half === true);
}

console.log(failures === 0 ? '\nALL CUSTOM GEOMETRY SIGMA TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
