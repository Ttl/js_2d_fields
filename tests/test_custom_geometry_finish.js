// Per-conductor surface finish in custom geometry: a conductor's own roughness (rq=) and
// its own plating material (plating_sigma=, plating_t=, plating_rq=), on both backends.
//
//   1. a value equal to the solver-wide one changes nothing
//   2. roughness on one conductor lands between all-smooth and all-rough
//   3. a pair with different finishes on its two traces leaves the half domain (the
//      fields mirror, the loss does not) and reads the mean of the two uniform pairs
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { parseAndEvaluate } from '../src/custom_geometry_text.js';

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
    return r.modes.map(m => ({ R: m.RLGC.R, L: m.RLGC.L, Z0: m.Z0?.re ?? m.Z0 }));
}
const TRI = { mesh_backend: 'triangular' };

// --- Text keys ---
{
    const g = parseAndEvaluate('units um\nsig+ x=0 y=0 w=100 h=35 rq=0.5 plating=top plating_sigma=1e7 plating_t=4 plating_rq=0.2');
    const r = g.rects[0];
    check('rq and plating keys evaluate to metres', g.errors.length === 0 && Math.abs(r.rq - 0.5e-6) < 1e-15
        && r.platingMaterial.sigma === 1e7 && Math.abs(r.platingMaterial.thickness - 4e-6) < 1e-15
        && Math.abs(r.platingMaterial.rq - 0.2e-6) < 1e-15);
    check('rq on a dielectric is an error', parseAndEvaluate('diel x=0 y=0 w=1 h=1 er=2 rq=1').errors.length === 1);
    check('negative rq is an error', parseAndEvaluate('sig+ x=0 y=0 w=1 h=1 rq=-1').errors.length === 1);
}

const HEAD = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4.3 tand=0.02\n';
const single = (extra) => HEAD + `sig+ x=-0.15 w=0.3 y=0.2 h=0.035 ${extra}\n`;
const pair = (l, r) => HEAD + `sig- x=-0.4 w=0.3 y=0.2 h=0.035 ${l}\nsig+ x=0.1 w=0.3 y=0.2 h=0.035 ${r}\n`;

// --- 1. Own value equal to the solver-wide one ---
for (const [name, extra] of [['QS', {}], ['full-wave', TRI]]) {
    const a = await solveR(single(''), { ...extra, rq: 0.5e-6 });
    const b = await solveR(single('rq=0.0005'), { ...extra, rq: 0.5e-6 });
    check(`${name}: rq equal to the solver-wide roughness is a no-op`, a[0].R === b[0].R && a[0].L === b[0].L, `R ${a[0].R.toFixed(3)}`);
}
{
    const stmt = build(HEAD + 'plating sigma=1.45e7 t=0.004 rq=0.0002\n' + 'sig+ x=-0.15 w=0.3 y=0.2 h=0.035 plating=top,sides\n');
    const own = build(single('plating=top,sides plating_sigma=1.45e7 plating_t=0.004 plating_rq=0.0002'));
    const key = c => JSON.stringify(Object.entries(c.plating || {}).sort());
    check('plating keys on the line build the same plating as the statement',
        key(stmt.conductors.find(c => c.is_signal)) === key(own.conductors.find(c => c.is_signal)));
    let msg = '';
    try { build(single('plating=top')); } catch (e) { msg = e.message; }
    check('plating without any material is rejected', /plating material/.test(msg), msg);
}

// --- 2. Roughness on the signal only ---
for (const [name, extra, tol] of [['QS', {}, 0], ['full-wave', TRI, 0]]) {
    const smooth = (await solveR(single(''), extra))[0].R;
    const rough = (await solveR(single(''), { ...extra, rq: 1e-6 }))[0].R;
    const sigOnly = (await solveR(single('rq=0.001'), extra))[0].R;
    check(`${name}: a rough trace over a smooth ground lands between smooth and rough`,
        sigOnly > smooth * 1.05 && sigOnly < rough * (1 - 0.01 + tol),
        `${smooth.toFixed(2)} < ${sigOnly.toFixed(2)} < ${rough.toFixed(2)} ohm/m`);
}

// --- 3. Differential pair with different finishes ---
{
    // A poorly conducting plating on one trace, a rough bare surface on the other.
    const A = 'plating=top,sides,bottom plating_sigma=2e6 plating_t=0.004', B = 'rq=0.001';
    const same = build(pair(A, A)), mixed = build(pair(A, B)), swapped = build(pair(B, A));
    check('identical finishes keep the half domain', same.sym_half === true && same.tri_symmetry !== false);
    check('different finishes solve on the full domain', mixed.sym_half === false);
    const triMixed = build(pair(A, B), TRI);
    await quiet(() => triMixed.solve_adaptive(SOLVE));
    check('full-wave: different finishes solve on the full domain', triMixed._triBackend.symmetry === false);

    for (const [name, extra, tol] of [['QS', {}, 0.03], ['full-wave', TRI, 0.04]]) {
        const rA = await solveR(pair(A, A), extra), rB = await solveR(pair(B, B), extra);
        const rM = await solveR(pair(A, B), extra), rS = await solveR(pair(B, A), extra);
        for (const [i, mode] of [[0, 'odd'], [1, 'even']]) {
            const mean = 0.5 * (rA[i].R + rB[i].R);
            check(`${name}: mixed pair ${mode}-mode R is the mean of the two uniform pairs`,
                rel(rM[i].R, mean) < tol && Math.abs(rA[i].R - rB[i].R) / mean > 0.05,
                `${rM[i].R.toFixed(2)} vs mean ${mean.toFixed(2)} of ${rA[i].R.toFixed(2)} / ${rB[i].R.toFixed(2)}`);
            check(`${name}: swapping the two finishes gives the same ${mode}-mode R`, rel(rM[i].R, rS[i].R) < 0.01,
                `${rM[i].R.toFixed(3)} vs ${rS[i].R.toFixed(3)}`);
        }
    }
}

console.log(failures === 0 ? '\nALL CUSTOM GEOMETRY FINISH TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
