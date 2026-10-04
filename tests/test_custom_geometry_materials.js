// Per-conductor metal and finish in custom geometry, on both backends: a conductor's own
// conductivity (sigma=), roughness (rq=) and plating material (plating_sigma=,
// plating_t=, plating_rq=).
//
//   1. a value equal to the solver-wide one changes nothing
//   2. sigma on every conductor = the solver-wide sigma set to that value
//   3. a poor trace or a rough trace only lands between the two uniform lines
//   4. low frequency: R tends to the series DC resistance of the two metals, and
//      plating several skin depths thick equals a trace made of the plating metal
//   5. a pair with a different metal or finish on its two traces leaves the half domain
//      (the fields mirror, the loss does not) and reads the mean of the two uniform
//      pairs, whichever trace carries which
//   6. full-wave, a rough ground reaching the open domain edge (the ideal-ground path
//      at low frequency): sigma on every conductor = the solver-wide sigma, L included
//   7. one net of two separate traces, one of a poor metal in its skin transition: QS
//      blends each trace on its own and agrees with full-wave
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { parseAndEvaluate, solverToGeometryText } from '../src/custom_geometry_text.js';
import { check, quiet, rel, APP, done } from './helpers.js';

const SOLVE = { ...APP, max_iters: 8 };
const build = (text, extra = {}) => new CustomGeometrySolver({ text, nx: 30, ny: 30, freq: 5e9, ...extra });
async function solveR(text, extra = {}) {
    const r = await quiet(() => build(text, extra).solve_adaptive(SOLVE));
    return r.modes.map(m => ({ R: m.RLGC.R, L: m.RLGC.L, via: m.lossVia }));
}
const BACKENDS = [['QS', {}], ['full-wave', { mesh_backend: 'triangular' }]];
const CU = 5.8e7, LOW = 1e7;

// Microstrip on a finite ground, open on every side: both conductors are rectangles.
const OPEN = 'units mm\nbounds open open open open\ndiel x=-1.5 w=3 y=0 h=0.2 er=4.3 tand=0.02\n';
const line = (sig, gnd) => OPEN + `gnd x=-1 w=2 y=-0.035 h=0.035 ${gnd}\nsig+ x=-0.15 w=0.3 y=0.2 h=0.035 ${sig}\n`;
const HEAD = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4.3 tand=0.02\n';
const single = extra => HEAD + `sig+ x=-0.15 w=0.3 y=0.2 h=0.035 ${extra}\n`;
const pair = (l, r) => HEAD + `sig- x=-0.4 w=0.3 y=0.2 h=0.035 ${l}\nsig+ x=0.1 w=0.3 y=0.2 h=0.035 ${r}\n`;

// Mixed pair A|B against the uniform pairs A|A and B|B, and swapped B|A.
async function mixedPair(name, what, A, B, extra, tol, spread) {
    const rA = await solveR(pair(A, A), extra), rB = await solveR(pair(B, B), extra);
    const rM = await solveR(pair(A, B), extra), rS = await solveR(pair(B, A), extra);
    for (const [i, mode] of [[0, 'odd'], [1, 'even']]) {
        const mean = 0.5 * (rA[i].R + rB[i].R);
        check(`${name}: ${what} pair ${mode}-mode R is the mean of the two uniform pairs`,
            rel(rM[i].R, mean) < tol && Math.abs(rA[i].R - rB[i].R) / mean > spread,
            `${rM[i].R.toFixed(2)} vs mean ${mean.toFixed(2)} of ${rA[i].R.toFixed(2)} / ${rB[i].R.toFixed(2)}`);
        check(`${name}: swapping the two traces of the ${what} pair gives the same ${mode}-mode R`, rel(rM[i].R, rS[i].R) < 0.01,
            `${rM[i].R.toFixed(3)} vs ${rS[i].R.toFixed(3)}`);
    }
}

// --- Text keys ---
{
    const g = parseAndEvaluate('units um\ns = 1e7\nsig+ x=0 y=0 w=100 h=35 sigma=s*2');
    check('sigma evaluates as a plain number', g.errors.length === 0 && g.rects[0].sigma === 2e7);
    const f = parseAndEvaluate('units um\nsig+ x=0 y=0 w=100 h=35 rq=0.5 plating=top plating_sigma=1e7 plating_t=4 plating_rq=0.2');
    const r = f.rects[0];
    check('rq and plating keys evaluate to metres', f.errors.length === 0 && Math.abs(r.rq - 0.5e-6) < 1e-15
        && r.platingMaterial.sigma === 1e7 && Math.abs(r.platingMaterial.thickness - 4e-6) < 1e-15
        && Math.abs(r.platingMaterial.rq - 0.2e-6) < 1e-15);
    for (const [what, text] of [['negative sigma on a dielectric', 'diel x=0 y=0 w=1 h=1 er=2 sigma=-1'],
        ['zero sigma', 'sig+ x=0 y=0 w=1 h=1 sigma=0'], ['rq on a dielectric', 'diel x=0 y=0 w=1 h=1 er=2 rq=1'],
        ['negative rq', 'sig+ x=0 y=0 w=1 h=1 rq=-1']]) {
        check(`${what} is an error`, parseAndEvaluate(text).errors.length === 1);
    }
    const s = build(HEAD + 'sig+ x=-0.1 w=0.2 y=0.2 h=0.035 sigma=1e7\n');
    check('sigma survives the conversion to text', build(solverToGeometryText(s)).conductors.find(c => c.is_signal).sigma === 1e7);
    const stmt = build(HEAD + 'plating sigma=1.45e7 t=0.004 rq=0.0002\n' + 'sig+ x=-0.15 w=0.3 y=0.2 h=0.035 plating=top,sides\n');
    const own = build(single('plating=top,sides plating_sigma=1.45e7 plating_t=0.004 plating_rq=0.0002'));
    const key = c => JSON.stringify(Object.entries(c.plating || {}).sort());
    check('plating keys on the line build the same plating as the statement',
        key(stmt.conductors.find(c => c.is_signal)) === key(own.conductors.find(c => c.is_signal)));
    let msg = '';
    try { build(single('plating=top')); } catch (e) { msg = e.message; }
    check('plating without any material is rejected', /plating material/.test(msg), msg);
}

for (const [name, extra] of BACKENDS) {
    const tolSame = extra.mesh_backend ? 0.01 : 1e-9;
    // --- 1 ---
    const cu = (await solveR(line('', ''), extra))[0];
    const own = (await solveR(line(`sigma=${CU}`, `sigma=${CU}`), extra))[0];
    check(`${name}: sigma equal to the solver-wide one is a no-op`, cu.R === own.R && cu.L === own.L, `R ${cu.R.toFixed(3)}`);
    if (extra.mesh_backend) check(`${name}: loss comes from the eddy-current solve`, cu.via === 'mqs', cu.via);
    const a = await solveR(single(''), { ...extra, rq: 0.5e-6 });
    const b = await solveR(single('rq=0.0005'), { ...extra, rq: 0.5e-6 });
    check(`${name}: rq equal to the solver-wide roughness is a no-op`, a[0].R === b[0].R && a[0].L === b[0].L, `R ${a[0].R.toFixed(3)}`);

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
    const smooth = (await solveR(single(''), extra))[0].R;
    const rough = (await solveR(single(''), { ...extra, rq: 1e-6 }))[0].R;
    const roughSig = (await solveR(single('rq=0.001'), extra))[0].R;
    check(`${name}: a rough trace over a smooth ground lands between smooth and rough`,
        roughSig > smooth * 1.05 && roughSig < rough * 0.99,
        `${smooth.toFixed(2)} < ${roughSig.toFixed(2)} < ${rough.toFixed(2)} ohm/m`);

    // --- 4 ---
    const dc = 1 / (LOW * 0.3e-3 * 35e-6) + 1 / (CU * 2e-3 * 35e-6);
    const lf = (await solveR(line(`sigma=${LOW}`, ''), { ...extra, freq: 1e4 }))[0];
    check(`${name}: low-frequency R is the series DC resistance`, rel(lf.R, dc) < 0.03, `${lf.R.toFixed(3)} vs ${dc.toFixed(3)} ohm/m`);
    // Plating several skin depths thick hides the bulk: the same line as a trace made
    // of the plating metal (delta = 1.2 um at 5 GHz and 3.8e7 S/m).
    const asMetal = (await solveR(line('sigma=3.8e7', ''), extra))[0];
    const asPlating = (await solveR(line('plating=top,sides,bottom plating_sigma=3.8e7 plating_t=0.01', ''), extra))[0];
    check(`${name}: thick plating reads as a trace of the plating metal`,
        rel(asMetal.R, asPlating.R) < 2e-3 && rel(asMetal.L, asPlating.L) < 1e-4 && asMetal.R > cu.R * 1.05,
        `R ${asPlating.R.toFixed(3)} vs ${asMetal.R.toFixed(3)}, bare copper ${cu.R.toFixed(3)}`);

    // --- 5 ---
    await mixedPair(name, 'two-metal', `sigma=${LOW}`, '', extra, 0.03, 0.2);
    // Solid plating on one trace is a trace of the plating metal, also below the skin
    // transition where the surface impedance cannot stand in for it.
    const solidOne = await solveR(pair('', 'plating=top,sides,bottom plating_sigma=1e7 plating_t=0.04'), { ...extra, freq: 1e6 });
    const ownOne = await solveR(pair('', 'sigma=1e7'), { ...extra, freq: 1e6 });
    check(`${name}: solid plating on one trace of a pair = that trace in the plating metal`,
        [0, 1].every(i => rel(solidOne[i].R, ownOne[i].R) < 0.01),
        solidOne.map((m, i) => `${m.R.toFixed(3)} vs ${ownOne[i].R.toFixed(3)}`).join(', '));
    // A poorly conducting plating on one trace, a rough bare surface on the other.
    await mixedPair(name, 'two-finish', 'plating=top,sides,bottom plating_sigma=2e6 plating_t=0.004', 'rq=0.001',
        extra, extra.mesh_backend ? 0.04 : 0.03, 0.05);
}
{
    // The perturbation loss (the fallback when the eddy-current solve cannot run)
    // evaluates each conductor's integral at its own sigma.
    const pert = async (text, extra = {}) => {
        const s = build(text, { mesh_backend: 'triangular', ...extra });
        s.tri_opts = { lossMethod: 'perturbation' };
        return (await quiet(() => s.solve_adaptive(SOLVE))).modes[0].RLGC.R;
    };
    const low = await pert(line('', ''), { sigma_cond: LOW });
    const lowOwn = await pert(line(`sigma=${LOW}`, `sigma=${LOW}`));
    check('full-wave perturbation: sigma on every conductor = the solver-wide sigma', rel(low, lowOwn) < 1e-9,
        `${lowOwn.toFixed(3)} vs ${low.toFixed(3)}`);
}
{
    check('traces of different metals solve on the full domain', build(pair(`sigma=${LOW}`, '')).sym_half === false
        && build(pair(`sigma=${LOW}`, `sigma=${LOW}`)).sym_half === true);
    const A = 'plating=top,sides,bottom plating_sigma=2e6 plating_t=0.004', B = 'rq=0.001';
    const same = build(pair(A, A));
    check('identical finishes keep the half domain', same.sym_half === true && same.tri_symmetry !== false);
    check('different finishes solve on the full domain', build(pair(A, B)).sym_half === false);
    const triMixed = build(pair(A, B), { mesh_backend: 'triangular' });
    await quiet(() => triMixed.solve_adaptive(SOLVE));
    check('full-wave: different finishes solve on the full domain', triMixed._triBackend.symmetry === false);
}

// --- 6. ideal ground of its own metal ---
{
    const geo = own => 'units mm\nbounds open open open open\ndiel x=-inf w=inf y=0 h=0.2 er=4\n'
        + `sig+ x=-0.15 w=0.3 y=0.2 h=0.035${own}\ngnd x=-inf w=inf y=-0.035 h=0.035 rq=2um${own}\n`;
    const run = async (own, sc) => {
        const s = new CustomGeometrySolver({ text: geo(own), nx: 30, ny: 30, freq: 3e6, sigma_cond: sc, mesh_backend: 'triangular' });
        s.use_causal_materials = false;
        return (await quiet(() => s.solve_adaptive(APP))).modes[0];
    };
    const a = await run('', 1e7), b = await run(' sigma=1e7', 5.8e7);
    check('full-wave ideal rough ground of its own sigma: R and L match the solver-wide sigma',
        rel(a.RLGC.R, b.RLGC.R) < 1e-6 && rel(a.L_internal, b.L_internal) < 1e-6,
        `R ${a.RLGC.R.toFixed(4)} / ${b.RLGC.R.toFixed(4)}, Lint ${(a.L_internal * 1e9).toFixed(3)} / ${(b.L_internal * 1e9).toFixed(3)} nH/m`);
}

// --- 7. a net of two metals in the skin transition ---
{
    // sigma = 1e6 at 1 GHz: delta = 16 um against t = 35 um, in the transition notch.
    const text = HEAD + 'sig+ x=-0.8 w=0.2 y=0.2 h=0.035 sigma=1e6\nsig+ x=0.6 w=0.2 y=0.2 h=0.035\n';
    const [q] = await solveR(text, { freq: 1e9 }), [f] = await solveR(text, { freq: 1e9, mesh_backend: 'triangular' });
    check('one net of a copper and a poor trace: QS R agrees with full-wave', rel(q.R, f.R) < 0.05,
        `QS ${q.R.toFixed(2)} vs full-wave ${f.R.toFixed(2)} ohm/m`);
}

done();
