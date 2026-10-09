// Custom geometries that are not mirror symmetric, on both backends:
//   1. QS conductor loss on an asymmetric grid takes the centred quadrature: R of a
//      geometry and its mirror image agree and match the full-wave R. The one-sided
//      rule read 286 / 280 ohm/m against full-wave 292.
//      The native offset broadside stripline takes the same rule.
//   2. dielectric paint order: mirrored overlapping dielectrics painted in mismatched
//      order are not symmetric, for the swap symmetry of a pair (which used to keep the
//      odd/even decomposition with C11 = C22 forced) as well as for the half domain,
//      on both mirror axes. Overlapping metals of a conductor follow the same rule.
//   3. a trace off centre by rounding of its typed numbers (x=-199.29 w=398.57 ends at
//      199.28) is not symmetric: QS solves it on the full domain and agrees with the
//      exactly centred trace. Its grid was once declared symmetric by the geometry test
//      and rejected by the mesher's stricter one, so the solve threw.
//   4. a layer whose rounded y + h misses the trace face by 3 nm is snapped onto it: the
//      film between the two lines was meshed by full-wave as slivers, which put C0 and
//      C 3.5x off (eps far above every eps_r) and the ground surface loss 2x off.
//   5. an asymmetric pair of unequal metal: the QS mode R is 1/2 i^T R i of its own
//      per-line matrix with the modal currents i. The skin-transition blend weighted the
//      DC floor of each trace by the +-1 currents of a symmetric pair (odd R 3x high).
//   6. decoupled lines (a ground plane between two stacked traces): [C0]^-1 [C] is
//      diagonal or nearly so, and its eigenvectors are the two traces. They were built
//      from the off-diagonal entry, which gave QS the first trace for both modes (exactly
//      0) and full-wave rounding noise (1e-21).
//
// Run: node tests/test_custom_geometry_asymmetry.js
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { conductorSwapSymmetric } from '../src/geometry_symmetry.js';
import { eig2x2 } from '../src/matrix.js';
import { check, quiet, rel, APP, done } from './helpers.js';

const TRI = { mesh_backend: 'triangular' };
const pct = v => `${(100 * v).toFixed(2)}%`;
async function solve(text, extra = {}, opts = APP) {
    const s = new CustomGeometrySolver({ text, nx: 30, ny: 30, freq: 1e9, ...extra });
    const r = await quiet(() => s.solve_adaptive(opts));
    return { s, r };
}

// --- 1. Loss quadrature ---
{
    // 20 um trace with a coplanar ground on one side, and its mirror image.
    const g = m => `units um\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=100 er=4\n` +
        `sig+ x=${-10 * m} w=${20 * m} y=100 h=5\ngnd x=${60 * m} w=${40 * m} y=100 h=5\n`;
    const sym = 'units um\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=100 er=4\nsig+ x=-10 w=20 y=100 h=5\n';
    // The centred rule follows from the asymmetric grid, the symmetric one keeps the default.
    const sa = new CustomGeometrySolver({ text: g(1) }), ss = new CustomGeometrySolver({ text: sym });
    check('an asymmetric geometry meshes an asymmetric grid, a symmetric one keeps the default rule',
        sa.mesher.symmetric === false && ss.mesher.symmetric === true && !ss.centred_loss_quadrature);
    const opts = { ...APP, energy_tol: 0.002, param_tol: 0.01, max_nodes: 80000, max_iters: 16 };
    const a = (await solve(g(1), {}, opts)).r.modes[0].RLGC.R;
    const b = (await solve(g(-1), {}, opts)).r.modes[0].RLGC.R;
    const fw = (await solve(g(1), TRI)).r.modes[0].RLGC.R;
    check('QS: R of an asymmetric geometry and of its mirror image agree', rel(a, b) < 0.005,
        `${a.toFixed(1)} / ${b.toFixed(1)} ohm/m`);
    check('QS: R of an asymmetric geometry matches full-wave', rel(a, fw) < 0.025,
        `QS ${a.toFixed(1)}, full-wave ${fw.toFixed(1)} ohm/m`);
}

// --- 2. Paint order ---
{
    // Left trace on er = 3, right trace on er = 2: the later rectangles override the
    // earlier ones, but the rectangle set is mirror symmetric.
    const order = (o) => `units um\nbounds open open open gnd\n` + o.map(([x, er]) =>
        `diel x=${x} w=60 y=0 h=100 er=${er}\n`).join('') +
        'sig- x=-30 w=20 y=100 h=5\nsig+ x=10 w=20 y=100 h=5\n';
    const bad = order([[-50, 2], [-10, 3], [-10, 2], [-50, 3]]);
    const good = order([[-50, 2], [-10, 2], [-50, 3], [-10, 3]]);
    const sBad = new CustomGeometrySolver({ text: bad }), sGood = new CustomGeometrySolver({ text: good });
    check('pair swap symmetry: a mismatched paint order is asymmetric, a mirrored one symmetric',
        conductorSwapSymmetric(sBad.conductors, sBad.dielectrics) === false
        && conductorSwapSymmetric(sGood.conductors, sGood.dielectrics) === true);

    // Broadside pair, the mirror plane is y = 50: the same order test on the y axis.
    const bs = o => 'units um\nbounds open open gnd gnd\n' + o.map(([y, er]) =>
        `diel x=-inf w=inf y=${y} h=60 er=${er}\n`).join('') +
        'sig+ x=-10 w=20 y=20 h=5\nsig- x=-10 w=20 y=75 h=5\n';
    const bsBad = new CustomGeometrySolver({ text: bs([[0, 2], [40, 3], [40, 2], [0, 3]]) });
    const bsGood = new CustomGeometrySolver({ text: bs([[0, 2], [40, 2], [0, 3], [40, 3]]) });
    check('broadside swap symmetry follows the paint order too',
        conductorSwapSymmetric(bsBad.conductors, bsBad.dielectrics) === false
        && conductorSwapSymmetric(bsGood.conductors, bsGood.dielectrics) === true);

    // Overlapping metals of one signal: the later one fills the overlap, so their mirror
    // images must overlap in the same order for the half domain.
    const metal = o => 'units um\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=100 er=4\n' +
        o.map(([x, sg]) => `sig+ x=${x} w=25 y=100 h=5${sg ? ' sigma=1e7' : ''}\n`).join('');
    const mBad = new CustomGeometrySolver({ text: metal([[-20, true], [-5, false], [-5, true], [-20, false]]), mesh_backend: 'triangular' });
    const mGood = new CustomGeometrySolver({ text: metal([[-20, true], [-5, true], [-20, false], [-5, false]]), mesh_backend: 'triangular' });
    check('overlapping metals painted in mismatched order solve on the full domain',
        mBad.tri_symmetry === false && mGood.tri_symmetry !== false);

    for (const [name, extra] of [['QS', {}], ['full-wave', TRI]]) {
        const { r } = await solve(bad, extra);
        const C = r.RLGC_matrix && r.RLGC_matrix.C;
        check(`${name}: a pair painted asymmetrically keeps its physical matrices`,
            !!r.physMatrix && !!C && rel(C[0][0], C[1][1]) > 0.05,
            C ? `C11 / C22 = ${(C[0][0] / C[1][1]).toFixed(3)}` : 'no RLGC_matrix');
    }
    const q = (await solve(bad)).r.modes, f = (await solve(bad, TRI)).r.modes;
    const d = Math.max(...q.map((m, i) => rel(m.eps_eff, f[i].eps_eff)));
    check('pair painted asymmetrically: the modes of the two backends agree', d < 0.02, `eps_eff ${pct(d)}`);
}

// --- 3. Nearly symmetric trace ---
{
    const g = x => 'units um\nbounds open open open open\ndiel x=-inf w=inf y=0 h=550 er=6.3\n' +
        `sig+ x=${x} w=398.57 y=550 h=50\ngnd x=-970 w=1940 y=-36 h=36\n`;
    let err = null, near = null;
    try { near = await solve(g(-199.29)); } catch (e) { err = e.message; }
    check('QS: a trace off centre by 0.01 um solves', !err, err || '');
    if (near) {
        const exact = await solve(g(-398.57 / 2));
        check('QS: it solves on the full domain', near.s.sym_half === false && exact.s.sym_half === true);
        const d = rel(near.r.modes[0].Z0, exact.r.modes[0].Z0);
        check('QS: and agrees with the exactly centred trace', d < 0.002, `Z0 ${pct(d)}`);
    }
}

// --- 4. Dielectric edge a hair off a conductor face ---
{
    const g = h => 'units um\nbounds open open open open\n' +
        `diel x=-inf w=inf y=0 h=228.07 er=8.7\ndiel x=-inf w=inf y=228.07 h=${h} er=6.8\n` +
        'sig+ x=-292.31 w=267.63 y=276.94 h=35.541\nsig- x=24.674 w=267.63 y=276.94 h=35.541\n' +
        'gnd x=-inf w=inf y=-34.497 h=34.497\n';
    const s = new CustomGeometrySolver({ text: g(48.866), nx: 30, ny: 30, freq: 1e8 });
    const top = Math.max(...s.dielectrics.filter(d => d.y_max < 300e-6).map(d => d.y_max));
    check('a layer top 3 nm below the trace face is snapped onto it', top === s.conductors.find(c => c.is_signal).y_min,
        `${(top * 1e6).toFixed(4)} um`);
    const near = (await solve(g(48.866), { ...TRI, freq: 1e8 })).r.modes;
    const exact = (await solve(g(276.94 - 228.07), { ...TRI, freq: 1e8 })).r.modes;
    const d = Math.max(...near.map((m, i) => rel(m.RLGC.C, exact[i].RLGC.C)));
    check('full-wave: the rounded layer gives the C of the exact one', d < 0.01, `C ${pct(d)}`);
}

// --- 5. Modal R of an asymmetric pair of unequal metal ---
{
    const text = 'units um\nbounds open open open open\n' +
        'diel x=-inf w=inf y=0 h=339.38 er=9.0\ndiel x=-inf w=inf y=339.38 h=226.48 er=6.6\n' +
        'diel x=-inf w=inf y=565.86 h=190.74 er=8.7\n' +
        'sig+ x=-131.66 w=51.854 y=565.86 h=13.779\nsig- x=79.803 w=51.854 y=565.86 h=13.779 sigma=1.1244e6\n' +
        'gnd x=-261.71 w=-inf y=565.86 h=13.779\ngnd x=-197.72 w=395.44 y=756.6 h=24.491\n';
    const { s, r } = await solve(text, { freq: 1.38e9, sigma_cond: 3.06e7 });
    for (const m of r.modes) {
        const I = s._modal_currents(m.mode), Rm = m.RLGC.Rm;
        const q = I && Rm ? 0.5 * (I[0] * I[0] * Rm[0] + 2 * I[0] * I[1] * Rm[1] + I[1] * I[1] * Rm[2]) : NaN;
        check(`QS: ${m.mode} R of an unequal-metal pair is 1/2 i^T R i of its per-line matrix`, rel(m.RLGC.R, q) < 0.1,
            `${m.RLGC.R.toFixed(1)} vs ${q.toFixed(1)}`);
    }
}

// --- 6. Eigenvectors of decoupled lines ---
{
    const axes = M => eig2x2(M).vecs.map(v => (Math.abs(v[0]) > Math.abs(v[1]) ? 0 : 1) + (Math.min(Math.abs(v[0]), Math.abs(v[1])) < 1e-6 ? '' : '?'));
    for (const [what, M, want] of [['diagonal', [[3, 0], [0, 2]], '0,1'], ['diagonal, larger second', [[2, 0], [0, 3]], '1,0'],
        ['off-diagonal at rounding level', [[3, 1e-21], [2e-21, 2]], '0,1'], ['equal diagonal', [[3, 0], [0, 3]], '0,1']]) {
        const got = axes(M).join(',');
        check(`eig2x2, ${what}: one eigenvector per trace`, got === want, got);
    }
    const v = eig2x2([[2.5, 0.3], [0.3, 2.5]]).vecs;
    check('eig2x2 keeps the +-1 vectors of a symmetric pair', Math.abs(v[0][0] - 1) < 1e-12 && Math.abs(v[0][1] - 1) < 1e-12
        && Math.abs(v[1][0] - 1) < 1e-12 && Math.abs(v[1][1] + 1) < 1e-12, JSON.stringify(v));
}

done();
