// Custom geometries that are not mirror symmetric, on both backends:
//   1. QS conductor loss on an asymmetric grid takes the centred quadrature: R of a
//      geometry and its mirror image agree and match the full-wave R. The one-sided
//      rule read 286 / 280 ohm/m against full-wave 292.
//      The native offset broadside stripline takes the same rule.
//   2. dielectric paint order: mirrored overlapping dielectrics painted in mismatched
//      order are not symmetric, for the swap symmetry of a pair (which used to keep the
//      odd/even decomposition with C11 = C22 forced) as well as for the half domain,
//      on both mirror axes. Overlapping metals of a conductor follow the same rule.
//
// Run: node tests/test_custom_geometry_asymmetry.js
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { conductorSwapSymmetric } from '../src/geometry_symmetry.js';
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

done();
