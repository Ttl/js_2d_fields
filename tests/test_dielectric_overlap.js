// Overlapping dielectrics: the later rectangle wins, on both backends.
//
// Each geometry is given twice, once as an overlay (a base rectangle with a later
// rectangle of another material on top of part of it) and once tiled from
// non-overlapping rectangles that describe the same materials. The tiles put their
// edges where the overlay has them, so the FDM grid and therefore the FDM results are
// bit-identical where the mesher lays out the same grid. The triangular mesh differs between the two (the tiles add internal
// edges), so there the material map of the overlay mesh is checked triangle by triangle
// against the tiled description, and the results to mesh-noise level.
//
// Reversing the order hides the inset under the base, which must read as the plain
// substrate. Runs with causal materials on, away from the 1 GHz reference, on a fresh
// solver each time (the Djordjevic-Sarkar correction mutates the painted arrays).
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { check, quiet, done } from './helpers.js';


// Coordinates are integer multiples of u = 2^-13 m (0.122 mm), so every edge, sum and
// difference is an exact double and the two descriptions share their grid lines exactly.
const U = 1 / 8192;
const HEAD = 'units m\nu = 1/8192\nbounds open open open gnd\ndomain -32*u 32*u 0 32*u\n';
const TRACE = 'sig+ x=-2*u w=4*u y=4*u h=u/4\n';
const BASE = 'diel x=-inf w=inf y=0 h=4*u er=4 tand=0.02\n';
const PLAIN = (x, y) => (y < 4 * U ? [4, 0.02] : [1, 0]);

// Inset inside the substrate, under the trace.
const INSET = {
    overlay: BASE + 'diel x=-4*u w=8*u y=u h=2*u er=9 tand=0.002\n',
    reversed: 'diel x=-4*u w=8*u y=u h=2*u er=9 tand=0.002\n' + BASE,
    tiled: 'diel x=-inf w=inf y=0 h=u er=4 tand=0.02\n' +
           'diel x=-4*u w=-inf y=u h=2*u er=4 tand=0.02\n' +
           'diel x=-4*u w=8*u y=u h=2*u er=9 tand=0.002\n' +
           'diel x=4*u w=inf y=u h=2*u er=4 tand=0.02\n' +
           'diel x=-inf w=inf y=3*u h=u er=4 tand=0.02\n',
    expect: (x, y) => ((Math.abs(x) < 4 * U && y > U && y < 3 * U) ? [9, 0.002] : PLAIN(x, y)),
    expectReversed: PLAIN,
    hidden: true,
};
// Blocks that cross the top of the substrate on both sides of the trace.
const BLOCKS = 'diel x=3*u w=5*u y=2*u h=4*u er=2.5 tand=0.01\ndiel x=-8*u w=5*u y=2*u h=4*u er=2.5 tand=0.01\n';
const inBlock = (x, y) => Math.abs(x) > 3 * U && Math.abs(x) < 8 * U && y > 2 * U && y < 6 * U;
const CROSS = {
    overlay: BASE + BLOCKS,
    reversed: BLOCKS + BASE,
    tiled: 'diel x=-inf w=inf y=0 h=2*u er=4 tand=0.02\n' +
           'diel x=-8*u w=-inf y=2*u h=2*u er=4 tand=0.02\n' +
           'diel x=-3*u w=6*u y=2*u h=2*u er=4 tand=0.02\n' +
           'diel x=8*u w=inf y=2*u h=2*u er=4 tand=0.02\n' + BLOCKS,
    expect: (x, y) => (inBlock(x, y) ? [2.5, 0.01] : PLAIN(x, y)),
    // Reversed, the base hides the part of the blocks inside the substrate.
    expectReversed: (x, y) => (y < 4 * U ? [4, 0.02] : inBlock(x, y) ? [2.5, 0.01] : [1, 0]),
    hidden: false,
};

// Low enough that the full-wave backend has no dispersion to set it apart from the FDM.
const FREQ = 1e8;
const FDM = { max_iters: 4, energy_tol: 0.01, param_tol: 0.05, max_nodes: 6000, min_converged_passes: 2 };
const TRI = { max_iters: 4, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };

function build(diels, extra = {}) {
    const s = new CustomGeometrySolver({ text: HEAD + diels + TRACE, freq: FREQ, nx: 30, ny: 30, ...extra });
    s.use_causal_materials = true;
    return s;
}
const numbers = r => r.modes.map(m => [m.Z0?.re ?? m.Z0, m.eps_eff, m.RLGC.R, m.RLGC.L, m.RLGC.G, m.RLGC.C]).flat();
const relMax = (a, b) => Math.max(...a.map((v, i) => Math.abs(v - b[i]) / (Math.abs(v) + 1e-300)));

// Nominal (not yet causal-corrected) cell materials against the analytic description.
function fdmCellMismatch(s, expect) {
    s.use_causal_materials = false;
    s._setup_geometry();
    const x0 = s.x_shift ?? 0;
    let bad = 0;
    for (let i = 0; i < s.y.length - 1; i++) {
        for (let j = 0; j < s.x.length - 1; j++) {
            const yc = 0.5 * (s.y[i] + s.y[i + 1]);
            if (yc < 0) continue;   // ground slab behind the gnd wall
            const [er, td] = expect(0.5 * (s.x[j] + s.x[j + 1]) + x0, yc);
            if (s.epsilon_cell[i][j] !== er || s.tand_cell[i][j] !== td) bad++;
        }
    }
    s.use_causal_materials = true;
    return bad;
}

const plain = await quiet(async () => numbers(await build(BASE).solve_adaptive(FDM)));

for (const [name, g] of [['inset', INSET], ['crossing block', CROSS]]) {
    // --- FDM ---
    const so = build(g.overlay), st = build(g.tiled);
    check(`${name}: overlay keeps the half domain`, so.sym_half === true && st.sym_half === true);
    const ro = numbers(await quiet(() => so.solve_adaptive(FDM)));
    const rt = numbers(await quiet(() => st.solve_adaptive(FDM)));
    // The mesher grades by the rectangle list, so tiles may add a far-field grid line.
    check(`${name}: FDM overlay = tiled`, g.hidden ? ro.every((v, i) => v === rt[i]) : relMax(ro, rt) < 1e-4,
        `Z0 ${ro[0].toFixed(3)} eps ${ro[1].toFixed(4)} G ${ro[4].toExponential(3)}`);
    check(`${name}: FDM cell materials follow the later rectangle`, fdmCellMismatch(so, g.expect) === 0);
    check(`${name}: overlay changes the answer`, relMax([ro[1]], [plain[1]]) > 5e-3,
        `eps ${ro[1].toFixed(4)} vs plain ${plain[1].toFixed(4)}`);

    const sf = build(g.overlay, { symmetry: false });
    const rf = numbers(await quiet(() => sf.solve_adaptive(FDM)));
    // Adaptive trajectories may part ways, same bands as test_qs_symmetry_half_full.
    const noR = v => [v[0], v[1], v[3], v[4], v[5]];
    check(`${name}: FDM half domain = full domain`, sf.sym_half === false
        && relMax(noR(ro), noR(rf)) < 1e-2 && relMax([ro[2]], [rf[2]]) < 5e-2,
        `max rel ${relMax(noR(ro), noR(rf)).toExponential(2)}, R ${relMax([ro[2]], [rf[2]]).toExponential(2)}`);

    const sr = build(g.reversed);
    const rr = numbers(await quiet(() => sr.solve_adaptive(FDM)));
    check(`${name}: FDM reversed order paints the base last`,
        fdmCellMismatch(sr, g.expectReversed) === 0);
    if (g.hidden) {
        check(`${name}: FDM hidden inset reads as the plain substrate`, relMax([rr[0], rr[1], rr[4]], [plain[0], plain[1], plain[4]]) < 5e-3,
            `max rel ${relMax([rr[0], rr[1], rr[4]], [plain[0], plain[1], plain[4]]).toExponential(2)}`);
    }

    // --- Triangular ---
    const to = build(g.overlay, { mesh_backend: 'triangular' }), tt = build(g.tiled, { mesh_backend: 'triangular' });
    const qo = numbers(await quiet(() => to.solve_adaptive(TRI)));
    const qt = numbers(await quiet(() => tt.solve_adaptive(TRI)));
    const mesh = to._triBackend.mesh;
    let bad = 0, n = 0;
    const x0 = to.x_shift ?? 0;
    for (let t = 0; t < mesh.nTris; t++) {
        const [a, b, c] = [mesh.tris[3 * t], mesh.tris[3 * t + 1], mesh.tris[3 * t + 2]];
        const xc = (mesh.nodes[2 * a] + mesh.nodes[2 * b] + mesh.nodes[2 * c]) / 3 + x0;
        const yc = (mesh.nodes[2 * a + 1] + mesh.nodes[2 * b + 1] + mesh.nodes[2 * c + 1]) / 3;
        // Conductor interiors carry whatever dielectric lies under them.
        if (yc < 0 || (Math.abs(xc) < 2 * U && yc > 4 * U && yc < 4.25 * U)) continue;
        n++;
        // Compare to the nominal value: causal correction only rescales it.
        const er = g.expect(xc, yc)[0];
        const got = mesh.epsMap[t].re;
        if (er === 1 ? Math.abs(got - 1) > 1e-9 : Math.abs(got / er - 1) > 0.1) bad++;
    }
    check(`${name}: tri material map follows the later rectangle`, bad === 0 && n > 100, `${bad} of ${n} triangles`);
    check(`${name}: tri overlay = tiled to mesh noise`, relMax([qo[0], qo[1], qo[4]], [qt[0], qt[1], qt[4]]) < 1e-2,
        `max rel ${relMax([qo[0], qo[1], qo[4]], [qt[0], qt[1], qt[4]]).toExponential(2)}`);
    check(`${name}: tri agrees with FDM`, relMax([qo[0], qo[1]], [ro[0], ro[1]]) < 2e-2,
        `Z0 ${qo[0].toFixed(3)} vs ${ro[0].toFixed(3)}, eps ${qo[1].toFixed(4)} vs ${ro[1].toFixed(4)}`);
    check(`${name}: tri dielectric loss agrees with FDM`, relMax([qo[4]], [ro[4]]) < 5e-2,
        `G ${qo[4].toExponential(3)} vs ${ro[4].toExponential(3)}`);
}

done();
