// A conductor drawn as several touching blocks, each with its own metal and finish.
//
//   1. blocks of one metal solve like the single rectangle they add up to
//   2. the same with a conductivity or a roughness of their own on every block
//   3. a thin top block of another metal is a plating: it has to agree with the
//      plating=top model, down to a block thinner than its skin depth. The full-wave
//      eddy-current solve meshes the block and is the reference; the quasi-static
//      surface integral takes the layered impedance for it.
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
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b));

const HEAD = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4.3 tand=0.02\n';
const S = (x, w, y, h, e = '') => `sig+ x=${x} w=${w} y=${y} h=${h} ${e}\n`;
const single = e => S(-0.15, 0.3, 0.2, 0.035, e);
const sideBySide = e => S(-0.15, 0.1, 0.2, 0.035, e) + S(-0.05, 0.2, 0.2, 0.035, e);
const stacked = e => S(-0.15, 0.3, 0.2, 0.0175, e) + S(-0.15, 0.3, 0.2175, 0.0175, e);
const topBlock = (t, e) => S(-0.15, 0.3, 0.2, 0.035 - t) + S(-0.15, 0.3, 0.235 - t, t, e);

const SOLVE = { max_iters: 8, energy_tol: 0.005, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };
async function solve(body, backend) {
    const s = new CustomGeometrySolver({ text: HEAD + body, nx: 30, ny: 30, freq: 5e9, mesh_backend: backend });
    if (backend === 'triangular') s.tri_opts = { lossMethod: 'auto' };
    const m = (await quiet(() => s.solve_adaptive(SOLVE))).modes[0];
    return { R: m.RLGC.R, Li: m.L_internal, Z0: m.Z0 };
}

const ref = {};
// The quasi-static R moves a few percent with the grid lines the extra block edges add.
for (const [name, backend, tol] of [['QS', 'rectilinear', 0.04], ['full-wave', 'triangular', 0.005]]) {
    for (const [what, e] of [['one metal', ''], ['sigma=1e7 on every block', 'sigma=1e7'], ['rq=1um on every block', 'rq=0.001']]) {
        const one = await solve(single(e), backend);
        for (const [shape, body] of [['side by side', sideBySide(e)], ['stacked', stacked(e)]]) {
            const b = await solve(body, backend);
            check(`${name}: ${shape} blocks, ${what} = single rectangle`,
                rel(b.R, one.R) < tol && rel(b.Li, one.Li) < tol && rel(b.Z0, one.Z0) < 0.005,
                `R ${b.R.toFixed(2)} vs ${one.R.toFixed(2)}, L_int ${(b.Li * 1e9).toFixed(3)} vs ${(one.Li * 1e9).toFixed(3)} nH/m`);
        }
    }
    // Top block of 1e7 S/m (skin depth 2.25 um): against the same block in copper, which
    // takes the grid effect out, and against the plating model on the single rectangle.
    const bare = await solve(single(''), backend);
    for (const t of [0.001, 0.005]) {
        const cu = await solve(topBlock(t, ''), backend);
        const blk = await solve(topBlock(t, 'sigma=1e7'), backend);
        const plt = await solve(single(`plating=top plating_sigma=1e7 plating_t=${t}`), backend);
        ref[`${backend}${t}`] = blk;
        // The 5 um block also puts the poor metal on the top corners' sides, which the
        // plating=top model of the full-wave solver leaves as copper.
        const tolP = t < 0.002 ? 0.015 : 0.06;
        check(`${name}: ${t * 1e3} um top block of another metal acts as a plating`,
            rel(blk.R / cu.R, plt.R / bare.R) < tolP && blk.R > cu.R * 1.03,
            `R x${(blk.R / cu.R).toFixed(3)} as a block, x${(plt.R / bare.R).toFixed(3)} as plating`);
    }
}
for (const t of [0.001, 0.005]) {
    const q = ref[`rectilinear${t}`], f = ref[`triangular${t}`];
    check(`${t * 1e3} um top block: quasi-static agrees with full-wave`, rel(q.R, f.R) < 0.05 && rel(q.Li, f.Li) < 0.06,
        `R ${q.R.toFixed(2)} vs ${f.R.toFixed(2)}, L_int ${(q.Li * 1e9).toFixed(3)} vs ${(f.Li * 1e9).toFixed(3)} nH/m`);
}
{
    let msg = '';
    try { new CustomGeometrySolver({ text: HEAD + S(-0.15, 0.15, 0.2, 0.035, 'plating=top plating_sigma=1e7 plating_t=0.002') + S(0, 0.15, 0.2, 0.035) }); }
    catch (e) { msg = e.message; }
    check('plating= on a block that touches another block of its conductor is rejected', /plating/i.test(msg), msg);
}

console.log(failures === 0 ? '\nALL CUSTOM GEOMETRY BLOCK TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
