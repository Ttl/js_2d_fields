// Solder-mask certification at app-default settings (rect backend).
//
// A 20 µm mask adds three thin dielectric regions per axis. Before the mesher
// fixes (thin-sheet bracket lines, per-region allocation floor handling,
// near-duplicate interface merging) those regions consumed the entire per-axis
// line target: the initial mesh started with the mask rasterized ~70% fat and
// the far field starved, and the default-budget certificate failed at 1.28%
// (true C error +0.83%). Post-fix the same case certifies at ~0.5% with true C
// error ~+0.24%. This pins:
//   1. bare default microstrip still certifies (and its mesh is untouched),
//   2. default microstrip + default solder mask certifies at the default
//      budget/tolerance,
//   3. the FDM answer agrees with the triangular backend on C.
import { MicrostripSolver } from '../src/microstrip.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}

const base = {
    trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
    nx: 30, ny: 30, boundaries: ['open', 'open', 'open', 'gnd'],
};
const sm = { use_sm: true, sm_t_sub: 20e-6, sm_t_trace: 20e-6, sm_t_side: 20e-6, sm_er: 3.5, sm_tand: 0.02 };
const APP = { max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2, certify: true };

async function solve(opts, backend = 'rectilinear') {
    const s = new MicrostripSolver({ ...base, ...opts, mesh_backend: backend });
    const log = console.log;
    console.log = () => {};
    try {
        const r = await s.solve_adaptive({ ...APP });
        return { s, m: r.modes[0], cert: s.certification };
    } finally {
        console.log = log;
    }
}

const bare = await solve({});
check('bare microstrip certificate passes at default budget',
    bare.cert && bare.cert.pass, bare.cert ? `err ${(100 * bare.cert.err).toFixed(3)}%` : 'no certificate');

const fdm = await solve(sm);
check('solder-mask certificate passes at default budget',
    fdm.cert && fdm.cert.pass, fdm.cert ? `err ${(100 * fdm.cert.err).toFixed(3)}%` : 'no certificate');
check('solder-mask certificate err within tolerance', fdm.cert && fdm.cert.err < 0.01,
    fdm.cert ? `${(100 * fdm.cert.err).toFixed(3)}% < 1%` : 'no certificate');

const tri = await solve(sm, 'triangular');
const dC = Math.abs(fdm.m.RLGC.C - tri.m.RLGC.C) / tri.m.RLGC.C;
check('FDM solder-mask C matches triangular backend (±1%)', dC < 0.01,
    `fdm ${(fdm.m.RLGC.C * 1e12).toFixed(3)} vs tri ${(tri.m.RLGC.C * 1e12).toFixed(3)} pF/m, ${(100 * dC).toFixed(2)}%`);

console.log(failures === 0 ? '\nALL SM CERTIFICATION TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
