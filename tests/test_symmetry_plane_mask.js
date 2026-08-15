// Symmetry-plane mesh integrity with a dielectric edge ON the plane.
//
// A differential pair whose solder-mask side strips meet exactly at x = 0 puts a
// thin dielectric sliver edge ON the half-domain symmetry plane. OCC's boolean
// fragmentation emits that sliver's boundary nodes tens of nanometres off the
// wall (Precision::Confusion = 1e-7 model units), which used to defeat every
// exact-coordinate classifier downstream: the refinement smoother saw FREE nodes
// and dragged the mesh boundary off the plane, opening a VOID over the odd
// mode's dominant field region — odd C came out ~23% low, non-convergent, with
// certificates still claiming sub-1% error. Fixed by snapping near-wall nodes
// onto the domain walls right after occ mesh extraction (occ_to_mesh.js).
//
// These checks are reference-free identities:
//   1. half-domain odd C == full-domain odd C (same backend, same geometry)
//   2. tri odd C == FDM odd C (independent discretizations; FDM is converged
//      to ~0.1% on this geometry)
//   3. even mode agrees too (it was never broken — guards against a fix that
//      trades one mode for the other)
import { MicrostripSolver } from '../src/microstrip.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};

// Mask side strips (14 µm) meet exactly at x = 0 between the traces (28 µm apart).
const GEOM = {
    trace_width: 0.27e-3, trace_spacing: 28e-6, substrate_height: 0.37e-3,
    trace_thickness: 17e-6, gnd_thickness: 14e-6, epsilon_r: 7.5, tan_delta: 0.01,
    sigma_cond: 4.8e7, freq: 1.3e9, rq: 0, boundaries: ['open', 'open', 'open', 'gnd'],
    use_sm: true, sm_t_sub: 29e-6, sm_t_trace: 12e-6, sm_t_side: 14e-6, sm_er: 3.4, sm_tand: 0.02,
    nx: 30, ny: 30,
};

async function solveC(backend, triOpts = null) {
    return quiet(async () => {
        const s = new MicrostripSolver({ ...GEOM, mesh_backend: backend });
        if (triOpts) s.tri_opts = triOpts;
        const r = await s.solve_adaptive({
            max_iters: 12, energy_tol: 0.01, param_tol: 0.05, max_nodes: 30000, min_converged_passes: 2,
        });
        return { odd: r.modes[0].RLGC.C, even: r.modes[1].RLGC.C };
    });
}

const half = await solveC('triangular');
const full = await solveC('triangular', { symmetry: false });
const fdm = await solveC('rectilinear');

check('half-domain odd C == full-domain odd C (±3%)', relDiff(half.odd, full.odd) < 0.03,
    `half ${(half.odd * 1e12).toFixed(2)} vs full ${(full.odd * 1e12).toFixed(2)} pF/m (was ~23% low)`);
check('tri odd C == FDM odd C (±4%)', relDiff(half.odd, fdm.odd) < 0.04,
    `tri ${(half.odd * 1e12).toFixed(2)} vs fdm ${(fdm.odd * 1e12).toFixed(2)} pF/m`);
check('even mode still agrees across backends (±3%)', relDiff(half.even, fdm.even) < 0.03,
    `tri ${(half.even * 1e12).toFixed(2)} vs fdm ${(fdm.even * 1e12).toFixed(2)} pF/m`);

console.log(failures === 0 ? '\nALL SYMMETRY-PLANE MASK TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
