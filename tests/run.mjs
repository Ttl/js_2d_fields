// Test runner: runs every test in a tier (no early stop), prints a one-line status
// per test, shows the output of failing tests only, and exits nonzero if any failed.
//
//   node tests/run.mjs [fast|slow|e2e|fuzz|all] [-j N]
//
// Tests run in PARALLEL (default jobs: min(8, cpus - 2); e2e defaults to 1 because its
// tests share the dev server on localhost:8731). Each test's output is buffered and
// only printed on failure, so parallel output stays readable; status lines appear in
// completion order. `-j 1` restores strictly sequential runs.
//
// Scheduling is longest-first using the per-entry `cost` hints (seconds, measured
// 2026-08-15 on a 16-core box) — with a worker pool the wall time is then bounded by
// max(longest single test, total/jobs). The hints only order the queue; a stale hint
// costs a little wall time, never correctness.
//
// Tiers:
//   fast — run before every commit (~30 s at -j8): geometry construction, eigensolver
//          pencil, enclosed-box analytic modes, closed-form parallel-plate, coax and
//          rectangular waveguide.
//   slow — full solver validation (~4 min at -j8, ~28 min sequential): tri-backend
//          correctness suite, mode viewer, per-side plating, plating behavior
//          invariants, the reference-free identity suite (test_invariants.js), a
//          pinned-seed QS-vs-fullwave fuzzer smoke, and the reference suite on BOTH
//          backends.
//   e2e  — browser tests; need the dev server on localhost:8731 and chromium.
//   fuzz — the full 40-case QS-vs-fullwave fuzzer run (~10 min, deterministic seed).
import { spawn } from 'child_process';
import os from 'os';

const TIERS = {
    fast: [
        { file: 'tests/test_geometry.js', cost: 1 },
        { file: 'src/tri_solver/tests/eigen_pencil_test.mjs', cost: 1 },
        { file: 'src/tri_solver/tests/box_modes_test.mjs', cost: 16 },
        { file: 'tests/test_parallel_plate.js', cost: 28 },
        { file: 'tests/test_coax.js', cost: 11 },
        { file: 'tests/test_rect_waveguide.js', cost: 7 },
        { file: 'tests/test_djordjevic_sarkar.js', cost: 1 },
        { file: 'tests/test_meshability.js', cost: 1 },
        { file: 'tests/test_causal_effect.js', cost: 2 },
        { file: 'tests/test_layered_roughness.js', cost: 1 },
        { file: 'tests/test_interpolating_sweep.js', cost: 2 },
        { file: 'tests/test_mtl_sparams.js', cost: 1 },
        { file: 'tests/test_open_boundary_warning.js', cost: 1 },
        { file: 'tests/test_broadside_proximity_warning.js', cost: 15 },
        { file: 'tests/test_gcpw_qs_continuity.js', cost: 9 },
    ],
    slow: [
        { file: 'tests/test_fullwave_correctness.js', cost: 98 },
        { file: 'tests/test_gcpw_mqs.js', cost: 82 },
        { file: 'tests/test_symmetry_plane_mask.js', cost: 61 },
        { file: 'tests/test_symmetry_half_full.js', cost: 177 },
        { file: 'src/tri_solver/tests/modes_test.mjs', cost: 15 },
        { file: 'src/tri_solver/tests/plating_perside_test.mjs', cost: 138 },
        { file: 'tests/test_vs_ref.js', cost: 44 },
        { file: 'tests/test_vs_ref.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_vs_ref.js [triangular]', cost: 122 },
        { file: 'tests/test_causal_effect.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_causal_effect.js [triangular]', cost: 6 },
        { file: 'tests/test_modal_continuity.js', cost: 46 },
        { file: 'tests/test_modal_continuity.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_modal_continuity.js [triangular]', cost: 226 },
        { file: 'tests/test_asym_sparam_export.js', cost: 13 },
        { file: 'tests/test_asym_sparam_export.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_asym_sparam_export.js [triangular]', cost: 52 },
        { file: 'tests/test_interp_asym_sparams.js', cost: 11 },
        { file: 'tests/test_interp_asym_sparams.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_interp_asym_sparams.js [triangular]', cost: 36 },
        { file: 'tests/test_thick_plating.js', cost: 12 },
        { file: 'tests/test_top_only_plating.js', cost: 27 },
        { file: 'tests/test_poor_side_plating.js', cost: 27 },
        { file: 'tests/test_plating_thicker_than_trace.js', cost: 8 },
        { file: 'tests/test_corner_plating.js', cost: 35 },
        { file: 'tests/test_corner_roughness.js', cost: 12 },
        { file: 'tests/test_corner_loss_regularization.js', cost: 43 },
        { file: 'tests/test_coax_plating.js', cost: 30 },
        { file: 'tests/test_surface_reactance.js', cost: 75 },
        { file: 'tests/test_certification.js', cost: 39 },
        { file: 'tests/test_invariants.js', cost: 110 },
        { file: 'tests/test_sm_certification.js', cost: 21 },
        { file: 'tests/fuzz_qs_vs_fullwave.js', args: ['6', '3', '15'], label: 'tests/fuzz_qs_vs_fullwave.js [smoke N=6 seed=3]', cost: 74 },
        { file: 'tests/fuzz_fullwave_interp.js', args: ['3', '1'], label: 'tests/fuzz_fullwave_interp.js [smoke N=3 seed=1]', cost: 176 },
    ],
    e2e: [
        'src/tri_solver/tests/e2e.mjs',
        'src/tri_solver/tests/e2e_occ.mjs',
        'src/tri_solver/tests/e2e_modes.mjs',
        'src/tri_solver/tests/e2e_waveguide.mjs',
    ],
    fuzz: [
        { file: 'tests/fuzz_qs_vs_fullwave.js', args: ['40', '1', '15'], label: 'tests/fuzz_qs_vs_fullwave.js [N=40 seed=1]', cost: 600 },
        { file: 'tests/fuzz_fullwave_interp.js', args: ['12', '1'], label: 'tests/fuzz_fullwave_interp.js [N=12 seed=1]', cost: 700 },
    ],
};

// ---- CLI ----
const argv = process.argv.slice(2);
let tier = 'fast', jobsArg = null;
// A malformed -j must fail. parseInt(undefined) is NaN and a NaN job count
// silently produced zero workers.
const parseJobs = (v, flag) => {
    const n = parseInt(v, 10);
    if (!Number.isFinite(n) || n < 1) {
        console.error(`Invalid ${flag} value "${v ?? ''}" — expected a positive integer.`);
        process.exit(2);
    }
    return n;
};
for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    if (a === '-j' || a === '--jobs') jobsArg = parseJobs(argv[++i], a);
    else if (/^-j/.test(a)) jobsArg = parseJobs(a.slice(2), '-j');
    else if (/^-/.test(a)) { console.error(`Unknown option "${a}".`); process.exit(2); }
    else tier = a;
}
const tests = tier === 'all' ? [...TIERS.fast, ...TIERS.slow] : TIERS[tier];
if (!tests) {
    console.error(`Unknown tier "${tier}" — use fast, slow, e2e or all.`);
    process.exit(2);
}
// e2e tests share one dev server / browser, so they stay sequential unless forced.
const defaultJobs = tier === 'e2e' ? 1 : Math.min(8, Math.max(1, os.cpus().length - 2));
const jobs = Math.max(1, jobsArg ?? defaultJobs);

const norm = t => (typeof t === 'string' ? { file: t } : t);
// Longest-first minimizes makespan for a worker pool; ties keep list order.
const queue = tests.map((t, i) => ({ ...norm(t), i }))
    .sort((a, b) => (b.cost ?? 30) - (a.cost ?? 30) || a.i - b.i);

function runOne(t) {
    const { file, args = [], env = {}, label = file } = t;
    return new Promise(resolve => {
        const t0 = Date.now();
        const p = spawn('node', [file, ...args], { env: { ...process.env, ...env } });
        let out = '';
        p.stdout.on('data', d => { out += d; });
        p.stderr.on('data', d => { out += d; });
        const killer = setTimeout(() => { p.kill('SIGKILL'); out += '\n[runner] TIMEOUT after 15 min — killed\n'; }, 15 * 60 * 1000);
        p.on('close', status => {
            clearTimeout(killer);
            const dt = ((Date.now() - t0) / 1000).toFixed(1);
            const ok = status === 0;
            console.log(`${ok ? '✓' : '✗'} ${label}  (${dt}s)`);
            if (!ok) console.log(out.split('\n').map(l => '    ' + l).join('\n'));
            resolve({ t: label, ok });
        });
        p.on('error', err => {
            clearTimeout(killer);
            console.log(`✗ ${label}  (spawn failed: ${err.message})`);
            resolve({ t: label, ok: false });
        });
    });
}

const t0 = Date.now();
if (jobs > 1) console.log(`running ${queue.length} tests with ${jobs} parallel jobs\n`);
const results = [];
await Promise.all(Array.from({ length: Math.min(jobs, queue.length) }, async () => {
    for (let t; (t = queue.shift()) !== undefined;) results.push(await runOne(t));
}));

const failed = results.filter(r => !r.ok);
console.log(`\n${results.length - failed.length}/${results.length} test files passed` +
    ` in ${((Date.now() - t0) / 1000).toFixed(0)}s` +
    (failed.length ? ` — FAILED: ${failed.map(f => f.t).join(', ')}` : ''));
process.exit(failed.length ? 1 : 0);
