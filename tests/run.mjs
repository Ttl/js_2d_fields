// Test runner: runs every test in a tier (no early stop), prints a one-line status
// per test, shows the output of failing tests only, and exits nonzero if any failed.
//
//   node tests/run.mjs [fast|slow|e2e|fuzz|all]     (default: fast)
//
// Tiers:
//   fast — run before every commit (~1.5 min): geometry construction, eigensolver
//          pencil, enclosed-box analytic modes.
//   slow — full solver validation: tri-backend correctness suite, mode viewer,
//          per-side plating, plating behavior invariants, a pinned-seed
//          QS-vs-fullwave fuzzer smoke, and the reference suite on BOTH backends
//          (the MESH_BACKEND=triangular pass is the longest single entry).
//   e2e  — browser tests; need the dev server on localhost:8731 and chromium.
//   fuzz — the full 40-case QS-vs-fullwave fuzzer run (~10 min, deterministic seed).
import { spawnSync } from 'child_process';

const TIERS = {
    fast: [
        'tests/test_geometry.js',
        'src/tri_solver/tests/eigen_pencil_test.mjs',
        'src/tri_solver/tests/box_modes_test.mjs',
        'tests/test_parallel_plate.js',
        'tests/test_djordjevic_sarkar.js',
        'tests/test_meshability.js',
        'tests/test_causal_effect.js',
        'tests/test_layered_roughness.js',
        'tests/test_interpolating_sweep.js',
        'tests/test_mtl_sparams.js',
        'tests/test_open_boundary_warning.js',
    ],
    slow: [
        'tests/test_fullwave_correctness.js',
        'src/tri_solver/tests/modes_test.mjs',
        'src/tri_solver/tests/plating_perside_test.mjs',
        'tests/test_vs_ref.js',
        { file: 'tests/test_vs_ref.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_vs_ref.js [triangular]' },
        { file: 'tests/test_causal_effect.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_causal_effect.js [triangular]' },
        'tests/test_modal_continuity.js',
        { file: 'tests/test_modal_continuity.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_modal_continuity.js [triangular]' },
        'tests/test_asym_sparam_export.js',
        { file: 'tests/test_asym_sparam_export.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_asym_sparam_export.js [triangular]' },
        'tests/test_interp_asym_sparams.js',
        { file: 'tests/test_interp_asym_sparams.js', env: { MESH_BACKEND: 'triangular' }, label: 'tests/test_interp_asym_sparams.js [triangular]' },
        'tests/test_thick_plating.js',
        'tests/test_top_only_plating.js',
        'tests/test_poor_side_plating.js',
        'tests/test_plating_thicker_than_trace.js',
        'tests/test_corner_plating.js',
        'tests/test_corner_roughness.js',
        { file: 'tests/fuzz_qs_vs_fullwave.js', args: ['6', '3', '15'], label: 'tests/fuzz_qs_vs_fullwave.js [smoke N=6 seed=3]' },
    ],
    e2e: [
        'src/tri_solver/tests/e2e.mjs',
        'src/tri_solver/tests/e2e_occ.mjs',
        'src/tri_solver/tests/e2e_modes.mjs',
    ],
    fuzz: [
        { file: 'tests/fuzz_qs_vs_fullwave.js', args: ['40', '1', '15'], label: 'tests/fuzz_qs_vs_fullwave.js [N=40 seed=1]' },
    ],
};

const arg = process.argv[2] || 'fast';
const tests = arg === 'all' ? [...TIERS.fast, ...TIERS.slow] : TIERS[arg];
if (!tests) {
    console.error(`Unknown tier "${arg}" — use fast, slow, e2e or all.`);
    process.exit(2);
}

const results = [];
for (const t of tests) {
    const { file, args = [], env = {}, label = file } = typeof t === 'string' ? { file: t } : t;
    const t0 = Date.now();
    const r = spawnSync('node', [file, ...args], {
        encoding: 'utf8', timeout: 15 * 60 * 1000, env: { ...process.env, ...env },
    });
    const dt = ((Date.now() - t0) / 1000).toFixed(1);
    const ok = r.status === 0;
    results.push({ t: label, ok });
    console.log(`${ok ? '✓' : '✗'} ${label}  (${dt}s)`);
    if (!ok) {
        const out = (r.stdout || '') + (r.stderr || '');
        console.log(out.split('\n').map(l => '    ' + l).join('\n'));
    }
}

const failed = results.filter(r => !r.ok);
console.log(`\n${results.length - failed.length}/${results.length} test files passed` +
    (failed.length ? ` — FAILED: ${failed.map(f => f.t).join(', ')}` : ''));
process.exit(failed.length ? 1 : 0);
