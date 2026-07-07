// Test runner: runs every test in a tier (no early stop), prints a one-line status
// per test, shows the output of failing tests only, and exits nonzero if any failed.
//
//   node tests/run.mjs [fast|slow|e2e|all]     (default: fast)
//
// Tiers:
//   fast — run before every commit (~1.5 min): geometry construction, eigensolver
//          pencil, enclosed-box analytic modes.
//   slow — full solver validation (~10 min): tri-backend correctness suite, mode
//          viewer, per-side plating, and the reference suite on the rectilinear
//          backend. For the triangular backend reference run (much slower):
//          MESH_BACKEND=triangular node tests/test_vs_ref.js
//   e2e  — browser tests; need the dev server on localhost:8731 and chromium.
import { spawnSync } from 'child_process';

const TIERS = {
    fast: [
        'tests/test_geometry.js',
        'src/tri_solver/tests/eigen_pencil_test.mjs',
        'src/tri_solver/tests/box_modes_test.mjs',
    ],
    slow: [
        'tests/test_fullwave_correctness.js',
        'src/tri_solver/tests/modes_test.mjs',
        'src/tri_solver/tests/plating_perside_test.mjs',
        'tests/test_vs_ref.js',
    ],
    e2e: [
        'src/tri_solver/tests/e2e.mjs',
        'src/tri_solver/tests/e2e_occ.mjs',
        'src/tri_solver/tests/e2e_modes.mjs',
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
    const t0 = Date.now();
    const r = spawnSync('node', [t], { encoding: 'utf8', timeout: 15 * 60 * 1000 });
    const dt = ((Date.now() - t0) / 1000).toFixed(1);
    const ok = r.status === 0;
    results.push({ t, ok });
    console.log(`${ok ? '✓' : '✗'} ${t}  (${dt}s)`);
    if (!ok) {
        const out = (r.stdout || '') + (r.stderr || '');
        console.log(out.split('\n').map(l => '    ' + l).join('\n'));
    }
}

const failed = results.filter(r => !r.ok);
console.log(`\n${results.length - failed.length}/${results.length} test files passed` +
    (failed.length ? ` — FAILED: ${failed.map(f => f.t).join(', ')}` : ''));
process.exit(failed.length ? 1 : 0);
