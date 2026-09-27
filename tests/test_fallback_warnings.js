// Every fallback of the loss model is reported as an accuracy warning:
//   1. quasi-static: the DC current solve of the traces failing, and the thin-sheet
//      ground solve failing or over its group cap, while the return current spreads
//   2. full-wave: the low-frequency solve with ideal grounds failing
//   3. ordinary lines take no fallback on either backend at low and high frequency
// The failures are forced by replacing the failing step on a real solver.
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { check, quiet, APP, done } from './helpers.js';

const SOLVE = { ...APP, max_iters: 12 };
const kinds = r => (r.warnings || []).map(w => w.reason || w.type);
// Warnings that report a fallback of the loss model.
const FALLBACKS = ['dc-inductance-failed', 'ground-sheet-failed', 'mqs-perturbation', 'sibc-failed', 'wg-loss-failed',
    'dc-inductance-approx', 'ideal-ground-failed', 'refine-eigen', 'refine-metric', 'refine-h-metric',
    'mqs-solve-failed', 'mqs-rejected', 'eigensolve', 'eigen-anchor'];

// A microstrip on a finite ground in open space: the ground current spreads at low frequency.
const OPEN = 'units mm\nbounds open open open open\ndiel x=-3 w=6 y=0 h=0.21 er=4.4\n' +
    'gnd x=-3 w=6 y=-0.035 h=0.035\nsig+ x=-0.175 w=0.35 y=0.21 h=0.035';
const open = B => new CustomGeometrySolver({ text: OPEN, sigma_cond: 5.8e7, freq: 1e6, nx: 10, ny: 10, ...B });

// 1. Quasi-static.
{
    const s = open({});
    s._dc_signal_inductance = async () => { throw new Error('forced'); };
    const r = await quiet(async () => s.computeAtFrequency(1e4, await s.solve_adaptive(SOLVE)));
    check('quasi-static: a failed DC current solve is reported', kinds(r).includes('dc-inductance-failed'), kinds(r).join(', '));
}
{
    const s = open({});
    s._ground_sheet_setup = async () => ({ failed: '2000 ground groups, over the 1500 of the dense solve' });
    const r = await quiet(async () => s.computeAtFrequency(1e4, await s.solve_adaptive(SOLVE)));
    check('quasi-static: the ground sheet solve over its group cap is reported', kinds(r).includes('ground-sheet-failed'),
        kinds(r).join(', '));
}

// 2. Full-wave: the ideal-ground solve of a GCPW, whose pours reach the domain edge.
{
    const s = new MicrostripSolver({ trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
        epsilon_r: 4.4, tan_delta: 0, sigma_cond: 5.8e7, freq: 1e9, nx: 10, ny: 10, use_coplanar_gnd: true, gap: 0.2e-3,
        via_gap: 0.5e-3, use_vias: true, boundaries: ['open', 'open', 'open', 'gnd'], mesh_backend: 'triangular' });
    const r0 = await quiet(() => s.solve_adaptive(SOLVE));
    // Fail the ideal-ground variant only: its assembly cache throws on first use.
    const tri = s._triBackend, solve = tri._mqsSolve.bind(tri);
    tri._mqsSolve = (mesh, cr, f, sigma, opts) => solve(mesh, cr, f, sigma, opts, { get mesh() { throw new Error('forced'); } });
    const r = await quiet(() => s.computeAtFrequency(1.1e3, r0));
    const k = kinds(r);
    check('full-wave: a failed ideal-ground solve falls back to the volume solve and is reported',
        k.includes('ideal-ground-failed') && r.modes[0].RLGC.R > 0, k.join(', '));
}

// 3. No fallback on ordinary lines.
const plain = {
    microstrip: B => new MicrostripSolver({ trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6,
        gnd_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, nx: 10, ny: 10,
        boundaries: ['open', 'open', 'open', 'gnd'], ...B }),
    GCPW: B => new MicrostripSolver({ trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
        gnd_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, nx: 10, ny: 10,
        use_coplanar_gnd: true, gap: 0.2e-3, via_gap: 0.5e-3, use_vias: true, boundaries: ['open', 'open', 'open', 'gnd'], ...B }),
    'finite ground in open space': open,
};
for (const [name, mk] of Object.entries(plain)) {
    for (const tri of [false, true]) {
        const s = mk(tri ? { mesh_backend: 'triangular' } : {});
        const r0 = await quiet(() => s.solve_adaptive(SOLVE));
        const got = [];
        for (const f of [0, 1e3, 1e6, 1e9]) got.push(...kinds(await quiet(() => s.computeAtFrequency(f, r0))));
        const fb = [...new Set(got.filter(k => FALLBACKS.includes(k)))];
        check(`${name}, ${tri ? 'full-wave' : 'quasi-static'}: no fallback warning`, fb.length === 0, fb.join(', '));
    }
}

done();

console.log('\nALL CHECKS PASSED');
