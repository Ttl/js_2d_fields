// Phase 1+4 integration: drive the solver exactly like app_solver does, with
// mesh_backend:'triangular'. Confirms dispatch, result shape, and field state.
import { MicrostripSolver } from '../../microstrip.js';

const s = new MicrostripSolver({
    trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7,
    mesh_backend: 'triangular', freq: 10e9,
});

// app flow: ensure_mesh (no-op for tri) → solve_adaptive (build+solve) → computeAtFrequency sweep
s.ensure_mesh();
let progressSeen = false;
const cached = await s.solve_adaptive({ onProgress: () => { progressSeen = true; }, max_iters: 6 });
console.log('solve_adaptive: Z0=%s eps_eff=%s (progress cb fired: %s)',
    cached.modes[0].Z0.toFixed(2), cached.modes[0].eps_eff.toFixed(3), progressSeen);

for (const f of [1e9, 10e9, 30e9]) {
    const r = await s.computeAtFrequency(f, cached);
    const m = r.modes[0];
    console.log(`  f=${(f/1e9).toFixed(0).padStart(2)}GHz  Z0=${m.Z0.toFixed(2)}  eps_eff=${m.eps_eff.toFixed(3)}  alpha_c=${m.alpha_c.toFixed(2)}  alpha_d=${m.alpha_d.toFixed(2)}  Zc=${m.Zc.re.toFixed(2)}${m.Zc.im>=0?'+':''}${m.Zc.im.toFixed(2)}j`);
}

console.log('solver fields set: x=%d y=%d Ex_modes=%d triMesh=%s solution_valid=%s',
    s.x.length, s.y.length, s.Ex.length, !!s.triMesh, s.solution_valid);

// differential through the same path
const d = new MicrostripSolver({
    trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6,
    trace_spacing: 0.1e-3, epsilon_r: 4.4, tan_delta: 0.02,
    mesh_backend: 'triangular', freq: 10e9,
});
d.ensure_mesh();
const dc = await d.solve_adaptive({});
const r = await d.computeAtFrequency(10e9, dc);
console.log(`diff: Z_diff=${r.Z_diff.toFixed(2)} Z_common=${r.Z_common.toFixed(2)} RLGC C11=${(r.RLGC_matrix.C[0][0]*1e12).toFixed(1)} C12=${(r.RLGC_matrix.C[0][1]*1e12).toFixed(1)} pF/m  Ex_modes=${d.Ex.length}`);

console.log('\nINTEGRATION OK');
