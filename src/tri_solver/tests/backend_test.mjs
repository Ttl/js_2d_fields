// Phase 3 acceptance: full tri backend result vs the FDM solver.
import { MicrostripSolver } from '../../microstrip.js';
import { initTriBackend, TriBackend } from '../tri_backend.js';

const ctx = await initTriBackend();
const F = 10e9;

async function fdm(opts) {
    const s = new MicrostripSolver({ ...opts, freq: F });
    s.ensure_mesh();
    const cached = await s.solve_adaptive({ max_iters: 7, energy_tol: 0.005, param_tol: 0.02, max_nodes: 40000 });
    const r = await s.computeAtFrequency(F, cached);
    return r;
}

async function tri(opts) {
    const s = new MicrostripSolver({ ...opts, freq: F });
    const b = new TriBackend(ctx, s);
    await b.buildMesh();
    return b.solveAt(F);
}

function pct(a, b) { return (100 * (a - b) / b).toFixed(1) + '%'; }

function cmpMode(name, t, f) {
    console.log(`${name}`);
    console.log(`  Z0      tri=${t.Z0.toFixed(2)}  fdm=${f.Z0.toFixed(2)}  Δ=${pct(t.Z0, f.Z0)}`);
    console.log(`  eps_eff tri=${t.eps_eff.toFixed(3)}  fdm=${f.eps_eff.toFixed(3)}  Δ=${pct(t.eps_eff, f.eps_eff)}`);
    console.log(`  alpha_c tri=${t.alpha_c.toFixed(3)}  fdm=${f.alpha_c.toFixed(3)}  Δ=${pct(t.alpha_c, f.alpha_c)}  dB/m`);
    console.log(`  alpha_d tri=${t.alpha_d.toFixed(4)}  fdm=${f.alpha_d.toFixed(4)}  Δ=${pct(t.alpha_d, f.alpha_d)}  dB/m`);
    console.log(`  R       tri=${t.RLGC.R.toFixed(2)}  fdm=${f.RLGC.R.toFixed(2)}  | L tri=${(t.RLGC.L*1e9).toFixed(1)} fdm=${(f.RLGC.L*1e9).toFixed(1)} nH/m | C tri=${(t.RLGC.C*1e12).toFixed(1)} fdm=${(f.RLGC.C*1e12).toFixed(1)} pF/m`);
}

// --- single-ended microstrip ---
const msOpts = { trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7 };
{
    const t = await tri(msOpts); const f = await fdm(msOpts);
    cmpMode('=== microstrip (single) ===', t.modes[0], f.modes[0]);
}

// --- microstrip with roughness ---
const rqOpts = { ...msOpts, rq: 1e-6 };
{
    const t = await tri(rqOpts); const f = await fdm(rqOpts);
    cmpMode('=== microstrip + roughness (Rq=1um) ===', t.modes[0], f.modes[0]);
}

// --- differential microstrip ---
const diffOpts = { trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6, trace_spacing: 0.1e-3, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7 };
{
    const t = await tri(diffOpts); const f = await fdm(diffOpts);
    console.log('=== differential microstrip ===');
    console.log(`  Z_diff   tri=${t.Z_diff.toFixed(2)}  fdm=${f.Z_diff.toFixed(2)}  Δ=${pct(t.Z_diff, f.Z_diff)}`);
    console.log(`  Z_common tri=${t.Z_common.toFixed(2)}  fdm=${f.Z_common.toFixed(2)}  Δ=${pct(t.Z_common, f.Z_common)}`);
    cmpMode('  -- odd --', t.modes.find(m=>m.mode==='odd'), f.modes.find(m=>m.mode==='odd'));
    cmpMode('  -- even --', t.modes.find(m=>m.mode==='even'), f.modes.find(m=>m.mode==='even'));
}

console.log('\nBACKEND TEST DONE');
