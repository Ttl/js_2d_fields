// Focused check of the tricky cases through the real backend.
import { MicrostripSolver } from '../../microstrip.js';
import { initTriBackend, TriBackend } from '../tri_backend.js';

const ctx = await initTriBackend();
function run(name, opts, extras = {}) {
    const s = new MicrostripSolver(opts);
    const b = new TriBackend(ctx, s, extras);
    const r = b.solveAt(opts.freq);
    const m = r.modes[0];
    console.log(`${name}: Z0=${m.Z0.toFixed(2)} eps_eff=${m.eps_eff.toFixed(3)} eps_mode=${m.eps_eff_mode.toFixed(3)} R=${m.RLGC.R.toFixed(2)} cond_loss=${m.alpha_c.toFixed(3)} diel=${m.alpha_d.toFixed(3)} C=${(m.RLGC.C*1e12).toFixed(1)}pF nTris=${b.mesh.nTris} sym=${b.symmetry}`);
}

run('50Ω (R ref 3.26, cl 0.285)', { substrate_height:1.6e-3, trace_width:3e-3, trace_thickness:35e-6, gnd_thickness:35e-6, epsilon_r:4.5, tan_delta:0.02, sigma_cond:5.8e7, freq:1e9, boundaries:['open','open','open','gnd'] });
run('1kHz (R ref 1.01)', { substrate_height:1.6e-3, trace_width:3e-3, trace_thickness:35e-6, gnd_thickness:35e-6, epsilon_r:4.5, tan_delta:0, sigma_cond:1e7, enclosure_width:50e-3, freq:1e3, boundaries:['open','open','open','gnd'] });
run('embed (eps ref 3.621)', { substrate_height:1.6e-3, trace_width:3e-3, trace_thickness:35e-6, gnd_thickness:35e-6, epsilon_r:4.5, tan_delta:0.02, sigma_cond:5.8e7, freq:1e9, top_diel_h:0.2e-3, top_diel_er:4.5, boundaries:['open','open','open','gnd'] });
run('embed hi-res', { substrate_height:1.6e-3, trace_width:3e-3, trace_thickness:35e-6, gnd_thickness:35e-6, epsilon_r:4.5, tan_delta:0.02, sigma_cond:5.8e7, freq:1e9, top_diel_h:0.2e-3, top_diel_er:4.5, boundaries:['open','open','open','gnd'] }, { maxNodes: 40000, maxRefineIters: 9 });
run('stripline (cl ref 3.1)', { substrate_height:0.2e-3, trace_width:0.15e-3, trace_thickness:35e-6, gnd_thickness:16e-6, epsilon_r:4.1, epsilon_r_top:4.1, tan_delta:0.02, tan_delta_top:0.02, enclosure_height:0.2e-3+35e-6, freq:1e9, boundaries:['open','open','gnd','gnd'] });
