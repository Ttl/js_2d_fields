// Modes viewer: solve the full-wave eigenproblem for several modes via the FieldSolver2D
// public API (solveModes/getModeField), and verify classification + field resampling.
import { MicrostripSolver } from '../../microstrip.js';

const s = new MicrostripSolver({ trace_width: 0.35e-3, substrate_height: 0.21e-3,
    trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, freq: 20e9 });
s.mesh_backend = 'triangular';

const f = 20e9;
const res = await s.solveModes(f, 6);
console.log('nconv=', res.nconv, 'N=', res.N, 'nTris=', res.nTris, 'modes=', res.modes.length);
if (res.error) console.log('ERROR:', res.error);

for (const m of res.modes) {
    const eeff = m.eps_eff != null ? m.eps_eff.toFixed(4) : '  –   ';
    console.log(`  #${m.idx}  ${m.status.padEnd(11)} eps_eff=${eeff}  g2=${m.g2Re.toExponential(2)}  ovl=${m.overlap.toFixed(3)}`);
}

// The quasi-TEM (fundamental) mode = highest overlap, must be propagating with eps_eff
// somewhere between 1 and eps_r.
let best = 0, bestOvl = -1;
res.modes.forEach((m, i) => { if (m.overlap > bestOvl) { bestOvl = m.overlap; best = i; } });
const q = res.modes[best];
console.log(`\nquasi-TEM = mode ${best}: status=${q.status} eps_eff=${q.eps_eff?.toFixed(4)} overlap=${q.overlap.toFixed(3)}`);

// Resample its field; verify the grid is populated and the field peaks near the trace.
const grid = s.getModeField(best);
let emax = 0, nonzero = 0, total = 0;
for (let j = 0; j < grid.E.length; j++) for (let i = 0; i < grid.E[0].length; i++) {
    const e = grid.E[j][i]; emax = Math.max(emax, e); total++;
    if (e > 1e-12) nonzero++;
}
console.log(`field grid ${grid.E[0].length}x${grid.E.length}  |Et|max=${emax.toExponential(2)}  nonzero=${nonzero}/${total}`);

const ok = res.modes.length > 0 && q.status === 'propagating' &&
    q.eps_eff > 1 && q.eps_eff < 4.4 && q.overlap > 0.5 &&
    emax > 0 && nonzero > total * 0.1;
console.log(ok ? '\nMODES OK' : '\nMODES SUSPECT');
process.exit(ok ? 0 : 1);
