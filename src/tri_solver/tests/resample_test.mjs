// Phase 4: verify the backend fills V/Ex/Ey grids + triMesh, and the field is sane.
import { MicrostripSolver } from '../../microstrip.js';
import { initTriBackend, TriBackend } from '../tri_backend.js';

const ctx = await initTriBackend();
const s = new MicrostripSolver({ trace_width: 1.5e-3, substrate_height: 0.5e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, freq: 5e9 });
const b = new TriBackend(ctx, s);
const res = b.solveAt(5e9);

console.log('result Z0=', res.modes[0].Z0.toFixed(2), 'eps_eff=', res.modes[0].eps_eff.toFixed(3));
console.log('solver.x len=', s.x.length, 'solver.y len=', s.y.length);
console.log('solver.Ex modes=', s.Ex.length, 'grid=', s.Ex[0].length, 'x', s.Ex[0][0].length);
console.log('triMesh tris=', s.triMesh.nTris, 'nodes=', s.triMesh.nodes.length / 2);

// Field sanity: V should range ~[0,1] (signal at 1V, ground 0), E nonzero near trace.
let vmin = 1e9, vmax = -1e9, emax = 0, nonzero = 0;
const V = s.V[0], Ex = s.Ex[0], Ey = s.Ey[0];
for (let j = 0; j < V.length; j++) for (let i = 0; i < V[0].length; i++) {
    vmin = Math.min(vmin, V[j][i]); vmax = Math.max(vmax, V[j][i]);
    const e = Math.hypot(Ex[j][i], Ey[j][i]); emax = Math.max(emax, e);
    if (Math.abs(V[j][i]) > 1e-6) nonzero++;
}
console.log(`V range [${vmin.toFixed(3)}, ${vmax.toFixed(3)}]  |E|max=${emax.toExponential(2)}  nonzeroV=${nonzero}/${V.length*V[0].length}`);

// V at a point just above the trace center (x=0, y just above substrate top 0.5mm) should be high
const xi = V[0].length >> 1;  // x≈0
let yi = 0; for (let j = 0; j < s.y.length; j++) if (s.y[j] <= 0.5e-3) yi = j;
console.log(`V near trace (x≈0, y≈0.5mm) = ${V[yi][xi].toFixed(3)} (expect close to 1)`);

console.log(vmin > -0.05 && vmax < 1.05 && vmax > 0.9 && emax > 0 ? '\nRESAMPLE OK' : '\nRESAMPLE SUSPECT');
