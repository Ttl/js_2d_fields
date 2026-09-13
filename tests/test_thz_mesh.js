// High-frequency meshing on the full-wave backend.
//
//   T1 near-field wavelength cap: at 300 GHz the main-solve mesh used to be sized by
//      geometry only (400 triangles), the eigensolve found no quasi-TEM candidate and
//      the reported eps_eff fell back to the static value (3.17 where dispersion puts
//      it near 4.1). The cap holds the whole domain at 3 cells per wavelength and the
//      region within three substrate-stack heights of the trace at 8, at the sweep's
//      top frequency, so the dispersed mode is resolved rather than merely found.
//   T2 the cap is inert where the geometric sizing already resolves the wavelength:
//      a 10 GHz solve builds the same mesh with the cap on and off.
//   T3 Modes tab domain shrink: on an auto-sized open domain the shrunken box gives the
//      same quasi-TEM eps_eff with fewer triangles at 100 GHz, and makes 300 GHz solvable
//      within the default budget (the full 6.3 mm domain is refused there).
//   T4 an enclosure is never shrunk.
//   T5 shift ladder: far above the quasi-TEM regime (0.15 mm trace on 0.1 mm FR4 at
//      0.4 to 1 THz) the eigenvalues near the static eps_eff are a dense cluster of
//      surface-wave modes and the quasi-TEM sits near max eps_r; the pick must walk its
//      shifts up to it instead of falling back to the static value (a 25% kink).
//   T6 Modes tab: with six modes requested the quasi-TEM is still listed at 800 GHz
//      (the hunt appends it when the modes near the shift do not overlap the drive).
//
// Run: node tests/test_thz_mesh.js
import { MicrostripSolver } from '../src/microstrip.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    console.log(`  ${cond ? 'PASS' : 'FAIL'}: ${name}${detail ? ` (${detail})` : ''}`);
    cond ? pass++ : fail++;
}
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};
const APP = { max_iters: 10, energy_tol: 0.01, param_tol: 0.05, max_nodes: 20000, min_converged_passes: 2 };
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
const MS = {
    trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
    epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, rq: 0, nx: 30, ny: 30,
    boundaries: ['open', 'open', 'open', 'gnd'],
};
function ms(opts, tri = {}) {
    const s = new MicrostripSolver({ ...MS, ...opts, mesh_backend: 'triangular' });
    s.use_causal_materials = false;
    s.tri_opts = { lossMethod: 'auto', ...tri };
    return s;
}
async function main(opts, tri = {}) {
    const s = ms(opts, tri);
    const t0 = Date.now();
    const r = await quiet(() => s.solve_adaptive({ ...APP }));
    const tb = await s._ensureTriBackend();
    return { s, r, tb, secs: (Date.now() - t0) / 1000, m: r.modes[0], warns: (s.modeWarnings || []).map(w => w.type) };
}
const MODES = { maxRefineIters: 10, refineTol: 0.01, maxNodes: 20000, minConvergedPasses: 2, certify: false, wavelengthDensity: 8 };
async function modes(opts, f, extra) {
    const s = ms(opts);
    const t0 = Date.now();
    try {
        const r = await quiet(() => s.solveModes(f, 8, null, { ...MODES, ...extra }));
        const tem = r.modes.filter(m => m.eps_eff != null && (m.overlap ?? 0) > 0.5).sort((a, b) => (b.overlap ?? 0) - (a.overlap ?? 0))[0];
        return { ok: true, nTris: r.nTris, eps: tem ? tem.eps_eff : NaN, overlap: tem ? tem.overlap : 0, secs: (Date.now() - t0) / 1000 };
    } catch (e) { return { ok: false, err: e.message.slice(0, 120) }; }
}

// T1
{
    console.log('T1 near-field wavelength cap at 300 GHz');
    const a = await main({ freq: 300e9 });
    check('cap applied', !!a.tb.nearField, a.tb.nearField ? `size ${(a.tb.nearField.size * 1e6).toFixed(1)} um, ${a.tb.nearField.nLambda} cells/lambda, dist ${(a.tb.nearField.dist * 1e3).toFixed(2)} mm` : 'null');
    check('eigensolve succeeds (no fallback warning)', !a.warns.includes('eigensolve'), a.warns.join(',') || 'none');
    check('dispersion resolved: eps_eff well above the static 3.17', a.m.eps_eff > 3.8, `eps_eff ${a.m.eps_eff.toFixed(4)}, ${a.tb.mesh.nTris} tris, ${a.secs.toFixed(1)} s`);
    // Without the cap the shift ladder still finds the mode on the geometric mesh;
    // the cap is resolution insurance for it, so the two must agree closely.
    const off = await main({ freq: 300e9 }, { nearFieldCap: false });
    check('eps_eff with and without the cap within 2%', rel(a.m.eps_eff, off.m.eps_eff) < 0.02, `${a.m.eps_eff.toFixed(4)} vs ${off.m.eps_eff.toFixed(4)}`);
    check('conductor loss unchanged by the cap (< 3%)', rel(a.m.RLGC.R, off.m.RLGC.R) < 0.03, `${a.m.RLGC.R.toFixed(2)} vs ${off.m.RLGC.R.toFixed(2)} ohm/m`);
}

// T2
{
    console.log('T2 cap inert at 10 GHz');
    const a = await main({ freq: 10e9 }), b = await main({ freq: 10e9 }, { nearFieldCap: false });
    check('no cap needed', a.tb.nearField === null);
    check('identical mesh with the option on and off', a.tb.mesh.nTris === b.tb.mesh.nTris && a.m.Z0 === b.m.Z0, `${a.tb.mesh.nTris} vs ${b.tb.mesh.nTris} tris`);
}

// T3
{
    console.log('T3 Modes tab domain shrink');
    const box = ms({ freq: 100e9 })._modes_domain_box();
    check('shrunken box smaller than the auto-sized domain', box && (box.x_max - box.x_min) < 0.6 * MS.trace_width * 18,
        box ? `${((box.x_max - box.x_min) * 1e3).toFixed(2)} x ${((box.y_max - box.y_min) * 1e3).toFixed(2)} mm` : 'null');
    const full = await modes({ freq: 100e9 }, 100e9, { shrinkDomain: false });
    const shrunk = await modes({ freq: 100e9 }, 100e9, { shrinkDomain: true });
    check('100 GHz: both solve', full.ok && shrunk.ok, (full.err || '') + (shrunk.err || ''));
    if (full.ok && shrunk.ok) {
        check('100 GHz: fewer triangles with the shrunken domain', shrunk.nTris < full.nTris, `${shrunk.nTris} vs ${full.nTris} (${shrunk.secs.toFixed(1)} s vs ${full.secs.toFixed(1)} s)`);
        check('100 GHz: quasi-TEM eps_eff within 1%', rel(shrunk.eps, full.eps) < 0.01, `${shrunk.eps.toFixed(4)} vs ${full.eps.toFixed(4)}`);
    }
    const refused = await modes({ freq: 300e9 }, 300e9, { shrinkDomain: false });
    const solved = await modes({ freq: 300e9 }, 300e9, { shrinkDomain: true });
    check('300 GHz: full domain refused at the default budget (control)', !refused.ok, refused.err || 'solved');
    check('300 GHz: shrunken domain solves', solved.ok, solved.err || `${solved.nTris} tris, ${solved.secs.toFixed(1)} s`);
    if (solved.ok) {
        const mainSolve = await main({ freq: 300e9 });
        check('300 GHz: modes quasi-TEM eps_eff within 5% of the main solve', rel(solved.eps, mainSolve.m.eps_eff) < 0.05,
            `${solved.eps.toFixed(4)} vs ${mainSolve.m.eps_eff.toFixed(4)} (overlap ${solved.overlap.toFixed(2)})`);
    }
}

// T4
{
    console.log('T4 enclosure is never shrunk');
    const enc = ms({ freq: 100e9, enclosure_width: 2e-3, enclosure_height: 1e-3, boundaries: ['gnd', 'gnd', 'gnd', 'gnd'] });
    check('no box for an enclosed line', enc._modes_domain_box() === null);
}

// T5
const THIN = { trace_width: 0.15e-3, substrate_height: 0.1e-3, tan_delta: 0.002 };
{
    console.log('T5 quasi-TEM pick far above the quasi-TEM regime');
    const a = await main({ ...THIN, freq: 1e12 });
    check('1 THz: eigensolve succeeds', !a.warns.includes('eigensolve'), a.warns.join(',') || 'none');
    check('1 THz: eps_eff dispersed toward eps_r (> 3.8)', a.m.eps_eff > 3.8, `eps_eff ${a.m.eps_eff.toFixed(4)}, ${a.tb.mesh.nTris} tris, ${a.secs.toFixed(1)} s`);
    const eps = [];
    let warned = false;
    for (const f of [200e9, 400e9, 600e9, 800e9]) {
        const rf = await quiet(() => a.s.computeAtFrequency(f, a.r));
        eps.push(rf.modes[0].eps_eff);
        if ((a.s.modeWarnings || []).some(w => w.type === 'eigensolve')) warned = true;
    }
    eps.push(a.m.eps_eff);
    check('200 GHz .. 1 THz: no eigensolve fallback', !warned);
    let monotone = true;
    for (let k = 1; k < eps.length; k++) if (eps[k] < eps[k - 1] * 0.99) monotone = false;
    check('eps_eff rises monotonically with frequency (no kink)', monotone, eps.map(v => v.toFixed(3)).join(' -> '));
}

// T6
{
    console.log('T6 Modes tab lists the quasi-TEM at 800 GHz with six modes');
    const s = ms({ ...THIN, freq: 800e9 });
    const t0 = Date.now();
    // 800 GHz needs a larger node budget even on the shrunken box (5.1k triangles at
    // 8 cells per wavelength against the 5k the default budget allows).
    const r = await quiet(() => s.solveModes(800e9, 6, null, { ...MODES, maxNodes: 40000, shrinkDomain: true }));
    const tem = r.modes.filter(m => m.status === 'propagating' && (m.overlap ?? 0) >= 0.5).sort((a, b) => b.overlap - a.overlap)[0];
    check('a propagating mode with overlap >= 0.5 is listed', !!tem,
        tem ? `eps_eff ${tem.eps_eff.toFixed(4)}, overlap ${tem.overlap.toFixed(3)}, ${r.modes.length} modes, ${r.nTris} tris, ${((Date.now() - t0) / 1000).toFixed(1)} s` : r.modes.map(m => `${m.eps_eff?.toFixed(3)}(${m.overlap.toFixed(2)})`).join(' '));
    if (tem) check('its eps_eff is well above the static value', tem.eps_eff > 3.8, tem.eps_eff.toFixed(4));
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
