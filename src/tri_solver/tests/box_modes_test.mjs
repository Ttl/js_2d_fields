// Enclosed-box mode-count test — the mode viewer against ANALYTIC waveguide truth.
//
// An all-PEC enclosed microstrip is a rectangular waveguide (plus a small trace
// perturbation): its higher-order modes have exact cutoffs kc = π·√((m/a)²+(n/b)²)
// and ε_eff = ε − (c·kc/2πf)². The assertion is a COUNT plus values: every mode
// above cutoff must be found — a tolerance-only check passes when modes are
// silently dropped, which is exactly what regression 5a34928 did (strict residual
// gate + fixed Krylov subspace lost all far-from-shift cavity modes; with εr=4.4
// only the quasi-TEM survived).
//
// Cavity dimensions (see occ_to_mesh._clipDomain: full-span boundary ground slabs
// — bottom, lid AND the side slabs microstrip.js builds for 'gnd' side boundaries,
// all t_gnd thick — are absorbed into PEC walls at their inner faces):
//   a = domain_width − 2·t_gnd            (side slab inner faces)
//   b = domain_height − t_gnd             (bottom face y=0 → lid inner face)
// Above cutoff, (m,n) with m,n ≥ 1 contributes TWO modes (TE + TM, degenerate in
// the empty box); (m,0)/(0,n) only TE.
//
// Run: node src/tri_solver/tests/box_modes_test.mjs
import { MicrostripSolver } from '../../microstrip.js';

const C0 = 299792458;
let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}

// Analytic mode list of an a×b PEC box, uniform fill eps_r, at frequency f.
// minEps skips modes too close to cutoff to assert robustly (trace perturbation).
function boxModes(a, b, epsr, f, minEps = 0.05) {
    const k0 = 2 * Math.PI * f / C0;
    const out = [];
    for (let m = 0; m <= 6; m++) for (let n = 0; n <= 6; n++) {
        if (m === 0 && n === 0) continue;
        const eps = epsr - (Math.PI * Math.hypot(m / a, n / b) / k0) ** 2;
        if (eps <= minEps) continue;
        for (let d = 0; d < (m >= 1 && n >= 1 ? 2 : 1); d++)
            out.push({ name: `${d ? 'TM' : 'TE'}${m}${n}`, eps });
    }
    return out.sort((p, q) => q.eps - p.eps);
}

function makeSolver(opts) {
    const s = new MicrostripSolver({
        substrate_height: 0.21e-3, freq: 10e9,
        enclosure_width: 4e-3, enclosure_height: 4e-3,
        boundaries: ['gnd', 'gnd', 'gnd', 'gnd'],
        ...opts,
    });
    s.mesh_backend = 'triangular';
    return s;
}

const quiet = async (fn) => {
    const log = console.log; console.log = () => {};
    try { return await fn(); } finally { console.log = log; }
};

// ---------- 1. Air-filled box + tiny trace: exact truth, count AND values ----------
{
    const f = 60e9;
    // Tiny trace so the box modes are barely perturbed (assertions get ~2% slack)
    const s = makeSolver({ trace_width: 0.1e-3, trace_thickness: 10e-6, epsilon_r: 1.0, tan_delta: 0 });
    const a = s.domain_width - 2 * s.t_gnd;
    const b = s.domain_height - s.t_gnd;
    const expected = boxModes(a, b, 1.0, f);
    console.log(`air box ${(a * 1e3).toFixed(3)}×${(b * 1e3).toFixed(3)} mm at ${f / 1e9} GHz — analytic: ` +
        expected.map(e => `${e.name}=${e.eps.toFixed(4)}`).join(' '));

    const res = await quiet(() => s.solveModes(f, expected.length + 5, null, { maxRefineIters: 1 }));
    check('air box: solve succeeded', !res.error, res.error || `nconv=${res.nconv}, nTris=${res.nTris}`);
    const prop = (res.modes || []).filter(m => m.status === 'propagating');
    console.log('  found propagating: ' + prop.map(m => m.eps_eff.toFixed(4)).join(' '));

    // quasi-TEM: highest overlap with the static drive; homogeneous fill → ε_eff = εr exactly
    let qi = -1, qOvl = -1;
    prop.forEach((m, i) => { if (m.overlap > qOvl) { qOvl = m.overlap; qi = i; } });
    const qtem = prop[qi];
    check('air box: quasi-TEM ε_eff = 1 (homogeneous fill invariant)',
        qtem && Math.abs(qtem.eps_eff - 1.0) < 0.01, qtem ? `eps=${qtem.eps_eff.toFixed(4)} ovl=${qOvl.toFixed(2)}` : 'missing');

    // COUNT: exactly the analytic modes + the quasi-TEM, nothing extra (spurious leak)
    const cavity = prop.filter((_, i) => i !== qi);
    check('air box: propagating mode COUNT matches analytic', cavity.length === expected.length,
        `found ${cavity.length}, expected ${expected.length}`);

    // VALUES: greedy 1-to-1 match, tolerance covers trace perturbation + mesh error
    const pool = cavity.slice();
    let maxErr = 0, allMatched = true;
    for (const e of expected) {
        let best = -1, bestErr = Infinity;
        pool.forEach((m, i) => { const err = Math.abs(m.eps_eff - e.eps); if (err < bestErr) { bestErr = err; best = i; } });
        const tol = Math.max(0.03 * e.eps, 0.02);
        if (best >= 0 && bestErr < tol) pool.splice(best, 1); else allMatched = false;
        maxErr = Math.max(maxErr, best >= 0 ? bestErr : Infinity);
        console.log(`  ${e.name}: analytic ${e.eps.toFixed(4)}, error ${isFinite(bestErr) ? bestErr.toFixed(4) : '—'}`);
    }
    check('air box: every analytic mode found at its ε_eff', allMatched,
        `max |Δε|=${isFinite(maxErr) ? maxErr.toFixed(4) : '∞'}`);
}

// ---------- 2. εr=4.4 substrate: the 5a34928 regression geometry ----------
// Dielectric loading only LOWERS cutoffs (variational), so every air-box mode is a
// lower bound and the air-box count is a minimum count. Pre-fix, this solve
// returned ONLY the quasi-TEM (shift at ε≈3.4, cavity modes at ε≈0.5 dropped).
{
    const f = 50e9;
    const s = makeSolver({ trace_width: 0.35e-3, trace_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02 });
    const a = s.domain_width - 2 * s.t_gnd;
    const b = s.domain_height - s.t_gnd;
    const airRef = boxModes(a, b, 1.0, f);
    console.log(`\nεr=4.4 box at ${f / 1e9} GHz — air-box lower bounds: ` +
        airRef.map(e => `${e.name}=${e.eps.toFixed(4)}`).join(' '));

    const res = await quiet(() => s.solveModes(f, 6, null, { maxRefineIters: 1 }));
    check('εr=4.4 box: solve succeeded', !res.error, res.error || `nconv=${res.nconv}, nTris=${res.nTris}`);
    const prop = (res.modes || []).filter(m => m.status === 'propagating');
    console.log('  found propagating: ' + prop.map(m => m.eps_eff.toFixed(4)).join(' '));

    let qi = -1, qOvl = -1;
    prop.forEach((m, i) => { if (m.overlap > qOvl) { qOvl = m.overlap; qi = i; } });
    const qtem = prop[qi];
    check('εr=4.4 box: quasi-TEM present with 1 < ε_eff < 4.4',
        qtem && qtem.eps_eff > 1 && qtem.eps_eff < 4.4, qtem ? `eps=${qtem.eps_eff.toFixed(4)}` : 'missing');

    // THE regression assertion: the far-from-shift cavity modes must be found
    const cavity = prop.filter((m, i) => i !== qi && m.eps_eff < 1);
    check('εr=4.4 box: cavity modes not dropped (regression 5a34928)',
        cavity.length >= airRef.length, `found ${cavity.length}, air-box minimum ${airRef.length}`);

    // each air-box mode has a loaded counterpart: ε_eff in [air value − slack, 1)
    const pool = cavity.slice();
    let allBounded = true;
    for (const e of airRef) {
        const i = pool.findIndex(m => m.eps_eff > e.eps - 0.02 && m.eps_eff < e.eps + 0.2);
        if (i >= 0) pool.splice(i, 1); else allBounded = false;
    }
    check('εr=4.4 box: each cavity mode within [air bound, air bound + loading]', allBounded);

    // classification sanity: nothing propagating above the densest medium
    check('εr=4.4 box: no propagating mode with ε_eff > εr',
        prop.every(m => m.eps_eff < 4.4 * 1.02));
}

console.log(failures === 0 ? '\nBOX MODES OK' : `\nBOX MODES: ${failures} FAILURE(S)`);
process.exit(failures === 0 ? 0 : 1);
