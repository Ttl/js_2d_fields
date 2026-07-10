// test_interp_asym_sparams.js — asymmetric S-parameters through the INTERPOLATING sweep.
//
// The interpolating sweep (the app default) replaces the discrete per-frequency results
// wholesale with InterpolatingSweep.buildResults() output. For an asymmetric pair those
// entries must carry the physical 2×2 {C, L, Tv} (physMatrix) — the export/plot path keys
// on it to use the MTL 4-port instead of the symmetric odd/even combination. A regression
// here is silent: the solver logs "Asymmetric pair" while the export quietly produces a
// symmetric S22 ≡ S11, SDC ≡ 0 file, and toggling the interpolating-sweep checkbox changes
// the exported S-parameters.
//
// Geometry: the broadside pair from test_asym_sparam_export.js — thick air middle, unequal
// top/bottom substrates, i.e. two decoupled DISSIMILAR lines, so the asymmetry in the
// 4-port is large and unambiguous.
//
// Properties checked (per output frequency):
//   1. CARRIER    — every buildResults() entry has physMatrix and RLGC_matrix.
//   2. MATCH      — s4p exported from interpolated results matches the s4p exported from
//                   discrete computeAtFrequency results at the same frequencies.
//   3. ASYMMETRY  — exported S11 ≠ S22 (the symmetric fallback makes them identical).
//   4. RECIPROCITY — exported 4-port satisfies S = Sᵀ (regression guard for the MTL
//                   ABCD→S operand-order fix, through the full export path).
//
// Backend: `MESH_BACKEND=triangular node tests/test_interp_asym_sparams.js`; default rectilinear.

import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';
import { InterpolatingSweep } from '../src/interpolating_sweep.js';
import { generateS4P } from '../src/snp_export.js';

const MESH_BACKEND = process.env.MESH_BACKEND || 'rectilinear';

// --- Geometry (asymmetric, thick air middle) --------------------------------
const w = 0.2e-3, t = 35e-6, sigma = 5.8e7;
const length = 10e-3, Z_ref = 50;
const F_MIN = 1e9, F_MAX = 1e10;
const OUT_FREQS = [1.5e9, 4e9, 8e9];   // deliberately between sample points

const TOL_MATCH = 0.02;    // interp-vs-discrete |ΔS| (interp tolerance is 1% on RLGC)
const TOL_RECIP = 1e-3;    // reciprocity after Touchstone round-trip
const MIN_ASYM = 0.02;     // |S11 − S22| floor: dissimilar lines reflect differently
const SOLVE_OPTS = { max_iters: 6, max_nodes: 20000, param_tol: 0.03 };

const solver = new BroadsideStriplineSolver({
    trace_width: w, trace_thickness: t, sigma_cond: sigma,
    h_bottom: 0.2e-3, er_bottom: 4.4, tand_bottom: 0.02,
    h_middle: 1.0e-3, er_middle: 1.0, tand_middle: 0.0,   // thick AIR gap
    h_top: 0.12e-3, er_top: 3.0, tand_top: 0.01,          // different h AND er
    x_offset: 0, freq: F_MAX, nx: 60, ny: 60,
    enclosure_width: 3e-3, boundaries: ['gnd', 'gnd', 'gnd', 'gnd'],
    mesh_backend: MESH_BACKEND,
});

const mute = () => { const l = console.log, w2 = console.warn; console.log = () => {}; console.warn = () => {}; return () => { console.log = l; console.warn = w2; }; };

// --- Touchstone (MA) parsing (same as test_asym_sparam_export.js) ------------
function parseTouchstone(text, nPorts) {
    const nums = text.split('\n')
        .filter(l => { const s = l.trim(); return s && !s.startsWith('!') && !s.startsWith('#'); })
        .join(' ').trim().split(/\s+/).map(Number);
    const rec = 1 + 2 * nPorts * nPorts;
    if (nums.length % rec !== 0) throw new Error(`token count ${nums.length} not a multiple of ${rec}`);
    const out = [];
    for (let o = 0; o < nums.length; o += rec) {
        const S = Array.from({ length: nPorts }, () => Array(nPorts));
        for (let i = 0; i < nPorts; i++) for (let j = 0; j < nPorts; j++) {
            const k = o + 1 + (i * nPorts + j) * 2;
            const mag = nums[k], ang = nums[k + 1] * Math.PI / 180;
            S[i][j] = { re: mag * Math.cos(ang), im: mag * Math.sin(ang) };
        }
        out.push({ freq: nums[o], S });
    }
    return out;
}
const cdiff = (a, b) => Math.hypot(a.re - b.re, a.im - b.im);

// --- Test -------------------------------------------------------------------
async function run() {
    let failed = 0;
    const check = (name, cond, detail = '') => { console.log(`  ${cond ? '✓' : '✗'} ${name}${detail && !cond ? ` — ${detail}` : ''}`); if (!cond) failed++; };

    console.log(`\n### Asymmetric S-parameters through interpolating sweep (mesh_backend = "${MESH_BACKEND}") ###`);

    const un = mute();
    let interpResults, discreteResults;
    try {
        const base = await solver.solve_adaptive(SOLVE_OPTS);
        if (!base.physMatrix) throw new Error('geometry did not produce an asymmetric physMatrix — test premise broken');

        const sweep = new InterpolatingSweep(solver, base, { tolerance: 0.01 });
        await sweep.run(F_MIN, F_MAX);
        interpResults = sweep.buildResults(OUT_FREQS);

        discreteResults = [];
        for (const f of OUT_FREQS) discreteResults.push({ freq: f, result: await solver.computeAtFrequency(f, base) });
    } finally { un(); }

    // (1) Every interpolated entry carries the asymmetric physical matrices.
    check('all interpolated entries carry physMatrix ({C, L})',
        interpResults.every(({ result }) => result.physMatrix && result.physMatrix.C && result.physMatrix.L));
    check('all interpolated entries carry RLGC_matrix',
        interpResults.every(({ result }) => result.RLGC_matrix && result.RLGC_matrix.R && result.RLGC_matrix.C));

    // Export both through the real Touchstone writer.
    const s4i = parseTouchstone(generateS4P(interpResults, length, Z_ref), 4);
    const s4d = parseTouchstone(generateS4P(discreteResults, length, Z_ref), 4);

    for (let fi = 0; fi < OUT_FREQS.length; fi++) {
        const fGHz = OUT_FREQS[fi] / 1e9;
        const Si = s4i[fi].S, Sd = s4d[fi].S;

        // (2) Interpolated export matches the discrete export.
        let dm = 0;
        for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) dm = Math.max(dm, cdiff(Si[i][j], Sd[i][j]));
        check(`${fGHz} GHz: interp export matches discrete export (|ΔS| < ${TOL_MATCH})`, dm < TOL_MATCH, `max |ΔS| ${dm.toFixed(4)}`);

        // (3) Asymmetry survives: the two near-end reflections differ. The symmetric
        //     odd/even fallback produces S11 ≡ S22 exactly.
        const asym = cdiff(Si[0][0], Si[1][1]);
        check(`${fGHz} GHz: asymmetry present (|S11 − S22| > ${MIN_ASYM})`, asym > MIN_ASYM, `|S11 − S22| ${asym.toFixed(4)}`);

        // (4) Reciprocity of the exported 4-port.
        let dr = 0;
        for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) dr = Math.max(dr, cdiff(Si[i][j], Si[j][i]));
        check(`${fGHz} GHz: exported 4-port is reciprocal (|S_ij − S_ji| < ${TOL_RECIP})`, dr < TOL_RECIP, `max ${dr.toExponential(2)}`);
    }

    console.log(failed ? `\nSUMMARY [${MESH_BACKEND}]: ${failed} check(s) failed` : `\nSUMMARY [${MESH_BACKEND}]: all checks passed`);
    if (failed) process.exitCode = 1;
}

run().catch(e => { console.error(e); process.exit(1); });
