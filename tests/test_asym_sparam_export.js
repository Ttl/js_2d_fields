// test_asym_sparam_export.js — asymmetric 4-port S-parameter export for the broadside pair.
//
// A broadside-coupled line with a THICK air middle layer and UNEQUAL top/bottom substrates
// (different thickness AND permittivity) is two weakly-coupled, dissimilar single lines:
//   - lower trace: microstrip over the bottom ground through er_bottom (substrate h_bottom),
//   - upper trace: the mirror microstrip over the top ground through er_top (substrate h_top).
// The air gap is thick, so coupling is low and the exported 4-port S-matrix must be
// block-diagonal — each 2-port diagonal block reproducing the corresponding standalone
// microstrip's exported 2-port S-parameters. This exercises the asymmetric export path
// (computeSParamsDiffAuto with the genuine physical 2×2 [C]/[L], giving S11≠S22), which the
// symmetric odd/even combination cannot represent.
//
// Properties checked (per frequency, on both backends):
//   1. EXPORT SHAPE   — generateS4P/generateS2P produce parseable 4-/2-port Touchstone files.
//   2. DECOUPLING     — the two diagonal 2-port blocks each match one of the two microstrips
//                       (distinct assignment), and cross-trace terms (isolation) are small.
//   3. x_offset       — laterally offsetting the upper trace leaves the 4-port S unchanged
//                       (with the gap thick the traces don't "see" each other's x position).
//
// Backend: `MESH_BACKEND=triangular node tests/test_asym_sparam_export.js`; default rectilinear.

import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { generateS4P, generateS2P } from '../src/snp_export.js';

const MESH_BACKEND = process.env.MESH_BACKEND || 'rectilinear';

// --- Geometry (asymmetric, thick air middle) --------------------------------
const w = 0.2e-3, t = 35e-6, sigma = 5.8e7;
const h_bottom = 0.2e-3, er_bottom = 4.4, tand_bottom = 0.02;   // bottom line
const h_middle = 1.0e-3, er_middle = 1.0, tand_middle = 0.0;    // thick AIR gap
const h_top = 0.12e-3, er_top = 3.0, tand_top = 0.01;           // top line — different h AND er
const length = 10e-3, Z_ref = 50;
const FREQS = [1e9, 5e9];
const X_OFFSET = 0.15e-3;

// Tolerances (empirically ~0.018 match / ~0.029 isolation / ~0.002 x-offset at 5 GHz).
const TOL_MATCH = 0.035;   // max |ΔS| between a diagonal block and its microstrip
const TOL_ISO = 0.045;     // max cross-trace |S| (low-coupling bound)
const TOL_XOFF = 0.01;     // max |ΔS| of the 4-port under a lateral offset
const SOLVE_OPTS = { max_iters: 6, max_nodes: 20000, param_tol: 0.03 };

function broadside(x_offset) {
    return new BroadsideStriplineSolver({
        trace_width: w, trace_thickness: t, sigma_cond: sigma,
        h_bottom, er_bottom, tand_bottom,
        h_middle, er_middle, tand_middle,
        h_top, er_top, tand_top,
        x_offset, freq: FREQS[0], nx: 60, ny: 60,
        enclosure_width: 3e-3, boundaries: ['gnd', 'gnd', 'gnd', 'gnd'],
        mesh_backend: MESH_BACKEND,
    });
}

// Equivalent microstrip: substrate (h, er) over its ground, air above, the opposing ground as a
// cover at the same height the broadside trace sees it (h + h_middle + h_other).
function microstrip(h, er, tand, cover_h) {
    return new MicrostripSolver({
        substrate_height: h, trace_width: w, trace_thickness: t, epsilon_r: er, tan_delta: tand,
        sigma_cond: sigma, epsilon_r_top: 1.0, freq: FREQS[0], nx: 60, ny: 60,
        enclosure_width: 3e-3, enclosure_height: cover_h,
        boundaries: ['gnd', 'gnd', 'gnd', 'gnd'], mesh_backend: MESH_BACKEND,
    });
}

// --- Solve helpers ----------------------------------------------------------
const mute = () => { const l = console.log, w = console.warn; console.log = () => {}; console.warn = () => {}; return () => { console.log = l; console.warn = w; }; };
async function solveSweep(solver) {
    // Solve once (adaptive mesh at FREQS[0]) then reuse the mesh for the rest of the sweep.
    const un = mute();
    try {
        const base = await solver.solve_adaptive(SOLVE_OPTS);
        const out = [{ freq: FREQS[0], result: base }];
        for (const f of FREQS.slice(1)) out.push({ freq: f, result: await solver.computeAtFrequency(f, base) });
        return out;
    } finally { un(); }
}

// --- Touchstone (MA) parsing ------------------------------------------------
// Strip comment (!)/option (#) lines, tokenize the rest, chunk into per-frequency records of
// `1 + 2·nPorts²` numbers (freq, then row-major (mag, angle) pairs). Returns [{freq, S}] with S
// an nPorts×nPorts array of {re, im}.
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

const cabs = c => Math.hypot(c.re, c.im);
const cdiff = (a, b) => Math.hypot(a.re - b.re, a.im - b.im);

// A 2-port block from the 4-port S at port indices (a, b). Port map: 1=in+,2=in-,3=out+,4=out-,
// so trace+ = ports {1,3} = indices {0,2} and trace- = ports {2,4} = indices {1,3}.
const block = (S, a, b) => ({ S11: S[a][a], S21: S[b][a], S12: S[a][b], S22: S[b][b] });
const blockDist = (blk, ms) => Math.max(cdiff(blk.S11, ms.S11), cdiff(blk.S21, ms.S21),
    cdiff(blk.S12, ms.S12), cdiff(blk.S22, ms.S22));

// --- Test -------------------------------------------------------------------
async function run() {
    let failed = 0;
    const check = (name, cond, detail = '') => { console.log(`  ${cond ? '✓' : '✗'} ${name}${detail && !cond ? ` — ${detail}` : ''}`); if (!cond) failed++; };

    console.log(`\n### Asymmetric S-parameter export (mesh_backend = "${MESH_BACKEND}") ###`);

    // Broadside 4-port export (x_offset = 0 and offset), plus the two microstrip 2-port exports.
    const bs = await solveSweep(broadside(0));
    const bsOff = await solveSweep(broadside(X_OFFSET));
    const msBot = await solveSweep(microstrip(h_bottom, er_bottom, tand_bottom, h_bottom + h_middle + h_top));
    const msTop = await solveSweep(microstrip(h_top, er_top, tand_top, h_top + h_middle + h_bottom));

    // Sanity: the differential result really carries the asymmetric physical matrix.
    check('broadside result has odd/even modes + physMatrix',
        !!bs[0].result.physMatrix && !!bs[0].result.modes.find(m => m.mode === 'odd') && !!bs[0].result.modes.find(m => m.mode === 'even'));

    // Export through the real Touchstone writers, then parse back.
    const s4 = parseTouchstone(generateS4P(bs, length, Z_ref), 4);
    const s4off = parseTouchstone(generateS4P(bsOff, length, Z_ref), 4);
    const s2bot = parseTouchstone(generateS2P(msBot, length, Z_ref), 2);
    const s2top = parseTouchstone(generateS2P(msTop, length, Z_ref), 2);
    check('generateS4P → 4-port, one record per frequency', s4.length === FREQS.length && s4[0].S.length === 4);
    check('generateS2P → 2-port, one record per frequency', s2bot.length === FREQS.length && s2bot[0].S.length === 2);

    for (let fi = 0; fi < FREQS.length; fi++) {
        const fGHz = FREQS[fi] / 1e9;
        const S = s4[fi].S;
        const msB = { S11: s2bot[fi].S[0][0], S21: s2bot[fi].S[1][0], S12: s2bot[fi].S[0][1], S22: s2bot[fi].S[1][1] };
        const msT = { S11: s2top[fi].S[0][0], S21: s2top[fi].S[1][0], S12: s2top[fi].S[0][1], S22: s2top[fi].S[1][1] };

        // (2) Decoupling: assign each diagonal block to its nearest microstrip; require a distinct
        //     assignment and both distances within tolerance.
        const b1 = block(S, 0, 2), b2 = block(S, 1, 3);
        const d1B = blockDist(b1, msB), d1T = blockDist(b1, msT);
        const d2B = blockDist(b2, msB), d2T = blockDist(b2, msT);
        const b1IsTop = d1T < d1B;                       // block-1 nearer the top microstrip?
        const dist1 = b1IsTop ? d1T : d1B;
        const dist2 = b1IsTop ? d2B : d2T;               // block-2 must take the OTHER microstrip
        check(`${fGHz} GHz: blocks match distinct microstrips`, b1IsTop !== (d2T < d2B),
            `both blocks matched the same line (d1B=${d1B.toFixed(3)} d1T=${d1T.toFixed(3)} d2B=${d2B.toFixed(3)} d2T=${d2T.toFixed(3)})`);
        check(`${fGHz} GHz: block↔microstrip |ΔS| < ${TOL_MATCH}`, Math.max(dist1, dist2) < TOL_MATCH,
            `max block distance ${Math.max(dist1, dist2).toFixed(4)}`);

        // (2) Isolation: cross-trace terms between {ports 1,3} and {ports 2,4}.
        const iso = Math.max(...[[0, 1], [0, 3], [2, 1], [2, 3]].map(([i, j]) => cabs(S[i][j])));
        check(`${fGHz} GHz: cross-trace isolation |S| < ${TOL_ISO}`, iso < TOL_ISO, `max cross |S| ${iso.toFixed(4)}`);

        // (3) x_offset invariance of the full 4-port.
        let xd = 0;
        for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) xd = Math.max(xd, cdiff(S[i][j], s4off[fi].S[i][j]));
        check(`${fGHz} GHz: x_offset leaves 4-port S unchanged (|ΔS| < ${TOL_XOFF})`, xd < TOL_XOFF, `max |ΔS| ${xd.toFixed(4)}`);
    }

    console.log(failed ? `\nSUMMARY [${MESH_BACKEND}]: ${failed} check(s) failed` : `\nSUMMARY [${MESH_BACKEND}]: all checks passed`);
    if (failed) process.exitCode = 1;
}

run().catch(e => { console.error(e); process.exit(1); });
