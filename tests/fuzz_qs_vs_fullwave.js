// Fuzzer: quasi-static (rectilinear FDM) vs full-wave (triangular FEM).
//
// Generates random transmission-line geometries across ALL UI line types and advanced
// features, solves each with BOTH backends, and flags:
//   • large discrepancies in the field-driven quantities (Z0, eps_eff, C) the two
//     methods should agree on,
//   • geometries whose triangular mesh comes out poor (sliver Qmax, degenerate, NaN),
//   • geometries that make ONE backend throw while the other succeeds (a backend bug).
//
// Coverage:
//   line types — microstrip, diff_microstrip, stripline, diff_stripline,
//                gcpw, diff_gcpw, broadside_stripline
//   features   — solder mask, top dielectric, ground cutout, enclosure (side/top
//                ground walls), surface plating, surface roughness, differential gaps
//
// Loss (R / alpha_c) is intentionally NOT compared — the two methods model conductor
// loss differently (full-wave MQS captures proximity the quasi-static perturbation
// under-predicts), so they legitimately diverge there. Geometries that BOTH backends
// reject (invalid random combos) are skipped, not flagged.
//
// Usage:   node tests/fuzz_qs_vs_fullwave.js [N] [seed] [threshold%]
//   e.g.   node tests/fuzz_qs_vs_fullwave.js 40 1 15
// Reproduce a flagged case: rerun with the same N and seed (generation is deterministic).

import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';

const N      = parseInt(process.argv[2]) || 40;
const SEED   = process.argv[3] !== undefined ? parseInt(process.argv[3]) : 1;
const THRESH = (parseFloat(process.argv[4]) || 15) / 100;
const BAD_Q  = 100;

// --- Seeded PRNG (deterministic, reproducible) ---
function createRng(seed) {
    let s = (seed | 0) || 1;
    const next = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; };
    return {
        f: (lo, hi) => lo + next() * (hi - lo),
        logf: (lo, hi) => lo * Math.pow(hi / lo, next()),
        pick: arr => arr[Math.floor(next() * arr.length)],
        bool: p => next() < p,
    };
}

const TL_TYPES = ['microstrip', 'diff_microstrip', 'stripline', 'diff_stripline',
                  'gcpw', 'diff_gcpw', 'broadside_stripline'];

// --- Random spec generator (geometry + which features are enabled) ---
function randomSpec(rng) {
    const tl = rng.pick(TL_TYPES);
    const spec = { tl, freq: rng.logf(0.5e9, 12e9), rq: rng.bool(0.3) ? rng.logf(0.1e-6, 1e-6) : 0 };

    if (tl === 'broadside_stripline') {
        Object.assign(spec, {
            bs_w: rng.logf(0.05e-3, 1e-3), bs_t: rng.logf(10e-6, 50e-6),
            bs_x_offset: rng.bool(0.5) ? rng.f(-0.15e-3, 0.15e-3) : 0,
            bs_h_bottom: rng.logf(0.05e-3, 0.5e-3), bs_er_bottom: rng.f(2.2, 6), bs_tand_bottom: rng.f(0.001, 0.02),
            bs_h_middle: rng.logf(0.05e-3, 0.5e-3), bs_er_middle: rng.f(2.2, 6), bs_tand_middle: rng.f(0.001, 0.02),
            bs_h_top: rng.logf(0.05e-3, 0.5e-3), bs_er_top: rng.f(2.2, 6), bs_tand_top: rng.f(0.001, 0.02),
            bs_sigma: 5.8e7,
        });
        spec.use_enclosure = rng.bool(0.3);
        spec.use_side_gnd = rng.bool(0.6);
        if (spec.use_enclosure) spec.enclosure_width = spec.bs_w * rng.f(3, 8);
        spec.use_plating = rng.bool(0.2);
        return spec;
    }

    const w = rng.logf(0.05e-3, 3e-3), h = rng.logf(0.05e-3, 1.5e-3), t = rng.logf(10e-6, 70e-6);
    Object.assign(spec, {
        w, h, t, er: rng.f(2.2, 10), tand: rng.f(0.001, 0.02), sigma: 5.8e7,
        gnd_thickness: rng.logf(10e-6, 35e-6),
    });
    const isDiff = tl.includes('diff');
    if (isDiff) spec.trace_spacing = rng.logf(0.02e-3, 0.6e-3);
    if (tl.includes('stripline')) { spec.er_top = spec.er; spec.tand_top = spec.tand; spec.stripline_top_h = rng.logf(0.05e-3, 1.5e-3); }
    if (tl.includes('gcpw')) { spec.gap = rng.logf(0.05e-3, 0.5e-3); spec.via_gap = rng.logf(0.05e-3, 0.5e-3); }

    // Advanced features (solder mask / top dielectric / ground cutout only apply to
    // the plain microstrip family in the app).
    if (tl === 'microstrip' || tl === 'diff_microstrip') {
        spec.use_sm = rng.bool(0.3);
        if (rng.bool(0.3)) { spec.use_top_diel = true; spec.top_diel_h = rng.logf(0.05e-3, 0.5e-3); }
        if (rng.bool(0.2)) { spec.use_gnd_cut = true; spec.gnd_cut_w = w * rng.f(0.2, 1.5); spec.gnd_cut_h = h * rng.f(0.1, 0.7); }
    }
    // Enclosure walls apply to any of these.
    if (rng.bool(0.3)) {
        spec.use_enclosure = true;
        const span = isDiff ? (2 * w + spec.trace_spacing) : w;
        spec.enclosure_width = span * rng.f(2.5, 6) + 2 * spec.gnd_thickness;
        spec.enclosure_height = (h + t) * rng.f(1.5, 4);
        spec.use_side_gnd = rng.bool(0.6);
        spec.use_top_gnd = rng.bool(0.4);
    }
    spec.use_plating = rng.bool(0.2);
    return spec;
}

// --- Build solver options for a spec on a backend (mirrors app_solver.updateGeometry) ---
function addCommon(o, spec) {
    if (spec.use_sm) { o.use_sm = true; o.sm_t_sub = 20e-6; o.sm_t_trace = 20e-6; o.sm_t_side = 20e-6; o.sm_er = 3.5; o.sm_tand = 0.02; }
    if (spec.use_top_diel) { o.top_diel_h = spec.top_diel_h; o.top_diel_er = 4.5; o.top_diel_tand = 0.02; }
    if (spec.use_gnd_cut) { o.gnd_cut_width = spec.gnd_cut_w; o.gnd_cut_sub_h = spec.gnd_cut_h; }
    if (spec.use_enclosure) {
        o.enclosure_width = spec.enclosure_width;
        if (o.enclosure_height === undefined) o.enclosure_height = spec.enclosure_height;
        const lr = spec.use_side_gnd ? 'gnd' : 'open';
        const top = spec.use_top_gnd ? 'gnd' : (o.boundaries ? o.boundaries[2] : 'open');
        const bot = o.boundaries ? o.boundaries[3] : 'gnd';
        o.boundaries = [lr, lr, top, bot];
    }
    if (spec.use_plating) o.plating = { sigma: 1e7, thickness: 4e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: true };
}

function buildSolver(spec, backend) {
    if (spec.tl === 'broadside_stripline') {
        const o = {
            trace_width: spec.bs_w, trace_thickness: spec.bs_t, x_offset: spec.bs_x_offset, sigma_cond: spec.bs_sigma,
            h_bottom: spec.bs_h_bottom, er_bottom: spec.bs_er_bottom, tand_bottom: spec.bs_tand_bottom,
            h_middle: spec.bs_h_middle, er_middle: spec.bs_er_middle, tand_middle: spec.bs_tand_middle,
            h_top: spec.bs_h_top, er_top: spec.bs_er_top, tand_top: spec.bs_tand_top,
            freq: spec.freq, rq: spec.rq, boundaries: ['open', 'open', 'gnd', 'gnd'], mesh_backend: backend,
        };
        if (spec.use_enclosure) { o.enclosure_width = spec.enclosure_width; if (spec.use_side_gnd) o.boundaries = ['gnd', 'gnd', 'gnd', 'gnd']; }
        if (spec.use_plating) o.plating = { sigma: 1e7, thickness: 4e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: true };
        const s = new BroadsideStriplineSolver(o);
        s.use_causal_materials = false;
        return s;
    }
    const base = {
        trace_width: spec.w, substrate_height: spec.h, trace_thickness: spec.t, gnd_thickness: spec.gnd_thickness,
        epsilon_r: spec.er, tan_delta: spec.tand, sigma_cond: spec.sigma, freq: spec.freq, rq: spec.rq,
        mesh_backend: backend,
    };
    const o = { ...base };
    if (spec.trace_spacing) o.trace_spacing = spec.trace_spacing;
    if (spec.tl.includes('gcpw')) { o.boundaries = ['open', 'open', 'open', 'gnd']; o.use_coplanar_gnd = true; o.gap = spec.gap; o.via_gap = spec.via_gap; o.use_vias = true; }
    else if (spec.tl.includes('stripline')) { o.epsilon_r_top = spec.er_top; o.tan_delta_top = spec.tand_top; o.enclosure_height = spec.stripline_top_h; o.boundaries = ['open', 'open', 'gnd', 'gnd']; }
    else o.boundaries = ['open', 'open', 'open', 'gnd'];
    addCommon(o, spec);
    const s = new MicrostripSolver(o);
    s.use_causal_materials = false;
    return s;
}

async function solveOn(spec, backend) {
    const log = console.log;
    console.log = () => {};
    try {
        const s = buildSolver(spec, backend);
        const r = await s.solve_adaptive();
        const m = r.modes[0];
        return {
            Z0: m.Z0, eps_eff: m.eps_eff, C: m.RLGC.C,
            maxQ: s.meshQuality ? s.meshQuality.maxQ : null,
            degenerate: s.meshQuality ? s.meshQuality.degenerateCount : 0,
            nan: s.meshQuality ? s.meshQuality.nanNodes : 0,
        };
    } catch (e) {
        return { error: e.message || String(e) };
    } finally {
        console.log = log;
    }
}

function fmtSpec(spec) {
    const u = (x, s = 1e3, d = 3) => (x * s).toFixed(d);
    let s = spec.tl;
    if (spec.tl === 'broadside_stripline') s += ` w=${u(spec.bs_w)}mm t=${u(spec.bs_t, 1e6, 1)}µm off=${u(spec.bs_x_offset)}mm hM=${u(spec.bs_h_middle)}mm`;
    else s += ` w=${u(spec.w)}mm h=${u(spec.h)}mm t=${u(spec.t, 1e6, 1)}µm er=${spec.er.toFixed(2)}`;
    s += ` f=${(spec.freq / 1e9).toFixed(2)}GHz`;
    if (spec.trace_spacing) s += ` gap=${u(spec.trace_spacing)}mm`;
    const feats = ['use_sm', 'use_top_diel', 'use_gnd_cut', 'use_enclosure', 'use_side_gnd', 'use_top_gnd', 'use_plating']
        .filter(k => spec[k]).map(k => k.replace('use_', ''));
    if (feats.length) s += ` [${feats.join(',')}]`;
    return s;
}

const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);

async function main() {
    console.log(`\n### Fuzz QS vs full-wave — N=${N} seed=${SEED} flag>${(THRESH * 100).toFixed(0)}% Qmax>${BAD_Q} ###`);
    console.log(`covers: ${TL_TYPES.join(', ')} + solder-mask/top-diel/gnd-cut/enclosure/plating/roughness\n`);
    const rng = createRng(SEED);
    const flagged = { discrepancy: [], badMesh: [], crash: [] };
    let ok = 0, skipped = 0;
    const byType = {};

    for (let i = 0; i < N; i++) {
        const spec = randomSpec(rng);
        byType[spec.tl] = (byType[spec.tl] || 0) + 1;
        const qs = await solveOn(spec, 'rectilinear');
        const fw = await solveOn(spec, 'triangular');

        // Both reject → invalid random geometry (validation is backend-independent); skip.
        if (qs.error && fw.error) { skipped++; continue; }
        // Only one backend fails → a real backend bug.
        if (qs.error || fw.error) {
            const who = qs.error ? `quasi-static threw (full-wave OK): ${qs.error}` : `full-wave threw (quasi-static OK): ${fw.error}`;
            flagged.crash.push({ spec, who });
            console.log(`[${i}] ✗ ONE-SIDED FAILURE — ${who}\n      ${fmtSpec(spec)}`);
            continue;
        }

        const badMesh = (fw.maxQ != null && fw.maxQ > BAD_Q) || fw.degenerate > 0 || fw.nan > 0;
        if (badMesh) {
            flagged.badMesh.push({ spec, maxQ: fw.maxQ });
            console.log(`[${i}] ✗ BAD MESH (Qmax=${fw.maxQ?.toFixed(0)} degen=${fw.degenerate} nan=${fw.nan})\n      ${fmtSpec(spec)}`);
        }

        const dZ = relDiff(qs.Z0, fw.Z0), dE = relDiff(qs.eps_eff, fw.eps_eff), dC = relDiff(qs.C, fw.C);
        const dMax = Math.max(dZ, dE, dC);
        if (dMax > THRESH) {
            flagged.discrepancy.push({ spec, dZ, dE, dC });
            console.log(`[${i}] ✗ DISCREPANCY ${(dMax * 100).toFixed(0)}% (Z0 ${(dZ * 100).toFixed(0)}% eps ${(dE * 100).toFixed(0)}% C ${(dC * 100).toFixed(0)}%)\n` +
                `      ${fmtSpec(spec)}\n` +
                `      qs: Z0=${qs.Z0.toFixed(2)} eps=${qs.eps_eff.toFixed(3)} C=${(qs.C * 1e12).toFixed(1)}pF | ` +
                `fw: Z0=${fw.Z0.toFixed(2)} eps=${fw.eps_eff.toFixed(3)} C=${(fw.C * 1e12).toFixed(1)}pF`);
        }
        if (!badMesh && dMax <= THRESH) ok++;
    }

    console.log(`\n${'='.repeat(72)}`);
    console.log(`coverage: ${Object.entries(byType).map(([k, v]) => `${k}:${v}`).join('  ')}`);
    console.log(`SUMMARY: ${ok}/${N} clean | ${flagged.discrepancy.length} discrepancy | ` +
        `${flagged.badMesh.length} bad-mesh | ${flagged.crash.length} one-sided-fail | ${skipped} skipped(invalid)`);
    console.log(`${'='.repeat(72)}`);
}

main();
