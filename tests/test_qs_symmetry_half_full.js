// Half-domain ≡ full-domain identity for the QS rectilinear FDM backend.
//
// Mirror-symmetric geometries auto-solve on the right half (x >= 0) with a
// magnetic wall (single/even) or electric wall (odd) at x = 0; {symmetry: false}
// forces the legacy full domain. This file demands the two agree, mirroring
// tests/test_symmetry_half_full.js (the tri-backend identity suite) with the
// same five feature-stacked families:
//   - single-ended: the plane cuts the TRACE itself (straddle x2 charge scaling,
//     plane-column contour/loss segments, plating conductor_id alignment)
//   - GCPW: pours + via slabs clipped at the mesher
//   - diff stripline (PEC vs PMC plane walls), enclosure (side wall conductors
//     clipped), ground cutout centered on the plane
//
// The strongest pin is the SAME-GRID identity: the full-domain solve on the
// exact mirror of a converged half grid must reproduce C, C0, R, G to solver
// roundoff — the half-domain path is the exact restriction of the calibrated
// full-domain discretization (PMC row doubling, centered plane quadrature),
// not merely a converging approximation of it.
import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
const relDiff = (a, b) => Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-30);
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};

const FAMILIES = [
    {
        name: 'microstrip + sm + plating', diff: false, cTol: 0.01, rTol: 0.03,
        geom: {
            trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
            epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
            boundaries: ['open', 'open', 'open', 'gnd'],
            use_sm: true, sm_t_sub: 20e-6, sm_t_trace: 20e-6, sm_t_side: 20e-6, sm_er: 3.5, sm_tand: 0.02,
            plating: { sigma: 1e7, thickness: 4e-6, rq: 0, top: true, sides: true, bottom: false, thick_corners: true },
        },
    },
    {
        name: 'gcpw + vias + sm', diff: false, cTol: 0.01, rTol: 0.03,
        geom: {
            trace_width: 0.3e-3, substrate_height: 0.254e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
            epsilon_r: 3.66, tan_delta: 0.003, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
            boundaries: ['open', 'open', 'open', 'gnd'],
            use_coplanar_gnd: true, gap: 0.15e-3, via_gap: 0.3e-3, use_vias: true,
            use_sm: true, sm_t_sub: 15e-6, sm_t_trace: 15e-6, sm_t_side: 15e-6, sm_er: 3.4, sm_tand: 0.02,
        },
    },
    {
        name: 'diff stripline', diff: true, cTol: 0.01, rTol: 0.08,
        geom: {
            trace_width: 0.15e-3, trace_spacing: 0.15e-3, substrate_height: 0.2e-3, trace_thickness: 17e-6,
            gnd_thickness: 17e-6, epsilon_r: 4.1, epsilon_r_top: 4.1, tan_delta: 0.02, tan_delta_top: 0.02,
            enclosure_height: 0.4e-3, sigma_cond: 5.8e7, freq: 2e9, rq: 0,
            boundaries: ['open', 'open', 'gnd', 'gnd'],
        },
    },
    {
        name: 'diff microstrip + enclosure', diff: true, cTol: 0.01, rTol: 0.08,
        geom: {
            trace_width: 0.2e-3, trace_spacing: 0.2e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
            gnd_thickness: 20e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
            enclosure_width: 2.4e-3, enclosure_height: 0.8e-3,
            boundaries: ['gnd', 'gnd', 'gnd', 'gnd'],
        },
    },
    {
        name: 'diff microstrip + gnd cutout', diff: true, cTol: 0.03, rTol: 0.08,
        geom: {
            trace_width: 0.2e-3, trace_spacing: 0.15e-3, substrate_height: 0.25e-3, trace_thickness: 35e-6,
            gnd_thickness: 20e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9, rq: 0,
            boundaries: ['open', 'open', 'open', 'gnd'],
            gnd_cut_width: 0.25e-3, gnd_cut_sub_h: 0.1e-3,
        },
    },
];

const SOLVE_OPTS = { max_iters: 12, energy_tol: 0.01, param_tol: 0.05, max_nodes: 30000, min_converged_passes: 2 };

async function solve(geom, extra = {}) {
    const s = new MicrostripSolver({ ...geom, nx: 30, ny: 30, ...extra });
    const r = await quiet(() => s.solve_adaptive({ ...SOLVE_OPTS }));
    return { s, r };
}

for (const fam of FAMILIES) {
    const half = await solve(fam.geom);
    const full = await solve(fam.geom, { symmetry: false });
    check(`${fam.name}: half domain auto-detected`, half.s.sym_half === true && full.s.sym_half === false);
    for (let m = 0; m < half.r.modes.length; m++) {
        const h = half.r.modes[m], f = full.r.modes[m];
        const tag = `${fam.name} [${h.mode}]`;
        check(`${tag}: half C == full C (±${100 * fam.cTol}%)`, relDiff(h.RLGC.C, f.RLGC.C) < fam.cTol,
            `${(h.RLGC.C * 1e12).toFixed(2)} vs ${(f.RLGC.C * 1e12).toFixed(2)} pF/m`);
        check(`${tag}: half Z0 == full Z0 (±${100 * fam.cTol}%)`, relDiff(h.Z0, f.Z0) < fam.cTol,
            `${h.Z0.toFixed(2)} vs ${f.Z0.toFixed(2)} Ω`);
        check(`${tag}: half R == full R (±${100 * fam.rTol}%)`, relDiff(h.RLGC.R, f.RLGC.R) < fam.rTol,
            `${h.RLGC.R.toFixed(2)} vs ${f.RLGC.R.toFixed(2)} Ω/m`);
        check(`${tag}: half alpha_d == full alpha_d (±3%)`, relDiff(h.alpha_d, f.alpha_d) < 0.03,
            `${h.alpha_d.toFixed(3)} vs ${f.alpha_d.toFixed(3)} dB/m`);
    }
    if (fam.diff) {
        check(`${fam.name}: Z_diff (±${100 * fam.cTol}%)`, relDiff(half.r.Z_diff, full.r.Z_diff) < fam.cTol,
            `${half.r.Z_diff.toFixed(2)} vs ${full.r.Z_diff.toFixed(2)} Ω`);
        check(`${fam.name}: Z_common (±${100 * fam.cTol}%)`, relDiff(half.r.Z_common, full.r.Z_common) < fam.cTol,
            `${half.r.Z_common.toFixed(2)} vs ${full.r.Z_common.toFixed(2)} Ω`);
    }
}

// ---- Same-grid exact identity: full solve on the mirror of the half grid ----
// Different matrices and elimination orders, so "exact" means solver roundoff,
// far below any physical tolerance. Two families pin the machinery from both
// sides: the straddling plated microstrip (plane cuts the trace: interior
// epsilon rule on the plane column, two-segment plating quadrature) and the
// diff stripline (PEC vs PMC plane walls).
for (const fam of [FAMILIES[0], FAMILIES[2]]) {
    const { s: sh, r: rh } = await solve(fam.geom);
    const sf = new MicrostripSolver({ ...fam.geom, nx: 30, ny: 30, symmetry: false });
    const xh = sh.x, n = xh.length;
    const xf = new Float64Array(2 * n - 1);
    for (let j = 0; j < n - 1; j++) xf[j] = -xh[n - 1 - j];
    for (let j = 0; j < n; j++) xf[n - 1 + j] = xh[j];
    sf.x = xf;
    sf.y = Float64Array.from(sh.y);
    sf.dx = new Float64Array(xf.length - 1);
    for (let j = 0; j < xf.length - 1; j++) sf.dx[j] = xf[j + 1] - xf[j];
    sf.dy = new Float64Array(sf.y.length - 1);
    for (let i = 0; i < sf.y.length - 1; i++) sf.dy[i] = sf.y[i + 1] - sf.y[i];
    sf._setup_geometry();
    sf.mesh_generated = true;
    const rf = await quiet(() => sf.solve_adaptive({ skip_mesh: true }));
    for (let m = 0; m < rh.modes.length; m++) {
        const h = rh.modes[m], f = rf.modes[m];
        for (const [k, hv, fv] of [['C', h.C, f.C], ['C0', h.C0, f.C0],
                                   ['R', h.RLGC.R, f.RLGC.R], ['G', h.RLGC.G, f.RLGC.G]]) {
            check(`${fam.name}: same-grid identity [${h.mode}] ${k} (±0.01%)`, relDiff(hv, fv) < 1e-4,
                `${hv.toExponential(8)} vs ${fv.toExponential(8)}`);
        }
    }

    // Fast sweep path reuses the cached half-domain fields — same scaling.
    const rh2 = await quiet(() => sh.computeAtFrequency(2 * fam.geom.freq, rh));
    const rf2 = await quiet(() => sf.computeAtFrequency(2 * fam.geom.freq, rf));
    check(`${fam.name}: same-grid identity holds at a second sweep frequency (±0.01%)`,
        rh2.modes.every((m2, k) => relDiff(m2.RLGC.R, rf2.modes[k].RLGC.R) < 1e-4));

    // Plot mirror: full-width grid with the correct per-mode parity.
    const pf = sh.getPlotFields();
    const nf = pf.x.length;
    const iy = Math.floor(pf.V[0].length * 0.5);
    const jq = Math.floor((nf - 1) / 2 + (nf - 1) / 8), jm = (nf - 1) - jq;
    check(`${fam.name}: getPlotFields mirrors the grid to the full width`,
        nf === 2 * n - 1 && pf.x[0] === -xh[n - 1] && pf.x[(nf - 1) / 2] === 0);
    check(`${fam.name}: solver-internal fields stay half-domain after getPlotFields`,
        sh.x.length === n);
    if (fam.diff) {
        // Odd mode: electric wall — the plane column must be exactly 0.
        let maxPlane = 0;
        for (let i = 0; i < sh.V[0].length; i++) maxPlane = Math.max(maxPlane, Math.abs(sh.V[0][i][0]));
        check('odd-mode V == 0 on the symmetry plane column', maxPlane === 0);
        check('mirrored odd V is antisymmetric, even V symmetric',
            pf.V[0][iy][jm] === -pf.V[0][iy][jq] && pf.V[1][iy][jm] === pf.V[1][iy][jq]);
    } else {
        check('mirrored single-mode V is symmetric', pf.V[0][iy][jm] === pf.V[0][iy][jq]);
    }
}

// ---- Detection / veto assertions ----
{
    const geom = FAMILIES[0].geom;
    check('opt-out {symmetry: false} disables the half domain',
        new MicrostripSolver({ ...geom, nx: 30, ny: 30, symmetry: false }).sym_half === false);
    const bs = new BroadsideStriplineSolver({
        trace_width: 0.2e-3, trace_thickness: 35e-6,
        h_bottom: 0.2e-3, h_middle: 0.2e-3, h_top: 0.2e-3,
        er_bottom: 4.4, er_middle: 4.4, er_top: 4.4,
        sigma_cond: 5.8e7, freq: 1e9,
    });
    check('broadside pair never takes the half domain', !bs.sym_half);
    const single = new MicrostripSolver({ ...geom, nx: 30, ny: 30 });
    check('single-ended straddling trace flagged for x2 charge scaling',
        single.sym_half === true && single._sym_signal_straddles === true);
    const pair = new MicrostripSolver({ ...FAMILIES[2].geom, nx: 30, ny: 30 });
    check('differential pair: no straddle scaling (full per-trace contour)',
        pair.sym_half === true && pair._sym_signal_straddles === false);
}

console.log(failures ? `\n${failures} FAILURE(S)` : '\nALL PASS');
process.exit(failures ? 1 : 0);
