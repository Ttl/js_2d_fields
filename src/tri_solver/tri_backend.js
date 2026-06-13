// Triangular full-wave FEM backend.
//
// Consumes the SAME backend-neutral geometry the rectilinear FDM solver uses
// (`solver.conductors` / `solver.dielectrics` from `_build_geometry_lists()`),
// meshes it with the generalized geo-kernel mesher, and returns the identical
// result-object shape as field_solver.js (`{ modes:[...], Z_diff?, Z_common?,
// RLGC_matrix? }`), so the app's UI/plotting/export paths are unchanged.
//
// Physics (the "better solver"):
//   - static P2 FEM (C/L/C0/eps_eff/Z0) — variational, gold-standard magnitude
//   - full-wave eigenmode (assembleTriFEM -> generalized eigensolve) — dispersive
//     eps_eff(f) and the mode fields, quasi-TEM mode picked by overlap with the
//     static drive (works for asymmetric geometry; no symmetry assumption)
//   - dielectric loss G from the lossy-permittivity energy integral
//   - conductor loss via the validated perturbation method (corner-corrected),
//     with surface roughness AND plating applied as a surface-impedance scaling
//     (gradient model, calculateZrough): plating = "roughness with a different
//     conductivity" on the plated conductor surfaces.
//
// Mesh + freedom maps are built once and reused across a frequency sweep; only
// the eigen-assembly (k²-dependent) and loss are recomputed per frequency.

import createModule from './eigen_solver.js';
import { createWasmHelpers } from './fem_core.js';
import { initGmsh } from './gmsh_mesh.js';
import { buildGmshMeshFromGeometry, tagMaterials } from './geom_to_mesh.js';
import { buildTriFreedomMap, solveTriStatic, computeTriEnergy, refineTriMesh,
         markTrianglesForRefinement, computeTriP2StaticMatrices,
         staticToEdgeDofs, assembleTriFEM } from './tri_fem.js';
import { staticConductorLoss, solveConductorLoss, computeHtZZMetric } from './conductor_loss.js';
import { mqsConductorLoss, refineSkinBand } from './mqs_loss.js';
import { analyzeTriMode } from './tri_ms_solver.js';

// Below this frequency the full-wave eigensolve degenerates near DC; use the
// static (variational) solve instead. Above it, use the full-wave eigenmode for
// the dispersive effective permittivity.
const F_STATIC_MAX = 100e6;
import { calculateZrough } from './surface_roughness.js';
import { resampleStatic } from './resample.js';
import { Complex } from '../complex.js';

const c0 = 299792458;
const eps0 = 8.854187817e-12;
const MU0 = 4 * Math.PI * 1e-7;
const NP_TO_DB = 8.685889638;

// ---- WASM context (lazy singleton) ----
let _ctxPromise = null;
export function initTriBackend() {
    if (!_ctxPromise) {
        _ctxPromise = (async () => {
            const M = await createModule();
            const helpers = createWasmHelpers(M);
            const G = await initGmsh();
            const wasmSolver = (N, csr, rhsArrays) => helpers.solveSparseMulti(N, csr, rhsArrays);
            return { M, helpers, G, wasmSolver };
        })();
    }
    return _ctxPromise;
}

// Build a loss-edge mask: conductor SURFACE edges + grounded-wall edges.
// (Generalizes conductor_loss.buildMicrostripLossEdges, which hardcodes y=0.)
function buildLossEdges(mesh, fm, condRect) {
    const { nodes, edges, nEdges, nTris, triEdges } = mesh;
    const TOL = 1e-12;
    const isLoss = new Uint8Array(nEdges);
    const touchesExt = new Uint8Array(nEdges);
    for (let t = 0; t < nTris; t++) {
        if (fm.faceF && fm.faceF[2 * t] < 0) continue;  // skip conductor-interior tris
        for (let k = 0; k < 3; k++) touchesExt[triEdges[3 * t + k]] = 1;
    }
    const w = condRect.wallPEC || {};
    const X0 = condRect.xmin_domain, X1 = condRect.xmax_domain;
    const Y0 = condRect.ymin_domain, Y1 = condRect.ymax_domain;
    const symX = condRect.symmetry > 1 ? X0 : null;  // symmetry plane: not a real surface
    const onWall = (x, y) => (
        (w.bottom && Math.abs(y - Y0) < TOL) || (w.top && Math.abs(y - Y1) < TOL) ||
        (w.left && Math.abs(x - X0) < TOL) || (w.right && Math.abs(x - X1) < TOL));
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2 * e], n1 = edges[2 * e + 1];
        const x0 = nodes[2 * n0], y0 = nodes[2 * n0 + 1];
        const x1 = nodes[2 * n1], y1 = nodes[2 * n1 + 1];
        if (symX !== null && Math.abs(x0 - symX) < TOL && Math.abs(x1 - symX) < TOL) continue;
        if (fm.isCondEdge && fm.isCondEdge[e] && touchesExt[e]) { isLoss[e] = 1; continue; }
        if (onWall(x0, y0) && onWall(x1, y1)) isLoss[e] = 1;
    }
    return isLoss;
}

// Effective conductor surface impedance σ/Rq: plating coats the current-carrying
// surface, so when any signal conductor is plated use its plating σ/Rq; else base.
function effectiveSurface(solver) {
    let sigma = solver.sigma_cond ?? 5.8e7;
    let rq = solver.rq ?? 0;
    for (const c of solver.conductors) {
        if (c.is_signal && c.plating && c.plating.sigma) {
            sigma = c.plating.sigma;
            rq = c.plating.rq ?? rq;
            break;
        }
    }
    return { sigma, rq };
}

// Is the geometry mirror-symmetric about x=0? (lets us mesh only the right half.)
function isXSymmetric(conductors, dielectrics, domainW) {
    const tol = domainW * 1e-6;
    const mirrorOf = (list, extra) => {
        // every rect must have a mirror partner (xmin↔-xmax) with same y + material
        const key = (r) => `${extra(r)}|${r.y_min.toFixed(12)}|${r.y_max.toFixed(12)}`;
        const buckets = new Map();
        for (const o of list) {
            const k = key(o);
            if (!buckets.has(k)) buckets.set(k, []);
            buckets.get(k).push([o.x_min, o.x_max]);
        }
        for (const spans of buckets.values()) {
            for (const [xmin, xmax] of spans) {
                const mlo = -xmax, mhi = -xmin;
                const found = spans.some(([a, b]) => Math.abs(a - mlo) < tol && Math.abs(b - mhi) < tol);
                if (!found) return false;
            }
        }
        return true;
    };
    const condOk = mirrorOf(conductors, c => (c.is_signal ? 's' + Math.abs(c.polarity || 0) : 'g'));
    const dielOk = mirrorOf(dielectrics, d => `${d.epsilon_r.toFixed(6)}:${(d.tan_delta || 0).toFixed(6)}`);
    return condOk && dielOk;
}

// Per-triangle electrostatic energy (refinement metric): ε·∫|∇φ|² over each tri.
function perElementEnergy(phi, mesh, epsMap) {
    const { nodes, tris, nTris, triEdges } = mesh;
    const m = new Float64Array(nTris);
    const lp = new Float64Array(6);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const Sz = computeTriP2StaticMatrices(nodes, v0, v1, v2).Sz;
        lp[0] = phi.phiVertex[v0]; lp[1] = phi.phiVertex[v1]; lp[2] = phi.phiVertex[v2];
        lp[3] = phi.phiEdge[triEdges[3 * t]]; lp[4] = phi.phiEdge[triEdges[3 * t + 1]]; lp[5] = phi.phiEdge[triEdges[3 * t + 2]];
        let e = 0;
        for (let li = 0; li < 6; li++) for (let lj = 0; lj < 6; lj++) e += Sz[li * 6 + lj] * lp[li] * lp[lj];
        m[t] = epsMap[t].re * Math.abs(e);
    }
    return m;
}

// Full-wave eigenmode at frequency f via the Lee-Jin eigensolve. Robust mode
// selection: among physical modes (γ²<0) with decent overlap to the static drive,
// pick the one whose eps_mode is closest to eps_static (the quasi-TEM mode).
// Returns { eps, vRe, vIm, g2Re, g2Im } for the chosen mode, or null.
function fullwaveMode(ctx, mesh, fm, abc, condRect, epsMap, f, phiEps, eps_static) {
    const k2 = (2 * Math.PI * f / 299792458) ** 2;
    let fem;
    try { fem = assembleTriFEM(mesh, fm, k2, epsMap, abc, condRect, null); }
    catch { return null; }
    const N = fem.N;
    const seed = staticToEdgeDofs(phiEps, mesh, fm);
    // 4 eigenpairs is ample: shift-invert is centered on the quasi-TEM eigenvalue
    // (-k2·eps_static), so the target mode converges first. (Was 8 — halving the
    // Krylov subspace markedly speeds the eigensolve with no change in the picked mode.)
    const nev = 4, ncv = Math.min(2 * nev + 1, N - 1);
    let res;
    try { res = ctx.helpers.solveGeneralized(N, fem.csrA, fem.csrB, [-k2 * eps_static, 0], nev, ncv, seed); }
    catch { return null; }
    if (!res || !res.nconv) return null;
    let bestIdx = -1, bestD = Infinity;
    for (let i = 0; i < res.nconv; i++) {
        if (res.evalsRe[i] >= 0) continue;                 // physical: γ² < 0
        const vR = res.evecsRe.subarray(i * N, (i + 1) * N);
        let dot = 0, nS = 0, nV = 0;
        for (let k = 0; k < fm.nFreeTransverse; k++) { dot += seed[k] * vR[k]; nS += seed[k] ** 2; nV += vR[k] ** 2; }
        const ovl = nS > 0 && nV > 0 ? Math.abs(dot) / Math.sqrt(nS * nV) : 0;
        if (ovl < 0.3) continue;                           // quasi-TEM overlaps the static drive
        const epsMode = -res.evalsRe[i] / k2;
        globalThis.__TRI_DEBUG__ && console.log(`    cand ${i}: eps=${epsMode.toFixed(4)} ovl=${ovl.toFixed(3)} (target ${eps_static.toFixed(4)})`);
        const d = Math.abs(epsMode - eps_static);
        if (d < bestD) { bestD = d; bestIdx = i; }
    }
    if (bestIdx < 0) return null;
    return {
        eps: -res.evalsRe[bestIdx] / k2,
        vRe: res.evecsRe.slice(bestIdx * N, (bestIdx + 1) * N),
        vIm: res.evecsIm.slice(bestIdx * N, (bestIdx + 1) * N),
        g2Re: res.evalsRe[bestIdx], g2Im: res.evalsIm[bestIdx],
    };
}

// Per-mode solve config: boundary condition at the symmetry plane and the
// energy→capacitance factor kC.
//   full domain : single kC=2, diff kC=1 (even=[1,1], odd=[1,-1] drives)
//   half domain : single kC=4, diff kC=2 (even/odd via abc.left; signal=1)
function modeConfig(mode, isDiff, symmetry) {
    const kC = (isDiff ? 1 : 2) * (symmetry ? 2 : 1);
    let abc = {};
    if (symmetry && mode !== 'odd') abc = { left: 'pmc' };  // even/single: PMC at x=0
    return { abc, kC };
}

function drivePotentials(condRect, mode, symmetry) {
    return condRect.rectRoles.map(r => {
        if (!r.is_signal) return 0;
        if (symmetry) return 1;            // one trace in the half domain; mode via abc
        if (mode === 'even') return 1;
        if (mode === 'odd') return r.polarity || 1;
        return 1;
    });
}

export class TriBackend {
    constructor(ctx, solver, opts = {}) {
        this.ctx = ctx;
        this.solver = solver;
        this.opts = opts;
        this.mesh = null;
        this.fm = null;
        this.domain = null;
        this._static = null;   // per-mode cached static solve + resampled fields
    }

    get modeNames() { return this.solver.is_differential ? ['odd', 'even'] : ['single']; }

    // Build mesh (with symmetry + adaptive refinement), then solve & cache the
    // (frequency-independent) static problem per mode and resample fields.
    // onProgress({iteration, max_iterations, energy_error, param_error, nodes_x,
    // nodes_y}) is called once per adaptive refinement pass (same shape the
    // rectilinear backend reports, so the UI shows real passes for both).
    async buildMesh(onProgress = null) {
        const s = this.solver;
        const dom = { x_min: -s.domain_width / 2, x_max: s.domain_width / 2,
                      y_min: -s.t_gnd, y_max: s.domain_height };
        this.domain = dom;
        // Use a half-domain symmetry solve when the geometry is mirror-symmetric.
        // BUT the x=0 plane can only separate the differential even/odd modes when
        // it lies BETWEEN the signal conductors. For broadside pairs (traces stacked
        // at x=0, straddling the plane) both modes are x-symmetric, so the plane
        // cannot distinguish them — fall back to the full domain there.
        const tolX = s.domain_width * 1e-6;
        const diffStraddles = s.is_differential && s.conductors.some(
            c => c.is_signal && c.x_min < -tolX && c.x_max > tolX);
        this.symmetry = (this.opts.symmetry ?? true) && !diffStraddles &&
            isXSymmetric(s.conductors, s.dielectrics, s.domain_width);

        const tAbs = Math.max(Math.abs(s.t ?? 35e-6), 1e-9);
        const wRef = Math.max(s.w ?? (dom.x_max - dom.x_min) / 10, 1e-9);
        // Start coarse; adaptive refinement adds resolution where the field needs it.
        // GRADED (fine at conductors, coarse in the bulk), not uniform: the MQS volume
        // eddy-current loss needs the conductor region resolved, so a uniform ultra-coarse
        // start (as the SIBC reference viewer used) starves it and gives garbage R. With
        // the graded start, hFine can be ~3× coarser than the original here with no
        // accuracy loss — the conductor SURFACE is re-resolved by refineSkinBand for the
        // loss solve anyway, and eps/Z0 converge at a few hundred triangles. The
        // rough-stripline roughness loss is the binding constraint (it sits ~7% high vs
        // ref on any mesh); coarsening further than this eats its margin.
        const hFine = this.opts.hFine ?? Math.min(tAbs, wRef / 4) * 1.5;
        const hCoarse = this.opts.hCoarse ?? (dom.y_max - dom.y_min) / 5;
        let mesh = buildGmshMeshFromGeometry(this.ctx.G, {
            conductors: s.conductors, dielectrics: s.dielectrics,
            domain: dom, boundaries: s.boundaries, hFine, hCoarse, symmetry: this.symmetry,
            meshConductorInterior: true,   // conductor interiors meshed for the MQS loss solve
        });
        this.condRect = mesh.condRect;

        // Loss routines read condRects[0].symmetry to scale the half-domain integral.
        this.condRect.rects.forEach(r => { r.symmetry = this.condRect.symmetry; });

        // ---- Adaptive mesh refinement ----
        // The band+hole mesher produces high-quality triangles (Q≈1), so error-driven
        // refinement converges rather than degrading the Nedelec eigensolve. Metric:
        // the Zienkiewicz–Zhu estimator on the projected H-field (refines the
        // conductor surfaces / corners that dominate conductor loss), blended with the
        // static field energy (dielectric interfaces) — as in the reference solver.
        // Driven by the full-wave eigensolve at a representative frequency.
        const maxIters = this.opts.maxRefineIters ??
            (typeof process !== 'undefined' && process.env && process.env.TRI_MAXREF !== undefined
                ? +process.env.TRI_MAXREF : 6);
        const maxNodes = this.opts.maxNodes ?? 18000;
        const refTol = this.opts.refineTol ?? 0.003;
        const refineFrac = this.opts.refineFrac ?? 0.3;
        const minRefineNodes = this.opts.minRefineNodes ?? 0;
        // Refine on ALL solved modes (differential: odd AND even) so each mode's
        // critical regions are resolved — crucially the inter-trace gap, where the
        // odd mode's field concentrates and the mutual resistance R[12] comes from.
        const refModes = this.modeNames;
        const fRef = Math.max(s.freq || 1e9, 1e9);               // full-wave freq for ZZ
        let prev = null, prevEnergy = null;
        for (let it = 0; it <= maxIters; it++) {
            let metric = null;
            const conv = [];   // convergence quantities (eps + loss surface integral per mode)
            const energy = []; // total field energy per mode (for the reported energy error)
            for (const rm of refModes) {
                const { abc } = modeConfig(rm, s.is_differential, this.symmetry);
                const pot = drivePotentials(this.condRect, rm, this.symmetry);
                const fm = buildTriFreedomMap(mesh, this.condRect, abc);
                const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, pot);
                const phiAir = solveTriStatic(mesh, fm, null, pot);
                const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
                const eps_static = W_eps / computeTriEnergy(phiAir, mesh, null);
                energy.push(W_eps);
                const fw = fullwaveMode(this.ctx, mesh, fm, abc, this.condRect, mesh.epsMap, fRef, phiEps, eps_static);
                conv.push(fw ? fw.eps : eps_static);
                const metricS = perElementEnergy(phiEps, mesh, mesh.epsMap);
                let zz = null;
                if (fw) { try { zz = computeHtZZMetric(mesh, fm, fw.vRe, fw.vIm, fw.g2Re, fw.g2Im, fRef, this.condRect.rects, null, null, this.ctx.wasmSolver); } catch { zz = null; } }
                let mS = 0, mZ = 0;
                for (let i = 0; i < metricS.length; i++) { if (metricS[i] > mS) mS = metricS[i]; if (zz && zz[i] > mZ) mZ = zz[i]; }
                if (!metric) metric = new Float64Array(metricS.length);
                for (let i = 0; i < metric.length; i++) {
                    const v = (mS > 0 ? metricS[i] / mS : 0) + (zz && mZ > 0 ? zz[i] / mZ : 0);
                    if (v > metric[i]) metric[i] = v;   // union across modes
                }
            }
            // converged only when ALL quantities (eps + loss integrals) are stable
            let maxRel = Infinity, energyRel = Infinity;
            if (prev && prev.length === conv.length) {
                maxRel = 0;
                for (let i = 0; i < conv.length; i++) maxRel = Math.max(maxRel, Math.abs(conv[i] - prev[i]) / Math.max(Math.abs(prev[i]), 1e-300));
            }
            if (prevEnergy && prevEnergy.length === energy.length) {
                energyRel = 0;
                for (let i = 0; i < energy.length; i++) energyRel = Math.max(energyRel, Math.abs(energy[i] - prevEnergy[i]) / Math.max(Math.abs(prevEnergy[i]), 1e-300));
            }
            prev = conv; prevEnergy = energy;
            // Report this pass (real triangle count + convergence) so the UI shows
            // the adaptive progress just like the rectilinear backend. When a progress
            // sink is present, yield to the event loop so the browser can paint each
            // pass live instead of all at once when buildMesh returns.
            if (onProgress) {
                onProgress({
                    iteration: it + 1, max_iterations: maxIters + 1,
                    energy_error: isFinite(energyRel) ? energyRel : 1,
                    param_error: isFinite(maxRel) ? maxRel : 1,
                    n_tris: mesh.nTris, nodes_x: mesh.nNodes, nodes_y: 0,
                });
                await new Promise(r => setTimeout(r, 0));
            }
            // The mutual resistance R[12] of a coupled pair converges slower than
            // eps (it needs the inter-trace gap resolved), so for differential pairs
            // keep refining past eps-convergence until a node floor is reached.
            const enoughNodes = mesh.nNodes >= minRefineNodes;
            if (it === maxIters || (maxRel < refTol && enoughNodes) || mesh.nNodes > maxNodes) break;
            const marked = markTrianglesForRefinement(metric, refineFrac);
            const refined = refineTriMesh(mesh, marked);
            refined.condRect = this.condRect;
            const mats = tagMaterials(refined, s.dielectrics);
            refined.epsMap = mats.epsMap; refined.lossMap = mats.lossMap;
            mesh = refined;
        }

        this.mesh = mesh;
        this._prepareStatic();
        return mesh;
    }

    _prepareStatic() {
        const { mesh } = this;
        const s = this.solver;
        this._static = {};
        for (const mode of this.modeNames) {
            const { abc, kC } = modeConfig(mode, s.is_differential, this.symmetry);
            const fm = buildTriFreedomMap(mesh, this.condRect, abc);
            const pot = drivePotentials(this.condRect, mode, this.symmetry);
            const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, pot);
            const phiAir = solveTriStatic(mesh, fm, null, pot);
            const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
            const W_air = computeTriEnergy(phiAir, mesh, null);
            const W_loss = computeTriEnergy(phiEps, mesh, mesh.lossMap);
            const C0 = kC * eps0 * W_air;                 // vacuum cap (geometric)
            const eps_eff_static = W_eps / W_air;
            const parity = this.symmetry ? (mode === 'odd' ? 'odd' : 'even') : null;
            const fields = resampleStatic(mesh, phiEps, this.domain,
                { resolution: this.opts.resolution, parity });
            this._static[mode] = {
                fm, abc, kC, phiEps, C0, W_loss, eps_eff_static,
                Z_static: 1 / (c0 * Math.sqrt(C0 * eps_eff_static * C0)),
                fields,
            };
        }
    }

    // Per-mode result at frequency f, reusing the cached static solve.
    // eps_eff comes from the full-wave eigenmode (f ≥ 100 MHz) or the static
    // solve (below); conductor loss is from the robust static-field method.
    _modeAtFreq(mode, f) {
        const { mesh } = this;
        const s = this.solver;
        const cr = this.condRect;
        const st = this._static[mode];
        const { fm, kC, phiEps, C0, eps_eff_static, W_loss } = st;
        const omega = 2 * Math.PI * f;

        // Full-wave eigenmode above 100 MHz (dispersive eps + the mode field used
        // for the perturbation conductor loss); static solve near DC.
        let eps_d = eps_eff_static, fw = null;
        if (f >= F_STATIC_MAX) {
            fw = fullwaveMode(this.ctx, mesh, fm, st.abc, cr, mesh.epsMap, f, phiEps, eps_eff_static);
            if (fw && fw.eps > 0) eps_d = fw.eps;
        }
        // Dielectric dispersion is a capacitance effect: C = eps_d·C0 (geometric
        // L_external = 1/(c²C0)). Reduces to the static result when eps_d = eps_static.
        const C = eps_d * C0;
        const L_external = 1 / (c0 * c0 * C0);
        const Z0 = 1 / (c0 * C0 * Math.sqrt(eps_d));

        // dielectric loss: G from the lossy-permittivity energy integral
        const G = omega * kC * eps0 * W_loss;
        const alpha_d = Z0 > 0 ? (G * Z0 / 2) * NP_TO_DB : 0;

        // Conductor loss → R_total and the internal inductance L_internal (which adds
        // the conductor-roughness/skin ΔL to the phase ε_eff). ε_eff stays anchored to
        // the (reliable) dielectric eps_d: L_internal is a clamped increment, so a bad
        // loss solve can never corrupt ε_eff.
        const { sigma, rq } = effectiveSurface(s);
        const omu = omega * MU0;
        const delta = f > 0 ? Math.sqrt(2 / (omu * sigma)) : Infinity;
        const Rs_smooth = 1 / (sigma * Math.min(delta, 1e30));
        const Zr = calculateZrough(f > 0 ? f : 1, sigma, rq);
        const fR = Rs_smooth > 0 ? Zr.re / Rs_smooth : 1, fL = Rs_smooth > 0 ? Zr.im / Rs_smooth : 1;
        // MQS (reference) applies when grounds are wall-absorbed (signal-only rects)
        // and the domain is symmetric (mode set by the BC). Else H-field perturbation.
        const hasGroundRect = cr.rectRoles.some(r => !r.is_signal);
        const useMQS = this.symmetry && !hasGroundRect && cr.rects.length > 0;
        let R_total = 0, L_internal = 0;
        if (f > 0 && useMQS) {
            let minDim = Infinity;
            for (const c of cr.rects) minDim = Math.min(minDim, c.xmax - c.xmin, c.ymax - c.ymin);
            // Skin-band target element size (× δ): the band must be resolved to a
            // fraction of δ for accurate R (matches ms2d's mqs_band_delta).
            const bandDelta = this.opts.mqsBandDelta ?? 0.5;
            // 2δ-wide skin band (was 3δ): the eddy current is ~85% within 2δ, so the
            // narrower band halves the refined-mesh size with negligible R change.
            const mqsBand = this.opts.mqsBand ?? 2;
            const mqsMesh = (delta < minDim) ? refineSkinBand(mesh, cr, delta, 12, mqsBand, bandDelta * delta) : mesh;
            let mqs = null;
            try {
                mqs = mqsConductorLoss(mqsMesh, cr, f, sigma, this.ctx.helpers.solveSparseMulti, 0,
                    { topGround: !!(cr.wallPEC && cr.wallPEC.top), Rq: rq, oddSymmetry: mode === 'odd' });
            } catch { mqs = null; }
            // Accept MQS only if physically sane.
            if (mqs && isFinite(mqs.R_total) && mqs.R_total > 0 && isFinite(mqs.L_loop) && mqs.L_loop > 0) {
                R_total = s.is_differential ? mqs.R_total / 2 : mqs.R_total;
                L_internal = Math.max(0, Math.min(mqs.L_loop - L_external, 3 * L_external));
            }
        }
        if (f > 0 && R_total === 0 && fw) {   // perturbation (full domain, or MQS fallback)
            const k2 = (omega / c0) ** 2;
            const an = analyzeTriMode(fw.vRe, fw.vIm, fw.g2Re, fw.g2Im, mesh, fm, k2, f, mesh.epsMap, cr);
            if (an && Math.abs(an.P) > 1e-30) {
                const loss = solveConductorLoss(cr.rects, f, sigma, mesh, fm, fw.vRe, fw.vIm,
                    fw.g2Re, fw.g2Im, Math.abs(an.P), Z0, mesh.epsMap, buildLossEdges(mesh, fm, cr), null, null, this.ctx.wasmSolver);
                const R_ac = loss.R_ac * fR;
                R_total = Math.sqrt(loss.R_dc * loss.R_dc + R_ac * R_ac);
                L_internal = Math.max(0, Math.min(loss.L_int * fL, 3 * L_external));
            }
        }
        if (f > 0 && R_total === 0) {   // static-field fallback (near DC / eigensolve failed)
            const loss = staticConductorLoss(cr.rects, f, sigma, mesh, fm, phiEps,
                Z0, eps_eff_static, eps_d, null, buildLossEdges(mesh, fm, cr));
            const R_ac = loss.R_ac * fR;
            R_total = Math.sqrt(loss.R_dc * loss.R_dc + R_ac * R_ac);
        }

        // assemble RLGC + Zc
        // Reported eps_eff = phase ε_eff = c²·L·C: dielectric dispersion (via eps_d)
        // plus the conductor-roughness ΔL increment (matches the FDM convention).
        const L_total = L_external + L_internal;
        const eps_eff = c0 * c0 * L_total * C;
        let Zc, alpha_c;
        if (f > 0) {
            const Znum = new Complex(R_total, omega * L_total);
            const Zden = new Complex(G, omega * C);
            Zc = Znum.div(Zden).sqrt();
            alpha_c = Zc.re > 0 ? NP_TO_DB * R_total / (2 * Zc.re) : 0;
        } else {
            Zc = new Complex(Z0, 0);
            alpha_c = 0;
        }

        return {
            mode, Z0, eps_eff, eps_eff_mode: eps_d, C, C0,
            RLGC: { R: R_total, L: L_total, G, C },
            Zc,
            alpha_c, alpha_d, alpha_total: alpha_c + alpha_d,
            L_internal, L_external,
            V: st.fields.V, Ex: st.fields.Ex, Ey: st.fields.Ey,
        };
    }

    // Public: solve at one frequency, return the unified result object and write
    // grid fields + triangle mesh back onto the solver for plotting/export.
    solveAt(f) {
        if (!this.mesh) this.buildMesh();
        const modes = this.modeNames.map(m => this._modeAtFreq(m, f));
        const result = { modes };
        if (this.solver.is_differential) {
            const odd = modes.find(m => m.mode === 'odd');
            const even = modes.find(m => m.mode === 'even');
            result.Z_diff = 2 * odd.Z0;
            result.Z_common = even.Z0 / 2;
            result.RLGC_matrix = modalToPhysicalRLGC(odd, even);
        }
        // Write plotting state onto the solver (mirrors the FDM backend).
        const f0 = this._static[this.modeNames[0]].fields;
        this.solver.x = f0.x;
        this.solver.y = f0.y;
        this.solver.V = modes.map(m => m.V);
        this.solver.Ex = modes.map(m => m.Ex);
        this.solver.Ey = modes.map(m => m.Ey);
        this.solver.triMesh = { nodes: this.mesh.nodes, tris: this.mesh.tris, nTris: this.mesh.nTris };
        this.solver.solution_valid = true;
        return result;
    }
}

// even/odd modal -> physical 2x2 (matches field_solver._modal_to_physical_rlgc)
function modalToPhysicalRLGC(odd, even) {
    const mk = (xo, xe) => [[(xo + xe) / 2, (xe - xo) / 2], [(xe - xo) / 2, (xo + xe) / 2]];
    return {
        R: mk(odd.RLGC.R, even.RLGC.R),
        L: mk(odd.RLGC.L, even.RLGC.L),
        G: mk(odd.RLGC.G, even.RLGC.G),
        C: mk(odd.RLGC.C, even.RLGC.C),
    };
}
