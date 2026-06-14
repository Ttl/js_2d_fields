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
//     (gradient model, calculate_Zrough): plating = "roughness with a different
//     conductivity" on the plated conductor surfaces.
//
// Mesh + freedom maps are built once and reused across a frequency sweep; only
// the eigen-assembly (k²-dependent) and loss are recomputed per frequency.

import createModule from '../wasm_solver/eigen_solver.js';
import { createWasmHelpers } from './fem_core.js';
import { initGmsh } from './gmsh_mesh.js';
import { buildOccMeshFromGeometry, tagMaterials } from './occ_to_mesh.js';
import { buildTriFreedomMap, solveTriStatic, computeTriEnergy, refineTriMesh,
         markTrianglesForRefinement, computeTriP2StaticMatrices,
         staticToEdgeDofs, assembleTriFEM } from './tri_fem.js';
import { staticConductorLoss, solveConductorLoss, computeHtZZMetric } from './conductor_loss.js';
import { mqsConductorLoss, refineSkinBand } from './mqs_loss.js';
import { analyzeTriMode } from './tri_ms_solver.js';
import { checkMeshQuality } from './tri_mesh.js';

// Below this frequency the full-wave eigensolve degenerates near DC; use the
// static (variational) solve instead. Above it, use the full-wave eigenmode for
// the dispersive effective permittivity.
const F_STATIC_MAX = 100e6;
import { calculate_Zrough, calculate_Zrough_layered } from '../surface_roughness.js';
import { resampleStatic, resampleModeField } from './resample.js';
import { Complex } from '../complex.js';
import { djordjevic_sarkar } from '../djordjevic_sarkar.js';

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

// Complex permittivity map for the mode-viewer EIGENSOLVE. The reference viewer folds the
// dielectric loss into ε (im = −ε_r·tanδ), and that small imaginary part DAMPS the discrete
// vector-FEM spurious modes — pushing them off the real-γ² propagating axis — so they no
// longer masquerade as physical propagating modes. tagMaterials stores ε real (im=0) with the
// loss kept separately (for the energy-based G), so reconstruct the complex map here. This is
// THE reason the OCC backend saw spurious modes the reference's (lossy) eigensolve did not —
// it was never the mesh.
function lossyEpsMap(mesh, tandFloor = 0.002) {
    const { epsMap, lossMap } = mesh;
    const out = new Array(epsMap.length);
    for (let t = 0; t < epsMap.length; t++) {
        const er = epsMap[t].re;
        // Floor a tiny loss so the spurious-mode damping still works for (near-)lossless
        // dielectrics; negligible vs a real tanδ and it barely moves ε_eff=−Re(γ²)/k².
        const loss = Math.max(lossMap ? lossMap[t].re : 0, er * tandFloor);
        out[t] = { re: er, im: -loss };
    }
    return out;
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

// Per-face plating model. Each conductor face (top / bottom / sides) gets its own
// surface impedance: plated faces use the layered (plating-over-bulk) gradient
// impedance — matching the FDM backend's per-face model; bare faces and grounded
// walls use the base metal. NOTE: thick-corner plating (geometric wrap of
// side/top plating around corners) is NOT modeled — disabled in the UI for the
// full-wave solver.
//
// Two consumers, sharing makePlatingZs:
//   • buildSurfaceGroups (PERTURBATION path) groups loss edges by surface
//     impedance so the surface integral is summed per group (loss is linear in
//     the per-edge surface resistance). `uniform` true ⇒ a single impedance
//     everywhere ⇒ the plain single-pass path.
//   • buildFaceZs (MQS path) returns Zs at a face midpoint, so the MQS volume
//     solve can weight each face's smooth current by its own plating impedance.
// Shared surface-impedance model: the bare-metal Zs plus a (rectIndex, face) → Zs
// resolver. Plated faces use the layered (plating-over-bulk) gradient impedance
// (cached per distinct plating); bare faces use the base metal. Used by both the
// edge-based grouping (perturbation path) and the point-based lookup (MQS path).
function makePlatingZs(solver, condRect, freq) {
    const roles = condRect.rectRoles || [];
    const sigmaBase = solver.sigma_cond ?? 5.8e7;
    const rqBase = solver.rq ?? 0;
    const Zbare = calculate_Zrough(freq, sigmaBase, rqBase);    // Complex (.re/.im)
    const layeredCache = new Map();
    const zForFace = (ri, face) => {
        const pl = roles[ri] && roles[ri].plating;
        if (!(pl && pl[face] && pl.sigma > 0)) return Zbare;
        const key = `${pl.sigma}|${pl.rq}|${pl.thickness}`;
        let z = layeredCache.get(key);
        if (!z) {
            z = calculate_Zrough_layered(freq, sigmaBase, pl.rq ?? 0, pl.sigma, pl.thickness ?? 0);
            layeredCache.set(key, z);
        }
        return z;
    };
    return { Zbare, zForFace };
}

function platingTol(condRect) {
    const diag = Math.hypot((condRect.xmax_domain ?? condRect.xmax) - (condRect.xmin_domain ?? condRect.xmin),
                            (condRect.ymax_domain ?? condRect.ymax) - (condRect.ymin_domain ?? condRect.ymin));
    return (diag > 0 ? diag : 1) * 1e-6;
}

// Point-based per-face surface impedance: Zs at a face midpoint (x, y) with edge
// orientation 'h' (top/bottom face) or 'v' (side face). Mesh-independent — used by
// the MQS loss to weight each face's smooth current by its plating impedance.
// Returns bare metal off any conductor face (e.g. grounded walls).
function buildFaceZs(solver, condRect, freq) {
    const rects = condRect.rects || [];
    const tol = platingTol(condRect);
    const { Zbare, zForFace } = makePlatingZs(solver, condRect, freq);
    return (x, y, orient) => {
        for (let ri = 0; ri < rects.length; ri++) {
            const r = rects[ri];
            if (x < r.xmin - tol || x > r.xmax + tol || y < r.ymin - tol || y > r.ymax + tol) continue;
            let face = null;
            if (orient === 'h') {
                if (Math.abs(y - r.ymax) < tol) face = 'top';
                else if (Math.abs(y - r.ymin) < tol) face = 'bottom';
            } else if (Math.abs(x - r.xmin) < tol || Math.abs(x - r.xmax) < tol) {
                face = 'sides';
            }
            if (!face) continue;
            return zForFace(ri, face);
        }
        return Zbare;
    };
}

function buildSurfaceGroups(solver, mesh, fm, condRect, baseMask, freq) {
    const { nodes, edges, nEdges } = mesh;
    const rects = condRect.rects || [];
    const tol = platingTol(condRect);
    const { Zbare, zForFace } = makePlatingZs(solver, condRect, freq);

    // Surface impedance for one loss edge. Conductor edges resolve to the face of
    // the conductor rect they lie on; everything else (grounded walls) is bare.
    const zOfEdge = (e) => {
        if (!(fm.isCondEdge && fm.isCondEdge[e])) return Zbare;
        const n0 = edges[2 * e], n1 = edges[2 * e + 1];
        const x0 = nodes[2 * n0], y0 = nodes[2 * n0 + 1];
        const x1 = nodes[2 * n1], y1 = nodes[2 * n1 + 1];
        for (let ri = 0; ri < rects.length; ri++) {
            const r = rects[ri];
            if (x0 < r.xmin - tol || x0 > r.xmax + tol || x1 < r.xmin - tol || x1 > r.xmax + tol) continue;
            if (y0 < r.ymin - tol || y0 > r.ymax + tol || y1 < r.ymin - tol || y1 > r.ymax + tol) continue;
            let face = null;
            if (Math.abs(y0 - r.ymax) < tol && Math.abs(y1 - r.ymax) < tol) face = 'top';
            else if (Math.abs(y0 - r.ymin) < tol && Math.abs(y1 - r.ymin) < tol) face = 'bottom';
            else if ((Math.abs(x0 - r.xmin) < tol && Math.abs(x1 - r.xmin) < tol) ||
                     (Math.abs(x0 - r.xmax) < tol && Math.abs(x1 - r.xmax) < tol)) face = 'sides';
            if (!face) continue;
            return zForFace(ri, face);
        }
        return Zbare;
    };

    const keyOf = (z) => `${z.re.toExponential(6)}|${z.im.toExponential(6)}`;
    const groupIdx = new Map();
    const groupZ = [];
    const edgeGroup = new Int32Array(nEdges).fill(-1);
    for (let e = 0; e < nEdges; e++) {
        if (!baseMask[e]) continue;
        const z = zOfEdge(e);
        const k = keyOf(z);
        let gi = groupIdx.get(k);
        if (gi === undefined) { gi = groupZ.length; groupIdx.set(k, gi); groupZ.push(z); }
        edgeGroup[e] = gi;
    }
    const groups = groupZ.map((z, gi) => {
        const mask = new Uint8Array(nEdges);
        for (let e = 0; e < nEdges; e++) if (edgeGroup[e] === gi) mask[e] = 1;
        return { Zs: z, mask };
    });
    return { uniform: groups.length <= 1, groups };
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
// selection: among physical modes (γ²<0, eps not near cutoff), pick the one with
// the HIGHEST overlap to the static drive (the static field is the quasi-TEM shape) —
// reliable even when several modes cluster at eps ≈ ε_r (homogeneous/stripline fills).
// Returns { eps, vRe, vIm, g2Re, g2Im } for the chosen mode, or null.
// k²-decomposition cache. The FEM matrices A(k²), B(k²) are exactly affine in k²
// (Lee–Jin form; the only √k² term, the radiating Robin ABC, is absent here — the walls
// are PEC/PMC), with fixed sparsity across frequency. Assemble at the first two distinct
// frequencies, store A0/A1 and B0/B1, then rebuild any frequency's CSR with a cheap
// a0+k²·a1 combine — skipping the per-frequency element re-assembly, which profiling
// showed costs MORE than the eigensolve itself. Exact to machine precision; falls back
// to per-frequency assembly if the two sparsity patterns ever differ.
function buildFemK2Cache(femA, kA, femB, kB) {
    const samePat = (x, y) => {
        if (x.colIdx.length !== y.colIdx.length || x.rowPtr.length !== y.rowPtr.length) return false;
        for (let i = 0; i < x.colIdx.length; i++) if (x.colIdx[i] !== y.colIdx[i]) return false;
        for (let i = 0; i < x.rowPtr.length; i++) if (x.rowPtr[i] !== y.rowPtr[i]) return false;
        return true;
    };
    if (!samePat(femA.csrA, femB.csrA) || !samePat(femA.csrB, femB.csrB)) return null;
    const inv = 1 / (kB - kA);
    const mk = (ca, cb) => {
        const n = ca.valRe.length, a0 = new Float64Array(n), a1 = new Float64Array(n);
        for (let i = 0; i < n; i++) { a1[i] = (cb.valRe[i] - ca.valRe[i]) * inv; a0[i] = ca.valRe[i] - kA * a1[i]; }
        return { rowPtr: ca.rowPtr, colIdx: ca.colIdx, a0, a1, zeros: new Float64Array(n) };
    };
    return { N: femA.N, A: mk(femA.csrA, femB.csrA), B: mk(femA.csrB, femB.csrB) };
}
function femFromK2Cache(cache, k2) {
    const build = (m) => {
        const n = m.a0.length, v = new Float64Array(n);
        for (let i = 0; i < n; i++) v[i] = m.a0[i] + k2 * m.a1[i];
        return { rowPtr: m.rowPtr, colIdx: m.colIdx, valRe: v, valIm: m.zeros };
    };
    return { csrA: build(cache.A), csrB: build(cache.B), N: cache.N };
}

function fullwaveMode(ctx, mesh, fm, abc, condRect, epsMap, f, phiEps, eps_static, seedVec = null, cacheBox = null) {
    const k2 = (2 * Math.PI * f / 299792458) ** 2;
    let fem;
    if (cacheBox && cacheBox.ready) {
        fem = femFromK2Cache(cacheBox.ready, k2);   // O(nnz) rebuild, no re-assembly
    } else {
        try { fem = assembleTriFEM(mesh, fm, k2, epsMap, abc, condRect, null); }
        catch { return null; }
        if (cacheBox) {   // lazily build the cache from the first two distinct frequencies
            if (cacheBox.prev && Math.abs(cacheBox.prev.k2 - k2) > 0) {
                cacheBox.ready = buildFemK2Cache(cacheBox.prev.fem, cacheBox.prev.k2, fem, k2);
                cacheBox.prev = null;
            } else if (!cacheBox.prev) {
                cacheBox.prev = { k2, fem };
            }
        }
    }
    if (!fem) return null;
    const N = fem.N;
    // Warm start the Arnoldi from the previous frequency's eigenvector when available
    // (same mesh → same N); falls back to the static-field projection. A near-converged
    // start cuts the iterations across a sweep. The quasi-TEM mode-pick below always
    // overlaps against the static drive (staticSeed), independent of the warm start.
    const staticSeed = staticToEdgeDofs(phiEps, mesh, fm);
    const seed = (seedVec && seedVec.length === N) ? seedVec : staticSeed;
    // 4 eigenpairs is ample: shift-invert is centered on the quasi-TEM eigenvalue
    // (-k2·eps_static), so the target mode converges first. (Was 8 — halving the
    // Krylov subspace markedly speeds the eigensolve with no change in the picked mode.)
    const nev = 4, ncv = Math.min(2 * nev + 1, N - 1);
    let res;
    try { res = ctx.helpers.solveGeneralized(N, fem.csrA, fem.csrB, [-k2 * eps_static, 0], nev, ncv, seed); }
    catch { return null; }
    if (!res || !res.nconv) return null;
    // Pick the quasi-TEM mode by MAXIMUM overlap with the static drive (the static
    // field IS the quasi-TEM shape). This is robust where several eigenmodes cluster
    // at eps ≈ ε_r — e.g. a homogeneous (stripline) fill in an open/parallel-plate
    // domain, where the odd quasi-TEM is near-degenerate with lateral resonances:
    // a "closest-eps to eps_static" pick then jumps between them across frequency,
    // producing the loss/eps_eff dips. A loose eps gate drops near-cutoff / continuum
    // modes (eps ≪ eps_static) that can otherwise score a deceptively high overlap.
    let bestIdx = -1, bestOvl = 0.3;                       // require overlap > 0.3
    for (let i = 0; i < res.nconv; i++) {
        if (res.evalsRe[i] >= 0) continue;                 // physical: γ² < 0
        const epsMode = -res.evalsRe[i] / k2;
        if (epsMode < 0.5 * eps_static) continue;          // drop near-cutoff / lateral-continuum modes
        const vR = res.evecsRe.subarray(i * N, (i + 1) * N);
        let dot = 0, nS = 0, nV = 0;
        for (let k = 0; k < fm.nFreeTransverse; k++) { dot += staticSeed[k] * vR[k]; nS += staticSeed[k] ** 2; nV += vR[k] ** 2; }
        const ovl = nS > 0 && nV > 0 ? Math.abs(dot) / Math.sqrt(nS * nV) : 0;
        globalThis.__TRI_DEBUG__ && console.log(`    cand ${i}: eps=${epsMode.toFixed(4)} ovl=${ovl.toFixed(3)} (target ${eps_static.toFixed(4)})`);
        if (ovl > bestOvl) { bestOvl = ovl; bestIdx = i; }
    }
    if (bestIdx < 0) return null;
    return {
        eps: -res.evalsRe[bestIdx] / k2,
        vRe: res.evecsRe.slice(bestIdx * N, (bestIdx + 1) * N),
        vIm: res.evecsIm.slice(bestIdx * N, (bestIdx + 1) * N),
        g2Re: res.evalsRe[bestIdx], g2Im: res.evalsIm[bestIdx],
    };
}

// Per-mode solve config: wall boundary conditions + energy→capacitance factor kC.
//   full domain : single kC=2, diff kC=1 (even=[1,1], odd=[1,-1] drives)
//   half domain : single kC=4, diff kC=2 (even/odd via abc.left; signal=1)
//
// Wall BCs are driven by the wall types so the triangular backend matches the
// rectilinear quasi-static solver: 'gnd' walls are PEC (Dirichlet V=0), 'open' walls
// are natural/Neumann (the FDM sets zero off-domain flux there). buildTriFreedomMap
// makes a wall PEC unless abc.{left,right,top,bottom} is truthy; 'pmc' selects the
// natural BC WITHOUT the full-wave radiation Robin term (which only fires on === true).
// Far 'open' walls give the same answer either way (the field has decayed), but a
// close wall (small enclosure) needs the matching natural BC.  wallPEC.X === true ⇒ gnd.
function modeConfig(mode, isDiff, symmetry, wallPEC) {
    const kC = (isDiff ? 1 : 2) * (symmetry ? 2 : 1);
    const abc = {};
    if (wallPEC) {
        if (!wallPEC.right) abc.right = 'pmc';
        if (!wallPEC.top) abc.top = 'pmc';
        if (!wallPEC.bottom) abc.bottom = 'pmc';
        if (!wallPEC.left) abc.left = 'pmc';
    }
    // Symmetry plane at x=0 overrides the left wall: PMC (even/single) or PEC (odd).
    if (symmetry) { if (mode === 'odd') delete abc.left; else abc.left = 'pmc'; }
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

    // Cap a bulk element size to a fraction of the shortest in-medium wavelength at the
    // mode-viewer frequency (opts.modesFreq), so high-frequency cavity / higher-order modes
    // are actually resolved. Without this, a genuine high-f mode is under-resolved on the
    // geometry-scale mesh, drifts between the fine and coarse solves, and the mesh-convergence
    // test mis-flags it as spurious. No effect in the quasi-TEM regime (λ/12 ≫ geometry hCoarse).
    _wavelengthCap(h) {
        const f = this.opts.modesFreq;
        if (!f) return h;
        const epsMax = Math.max(1, ...this.solver.dielectrics.map(d => d.epsilon_r || 1));
        const lamMin = c0 / (f * Math.sqrt(epsMax));
        return Math.min(h, lamMin / 12);
    }

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
        const hCoarse = this._wavelengthCap(this.opts.hCoarse ?? (dom.y_max - dom.y_min) / 5);
        let mesh = buildOccMeshFromGeometry(this.ctx.G, {
            conductors: s.conductors, dielectrics: s.dielectrics,
            domain: dom, boundaries: s.boundaries, hFine, hCoarse, symmetry: this.symmetry,
            meshConductorInterior: true,   // conductor interiors meshed for the MQS loss solve
            gmshOptions: this.opts.gmshOptions,
        });
        this.condRect = mesh.condRect;

        // Loss routines read condRects[0].symmetry to scale the half-domain integral.
        this.condRect.rects.forEach(r => { r.symmetry = this.condRect.symmetry; });

        // ---- Adaptive mesh refinement ----
        // The OCC field-graded mesh produces high-quality triangles (Q≈2), so error-driven
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
        const refineFrac = this.opts.refineFrac ?? 0.15;   // matches the FDM backend's fraction
        const minRefineNodes = this.opts.minRefineNodes ?? 0;
        // Require the convergence test to hold for this many CONSECUTIVE passes before
        // stopping (matches the FDM backend's min_converged_passes / UI "Min Converged
        // Passes"): guards against a premature stop on a single lucky pass.
        const minConvergedPasses = Math.max(1, this.opts.minConvergedPasses ?? 1);
        let convergedCount = 0;
        // Refine on ALL solved modes (differential: odd AND even) so each mode's
        // critical regions are resolved — crucially the inter-trace gap, where the
        // odd mode's field concentrates and the mutual resistance R[12] comes from.
        const refModes = this.modeNames;
        // Drive mesh refinement with the cheaper static field-energy metric (default):
        // it refines the same critical regions as the full-wave H-field ZZ estimator
        // (verified identical Z0/eps/loss across the test suite) while skipping the
        // per-iteration eigensolve, roughly halving buildMesh time. Set true to restore
        // the full-wave ZZ-driven refinement.
        const refineFullwave = this.opts.refineFullwave ?? false;
        const fRef = Math.max(s.freq || 1e9, 1e9);               // full-wave freq for ZZ
        let prev = null, prevEnergy = null;
        for (let it = 0; it <= maxIters; it++) {
            let metric = null;
            const conv = [];   // convergence quantities (eps + loss surface integral per mode)
            const energy = []; // total field energy per mode (for the reported energy error)
            for (const rm of refModes) {
                const { abc } = modeConfig(rm, s.is_differential, this.symmetry, this.condRect.wallPEC);
                const pot = drivePotentials(this.condRect, rm, this.symmetry);
                const fm = buildTriFreedomMap(mesh, this.condRect, abc);
                const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, pot);
                const phiAir = solveTriStatic(mesh, fm, null, pot);
                const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
                const W_air = computeTriEnergy(phiAir, mesh, null);
                const eps_static = W_eps / W_air;
                energy.push(W_eps);
                const fw = refineFullwave
                    ? fullwaveMode(this.ctx, mesh, fm, abc, this.condRect, mesh.epsMap, fRef, phiEps, eps_static)
                    : null;
                // Drive convergence on the static characteristic impedance Z = sqrt(L/C)
                // ∝ 1/sqrt(W_eps·W_air) (constant factors cancel in the relative change),
                // matching the reference adaptive driver. Z keeps refining the conductor /
                // inter-trace-gap regions (which set L) past the point where eps — a pure
                // C-ratio — has already stabilised, so the adaptive does meaningful work and
                // the loss / mutual-R cases gain margin.
                conv.push(1 / Math.sqrt(Math.max(W_eps * W_air, 1e-300)));
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
            // Count consecutive converged passes; only stop once minConvergedPasses in a
            // row meet the tolerance (reset the streak on any non-converged pass).
            if (maxRel < refTol && enoughNodes) convergedCount++; else convergedCount = 0;
            if (it === maxIters || convergedCount >= minConvergedPasses || mesh.nNodes > maxNodes) break;
            const marked = markTrianglesForRefinement(metric, refineFrac);
            const refined = refineTriMesh(mesh, marked);
            refined.condRect = this.condRect;
            const mats = tagMaterials(refined, s.dielectrics);
            refined.epsMap = mats.epsMap; refined.lossMap = mats.lossMap;
            mesh = refined;
        }

        this.mesh = mesh;
        // Mesh quality of the final (post-refinement) mesh. Q = circumradius/(2·inradius),
        // ideal 1; a very high Qmax means a sliver that ill-conditions the FEM solve.
        // Surfaced to the UI as a warning (empty constraint args → quality metrics only).
        try {
            const mq = checkMeshQuality(mesh, [], []);
            this.meshQuality = mq.metrics;
            if (this.solver) this.solver.meshQuality = mq.metrics;
        } catch { this.meshQuality = null; }
        this._prepareStatic();
        return mesh;
    }

    _prepareStatic() {
        const { mesh } = this;
        const s = this.solver;
        this._static = {};
        // Resample fields onto the FDM mesher's graded grid (fine near conductors and
        // the ground surface) so contour plots match the rectilinear backend. A uniform
        // grid aliases sub-grid-thickness traces/grounds into spurious contour bands and
        // a step line above the ground. Fall back to a uniform grid if no mesher exists.
        let grid = null;
        try {
            if (s.mesher && typeof s.mesher.generate_mesh === 'function') {
                const [gx, gy] = s.mesher.generate_mesh();
                if (gx && gy && gx.length && gy.length) grid = { x: gx, y: gy };
            }
        } catch { grid = null; }
        for (const mode of this.modeNames) {
            const { abc, kC } = modeConfig(mode, s.is_differential, this.symmetry, this.condRect.wallPEC);
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
                { resolution: this.opts.resolution, parity, grid });
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
        // Conductor-loss method (UI "Solver" dropdown for the full-wave backend):
        //   'auto'         — MQS volume eddy-current where applicable (symmetric, no
        //                    ground rect), else the blended perturbation. Most
        //                    accurate; default ("Full-wave (MQS)").
        //   'perturbation' — blended perturbation (static-field in the skin transition,
        //                    SIBC eigenmode at deep skin); much cheaper than MQS per
        //                    point, smooth across frequency. ("Full-wave (perturbation)").
        //   'static'       — static-field perturbation only (no SIBC blend).
        const lossMethod = this.opts.lossMethod ?? 'auto';
        let eps_d = eps_eff_static, fw = null;
        if (f >= F_STATIC_MAX) {
            // Causal materials change epsMap per frequency, so the k2-interpolation cache
            // (which assumes a constant epsMap) must be bypassed — re-assemble each freq.
            let cacheBox = null;
            if (!s.use_causal_materials) cacheBox = st.femCacheBox ??= { ready: null, prev: null };
            // Seed the eigensolve with the STATIC field (frequency-deterministic), NOT a
            // cross-frequency eigenvector warm-start. The static seed is the quasi-TEM
            // shape (a good initial guess) AND it makes the picked eigenvector independent
            // of evaluation ORDER. With a cross-frequency warm-start, a near-degenerate
            // cluster (e.g. the open homogeneous stripline odd mode + lateral resonances)
            // converges to a different in-subspace mixture depending on which frequency
            // ran last — so the out-of-order interpolating sweep produced order-dependent
            // conductor loss → ripple + non-convergent interpolant (excess solves).
            fw = fullwaveMode(this.ctx, mesh, fm, st.abc, cr, mesh.epsMap, f, phiEps, eps_eff_static,
                null, cacheBox);
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
        const Zr = calculate_Zrough(f > 0 ? f : 1, sigma, rq);
        const fR = Rs_smooth > 0 ? Zr.re / Rs_smooth : 1, fL = Rs_smooth > 0 ? Zr.im / Rs_smooth : 1;

        // Per-face plating: split the loss surface into surface-impedance groups.
        // Used by the perturbation path below (when MQS doesn't apply) to evaluate
        // and sum the loss per group; >1 group (e.g. top/sides plated, bottom bare)
        // means a face-dependent impedance. (MQS handles per-face plating itself via
        // surfaceZs, so plating no longer forces the perturbation path.)
        const lossEdgeMask = f > 0 ? buildLossEdges(mesh, fm, cr) : null;
        const surf = lossEdgeMask ? buildSurfaceGroups(s, mesh, fm, cr, lossEdgeMask, f) : { uniform: true, groups: [] };
        const platingPerSide = !surf.uniform;

        // MQS (reference) applies when grounds are wall-absorbed (signal-only rects)
        // and the domain is symmetric (mode set by the BC). Else H-field perturbation.
        // Per-face plating is handled INSIDE MQS (surfaceZs weights each face's smooth
        // current by its own impedance), so plating no longer forces perturbation;
        // when MQS doesn't apply, the perturbation path handles plating per-group.
        const hasGroundRect = cr.rectRoles.some(r => !r.is_signal);
        const anyPlating = cr.rectRoles.some(r => r.plating && r.plating.sigma > 0
            && (r.plating.top || r.plating.sides || r.plating.bottom));
        const useMQS = (lossMethod === 'mqs')
            || (lossMethod === 'auto' && this.symmetry && !hasGroundRect && cr.rects.length > 0);
        let R_total = 0, L_internal = 0;
        if (f > 0 && useMQS) {
            // The volume eddy solve runs the conductor BODY at the bulk metal σ;
            // plating is a SURFACE effect applied only through surfaceZs (relative to
            // the bulk Rs). So with plating, use the bulk σ here — NOT effectiveSurface's
            // plating σ (which would bake plating into the body and make surfaceZs a
            // no-op). The skin band is sized by the matching bulk skin depth.
            const mqsSigma = anyPlating ? (s.sigma_cond ?? 5.8e7) : sigma;
            const mqsDelta = anyPlating ? Math.sqrt(2 / (omu * mqsSigma)) : delta;
            let minDim = Infinity;
            for (const c of cr.rects) minDim = Math.min(minDim, c.xmax - c.xmin, c.ymax - c.ymin);
            // Skin-band target element size (× δ): the band must be resolved to a
            // fraction of δ for accurate R (matches ms2d's mqs_band_delta).
            // Skin-band target element size (×δ) and band width (×δ): resolve the skin
            // layer to bandDelta·δ within mqsBand·δ of each conductor surface. Tuned down
            // from (2, 0.5): a 1.5δ band with 0.6δ elements roughly halves the refined
            // mesh on tight/high-frequency geometries with ≤1% change in R (validated
            // against the reference suite), keeping the reused mesh small and the solve
            // fast.
            const bandDelta = this.opts.mqsBandDelta ?? 0.6;
            const mqsBand = this.opts.mqsBand ?? 1.5;
            // DEFAULT: f_max-reuse. Build the skin mesh at the HIGHEST frequency seen
            // (finest target, narrowest band) and reuse it for all lower frequencies.
            // The fine near-surface band plus the always-meshed conductor interior
            // resolve the lower-frequency (more uniform) current accurately — validated
            // against the reference suite to a 100× frequency mismatch. One mesh ⇒ smooth
            // R(f) with no per-frequency-remesh "dip", at moderate size (the high-f band
            // is narrow), unlike a whole-range mesh which over-resolves the low-f band.
            // opts.mqsCacheMesh === false falls back to a per-frequency remesh (marginally
            // faster on wide sweeps, but can show a small non-monotonic wiggle).
            let mqsMesh;
            if (mqsDelta >= minDim) {
                mqsMesh = mesh;
            } else if (this.opts.mqsCacheMesh === false) {
                mqsMesh = refineSkinBand(mesh, cr, mqsDelta, 12, mqsBand, bandDelta * mqsDelta);
            } else {
                const sc = st.skinCache || (st.skinCache = { dB: Infinity, mesh: null });
                if (!sc.mesh || mqsDelta < sc.dB * (1 - 1e-9)) {   // higher freq appeared → rebuild finer
                    sc.dB = mqsDelta;
                    sc.mesh = refineSkinBand(mesh, cr, mqsDelta, 12, mqsBand, bandDelta * mqsDelta);
                }
                mqsMesh = sc.mesh;
            }
            // Per-face plating ⇒ weight each face's smooth current by its own surface
            // impedance (surfaceZs); otherwise the uniform roughness factor (Rq).
            const mqsOpts = { topGround: !!(cr.wallPEC && cr.wallPEC.top), oddSymmetry: mode === 'odd' };
            if (anyPlating) mqsOpts.surfaceZs = buildFaceZs(s, cr, f);
            else mqsOpts.Rq = rq;
            let mqs = null;
            try {
                mqs = mqsConductorLoss(mqsMesh, cr, f, mqsSigma, this.ctx.helpers.solveSparseMulti, 0, mqsOpts);
            } catch { mqs = null; }
            // Accept MQS only if physically sane.
            if (mqs && isFinite(mqs.R_total) && mqs.R_total > 0 && isFinite(mqs.L_loop) && mqs.L_loop > 0) {
                R_total = s.is_differential ? mqs.R_total / 2 : mqs.R_total;
                L_internal = Math.max(0, Math.min(mqs.L_loop - L_external, 3 * L_external));
            }
        }
        // Perturbation conductor loss (used for the 'perturbation' option and as the
        // 'auto' fallback when MQS doesn't apply). Two estimators are blended by the
        // skin-depth / conductor-thickness ratio:
        //   • static-field perturbation — smooth and robust across the skin transition
        //     (δ ≳ thickness), where it matches MQS; used at low frequency.
        //   • SIBC eigenmode perturbation — more accurate at deep skin (δ ≪ thickness,
        //     e.g. rough conductors at high f), but its corner-singularity model is
        //     non-monotonic for δ ≳ thickness (a spurious jump/flat band). Used at high
        //     frequency only.
        // The crossover (δ/t ≈ 0.08–0.15) is where the two agree, so the blend is
        // continuous — fixing the old discontinuity at F_STATIC_MAX. 'static' forces
        // the static estimator everywhere.
        if (f > 0 && R_total === 0) {
            const edges = lossEdgeMask;
            // Per-edge surface groups (one group = one surface impedance). Without
            // per-face plating there is a single group whose Zs equals the effective
            // surface impedance — identical to the original single-pass behaviour.
            // The static/SIBC loss routines return R_ac ∝ Rs(sigmaRef); scaling each
            // group by Zs_group.re/Rs(sigmaRef) recovers the per-face surface
            // resistance, and the groups sum (the loss is linear over edges).
            const groups = platingPerSide
                ? surf.groups
                : [{ Zs: { re: Zr.re, im: Zr.im }, mask: edges }];
            const sigmaRef = platingPerSide ? (s.sigma_cond ?? 5.8e7) : sigma;
            const deltaRef = Math.sqrt(2 / (omu * sigmaRef));
            const RsRef = 1 / (sigmaRef * Math.min(deltaRef, 1e30));

            let minDim = Infinity;
            for (const c of cr.rects) minDim = Math.min(minDim, c.xmax - c.xmin, c.ymax - c.ymin);
            const useFW = !!fw && lossMethod !== 'static';
            let an = null;
            if (useFW) {
                const k2 = (omega / c0) ** 2;
                an = analyzeTriMode(fw.vRe, fw.vIm, fw.g2Re, fw.g2Im, mesh, fm, k2, f, mesh.epsMap, cr);
            }
            const haveFW = useFW && an && Math.abs(an.P) > 1e-30 && minDim > 0;
            const wFW = haveFW ? Math.max(0, Math.min(1, (0.15 - deltaRef / minDim) / (0.15 - 0.08))) : 0;

            let R_ac = 0, R_dc = 0;
            for (const g of groups) {
                const fRg = RsRef > 0 ? g.Zs.re / RsRef : 1;
                const lossS = staticConductorLoss(cr.rects, f, sigmaRef, mesh, fm, phiEps,
                    Z0, eps_eff_static, eps_d, null, g.mask);
                R_dc = lossS.R_dc;   // mask-independent (full cross-section at sigmaRef)
                let racg = lossS.R_ac * fRg;
                if (haveFW && wFW > 0) {
                    const lossW = solveConductorLoss(cr.rects, f, sigmaRef, mesh, fm, fw.vRe, fw.vIm,
                        fw.g2Re, fw.g2Im, Math.abs(an.P), Z0, mesh.epsMap, g.mask, null, null, this.ctx.wasmSolver);
                    racg = (1 - wFW) * racg + wFW * (lossW.R_ac * fRg);
                }
                R_ac += racg;
            }
            R_total = Math.sqrt(R_dc * R_dc + R_ac * R_ac);
            // Internal (skin) inductance: in the skin regime Zs = Rs(1+j), so the
            // internal reactance equals the AC resistance → L_int ≈ R_ac/ω. Neither the
            // static-field nor the SIBC routine returns a reliable L_int across the
            // transition; this closed form keeps eps_eff (= c²·L·C) and Re(Zc) consistent
            // with the MQS/FDM backends, which include the internal inductance. The √f
            // rolloff (R_ac∝√f) makes it vanish at high f and the cap bounds it near DC.
            L_internal = (omega > 0) ? Math.max(0, Math.min(R_ac / omega, 0.5 * L_external)) : 0;
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
    // Apply the Djordjevic–Sarkar causal dielectric model at frequency f: rebuild the
    // per-triangle eps/loss maps from the frequency-dependent permittivity and re-solve
    // the static field so eps_eff (below F_STATIC_MAX), the full-wave eigensolve (which
    // reads mesh.epsMap), C, and the dielectric-loss energy W_loss all track the causal
    // model. The k2 eigensolve cache is bypassed separately in _modeAtFreq.
    _applyCausal(f) {
        const s = this.solver;
        const fref = this.opts.causalFref ?? 1e9;
        const causalDiel = s.dielectrics.map(d => {
            const er = d.epsilon_r, td = d.tan_delta || 0;
            const rect = { x_min: d.x_min, x_max: d.x_max, y_min: d.y_min, y_max: d.y_max };
            if (Math.abs(er - 1) < 1e-6 || Math.abs(td) < 1e-10) return { ...rect, epsilon_r: er, tan_delta: td };
            const { eps_real, tand_actual } = djordjevic_sarkar(f, er, td, fref);
            return { ...rect, epsilon_r: eps_real, tan_delta: tand_actual };
        });
        const { epsMap, lossMap } = tagMaterials(this.mesh, causalDiel);
        this.mesh.epsMap = epsMap; this.mesh.lossMap = lossMap;
        for (const mode of this.modeNames) {
            const st = this._static[mode];
            const pot = drivePotentials(this.condRect, mode, this.symmetry);
            const phiEps = solveTriStatic(this.mesh, st.fm, epsMap, pot);
            const W_air = st.C0 / (st.kC * eps0);   // geometric, eps-independent
            st.phiEps = phiEps;
            st.eps_eff_static = computeTriEnergy(phiEps, this.mesh, epsMap) / W_air;
            st.W_loss = computeTriEnergy(phiEps, this.mesh, lossMap);
        }
    }

    solveAt(f) {
        if (!this.mesh) this.buildMesh();
        if (this.solver.use_causal_materials) this._applyCausal(f);
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
        // The resampled grid (solver.x/y/Ex/Ey) is now valid — mark the mesh ready so
        // plotting paths gated on mesh_generated (e.g. the geometry-view E-field contour
        // overlay) render. ensure_mesh() is a no-op for the triangular backend, so this
        // flag would otherwise stay false.
        this.solver.mesh_generated = true;
        return result;
    }

    // ---- Mode viewer ----------------------------------------------------------
    // Solve the full-wave eigenproblem at frequency f for the lowest `nev` modes and
    // return a classified list (propagating / evanescent / spurious / null-space), the
    // same way microstrip_viewer_gmsh.html does. Unlike solveAt (which picks only the
    // single quasi-TEM mode), this exposes EVERY converged mode so the UI can list them
    // and plot each one's field. Uses the full domain (the modes backend is built with
    // symmetry:false) and the natural wall BCs — no differential drive, since the
    // eigenmodes are a property of the geometry, not the excitation.
    // Returns { modes:[{idx, g2Re, g2Im, eps_eff|null, overlap, status}], nconv, N, nTris }.
    // Per-mode eigenvectors are cached for getModeField(); eps_eff is -γ²/k² for
    // propagating modes (γ²<0), null otherwise.
    solveModes(f, nev = 4) {
        if (!this.mesh) this.buildMesh();
        if (this.solver.use_causal_materials) this._applyCausal(f);
        const s = this.solver, cr = this.condRect, mesh = this.mesh;
        const st = this._static[this.modeNames[0]];   // any cached static solve seeds the shift
        const { fm, phiEps, eps_eff_static } = st;
        const k2 = (2 * Math.PI * f / c0) ** 2;
        // Radiating first-order ABC on the OPEN walls (=== true → the Robin radiation term
        // in assembleTriFEM), exactly like the reference viewer (microstrip_viewer_gmsh.html:
        // abc={top,left,right:true}). This absorbs outgoing waves so the open structure's
        // bound quasi-TEM mode is found ONCE and cleanly, instead of the spurious
        // near-degenerate duplicates a closed PMC box produces. 'gnd' walls stay PEC. The
        // freedom map (built with 'pmc', same truthiness → same free DOFs) is reused as-is;
        // only the assembled matrices differ (ABC adds an imaginary boundary term).
        const w = cr.wallPEC || {};
        const abc = {};
        if (!w.left) abc.left = true;
        if (!w.right) abc.right = true;
        if (!w.top) abc.top = true;
        if (!w.bottom) abc.bottom = true;
        let fem;
        try { fem = assembleTriFEM(mesh, fm, k2, lossyEpsMap(mesh), abc, cr, null); }
        catch (e) { return { modes: [], nconv: 0, error: String(e && e.message || e) }; }
        const N = fem.N;
        // Static-field "drive" seeds for every physical conductor excitation (single-ended:
        // one; differential: even AND odd). A genuine quasi-TEM mode matches ONE of these
        // (high overlap); a spurious discrete gradient mode matches none. All modes share
        // this fm (symmetry off), so each drive's phiEps projects into the same DOF space.
        const seeds = this.modeNames.map(m => staticToEdgeDofs(this._static[m].phiEps, mesh, fm));
        const seed = seeds[0];
        // Shift-invert centered on the quasi-TEM eigenvalue (-k²·eps_static) converges the
        // fundamental and the nearby higher-order modes first. Widen the Krylov subspace so
        // several modes resolve at once.
        const ncv = Math.min(Math.max(2 * nev + 1, 20), N - 1);
        let res;
        try { res = this.ctx.helpers.solveGeneralized(N, fem.csrA, fem.csrB, [-k2 * eps_eff_static, 0], nev, ncv, seed); }
        catch (e) { return { modes: [], nconv: 0, error: String(e && e.message || e) }; }
        if (!res || !res.nconv) return { modes: [], nconv: 0 };

        const epsMax = Math.max(1, ...s.dielectrics.map(d => d.epsilon_r || 1));
        const nullThresh = k2 * 1e-6;
        // Classification — the reference viewer's rule, which now works because the lossy
        // eigensolve (lossyEpsMap) damps the discrete vector-FEM spurious modes out of the
        // propagating range, so what remains in 0 < ε_eff < ε_r are the genuine quasi-TEM +
        // cavity/higher-order modes (verified to match the reference at 10/50/100 GHz):
        //   • ε_eff > ε_r           → spurious (non-physical, faster-than-the-densest-medium)
        //   • |γ²| ≈ 0              → null space
        //   • γ² < 0                → propagating (ε_eff = −γ²/k²)
        //   • γ² > 0                → evanescent / below cutoff
        // Plus: an OPEN domain (any radiating ABC wall) has no bound modes below the light
        // line — ε_eff < 1 there is the radiation continuum of the truncated domain, not a
        // guided mode — so flag those spurious. An enclosed (all-PEC) box does have them
        // (cavity modes above their cutoff), so they stay propagating.
        const isEnclosed = !(abc.left || abc.right || abc.top || abc.bottom);
        const list = [];
        for (let i = 0; i < res.nconv; i++) {
            const g2Re = res.evalsRe[i], g2Im = res.evalsIm[i];
            const vRe = res.evecsRe.slice(i * N, (i + 1) * N);
            const vIm = res.evecsIm.slice(i * N, (i + 1) * N);
            // overlap of the transverse part with the BEST-matching conductor drive (for
            // display + auto-selecting the quasi-TEM mode; no longer used for classification).
            let overlap = 0;
            for (const sd of seeds) {
                let dot = 0, nS = 0, nV = 0;
                for (let k = 0; k < fm.nFreeTransverse; k++) { dot += sd[k] * vRe[k]; nS += sd[k] ** 2; nV += vRe[k] ** 2; }
                const o = nS > 0 && nV > 0 ? Math.abs(dot) / Math.sqrt(nS * nV) : 0;
                if (o > overlap) overlap = o;
            }
            let status, eps_eff = null;
            if (g2Re < -k2 * epsMax * 1.01) status = 'spurious';
            else if (Math.abs(g2Re) < nullThresh) status = 'nullspace';
            else if (g2Re < 0) {
                eps_eff = -g2Re / k2;
                status = (!isEnclosed && eps_eff < 1.0) ? 'spurious' : 'propagating';
            } else status = 'evanescent';
            list.push({ idx: i, g2Re, g2Im, eps_eff, overlap, status, vRe, vIm });
        }
        list.sort((a, b) => a.g2Re - b.g2Re);
        this._modesState = { fm, list, f, k2 };
        return {
            modes: list.map(({ vRe, vIm, ...m }) => m),   // strip eigenvectors from the summary
            nconv: res.nconv, N, nTris: mesh.nTris,
        };
    }

    // Resample the transverse E-field of the sortedIdx-th mode from the last solveModes()
    // onto a grid for plotting. Returns { x, y, E, Ex, Ey } (see resampleModeField) or null.
    getModeField(sortedIdx) {
        const ms = this._modesState;
        if (!ms || !ms.list[sortedIdx]) return null;
        const m = ms.list[sortedIdx];
        let grid = null;
        try {
            const mr = this.solver.mesher;
            if (mr && typeof mr.generate_mesh === 'function') {
                const [gx, gy] = mr.generate_mesh();
                if (gx && gy && gx.length && gy.length) grid = { x: gx, y: gy };
            }
        } catch { grid = null; }
        return resampleModeField(this.mesh, ms.fm, m.vRe, m.vIm, this.domain,
            { resolution: this.opts.resolution, grid });
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
