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
import { buildOccMeshFromGeometry, estimateOccTriCount, tagMaterials, validateTriMesh, _clipDomain } from './occ_to_mesh.js';
import { shapeArea, shapeBBox, shapeSignedDist, isComplement } from '../shapes.js';
import { buildTriFreedomMap, solveTriStatic, computeTriEnergy, refineTriMesh, refineTriMeshNested,
         markTrianglesForRefinement, triP2Stiffness,
         triCoefficients, lvGrad, leGrad,
         staticToEdgeDofs, analyticSeedDofs, assembleTriFEM, assembleTriFEMDecomposed,
         femFromDecomposition, decompositionClassEps, buildTriGauge, gaugeSeed, gaugeExpand } from './tri_fem.js';
import { staticConductorLoss, solveConductorLoss, computeHtZZMetric,
         projectH, computePoyntingFromProjectedH } from './conductor_loss.js';
import { csqrt } from './fem_core.js';
import { mqsConductorLoss, mqsPecInductance, refineSkinBand } from './mqs_loss.js';
import { checkMeshQuality } from './tri_mesh.js';

// Below this frequency use the static solve. Above it, the full-wave eigenmode
// for the dispersive effective permittivity, anchored to the static solve at
// this frequency (see TriBackend._eigenBias) so the two agree at the switch.
const F_STATIC_MAX = 100e6;
import { calculate_Zrough, calculate_Zrough_layered } from '../surface_roughness.js';
import { resampleStatic, resampleModeField, buildGridFromMesh } from './resample.js';
import { Complex } from '../complex.js';
import { classifyModalDecomposition, halfDomainSymmetry } from '../geometry_symmetry.js';
import { buildPhysicalRLGC } from '../sparameters.js';
import { djordjevic_sarkar } from '../djordjevic_sarkar.js';

const c0 = 299792458;
const eps0 = 8.854187817e-12;
const MU0 = 4 * Math.PI * 1e-7;
const NP_TO_DB = 8.685889638;


// ---- Mesh-size budget (memory parity with the FDM backend) ----
// The full-wave eigensolve is far heavier per mesh entity than the quasi-static FDM
// Poisson solve: it factorizes a COMPLEX matrix over EDGE DOFs (~1.5 per triangle)
// instead of a real one over node DOFs. Per triangle that is roughly
//   1.5 (edges/tri) × 2 (complex doubles) × ~1.4 (Nédélec vs 5-point sparsity) ≈ 4×
// the bytes the FDM solver spends per node. So a triangle budget of maxNodes/4 keeps the
// full-wave memory footprint ABOUT THE SAME as the FDM backend at the same Max Nodes
// setting — letting both backends accept the same max_nodes value (UI: Max Nodes
// (thousands)). Without this the eigensolve can be handed an unbounded mesh (seen: a
// 62k-triangle modes mesh) and the eigensolver heap-abort()s with an opaque assertion.
const FW_NODES_PER_TRI = 4;
function maxTrisForBudget(maxNodes) {
    return Math.max(800, Math.floor((maxNodes || 18000) / FW_NODES_PER_TRI));
}

// Absolute backend ceiling, mirroring the FDM solver's "Problem too large" guard
// (field_solver.js throws above 1 GB). Refuse a solve whose complex matrices would not
// fit the eigensolver's WASM heap (2 GB hard cap) with a clean error instead of an opaque
// eigensolver abort(). bytes ≈ complex (2×) sparse-factorization estimate, same shape as
// the FDM estimate 10·(12·nnz + 20·N). With the triangle budget above this never trips in
// normal use; it backstops direct/programmatic misuse and extreme Max Nodes settings.
const MAX_SOLVE_BYTES = 1e9;
// The eigensolver factorizes C = A - σB, whose sparsity pattern is the union of
// the A and B patterns, and B (with the g2 coupling blocks) is typically ~2.6x
// denser than A. Size the estimate on the denser of the two patterns.
// nnzScale: the tree-cotree gauged pencil (fullwaveMode) carries ~20% more B entries
// than the plain one for the same factor size (measured L+U fill 21.9M vs 22.1M on a
// 12k-triangle mesh), so its estimate is scaled back to keep the guard at parity.
function eigenSolveBytes(N, fem, nnzScale = 1) {
    const nnz = nnzScale * Math.max(fem.csrA.colIdx.length, fem.csrB.colIdx.length);
    return 2 * 10 * (12 * nnz + 20 * N);
}
const GAUGED_NNZ_SCALE = 1 / 1.2;

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
// surface, so when any signal conductor is ACTIVELY plated (σ set AND at least one
// face selected — the top/sides/bottom checkboxes can all be off, which means bare
// metal, matching the FDM backend) use its plating σ/Rq; else base.
function effectiveSurface(solver) {
    let sigma = solver.sigma_cond ?? 5.8e7;
    let rq = solver.rq ?? 0;
    for (const c of solver.conductors) {
        if (c.is_signal && c.plating && c.plating.sigma > 0 &&
            (c.plating.top || c.plating.sides || c.plating.bottom)) {
            sigma = c.plating.sigma;
            rq = c.plating.rq ?? rq;
            break;
        }
    }
    return { sigma, rq };
}

// Per-face plating model. Each conductor face (top / bottom / sides) gets its
// own surface impedance: plated faces use the layered (plating-over-bulk)
// gradient impedance matching the FDM backend's per-face model and bare faces
// and grounded walls use the base metal. Thick-corner plating (geometric wrap
// of side plating onto the bottom face within one plating thickness of each
// corner) is additionally modeled position-aware via zAt. See zAt for FDM
// getZsurf differences.
//
// Two consumers, sharing makePlatingZs:
//   * buildSurfaceGroups (pertrubation path) groups loss edges by surface
//     impedance so the surface integral is summed per group (loss is linear in
//     the per-edge surface resistance). `uniform` true, a single impedance
//     everywhere, the plain single-pass path.
//   * buildFaceZs (MQS path) returns Zs at a face midpoint, so the MQS volume
//     solve can weight each face's smooth current by its own plating impedance.
// Shared surface-impedance model: the bare-metal Zs plus a (rectIndex, face) -> Zs
// resolver. Plated faces use the layered (plating-over-bulk) gradient impedance
// (cached per distinct plating), bare faces use the base metal. Used by both the
// edge-based grouping (perturbation path) and the point-based lookup (MQS path).
function makePlatingZs(solver, condRect, freq) {
    const roles = condRect.rectRoles || [];
    const rects = condRect.rects || [];
    const sigmaBase = solver.sigma_cond ?? 5.8e7;
    const rqBase = solver.rq ?? 0;
    const Zbare = calculate_Zrough(freq, sigmaBase, rqBase);    // Complex (.re/.im)
    const layeredCache = new Map();
    const singleCache = new Map();
    // face is 'top' | 'bottom' | 'sides' | 'all'. 'all' is the shaped-conductor case:
    // a circle has ONE continuous surface, so there is nothing to select between and
    // the plating covers the whole boundary (CoaxSolver sets pl.all).
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
    // Single-layer plating impedance (no bulk underneath) for wrap-around regions.
    const zSingle = (sigma, rq) => {
        const key = `${sigma}|${rq}`;
        let z = singleCache.get(key);
        if (!z) { z = calculate_Zrough(freq, sigma, rq); singleCache.set(key, z); }
        return z;
    };
    // Position-aware resolver: same per-face Zs as zForFace, plus the FDM backend's
    // thick-corner wrap (field_solver.getZsurf): side
    // plating physically wraps onto the bottom face within one plating thickness
    // of each corner (sides plated, bottom not, thick_corners on). Single-layer
    // plating σ with the bulk rq (the surface prep is the original bottom).
    //
    // Two deliberate differences from the FDM getZsurf:
    //   * The FDM's other coverage rule: top-only plating extending down the
    //     sides by one thickness is not implemented here.
    //   * The wrap is a hard step at the query point (edge midpoint / face
    //     midpoint), where the FDM blends fractional cell coverage. Under a
    //     plating-thickness sweep the tri answer therefore steps where the FDM
    //     ramps, the two agree once the wrap is resolved by more than one
    //     element, which the skin-band refinement gives at the corners.
    const zAt = (ri, face, x, y) => {
        const pl = roles[ri] && roles[ri].plating;
        const r = rects[ri];
        if (pl && pl.sigma > 0 && r && !r.shape && (pl.thickness ?? 0) > 0 &&
            face === 'bottom' && pl.sides && !pl.bottom && pl.thick_corners) {
            const d = Math.min(x - r.xmin, r.xmax - x);
            if (d <= pl.thickness) return zSingle(pl.sigma, rqBase);
        }
        return zForFace(ri, face);
    };
    return { Zbare, zForFace, zAt };
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
    const { Zbare, zForFace, zAt } = makePlatingZs(solver, condRect, freq);
    return (x, y, orient) => {
        for (let ri = 0; ri < rects.length; ri++) {
            const r = rects[ri];
            if (r.shape) {
                // Curved surface: orientation carries no face information. A point on
                // the boundary is on THE surface.
                if (Math.abs(shapeSignedDist(r.shape, x, y)) > tol) continue;
                return zForFace(ri, 'all');
            }
            if (x < r.xmin - tol || x > r.xmax + tol || y < r.ymin - tol || y > r.ymax + tol) continue;
            let face = null;
            if (orient === 'h') {
                if (Math.abs(y - r.ymax) < tol) face = 'top';
                else if (Math.abs(y - r.ymin) < tol) face = 'bottom';
            } else if (Math.abs(x - r.xmin) < tol || Math.abs(x - r.xmax) < tol) {
                face = 'sides';
            }
            if (!face) continue;
            return zAt(ri, face, x, y);
        }
        return Zbare;
    };
}

function buildSurfaceGroups(solver, mesh, fm, condRect, baseMask, freq, cache = null) {
    const { nodes, edges, nEdges } = mesh;
    const rects = condRect.rects || [];
    const tol = platingTol(condRect);
    const { Zbare, zForFace, zAt } = makePlatingZs(solver, condRect, freq);
    // 4 slots, not 3: shaped conductors add 'all' (one continuous curved surface).
    const FACES = ['top', 'bottom', 'sides', 'all'];
    const NFACE = FACES.length;

    // Geometric classification: loss edge → (rectIndex, face) or bare (−1).
    // Frequency-INVARIANT, so it's cached per (mesh, fm, baseMask); only the
    // Zs values (and thus the grouping) are re-evaluated per frequency.
    let edgeFace = (cache && cache.clsMesh === mesh && cache.clsFm === fm &&
                    cache.clsMask === baseMask) ? cache.edgeFace : null;
    if (!edgeFace) {
        edgeFace = new Int32Array(nEdges).fill(-1);
        for (let e = 0; e < nEdges; e++) {
            if (!baseMask[e]) continue;
            if (!(fm.isCondEdge && fm.isCondEdge[e])) continue;   // grounded walls → bare
            const n0 = edges[2 * e], n1 = edges[2 * e + 1];
            const x0 = nodes[2 * n0], y0 = nodes[2 * n0 + 1];
            const x1 = nodes[2 * n1], y1 = nodes[2 * n1 + 1];
            for (let ri = 0; ri < rects.length; ri++) {
                const r = rects[ri];
                if (r.shape) {
                    // Both endpoints on the curved boundary => a surface edge. There is
                    // only one face, so no orientation test.
                    if (Math.abs(shapeSignedDist(r.shape, x0, y0)) > tol) continue;
                    if (Math.abs(shapeSignedDist(r.shape, x1, y1)) > tol) continue;
                    edgeFace[e] = ri * NFACE + 3;                                        // all
                    break;
                }
                if (x0 < r.xmin - tol || x0 > r.xmax + tol || x1 < r.xmin - tol || x1 > r.xmax + tol) continue;
                if (y0 < r.ymin - tol || y0 > r.ymax + tol || y1 < r.ymin - tol || y1 > r.ymax + tol) continue;
                let face = -1;
                if (Math.abs(y0 - r.ymax) < tol && Math.abs(y1 - r.ymax) < tol) face = 0;         // top
                else if (Math.abs(y0 - r.ymin) < tol && Math.abs(y1 - r.ymin) < tol) face = 1;    // bottom
                else if ((Math.abs(x0 - r.xmin) < tol && Math.abs(x1 - r.xmin) < tol) ||
                         (Math.abs(x0 - r.xmax) < tol && Math.abs(x1 - r.xmax) < tol)) face = 2;  // sides
                if (face < 0) continue;
                edgeFace[e] = ri * NFACE + face;
                break;
            }
        }
        if (cache) { cache.clsMesh = mesh; cache.clsFm = fm; cache.clsMask = baseMask; cache.edgeFace = edgeFace; }
    }

    const keyOf = (z) => `${z.re.toExponential(6)}|${z.im.toExponential(6)}`;
    const groupIdx = new Map();
    const groupZ = [];
    const edgeGroup = new Int32Array(nEdges).fill(-1);
    for (let e = 0; e < nEdges; e++) {
        if (!baseMask[e]) continue;
        const fid = edgeFace[e];
        // Position-aware Zs (zAt): the thick-corner / top-only wrap regions get the
        // single-layer plating impedance, resolved at the edge midpoint. 'all'
        // (shaped) faces have no wrap geometry, zAt falls through to zForFace.
        const n0 = edges[2 * e], n1 = edges[2 * e + 1];
        const z = fid < 0 ? Zbare
            : zAt((fid / NFACE) | 0, FACES[fid % NFACE],
                  (nodes[2 * n0] + nodes[2 * n1]) / 2, (nodes[2 * n0 + 1] + nodes[2 * n1 + 1]) / 2);
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

// Per-triangle electrostatic energy (legacy refinement metric): ε·∫|∇φ|² over
// each tri.  Kept as the fallback for perElementKelly below. It refines where
// the field is strong, not where the discretization error is.
function perElementEnergy(phi, mesh, epsMap) {
    const { nodes, tris, nTris, triEdges } = mesh;
    const m = new Float64Array(nTris);
    const lp = new Float64Array(6);
    const Sz = new Float64Array(36);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        triP2Stiffness(nodes, v0, v1, v2, Sz);
        lp[0] = phi.phiVertex[v0]; lp[1] = phi.phiVertex[v1]; lp[2] = phi.phiVertex[v2];
        lp[3] = phi.phiEdge[triEdges[3 * t]]; lp[4] = phi.phiEdge[triEdges[3 * t + 1]]; lp[5] = phi.phiEdge[triEdges[3 * t + 2]];
        let e = 0;
        for (let li = 0; li < 6; li++) for (let lj = 0; lj < 6; lj++) e += Sz[li * 6 + lj] * lp[li] * lp[lj];
        m[t] = epsMap[t].re * Math.abs(e);
    }
    return m;
}

// Per-triangle Kelly (flux-jump) error indicator for the P2 static solve:
//   η_T² = Σ_{edges e of T} w_e · h_e · ∫_e [[ε ∂φ/∂n]]² ds
// The exact solution has a continuous normal flux ε∂φ/∂n across element interfaces, so
// the discrete jump measures the local discretization error. Interior
// edges split the contribution between their two triangles (w = 0.5). Edges
// whose midpoint DOF is constrained (fm.edgeNodeF < 0: conductor surfaces and
// interiors, PEC walls) are skipped because the flux jump at a Dirichlet
// surface is the physical surface charge, not an error.
function perElementKelly(phi, mesh, epsMap, fm) {
    const { nodes, tris, triEdges, edges, nTris, nEdges } = mesh;
    const { phiVertex, phiEdge } = phi;
    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    // edge → adjacent triangles
    const edgeT = new Int32Array(2 * nEdges).fill(-1);
    for (let t = 0; t < nTris; t++)
        for (let k = 0; k < 3; k++) {
            const e = triEdges[3 * t + k];
            if (edgeT[2 * e] < 0) edgeT[2 * e] = t; else edgeT[2 * e + 1] = t;
        }

    // Per-triangle barycentric coefficients, computed once (gradAt is called up to
    // 4× per edge below).
    const coeffs = new Array(nTris);
    const coeffOf = t => coeffs[t] ??
        (coeffs[t] = triCoefficients(nodes, tris[3 * t], tris[3 * t + 1], tris[3 * t + 2]).coeff);

    // ∇φ of triangle t's P2 solution at (x, y)
    const gradAt = (t, x, y) => {
        const coeff = coeffOf(t);
        let gx = 0, gy = 0;
        for (let k = 0; k < 3; k++) {
            const p = phiVertex[tris[3 * t + k]];
            if (p !== 0) { const g = lvGrad(coeff, k, x, y); gx += p * g[0]; gy += p * g[1]; }
        }
        for (let k = 0; k < 3; k++) {
            const p = phiEdge[triEdges[3 * t + k]];
            if (p !== 0) {
                const [i, j] = edgeVerts[k];
                const g = leGrad(coeff, i, j, x, y); gx += p * g[0]; gy += p * g[1];
            }
        }
        return [gx, gy];
    };

    // 2-point Gauss on the edge (P2, the normal flux is linear along the edge, so this
    // integrates J^2 exactly on straight edges).
    const GP = [0.5 - 0.5 / Math.sqrt(3), 0.5 + 0.5 / Math.sqrt(3)];
    const m = new Float64Array(nTris);
    for (let e = 0; e < nEdges; e++) {
        if (fm.edgeNodeF[e] < 0) continue;           // Dirichlet-constrained: physical flux, skip
        const t1 = edgeT[2 * e], t2 = edgeT[2 * e + 1];
        if (t1 < 0) continue;
        const n0 = edges[2 * e], n1 = edges[2 * e + 1];
        const x0 = nodes[2 * n0], y0 = nodes[2 * n0 + 1];
        const ex = nodes[2 * n1] - x0, ey = nodes[2 * n1 + 1] - y0;
        const h = Math.hypot(ex, ey);
        if (!(h > 0)) continue;
        const nxu = ey / h, nyu = -ex / h;           // unit normal (orientation cancels in J²)
        const e1 = epsMap ? epsMap[t1].re : 1;
        const e2 = (t2 >= 0 && epsMap) ? epsMap[t2].re : 1;
        let I = 0;
        for (const s of GP) {
            const x = x0 + s * ex, y = y0 + s * ey;
            const g1 = gradAt(t1, x, y);
            let J = e1 * (g1[0] * nxu + g1[1] * nyu);
            if (t2 >= 0) {
                const g2 = gradAt(t2, x, y);
                J -= e2 * (g2[0] * nxu + g2[1] * nyu);
            }
            I += 0.5 * J * J;                        // Gauss weights ½ each
        }
        const eta2 = h * h * I;                      // h_e · ∫ J² ds
        if (t2 >= 0) { m[t1] += 0.5 * eta2; m[t2] += 0.5 * eta2; }
        else m[t1] += eta2;
    }
    return m;
}

// Full-wave eigenmode at frequency f via the Lee-Jin eigensolve. Robust mode
// selection: among physical modes (γ²<0, eps not near cutoff), pick the one with
// the HIGHEST overlap to the static drive (the static field is the quasi-TEM shape) —
// reliable even when several modes cluster at eps ≈ ε_r (homogeneous/stripline fills).
// Returns { eps, vRe, vIm, g2Re, g2Im } for the chosen mode, or null.
// Throws on assembly/eigensolver failure (callers catch and surface a warning);
// returns null when the eigensolve converged but no quasi-TEM candidate passed
// the physicality gates.
function fullwaveMode(ctx, mesh, fm, abc, condRect, epsMap, f, phiEps, eps_static, cache = null) {
    const k2 = (2 * Math.PI * f / c0) ** 2;
    // The system is affine in k² with the ABC term linear in k₀, so the expensive
    // quadrature assembly is done ONCE per (mesh, fm, epsMap, abc) as a
    // decomposition A0 + k²·A1 + j·k0·Ar and combined per frequency in O(nnz).
    // The cache self-validates on object identity (a refinement pass or a causal
    // epsMap rebuild changes the objects → clean recompute).
    // The pencil is assembled in the tree-cotree gauge (buildTriGauge): the gradient
    // part of e_t is an explicit nodal variable, so the O(k^2) terms the quasi-TEM lives
    // on never sit inside the roundoff of the 1/h^2 curl-curl entries. Without it the
    // solve on a finely graded mesh returned nothing below ~200 MHz (the true mode's
    // Ritz residual at the gate) and fell back from LDL^T to LU at every frequency.
    // A rebuilt material map (the causal model re-tags the same layout with new
    // values at every frequency) keeps the decomposition and only changes the
    // per-class permittivities it is combined with (decompositionClassEps).
    let dec = null, gauge, classEps = null;
    const abcKey = JSON.stringify(abc);
    if (cache && cache.dec && cache.mesh === mesh && cache.fm === fm && cache.abcKey === abcKey) {
        classEps = cache.epsMap === epsMap ? cache.classEps : decompositionClassEps(cache.dec, epsMap);
        if (classEps) { dec = cache.dec; gauge = cache.gauge; cache.epsMap = epsMap; cache.classEps = classEps; }
    }
    if (!dec) {
        gauge = (cache && cache.mesh === mesh && cache.fm === fm && cache.gauge)
            ? cache.gauge : buildTriGauge(mesh, fm);
        dec = assembleTriFEMDecomposed(mesh, fm, epsMap, abc, condRect, gauge);
        classEps = dec.classes;
        if (cache) {
            cache.mesh = mesh; cache.fm = fm; cache.epsMap = epsMap; cache.classEps = classEps;
            cache.abcKey = abcKey; cache.dec = dec; cache.gauge = gauge;
        }
    }
    const fem = femFromDecomposition(dec, k2, classEps);
    if (!fem) return null;
    const N = fem.N;
    // Seed the Arnoldi from the static-field projection (the quasi-TEM shape), in the
    // gauged layout. The mode-pick below overlaps every candidate against the same
    // static drive in the original layout (staticSeed).
    const staticSeed = staticToEdgeDofs(phiEps, mesh, fm);
    const seed = gaugeSeed(gauge, phiEps, mesh, fm);
    // Request enough eigenpairs that the quasi-TEM is in the returned set even when it has
    // DISPERSED away from the shift. Shift-invert is centered on -k2·eps_static, so the
    // solver returns the `nev` modes nearest that shift. At high frequency on an inhomogeneous
    // line the quasi-TEM ε_eff disperses well ABOVE eps_static (e.g. a wide GCPW at 40 GHz:
    // eps_static 3.6 but the quasi-TEM has risen to ~5), so with too small an `nev` the only
    // modes near the shift are spurious/lower ones and the pick mode-hops to a wrong, lower ε.
    // nev=8 keeps the dispersed quasi-TEM in the candidate set; the overlap-with-static-drive
    // selection below then picks it (it carries the highest overlap, ~0.7, vs ~0.6 for the
    // spurious neighbours). (Was 4, which mode-hopped a wide-GCPW case at 40 GHz to ε≈2.7.)
    // For the radiating-ABC pick path ncv is already floored at 20 (the non-Hermitian matrix
    // needs the wider Krylov subspace), so raising nev 4→8 there costs nothing — the same
    // 20-vector subspace, just more Ritz pairs extracted. The closed 'pmc' assembly
    // (refinement metric, enclosed structures) pays the modest ncv 9→17.
    const hasRadiatingAbc = abc && Object.values(abc).some(v => v === true);
    // Clamp nev to the problem size (a tiny mesh can have N ≤ 8).
    const nev = Math.max(1, Math.min(8, N - 1));
    const ncv = Math.min(hasRadiatingAbc ? Math.max(2 * nev + 1, 20) : 2 * nev + 1, N - 1);
    // Absolute ceiling (mirrors the FDM "Problem too large" guard): refuse a solve that
    // would overrun the eigensolver heap with a clear error rather than an opaque abort().
    // The triangle budget keeps normal solves well under this; this backstops misuse.
    if (eigenSolveBytes(N, fem, GAUGED_NNZ_SCALE) > MAX_SOLVE_BYTES)
        throw new Error(`Problem too large for the full-wave solver (~${(eigenSolveBytes(N, fem, GAUGED_NNZ_SCALE) / 1e9).toFixed(1)} GB). Reduce Max Nodes or the domain/frequency.`);
    // solveGeneralized throws on solver failure (negative return code) — let it
    // propagate so the caller can warn. nconv === 0 (nothing converged) → null.
    const res = ctx.helpers.solveGeneralized(N, fem.csrA, fem.csrB, [-k2 * eps_static, 0], nev, ncv, seed);
    if (!res || res.nconv <= 0) return null;
    // Pick the quasi-TEM mode by overlap with the static drive (the static field IS the
    // quasi-TEM shape), with a tie-break by proximity to eps_static. A loose eps gate drops
    // near-cutoff / continuum modes (eps ≪ eps_static) that can score a deceptively high
    // overlap. Collect every physical candidate, then select.
    // Physical upper bound for a guided quasi-TEM mode: its phase velocity cannot be slower
    // than light in the densest dielectric, so ε_eff ≤ max(ε_r). A mode above that is a
    // spurious discrete-FEM eigenvector (same rule solveModes uses to classify spurious). This
    // gate matters specifically with the larger nev needed to catch the dispersed quasi-TEM:
    // the wider candidate set can surface a spurious high-ε mode with deceptively high overlap
    // (e.g. a near-degenerate broadside stripline produced ε≈15 with εr≤6) that would otherwise
    // win the pick. The 1.01 slack mirrors solveModes.
    let epsRMax = 1;
    for (let t = 0; t < epsMap.length; t++) if (epsMap[t].re > epsRMax) epsRMax = epsMap[t].re;
    const cands = [];
    for (let i = 0; i < res.nconv; i++) {
        if (res.evalsRe[i] >= 0) continue;                 // physical: γ² < 0
        const epsMode = -res.evalsRe[i] / k2;
        if (epsMode < 0.5 * eps_static) continue;          // drop near-cutoff / lateral-continuum modes
        if (epsMode > epsRMax * 1.01) continue;            // drop spurious super-dense modes (ε_eff > max εr)
        const vR = gaugeExpand(gauge, res.evecsRe.subarray(i * N, (i + 1) * N));
        let dot = 0, nS = 0, nV = 0;
        for (let k = 0; k < fm.nFreeTransverse; k++) { dot += staticSeed[k] * vR[k]; nS += staticSeed[k] ** 2; nV += vR[k] ** 2; }
        const ovl = nS > 0 && nV > 0 ? Math.abs(dot) / Math.sqrt(nS * nV) : 0;
        globalThis.__TRI_DEBUG__ && console.log(`    cand ${i}: eps=${epsMode.toFixed(4)} ovl=${ovl.toFixed(3)} (target ${eps_static.toFixed(4)})`);
        if (ovl > 0.3) cands.push({ i, eps: epsMode, ovl, vR });   // require overlap > 0.3
    }
    if (!cands.length) return null;
    // Among the overlap-COMPETITIVE candidates (within 0.05 of the best overlap), pick the one
    // CLOSEST to eps_static. The quasi-TEM disperses smoothly from eps_static, so it stays near
    // it, whereas a box/cavity resonance can transiently MATCH the overlap (a near-tie) while
    // sitting far from eps_static. Pure max-overlap then coin-flips between them right at the
    // crossing — so the reported eps jumps to the resonance depending on the mesh. The tie-break
    // tracks the quasi-TEM robustly, and reduces to max-overlap whenever there is no tie.
    const maxOvl = Math.max(...cands.map(c => c.ovl));
    const competitive = cands.filter(c => c.ovl >= maxOvl - 0.05);
    let pick = competitive[0];
    for (const c of competitive) if (Math.abs(c.eps - eps_static) < Math.abs(pick.eps - eps_static)) pick = c;
    const bestIdx = pick.i, bestEps = pick.eps, bestOvl = pick.ovl;
    // Ambiguity: a competitive runner-up carries comparable overlap but a very different ε
    // (near-degenerate modes). The pick is now robust, but the closeness is still worth flagging.
    const runner = cands.filter(c => c.i !== bestIdx).sort((a, b) => b.ovl - a.ovl)[0] || { ovl: 0, eps: bestEps };
    const altOvl = runner.ovl, altEps = runner.eps;
    // Ambiguity: a runner-up mode carries comparable overlap with the static drive but a
    // very different ε. That means the quasi-TEM has fragmented across degenerate modes
    // (inhomogeneous fill, esp. two-ground / differential geometries) and the single-mode
    // pick is unreliable — the eigenvalue is one fragment, not the true quasi-TEM ε. See
    // the broadside-stripline open-side case. Distinct from a clean near-degenerate pair
    // that AGREES on ε (e.g. embedded microstrip), which is harmless.
    const ambiguous = altOvl >= 0.4 && altOvl >= 0.6 * bestOvl
        && Math.abs(altEps - bestEps) / Math.max(bestEps, 1e-9) >= 0.10;
    return {
        eps: -res.evalsRe[bestIdx] / k2,
        vRe: pick.vR,
        vIm: gaugeExpand(gauge, res.evecsIm.subarray(bestIdx * N, (bestIdx + 1) * N)),
        g2Re: res.evalsRe[bestIdx], g2Im: res.evalsIm[bestIdx],
        ambiguous, bestOvl, altOvl, bestEps, altEps,
    };
}

// Rectangular waveguide
//
// A hollow guide has no driven conductor, so there is no static solve to seed the Arnoldi
// from, no eps_eff_static to centre the shift on, and no quasi-TEM mode to overlap
// against. What replaces all three is the analytic cutoff wavenumber: kc is purely
// geometric, so the target eigenvalue γ² = kc² - k0²*eps_r is known before the solve. A
// shift placed exactly there converges in a handful of Arnoldi steps and cannot mode-hop.

// Largest real permittivity on the mesh.
function maxEpsRe(epsMap) {
    let m = 1;
    for (let t = 0; t < epsMap.length; t++) if (epsMap[t].re > m) m = epsMap[t].re;
    return m;
}

// Loss tangent of the (homogeneous) fill, lossMap carries eps_r*tand per triangle.
function wgFillTand(mesh) {
    let t = 0;
    for (let i = 0; i < mesh.nTris; i++) {
        const e = mesh.epsMap[i].re, l = mesh.lossMap[i].re;
        if (e > 0 && l / e > t) t = l / e;
    }
    return t;
}

// Guide dimensions and the Zpv/Z_TE ratio, taken from the meshed domain so they can never
// disagree with what was actually solved.
//
// κ = Zpv/Z_TE = 2*(extent along E)/(extent along the sinusoid) = 2*min/max, for either
// orientation: TE10 has E along y over b and varies along x over a (κ = 2b/a). TE01 has E
// along x over a and varies along y over b (κ = 2a/b).
function wgGeom(condRect) {
    const a = condRect.xmax_domain - condRect.xmin_domain;
    const b = condRect.ymax_domain - condRect.ymin_domain;
    return { a, b, kappa: 2 * Math.min(a, b) / Math.max(a, b) };
}

// Analytic fundamental-mode field of the meshed box, for the Arnoldi seed.
//   broad wall along x -> TE10, E = ŷ·sin(π(x−X0)/a)
//   broad wall along y -> TE01, E = x̂·sin(π(y−Y0)/b)
function wgFundamentalEfun(condRect) {
    const X0 = condRect.xmin_domain, Y0 = condRect.ymin_domain;
    const { a, b } = wgGeom(condRect);
    if (a >= b) return (x) => [0, Math.sin(Math.PI * (x - X0) / a)];
    return (x, y) => [Math.sin(Math.PI * (y - Y0) / b), 0];
}

// Fundamental guided mode of an enclosed PEC cross-section at frequency f.
//
// Unlike fullwaveMode the pick is not an overlap test but an ordering: γ² = kc² - k0²*eps_r
// is monotone increasing in kc, so the smallest γ² is always the lower cutoff, the
// fundamental at any frequency, including above the second cutoff where several modes
// propagate. That is the whole of "solve only the lowest mode", with no per-frequency
// special case. (solveModes sorts on the same key.)
//
// Returns { kc, vRe, vIm, g2Re, g2Im } or null when nothing physical converged.
//
// No assembly cache, unlike fullwaveMode. This medium runs exactly one eigensolve per
// (mesh, freedom map), one per refinement pass, then one on the final mesh and every
// call therefore arrives with a freshly built fm, so there is nothing a cache could hit.
// The sweep itself never comes back here at all: _waveguideModeAtFreq propagates beta(f)
// analytically from the single cached eigenvector.
function waveguideEigen(ctx, mesh, fm, abc, condRect, epsMap, f, seed, kcAnalytic) {
    const k2 = (2 * Math.PI * f / c0) ** 2;
    const erMax = maxEpsRe(epsMap);
    // Same affine-in-k² decomposition as fullwaveMode.
    const dec = assembleTriFEMDecomposed(mesh, fm, epsMap, abc, condRect);
    const fem = femFromDecomposition(dec, k2);
    if (!fem) return null;
    const N = fem.N;
    if (eigenSolveBytes(N, fem) > MAX_SOLVE_BYTES)
        throw new Error(`Problem too large for the full-wave solver (~${(eigenSolveBytes(N, fem) / 1e9).toFixed(1)} GB). Reduce Max Nodes or the guide dimensions.`);
    // The shift IS the analytic target eigenvalue for a homogeneous fill.
    const sigma = kcAnalytic * kcAnalytic - k2 * erMax;
    const nev = Math.max(1, Math.min(6, N - 1));
    const ncv = Math.min(2 * nev + 1, N - 1);
    const res = ctx.helpers.solveGeneralized(N, fem.csrA, fem.csrB, [sigma, 0], nev, ncv, seed || null);
    if (!res || res.nconv <= 0) return null;
    // Physicality gates, same thresholds solveModes uses: drop the discrete null space and
    // anything faster than light in the densest medium.
    const nullThresh = k2 * 1e-6;
    let best = -1, bestG2 = Infinity;
    for (let i = 0; i < res.nconv; i++) {
        const g2 = res.evalsRe[i];
        if (Math.abs(g2) < nullThresh) continue;
        if (g2 < -k2 * erMax * 1.01) continue;
        if (g2 < bestG2) { bestG2 = g2; best = i; }
    }
    if (best < 0) return null;
    return {
        kc: Math.sqrt(Math.max(bestG2 + k2 * erMax, 0)),
        vRe: res.evecsRe.slice(best * N, (best + 1) * N),
        vIm: res.evecsIm.slice(best * N, (best + 1) * N),
        g2Re: bestG2, g2Im: res.evalsIm[best],
    };
}

// ---- Anchor caches (interpolation between exact solves along a sweep) ----
// Three per-mode scalars are smooth, slowly varying functions of ln f on the MQS
// loss path: the full-wave eps_d (quasi-TEM dispersion) and the MQS R and
// L_internal. Sweep points between exact anchors interpolate them with a
// monotone piecewise cubic (PCHIP, Fritsch-Carlson slopes, no overshoot) over
// x = ln f once a leave-one-out check says the local interpolation is within
// tolerance. A cache with `log` set stores ln(value), so power laws such as the
// skin-regime R ~ sqrt(f) are straight lines the cubic reproduces exactly and
// the check only reacts to genuine curvature (the skin transition).
function pchipEval(xs, ys, x) {
    const n = xs.length;
    if (n < 2 || x < xs[0] || x > xs[n - 1]) return null;
    // slopes
    const h = new Array(n - 1), del = new Array(n - 1);
    for (let i = 0; i < n - 1; i++) { h[i] = xs[i + 1] - xs[i]; del[i] = (ys[i + 1] - ys[i]) / h[i]; }
    const d = new Array(n).fill(0);
    for (let i = 1; i < n - 1; i++) {
        if (del[i - 1] * del[i] > 0) {
            const w1 = 2 * h[i] + h[i - 1], w2 = h[i] + 2 * h[i - 1];
            d[i] = (w1 + w2) / (w1 / del[i - 1] + w2 / del[i]);
        }
    }
    const endSlope = (h0, h1, d0, d1) => {
        let dd = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if (dd * d0 <= 0) dd = 0;
        else if (d0 * d1 <= 0 && Math.abs(dd) > 3 * Math.abs(d0)) dd = 3 * d0;
        return dd;
    };
    if (n >= 3) {
        d[0] = endSlope(h[0], h[1], del[0], del[1]);
        d[n - 1] = endSlope(h[n - 2], h[n - 3], del[n - 2], del[n - 3]);
    } else { d[0] = d[1] = del[0]; }
    // bracket + cubic Hermite
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
    const hh = xs[lo + 1] - xs[lo], t = (x - xs[lo]) / hh;
    const t2 = t * t, t3 = t2 * t;
    return ys[lo] * (2 * t3 - 3 * t2 + 1) + hh * d[lo] * (t3 - 2 * t2 + t)
         + ys[lo + 1] * (-2 * t3 + 3 * t2) + hh * d[lo + 1] * (t3 - t2);
}

// Insert an exact anchor (keeps xs ascending; ignores duplicate frequencies and
// values a log cache cannot hold).
function dispersionInsert(dc, f, v) {
    const y = dc.log ? Math.log(v) : v;
    if (!Number.isFinite(y)) return;
    const x = Math.log(f);
    let lo = 0, hi = dc.xs.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (dc.xs[m] < x) lo = m + 1; else hi = m; }
    if (lo < dc.xs.length && Math.abs(dc.xs[lo] - x) < 1e-12) return;
    dc.xs.splice(lo, 0, x); dc.ys.splice(lo, 0, y);
}

// Interpolate at f, or null when not (yet) trustworthy. Trust requires >= 5
// anchors, f strictly interior, and a leave-one-out check on each interior
// anchor bracketing f: removing that anchor and predicting it from the rest must
// land within tol (relative). The check only says something about the gap it
// is made across when the removed anchor sits well inside that gap: a cubic is
// always accurate next to its own end points, so an anchor at the edge of a wide
// gap (a frequency list ascending towards an anchor solved first at f_max)
// would pass while the interpolant is wrong in the middle. The two intervals
// meeting at the removed anchor may differ by at most 2x, which is what
// bisection produces; anything more lopsided fails outright, which sends the
// point to an exact solve and narrows the gap.
function dispersionInterp(dc, f, tol) {
    const n = dc.xs.length;
    if (n < 5) return null;
    const x = Math.log(f);
    if (x <= dc.xs[0] || x >= dc.xs[n - 1]) return null;
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (dc.xs[m] <= x) lo = m; else hi = m; }
    for (const j of new Set([Math.max(1, lo), Math.min(n - 2, hi)])) {
        const tj = (dc.xs[j] - dc.xs[j - 1]) / (dc.xs[j + 1] - dc.xs[j - 1]);
        if (tj < 1 / 3 - 1e-9 || tj > 2 / 3 + 1e-9) return null;
        const xsL = dc.xs.slice(0, j).concat(dc.xs.slice(j + 1));
        const ysL = dc.ys.slice(0, j).concat(dc.ys.slice(j + 1));
        const pred = pchipEval(xsL, ysL, dc.xs[j]);
        if (pred === null) return null;
        const err = dc.log ? Math.abs(pred - dc.ys[j]) : Math.abs(pred - dc.ys[j]) / Math.abs(dc.ys[j]);
        if (!(err <= tol)) return null;
    }
    const y = pchipEval(dc.xs, dc.ys, x);
    return dc.log ? Math.exp(y) : y;
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

// --- Helpers for general (asymmetric) two-conductor modal analysis -----------
// Linear combination a·A + b·B of two static field objects (every component is
// linear in the conductor drive, so the field for a combined drive is the same
// combination of the per-trace fields). Used to synthesise the modal-eigenvector
// field from the per-trace solves without re-solving.
function combineStatic(a, A, b, B) {
    const comb = (x, y) => {
        if (!x || !y) return x || y || null;
        const out = new Float64Array(x.length);
        for (let i = 0; i < x.length; i++) out[i] = a * x[i] + b * y[i];
        return out;
    };
    return {
        phiVertex: comb(A.phiVertex, B.phiVertex),
        phiEdge: comb(A.phiEdge, B.phiEdge),
    };
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
        // Source-free waveguide medium: no driven conductor, so no static solve exists.
        // Every stage that would normally consume _prepareStatic takes a waveguide branch.
        this._isWG = !!(solver && solver.mode_type === 'waveguide');
        this._wg = null;       // cached fundamental eigenmode (see _prepareWaveguide)
    }

    get modeNames() { return this.solver.is_differential ? ['odd', 'even'] : ['single']; }

    // Cap a bulk element size to a fraction of the shortest in-medium wavelength at the
    // mode-viewer frequency (opts.modesFreq), so high-frequency cavity / higher-order modes
    // are actually resolved. Without this, a genuine high-f mode is under-resolved on the
    // geometry-scale mesh, drifts between the fine and coarse solves, and the mesh-convergence
    // test mis-flags it as spurious. No effect in the quasi-TEM regime (λ/N ≫ geometry hCoarse).
    // opts.wavelengthDensity sets the cells-per-wavelength N (default 12); lower values trade
    // higher-order-mode reliability for a much smaller mesh (nTris ∝ N²).
    _wavelengthCap(h) {
        const f = this.opts.modesFreq;
        if (!f) return h;
        const nLambda = this.opts.wavelengthDensity > 0 ? this.opts.wavelengthDensity : 12;
        const epsMax = Math.max(1, ...this.solver.dielectrics.map(d => d.epsilon_r || 1));
        const lamMin = c0 / (f * Math.sqrt(epsMax));
        return Math.min(h, lamMin / nLambda);
    }

    // Build mesh (with symmetry + adaptive refinement), then solve & cache the
    // (frequency-independent) static problem per mode and resample fields.
    // onProgress({iteration, max_iterations, energy_error, param_error, nodes_x,
    // nodes_y}) is called once per adaptive refinement pass (same shape the
    // rectilinear backend reports, so the UI shows real passes for both).
    // shouldStop: optional callback polled between refinement passes (the UI Stop
    // button). Like the FDM backend, a stop request ends refinement gracefully —
    // the current mesh is kept and the solve proceeds on it.
    async buildMesh(onProgress = null, shouldStop = null) {
        const s = this.solver;
        const dom = { x_min: -s.domain_width / 2, x_max: s.domain_width / 2,
                      y_min: -s.t_gnd, y_max: s.domain_height };
        this.domain = dom;
        // Use a half-domain symmetry solve when the geometry is mirror-symmetric (the
        // shared halfDomainSymmetry disables symmetry for broadside).
        // A medium may also veto the half domain outright. A rectangular waveguide is
        // mirror-symmetric and its fundamental has a magnetic wall at x=0, so the half
        // solve would be valid, but resampleModeField has no parity/mirroring support the
        // way resampleStatic does, so the plotted mode field would be blank for x < 0.
        // (Routed through the solver, not tri_opts: app_solver overwrites tri_opts wholesale
        // after construction. MicrostripSolver also sets it for {symmetry: false}.)
        this.symmetry = (this.opts.symmetry ?? true) && s.tri_symmetry !== false &&
            halfDomainSymmetry(s.conductors, s.dielectrics, s.domain_width, s.is_differential).ok;

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
        //
        // The thickness term uses 2t, not t, when the MQS loss path applies:
        // the base mesh then only needs ~one element across the conductor
        // thickness, the skin-depth resolution for the loss solve comes from
        // refineSkinBand, not from here. (Coarser than 2t the sweep gets
        // slower, skin-band refinement from a too-coarse base grows a bigger
        // skin mesh.) The perturbation/SIBC path integrates the eigenmode's
        // surface |H|^2 on this mesh, and its corner-dominated R is
        // mesh-sensitive, so replicate the useMQS gate (see _modeAtFreq)
        // pre-mesh: wall-absorption via _clipDomain instead of cr.rectRoles,
        // and keep the proven 1.5t sizing when MQS won't apply.
        //
        // On the MQS path the thickness term is also floored at a quarter of
        // the smaller of the trace width and twice the clearance to the nearest
        // other conductor or PEC wall (gapRef). A thin conductor otherwise gets
        // its whole perimeter meshed at the thickness scale for nothing: the
        // static and eigen solves need the edge singularity, which the adaptive
        // passes refine, and the loss solve re-resolves the conductor through
        // the skin band, which copes with the flat sliver elements gmsh puts
        // inside a rect thinner than its surface size (a 1 um trace on the
        // 35 um sizing: R within 0.1% of the thickness-sized mesh from 100 MHz
        // to 100 GHz on a tenth of the base triangles). The floor tracks the
        // thick-copper sizing, where 2t is near w/4 anyway, so the base mesh of
        // a thin conductor costs about what the same trace in thick copper does.
        // The perturbation path keeps the thickness sizing: its surface |H|^2
        // integral is only as good as the base mesh's edge resolution.
        const lmPre = this.opts.lossMethod ?? 'auto';
        const clip = _clipDomain(dom, s.conductors, s.boundaries, s.domain_width * 1e-9);
        const survives = c => Math.min(c.x_max, clip.X1) - Math.max(c.x_min, clip.X0) > 0
                          && Math.min(c.y_max, clip.Y1) - Math.max(c.y_min, clip.Y0) > 0;
        // MQS handles explicit ground rects as passive (C = 0) return
        // conductors, only shaped conductors do (coax), and
        // a surviving signal rect must exist to drive.
        const shapedPre = s.conductors.some(c => c.shape);
        const sigCondPre = s.conductors.some(c => c.is_signal && survives(c));
        // A differential pair without the symmetry-plane mode walls runs the
        // per-conductor-drive MQS (modeCurrents) when both polarity groups
        // are in the meshed domain, broadside/asymmetric pairs and
        // symmetry-disabled solves included. Only when a polarity group is
        // missing does it fall back to perturbation, whose corner-dominated R
        // needs the finer 1.5t sizing.
        const diffMultiPre = s.is_differential
            && s.conductors.some(c => c.is_signal && (c.polarity || 1) > 0 && survives(c))
            && s.conductors.some(c => c.is_signal && (c.polarity || 1) < 0 && survives(c));
        const mqsPre = ((lmPre === 'mqs' && !shapedPre && sigCondPre)
            || (lmPre === 'auto' && (this.symmetry || diffMultiPre) && !shapedPre && sigCondPre))
            && (!s.is_differential || this.symmetry || diffMultiPre);
        // A medium may supply its own base sizing (coax: derived from the conductor
        // radii, since w/t have no meaning for a round conductor).
        const hints = s.tri_mesh_hints || {};
        // Clearance from the signal rects to the nearest other conductor rect or
        // PEC wall. Absorbed wall slabs do not survive the clip, the wall they
        // became does. Touching conductors give 0, which disables the floor.
        let gapRef = Infinity;
        const rectsPre = s.conductors.filter(c => !c.shape && survives(c));
        for (const a of rectsPre) {
            if (!a.is_signal) continue;
            for (const b of rectsPre) {
                if (b === a) continue;
                const dx = Math.max(0, Math.max(a.x_min, b.x_min) - Math.min(a.x_max, b.x_max));
                const dy = Math.max(0, Math.max(a.y_min, b.y_min) - Math.min(a.y_max, b.y_max));
                gapRef = Math.min(gapRef, Math.hypot(dx, dy));
            }
            if (clip.wallPEC.left) gapRef = Math.min(gapRef, a.x_min - clip.X0);
            if (clip.wallPEC.right) gapRef = Math.min(gapRef, clip.X1 - a.x_max);
            if (clip.wallPEC.bottom) gapRef = Math.min(gapRef, a.y_min - clip.Y0);
            if (clip.wallPEC.top) gapRef = Math.min(gapRef, clip.Y1 - a.y_max);
        }
        const hFineFloor = (mqsPre && gapRef > 0) ? Math.min(wRef, 2 * gapRef) / 4 : 0;
        let hFine = this.opts.hFine ?? hints.hFine
            ?? Math.min(Math.max((mqsPre ? 2 : 1) * tAbs, hFineFloor), wRef / 4) * 1.5;
        let hCoarse = this._wavelengthCap(this.opts.hCoarse ?? hints.hCoarse ?? (dom.y_max - dom.y_min) / 5);
        // Triangle budget derived from the Max Nodes setting (memory parity with the FDM
        // backend — see maxTrisForBudget). Cap the INITIAL mesh here: element count ∝ 1/h²,
        // so if the feature-/wavelength-sized start overshoots (thin trace → tiny hFine, or a
        // high modesFreq → fine bulk), coarsen hFine/hCoarse by √(nTris/budget) and rebuild.
        // The refinement loop below also stops at this budget; capping the initial mesh too is
        // what prevents the very first eigensolve from running on an unbounded mesh.
        const maxTris = maxTrisForBudget(this.opts.maxNodes ?? 18000);
        let mesh, prevAtt = null;   // previous attempt's (h, nTris), for the scaling fit
        // Passive-ground surface-size relaxation: when the initial mesh
        // overruns the budget, first coarsen the size painted on passive-ground
        // curves (gndScale, see occ_to_mesh gndSizeScale) and only fall back to
        // global hFine/hCoarse coarsening once that stops paying. Rationale:
        // a large interior ground can cost more triangles than the entire
        // signal region. Coarsening hFine globally to absorb it left the trace
        // with elements bigger than the trace. Signal resolution is what the
        // solved quantities live on, so it is the last thing the budget may
        // take.
        //
        // The two mechanisms get separate pass budgets. Sharing one counter let a
        // geometry spend the whole budget relaxing grounds and then hard-stop with
        // an over-budget mesh, never having coarsened globally at all and maxTris
        // is a memory guard not a preference. Worst case is now
        // COARSEN_MAX_PASSES + GND_MAX_PASSES + 1 gmsh runs, and a geometry with no
        // relaxable ground curves builds exactly the same meshes as before.
        const COARSEN_MAX_PASSES = 4, GND_MAX_PASSES = 3;
        let gndScale = 1, gndStalled = false, gndPrevN = null;
        let coarsenPasses = 0, gndPasses = 0;
        // Accept a near-budget mesh rather than pay a rebuild to shave it (see the loop).
        const OVERSHOOT_OK = 1.3;
        const surfScaleB = this.opts.occSurfScale ?? 0.35;
        const occBase = {
            conductors: s.conductors, dielectrics: s.dielectrics,
            domain: dom, boundaries: s.boundaries, symmetry: this.symmetry,
            // Non-rectangular meshed domain (coax: the dielectric disk itself).
            domainShape: s.domain_shape || null,
            // Conductor interiors are meshed only for the MQS volume eddy-current
            // solve. mqsPre is the same applicability test _modeAtFreq uses. On
            // the perturbation path those elements carry no DOFs, every
            // node/edge inside metal is PEC, so meshing them is pure cost. Coax
            // always lands here: its shaped (non-rectangular) conductors rule
            // MQS out.
            meshConductorInterior: mqsPre,
            occSurfScale: this.opts.occSurfScale, gradeRate: this.opts.gradeRate,
            gmshOptions: this.opts.gmshOptions,
        };
        // Pre-mesh budget sizing
        //
        // Predict the triangle count of the feature-sized start from the size-field
        // statistics (estimateOccTriCount: the OCC model is built and classified, but
        // not meshed, ~0.1 s) and apply the loop's own relaxation policy analytically,
        // passive grounds first, then global, before the first gmsh run. The estimate
        // is good to ~±30%, so it only takes over when the start is far enough over
        // budget (PRESIZE_MIN_OVERSHOOT) that the loop would need two or more extra
        // builds. A moderately over-budget start keeps the loop's trial coarsening
        // and its exact meshes (mesh-sensitive quantities such as the perturbation
        // loss were tuned on them), and a start that fits is untouched. The loop
        // stays as the safety net for a misestimate either way.
        // opts.preMeshSizing = false disables it (diagnostics).
        const PRESIZE_MIN_OVERSHOOT = 2.5;
        if (this.opts.preMeshSizing !== false) try {
            const pre = buildOccMeshFromGeometry(this.ctx.G, { ...occBase, hFine, hCoarse, statsOnly: true });
            const est = (hF, hC, gS) => estimateOccTriCount(pre.sizeStats, {
                sizeMin: hF * surfScaleB, sizeMinGnd: Math.min(hF * surfScaleB * Math.max(1, gS), hC),
                hCoarse: hC, gradeRate: this.opts.gradeRate ?? 0.35 });
            // Aim at the budget itself: the loop accepts up to OVERSHOOT_OK×, and an
            // over-estimate (the common error direction) then lands under it.
            const target = maxTris;
            const bisect = (f) => {   // largest x in [0, 1] with f(x) > target, log-scale
                let lo = 0, hi = 1;
                for (let k = 0; k < 40; k++) { const mid = (lo + hi) / 2; if (f(mid) > target) lo = mid; else hi = mid; }
                return hi;
            };
            const est0 = est(hFine, hCoarse, 1);
            if (est0 > PRESIZE_MIN_OVERSHOOT * maxTris) {
                // Ground relaxation up to the loop's cap (the ground field is inert at
                // hCoarse): the smallest scale that fits, or the cap.
                const gndCap = hCoarse / (hFine * surfScaleB);
                if (pre.nGndCurves > 0 && gndCap > 1) {
                    if (est(hFine, hCoarse, gndCap) > target) gndScale = gndCap;
                    else gndScale = Math.exp(bisect(x => est(hFine, hCoarse, Math.exp(x * Math.log(gndCap)))) * Math.log(gndCap));
                }
                // Global coarsening (hFine and hCoarse together, as the loop does),
                // bounded at 64x: beyond that the geometry is so far over budget that
                // the loop's measured exponent is the better guide.
                if (est(hFine, hCoarse, gndScale) > target) {
                    const lnMax = Math.log(64);
                    const x = est(hFine * 64, hCoarse * 64, gndScale) > target ? 1
                        : bisect(x => est(hFine * Math.exp(x * lnMax), hCoarse * Math.exp(x * lnMax), gndScale));
                    const f = Math.exp(x * lnMax);
                    hFine *= f; hCoarse *= f;
                }
                if (globalThis.__OCC_DEBUG__)
                    console.error('[tri] pre-mesh sizing', { est0, est: est(hFine, hCoarse, gndScale), maxTris, hFine, hCoarse, gndScale });
            }
        } catch (e) {
            // Estimate unavailable: the loop below sizes by trial as before.
            globalThis.__OCC_DEBUG__ && console.error('[tri] pre-mesh sizing failed', e);
        }
        for (;;) {
            mesh = buildOccMeshFromGeometry(this.ctx.G, { ...occBase, hFine, hCoarse, gndSizeScale: gndScale });
            if (globalThis.__OCC_DEBUG__)
                console.error('[tri] initial mesh attempt', { nTris: mesh.nTris, maxTris, hFine, hCoarse, gndScale, coarsenPasses, gndPasses,
                    est: estimateOccTriCount(mesh.sizeStats, mesh.sizeStats) });
            if (mesh.nTris <= maxTris || coarsenPasses >= COARSEN_MAX_PASSES) break;
            // Accept a near-budget mesh before any relaxation/coarsening. This is the
            // same OVERSHOOT_OK acceptance as below, hoisted so a mesh that today gets
            // accepted at, say, 1.2x budget is still accepted identically rather than
            // paying a rebuild with relaxed grounds (keeps the meshes of every family
            // that never over-runs badly bit-for-bit unchanged).
            if (mesh.nTris < maxTris * OVERSHOOT_OK) break;
            // Ground relaxation first (only when there are relaxable ground curves and
            // headroom below the hCoarse cap). Element count along a painted curve is
            // ~1/size, so scale by the overshoot, +5% to land under. If a relaxation
            // pass cut the mesh by less than 10%, the grounds were not the cost, stop
            // relaxing and let the global path take over.
            const gndCap = hCoarse / (hFine * surfScaleB);
            if (!gndStalled && gndPasses < GND_MAX_PASSES && mesh.nGndCurves > 0 && gndScale < gndCap) {
                if (gndPrevN !== null && mesh.nTris > gndPrevN * 0.9) {
                    gndStalled = true;
                } else {
                    gndPrevN = mesh.nTris;
                    gndScale = Math.min(gndScale * (mesh.nTris / maxTris) * 1.05, gndCap);
                    gndPasses++;
                    continue;
                }
            }
            // Element count vs element size is NOT the uniform-fill nTris ∝ 1/h². The
            // conductor-surface size field makes the mesh a graded band around the
            // metal, so the true exponent is shallower — measured ~0.75 on coax, where
            // assuming 2 under-coarsens so badly that the loop burns every coarsening
            // pass, never reaches the budget, and ends up SLOWER with a SMALL "Max Nodes"
            // than a large one (7.6 s vs 4.1 s). Learn the exponent from the last two
            // attempts instead; the first step has no data and keeps the 1/h² guess.
            let p = 2;
            if (prevAtt) {
                const e = Math.log(prevAtt.n / mesh.nTris) / Math.log(hFine / prevAtt.h);
                if (Number.isFinite(e) && e > 0.3) p = Math.min(e, 3);
            }
            const factor = Math.pow(mesh.nTris / maxTris, 1 / p) * 1.05;   // +5% to land under
            // Backstop: give up if coarsening has stopped paying (some geometries floor
            // out). The companion OVERSHOOT_OK acceptance, accept a
            // mesh already near the budget rather than pay a whole gmsh run to shave
            // it is hoisted above the ground relaxation. Its test is on the
            // overshoot, not on `factor`: factor is the p-th root, so the same cutoff
            // there would mean a different budget overrun for every fitted exponent (a
            // `factor < 1.2` bail lets p=3 through at 1.5x budget while p=2 stops at
            // 1.3x). maxTris is a memory guard, so what it may be exceeded by has to
            // be stated in its own units.
            if (mesh.nTris > (prevAtt ? prevAtt.n : Infinity) * 0.9) break;
            prevAtt = { h: hFine, n: mesh.nTris };
            hFine *= factor; hCoarse *= factor;
            coarsenPasses++;
        }
        this.condRect = mesh.condRect;
        // A shaped domain does not fill the solver's bounding box, so tighten the
        // plotting/resample box onto the actual body. Taken from the FULL polygon, not
        // mesh.meshedDomain — under symmetry the latter is only the x >= 0 half, while
        // buildGridFromMesh(mirrorX) reconstructs the whole cross-section.
        if (s.domain_shape) {
            const bb = shapeBBox(s.domain_shape);
            this.domain = { x_min: bb.xmin, x_max: bb.xmax, y_min: bb.ymin, y_max: bb.ymax };
        }

        // Loss routines read per-rect metadata: .symmetry scales the half-domain
        // integral, .xmin_domain locates the symmetry plane (corner/edge exclusion —
        // rects[0] is what the loss routines receive, and without this the plane
        // test compared against undefined and never fired), .is_signal separates
        // signal from ground rects in the DC-resistance area.
        this.condRect.rects.forEach((r, i) => {
            r.symmetry = this.condRect.symmetry;
            r.xmin_domain = this.condRect.xmin_domain;
            r.ymin_domain = this.condRect.ymin_domain;
            const role = this.condRect.rectRoles && this.condRect.rectRoles[i];
            r.is_signal = role ? !!role.is_signal : true;
        });

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
        const refTol = this.opts.refineTol ?? 0.003;
        const refineFrac = this.opts.refineFrac ?? 0.15;   // matches the FDM backend's fraction
        const minRefineNodes = this.opts.minRefineNodes ?? 0;
        // Certify the tolerance against a controlled uniform-refinement solve before
        // accepting the (optimistic) pass-to-pass convergence gate. See _certifyStatic.
        const certify = this.opts.certify ?? true;
        this.certification = null;
        this._certWarn = null;
        // Require the convergence test to hold for this many CONSECUTIVE passes before
        // stopping (matches the FDM backend's min_converged_passes / UI "Min Converged
        // Passes"): guards against a premature stop on a single lucky pass.
        const minConvergedPasses = Math.max(1, this.opts.minConvergedPasses ?? 1);
        let convergedCount = 0;
        // Refine on ALL solved modes (differential: odd AND even) so each mode's
        // critical regions are resolved — crucially the inter-trace gap, where the
        // odd mode's field concentrates and the mutual resistance R[12] comes from.
        const refModes = this.modeNames;
        // Run the full-wave eigensolve EVERY refinement pass and converge on the actually
        // REPORTED quantities (the dispersive ε_eff from the eigenmode and the static Z0),
        // not on a cheaper static-only proxy. The static field-energy metric stabilises a
        // mesh or two before the full-wave eigenmode does, so a static-only convergence test
        // signs off while the reported ε_eff is still drifting (e.g. 4.305 → 4.416 on an open
        // asymmetric stripline). The per-pass eigensolve costs more, but an adaptive mesh
        // whose entire purpose is an accurate ε_eff/Z0/loss should converge on those numbers.
        const fRef = Math.max(s.freq || 1e9, 1e9);               // full-wave eval frequency
        let prev = null, prevEnergy = null;
        for (let it = 0; it <= maxIters; it++) {
            let metric = null;
            const conv = [];   // convergence quantities (eps + loss surface integral per mode)
            const energy = []; // total field energy per mode (for the reported energy error)
            const certQ0 = []; // this pass's (W_eps, W_air) per mode, for the certificate
            this._passEnergies = null;
            // Per-mode work of this pass (freedom map, statics, pencil decomposition,
            // eigenmode at fRef). When this turns out to be the final mesh,
            // _prepareStatic and the first solveAt would redo all of it on the same
            // objects, so it is handed over instead (see this._passWork).
            const passWork = { mesh, byMode: {} };
            // Source-free waveguide: there is no static solve, so the whole per-mode body
            // in the else branch is replaced, see _wgRefinePass.
            if (this._isWG) {
                const p = this._wgRefinePass(mesh);
                metric = p.metric; conv.push(...p.conv); energy.push(...p.energy);
            } else for (const rm of refModes) {
                const { abc } = modeConfig(rm, s.is_differential, this.symmetry, this.condRect.wallPEC);
                const pot = drivePotentials(this.condRect, rm, this.symmetry);
                const fm = buildTriFreedomMap(mesh, this.condRect, abc);
                const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, pot, this.ctx.wasmSolver);
                const phiAir = solveTriStatic(mesh, fm, null, pot, this.ctx.wasmSolver);
                const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
                const W_air = computeTriEnergy(phiAir, mesh, null);
                const eps_static = W_eps / W_air;
                energy.push(W_eps);
                // Same (mesh, mode, drive) pair _staticEnergies would recompute as the
                // certificate's level-0 solve, hand it over instead (see _certifyStatic).
                certQ0.push(W_eps, W_air);
                // Full-wave eigenmode at fRef: source of the reported dispersive ε_eff
                // (fw.eps) and the H-field for the ZZ refinement metric. Computed every pass.
                let fw = null;
                const femCache = {};
                try { fw = fullwaveMode(this.ctx, mesh, fm, abc, this.condRect, mesh.epsMap, fRef, phiEps, eps_static, femCache); } catch { fw = null; }
                passWork.byMode[rm] = { fm, phiEps, phiAir, W_eps, W_air, femCache, fw, fRef, epsMap: mesh.epsMap };
                // Converge on the REPORTED quantities, not a static proxy:
                //   • static characteristic impedance Z0 = sqrt(L/C) ∝ 1/sqrt(W_eps·W_air)
                //   • full-wave effective permittivity ε_eff = fw.eps (what _modeAtFreq reports)
                // The eigenmode ε_eff lags the static energy, so including it is what makes the
                // adaptive keep refining until the number the user sees has settled. Fall back
                // to the static C-ratio when the eigensolve fails on a pass, so the convergence
                // vector keeps a consistent length across passes.
                conv.push(1 / Math.sqrt(Math.max(W_eps * W_air, 1e-300)));
                conv.push(fw && fw.eps > 0 ? fw.eps : eps_static);
                // Static marking metric: Kelly flux-jump error indicator (refines where
                // the discretization error is), falling back to the legacy field-energy
                // metric if it degenerates (e.g. every edge constrained).
                let metricS = null;
                try { metricS = perElementKelly(phiEps, mesh, mesh.epsMap, fm); } catch { metricS = null; }
                if (!metricS || !metricS.some(v => v > 0)) metricS = perElementEnergy(phiEps, mesh, mesh.epsMap);
                let zz = null;
                if (fw) { try { zz = computeHtZZMetric(mesh, fm, fw.vRe, fw.vIm, fw.g2Re, fw.g2Im, fRef, this.condRect.rects, null, this.ctx.wasmSolver); } catch { zz = null; } }
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
            if (certQ0.length) this._passEnergies = { mesh, q: certQ0 };
            this._passWork = this._isWG ? null : passWork;
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
            // Stop before exceeding the triangle budget. Refining marks refineFrac of the
            // triangles and splits each ~1→4 (≈ +3·refineFrac growth), so project the next
            // size and stop now if it would blow the budget — keeps every SOLVED mesh within
            // the memory-parity budget rather than overshooting by a full pass.
            const projTris = Math.round(mesh.nTris * (1 + 3 * refineFrac));
            if (shouldStop && shouldStop()) {
                console.log('Adaptive refinement stopped by user');
                break;
            }
            // Verification certificate: the pass-to-pass gate above measures
            // the rate of approach, not the absolute error. Each pass refines
            // only refineFrac of the triangles, so with a per-pass error-decay
            // ratio r near 1 the remaining error is (observed change)*r/(1−r).
            // A tripped gate is a candidate stop. Certify the actual error
            // first (see _certifyStatic) and keep refining if the certificate
            // fails. Waveguides have no static solve to certify with, and
            // a mesh past certifyMaxTris (cert === null) keeps the legacy
            // behavior.
            if (convergedCount >= minConvergedPasses) {
                let cert = null;
                if (certify && !this._isWG) {
                    try { cert = this._certifyStatic(mesh, refTol); } catch { cert = null; }
                    if (cert) this.certification = cert;
                }
                if (!certify || this._isWG || !cert || cert.pass) break;
                convergedCount = 0;   // gate was optimistic — keep refining
            }
            if (it === maxIters || projTris > maxTris) break;
            const marked = markTrianglesForRefinement(metric, refineFrac);
            const refined = refineTriMesh(mesh, marked);
            refined.condRect = this.condRect;
            const mats = tagMaterials(refined, s.dielectrics);
            refined.epsMap = mats.epsMap; refined.lossMap = mats.lossMap;
            mesh = refined;
        }

        this.mesh = mesh;
        // Accuracy report: if refinement ended without a passing certificate
        // (iteration cap, triangle budget), measure the final mesh once so the
        // user gets an estimated error. A certificate that already covered this
        // mesh (pass or fail) is not recomputed.
        if (certify && !this._isWG) {
            let cert = this.certification;
            if (!cert || (!cert.pass && cert.tris !== mesh.nTris)) {
                try { cert = this._certifyStatic(mesh, refTol); this.certification = cert; }
                catch { /* keep whatever we had */ }
            }
            if (cert && !cert.pass) {
                // Three failure regimes need three phrasings: pre-asymptotic (estimate
                // is only a lower bound), estimate over the tolerance, and estimate
                // under the tolerance but without the safety margin the certificate
                // requires (saying "stopped before reaching 1%: error is 0.7%" would
                // read as a contradiction).
                const pct = (x) => (100 * x).toFixed(2);
                const estStr = cert.preAsymptotic
                    ? `at least ${pct(cert.err)}% (mesh still pre-asymptotic — the estimate is a lower bound)`
                    : cert.err < refTol
                        ? `about ${pct(cert.err)}%, within the tolerance but without the ` +
                          `×${cert.safety ?? 1.5} margin needed to certify it`
                        : `about ${pct(cert.err)}%`;
                // `reason` matches the rectilinear backend's certificate warning: machine
                // consumers (the fuzzer's R-gate relaxation, the UI's warning dedup) key
                // on it to tell a certificate note from a loss-accuracy note.
                this._certWarn = { type: 'accuracy', reason: 'certificate', mode: 'all', message:
                    `Full-wave mesh refinement could not verify the requested tolerance ` +
                    `(${pct(refTol)}%): estimated remaining error is ${estStr}. ` +
                    `Increase Max Nodes / Max Iterations, or relax Tolerance.` };
                globalThis.__TRI_DEBUG__ && console.warn('[tri] ' + this._certWarn.message);
            }
        }
        this.solver.certification = this.certification;
        // Mesh quality of the final (post-refinement) mesh. Q = circumradius/(2·inradius),
        // ideal 1; a very high Qmax means a sliver that ill-conditions the FEM solve.
        // Surfaced to the UI as a warning (empty constraint args → quality metrics only).
        // Conductor-interior triangles are left out: they carry no field DOFs (every
        // node and edge inside metal is PEC in the static and eigen solves) and the
        // skin band re-refines them before the MQS solve, so their shape says nothing
        // about this mesh's accuracy. A rect thinner than its surface element size
        // (a 1 um trace on the thick-copper sizing) holds flat slivers of Q in the
        // hundreds that would otherwise trip the UI warning. Shaped conductors are
        // never interior-meshed and their bbox can span the domain (coax shield), so
        // only rects are tested.
        try {
            const skip = new Uint8Array(mesh.nTris);
            const rects = (this.condRect.rects || []).filter(r => !r.shape);
            if (rects.length) for (let t = 0; t < mesh.nTris; t++) {
                const v0 = mesh.tris[3 * t], v1 = mesh.tris[3 * t + 1], v2 = mesh.tris[3 * t + 2];
                const xc = (mesh.nodes[2 * v0] + mesh.nodes[2 * v1] + mesh.nodes[2 * v2]) / 3;
                const yc = (mesh.nodes[2 * v0 + 1] + mesh.nodes[2 * v1 + 1] + mesh.nodes[2 * v2 + 1]) / 3;
                for (const r of rects) {
                    if (xc > r.xmin && xc < r.xmax && yc > r.ymin && yc < r.ymax) { skip[t] = 1; break; }
                }
            }
            const mq = checkMeshQuality(mesh, [], [], { skip });
            this.meshQuality = mq.metrics;
            if (this.solver) this.solver.meshQuality = mq.metrics;
        } catch { this.meshQuality = null; }
        // Reject a degenerate/collapsed mesh up front (the full-wave analogue of the
        // rectilinear validate_laplace_inputs guard) rather than producing NaN/Inf in the
        // eigensolve.
        const meshErrors = validateTriMesh(mesh, this.condRect);
        if (meshErrors.length)
            throw new Error('Full-wave mesh validation failed:\n' + meshErrors.map(e => ' - ' + e).join('\n'));
        if (this._isWG) this._prepareWaveguide(); else this._prepareStatic();
        return mesh;
    }

    // Verification certificate
    // Refine the mesh uniformly and nested (refineTriMeshNested: red split,
    // nodes fixed) and re-solve the statics:
    //   d1 = max rel change of (W_eps, W_air) per mode, mesh -> level 1
    //   certified error = d1/(1-r),  r = level-to-level error ratio
    // On a nested P2 hierarchy r is set by the strongest field singularity: h halves per
    // level and the energy error goes as h^(2λ), so a conductor corner (λ = 2/3 for its
    // 270° exterior angle) gives r = ~0.40 and a knife edge (λ = 1/2) is the worst case at
    // r = 0.5. R_BOUND = 0.5 is therefore the one-level estimate.
    // A second level, which measures r, costs 16x the base mesh
    // and only changes the verdict when it is borderline (d1 alone passes, 2*d1 does
    // not), so it runs only then.
    // Every reported static quantity (C0, eps_static, Z0) is a combination of W_eps and
    // W_air, and the eigensolve part of the reported eps_eff is anchored to the static
    // solve (_eigenBias), so the static error is the reported error. Gates:
    //   * d1 < noise floor  -> pass outright: below the linear-solver reproducibility.
    //   * d1*safety ≥ tol   -> already failed, no second level (cost control).
    //   * 2*d1*safety < tol -> pass on level 1 (rSource 'bound').
    //   * otherwise level 2 when it fits certifyMaxTris: r ≥ certifyRMax -> pre-asymptotic,
    //     cannot certify, else err = d1/(1−r). Over the cap -> fail on the level-1 bound.
    //     A base mesh whose level 1 is over the cap returns null.
    // The pass decision applies certifySafety (x1.5) on top of the estimate. `err` stays
    // the un-inflated best estimate (it is what warnings report).
    //
    // Result shape is the same as FieldSolver2D._certifyStatic so the UI and tests can
    // read either backend's `solver.certification` the same way:
    //   { pass, err, d1, r, rSource, levels, safety, preAsymptotic?, <size> }
    // rSource is 'bound' (r = R_BOUND), 'measured' (level-2 ratio) or null (no r, the
    // d1-alone failure branch). The size field is deliberately named in each backend's
    // own unit (tris or nodes).
    _nestedRefine(mesh) {
        const refined = refineTriMeshNested(mesh);
        refined.condRect = this.condRect;
        // Inherit the parent's material map (children 4t..4t+3 of parent t): the
        // certificate must solve the same problem on nested spaces.
        const epsMap = new Array(refined.nTris), lossMap = new Array(refined.nTris);
        for (let t = 0; t < mesh.nTris; t++) {
            for (let k = 0; k < 4; k++) { epsMap[4 * t + k] = mesh.epsMap[t]; lossMap[4 * t + k] = mesh.lossMap[t]; }
        }
        refined.epsMap = epsMap; refined.lossMap = lossMap;
        return refined;
    }

    _staticEnergies(mesh) {
        const s = this.solver;
        const out = [];
        for (const rm of this.modeNames) {
            const { abc } = modeConfig(rm, s.is_differential, this.symmetry, this.condRect.wallPEC);
            const pot = drivePotentials(this.condRect, rm, this.symmetry);
            const fm = buildTriFreedomMap(mesh, this.condRect, abc);
            const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, pot, this.ctx.wasmSolver);
            const phiAir = solveTriStatic(mesh, fm, null, pot, this.ctx.wasmSolver);
            out.push(computeTriEnergy(phiEps, mesh, mesh.epsMap), computeTriEnergy(phiAir, mesh, null));
        }
        return out;
    }

    _certifyStatic(mesh, tol) {
        const maxDiff = (a, b) => {
            let m = 0;
            for (let i = 0; i < a.length; i++)
                m = Math.max(m, Math.abs(a[i] - b[i]) / Math.max(Math.abs(b[i]), 1e-300));
            return m;
        };
        const GROW = 4;                                   // red split: 4 children per triangle
        const cap = this.opts.certifyMaxTris ?? 150000;
        if (mesh.nTris * GROW > cap) return null;
        const rMax = this.opts.certifyRMax ?? 0.7;
        const safety = this.opts.certifySafety ?? 1.5;
        const R_BOUND = 0.5;                              // knife-edge rate, worst singularity
        // Level 0 is exactly what the refinement pass on this mesh already solved
        // (same freedom map, same drives), reuse it rather than re-solving.
        const pe = this._passEnergies;
        const q0 = (pe && pe.mesh === mesh) ? pe.q : this._staticEnergies(mesh);
        const m1 = this._nestedRefine(mesh);
        const q1 = this._staticEnergies(m1);
        const d1 = maxDiff(q0, q1);
        const base = { tris: mesh.nTris, d1, safety };
        if (d1 < 5e-5) return { ...base, pass: true, err: d1, r: 0, rSource: null, levels: 1 };
        if (d1 * safety >= tol)
            return { ...base, pass: false, err: d1, r: null, rSource: null, levels: 1 };
        const errBound = d1 / (1 - R_BOUND);
        if (errBound * safety < tol)
            return { ...base, pass: true, err: errBound, r: R_BOUND, rSource: 'bound', levels: 1 };
        // Borderline: only a measured ratio below the bound can still certify this mesh.
        if (m1.nTris * GROW > cap)
            return { ...base, pass: false, err: errBound, r: R_BOUND, rSource: 'bound', levels: 1 };
        const m2 = this._nestedRefine(m1);
        const q2 = this._staticEnergies(m2);
        const r = maxDiff(q1, q2) / d1;
        if (r >= rMax)
            return { ...base, pass: false, err: d1, r, rSource: 'measured', levels: 2, preAsymptotic: true };
        const err = d1 / (1 - Math.min(r, rMax));
        return { ...base, pass: err * safety < tol, err, r, rSource: 'measured', levels: 2 };
    }

    // Rectangular waveguide path
    //
    // Reference frequency for the one eigensolve this medium runs: 1.5x cutoff.
    //
    // kc² = γ² + k0²·εr is a cancellation of two large numbers, so the relative error in
    // γ² is amplified by (f/fc)² into kc², x100 at 10*fc, only x2.25 at 1.5*fc. Pinning
    // it to a multiple of fc rather than the user's frequency also makes the refined mesh
    // independent of the sweep range.
    wgRefFreq() { return 1.5 * this.solver.fc; }

    // One adaptive-refinement pass for the waveguide. Replaces the static energy /
    // quasi-TEM eps_eff pair with the two quantities that actually have an analytic answer
    // here, and refines on the projected-Ht ZZ metric alone (there is no static field to
    // blend in, and the walls are where the loss integral lives anyway).
    _wgRefinePass(mesh) {
        const s = this.solver;
        const f = this.wgRefFreq();
        const { abc } = modeConfig('single', false, this.symmetry, this.condRect.wallPEC);
        const fm = buildTriFreedomMap(mesh, this.condRect, abc);
        const seed = analyticSeedDofs(mesh, fm, wgFundamentalEfun(this.condRect));
        let wg = null;
        try {
            wg = waveguideEigen(this.ctx, mesh, fm, abc, this.condRect, mesh.epsMap, f,
                seed, s.kc_analytic);
        } catch { wg = null; }
        // Keep the convergence vector a consistent length across passes when a pass fails.
        if (!wg) return { metric: null, conv: [s.kc_analytic, 0], energy: [s.kc_analytic] };

        // Converge on kc and on alpha_c. The two quantities with an exact analytic answer
        // here, both evaluated at the same fRef every pass so this is a relative test.
        //
        // kc is the binding one. alpha_c is a pure SURFACE integral, which is what makes it
        // converge last on every OTHER medium (see CoaxSolver._hFine), but
        // a waveguide has straight walls and a smooth sinusoidal mode, no
        // corner singularities or polygonized curve, so measured it is
        // converged to five figures on the coarsest usable mesh while kc is
        // still moving. It stays in the vector anyway: it is cheap next to the
        // eigensolve, and it is the only thing that would notice a mesh
        // resolving the bulk field but not the wall currents.
        const beta = Math.sqrt(Math.max(-wg.g2Re, 0));
        const { kappa } = wgGeom(this.condRect);
        const Zpv = beta > 0 ? kappa * 2 * Math.PI * f * MU0 / beta : 1;
        const wgState = { fm, vRe: wg.vRe, vIm: wg.vIm, projCache: {},
            lossMask: buildLossEdges(mesh, fm, this.condRect) };
        const conv = [wg.kc, this._wgConductorLoss(mesh, wgState, f, beta, Zpv).alpha_c_np];

        let zz = null;
        try {
            zz = computeHtZZMetric(mesh, fm, wg.vRe, wg.vIm, wg.g2Re, wg.g2Im, f,
                this.condRect.rects, null, this.ctx.wasmSolver);
        } catch { zz = null; }
        let metric = null;
        if (zz) {
            let mZ = 0;
            for (let i = 0; i < zz.length; i++) if (zz[i] > mZ) mZ = zz[i];
            metric = new Float64Array(zz.length);
            if (mZ > 0) for (let i = 0; i < zz.length; i++) metric[i] = zz[i] / mZ;
        }
        return { metric, conv, energy: [wg.kc] };
    }

    // Waveguide substitute for _prepareStatic. What the rest of the backend actually needs
    // from that method is a freedom map, a resampled field for plotting, and something to
    // centre the eigen shift on; all three come from ONE eigensolve of the fundamental.
    //
    // The eigenvector is then reused at every sweep frequency, exactly rather than
    // approximately: kc is geometric, and projectH builds Ht from (Et, ∇Ez, γ) with
    // Ez ≡ 0 for a TE mode, so the only frequency dependence is through γ and ω, both
    // passed as arguments. The eigenvector's arbitrary scale cancels in h2dl/P.
    _prepareWaveguide() {
        const s = this.solver, mesh = this.mesh;
        const { abc } = modeConfig('single', false, this.symmetry, this.condRect.wallPEC);
        const fm = buildTriFreedomMap(mesh, this.condRect, abc);
        const seed = analyticSeedDofs(mesh, fm, wgFundamentalEfun(this.condRect));
        const wg = waveguideEigen(this.ctx, mesh, fm, abc, this.condRect, mesh.epsMap,
            this.wgRefFreq(), seed, s.kc_analytic);
        if (!wg) throw new Error(
            'Waveguide eigensolve failed: no physical mode converged at the reference ' +
            'frequency. Check the guide dimensions and the Max Nodes budget.');

        // Sanity gate against the closed form. Outside it the mesh, the boundary conditions
        // or the pick went wrong, and the analytic kc is the honest fallback, unlike the
        // quasi-static eps the quasi-TEM path falls back to, which is meaningless here.
        const gated = Math.abs(wg.kc / s.kc_analytic - 1) > 0.02;
        const kc = gated ? s.kc_analytic : wg.kc;

        this._staticGrid = buildGridFromMesh(mesh, this.domain, { resolution: this.opts.resolution });
        this._plotRects = [];   // no conductors to clamp the resampler's E baseline against
        const fields = resampleModeField(mesh, fm, wg.vRe, wg.vIm, this.domain,
            { resolution: this.opts.resolution, grid: this._staticGrid });
        this._wg = { fm, abc, seed, kc, kcRaw: wg.kc, gated,
            vRe: wg.vRe, vIm: wg.vIm, g2Re: wg.g2Re, g2Im: wg.g2Im,
            fields, projCache: {}, lossMask: buildLossEdges(mesh, fm, this.condRect) };
        // solveModes and solveAt's plotting writeback both index _static by mode name;
        // populating it here keeps them branch-free. modeNames is ['single'] (a waveguide
        // is never differential), and `fields` is the SAME object as _wg.fields.
        this._static = { single: { fm, abc, kC: 2, phiEps: null, C0: null, W_loss: null,
            eps_eff_static: null, Z_static: null, fields } };
    }

    // Conductor loss of a waveguide mode at f. Returns { alpha_c_np, R_ac }.
    //
    // solveConductorLoss needs NO waveguide special-casing: with an empty rect list
    // signal_area is 0, so R_dc = 0, R_combined = R_ac, and alpha_c = R_ac/(2·Z0) is
    // exactly the Z0-independent alpha_c_ac = Rs*h2dl/(4*P), the Z0 passed in cancels.
    // Passing Zpv makes the R_ac it returns the series resistance this medium's equivalent
    // circuit wants, so one call yields both numbers.
    //
    // Roughness and plating are applied HERE rather than through buildSurfaceGroups: that
    // path resolves a surface impedance per (rect, face), and a waveguide has no rects at
    // all, so every wall edge would silently fall through to bare metal. The wall is ONE
    // continuous surface, so a single scalar Zs.re/Rs scaling is exact, the same
    // convention the per-face groups use on the quasi-TEM path.
    _wgConductorLoss(mesh, wg, f, beta, Zpv) {
        const s = this.solver;
        const sigma = s.sigma_cond ?? 5.8e7;
        const omega = 2 * Math.PI * f;
        if (!(beta > 0)) return { alpha_c_np: 0, R_ac: 0, X_ac: 0 };
        let projH, P;
        try {
            projH = projectH(mesh, wg.fm, wg.vRe, wg.vIm, { re: 0, im: beta }, f,
                this.ctx.wasmSolver, wg.projCache);
            P = Math.abs(computePoyntingFromProjectedH(mesh, wg.fm, wg.vRe, wg.vIm,
                projH.htRe, projH.htIm, omega * MU0, projH.hDofs));
        } catch { return { alpha_c_np: 0, R_ac: 0, X_ac: 0 }; }
        if (!(P > 1e-30)) return { alpha_c_np: 0, R_ac: 0, X_ac: 0 };
        const loss = solveConductorLoss(this.condRect.rects, f, sigma, mesh, wg.fm,
            wg.vRe, wg.vIm, -beta * beta, 0, P, Zpv, mesh.epsMap, wg.lossMask,
            projH, this.ctx.wasmSolver);
        // Surface finish: smooth-metal Rs is what solveConductorLoss assumed.
        const delta = Math.sqrt(2 / (omega * MU0 * sigma));
        const Rs = 1 / (sigma * delta);
        const pl = s.plating;
        const Zs = (pl && pl.sigma > 0)
            ? calculate_Zrough_layered(f, sigma, pl.rq ?? s.rq ?? 0, pl.sigma, pl.thickness ?? 0)
            : calculate_Zrough(f, sigma, s.rq ?? 0);
        // Re(Zs) carries the loss and Im(Zs) the internal inductance; solveConductorLoss
        // evaluated BOTH at the smooth Rs, so each is recovered by its own scaling of the
        // same integral. Smooth metal has Zs = Rs(1+j) => fX = fR, which is why a bare
        // guide is unchanged by carrying X_ac.
        const fR = Rs > 0 ? Zs.re / Rs : 1;
        const fX = Rs > 0 ? Zs.im / Rs : 1;
        return { alpha_c_np: loss.alpha_c * fR, R_ac: loss.R_ac * fR, X_ac: loss.R_ac * fX };
    }

    // Per-frequency result for the waveguide fundamental. Everything here is analytic in f
    // on top of the single cached eigensolve: kc is geometric, so beta, the impedances and
    // alpha_d follow in closed form, and alpha_c reuses the cached eigenvector.
    //
    // The reported impedance is the POWER-VOLTAGE definition Zpv = kappa*Z_TE, which for
    // TE10 follows exactly from V = E0*b and P = E0²ab/(4*Z_TE). Waveguide impedance is not
    // unique, the choice is a normalization, and every normalization-independent quantity
    // here (gamma, eps_eff, alpha) is unaffected by it.
    _waveguideModeAtFreq(f) {
        const s = this.solver, mesh = this.mesh, wg = this._wg;
        const kc = wg.kc;
        const omega = 2 * Math.PI * f;
        const k0 = omega / c0;
        // Read the fill from the MESH so causal materials are tracked automatically.
        const er = maxEpsRe(mesh.epsMap);
        const tand = wgFillTand(mesh);
        const { kappa } = wgGeom(this.condRect);
        const fc = c0 * kc / (2 * Math.PI * Math.sqrt(er));
        const fc2 = s.fc2;
        const warn = (type, message) => {
            if (this._modeWarnings && !this._modeWarnings.some(w => w.type === type))
                this._modeWarnings.push({ type, mode: 'single', freq: f, message });
        };
        if (wg.gated) {
            warn('wg-kc-gate',
                `Waveguide cutoff from the field solve (${(wg.kcRaw / (2 * Math.PI) * c0 / 1e9).toFixed(3)} GHz) ` +
                `disagrees with the closed form by more than 2%, the mesh is too coarse. ` +
                `Falling back to the analytic cutoff. Increase Max Nodes.`);
        }
        if (fc2 && f >= fc2) {
            warn('wg-overmoded',
                `${(f / 1e9).toFixed(3)} GHz is at or above the ${(fc2 / 1e9).toFixed(3)} GHz second ` +
                `cutoff, the guide is over-moded. Only the fundamental mode is reported. ` +
                `Use the Modes tab to inspect the higher-order modes.`);
        }

        const g2 = kc * kc - k0 * k0 * er;
        if (g2 >= 0) {
            // Below cutoff. Only the attenuation is meaningful. Everything else is reported
            // as NaN so the Results traces break cleanly at cutoff instead of autoscaling
            // around the (correct, but large and negative) evanescent shunt capacitance.
            warn('wg-cutoff',
                `${(f / 1e9).toFixed(3)} GHz is below the ${(fc / 1e9).toFixed(3)} GHz cutoff. ` +
                `The mode is evanescent, so only the attenuation is reported.`);
            return {
                mode: 'single', Z0: NaN, eps_eff: NaN, eps_eff_mode: NaN, C: NaN, C0: NaN,
                RLGC: { R: NaN, L: NaN, G: NaN, C: NaN },
                Zc: new Complex(NaN, NaN),
                alpha_c: NaN, alpha_d: NaN, alpha_total: NP_TO_DB * Math.sqrt(g2),
                L_internal: NaN, L_external: NaN,
                kc, fc, beta: 0, Z_TE: NaN, Zpv: NaN, belowCutoff: true,
                self_referenced: true,
                V: null, Ex: wg.fields.Ex, Ey: wg.fields.Ey,
            };
        }

        const beta = Math.sqrt(-g2);
        const Z_TE = omega * MU0 / beta;          // = k0·eta0/beta, independent of er
        const Zpv = kappa * Z_TE;
        const { alpha_c_np, R_ac, X_ac } = this._wgConductorLoss(mesh, wg, f, beta, Zpv);
        const alpha_d_np = k0 * k0 * er * tand / (2 * beta);

        // Exact per-unit-length equivalent circuit of the mode, normalized to Zpv:
        //   series Z = j*omega*mu0*kappa  (+ the wall surface resistance)
        //   shunt  Y = (omega*eps0*er*tand + j*omega*eps0*er − j*kc²/(omega*mu0)) / kappa
        // It reproduces gamma and Zc exactly, and eps_eff = c0²*L*C = er - (kc/k0)² falls
        // out with kappa cancelling. The R/G split is normalization-dependent, gamma is not.
        const L_external = kappa * MU0;
        // The wall's internal inductance is its surface REACTANCE over omega, not its
        // resistance: those coincide only while Zs = Rs(1+j), which a rough or plated
        // wall breaks (Im(Zs)/Re(Zs) reaches ~5 at 1 um rms). Bare walls are unchanged.
        const L_internal = omega > 0 ? Math.min(X_ac / omega, 0.5 * L_external) : 0;
        const L = L_external + L_internal;
        const C = (eps0 * er - kc * kc / (omega * omega * MU0)) / kappa;
        const G = omega * eps0 * er * tand / kappa;
        const R = R_ac;
        const Zc = new Complex(R, omega * L).div(new Complex(G, omega * C)).sqrt();
        return {
            mode: 'single',
            Z0: Zpv,
            eps_eff: c0 * c0 * L * C, eps_eff_mode: er - (kc / k0) ** 2,
            C, C0: eps0,
            RLGC: { R, L, G, C },
            Zc,
            alpha_c: NP_TO_DB * alpha_c_np, alpha_d: NP_TO_DB * alpha_d_np,
            alpha_total: NP_TO_DB * (alpha_c_np + alpha_d_np),
            L_internal, L_external,
            kc, fc, beta, Z_TE, Zpv, belowCutoff: false,
            // There is no external TEM reference for a guide: the S-parameter paths
            // normalize each point to Zc, giving S11 = 0 and S21 = exp(-gamma*L).
            self_referenced: true,
            V: null, Ex: wg.fields.Ex, Ey: wg.fields.Ey,
        };
    }

    _prepareStatic() {
        const { mesh } = this;
        const s = this.solver;
        this._static = {};
        // Resample fields onto a graded grid derived from the triangle mesh itself:
        // grid-line density follows the adaptively refined mesh (fine at corners,
        // surfaces and gaps — exactly where the solution needed resolution), with the
        // conductor/dielectric interface lines forced in so sub-grid-thickness features
        // keep crisp plateaus instead of aliasing into spurious contour bands. This
        // replaces the FDM mesher's geometry-heuristic grid, which knew nothing about
        // where the field actually concentrates.
        const forcedX = [], forcedY = [];
        for (const r of [...(s.conductors || []), ...(s.dielectrics || [])]) {
            // A shaped body has no axis-aligned interfaces to force lines onto, and a
            // complement's bbox lies outside the meshed domain entirely.
            if (r.shape) continue;
            forcedX.push(r.x_min, r.x_max);
            forcedY.push(r.y_min, r.y_max);
        }
        // Near-surface companion lines (the FDM mesher's "boundary line" trick): one
        // grid line a conductor-dimension/20 outside each face. The row ON a face
        // carries the one-sided surface gradient; without a nearby second row the
        // segment down to it spans the whole first (coarse-mesh) cell and any surface
        // bias shows as a visible kink where contours meet the ground. With the
        // companion line the segment is ~µm-scale and the join renders clean.
        for (const r of (s.conductors || [])) {
            if (r.shape) continue;   // no flat faces to place a companion line beside
            const off = Math.min(r.x_max - r.x_min, r.y_max - r.y_min) / 20;
            if (!(off > 0)) continue;
            forcedX.push(r.x_min - off, r.x_max + off);
            forcedY.push(r.y_min - off, r.y_max + off);
        }
        const grid = buildGridFromMesh(mesh, this.domain, {
            resolution: this.opts.resolution, forcedX, forcedY, mirrorX: this.symmetry,
        });
        // Cached for the per-frequency causal-materials rebuild (_applyCausal), which re-runs
        // the same static preparation under the updated permittivity.
        this._staticGrid = grid;
        // Conductor rects in full-domain coordinates for the resampler's E-baseline
        // clamping — taken from the solver geometry because it includes the ground
        // planes (mesh.condRect.rects only carries the non-wall conductors).
        // A complement conductor (coax shield) is EXCLUDED: its bounding box spans the
        // whole domain, so the resampler would treat every grid point as inside metal
        // and clamp |E| to zero across the entire plot. Its surface is the domain
        // outline, which the resampler already handles as a boundary.
        this._plotRects = (s.conductors || [])
            .filter(c => !(c.shape && isComplement(c.shape)))
            .map(c => ({ xmin: c.x_min, xmax: c.x_max, ymin: c.y_min, ymax: c.y_max, shape: c.shape || null }));
        // Asymmetric differential pair on the full domain: the odd/even excitation basis is
        // NOT the modal basis (the two traces sit in different environments), so driving
        // odd/even and picking one eigenmode per drive can return the SAME eigenmode twice
        // (e.g. broadside stripline with an unequal er stack). Use the genuine line modes
        // instead — the eigenvectors of [C0]⁻¹[C] from per-trace drives.
        // The last refinement pass solved exactly these statics (same mesh, freedom
        // map, drives and material map) and assembled the pencil for its eigensolve
        // at fRef; reuse them rather than recompute (this._passWork, buildMesh). The
        // modal path below derives its own basis and does not use them.
        const pw = (this._passWork && this._passWork.mesh === mesh) ? this._passWork.byMode : null;
        this._passWork = null;
        if (s.is_differential && !this.symmetry &&
            this._prepareStaticModal(grid, s.use_causal_materials ? 'defer' : 'build')) return;
        for (const mode of this.modeNames) {
            const { abc, kC } = modeConfig(mode, s.is_differential, this.symmetry, this.condRect.wallPEC);
            const w = (pw && pw[mode] && pw[mode].epsMap === mesh.epsMap) ? pw[mode] : null;
            const fm = w ? w.fm : buildTriFreedomMap(mesh, this.condRect, abc);
            const pot = drivePotentials(this.condRect, mode, this.symmetry);
            const phiEps = w ? w.phiEps : solveTriStatic(mesh, fm, mesh.epsMap, pot, this.ctx.wasmSolver);
            const phiAir = w ? w.phiAir : solveTriStatic(mesh, fm, null, pot, this.ctx.wasmSolver);
            const W_eps = w ? w.W_eps : computeTriEnergy(phiEps, mesh, mesh.epsMap);
            const W_air = w ? w.W_air : computeTriEnergy(phiAir, mesh, null);
            const W_loss = computeTriEnergy(phiEps, mesh, mesh.lossMap);
            const C0 = kC * eps0 * W_air;                 // vacuum cap (geometric)
            const eps_eff_static = W_eps / W_air;
            const parity = this.symmetry ? (mode === 'odd' ? 'odd' : 'even') : null;
            // With causal materials the first solveAt immediately re-solves the statics
            // under eps(f) and resamples again, so resampling here would be thrown away.
            // _applyCausal fills st.fields in whenever it is still missing.
            const fields = s.use_causal_materials ? null : resampleStatic(mesh, phiEps, this.domain,
                { resolution: this.opts.resolution, parity, grid, rects: this._plotRects });
            this._static[mode] = {
                fm, abc, kC, phiEps, C0, W_loss, eps_eff_static,
                Z_static: 1 / (c0 * Math.sqrt(C0 * eps_eff_static * C0)),
                fields,
                // Closed-pick pencil decomposition from the pass (fullwaveMode's own
                // identity-validated cache) and its eigenmode at fRef, which the first
                // solveAt asks for again (_eigenPick).
                femCache: w ? w.femCache : undefined,
                fwSeed: (w && w.fw) ? { f: w.fRef, fw: w.fw, epsMap: w.epsMap, phiEps: w.phiEps } : undefined,
            };
        }
    }

    // General two-conductor modal analysis for an asymmetric differential pair (full domain).
    // Drives each trace independently to assemble the full 2×2 capacitance matrices [C]
    // (dielectric) and [C0] (vacuum), diagonalises [C0]⁻¹[C] to get the two genuine line
    // modes (eps_eff = eigenvalues, modal voltage ratios = eigenvectors), and builds each
    // mode's static field as the eigenvector combination of the per-trace fields. The two
    // modes are labelled odd/even by differential character (opposite-sign components → odd)
    // so the downstream result assembly is unchanged. For a symmetric pair the eigenvectors
    // come out as [1,∓1], reproducing the existing odd/even basis exactly.
    // fieldsMode: 'build' resamples the plot fields, 'defer' leaves them null (the
    // causal re-solve in _applyCausal supersedes them immediately), 'reuse' keeps the
    // ones already there (mid-sweep: nothing displays them, see _applyCausal).
    _prepareStaticModal(grid, fieldsMode = 'build') {
        const { mesh } = this;
        const s = this.solver;
        const roles = this.condRect.rectRoles;
        const { abc, kC } = modeConfig('odd', true, false, this.condRect.wallPEC);
        const fm = buildTriFreedomMap(mesh, this.condRect, abc);
        // Per-trace drives (positive- and negative-polarity signal conductors).
        const potA = roles.map(r => (r.is_signal && (r.polarity || 1) > 0) ? 1 : 0);
        const potB = roles.map(r => (r.is_signal && (r.polarity || 1) < 0) ? 1 : 0);
        // Both polarities must actually be present in the meshed domain — a missing
        // group (both traces tagged +, or one clipped out) makes phiB ≡ 0 and C0m
        // singular, producing a silent NaN cascade in the classifier. Fall back to
        // the odd/even drive instead.
        if (!potA.some(v => v) || !potB.some(v => v)) { this._modalPhys = null; return false; }
        const phiA = solveTriStatic(mesh, fm, mesh.epsMap, potA, this.ctx.wasmSolver);
        const phiB = solveTriStatic(mesh, fm, mesh.epsMap, potB, this.ctx.wasmSolver);
        const phiAa = solveTriStatic(mesh, fm, null, potA, this.ctx.wasmSolver);
        const phiBa = solveTriStatic(mesh, fm, null, potB, this.ctx.wasmSolver);
        // 2×2 capacitance matrices via the energy form W = ½·vᵀ·C·v.
        const We = X => computeTriEnergy(X, mesh, mesh.epsMap);
        const Wa = X => computeTriEnergy(X, mesh, null);
        const eaa = We(phiA), ebb = We(phiB), eab = We(combineStatic(1, phiA, 1, phiB));
        const aaa = Wa(phiAa), abb = Wa(phiBa), aab = Wa(combineStatic(1, phiAa, 1, phiBa));
        const C = [[2 * eaa, eab - eaa - ebb], [eab - eaa - ebb, 2 * ebb]];
        const C0m = [[2 * aaa, aab - aaa - abb], [aab - aaa - abb, 2 * abb]];
        // Shared symmetric/degenerate/modal decision (thresholds, ordering — see
        // classifyModalDecomposition; the rectilinear backend uses the identical guard).
        // The energy-form matrices are scaled to physical F/m by kC·ε0 first (same factor
        // as the per-mode C0 = kC·ε0·W_air), so physMatrix comes out in SI units.
        // Symmetric: return false and let _prepareStatic do the odd/even prep, with no
        // physMatrix. Degenerate: physMatrix without Tv still drives the asymmetric
        // 4-port S-matrix while _prepareStatic builds the odd/even per-mode static.
        const ksc = kC * eps0;
        const sc2 = (M, s) => [[M[0][0] * s, M[0][1] * s], [M[1][0] * s, M[1][1] * s]];
        const { physMatrix, modalVecs } = classifyModalDecomposition(
            sc2(C, ksc), sc2(C0m, ksc), s.conductors, s.dielectrics);
        this._modalPhys = physMatrix;
        if (!modalVecs) return false;
        ['odd', 'even'].forEach((label, li) => {
            const v = modalVecs[li];
            const phiEps = combineStatic(v[0], phiA, v[1], phiB);
            const phiAir = combineStatic(v[0], phiAa, v[1], phiBa);
            const W_air = Wa(phiAir);
            const C0 = kC * eps0 * W_air;
            const eps_eff_static = We(phiEps) / W_air;
            const prevFields = this._static[label] && this._static[label].fields;
            const fields = fieldsMode === 'defer' ? null
                : (fieldsMode === 'reuse' && prevFields) ? prevFields
                : resampleStatic(mesh, phiEps, this.domain,
                    { resolution: this.opts.resolution, parity: null, grid, rects: this._plotRects });
            this._static[label] = {
                fm, abc, kC, phiEps, C0, W_loss: computeTriEnergy(phiEps, mesh, mesh.lossMap),
                eps_eff_static,
                Z_static: 1 / (c0 * Math.sqrt(C0 * eps_eff_static * C0)),
                fields, modalVec: v,
            };
        });
        return true;
    }

    // Target net currents for the full-domain (per-conductor-drive) MQS solve,
    // ordered [positive-polarity group, negative-polarity group] to match
    // mqsPrecompute's group order. Symmetric pair: the odd/even unit vectors.
    // Asymmetric pair with a modal decomposition (_prepareStaticModal): the
    // genuine modal current vector i \propto C*v (v = the modal voltage eigenvector,
    // C = the SI capacitance matrix, γI = YV makes i an eigenvector of YZ),
    // normalized to Σ|I|^2 = 2 so the per-line power convention matches +-1.
    _mqsModeCurrents(mode, st) {
        const v = st.modalVec, Cm = this._modalPhys && this._modalPhys.C;
        if (v && Cm) {
            const i0 = Cm[0][0] * v[0] + Cm[0][1] * v[1];
            const i1 = Cm[1][0] * v[0] + Cm[1][1] * v[1];
            const n = Math.sqrt(2 / (i0 * i0 + i1 * i1));
            if (isFinite(n) && n > 0) return [i0 * n, i1 * n];
        }
        return mode === 'odd' ? [1, -1] : [1, 1];
    }

    // One exact quasi-TEM eigensolve at f for a prepared mode: the closed-wall pick
    // first, escalating to the radiating-ABC pick when that fails or is ambiguous (the
    // two-stage pick is explained in _modeAtFreq). Returns { fw, fwErr }.
    _eigenPick(st, f, phiEps, eps_eff_static) {
        const { mesh } = this, cr = this.condRect, fm = st.fm;
        let fw = null, fwErr = null;
        // The closed pick at this exact frequency, material map and static field may
        // already exist from the final refinement pass (st.fwSeed, see _prepareStatic).
        const seed = st.fwSeed;
        if (seed && seed.f === f && seed.epsMap === mesh.epsMap && seed.phiEps === phiEps) {
            fw = seed.fw;
        } else {
            try {
                fw = fullwaveMode(this.ctx, mesh, fm, st.abc, cr, mesh.epsMap, f, phiEps, eps_eff_static,
                    st.femCache || (st.femCache = {}));
            } catch (e) { fw = null; fwErr = e; }
        }
        if (!fw || fw.ambiguous) {
            const pickAbc = {};
            let radiates = false;
            for (const k of ['left', 'right', 'top', 'bottom']) {
                const v = st.abc[k];
                if (v === undefined) continue;            // PEC ground wall → leave absent (Dirichlet)
                const rad = !(this.symmetry && k === 'left');
                pickAbc[k] = rad ? true : 'pmc';
                if (rad) radiates = true;
            }
            if (radiates) {
                try {
                    const fw2 = fullwaveMode(this.ctx, mesh, fm, pickAbc, cr, mesh.epsMap, f,
                        phiEps, eps_eff_static, st.femCacheAbc || (st.femCacheAbc = {}));
                    if (fw2) { fw = fw2; fwErr = null; }
                } catch (e) { if (!fw) fwErr = e; }
            }
        }
        return { fw, fwErr };
    }

    // Static anchoring of the eigensolve. On a given mesh the quasi-TEM
    // eigenvalue sits above the variational static solve of the same mesh by
    // a discretization bias that is flat in frequency below a few GHz, while
    // the dispersion ratio eps_fw(f)/eps_fw(F_STATIC_MAX) is mesh-independent
    // to 1e-4. So the reported eps_d is
    //     eps_static * eps_fw(f) / eps_fw(F_STATIC_MAX)
    // the static energies (which the verification certificate covers) set the level and
    // the eigensolve contributes only the dispersion. It also makes eps_eff continuous at
    // the static/full-wave switch. One extra eigensolve per (mesh, mode), cached on the
    // backend keyed by mesh identity (survives the per-frequency causal re-solves of
    // this._static[*], whose material changes move the bias only in the noise). When no
    // anchor frequency yields a usable solve, or the candidate is implausible (|bias - 1|
    // > 5%: a different mode was picked), the raw eigenvalue stays in place and a warning
    // is surfaced once. `fwAtAnchor` lets the caller hand over a solve it already did at
    // exactly F_STATIC_MAX.
    _eigenBias(mode, st, phiEps, eps_eff_static, fwAtAnchor = null) {
        if (!this._eigenBiasCache || this._eigenBiasCache.mesh !== this.mesh)
            this._eigenBiasCache = { mesh: this.mesh, byMode: {} };
        const cache = this._eigenBiasCache.byMode;
        if (cache[mode] !== undefined) return cache[mode];
        // Safety net: walk up to 1 GHz if the anchor solve returns no quasi-TEM candidate.
        // This shouldn't be needed with tree-cotree cauge.
        const ANCHORS = [F_STATIC_MAX, 2e8, 5e8, 1e9];
        let bias = 1, ok = false;
        const reasons = [];
        const accept = (fw, fwErr, fa) => {
            const label = `${(fa / 1e6).toFixed(0)} MHz`;
            if (!fw) { reasons.push(`${label}: ${fwErr ? `eigensolver error: ${fwErr.message || fwErr}` : 'no converged quasi-TEM candidate'}`); return false; }
            if (!(fw.eps > 0 && eps_eff_static > 0)) { reasons.push(`${label}: non-physical eigenvalue`); return false; }
            const b = fw.eps / eps_eff_static;
            if (Math.abs(b - 1) > 0.05) {
                reasons.push(`${label}: ε_eff=${fw.eps.toFixed(3)} is ${((b - 1) * 100).toFixed(1)}% off the static ε_eff=${eps_eff_static.toFixed(3)}, not a discretization bias`);
                return false;
            }
            bias = b; return true;
        };
        if (fwAtAnchor) ok = accept(fwAtAnchor, null, F_STATIC_MAX);
        for (const fa of (fwAtAnchor ? ANCHORS.slice(1) : ANCHORS)) {
            if (ok) break;
            const { fw, fwErr } = this._eigenPick(st, fa, phiEps, eps_eff_static);
            ok = accept(fw, fwErr, fa);
        }
        if (!ok && this._modeWarnings && !this._modeWarnings.some(w => w.mode === mode && w.type === 'eigen-anchor')) {
            const msg = `${mode} mode: no usable full-wave anchor eigensolve (${reasons.join('; ')}) — reported `
                + `ε_eff keeps the raw eigenvalue, which can sit up to ~1% above the static ε_eff on a coarse mesh.`;
            this._modeWarnings.push({ type: 'eigen-anchor', mode, freq: F_STATIC_MAX, message: msg });
            globalThis.__TRI_DEBUG__ && console.warn('[tri full-wave] ' + msg);
        }
        cache[mode] = bias;
        return bias;
    }

    // Per-mode result at frequency f, reusing the cached static solve.
    // eps_eff comes from the full-wave eigenmode (f ≥ 100 MHz), anchored to the
    // static solve (_eigenBias), or the static solve alone (below), conductor loss is
    // from the robust static-field method.
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
        //   'auto'         — MQS volume eddy-current where applicable (symmetric,
        //                    rectangular conductors), else the blended perturbation.
        //                    Most accurate; default ("Full-wave (MQS)").
        //   'perturbation' — blended perturbation (static-field in the skin transition,
        //                    SIBC eigenmode at deep skin); much cheaper than MQS per
        //                    point, smooth across frequency. ("Full-wave (perturbation)").
        //   'static'       — static-field perturbation only (no SIBC blend).
        const lossMethod = this.opts.lossMethod ?? 'auto';
        // MQS applicability, hoisted ABOVE the eigensolve: the MQS loss path consumes
        // only the scalar eps_d from the full-wave solve (never the eigenvector), which
        // is what makes the dispersion-cache fast path below safe on that path. MQS
        // applies to rectangular conductors, explicit ground rects (coplanar GCPW
        // grounds, via slabs, cutout remnants) are handled inside the solve as passive
        // C = 0 return conductors, and on 'auto' a symmetric domain (mode set by
        // the BC). Shaped conductors (coax) have no rect-based skin-band/
        // classification support, so they fall back to the H-field perturbation.
        // Per-face plating is handled inside MQS (surfaceZs weights each face's smooth
        // current by its own impedance), so plating doesn't force perturbation.
        const mqsOk = cr.rects.length > 0 && !cr.rects.some(r => r.shape)
            && cr.rectRoles.some(r => r.is_signal);
        const anyPlating = cr.rectRoles.some(r => r.plating && r.plating.sigma > 0
            && (r.plating.top || r.plating.sides || r.plating.bottom));
        // Refuse a forced 'mqs' override where it cannot apply (with a warning)
        // rather than produce garbage.
        if (lossMethod === 'mqs' && !mqsOk && this._modeWarnings
            && !this._modeWarnings.some(w => w.type === 'mqs-shape')) {
            this._modeWarnings.push({ type: 'mqs-shape', mode, freq: f,
                message: 'MQS conductor loss is not applicable to this geometry ' +
                         '(shaped conductors or no signal conductor), using the perturbation method instead.' });
        }
        // A differential pair without the symmetry-plane mode walls uses the
        // per-conductor-drive MQS: one unit solve per polarity group and an NxN
        // current-matrix solve for the drives realizing the target mode currents
        // (mqs_loss opts.modeCurrents). A single shared drive there would make
        // the two traces one parallel conductor.
        // The multi-drive path needs both polarity groups present in the meshed
        // domain. If one was clipped away, refuse rather than mis-normalize.
        const hasPol = p => cr.rectRoles.some(r => r.is_signal && (r.polarity || 1) * p > 0);
        const mqsMulti = this.solver.is_differential && !this.symmetry
            && mqsOk && hasPol(1) && hasPol(-1);
        const mqsModeSelectable = !this.solver.is_differential || this.symmetry || mqsMulti;
        if (lossMethod === 'mqs' && mqsOk && !mqsModeSelectable && this._modeWarnings
            && !this._modeWarnings.some(w => w.type === 'mqs-asym-pair')) {
            this._modeWarnings.push({ type: 'mqs-asym-pair', mode, freq: f,
                message: 'MQS conductor loss on a full-domain differential pair needs both ' +
                         'traces present in the solved domain, but one polarity group is ' +
                         'missing. Using the perturbation method instead.' });
        }
        // The volume eddy-current solve reads the elements INSIDE the metal, which the
        // mesher only emits when buildMesh's own applicability test (mqsPre) said MQS
        // could apply. That test is derived from the pre-mesh geometry and this one from
        // the built condRect, so they are computed independently and could in principle
        // disagree, on a mesh with no conductor interior MQS would return a plausible
        // but wrong R rather than fail, so fall back to perturbation and say so.
        let useMQS = ((lossMethod === 'mqs' && mqsOk)
            || (lossMethod === 'auto' && (this.symmetry || mqsMulti) && mqsOk))
            && mqsModeSelectable;
        // Both modes need to use MQS if one uses it. Otherwise R12 could be
        // non-physical.
        if (this._noMqsThisSolve) useMQS = false;
        if (useMQS && mesh.condInteriorMeshed === false) {
            useMQS = false;
            if (this._modeWarnings && !this._modeWarnings.some(w => w.type === 'mqs-no-interior')) {
                this._modeWarnings.push({ type: 'mqs-no-interior', mode, freq: f,
                    message: 'MQS conductor loss needs the conductor interiors meshed, but this ' +
                             'mesh was built without it. Using the perturbation method instead.' });
            }
        }
        let eps_d = eps_eff_static, fw = null, eigen_bias = 1;
        if (f >= F_STATIC_MAX) {
            // Seed the eigensolve with the STATIC field (frequency-deterministic), NOT a
            // cross-frequency eigenvector warm-start. The static seed is the quasi-TEM
            // shape (a good initial guess) AND it makes the picked eigenvector independent
            // of evaluation ORDER. With a cross-frequency warm-start, a near-degenerate
            // cluster (e.g. the open homogeneous stripline odd mode + lateral resonances)
            // converges to a different in-subspace mixture depending on which frequency
            // ran last — so the out-of-order interpolating sweep produced order-dependent
            // conductor loss → ripple + non-convergent interpolant (excess solves).
            //
            // TWO-STAGE quasi-TEM pick:
            //   1. Closed-wall pick ('pmc' naturals, st.abc): all-real matrices, so the
            //      WASM real-arithmetic Arnoldi fast path applies — the complex SparseLU
            //      factorization dominates the sweep at these problem sizes, and the real
            //      factorization is several times cheaper.
            //   2. Only when the closed pick FAILS or is AMBIGUOUS (a competing mode with
            //      comparable overlap but different ε — the closed box manufacturing
            //      near-degenerate cavity modes next to the quasi-TEM), ESCALATE to the
            //      radiating first-order ABC pick (=== true on open walls, complex): the
            //      radiation term absorbs outgoing waves, so the open structure's bound
            //      quasi-TEM is found without box-mode competition. Reference output was
            //      validated byte-identical between the closed and the ABC pick across the
            //      whole suite — escalation only changes WHICH solve resolves the ambiguous
            //      minority of points. PEC ('gnd') walls stay PEC either way; a fully
            //      enclosed structure has no open wall and never escalates. The symmetry
            //      plane stays 'pmc' (even/single) / PEC (odd, absent).
            // The FEM is NOT re-assembled per frequency: the system is affine in k² and
            // the ABC Robin term is linear in k₀ = √k², so fullwaveMode combines a cached
            // decomposition A0 + k²·A1 + j·k0·Ar per point (st.femCache / st.femCacheAbc —
            // one slot per BC set so escalation doesn't thrash the closed-pick cache).
            //
            // DISPERSION CACHE (MQS loss path only): there the eigensolve contributes just
            // the scalar eps_d(f) — a very smooth dispersion curve — while the eigenvector
            // is unused. Exact eigensolves fill st.disp with (f, eps_d) anchors; once a
            // leave-one-out check on the neighbouring anchors proves the local PCHIP
            // interpolation error ≤ dispTol, intermediate sweep points interpolate eps_d
            // instead of paying a full eigensolve. Anchor points stay exact, and the
            // interpolating sweep's own RLGC error control still verifies every reported
            // quantity downstream. The perturbation path needs the eigenvector (SIBC
            // projection), so it never interpolates.
            const dispTol = this.opts.dispTol ?? 1e-3;
            const dc = useMQS ? (st.disp || (st.disp = { xs: [], ys: [] })) : null;
            const epsI = dc ? dispersionInterp(dc, f, dispTol) : null;
            // Raw eigenvalue (exact or interpolated); anchored to the static solve below.
            let haveEigen = false;
            if (epsI !== null) {
                eps_d = epsI; haveEigen = true;
            } else {
            let fwErr = null;
            ({ fw, fwErr } = this._eigenPick(st, f, phiEps, eps_eff_static));
            if (fw && fw.eps > 0) {
                eps_d = fw.eps; haveEigen = true;
                if (dc) dispersionInsert(dc, f, fw.eps);
            }
            // The eigensolve failed (or converged to no physical quasi-TEM) at this
            // frequency: the point silently degrades to the quasi-static ε, which
            // kinks the sweep. Surface it — one warning per mode, like the
            // ambiguity warning below.
            if (!fw && this._modeWarnings) {
                const reason = fwErr ? `eigensolver error: ${fwErr.message || fwErr}`
                                     : 'no converged quasi-TEM candidate';
                const msg = `${mode} mode: full-wave eigensolve failed at ${(f / 1e9).toFixed(2)} GHz `
                    + `(${reason}) — falling back to the quasi-static ε_eff=${eps_eff_static.toFixed(3)} `
                    + `for this point; the dispersion curve may show a kink here.`;
                if (!this._modeWarnings.some(w => w.mode === mode && w.type === 'eigensolve')) {
                    this._modeWarnings.push({ type: 'eigensolve', mode, freq: f, message: msg });
                    globalThis.__TRI_DEBUG__ && console.warn('[tri full-wave] ' + msg);
                }
            }
            if (fw && fw.ambiguous && this._modeWarnings) {
                const msg = `${mode} mode: full-wave quasi-TEM pick is ambiguous at ${(f / 1e9).toFixed(2)} GHz — `
                    + `picked ε_eff=${fw.bestEps.toFixed(3)} (overlap ${fw.bestOvl.toFixed(2)}) but a competing mode `
                    + `ε_eff=${fw.altEps.toFixed(3)} (overlap ${fw.altOvl.toFixed(2)}) carries comparable weight; `
                    + `static ε_eff=${eps_eff_static.toFixed(3)}. The quasi-TEM has likely fragmented across degenerate `
                    + `modes (inhomogeneous fill) — reported ε_eff/Z0 may be unreliable.`;
                if (!this._modeWarnings.some(w => w.mode === mode)) {
                    // Surfaced to the UI via result.warnings / solver.modeWarnings; the console
                    // line is dev-only (avoid duplicate user-facing noise).
                    this._modeWarnings.push({ mode, freq: f, bestEps: fw.bestEps, altEps: fw.altEps,
                        bestOvl: fw.bestOvl, altOvl: fw.altOvl, staticEps: eps_eff_static, message: msg });
                    globalThis.__TRI_DEBUG__ && console.warn('[tri full-wave] ' + msg);
                }
            }
            }   // end exact-eigensolve branch (dispersion-cache miss)
            if (haveEigen) {
                eigen_bias = this._eigenBias(mode, st, phiEps, eps_eff_static,
                    (f === F_STATIC_MAX && fw && fw.eps > 0) ? fw : null);
                eps_d /= eigen_bias;
            }
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
        const Zr = calculate_Zrough(f > 0 ? f : 1, sigma, rq);

        // Per-face plating: split the loss surface into surface-impedance groups.
        // Used by the perturbation path below (when MQS doesn't apply) to evaluate
        // and sum the loss per group; >1 group (e.g. top/sides plated, bottom bare)
        // means a face-dependent impedance. (MQS handles per-face plating itself via
        // surfaceZs, so plating no longer forces the perturbation path.)
        // The loss-edge mask and the surface-group face classification are purely
        // geometric — cache them per mode (validated against mesh/fm identity);
        // only the per-frequency Zs grouping is re-evaluated each point.
        const lc = st.lossCache || (st.lossCache = {});
        if (lc.mesh !== mesh || lc.fm !== fm) {
            lc.mesh = mesh; lc.fm = fm;
            lc.mask = buildLossEdges(mesh, fm, cr);
            lc.clsMask = null;   // invalidate the surface-group classification too
        }
        const lossEdgeMask = f > 0 ? lc.mask : null;

        // (mqsOk / anyPlating / useMQS are computed above the eigensolve —
        // the dispersion-cache fast path there needs the MQS applicability.)
        let R_total = 0, L_internal = 0;
        // 'mqs' only when R_total came from the eddy solve (directly or from its
        // interpolated anchors), solveAt compares this across modes.
        let lossVia = 'perturbation';
        // Smallest conductor cross-section dimension (gates the skin-band remesh and
        // the static/SIBC blend below).
        let minDim = Infinity;
        for (const c of cr.rects) {
            // A complement conductor (coax shield) is a zero-thickness PEC shell whose
            // bbox is the whole cavity, it has no cross-section to gate anything on.
            if (c.shape && isComplement(c.shape)) continue;
            minDim = Math.min(minDim, c.xmax - c.xmin, c.ymax - c.ymin);
        }
        if (f > 0 && useMQS) {
            // The volume eddy solve runs the conductor BODY at the bulk metal σ;
            // plating is a SURFACE effect applied only through surfaceZs (relative to
            // the bulk Rs). So with plating, use the bulk σ here — NOT effectiveSurface's
            // plating σ (which would bake plating into the body and make surfaceZs a
            // no-op). The skin band is sized by the matching bulk skin depth.
            const mqsSigma = anyPlating ? (s.sigma_cond ?? 5.8e7) : sigma;
            const mqsDelta = anyPlating ? Math.sqrt(2 / (omu * mqsSigma)) : delta;
            // Skin-band element size at the conductor surface (xδ) and band width (xδ).
            // Resolve the skin layer to bandDelta*δ within mqsBand*δ of each surface,
            // with the target relaxing away from the surface (mqsBandDepthSlope, see
            // buildSkin / refineSkinBand).
            //
            // The surface target is deliberately ~δ rather than a fraction of it: with
            // the depth grading in place, R converges on the grading , not on the surface
            // element size, a P2 element about a skin depth across already carries
            // exp(-d/δ) well. The previous uniform-target band needed 0.6δ elements
            // everywhere inside 1.5δ and still ran into its triangle cap (leaving R
            // +0.6% at 10 GHz and +1.3% at 40 GHz on the reference differential GCPW).
            // The graded band converges inside the same budget at ~0.1%, with a smaller
            // mesh and a faster MQS solve.
            const bandDelta = this.opts.mqsBandDelta ?? 1.25;
            // Skin-band triangle budget: how many triangles the size-aware band
            // refinement may add on top of the base mesh, a memory guard for
            // large-perimeter/small-δ geometries (hitting it leaves the band
            // partially resolved rather than aborting, see mqs-band-capped below).
            // Not an absolute cap, otherwise base mesh near the limit would not
            // get band at all.
            const mqsMaxTris = this.opts.mqsMaxTris ?? 80000;
            const mqsBand = this.opts.mqsBand ?? 1.5;
            // Skin-band refinement runs in two passes, trace first: the signal
            // band may add mqsMaxTris triangles to the base mesh, then the passive
            // ground rects (GCPW coplanar grounds, via slabs) refine on top with
            // their own additional budget and a distance-graded target, full
            // bandDelta*δ resolution only within ~Dfine of the signal, where
            // the return current concentrates at the slot edge, growing
            // linearly beyond (the ground surface current decays away from the
            // slot, and past the via fence it is negligible). The grading
            // keeps the band cost independent of how far the ground pour
            // extends, the pass order keeps a large ground perimeter from
            // eating the trace's budget (one shared budget would abort the
            // whole refinement with the trace still coarse, worse than no
            // ground band at all).
            const sigRects = [], gndRects = [];
            cr.rects.forEach((r, i) => (cr.rectRoles[i].is_signal ? sigRects : gndRects).push(r));
            let gndGrading = null;
            if (gndRects.length) {
                // Dfine keys off the signal-ground clearance (the GCPW slot width):
                // points on the ground's inner edge sit exactly gapMin from the
                // signal, so Dfine = 1.5*gapMin keeps ~half a slot-width of the
                // ground at full resolution beyond its inner edge.
                let gapMin = Infinity;
                for (const gr of gndRects) for (const sr of sigRects) {
                    const dx = Math.max(0, Math.max(gr.xmin, sr.xmin) - Math.min(gr.xmax, sr.xmax));
                    const dy = Math.max(0, Math.max(gr.ymin, sr.ymin) - Math.min(gr.ymax, sr.ymax));
                    gapMin = Math.min(gapMin, Math.hypot(dx, dy));
                }
                if (!isFinite(gapMin)) gapMin = 0;
                gndGrading = { sigRects,
                    Dfine: (this.opts.mqsGndFine ?? 1.5) * gapMin,
                    slope: this.opts.mqsGndSlope ?? 0.5 };
            }
            const gndBudget = this.opts.mqsGndMaxTris ?? Math.floor(mqsMaxTris / 2);
            // Depth grading across the band (see refineSkinBand). The target size
            // relaxes with distance from the metal surface.
            const depthSlope = this.opts.mqsBandDepthSlope ?? 3;
            // Bisection-pass budget for the band refinement. refineSkinBand
            // stops on its own as soon as nothing is left above the size target
            // (or the triangle cap trips), so an unused pass costs nothing. The
            // triangle cap is the real bound on runaway growth.
            const bandPasses = this.opts.mqsBandPasses ?? 20;
            const buildSkin = (base, dlt) => {
                let m = refineSkinBand(base, { rects: sigRects }, dlt, bandPasses, mqsBand, bandDelta * dlt,
                    base.nTris + mqsMaxTris, null, depthSlope);
                // Each refineSkinBand stamps its own bandTrunc on the mesh it returns,
                // so collect the signal band's before the ground pass overwrites it,
                // and leave the merged list (or null) on the mesh that gets cached.
                const trunc = [];
                if (m.bandTrunc) trunc.push({ band: 'signal', ...m.bandTrunc });
                if (gndRects.length) {
                    m = refineSkinBand(m, { rects: gndRects }, dlt, bandPasses, mqsBand,
                        bandDelta * dlt, m.nTris + gndBudget, gndGrading, depthSlope);
                    if (m.bandTrunc) trunc.push({ band: 'ground', ...m.bandTrunc });
                }
                m.bandTrunc = trunc.length ? trunc : null;
                return m;
            };
            // DEFAULT: f_max-reuse. Build the skin mesh at the HIGHEST frequency seen
            // (finest target, narrowest band) and reuse it for all lower frequencies.
            // The fine near-surface band plus the conductor interior (always meshed on
            // this path, see the guard on condInteriorMeshed above) resolve the
            // lower-frequency (more uniform) current accurately, validated
            // against the reference suite to a 100× frequency mismatch. One mesh ⇒ smooth
            // R(f) with no per-frequency-remesh "dip", at moderate size (the high-f band
            // is narrow), unlike a whole-range mesh which over-resolves the low-f band.
            // opts.mqsCacheMesh === false falls back to a per-frequency remesh (marginally
            // faster on wide sweeps, but can show a small non-monotonic wiggle).
            // Band-sizing skin depth: when a sweep announces its maximum frequency
            // (solver._sweepFmax, set by solve_sweep / InterpolatingSweep), size the
            // band for the WHOLE sweep up front — an ascending discrete sweep then
            // builds ONE skin mesh at the first point and every later point reuses
            // it (plus its cached MQS assembly), instead of re-refining per point.
            // Without the hint this reduces to the current frequency (old behavior:
            // rebuild whenever a higher frequency appears).
            const fBand = Math.max(f, this.solver._sweepFmax || 0);
            const deltaBand = fBand > f ? Math.sqrt(2 / (2 * Math.PI * fBand * MU0 * mqsSigma)) : mqsDelta;
            // The band is built at every skin depth. refineSkinBand is size-aware
            // and leaves the base mesh untouched where it already resolves
            // bandDelta*delta, so at delta well above the conductor dimensions it
            // costs one marking scan. Skipping it for delta >= thickness would be
            // wrong for a thin conductor on a coarse base mesh: the current is
            // uniform across the thickness but still crowds toward the edges over
            // a lateral scale of order delta^2/t, which the sliver elements inside
            // the conductor cannot follow (1 um trace at 3 GHz: R +12% without
            // the band).
            let mqsMesh, bandTrunc = null;
            if (this.opts.mqsCacheMesh === false) {
                mqsMesh = buildSkin(mesh, mqsDelta);
                bandTrunc = mqsMesh.bandTrunc;
            } else {
                // The skin mesh depends only on the geometry and the band skin depth —
                // NOT on the mode (the odd/even BC enters the MQS solve, not the mesh) —
                // so one cache on the backend serves both modes: the band refinement
                // runs once per sweep instead of once per mode.
                const sc = (this._skinCache && this._skinCache.base === mesh) ? this._skinCache
                    : (this._skinCache = { base: mesh, dB: Infinity, mesh: null });
                if (!sc.mesh || deltaBand < sc.dB * (1 - 1e-9)) {   // higher freq appeared → rebuild finer
                    sc.dB = deltaBand;
                    sc.mesh = buildSkin(mesh, deltaBand);
                }
                mqsMesh = sc.mesh;
                bandTrunc = mqsMesh.bandTrunc;
            }
            // A band that stopped short of bandDelta*δ leaves the loss integral
            // on an under-resolved skin layer, and it is invisible from the
            // outside. The answer looks converged. Read the flag off the mesh,
            // the skin mesh is cached across modes and sweep points, and
            // a point that reuses it must warn just the same.
            if (bandTrunc && this._modeWarnings
                && !this._modeWarnings.some(w => w.type === 'mqs-band-capped')) {
                // The message goes straight to the UI log, so it says what is
                // degraded and what to use instead. Neither mqsMaxTris nor
                // mqsBandPasses is reachable from the UI, so naming them would
                // only be noise.
                this._modeWarnings.push({ type: 'mqs-band-capped', mode, freq: f,
                    detail: bandTrunc, message:
                    `Conductor-loss accuracy is reduced. Fully resolving the ` +
                    `${(deltaBand * 1e6).toFixed(2)} µm skin layer on this geometry needs a ` +
                    `finer mesh than the full-wave solver's budget allows. ` +
                    `Conductor loss can be several percent high.` });
            }
            // R(f) and L_internal(f) anchor caches: the sweep consumes only these
            // two scalars from the solve, so past the first few anchor frequencies
            // interpolate them the same way the eigensolve's eps_d dispersion cache
            // does (log-space PCHIP + leave-one-out gate) and skip the MQS
            // factorization entirely. Keyed to the skin mesh so a finer rebuild
            // (higher f seen) discards anchors from the coarser mesh.
            //
            // Each cache is gated on its own leave-one-out check and both must
            // pass: a hit on R alone would otherwise force L_internal to be
            // reconstructed from it, which is exactly the R/omega assumption this
            // path must avoid.
            const rTol = this.opts.mqsInterpTol ?? 2e-3;
            const rc = (st.mqsR && st.mqsR.mesh === mqsMesh) ? st.mqsR
                : (st.mqsR = { mesh: mqsMesh, xs: [], ys: [], log: true });
            const lc = (st.mqsL && st.mqsL.mesh === mqsMesh) ? st.mqsL
                : (st.mqsL = { mesh: mqsMesh, xs: [], ys: [], log: true });
            const rInterp = dispersionInterp(rc, f, rTol);
            const lInterp = dispersionInterp(lc, f, rTol);
            if (rInterp !== null && rInterp > 0 && lInterp !== null && lInterp > 0) {
                R_total = rInterp;
                L_internal = lInterp;
                lossVia = 'mqs';                      // interpolated MQS anchors
            } else {
            // Per-face plating ⇒ weight each face's smooth current by its own surface
            // impedance (surfaceZs); otherwise the uniform roughness factor (Rq).
            // oddSymmetry (the A = 0 wall at x = 0) only exists on a half-domain
            // solve, the full-domain differential path selects the mode via
            // modeCurrents instead, and shares one cache across both modes. The
            // assembly and the per-frequency unit solves are mode-independent
            // there, so the second mode at each frequency skips the factorization.
            const mqsOpts = { wallPEC: cr.wallPEC || null, wallThick: cr.wallThick || null,
                              topGround: !!(cr.wallPEC && cr.wallPEC.top),   // legacy fallback
                              oddSymmetry: this.symmetry && mode === 'odd',
                              diffPair: !!s.is_differential,
                              // Frequency-invariant assembly cache (validated
                              // against the mesh object inside mqsConductorLoss, so a
                              // skin-mesh rebuild invalidates it automatically).
                              cache: mqsMulti ? (this._mqsMultiCache || (this._mqsMultiCache = {}))
                                              : (st.mqsCache || (st.mqsCache = {})) };
            if (mqsMulti) mqsOpts.modeCurrents = this._mqsModeCurrents(mode, st);
            if (anyPlating) mqsOpts.surfaceZs = buildFaceZs(s, cr, f);
            else mqsOpts.Rq = rq;
            let mqs = null;
            try {
                mqs = mqsConductorLoss(mqsMesh, cr, f, mqsSigma, this.ctx.helpers.solveComplexSymmetric, 0, mqsOpts);
            } catch (e) {
                // Surface the downgrade instead of silently falling back: an MQS solve
                // failure here is almost always the factorization exhausting the WASM heap
                // on a very large skin mesh, and the perturbation estimate it falls
                // back to is markedly less accurate (corner-sensitive).
                mqs = null;
                if (this._modeWarnings) this._modeWarnings.push({ type: 'mqs-solve-failed', mode, freq: f, message:
                    `MQS conductor-loss solve failed (${String(e && e.message || e).slice(0, 120)}); ` +
                    `using the perturbation loss estimate instead. Reduce Max Nodes if this persists.` });
            }
            // Accept MQS only if physically sane.
            if (mqs && isFinite(mqs.R_total) && mqs.R_total > 0 && isFinite(mqs.L_loop) && mqs.L_loop > 0
                && isFinite(mqs.X_total) && mqs.X_total > 0) {
                R_total = s.is_differential ? mqs.R_total / 2 : mqs.R_total;
                lossVia = 'mqs';
                // Internal inductance from the volume solve: the MQS loop inductance
                // minus the loop inductance of the same box with perfect conductors
                // (mqsPecInductance, one real solve per skin mesh), plus the skin
                // reactance of the A = 0 walls. The difference is taken on one mesh
                // with one set of walls, so the external part cancels exactly and
                // what remains holds at every skin depth. The surface reactance
                // X_total/omega would only be right for delta well below the
                // conductor thickness (Zs = Rs(1+j)); a film thinner than delta has
                // R/omega orders of magnitude above its true internal inductance.
                const pc = mqsOpts.cache;
                if (!(pc.pecMesh === mqsMesh && pc.pecMode === mode)) {
                    pc.Lpec = mqsPecInductance(mqsMesh, cr, this.ctx.helpers.solveSparseMulti, mqsOpts);
                    pc.pecMesh = mqsMesh; pc.pecMode = mode;
                }
                L_internal = Math.max(0, mqs.L_loop - pc.Lpec + mqs.L_wall);
                dispersionInsert(rc, f, R_total);
                dispersionInsert(lc, f, L_internal);
            }
            }   // end R-cache miss (exact MQS solve)
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
        // DC resistance from the geometry lists — the SAME convention as the FDM
        // backend (field_solver.calculate_conductor_loss): current flows through the
        // signal cross-section and returns through the ground cross-section in
        // series. cr.rects can NOT be used for this: it contains coplanar ground
        // rects (which are return path, not parallel signal metal) and omits
        // wall-absorbed grounds entirely.
        const sigmaDC = s.sigma_cond ?? 5.8e7;
        let sigArea = 0, gndArea = 0;
        for (const c of s.conductors) {
            // shapeArea is the bbox product for a plain rect (unchanged) but
            // the true cross-section for a shaped one. A complement shell
            // returns 0: the coax shield is modelled as infinitely thick, so it
            // carries no DC resistance.
            const a = shapeArea(c);
            if (c.is_signal) sigArea += a; else gndArea += a;
        }
        // Mode-aware per-line convention for a differential pair (mirrors
        // field_solver.calculate_conductor_loss): sigArea sums both traces, the odd
        // mode's DC return is the partner trace (no net ground current), the even
        // mode returns 2I through the ground.
        let R_dc;
        if (s.is_differential && (mode === 'odd' || mode === 'even')) {
            const R_trace = sigArea > 0 ? 2 / (sigmaDC * sigArea) : 0;
            R_dc = mode === 'odd' ? R_trace
                 : R_trace + (gndArea > 0 ? 2 / (sigmaDC * gndArea) : 0);
        } else {
            R_dc = (sigArea > 0 ? 1 / (sigmaDC * sigArea) : 0)
                 + (gndArea > 0 ? 1 / (sigmaDC * gndArea) : 0);
        }
        if (f === 0) {
            // Pure DC point: no skin effect, R is the geometric DC resistance
            // (the FDM backend returns the same; this used to return R = 0).
            R_total = R_dc;
            L_internal = 0;
        }
        if (f > 0 && R_total === 0) {
            // Per-edge surface groups (one group = one surface impedance), evaluated
            // only HERE — the MQS path above never reads them (it handles per-face
            // plating itself via surfaceZs), so grouping per sweep point up front
            // wasted an O(nEdges) classification plus surface-impedance evaluations.
            // The groups from buildSurfaceGroups are used even when uniform: a
            // fully-plated surface then correctly keeps the LAYERED plating-over-bulk
            // impedance (substituting effectiveSurface's solid-plating-metal Zs here
            // used to lose the bulk underneath).
            const surf = buildSurfaceGroups(s, mesh, fm, cr, lossEdgeMask, f, lc);
            const groups = surf.groups.length ? surf.groups
                : [{ Zs: { re: Zr.re, im: Zr.im }, mask: lossEdgeMask }];
            // The loss routines run at the base bulk σ; plating/roughness enter only
            // through the per-group fRg = Re(Zs)/Rs(σ_base) scaling (the group Zs from
            // makePlatingZs is defined relative to the base metal).
            const sigmaRef = s.sigma_cond ?? 5.8e7;
            const deltaRef = Math.sqrt(2 / (omu * sigmaRef));
            const RsRef = 1 / (sigmaRef * Math.min(deltaRef, 1e30));

            const useFW = !!fw && lossMethod !== 'static';
            // Galerkin H-projection of the eigenmode, computed ONCE per frequency and
            // shared by every plating group's SIBC evaluation (solveConductorLoss used
            // to re-project per group). Replaces the analyzeTriMode call whose only
            // consumed output was the Poynting power — and whose voltage-path contour
            // walk could crash on geometries without a node at its hardcoded ground.
            let projH = null, Pfw = 0;
            if (useFW) {
                let gamma = csqrt(fw.g2Re, fw.g2Im);
                if (gamma.im < 0) gamma = { re: -gamma.re, im: -gamma.im };
                try {
                    projH = projectH(mesh, fm, fw.vRe, fw.vIm, gamma, f, this.ctx.wasmSolver,
                        st.projCache || (st.projCache = {}));
                    Pfw = Math.abs(computePoyntingFromProjectedH(mesh, fm, fw.vRe, fw.vIm,
                        projH.htRe, projH.htIm, omu, projH.hDofs));
                } catch { projH = null; Pfw = 0; }
            }
            const haveFW = useFW && projH && Pfw > 1e-30 && minDim > 0;
            const wFW = haveFW ? Math.max(0, Math.min(1, (0.15 - deltaRef / minDim) / (0.15 - 0.08))) : 0;

            // A differential drive carries TWICE the transmitted power of the
            // single-conductor normalization the static estimator uses (two ±1 V
            // traces), while its loss integral covers both traces — so its R_ac is
            // reported per drive pair and must be halved to the per-trace mode R
            // (same convention as the MQS mqs.R_total/2 above and the FDM backend).
            // The SIBC estimator self-normalizes by the eigenmode's own Poynting
            // power, so it needs no correction.
            const driveScale = s.is_differential ? 0.5 : 1;
            // Both halves of the surface impedance, out of one loss integral per group.
            // The perturbation integral is ∮|K|²dl/(2P) evaluated at the reference Rs and
            // is linear in the surface impedance, so the same integral scaled by
            // Re(Zs)/Rs is the resistance and scaled by Im(Zs)/Rs is the surface
            // reactance. Carrying the second scaling costs nothing, which is
            // why the reactance is accumulated here rather than reconstructed
            // from R_ac afterwards (see L_internal below).
            let R_ac = 0, X_ac = 0;
            for (const g of groups) {
                const lossS = staticConductorLoss(cr.rects, f, sigmaRef, mesh, fm, phiEps,
                    Z0, eps_eff_static, eps_d, g.mask);
                let baseg = lossS.R_ac * driveScale;   // driveScale: static estimator only
                if (haveFW && wFW > 0) {
                    const lossW = solveConductorLoss(cr.rects, f, sigmaRef, mesh, fm, fw.vRe, fw.vIm,
                        fw.g2Re, fw.g2Im, Pfw, Z0, mesh.epsMap, g.mask, projH, this.ctx.wasmSolver);
                    baseg = (1 - wFW) * baseg + wFW * lossW.R_ac;
                }
                R_ac += baseg * (RsRef > 0 ? g.Zs.re / RsRef : 1);
                X_ac += baseg * (RsRef > 0 ? g.Zs.im / RsRef : 1);
            }
            R_total = Math.sqrt(R_dc * R_dc + R_ac * R_ac);
            // Internal (skin) inductance is the surface reactance over ω.
            // Smooth metal in the skin regime has Zs = Rs(1+j), which makes
            // that exactly R_ac/ω.  Roughness and plating tilt Zs off 45°, and
            // only X_ac follows: at 2 µm rms Im(Zs)/Re(Zs) reaches ~5.5, so
            // R_ac/ω recovered under a fifth of the true increment and the line
            // came out measurably fast.  Neither the static-field nor the SIBC
            // routine returns a usable L_int of its own across the skin
            // transition, which is why it is rebuilt here at all. The cap
            // bounds it near dc, below the skin regime where no surface model
            // holds.
            L_internal = (omega > 0) ? Math.max(0, Math.min(X_ac / omega, 0.5 * L_external)) : 0;
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
            // The static-anchoring factor divided out of the eigensolve's eps (1 below
            // F_STATIC_MAX or when the anchor was unusable): eps_eff_mode * eigen_bias is
            // the raw eigenvalue, kept visible for diagnostics of the eigensolve itself.
            eigen_bias,
            RLGC: { R: R_total, L: L_total, G, C },
            Zc,
            alpha_c, alpha_d, alpha_total: alpha_c + alpha_d,
            L_internal, L_external, lossVia,
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
    _applyCausal(f, skipFields = false) {
        const s = this.solver;
        const fref = this.opts.causalFref ?? 1e9;
        const causalDiel = s.dielectrics.map(d => {
            const er = d.epsilon_r, td = d.tan_delta || 0;
            // `shape` must ride along: without it the rebuilt dielectric is only a
            // bounding box, and the causal re-tag would label everything outside the
            // real body (the polygon's vertex "horns") as air.
            const rect = { x_min: d.x_min, x_max: d.x_max, y_min: d.y_min, y_max: d.y_max, shape: d.shape || null };
            if (Math.abs(er - 1) < 1e-6 || Math.abs(td) < 1e-10) return { ...rect, epsilon_r: er, tan_delta: td };
            const { eps_real, tand_actual } = djordjevic_sarkar(f, er, td, fref);
            return { ...rect, epsilon_r: eps_real, tan_delta: tand_actual };
        });
        const { epsMap, lossMap } = tagMaterials(this.mesh, causalDiel);
        this.mesh.epsMap = epsMap; this.mesh.lossMap = lossMap;
        // A waveguide has no static solve to redo and no field to resample: kc is purely
        // geometric, so a causal shift in εr moves beta and the losses (both read the maps
        // above per frequency) without touching the mode pattern or the cutoff.
        if (this._isWG) return;
        // Asymmetric differential pair: the odd/even drive is NOT the modal basis, so re-derive
        // the genuine line modes (and the physical [C]/[L] matrices that drive the 4-port
        // S-parameters) under the causal permittivity — exactly as the initial build does.
        // _prepareStaticModal rebuilds this._static[*] (fields, eps_eff_static, W_loss) and
        // this._modalPhys from the now-causal mesh.epsMap. It returns false for an electrically
        // symmetric or velocity-degenerate pair, in which case the odd/even drive below is the
        // correct basis (and any _modalPhys it set is already causal-consistent).
        if (s.is_differential && !this.symmetry &&
            this._prepareStaticModal(this._staticGrid, skipFields ? 'reuse' : 'build')) return;
        for (const mode of this.modeNames) {
            const st = this._static[mode];
            const pot = drivePotentials(this.condRect, mode, this.symmetry);
            const phiEps = solveTriStatic(this.mesh, st.fm, epsMap, pot, this.ctx.wasmSolver);
            const W_air = st.C0 / (st.kC * eps0);   // geometric, eps-independent
            st.phiEps = phiEps;
            st.eps_eff_static = computeTriEnergy(phiEps, this.mesh, epsMap) / W_air;
            st.W_loss = computeTriEnergy(phiEps, this.mesh, lossMap);
            // Refresh the resampled plot fields so contour plots track the causal ε.
            // Skipped on mid-sweep calls (skipFields): each sweep point's resample
            // is immediately overwritten by the next, and the plot the user sees
            // comes from the main solve — the physics above never depends on it.
            // (Unless nothing has resampled yet: _prepareStatic leaves fields null on
            // the causal path precisely because this call supersedes it.)
            if (!skipFields || !st.fields) {
                const parity = this.symmetry ? (mode === 'odd' ? 'odd' : 'even') : null;
                st.fields = resampleStatic(this.mesh, phiEps, this.domain,
                    { resolution: this.opts.resolution, parity, grid: this._staticGrid, rects: this._plotRects });
            }
        }
    }

    solveAt(f, opts = {}) {
        if (!this.mesh) throw new Error('TriBackend: buildMesh() must be awaited before solving (mesh not built).');
        if (this.solver.use_causal_materials) this._applyCausal(f, opts.skipFieldResample === true);
        this._modeWarnings = [];
        // Surface the buildMesh-time accuracy warning (failed verification certificate)
        // through the same per-solve channel as the mode warnings, so the UI logs it.
        if (this._certWarn) this._modeWarnings.push(this._certWarn);
        // A waveguide never reaches _modeAtFreq, so its quasi-static fallback (which would
        // be physically meaningless below cutoff) is unreachable for this medium.
        let modes = this._isWG ? [this._waveguideModeAtFreq(f)]
            : this.modeNames.map(m => this._modeAtFreq(m, f));
        // Use the same loss (MQS or pertrubation) for both modes. Keeps R12
        // physical. If MQS fails for one mode, then use pertrubation for both.
        if (!this._isWG && modes.length > 1 && !this._noMqsThisSolve) {
            const via = modes.map(m => m.lossVia);
            if (via.some(v => v === 'mqs') && via.some(v => v !== 'mqs')) {
                // Drop the MQS R/L anchors the succeeding mode already cached
                // so that they aren't used in frequency sweep other points.
                for (const m of this.modeNames) {
                    const st = this._static && this._static[m];
                    if (st) { st.mqsR = null; st.mqsL = null; }
                }
                this._noMqsThisSolve = true;
                try {
                    modes = this.modeNames.map(m => this._modeAtFreq(m, f));
                } finally {
                    this._noMqsThisSolve = false;
                }
                if (this._modeWarnings && !this._modeWarnings.some(w => w.type === 'mqs-mixed-modes')) {
                    this._modeWarnings.push({ type: 'mqs-mixed-modes', freq: f, message:
                        `The MQS conductor-loss solve succeeded for some modes of this pair and ` +
                        `failed for others (${this.modeNames.map((m, i) => `${m}:${via[i]}`).join(', ')}). ` +
                        `Mixing the two estimators corrupts the mutual resistance R12, so all modes ` +
                        `were recomputed with the perturbation estimate instead.` });
                }
            }
        }
        const result = { modes };
        if (this._modeWarnings.length) result.warnings = this._modeWarnings;
        this.solver.modeWarnings = this._modeWarnings;
        if (this.solver.is_differential) {
            const odd = modes.find(m => m.mode === 'odd');
            const even = modes.find(m => m.mode === 'even');
            result.Z_diff = 2 * odd.Z0;
            result.Z_common = even.Z0 / 2;
            // Genuine asymmetric [C]/[L] (with Tv-reconstructed [R]/[G]) when an asymmetric
            // physMatrix was computed for this pair, else the symmetric odd/even reconstruction.
            // Shared with the S-parameter path so RLGC_matrix matches what drives the S-params.
            result.RLGC_matrix = buildPhysicalRLGC(odd.RLGC, even.RLGC, this._modalPhys);
            // True physical 2×2 [C]/[L] for the asymmetric MTL 4-port S-parameter path.
            if (this._modalPhys) result.physMatrix = this._modalPhys;
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
        if (!this.mesh) throw new Error('TriBackend: buildMesh() must be awaited before solving (mesh not built).');
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
        // Absolute ceiling (mirrors the FDM "Problem too large" guard) — return a clean
        // error instead of letting the eigensolver heap-abort(). The mesh is already held to
        // the triangle budget, so this only trips on extreme settings.
        if (eigenSolveBytes(N, fem) > MAX_SOLVE_BYTES)
            return { modes: [], nconv: 0, error: `Problem too large for the full-wave solver (~${(eigenSolveBytes(N, fem) / 1e9).toFixed(1)} GB). Reduce Max Nodes or the domain/frequency.` };
        // Static-field "drive" seeds for every physical conductor excitation (single-ended:
        // one; differential: even AND odd). A genuine quasi-TEM mode matches ONE of these
        // (high overlap); a spurious discrete gradient mode matches none. All modes share
        // this fm (symmetry off), so each drive's phiEps projects into the same DOF space.
        // A waveguide has no conductor drive: the analytic fundamental takes the seed's
        // place, so the `overlap` column becomes "overlap with the fundamental" instead of
        // "overlap with the static quasi-TEM field", the same thing it means elsewhere.
        const seeds = this._isWG
            ? [this._wg.seed]
            : this.modeNames.map(m => staticToEdgeDofs(this._static[m].phiEps, mesh, fm));
        const seed = seeds[0];
        // Shift-invert centered on the quasi-TEM eigenvalue (-k²·eps_static) converges the
        // fundamental and the nearby higher-order modes first. Widen the Krylov subspace so
        // several modes resolve at once. ncvMax lets the eigensolver keep DOUBLING the
        // subspace (reusing its factorization) until nev pairs truly converge: cavity /
        // higher-order modes sit far from the quasi-TEM shift (e.g. ε_eff≈0.5 vs a shift at
        // 3.4 on an enclosed εr=4.4 line) and need a much larger subspace than the
        // fixed 20-vector pass, whose strict residual gate silently dropped them.
        // For a waveguide the equivalent centre is the fundamental's own eigenvalue,
        // γ² = kc² - k²*εr, which is known exactly rather than estimated from a static solve.
        const shift = this._isWG
            ? this._wg.kc * this._wg.kc - k2 * maxEpsRe(mesh.epsMap)
            : -k2 * eps_eff_static;
        const ncv = Math.min(Math.max(2 * nev + 1, 20), N - 1);
        const ncvMax = Math.min(320, N - 1);
        let res;
        try { res = this.ctx.helpers.solveGeneralized(N, fem.csrA, fem.csrB, [shift, 0], nev, ncv, seed, ncvMax); }
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
        // Refinements on the γ²<0 branch:
        //   • ε_eff ≈ 0 (γ² negligible vs the k²·ε mode scale, only escaped the absolute
        //     nullThresh band) → null space, not a guided mode.
        //   • OPEN domain (any radiating ABC wall) has no bound modes below the light line —
        //     ε_eff < 1 there is the radiation continuum of the truncated domain → spurious.
        //   • enclosed (all-PEC) box keeps sub-light-line modes (cavity modes above cutoff).
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
                const epsCand = -g2Re / k2;
                // ε_eff ≈ 0 → null-space residue: |γ²| is negligible against the k²·ε ~ 10⁵
                // mode scale (it only escaped the absolute nullThresh band), and such a mode
                // is well separated from the genuine guided/cavity modes — classify it as
                // null space rather than a real propagating mode. Above that, an OPEN domain
                // still has no bound mode below the light line (ε_eff < 1 = radiation
                // continuum) → spurious; an enclosed box keeps cavity modes as propagating.
                if (epsCand < 1e-2) status = 'nullspace';
                else { eps_eff = epsCand; status = (!isEnclosed && eps_eff < 1.0) ? 'spurious' : 'propagating'; }
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
    //
    // The grid is a fine UNIFORM grid whose sample points are CELL CENTERS spanning exactly
    // the domain. We deliberately do NOT reuse the FDM mesher's graded grid here: its bulk
    // cells are coarse (~0.5 mm wide at the side walls), and Plotly's heatmap centers each
    // cell on its coordinate — so a wall cell whose centre sits on the domain edge renders
    // half its width (~0.28 mm) PAST the geometry, drawing field over the gridlines beyond
    // the structure. Cell-centred samples (x_min+dx/2, x_max-dx/2) make the heatmap fill
    // the domain exactly, with no spill beyond the walls.
    getModeField(sortedIdx) {
        const ms = this._modesState;
        if (!ms || !ms.list[sortedIdx]) return null;
        const m = ms.list[sortedIdx];
        const d = this.domain;
        const res = this.opts.resolution || 320;
        const ax = d.x_max - d.x_min, ay = d.y_max - d.y_min;
        const aspect = ax / (ay || 1);
        const nx = aspect >= 1 ? res : Math.max(80, Math.round(res * aspect));
        const ny = aspect >= 1 ? Math.max(80, Math.round(res / aspect)) : res;
        const dx = ax / nx, dy = ay / ny;
        const x = new Float64Array(nx), y = new Float64Array(ny);
        for (let i = 0; i < nx; i++) x[i] = d.x_min + dx * (i + 0.5);
        for (let j = 0; j < ny; j++) y[j] = d.y_min + dy * (j + 0.5);
        return resampleModeField(this.mesh, ms.fm, m.vRe, m.vIm, this.domain,
            { resolution: this.opts.resolution, grid: { x, y } });
    }

}

