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
import { buildOccMeshFromGeometry, tagMaterials, validateTriMesh } from './occ_to_mesh.js';
import { buildTriFreedomMap, solveTriStatic, computeTriEnergy, refineTriMesh,
         markTrianglesForRefinement, computeTriP2StaticMatrices,
         staticToEdgeDofs, assembleTriFEM, assembleTriFEMDecomposed,
         femFromDecomposition } from './tri_fem.js';
import { staticConductorLoss, solveConductorLoss, computeHtZZMetric,
         projectH, computePoyntingFromProjectedH } from './conductor_loss.js';
import { csqrt } from './fem_core.js';
import { mqsConductorLoss, refineSkinBand } from './mqs_loss.js';
import { checkMeshQuality } from './tri_mesh.js';

// Below this frequency the full-wave eigensolve degenerates near DC; use the
// static (variational) solve instead. Above it, use the full-wave eigenmode for
// the dispersive effective permittivity.
const F_STATIC_MAX = 100e6;
import { calculate_Zrough, calculate_Zrough_layered } from '../surface_roughness.js';
import { resampleStatic, resampleModeField, buildGridFromMesh } from './resample.js';
import { Complex } from '../complex.js';
import { classifyModalDecomposition } from '../geometry_symmetry.js';
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
function eigenSolveBytes(N, nnz) {
    return 2 * 10 * (12 * nnz + 20 * N);
}

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

function buildSurfaceGroups(solver, mesh, fm, condRect, baseMask, freq, cache = null) {
    const { nodes, edges, nEdges } = mesh;
    const rects = condRect.rects || [];
    const tol = platingTol(condRect);
    const { Zbare, zForFace } = makePlatingZs(solver, condRect, freq);
    const FACES = ['top', 'bottom', 'sides'];

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
                if (x0 < r.xmin - tol || x0 > r.xmax + tol || x1 < r.xmin - tol || x1 > r.xmax + tol) continue;
                if (y0 < r.ymin - tol || y0 > r.ymax + tol || y1 < r.ymin - tol || y1 > r.ymax + tol) continue;
                let face = -1;
                if (Math.abs(y0 - r.ymax) < tol && Math.abs(y1 - r.ymax) < tol) face = 0;         // top
                else if (Math.abs(y0 - r.ymin) < tol && Math.abs(y1 - r.ymin) < tol) face = 1;    // bottom
                else if ((Math.abs(x0 - r.xmin) < tol && Math.abs(x1 - r.xmin) < tol) ||
                         (Math.abs(x0 - r.xmax) < tol && Math.abs(x1 - r.xmax) < tol)) face = 2;  // sides
                if (face < 0) continue;
                edgeFace[e] = ri * 3 + face;
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
        const z = fid < 0 ? Zbare : zForFace((fid / 3) | 0, FACES[fid % 3]);
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
    let dec;
    const abcKey = JSON.stringify(abc);
    if (cache && cache.mesh === mesh && cache.fm === fm && cache.epsMap === epsMap &&
        cache.abcKey === abcKey) {
        dec = cache.dec;
    } else {
        dec = assembleTriFEMDecomposed(mesh, fm, epsMap, abc, condRect);
        if (cache) {
            cache.mesh = mesh; cache.fm = fm; cache.epsMap = epsMap;
            cache.abcKey = abcKey; cache.dec = dec;
        }
    }
    const fem = femFromDecomposition(dec, k2);
    if (!fem) return null;
    const N = fem.N;
    // Seed the Arnoldi from the static-field projection (the quasi-TEM shape). The
    // mode-pick below also overlaps every candidate against this static drive.
    const staticSeed = staticToEdgeDofs(phiEps, mesh, fm);
    const seed = staticSeed;
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
    if (eigenSolveBytes(N, fem.csrA.colIdx.length) > MAX_SOLVE_BYTES)
        throw new Error(`Problem too large for the full-wave solver (~${(eigenSolveBytes(N, fem.csrA.colIdx.length) / 1e9).toFixed(1)} GB). Reduce Max Nodes or the domain/frequency.`);
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
        const vR = res.evecsRe.subarray(i * N, (i + 1) * N);
        let dot = 0, nS = 0, nV = 0;
        for (let k = 0; k < fm.nFreeTransverse; k++) { dot += staticSeed[k] * vR[k]; nS += staticSeed[k] ** 2; nV += vR[k] ** 2; }
        const ovl = nS > 0 && nV > 0 ? Math.abs(dot) / Math.sqrt(nS * nV) : 0;
        globalThis.__TRI_DEBUG__ && console.log(`    cand ${i}: eps=${epsMode.toFixed(4)} ovl=${ovl.toFixed(3)} (target ${eps_static.toFixed(4)})`);
        if (ovl > 0.3) cands.push({ i, eps: epsMode, ovl });   // require overlap > 0.3
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
        vRe: res.evecsRe.slice(bestIdx * N, (bestIdx + 1) * N),
        vIm: res.evecsIm.slice(bestIdx * N, (bestIdx + 1) * N),
        g2Re: res.evalsRe[bestIdx], g2Im: res.evalsIm[bestIdx],
        ambiguous, bestOvl, altOvl, bestEps, altEps,
    };
}

// ---- Dispersion cache (eps_d(f) interpolation between exact eigensolve anchors) ----
// The full-wave eps_d(f) is a smooth, slowly-varying scalar (quasi-TEM dispersion).
// On the MQS loss path the eigensolve contributes ONLY this scalar, so sweep points
// between exact anchors can interpolate it — monotone piecewise-cubic (PCHIP, no
// overshoot) over x = ln f — once a leave-one-out check proves the local
// interpolation error is within tolerance. Fritsch–Carlson slopes.
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

// Insert an exact anchor (keeps xs ascending; ignores duplicate frequencies).
function dispersionInsert(dc, f, eps) {
    const x = Math.log(f);
    let lo = 0, hi = dc.xs.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (dc.xs[m] < x) lo = m + 1; else hi = m; }
    if (lo < dc.xs.length && Math.abs(dc.xs[lo] - x) < 1e-12) return;
    dc.xs.splice(lo, 0, x); dc.eps.splice(lo, 0, eps);
}

// Interpolate eps_d at f, or null when not (yet) trustworthy. Trust requires:
// ≥5 anchors, f strictly interior, and a leave-one-out check on the interior
// anchors bracketing f: removing that anchor and predicting it from the rest
// must land within tol (relative). The LOO gap (~2 intervals) is wider than any
// interpolation gap actually bridged, so passing it is a conservative bound.
function dispersionInterp(dc, f, tol) {
    const n = dc.xs.length;
    if (n < 5) return null;
    const x = Math.log(f);
    if (x <= dc.xs[0] || x >= dc.xs[n - 1]) return null;
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (dc.xs[m] <= x) lo = m; else hi = m; }
    for (const j of new Set([Math.max(1, lo), Math.min(n - 2, hi)])) {
        const xsL = dc.xs.slice(0, j).concat(dc.xs.slice(j + 1));
        const epL = dc.eps.slice(0, j).concat(dc.eps.slice(j + 1));
        const pred = pchipEval(xsL, epL, dc.xs[j]);
        if (pred === null || Math.abs(pred - dc.eps[j]) > tol * Math.abs(dc.eps[j])) return null;
    }
    return pchipEval(dc.xs, dc.eps, x);
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
        //
        // The thickness term uses 2t, not t: the base mesh only needs ~one
        // element across the conductor thickness, the skin-depth resolution for
        // the loss solve comes from refineSkinBand, not from here. Coarser than
        // that the sweep gets slower due to skin-band refinement from
        // a too-coarse base.
        let hFine = this.opts.hFine ?? Math.min(2 * tAbs, wRef / 4) * 1.5;
        let hCoarse = this._wavelengthCap(this.opts.hCoarse ?? (dom.y_max - dom.y_min) / 5);
        // Triangle budget derived from the Max Nodes setting (memory parity with the FDM
        // backend — see maxTrisForBudget). Cap the INITIAL mesh here: element count ∝ 1/h²,
        // so if the feature-/wavelength-sized start overshoots (thin trace → tiny hFine, or a
        // high modesFreq → fine bulk), coarsen hFine/hCoarse by √(nTris/budget) and rebuild.
        // The refinement loop below also stops at this budget; capping the initial mesh too is
        // what prevents the very first eigensolve from running on an unbounded mesh.
        const maxTris = maxTrisForBudget(this.opts.maxNodes ?? 18000);
        let mesh;
        for (let attempt = 0; ; attempt++) {
            mesh = buildOccMeshFromGeometry(this.ctx.G, {
                conductors: s.conductors, dielectrics: s.dielectrics,
                domain: dom, boundaries: s.boundaries, hFine, hCoarse, symmetry: this.symmetry,
                meshConductorInterior: true,   // conductor interiors meshed for the MQS loss solve
                occSurfScale: this.opts.occSurfScale, gradeRate: this.opts.gradeRate,
                gmshOptions: this.opts.gmshOptions,
            });
            if (mesh.nTris <= maxTris || attempt >= 4) break;
            const factor = Math.sqrt(mesh.nTris / maxTris) * 1.05;   // +5% margin to land under
            hFine *= factor; hCoarse *= factor;
        }
        this.condRect = mesh.condRect;

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
            for (const rm of refModes) {
                const { abc } = modeConfig(rm, s.is_differential, this.symmetry, this.condRect.wallPEC);
                const pot = drivePotentials(this.condRect, rm, this.symmetry);
                const fm = buildTriFreedomMap(mesh, this.condRect, abc);
                const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, pot, this.ctx.wasmSolver);
                const phiAir = solveTriStatic(mesh, fm, null, pot, this.ctx.wasmSolver);
                const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
                const W_air = computeTriEnergy(phiAir, mesh, null);
                const eps_static = W_eps / W_air;
                energy.push(W_eps);
                // Full-wave eigenmode at fRef: source of the reported dispersive ε_eff
                // (fw.eps) and the H-field for the ZZ refinement metric. Computed every pass.
                let fw = null;
                try { fw = fullwaveMode(this.ctx, mesh, fm, abc, this.condRect, mesh.epsMap, fRef, phiEps, eps_static); } catch { fw = null; }
                // Converge on the REPORTED quantities, not a static proxy:
                //   • static characteristic impedance Z0 = sqrt(L/C) ∝ 1/sqrt(W_eps·W_air)
                //   • full-wave effective permittivity ε_eff = fw.eps (what _modeAtFreq reports)
                // The eigenmode ε_eff lags the static energy, so including it is what makes the
                // adaptive keep refining until the number the user sees has settled. Fall back
                // to the static C-ratio when the eigensolve fails on a pass, so the convergence
                // vector keeps a consistent length across passes.
                conv.push(1 / Math.sqrt(Math.max(W_eps * W_air, 1e-300)));
                conv.push(fw && fw.eps > 0 ? fw.eps : eps_static);
                const metricS = perElementEnergy(phiEps, mesh, mesh.epsMap);
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
            if (it === maxIters || convergedCount >= minConvergedPasses || projTris > maxTris) break;
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
        // Reject a degenerate/collapsed mesh up front (the full-wave analogue of the
        // rectilinear validate_laplace_inputs guard) rather than producing NaN/Inf in the
        // eigensolve.
        const meshErrors = validateTriMesh(mesh, this.condRect);
        if (meshErrors.length)
            throw new Error('Full-wave mesh validation failed:\n' + meshErrors.map(e => ' - ' + e).join('\n'));
        this._prepareStatic();
        return mesh;
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
        this._plotRects = (s.conductors || []).map(c =>
            ({ xmin: c.x_min, xmax: c.x_max, ymin: c.y_min, ymax: c.y_max }));
        // Asymmetric differential pair on the full domain: the odd/even excitation basis is
        // NOT the modal basis (the two traces sit in different environments), so driving
        // odd/even and picking one eigenmode per drive can return the SAME eigenmode twice
        // (e.g. broadside stripline with an unequal er stack). Use the genuine line modes
        // instead — the eigenvectors of [C0]⁻¹[C] from per-trace drives.
        if (s.is_differential && !this.symmetry && this._prepareStaticModal(grid)) return;
        for (const mode of this.modeNames) {
            const { abc, kC } = modeConfig(mode, s.is_differential, this.symmetry, this.condRect.wallPEC);
            const fm = buildTriFreedomMap(mesh, this.condRect, abc);
            const pot = drivePotentials(this.condRect, mode, this.symmetry);
            const phiEps = solveTriStatic(mesh, fm, mesh.epsMap, pot, this.ctx.wasmSolver);
            const phiAir = solveTriStatic(mesh, fm, null, pot, this.ctx.wasmSolver);
            const W_eps = computeTriEnergy(phiEps, mesh, mesh.epsMap);
            const W_air = computeTriEnergy(phiAir, mesh, null);
            const W_loss = computeTriEnergy(phiEps, mesh, mesh.lossMap);
            const C0 = kC * eps0 * W_air;                 // vacuum cap (geometric)
            const eps_eff_static = W_eps / W_air;
            const parity = this.symmetry ? (mode === 'odd' ? 'odd' : 'even') : null;
            const fields = resampleStatic(mesh, phiEps, this.domain,
                { resolution: this.opts.resolution, parity, grid, rects: this._plotRects });
            this._static[mode] = {
                fm, abc, kC, phiEps, C0, W_loss, eps_eff_static,
                Z_static: 1 / (c0 * Math.sqrt(C0 * eps_eff_static * C0)),
                fields,
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
    _prepareStaticModal(grid) {
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
            const fields = resampleStatic(mesh, phiEps, this.domain,
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
        // MQS applicability, hoisted ABOVE the eigensolve: the MQS loss path consumes
        // only the scalar eps_d from the full-wave solve (never the eigenvector), which
        // is what makes the dispersion-cache fast path below safe on that path. MQS
        // (reference) applies when grounds are wall-absorbed (signal-only rects) and
        // the domain is symmetric (mode set by the BC); else H-field perturbation.
        // Per-face plating is handled INSIDE MQS (surfaceZs weights each face's smooth
        // current by its own impedance), so plating doesn't force perturbation.
        const hasGroundRect = cr.rectRoles.some(r => !r.is_signal);
        const anyPlating = cr.rectRoles.some(r => r.plating && r.plating.sigma > 0
            && (r.plating.top || r.plating.sides || r.plating.bottom));
        // MQS drives every meshed conductor rect with the same source current, so
        // explicit ground rects (coplanar GCPW grounds etc.) would carry FORWARD
        // current instead of return current — nonsense R. Refuse the forced 'mqs'
        // override there (with a warning) rather than produce garbage.
        if (lossMethod === 'mqs' && hasGroundRect && this._modeWarnings
            && !this._modeWarnings.some(w => w.type === 'mqs-grounds')) {
            this._modeWarnings.push({ type: 'mqs-grounds', mode, freq: f,
                message: 'MQS conductor loss is not applicable with explicit ground conductors ' +
                         '(they would be driven as signal); using the perturbation method instead.' });
        }
        const useMQS = (lossMethod === 'mqs' && !hasGroundRect)
            || (lossMethod === 'auto' && this.symmetry && !hasGroundRect && cr.rects.length > 0);
        let eps_d = eps_eff_static, fw = null;
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
            const dc = useMQS ? (st.disp || (st.disp = { xs: [], eps: [] })) : null;
            const epsI = dc ? dispersionInterp(dc, f, dispTol) : null;
            if (epsI !== null) {
                eps_d = epsI;
            } else {
            let fwErr = null;
            try {
                fw = fullwaveMode(this.ctx, mesh, fm, st.abc, cr, mesh.epsMap, f, phiEps, eps_eff_static,
                    st.femCache || (st.femCache = {}));
            } catch (e) { fw = null; fwErr = e; }
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
            if (fw && fw.eps > 0) {
                eps_d = fw.eps;
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

        // (hasGroundRect / anyPlating / useMQS are computed above the eigensolve —
        // the dispersion-cache fast path there needs the MQS applicability.)
        let R_total = 0, L_internal = 0;
        // Smallest conductor cross-section dimension (gates the skin-band remesh and
        // the static/SIBC blend below).
        let minDim = Infinity;
        for (const c of cr.rects) minDim = Math.min(minDim, c.xmax - c.xmin, c.ymax - c.ymin);
        if (f > 0 && useMQS) {
            // The volume eddy solve runs the conductor BODY at the bulk metal σ;
            // plating is a SURFACE effect applied only through surfaceZs (relative to
            // the bulk Rs). So with plating, use the bulk σ here — NOT effectiveSurface's
            // plating σ (which would bake plating into the body and make surfaceZs a
            // no-op). The skin band is sized by the matching bulk skin depth.
            const mqsSigma = anyPlating ? (s.sigma_cond ?? 5.8e7) : sigma;
            const mqsDelta = anyPlating ? Math.sqrt(2 / (omu * mqsSigma)) : delta;
            // Skin-band target element size (× δ): the band must be resolved to a
            // fraction of δ for accurate R (matches ms2d's mqs_band_delta).
            // Skin-band target element size (×δ) and band width (×δ): resolve the skin
            // layer to bandDelta·δ within mqsBand·δ of each conductor surface. Tuned down
            // from (2, 0.5): a 1.5δ band with 0.6δ elements roughly halves the refined
            // mesh on tight/high-frequency geometries with ≤1% change in R (validated
            // against the reference suite), keeping the reused mesh small and the solve
            // fast.
            const bandDelta = this.opts.mqsBandDelta ?? 0.6;
            // Skin-mesh triangle budget: bounds the size-aware band refinement on
            // large-perimeter/small-δ geometries (memory guard; hitting it leaves the
            // band partially resolved rather than aborting).
            const mqsMaxTris = this.opts.mqsMaxTris ?? 40000;
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
            // Band-sizing skin depth: when a sweep announces its maximum frequency
            // (solver._sweepFmax, set by solve_sweep / InterpolatingSweep), size the
            // band for the WHOLE sweep up front — an ascending discrete sweep then
            // builds ONE skin mesh at the first point and every later point reuses
            // it (plus its cached MQS assembly), instead of re-refining per point.
            // Without the hint this reduces to the current frequency (old behavior:
            // rebuild whenever a higher frequency appears).
            const fBand = Math.max(f, this.solver._sweepFmax || 0);
            const deltaBand = fBand > f ? Math.sqrt(2 / (2 * Math.PI * fBand * MU0 * mqsSigma)) : mqsDelta;
            let mqsMesh;
            if (deltaBand >= minDim) {
                mqsMesh = mesh;
            } else if (this.opts.mqsCacheMesh === false) {
                mqsMesh = refineSkinBand(mesh, cr, mqsDelta, 12, mqsBand, bandDelta * mqsDelta, mqsMaxTris);
            } else {
                // The skin mesh depends only on the geometry and the band skin depth —
                // NOT on the mode (the odd/even BC enters the MQS solve, not the mesh) —
                // so one cache on the backend serves both modes: the band refinement
                // runs once per sweep instead of once per mode.
                const sc = (this._skinCache && this._skinCache.base === mesh) ? this._skinCache
                    : (this._skinCache = { base: mesh, dB: Infinity, mesh: null });
                if (!sc.mesh || deltaBand < sc.dB * (1 - 1e-9)) {   // higher freq appeared → rebuild finer
                    sc.dB = deltaBand;
                    sc.mesh = refineSkinBand(mesh, cr, deltaBand, 12, mqsBand, bandDelta * deltaBand, mqsMaxTris);
                }
                mqsMesh = sc.mesh;
            }
            // Per-face plating ⇒ weight each face's smooth current by its own surface
            // impedance (surfaceZs); otherwise the uniform roughness factor (Rq).
            const mqsOpts = { topGround: !!(cr.wallPEC && cr.wallPEC.top), oddSymmetry: mode === 'odd',
                              diffPair: !!s.is_differential,
                              // Frequency-invariant assembly cache (per mode; validated
                              // against the mesh object inside mqsConductorLoss, so a
                              // skin-mesh rebuild invalidates it automatically).
                              cache: st.mqsCache || (st.mqsCache = {}) };
            if (anyPlating) mqsOpts.surfaceZs = buildFaceZs(s, cr, f);
            else mqsOpts.Rq = rq;
            let mqs = null;
            try {
                mqs = mqsConductorLoss(mqsMesh, cr, f, mqsSigma, this.ctx.helpers.solveSparseMulti, 0, mqsOpts);
            } catch { mqs = null; }
            // Accept MQS only if physically sane.
            if (mqs && isFinite(mqs.R_total) && mqs.R_total > 0 && isFinite(mqs.L_loop) && mqs.L_loop > 0) {
                R_total = s.is_differential ? mqs.R_total / 2 : mqs.R_total;
                // Internal (skin) inductance from the closed form L_int = R_ac/ω, NOT the
                // difference (mqs.L_loop − L_external): that subtracts two large near-equal
                // numbers (~L_external each, computed on different bases — MQS magnetic
                // energy vs analytic 1/c²C0) so the small ~1 nH internal part is lost to
                // noise and clamps to 0 at high f, undershooting eps_eff vs the FDM/static
                // backends. The skin-regime closed form (Zs=Rs(1+j) ⇒ X_int=R_ac) is well
                // conditioned and matches the perturbation path below + the FDM backend.
                L_internal = (omega > 0) ? Math.max(0, Math.min(R_total / omega, 0.5 * L_external)) : 0;
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
        // DC resistance from the geometry lists — the SAME convention as the FDM
        // backend (field_solver.calculate_conductor_loss): current flows through the
        // signal cross-section and returns through the ground cross-section in
        // series. cr.rects can NOT be used for this: it contains coplanar ground
        // rects (which are return path, not parallel signal metal) and omits
        // wall-absorbed grounds entirely.
        const sigmaDC = s.sigma_cond ?? 5.8e7;
        let sigArea = 0, gndArea = 0;
        for (const c of s.conductors) {
            const a = Math.abs((c.width ?? (c.x_max - c.x_min)) * (c.height ?? (c.y_max - c.y_min)));
            if (c.is_signal) sigArea += a; else gndArea += a;
        }
        const R_dc = (sigArea > 0 ? 1 / (sigmaDC * sigArea) : 0)
                   + (gndArea > 0 ? 1 / (sigmaDC * gndArea) : 0);
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
            let R_ac = 0;
            for (const g of groups) {
                const fRg = RsRef > 0 ? g.Zs.re / RsRef : 1;
                const lossS = staticConductorLoss(cr.rects, f, sigmaRef, mesh, fm, phiEps,
                    Z0, eps_eff_static, eps_d, g.mask);
                let racg = lossS.R_ac * driveScale * fRg;
                if (haveFW && wFW > 0) {
                    const lossW = solveConductorLoss(cr.rects, f, sigmaRef, mesh, fm, fw.vRe, fw.vIm,
                        fw.g2Re, fw.g2Im, Pfw, Z0, mesh.epsMap, g.mask, projH, this.ctx.wasmSolver);
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
    _applyCausal(f, skipFields = false) {
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
        // Asymmetric differential pair: the odd/even drive is NOT the modal basis, so re-derive
        // the genuine line modes (and the physical [C]/[L] matrices that drive the 4-port
        // S-parameters) under the causal permittivity — exactly as the initial build does.
        // _prepareStaticModal rebuilds this._static[*] (fields, eps_eff_static, W_loss) and
        // this._modalPhys from the now-causal mesh.epsMap. It returns false for an electrically
        // symmetric or velocity-degenerate pair, in which case the odd/even drive below is the
        // correct basis (and any _modalPhys it set is already causal-consistent).
        if (s.is_differential && !this.symmetry && this._prepareStaticModal(this._staticGrid)) return;
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
            if (!skipFields) {
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
        const modes = this.modeNames.map(m => this._modeAtFreq(m, f));
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
        if (eigenSolveBytes(N, fem.csrA.colIdx.length) > MAX_SOLVE_BYTES)
            return { modes: [], nconv: 0, error: `Problem too large for the full-wave solver (~${(eigenSolveBytes(N, fem.csrA.colIdx.length) / 1e9).toFixed(1)} GB). Reduce Max Nodes or the domain/frequency.` };
        // Static-field "drive" seeds for every physical conductor excitation (single-ended:
        // one; differential: even AND odd). A genuine quasi-TEM mode matches ONE of these
        // (high overlap); a spurious discrete gradient mode matches none. All modes share
        // this fm (symmetry off), so each drive's phiEps projects into the same DOF space.
        const seeds = this.modeNames.map(m => staticToEdgeDofs(this._static[m].phiEps, mesh, fm));
        const seed = seeds[0];
        // Shift-invert centered on the quasi-TEM eigenvalue (-k²·eps_static) converges the
        // fundamental and the nearby higher-order modes first. Widen the Krylov subspace so
        // several modes resolve at once. ncvMax lets the eigensolver keep DOUBLING the
        // subspace (reusing its factorization) until nev pairs truly converge: cavity /
        // higher-order modes sit far from the quasi-TEM shift (e.g. ε_eff≈0.5 vs a shift at
        // 3.4 on an enclosed εr=4.4 line) and need a much larger subspace than the
        // fixed 20-vector pass, whose strict residual gate silently dropped them.
        const ncv = Math.min(Math.max(2 * nev + 1, 20), N - 1);
        const ncvMax = Math.min(320, N - 1);
        let res;
        try { res = this.ctx.helpers.solveGeneralized(N, fem.csrA, fem.csrB, [-k2 * eps_eff_static, 0], nev, ncv, seed, ncvMax); }
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

