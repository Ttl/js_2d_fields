// Magneto-quasi-static (MQS) volume eddy-current conductor loss.
//
// Solves the 2D skin-effect problem with the conductor interior meshed at finite σ
// — the same physics as Ansys 2D Extractor's RL solve. No Leontovich/SIBC
// assumption, no PEC-field singular boundary integral: corners are handled by the
// actual current distribution. Requires a mesh generated with
// generateGmshMesh({..., meshConductorInterior: true}).
//
// Formulation (A_z, e^{jωt}, quasi-TEM so the H pattern is ε-independent):
//   signal:     ∇²A − jωμ₀σ·A + μ₀σ·C = 0     J = σ(C − jωA), C = −dV/dz (V/m)
//   ground rect: ∇²A − jωμ₀σ·A = 0            J = −jωσA (passive return/eddy)
//   dielectric: ∇²A = 0
//   ground y=0 and outer walls: A = 0 (perfect ground; its loss is added
//   perturbatively with the flat-surface skin formula — exact for a plane).
//   Symmetry plane (condRect.symmetry > 1): natural BC at x = xmin_domain.
//
// Explicit ground rects (coplanar GCPW grounds, via-fence slabs, ground-cutout
// remnants) are return conductors tied to the reference at both line ends, so
// their per-unit-length voltage gradient is C = 0: they get finite σ (the jωμ₀σ
// mass term) but no source term, and the induced return current splits between
// them and the PEC boundary walls by the field solution itself. The roles come
// from condRect.rectRoles (is_signal), without rectRoles every rect is driven.
//
// The problem is linear in the single drive C: solve K*A1 = μ₀σ*Fc once with
// C = 1 (K = S + jωμ₀σ*M, complex symmetric, Fc spans the signal rects only),
// then scale C so the total signal current equals I. The complex system is
// solved as the real symmetric indefinite block
// [[S, −βM],[−βM, −S]]·[Ar; Ai] = [μ₀σFc; 0], symmetry is required because
// the WASM solver's LDLT fast path assumes a symmetric matrix.
//
// Per-conductor drives (opts.modeCurrents): a differential pair with no
// symmetry walls, or asymmetric pair, cannot use the single shared drive: both
// traces would act as a parallel conductor. Instead the signal rects are
// grouped by polarity (one group per net), each group gets its own drive C_k,
// and linearity gives A = Σ C_k*A_k from one unit solve per group (same
// factorization, solveSparseMulti takes all RHS at once). The small NxN current
// matrix D_jk = σ(δ_jk*area_j - jω*Fc_j*A_k) maps drives to net currents,
// solving D*C = I_target realizes the requested mode currents (+1/-1 odd, +1/+1
// even, or the modal current vector for an asymmetric pair).

import { tripletsToCSRMulti, GL3p, GL3w } from './fem_core.js';
import { triCoefficients, lv, le, lvGrad, leGrad, QW, QL1, QL2, QL3, NQ,
         triP2Stiffness, P2_MASS, P2_LOAD, refineTriMesh } from './tri_fem.js';
import { calculate_Zrough } from '../surface_roughness.js';

const MU0 = 4 * Math.PI * 1e-7;
const edgeVerts = [[0,1],[1,2],[2,0]];

function inAnyRect(rects, x, y, tol) {
    for (const r of rects) {
        if (x > r.xmin - tol && x < r.xmax + tol && y > r.ymin - tol && y < r.ymax + tol) return true;
    }
    return false;
}

// Refine triangles near any conductor surface (inside or outside, within `band`
// skin depths) up to `passes` times. A triangle is marked when any of its
// vertices (or its centroid) lies within the band — vertex-based marking makes
// the refinement engage even when elements are much larger than the band
// (coarse initial meshes: triangles touching the surface always qualify).
// Element size roughly halves per pass. If `targetH` is given, refinement stops
// early once the smallest conductor-surface edge is ≤ targetH.
//
// `grading` (optional) relaxes the target size with distance from the signal
// conductors: target(x,y) = max(targetH, slope*(d_sig - Dfine)) where d_sig is
// the distance to the nearest rect in grading.sigRects. Used for the ground-rect
// band (GCPW coplanar grounds, via slabs): their surface current decays away
// from the slot, so full δ-resolution is only needed within ~Dfine of the signal
// and the band cost stays independent of how far the ground pour extends.
// Refinement stops naturally once the graded target exceeds the base mesh size.
export function refineSkinBand(mesh, condRect, delta, passes, band = 3, targetH = 0, maxTris = Infinity, grading = null, depthSlope = 0) {
    const rects = condRect.rects || [condRect];
    const bw = band * delta;
    function distToRectBoundary(r, x, y) {
        const dx = Math.max(r.xmin - x, 0, x - r.xmax);
        const dy = Math.max(r.ymin - y, 0, y - r.ymax);
        const outside = Math.hypot(dx, dy);
        if (outside > 0) return outside;
        return -Math.min(x - r.xmin, r.xmax - x, y - r.ymin, r.ymax - y);
    }
    function nearSurface(x, y) {
        for (const r of rects) if (Math.abs(distToRectBoundary(r, x, y)) < bw) return true;
        return false;
    }
    // Graded target size at a point: outside distance to the nearest signal rect.
    function targetAt(x, y) {
        if (!grading) return targetH;
        let d = Infinity;
        for (const r of grading.sigRects) {
            const dx = Math.max(r.xmin - x, 0, x - r.xmax);
            const dy = Math.max(r.ymin - y, 0, y - r.ymax);
            d = Math.min(d, Math.hypot(dx, dy));
        }
        return Math.max(targetH, grading.slope * (d - grading.Dfine));
    }
    // Depth grading (across the band, not along the surface). The eddy current
    // decays as exp(-d/δ), so the outer part of the
    // band carries a small, smooth fraction of the total and does not need the
    // surface's element size. Isotropic bisection means an element's tangential 
    // extent shrinks with its target too, so a layer at 2*targetH costs 4x fewer
    // triangles.
    //     target(d) = targetH*(1 + depthSlope*|d|/δ)
    // depthSlope = 0 for uniform meshing.
    const depthTarget = (dSurf) => targetH * (1 + depthSlope * dSurf / delta);
    // Size-aware marking: refine a band triangle only while it is still larger
    // than targetH, and stop when nothing needs refining. (The previous early-stop
    // keyed on the MINIMUM conductor-surface edge, which the main mesh's
    // corner-concentrated adaptive refinement satisfies passes early — aborting
    // while the face centers were still ≫ δ — and conversely kept refining
    // already-fine corner elements on every pass.)
    // maxTris bounds the growth: a large conductor perimeter at a tiny skin depth
    // would otherwise refine without limit (band-tris ∝ perimeter/δ) and exhaust
    // memory. Hitting the cap degrades gracefully to a partially-resolved band —
    // the same accuracy compromise the old early-stop made.
    //
    // "Gracefully" still costs accuracy. So record a band that stopped before
    // everything reached targetH, by the triangle cap or by running out of
    // passes, on the returned mesh as `bandTrunc`, for the caller to surface as
    // a warning. null = the band converged (a final marking scan found nothing
    // left to refine).
    let trunc = null;
    for (let p = 0; ; p++) {
        const marked = new Uint8Array(mesh.nTris);
        let any = false, nMarked = 0;
        for (let t = 0; t < mesh.nTris; t++) {
            const v0 = mesh.tris[3*t], v1 = mesh.tris[3*t+1], v2 = mesh.tris[3*t+2];
            const x0 = mesh.nodes[2*v0], y0 = mesh.nodes[2*v0+1];
            const x1 = mesh.nodes[2*v1], y1 = mesh.nodes[2*v1+1];
            const x2 = mesh.nodes[2*v2], y2 = mesh.nodes[2*v2+1];
            const xc = (x0+x1+x2)/3, yc = (y0+y1+y2)/3;
            if (!(nearSurface(xc, yc) || nearSurface(x0, y0) ||
                  nearSurface(x1, y1) || nearSurface(x2, y2))) continue;
            let tgt = grading ? targetAt(xc, yc) : targetH;
            if (depthSlope > 0 && targetH > 0) {
                // Distance of the element's closest point (of centroid + vertices) to
                // metal. A triangle straddling the surface must get the surface target.
                let dS = Infinity;
                for (const r of rects) {
                    dS = Math.min(dS, Math.abs(distToRectBoundary(r, xc, yc)),
                                      Math.abs(distToRectBoundary(r, x0, y0)),
                                      Math.abs(distToRectBoundary(r, x1, y1)),
                                      Math.abs(distToRectBoundary(r, x2, y2)));
                }
                tgt = Math.max(tgt, depthTarget(dS));
            }
            if (tgt > 0) {
                const hMax = Math.max(Math.hypot(x1-x0, y1-y0),
                                      Math.hypot(x2-x1, y2-y1),
                                      Math.hypot(x0-x2, y0-y2));
                if (hMax <= tgt) continue;   // this element already resolves the band
            }
            marked[t] = 1; any = true; nMarked++;
        }
        if (!any) break;                     // target reached everywhere
        // The pass counter is checked here rather than in the loop header so the
        // scan above runs once more after the last refinement: that is what tells
        // converged-on-the-final-pass apart from out-of-passes. The number of
        // refinements is unchanged.
        if (p >= passes) { trunc = { reason: 'passes', passes, nTris: mesh.nTris, nMarked }; break; }
        // Refining splits each marked triangle ~1→4 (plus conformity closure);
        // stop before the projected size blows the budget.
        if (mesh.nTris + 3 * nMarked > maxTris) {
            trunc = { reason: 'maxTris', maxTris, nTris: mesh.nTris, nMarked };
            break;
        }
        mesh = refineTriMesh(mesh, marked);
    }
    mesh.bandTrunc = trunc;
    return mesh;
}

// Compute conductor loss by volume eddy-current solve.
//   mesh        — mesh WITH conductor interior triangles
//   condRect    — conductor geometry object ({rects: [...]} for multi-rect signal)
//   freq, sigma — frequency (Hz), conductivity (S/m)
//   solveSparseMulti — WASM helper from createWasmHelpers
//   Z0          — (optional) line impedance for α_c conversion
//   opts.wallPEC — {left,right,top,bottom} flags saying which domain walls are metal
//   (from the mesher. `left` is already cleared on a half domain so the symmetry
//   plane is never counted). Only these walls contribute to the wall loss, every
//   wall is Dirichlet in the solve regardless. Omit it and the legacy behaviour
//   applies: bottom always, top per opts.topGround, sides never.
//   opts.topGround — legacy fallback for opts.wallPEC.top (stripline)
//   opts.oddSymmetry — odd-mode symmetry: A = 0 (Dirichlet) on the symmetry plane
//   x = xmin_domain instead of the natural (even-mode) BC. For a centered
//   differential pair meshed as a half domain with one full trace, this solves
//   the odd mode; the default natural BC solves the even mode.
//   opts.modeCurrents — per-conductor-drive path (full domain only): array of
//   target net currents, one per signal polarity group (positive-polarity group
//   first, pre.groups order). Use [1, -1] for the odd mode, [1, 1] for even, or
//   the modal current vector (normalized to Σ|I|^2 = 2) for an asymmetric pair.
//   Replaces the symmetry-plane mode selection: broadside/asymmetric pairs and
//   symmetry-disabled solves get a per-mode MQS answer instead of the
//   perturbation fallback. R_total/X_total keep the differential convention
//   (both traces, per-line mode R = R_total/2), L_loop is per-line.
//   opts.Rq — RMS surface roughness (m), gradient model (single layer, same Rq on
//   all surfaces). The dissipation is scaled by Ψ_R(f) = Re(Z_rough)/Rs and the
//   loop inductance gets the matching surface-reactance increment
//   ΔL = (Im(Z_rough) − Rs)/ω · ∮|K|² = R_smooth·(Im(Z_rough)/Rs − 1)/ω, keeping
//   R(f)/L(f) causal. With a single Rq the per-surface factor is uniform, so
//   scaling the smooth-σ solution is exact in the skin regime (t ≫ δ) and
//   degrades gracefully to the correct DC limit (Ψ → 1 as δ grows).
//   opts.surfaceZs(x, y, orient) — OPTIONAL per-face surface impedance {re,im} at a
//   face midpoint (orient 'h' = top/bottom, 'v' = side), for per-side plating. The
//   smooth trace (and ground) loss is scaled by the |K|²-weighted average of
//   Re(Zs_face)/Rs over the faces — generalizing the uniform Ψ_R. The weights come
//   from the MQS smooth current distribution (corner-regularized, unlike a
//   perturbation surface integral). Reduces exactly to the uniform factor when
//   every face shares one Zs. Takes precedence over Rq when provided.
// Returns { R_trace, R_gnd, R_total, X_total, L_loop, alpha_c, alpha_c_dBm, delta, nDofs }
// L_loop is the series inductance from Im(Z_pul) = ωL — includes trace internal
// inductance and ground-current spreading (ground itself is PEC here).
// X_total is the surface REACTANCE per unit length, the twin of R_total: the same
// |K|²-weighted surface integral taken against Im(Zs) where R_total takes Re(Zs). It is
// what the caller builds the internal inductance from (X_total/ω), because L_loop −
// L_external cancels two large near-equal numbers and loses it. Smooth metal has
// Zs = Rs(1 + j), so X_total === R_total there and only rough/plated surfaces separate
// them (Im(Zs)/Re(Zs) reaches ~5 at 1 µm rms). Same differential convention as R_total.
// Frequency-INVARIANT part of the MQS solve on a given mesh: conductor
// classification, DOF map, the S (everywhere) / M, Fc (conductor) assembly, and
// the symbolic block-CSR with separate S / M value templates. The block system
// is [[S, −βM], [−βM, −S]] with β = ωμ₀σ the ONLY frequency-dependent piece, so
// per frequency the values are just valS + β·valM on a fixed pattern. Cached by
// the caller (opts.cache) per (mesh, oddSymmetry) — the skin mesh is reused
// across sweep points (f_max-reuse), and re-assembling it every point used to
// dominate the MQS path's JS time.
export function mqsPrecompute(mesh, condRect, opts = {}) {
    const { nodes, edges, tris, triEdges, nNodes, nEdges, nTris } = mesh;
    const rects = condRect.rects || [condRect];
    const roles = condRect.rectRoles || null;
    const sym = condRect.symmetry > 1 ? 2 : 1;
    const TOL = 1e-12;

    // Conductor triangles by centroid: 1 = signal (driven, C), 2 = ground rect
    // (passive, C = 0). Without rectRoles every rect is driven.
    // Signal wins where rects overlap; ground rects may overlap each other (GCPW
    // via slab under the coplanar ground), same class either way.
    const sigRects = [], sigPol = [];
    rects.forEach((r, i) => {
        if (!roles || roles[i].is_signal) { sigRects.push(r); sigPol.push(roles ? (roles[i].polarity || 0) : 0); }
    });
    const gndRects = roles ? rects.filter((_, i) => !roles[i].is_signal) : [];
    // Drive groups: one per distinct signal polarity (positive first). A
    // single-ended line has one group; a differential pair two (+/-). Used by
    // the per-conductor-drive path (opts.modeCurrents), the single-drive path
    // ignores the grouping entirely.
    const groups = [...new Set(sigPol)].sort((a, b) => b - a);
    const groupOfSig = sigPol.map(p => groups.indexOf(p));
    const isCondTri = new Uint8Array(nTris);
    const triGroup = new Int32Array(nTris).fill(-1);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const xc = (nodes[2*v0]+nodes[2*v1]+nodes[2*v2])/3;
        const yc = (nodes[2*v0+1]+nodes[2*v1+1]+nodes[2*v2+1])/3;
        // Same containment test as inAnyRect, but keeping the matching rect's
        // index so the triangle lands in its polarity group (first match wins,
        // identical iteration order = identical classification).
        let si = -1;
        for (let i = 0; i < sigRects.length; i++) {
            const r = sigRects[i];
            if (xc > r.xmin - TOL && xc < r.xmax + TOL && yc > r.ymin - TOL && yc < r.ymax + TOL) { si = i; break; }
        }
        if (si >= 0) { isCondTri[t] = 1; triGroup[t] = groupOfSig[si]; }
        else if (gndRects.length && inAnyRect(gndRects, xc, yc, TOL)) isCondTri[t] = 2;
    }

    // DOFs: vertices + edge midpoints, Dirichlet at ground/outer walls.
    // The symmetry plane (x = xmin_domain when symmetry) gets the natural BC.
    // The bottom ground wall is the DOMAIN bottom (ymin_domain — usually 0 after
    // the ground slab is wall-absorbed, but not by construction).
    const xmin_d = condRect.xmin_domain, xmax_d = condRect.xmax_domain;
    const ymax_d = condRect.ymax_domain;
    const ymin_d = condRect.ymin_domain ?? 0;
    function isDirichletPt(x, y) {
        if (Math.abs(y - ymin_d) < 1e-9 || Math.abs(y - ymax_d) < 1e-9) return true;
        if (Math.abs(x - xmax_d) < 1e-9) return true;
        if ((sym === 1 || opts.oddSymmetry) && Math.abs(x - xmin_d) < 1e-9) return true;
        return false;
    }
    const dofOf = new Int32Array(nNodes + nEdges).fill(-1);
    let nF = 0;
    for (let n = 0; n < nNodes; n++) {
        if (!isDirichletPt(nodes[2*n], nodes[2*n+1])) dofOf[n] = nF++;
    }
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const xm = (nodes[2*n0]+nodes[2*n1])/2, ym = (nodes[2*n0+1]+nodes[2*n1+1])/2;
        if (!isDirichletPt(xm, ym)) dofOf[nNodes + e] = nF++;
    }

    // Assemble S (everywhere), M (all metal: signal + passive grounds) and Fc
    // (signal only, grounds have C = 0, so no source term). condArea is the
    // signal cross-section: it is only used to normalize the signal current.
    // The three element integrals are closed-form on a straight-sided triangle.
    // The P2 stiffness comes from triP2Stiffness and the mass / load are Area
    // x the reference P2_MASS / P2_LOAD. That replaces a 6-point quadrature
    // loop with six basis evaluations per point, the assembly used to dominate
    // the MQS path's JS time on the (skin-refined, hence large) conductor mesh.
    const Nb = 2 * nF;
    const maxCoo = 4 * 36 * nTris;   // ≤ 2 block entries per S and per M local entry
    const R = new Int32Array(maxCoo), Cc = new Int32Array(maxCoo);
    const VS = new Float64Array(maxCoo), VM = new Float64Array(maxCoo);
    let nCoo = 0;
    const Fc = new Float64Array(nF);
    // Per-group source vectors / areas for the per-conductor-drive path. The
    // combined Fc/condArea keep their own accumulation (not a sum of these) so
    // the single-drive path's floating-point stream is unchanged.
    const FcG = groups.map(() => new Float64Array(nF));
    const areaG = new Float64Array(groups.length);
    let condArea = 0;
    const lg = new Int32Array(6);
    const Sl = new Float64Array(36);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const Area = triP2Stiffness(nodes, v0, v1, v2, Sl);
        lg[0] = dofOf[v0]; lg[1] = dofOf[v1]; lg[2] = dofOf[v2];
        for (let k = 0; k < 3; k++) lg[3+k] = dofOf[nNodes + triEdges[3*t+k]];
        const cond = isCondTri[t];
        const driven = cond === 1;
        if (driven) { condArea += Area; areaG[triGroup[t]] += Area; }
        const FcGt = driven ? FcG[triGroup[t]] : null;
        for (let i = 0; i < 6; i++) {
            const gi = lg[i]; if (gi < 0) continue;
            if (driven) { const f = Area * P2_LOAD[i]; Fc[gi] += f; FcGt[gi] += f; }
            for (let j = 0; j < 6; j++) {
                const gj = lg[j]; if (gj < 0) continue;
                // Block system [[S, -beta * M], [-beta * M, -S]] with beta
                // = w * mu0 * sigma, the only frequency-dependent factor. 
                // Put each local entry into both of its blocks with separate
                // S / M values on one shared index stream.
                const sv = Sl[6*i+j];
                if (sv !== 0) {
                    R[nCoo] = gi; Cc[nCoo] = gj; VS[nCoo] = sv; VM[nCoo] = 0; nCoo++;
                    R[nCoo] = nF + gi; Cc[nCoo] = nF + gj; VS[nCoo] = -sv; VM[nCoo] = 0; nCoo++;
                }
                const mv = cond ? Area * P2_MASS[6*i+j] : 0;
                if (mv !== 0) {
                    R[nCoo] = gi; Cc[nCoo] = nF + gj; VS[nCoo] = 0; VM[nCoo] = -mv; nCoo++;
                    R[nCoo] = nF + gi; Cc[nCoo] = gj; VS[nCoo] = 0; VM[nCoo] = -mv; nCoo++;
                }
            }
        }
    }

    // Symbolic block system with separate S / M value templates on one pattern:
    // valS and valM line up entry-for-entry, so a frequency only needs
    // valS + beta * valM.
    const csr2 = tripletsToCSRMulti(R.subarray(0, nCoo), Cc.subarray(0, nCoo), Nb,
                                    [VS.subarray(0, nCoo), VM.subarray(0, nCoo)]);
    const csrS = { rowPtr: csr2.rowPtr, colIdx: csr2.colIdx, valRe: csr2.vals[0],
                   valIm: new Float64Array(csr2.colIdx.length) };
    const csrM = { valRe: csr2.vals[1] };

    // edge → one adjacent triangle (for the ground-loss surface integral)
    const edgeToTri = new Int32Array(nEdges).fill(-1);
    for (let t = 0; t < nTris; t++)
        for (let k = 0; k < 3; k++)
            if (edgeToTri[triEdges[3*t+k]] === -1) edgeToTri[triEdges[3*t+k]] = t;

    return {
        isCondTri, dofOf, nF, Nb, Fc, condArea, edgeToTri,
        groups, triGroup, FcG, areaG,
        rowPtr: csrS.rowPtr, colIdx: csrS.colIdx,
        valS: csrS.valRe, valM: csrM.valRe, valIm: csrS.valIm,   // valIm stays zero
    };
}

// Dense complex NxN solve (Gauss-Jordan, partial pivot) for the drive-to-current
// matrix, N is the number of driven nets (2 for a differential pair).
// A: rows of {re, im}; b: array of {re, im}. Returns array of {re, im}.
function solveComplexLinear(A, b) {
    const n = b.length;
    const M = A.map((row, i) => row.map(z => [z.re, z.im]).concat([[b[i].re, b[i].im]]));
    for (let c = 0; c < n; c++) {
        let p = c, best = M[c][c][0] ** 2 + M[c][c][1] ** 2;
        for (let r = c + 1; r < n; r++) {
            const m = M[r][c][0] ** 2 + M[r][c][1] ** 2;
            if (m > best) { best = m; p = r; }
        }
        if (!(best > 0)) throw new Error('mqs: singular drive-current matrix');
        const tmp = M[c]; M[c] = M[p]; M[p] = tmp;
        const ar = M[c][c][0], ai = M[c][c][1], den = ar * ar + ai * ai;
        for (let r = 0; r < n; r++) {
            if (r === c) continue;
            const br = M[r][c][0], bi = M[r][c][1];
            const fr = (br * ar + bi * ai) / den, fi = (bi * ar - br * ai) / den;
            if (fr === 0 && fi === 0) continue;
            for (let k = c; k <= n; k++) {
                const mr = M[c][k][0], mi = M[c][k][1];
                M[r][k][0] -= fr * mr - fi * mi;
                M[r][k][1] -= fr * mi + fi * mr;
            }
        }
    }
    return M.map((row, i) => {
        const ar = row[i][0], ai = row[i][1], den = ar * ar + ai * ai;
        return { re: (row[n][0] * ar + row[n][1] * ai) / den,
                 im: (row[n][1] * ar - row[n][0] * ai) / den };
    });
}

export function mqsConductorLoss(mesh, condRect, freq, sigma, solveSparseMulti, Z0 = 0, opts = {}) {
    const { nodes, edges, tris, triEdges, nNodes, nEdges, nTris } = mesh;
    const rects = condRect.rects || [condRect];
    const sym = condRect.symmetry > 1 ? 2 : 1;
    const omega = 2 * Math.PI * freq;
    const delta = Math.sqrt(2 / (omega * MU0 * sigma));
    const Rs = 1 / (sigma * delta);
    const TOL = 1e-12;
    const xmin_d = condRect.xmin_domain;
    const xmax_d = condRect.xmax_domain;
    const ymax_d = condRect.ymax_domain;
    const ymin_d = condRect.ymin_domain ?? 0;

    // Frequency-invariant assembly: reuse the caller's cache when it matches
    // this exact mesh (identity) and symmetry mode; else recompute (and store).
    const cc = opts.cache;
    let pre = (cc && cc.mesh === mesh && cc.odd === !!opts.oddSymmetry) ? cc.pre : null;
    if (!pre) {
        pre = mqsPrecompute(mesh, condRect, opts);
        if (cc) { cc.mesh = mesh; cc.odd = !!opts.oddSymmetry; cc.pre = pre; }
    }
    const { isCondTri, dofOf, nF, Nb, Fc, condArea, edgeToTri, triGroup } = pre;
    const lg = new Int32Array(6);

    // Per-frequency system: values = valS + β·valM on the cached pattern.
    const beta = omega * MU0 * sigma;
    const val = new Float64Array(pre.valS.length);
    for (let k = 0; k < val.length; k++) val[k] = pre.valS[k] + beta * pre.valM[k];
    const csr = { rowPtr: pre.rowPtr, colIdx: pre.colIdx, valRe: val, valIm: pre.valIm };

    // Per-conductor-drive path (opts.modeCurrents, full domain only): one unit
    // solve per polarity group, then the NxN current-matrix solve for the drive
    // constants that realize the target net currents.
    const multiI = (opts.modeCurrents && sym === 1 && pre.groups.length > 1
        && opts.modeCurrents.length === pre.groups.length) ? opts.modeCurrents : null;
    if (opts.modeCurrents && !multiI) throw new Error('mqs: modeCurrents needs a full-domain mesh with one signal group per entry');
    let sol, Cr, Ci, Zmode = null, CgR = null, CgI = null;
    if (multiI) {
        const nG = pre.groups.length;
        // Unit solutions A_k are mode-independent, cache them per (mesh, beta) so
        // the second mode of the pair at the same frequency skips the block LU.
        let sols = (cc && cc.sols && cc.solsMesh === mesh && cc.solsBeta === beta) ? cc.sols : null;
        if (!sols) {
            const rhsList = pre.FcG.map(F => {
                const r = new Float64Array(Nb);
                for (let i = 0; i < nF; i++) r[i] = MU0 * sigma * F[i];
                return r;
            });
            sols = solveSparseMulti(Nb, csr, rhsList);
            if (cc) { cc.sols = sols; cc.solsMesh = mesh; cc.solsBeta = beta; }
        }
        // Net current in group j for unit drive on group k:
        // D_jk = σ(δ_jk·area_j − jω·Fc_j·A_k)
        const D = [];
        for (let j = 0; j < nG; j++) {
            const row = [];
            const F = pre.FcG[j];
            for (let k = 0; k < nG; k++) {
                const A = sols[k];
                let fr = 0, fi = 0;
                for (let i = 0; i < nF; i++) { fr += F[i] * A[i]; fi += F[i] * A[nF + i]; }
                row.push({ re: sigma * ((j === k ? pre.areaG[j] : 0) + omega * fi),
                           im: -sigma * omega * fr });
            }
            D.push(row);
        }
        const Cg = solveComplexLinear(D, multiI.map(v => ({ re: v, im: 0 })));
        // Combined physical field A = Σ C_k*A_k, downstream integrals then run
        // with unit scale (Cr = 1) and the per-group complex drive in CgR/CgI.
        const comb = new Float64Array(Nb);
        for (let k = 0; k < nG; k++) {
            const A = sols[k], cRk = Cg[k].re, cIk = Cg[k].im;
            for (let i = 0; i < nF; i++) {
                comb[i] += cRk * A[i] - cIk * A[nF + i];
                comb[nF + i] += cRk * A[nF + i] + cIk * A[i];
            }
        }
        sol = comb;
        CgR = Cg.map(z => z.re); CgI = Cg.map(z => z.im);
        Cr = 1; Ci = 0;
        // Per-line mode impedance Z = Σ C_k*conj(I_k) / Σ|I_k|^2 (C = −dV/dz): the
        // multi-net twin of the single-drive Z = C/I. Re(Z) is the per-line mode
        // R (power balance), Im(Z)/ω the per-line loop L.
        let numR = 0, numI = 0, den = 0;
        for (let k = 0; k < nG; k++) {
            numR += Cg[k].re * multiI[k]; numI += Cg[k].im * multiI[k];
            den += multiI[k] * multiI[k];
        }
        Zmode = { re: numR / den, im: numI / den };
    } else {
        const rhs = new Float64Array(Nb);
        for (let i = 0; i < nF; i++) rhs[i] = MU0 * sigma * Fc[i];
        [sol] = solveSparseMulti(Nb, csr, [rhs]);

        // Rescale C for trace current 1 A. When the signal conductor straddles the
        // symmetry plane (single-ended half mesh), the meshed part carries half the
        // current. When it lies entirely inside the half domain (differential pair,
        // one full trace meshed), it carries the full per-trace current. Ground rects
        // are excluded. They carry return current, not the normalized drive current.
        const rolesCR = condRect.rectRoles || null;
        const straddles = rects.some((r, i) =>
            (!rolesCR || rolesCR[i].is_signal) && r.xmin <= xmin_d + 1e-12);
        const I_mesh = (sym === 2 && straddles) ? 0.5 : 1;
        let fr = 0, fi = 0;
        for (let i = 0; i < nF; i++) { fr += Fc[i] * sol[i]; fi += Fc[i] * sol[nF + i]; }
        const dR = sigma * (condArea + omega * fi);
        const dI = -sigma * omega * fr;
        const dMag2 = dR*dR + dI*dI;
        Cr = I_mesh * dR / dMag2; Ci = -I_mesh * dI / dMag2;
    }
    const Cmag2 = Cr*Cr + Ci*Ci;

    // Conductor dissipation: J/σ = C*(u - jωA1) with the drive u = 1 in the
    // signal (class 1) and u = 0 in passive ground rects (class 2, pure eddy /
    // return current). Signal and ground-rect dissipation accumulate separately
    // so the ground share reports (and plating-scales) with R_gnd, not R_trace.
    let Psig = 0, PgndRect = 0;
    for (let t = 0; t < nTris; t++) {
        const cls = isCondTri[t];
        if (!cls) continue;
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
        lg[0] = dofOf[v0]; lg[1] = dofOf[v1]; lg[2] = dofOf[v2];
        for (let k = 0; k < 3; k++) lg[3+k] = dofOf[nNodes + triEdges[3*t+k]];
        const xs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
        const ys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];
        // Drive term: 1 (single) or the group's complex C_k (multi) in the
        // signal, 0 in passive ground rects (pure eddy / return current).
        let dvR = 0, dvI = 0;
        if (cls === 1) { if (CgR) { const g = triGroup[t]; dvR = CgR[g]; dvI = CgI[g]; } else dvR = 1; }
        let Ptri = 0;
        for (let q = 0; q < NQ; q++) {
            const w = QW[q] * Area;
            const xq = xs[0]*QL1[q] + xs[1]*QL2[q] + xs[2]*QL3[q];
            const yq = ys[0]*QL1[q] + ys[1]*QL2[q] + ys[2]*QL3[q];
            let aR = 0, aI = 0;
            for (let k = 0; k < 6; k++) {
                const g = lg[k]; if (g < 0) continue;
                const Nk = k < 3 ? lv(coeff, k, xq, yq) : le(coeff, edgeVerts[k-3][0], edgeVerts[k-3][1], xq, yq);
                aR += Nk * sol[g]; aI += Nk * sol[nF + g];
            }
            const uR = dvR + omega * aI, uI = dvI - omega * aR;
            const eR = Cr * uR - Ci * uI, eI = Cr * uI + Ci * uR;
            Ptri += 0.5 * sigma * (eR*eR + eI*eI) * w;
        }
        if (cls === 1) Psig += Ptri; else PgndRect += Ptri;
    }

    // Domain-wall loss
    //
    // Flat-surface skin formula on the tangential H at each wall
    // that is actually metal. |Hx| = |∂A/∂y|/μ₀ on a horizontal wall, |Hy| =
    // |∂A/∂x|/μ₀ on a vertical one. (edgeToTri comes precomputed from mqsPrecompute.)
    //
    // Every domain wall is Dirichlet in the MQS solve (A = 0 confines the flux), but
    // only the metal ones dissipate. An 'open' far-field truncation and the symmetry
    // plane are boundary conditions, not surfaces. opts.wallPEC carries that
    // distinction from the mesher (which clears `left` on a half domain, so the
    // symmetry plane can never be mistaken for metal).
    let Pgnd = 0;
    // Per-face plating weights (∮|K|²dl and Σ Re/Im(Zs)·|K|²dl) over the ground and
    // conductor surfaces, used below to scale the smooth loss per face.
    let gndS = 0, gndZreS = 0, gndZimS = 0;
    // Dallback for a caller with no wall map: bottom is ground on every
    // geometry the mesher builds, top came from opts.topGround, sides skipped.
    const wp = opts.wallPEC || { bottom: true, top: !!opts.topGround };
    const onLine = (a, b, v) => Math.abs(a - v) < 1e-9 && Math.abs(b - v) < 1e-9;
    // Orientation of the metal wall this edge lies on ('h' | 'v'), else null.
    function wallOrient(x0, y0, x1, y1) {
        if (wp.bottom && onLine(y0, y1, ymin_d)) return 'h';
        if (wp.top && onLine(y0, y1, ymax_d)) return 'h';
        if (wp.left && onLine(x0, x1, xmin_d)) return 'v';
        if (wp.right && onLine(x0, x1, xmax_d)) return 'v';
        return null;
    }
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const x0 = nodes[2*n0], y0 = nodes[2*n0+1];
        const x1 = nodes[2*n1], y1 = nodes[2*n1+1];
        const orient = wallOrient(x0, y0, x1, y1);
        if (!orient) continue;
        const adj = edgeToTri[e];
        // A meshed conductor sitting on the wall (a GCPW pour or via slab reaching
        // it) bonds to it. There is no exposed wall surface there, and the adjacent
        // triangle is metal, so ∂A/∂n is an interior gradient, not the exterior H.
        // That segment's dissipation belongs to the volume eddy term, not here.
        if (adj < 0 || isCondTri[adj]) continue;
        const v0 = tris[3*adj], v1 = tris[3*adj+1], v2 = tris[3*adj+2];
        const { coeff } = triCoefficients(nodes, v0, v1, v2);
        lg[0] = dofOf[v0]; lg[1] = dofOf[v1]; lg[2] = dofOf[v2];
        for (let k = 0; k < 3; k++) lg[3+k] = dofOf[nNodes + triEdges[3*adj+k]];
        const L = Math.hypot(x1 - x0, y1 - y0);
        // Tangential H comes from the gradient component along the wall normal.
        const comp = orient === 'h' ? 1 : 0;
        // Walls are bare metal (never plated). Evaluate Zs once at the edge midpoint.
        const Zg = opts.surfaceZs ? opts.surfaceZs((x0 + x1) / 2, (y0 + y1) / 2, orient) : null;
        for (let q = 0; q < 3; q++) {
            const xq = x0 + GL3p[q] * (x1 - x0), yq = y0 + GL3p[q] * (y1 - y0);
            let gR = 0, gI = 0;
            for (let k = 0; k < 6; k++) {
                const g = lg[k]; if (g < 0) continue;
                const gr = k < 3 ? lvGrad(coeff, k, xq, yq) : leGrad(coeff, edgeVerts[k-3][0], edgeVerts[k-3][1], xq, yq);
                gR += gr[comp] * sol[g]; gI += gr[comp] * sol[nF + g];
            }
            const K2 = Cmag2 * (gR*gR + gI*gI) / (MU0*MU0);
            Pgnd += 0.5 * Rs * K2 * GL3w[q] * L;
            if (Zg) {
                const Sseg = (gR*gR + gI*gI) * GL3w[q] * L;   // ∝ |K|² (global factors cancel in the ratio)
                gndS += Sseg; gndZreS += Zg.re * Sseg; gndZimS += Zg.im * Sseg;
            }
        }
    }

    // Conductor-surface weights: the trace loss is a volume integral, so to apply a
    // per-face impedance we weight each face by its surface current ∮|K|²dl, taken
    // on the EXTERIOR (dielectric) side of the face — like the ground above. K is
    // the tangential H: ∂A/∂y for a horizontal (top/bottom) face, ∂A/∂x for a side.
    // Signal faces and ground-rect faces accumulate separate buckets: each scales
    // its own volume loss (a plated trace next to a bare ground must not dilute).
    let trS = 0, trZreS = 0, trZimS = 0;
    let grS = 0, grZreS = 0, grZimS = 0;
    if (opts.surfaceZs) {
        const eA = new Int32Array(2 * nEdges).fill(-1);
        for (let t = 0; t < nTris; t++) for (let k = 0; k < 3; k++) {
            const e = triEdges[3*t+k];
            if (eA[2*e] < 0) eA[2*e] = t; else eA[2*e+1] = t;
        }
        const lge = new Int32Array(6);
        for (let e = 0; e < nEdges; e++) {
            const ta = eA[2*e], tb = eA[2*e+1];
            if (ta < 0 || tb < 0) continue;          // boundary edge, not a cond/dielectric interface
            const aMetal = isCondTri[ta] > 0, bMetal = isCondTri[tb] > 0;
            if (aMetal === bMetal) continue;         // both metal or both dielectric → not a surface
            const ext = aMetal ? tb : ta;            // exterior (dielectric) triangle
            const cnd = aMetal ? ta : tb;            // conductor-interior triangle
            const n0 = edges[2*e], n1 = edges[2*e+1];
            const x0 = nodes[2*n0], y0 = nodes[2*n0+1], x1 = nodes[2*n1], y1 = nodes[2*n1+1];
            if (Math.abs(x0 - xmin_d) < 1e-9 && Math.abs(x1 - xmin_d) < 1e-9) continue;  // symmetry-plane cut
            const horiz = Math.abs(y1 - y0) < Math.abs(x1 - x0);
            // Classify against the conductor this edge borders, then SNAP the query
            // point onto that face. The skin-band refinement makes the centroid-based
            // interface wander ~an element size off the nominal face, so a raw midpoint
            // can miss the face's exact-coordinate test (esp. the short side faces) and
            // be mis-read as bare. Snapping uses the conductor-side triangle's rect.
            let qx = (x0 + x1) / 2, qy = (y0 + y1) / 2;
            const ccx = (nodes[2*tris[3*cnd]] + nodes[2*tris[3*cnd+1]] + nodes[2*tris[3*cnd+2]]) / 3;
            const ccy = (nodes[2*tris[3*cnd]+1] + nodes[2*tris[3*cnd+1]+1] + nodes[2*tris[3*cnd+2]+1]) / 3;
            for (const r of rects) {
                if (ccx <= r.xmin - TOL || ccx >= r.xmax + TOL || ccy <= r.ymin - TOL || ccy >= r.ymax + TOL) continue;
                if (horiz) qy = Math.abs(qy - r.ymax) < Math.abs(qy - r.ymin) ? r.ymax : r.ymin;
                else qx = Math.abs(qx - r.xmax) < Math.abs(qx - r.xmin) ? r.xmax : r.xmin;
                break;
            }
            const Zs = opts.surfaceZs(qx, qy, horiz ? 'h' : 'v');
            const v0 = tris[3*ext], v1 = tris[3*ext+1], v2 = tris[3*ext+2];
            const { coeff } = triCoefficients(nodes, v0, v1, v2);
            lge[0] = dofOf[v0]; lge[1] = dofOf[v1]; lge[2] = dofOf[v2];
            for (let k = 0; k < 3; k++) lge[3+k] = dofOf[nNodes + triEdges[3*ext+k]];
            const L = Math.hypot(x1 - x0, y1 - y0);
            let Sseg = 0;
            for (let q = 0; q < 3; q++) {
                const xq = x0 + GL3p[q]*(x1-x0), yq = y0 + GL3p[q]*(y1-y0);
                let gR = 0, gI = 0;
                for (let k = 0; k < 6; k++) {
                    const g = lge[k]; if (g < 0) continue;
                    const gr = k < 3 ? lvGrad(coeff, k, xq, yq) : leGrad(coeff, edgeVerts[k-3][0], edgeVerts[k-3][1], xq, yq);
                    const comp = horiz ? gr[1] : gr[0];   // tangential-H-producing gradient component
                    gR += comp * sol[g]; gI += comp * sol[nF + g];
                }
                Sseg += (gR*gR + gI*gI) * GL3w[q] * L;
            }
            if (isCondTri[cnd] === 1) { trS += Sseg; trZreS += Zs.re * Sseg; trZimS += Zs.im * Sseg; }
            else { grS += Sseg; grZreS += Zs.re * Sseg; grZimS += Zs.im * Sseg; }
        }
    }

    // Totals: R = 2P/|I|². Single-ended (straddling): R_total is the line R for
    // |I| = 1 A. Differential half mesh (full trace meshed, per-trace |I| = 1):
    // R_total covers BOTH traces (mirror included) — per-trace mode R is
    // R_total/2, and L_loop is already per-trace (from the drive field C/I₁).
    // Cross-check: Re(C) = per-trace dissipation (power balance of the drive
    // field: signal + ground-rect volume loss. The PEC boundary walls add
    // nothing there, their loss is the perturbative Pgnd term).
    // R_gnd has two parts: the meshed passive ground rects (volume eddy loss)
    // and the PEC boundary walls (flat-surface skin formula).
    let R_trace = 2 * sym * Psig;
    let R_gr = 2 * sym * PgndRect;
    let R_gw = 2 * sym * Pgnd;
    // Z_pul = C/I = R + jωL (trace internal L included). On the multi-drive path
    // the per-line mode impedance Zmode plays the role of C (per-line current
    // normalized to the target vector), same convention.
    let L_loop = (Zmode ? Zmode.im : Ci) / omega;

    // Surface roughness / plating post-processing. Scale the smooth-σ loss by the
    // effective surface impedance and add the matching surface-reactance increment
    // to L_loop (keeping R(f)/L(f) causal). Either:
    //   • per-face (opts.surfaceZs): the |K|²-weighted Re/Im(Zs)/Rs over each face's
    //     surface current — trace and ground weighted independently; or
    //   • uniform (opts.Rq): a single Ψ_R = Re(Z_rough)/Rs for the whole surface.
    let PsiR = 1;
    const Rq = opts.Rq || 0;
    // For a differential pair (opts.diffPair) R_trace/R_gnd cover BOTH traces
    // (the caller halves R_total to the per-trace mode R), while L_loop is
    // per-trace — so the surface-reactance increment must use the per-trace R.
    const perTrace = opts.diffPair ? 0.5 : 1;
    // X_* are the REACTANCE twins of the R components: the same smooth-metal
    // loss weighted by Im(Zs)/Rs where R is weighted by Re(Zs)/Rs. They are what the
    // internal inductance is built from (X_total/omega), and they carry the same units
    // and the same differential convention as R, so the caller halves them alike.
    // Smooth metal has Im(Zs) = Re(Zs) = Rs, so these stay equal to R, the seed
    // value below is that smooth case, and each branch overwrites it from the
    // SMOOTH R before R is scaled by psiR. Each component scales by its own
    // face bucket: signal faces -> R_trace, ground-rect faces -> R_gr, boundary
    // walls -> R_gw (walls are bare metal unless surfaceZs says otherwise).
    let X_trace = R_trace, X_gr = R_gr, X_gw = R_gw;
    const R_smooth_total = R_trace + R_gr + R_gw;   // before psiR rewrites any part
    if (opts.surfaceZs) {
        if (trS > 0) {
            const psiR = trZreS / (Rs * trS), psiX = trZimS / (Rs * trS);
            X_trace = R_trace * psiX;
            R_trace *= psiR;
            PsiR = psiR;
        }
        if (grS > 0) {
            const psiR = grZreS / (Rs * grS), psiX = grZimS / (Rs * grS);
            X_gr = R_gr * psiX;
            R_gr *= psiR;
        }
        if (gndS > 0) {
            const psiR = gndZreS / (Rs * gndS), psiX = gndZimS / (Rs * gndS);
            X_gw = R_gw * psiX;
            R_gw *= psiR;
        }
    } else if (Rq > 0) {
        const Zs = calculate_Zrough(freq, sigma, Rq);
        PsiR = Zs.re / Rs;
        const PsiX = Zs.im / Rs;
        X_trace = R_trace * PsiX;
        X_gr = R_gr * PsiX;
        X_gw = R_gw * PsiX;
        R_trace *= PsiR;
        R_gr *= PsiR;
        R_gw *= PsiR;
    }

    const R_gnd = R_gr + R_gw;
    const X_gnd = X_gr + X_gw;
    const R_total = R_trace + R_gnd;
    const X_total = X_trace + X_gnd;
    // The surface-reactance increment on the loop inductance is (X − R_smooth)/ω by
    // construction, so it is formed ONCE from X_total here rather than accumulated
    // branch by branch. That keeps L_loop and X_total from ever describing different
    // surfaces, and it is identically zero for smooth metal (X_total === R_smooth_total).
    L_loop += perTrace * (X_total - R_smooth_total) / omega;
    const alpha_c = Z0 > 0 ? R_total / (2 * Z0) : NaN;
    return {
        R_trace, R_gnd, R_total, X_total, L_loop, PsiR,
        alpha_c, alpha_c_dBm: alpha_c * 8.686,
        Rs, delta, nDofs: Nb, modeZ: Zmode,
    };
}
