// Magneto-quasi-static (MQS) volume eddy-current conductor loss.
//
// Solves the 2D skin-effect problem with the conductor interior meshed at finite σ
// — the same physics as Ansys 2D Extractor's RL solve. No Leontovich/SIBC
// assumption, no PEC-field singular boundary integral: corners are handled by the
// actual current distribution. Requires a mesh generated with
// generateGmshMesh({..., meshConductorInterior: true}).
//
// Formulation (A_z, e^{jωt}, quasi-TEM so the H pattern is ε-independent):
//   conductor:  ∇²A − jωμ₀σ·A + μ₀σ·C = 0     J = σ(C − jωA), C = −dV/dz (V/m)
//   dielectric: ∇²A = 0
//   ground y=0 and outer walls: A = 0 (perfect ground; its loss is added
//   perturbatively with the flat-surface skin formula — exact for a plane).
//   Symmetry plane (condRect.symmetry > 1): natural BC at x = xmin_domain.
//
// The problem is linear in C: solve K·A1 = μ₀σ·Fc once with C = 1
// (K = S + jωμ₀σ·Mc, complex symmetric), then scale C so the total conductor
// current equals I. The complex system is solved as the REAL SYMMETRIC indefinite
// block [[S, −βM],[−βM, −S]]·[Ar; Ai] = [μ₀σFc; 0] — symmetry is required because
// the WASM solver's LDLT fast path assumes a symmetric matrix.

import { tripletsToCSR, GL3p, GL3w } from './fem_core.js';
import { triCoefficients, lv, le, lvGrad, leGrad, QW, QL1, QL2, QL3, NQ,
         refineTriMesh } from './tri_fem.js';
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
export function refineSkinBand(mesh, condRect, delta, passes, band = 3, targetH = 0, maxTris = Infinity) {
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
    for (let p = 0; p < passes; p++) {
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
            if (targetH > 0) {
                const hMax = Math.max(Math.hypot(x1-x0, y1-y0),
                                      Math.hypot(x2-x1, y2-y1),
                                      Math.hypot(x0-x2, y0-y2));
                if (hMax <= targetH) continue;   // this element already resolves the band
            }
            marked[t] = 1; any = true; nMarked++;
        }
        if (!any) break;
        // Refining splits each marked triangle ~1→4 (plus conformity closure);
        // stop before the projected size blows the budget.
        if (mesh.nTris + 3 * nMarked > maxTris) break;
        mesh = refineTriMesh(mesh, marked);
    }
    return mesh;
}

// Compute conductor loss by volume eddy-current solve.
//   mesh        — mesh WITH conductor interior triangles
//   condRect    — conductor geometry object ({rects: [...]} for multi-rect signal)
//   freq, sigma — frequency (Hz), conductivity (S/m)
//   solveSparseMulti — WASM helper from createWasmHelpers
//   Z0          — (optional) line impedance for α_c conversion
//   opts.topGround — include the y=ymax_domain wall in the ground loss (stripline)
//   opts.oddSymmetry — odd-mode symmetry: A = 0 (Dirichlet) on the symmetry plane
//   x = xmin_domain instead of the natural (even-mode) BC. For a centered
//   differential pair meshed as a half domain with one full trace, this solves
//   the odd mode; the default natural BC solves the even mode.
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
    const sym = condRect.symmetry > 1 ? 2 : 1;
    const TOL = 1e-12;

    // Conductor triangles by centroid
    const isCondTri = new Uint8Array(nTris);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const xc = (nodes[2*v0]+nodes[2*v1]+nodes[2*v2])/3;
        const yc = (nodes[2*v0+1]+nodes[2*v1+1]+nodes[2*v2+1])/3;
        if (inAnyRect(rects, xc, yc, TOL)) isCondTri[t] = 1;
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

    // Assemble S (everywhere), Mc and Fc (conductor only)
    const sR = [], sC = [], sV = [], mR = [], mC = [], mV = [];
    const Fc = new Float64Array(nF);
    let condArea = 0;
    const lg = new Int32Array(6);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
        lg[0] = dofOf[v0]; lg[1] = dofOf[v1]; lg[2] = dofOf[v2];
        for (let k = 0; k < 3; k++) lg[3+k] = dofOf[nNodes + triEdges[3*t+k]];
        const xs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
        const ys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];
        const cond = isCondTri[t];
        if (cond) condArea += Area;
        const Sl = new Float64Array(36), Ml = new Float64Array(36), Fl = new Float64Array(6);
        for (let q = 0; q < NQ; q++) {
            const w = QW[q] * Area;
            const xq = xs[0]*QL1[q] + xs[1]*QL2[q] + xs[2]*QL3[q];
            const yq = ys[0]*QL1[q] + ys[1]*QL2[q] + ys[2]*QL3[q];
            const N = new Float64Array(6), Gx = new Float64Array(6), Gy = new Float64Array(6);
            for (let k = 0; k < 3; k++) {
                N[k] = lv(coeff, k, xq, yq);
                const g = lvGrad(coeff, k, xq, yq); Gx[k] = g[0]; Gy[k] = g[1];
                const [p, qq] = edgeVerts[k];
                N[3+k] = le(coeff, p, qq, xq, yq);
                const ge = leGrad(coeff, p, qq, xq, yq); Gx[3+k] = ge[0]; Gy[3+k] = ge[1];
            }
            for (let i = 0; i < 6; i++) {
                if (cond) Fl[i] += w * N[i];
                for (let j = 0; j < 6; j++) {
                    Sl[6*i+j] += w * (Gx[i]*Gx[j] + Gy[i]*Gy[j]);
                    if (cond) Ml[6*i+j] += w * N[i]*N[j];
                }
            }
        }
        for (let i = 0; i < 6; i++) {
            const gi = lg[i]; if (gi < 0) continue;
            if (cond) Fc[gi] += Fl[i];
            for (let j = 0; j < 6; j++) {
                const gj = lg[j]; if (gj < 0) continue;
                if (Sl[6*i+j] !== 0) { sR.push(gi); sC.push(gj); sV.push(Sl[6*i+j]); }
                if (cond && Ml[6*i+j] !== 0) { mR.push(gi); mC.push(gj); mV.push(Ml[6*i+j]); }
            }
        }
    }

    // Symbolic block system with separate S / M value templates. The SAME
    // (row, col) sequence is converted twice, so the dedup produces an identical
    // pattern for both — valS and valM line up entry-for-entry.
    const Nb = 2 * nF;
    const R = [], Cc = [], VS = [], VM = [];
    for (let k = 0; k < sR.length; k++) {
        R.push(sR[k]); Cc.push(sC[k]); VS.push(sV[k]); VM.push(0);
        R.push(nF + sR[k]); Cc.push(nF + sC[k]); VS.push(-sV[k]); VM.push(0);
    }
    for (let k = 0; k < mR.length; k++) {
        R.push(mR[k]); Cc.push(nF + mC[k]); VS.push(0); VM.push(-mV[k]);
        R.push(nF + mR[k]); Cc.push(mC[k]); VS.push(0); VM.push(-mV[k]);
    }
    const csrS = tripletsToCSR(R, Cc, VS, Nb);
    const csrM = tripletsToCSR(R, Cc, VM, Nb);

    // edge → one adjacent triangle (for the ground-loss surface integral)
    const edgeToTri = new Int32Array(nEdges).fill(-1);
    for (let t = 0; t < nTris; t++)
        for (let k = 0; k < 3; k++)
            if (edgeToTri[triEdges[3*t+k]] === -1) edgeToTri[triEdges[3*t+k]] = t;

    return {
        isCondTri, dofOf, nF, Nb, Fc, condArea, edgeToTri,
        rowPtr: csrS.rowPtr, colIdx: csrS.colIdx,
        valS: csrS.valRe, valM: csrM.valRe, valIm: csrS.valIm,   // valIm stays zero
    };
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
    const { isCondTri, dofOf, nF, Nb, Fc, condArea, edgeToTri } = pre;
    const lg = new Int32Array(6);

    // Per-frequency system: values = valS + β·valM on the cached pattern.
    const beta = omega * MU0 * sigma;
    const val = new Float64Array(pre.valS.length);
    for (let k = 0; k < val.length; k++) val[k] = pre.valS[k] + beta * pre.valM[k];
    const csr = { rowPtr: pre.rowPtr, colIdx: pre.colIdx, valRe: val, valIm: pre.valIm };
    const rhs = new Float64Array(Nb);
    for (let i = 0; i < nF; i++) rhs[i] = MU0 * sigma * Fc[i];
    const [sol] = solveSparseMulti(Nb, csr, [rhs]);

    // Rescale C for trace current 1 A. When the conductor straddles the symmetry
    // plane (single-ended half mesh), the meshed part carries half the current;
    // when the conductor lies entirely inside the half domain (differential pair,
    // one full trace meshed), it carries the full per-trace current.
    const straddles = rects.some(r => r.xmin <= xmin_d + 1e-12);
    const I_mesh = (sym === 2 && straddles) ? 0.5 : 1;
    let fr = 0, fi = 0;
    for (let i = 0; i < nF; i++) { fr += Fc[i] * sol[i]; fi += Fc[i] * sol[nF + i]; }
    const dR = sigma * (condArea + omega * fi);
    const dI = -sigma * omega * fr;
    const dMag2 = dR*dR + dI*dI;
    const Cr = I_mesh * dR / dMag2, Ci = -I_mesh * dI / dMag2;
    const Cmag2 = Cr*Cr + Ci*Ci;

    // Conductor dissipation: J/σ = C·(1 − jωA1)
    let Pcond = 0;
    for (let t = 0; t < nTris; t++) {
        if (!isCondTri[t]) continue;
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
        lg[0] = dofOf[v0]; lg[1] = dofOf[v1]; lg[2] = dofOf[v2];
        for (let k = 0; k < 3; k++) lg[3+k] = dofOf[nNodes + triEdges[3*t+k]];
        const xs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
        const ys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];
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
            const uR = 1 + omega * aI, uI = -omega * aR;
            const eR = Cr * uR - Ci * uI, eI = Cr * uI + Ci * uR;
            Pcond += 0.5 * sigma * (eR*eR + eI*eI) * w;
        }
    }

    // Ground loss: flat-surface skin formula on |Hx(y=0)| = |∂A/∂y|/μ₀
    // (edgeToTri comes precomputed from mqsPrecompute)
    let Pgnd = 0;
    // Per-face plating weights (∮|K|²dl and Σ Re/Im(Zs)·|K|²dl) over the ground and
    // conductor surfaces, used below to scale the smooth loss per face.
    let gndS = 0, gndZreS = 0, gndZimS = 0;
    function isGroundY(y) {
        if (Math.abs(y - ymin_d) < 1e-9) return true;
        if (opts.topGround && Math.abs(y - ymax_d) < 1e-9) return true;
        return false;
    }
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        if (!isGroundY(nodes[2*n0+1]) || !isGroundY(nodes[2*n1+1])) continue;
        const adj = edgeToTri[e];
        if (adj < 0) continue;
        const v0 = tris[3*adj], v1 = tris[3*adj+1], v2 = tris[3*adj+2];
        const { coeff } = triCoefficients(nodes, v0, v1, v2);
        lg[0] = dofOf[v0]; lg[1] = dofOf[v1]; lg[2] = dofOf[v2];
        for (let k = 0; k < 3; k++) lg[3+k] = dofOf[nNodes + triEdges[3*adj+k]];
        const x0 = nodes[2*n0], x1 = nodes[2*n1];
        const yWall = nodes[2*n0+1];
        const L = Math.abs(x1 - x0);
        // Ground is bare metal (never plated); evaluate Zs once at the edge midpoint.
        const Zg = opts.surfaceZs ? opts.surfaceZs((x0 + x1) / 2, yWall, 'h') : null;
        for (let q = 0; q < 3; q++) {
            const xq = x0 + GL3p[q] * (x1 - x0);
            let gyR = 0, gyI = 0;
            for (let k = 0; k < 6; k++) {
                const g = lg[k]; if (g < 0) continue;
                const gr = k < 3 ? lvGrad(coeff, k, xq, yWall) : leGrad(coeff, edgeVerts[k-3][0], edgeVerts[k-3][1], xq, yWall);
                gyR += gr[1] * sol[g]; gyI += gr[1] * sol[nF + g];
            }
            const K2 = Cmag2 * (gyR*gyR + gyI*gyI) / (MU0*MU0);
            Pgnd += 0.5 * Rs * K2 * GL3w[q] * L;
            if (Zg) {
                const Sseg = (gyR*gyR + gyI*gyI) * GL3w[q] * L;   // ∝ |K|² (global factors cancel in the ratio)
                gndS += Sseg; gndZreS += Zg.re * Sseg; gndZimS += Zg.im * Sseg;
            }
        }
    }

    // Conductor-surface weights: the trace loss is a volume integral, so to apply a
    // per-face impedance we weight each face by its surface current ∮|K|²dl, taken
    // on the EXTERIOR (dielectric) side of the face — like the ground above. K is
    // the tangential H: ∂A/∂y for a horizontal (top/bottom) face, ∂A/∂x for a side.
    let trS = 0, trZreS = 0, trZimS = 0;
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
            if (isCondTri[ta] === isCondTri[tb]) continue;   // both in or both out → not a surface
            const ext = isCondTri[ta] ? tb : ta;     // exterior (dielectric) triangle
            const cnd = isCondTri[ta] ? ta : tb;      // conductor-interior triangle
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
            trS += Sseg; trZreS += Zs.re * Sseg; trZimS += Zs.im * Sseg;
        }
    }

    // Totals: R = 2P/|I|². Single-ended (straddling): R_total is the line R for
    // |I| = 1 A. Differential half mesh (full trace meshed, per-trace |I| = 1):
    // R_total covers BOTH traces (mirror included) — per-trace mode R is
    // R_total/2, and L_loop is already per-trace (from the drive field C/I₁).
    // Cross-check: Re(C) = per-trace dissipation (power balance of the drive
    // field; ground is PEC inside the solve so it adds nothing there).
    let R_trace = 2 * sym * Pcond;
    let R_gnd = 2 * sym * Pgnd;
    let L_loop = Ci / omega;  // Z_pul = C/I = R + jωL (trace internal L included)

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
    // X_trace/X_gnd are the REACTANCE twins of R_trace/R_gnd: the same smooth-metal
    // loss weighted by Im(Zs)/Rs where R is weighted by Re(Zs)/Rs. They are what the
    // internal inductance is built from (X_total/omega), and they carry the same units
    // and the same differential convention as R, so the caller halves them alike.
    // Smooth metal has Im(Zs) = Re(Zs) = Rs, so these stay equal to R, the seed
    // value below is that smooth case, and each branch overwrites it from the
    // SMOOTH R before R is scaled by psiR.
    let X_trace = R_trace, X_gnd = R_gnd;
    const R_smooth_total = R_trace + R_gnd;   // before psiR rewrites either one
    if (opts.surfaceZs) {
        if (trS > 0) {
            const psiR = trZreS / (Rs * trS), psiX = trZimS / (Rs * trS);
            X_trace = R_trace * psiX;
            R_trace *= psiR;
            PsiR = psiR;
        }
        if (gndS > 0) {
            const psiR = gndZreS / (Rs * gndS), psiX = gndZimS / (Rs * gndS);
            X_gnd = R_gnd * psiX;
            R_gnd *= psiR;
        }
    } else if (Rq > 0) {
        const Zs = calculate_Zrough(freq, sigma, Rq);
        PsiR = Zs.re / Rs;
        const PsiX = Zs.im / Rs;
        X_trace = R_trace * PsiX;
        X_gnd = R_gnd * PsiX;
        R_trace *= PsiR;
        R_gnd *= PsiR;
    }

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
        Rs, delta, nDofs: Nb,
    };
}
