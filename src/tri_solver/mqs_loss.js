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

import { tripletsToCSR } from './fem_core.js';
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
export function refineSkinBand(mesh, condRect, delta, passes, band = 3, targetH = 0) {
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
    function onSurface(x, y) {
        for (const r of rects) {
            const T = 1e-12;
            const onX = (Math.abs(x-r.xmin)<T || Math.abs(x-r.xmax)<T) && y>r.ymin-T && y<r.ymax+T;
            const onY = (Math.abs(y-r.ymin)<T || Math.abs(y-r.ymax)<T) && x>r.xmin-T && x<r.xmax+T;
            if (onX || onY) return true;
        }
        return false;
    }
    function minSurfaceEdge() {
        let h = Infinity;
        for (let e = 0; e < mesh.nEdges; e++) {
            const n0 = mesh.edges[2*e], n1 = mesh.edges[2*e+1];
            if (!onSurface(mesh.nodes[2*n0], mesh.nodes[2*n0+1])) continue;
            if (!onSurface(mesh.nodes[2*n1], mesh.nodes[2*n1+1])) continue;
            const L = Math.hypot(mesh.nodes[2*n1]-mesh.nodes[2*n0], mesh.nodes[2*n1+1]-mesh.nodes[2*n0+1]);
            if (L < h) h = L;
        }
        return h;
    }
    for (let p = 0; p < passes; p++) {
        if (targetH > 0 && minSurfaceEdge() <= targetH) break;
        const marked = new Uint8Array(mesh.nTris);
        for (let t = 0; t < mesh.nTris; t++) {
            const v0 = mesh.tris[3*t], v1 = mesh.tris[3*t+1], v2 = mesh.tris[3*t+2];
            const xc = (mesh.nodes[2*v0]+mesh.nodes[2*v1]+mesh.nodes[2*v2])/3;
            const yc = (mesh.nodes[2*v0+1]+mesh.nodes[2*v1+1]+mesh.nodes[2*v2+1])/3;
            if (nearSurface(xc, yc) ||
                nearSurface(mesh.nodes[2*v0], mesh.nodes[2*v0+1]) ||
                nearSurface(mesh.nodes[2*v1], mesh.nodes[2*v1+1]) ||
                nearSurface(mesh.nodes[2*v2], mesh.nodes[2*v2+1])) marked[t] = 1;
        }
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
// Returns { R_trace, R_gnd, R_total, L_loop, alpha_c, alpha_c_dBm, delta, nDofs }
// L_loop is the series inductance from Im(Z_pul) = ωL — includes trace internal
// inductance and ground-current spreading (ground itself is PEC here).
export function mqsConductorLoss(mesh, condRect, freq, sigma, solveSparseMulti, Z0 = 0, opts = {}) {
    const { nodes, edges, tris, triEdges, nNodes, nEdges, nTris } = mesh;
    const rects = condRect.rects || [condRect];
    const sym = condRect.symmetry > 1 ? 2 : 1;
    const omega = 2 * Math.PI * freq;
    const delta = Math.sqrt(2 / (omega * MU0 * sigma));
    const Rs = 1 / (sigma * delta);
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
    const xmin_d = condRect.xmin_domain, xmax_d = condRect.xmax_domain;
    const ymax_d = condRect.ymax_domain;
    function isDirichletPt(x, y) {
        if (Math.abs(y) < 1e-9 || Math.abs(y - ymax_d) < 1e-9) return true;
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

    // Symmetric real block system, C=1 solve
    const beta = omega * MU0 * sigma;
    const Nb = 2 * nF;
    const R = [], Cc = [], V = [];
    for (let k = 0; k < sR.length; k++) {
        R.push(sR[k]); Cc.push(sC[k]); V.push(sV[k]);
        R.push(nF + sR[k]); Cc.push(nF + sC[k]); V.push(-sV[k]);
    }
    for (let k = 0; k < mR.length; k++) {
        R.push(mR[k]); Cc.push(nF + mC[k]); V.push(-beta * mV[k]);
        R.push(nF + mR[k]); Cc.push(mC[k]); V.push(-beta * mV[k]);
    }
    const rhs = new Float64Array(Nb);
    for (let i = 0; i < nF; i++) rhs[i] = MU0 * sigma * Fc[i];
    const csr = tripletsToCSR(R, Cc, V, Nb);
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
    const edgeToTri = new Array(nEdges);
    for (let e = 0; e < nEdges; e++) edgeToTri[e] = -1;
    for (let t = 0; t < nTris; t++)
        for (let k = 0; k < 3; k++)
            if (edgeToTri[triEdges[3*t+k]] === -1) edgeToTri[triEdges[3*t+k]] = t;
    let Pgnd = 0;
    const GL3p = [0.11270, 0.5, 0.88730], GL3w = [0.27778, 0.44444, 0.27778];
    function isGroundY(y) {
        if (Math.abs(y) < 1e-9) return true;
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

    // Surface roughness (gradient model, single Rq on all surfaces)
    const Rq = opts.Rq || 0;
    let PsiR = 1;
    if (Rq > 0) {
        const Zs = calculate_Zrough(freq, sigma, Rq);
        PsiR = Zs.re / Rs;
        const R_smooth = R_trace + R_gnd;
        L_loop += R_smooth * (Zs.im / Rs - 1) / omega;
        R_trace *= PsiR;
        R_gnd *= PsiR;
    }

    const R_total = R_trace + R_gnd;
    const alpha_c = Z0 > 0 ? R_total / (2 * Z0) : NaN;
    return {
        R_trace, R_gnd, R_total, L_loop, PsiR,
        alpha_c, alpha_c_dBm: alpha_c * 8.686,
        Rs, delta, nDofs: Nb,
    };
}
