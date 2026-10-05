// Magneto-quasi-static (MQS) eddy-current conductor loss.
//
// Solves the 2D skin-effect problem for A_z with the conductor interiors meshed at
// finite sigma, so the current distribution (corners included) comes out of the
// solve rather than from a surface-impedance model. The mesh must include the
// conductor interiors (buildOccMeshFromGeometry, meshConductorInterior on).
//
// Formulation (e^{jwt}, quasi-TEM so the H pattern does not depend on eps):
//   signal rect:  lap(A) - jw*mu0*sigma*A + mu0*sigma*C = 0   J = sigma*(C - jw*A)
//   ground rect:  lap(A) - jw*mu0*sigma*A = 0                 J = -jw*sigma*A
//   dielectric:   lap(A) = 0
// C = -dV/dz is the drive. Ground rects (coplanar grounds, via slabs) are tied to
// the reference at both line ends, so their C is 0 and the return current they
// carry is induced. The roles come from condRect.rectRoles, without them every
// rect is driven.
//
// Boundaries: A = 0 on the metal walls (opts.wallPEC) and on the symmetry plane of
// an odd mode. Open walls and the even-mode symmetry plane keep the natural BC.
// Without a wall map, or with neither metal walls nor ground rects, every outer
// wall is A = 0. Ideal grounds (opts.idealRects) are held at A = 0 too. The loss of
// the walls and ideal grounds is a surface term from the tangential H, with the
// slab impedance of the wall thickness (opts.wallThick).
//
// The system is linear in C: one solve of K*A1 = mu0*sigma*Fc with C = 1, then C is
// scaled so the signal current equals I. K = S + jw*mu0*sigma*M is complex symmetric
// (K^T = K, not Hermitian) and is factorized by the WASM LDL^T
// (solveComplexSymmetric), which falls back to a pivoted LU when its residual check
// fails. Vectors are stored as [Ar; Ai], the N real parts then the N imaginary parts.
//
// A differential pair without a symmetry plane (and any asymmetric pair) needs one
// drive per net, opts.modeCurrents: the signal rects are grouped by polarity, each
// group gets a unit solve on the same factorization, and the small current matrix
// D_jk = sigma*(delta_jk*area_j - jw*Fc_j*A_k) gives the drives C_k for the target
// net currents (+1/-1 odd, +1/+1 even, or the modal currents).
//
// Roughness and plating scale the smooth loss per face by the |K|^2-weighted
// Re(Zs)/Rs (opts.Rq, opts.surfaceZs, opts.surfaceDZ). mqsPecInductance gives the
// sigma -> infinity inductance on the same mesh, the reference the caller subtracts
// to get the internal inductance.

import { tripletsToCSRMulti, GL3p, GL3w } from './fem_core.js';
import { triCoefficients, lvGrad, leGrad, QW, QL1, QL2, QL3, NQ,
         triP2Stiffness, P2_MASS, P2_LOAD, refineTriMesh } from './tri_fem.js';
import { calculate_Zrough, wallSpreadFactor, slabCoth } from '../surface_roughness.js';
import { shapeContains, shapeSignedDist, shapeFaceAt, insideRingHole } from '../shapes.js';

const MU0 = 4 * Math.PI * 1e-7;
// Thin dimension of a conductor entry: the smaller side of a rect, the thickness a
// shape declares (a ring's wall), else the smaller side of its bounding box.
export function condThinDim(r) {
    if (r.shape && r.shape.thickness > 0) return r.shape.thickness;
    return Math.min(r.xmax - r.xmin, r.ymax - r.ymin);
}
const edgeVerts = [[0,1],[1,2],[2,0]];
// P2 shape functions at the quadrature points, P2_AT_Q[6*q + k]: vertices 2λ^2 - λ,
// edges (edgeVerts) 4 λ_p λ_q. Constant on every straight-sided triangle.
const P2_AT_Q = (() => {
    const out = new Float64Array(6 * NQ);
    for (let q = 0; q < NQ; q++) {
        const l = [QL1[q], QL2[q], QL3[q]];
        for (let k = 0; k < 3; k++) out[6*q + k] = 2*l[k]*l[k] - l[k];
        for (let k = 0; k < 3; k++) out[6*q + 3 + k] = 4*l[edgeVerts[k][0]]*l[edgeVerts[k][1]];
    }
    return out;
})();

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
//
// condRect.symX (optional) is the x of a symmetry plane. A rect cut by it is measured
// as its mirrored whole, so the cut face gets no band. Shaped conductors keep their
// full polygon in a half-domain solve and need no mirroring.
//
// `aniso` (optional, { cornerGrade, maxAspect, cornerTurn }) refines to a metric aligned with the
// surface instead of isotropically: the size across the band is the depth-graded
// target above, the size along it only shrinks towards corners, max(across,
// cornerGrade * distance to the nearest corner), at most maxAspect times the size
// across. The current along a face varies on the scale of the face, not of delta, so
// the band elements stretch along it. Corners are the rect corners and the polygon
// vertices turning by more than cornerTurn (CORNER_TURN), the vertices of an n-gon standing in for
// a circle are not. Each triangle takes one metric, at its point closest to the
// surface, and bisects its longest edge in it: longest-edge bisection in a fixed
// affine frame keeps the triangles' shape in the metric bounded. (Edges measured each
// in their own metric let a triangle keep splitting a short edge where the normal
// turns, down to a degenerate needle.)
const CORNER_TURN = 20 * Math.PI / 180;
export function refineSkinBand(mesh, condRect, delta, passes, band = 3, targetH = 0, maxTris = Infinity, grading = null, depthSlope = 0, aniso = null) {
    const symX = condRect.symX ?? null;
    const rects = (condRect.rects || [condRect]).map(r =>
        symX !== null && !r.shape && Math.abs(r.xmin - symX) < 1e-12 ? { ...r, xmin: 2 * symX - r.xmax } : r);
    const bw = band * delta;
    // A ground ring around every signal conductor is a shield: its current runs on the
    // hole side, and its outer surface needs no band.
    const shields = new Set(grading ? rects.filter(r => r.shape && r.shape.type === 'ring'
        && grading.sigRects.every(sr => insideRingHole(r.shape, sr))) : []);
    function distToRectBoundary(r, x, y) {
        if (shields.has(r)) return -shapeSignedDist({ type: 'polygon', poly: r.shape.hole }, x, y);
        // A polygon's signed distance is exact inside and never too large outside.
        if (r.shape) return shapeSignedDist(r.shape, x, y);
        const dx = Math.max(r.xmin - x, 0, x - r.xmax);
        const dy = Math.max(r.ymin - y, 0, y - r.ymax);
        const outside = Math.hypot(dx, dy);
        if (outside > 0) return outside;
        return -Math.min(x - r.xmin, r.xmax - x, y - r.ymin, r.ymax - y);
    }
    // Graded target size at a point: outside distance to the nearest signal rect.
    function targetAt(x, y) {
        if (!grading) return targetH;
        let d = Infinity;
        for (const r of grading.sigRects) {
            if (r.shape) { d = Math.min(d, Math.max(0, shapeSignedDist(r.shape, x, y))); continue; }
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
    // Per-rect boxes grown by the band width: a triangle whose bounding box misses
    // every one of them has no point within bw of any surface (a point within bw of a
    // rect's boundary lies inside the grown box), so the far field, most of the mesh
    // on every pass, is rejected on four comparisons per rect instead of the
    // distance evaluations below.
    const grown = rects.map(r => ({ x0: r.xmin - bw, x1: r.xmax + bw, y0: r.ymin - bw, y1: r.ymax + bw }));
    // Distance of a point to the nearest conductor surface.
    function surfDist(x, y) {
        let d = Infinity;
        for (const r of rects) d = Math.min(d, Math.abs(distToRectBoundary(r, x, y)));
        return d;
    }
    const metric = aniso && targetH > 0 ? bandMetric() : null;
    for (let p = 0; ; p++) {
        if (metric) {
            const r = metric.pass(mesh);
            if (!r) break;
            if (p >= passes) { trunc = { reason: 'passes', passes, nTris: mesh.nTris, nMarked: r.nMarked }; break; }
            if (mesh.nTris + 2 * r.nMarked > maxTris) {
                trunc = { reason: 'maxTris', maxTris, nTris: mesh.nTris, nMarked: r.nMarked };
                break;
            }
            mesh = refineTriMesh(mesh, r.marked, { longest: r.longest, smooth: false });
            continue;
        }
        const marked = new Uint8Array(mesh.nTris);
        let any = false, nMarked = 0;
        // Vertex distances, evaluated once per node on first use (NaN = not yet).
        const nodeDist = new Float64Array(mesh.nNodes).fill(NaN);
        const vertDist = (v) => {
            let d = nodeDist[v];
            if (d !== d) d = nodeDist[v] = surfDist(mesh.nodes[2*v], mesh.nodes[2*v+1]);
            return d;
        };
        for (let t = 0; t < mesh.nTris; t++) {
            const v0 = mesh.tris[3*t], v1 = mesh.tris[3*t+1], v2 = mesh.tris[3*t+2];
            const x0 = mesh.nodes[2*v0], y0 = mesh.nodes[2*v0+1];
            const x1 = mesh.nodes[2*v1], y1 = mesh.nodes[2*v1+1];
            const x2 = mesh.nodes[2*v2], y2 = mesh.nodes[2*v2+1];
            const bx0 = Math.min(x0, x1, x2), bx1 = Math.max(x0, x1, x2);
            const by0 = Math.min(y0, y1, y2), by1 = Math.max(y0, y1, y2);
            let candidate = false;
            for (const g of grown) {
                if (bx1 < g.x0 || bx0 > g.x1 || by1 < g.y0 || by0 > g.y1) continue;
                candidate = true; break;
            }
            if (!candidate) continue;
            const xc = (x0+x1+x2)/3, yc = (y0+y1+y2)/3;
            // Distance of the element's closest point (of centroid + vertices) to
            // metal. Within the band iff dS < bw (the former per-point nearSurface
            // test, evaluated once), a triangle straddling the surface must get the
            // surface target.
            const dS = Math.min(surfDist(xc, yc), vertDist(v0), vertDist(v1), vertDist(v2));
            if (!(dS < bw)) continue;
            let tgt = grading ? targetAt(xc, yc) : targetH;
            if (depthSlope > 0 && targetH > 0) tgt = Math.max(tgt, depthTarget(dS));
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

    // The surface-aligned metric of `aniso`. pass(mesh) marks the band triangles with an
    // edge longer than 1 in their metric and returns { marked, nMarked, longest } (the
    // local edge each band triangle bisects, -1 off the band), or null when every band
    // triangle fits.
    function bandMetric() {
        const corners = [];
        // A vertex is a corner when the outline turns by more than CORNER_TURN within
        // the band width of it: a rounding much smaller than delta is a corner at the
        // skin depth's scale, a circle's n-gon is not.
        const cornerTurn = aniso.cornerTurn ?? CORNER_TURN;
        const addLoop = (poly) => {
            const n = poly.length / 2;
            const turn = new Float64Array(n), seg = new Float64Array(n);
            for (let i = 0; i < n; i++) {
                const j = (i + n - 1) % n, k = (i + 1) % n;
                const ax = poly[2*i] - poly[2*j], ay = poly[2*i+1] - poly[2*j+1];
                const bx = poly[2*k] - poly[2*i], by = poly[2*k+1] - poly[2*i+1];
                turn[i] = Math.atan2(ax * by - ay * bx, ax * bx + ay * by);
                seg[i] = Math.hypot(bx, by);   // vertex i to i+1
            }
            for (let i = 0; i < n; i++) {
                let sum = turn[i];
                for (let s = 0, k = i; ; ) {        // forward
                    s += seg[k]; k = (k + 1) % n;
                    if (s > bw || k === i) break;
                    sum += turn[k];
                }
                for (let s = 0, k = i; ; ) {        // backward
                    k = (k + n - 1) % n; s += seg[k];
                    if (s > bw || k === i) break;
                    sum += turn[k];
                }
                if (Math.abs(sum) > cornerTurn) corners.push(poly[2*i], poly[2*i+1]);
            }
        };
        for (const r of rects) {
            if (!r.shape) { corners.push(r.xmin, r.ymin, r.xmax, r.ymin, r.xmax, r.ymax, r.xmin, r.ymax); continue; }
            if (!shields.has(r) && r.shape.poly) addLoop(r.shape.poly);
            if (r.shape.hole) addLoop(r.shape.hole);
        }
        const cornerDist = (x, y) => {
            let d2 = Infinity;
            for (let i = 0; i < corners.length; i += 2) d2 = Math.min(d2, (x - corners[i]) ** 2 + (y - corners[i+1]) ** 2);
            return Math.sqrt(d2);
        };
        const grade = aniso.cornerGrade ?? 0.25;
        const maxAspect = aniso.maxAspect ?? 64;
        // airSide 'corners' (default) refines the air side of the band near corners
        // only, 'all' everywhere within the band width.
        const airCorners = (aniso.airSide ?? 'corners') === 'corners';
        // Unit normal of the nearest surface at (x, y), the gradient of its signed
        // distance, or null on a medial axis (the middle of a trace, the diagonal of a
        // corner) where the gradient collapses and there is no normal.
        const normalAt = (x, y) => {
            let best = null, bd = Infinity;
            for (const r of rects) {
                const d = Math.abs(distToRectBoundary(r, x, y));
                if (d < bd) { bd = d; best = r; }
            }
            const h = 1e-3 * targetH;
            const gx = distToRectBoundary(best, x + h, y) - distToRectBoundary(best, x - h, y);
            const gy = distToRectBoundary(best, x, y + h) - distToRectBoundary(best, x, y - h);
            const g = Math.hypot(gx, gy);
            return g > h ? [gx / g, gy / g] : null;
        };
        return { pass(m) {
            const { nodes, tris, nTris } = m;
            const nodeDist = new Float64Array(m.nNodes).fill(NaN);
            const vertDist = (v) => {
                let d = nodeDist[v];
                if (d !== d) d = nodeDist[v] = surfDist(nodes[2*v], nodes[2*v+1]);
                return d;
            };
            const marked = new Uint8Array(nTris);
            const longest = new Int8Array(nTris).fill(-1);
            let nMarked = 0;
            const P = new Float64Array(6), L2 = new Float64Array(3);
            for (let t = 0; t < nTris; t++) {
                for (let k = 0; k < 3; k++) { const v = tris[3*t+k]; P[2*k] = nodes[2*v]; P[2*k+1] = nodes[2*v+1]; }
                const bx0 = Math.min(P[0], P[2], P[4]), bx1 = Math.max(P[0], P[2], P[4]);
                const by0 = Math.min(P[1], P[3], P[5]), by1 = Math.max(P[1], P[3], P[5]);
                if (!grown.some(g => !(bx1 < g.x0 || bx0 > g.x1 || by1 < g.y0 || by0 > g.y1))) continue;
                const xc = (P[0] + P[2] + P[4]) / 3, yc = (P[1] + P[3] + P[5]) / 3;
                // The closest point to the surface, of the centroid and the vertices.
                let dS = surfDist(xc, yc), px = xc, py = yc;
                for (let k = 0; k < 3; k++) {
                    const d = vertDist(tris[3*t+k]);
                    if (d < dS) { dS = d; px = P[2*k]; py = P[2*k+1]; }
                }
                if (!(dS < bw)) continue;
                const air = airCorners && !rects.some(r => distToRectBoundary(r, xc, yc) < 0);
                let hn = grading ? targetAt(xc, yc) : targetH;
                if (depthSlope > 0) hn = Math.max(hn, depthTarget(dS));
                // The normal at the closest point, or at the centroid when that point is
                // on the surface of two faces (a corner vertex).
                const n = normalAt(px, py) || normalAt(xc, yc);
                // Corner distance of the triangle's nearest point to a corner (of the
                // centroid and the vertices).
                const dc = Math.min(cornerDist(xc, yc), cornerDist(P[0], P[1]), cornerDist(P[2], P[3]), cornerDist(P[4], P[5]));
                // Outside the metal the field is smooth along a face, the band only needs
                // the air side where the corner grading still holds the elements short.
                if (air && grade * dc >= maxAspect * hn) continue;
                const ht = n ? Math.min(maxAspect * hn, Math.max(hn, grade * dc)) : hn;
                let lMax = 0, kMax = 0;
                for (let k = 0; k < 3; k++) {
                    const a = k, b = (k + 1) % 3;
                    const ex = P[2*b] - P[2*a], ey = P[2*b+1] - P[2*a+1];
                    const l2 = ex * ex + ey * ey, en = n ? ex * n[0] + ey * n[1] : 0;
                    L2[k] = n ? en * en / (hn * hn) + Math.max(0, l2 - en * en) / (ht * ht) : l2 / (hn * hn);
                    if (L2[k] > lMax) { lMax = L2[k]; kMax = k; }
                }
                longest[t] = kMax;
                if (lMax > 1) { marked[t] = 1; nMarked++; }
            }
            return nMarked ? { marked, nMarked, longest } : null;
        } };
    }
}

// Conductor loss by volume eddy-current solve.
//   mesh        - mesh with the conductor interiors meshed
//   condRect    - conductor geometry ({rects: [...]}, rectRoles, domain bounds)
//   freq, sigma - frequency (Hz) and reference conductivity (S/m)
//   solveComplexSymmetric - WASM helper from createWasmHelpers ([re; im] layout)
//   Z0          - optional line impedance for the alpha_c conversion
//   opts.wallPEC - {left, right, top, bottom}: which domain walls are metal (from the
//     mesher, `left` already cleared on a half domain). These walls are A = 0 and
//     carry the wall loss. Without it every outer wall is A = 0 and only the bottom
//     one dissipates.
//   opts.wallThick, opts.wallSigma - thickness of each wall's metal (slab impedance,
//     default infinite) and its conductivity (default sigma).
//   opts.oddSymmetry - A = 0 on the symmetry plane x = xmin_domain (odd mode) instead
//     of the natural BC (even mode), for a half-domain pair.
//   opts.modeCurrents - one drive per signal polarity group (positive group first,
//     pre.groups order), full domain only: the target net currents, [1, -1] odd,
//     [1, 1] even, or the modal currents (normalized to sum |I|^2 = 2) of an
//     asymmetric pair. R_total/X_total keep the differential convention (both
//     traces, per-line mode R = R_total/2), L_loop is per line.
//   opts.diffPair - the mesh holds one full trace of a differential pair: on a half
//     domain it carries the full trace current, and R_trace/R_gnd cover both traces.
//   opts.idealRects - per rect, true for a perfect ground held at A = 0, its loss a
//     surface term like the walls.
//   opts.Rq - rms surface roughness (m), gradient model, the same on every surface.
//     The dissipation is scaled by Re(Z_rough)/Rs and the inductance gets the
//     matching reactance increment, which tends to the DC limit as delta grows.
//   opts.surfaceZs(x, y, orient) - per-face surface impedance {re, im} at a face
//     midpoint (orient 'h' top/bottom, 'v' side), for per-side plating. The smooth
//     loss is scaled by the |K|^2-weighted Re(Zs)/Rs over the faces. Overrides Rq.
//   opts.surfaceDZ(x, y, orient) - impedance {re, im} or null added on top of
//     surfaceZs over a face, R += Re(dZ)*|K|^2 and X += Im(dZ)*|K|^2 (meshed plating:
//     the roughness of the plating/bulk interface the mesh holds smooth).
//   opts.rectSigmaRel - sigma_rect / sigma per rect for conductors of different
//     metals. The mass matrix, drive, net current and dissipation carry the ratio per
//     triangle, and the surface scaling is taken per rect against its own Rs.
//   opts.cache - per-caller cache of the frequency-independent assembly (mqsPrecompute).
//   opts.fieldOut - object that receives the solved field, for current density plots.
// Returns { R_trace, R_gnd, R_total, X_total, L_loop, alpha_c, alpha_c_dBm, delta, nDofs }.
// L_loop is the series inductance from Im(Z) = w*L of the solve. X_total is the surface
// reactance per unit length: the |K|^2-weighted integral that gives R_total, taken
// against Im(Zs) instead of Re(Zs). Smooth metal has Zs = Rs*(1 + j), so X_total equals
// R_total there and only rough or plated surfaces separate them.
// Frequency-independent part of the MQS solve on a given mesh: conductor
// classification, DOF map, the S (everywhere), M and Fc (conductor) assembly, and
// the CSR pattern with separate S and M value templates. In K = S + j*beta*M only
// beta = w*mu0*sigma depends on frequency, so each frequency just sets valRe = valS,
// valIm = beta*valM. Cached by the caller (opts.cache), the skin mesh is reused
// across sweep points.
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
    const sigRects = [], sigPol = [], sigIdx = [], gndIdx = [];
    rects.forEach((r, i) => {
        if (!roles || roles[i].is_signal) { sigRects.push(r); sigIdx.push(i); sigPol.push(roles ? (roles[i].polarity || 0) : 0); }
        else gndIdx.push(i);
    });
    const gndRects = gndIdx.map(i => rects[i]);
    // Conductivity of each rect over the reference sigma of the solve.
    const sRelRect = opts.rectSigmaRel || null;
    // Drive groups: one per distinct signal polarity (positive first). A
    // single-ended line has one group; a differential pair two (+/-). Used by
    // the per-conductor-drive path (opts.modeCurrents), the single-drive path
    // ignores the grouping entirely.
    const groups = [...new Set(sigPol)].sort((a, b) => b - a);
    const groupOfSig = sigPol.map(p => groups.indexOf(p));
    const isCondTri = new Uint8Array(nTris);
    const triGroup = new Int32Array(nTris).fill(-1);
    // Rect (index into rects) each conductor triangle belongs to.
    const triRect = new Int32Array(nTris).fill(-1);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const xc = (nodes[2*v0]+nodes[2*v1]+nodes[2*v2])/3;
        const yc = (nodes[2*v0+1]+nodes[2*v1+1]+nodes[2*v2+1])/3;
        // Containment in a signal rect (or shape), keeping the matching rect's
        // index so the triangle lands in its polarity group. Where rects of one kind
        // overlap the later one wins, as for dielectrics.
        let si = -1;
        for (let i = sigRects.length - 1; i >= 0; i--) {
            const r = sigRects[i];
            if (r.shape ? shapeContains(r, xc, yc, TOL)
                : (xc > r.xmin - TOL && xc < r.xmax + TOL && yc > r.ymin - TOL && yc < r.ymax + TOL)) { si = i; break; }
        }
        if (si >= 0) { isCondTri[t] = 1; triGroup[t] = groupOfSig[si]; triRect[t] = sigIdx[si]; }
        else {
            for (let i = gndRects.length - 1; i >= 0; i--) {
                const r = gndRects[i];
                if (r.shape ? shapeContains(r, xc, yc, TOL)
                    : (xc > r.xmin - TOL && xc < r.xmax + TOL && yc > r.ymin - TOL && yc < r.ymax + TOL)) {
                    isCondTri[t] = 2; triRect[t] = gndIdx[i]; break;
                }
            }
        }
    }
    // Ideal grounds (opts.idealRects, per rect): perfect conductors held at A = 0 like
    // the walls, class 3, their loss a surface term in mqsConductorLoss.
    const ideal = opts.idealRects || null;
    const fixedDof = ideal ? new Uint8Array(nNodes + nEdges) : null;
    if (ideal) {
        for (let t = 0; t < nTris; t++) {
            if (isCondTri[t] !== 2 || !ideal[triRect[t]]) continue;
            isCondTri[t] = 3;
            for (let k = 0; k < 3; k++) { fixedDof[tris[3*t+k]] = 1; fixedDof[nNodes + triEdges[3*t+k]] = 1; }
        }
    }

    // DOFs: vertices + edge midpoints, Dirichlet at ground/outer walls.
    // The symmetry plane (x = xmin_domain when symmetry) gets the natural BC.
    // The bottom ground wall is the DOMAIN bottom (ymin_domain — usually 0 after
    // the ground slab is wall-absorbed, but not by construction).
    const xmin_d = condRect.xmin_domain, xmax_d = condRect.xmax_domain;
    const ymax_d = condRect.ymax_domain;
    const ymin_d = condRect.ymin_domain ?? 0;
    // A = 0 on the metal walls (opts.wallPEC) and on the symmetry plane of an odd
    // mode. An open wall is a field truncation and keeps the natural BC, the dual of
    // the static solve's: no current can close through it, so the return current
    // stays in the metal whatever the domain size. With no metal wall the conductor
    // mass term keeps the system regular, which needs a ground rect; without a wall
    // map, or with neither walls nor ground rects, every outer wall is A = 0.
    const wpec = opts.wallPEC;
    const anyWall = !!wpec && (wpec.left || wpec.right || wpec.top || wpec.bottom);
    const perWall = !!wpec && (anyWall || gndRects.length > 0);
    function isDirichletPt(x, y) {
        const onPlane = Math.abs(x - xmin_d) < 1e-9;
        if (perWall) {
            if ((wpec.bottom && Math.abs(y - ymin_d) < 1e-9) || (wpec.top && Math.abs(y - ymax_d) < 1e-9)
                || (wpec.right && Math.abs(x - xmax_d) < 1e-9)) return true;
            return onPlane && (sym === 2 ? !!opts.oddSymmetry : !!wpec.left);
        }
        if (Math.abs(y - ymin_d) < 1e-9 || Math.abs(y - ymax_d) < 1e-9) return true;
        if (Math.abs(x - xmax_d) < 1e-9) return true;
        if ((sym === 1 || opts.oddSymmetry) && onPlane) return true;
        return false;
    }
    const dofOf = new Int32Array(nNodes + nEdges).fill(-1);
    let nF = 0;
    for (let n = 0; n < nNodes; n++) {
        if (!isDirichletPt(nodes[2*n], nodes[2*n+1]) && !(fixedDof && fixedDof[n])) dofOf[n] = nF++;
    }
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const xm = (nodes[2*n0]+nodes[2*n1])/2, ym = (nodes[2*n0+1]+nodes[2*n1+1])/2;
        if (!isDirichletPt(xm, ym) && !(fixedDof && fixedDof[nNodes + e])) dofOf[nNodes + e] = nF++;
    }

    // Assemble S (everywhere), M (all metal: signal + passive grounds) and Fc
    // (signal only, grounds have C = 0, so no source term). condArea is the
    // signal cross-section: it is only used to normalize the signal current.
    // The three element integrals are closed-form on a straight-sided triangle.
    // The P2 stiffness comes from triP2Stiffness and the mass / load are Area
    // x the reference P2_MASS / P2_LOAD. That replaces a 6-point quadrature
    // loop with six basis evaluations per point, the assembly used to dominate
    // the MQS path's JS time on the (skin-refined, hence large) conductor mesh.
    const Nb = 2 * nF;               // length of a [re; im] state vector
    const maxCoo = 36 * nTris;       // one entry per local (i, j) pair
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
        const cond = isCondTri[t] === 1 || isCondTri[t] === 2 ? isCondTri[t] : 0;
        const driven = cond === 1;
        // Area weighted by the relative conductivity: the mass, the drive and the net
        // current all carry sigma per triangle.
        const AreaS = (cond && sRelRect) ? Area * sRelRect[triRect[t]] : Area;
        if (driven) { condArea += AreaS; areaG[triGroup[t]] += AreaS; }
        const FcGt = driven ? FcG[triGroup[t]] : null;
        for (let i = 0; i < 6; i++) {
            const gi = lg[i]; if (gi < 0) continue;
            if (driven) { const f = AreaS * P2_LOAD[i]; Fc[gi] += f; FcGt[gi] += f; }
            for (let j = 0; j < 6; j++) {
                const gj = lg[j]; if (gj < 0) continue;
                // K = S + jβM with β = ωμ₀σ the only frequency-dependent factor:
                // the stiffness and the (metal-only) mass go into separate value
                // templates on one shared index stream.
                // An entry is dropped only when its transpose is zero too: rounding
                // can leave one of S_ij, S_ji exactly zero (right-angle triangles),
                // and a nonsymmetric pattern keeps the solvers off LDL^T.
                const sv = Sl[6*i+j];
                const mv = cond ? AreaS * P2_MASS[6*i+j] : 0;
                if (sv === 0 && Sl[6*j+i] === 0 && mv === 0) continue;
                R[nCoo] = gi; Cc[nCoo] = gj; VS[nCoo] = sv; VM[nCoo] = mv; nCoo++;
            }
        }
    }

    // Symbolic system with separate S / M value templates on one pattern:
    // valS and valM line up entry-for-entry, so a frequency only needs
    // valIm = β·valM next to valRe = valS.
    const csr2 = tripletsToCSRMulti(R.subarray(0, nCoo), Cc.subarray(0, nCoo), nF,
                                    [VS.subarray(0, nCoo), VM.subarray(0, nCoo)]);

    // edge → one adjacent triangle (for the ground-loss surface integral)
    const edgeToTri = new Int32Array(nEdges).fill(-1);
    for (let t = 0; t < nTris; t++)
        for (let k = 0; k < 3; k++)
            if (edgeToTri[triEdges[3*t+k]] === -1) edgeToTri[triEdges[3*t+k]] = t;

    return {
        isCondTri, dofOf, nF, Nb, Fc, condArea, edgeToTri,
        groups, triGroup, triRect, sRelRect, FcG, areaG,
        rowPtr: csr2.rowPtr, colIdx: csr2.colIdx,
        valS: csr2.vals[0], valM: csr2.vals[1],
    };
}

// mqsPrecompute through the caller's cache (opts.cache), reused while the mesh
// (identity), the symmetry mode and the ideal-ground set match.
function cachedPrecompute(mesh, condRect, opts) {
    const cc = opts.cache;
    const idealKey = opts.idealRects ? Array.from(opts.idealRects, v => (v ? 1 : 0)).join('') : '';
    if (cc && cc.mesh === mesh && cc.odd === !!opts.oddSymmetry && (cc.ideal ?? '') === idealKey) return cc.pre;
    const pre = mqsPrecompute(mesh, condRect, opts);
    if (cc) { cc.mesh = mesh; cc.odd = !!opts.oddSymmetry; cc.ideal = idealKey; cc.pre = pre; }
    return pre;
}

// Loop inductance of the same problem with every conductor perfect, on the same
// mesh with the same walls and drive convention. L_loop(f) minus this is the
// internal inductance: the volume skin/proximity part plus the excess of the
// finite-sigma current distribution over the surface (TEM) one, which for a film
// thinner than the skin depth is the dominant term (a 100 nm microstrip: ~20 nH/m
// at low frequency, where the deep-skin identity L_int = R/omega would give
// microhenries). Referencing against the static L_external instead would not
// work: that one comes from the open-boundary static solve, this box has A = 0
// on every wall.
//
// In the sigma -> infinity limit A is constant on each conductor, so the exterior
// problem is Laplace's equation with phi = 1 on the driven group, 0 on the other
// metal and on the walls, and the inductance is mu0 over the Dirichlet energy
// Phi = phi^T S phi (the magnetostatic twin of C = eps0 * energy). Conductor
// interiors carry a constant and add nothing to Phi, so the full S of the MQS
// assembly serves as is. Several driven groups (modeCurrents) give the matrix
// Phi_jk from unit solves and L = mu0 * Phi^-1, contracted with the target
// currents the way the MQS Zmode is. Real solve, frequency-independent, cached
// per skin mesh by the caller. Same I_mesh normalization as mqsConductorLoss.
export function mqsPecInductance(mesh, condRect, solveSparseMulti, opts = {}) {
    const { tris, triEdges, nNodes, nTris } = mesh;
    const sym = condRect.symmetry > 1 ? 2 : 1;
    const pre = cachedPrecompute(mesh, condRect, opts);
    const { isCondTri, dofOf, nF, triGroup, rowPtr, colIdx, valS } = pre;
    const nG = pre.groups.length;
    // Metal DOF classes: 1 + group for signal metal, -1 for passive ground metal,
    // 0 free. Signal wins where a DOF touches both (touching rects).
    const cls = new Int32Array(nF);
    for (let t = 0; t < nTris; t++) {
        const c = isCondTri[t]; if (!c) continue;
        const v = c === 1 ? 1 + triGroup[t] : -1;
        for (let k = 0; k < 3; k++) {
            const gv = dofOf[tris[3*t+k]], ge = dofOf[nNodes + triEdges[3*t+k]];
            if (gv >= 0 && cls[gv] <= 0 && (v > 0 || cls[gv] === 0)) cls[gv] = v;
            if (ge >= 0 && cls[ge] <= 0 && (v > 0 || cls[ge] === 0)) cls[ge] = v;
        }
    }
    // Reduced system on the free DOFs, one right-hand side per driven group.
    const red = new Int32Array(nF).fill(-1);
    let nR = 0;
    for (let i = 0; i < nF; i++) if (cls[i] === 0) red[i] = nR++;
    const rp = new Int32Array(nR + 1), ci = [], cv = [];
    const rhs = Array.from({ length: nG }, () => new Float64Array(nR));
    for (let i = 0, r = 0; i < nF; i++) {
        if (cls[i] !== 0) continue;
        for (let k = rowPtr[i]; k < rowPtr[i+1]; k++) {
            const j = colIdx[k], sv = valS[k];
            if (cls[j] === 0) { ci.push(red[j]); cv.push(sv); }
            else if (cls[j] > 0) rhs[cls[j] - 1][r] -= sv;
        }
        rp[++r] = ci.length;
    }
    const sols = nR > 0 ? solveSparseMulti(nR, { rowPtr: rp, colIdx: new Int32Array(ci), valRe: new Float64Array(cv) }, rhs) : [];
    // Full potentials (1 on the group, solution outside, 0 elsewhere) and the
    // energy matrix Phi_jk = phi_j^T S phi_k.
    const phi = [];
    for (let g = 0; g < nG; g++) {
        const v = new Float64Array(nF);
        for (let i = 0; i < nF; i++) v[i] = cls[i] === 0 ? sols[g][red[i]] : (cls[i] === 1 + g ? 1 : 0);
        phi.push(v);
    }
    const Phi = [];
    for (let j = 0; j < nG; j++) {
        Phi.push([]);
        const Sp = new Float64Array(nF);
        for (let i = 0; i < nF; i++) {
            let a = 0;
            for (let k = rowPtr[i]; k < rowPtr[i+1]; k++) a += valS[k] * phi[j][colIdx[k]];
            Sp[i] = a;
        }
        for (let k = 0; k < nG; k++) { let a = 0; for (let i = 0; i < nF; i++) a += Sp[i] * phi[k][i]; Phi[j].push(a); }
    }
    // Same current normalization as the MQS solve.
    const I_mesh = (sym === 2 && !opts.diffPair) ? 0.5 : 1;
    const I = (opts.modeCurrents && sym === 1 && nG > 1 && opts.modeCurrents.length === nG) ? opts.modeCurrents : null;
    if (!I) return I_mesh * MU0 / Phi[0][0];
    // L = mu0 * Phi^-1 contracted with the mode currents: (I^T L I) / (I^T I).
    const M = Phi.map(row => row.slice());
    const B = I.slice();
    for (let c = 0; c < nG; c++) {   // Gaussian elimination with partial pivoting
        let p = c; for (let r = c + 1; r < nG; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
        [M[c], M[p]] = [M[p], M[c]]; [B[c], B[p]] = [B[p], B[c]];
        for (let r = 0; r < nG; r++) {
            if (r === c) continue;
            const f = M[r][c] / M[c][c]; if (f === 0) continue;
            for (let k = c; k < nG; k++) M[r][k] -= f * M[c][k];
            B[r] -= f * B[c];
        }
    }
    let num = 0, den = 0;
    for (let k = 0; k < nG; k++) { num += (B[k] / M[k][k]) * I[k]; den += I[k] * I[k]; }
    return MU0 * num / den;
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

export function mqsConductorLoss(mesh, condRect, freq, sigma, solveComplexSymmetric, Z0 = 0, opts = {}) {
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

    // Frequency-invariant assembly, and the per-frequency unit solves below, through
    // the caller's cache.
    const cc = opts.cache;
    const pre = cachedPrecompute(mesh, condRect, opts);
    const { isCondTri, dofOf, nF, Nb, Fc, condArea, edgeToTri, triGroup } = pre;
    const lg = new Int32Array(6);

    // Per-frequency system K = S + jβM on the cached pattern: only the imaginary
    // template depends on frequency. Vectors are [re(nF); im(nF)] = [Ar; Ai].
    const beta = omega * MU0 * sigma;
    const valIm = new Float64Array(pre.valM.length);
    for (let k = 0; k < valIm.length; k++) valIm[k] = beta * pre.valM[k];
    const csr = { rowPtr: pre.rowPtr, colIdx: pre.colIdx, valRe: pre.valS, valIm };

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
        // the second mode of the pair at the same frequency skips the factorization.
        let sols = (cc && cc.sols && cc.solsMesh === mesh && cc.solsBeta === beta) ? cc.sols : null;
        if (!sols) {
            const rhsList = pre.FcG.map(F => {
                const r = new Float64Array(Nb);
                for (let i = 0; i < nF; i++) r[i] = MU0 * sigma * F[i];
                return r;
            });
            sols = solveComplexSymmetric(nF, csr, rhsList);
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
        [sol] = solveComplexSymmetric(nF, csr, [rhs]);

        // Rescale C for trace current 1 A. A single-ended half mesh carries half the
        // line current, whether its signal straddles the symmetry plane or is one of a
        // mirrored pair of traces. A differential half mesh holds one full trace and
        // carries the full per-trace current. Ground rects are excluded, they carry
        // return current, not the normalized drive current.
        const I_mesh = (sym === 2 && !opts.diffPair) ? 0.5 : 1;
        let fr = 0, fi = 0;
        for (let i = 0; i < nF; i++) { fr += Fc[i] * sol[i]; fi += Fc[i] * sol[nF + i]; }
        const dR = sigma * (condArea + omega * fi);
        const dI = -sigma * omega * fr;
        const dMag2 = dR*dR + dI*dI;
        Cr = I_mesh * dR / dMag2; Ci = -I_mesh * dI / dMag2;
    }
    const Cmag2 = Cr*Cr + Ci*Ci;
    // The solved field for a current density plot (sampleMqsCurrent): J/σ = C(u − jωA)
    // with the drive u of the triangle's class and group.
    if (opts.fieldOut) Object.assign(opts.fieldOut, {
        mesh, sol, nF, dofOf, isCondTri, triGroup, triRect: pre.triRect, sRel: pre.sRelRect,
        sigma, omega, Cr, Ci, CgR, CgI,
    });

    // Conductor dissipation: J/σ = C*(u - jωA1) with the drive u = 1 in the
    // signal (class 1) and u = 0 in passive ground rects (class 2, pure eddy /
    // return current). Signal and ground-rect dissipation accumulate separately
    // so the ground share reports (and plating-scales) with R_gnd, not R_trace.
    let Psig = 0, PgndRect = 0;
    const sRel = pre.sRelRect;
    const Prect = new Float64Array(rects.length);
    for (let t = 0; t < nTris; t++) {
        const cls = isCondTri[t];
        if (!cls || cls === 3) continue;
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const x0 = nodes[2*v0], y0 = nodes[2*v0+1];
        const Area = 0.5 * Math.abs((nodes[2*v1] - x0) * (nodes[2*v2+1] - y0) - (nodes[2*v2] - x0) * (nodes[2*v1+1] - y0));
        lg[0] = dofOf[v0]; lg[1] = dofOf[v1]; lg[2] = dofOf[v2];
        for (let k = 0; k < 3; k++) lg[3+k] = dofOf[nNodes + triEdges[3*t+k]];
        // Drive term: 1 (single) or the group's complex C_k (multi) in the
        // signal, 0 in passive ground rects (pure eddy / return current).
        let dvR = 0, dvI = 0;
        if (cls === 1) { if (CgR) { const g = triGroup[t]; dvR = CgR[g]; dvI = CgI[g]; } else dvR = 1; }
        let Ptri = 0;
        for (let q = 0; q < NQ; q++) {
            const w = QW[q] * Area;
            let aR = 0, aI = 0;
            for (let k = 0; k < 6; k++) {
                const g = lg[k]; if (g < 0) continue;
                const Nk = P2_AT_Q[6*q + k];
                aR += Nk * sol[g]; aI += Nk * sol[nF + g];
            }
            const uR = dvR + omega * aI, uI = dvI - omega * aR;
            const eR = Cr * uR - Ci * uI, eI = Cr * uI + Ci * uR;
            Ptri += 0.5 * sigma * (eR*eR + eI*eI) * w;
        }
        if (sRel) Ptri *= sRel[pre.triRect[t]];
        Prect[pre.triRect[t]] += Ptri;
        if (cls === 1) Psig += Ptri; else PgndRect += Ptri;
    }

    // Domain-wall loss
    //
    // Flat-surface skin formula on the tangential H at each wall
    // that is actually metal. |Hx| = |∂A/∂y|/μ₀ on a horizontal wall, |Hy| =
    // |∂A/∂x|/μ₀ on a vertical one. (edgeToTri comes precomputed from mqsPrecompute.)
    //
    // The metal walls are A = 0 in the MQS solve and dissipate. An 'open' far-field
    // truncation and the symmetry plane are boundary conditions, not surfaces.
    // opts.wallPEC carries that distinction from the mesher (which clears `left` on a
    // half domain, so the symmetry plane can never be mistaken for metal).
    let Pgnd = 0, Xgnd = 0;
    // Half the |K|^2 integral over the walls and ideal grounds, unweighted: the
    // surface-layer increment of roughness or plating on them is (Zs - Zs_smooth) times it.
    let wallS2 = 0;
    // Per-face plating weights (∮|K|²dl and Σ Re/Im(Zs)·|K|²dl) over the ground and
    // conductor surfaces, used below to scale the smooth loss per face.
    // gndDXS: the reactance increment Im(Zs) - Rs of each segment against its own metal.
    let gndS = 0, gndZreS = 0, gndDXS = 0;
    // Without a wall map only the bottom, ground on every geometry the mesher builds.
    const wp = opts.wallPEC || { bottom: true };
    const onLine = (a, b, v) => Math.abs(a - v) < 1e-9 && Math.abs(b - v) < 1e-9;
    // The metal wall this edge lies on, else null.
    function wallOf(x0, y0, x1, y1) {
        if (wp.bottom && onLine(y0, y1, ymin_d)) return 'bottom';
        if (wp.top && onLine(y0, y1, ymax_d)) return 'top';
        if (wp.left && onLine(x0, x1, xmin_d)) return 'left';
        if (wp.right && onLine(x0, x1, xmax_d)) return 'right';
        return null;
    }
    // Smooth surface impedance of a wall: the slab formula (1+j)Rs*coth((1+j)d/delta)
    // for an absorbed slab of thickness d with a field-free back, which is Rs(1+j)
    // for d well above delta and tends to the sheet resistance 1/(sigma*d) with a
    // vanishing reactance as delta grows past d. A 'gnd' boundary is infinitely thick.
    // The walls are the bulk ground metal: opts.wallSigma differs from the solve's
    // sigma when the meshed conductors run at another metal (solid plating).
    const sigmaW = opts.wallSigma ?? sigma;
    const deltaW = Math.sqrt(2 / (omega * MU0 * sigmaW)), RsW = 1 / (sigmaW * deltaW);
    const wt = opts.wallThick || {};
    const zWall = {}, wS1 = {}, wS2 = {};
    for (const w of ['bottom', 'top', 'left', 'right']) {
        wS1[w] = 0; wS2[w] = 0;
        const z = slabCoth((wt[w] ?? Infinity) / deltaW);
        zWall[w] = { re: RsW * z.re, im: RsW * z.im };
    }
    // |dA/dn|^2 at the 3 Gauss points of the edge (x0, y0)-(x1, y1), from the P2 field
    // of triangle t beside it, with n the edge normal (exactly an axis on the faces of
    // a rect, slanted on a polygon side).
    const lge = new Int32Array(6), g2 = new Float64Array(3);
    const edgeGrad2 = (t, x0, y0, x1, y1) => {
        const len = Math.hypot(x1 - x0, y1 - y0);
        const nx = (y0 - y1) / len, ny = (x1 - x0) / len;
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const { coeff } = triCoefficients(nodes, v0, v1, v2);
        lge[0] = dofOf[v0]; lge[1] = dofOf[v1]; lge[2] = dofOf[v2];
        for (let k = 0; k < 3; k++) lge[3+k] = dofOf[nNodes + triEdges[3*t+k]];
        for (let q = 0; q < 3; q++) {
            const xq = x0 + GL3p[q] * (x1 - x0), yq = y0 + GL3p[q] * (y1 - y0);
            let gR = 0, gI = 0;
            for (let k = 0; k < 6; k++) {
                const g = lge[k]; if (g < 0) continue;
                const gr = k < 3 ? lvGrad(coeff, k, xq, yq) : leGrad(coeff, edgeVerts[k-3][0], edgeVerts[k-3][1], xq, yq);
                const gn = gr[0] * nx + gr[1] * ny;
                gR += gn * sol[g]; gI += gn * sol[nF + g];
            }
            g2[q] = gR*gR + gI*gI;
        }
        return g2;
    };
    // The two triangles of each interior edge, eA[2e] and eA[2e+1] (-1 on the boundary).
    let eA = null;
    if (opts.idealRects || opts.surfaceZs) {
        eA = new Int32Array(2 * nEdges).fill(-1);
        for (let t = 0; t < nTris; t++) for (let k = 0; k < 3; k++) {
            const e = triEdges[3*t+k];
            if (eA[2*e] < 0) eA[2*e] = t; else eA[2*e+1] = t;
        }
    }
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const x0 = nodes[2*n0], y0 = nodes[2*n0+1];
        const x1 = nodes[2*n1], y1 = nodes[2*n1+1];
        const wall = wallOf(x0, y0, x1, y1);
        if (!wall) continue;
        const orient = (wall === 'bottom' || wall === 'top') ? 'h' : 'v';
        const adj = edgeToTri[e];
        // A meshed conductor sitting on the wall (a GCPW pour or via slab reaching
        // it) bonds to it. There is no exposed wall surface there, and the adjacent
        // triangle is metal, so ∂A/∂n is an interior gradient, not the exterior H.
        // That segment's dissipation belongs to the volume eddy term, not here.
        if (adj < 0 || isCondTri[adj]) continue;
        const L = Math.hypot(x1 - x0, y1 - y0);
        // Tangential H comes from the gradient component along the wall normal.
        const G2 = edgeGrad2(adj, x0, y0, x1, y1);
        // Walls are bare metal (never plated). Evaluate Zs once at the edge midpoint.
        const Zg = opts.surfaceZs ? opts.surfaceZs((x0 + x1) / 2, (y0 + y1) / 2, orient) : null;
        for (let q = 0; q < 3; q++) {
            const K2 = Cmag2 * G2[q] / (MU0*MU0);
            wS1[wall] += Math.sqrt(K2) * GL3w[q] * L;
            wS2[wall] += K2 * GL3w[q] * L;
            if (Zg) {
                const Sseg = G2[q] * GL3w[q] * L;   // ∝ |K|² (global factors cancel in the ratio)
                gndS += Sseg; gndZreS += Zg.re * Sseg; gndDXS += (Zg.im - RsW) * Sseg;
            }
        }
    }

    // Per-wall totals: slab impedance times the PEC-distribution integral, the
    // resistance reduced by the lateral spreading of a thin wall (wallSpreadFactor). A
    // wall is a plane of unlimited width, so the spreading has no floor and the wall
    // resistance vanishes towards DC, where the wall is an ideal return. The reactance
    // keeps the slab value of the confined current (the DC convention of
    // dcLineParameters). A wall of infinite thickness has no spreading. The half-domain
    // integrals cover half the wall; the width ratio is the same on either.
    for (const w of ['bottom', 'top', 'left', 'right']) {
        const S2 = wS2[w];
        if (!(S2 > 0)) continue;
        const d = wt[w] ?? Infinity;
        let g = 1;
        if (d < Infinity) {
            const Wk = sym * wS1[w] * wS1[w] / S2;
            g = wallSpreadFactor(2 * Math.PI * (deltaW * deltaW / d) / Wk);
        }
        Pgnd += 0.5 * zWall[w].re * S2 * g;
        Xgnd += 0.5 * zWall[w].im * S2;
        wallS2 += 0.5 * S2;
    }

    // Ideal grounds (class 3, opts.idealRects): perfect conductors in the solve, like
    // the walls, with the same surface term per rect: the slab impedance of its
    // thickness (the thin side) at its own sigma times the perfect-conductor surface
    // current on its faces, the resistance reduced by the unlimited lateral spreading.
    // uIdeal: the largest spreading parameter delta^2 / (d W_K), for the caller's blend.
    let uIdeal = 0;
    if (opts.idealRects) {
        const nR = rects.length, iS1 = new Float64Array(nR), iS2 = new Float64Array(nR);
        for (let e = 0; e < nEdges; e++) {
            const ta = eA[2*e], tb = eA[2*e+1];
            if (ta < 0 || tb < 0) continue;
            const aI = isCondTri[ta] === 3, bI = isCondTri[tb] === 3;
            if (aI === bI) continue;
            const ext = aI ? tb : ta, cnd = aI ? ta : tb;
            if (isCondTri[ext]) continue;            // bonded to other metal, no exposed face
            const n0 = edges[2*e], n1 = edges[2*e+1];
            const x0 = nodes[2*n0], y0 = nodes[2*n0+1], x1 = nodes[2*n1], y1 = nodes[2*n1+1];
            const horiz = Math.abs(y1 - y0) < Math.abs(x1 - x0);
            const L = Math.hypot(x1 - x0, y1 - y0);
            const ri = pre.triRect[cnd];
            const G2 = edgeGrad2(ext, x0, y0, x1, y1);
            let Sseg = 0;
            for (let q = 0; q < 3; q++) {
                const K2 = Cmag2 * G2[q] / (MU0*MU0);
                iS1[ri] += Math.sqrt(K2) * GL3w[q] * L;
                iS2[ri] += K2 * GL3w[q] * L;
                Sseg += G2[q] * GL3w[q] * L;
            }
            if (opts.surfaceZs) {
                // Weighted against the wall metal like the walls: Zs relative to the rect's
                // own Rs for R. The reactance increment is against its own Rs.
                const sR = pre.sRelRect ? pre.sRelRect[ri] : 1;
                const Zs = opts.surfaceZs((x0 + x1) / 2, (y0 + y1) / 2, horiz ? 'h' : 'v');
                const RsR = Rs / Math.sqrt(sR), k = RsW / RsR;
                gndS += Sseg; gndZreS += Zs.re * k * Sseg;
                gndDXS += (Zs.im - RsR) * Sseg;
            }
        }
        for (let ri = 0; ri < nR; ri++) {
            const S2 = iS2[ri];
            if (!(S2 > 0)) continue;
            const r = rects[ri];
            const sigmaR = sigma * (pre.sRelRect ? pre.sRelRect[ri] : 1);
            const deltaR = Math.sqrt(2 / (omega * MU0 * sigmaR)), RsR = 1 / (sigmaR * deltaR);
            const d = condThinDim(r);
            const z = slabCoth(d / deltaR), zr = RsR * z.re, zi = RsR * z.im;
            // A rect on the symmetry plane is seen by half.
            const straddles = sym === 2 && r.xmin <= xmin_d + 1e-12;
            const Wk = (straddles ? 2 : 1) * iS1[ri] * iS1[ri] / S2;
            const u = deltaR * deltaR / d / Wk;
            uIdeal = Math.max(uIdeal, u);
            Pgnd += 0.5 * zr * S2 * wallSpreadFactor(2 * Math.PI * u);
            Xgnd += 0.5 * zi * S2;
            wallS2 += 0.5 * S2;
        }
    }

    // Conductor-surface weights: the trace loss is a volume integral, so to apply a
    // per-face impedance we weight each face by its surface current ∮|K|²dl, taken
    // on the EXTERIOR (dielectric) side of the face — like the ground above. K is
    // the tangential H: ∂A/∂n along the face normal.
    // Signal faces and ground-rect faces accumulate separate buckets: each scales
    // its own volume loss (a plated trace next to a bare ground must not dilute).
    let trS = 0, trZreS = 0, trZimS = 0;
    let grS = 0, grZreS = 0, grZimS = 0;
    // The same weights per rect, for conductors of different metals.
    const rS = new Float64Array(rects.length), rZreS = new Float64Array(rects.length),
        rZimS = new Float64Array(rects.length);
    // Additive surface terms (opts.surfaceDZ), 0.5 * dZ * |K|^2 summed like the walls.
    let dRtr = 0, dXtr = 0, dRgr = 0, dXgr = 0;
    if (opts.surfaceZs) {
        for (let e = 0; e < nEdges; e++) {
            const ta = eA[2*e], tb = eA[2*e+1];
            if (ta < 0 || tb < 0) continue;          // boundary edge, not a cond/dielectric interface
            if (isCondTri[ta] === 3 || isCondTri[tb] === 3) continue;   // ideal ground, below
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
            // Later rects first: where rects overlap the later one is the metal.
            for (let ri = rects.length - 1; ri >= 0; ri--) {
                const r = rects[ri];
                if (r.shape) {
                    // A polygon side: snap onto the nearest point of its boundary.
                    if (!shapeContains(r, ccx, ccy, TOL)) continue;
                    ({ x: qx, y: qy } = shapeFaceAt(r.shape, qx, qy));
                    break;
                }
                if (ccx <= r.xmin - TOL || ccx >= r.xmax + TOL || ccy <= r.ymin - TOL || ccy >= r.ymax + TOL) continue;
                if (horiz) qy = Math.abs(qy - r.ymax) < Math.abs(qy - r.ymin) ? r.ymax : r.ymin;
                else qx = Math.abs(qx - r.xmax) < Math.abs(qx - r.xmin) ? r.xmax : r.xmin;
                break;
            }
            const Zs = opts.surfaceZs(qx, qy, horiz ? 'h' : 'v');
            const L = Math.hypot(x1 - x0, y1 - y0);
            // The normal gradient, which gives the tangential H.
            const G2 = edgeGrad2(ext, x0, y0, x1, y1);
            let Sseg = 0;
            for (let q = 0; q < 3; q++) Sseg += G2[q] * GL3w[q] * L;
            if (isCondTri[cnd] === 1) { trS += Sseg; trZreS += Zs.re * Sseg; trZimS += Zs.im * Sseg; }
            else { grS += Sseg; grZreS += Zs.re * Sseg; grZimS += Zs.im * Sseg; }
            const ri = pre.triRect[cnd];
            rS[ri] += Sseg; rZreS[ri] += Zs.re * Sseg; rZimS[ri] += Zs.im * Sseg;
            const dZ = opts.surfaceDZ ? opts.surfaceDZ(qx, qy, horiz ? 'h' : 'v') : null;
            if (dZ) {
                const K2 = 0.5 * Cmag2 * Sseg / (MU0 * MU0);
                if (isCondTri[cnd] === 1) { dRtr += dZ.re * K2; dXtr += dZ.im * K2; }
                else { dRgr += dZ.re * K2; dXgr += dZ.im * K2; }
            }
        }
    }

    // Totals: R = 2P/|I|². Single-ended: R_total is the line R for
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
    const X_gw_smooth = 2 * sym * Xgnd;
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
    // loss weighted by Im(Zs)/Rs where R is weighted by Re(Zs)/Rs. They carry the
    // same units and the same differential convention as R, so the caller halves
    // them alike. Smooth eddy metal has Im(Zs) = Re(Zs) = Rs, so X_trace and X_gr
    // seed equal to R; the walls seed from their slab reactance. Each branch
    // overwrites its X from the SMOOTH value before R is scaled by psiR, per face
    // bucket: signal faces -> R_trace, ground-rect faces -> R_gr, boundary walls
    // -> R_gw (walls are bare metal unless surfaceZs says otherwise).
    // The reactance increment of a meshed conductor is (Im(Zs) - Rs) * |K|^2 over its
    // faces. The volume R gives that integral as R / Rs only in the skin regime: a
    // conductor thinner than delta has R near R_dc, and R * (Im(Zs)/Rs - 1) / omega
    // would grow as 1/sqrt(f). Dividing by the slab resistance factor
    // Re[(1+j) coth((1+j) d/delta)] recovers |K|^2 at every delta (it is 1 for a thick
    // conductor and delta/d for a thin one, where R = |K|^2 / (sigma d)). A signal
    // trace carries current on both faces, so each face sees half the thickness;
    // ground rects are one-sided slabs. The walls and ideal grounds add the same
    // increment, (Zs - Zs_smooth) * |K|^2 (wallS2), to their slab reactance.
    const slabR = d => slabCoth(d / delta).re;
    const rolesX = condRect.rectRoles || null;
    let dSig = Infinity, dGr = Infinity;
    rects.forEach((r, i) => {
        const d = condThinDim(r);
        if (!rolesX || rolesX[i].is_signal) dSig = Math.min(dSig, d / 2); else dGr = Math.min(dGr, d);
    });
    const kTrace = 1 / slabR(dSig), kGr = 1 / slabR(dGr);
    let X_trace = R_trace, X_gr = R_gr, X_gw = X_gw_smooth;
    const X_smooth_total = R_trace + R_gr + X_gw_smooth;   // before psi rewrites any part
    if (sRel) {
        // Per rect: its own volume loss, scaled by its faces' impedance over its own
        // smooth Rs = Rs / sqrt(sigma_rect / sigma), with its own slab factor.
        R_trace = 0; R_gr = 0; X_trace = 0; X_gr = 0;
        let sigS = 0, sigPsi = 0;
        rects.forEach((r, i) => {
            const Ri = 2 * sym * Prect[i];
            const RsI = Rs / Math.sqrt(sRel[i]), deltaI = delta / Math.sqrt(sRel[i]);
            const isSig = !rolesX || rolesX[i].is_signal;
            let psiR = 1, psiX = 1;
            if (opts.surfaceZs && rS[i] > 0) {
                psiR = rZreS[i] / (RsI * rS[i]); psiX = rZimS[i] / (RsI * rS[i]);
            } else if (!opts.surfaceZs && Rq > 0) {
                const Zs = calculate_Zrough(freq, sigma * sRel[i], Rq);
                psiR = Zs.re / RsI; psiX = Zs.im / RsI;
            }
            const d = condThinDim(r);
            const k = 1 / slabCoth((isSig ? d / 2 : d) / deltaI).re;
            const Xi = Ri * (1 + k * (psiX - 1));
            if (isSig) { R_trace += Ri * psiR; X_trace += Xi; sigS += Ri; sigPsi += Ri * psiR; }
            else { R_gr += Ri * psiR; X_gr += Xi; }
        });
        if (sigS > 0) PsiR = sigPsi / sigS;
        if (gndS > 0) {
            X_gw = X_gw_smooth + gndDXS / gndS * 2 * sym * wallS2;
            R_gw *= gndZreS / (RsW * gndS);
        } else if (!opts.surfaceZs && Rq > 0) {
            const Zs = calculate_Zrough(freq, sigmaW, Rq);
            X_gw = X_gw_smooth + (Zs.im - RsW) * 2 * sym * wallS2;
            R_gw *= Zs.re / RsW;
        }
    } else if (opts.surfaceZs) {
        if (trS > 0) {
            const psiR = trZreS / (Rs * trS), psiX = trZimS / (Rs * trS);
            X_trace = R_trace * (1 + kTrace * (psiX - 1));
            R_trace *= psiR;
            PsiR = psiR;
        }
        if (grS > 0) {
            const psiR = grZreS / (Rs * grS), psiX = grZimS / (Rs * grS);
            X_gr = R_gr * (1 + kGr * (psiX - 1));
            R_gr *= psiR;
        }
        if (gndS > 0) {
            const psiR = gndZreS / (RsW * gndS);
            X_gw = X_gw_smooth + gndDXS / gndS * 2 * sym * wallS2;
            R_gw *= psiR;
        }
    } else if (Rq > 0) {
        const Zs = calculate_Zrough(freq, sigma, Rq);
        PsiR = Zs.re / Rs;
        const PsiX = Zs.im / Rs;
        X_trace = R_trace * (1 + kTrace * (PsiX - 1));
        X_gr = R_gr * (1 + kGr * (PsiX - 1));
        X_gw = X_gw_smooth + RsW * (PsiX - 1) * 2 * sym * wallS2;
        R_trace *= PsiR;
        R_gr *= PsiR;
        R_gw *= PsiR;
    }

    R_trace += 2 * sym * dRtr; X_trace += 2 * sym * dXtr;
    R_gr += 2 * sym * dRgr; X_gr += 2 * sym * dXgr;
    const R_gnd = R_gr + R_gw;
    const X_gnd = X_gr + X_gw;
    const R_total = R_trace + R_gnd;
    const X_total = X_trace + X_gnd;
    // The surface-reactance increment on the loop inductance is (X − X_smooth)/ω by
    // construction, so it is formed ONCE from X_total here rather than accumulated
    // branch by branch. That keeps L_loop and X_total from ever describing different
    // surfaces, and it is identically zero for smooth metal (X_total === X_smooth_total).
    const L_surface = perTrace * (X_total - X_smooth_total) / omega;
    L_loop += L_surface;
    // Internal inductance of the domain walls. They are A = 0 in the solve, so
    // their skin layer is not in L_loop; the smooth slab reactance over omega
    // supplies it (their rough/plated increment is in L_loop above).
    const L_wall = perTrace * X_gw_smooth / omega;
    const alpha_c = Z0 > 0 ? R_total / (2 * Z0) : NaN;
    return {
        R_trace, R_gnd, R_total, X_total, L_loop, L_wall, L_surface, PsiR, uIdeal,
        alpha_c, alpha_c_dBm: alpha_c * 8.686,
        Rs, delta, nDofs: nF, modeZ: Zmode,
    };
}
