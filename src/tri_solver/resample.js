// Resample a triangular-mesh static FEM solution (P2 scalar potential) onto a
// regular grid, producing the SAME { x, y, V[ny][nx], Ex[ny][nx], Ey[ny][nx] }
// shape the rectilinear FDM solver exposes — so plot.js (heatmap/contour) and
// streamlines.js work unchanged on the triangular backend's output.
//
// We resample the static potential φ and E = −∇φ (real-valued, always available,
// and the same quasi-TEM field the FDM plots), not the complex full-wave
// eigenvector. The triangle mesh itself is also returned for the mesh overlay.

import { triCoefficients, lv, le, lvGrad, leGrad } from './tri_fem.js';
import { shapeContains, distToShapeBoundary } from '../shapes.js';
import { evalFieldsAtPoint } from './tri_ms_solver.js';

const EDGE_VERTS = [[0, 1], [1, 2], [2, 0]];

// Build the regular sample grid for `domain`, honoring an optional caller-supplied
// grid (normally buildGridFromMesh's mesh-derived graded grid). Shared by
// resampleStatic / resampleModeField.
function makeGrid(domain, opts) {
    const { x_min, x_max, y_min, y_max } = domain;
    if (opts.grid && opts.grid.x && opts.grid.y) return { x: opts.grid.x, y: opts.grid.y };
    const aspect = (x_max - x_min) / (y_max - y_min || 1);
    const nLong = opts.resolution || 240;
    const nx = aspect >= 1 ? nLong : Math.max(40, Math.round(nLong * aspect));
    const ny = aspect >= 1 ? Math.max(40, Math.round(nLong / aspect)) : nLong;
    const x = new Float64Array(nx), y = new Float64Array(ny);
    for (let i = 0; i < nx; i++) x[i] = x_min + (x_max - x_min) * i / (nx - 1);
    for (let j = 0; j < ny; j++) y[j] = y_min + (y_max - y_min) * j / (ny - 1);
    return { x, y };
}

// Build a graded rectilinear sample grid from the triangle mesh itself. Grid-line
// density follows the mesh's node density via quantile decimation, so the plot grid
// is fine exactly where the (adaptively refined) mesh is fine — corners, surfaces,
// gaps — with no dependency on the FDM mesher's geometry heuristics (which know
// nothing about where the solution actually needed resolution). Conductor and
// dielectric interface coordinates (opts.forcedX/forcedY) are inserted exactly so
// thin features keep crisp plateaus. opts.mirrorX mirrors the node x-coordinates
// across x=0 for a half-domain symmetry solve, giving a symmetric full-domain grid.
export function buildGridFromMesh(mesh, domain, opts = {}) {
    const res = Math.max(60, opts.resolution || 300);
    const { nodes, nNodes } = mesh;
    const mirrorX = !!opts.mirrorX;
    const xs = new Float64Array(mirrorX ? 2 * nNodes : nNodes);
    const ys = new Float64Array(nNodes);
    for (let i = 0; i < nNodes; i++) {
        xs[i] = nodes[2 * i];
        ys[i] = nodes[2 * i + 1];
        if (mirrorX) xs[nNodes + i] = -xs[i];
    }
    return {
        x: quantileAxis(xs, domain.x_min, domain.x_max, res, opts.forcedX || []),
        y: quantileAxis(ys, domain.y_min, domain.y_max, res, opts.forcedY || []),
    };
}

// One axis of buildGridFromMesh: decimate the sorted coordinate multiset to ~n
// quantile picks (line density ∝ node density), then union with the forced
// interface lines and the domain endpoints. Picks within `tol` of a forced line
// or of the previous kept pick are dropped: tol is a fraction of the target
// pitch, so the grading can go ~20× finer than uniform where the mesh is dense,
// but skin-band-scale node spacing (nm at GHz) doesn't leak into the plot grid.
function quantileAxis(vals, lo, hi, n, forced) {
    vals.sort();
    const tol = (hi - lo) / (n * 20);
    const lines = [];
    for (const v of [lo, hi, ...forced].sort((a, b) => a - b)) {
        if (v < lo - tol || v > hi + tol) continue;
        if (!lines.length || v - lines[lines.length - 1] > tol) lines.push(v);
    }
    const N = vals.length;
    const picks = [];
    for (let i = 0; i < n; i++) picks.push(vals[Math.round(i * (N - 1) / Math.max(n - 1, 1))]);
    picks.sort((a, b) => a - b);
    const out = lines.slice();
    let j = 0, lastPick = -Infinity;
    for (const p of picks) {
        if (p <= lo + tol || p >= hi - tol || p - lastPick <= tol) continue;
        while (j < lines.length && lines[j] < p) j++;
        const dNext = j < lines.length ? lines[j] - p : Infinity;
        const dPrev = j > 0 ? p - lines[j - 1] : Infinity;
        if (Math.min(dNext, dPrev) <= tol) continue;
        out.push(p);
        lastPick = p;
    }
    return Float64Array.from(out.sort((a, b) => a - b));
}

// Bucketed point-in-triangle locator over the mesh.
function buildLocator(mesh) {
    const { nodes, tris, nTris } = mesh;
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (let i = 0; i < mesh.nNodes; i++) {
        const x = nodes[2 * i], y = nodes[2 * i + 1];
        if (x < xmin) xmin = x; if (x > xmax) xmax = x;
        if (y < ymin) ymin = y; if (y > ymax) ymax = y;
    }
    const nb = Math.max(1, Math.round(Math.sqrt(nTris / 2)));
    const dx = (xmax - xmin) / nb || 1, dy = (ymax - ymin) / nb || 1;
    const buckets = Array.from({ length: nb * nb }, () => []);
    const bx = (x) => Math.min(nb - 1, Math.max(0, Math.floor((x - xmin) / dx)));
    const by = (y) => Math.min(nb - 1, Math.max(0, Math.floor((y - ymin) / dy)));
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const txmin = Math.min(nodes[2*v0], nodes[2*v1], nodes[2*v2]);
        const txmax = Math.max(nodes[2*v0], nodes[2*v1], nodes[2*v2]);
        const tymin = Math.min(nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]);
        const tymax = Math.max(nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]);
        for (let i = bx(txmin); i <= bx(txmax); i++)
            for (let j = by(tymin); j <= by(tymax); j++)
                buckets[i * nb + j].push(t);
    }
    // cache per-triangle coefficients lazily
    const coeffCache = new Array(nTris).fill(null);
    function coeffOf(t) {
        if (!coeffCache[t]) {
            const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
            coeffCache[t] = triCoefficients(nodes, v0, v1, v2).coeff;
        }
        return coeffCache[t];
    }
    function locate(x, y) {
        if (x < xmin || x > xmax || y < ymin || y > ymax) return -1;
        const list = buckets[bx(x) * nb + by(y)];
        for (const t of list) {
            const c = coeffOf(t);
            const l0 = c[0][0] + c[0][1] * x + c[0][2] * y;
            const l1 = c[1][0] + c[1][1] * x + c[1][2] * y;
            const l2 = c[2][0] + c[2][1] * x + c[2][2] * y;
            if (l0 >= -1e-9 && l1 >= -1e-9 && l2 >= -1e-9) return t;
        }
        return -1;
    }
    return { locate, coeffOf };
}

// Evaluate φ and ∇φ of a P2 static solution at point (x,y) in triangle t.
function evalPhi(phi, mesh, coeff, t, x, y) {
    const { tris, triEdges } = mesh;
    const v = [tris[3 * t], tris[3 * t + 1], tris[3 * t + 2]];
    let val = 0, gx = 0, gy = 0;
    for (let k = 0; k < 3; k++) {
        const pv = phi.phiVertex[v[k]];
        val += pv * lv(coeff, k, x, y);
        const g = lvGrad(coeff, k, x, y);
        gx += pv * g[0]; gy += pv * g[1];
    }
    for (let k = 0; k < 3; k++) {
        const pe = phi.phiEdge[triEdges[3 * t + k]];
        const [i, j] = EDGE_VERTS[k];
        val += pe * le(coeff, i, j, x, y);
        const g = leGrad(coeff, i, j, x, y);
        gx += pe * g[0]; gy += pe * g[1];
    }
    return { V: val, Ex: -gx, Ey: -gy };
}

// Resample a static solution onto a regular grid spanning `domain`.
// Returns { x:Float64Array(nx), y:Float64Array(ny), V, Ex, Ey } with V/Ex/Ey as
// [ny][nx] arrays (row index = y, matching plot.js / streamlines.js). Points
// outside the mesh (inside conductors, below ground) are left at 0.
export function resampleStatic(mesh, phi, domain, opts = {}) {
    // Sample on the caller-supplied grid (normally the mesh-derived graded grid from
    // buildGridFromMesh) when given, so the contour plot resolves thin conductors and
    // the ground surface; otherwise a uniform grid spanning the domain.
    const { x, y } = makeGrid(domain, opts);
    const nx = x.length, ny = y.length;

    // For a half-domain (symmetry) solve, mirror x<0 from the meshed x>0 half.
    // parity: 'even' → V even; 'odd' → V odd. E-field parity follows automatically
    // from differentiating the mirrored V grid below.
    const parity = opts.parity || null;
    const { locate, coeffOf } = buildLocator(mesh);
    // Sample only V (the P2 potential — continuous, and conductor interiors carry
    // their exact Dirichlet potential since they are meshed), then compute E = −∇V
    // by central differences. Evaluating ∇φ per P2 element instead (tried) leaves
    // the element-boundary gradient discontinuities in the data — faceted/kinked
    // |E| contours on coarse elements, and a hard one-sided jump along
    // material-interface rows where a centered stencil blends the two sides.
    const V = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ex = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ey = Array.from({ length: ny }, () => new Float64Array(nx));
    // Local element size at each sample — sets the differencing baseline below.
    const hT = Array.from({ length: ny }, () => new Float32Array(nx));
    // Deterministic side selection for on-edge samples (grid rows sit exactly on
    // mesh lines): nudge the LOCATE query by a tiny NE bias; V is continuous so
    // the evaluated value is side-independent, this only makes degenerate
    // on-vertex lookups deterministic.
    const eps = 1e-9 * Math.hypot(
        (domain.x_max ?? 0) - (domain.x_min ?? 0),
        (domain.y_max ?? 0) - (domain.y_min ?? 0)) || 1e-15;
    const { nodes, tris } = mesh;
    // CONTINUOUS local element-size field: node-averaged triangle size, interpolated
    // linearly within each triangle. Using the containing element's own size instead
    // (tried) makes the differencing baseline below jump at every element boundary,
    // which itself puts steps into E where the field has curvature.
    const nodeH = new Float64Array(mesh.nNodes), nodeW = new Float64Array(mesh.nNodes);
    for (let t = 0; t < mesh.nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const area = Math.abs(
            (nodes[2*v1] - nodes[2*v0]) * (nodes[2*v2+1] - nodes[2*v0+1]) -
            (nodes[2*v2] - nodes[2*v0]) * (nodes[2*v1+1] - nodes[2*v0+1])) / 2;
        const sz = Math.sqrt(2 * area) || 1e-30;
        for (const v of [v0, v1, v2]) { nodeH[v] += sz; nodeW[v]++; }
    }
    for (let i = 0; i < mesh.nNodes; i++) if (nodeW[i] > 0) nodeH[i] /= nodeW[i];
    const sizeAt = (t, qx, qy) => {
        const c = coeffOf(t);
        const l0 = c[0][0] + c[0][1] * qx + c[0][2] * qy;
        const l1 = c[1][0] + c[1][1] * qx + c[1][2] * qy;
        const l2 = c[2][0] + c[2][1] * qx + c[2][2] * qy;
        return l0 * nodeH[tris[3 * t]] + l1 * nodeH[tris[3 * t + 1]] + l2 * nodeH[tris[3 * t + 2]];
    };
    // Point sample of V anywhere in the (mirrored) domain; NaN outside the mesh.
    const sampleV = (qx, qy) => {
        let s = 1;
        if (parity && qx < 0) {
            qx = -qx;
            if (parity === 'odd') s = -1;
        }
        let t = locate(qx + eps, qy + eps);
        if (t < 0) t = locate(qx, qy);
        if (t < 0) return NaN;
        return s * evalPhi(phi, mesh, coeffOf(t), t, qx, qy).V;
    };
    for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
            let qx = x[i], sV = 1;
            if (parity && qx < 0) {
                qx = -qx;
                if (parity === 'odd') sV = -1;
            }
            let t = locate(qx + eps, y[j] + eps);
            if (t < 0) t = locate(qx, y[j]);   // domain edge: fall back to the exact point
            if (t < 0) continue;
            const r = evalPhi(phi, mesh, coeffOf(t), t, qx, y[j]);
            V[j][i] = sV * r.V;
            hT[j][i] = sizeAt(t, qx, y[j]);
        }
    }
    // E = −∇V by central differences over an ELEMENT-SIZE-AWARE baseline: half-width
    // max(local grid pitch, half the containing element's size). Where the mesh is
    // finer than the grid (near conductors) this is the plain grid-pitch central
    // difference. Where the mesh is coarser (the far field — the tensor grid carries
    // the trace region's µm-fine lines all the way out), a grid-pitch baseline would
    // sample the P2 gradient jump across each element edge as a sharp step, and the
    // |E| contours come out jagged wherever they cross element boundaries; widening
    // the baseline to the element scale averages the jump away, which is exactly the
    // resolution the FEM solution actually has there. The baseline SHRINKS
    // CONTINUOUSLY as a conductor surface or domain edge approaches — clamped to the
    // distance so an endpoint lands at most ON the surface (a valid Dirichlet
    // sample), never across it into the interior V plateau (which would bias the
    // difference low right where the field is largest). A continuous clamp, not a
    // switch to another estimator: a hard hand-off to the grid stencil near
    // surfaces (tried) leaves a small systematic |E| step along the hand-off band,
    // which drags contour lines sideways just before they reach the ground plane
    // instead of letting them meet it straight. Only when the clamp collapses the
    // baseline below the local grid pitch (points ON a surface line, or inside a
    // conductor next to its face) does it fall back to the plain non-uniform grid
    // stencil. Boundary rows/columns stay 0.
    const x0 = domain.x_min ?? x[0], x1 = domain.x_max ?? x[nx - 1];
    const y0 = domain.y_min ?? y[0], y1 = domain.y_max ?? y[ny - 1];
    // Conductor rects in PLOT (full-domain) coordinates. Prefer the caller-supplied
    // list (tri_backend passes the solver's conductor geometry — it includes the
    // ground planes, which mesh.condRect.rects does NOT: grounds are handled as PEC
    // walls there, and without them the baseline would reach into the ground metal's
    // V≡0 region and bias |E| low across the whole near-ground band). Fall back to
    // the mesh rects (mirrored for a half-domain solve) when none are given.
    let rects = opts.rects;
    if (!rects) {
        rects = [];
        for (const c of (mesh.condRect && mesh.condRect.rects) || []) {
            rects.push(c);
            if (parity) rects.push({ xmin: -c.xmax, xmax: -c.xmin, ymin: c.ymin, ymax: c.ymax });
        }
    }
    const tolC = eps;
    for (let j = 1; j < ny - 1; j++) {
        const dyd = y[j] - y[j - 1], dyu = y[j + 1] - y[j];
        for (let i = 1; i < nx - 1; i++) {
            const h = hT[j][i];
            if (!h) continue;   // outside the mesh
            const px = x[i], py = y[j];
            const dxl = px - x[i - 1], dxr = x[i + 1] - px;
            let bx = Math.max((dxl + dxr) / 2, 0.5 * h);
            let by = Math.max((dyd + dyu) / 2, 0.5 * h);
            bx = Math.min(bx, px - x0, x1 - px);
            by = Math.min(by, py - y0, y1 - py);
            // Baseline vs conductor faces. Center inside a conductor: clamp the
            // baseline to stay inside (interior plateau, E→0). Otherwise, when the
            // baseline would cross the NEAREST face, don't clamp — MIRROR the sample
            // across it (image theory: |E| is even about a flat Dirichlet face, so
            // the odd extension Ṽ(face−d) = 2·V_face − V(face+d) is the correct
            // continuation). This keeps the SAME baseline width for every row near
            // the face, so the E estimate stays consistent from the bulk all the way
            // onto the surface row and contours meet the ground at 90° — clamping or
            // switching estimators near the face (tried) leaves a percent-level |E|
            // step across the hand-off band that visibly bends the contours sideways
            // just above the plane.
            let inside = false;
            let mxPlane = null, mxSide = 0, mxDist = Infinity;   // nearest x-face the baseline crosses
            let myPlane = null, mySide = 0, myDist = Infinity;   // nearest y-face the baseline crosses
            for (const c of rects) {
                if (c.shape) {
                    // Curved conductor: clamp the baseline inside it as for a rect, but
                    // do NOT set up a mirror plane. The image-theory extension below is
                    // derived for a FLAT Dirichlet face and does not hold on a curved
                    // one; near a coax surface the grid is mesh-quantile dense anyway,
                    // so the plain non-uniform stencil is accurate there.
                    if (shapeContains(c, px, py, -tolC)) {
                        const d = distToShapeBoundary(c, px, py);
                        bx = Math.min(bx, d); by = Math.min(by, d);
                        inside = true;
                    }
                    continue;
                }
                const inX = px > c.xmin + tolC && px < c.xmax - tolC;
                const inY = py > c.ymin + tolC && py < c.ymax - tolC;
                if (inX && inY) {
                    bx = Math.min(bx, px - c.xmin, c.xmax - px);
                    by = Math.min(by, py - c.ymin, c.ymax - py);
                    inside = true;
                } else if (inY) {
                    // conductor beside the point; side = which side of the face the metal is on
                    const left = px <= c.xmin + tolC;   // metal to the right of its xmin face
                    const plane = left ? c.xmin : c.xmax;
                    const d = Math.abs(px - plane);
                    if (d < bx && d < mxDist) { mxDist = d; mxPlane = plane; mxSide = left ? 1 : -1; }
                } else if (inX) {
                    const below = py >= c.ymax - tolC;  // metal below its ymax face
                    const plane = below ? c.ymax : c.ymin;
                    const d = Math.abs(py - plane);
                    if (d < by && d < myDist) { myDist = d; myPlane = plane; mySide = below ? -1 : 1; }
                }
            }
            // Mirror any sample that falls on the METAL side of the nearest face
            // (including when the center itself sits exactly on the face).
            const sampleX = (qx) => {
                if (!inside && mxPlane !== null && (qx - mxPlane) * mxSide > 0) {
                    return 2 * sampleV(mxPlane, py) - sampleV(2 * mxPlane - qx, py);
                }
                return sampleV(qx, py);
            };
            const sampleY = (qy) => {
                if (!inside && myPlane !== null && (qy - myPlane) * mySide > 0) {
                    return 2 * sampleV(px, myPlane) - sampleV(px, 2 * myPlane - qy);
                }
                return sampleV(px, qy);
            };
            let e;
            if (bx > 0.5 * Math.min(dxl, dxr) &&
                isFinite(e = -(sampleX(px + bx) - sampleX(px - bx)) / (2 * bx))) {
                Ex[j][i] = e;
            } else {
                Ex[j][i] = -(
                    (dxl / (dxr * (dxl + dxr))) * V[j][i + 1] +
                    ((dxr - dxl) / (dxl * dxr)) * V[j][i] -
                    (dxr / (dxl * (dxl + dxr))) * V[j][i - 1]
                );
            }
            if (by > 0.5 * Math.min(dyd, dyu) &&
                isFinite(e = -(sampleY(py + by) - sampleY(py - by)) / (2 * by))) {
                Ey[j][i] = e;
            } else {
                Ey[j][i] = -(
                    (dyd / (dyu * (dyd + dyu))) * V[j + 1][i] +
                    ((dyu - dyd) / (dyd * dyu)) * V[j][i] -
                    (dyu / (dyd * (dyd + dyu))) * V[j - 1][i]
                );
            }
        }
    }
    return { x, y, V, Ex, Ey };
}

// Per-triangle material region for the nodal recovery. Averaging must never cross a
// dielectric interface (the normal E component genuinely jumps by the permittivity
// ratio) or a conductor surface (the field is identically zero inside the metal):
// cross-region averaging smears the physical discontinuity over one element, which at
// wavelength-scale bulk element sizes renders as jagged element-shaped blobs along the
// substrate and around the trace. Triangles are grouped by their epsMap permittivity,
// with conductor interiors (centroid inside a conductor rect) as their own region.
function buildTriRegions(mesh) {
    const { nodes, tris, nTris, epsMap, condRect } = mesh;
    const regionOf = new Int32Array(nTris);
    if (!epsMap || epsMap.length !== nTris) return { regionOf, nRegions: 1 };
    const rects = (condRect && condRect.rects) || [];
    const ids = new Map();
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const xc = (nodes[2*v0] + nodes[2*v1] + nodes[2*v2]) / 3;
        const yc = (nodes[2*v0+1] + nodes[2*v1+1] + nodes[2*v2+1]) / 3;
        let key = 'cond';
        let inCond = false;
        for (const r of rects) {
            if (shapeContains(r, xc, yc, 0)) { inCond = true; break; }
        }
        if (!inCond) {
            const e = epsMap[t];
            key = e.re.toPrecision(6) + ',' + (e.im ? e.im.toPrecision(6) : '0');
        }
        let id = ids.get(key);
        if (id === undefined) { id = ids.size; ids.set(key, id); }
        regionOf[t] = id;
    }
    return { regionOf, nRegions: ids.size };
}

// Recover a nodal transverse field, CONTINUOUS within each material region, by
// area-weighted averaging of the per-element Nedelec field at each node. The Lee–Jin
// edge-element field e_t is DISCONTINUOUS across triangle edges, so sampling it directly
// at grid points gives faceted/discontinuity artifacts (worst near the trace, where the
// mesh is fine and the field is strong). Averaging the four complex components
// (Ex, Ey × re, im) to nodes and interpolating linearly within each triangle yields a
// smooth field — the same recovery resampleStatic used for the static E-field before it
// switched to FD-of-V. The averages are kept PER REGION (node slot = node·nRegions +
// region), so the genuine field jumps at dielectric interfaces and conductor surfaces
// stay sharp; a triangle always contributes to its own region's slots, so every
// (vertex, region-of-containing-triangle) slot an interpolation reads is populated.
function recoverNodalModeField(mesh, fm, vRe, vIm) {
    const { nodes, tris, nTris, nNodes } = mesh;
    const { regionOf, nRegions } = buildTriRegions(mesh);
    const nSlots = nNodes * nRegions;
    const exr = new Float64Array(nSlots), exi = new Float64Array(nSlots);
    const eyr = new Float64Array(nSlots), eyi = new Float64Array(nSlots);
    const w = new Float64Array(nSlots);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const ax = nodes[2*v0], ay = nodes[2*v0+1], bx = nodes[2*v1], by = nodes[2*v1+1], cx = nodes[2*v2], cy = nodes[2*v2+1];
        const area = Math.abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) / 2;
        if (!(area > 0)) continue;
        const reg = regionOf[t];
        for (const vk of [v0, v1, v2]) {
            const f = evalFieldsAtPoint(t, nodes[2 * vk], nodes[2 * vk + 1], mesh, fm, vRe, vIm);
            const s = vk * nRegions + reg;
            exr[s] += area * f.exr; exi[s] += area * f.exi;
            eyr[s] += area * f.eyr; eyi[s] += area * f.eyi; w[s] += area;
        }
    }
    for (let i = 0; i < nSlots; i++) {
        if (w[i] > 0) { exr[i] /= w[i]; exi[i] /= w[i]; eyr[i] /= w[i]; eyi[i] /= w[i]; }
    }
    return { exr, exi, eyr, eyi, regionOf, nRegions };
}

// Resample a full-wave eigenmode's transverse E-field onto a regular grid.
// The Lee–Jin eigenvector stores e_t = γ·Et, so the spatial pattern of |e_t| is the
// transverse field pattern (the constant γ only scales the whole map), which is what a
// mode plot shows. Returns { x, y, E[ny][nx], Ex[ny][nx], Ey[ny][nx] }: E is the
// transverse magnitude |Et| (for the heatmap); Ex/Ey are the real-part components (for
// an optional quiver/streamline overlay). Points outside the mesh are left at 0.
export function resampleModeField(mesh, fm, vRe, vIm, domain, opts = {}) {
    const { x, y } = makeGrid(domain, opts);
    const nx = x.length, ny = y.length;
    const { locate, coeffOf } = buildLocator(mesh);
    const nodal = recoverNodalModeField(mesh, fm, vRe, vIm);   // continuous field
    const { tris } = mesh;
    const E = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ex = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ey = Array.from({ length: ny }, () => new Float64Array(nx));
    for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
            const t = locate(x[i], y[j]);
            if (t < 0) continue;
            const coeff = coeffOf(t);
            const l0 = coeff[0][0] + coeff[0][1] * x[i] + coeff[0][2] * y[j];
            const l1 = coeff[1][0] + coeff[1][1] * x[i] + coeff[1][2] * y[j];
            const l2 = coeff[2][0] + coeff[2][1] * x[i] + coeff[2][2] * y[j];
            // Barycentric interpolation of the recovered nodal field, reading the nodal
            // slots of the containing triangle's material region → continuous within a
            // region, sharp at interfaces.
            const nr = nodal.nRegions, reg = nodal.regionOf[t];
            const s0 = tris[3 * t] * nr + reg, s1 = tris[3 * t + 1] * nr + reg, s2 = tris[3 * t + 2] * nr + reg;
            const exr = l0 * nodal.exr[s0] + l1 * nodal.exr[s1] + l2 * nodal.exr[s2];
            const exi = l0 * nodal.exi[s0] + l1 * nodal.exi[s1] + l2 * nodal.exi[s2];
            const eyr = l0 * nodal.eyr[s0] + l1 * nodal.eyr[s1] + l2 * nodal.eyr[s2];
            const eyi = l0 * nodal.eyi[s0] + l1 * nodal.eyi[s1] + l2 * nodal.eyi[s2];
            E[j][i] = Math.hypot(Math.hypot(exr, exi), Math.hypot(eyr, eyi));
            Ex[j][i] = exr; Ey[j][i] = eyr;
        }
    }
    return { x, y, E, Ex, Ey };
}
