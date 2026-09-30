// Resample a triangular-mesh static FEM solution (P2 scalar potential) onto a
// regular grid, producing the SAME { x, y, V[ny][nx], Ex[ny][nx], Ey[ny][nx] }
// shape the rectilinear FDM solver exposes — so plot.js (heatmap/contour) and
// streamlines.js work unchanged on the triangular backend's output.
//
// Also the field on the triangles themselves for the |E| image (static, full-wave mode
// and Modes tab), the full-wave mode field on a grid, and the MQS current density and
// surface current plots.

import { triCoefficients } from './tri_fem.js';
import { shapeContains, distToShapeBoundary } from '../shapes.js';
import { evalFieldsAtPoint } from './tri_ms_solver.js';
import { segmentBuffer } from '../surface_segments.js';


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
    // The domain endpoints are always lines. A forced line outside the domain (the
    // companion line of a conductor face on a wall) or within tol of an endpoint is
    // dropped: a sample outside the mesh reads 0.
    const lines = [lo];
    for (const v of forced.slice().sort((a, b) => a - b)) {
        if (v <= lo + tol || v >= hi - tol) continue;
        if (v - lines[lines.length - 1] > tol) lines.push(v);
    }
    lines.push(hi);
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
    // The adaptively refined triangles crowd the few buckets around the conductors,
    // exactly where the graded plot grid samples most. A bucket holding many triangles
    // gets a sub-grid, each sub-cell listing (in the same ascending order) the bucket's
    // triangles whose bounding box, widened past the containment tolerance below,
    // reaches it. Every triangle that can pass the test at a point is in that point's
    // sub-cell, so the first hit is the same triangle the whole bucket would give.
    const SUB_MIN = 24;
    const subs = new Array(nb * nb).fill(null);
    for (let k = 0; k < nb * nb; k++) {
        const list = buckets[k];
        if (list.length <= SUB_MIN) continue;
        const i = Math.floor(k / nb), j = k % nb;
        const cx0 = xmin + i * dx, cy0 = ymin + j * dy;
        const m = Math.ceil(Math.sqrt(list.length / 4));
        const sdx = dx / m, sdy = dy / m;
        const cells = Array.from({ length: m * m }, () => []);
        const sx = v => Math.min(m - 1, Math.max(0, Math.floor((v - cx0) / sdx)));
        const sy = v => Math.min(m - 1, Math.max(0, Math.floor((v - cy0) / sdy)));
        for (const t of list) {
            const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
            const x0 = Math.min(nodes[2*v0], nodes[2*v1], nodes[2*v2]), x1 = Math.max(nodes[2*v0], nodes[2*v1], nodes[2*v2]);
            const y0 = Math.min(nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]), y1 = Math.max(nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]);
            const pad = 1e-6 * Math.hypot(x1 - x0, y1 - y0);
            for (let a = sx(x0 - pad); a <= sx(x1 + pad); a++)
                for (let b = sy(y0 - pad); b <= sy(y1 + pad); b++) cells[a * m + b].push(t);
        }
        subs[k] = { cells, sx, sy, m };
    }
    function locate(x, y) {
        if (x < xmin || x > xmax || y < ymin || y > ymax) return -1;
        const k = bx(x) * nb + by(y);
        const sub = subs[k];
        const list = sub ? sub.cells[sub.sx(x) * sub.m + sub.sy(y)] : buckets[k];
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

// φ of a P2 static solution at (x, y) in triangle t: the vertex bases -λ + 2λ² and the
// edge bases 4 λi λj, without the gradient (E comes from differences of V). The
// resampler samples V several times per grid point, so this is its inner loop.
function evalPhiValue(phi, mesh, coeff, t, x, y) {
    const { tris, triEdges } = mesh;
    const c0 = coeff[0], c1 = coeff[1], c2 = coeff[2];
    const l0 = c0[0] + c0[1]*x + c0[2]*y;
    const l1 = c1[0] + c1[1]*x + c1[2]*y;
    const l2 = c2[0] + c2[1]*x + c2[2]*y;
    const pv = phi.phiVertex, pe = phi.phiEdge, e = 3 * t;
    let val = 0;
    val += pv[tris[e]] * (-l0 + 2*l0*l0);
    val += pv[tris[e + 1]] * (-l1 + 2*l1*l1);
    val += pv[tris[e + 2]] * (-l2 + 2*l2*l2);
    val += pe[triEdges[e]] * (4*l0*l1);
    val += pe[triEdges[e + 1]] * (4*l1*l2);
    val += pe[triEdges[e + 2]] * (4*l2*l0);
    return val;
}

// -dV/dx on the non-uniform grid stencil (vm, v0, vp at spacings dl, dr). A neighbor
// outside the mesh (an unmeshed conductor interior) holds no potential, so the
// difference goes one-sided away from it; with both neighbors outside it is 0.
function gridDiff(vm, v0, vp, dl, dr, okM, okP) {
    if (okM && okP) return -((dl / (dr * (dl + dr))) * vp + ((dr - dl) / (dl * dr)) * v0 - (dr / (dl * (dl + dr))) * vm);
    if (okP) return -(vp - v0) / dr;
    if (okM) return -(v0 - vm) / dl;
    return 0;
}

// Point evaluator of the smoothed static field E = -grad V of a P2 solve `phi` (plot
// coordinates, a half-domain solve mirrored by opts.parity as in resampleStatic):
// evalE(px, py, h, dxl, dxr, dyd, dyu) differences V over an element-size baseline
// (h the local element size, d* the grid pitches around the point, 0 off a grid) and
// returns [Ex, Ey], a component NaN where no baseline fits (the caller falls back).
// opts.oneSidedEdges: at the domain edge the stencil turns one-sided instead of
// shrinking (a mesh point on the edge has no grid neighbour to fall back to).
function staticFieldEvaluator(mesh, phi, domain, opts = {}) {
    const parity = opts.parity || null;
    const { locate, coeffOf } = buildLocator(mesh);
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
    // Material region of the point (buildTriRegions), -1 outside the mesh.
    const { regionOf, condRegion } = buildTriRegions(mesh);
    const regionAt = (qx, qy) => {
        if (parity && qx < 0) qx = -qx;
        let t = locate(qx + eps, qy + eps);
        if (t < 0) t = locate(qx, qy);
        return t < 0 ? -1 : regionOf[t];
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
        return s * evalPhiValue(phi, mesh, coeffOf(t), t, qx, qy);
    };
    const x0 = domain.x_min, x1 = domain.x_max;
    const y0 = domain.y_min, y1 = domain.y_max;
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
    // The end of the stencil from p toward q, pulled back to the first point where the
    // line leaves dielectric region reg for another dielectric, and with toMetal also
    // for a conductor or the edge of the mesh (a shaped conductor, a shield): a scan in
    // eight steps, then bisection.
    const interfaceToward = (regionOn, p, q, reg, toMetal) => {
        const other = r => r !== reg && (toMetal || (r >= 0 && r !== condRegion));
        let a = p;
        for (let k = 1; k <= 8; k++) {
            let b = p + (q - p) * k / 8;
            if (other(regionOn(b))) {
                for (let it = 0; it < 30; it++) {
                    const m = (a + b) / 2;
                    if (other(regionOn(m))) b = m; else a = m;
                }
                return a;
            }
            a = b;
        }
        return q;
    };
    const oneSided = !!opts.oneSidedEdges;
    const evalE = (px, py, h, dxl, dxr, dyd, dyu) => {
        let bx = Math.max((dxl + dxr) / 2, 0.5 * h);
        let by = Math.max((dyd + dyu) / 2, 0.5 * h);
        if (!oneSided) {
            bx = Math.min(bx, px - x0, x1 - px);
            by = Math.min(by, py - y0, y1 - py);
        }
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
        // A point the mesh puts in a dielectric is outside every shaped conductor: the
        // exact polygon test (a shield ring's loops) only runs for metal and off-mesh points.
        const regP = regionAt(px, py);
        const maybeInShape = regP < 0 || regP === condRegion;
        for (const c of rects) {
            if (c.shape) {
                if (!maybeInShape) continue;
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
        // The normal E jumps at a dielectric interface, so the baseline must not blend
        // the two sides: an end in another dielectric is pulled back to the interface
        // (V is continuous there) and the derivative taken from the quadratic through
        // the two ends and their midpoint, all on the sample's own side. An end in a
        // conductor with no mirror plane (a curved or slanted surface) is pulled back
        // to its surface the same way, the metal's constant V would bias E low.
        const reg = inside ? -1 : regP;
        const derivative = (sample, p, b, at, pitch, mirrored, edgeLo, edgeHi) => {
            let lo = p - b, hi = p + b;
            if (oneSided) { lo = Math.max(lo, edgeLo); hi = Math.min(hi, edgeHi); }
            if (reg >= 0 && reg !== condRegion) {
                lo = interfaceToward(at, p, lo, reg, !mirrored);
                hi = interfaceToward(at, p, hi, reg, !mirrored);
            }
            if (lo === p - b && hi === p + b) return -(sample(hi) - sample(lo)) / (2 * b);
            if (hi - lo < 0.5 * pitch || !(hi - lo > 1e-6 * b)) return NaN;
            const m = (lo + hi) / 2;
            return -(sample(lo) * (2 * p - m - hi) / ((lo - m) * (lo - hi))
                + sample(m) * (2 * p - lo - hi) / ((m - lo) * (m - hi))
                + sample(hi) * (2 * p - lo - m) / ((hi - lo) * (hi - m)));
        };
        const pitchX = Math.min(dxl, dxr), pitchY = Math.min(dyd, dyu);
        const ex = bx > 0.5 * pitchX
            ? derivative(sampleX, px, bx, q => regionAt(q, py), pitchX, mxPlane !== null, x0, x1) : NaN;
        const ey = by > 0.5 * pitchY
            ? derivative(sampleY, py, by, q => regionAt(px, q), pitchY, myPlane !== null, y0, y1) : NaN;
        return [ex, ey];
    };
    return { locate, coeffOf, sampleV, sizeAt, regionOf, condRegion, regionAt, evalE, eps };
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
    const parity = opts.parity || null;
    const ev = staticFieldEvaluator(mesh, phi, { x_min: domain.x_min ?? x[0], x_max: domain.x_max ?? x[nx - 1],
        y_min: domain.y_min ?? y[0], y_max: domain.y_max ?? y[ny - 1] }, opts);
    const { locate, coeffOf, sizeAt, eps } = ev;
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
            V[j][i] = sV * evalPhiValue(phi, mesh, coeffOf(t), t, qx, y[j]);
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
    for (let j = 1; j < ny - 1; j++) {
        const dyd = y[j] - y[j - 1], dyu = y[j + 1] - y[j];
        for (let i = 1; i < nx - 1; i++) {
            const h = hT[j][i];
            if (!h) continue;   // outside the mesh
            const dxl = x[i] - x[i - 1], dxr = x[i + 1] - x[i];
            const [ex, ey] = ev.evalE(x[i], y[j], h, dxl, dxr, dyd, dyu);
            Ex[j][i] = isFinite(ex) ? ex
                : gridDiff(V[j][i - 1], V[j][i], V[j][i + 1], dxl, dxr, hT[j][i - 1] > 0, hT[j][i + 1] > 0);
            Ey[j][i] = isFinite(ey) ? ey
                : gridDiff(V[j - 1][i], V[j][i], V[j + 1][i], dyd, dyu, hT[j - 1][i] > 0, hT[j + 1][i] > 0);
        }
    }
    return { x, y, V, Ex, Ey };
}

// -grad of the P2 potential phi in triangle t at (x, y).
function evalPhiField(phi, mesh, coeff, t, x, y) {
    const { tris, triEdges } = mesh;
    const pv = phi.phiVertex, pe = phi.phiEdge;
    const l = [0, 1, 2].map(a => coeff[a][0] + coeff[a][1] * x + coeff[a][2] * y);
    let gx = 0, gy = 0;
    for (let a = 0; a < 3; a++) {
        const b = (a + 1) % 3;
        const fv = pv[tris[3 * t + a]] * (4 * l[a] - 1), fe = 4 * pe[triEdges[3 * t + a]];
        gx += fv * coeff[a][1] + fe * (l[a] * coeff[b][1] + l[b] * coeff[a][1]);
        gy += fv * coeff[a][2] + fe * (l[a] * coeff[b][2] + l[b] * coeff[a][2]);
    }
    return [-gx, -gy];
}

// The smoothed static field (staticFieldEvaluator, as on the plot grid) at the three
// vertices and three edge midpoints of every dielectric triangle, each evaluated on the
// triangle's own side of an interface. Point order per triangle: v0, v1, v2, then the
// midpoints of v0-v1, v1-v2, v2-v0. Returns { tri, Ex, Ey } (6 values per triangle).
export function staticFieldOnMesh(mesh, phi, domain, opts = {}) {
    const ev = staticFieldEvaluator(mesh, phi, domain, { ...opts, oneSidedEdges: true });
    const { nodes, tris, nTris } = mesh;
    const tri = [];
    for (let t = 0; t < nTris; t++) if (ev.regionOf[t] !== ev.condRegion) tri.push(t);
    const Ex = new Float64Array(6 * tri.length), Ey = new Float64Array(6 * tri.length);
    const px = new Float64Array(6), py = new Float64Array(6);
    tri.forEach((t, k) => {
        const v = [tris[3 * t], tris[3 * t + 1], tris[3 * t + 2]];
        for (let a = 0; a < 3; a++) {
            const b = (a + 1) % 3;
            px[a] = nodes[2 * v[a]]; py[a] = nodes[2 * v[a] + 1];
            px[3 + a] = (nodes[2 * v[a]] + nodes[2 * v[b]]) / 2;
            py[3 + a] = (nodes[2 * v[a] + 1] + nodes[2 * v[b] + 1]) / 2;
        }
        const cx = (px[0] + px[1] + px[2]) / 3, cy = (py[0] + py[1] + py[2]) / 3;
        const coeff = ev.coeffOf(t);
        for (let q = 0; q < 6; q++) {
            // A hair toward the centroid, so a point on an interface reads this side.
            const qx = px[q] + 1e-4 * (cx - px[q]), qy = py[q] + 1e-4 * (cy - py[q]);
            let [ex, ey] = ev.evalE(qx, qy, ev.sizeAt(t, qx, qy), 0, 0, 0, 0);
            if (!isFinite(ex) || !isFinite(ey)) {
                const g = evalPhiField(phi, mesh, coeff, t, qx, qy);
                if (!isFinite(ex)) ex = g[0];
                if (!isFinite(ey)) ey = g[1];
            }
            Ex[6 * k + q] = ex; Ey[6 * k + q] = ey;
        }
    });
    return { tri: Int32Array.from(tri), Ex, Ey };
}

// |E| on the dielectric triangles for a plot: each triangle of staticFieldOnMesh split in
// four at its edge midpoints, { tris: [x0, y0, x1, y1, x2, y2, ...], E: [E0, E1, E2, ...] }.
// combine(t, q, Ex, Ey) maps the static field at point q of triangle t to the plotted one
// ([Ex, Ey]); mirror adds the image about x = 0 of a half-domain solve.
//
// Each triangle's field is its own, so a vertex or edge midpoint shared by triangles would
// read a different |E| in each (by tens of percent at a conductor corner, where the field
// is singular) and a contour line would end on the shared edge. The plotted value is the
// mean over the triangles of the same material sharing the point: continuous within a
// material, sharp across an interface.
export function meshFieldBlock(mesh, sf, combine = null, mirror = false) {
    const { tris, triEdges, nNodes, nEdges } = mesh;
    const { regionOf, nRegions: nR } = buildTriRegions(mesh);
    const n = sf.tri.length;
    const sum = new Float64Array((nNodes + nEdges) * nR), cnt = new Uint32Array((nNodes + nEdges) * nR);
    // Slot of point q (vertices 0-2, then edge midpoints) of triangle t.
    const slot = (t, q) => (q < 3 ? tris[3 * t + q] : nNodes + triEdges[3 * t + q - 3]) * nR + regionOf[t];
    for (let k = 0; k < n; k++) {
        const t = sf.tri[k];
        for (let q = 0; q < 6; q++) {
            const e = combine ? combine(t, q, sf.Ex[6 * k + q], sf.Ey[6 * k + q]) : [sf.Ex[6 * k + q], sf.Ey[6 * k + q]];
            const i = slot(t, q);
            sum[i] += Math.hypot(e[0], e[1]); cnt[i]++;
        }
    }
    return splitTriBlock(mesh, sf.tri, (k, t, mag) => {
        for (let q = 0; q < 6; q++) { const i = slot(t, q); mag[q] = sum[i] / cnt[i]; }
    }, mirror);
}

// Triangles `tri` of the mesh each split in four at the edge midpoints, with a value at
// the three vertices and three edge midpoints (order v0, v1, v2, then v0-v1, v1-v2,
// v2-v0) written by fill(k, t, mag) for the k-th triangle t. { tris: [x0, y0, x1, y1,
// x2, y2, ...], E: [E0, E1, E2, ...] }; mirror adds the image about x = 0.
function splitTriBlock(mesh, tri, fill, mirror = false) {
    const { nodes, tris } = mesh;
    const SUB = [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]];
    const n = tri.length, copies = mirror ? 2 : 1;
    const T = new Float64Array(24 * n * copies), E = new Float64Array(12 * n * copies);
    const px = new Float64Array(6), py = new Float64Array(6), mag = new Float64Array(6);
    for (let k = 0; k < n; k++) {
        const t = tri[k];
        for (let a = 0; a < 3; a++) {
            const va = tris[3 * t + a], vb = tris[3 * t + (a + 1) % 3];
            px[a] = nodes[2 * va]; py[a] = nodes[2 * va + 1];
            px[3 + a] = (nodes[2 * va] + nodes[2 * vb]) / 2; py[3 + a] = (nodes[2 * va + 1] + nodes[2 * vb + 1]) / 2;
        }
        fill(k, t, mag);
        for (let c = 0; c < copies; c++) {
            const sx = c ? -1 : 1, base = (c * n + k) * 4;
            SUB.forEach((sub, u) => {
                for (let a = 0; a < 3; a++) {
                    T[6 * (base + u) + 2 * a] = sx * px[sub[a]];
                    T[6 * (base + u) + 2 * a + 1] = py[sub[a]];
                    E[3 * (base + u) + a] = mag[sub[a]];
                }
            });
        }
    }
    return { tris: T, E };
}

// Per-triangle material region for the nodal recovery. Averaging must never cross a
// dielectric interface (the normal E component genuinely jumps by the permittivity
// ratio) or a conductor surface (the field is identically zero inside the metal):
// cross-region averaging smears the physical discontinuity over one element, which at
// wavelength-scale bulk element sizes renders as jagged element-shaped blobs along the
// substrate and around the trace. Triangles are grouped by their epsMap permittivity,
// with conductor interiors (centroid inside a conductor rect) as their own region.
const regionCache = new WeakMap();
export function buildTriRegions(mesh) {
    // Cached per mesh and material map: the conductor test is a polygon test per
    // triangle, and every plot field of a solve needs the same regions.
    const hit = regionCache.get(mesh);
    if (hit && hit.epsMap === mesh.epsMap) return hit.value;
    const value = triRegions(mesh);
    regionCache.set(mesh, { epsMap: mesh.epsMap, value });
    return value;
}

function triRegions(mesh) {
    const { nodes, tris, nTris, epsMap, condRect } = mesh;
    const regionOf = new Int32Array(nTris);
    if (!epsMap || epsMap.length !== nTris) return { regionOf, nRegions: 1, condRegion: -1 };
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
    return { regionOf, nRegions: ids.size, condRegion: ids.has('cond') ? ids.get('cond') : -1 };
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
const nodalCache = new WeakMap();
export function recoverNodalModeField(mesh, fm, vRe, vIm) {
    // Cached per eigenvector: the E plot recovers each mode field on the grid and on
    // the mesh.
    const hit = nodalCache.get(vRe);
    if (hit && hit.mesh === mesh && hit.vIm === vIm) return hit.value;
    const value = nodalModeField(mesh, fm, vRe, vIm);
    nodalCache.set(vRe, { mesh, vIm, value });
    return value;
}

function nodalModeField(mesh, fm, vRe, vIm) {
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

// |E_t| of an eigenmode on the dielectric triangles for a plot, in meshFieldBlock's
// format { tris, E }: the recovered nodal field (recoverNodalModeField) at the vertices
// and the element field averaged per material region at the edge midpoints, each
// triangle split in four at the midpoints. The magnitude is complex (|Re|^2 + |Im|^2),
// so the eigenvector's arbitrary phase drops out.
export function modeFieldMeshBlock(mesh, fm, vRe, vIm) {
    const { nodes, tris, nTris, triEdges, nEdges } = mesh;
    const nod = recoverNodalModeField(mesh, fm, vRe, vIm);
    const { regionOf, nRegions: nR } = nod;
    const condRegion = buildTriRegions(mesh).condRegion;
    // Edge midpoint field per (edge, region), averaged over the triangles sharing it.
    const mid = new Float64Array(5 * nEdges * nR);
    for (let t = 0; t < nTris; t++) {
        const reg = regionOf[t];
        if (reg === condRegion) continue;
        for (let a = 0; a < 3; a++) {
            const va = tris[3 * t + a], vb = tris[3 * t + (a + 1) % 3];
            const f = evalFieldsAtPoint(t, (nodes[2 * va] + nodes[2 * vb]) / 2,
                (nodes[2 * va + 1] + nodes[2 * vb + 1]) / 2, mesh, fm, vRe, vIm);
            const s = 5 * (triEdges[3 * t + a] * nR + reg);
            mid[s] += f.exr; mid[s + 1] += f.exi; mid[s + 2] += f.eyr; mid[s + 3] += f.eyi; mid[s + 4]++;
        }
    }
    const tri = [];
    for (let t = 0; t < nTris; t++) if (regionOf[t] !== condRegion) tri.push(t);
    return splitTriBlock(mesh, tri, (k, t, mag) => {
        const reg = regionOf[t];
        for (let a = 0; a < 3; a++) {
            const s = tris[3 * t + a] * nR + reg;
            mag[a] = Math.hypot(nod.exr[s], nod.exi[s], nod.eyr[s], nod.eyi[s]);
            const m = 5 * (triEdges[3 * t + a] * nR + reg), c = mid[m + 4] || 1;
            mag[3 + a] = Math.hypot(mid[m], mid[m + 1], mid[m + 2], mid[m + 3]) / c;
        }
    });
}

// Resample a full-wave eigenmode's transverse E-field onto a regular grid.
// The Lee–Jin eigenvector stores e_t = γ·Et, so the spatial pattern of |e_t| is the
// transverse field pattern (the constant γ only scales the whole map), which is what a
// mode plot shows. Returns { x, y, E[ny][nx], Ex[ny][nx], Ey[ny][nx], ExIm, EyIm }: E is
// the transverse magnitude |Et| (for the heatmap); Ex/Ey are the real-part components (for
// an optional quiver/streamline overlay), ExIm/EyIm the imaginary parts. Points outside
// the mesh are left at 0. opts.parity mirrors a half-domain solve onto x < 0 as in
// resampleStatic: 'even' (V even) flips Ex, 'odd' flips Ey.
export function resampleModeField(mesh, fm, vRe, vIm, domain, opts = {}) {
    const { x, y } = makeGrid(domain, opts);
    const nx = x.length, ny = y.length;
    const { locate, coeffOf } = buildLocator(mesh);
    const nodal = recoverNodalModeField(mesh, fm, vRe, vIm);   // continuous field
    const { tris } = mesh;
    const parity = opts.parity || null;
    // A sample on a mesh line reads the triangle to its NE, the side resampleStatic reads:
    // the full-wave plot adds the two, and on an interface row the sides differ.
    const eps = 1e-9 * Math.hypot((domain.x_max ?? 0) - (domain.x_min ?? 0), (domain.y_max ?? 0) - (domain.y_min ?? 0)) || 1e-15;
    const E = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ex = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ey = Array.from({ length: ny }, () => new Float64Array(nx));
    const ExIm = Array.from({ length: ny }, () => new Float64Array(nx));
    const EyIm = Array.from({ length: ny }, () => new Float64Array(nx));
    for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
            let qx = x[i], sx = 1, sy = 1;
            if (parity && qx < 0) {
                qx = -qx;
                if (parity === 'odd') sy = -1; else sx = -1;
            }
            let t = locate(qx + eps, y[j] + eps);
            if (t < 0) t = locate(qx, y[j]);
            if (t < 0) continue;
            const coeff = coeffOf(t);
            const l0 = coeff[0][0] + coeff[0][1] * qx + coeff[0][2] * y[j];
            const l1 = coeff[1][0] + coeff[1][1] * qx + coeff[1][2] * y[j];
            const l2 = coeff[2][0] + coeff[2][1] * qx + coeff[2][2] * y[j];
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
            Ex[j][i] = sx * exr; Ey[j][i] = sy * eyr;
            ExIm[j][i] = sx * exi; EyIm[j][i] = sy * eyi;
        }
    }
    return { x, y, E, Ex, Ey, ExIm, EyIm };
}

// |J| of an MQS eddy-current solve (mqsConductorLoss opts.fieldOut, per 1 A in each
// trace), one block per conductor of `rects` (solve coordinates) carrying current. A
// rectangle is sampled on a graded grid whose lines follow the skin-refined mesh nodes
// inside it, n per axis at most, the faces sampled just inside the metal:
// { x, y, J[ny][nx] }, null outside the metal. A shaped conductor (polygon,
// ring) is its own triangles with |J| at each centroid and each vertex, { tris: [x0, y0,
// x1, y1, x2, y2, ...], J, Jv: [J0, J1, J2, ...] }: a grid cannot follow a slanted or
// curved skin layer. symX mirrors a
// half-domain solve (|J| is even about the plane in either mode).
export function sampleMqsCurrent(F, rects, n = 160, symX = null) {
    const { nodes, tris, nTris } = F.mesh;
    const ev = mqsFieldEval(F);
    // |J| of conductor k at (x, y), null outside its metal. The bounding box of a ring
    // holds the conductors inside it, which have their own blocks.
    const evalJ = (k, x, y) => {
        const t = ev.locate(x, y);
        if (t < 0 || F.triRect[t] !== k) return null;
        const J = ev.J(t, x, y);
        return J ? Math.hypot(J[0], J[1]) : null;
    };
    const out = [];
    rects.forEach((r, k) => {
        if (r.shape) {
            const xy = [], Jt = [], Jv = [];
            const mag = (t, x, y) => { const J = ev.J(t, x, y); return J ? Math.hypot(J[0], J[1]) : 0; };
            for (let t = 0; t < nTris; t++) {
                if (F.triRect[t] !== k) continue;
                let cx = 0, cy = 0;
                for (let a = 0; a < 3; a++) {
                    const v = tris[3 * t + a], x = nodes[2 * v], y = nodes[2 * v + 1];
                    xy.push(x, y);
                    Jv.push(mag(t, x, y));
                    cx += x / 3; cy += y / 3;
                }
                Jt.push(mag(t, cx, cy));
            }
            if (!Jt.some(v => v > 0)) return;
            const block = { tris: Float64Array.from(xy), J: Float64Array.from(Jt), Jv: Float64Array.from(Jv) };
            out.push(block);
            if (symX !== null) out.push({ ...block, tris: Float64Array.from(xy, (v, i) => i % 2 ? v : 2 * symX - v) });
            return;
        }
        // The solved part of the box: a half-domain solve holds x >= symX only.
        const xmin = symX !== null ? Math.max(r.xmin, symX) : r.xmin;
        if (!(r.xmax > xmin)) return;
        // Grid lines follow the nodes of the conductor's own triangles.
        const xs = [], ys = [];
        for (let t = 0; t < nTris; t++) {
            if (F.triRect[t] !== k) continue;
            for (let a = 0; a < 3; a++) { const v = tris[3 * t + a]; xs.push(nodes[2 * v]); ys.push(nodes[2 * v + 1]); }
        }
        if (!xs.length) return;
        const gx = quantileAxis(Float64Array.from(xs), xmin, r.xmax, n, []);
        const gy = quantileAxis(Float64Array.from(ys), r.ymin, r.ymax, n, []);
        // Samples on a face would land in the dielectric beside it.
        const ex = 1e-6 * (r.xmax - xmin), ey = 1e-6 * (r.ymax - r.ymin);
        const J = Array.from(gy, y => Array.from(gx, x =>
            evalJ(k, Math.min(Math.max(x, xmin + ex), r.xmax - ex), Math.min(Math.max(y, r.ymin + ey), r.ymax - ey))));
        // No current in it (an ideal ground): the plot keeps its metal fill.
        if (!J.some(row => row.some(v => v !== null))) return;
        out.push({ x: gx, y: gy, J });
        if (symX !== null) {
            out.push({ x: Float64Array.from(gx, x => 2 * symX - x).reverse(), y: gy, J: J.map(row => row.slice().reverse()) });
        }
    });
    return out;
}

// Point evaluation of an MQS solve (mqsConductorLoss opts.fieldOut): the complex J in a
// conductor triangle, null elsewhere, and the complex gradient of A in any triangle.
const fieldEvalCache = new WeakMap();
function mqsFieldEval(F) {
    // One locator per solve: the |J| and |K| plots both evaluate it.
    let ev = fieldEvalCache.get(F);
    if (!ev) fieldEvalCache.set(F, ev = mqsFieldEvaluator(F));
    return ev;
}

function mqsFieldEvaluator(F) {
    const { mesh, sol, nF, dofOf, isCondTri, triGroup, triRect, sRel, sigma, omega, Cr, Ci, CgR, CgI } = F;
    const { nodes, tris, triEdges, nNodes } = mesh;
    const { locate, coeffOf } = buildLocator(mesh);
    const dofs = t => [dofOf[tris[3 * t]], dofOf[tris[3 * t + 1]], dofOf[tris[3 * t + 2]],
        dofOf[nNodes + triEdges[3 * t]], dofOf[nNodes + triEdges[3 * t + 1]], dofOf[nNodes + triEdges[3 * t + 2]]];
    const bary = (t, x, y) => { const c = coeffOf(t); return [c, [0, 1, 2].map(k => c[k][0] + c[k][1] * x + c[k][2] * y)]; };
    const J = (t, x, y) => {
        if (t < 0) return null;
        const cls = isCondTri[t];
        if (cls !== 1 && cls !== 2) return null;
        const [, l] = bary(t, x, y);
        const N = [2 * l[0] * l[0] - l[0], 2 * l[1] * l[1] - l[1], 2 * l[2] * l[2] - l[2],
            4 * l[0] * l[1], 4 * l[1] * l[2], 4 * l[2] * l[0]];
        const g = dofs(t);
        let aR = 0, aI = 0;
        for (let k = 0; k < 6; k++) if (g[k] >= 0) { aR += N[k] * sol[g[k]]; aI += N[k] * sol[nF + g[k]]; }
        let dR = 0, dI = 0;
        if (cls === 1) { if (CgR) { dR = CgR[triGroup[t]]; dI = CgI[triGroup[t]]; } else dR = 1; }
        const uR = dR + omega * aI, uI = dI - omega * aR;
        const s = sigma * (sRel ? sRel[triRect[t]] : 1);
        return [s * (Cr * uR - Ci * uI), s * (Cr * uI + Ci * uR)];
    };
    // [dA/dx, dA/dy] as [re, im] pairs.
    const gradA = (t, x, y) => {
        const [c, l] = bary(t, x, y);
        const g = dofs(t);
        const out = [0, 0, 0, 0];
        for (let a = 0; a < 3; a++) {
            const b = (a + 1) % 3;
            const wv = 4 * l[a] - 1;
            const bases = [[g[a], wv * c[a][1], wv * c[a][2]],
                [g[3 + a], 4 * (l[a] * c[b][1] + l[b] * c[a][1]), 4 * (l[a] * c[b][2] + l[b] * c[a][2])]];
            for (const [dof, dx, dy] of bases) {
                if (dof < 0) continue;
                out[0] += dx * sol[dof]; out[1] += dx * sol[nF + dof];
                out[2] += dy * sol[dof]; out[3] += dy * sol[nF + dof];
            }
        }
        return out;
    };
    return { locate, J, gradA, isCondTri };
}

// Surface current density per 1 A in each trace from an MQS solve (mqsConductorLoss
// opts.fieldOut): the tangential magnetic field |n x H| = |dA/dn| / mu0 on every metal
// surface, conductor faces and the perfect conductors of the solve (metal walls, ideal
// grounds, A = 0) alike. It is the surface current K in the skin regime and its
// contour integral around a conductor is the conductor's current at any frequency; once
// the current fills the metal it is the surface field, not a current sheet. symX
// mirrors a half-domain solve. Returns segments { x0, y0, x1, y1, K }.
export function mqsSurfaceCurrent(F, maxLen, symX = null) {
    const MU0 = 4 * Math.PI * 1e-7;
    const { mesh, Cr, Ci } = F;
    const { nodes, edges, triEdges, nTris, nEdges } = mesh;
    const ev = mqsFieldEval(F);
    const cls = ev.isCondTri;
    const edgeTris = new Int32Array(2 * nEdges).fill(-1);
    for (let t = 0; t < nTris; t++) for (let k = 0; k < 3; k++) {
        const e = triEdges[3 * t + k];
        if (edgeTris[2 * e] < 0) edgeTris[2 * e] = t; else edgeTris[2 * e + 1] = t;
    }
    const TOL = 1e-12;
    const Cmag = Math.hypot(Cr, Ci);
    // |dA/dn| / mu0 times the drive, (nx, ny) the unit normal of the face. Evaluated in
    // the dielectric triangle: the metal side of a perfect conductor has no field.
    const surfaceH = (t, x, y, nx, ny) => {
        const g = ev.gradA(t, x, y);
        return Cmag * Math.hypot(g[0] * nx + g[2] * ny, g[1] * nx + g[3] * ny) / MU0;
    };
    // A boundary edge with A = 0 on all its unknowns is metal (a wall, a coax shield),
    // except on the symmetry plane of an odd mode.
    const { dofOf } = F, nNodes = mesh.nNodes;
    const metalBoundary = (e, n0, n1) => dofOf[n0] < 0 && dofOf[n1] < 0 && dofOf[nNodes + e] < 0
        && !(symX !== null && Math.abs(nodes[2 * n0] - symX) < 1e-12 && Math.abs(nodes[2 * n1] - symX) < 1e-12);
    const seg = segmentBuffer();
    for (let e = 0; e < nEdges; e++) {
        const ta = edgeTris[2 * e], tb = edgeTris[2 * e + 1];
        const n0 = edges[2 * e], n1 = edges[2 * e + 1];
        const x0 = nodes[2 * n0], y0 = nodes[2 * n0 + 1], dx = nodes[2 * n1] - x0, dy = nodes[2 * n1 + 1] - y0;
        // The dielectric triangle beside a metal face (driven, return or ideal metal).
        let td;
        if (tb < 0) {
            if (cls[ta] === 0 && metalBoundary(e, n0, n1)) td = ta;
            else continue;
        } else if (cls[ta] && !cls[tb]) td = tb;
        else if (cls[tb] && !cls[ta]) td = ta;
        else continue;
        const L = Math.hypot(dx, dy);
        const n = Math.max(1, Math.ceil(L / maxLen));
        for (let k = 0; k < n; k++) {
            const px = x0 + dx * (k + 0.5) / n, py = y0 + dy * (k + 0.5) / n;
            const K = surfaceH(td, px, py, dy / L, -dx / L);
            const xa = x0 + dx * k / n, ya = y0 + dy * k / n, xb = x0 + dx * (k + 1) / n, yb = y0 + dy * (k + 1) / n;
            seg.push(xa, ya, xb, yb, K);
            if (symX !== null && px > symX + TOL) seg.push(2 * symX - xa, ya, 2 * symX - xb, yb, K);
        }
    }
    return seg.out();
}

// Surface current density of the perfect-conductor limit on the conductor surfaces and
// metal walls (the loss edges), per ampere of signal current. The vacuum static solve
// phiAir has the surface charge distribution of the TEM mode in the lossless limit, and
// K = |H_t| is proportional to the normal field on the surface, so normalizing by the
// charge on the driven conductor gives K / I. Each edge is sampled every maxLen at most.
// symX mirrors a half-domain solve onto x < symX; a driven conductor cut by the plane
// carries twice the half-domain current. nets: the number of separate nets at the
// driving potential (the current is per net). Returns segments { x0, y0, x1, y1, K }.
export function surfaceCurrentPoints(mesh, fm, phi, lossMask, maxLen, symX = null, nets = 1) {
    const { nodes, tris, edges, triEdges, nTris, nEdges } = mesh;
    const pv = phi.phiVertex;
    const edgeTri = new Int32Array(nEdges).fill(-1);
    for (let t = 0; t < nTris; t++) {
        if (fm.faceF && fm.faceF[2 * t] < 0) continue;   // conductor interior
        for (let k = 0; k < 3; k++) edgeTri[triEdges[3 * t + k]] = t;
    }
    let pMax = -Infinity;
    for (let e = 0; e < nEdges; e++) {
        if (!lossMask[e]) continue;
        pMax = Math.max(pMax, pv[edges[2 * e]], pv[edges[2 * e + 1]]);
    }
    if (!(pMax > 0)) return null;
    const onDriven = (n) => Math.abs(pv[n] - pMax) < 1e-9 * pMax;
    const TOL = 1e-12;
    const seg = segmentBuffer();
    let I = 0, cutByPlane = false;
    for (let e = 0; e < nEdges; e++) {
        const t = lossMask[e] ? edgeTri[e] : -1;
        if (t < 0) continue;
        const n0 = edges[2 * e], n1 = edges[2 * e + 1];
        const x0 = nodes[2 * n0], y0 = nodes[2 * n0 + 1];
        const dx = nodes[2 * n1] - x0, dy = nodes[2 * n1 + 1] - y0;
        const L = Math.hypot(dx, dy);
        const driven = onDriven(n0) && onDriven(n1);
        if (driven && symX !== null && (Math.abs(x0 - symX) < TOL || Math.abs(x0 + dx - symX) < TOL)) cutByPlane = true;
        const c = triCoefficients(nodes, tris[3 * t], tris[3 * t + 1], tris[3 * t + 2]).coeff;
        const n = Math.max(1, Math.ceil(L / maxLen));
        for (let k = 0; k < n; k++) {
            const px = x0 + dx * (k + 0.5) / n, py = y0 + dy * (k + 0.5) / n;
            const E = Math.hypot(...evalPhiField(phi, mesh, c, t, px, py));
            if (driven) I += E * L / n;
            const xa = x0 + dx * k / n, ya = y0 + dy * k / n, xb = x0 + dx * (k + 1) / n, yb = y0 + dy * (k + 1) / n;
            seg.push(xa, ya, xb, yb, E);
            if (symX !== null && px > symX + TOL) seg.push(2 * symX - xa, ya, 2 * symX - xb, yb, E);
        }
    }
    if (cutByPlane) I *= 2;
    I /= nets;
    if (!(I > 0)) return null;
    return seg.out(1 / I);
}
