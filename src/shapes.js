// Non-rectangular geometry primitives for the full-wave backend.
//
// The whole solver's geometry contract has historically been "axis-aligned rectangle
// exposing x_min/x_max/y_min/y_max", and roughly a dozen places do a point-in-box test
// against it. Coax is the first medium that cannot be expressed that way, so
// Dielectric/Conductor gained an optional `shape` descriptor and those tests route
// through shapeContains() here.
//
// THE INVARIANT THIS FILE EXISTS TO PRESERVE: when an object carries no `shape`,
// shapeContains() evaluates the literal legacy bbox expression, so every existing
// rectangular medium is bit-identical. Nothing in this module runs for them.
//
// Circles are materialized as CONVEX POLYGONS (regular n-gons), not treated as ideal
// circles, and that is deliberate rather than a compromise:
//
//   • The mesher (gmsh via OCC) has no disk primitive exported, but more importantly
//     refineTriMesh() is pure longest-edge bisection that never consults the geometry.
//     An ideal circle would be discarded at the first mesh extraction anyway.
//   • buildTriFreedomMap classifies a PEC edge by its MIDPOINT. A mesh edge on a
//     polygon side is a sub-segment of that side, so its midpoint lies exactly ON the
//     boundary at any refinement depth. Against an ideal circle those edges are chords
//     whose midpoints sit a sagitta r(1-cos(pi/n)) inside — harmless for a solid disk
//     (still inside the metal) but fatal for an `outside_circle` shield, where every
//     boundary edge would fail the test and the shield would leak.
//
// So the polygon IS the geometry, not an approximation of it, and every predicate here
// is exact with respect to that polygon.

// Relative tolerance for "is this point on the boundary" tests, scaled by the caller's
// domain diagonal (or by the local feature size where that is the meaningful scale).
// Boundary points must classify as inside for both polarities so a node sitting exactly
// on a conductor surface is PEC.
//
// Every stage that decides "is this node/edge/curve ON a shape boundary" must agree, or a
// node classifies as metal in one and as free space in the next and the PEC surface
// develops holes. So this is the single definition, imported by the mesher (its geometric
// tolerance), the freedom map and the refinement smoother's segment pinning, rather than
// each re-deriving its own 1e-9.
export const REL_SHAPE_TOL = 1e-9;

const TWO_PI = 2 * Math.PI;

// --- Polygon radius -------------------------------------------------------------
// Radius of the regular n-gon that encloses exactly the same area as a circle of
// radius r: n/2 * R^2 * sin(2pi/n) = pi r^2  =>  R = r * sqrt((2pi/n)/sin(2pi/n)).
//
// Area matching (rather than inscribing at R = r) is what keeps the electrostatics
// right: capacitance is set by the shape's logarithmic capacity, and for a convex
// polygon the area-equivalent radius tracks that far better than the circumradius —
// it removes the bulk of the leading O(1/n^2) error. At n = 64 the residual effective
// -radius error is ~1e-4, i.e. ~0.02% on a coax Z0, two orders below the FEM
// discretization error. The perimeter is 2*pi*r*(1 + pi^2/(6n^2)), so the conductor
// -loss line integral inherits a matching ~4e-4 bias.
export function polyRadiusForArea(r, n) {
    const th = TWO_PI / n;
    return r * Math.sqrt(th / Math.sin(th));
}

// Vertices of the area-matched regular n-gon, CCW, vertex k at angle phase + 2pi k/n.
// Pure and deterministic: every consumer (OCC geometry, material tagging, freedom map,
// loss, plotting) must materialize the SAME polygon or they silently disagree about
// where the metal is, so this must never depend on mesh state or call order.
export function circlePolygon(cx, cy, r, n, phase = 0) {
    const R = polyRadiusForArea(r, n);
    const poly = new Float64Array(2 * n);
    for (let k = 0; k < n; k++) {
        const th = phase + TWO_PI * k / n;
        poly[2 * k] = cx + R * Math.cos(th);
        poly[2 * k + 1] = cy + R * Math.sin(th);
    }
    return poly;
}

// The x >= 0 half of a circle polygon, as a closed convex polygon: the arc from
// angle -pi/2 up through 0 to +pi/2, closed by the chord back down the y axis.
//
// Requires n % 4 === 0 and phase === 0 so that vertices land EXACTLY on +-90 degrees.
// Without that the half is not a half: the chord would cut through a side, the two
// halves would not tile the full polygon, and the half-domain symmetry solve would be
// integrating a slightly different body than the full-domain one it is validated
// against. Callers building symmetric geometry must clamp n to a multiple of 4.
function halfCirclePolygon(cx, cy, r, n, phase = 0) {
    if (n % 4 !== 0 || phase !== 0) {
        throw new Error(`halfCirclePolygon needs n % 4 === 0 and phase === 0 (got n=${n}, phase=${phase})`);
    }
    const R = polyRadiusForArea(r, n);
    const q = n / 4;                       // vertex index of +90 degrees
    // Arc vertices from -90 degrees (index 3n/4) CCW through 0 to +90 degrees (index n/4).
    const nArc = 2 * q + 1;
    const poly = new Float64Array(2 * nArc);
    for (let i = 0; i < nArc; i++) {
        const k = 3 * q + i;               // wraps past n back through 0
        const th = TWO_PI * (k % n) / n;
        poly[2 * i] = cx + R * Math.cos(th);
        poly[2 * i + 1] = cy + R * Math.sin(th);
    }
    // The chord from (cx, cy+R) back to (cx, cy-R) closes the loop implicitly (the
    // last vertex connects to the first), so no explicit chord vertices are needed.
    return poly;
}

// --- Shape descriptors ----------------------------------------------------------
//
//   { type: 'circle',          cx, cy, r, n, phase, xSymmetric }
//   { type: 'outside_circle',  cx, cy, r, n, phase, xSymmetric }
//   { type: 'polygon',         poly: Float64Array }     CONVEX, CCW
//   { type: 'outside_polygon', poly: Float64Array }     CONVEX, CCW
//   { type: 'ring',            poly, hole }             convex outer and hole loops, CCW
//
// Custom geometry primitives (trapezoid, n-gon, n-gon ring) are 'polygon' and 'ring'
// shapes with a few optional fields:
//   faces     - face name per polygon edge (edge i runs from vertex i to i + 1): 'top',
//               'sides' or 'bottom' on a trapezoid, for per-face plating. Without it
//               every edge is 'all'.
//   thickness - the conductor's thin dimension (slab reactance, skin band gating)
//   round     - a round wire (n-gon), whose internal inductance follows the wire rule
//
// An `outside_*` shape is the COMPLEMENT of its body: the coax shield is "everything at
// radius >= b", which has zero meshed area (the meshed domain stops at the boundary)
// but still owns every node and edge on that boundary. That is what makes the outer
// boundary PEC and gives it loss edges, without meshing any shield metal or leaving
// dead air cavities in the corners of a bounding box.

export function isComplement(shape) {
    return shape.type === 'outside_circle' || shape.type === 'outside_polygon';
}

// Polygon and ring shapes (custom geometry primitives), as opposed to the circles of
// the coax model.
export const isPolyShape = shape => !!shape && (shape.type === 'polygon' || shape.type === 'ring');

// Bounds of a rect { x_min.. } or { xmin.. } as { xmin, xmax, ymin, ymax }.
export const rectOf = o => ({ xmin: o.xmin ?? o.x_min, xmax: o.xmax ?? o.x_max, ymin: o.ymin ?? o.y_min, ymax: o.ymax ?? o.y_max });

function isCircular(shape) {
    return shape.type === 'circle' || shape.type === 'outside_circle';
}

// Materialized polygon for a shape, memoized on the shape object.
// `half` returns the x >= 0 half (symmetry solves); it is a DIFFERENT body than the
// full polygon and is cached separately. Only the mesher and the constraint segments
// use the half — containment always tests against the FULL polygon (see shapeContains).
export function shapePoly(shape, { half = false } = {}) {
    const key = half ? '_polyHalf' : '_poly';
    let poly = shape[key];
    if (poly) return poly;
    if (isCircular(shape)) {
        poly = half
            ? halfCirclePolygon(shape.cx, shape.cy, shape.r, shape.n, shape.phase || 0)
            : circlePolygon(shape.cx, shape.cy, shape.r, shape.n, shape.phase || 0);
    } else if (half) {
        // The part at x >= 0 of any other convex polygon, empty when it lies left of
        // the plane. A ring's poly is its outer loop.
        if (shape.type === 'outside_polygon') throw new Error('shapePoly: {half} is not supported for outside_polygon');
        poly = clipPolyX(shape.poly, 0);
    } else {
        poly = shape.poly;
    }
    // Non-enumerable so a shape object still serializes/spreads cleanly.
    Object.defineProperty(shape, key, { value: poly, writable: true, configurable: true });
    return poly;
}

// Bounding box of the shape's positive BODY (for `outside_*`, the hole it surrounds).
export function shapeBBox(shape, opts) {
    const poly = shapePoly(shape, opts);
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (let i = 0; i < poly.length; i += 2) {
        if (poly[i] < xmin) xmin = poly[i];
        if (poly[i] > xmax) xmax = poly[i];
        if (poly[i + 1] < ymin) ymin = poly[i + 1];
        if (poly[i + 1] > ymax) ymax = poly[i + 1];
    }
    return { xmin, xmax, ymin, ymax };
}

// Signed distance from (x, y) to the polygon boundary: > 0 outside the body,
// < 0 strictly inside, 0 on the boundary. Convexity makes this the max over the
// per-edge half-plane distances, which is exact (no ray casting, no on-boundary
// ambiguity) — the reason shapes are required to be convex.
//
// This is on the hot path: buildTriFreedomMap calls shapeContains ~3x per node per
// mode per refinement pass. Circle-derived shapes get a radial fast path — a point
// inside the inradius or outside the circumradius is decided in O(1) — so the O(n)
// half-plane loop only runs for the thin annulus between them.
export function shapeSignedDist(shape, x, y) {
    if (isCircular(shape)) {
        let rad = shape._radii;
        if (!rad) {
            const R = polyRadiusForArea(shape.r, shape.n);          // circumradius
            rad = { R, rIn: R * Math.cos(Math.PI / shape.n) };      // inradius (apothem)
            Object.defineProperty(shape, '_radii', { value: rad, writable: true, configurable: true });
        }
        const d = Math.hypot(x - shape.cx, y - shape.cy);
        if (d <= rad.rIn) return d - rad.rIn;                       // strictly inside
        if (d >= rad.R) return d - rad.R;                           // strictly outside
        // In the annulus between them: fall through to the exact half-plane test.
    }
    // A regular n-gon (or ring of two): inside the inscribed circle or outside the
    // circumscribed one the radial distance decides, with the same sign and never more
    // than the true distance. Only the band between the circles takes the edge loop.
    const rad = shape.radial;
    if (rad) {
        const d = Math.hypot(x - rad.cx, y - rad.cy);
        const outer = d <= rad.rIn ? d - rad.rIn : d >= rad.rOut ? d - rad.rOut : null;
        if (shape.type !== 'ring') return outer ?? convexSignedDist(shape.poly, x, y);
        const hole = d <= rad.holeIn ? d - rad.holeIn : d >= rad.holeOut ? d - rad.holeOut : null;
        return Math.max(outer ?? convexSignedDist(shape.poly, x, y), -(hole ?? convexSignedDist(shape.hole, x, y)));
    }
    // A ring is the outer body minus the hole: outside the ring is outside the outer
    // loop or inside the hole.
    if (shape.type === 'ring') return Math.max(convexSignedDist(shape.poly, x, y), -convexSignedDist(shape.hole, x, y));
    return convexSignedDist(shapePoly(shape), x, y);
}

function convexSignedDist(poly, x, y) {
    const n = poly.length >> 1;
    let best = -Infinity;
    for (let i = 0; i < n; i++) {
        const j = (i + 1) % n;
        const ax = poly[2 * i], ay = poly[2 * i + 1];
        const ex = poly[2 * j] - ax, ey = poly[2 * j + 1] - ay;
        const len = Math.hypot(ex, ey);
        if (!(len > 0)) continue;
        // CCW polygon: outward normal is (ey, -ex)/len, so this is > 0 outside the edge.
        const d = ((x - ax) * ey - (y - ay) * ex) / len;
        if (d > best) best = d;
    }
    return best;
}

// Is (x, y) inside the object, within tol?
//
// No shape  -> the LEGACY bbox test, character for character. This is the fallback that
//              keeps every rectangular medium bit-identical; do not "simplify" it.
// body      -> s <= +tol  (true ON the boundary)
// complement-> s >= -tol  (true ON the boundary)
//
// Both polarities include their boundary, which is what makes a node sitting exactly on
// a conductor surface PEC. It also means the boundary belongs to both a disk and its
// complement, but those are never both present as separate conductors in one geometry.
//
// `o` is the geometry object (Conductor/Dielectric or a mesher condRect), accepting
// either the x_min/x_max naming or the xmin/xmax naming used inside the mesher.
export function shapeContains(o, x, y, tol = 0) {
    const shape = o.shape;
    if (!shape) {
        const xmin = o.xmin !== undefined ? o.xmin : o.x_min;
        const xmax = o.xmax !== undefined ? o.xmax : o.x_max;
        const ymin = o.ymin !== undefined ? o.ymin : o.y_min;
        const ymax = o.ymax !== undefined ? o.ymax : o.y_max;
        return x >= xmin - tol && x <= xmax + tol && y >= ymin - tol && y <= ymax + tol;
    }
    const s = shapeSignedDist(shape, x, y);
    return isComplement(shape) ? (s >= -tol) : (s <= tol);
}

// Cross-sectional area of the object.
//
// Feeds the DC resistance (R_dc = 1/(sigma*A)), where a bbox would be badly wrong:
// a disk's bbox area is 4a^2 vs the true pi*a^2, a +27% error on R_dc.
// A complement shape has no finite area — it is a zero-thickness PEC shell in this
// model — and returns 0 so it contributes nothing to any area sum.
export function shapeArea(o, opts) {
    const shape = o.shape;
    if (!shape) {
        const w = o.width !== undefined ? o.width
            : ((o.xmax !== undefined ? o.xmax : o.x_max) - (o.xmin !== undefined ? o.xmin : o.x_min));
        const h = o.height !== undefined ? o.height
            : ((o.ymax !== undefined ? o.ymax : o.y_max) - (o.ymin !== undefined ? o.ymin : o.y_min));
        return Math.abs(w * h);
    }
    if (isComplement(shape)) return 0;
    if (shape.type === 'ring') {
        const hole = (opts && opts.half) ? clipPolyX(shape.hole, 0) : shape.hole;
        return polyArea(shapePoly(shape, opts)) - polyArea(hole);
    }
    return polyArea(shapePoly(shape, opts));
}

function polyArea(poly) {
    const n = poly.length >> 1;
    let a2 = 0;
    for (let i = 0; i < n; i++) {
        const j = (i + 1) % n;
        a2 += poly[2 * i] * poly[2 * j + 1] - poly[2 * j] * poly[2 * i + 1];
    }
    return Math.abs(a2) / 2;
}

// Boundary segments [{x0,y0,x1,y1}, ...] of a shaped object; [] for a bare rect
// (rects use the existing constraintXRanges/constraintYRanges mechanism instead).
//
// These become mesh constraint segments: a node lying on one is PINNED during the
// refinement smoother. Without that, Laplacian smoothing classifies polygon-boundary
// nodes as free (they are not on the domain bbox walls) and pulls them inward every
// pass, shrinking the shield until the complement test stops matching and the PEC
// develops holes. Bisection midpoints of a straight segment land exactly on it, so
// pinning costs nothing geometrically.
//
// For a complement shape the segments are its inner boundary — the same polygon.
//
// In the {half} case the CLOSING edge (last vertex back to first) is the chord down
// x = cx, i.e. the symmetry plane. It is deliberately omitted: nodes there must stay
// free to slide ALONG the plane (the mesher constrains them via constraintXRanges),
// and pinning a whole plane of nodes would strand slivers the smoother could not fix.
export function shapeSegments(o, opts) {
    const shape = o.shape;
    if (!shape) return [];
    if (!isCircular(shape)) {
        // Polygon and ring loops. On a half domain the edges lying on the plane are
        // the cut, not a surface, and are left out like the circle's chord.
        const half = !!(opts && opts.half);
        const segs = [];
        const add = (poly, reverse) => {
            const n = poly.length >> 1;
            for (let i = 0; i < n; i++) {
                const j = (i + 1) % n;
                const [a, b] = reverse ? [j, i] : [i, j];
                const seg = { x0: poly[2 * a], y0: poly[2 * a + 1], x1: poly[2 * b], y1: poly[2 * b + 1] };
                if (half && seg.x0 === 0 && seg.x1 === 0) continue;
                segs.push(seg);
            }
        };
        add(shapePoly(shape, opts), false);
        // The hole runs clockwise so the outward normal of every segment points away
        // from the metal.
        if (shape.type === 'ring') add(half ? clipPolyX(shape.hole, 0) : shape.hole, true);
        return segs;
    }
    const poly = shapePoly(shape, opts);
    const n = poly.length >> 1;
    const last = (opts && opts.half) ? n - 1 : n;   // skip the closing chord on a half
    const segs = [];
    for (let i = 0; i < last; i++) {
        const j = (i + 1) % n;
        segs.push({ x0: poly[2 * i], y0: poly[2 * i + 1], x1: poly[2 * j], y1: poly[2 * j + 1] });
    }
    return segs;
}

// Perimeter of a shaped object's boundary (0 for a bare rect — callers use 2(w+h)).
export function shapePerimeter(o) {
    if (!o.shape) return 0;
    let p = 0;
    for (const s of shapeSegments(o)) p += Math.hypot(s.x1 - s.x0, s.y1 - s.y0);
    return p;
}

// Walk a shaped object's boundary by ARC LENGTH: t in [0,1) → a point on the boundary
// plus the unit OUTWARD normal there. Streamline seeding uses it to place seeds evenly
// around a curved conductor and to weight them by the normal field |E·n|.
// For a complement shape the normal points into the hole (toward the field region),
// which is still "away from the metal".
export function shapePerimeterPoint(o, t) {
    const segs = shapeSegments(o);
    const total = shapePerimeter(o);
    let target = ((t % 1) + 1) % 1 * total;
    for (const s of segs) {
        const ex = s.x1 - s.x0, ey = s.y1 - s.y0;
        const len = Math.hypot(ex, ey);
        if (target > len && len > 0) { target -= len; continue; }
        const u = len > 0 ? target / len : 0;
        // CCW polygon ⇒ outward normal of edge (ex, ey) is (ey, -ex)/len.
        const sgn = isComplement(o.shape) ? -1 : 1;
        return { x: s.x0 + u * ex, y: s.y0 + u * ey,
                 nx: sgn * ey / len, ny: -sgn * ex / len };
    }
    const s = segs[segs.length - 1];
    const ex = s.x1 - s.x0, ey = s.y1 - s.y0, len = Math.hypot(ex, ey);
    const sgn = isComplement(o.shape) ? -1 : 1;
    return { x: s.x1, y: s.y1, nx: sgn * ey / len, ny: -sgn * ex / len };
}

// Distance from (x, y) to the object's boundary, unsigned. Streamline stepping uses
// it to size steps near metal.
export function distToShapeBoundary(o, x, y) {
    if (!o.shape) {
        const { xmin, xmax, ymin, ymax } = rectOf(o);
        const dx = Math.max(xmin - x, 0, x - xmax);
        const dy = Math.max(ymin - y, 0, y - ymax);
        if (dx > 0 || dy > 0) return Math.hypot(dx, dy);
        return Math.min(x - xmin, xmax - x, y - ymin, ymax - y);
    }
    return Math.abs(shapeSignedDist(o.shape, x, y));
}

// Where does the segment (x1,y1)->(x2,y2) first enter the object? Returns {x, y} or
// null. Used by streamline integration to stop a line exactly at a conductor surface.
// Bisection on the signed distance: the segment is short (one integration step) and
// the shape is convex, so the sign change is unique and 40 iterations is exact to
// double precision.
export function segShapeBoundaryHit(x1, y1, x2, y2, o) {
    if (!o.shape) return null;
    const inside = (t) => shapeContains(o, x1 + t * (x2 - x1), y1 + t * (y2 - y1), 0);
    if (inside(0) || !inside(1)) return null;      // started inside, or never enters
    let lo = 0, hi = 1;
    for (let i = 0; i < 40; i++) {
        const mid = (lo + hi) / 2;
        if (inside(mid)) hi = mid; else lo = mid;
    }
    return { x: x1 + hi * (x2 - x1), y: y1 + hi * (y2 - y1), conductor: o };
}

// SVG path for an annular ring, as two closed polygonal loops for evenodd filling.
// Plotly's layout.shapes[].path grammar allows only M/L/H/V/Q/C/T/S/Z — there is no
// elliptical-arc command — so the ring must be drawn as polygons rather than arcs.
export function svgRingPath(cx, cy, rIn, rOut, n = 180) {
    const loop = (r, dir) => {
        let d = '';
        for (let k = 0; k < n; k++) {
            const th = dir * TWO_PI * k / n;
            const x = cx + r * Math.cos(th), y = cy + r * Math.sin(th);
            d += (k === 0 ? 'M' : 'L') + x.toFixed(6) + ',' + y.toFixed(6);
        }
        return d + 'Z';
    };
    // Opposite winding for the hole is not required by evenodd, but keeps the path
    // valid under nonzero filling too.
    return loop(rOut, 1) + loop(rIn, -1);
}

// --- Polygon primitives (custom geometry) -----------------------------------------

// Part of the convex CCW polygon `poly` at x >= x0 (one Sutherland-Hodgman stage),
// still CCW, empty when nothing is right of x0. Vertices on the plane stay exact, and
// a crossing edge gets its new vertex at exactly x = x0.
function clipPolyX(poly, x0) {
    const n = poly.length >> 1;
    // A vertex meant to lie on the plane (a rotated n-gon's, an inset core's) can come
    // out of the arithmetic an ulp off it: it is put on the plane, or the cut would add
    // a second vertex beside it and leave a zero-length side.
    let ext = 0;
    for (let i = 0; i < poly.length; i++) ext = Math.max(ext, Math.abs(poly[i] - (i % 2 === 0 ? x0 : 0)));
    const tol = ext * 1e-12;
    const px = i => (Math.abs(poly[2 * i] - x0) <= tol ? x0 : poly[2 * i]);
    const out = [];
    for (let i = 0; i < n; i++) {
        const j = (i + 1) % n;
        const ax = px(i), ay = poly[2 * i + 1], bx = px(j), by = poly[2 * j + 1];
        const aIn = ax >= x0, bIn = bx >= x0;
        if (aIn) out.push(ax, ay);
        if (aIn !== bIn && ax !== x0 && bx !== x0) {
            out.push(x0, ay + (by - ay) * (x0 - ax) / (bx - ax));
        }
    }
    // Coincident neighbours (a cut through a vertex) leave a zero-length side.
    const kept = [];
    const m = out.length >> 1;
    for (let i = 0; i < m; i++) {
        const j = (i + 1) % m;
        if (Math.hypot(out[2 * j] - out[2 * i], out[2 * j + 1] - out[2 * i + 1]) > tol) kept.push(out[2 * i], out[2 * i + 1]);
    }
    // A polygon touching the plane from the left leaves a point or a segment.
    let right = false;
    for (let i = 0; i < kept.length; i += 2) if (kept[i] > x0) { right = true; break; }
    return right && kept.length >= 6 ? new Float64Array(kept) : new Float64Array(0);
}

// Closed loops of a shape for the mesher: [outer, hole] for a ring, [polygon] for the
// others. On a half domain (x >= 0) a ring cut by the plane is one C-shaped loop:
// the outer arc up, down the plane to the hole, the hole arc back down, and the plane
// again to the start.
export function shapeLoops(shape, { half = false } = {}) {
    if (shape.type !== 'ring') {
        const poly = shapePoly(shape, { half });
        return poly.length ? [poly] : [];
    }
    if (!half) return [shape.poly, shape.hole];
    const outer = clipPolyX(shape.poly, 0);
    if (!outer.length) return [];
    if (outer.length === shape.poly.length && !outer.some((v, i) => v !== shape.poly[i])) return [shape.poly, shape.hole];
    const hole = clipPolyX(shape.hole, 0);
    if (!hole.length) return [outer];
    // Arc of a clipped CCW loop from its lower plane vertex to its upper one.
    const arc = (poly) => {
        const n = poly.length >> 1;
        for (let i = 0; i < n; i++) {
            const j = (i + 1) % n;
            // The cut runs down the plane, from the upper plane vertex i to the lower j.
            if (poly[2 * i] === 0 && poly[2 * j] === 0 && poly[2 * i + 1] > poly[2 * j + 1]) {
                const pts = [];
                for (let k = 0; k < n; k++) { const m = (j + k) % n; pts.push(poly[2 * m], poly[2 * m + 1]); }
                return pts;
            }
        }
        return null;
    };
    const a = arc(outer), b = arc(hole);
    if (!a || !b) throw new Error('shapeLoops: a ring cut by the symmetry plane must be centred on it');
    const loop = [...a];
    for (let k = (b.length >> 1) - 1; k >= 0; k--) loop.push(b[2 * k], b[2 * k + 1]);
    return [new Float64Array(loop)];
}

// Distance from (x, y) to the segment (ax, ay)-(bx, by).
function pointSegDist(x, y, ax, ay, bx, by) {
    const ex = bx - ax, ey = by - ay;
    const l2 = ex * ex + ey * ey;
    let t = l2 > 0 ? ((x - ax) * ex + (y - ay) * ey) / l2 : 0;
    t = Math.max(0, Math.min(1, t));
    return Math.hypot(x - ax - t * ex, y - ay - t * ey);
}

// Face name and nearest boundary point of a shape at (x, y): the edge closest to the
// point decides. Polygons name their edges through shape.faces, anything else is 'all'.
export function shapeFaceAt(shape, x, y) {
    const loops = shape.type === 'ring' ? [shape.poly, shape.hole] : [shapePoly(shape)];
    let best = Infinity, face = 'all', px = x, py = y;
    loops.forEach((poly, li) => {
        const n = poly.length >> 1;
        for (let i = 0; i < n; i++) {
            const j = (i + 1) % n;
            const ax = poly[2 * i], ay = poly[2 * i + 1], bx = poly[2 * j], by = poly[2 * j + 1];
            const d = pointSegDist(x, y, ax, ay, bx, by);
            if (d < best) {
                best = d;
                face = (li === 0 && shape.faces) ? shape.faces[i] : 'all';
                const ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey;
                const t = l2 > 0 ? Math.max(0, Math.min(1, ((x - ax) * ex + (y - ay) * ey) / l2)) : 0;
                px = ax + t * ex; py = ay + t * ey;
            }
        }
    });
    return { face, x: px, y: py, dist: best };
}

// Boundary loops of a conductor or dielectric: its rectangle, or the loops of its
// shape. null for a complement shape, which has no finite body.
function bodyLoops(o) {
    if (o.shape) {
        if (isComplement(o.shape)) return null;
        return o.shape.type === 'ring' ? [o.shape.poly, o.shape.hole] : [shapePoly(o.shape)];
    }
    const { xmin: x0, xmax: x1, ymin: y0, ymax: y1 } = rectOf(o);
    return [new Float64Array([x0, y0, x1, y0, x1, y1, x0, y1])];
}

function segmentsCross(ax, ay, bx, by, cx, cy, dx, dy) {
    const o = (px, py, qx, qy, rx, ry) => Math.sign((qx - px) * (ry - py) - (qy - py) * (rx - px));
    const o1 = o(ax, ay, bx, by, cx, cy), o2 = o(ax, ay, bx, by, dx, dy);
    const o3 = o(cx, cy, dx, dy, ax, ay), o4 = o(cx, cy, dx, dy, bx, by);
    return o1 * o2 < 0 && o3 * o4 < 0;
}

// Shortest distance between two bodies (rectangles or shaped objects), 0 when they
// overlap or one contains the other. Between two separate bodies the closest pair of
// points always includes a vertex, so vertex-to-edge distances are exact.
export function bodyDistance(a, b) {
    const la = bodyLoops(a), lb = bodyLoops(b);
    if (!la || !lb) return Infinity;
    const inside = (o, x, y) => shapeContains(o, x, y, 0);
    let best = Infinity;
    const scan = (loops, other, otherLoops) => {
        for (const p of loops) {
            const n = p.length >> 1;
            for (let i = 0; i < n; i++) {
                const x = p[2 * i], y = p[2 * i + 1];
                if (inside(other, x, y)) { best = 0; return; }
                for (const q of otherLoops) {
                    const m = q.length >> 1;
                    for (let k = 0; k < m; k++) {
                        const l = (k + 1) % m;
                        best = Math.min(best, pointSegDist(x, y, q[2 * k], q[2 * k + 1], q[2 * l], q[2 * l + 1]));
                    }
                }
            }
        }
    };
    scan(la, b, lb);
    if (best > 0) scan(lb, a, la);
    if (best === 0) return 0;
    // Crossing edges without a vertex inside (two bars forming a cross).
    for (const p of la) for (const q of lb) {
        const n = p.length >> 1, m = q.length >> 1;
        for (let i = 0; i < n; i++) {
            const j = (i + 1) % n;
            for (let k = 0; k < m; k++) {
                const l = (k + 1) % m;
                if (segmentsCross(p[2 * i], p[2 * i + 1], p[2 * j], p[2 * j + 1],
                                  q[2 * k], q[2 * k + 1], q[2 * l], q[2 * l + 1])) return 0;
            }
        }
    }
    return best;
}

// Mirror image of a polygon about x = 0, still CCW. Edge j of the image is the image
// of edge n-2-j of the original, which is where its face name goes.
function mirrorPoly(poly) {
    const n = poly.length >> 1;
    const out = new Float64Array(2 * n);
    for (let j = 0; j < n; j++) {
        const k = n - 1 - j;
        out[2 * j] = -poly[2 * k];
        out[2 * j + 1] = poly[2 * k + 1];
    }
    return out;
}

// A polygon or ring shape mirrored about x = 0.
export function mirrorShapeX(shape) {
    const out = { ...shape, poly: mirrorPoly(shape.poly) };
    if (shape.hole) out.hole = mirrorPoly(shape.hole);
    if (shape.radial) out.radial = { ...shape.radial, cx: -shape.radial.cx };
    if (shape.faces) {
        const n = shape.faces.length;
        out.faces = shape.faces.map((_, j) => shape.faces[((n - 2 - j) % n + n) % n]);
    }
    return out;
}

// A polygon or ring shape moved by dx along x.
export function translateShapeX(shape, dx) {
    const move = poly => poly.map((v, i) => (i % 2 === 0 ? v + dx : v));
    const out = { ...shape, poly: move(shape.poly) };
    if (shape.hole) out.hole = move(shape.hole);
    if (shape.radial) out.radial = { ...shape.radial, cx: shape.radial.cx + dx };
    return out;
}

// Same vertex set within tol, whatever the starting vertex.
function samePoly(a, b, tol) {
    if (a.length !== b.length) return false;
    const n = a.length >> 1;
    for (let i = 0; i < n; i++) {
        let found = false;
        for (let k = 0; k < n && !found; k++) {
            found = Math.abs(a[2 * i] - b[2 * k]) <= tol && Math.abs(a[2 * i + 1] - b[2 * k + 1]) <= tol;
        }
        if (!found) return false;
    }
    return true;
}

// Is shape b the mirror image of shape a about x = 0?
export function isMirrorShape(a, b, tol) {
    if (!a || !b || a.type !== b.type || (a.type !== 'polygon' && a.type !== 'ring')) return false;
    if (!samePoly(mirrorPoly(a.poly), b.poly, tol)) return false;
    return a.type !== 'ring' || samePoly(mirrorPoly(a.hole), b.hole, tol);
}

// Closed path of a polygon or ring shape for a Plotly path shape, in mm.
export function svgShapePath(shape) {
    const loops = shape.type === 'ring' ? [shape.poly, shape.hole] : [shapePoly(shape)];
    return loops.map(poly => {
        let d = '';
        for (let i = 0; i < poly.length; i += 2) {
            d += (i === 0 ? 'M' : 'L') + (poly[i] * 1000).toFixed(9) + ',' + (poly[i + 1] * 1000).toFixed(9);
        }
        return d + 'Z';
    }).join(' ');
}

// Intersection of the half-planes of a convex CCW polygon's edges, each edge moved
// inward by offsets[i] (outward for a negative offset): the polygon with some faces
// inset (a plating layer's inside) or grown (a hole widened by a plating layer). Empty
// when nothing is left.
function offsetConvex(poly, offsets) {
    const n = poly.length >> 1;
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (let i = 0; i < n; i++) {
        xmin = Math.min(xmin, poly[2 * i]); xmax = Math.max(xmax, poly[2 * i]);
        ymin = Math.min(ymin, poly[2 * i + 1]); ymax = Math.max(ymax, poly[2 * i + 1]);
    }
    const grow = Math.max(0, ...offsets.map(d => -d)) * 2 + (xmax - xmin + ymax - ymin) * 1e-6;
    let q = [xmin - grow, ymin - grow, xmax + grow, ymin - grow, xmax + grow, ymax + grow, xmin - grow, ymax + grow];
    for (let i = 0; i < n && q.length >= 6; i++) {
        const j = (i + 1) % n;
        const ax = poly[2 * i], ay = poly[2 * i + 1];
        const ex = poly[2 * j] - ax, ey = poly[2 * j + 1] - ay;
        const l = Math.hypot(ex, ey);
        if (!(l > 0)) continue;
        const nx = -ey / l, ny = ex / l, d = offsets[i];   // inward normal of a CCW edge
        const side = (x, y) => (x - ax) * nx + (y - ay) * ny - d;
        const out = [];
        const m = q.length >> 1;
        for (let k = 0; k < m; k++) {
            const k2 = (k + 1) % m;
            const px = q[2 * k], py = q[2 * k + 1], rx = q[2 * k2], ry = q[2 * k2 + 1];
            const sp = side(px, py), sr = side(rx, ry);
            if (sp >= 0) out.push(px, py);
            if ((sp >= 0) !== (sr >= 0)) {
                const t = sp / (sp - sr);
                out.push(px + t * (rx - px), py + t * (ry - py));
            }
        }
        q = out;
    }
    return q.length >= 6 && polyArea(q) > 0 ? new Float64Array(q) : new Float64Array(0);
}

// Two conductors or dielectrics { x_min.., shape? } share some area. Touching blocks
// whose shared edge differs by rounding do not overlap. Shapes are tested on a grid of
// points over the common box.
export function bodiesOverlap(a, b) {
    if ((a.shape && isComplement(a.shape)) || (b.shape && isComplement(b.shape))) return false;
    const w = Math.min(a.x_max, b.x_max) - Math.max(a.x_min, b.x_min);
    const h = Math.min(a.y_max, b.y_max) - Math.max(a.y_min, b.y_min);
    const tol = 1e-9 * Math.max(a.x_max - a.x_min, a.y_max - a.y_min, b.x_max - b.x_min, b.y_max - b.y_min);
    if (!(w > tol && h > tol)) return false;
    if (!a.shape && !b.shape) return true;
    const n = 16;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
        const x = Math.max(a.x_min, b.x_min) + w * (i + 0.5) / n, y = Math.max(a.y_min, b.y_min) + h * (j + 0.5) / n;
        if (shapeContains(a, x, y, 0) && shapeContains(b, x, y, 0)) return true;
    }
    return false;
}

// Cross-section each conductor contributes where conductors of one kind (positive,
// negative, ground) overlap: an overlapped area counts once, for the later conductor,
// whose metal fills it. A Map from conductor index, only for the overlapping ones; null
// when nothing overlaps. Axis-aligned rects are cut into the cells of their edge lines,
// which is exact; shapes add a fine grid over their bounding box.
export function visibleAreas(conductors) {
    const kind = c => (c.is_signal ? (c.polarity < 0 ? -1 : 1) : 0);
    const box = c => ({ xmin: c.x_min, xmax: c.x_max, ymin: c.y_min, ymax: c.y_max });
    const cs = conductors;
    const inCluster = new Set();
    for (let i = 0; i < cs.length; i++) for (let j = i + 1; j < cs.length; j++) {
        if (kind(cs[i]) === kind(cs[j]) && bodiesOverlap(cs[i], cs[j])) { inCluster.add(i); inCluster.add(j); }
    }
    if (!inCluster.size) return null;
    const out = new Map();
    for (const k of [-1, 0, 1]) {
        // In list order: the later conductor wins.
        const idx = [...inCluster].filter(i => kind(cs[i]) === k).sort((a, b) => a - b);
        if (!idx.length) continue;
        const xs = new Set(), ys = new Set();
        let bx0 = Infinity, bx1 = -Infinity, by0 = Infinity, by1 = -Infinity;
        for (const i of idx) {
            const b = box(cs[i]);
            xs.add(b.xmin); xs.add(b.xmax); ys.add(b.ymin); ys.add(b.ymax);
            bx0 = Math.min(bx0, b.xmin); bx1 = Math.max(bx1, b.xmax); by0 = Math.min(by0, b.ymin); by1 = Math.max(by1, b.ymax);
        }
        if (idx.some(i => cs[i].shape)) {
            const N = 400;
            for (let m = 1; m < N; m++) { xs.add(bx0 + (bx1 - bx0) * m / N); ys.add(by0 + (by1 - by0) * m / N); }
        }
        const X = [...xs].sort((a, b) => a - b), Y = [...ys].sort((a, b) => a - b);
        for (const i of idx) out.set(i, 0);
        for (let a = 0; a + 1 < X.length; a++) {
            for (let b = 0; b + 1 < Y.length; b++) {
                const x = (X[a] + X[a + 1]) / 2, y = (Y[b] + Y[b + 1]) / 2;
                for (let q = idx.length - 1; q >= 0; q--) {
                    const c = cs[idx[q]];
                    if (shapeContains(c, x, y, 0)) { out.set(idx[q], out.get(idx[q]) + (X[a + 1] - X[a]) * (Y[b + 1] - Y[b])); break; }
                }
            }
        }
    }
    return out;
}

// --- Plating geometry -------------------------------------------------------------
// A plating layer lies inside its conductor's outline. These give its cross-section
// and what is left inside it, for a rect { x_min.. } / { xmin.. } or a shaped object.

const isPlated = pl => !!(pl && pl.sigma > 0 && (pl.thickness ?? 0) > 0 && (pl.top || pl.sides || pl.bottom || pl.all));

// Whether the edges of a shape carry plating `pl`: per face name, every edge for a
// round shape or when the plating is all around.
function platedEdge(shape, pl, i) {
    return !!(pl.all || !shape.faces || pl[shape.faces[i]]);
}

// Inside of the plating layer `pl` of object o: its outline with the plated faces moved
// in by the plating thickness, as { rects: [rect] } or { shape }, null when nothing is
// left (the conductor is plating metal through), undefined without plating or for a
// shape the layer cannot be built for (a complement).
export function platingCoreOf(o, pl = o.plating) {
    if (!isPlated(pl)) return undefined;
    const t = pl.thickness;
    const sh = o.shape;
    if (sh) {
        if (sh.type !== 'polygon' && sh.type !== 'ring') return undefined;
        const n = sh.poly.length >> 1;
        const poly = offsetConvex(sh.poly, Array.from({ length: n }, (_, i) => (platedEdge(sh, pl, i) ? t : 0)));
        if (!poly.length) return null;
        if (sh.type === 'polygon') return { shape: { type: 'polygon', prim: sh.prim, poly, faces: sh.faces && sh.faces.length === n ? sh.faces : undefined } };
        // A ring is plated on both surfaces: the hole grows by the layer.
        const hole = offsetConvex(sh.hole, new Array(sh.hole.length >> 1).fill(-t));
        for (let i = 0; i < hole.length; i += 2) {
            if (!shapeContains({ shape: { type: 'polygon', poly } }, hole[i], hole[i + 1], -t * 1e-3)) return null;
        }
        return { shape: { type: 'ring', prim: sh.prim, poly, hole } };
    }
    const r = rectOf(o);
    const core = { xmin: r.xmin + (pl.sides ? t : 0), xmax: r.xmax - (pl.sides ? t : 0),
                   ymin: r.ymin + (pl.bottom ? t : 0), ymax: r.ymax - (pl.top ? t : 0) };
    // A plating at least as thick as the conductor is plating through, whatever faces.
    if (t >= r.ymax - r.ymin || !(core.xmax > core.xmin && core.ymax > core.ymin)) return null;
    return { rects: [core] };
}

// Plating `pl` fills the whole cross-section of o: it is plating metal through.
export function platedThrough(o, pl = o.plating) {
    return isPlated(pl) && platingCoreOf(o, pl) === null;
}

// Cross-section of the plating layer of o (0 without plating, the whole area when it
// is plated through): the plated faces times the thickness, corners counted once.
export function platingArea(o, pl = o.plating) {
    if (!isPlated(pl)) return 0;
    const core = platingCoreOf(o, pl);
    const total = shapeArea(o);
    if (core === null) return total;
    if (core === undefined) return 0;
    const inner = core.shape ? shapeArea({ shape: core.shape }) : core.rects.reduce((a, k) => a + (k.xmax - k.xmin) * (k.ymax - k.ymin), 0);
    return Math.max(0, total - inner);
}

// Every vertex of object b (its loops, or its rect corners) inside the hole of the
// ring shape `ring`: b lies in the cavity the ring shields.
export function insideRingHole(ring, b) {
    const loops = bodyLoops(b);
    if (!loops) return false;
    const hole = { shape: { type: 'polygon', poly: ring.hole } };
    return loops.every(p => {
        for (let i = 0; i < p.length; i += 2) if (!shapeContains(hole, p[i], p[i + 1], 0)) return false;
        return true;
    });
}
