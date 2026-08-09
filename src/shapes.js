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
// domain diagonal. Boundary points must classify as INSIDE for both polarities so a
// node sitting exactly on a conductor surface is PEC.
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
export function halfCirclePolygon(cx, cy, r, n, phase = 0) {
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
//
// An `outside_*` shape is the COMPLEMENT of its body: the coax shield is "everything at
// radius >= b", which has zero meshed area (the meshed domain stops at the boundary)
// but still owns every node and edge on that boundary. That is what makes the outer
// boundary PEC and gives it loss edges, without meshing any shield metal or leaving
// dead air cavities in the corners of a bounding box.

export function isComplement(shape) {
    return shape.type === 'outside_circle' || shape.type === 'outside_polygon';
}

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
    } else {
        if (half) throw new Error('shapePoly: {half} is only supported for circle shapes');
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
    const poly = shapePoly(shape);
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
    const poly = shapePoly(shape, opts);
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
        const xmin = o.xmin !== undefined ? o.xmin : o.x_min;
        const xmax = o.xmax !== undefined ? o.xmax : o.x_max;
        const ymin = o.ymin !== undefined ? o.ymin : o.y_min;
        const ymax = o.ymax !== undefined ? o.ymax : o.y_max;
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
