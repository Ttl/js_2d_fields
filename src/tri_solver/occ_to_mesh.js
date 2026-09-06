// OCC (OpenCASCADE) boolean-fragment mesher — alternative to the GEO-kernel band/hole
// mesher in geom_to_mesh.js. Builds the cross-section as overlapping axis-aligned
// rectangles (a background "air" rectangle + every dielectric + every conductor) and
// lets OCC `fragment` compute a fully CONFORMING arrangement: one face per material
// region, with mesh-aligned edges at EVERY material boundary — horizontal AND vertical,
// and only where materials actually change (local interfaces, not whole-domain lines).
// This fixes the band mesher's two structural weaknesses: solder-mask side edges weren't
// conforming, and a local horizontal interface (e.g. SM-trace-top) was forced to span the
// whole domain. Element sizing is a Distance→Threshold background field (fine at
// conductors, graded to coarse), which gmsh meshes robustly.
//
// Output is the standard solver mesh object (nodes/tris/edges/condRect/epsMap/...) that
// the rest of the tri backend consumes. Material/role tagging stays centroid-based
// (tagMaterials + condRect.rects), so no per-face bookkeeping is needed.
//
// SHAPED GEOMETRY (coax): a conductor/dielectric may carry a `shape` descriptor (see
// shapes.js) instead of being its own bounding box. Such objects are emitted as convex
// polygons via AddPoint/AddLine/AddCurveLoop/AddPlaneSurface rather than AddRectangle,
// and every containment test routes through shapeContains(). `opts.domainShape` further
// replaces the background rectangle with a polygon, so the meshed domain can itself be
// a disk. A complement shape ('outside_circle') emits no OCC geometry at all. Its
// boundary is the domain outline, but still gets a condRects entry so the freedom map
// makes that outline PEC and the loss integral finds its surface edges.

import { shapeContains, shapePoly, shapeBBox, shapeArea, shapeSegments, shapeSignedDist,
         isComplement, REL_SHAPE_TOL } from '../shapes.js';

// Domain-diagonal-relative geometric tolerance. Shared with the freedom map and the
// refinement smoother (see REL_SHAPE_TOL) so all three agree on where a boundary is.
const REL_TOL = REL_SHAPE_TOL;

// OCC model units. OpenCASCADE works to a fixed Precision::Confusion of 1e-7
// model units, so a cross-section built in metres falls apart for sub-micron
// features. Every coordinate and size handed to gmsh is multiplied by OCC_SCALE
// and every coordinate read back divided by it.
const OCC_SCALE = 1e6;

// Per-triangle ε / loss maps by centroid (last-match-wins over the dielectric list,
// matching _setup_geometry). Reusable: call again after refineTriMesh.
export function tagMaterials(mesh, dielectrics, tol) {
    const { nodes, tris, nTris } = mesh;
    if (tol === undefined) tol = 1e-12;
    const epsMap = new Array(nTris), lossMap = new Array(nTris);
    // shapeContains() falls back to the literal bbox test when d.shape is absent, so
    // shapeless dielectrics behave exactly as before.
    const defs = dielectrics.map(d => ({ o: d, er: d.epsilon_r, tand: d.tan_delta || 0 }));
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const xc = (nodes[2 * v0] + nodes[2 * v1] + nodes[2 * v2]) / 3;
        const yc = (nodes[2 * v0 + 1] + nodes[2 * v1 + 1] + nodes[2 * v2 + 1]) / 3;
        let er = 1.0, tand = 0.0;
        for (const d of defs) {
            if (shapeContains(d.o, xc, yc, tol)) { er = d.er; tand = d.tand; }
        }
        epsMap[t] = { re: er, im: 0 };
        lossMap[t] = { re: er * tand, im: 0 };
    }
    return { epsMap, lossMap };
}

// Validate a built triangular mesh before the FEM solve — the full-wave analogue of
// field_solver.js's validate_laplace_inputs(). Returns an array of human-readable error
// strings (empty ⇒ OK). Catches the degenerate-geometry / collapsed-mesh cases that would
// otherwise surface as opaque NaN/Inf (Area=0 → division) or a heap abort deep in the
// eigensolve, mirroring how the rectilinear backend rejects a bad Laplace input up front.
export function validateTriMesh(mesh, condRect) {
    const errors = [];
    if (!mesh || !mesh.nNodes || !mesh.nTris) {
        errors.push('Mesh is empty (no nodes or triangles).');
        return errors;
    }
    const { nodes, tris, nNodes, nTris, epsMap } = mesh;
    for (let i = 0; i < 2 * nNodes; i++) {
        if (!Number.isFinite(nodes[i])) { errors.push(`Non-finite node coordinate at index ${i}.`); break; }
    }
    let degenerate = 0, badIdx = false;
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        if (v0 < 0 || v1 < 0 || v2 < 0 || v0 >= nNodes || v1 >= nNodes || v2 >= nNodes) { badIdx = true; continue; }
        const ax = nodes[2 * v0], ay = nodes[2 * v0 + 1];
        const bx = nodes[2 * v1], by = nodes[2 * v1 + 1];
        const cx = nodes[2 * v2], cy = nodes[2 * v2 + 1];
        if (!(Math.abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) > 0)) degenerate++;
    }
    if (badIdx) errors.push('A triangle references an out-of-range node index.');
    if (degenerate) errors.push(`${degenerate} degenerate (zero-area) triangle(s) in the mesh.`);
    if (!epsMap || epsMap.length !== nTris) {
        errors.push(`Permittivity map length ${epsMap ? epsMap.length : 'null'} does not match triangle count ${nTris}.`);
    } else {
        for (let t = 0; t < nTris; t++) {
            const e = epsMap[t];
            if (!e || !Number.isFinite(e.re) || !(e.re > 0)) { errors.push(`Invalid permittivity at triangle ${t} (re=${e && e.re}).`); break; }
        }
    }
    const roles = condRect && condRect.rectRoles;
    // A fully enclosed all-PEC domain guides modes with no conductor meshed inside it at
    // all: A hollow waveguide, where the walls themselves are the only metal and they
    // live in the boundary conditions rather than in rectRoles. Every other geometry
    // without a signal conductor is degenerate, so the check stays for them.
    // Under half-domain symmetry the left wall is the symmetry plane, not
    // metal, so it is exempt from the test.
    const w = condRect && condRect.wallPEC;
    const enclosed = !!(w && (w.left || condRect.symmetry > 1) && w.right && w.top && w.bottom);
    if ((!roles || !roles.some(r => r.is_signal)) && !enclosed)
        errors.push('No signal conductor found in the geometry.');
    return errors;
}

function _cstr(G, str) {
    const len = G.lengthBytesUTF8(str) + 1;
    const ptr = G.stackAlloc(len);
    G.stringToUTF8(str, ptr, len);
    return ptr;
}
function _readDoubleArray(G, ptrPtr, nPtr) {
    const ptr = G.getValue(ptrPtr, 'i32'); const n = G.getValue(nPtr, 'i32');
    const arr = new Float64Array(n);
    for (let i = 0; i < n; i++) arr[i] = G.getValue(ptr + i * 8, 'double');
    G._gmshFree(ptr); return arr;
}
function _readIntArray(G, ptrPtr, nPtr) {
    const ptr = G.getValue(ptrPtr, 'i32'); const n = G.getValue(nPtr, 'i32');
    const arr = new Int32Array(n);
    for (let i = 0; i < n; i++) arr[i] = G.getValue(ptr + i * 4, 'i32');
    G._gmshFree(ptr); return arr;
}
function _rectOf(o) {
    if (o.xmin !== undefined) return { xmin: o.xmin, xmax: o.xmax, ymin: o.ymin, ymax: o.ymax };
    return { xmin: o.x_min, xmax: o.x_max, ymin: o.y_min, ymax: o.y_max };
}

// Absorb full-span boundary ground slabs into wall PEC BCs + clip the meshed domain
// (identical logic to geom_to_mesh._clipDomain). Exported for tests/test_geometry.js —
// the wall-absorption rule defines the effective cavity that the analytic mode tests
// (box_modes_test.mjs) compute their truth from.
export function _clipDomain(domain, conductors, boundaries, tol) {
    let { x_min: X0, x_max: X1, y_min: Y0, y_max: Y1 } = domain;
    const b = boundaries || ['open', 'open', 'open', 'gnd'];
    const wallPEC = { left: b[0] === 'gnd', right: b[1] === 'gnd', top: b[2] === 'gnd', bottom: b[3] === 'gnd' };
    // Metal thickness of each PEC wall: the absorbed slab's thickness (stacked
    // slabs add up), Infinity for a bare 'gnd' boundary. A slab against a 'gnd'
    // boundary is the ground plane itself, so its thickness is the slab's.
    const wallThick = { left: Infinity, right: Infinity, top: Infinity, bottom: Infinity };
    const absorb = (side, t) => {
        wallThick[side] = (Number.isFinite(wallThick[side]) ? wallThick[side] : 0) + t;
        wallPEC[side] = true;
    };
    let changed = true;
    while (changed) {
        changed = false;
        for (const c of conductors) {
            if (c.is_signal) continue;
            // Shaped grounds (e.g. a coax shield) are not full-span slabs and must not
            // be absorbed into a wall: their bounding box spans the domain but their
            // actual body does not fill it.
            if (c.shape) continue;
            const r = _rectOf(c);
            const touchesL = r.xmin <= X0 + tol, touchesR = r.xmax >= X1 - tol;
            const touchesB = r.ymin <= Y0 + tol, touchesT = r.ymax >= Y1 - tol;
            if (touchesB && touchesL && touchesR && r.ymax < Y1 - tol && r.ymax > Y0 + tol) { absorb('bottom', r.ymax - Y0); Y0 = r.ymax; changed = true; continue; }
            if (touchesT && touchesL && touchesR && r.ymin > Y0 + tol && r.ymin < Y1 - tol) { absorb('top', Y1 - r.ymin); Y1 = r.ymin; changed = true; continue; }
            if (touchesL && touchesB && touchesT && r.xmax < X1 - tol && r.xmax > X0 + tol) { absorb('left', r.xmax - X0); X0 = r.xmax; changed = true; continue; }
            if (touchesR && touchesB && touchesT && r.xmin > X0 + tol && r.xmin < X1 - tol) { absorb('right', X1 - r.xmin); X1 = r.xmin; changed = true; continue; }
        }
    }
    return { X0, X1, Y0, Y1, wallPEC, wallThick };
}

// Pre-mesh triangle-count estimate for the Distance->Threshold size field that
// buildOccMeshFromGeometry paints (see the field setup there). Triangles concentrate
// in the graded band around each painted curve: with surface size s and the size
// growing at the rate g (gradeRate) with the distance d beyond d = s, one side of the
// band holds L * integral(dd / (K * h(d)^2)) triangles, K being the mean triangle
// area in units of h^2. The exterior side integrates to L * (1 + 1/g) / (K * s). A
// meshed conductor interior is truncated at half the rect's thinner dimension and
// bounded by the rect area. The bulk adds area / (K * hCoarse^2).
//
// K = 0.33 was fitted on gmsh's Delaunay output for these fields (six meshes of a
// solder-masked differential GCPW, 4.7k..35k triangles: within -5..+20% wherever the
// Distance field's 200-point curve sampling resolves s, over-predicting only when s
// is far below the sample spacing of a long curve, which is the safe direction for
// a budget check). Curves on the domain outline get no exterior band (sigLenExt /
// gndLenExt). Bands of neighbouring curves overlap and are counted twice, so the
// estimate is biased high for closely stacked conductors (a cutout's ground slabs
// under a pair: +30%). `st` is the sizeStats object from a buildOccMeshFromGeometry
// call (statsOnly or not), `sz` the sizing {sizeMin, sizeMinGnd, hCoarse, gradeRate}.
export function estimateOccTriCount(st, sz) {
    const K = 0.33;
    const g = sz.gradeRate ?? 0.35;
    // Triangles per unit curve length on one side of the band, depth D available.
    const side = (s, D) => {
        if (!(s > 0)) return 0;
        if (D <= s) return D / (s * s * K);
        return (1 / s + (1 / s - 1 / (s + g * (D - s))) / g) / K;
    };
    const cls = (len, lenExt, area, thin, s) => {
        if (!(len > 0) || !(s > 0)) return 0;
        let n = lenExt * side(s, Infinity);
        if (st.interior && area > 0) {
            const D = Number.isFinite(thin) ? thin / 2 : Infinity;
            n += Math.min(len * side(s, D), area / (K * s * s));
        }
        return n;
    };
    const bulk = sz.hCoarse > 0 ? st.area / (K * sz.hCoarse * sz.hCoarse) : 0;
    return cls(st.sigLen, st.sigLenExt ?? st.sigLen, st.sigArea, st.sigThin, sz.sizeMin)
         + cls(st.gndLen, st.gndLenExt ?? st.gndLen, st.gndArea, st.gndThin, sz.sizeMinGnd) + bulk;
}

export function buildOccMeshFromGeometry(G, opts) {
    const conductors = opts.conductors || [];
    const dielectrics = opts.dielectrics || [];
    const domain = opts.domain;
    const boundaries = opts.boundaries;

    const diag = Math.hypot(domain.x_max - domain.x_min, domain.y_max - domain.y_min);
    const tol = diag * REL_TOL;
    const hFine = opts.hFine;
    const hCoarse = opts.hCoarse;
    const gradeRate = opts.gradeRate ?? 0.35;

    const symmetry = !!opts.symmetry;
    const domainShape = opts.domainShape || null;
    // Half-domain meshing of a shaped domain uses the x >= 0 half of every polygon.
    // Containment, however, always tests the FULL polygon (see condRects below).
    const meshOpts = { half: symmetry };

    let X0, X1, Y0, Y1, wallPEC, wallThick;
    if (domainShape) {
        // The meshed domain IS the shape, so there is no full-span slab to absorb and
        // no rectangle to clip: the bounds are simply the materialized polygon's bbox.
        ({ xmin: X0, xmax: X1, ymin: Y0, ymax: Y1 } = shapeBBox(domainShape, meshOpts));
        const b = boundaries || ['gnd', 'gnd', 'gnd', 'gnd'];
        wallPEC = { left: b[0] === 'gnd', right: b[1] === 'gnd', top: b[2] === 'gnd', bottom: b[3] === 'gnd' };
        wallThick = { left: Infinity, right: Infinity, top: Infinity, bottom: Infinity };
        if (symmetry) { X0 = 0; wallPEC.left = false; }
    } else {
        ({ X0, X1, Y0, Y1, wallPEC, wallThick } = _clipDomain(domain, conductors, boundaries, tol));
        if (symmetry) { X0 = 0; wallPEC = { ...wallPEC, left: false }; }
    }

    // Conductor rects clipped to the meshed domain (absorbed/outside ones dropped).
    const condRects = [], condRoles = [];
    for (const c of conductors) {
        if (c.shape) {
            // Bounds come from the shape's positive body. For a complement that is the
            // hole it surrounds, which is exactly the extent of its boundary in the mesh
            // (the shield itself has no meshed area). Downstream bbox prefilters stay
            // meaningful, while every real test goes through `shape`.
            //
            // `shape` is the full polygon even in a half-domain solve: using the half
            // would put the x=0 chord on the shield's boundary and make every symmetry
            // -plane node PEC, shorting the plane. `meshArea` carries the half instead,
            // because that is what the loss code multiplies back up by `symmetry`.
            const bb = shapeBBox(c.shape);
            condRects.push({
                xmin: bb.xmin, xmax: bb.xmax, ymin: bb.ymin, ymax: bb.ymax,
                shape: c.shape, meshArea: shapeArea(c, meshOpts),
            });
            condRoles.push({ is_signal: !!c.is_signal, polarity: c.polarity || 0, plating: c.plating || null });
            continue;
        }
        const r = _rectOf(c);
        const xmin = Math.max(r.xmin, X0), xmax = Math.min(r.xmax, X1);
        const ymin = Math.max(r.ymin, Y0), ymax = Math.min(r.ymax, Y1);
        if (xmax - xmin <= tol || ymax - ymin <= tol) continue;
        condRects.push({ xmin, xmax, ymin, ymax, meshArea: (xmax - xmin) * (ymax - ymin) });
        condRoles.push({ is_signal: !!c.is_signal, polarity: c.polarity || 0, plating: c.plating || null });
    }

    const stack = G.stackSave();
    const ierr = G.stackAlloc(4);
    const check = (label) => { const e = G.getValue(ierr, 'i32'); if (e) throw new Error(`occ: ${label} failed (err=${e})`); };

    G._gmshClear(ierr); check('Clear');
    G._gmshModelAdd(_cstr(G, 'occmesh'), ierr); check('ModelAdd');

    const addRect = (r) => {
        const t = G._gmshModelOccAddRectangle(r.xmin * OCC_SCALE, r.ymin * OCC_SCALE, 0,
            (r.xmax - r.xmin) * OCC_SCALE, (r.ymax - r.ymin) * OCC_SCALE, -1, 0, ierr);
        check('AddRectangle'); return t;
    };
    // Convex polygon as an OCC plane surface. The trimmed gmsh build exports no
    // AddDisk/AddCircle, but a circle is meshed as a polygon regardless (gmsh
    // discretizes curves, and refineTriMesh never revisits the CAD geometry), so the
    // polygon is the geometry rather than an approximation of it (see shapes.js).
    // meshSize 0 on every point: Mesh.MeshSizeFromPoints is disabled below, so sizing
    // comes exclusively from the background field.
    const addPolygon = (poly) => {
        const n = poly.length >> 1;
        const pts = new Array(n), lines = new Array(n);
        for (let i = 0; i < n; i++) {
            pts[i] = G._gmshModelOccAddPoint(poly[2 * i] * OCC_SCALE, poly[2 * i + 1] * OCC_SCALE, 0, 0, -1, ierr);
            check('AddPoint');
        }
        for (let i = 0; i < n; i++) {
            lines[i] = G._gmshModelOccAddLine(pts[i], pts[(i + 1) % n], -1, ierr);
            check('AddLine');
        }
        const lb = G.stackAlloc(n * 4);
        for (let i = 0; i < n; i++) G.setValue(lb + i * 4, lines[i], 'i32');
        const loop = G._gmshModelOccAddCurveLoop(lb, n, -1, ierr); check('AddCurveLoop');
        const wb = G.stackAlloc(4); G.setValue(wb, loop, 'i32');
        const face = G._gmshModelOccAddPlaneSurface(wb, 1, -1, ierr); check('AddPlaneSurface');
        return face;
    };
    const clipToDomain = (r) => ({
        xmin: Math.max(r.xmin, X0), xmax: Math.min(r.xmax, X1),
        ymin: Math.max(r.ymin, Y0), ymax: Math.min(r.ymax, Y1),
    });

    // Background "air" region (object) + all dielectric & conductor bodies (tools).
    const domTag = domainShape
        ? addPolygon(shapePoly(domainShape, meshOpts))
        : addRect({ xmin: X0, xmax: X1, ymin: Y0, ymax: Y1 });
    const toolTags = [];
    for (const d of dielectrics) {
        if (d.shape) {
            // A dielectric that IS the domain needs no tool — the background face
            // already covers it, and fragmenting a face against itself is degenerate.
            if (domainShape && d.shape === domainShape) continue;
            toolTags.push(addPolygon(shapePoly(d.shape, meshOpts)));
            continue;
        }
        const r = clipToDomain(_rectOf(d));
        if (r.xmax - r.xmin > tol && r.ymax - r.ymin > tol) toolTags.push(addRect(r));
    }
    for (const c of condRects) {
        // A complement conductor is a zero-area PEC shell whose boundary is already the
        // domain outline: it contributes no OCC geometry, only a condRects entry.
        if (c.shape && isComplement(c.shape)) continue;
        toolTags.push(c.shape ? addPolygon(shapePoly(c.shape, meshOpts)) : addRect(c));
    }

    // fragment([(2,domTag)], [(2,tool)...]) → conforming arrangement (remove originals).
    const obj = G.stackAlloc(8); G.setValue(obj, 2, 'i32'); G.setValue(obj + 4, domTag, 'i32');
    const toolBuf = G.stackAlloc(toolTags.length * 8);
    toolTags.forEach((t, i) => { G.setValue(toolBuf + i * 8, 2, 'i32'); G.setValue(toolBuf + i * 8 + 4, t, 'i32'); });
    const oDT = G.stackAlloc(4), oDTn = G.stackAlloc(4), oMap = G.stackAlloc(4), oMapN = G.stackAlloc(4), oMapNN = G.stackAlloc(4);
    // NB: the *_n args are the array length in INTS (2 per (dim,tag) pair), not pair count.
    G._gmshModelOccFragment(obj, 2, toolBuf, toolTags.length * 2, oDT, oDTn, oMap, oMapN, oMapNN, -1, 1, 1, ierr);
    check('Fragment');
    // free fragment outputs (outDimTags + the nested map)
    { const p = G.getValue(oDT, 'i32'); if (p) G._gmshFree(p); }
    { const mapPtr = G.getValue(oMap, 'i32'); const mapNPtr = G.getValue(oMapN, 'i32'); const nn = G.getValue(oMapNN, 'i32');
      for (let i = 0; i < nn; i++) { const sub = G.getValue(mapPtr + i * 4, 'i32'); if (sub) G._gmshFree(sub); }
      if (mapPtr) G._gmshFree(mapPtr); if (mapNPtr) G._gmshFree(mapNPtr); }

    G._gmshModelOccSynchronize(ierr); check('Synchronize');

    // ---- Size field: distance to conductor boundary curves → threshold ----
    // Identify conductor faces by enumerating all fragment faces and matching each face's
    // bounding-box centre to a conductor rect; collect their boundary curves for the field.
    // Signal and passive-ground conductor curves feed separate size fields. The
    // budget-coarsening loop in tri_backend.buildMesh relaxes the ground painting
    // (gndSizeScale) before it touches the signal hFine. With one union field a
    // large interior ground (e.g. the split ground planes of a cutout, ~the whole
    // domain width of curves painted at sizeMin) blows the initial-mesh budget and
    // the only lever was global coarsening, which left the signal trace with
    // elements larger than the trace itself.
    const sigCurves = new Set(), gndCurves = new Set();
    // Size-field statistics for the pre-mesh triangle-count estimate
    // (estimateOccTriCount): painted curve length, meshed conductor area and the
    // thinnest rect dimension per class, plus the meshed domain area.
    const sizeStats = { sigLen: 0, gndLen: 0, sigLenExt: 0, gndLenExt: 0,
                        sigArea: 0, gndArea: 0, sigThin: Infinity, gndThin: Infinity,
                        area: domainShape ? shapeArea({ shape: domainShape }, meshOpts) : (X1 - X0) * (Y1 - Y0),
                        interior: opts.meshConductorInterior !== false };
    {
        const fe = G.stackAlloc(4), feN = G.stackAlloc(4);
        G._gmshModelGetEntities(fe, feN, 2, ierr); check('GetEntities(2)');
        const faces = _readIntArray(G, fe, feN);   // (dim,tag) pairs
        const bb = [G.stackAlloc(8), G.stackAlloc(8), G.stackAlloc(8), G.stackAlloc(8), G.stackAlloc(8), G.stackAlloc(8)];
        // Bounding box of a model entity, in metres.
        const readBB = (dim, tag, label) => {
            G._gmshModelGetBoundingBox(dim, tag, bb[0], bb[1], bb[2], bb[3], bb[4], bb[5], ierr); check(label);
            return { x0: G.getValue(bb[0], 'double') / OCC_SCALE, y0: G.getValue(bb[1], 'double') / OCC_SCALE,
                     x1: G.getValue(bb[3], 'double') / OCC_SCALE, y1: G.getValue(bb[4], 'double') / OCC_SCALE };
        };
        // Only shapeless (non-rect) conductors take the face-bbox route below. A shaped
        // conductor's bounding box is far larger than its body (a complement shield's
        // spans the whole domain), so it would match the air face and collect the
        // entire domain outline as "conductor" curves, refining everything at sizeMin.
        // Roles ride along (condRoles is index-parallel to condRects) so each face's
        // curves land in the signal or ground set.
        const bboxRects = [], shapedRects = [];
        condRects.forEach((c, i) => {
            const sig = !!(condRoles[i] && condRoles[i].is_signal);
            if (c.shape) shapedRects.push(c); else bboxRects.push({ c, sig });
        });
        // A conductor sub-face's bounding box is contained in its conductor rect.
        // Testing only the bbox CENTER misclassifies a face-with-hole: the air
        // face surrounding a centered conductor (e.g. symmetric stripline) has
        // its bbox center inside the conductor rect, so the whole domain outline
        // would be collected as "conductor" curves and refined at sizeMin.
        // OCC enlarges bounding boxes by Precision::Confusion, the containment
        // pad must absorb that, with a 2x margin.
        const OCC_BBOX_PAD = 2e-7 / OCC_SCALE;
        const pad = Math.max(tol, OCC_BBOX_PAD);
        // Is the bbox (x0,y0)-(x1,y1) contained in this rect, within the OCC pad? Shared
        // with the conductor-interior removal below so both use one padding rule.
        const bboxInRect = (c, x0, y0, x1, y1) => {
            return x0 > c.xmin - pad && x1 < c.xmax + pad &&
                   y0 > c.ymin - pad && y1 < c.ymax + pad;
        };
        const inCond = (x0, y0, x1, y1) => bboxRects.find(e => bboxInRect(e.c, x0, y0, x1, y1)) || null;
        for (let i = 0; i < faces.length; i += 2) {
            const ftag = faces[i + 1];
            const b = readBB(2, ftag, 'GetBoundingBox');
            const hit = inCond(b.x0, b.y0, b.x1, b.y1);
            if (!hit) continue;
            const fbuf = G.stackAlloc(8); G.setValue(fbuf, 2, 'i32'); G.setValue(fbuf + 4, ftag, 'i32');
            const cb = G.stackAlloc(4), cbN = G.stackAlloc(4);
            G._gmshModelGetBoundary(fbuf, 2, cb, cbN, 0, 0, 0, ierr); check('GetBoundary');
            const curves = _readIntArray(G, cb, cbN);
            const set = hit.sig ? sigCurves : gndCurves;
            for (let k = 1; k < curves.length; k += 2) set.add(Math.abs(curves[k]));
        }

        // Shaped conductors are classified per CURVE instead of per face. This is what
        // makes a complement shield's surface get refined at all: it has zero meshed
        // area, so no face-based rule can ever find it, yet it is a real PEC surface
        // carrying half the conductor loss (alpha_c goes as 1/a + 1/b).
        //
        // Every curve here is a straight sub-segment of a polygon side, so its bounding
        // -box CENTRE is the segment midpoint and lies exactly on the shape boundary.
        // (Testing bbox corners would not work — for a non-axis-aligned segment they
        // are off the line.)
        if (shapedRects.length) {
            const ce = G.stackAlloc(4), ceN = G.stackAlloc(4);
            G._gmshModelGetEntities(ce, ceN, 1, ierr); check('GetEntities(1)');
            const curves = _readIntArray(G, ce, ceN);   // (dim,tag) pairs
            for (let i = 0; i < curves.length; i += 2) {
                const ctag = curves[i + 1];
                const b = readBB(1, ctag, 'GetBoundingBox(1)');
                const cx = (b.x0 + b.x1) / 2, cy = (b.y0 + b.y1) / 2;
                // Shaped conductors always take the fine (signal) field: a coax
                // shield is the return conductor and its surface bounds the whole
                // field region, so it never qualifies for the ground relaxation.
                for (const c of shapedRects) {
                    if (Math.abs(shapeSignedDist(c.shape, cx, cy)) < tol) { sigCurves.add(ctag); break; }
                }
            }
        }

        // Curve lengths per class (straight segments: the bbox diagonal), for the
        // pre-mesh triangle-count estimate (estimateOccTriCount). A curve lying on
        // the meshed-domain outline (a ground slab's bottom face on y = Y0, a via
        // slab's outer face on x = X1) has no exterior side to mesh, so it counts
        // toward the interior term only (`*LenExt` excludes it).
        const bTol = OCC_BBOX_PAD + tol;
        const curveLen = (ctag) => {
            const { x0, y0, x1, y1 } = readBB(1, ctag, 'GetBoundingBox(1) len');
            const onOutline = !domainShape && (
                (x1 - x0 < bTol && (Math.abs(x0 - X0) < bTol || Math.abs(x1 - X1) < bTol)) ||
                (y1 - y0 < bTol && (Math.abs(y0 - Y0) < bTol || Math.abs(y1 - Y1) < bTol)));
            return { len: Math.hypot(x1 - x0, y1 - y0), onOutline };
        };
        for (const c of sigCurves) { const { len, onOutline } = curveLen(c); sizeStats.sigLen += len; if (!onOutline) sizeStats.sigLenExt += len; }
        for (const c of gndCurves) { const { len, onOutline } = curveLen(c); sizeStats.gndLen += len; if (!onOutline) sizeStats.gndLenExt += len; }
        condRects.forEach((c, i) => {
            if (c.shape && isComplement(c.shape)) return;
            const thin = c.shape ? Infinity : Math.min(c.xmax - c.xmin, c.ymax - c.ymin);
            if (condRoles[i] && condRoles[i].is_signal) {
                sizeStats.sigArea += c.meshArea; sizeStats.sigThin = Math.min(sizeStats.sigThin, thin);
            } else {
                sizeStats.gndArea += c.meshArea; sizeStats.gndThin = Math.min(sizeStats.gndThin, thin);
            }
        });
        // Statistics-only call (tri_backend's pre-mesh budget sizing): the OCC model
        // is built and classified exactly as for a real mesh, but nothing is meshed.
        if (opts.statsOnly) {
            G.stackRestore(stack);
            return { sizeStats, nGndCurves: gndCurves.size };
        }

        // ---- Optionally drop the conductor INTERIORS from the meshed model ----
        // Only the MQS volume eddy-current solve ever looks inside the metal. On the
        // perturbation path those elements carry NO degrees of freedom at all — every
        // node and edge inside a conductor is PEC, and their triangles get no face DOFs
        // — so they cost mesh generation, storage and refinement passes for nothing.
        //
        // Removing the face (recursive = 0) keeps its boundary curves, which the
        // surrounding face already shares, so the conductor simply becomes a hole and
        // the surface edges the loss integral needs are untouched.
        //
        // Must run AFTER the curve collection above (which finds conductor curves via
        // their faces) and complements are skipped — they have no face to remove, and
        // their bounding box would match the surrounding dielectric.
        if (opts.meshConductorInterior === false && condRects.length) {
            const holes = [];
            for (let i = 0; i < faces.length; i += 2) {
                const ftag = faces[i + 1];
                const b = readBB(2, ftag, 'GetBoundingBox');
                const hit = condRects.some(c =>
                    !(c.shape && isComplement(c.shape)) && bboxInRect(c, b.x0, b.y0, b.x1, b.y1));
                if (hit) holes.push(ftag);
            }
            let dangling = 0;
            if (holes.length) {
                const hb = G.stackAlloc(holes.length * 8);
                holes.forEach((t, i) => { G.setValue(hb + i * 8, 2, 'i32'); G.setValue(hb + i * 8 + 4, t, 'i32'); });
                G._gmshModelOccRemove(hb, holes.length * 2, 0, ierr); check('OccRemove(conductor faces)');
                G._gmshModelOccSynchronize(ierr); check('Synchronize(after remove)');

                // Removing a face non-recursively keeps ALL its curves, including any
                // that ran through its interior — the fragment splits a line crossing a
                // conductor into three, and the middle piece is now bounded by nothing.
                // gmsh still meshes such a dangling curve, seeding stray nodes inside the
                // metal that belong to no triangle (harmless for the solve — they get no
                // DOFs — but they inflate the node count and show up in the mesh overlay).
                // Recursive removal is not an option: these curves are shared with the
                // surrounding face. Delete just the ones strictly inside a conductor.
                const ce2 = G.stackAlloc(4), ceN2 = G.stackAlloc(4);
                G._gmshModelGetEntities(ce2, ceN2, 1, ierr); check('GetEntities(1) post-remove');
                const cs = _readIntArray(G, ce2, ceN2);
                const kill = [];
                for (let i = 0; i < cs.length; i += 2) {
                    const ctag = cs[i + 1];
                    const b = readBB(1, ctag, 'GetBoundingBox(1) post-remove');
                    const cx = (b.x0 + b.x1) / 2, cy = (b.y0 + b.y1) / 2;
                    for (const c of condRects) {
                        if (c.shape && isComplement(c.shape)) continue;
                        // STRICTLY inside (−tol): a curve ON the conductor boundary is a
                        // real surface the loss integral needs, and must be kept.
                        if (shapeContains(c, cx, cy, -tol)) { kill.push(ctag); break; }
                    }
                }
                if (kill.length) {
                    const kb = G.stackAlloc(kill.length * 8);
                    kill.forEach((t, i) => { G.setValue(kb + i * 8, 1, 'i32'); G.setValue(kb + i * 8 + 4, t, 'i32'); });
                    G._gmshModelOccRemove(kb, kill.length * 2, 0, ierr); check('OccRemove(dangling curves)');
                    G._gmshModelOccSynchronize(ierr); check('Synchronize(after curve remove)');
                }
                dangling = kill.length;
            }
            if (globalThis.__OCC_DEBUG__)
                console.error('[occ] removed', holes.length, 'conductor face(s),', dangling, 'dangling curve(s)');
        }
    }

    // Conductor-surface element size — finer than global hFine (see field setup below).
    const surfScale = opts.occSurfScale ?? 0.35;
    const sizeMin = hFine * surfScale;
    // Passive-ground surface size: sizeMin relaxed by gndSizeScale (default 1 = the
    // grounds mesh exactly as fine as the signals, tri_backend raises it only when the
    // initial mesh overruns its triangle budget). At the hCoarse cap the ground field
    // is inert, nearby signal grading (the signal Threshold is a background field
    // evaluated everywhere, ground surfaces included) and thin-rect conformity then
    // set the ground resolution, and the adaptive passes sharpen what matters.
    const gndScale = Math.max(1, opts.gndSizeScale ?? 1);
    const sizeMinGnd = Math.min(sizeMin * gndScale, hCoarse);
    Object.assign(sizeStats, { sizeMin, sizeMinGnd, hFine, hCoarse, gradeRate });
    if (globalThis.__OCC_DEBUG__) console.error('[occ] condRects', condRects.length, 'sigCurves', sigCurves.size, 'gndCurves', gndCurves.size, 'hFine', hFine, 'sizeMin', sizeMin, 'sizeMinGnd', sizeMinGnd, 'hCoarse', hCoarse, 'stats', JSON.stringify(sizeStats));
    {
        // sizeMin (conductor-surface element size) is finer than the global hFine so the
        // current distribution (skin/proximity loss) and a narrow inter-trace gap (only
        // ~1·hFine wide) are resolved WITHOUT globally refining the far field. Adaptive ZZ
        // passes sharpen it further.
        const addCurveField = (curveSet, sMin) => {
            const fd = G._gmshModelMeshFieldAdd(_cstr(G, 'Distance'), -1, ierr); check('FieldAdd Distance');
            const cl = [...curveSet];
            const clBuf = G.stackAlloc(cl.length * 8);
            cl.forEach((v, i) => G.setValue(clBuf + i * 8, v, 'double'));
            G._gmshModelMeshFieldSetNumbers(fd, _cstr(G, 'CurvesList'), clBuf, cl.length, ierr); check('CurvesList');
            G._gmshModelMeshFieldSetNumber(fd, _cstr(G, 'Sampling'), 200, ierr);
            const ft = G._gmshModelMeshFieldAdd(_cstr(G, 'Threshold'), -1, ierr); check('FieldAdd Threshold');
            G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'InField'), fd, ierr);
            G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'SizeMin'), sMin * OCC_SCALE, ierr);
            G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'SizeMax'), hCoarse * OCC_SCALE, ierr);
            G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'DistMin'), sMin * OCC_SCALE, ierr);
            G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'DistMax'), (hFine + (hCoarse - hFine) / gradeRate) * OCC_SCALE, ierr);
            return ft;
        };
        const fields = [];
        if (sigCurves.size > 0) fields.push(addCurveField(sigCurves, sizeMin));
        // A fully relaxed ground field would be constant hCoarse — skip it.
        if (gndCurves.size > 0 && sizeMinGnd < hCoarse) fields.push(addCurveField(gndCurves, sizeMinGnd));
        if (fields.length === 1) {
            G._gmshModelMeshFieldSetAsBackgroundMesh(fields[0], ierr); check('SetBackground');
        } else if (fields.length > 1) {
            // With gndSizeScale = 1 both Thresholds have identical parameters, so
            // Min over the two curve subsets equals the old single union field.
            const fm = G._gmshModelMeshFieldAdd(_cstr(G, 'Min'), -1, ierr); check('FieldAdd Min');
            const fBuf = G.stackAlloc(fields.length * 8);
            fields.forEach((v, i) => G.setValue(fBuf + i * 8, v, 'double'));
            G._gmshModelMeshFieldSetNumbers(fm, _cstr(G, 'FieldsList'), fBuf, fields.length, ierr); check('FieldsList');
            G._gmshModelMeshFieldSetAsBackgroundMesh(fm, ierr); check('SetBackground');
        }
    }

    // Field-only sizing (don't let point/curvature/boundary sizing fight it).
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeExtendFromBoundary'), 0, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeFromPoints'), 0, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeFromCurvature'), 0, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeMin'), Math.min(hFine * 0.5, sizeMin * 0.5) * OCC_SCALE, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeMax'), hCoarse * OCC_SCALE, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.Algorithm'), 5, ierr);
    // One Laplacian smoothing pass (gmsh's default). Ten passes cost ~0.2 ms per
    // triangle in this WASM build (1.1 s of a 1.6 s build at 5k triangles) and
    // measured no quality gain over one: same maxQ/avgQ on a 1 µm-trace GCPW,
    // slightly better on a 35 µm one, results within 0.01%.
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.Smoothing'), 1, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.RandomSeed'), 1, ierr);
    // Caller overrides (experimental tuning of the gmsh mesher, e.g. Mesh.Optimize).
    if (opts.gmshOptions) {
        for (const [k, v] of Object.entries(opts.gmshOptions)) {
            G._gmshOptionSetNumber(_cstr(G, k), v, ierr); check('opt ' + k);
        }
    }

    G._gmshModelMeshGenerate(2, ierr); check('MeshGenerate');

    // ---- Extract nodes ----
    const nodePtrPtr = G.stackAlloc(4), nodeNPtr = G.stackAlloc(4);
    const coordPtrPtr = G.stackAlloc(4), coordNPtr = G.stackAlloc(4);
    const paramPtrPtr = G.stackAlloc(4), paramNPtr = G.stackAlloc(4);
    G._gmshModelMeshGetNodes(nodePtrPtr, nodeNPtr, coordPtrPtr, coordNPtr, paramPtrPtr, paramNPtr, -1, -1, 0, 0, ierr);
    check('GetNodes');
    const gmshTags = _readIntArray(G, nodePtrPtr, nodeNPtr);
    const coords = _readDoubleArray(G, coordPtrPtr, coordNPtr);
    const paramPtr = G.getValue(paramPtrPtr, 'i32'); if (paramPtr) G._gmshFree(paramPtr);

    const nNodes = gmshTags.length;
    const tagToIdx = new Map();
    for (let i = 0; i < nNodes; i++) tagToIdx.set(gmshTags[i], i);
    const nodes = new Float64Array(2 * nNodes);
    for (let i = 0; i < nNodes; i++) { nodes[2 * i] = coords[3 * i] / OCC_SCALE; nodes[2 * i + 1] = coords[3 * i + 1] / OCC_SCALE; }

    // Snap near-wall nodes exactly onto the domain walls
    //
    // OCC's boolean fragmentation works to Precision::Confusion() (1e-7 model
    // units, 1e-7 / OCC_SCALE in metres), so a nearly-degenerate fragment (for
    // example solder mask) can emit nodes slightly off the wall. Adaptive
    // passes can then drag the node causing it to cause a void to appear near
    // the symmetry plane. Shaped domains (coax) are excluded.
    if (!domainShape) {
        const snapTol = Math.max(tol, 2e-7 / OCC_SCALE);
        for (let i = 0; i < nNodes; i++) {
            const x = nodes[2 * i], y = nodes[2 * i + 1];
            if (Math.abs(x - X0) < snapTol) nodes[2 * i] = X0;
            else if (Math.abs(x - X1) < snapTol) nodes[2 * i] = X1;
            if (Math.abs(y - Y0) < snapTol) nodes[2 * i + 1] = Y0;
            else if (Math.abs(y - Y1) < snapTol) nodes[2 * i + 1] = Y1;
        }
    }

    // ---- Extract triangles ----
    const etPtrPtr = G.stackAlloc(4), etNPtr = G.stackAlloc(4);
    const enPtrPtr = G.stackAlloc(4), enNPtr = G.stackAlloc(4);
    G._gmshModelMeshGetElementsByType(2, etPtrPtr, etNPtr, enPtrPtr, enNPtr, -1, 0, 1, ierr);
    check('GetElementsByType');
    const nTris = G.getValue(etNPtr, 'i32');
    const enPtr = G.getValue(enPtrPtr, 'i32');
    G._gmshFree(G.getValue(etPtrPtr, 'i32'));
    const tris = new Int32Array(3 * nTris);
    // tagToIdx.get() returns undefined for an unknown tag, which the Int32Array would
    // silently coerce to node 0 (corrupt geometry). Flag it and fail loudly below.
    let badTag = false;
    for (let t = 0; t < nTris; t++)
        for (let k = 0; k < 3; k++) {
            const idx = tagToIdx.get(G.getValue(enPtr + (3 * t + k) * 4, 'i32'));
            if (idx === undefined) badTag = true;
            tris[3 * t + k] = idx | 0;
        }
    G._gmshFree(enPtr);

    G.stackRestore(stack);

    // Fail loudly on a collapsed/empty or corrupt mesh rather than feeding empty/garbage
    // arrays into the FEM solve (where they surface as opaque NaN/Inf or a heap abort).
    if (nTris === 0)
        throw new Error('Full-wave mesher produced no triangles — the geometry may be degenerate or the OCC fragment collapsed.');
    if (badTag)
        throw new Error('Full-wave mesh extraction referenced an unknown node tag (corrupt mesh).');

    // ---- CCW winding ----
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const ax = nodes[2 * v0], ay = nodes[2 * v0 + 1];
        const bx = nodes[2 * v1], by = nodes[2 * v1 + 1];
        const cx = nodes[2 * v2], cy = nodes[2 * v2 + 1];
        if ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay) < 0) { tris[3 * t + 1] = v2; tris[3 * t + 2] = v1; }
    }

    // ---- Edges, triEdges, triSigns ----
    const edgeMap = new Map(), edgeList = [];
    const triEdges = new Int32Array(3 * nTris), triSigns = new Int8Array(3 * nTris);
    const localEdges = [[0, 1], [1, 2], [2, 0]];
    for (let t = 0; t < nTris; t++) {
        for (let le = 0; le < 3; le++) {
            const na = tris[3 * t + localEdges[le][0]], nb = tris[3 * t + localEdges[le][1]];
            const n0 = na < nb ? na : nb, n1 = na < nb ? nb : na;
            const key = n0 * (nNodes + 1) + n1;
            let eIdx = edgeMap.get(key);
            if (eIdx === undefined) { eIdx = edgeList.length >> 1; edgeList.push(n0, n1); edgeMap.set(key, eIdx); }
            triEdges[3 * t + le] = eIdx;
            triSigns[3 * t + le] = (na === n0) ? 1 : -1;
        }
    }
    const nEdges = edgeList.length >> 1;
    const edges = new Int32Array(edgeList);

    // ---- Per-triangle material maps (centroid, last-match-wins over dielectrics) ----
    const { epsMap, lossMap } = tagMaterials({ nodes, tris, nTris }, dielectrics, tol);

    // ---- condRect for the FEM freedom map ----
    const condRect = {
        rects: condRects, rectRoles: condRoles,
        xmin_domain: X0, xmax_domain: X1, ymin_domain: Y0, ymax_domain: Y1,
        // Absolute geometric tolerance (domain-diagonal relative). Shaped containment
        // tests need it: an exact-zero tolerance would reject boundary nodes that are
        // one ULP off the polygon side after the mesher's own arithmetic.
        geomTol: tol,
        wallPEC, wallThick, symmetry: symmetry ? 2 : 1,
        xmin: condRects.length ? Math.min(...condRects.map(c => c.xmin)) : X0,
        xmax: condRects.length ? Math.max(...condRects.map(c => c.xmax)) : X1,
        ymin: condRects.length ? Math.min(...condRects.map(c => c.ymin)) : Y0,
        ymax: condRects.length ? Math.max(...condRects.map(c => c.ymax)) : Y1,
    };
    condRect.rects.forEach(r => { r.symmetry = condRect.symmetry; });

    // Constraint metadata (checkMeshQuality / refinement smoothing): a node on
    // one of these lines may only move ALONG it. Constrained lines are the
    // CONDUCTOR faces (all four sides, over the rect's own extent — the freedom
    // map classifies PEC nodes/edges by the exact rect coordinates, and the
    // loss/MQS surface integrals live there) plus the symmetry plane.
    // Dielectric interfaces are deliberately NOT constrained: epsMap is re-tagged
    // by centroid after every refinement pass, so a node sliding off a dielectric
    // line only displaces that interface locally by ~one element (second-order),
    // while forbidding the move can leave unfixable slivers that ill-condition
    // the whole solve (fuzzer: Qmax 183 on solder-mask cases when constrained).
    const constraintYRanges = {};
    const addYR = (y, lo, hi) => { const e = constraintYRanges[y]; constraintYRanges[y] = e ? [Math.min(e[0], lo), Math.max(e[1], hi)] : [lo, hi]; };
    const constraintXRanges = {};
    const addXR = (x, lo, hi) => { const e = constraintXRanges[x]; constraintXRanges[x] = e ? [Math.min(e[0], lo), Math.max(e[1], hi)] : [lo, hi]; };
    // Shaped conductors constrain arbitrary line SEGMENTS (their polygon sides) rather
    // than axis-aligned lines. They must contribute no addXR/addYR: their bounding box
    // is not a real surface, and for a complement it spans the whole domain, which
    // would pin the domain walls for no reason.
    const constraintSegments = [];
    for (const c of condRects) {
        if (c.shape) { constraintSegments.push(...shapeSegments(c, meshOpts)); continue; }
        addYR(c.ymin, c.xmin, c.xmax); addYR(c.ymax, c.xmin, c.xmax);
        addXR(c.xmin, c.ymin, c.ymax); addXR(c.xmax, c.ymin, c.ymax);
    }
    if (symmetry) addXR(X0, Y0, Y1);
    const constraintXs = Object.keys(constraintXRanges).map(Number);

    return {
        nodes, tris, edges, triEdges, triSigns, nNodes, nTris, nEdges,
        constraintYRanges, constraintXRanges, constraintXs, constraintSegments,
        condRect, epsMap, lossMap,
        // Whether conductor interiors carry triangles. The MQS volume eddy solve is the
        // only consumer of them and would silently return a wrong R on a mesh built
        // without them, so it is recorded on the mesh rather than left implicit in the
        // caller's (separately derived) applicability test. See TriBackend._modeAtFreq.
        condInteriorMeshed: opts.meshConductorInterior !== false,
        // How many passive-ground boundary curves fed the (relaxable) ground size
        // field — tri_backend's budget loop only tries gndSizeScale when nonzero.
        nGndCurves: gndCurves.size,
        // Size-field statistics + the sizing this mesh was built with, for
        // estimateOccTriCount (tri_backend's budget pre-sizing / diagnostics).
        sizeStats,
        meshedDomain: { X0, X1, Y0, Y1, wallPEC, wallThick },
    };
}
