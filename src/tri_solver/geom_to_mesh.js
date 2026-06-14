// Generalized geo-kernel gmsh mesher for arbitrary 2D transmission-line
// cross-sections, built from the backend-neutral geometry the app already
// produces: a list of axis-aligned Conductor and Dielectric rectangles
// (microstrip.js / broadside_stripline.js `_build_geometry_lists()`).
//
// The vendored gmsh WASM is a minimal geo-kernel build: no OpenCASCADE (no
// boolean fragment), no Mesh.Field, no mesh `embed`, no physical groups. To make
// the mesh CONFORM to every conductor and dielectric boundary, each such boundary
// must be a real surface edge. We therefore build a rectilinear ARRANGEMENT:
//   - collect every distinct x and y edge coordinate (conductors + dielectrics +
//     domain), giving a non-uniform grid of cells;
//   - emit a shared point/line grid (lines reused across cells) and one plane
//     surface per non-conductor cell;
//   - conductor cells are left as holes (their boundary lines still exist, shared
//     by neighbouring air cells), unless meshConductorInterior is set.
// gmsh conforms to all shared lines, so ε jumps and PEC faces are mesh-aligned.
//
// Full-span boundary ground slabs (a ground rectangle covering an entire domain
// edge, e.g. the bottom ground) are absorbed into a wall PEC boundary condition
// and the meshed domain is clipped to their inner face — this avoids meshing the
// ground interior and matches the ms2d "ground = wall" convention.
//
// Output: the standard solver mesh object plus `condRect` (with per-rect roles)
// and per-triangle `epsMap`/`lossMap`, ready for buildTriFreedomMap / assembleTriFEM.

import { generateGmshMesh } from './gmsh_mesh.js';  // re-export initGmsh via this module's consumer

const REL_TOL = 1e-9;

function _cstr(G, str) {
    const len = G.lengthBytesUTF8(str) + 1;
    const ptr = G.stackAlloc(len);
    G.stringToUTF8(str, ptr, len);
    return ptr;
}
function _readDoubleArray(G, ptrPtr, nPtr) {
    const ptr = G.getValue(ptrPtr, 'i32');
    const n = G.getValue(nPtr, 'i32');
    const arr = new Float64Array(n);
    for (let i = 0; i < n; i++) arr[i] = G.getValue(ptr + i * 8, 'double');
    G._gmshFree(ptr);
    return arr;
}
function _readIntArray(G, ptrPtr, nPtr) {
    const ptr = G.getValue(ptrPtr, 'i32');
    const n = G.getValue(nPtr, 'i32');
    const arr = new Int32Array(n);
    for (let i = 0; i < n; i++) arr[i] = G.getValue(ptr + i * 4, 'i32');
    G._gmshFree(ptr);
    return arr;
}

// Normalize a Conductor/Dielectric-like object (has x_min/x_max/y_min/y_max
// getters) or a plain {xmin,xmax,ymin,ymax} into a plain rect.
function _rectOf(o) {
    if (o.xmin !== undefined) return { xmin: o.xmin, xmax: o.xmax, ymin: o.ymin, ymax: o.ymax };
    return { xmin: o.x_min, xmax: o.x_max, ymin: o.y_min, ymax: o.y_max };
}

// Sorted unique coordinate list with absolute-tolerance dedup.
function _uniqSorted(vals, tol) {
    const s = [...vals].sort((a, b) => a - b);
    const out = [];
    for (const v of s) {
        if (out.length === 0 || v - out[out.length - 1] > tol) out.push(v);
    }
    return out;
}

// Determine which full domain walls are grounded, and clip the meshed domain to
// absorb full-span boundary ground slabs. boundaries = [left,right,top,bottom]
// each 'open' | 'gnd'. Returns { X0,X1,Y0,Y1, wallPEC:{left,right,top,bottom} }.
function _clipDomain(domain, conductors, boundaries, tol) {
    let { x_min: X0, x_max: X1, y_min: Y0, y_max: Y1 } = domain;
    const b = boundaries || ['open', 'open', 'open', 'gnd'];
    const wallPEC = { left: b[0] === 'gnd', right: b[1] === 'gnd', top: b[2] === 'gnd', bottom: b[3] === 'gnd' };

    // Iteratively absorb ground rects that span the full current domain extent on
    // a wall (and are thin against the perpendicular extent). Repeat because
    // absorbing one slab can expose another (e.g. corner stacks).
    let changed = true;
    while (changed) {
        changed = false;
        for (const c of conductors) {
            if (c.is_signal) continue;
            const r = _rectOf(c);
            // "touches" = reaches the (clipped) wall or extends beyond it — so a
            // ground that runs past an already-absorbed wall still counts (e.g. a
            // side ground whose ymin = −t_gnd sits below an absorbed bottom ground).
            const touchesL = r.xmin <= X0 + tol, touchesR = r.xmax >= X1 - tol;
            const touchesB = r.ymin <= Y0 + tol, touchesT = r.ymax >= Y1 - tol;
            // Bottom slab: full width, sits at bottom, thin in y
            if (touchesB && touchesL && touchesR && r.ymax < Y1 - tol && r.ymax > Y0 + tol) {
                Y0 = r.ymax; wallPEC.bottom = true; changed = true; continue;
            }
            if (touchesT && touchesL && touchesR && r.ymin > Y0 + tol && r.ymin < Y1 - tol) {
                Y1 = r.ymin; wallPEC.top = true; changed = true; continue;
            }
            if (touchesL && touchesB && touchesT && r.xmax < X1 - tol && r.xmax > X0 + tol) {
                X0 = r.xmax; wallPEC.left = true; changed = true; continue;
            }
            if (touchesR && touchesB && touchesT && r.xmin > X0 + tol && r.xmin < X1 - tol) {
                X1 = r.xmin; wallPEC.right = true; changed = true; continue;
            }
        }
    }
    return { X0, X1, Y0, Y1, wallPEC };
}

function pointInPoly(x, y, poly) {
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
        const xi = poly[i][0], yi = poly[i][1], xj = poly[j][0], yj = poly[j][1];
        if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) inside = !inside;
    }
    return inside;
}

// Outer boundary polygon(s) of the air region = band rect [X0,X1]×[ya,yb] minus
// the wall-flush conductor rectangles (touchers). Traced from a local cell mask
// over the touchers' edge coordinates, so arbitrary notches (gnd cut, coplanar
// grounds, vias) yield clean rectilinear polygons. Returns an array of CCW rings.
// forcedX: x-coordinates at which the horizontal band-interface edges (y=ya, y=yb)
// must be split so that ADJACENT bands share identical lines (conformity). These
// are the x-edges of every conductor touching either interface, from either side.
function traceAirBoundary(X0, X1, ya, yb, touchers, tol, forcedXbottom = [], forcedXtop = []) {
    // A collinear vertex on the bottom (top) interface is kept only if it lies at a
    // forced split for THAT interface — so each interface is split to match only its
    // own adjacent band, never the opposite one.
    const inSet = (x, set) => set.some(v => Math.abs(v - x) < tol);
    const keepIface = (x, y) => (Math.abs(y - ya) < tol && inSet(x, forcedXbottom)) ||
                                (Math.abs(y - yb) < tol && inSet(x, forcedXtop));
    if (touchers.length === 0 && forcedXbottom.length === 0 && forcedXtop.length === 0)
        return [[[X0, ya], [X1, ya], [X1, yb], [X0, yb]]];
    const gx = _uniqSorted([X0, X1,
        ...forcedXbottom.filter(x => x > X0 - tol && x < X1 + tol),
        ...forcedXtop.filter(x => x > X0 - tol && x < X1 + tol),
        ...touchers.flatMap(t => [Math.max(t.xmin, X0), Math.min(t.xmax, X1)])], tol);
    const gy = _uniqSorted([ya, yb, ...touchers.flatMap(t => [Math.max(t.ymin, ya), Math.min(t.ymax, yb)])], tol);
    const nx = gx.length - 1, ny = gy.length - 1;
    const isAir = (i, j) => {
        if (i < 0 || j < 0 || i >= nx || j >= ny) return false;
        const xc = (gx[i] + gx[i + 1]) / 2, yc = (gy[j] + gy[j + 1]) / 2;
        for (const t of touchers)
            if (xc > t.xmin - tol && xc < t.xmax + tol && yc > t.ymin - tol && yc < t.ymax + tol) return false;
        return true;
    };
    // Directed boundary segments (air on the left → CCW outer rings).
    const segs = new Map();  // "x,y" start -> [x,y] end
    const k = (x, y) => `${x},${y}`;
    const addSeg = (x1, y1, x2, y2) => segs.set(k(x1, y1), [x2, y2]);
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
        if (!isAir(i, j)) continue;
        if (!isAir(i, j - 1)) addSeg(gx[i], gy[j], gx[i + 1], gy[j]);       // bottom +x
        if (!isAir(i, j + 1)) addSeg(gx[i + 1], gy[j + 1], gx[i], gy[j + 1]); // top -x
        if (!isAir(i - 1, j)) addSeg(gx[i], gy[j + 1], gx[i], gy[j]);       // left -y
        if (!isAir(i + 1, j)) addSeg(gx[i + 1], gy[j], gx[i + 1], gy[j + 1]); // right +y
    }
    const rings = [];
    const seen = new Set();
    for (const start of segs.keys()) {
        if (seen.has(start)) continue;
        const ring = [];
        let cur = start;
        while (cur && !seen.has(cur)) {
            seen.add(cur);
            const [px, py] = cur.split(',').map(Number);
            ring.push([px, py]);
            const nxt = segs.get(cur);
            if (!nxt) break;
            cur = k(nxt[0], nxt[1]);
        }
        // drop collinear points, EXCEPT those on a band interface (y=ya/yb): those
        // splits must be kept so adjacent bands share identical interface lines.
        const simp = [];
        for (let i = 0; i < ring.length; i++) {
            const a = ring[(i - 1 + ring.length) % ring.length], b = ring[i], c = ring[(i + 1) % ring.length];
            const cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            if (Math.abs(cross) > tol * tol || keepIface(b[0], b[1])) simp.push(b);
        }
        if (simp.length >= 3) rings.push(simp);
    }
    return rings.length ? rings : [[[X0, ya], [X1, ya], [X1, yb], [X0, yb]]];
}

// Per-triangle ε / loss maps by centroid (last-match-wins over the dielectric
// list, matching _setup_geometry). Reusable: call again after refineTriMesh.
export function tagMaterials(mesh, dielectrics, tol) {
    const { nodes, tris, nTris } = mesh;
    if (tol === undefined) tol = 1e-12;
    const epsMap = new Array(nTris), lossMap = new Array(nTris);
    const rects = dielectrics.map(d => ({ r: _rectOf(d), er: d.epsilon_r, tand: d.tan_delta || 0 }));
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const xc = (nodes[2 * v0] + nodes[2 * v1] + nodes[2 * v2]) / 3;
        const yc = (nodes[2 * v0 + 1] + nodes[2 * v1 + 1] + nodes[2 * v2 + 1]) / 3;
        let er = 1.0, tand = 0.0;
        for (const d of rects) {
            if (xc >= d.r.xmin - tol && xc <= d.r.xmax + tol && yc >= d.r.ymin - tol && yc <= d.r.ymax + tol) {
                er = d.er; tand = d.tand;
            }
        }
        epsMap[t] = { re: er, im: 0 };
        lossMap[t] = { re: er * tand, im: 0 };
    }
    return { epsMap, lossMap };
}

// Build a generalized triangular mesh from geometry rectangle lists.
//
// opts:
//   conductors  — array of Conductor-like rects (x_min/x_max/y_min/y_max, is_signal, polarity, plating)
//   dielectrics — array of Dielectric-like rects (+ epsilon_r, tan_delta), ordered low→high priority
//                 (later entries win on overlap, matching _setup_geometry)
//   domain      — { x_min, x_max, y_min, y_max }
//   boundaries  — [left,right,top,bottom] of 'open'|'gnd'
//   hFine, hCoarse — characteristic mesh sizes near conductors / in the far field
//   meshConductorInterior — also mesh conductor cells (for volume solves); default false
//
// Returns: { ...meshObject, condRect, epsMap, lossMap, meshedDomain }
export function buildGmshMeshFromGeometry(G, opts) {
    const conductors = opts.conductors || [];
    const dielectrics = opts.dielectrics || [];
    const domain = opts.domain;
    const boundaries = opts.boundaries;
    const meshInterior = !!opts.meshConductorInterior;

    const diag = Math.hypot(domain.x_max - domain.x_min, domain.y_max - domain.y_min);
    const tol = diag * REL_TOL;
    const hFine = opts.hFine;
    const hCoarse = opts.hCoarse;

    let { X0, X1, Y0, Y1, wallPEC } = _clipDomain(domain, conductors, boundaries, tol);

    // Symmetry: mesh only the right half [0, X1]. The plane x=0 is a symmetry
    // boundary (PMC for even/single, PEC for odd) applied per-mode via the freedom
    // map's abc.left — so it is NOT a fixed PEC wall here.
    const symmetry = !!opts.symmetry;
    if (symmetry) { X0 = 0; wallPEC = { ...wallPEC, left: false }; }

    // Conductor rects clipped to the meshed domain (those fully outside are dropped;
    // absorbed wall slabs collapse to zero extent and are dropped).
    const condRects = [];
    const condRoles = [];
    for (const c of conductors) {
        const r = _rectOf(c);
        const xmin = Math.max(r.xmin, X0), xmax = Math.min(r.xmax, X1);
        const ymin = Math.max(r.ymin, Y0), ymax = Math.min(r.ymax, Y1);
        if (xmax - xmin <= tol || ymax - ymin <= tol) continue;  // absorbed/outside
        condRects.push({ xmin, xmax, ymin, ymax });
        condRoles.push({ is_signal: !!c.is_signal, polarity: c.polarity || 0, plating: c.plating || null });
    }

    // ---- Arrangement coordinates ----
    const xsAll = [X0, X1];
    const ysAll = [Y0, Y1];
    const within = (v, lo, hi) => v > lo + tol && v < hi - tol;
    for (const c of condRects) {
        if (within(c.xmin, X0, X1)) xsAll.push(c.xmin);
        if (within(c.xmax, X0, X1)) xsAll.push(c.xmax);
        if (within(c.ymin, Y0, Y1)) ysAll.push(c.ymin);
        if (within(c.ymax, Y0, Y1)) ysAll.push(c.ymax);
    }
    for (const d of dielectrics) {
        const r = _rectOf(d);
        if (within(r.xmin, X0, X1)) xsAll.push(r.xmin);
        if (within(r.xmax, X0, X1)) xsAll.push(r.xmax);
        if (within(r.ymin, Y0, Y1)) ysAll.push(r.ymin);
        if (within(r.ymax, Y0, Y1)) ysAll.push(r.ymax);
    }
    const xs = _uniqSorted(xsAll, tol);
    const ys = _uniqSorted(ysAll, tol);
    const nx = xs.length, ny = ys.length;

    // Snap conductor rect coords to the arrangement grid (so membership tests
    // and the FEM freedom map align exactly with mesh lines).
    const snap = (v, arr) => {
        let best = arr[0], bd = Math.abs(v - arr[0]);
        for (const a of arr) { const d = Math.abs(v - a); if (d < bd) { bd = d; best = a; } }
        return best;
    };
    for (const c of condRects) {
        c.xmin = snap(c.xmin, xs); c.xmax = snap(c.xmax, xs);
        c.ymin = snap(c.ymin, ys); c.ymax = snap(c.ymax, ys);
    }

    // Point on a conductor boundary? (for fine sizing)
    const onConductor = (x, y) => {
        for (const c of condRects) {
            if (x >= c.xmin - tol && x <= c.xmax + tol && y >= c.ymin - tol && y <= c.ymax + tol) {
                const onV = (Math.abs(x - c.xmin) < tol || Math.abs(x - c.xmax) < tol);
                const onH = (Math.abs(y - c.ymin) < tol || Math.abs(y - c.ymax) < tol);
                if (onV || onH) return true;
            }
        }
        return false;
    };
    // Cell centroid inside a conductor?
    const cellInConductor = (xc, yc) => {
        for (const c of condRects) {
            if (xc > c.xmin - tol && xc < c.xmax + tol && yc > c.ymin - tol && yc < c.ymax + tol) return true;
        }
        return false;
    };

    // ---- gmsh model ----
    const stack = G.stackSave();
    const ierr = G.stackAlloc(4);
    const check = (label) => { const e = G.getValue(ierr, 'i32'); if (e) throw new Error(`gmsh: ${label} failed (err=${e})`); };

    G._gmshClear(ierr); check('Clear');
    G._gmshModelAdd(_cstr(G, 'mesh'), ierr); check('ModelAdd');

    // Feature-size at a point: conductor corners get hFine; the size then grades
    // with DISTANCE to the nearest conductor (so the smooth field region around and
    // BETWEEN conductors — the inter-trace gap that sets the mutual R — is resolved,
    // not just the corners); and is capped by the local dielectric layer thickness.
    const sizeFactor = opts.sizeFactor ?? 0.6;
    const gradeRate = opts.gradeRate ?? 0.35;   // size growth per unit distance from a conductor
    const dielRects = dielectrics.map(d => _rectOf(d));
    const distToRect = (x, y, r) => {
        const dx = Math.max(r.xmin - x, 0, x - r.xmax);
        const dy = Math.max(r.ymin - y, 0, y - r.ymax);
        return Math.hypot(dx, dy);
    };
    // Dielectric band y-levels (shared by element sizing below and band construction
    // later). Interfaces closer than ~hFine are merged: a sub-element-thin strip (e.g. a
    // solder-mask layer whose top nearly coincides with the trace top when SM thickness
    // ≈ trace thickness) would otherwise be unmeshable. Real layers (≳hFine) survive.
    const bandYset = [Y0, Y1];
    for (const d of dielectrics) {
        const r = _rectOf(d);
        if (r.ymin > Y0 + tol && r.ymin < Y1 - tol) bandYset.push(r.ymin);
        if (r.ymax > Y0 + tol && r.ymax < Y1 - tol) bandYset.push(r.ymax);
    }
    const bandYs = _uniqSorted(bandYset, Math.max(tol, 0.25 * hFine));
    const localSize = (x, y) => {
        if (onConductor(x, y)) return hFine;
        let h = hCoarse;
        // grade from nearest conductor
        let dmin = Infinity;
        for (const c of condRects) dmin = Math.min(dmin, distToRect(x, y, c));
        if (isFinite(dmin)) h = Math.min(h, hFine + gradeRate * dmin);
        // cap by containing dielectric layer thickness
        for (const r of dielRects) {
            if (x >= r.xmin - tol && x <= r.xmax + tol && y >= r.ymin - tol && y <= r.ymax + tol) {
                h = Math.min(h, sizeFactor * Math.min(r.xmax - r.xmin, r.ymax - r.ymin));
            }
        }
        // cap by the local band spacing — resolves thin SUB-bands of a large region
        // (e.g. air slivers left by full-width solder-mask interface splits), which the
        // per-rect cap above misses because their containing rect (air) is huge.
        let bi = 0;
        while (bi < bandYs.length - 1 && bandYs[bi + 1] < y - tol) bi++;
        if (bi < bandYs.length - 1) h = Math.min(h, sizeFactor * (bandYs[bi + 1] - bandYs[bi]));
        return Math.max(hFine, Math.min(hCoarse, h));
    };
    // Coordinate-keyed point + line caches (shared across bands → conforming interfaces).
    const ptCache = new Map();
    const keyOf = (x, y) => `${Math.round((x - X0) / (diag * 1e-9))},${Math.round((y - Y0) / (diag * 1e-9))}`;
    const getPt = (x, y) => {
        const k = keyOf(x, y);
        let t = ptCache.get(k);
        if (t === undefined) { t = G._gmshModelGeoAddPoint(x, y, 0, localSize(x, y), -1, ierr); check('pt'); ptCache.set(k, t); }
        return t;
    };
    const lineCache = new Map();
    const getLine = (x1, y1, x2, y2) => {
        const ka = keyOf(x1, y1), kb = keyOf(x2, y2);
        const fwd = ka + '|' + kb, rev = kb + '|' + ka;
        if (lineCache.has(fwd)) return lineCache.get(fwd);
        if (lineCache.has(rev)) return -lineCache.get(rev);
        const l = G._gmshModelGeoAddLine(getPt(x1, y1), getPt(x2, y2), -1, ierr); check('line');
        lineCache.set(fwd, l); return l;
    };
    // CCW polygon (array of [x,y]) → curve loop tag.
    const makeLoop = (poly, label) => {
        const tags = [];
        for (let i = 0; i < poly.length; i++) {
            const a = poly[i], b = poly[(i + 1) % poly.length];
            tags.push(getLine(a[0], a[1], b[0], b[1]));
        }
        const buf = G.stackAlloc(tags.length * 4);
        tags.forEach((v, k) => G.setValue(buf + k * 4, v, 'i32'));
        const loop = G._gmshModelGeoAddCurveLoop(buf, tags.length, -1, 0, ierr); check(label);
        return loop;
    };

    // ---- Build dielectric bands (split only at dielectric interface y-values) ----
    // Conductors become HOLES (interior) or boundary notches (wall-touching), so
    // conductor edges never create thin full-width strips → no sliver triangles.
    // bandYs is computed above (shared with localSize).

    for (let bi = 0; bi < bandYs.length - 1; bi++) {
        const ya = bandYs[bi], yb = bandYs[bi + 1];
        // Conductors clipped to this band (always carved from the air region so the
        // conductor boundary is a conforming mesh edge; with meshInterior they are
        // also meshed as their own surfaces at finite σ for the MQS loss solve).
        const inBand = [];
        for (const c of condRects) {
            const xmin = Math.max(c.xmin, X0), xmax = Math.min(c.xmax, X1);
            const cymin = Math.max(c.ymin, ya), cymax = Math.min(c.ymax, yb);
            if (xmax - xmin > tol && cymax - cymin > tol) inBand.push({ xmin, xmax, ymin: cymin, ymax: cymax });
        }
        const holes = [], touchers = [];
        for (const c of inBand) {
            const touch = Math.abs(c.xmin - X0) < tol || Math.abs(c.xmax - X1) < tol ||
                          Math.abs(c.ymin - ya) < tol || Math.abs(c.ymax - yb) < tol;
            (touch ? touchers : holes).push(c);
        }
        // Conformity: split the bottom interface (ya) at the x-edges of conductors
        // touching ya, and the top interface (yb) at conductors touching yb — so each
        // shared interface line is split identically in both adjacent bands.
        const forcedXbottom = [], forcedXtop = [];
        for (const c of condRects) {
            if (Math.abs(c.ymin - ya) < tol || Math.abs(c.ymax - ya) < tol) forcedXbottom.push(c.xmin, c.xmax);
            if (Math.abs(c.ymin - yb) < tol || Math.abs(c.ymax - yb) < tol) forcedXtop.push(c.xmin, c.xmax);
        }
        // Outer air boundary = band rect minus wall-touching conductors. Traced from
        // a local cell mask over the touchers' edge coordinates (robust for notches).
        const outers = traceAirBoundary(X0, X1, ya, yb, touchers, tol, forcedXbottom, forcedXtop);
        for (const outer of outers) {
            const loopTags = [makeLoop(outer, 'bandOuter')];
            // Interior holes whose centroid lies inside this outer ring become holes.
            for (const h of holes) {
                const hx = (h.xmin + h.xmax) / 2, hy = (h.ymin + h.ymax) / 2;
                if (pointInPoly(hx, hy, outer)) {
                    loopTags.push(makeLoop([[h.xmin, h.ymin], [h.xmax, h.ymin], [h.xmax, h.ymax], [h.xmin, h.ymax]], 'hole'));
                }
            }
            const wbuf = G.stackAlloc(loopTags.length * 4);
            loopTags.forEach((v, k) => G.setValue(wbuf + k * 4, v, 'i32'));
            G._gmshModelGeoAddPlaneSurface(wbuf, loopTags.length, -1, ierr); check('bandSurf');
        }
        // Conductor interior surfaces (for the MQS volume-eddy loss solve). They
        // share boundary lines with the air region (cached) → conforming.
        if (meshInterior) {
            for (const c of inBand) {
                const loop = makeLoop([[c.xmin, c.ymin], [c.xmax, c.ymin], [c.xmax, c.ymax], [c.xmin, c.ymax]], 'condSurf');
                const wbuf = G.stackAlloc(4);
                G.setValue(wbuf, loop, 'i32');
                G._gmshModelGeoAddPlaneSurface(wbuf, 1, -1, ierr); check('condSurf');
            }
        }
    }

    // ---- Mesh options ----
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.Algorithm'), 5, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.Smoothing'), 10, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.RandomSeed'), 1, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.CharacteristicLengthMin'), Math.min(hFine, hCoarse) * 0.25, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.CharacteristicLengthMax'), hCoarse, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.CharacteristicLengthExtendFromBoundary'), 1, ierr);

    G._gmshModelGeoSynchronize(ierr); check('Synchronize');
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
    for (let i = 0; i < nNodes; i++) { nodes[2 * i] = coords[3 * i]; nodes[2 * i + 1] = coords[3 * i + 1]; }

    // ---- Extract triangles ----
    const etPtrPtr = G.stackAlloc(4), etNPtr = G.stackAlloc(4);
    const enPtrPtr = G.stackAlloc(4), enNPtr = G.stackAlloc(4);
    G._gmshModelMeshGetElementsByType(2, etPtrPtr, etNPtr, enPtrPtr, enNPtr, -1, 0, 1, ierr);
    check('GetElementsByType');
    const nTris = G.getValue(etNPtr, 'i32');
    const enPtr = G.getValue(enPtrPtr, 'i32');
    G._gmshFree(G.getValue(etPtrPtr, 'i32'));
    const tris = new Int32Array(3 * nTris);
    for (let t = 0; t < nTris; t++)
        for (let k = 0; k < 3; k++) tris[3 * t + k] = tagToIdx.get(G.getValue(enPtr + (3 * t + k) * 4, 'i32'));
    G._gmshFree(enPtr);

    G.stackRestore(stack);

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
    const meshObj0 = { nodes, tris, nTris };
    const { epsMap, lossMap } = tagMaterials(meshObj0, dielectrics, tol);

    // ---- condRect for the FEM freedom map ----
    const condRect = {
        rects: condRects,
        rectRoles: condRoles,
        xmin_domain: X0, xmax_domain: X1, ymin_domain: Y0, ymax_domain: Y1,
        wallPEC, symmetry: symmetry ? 2 : 1,
        // legacy single-rect bbox (some callers read these)
        xmin: condRects.length ? Math.min(...condRects.map(c => c.xmin)) : X0,
        xmax: condRects.length ? Math.max(...condRects.map(c => c.xmax)) : X1,
        ymin: condRects.length ? Math.min(...condRects.map(c => c.ymin)) : Y0,
        ymax: condRects.length ? Math.max(...condRects.map(c => c.ymax)) : Y1,
    };

    // constraint metadata (used by refineTriMesh + checkMeshQuality): each
    // conductor edge spans only the conductor's own x-range (union if shared), plus
    // dielectric interfaces span the full width.
    const constraintYRanges = {};
    const addYR = (y, lo, hi) => {
        const e = constraintYRanges[y];
        constraintYRanges[y] = e ? [Math.min(e[0], lo), Math.max(e[1], hi)] : [lo, hi];
    };
    for (const c of condRects) { addYR(c.ymin, c.xmin, c.xmax); addYR(c.ymax, c.xmin, c.xmax); }
    for (const d of dielectrics) {
        const r = _rectOf(d);
        if (r.ymin > Y0 + tol && r.ymin < Y1 - tol) addYR(r.ymin, X0, X1);
        if (r.ymax > Y0 + tol && r.ymax < Y1 - tol) addYR(r.ymax, X0, X1);
    }
    const constraintXs = _uniqSorted(
        [...(symmetry ? [X0] : []), ...condRects.flatMap(c => [c.xmin, c.xmax])], tol);

    return {
        nodes, tris, edges, triEdges, triSigns, nNodes, nTris, nEdges,
        constraintYRanges, constraintXs,
        condRect, epsMap, lossMap,
        meshedDomain: { X0, X1, Y0, Y1, wallPEC },
    };
}

export { generateGmshMesh };
