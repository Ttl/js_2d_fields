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

const REL_TOL = 1e-9;

// Per-triangle ε / loss maps by centroid (last-match-wins over the dielectric list,
// matching _setup_geometry). Reusable: call again after refineTriMesh.
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
    if (!roles || !roles.some(r => r.is_signal)) errors.push('No signal conductor found in the geometry.');
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
    let changed = true;
    while (changed) {
        changed = false;
        for (const c of conductors) {
            if (c.is_signal) continue;
            const r = _rectOf(c);
            const touchesL = r.xmin <= X0 + tol, touchesR = r.xmax >= X1 - tol;
            const touchesB = r.ymin <= Y0 + tol, touchesT = r.ymax >= Y1 - tol;
            if (touchesB && touchesL && touchesR && r.ymax < Y1 - tol && r.ymax > Y0 + tol) { Y0 = r.ymax; wallPEC.bottom = true; changed = true; continue; }
            if (touchesT && touchesL && touchesR && r.ymin > Y0 + tol && r.ymin < Y1 - tol) { Y1 = r.ymin; wallPEC.top = true; changed = true; continue; }
            if (touchesL && touchesB && touchesT && r.xmax < X1 - tol && r.xmax > X0 + tol) { X0 = r.xmax; wallPEC.left = true; changed = true; continue; }
            if (touchesR && touchesB && touchesT && r.xmin > X0 + tol && r.xmin < X1 - tol) { X1 = r.xmin; wallPEC.right = true; changed = true; continue; }
        }
    }
    return { X0, X1, Y0, Y1, wallPEC };
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

    let { X0, X1, Y0, Y1, wallPEC } = _clipDomain(domain, conductors, boundaries, tol);
    const symmetry = !!opts.symmetry;
    if (symmetry) { X0 = 0; wallPEC = { ...wallPEC, left: false }; }

    // Conductor rects clipped to the meshed domain (absorbed/outside ones dropped).
    const condRects = [], condRoles = [];
    for (const c of conductors) {
        const r = _rectOf(c);
        const xmin = Math.max(r.xmin, X0), xmax = Math.min(r.xmax, X1);
        const ymin = Math.max(r.ymin, Y0), ymax = Math.min(r.ymax, Y1);
        if (xmax - xmin <= tol || ymax - ymin <= tol) continue;
        condRects.push({ xmin, xmax, ymin, ymax });
        condRoles.push({ is_signal: !!c.is_signal, polarity: c.polarity || 0, plating: c.plating || null });
    }

    const stack = G.stackSave();
    const ierr = G.stackAlloc(4);
    const check = (label) => { const e = G.getValue(ierr, 'i32'); if (e) throw new Error(`occ: ${label} failed (err=${e})`); };

    G._gmshClear(ierr); check('Clear');
    G._gmshModelAdd(_cstr(G, 'occmesh'), ierr); check('ModelAdd');

    const addRect = (r) => {
        const t = G._gmshModelOccAddRectangle(r.xmin, r.ymin, 0, r.xmax - r.xmin, r.ymax - r.ymin, -1, 0, ierr);
        check('AddRectangle'); return t;
    };
    const clipToDomain = (r) => ({
        xmin: Math.max(r.xmin, X0), xmax: Math.min(r.xmax, X1),
        ymin: Math.max(r.ymin, Y0), ymax: Math.min(r.ymax, Y1),
    });

    // Background "air" rectangle (object) + all dielectric & conductor rects (tools).
    const domTag = addRect({ xmin: X0, xmax: X1, ymin: Y0, ymax: Y1 });
    const toolTags = [];
    for (const d of dielectrics) {
        const r = clipToDomain(_rectOf(d));
        if (r.xmax - r.xmin > tol && r.ymax - r.ymin > tol) toolTags.push(addRect(r));
    }
    for (const c of condRects) toolTags.push(addRect(c));

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
    const condCurves = new Set();
    {
        const fe = G.stackAlloc(4), feN = G.stackAlloc(4);
        G._gmshModelGetEntities(fe, feN, 2, ierr); check('GetEntities(2)');
        const faces = _readIntArray(G, fe, feN);   // (dim,tag) pairs
        const bb = [G.stackAlloc(8), G.stackAlloc(8), G.stackAlloc(8), G.stackAlloc(8), G.stackAlloc(8), G.stackAlloc(8)];
        // A conductor sub-face's bounding box is CONTAINED in its conductor rect.
        // Testing only the bbox CENTER misclassifies a face-with-hole: the air
        // face surrounding a centered conductor (e.g. symmetric stripline) has
        // its bbox center inside the conductor rect, so the whole domain outline
        // would be collected as "conductor" curves and refined at sizeMin.
        // OCC enlarges bounding boxes by Precision::Confusion (1e-7 in model
        // units, verified empirically) — the containment pad must absorb that,
        // capped per rect so a very thin conductor can't false-positive.
        const OCC_BBOX_PAD = 2e-7;
        const inCond = (x0, y0, x1, y1) => condRects.some(c => {
            const pad = Math.max(tol, Math.min(OCC_BBOX_PAD,
                0.4 * Math.min(c.xmax - c.xmin, c.ymax - c.ymin)));
            return x0 > c.xmin - pad && x1 < c.xmax + pad &&
                   y0 > c.ymin - pad && y1 < c.ymax + pad;
        });
        for (let i = 0; i < faces.length; i += 2) {
            const ftag = faces[i + 1];
            G._gmshModelGetBoundingBox(2, ftag, bb[0], bb[1], bb[2], bb[3], bb[4], bb[5], ierr); check('GetBoundingBox');
            const bx0 = G.getValue(bb[0], 'double'), by0 = G.getValue(bb[1], 'double');
            const bx1 = G.getValue(bb[3], 'double'), by1 = G.getValue(bb[4], 'double');
            if (!inCond(bx0, by0, bx1, by1)) continue;
            const fbuf = G.stackAlloc(8); G.setValue(fbuf, 2, 'i32'); G.setValue(fbuf + 4, ftag, 'i32');
            const cb = G.stackAlloc(4), cbN = G.stackAlloc(4);
            G._gmshModelGetBoundary(fbuf, 2, cb, cbN, 0, 0, 0, ierr); check('GetBoundary');
            const curves = _readIntArray(G, cb, cbN);
            for (let k = 1; k < curves.length; k += 2) condCurves.add(Math.abs(curves[k]));
        }
    }

    // Conductor-surface element size — finer than global hFine (see field setup below).
    const surfScale = opts.occSurfScale ?? 0.35;
    const sizeMin = hFine * surfScale;
    if (globalThis.__OCC_DEBUG__) console.error('[occ] condRects', condRects.length, 'condCurves', condCurves.size, 'hFine', hFine, 'sizeMin', sizeMin, 'hCoarse', hCoarse);
    if (condCurves.size > 0) {
        const fd = G._gmshModelMeshFieldAdd(_cstr(G, 'Distance'), -1, ierr); check('FieldAdd Distance');
        const cl = [...condCurves];
        const clBuf = G.stackAlloc(cl.length * 8);
        cl.forEach((v, i) => G.setValue(clBuf + i * 8, v, 'double'));
        G._gmshModelMeshFieldSetNumbers(fd, _cstr(G, 'CurvesList'), clBuf, cl.length, ierr); check('CurvesList');
        G._gmshModelMeshFieldSetNumber(fd, _cstr(G, 'Sampling'), 200, ierr);
        // sizeMin (conductor-surface element size) is finer than the global hFine so the
        // current distribution (skin/proximity loss) and a narrow inter-trace gap (only
        // ~1·hFine wide) are resolved WITHOUT globally refining the far field. Adaptive ZZ
        // passes sharpen it further.
        const ft = G._gmshModelMeshFieldAdd(_cstr(G, 'Threshold'), -1, ierr); check('FieldAdd Threshold');
        G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'InField'), fd, ierr);
        G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'SizeMin'), sizeMin, ierr);
        G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'SizeMax'), hCoarse, ierr);
        G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'DistMin'), sizeMin, ierr);
        G._gmshModelMeshFieldSetNumber(ft, _cstr(G, 'DistMax'), hFine + (hCoarse - hFine) / gradeRate, ierr);
        G._gmshModelMeshFieldSetAsBackgroundMesh(ft, ierr); check('SetBackground');
    }

    // Field-only sizing (don't let point/curvature/boundary sizing fight it).
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeExtendFromBoundary'), 0, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeFromPoints'), 0, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeFromCurvature'), 0, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeMin'), Math.min(hFine * 0.5, sizeMin * 0.5), ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.MeshSizeMax'), hCoarse, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.Algorithm'), 5, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'Mesh.Smoothing'), 10, ierr);
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
        wallPEC, symmetry: symmetry ? 2 : 1,
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
    for (const c of condRects) {
        addYR(c.ymin, c.xmin, c.xmax); addYR(c.ymax, c.xmin, c.xmax);
        addXR(c.xmin, c.ymin, c.ymax); addXR(c.xmax, c.ymin, c.ymax);
    }
    if (symmetry) addXR(X0, Y0, Y1);
    const constraintXs = Object.keys(constraintXRanges).map(Number);

    return {
        nodes, tris, edges, triEdges, triSigns, nNodes, nTris, nEdges,
        constraintYRanges, constraintXRanges, constraintXs,
        condRect, epsMap, lossMap,
        meshedDomain: { X0, X1, Y0, Y1, wallPEC },
    };
}
