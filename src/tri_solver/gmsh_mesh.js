// gmsh WASM mesh generator for 2D transmission line cross-sections.
// Builds geometry from domain bounds, dielectric interface, and PEC conductor rectangles.
// Supports symmetry (half-domain with PMC at x=0).
//
// Usage:
//   const G = await initGmsh();
//   const mesh = generateGmshMesh(G, { domainW, domainH, conductors, interfaceY, hFine, hCoarse });

// ---- Internal helpers ----

function _cstr(G, str) {
    const len = G.lengthBytesUTF8(str) + 1;
    const ptr = G.stackAlloc(len);
    G.stringToUTF8(str, ptr, len);
    return ptr;
}

function _readDoubleArray(G, ptrPtr, nPtr) {
    const ptr = G.getValue(ptrPtr, 'i32');
    const n   = G.getValue(nPtr,   'i32');
    const arr = new Float64Array(n);
    for (let i = 0; i < n; i++) arr[i] = G.getValue(ptr + i * 8, 'double');
    G._gmshFree(ptr);
    return arr;
}

function _readIntArray(G, ptrPtr, nPtr) {
    const ptr = G.getValue(ptrPtr, 'i32');
    const n   = G.getValue(nPtr,   'i32');
    const arr = new Int32Array(n);
    for (let i = 0; i < n; i++) arr[i] = G.getValue(ptr + i * 4, 'i32');
    G._gmshFree(ptr);
    return arr;
}

// ---- Public API ----

// Load and initialize the gmsh WASM module. Call once; pass the returned G
// to generateGmshMesh for every mesh.
export async function initGmsh() {
    const gmshUrl = new URL('./gmsh/build-wasm/gmsh.js', import.meta.url).href;
    const G = await (await import(gmshUrl)).default();
    const stack = G.stackSave();
    const ierr  = G.stackAlloc(4);
    G._gmshInitialize(0, 0, 1, 0, ierr);
    G._gmshOptionSetNumber(_cstr(G, 'General.Verbosity'), 0, ierr);
    G.stackRestore(stack);
    // Store a persistent ierr slot on the module object for convenience
    return G;
}

// Generate a solver-compatible triangular mesh for a 2D cross-section.
//
// Parameters:
//   domainW    — full domain width; domain spans x ∈ [-domainW/2, +domainW/2]
//                (or x ∈ [0, domainW/2] when symmetry is true)
//   domainH    — full domain height; y ∈ [0, domainH]
//   conductors — [{xmin, xmax, ymin, ymax}, ...], sorted by xmin (or will be sorted)
//   interfaceY — y-coordinate of the dielectric interface (e.g. substrate height).
//                If omitted, derived from conductors[0].ymin for backward compat.
//   hFine      — mesh size at conductor corners
//   hCoarse    — mesh size at domain corners
//   hMed       — (optional) mesh size at interface/domain-wall endpoints;
//                defaults to (hFine + hCoarse) / 2
//   symmetry   — (optional) if true, mesh only the right half [0, domainW/2]
//                with a symmetry plane at x=0
//   meshConductorInterior — (optional) if true, conductor interiors are meshed
//                as their own surfaces (for volume eddy-current solves) instead
//                of being left as holes. The FEM freedom maps already skip
//                conductor-interior triangles, so eigensolves are unaffected.
//
// Returns the standard solver mesh object:
//   { nodes, tris, edges, triEdges, triSigns, nNodes, nTris, nEdges,
//     constraintYRanges, constraintXs }
export function generateGmshMesh(G, { domainW, domainH, conductors, interfaceY, hFine, hCoarse, hMed, symmetry, meshConductorInterior }) {
    if (!hMed) hMed = (hFine + hCoarse) / 2;

    // Sort conductors left to right
    conductors = [...conductors].sort((a, b) => a.xmin - b.xmin);
    const N = conductors.length;
    const h_sub = interfaceY !== undefined ? interfaceY : conductors[0].ymin;

    const stack = G.stackSave();
    const ierr  = G.stackAlloc(4);

    function check(label) {
        const e = G.getValue(ierr, 'i32');
        if (e) throw new Error(`gmsh: ${label} failed (err=${e})`);
    }
    function cstr(s) { return _cstr(G, s); }

    // Reuse module: clear previous model and add a fresh one
    G._gmshClear(ierr); check('Clear');
    G._gmshModelAdd(cstr('mesh'), ierr); check('ModelAdd');

    const xL = symmetry ? 0 : -domainW / 2, xR = domainW / 2;

    const TOL = 1e-10;
    // Detect conductors whose edges coincide with domain walls.
    // Wall-touching conductors become boundary notches instead of holes.
    const leftTouch = N > 0 && Math.abs(conductors[0].xmin - xL) < TOL;

    // ---- Points ----
    // Domain corners
    const pBL = G._gmshModelGeoAddPoint(xL, 0,      0, hCoarse, -1, ierr); check('pBL');
    const pBR = G._gmshModelGeoAddPoint(xR, 0,      0, hCoarse, -1, ierr); check('pBR');
    const pTR = G._gmshModelGeoAddPoint(xR, domainH, 0, hCoarse, -1, ierr); check('pTR');
    const pTL = G._gmshModelGeoAddPoint(xL, domainH, 0, hCoarse, -1, ierr); check('pTL');

    // Interface endpoints on domain walls
    const pIL = G._gmshModelGeoAddPoint(xL, h_sub, 0, leftTouch ? hFine : hMed, -1, ierr); check('pIL');
    const pIR = G._gmshModelGeoAddPoint(xR, h_sub, 0, hMed, -1, ierr); check('pIR');

    // Per-conductor points on h_sub and at top
    const pCbot = [];  // [i] = [BL tag, BR tag]
    const pCtop = [];  // [i] = [TL tag, TR tag]
    for (let i = 0; i < N; i++) {
        const c = conductors[i];
        let bl, tl;
        if (i === 0 && leftTouch) {
            bl = pIL;  // reuse interface-left point (same location)
            tl = G._gmshModelGeoAddPoint(xL, c.ymax, 0, hFine, -1, ierr); check(`c${i}TL`);
        } else {
            bl = G._gmshModelGeoAddPoint(c.xmin, h_sub,  0, hFine, -1, ierr); check(`c${i}BL`);
            tl = G._gmshModelGeoAddPoint(c.xmin, c.ymax, 0, hFine, -1, ierr); check(`c${i}TL`);
        }
        const br = G._gmshModelGeoAddPoint(c.xmax, h_sub,  0, hFine, -1, ierr); check(`c${i}BR`);
        const tr = G._gmshModelGeoAddPoint(c.xmax, c.ymax, 0, hFine, -1, ierr); check(`c${i}TR`);
        pCbot.push([bl, br]);
        pCtop.push([tl, tr]);
    }

    // ---- Lines ----
    // Domain perimeter (split at h_sub and at wall-touching conductor tops)
    const lBottom      = G._gmshModelGeoAddLine(pBL, pBR, -1, ierr); check('lBottom');
    const lRightLower  = G._gmshModelGeoAddLine(pBR, pIR, -1, ierr); check('lRightLower');
    const lRightUpper  = G._gmshModelGeoAddLine(pIR, pTR, -1, ierr); check('lRightUpper');
    const lTop         = G._gmshModelGeoAddLine(pTR, pTL, -1, ierr); check('lTop');

    // Left wall: when leftTouch, split into upper-upper (pTL→c0TL) + lower (pIL→pBL).
    // The segment c0TL→pIL is the conductor's left face on the wall — not a separate line,
    // since the air boundary wraps around the conductor notch instead.
    let lLeftUpper = null;       // pTL → pIL (normal case)
    let lLeftUpperUpper = null;  // pTL → c0TL (leftTouch case)
    if (leftTouch) {
        lLeftUpperUpper = G._gmshModelGeoAddLine(pTL, pCtop[0][0], -1, ierr); check('lLeftUpperUpper');
    } else {
        lLeftUpper = G._gmshModelGeoAddLine(pTL, pIL, -1, ierr); check('lLeftUpper');
    }
    const lLeftLower = G._gmshModelGeoAddLine(pIL, pBL, -1, ierr); check('lLeftLower');

    // Interface segments: gaps between/around conductors + conductor-bottom segments.
    // A gap is null when the conductor reuses a wall point (zero-length).
    const lIface = new Array(N + 1);
    const lCbot  = new Array(N);

    const ifacePts = [pIL, ...pCbot.flatMap(([bl, br]) => [bl, br]), pIR];
    for (let i = 0; i < N; i++) {
        const fromPt = ifacePts[2*i], toPt = ifacePts[2*i+1];
        if (fromPt === toPt) {
            lIface[i] = null;  // conductor touches wall — zero-length gap
        } else {
            lIface[i] = G._gmshModelGeoAddLine(fromPt, toPt, -1, ierr); check(`lIface${i}`);
        }
        lCbot[i] = G._gmshModelGeoAddLine(ifacePts[2*i+1], ifacePts[2*i+2], -1, ierr); check(`lCbot${i}`);
    }
    if (ifacePts[2*N] === pIR) {
        lIface[N] = null;
    } else {
        lIface[N] = G._gmshModelGeoAddLine(ifacePts[2*N], pIR, -1, ierr); check(`lIfaceN`);
    }

    // Per-conductor walls and top (no left line for wall-touching conductor)
    const lCLeft  = new Array(N);
    const lCTop   = new Array(N);
    const lCRight = new Array(N);
    for (let i = 0; i < N; i++) {
        const [bl, br] = pCbot[i];
        const [tl, tr] = pCtop[i];
        if (i === 0 && leftTouch) {
            lCLeft[i] = null;  // left face is on domain wall
        } else {
            lCLeft[i] = G._gmshModelGeoAddLine(bl, tl, -1, ierr); check(`lCLeft${i}`);
        }
        lCTop[i]   = G._gmshModelGeoAddLine(tl, tr, -1, ierr); check(`lCTop${i}`);
        lCRight[i] = G._gmshModelGeoAddLine(tr, br, -1, ierr); check(`lCRight${i}`);
    }

    // ---- Curve loops ----
    function makeCurveLoop(tags, label) {
        const buf = G.stackAlloc(tags.length * 4);
        tags.forEach((v, k) => G.setValue(buf + k*4, v, 'i32'));
        const loop = G._gmshModelGeoAddCurveLoop(buf, tags.length, -1, 0, ierr); check(label);
        return loop;
    }

    // Substrate loop (CCW): bottom → right-lower → interface(R→L) → left-lower
    {
        const subTags = [lBottom, lRightLower];
        for (let i = N; i >= 0; i--) {
            if (lIface[i] !== null) subTags.push(-lIface[i]);
            if (i > 0) subTags.push(-lCbot[i-1]);
        }
        subTags.push(lLeftLower);
        const loop = makeCurveLoop(subTags, 'subLoop');
        const wbuf = G.stackAlloc(4);
        G.setValue(wbuf, loop, 'i32');
        G._gmshModelGeoAddPlaneSurface(wbuf, 1, -1, ierr); check('subSurf');
    }

    // Air outer loop (CCW).
    // When leftTouch, conductor[0] is a notch in the boundary rather than a hole:
    //   c0BR → [interface gaps] → pIR → pTR → pTL → c0TL → c0TR → c0BR
    // Otherwise: pIL → [interface L→R] → pIR → pTR → pTL → pIL
    const airOutTags = [];
    if (leftTouch) {
        // Interface segments after conductor 0
        for (let i = 1; i < N; i++) {
            if (lIface[i] !== null) airOutTags.push(lIface[i]);
            airOutTags.push(lCbot[i]);
        }
        if (lIface[N] !== null) airOutTags.push(lIface[N]);
        airOutTags.push(lRightUpper, lTop, lLeftUpperUpper, lCTop[0], lCRight[0]);
    } else {
        for (let i = 0; i < N; i++) {
            if (lIface[i] !== null) airOutTags.push(lIface[i]);
            airOutTags.push(lCbot[i]);
        }
        if (lIface[N] !== null) airOutTags.push(lIface[N]);
        airOutTags.push(lRightUpper, lTop, lLeftUpper);
    }
    const airOutLoop = makeCurveLoop(airOutTags, 'airOutLoop');

    // Conductor hole loops (skip wall-touching conductors — they're notches, not holes)
    const condLoops = [];
    for (let i = 0; i < N; i++) {
        if (i === 0 && leftTouch) continue;
        const tags = [lCLeft[i], lCTop[i], lCRight[i], -lCbot[i]];
        const loop = makeCurveLoop(tags, `condLoop${i}`);
        condLoops.push(loop);
    }

    // Air surface: outer + holes
    const airWireTags = [airOutLoop, ...condLoops];
    const airWireBuf = G.stackAlloc(airWireTags.length * 4);
    airWireTags.forEach((v, k) => G.setValue(airWireBuf + k*4, v, 'i32'));
    G._gmshModelGeoAddPlaneSurface(airWireBuf, airWireTags.length, -1, ierr); check('airSurf');

    // Conductor interior surfaces (for volume eddy-current solves).
    // The element export collects type-2 elements from all surfaces, so these
    // triangles appear in the output mesh automatically.
    if (meshConductorInterior) {
        for (let i = 0; i < N; i++) {
            let loopTags;
            if (i === 0 && leftTouch) {
                // Notch conductor: its left face lies on the domain wall and has no
                // curve yet — add one (used only by this surface).
                const [bl] = pCbot[0];
                const [tl] = pCtop[0];
                const lWall = G._gmshModelGeoAddLine(tl, bl, -1, ierr); check('lCWall0');
                loopTags = [lWall, lCbot[0], -lCRight[0], -lCTop[0]];
            } else {
                loopTags = [lCLeft[i], lCTop[i], lCRight[i], -lCbot[i]];
            }
            const loop = makeCurveLoop(loopTags, `condSurfLoop${i}`);
            const wbuf = G.stackAlloc(4);
            G.setValue(wbuf, loop, 'i32');
            G._gmshModelGeoAddPlaneSurface(wbuf, 1, -1, ierr); check(`condSurf${i}`);
        }
    }

    // ---- Mesh options ----
    G._gmshOptionSetNumber(cstr('Mesh.Algorithm'), 5, ierr);   // Delaunay
    G._gmshOptionSetNumber(cstr('Mesh.Smoothing'), 10, ierr);
    G._gmshOptionSetNumber(cstr('Mesh.RandomSeed'), 1, ierr);  // deterministic
    G._gmshOptionSetNumber(cstr('Mesh.CharacteristicLengthMin'), Math.min(hFine, hCoarse) * 0.5, ierr);
    G._gmshOptionSetNumber(cstr('Mesh.CharacteristicLengthMax'), hCoarse, ierr);
    G._gmshOptionSetNumber(cstr('Mesh.CharacteristicLengthExtendFromBoundary'), 1, ierr);

    G._gmshModelGeoSynchronize(ierr); check('Synchronize');
    G._gmshModelMeshGenerate(2, ierr); check('MeshGenerate');

    // ---- Extract nodes ----
    const nodePtrPtr  = G.stackAlloc(4), nodeNPtr  = G.stackAlloc(4);
    const coordPtrPtr = G.stackAlloc(4), coordNPtr = G.stackAlloc(4);
    const paramPtrPtr = G.stackAlloc(4), paramNPtr = G.stackAlloc(4);
    G._gmshModelMeshGetNodes(
        nodePtrPtr, nodeNPtr, coordPtrPtr, coordNPtr, paramPtrPtr, paramNPtr,
        -1, -1, 0, 0, ierr);
    check('GetNodes');

    const gmshTags = _readIntArray(G, nodePtrPtr, nodeNPtr);
    const coords   = _readDoubleArray(G, coordPtrPtr, coordNPtr);
    const paramPtr = G.getValue(paramPtrPtr, 'i32');
    if (paramPtr) G._gmshFree(paramPtr);

    const nNodes = gmshTags.length;
    const tagToIdx = new Map();
    for (let i = 0; i < nNodes; i++) tagToIdx.set(gmshTags[i], i);

    const nodes = new Float64Array(2 * nNodes);
    for (let i = 0; i < nNodes; i++) {
        nodes[2*i]   = coords[3*i];
        nodes[2*i+1] = coords[3*i+1];
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
    for (let t = 0; t < nTris; t++) {
        for (let k = 0; k < 3; k++) {
            tris[3*t+k] = tagToIdx.get(G.getValue(enPtr + (3*t+k)*4, 'i32'));
        }
    }
    G._gmshFree(enPtr);

    G.stackRestore(stack);

    // ---- Enforce CCW winding ----
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const ax = nodes[2*v0], ay = nodes[2*v0+1];
        const bx = nodes[2*v1], by = nodes[2*v1+1];
        const cx = nodes[2*v2], cy = nodes[2*v2+1];
        if ((bx-ax)*(cy-ay) - (cx-ax)*(by-ay) < 0) {
            tris[3*t+1] = v2; tris[3*t+2] = v1;
        }
    }

    // ---- Build edge list, triEdges, triSigns ----
    const edgeMap      = new Map();
    const edgeList     = [];
    const triEdgesArr  = new Int32Array(3 * nTris);
    const triSignsArr  = new Int8Array(3 * nTris);
    const localEdges   = [[0,1],[1,2],[2,0]];

    for (let t = 0; t < nTris; t++) {
        for (let le = 0; le < 3; le++) {
            const na = tris[3*t + localEdges[le][0]];
            const nb = tris[3*t + localEdges[le][1]];
            const n0 = na < nb ? na : nb;
            const n1 = na < nb ? nb : na;
            const key = n0 * (nNodes + 1) + n1;
            let eIdx = edgeMap.get(key);
            if (eIdx === undefined) {
                eIdx = edgeList.length >> 1;
                edgeList.push(n0, n1);
                edgeMap.set(key, eIdx);
            }
            triEdgesArr[3*t+le] = eIdx;
            triSignsArr[3*t+le] = (na === n0) ? 1 : -1;
        }
    }
    const nEdges = edgeList.length >> 1;
    const edges  = new Int32Array(edgeList);

    // ---- Constraint metadata for refineTriMesh ----
    const constraintYRanges = {};
    for (const c of conductors) {
        constraintYRanges[c.ymin] = [xL, xR];
        constraintYRanges[c.ymax] = [xL, xR];
    }
    const constraintXs = [...new Set(conductors.flatMap(c => [c.xmin, c.xmax]))].sort((a, b) => a - b);

    return { nodes, tris, edges, triEdges: triEdgesArr, triSigns: triSignsArr,
             nNodes, nTris, nEdges, constraintYRanges, constraintXs };
}
