// Mesh utilities for 2D FEM: element order assignment and quality checking

// --- Assign element polynomial orders for hp-refinement ---
// Elements near PEC corners get P3, others stay P2.
// corners: array of {x, y} corner points
// radius: elements within this distance from a corner get P3
export function assignElementOrders(mesh, corners, radius) {
    const { nodes, tris, nTris } = mesh;
    const elemOrder = new Uint8Array(nTris).fill(2);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const xc = (nodes[2*v0] + nodes[2*v1] + nodes[2*v2]) / 3;
        const yc = (nodes[2*v0+1] + nodes[2*v1+1] + nodes[2*v2+1]) / 3;
        for (const c of corners) {
            const d = Math.sqrt((xc - c.x) ** 2 + (yc - c.y) ** 2);
            if (d < radius) { elemOrder[t] = 3; break; }
        }
    }
    return elemOrder;
}

// --- Mesh quality check ---
// Validates mesh quality before FEM solve. Returns { ok, warnings, errors, metrics }.
// constraintYs/constraintXs: arrays of y/x values where material interfaces exist.
export function checkMeshQuality(mesh, constraintYs, constraintXs) {
    const { nodes, tris, edges, nTris, nEdges, nNodes } = mesh;
    const warnings = [], errors = [];

    // --- Triangle quality (Q = circumradius / (2 * inradius), ideal = 1) ---
    let maxQ = 0, sumQ = 0, badCount = 0, degenerateCount = 0, worstTri = -1;
    for (let t = 0; t < nTris; t++) {
        const ax = nodes[2*tris[3*t]], ay = nodes[2*tris[3*t]+1];
        const bx = nodes[2*tris[3*t+1]], by = nodes[2*tris[3*t+1]+1];
        const cx = nodes[2*tris[3*t+2]], cy = nodes[2*tris[3*t+2]+1];
        const al = Math.sqrt((bx-cx)**2+(by-cy)**2);
        const bl = Math.sqrt((ax-cx)**2+(ay-cy)**2);
        const cl = Math.sqrt((ax-bx)**2+(ay-by)**2);
        const s = (al+bl+cl)/2;
        const area = Math.abs((bx-ax)*(cy-ay)-(cx-ax)*(by-ay))/2;
        if (area < 1e-30) { degenerateCount++; continue; }
        const q = al*bl*cl/(8*area*(area/s));
        if (q > maxQ) { maxQ = q; worstTri = t; }
        sumQ += q;
        if (q > 5) badCount++;
    }

    // --- Constraint crossings ---
    let crossings = 0;
    const cYs = constraintYs || [], cXs = constraintXs || [];
    const cyR = mesh.constraintYRanges || {};
    for (let t = 0; t < nTris; t++) {
        const ys = [nodes[2*tris[3*t]+1], nodes[2*tris[3*t+1]+1], nodes[2*tris[3*t+2]+1]];
        const xs = [nodes[2*tris[3*t]], nodes[2*tris[3*t+1]], nodes[2*tris[3*t+2]]];
        const yMin = Math.min(...ys), yMax = Math.max(...ys);
        const xMin = Math.min(...xs), xMax = Math.max(...xs);
        for (const cy of cYs) {
            if (yMin < cy - 1e-10 && yMax > cy + 1e-10) {
                const r = cyR[cy];
                if (!r || (xMax > r[0] + 1e-10 && xMin < r[1] - 1e-10)) { crossings++; break; }
            }
        }
        for (const cx of cXs) {
            if (xMin < cx - 1e-10 && xMax > cx + 1e-10) { crossings++; break; }
        }
    }

    // --- Missing constraint edges ---
    let missingEdges = 0;
    const edgeSet = new Set();
    for (let e = 0; e < nEdges; e++) {
        const a = edges[2*e], b = edges[2*e+1];
        edgeSet.add(Math.min(a,b)+','+Math.max(a,b));
    }
    function checkLine(axis, val, lo, hi) {
        const pts = [];
        for (let n = 0; n < nNodes; n++) {
            const coord = axis === 'y' ? nodes[2*n+1] : nodes[2*n];
            const pos = axis === 'y' ? nodes[2*n] : nodes[2*n+1];
            if (Math.abs(coord - val) < 1e-10 && pos >= lo - 1e-10 && pos <= hi + 1e-10) {
                pts.push({ pos, n });
            }
        }
        pts.sort((a, b) => a.pos - b.pos);
        for (let i = 0; i < pts.length - 1; i++) {
            const a = pts[i].n, b = pts[i+1].n;
            if (!edgeSet.has(Math.min(a,b)+','+Math.max(a,b))) missingEdges++;
        }
    }
    const xRange = [Infinity, -Infinity], yRange = [Infinity, -Infinity];
    for (let n = 0; n < nNodes; n++) {
        xRange[0] = Math.min(xRange[0], nodes[2*n]);
        xRange[1] = Math.max(xRange[1], nodes[2*n]);
        yRange[0] = Math.min(yRange[0], nodes[2*n+1]);
        yRange[1] = Math.max(yRange[1], nodes[2*n+1]);
    }
    for (const cy of cYs) {
        const r = cyR[cy];
        checkLine('y', cy, r ? r[0] : xRange[0], r ? r[1] : xRange[1]);
    }
    for (const cx of cXs) checkLine('x', cx, yRange[0], yRange[1]);

    // --- Extreme area ratio (indicates ill-conditioned FEM matrices) ---
    let minArea = Infinity, maxArea = 0;
    for (let t = 0; t < nTris; t++) {
        const ax = nodes[2*tris[3*t]], ay = nodes[2*tris[3*t]+1];
        const bx = nodes[2*tris[3*t+1]], by = nodes[2*tris[3*t+1]+1];
        const cx = nodes[2*tris[3*t+2]], cy = nodes[2*tris[3*t+2]+1];
        const area = Math.abs((bx-ax)*(cy-ay)-(cx-ax)*(by-ay))/2;
        if (area > 1e-30 && area < minArea) minArea = area;
        if (area > maxArea) maxArea = area;
    }
    const areaRatio = minArea > 0 ? maxArea / minArea : Infinity;

    // --- NaN/Inf check on node coordinates ---
    let nanNodes = 0;
    for (let n = 0; n < nNodes; n++) {
        if (!isFinite(nodes[2*n]) || !isFinite(nodes[2*n+1])) nanNodes++;
    }

    // --- Build result ---
    const badFraction = nTris > 0 ? badCount / nTris : 0;
    const metrics = { maxQ, avgQ: nTris > 0 ? sumQ / nTris : 0, badCount, badFraction,
                      degenerateCount, crossings, missingEdges, areaRatio, minArea, nanNodes };

    if (nanNodes > 0) errors.push(`${nanNodes} nodes with NaN/Inf coordinates`);
    if (crossings > 0) errors.push(`${crossings} triangles cross constraint lines`);
    if (degenerateCount > 0) errors.push(`${degenerateCount} degenerate triangles (zero area)`);
    if (missingEdges > 0) errors.push(`${missingEdges} missing constraint edges`);
    if (areaRatio > 1e6) errors.push(`area ratio ${areaRatio.toExponential(1)} (extreme element size variation, will cause ill-conditioning)`);
    if (maxQ > 10) warnings.push(`max Q=${maxQ.toFixed(1)} (poor quality triangle)`);
    if (badFraction > 0.05) warnings.push(`${badCount}/${nTris} (${(badFraction*100).toFixed(1)}%) triangles with Q>5`);
    if (areaRatio > 1e4 && areaRatio <= 1e6) warnings.push(`area ratio ${areaRatio.toExponential(1)} (large element size variation)`);

    return { ok: errors.length === 0, warnings, errors, metrics };
}
