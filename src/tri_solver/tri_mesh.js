// Mesh utilities for 2D FEM: quality checking

import { triQuality } from './fem_core.js';

// --- Mesh quality check ---
// Validates mesh quality before FEM solve. Returns { ok, warnings, errors, metrics }.
// opts.skip: per-triangle mask (truthy = leave the triangle out of the shape
// statistics: Q, badFraction, area ratio). The NaN node check always covers the
// whole mesh.
export function checkMeshQuality(mesh, opts = {}) {
    const { nodes, tris, nTris, nNodes } = mesh;
    const warnings = [], errors = [];
    const skip = opts.skip || null;

    // --- Triangle quality (Q = circumradius / (2 * inradius), ideal = 1) ---
    let maxQ = 0, sumQ = 0, badCount = 0, degenerateCount = 0, worstTri = -1, nRated = 0;
    for (let t = 0; t < nTris; t++) {
        if (skip && skip[t]) continue;
        nRated++;
        const ax = nodes[2*tris[3*t]], ay = nodes[2*tris[3*t]+1];
        const bx = nodes[2*tris[3*t+1]], by = nodes[2*tris[3*t+1]+1];
        const cx = nodes[2*tris[3*t+2]], cy = nodes[2*tris[3*t+2]+1];
        const q = triQuality(ax, ay, bx, by, cx, cy);
        if (q >= 1e10) { degenerateCount++; continue; }
        if (q > maxQ) { maxQ = q; worstTri = t; }
        sumQ += q;
        if (q > 5) badCount++;
    }

    // --- Extreme area ratio (indicates ill-conditioned FEM matrices) ---
    let minArea = Infinity, maxArea = 0;
    for (let t = 0; t < nTris; t++) {
        if (skip && skip[t]) continue;
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
    const badFraction = nRated > 0 ? badCount / nRated : 0;
    const metrics = { maxQ, avgQ: nRated > 0 ? sumQ / nRated : 0, badCount, badFraction,
                      degenerateCount, areaRatio, minArea, nanNodes };

    if (nanNodes > 0) errors.push(`${nanNodes} nodes with NaN/Inf coordinates`);
    if (degenerateCount > 0) errors.push(`${degenerateCount} degenerate triangles (zero area)`);
    if (areaRatio > 1e6) errors.push(`area ratio ${areaRatio.toExponential(1)} (extreme element size variation, will cause ill-conditioning)`);
    if (maxQ > 10) warnings.push(`max Q=${maxQ.toFixed(1)} (poor quality triangle)`);
    if (badFraction > 0.05) warnings.push(`${badCount}/${nRated} (${(badFraction*100).toFixed(1)}%) triangles with Q>5`);
    if (areaRatio > 1e4 && areaRatio <= 1e6) warnings.push(`area ratio ${areaRatio.toExponential(1)} (large element size variation)`);

    return { ok: errors.length === 0, warnings, errors, metrics };
}
