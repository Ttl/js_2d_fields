// Microstrip analysis functions for triangular P2 Nedelec FEM

import { csqrt } from './fem_core.js';
import { ne1, ne2, ne3, nf1, nf2, nf3, nf4, nf5, nf6,
         lvGrad, leGrad, le3Grad, lf3Grad, triCoefficients,
         ne1Curl, ne2Curl, ne3Curl, nf1Curl, nf2Curl, nf3Curl, nf4Curl, nf5Curl, nf6Curl,
         QW, QL1, QL2, QL3, NQ, QW12, QL1_12, QL2_12, QL3_12, NQ12,
         getLzOffsets } from './tri_fem.js';

// --- Evaluate e_t and grad(ez) at point (px,py) inside triangle t ---
// Returns the Lee-Jin eigenvector fields: e_t (= γ·Et) and ∇ez (= ∇Ez).
// H can be computed from these as: jωμ₀Ht = ẑ×∇ez + ẑ×e_t (no γ multiply).
// Optional enrichment: { corners, evalFn, vecRe, vecIm, nStd, nCorners }
//   Adds ∇φ_ci·coeff to both e_t and ∇ez (de Rham compatible enrichment).
export function evalFieldsAtPoint(t, px, py, mesh, fm, vecRe, vecIm, enrichment) {
    const { tris, triEdges, triSigns, nodes } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF,
            elemOrder, edgeF3, faceF3, edgeNodeF3, faceNodeF } = fm;
    const { lzOff, lzEdgeMidOff, lzEdge3Off, lzFaceNodeOff } = getLzOffsets(fm);
    const _edgeVerts = [[0, 1], [1, 2], [2, 0]];

    const vv0 = tris[3*t], vv1 = tris[3*t+1], vv2 = tris[3*t+2];
    const { coeff: cf } = triCoefficients(nodes, vv0, vv1, vv2);
    const order = elemOrder ? elemOrder[t] : 2;
    const nNed = order >= 3 ? 15 : 8;
    const nLag = order >= 3 ? 10 : 6;

    // Gather transverse DOFs
    const eR = new Float64Array(nNed), eI = new Float64Array(nNed);
    for (let k = 0; k < 3; k++) {
        const eIdx = triEdges[3*t+k], s = triSigns[3*t+k];
        const ef1 = edgeF[2*eIdx]; if (ef1 >= 0) { eR[k] = s*vecRe[ef1]; eI[k] = s*vecIm[ef1]; }
        const ef2 = edgeF[2*eIdx+1]; if (ef2 >= 0) { eR[k+4] = vecRe[ef2]; eI[k+4] = vecIm[ef2]; }
        if (order >= 3 && edgeF3) {
            const ef3 = edgeF3[eIdx]; if (ef3 >= 0) { eR[k+8] = s*vecRe[ef3]; eI[k+8] = s*vecIm[ef3]; }
        }
    }
    const ff1 = faceF[2*t]; if (ff1 >= 0) { eR[3] = vecRe[ff1]; eI[3] = vecIm[ff1]; }
    const ff2 = faceF[2*t+1]; if (ff2 >= 0) { eR[7] = vecRe[ff2]; eI[7] = vecIm[ff2]; }
    if (order >= 3 && faceF3) {
        for (let k = 0; k < 4; k++) { const ff = faceF3[4*t+k]; if (ff >= 0) { eR[11+k] = vecRe[ff]; eI[11+k] = vecIm[ff]; } }
    }

    // Evaluate e_t from all Nedelec DOFs
    let exr = 0, exi = 0, eyr = 0, eyi = 0;
    for (let k = 0; k < 3; k++) {
        const [p, q] = _edgeVerts[k];
        {const [wx, wy] = ne1(cf, p, q, px, py); exr+=wx*eR[k]; exi+=wx*eI[k]; eyr+=wy*eR[k]; eyi+=wy*eI[k];}
        {const [wx, wy] = ne2(cf, p, q, px, py); exr+=wx*eR[k+4]; exi+=wx*eI[k+4]; eyr+=wy*eR[k+4]; eyi+=wy*eI[k+4];}
        if (order >= 3) {const [wx, wy] = ne3(cf, p, q, px, py); exr+=wx*eR[k+8]; exi+=wx*eI[k+8]; eyr+=wy*eR[k+8]; eyi+=wy*eI[k+8];}
    }
    { const [wx,wy] = nf1(cf, px, py); exr+=wx*eR[3]; exi+=wx*eI[3]; eyr+=wy*eR[3]; eyi+=wy*eI[3]; }
    { const [wx,wy] = nf2(cf, px, py); exr+=wx*eR[7]; exi+=wx*eI[7]; eyr+=wy*eR[7]; eyi+=wy*eI[7]; }
    if (order >= 3) {
        { const [wx,wy] = nf3(cf, px, py); exr+=wx*eR[11]; exi+=wx*eI[11]; eyr+=wy*eR[11]; eyi+=wy*eI[11]; }
        { const [wx,wy] = nf4(cf, px, py); exr+=wx*eR[12]; exi+=wx*eI[12]; eyr+=wy*eR[12]; eyi+=wy*eI[12]; }
        { const [wx,wy] = nf5(cf, px, py); exr+=wx*eR[13]; exi+=wx*eI[13]; eyr+=wy*eR[13]; eyi+=wy*eI[13]; }
        { const [wx,wy] = nf6(cf, px, py); exr+=wx*eR[14]; exi+=wx*eI[14]; eyr+=wy*eR[14]; eyi+=wy*eI[14]; }
    }

    // Gather longitudinal DOFs and evaluate grad(Ez)
    const nDR = new Float64Array(nLag), nDI = new Float64Array(nLag);
    const vts = [vv0, vv1, vv2];
    for (let k = 0; k < 3; k++) { const nf = nodeF[vts[k]]; if (nf >= 0) { nDR[k] = vecRe[lzOff+nf]; nDI[k] = vecIm[lzOff+nf]; } }
    for (let k = 0; k < 3; k++) { const enf = edgeNodeF[triEdges[3*t+k]]; if (enf >= 0) { nDR[k+3] = vecRe[lzEdgeMidOff+enf]; nDI[k+3] = vecIm[lzEdgeMidOff+enf]; } }
    if (order >= 3) {
        for (let k = 0; k < 3; k++) {
            const enf3 = edgeNodeF3 ? edgeNodeF3[triEdges[3*t+k]] : -1;
            if (enf3 >= 0) { const s = triSigns[3*t+k]; nDR[k+6] = s*vecRe[lzEdge3Off+enf3]; nDI[k+6] = s*vecIm[lzEdge3Off+enf3]; }
        }
        if (faceNodeF) { const fnf = faceNodeF[t]; if (fnf >= 0) { nDR[9] = vecRe[lzFaceNodeOff+fnf]; nDI[9] = vecIm[lzFaceNodeOff+fnf]; } }
    }

    let dezdxr = 0, dezdxi = 0, dezdyr = 0, dezdyi = 0;
    for (let k = 0; k < 3; k++) { const [gx,gy] = lvGrad(cf, k, px, py); dezdxr += gx*nDR[k]; dezdxi += gx*nDI[k]; dezdyr += gy*nDR[k]; dezdyi += gy*nDI[k]; }
    for (let k = 0; k < 3; k++) { const [p,q] = _edgeVerts[k]; const [gx,gy] = leGrad(cf, p, q, px, py); dezdxr += gx*nDR[k+3]; dezdxi += gx*nDI[k+3]; dezdyr += gy*nDR[k+3]; dezdyi += gy*nDI[k+3]; }
    if (order >= 3) {
        for (let k = 0; k < 3; k++) { const [p,q] = _edgeVerts[k]; const [gx,gy] = le3Grad(cf, p, q, px, py); dezdxr += gx*nDR[k+6]; dezdxi += gx*nDI[k+6]; dezdyr += gy*nDR[k+6]; dezdyi += gy*nDI[k+6]; }
        { const [gx,gy] = lf3Grad(cf, px, py); dezdxr += gx*nDR[9]; dezdxi += gx*nDI[9]; dezdyr += gy*nDR[9]; dezdyi += gy*nDI[9]; }
    }

    // Add enrichment contributions to e_t and ∇ez
    if (enrichment) {
        const { corners, evalFn, vecRe: evR, vecIm: evI, nStd, nCorners } = enrichment;
        for (let ci = 0; ci < nCorners; ci++) {
            const { dphidx, dphidy } = evalFn(px, py, corners[ci]);
            if (dphidx === 0 && dphidy === 0) continue;
            const ctR = evR[nStd + ci], ctI = evI[nStd + ci];
            exr += dphidx * ctR; exi += dphidx * ctI;
            eyr += dphidy * ctR; eyi += dphidy * ctI;
            const czR = evR[nStd + nCorners + ci], czI = evI[nStd + nCorners + ci];
            dezdxr += dphidx * czR; dezdxi += dphidx * czI;
            dezdyr += dphidy * czR; dezdyi += dphidy * czI;
        }
    }

    // curl(Et)_z = ∂Ey/∂x - ∂Ex/∂y — computed from analytical Nedelec curl.
    let curlEtRe = 0, curlEtIm = 0;
    for (let k = 0; k < 3; k++) {
        const [p, q] = _edgeVerts[k];
        {const c = ne1Curl(cf, p, q); curlEtRe += c*eR[k]; curlEtIm += c*eI[k];}
        {const c = ne2Curl(cf, p, q, px, py); curlEtRe += c*eR[k+4]; curlEtIm += c*eI[k+4];}
        if (order >= 3) {const c = ne3Curl(cf, p, q, px, py); curlEtRe += c*eR[k+8]; curlEtIm += c*eI[k+8];}
    }
    { const c = nf1Curl(cf, px, py); curlEtRe += c * eR[3]; curlEtIm += c * eI[3]; }
    { const c = nf2Curl(cf, px, py); curlEtRe += c * eR[7]; curlEtIm += c * eI[7]; }
    if (order >= 3) {
        { const c = nf3Curl(cf, px, py); curlEtRe += c*eR[11]; curlEtIm += c*eI[11]; }
        { const c = nf4Curl(cf, px, py); curlEtRe += c*eR[12]; curlEtIm += c*eI[12]; }
        { const c = nf5Curl(cf, px, py); curlEtRe += c*eR[13]; curlEtIm += c*eI[13]; }
        { const c = nf6Curl(cf, px, py); curlEtRe += c*eR[14]; curlEtIm += c*eI[14]; }
    }

    return { exr, exi, eyr, eyi, dezdxr, dezdxi, dezdyr, dezdyi, curlEtRe, curlEtIm };
}

// --- Element-interior contour for I computation ---
// Finds the ring of mesh edges one layer out from the conductor. For each edge in the
// ring, evaluates H at a point INSIDE one of the adjacent triangles (offset from the
// edge midpoint toward the triangle centroid). This ensures every H evaluation is inside
// a single element and the integration path follows actual mesh edges.

export function buildElementContour(mesh, fm, condRect, distFrac = 1/3) {
    const { nodes, tris, edges, triEdges, nTris, nEdges, nNodes } = mesh;
    const { isCondNode, isCondEdge } = fm;
    const TOL = 1e-12;

    // Build edge→triangle adjacency
    const edgeToTriList = new Array(nEdges);
    for (let e = 0; e < nEdges; e++) edgeToTriList[e] = [];
    for (let t = 0; t < nTris; t++) {
        for (let le = 0; le < 3; le++) {
            edgeToTriList[triEdges[3*t+le]].push(t);
        }
    }

    // Check which triangles are inside conductor
    const isCondTri = new Uint8Array(nTris);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const xc = (nodes[2*v0] + nodes[2*v1] + nodes[2*v2]) / 3;
        const yc = (nodes[2*v0+1] + nodes[2*v1+1] + nodes[2*v2+1]) / 3;
        if (xc >= condRect.xmin - TOL && xc <= condRect.xmax + TOL &&
            yc >= condRect.ymin - TOL && yc <= condRect.ymax + TOL)
            isCondTri[t] = 1;
    }

    // Compute target distance: fraction of min clearance to nearest boundary.
    const clearGnd = condRect.ymin;
    const clearTop = condRect.ymax_domain - condRect.ymax;
    const clearL = condRect.xmin - condRect.xmin_domain;
    const clearR = condRect.xmax_domain - condRect.xmax;
    const minClear = Math.min(clearGnd, clearTop, clearL, clearR);
    const targetDist = minClear * distFrac;

    // BFS outward from conductor, layer by layer. For each layer, check if the
    // outer boundary forms a single closed loop at sufficient distance.
    const allExpanded = new Set();
    // Seed: conductor-adjacent non-conductor triangles
    for (let e = 0; e < nEdges; e++) {
        if (!isCondEdge[e]) continue;
        const n0 = edges[2*e], n1 = edges[2*e+1];
        if (Math.abs(nodes[2*n0+1]) < TOL && Math.abs(nodes[2*n1+1]) < TOL) continue;
        for (const t of edgeToTriList[e]) {
            if (!isCondTri[t]) allExpanded.add(t);
        }
    }
    if (allExpanded.size === 0) return null;

    // Find outer boundary of current expanded region, check closure, expand if needed
    let ringEdges, ringEdgeTri, closedLoop = false;

    for (let layer = 0; layer < 30; layer++) {
        // Find boundary: interior edges where one side is expanded, other is not
        ringEdges = [];
        ringEdgeTri = [];
        for (let e = 0; e < nEdges; e++) {
            const tl = edgeToTriList[e];
            if (tl.length !== 2) continue; // skip domain boundary edges
            const in0 = allExpanded.has(tl[0]), in1 = allExpanded.has(tl[1]);
            const c0 = isCondTri[tl[0]], c1 = isCondTri[tl[1]];
            // One side expanded, other side is non-conductor and not expanded
            if (in0 && !in1 && !c1) {
                ringEdges.push(e); ringEdgeTri.push(tl[0]);
            } else if (in1 && !in0 && !c0) {
                ringEdges.push(e); ringEdgeTri.push(tl[1]);
            }
        }

        if (ringEdges.length === 0) break;

        // Try to order into a closed loop
        const vertToIdx = new Map();
        for (let i = 0; i < ringEdges.length; i++) {
            const e = ringEdges[i];
            const n0 = edges[2*e], n1 = edges[2*e+1];
            if (!vertToIdx.has(n0)) vertToIdx.set(n0, []);
            if (!vertToIdx.has(n1)) vertToIdx.set(n1, []);
            vertToIdx.get(n0).push(i);
            vertToIdx.get(n1).push(i);
        }

        const used = new Uint8Array(ringEdges.length);
        used[0] = 1;
        const ordered = [0];
        const firstEdge = ringEdges[0];
        const startVert = edges[2*firstEdge];
        let curVert = edges[2*firstEdge+1];

        for (let step = 0; step < ringEdges.length; step++) {
            const cands = vertToIdx.get(curVert);
            if (!cands) break;
            let found = false;
            for (const idx of cands) {
                if (used[idx]) continue;
                used[idx] = 1;
                ordered.push(idx);
                const e = ringEdges[idx];
                curVert = (edges[2*e] === curVert) ? edges[2*e+1] : edges[2*e];
                found = true;
                break;
            }
            if (!found) break;
        }

        closedLoop = (ordered.length === ringEdges.length && curVert === startVert);

        // Check average distance
        let avgDist = 0;
        for (const e of ringEdges) {
            const n0 = edges[2*e], n1 = edges[2*e+1];
            const mx = (nodes[2*n0] + nodes[2*n1]) / 2;
            const my = (nodes[2*n0+1] + nodes[2*n1+1]) / 2;
            const cx = Math.max(condRect.xmin, Math.min(condRect.xmax, mx));
            const cy = Math.max(condRect.ymin, Math.min(condRect.ymax, my));
            avgDist += Math.sqrt((mx-cx)**2 + (my-cy)**2);
        }
        avgDist /= ringEdges.length;

        if (closedLoop && avgDist >= targetDist) break;

        // Expand one more layer
        const frontier = [];
        for (const e of ringEdges) {
            for (const t of edgeToTriList[e]) {
                if (!allExpanded.has(t) && !isCondTri[t]) frontier.push(t);
            }
        }
        if (frontier.length === 0) break;
        for (const t of frontier) allExpanded.add(t);
    }

    if (!closedLoop || ringEdges.length === 0) return null;

    // Order ring edges into a closed loop using vertex connectivity
    const vertToRingEdge = new Map();
    for (let i = 0; i < ringEdges.length; i++) {
        const e = ringEdges[i];
        const n0 = edges[2*e], n1 = edges[2*e+1];
        if (!vertToRingEdge.has(n0)) vertToRingEdge.set(n0, []);
        if (!vertToRingEdge.has(n1)) vertToRingEdge.set(n1, []);
        vertToRingEdge.get(n0).push(i);
        vertToRingEdge.get(n1).push(i);
    }

    const ordered = [];
    const used = new Uint8Array(ringEdges.length);
    used[0] = 1;
    ordered.push(0);
    let curVert = edges[2*ringEdges[0]+1]; // start from second vertex of first edge

    for (let step = 0; step < ringEdges.length; step++) {
        const candidates = vertToRingEdge.get(curVert);
        if (!candidates) break;
        let found = false;
        for (const idx of candidates) {
            if (used[idx]) continue;
            used[idx] = 1;
            ordered.push(idx);
            const e = ringEdges[idx];
            const n0 = edges[2*e], n1 = edges[2*e+1];
            curVert = (n0 === curVert) ? n1 : n0;
            found = true;
            break;
        }
        if (!found) break;
    }

    // Measure distance from ring edges to conductor surface
    let dMin = Infinity, dMax = 0, dSum = 0;
    for (let i = 0; i < ringEdges.length; i++) {
        const e = ringEdges[i];
        const n0 = edges[2*e], n1 = edges[2*e+1];
        const mx = (nodes[2*n0] + nodes[2*n1]) / 2;
        const my = (nodes[2*n0+1] + nodes[2*n1+1]) / 2;
        // Distance to nearest conductor boundary point
        const cx = Math.max(condRect.xmin, Math.min(condRect.xmax, mx));
        const cy = Math.max(condRect.ymin, Math.min(condRect.ymax, my));
        const d = Math.sqrt((mx-cx)**2 + (my-cy)**2);
        dMin = Math.min(dMin, d); dMax = Math.max(dMax, d); dSum += d;
    }
    // (diagnostics printed after filtering)

    // Build segments: for each ordered ring edge, evaluate H inside the ring triangle
    // at a point offset from edge midpoint toward the triangle centroid.
    const alpha = 0.3; // how far from edge midpoint toward centroid
    const segs = [];
    let nSkipped = 0;
    for (let i = 0; i < ordered.length; i++) {
        const idx = ordered[i];
        const e = ringEdges[idx];
        const t = ringEdgeTri[idx];
        const n0 = edges[2*e], n1 = edges[2*e+1];

        // Edge midpoint
        const emx = (nodes[2*n0] + nodes[2*n1]) / 2;
        const emy = (nodes[2*n0+1] + nodes[2*n1+1]) / 2;

        // Distance to conductor surface
        const cpx = Math.max(condRect.xmin, Math.min(condRect.xmax, emx));
        const cpy = Math.max(condRect.ymin, Math.min(condRect.ymax, emy));
        const dCond = Math.sqrt((emx-cpx)**2 + (emy-cpy)**2);

        // Skip edges too close to conductor (corner-adjacent, poor field quality)
        if (dCond < targetDist * 0.5) { nSkipped++; continue; }

        // Triangle centroid
        const tv0 = tris[3*t], tv1 = tris[3*t+1], tv2 = tris[3*t+2];
        const tcx = (nodes[2*tv0] + nodes[2*tv1] + nodes[2*tv2]) / 3;
        const tcy = (nodes[2*tv0+1] + nodes[2*tv1+1] + nodes[2*tv2+1]) / 3;

        // Evaluation point: offset from edge midpoint toward centroid
        const px = emx + alpha * (tcx - emx);
        const py = emy + alpha * (tcy - emy);

        // Determine CCW direction: interior (expanded region) on the LEFT.
        let oppV = -1;
        if (tv0 !== n0 && tv0 !== n1) oppV = tv0;
        else if (tv1 !== n0 && tv1 !== n1) oppV = tv1;
        else oppV = tv2;
        const oppX = nodes[2*oppV], oppY = nodes[2*oppV+1];

        // Edge vector n0→n1
        const tx = nodes[2*n1] - nodes[2*n0], ty = nodes[2*n1+1] - nodes[2*n0+1];
        // Cross product: if oppV is on LEFT of n0→n1, then n0→n1 is CCW
        const cross = tx * (oppY - nodes[2*n0+1]) - ty * (oppX - nodes[2*n0]);
        const sign = cross > 0 ? 1 : -1;

        segs.push({
            xm: px, ym: py, tri: t,
            dx: sign * tx, dy: sign * ty,
        });
    }
    if (nSkipped > 0) globalThis.__TRI_DEBUG__ && console.log('  Ring: ' + segs.length + ' used, ' + nSkipped + ' skipped (too close)');

    return segs;
}

// --- Mode analysis ---
// Note: enrichment DOFs are NOT included here. The enrichment basis has compact support
// near PEC corners with small coefficients (~1e-4), so its contribution to full-domain
// integrals (Poynting power, voltage, ∮H·dl) is negligible for impedance computation.
// Conductor loss handles enrichment separately via the Galerkin H projection in
// conductor_loss.js::projectH, where the singular corner fields matter most.

export function analyzeTriMode(vecRe, vecIm, gamma2Re, gamma2Im, mesh, fm, k2, f, epsMap, condRect, distFrac) {
    const { nodes, tris, edges, triEdges, triSigns, nTris, nEdges } = mesh;
    const { edgeF, faceF, nodeF, edgeNodeF, isCondNode,
            elemOrder, edgeF3, faceF3, edgeNodeF3, faceNodeF } = fm;
    const nNodes = mesh.nNodes;
    const omega = 2 * Math.PI * f;
    const mu0 = 4 * Math.PI * 1e-7;

    const omu = omega * mu0;
    const { lzOff, lzEdgeMidOff, lzEdge3Off, lzFaceNodeOff } = getLzOffsets(fm);

    let gamma = csqrt(gamma2Re, gamma2Im);
    if (gamma.im < 0) gamma = { re: -gamma.re, im: -gamma.im };
    const gMag2 = gamma.re*gamma.re + gamma.im*gamma.im;

    const edgeVerts = [[0, 1], [1, 2], [2, 0]];

    let SetR = 0, SetI = 0;
    let Wm = 0;
    let Pdiel_et = 0;

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const txs = [nodes[2*v0], nodes[2*v1], nodes[2*v2]];
        const tys = [nodes[2*v0+1], nodes[2*v1+1], nodes[2*v2+1]];
        const { coeff, Area } = triCoefficients(nodes, v0, v1, v2);
        const order = elemOrder ? elemOrder[t] : 2;
        const nNed = order >= 3 ? 15 : 8;
        const nLag = order >= 3 ? 10 : 6;

        // Gather transverse DOFs
        const eDofRe = new Float64Array(nNed), eDofIm = new Float64Array(nNed);
        for (let k = 0; k < 3; k++) {
            const eIdx = triEdges[3*t+k], s = triSigns[3*t+k];
            const ef1 = edgeF[2*eIdx]; if (ef1 >= 0) { eDofRe[k] = s*vecRe[ef1]; eDofIm[k] = s*vecIm[ef1]; }
            const ef2 = edgeF[2*eIdx+1]; if (ef2 >= 0) { eDofRe[k+4] = vecRe[ef2]; eDofIm[k+4] = vecIm[ef2]; }
            if (order >= 3 && edgeF3) { const ef3 = edgeF3[eIdx]; if (ef3 >= 0) { eDofRe[k+8] = s*vecRe[ef3]; eDofIm[k+8] = s*vecIm[ef3]; } }
        }
        { const ff = faceF[2*t]; if (ff >= 0) { eDofRe[3] = vecRe[ff]; eDofIm[3] = vecIm[ff]; } }
        { const ff = faceF[2*t+1]; if (ff >= 0) { eDofRe[7] = vecRe[ff]; eDofIm[7] = vecIm[ff]; } }
        if (order >= 3 && faceF3) { for (let k = 0; k < 4; k++) { const ff = faceF3[4*t+k]; if (ff >= 0) { eDofRe[11+k] = vecRe[ff]; eDofIm[11+k] = vecIm[ff]; } } }

        // Gather longitudinal DOFs
        const nDofRe = new Float64Array(nLag), nDofIm = new Float64Array(nLag);
        const verts = [v0, v1, v2];
        for (let k = 0; k < 3; k++) { const nf = nodeF[verts[k]]; if (nf >= 0) { nDofRe[k] = vecRe[lzOff+nf]; nDofIm[k] = vecIm[lzOff+nf]; } }
        for (let k = 0; k < 3; k++) { const enf = edgeNodeF[triEdges[3*t+k]]; if (enf >= 0) { nDofRe[k+3] = vecRe[lzEdgeMidOff+enf]; nDofIm[k+3] = vecIm[lzEdgeMidOff+enf]; } }
        if (order >= 3) {
            for (let k = 0; k < 3; k++) { const enf3 = edgeNodeF3 ? edgeNodeF3[triEdges[3*t+k]] : -1; if (enf3 >= 0) { const s = triSigns[3*t+k]; nDofRe[k+6] = s*vecRe[lzEdge3Off+enf3]; nDofIm[k+6] = s*vecIm[lzEdge3Off+enf3]; } }
            if (faceNodeF) { const fnf = faceNodeF[t]; if (fnf >= 0) { nDofRe[9] = vecRe[lzFaceNodeOff+fnf]; nDofIm[9] = vecIm[lzFaceNodeOff+fnf]; } }
        }

        const useQ12 = order >= 3;
        const nqp = useQ12 ? NQ12 : NQ;
        const qw = useQ12 ? QW12 : QW;
        const ql1 = useQ12 ? QL1_12 : QL1, ql2 = useQ12 ? QL2_12 : QL2, ql3 = useQ12 ? QL3_12 : QL3;

        for (let q = 0; q < nqp; q++) {
            const xq = txs[0]*ql1[q] + txs[1]*ql2[q] + txs[2]*ql3[q];
            const yq = tys[0]*ql1[q] + tys[1]*ql2[q] + tys[2]*ql3[q];

            // Evaluate Et
            let exr = 0, exi = 0, eyr = 0, eyi = 0;
            for (let k = 0; k < 3; k++) {
                const [p, qq] = edgeVerts[k];
                {const [wx,wy]=ne1(coeff,p,qq,xq,yq); exr+=wx*eDofRe[k]; exi+=wx*eDofIm[k]; eyr+=wy*eDofRe[k]; eyi+=wy*eDofIm[k];}
                {const [wx,wy]=ne2(coeff,p,qq,xq,yq); exr+=wx*eDofRe[k+4]; exi+=wx*eDofIm[k+4]; eyr+=wy*eDofRe[k+4]; eyi+=wy*eDofIm[k+4];}
                if (order >= 3) {const [wx,wy]=ne3(coeff,p,qq,xq,yq); exr+=wx*eDofRe[k+8]; exi+=wx*eDofIm[k+8]; eyr+=wy*eDofRe[k+8]; eyi+=wy*eDofIm[k+8];}
            }
            {const [wx,wy]=nf1(coeff,xq,yq); exr+=wx*eDofRe[3]; exi+=wx*eDofIm[3]; eyr+=wy*eDofRe[3]; eyi+=wy*eDofIm[3];}
            {const [wx,wy]=nf2(coeff,xq,yq); exr+=wx*eDofRe[7]; exi+=wx*eDofIm[7]; eyr+=wy*eDofRe[7]; eyi+=wy*eDofIm[7];}
            if (order >= 3) {
                {const [wx,wy]=nf3(coeff,xq,yq); exr+=wx*eDofRe[11]; exi+=wx*eDofIm[11]; eyr+=wy*eDofRe[11]; eyi+=wy*eDofIm[11];}
                {const [wx,wy]=nf4(coeff,xq,yq); exr+=wx*eDofRe[12]; exi+=wx*eDofIm[12]; eyr+=wy*eDofRe[12]; eyi+=wy*eDofIm[12];}
                {const [wx,wy]=nf5(coeff,xq,yq); exr+=wx*eDofRe[13]; exi+=wx*eDofIm[13]; eyr+=wy*eDofRe[13]; eyi+=wy*eDofIm[13];}
                {const [wx,wy]=nf6(coeff,xq,yq); exr+=wx*eDofRe[14]; exi+=wx*eDofIm[14]; eyr+=wy*eDofRe[14]; eyi+=wy*eDofIm[14];}
            }

            // Evaluate grad(Ez)
            let dezdxr = 0, dezdxi = 0, dezdyr = 0, dezdyi = 0;
            for (let k = 0; k < 3; k++) { const [gx,gy]=lvGrad(coeff,k,xq,yq); dezdxr+=gx*nDofRe[k]; dezdxi+=gx*nDofIm[k]; dezdyr+=gy*nDofRe[k]; dezdyi+=gy*nDofIm[k]; }
            for (let k = 0; k < 3; k++) { const [p,qq]=edgeVerts[k]; const [gx,gy]=leGrad(coeff,p,qq,xq,yq); dezdxr+=gx*nDofRe[k+3]; dezdxi+=gx*nDofIm[k+3]; dezdyr+=gy*nDofRe[k+3]; dezdyi+=gy*nDofIm[k+3]; }
            if (order >= 3) {
                for (let k = 0; k < 3; k++) { const [p,qq]=edgeVerts[k]; const [gx,gy]=le3Grad(coeff,p,qq,xq,yq); dezdxr+=gx*nDofRe[k+6]; dezdxi+=gx*nDofIm[k+6]; dezdyr+=gy*nDofRe[k+6]; dezdyi+=gy*nDofIm[k+6]; }
                { const [gx,gy]=lf3Grad(coeff,xq,yq); dezdxr+=gx*nDofRe[9]; dezdxi+=gx*nDofIm[9]; dezdyr+=gy*nDofRe[9]; dezdyi+=gy*nDofIm[9]; }
            }

            const nHxr = -(dezdyr + eyr), nHxi = -(dezdyi + eyi);
            const nHyr = exr + dezdxr, nHyi = exi + dezdxi;
            const hxr = nHxi/omu, hxi = -nHxr/omu;
            const hyr = nHyi/omu, hyi = -nHyr/omu;

            const szR = (exr*hyr + exi*hyi) - (eyr*hxr + eyi*hxi);
            const szI = (exi*hyr - exr*hyi) - (eyi*hxr - eyr*hxi);
            SetR += qw[q] * Area * szR;
            SetI += qw[q] * Area * szI;

            const hMag2 = hxr*hxr + hxi*hxi + hyr*hyr + hyi*hyi;
            Wm += qw[q] * Area * 0.5 * mu0 * hMag2;

            const epsIm = epsMap ? Math.abs(epsMap[t].im) : 0;
            if (epsIm > 0) {
                const et2 = exr*exr + exi*exi + eyr*eyr + eyi*eyi;
                Pdiel_et += qw[q] * Area * 0.5 * omega * 8.854187817e-12 * epsIm * et2;
            }
        }
    }

    // Convert accumulated e_t integrals to physical values via single γ division
    // P = 0.5·Re(S_et / γ) where S_et = ∫(e_t × Ht*) dA
    const Ppoynting = 0.5 * (SetR * gamma.re + SetI * gamma.im) / gMag2;
    // Pdiel_phys = Pdiel_et / |γ|² (since |Et|² = |e_t|²/|γ|²)
    const Pdiel = Pdiel_et / gMag2;

    // --- Voltage: path integral from ground to conductor ---
    const xc = (condRect.xmin + condRect.xmax) / 2;
    const TOL2 = (condRect.xmax - condRect.xmin) * 0.6;

    const nodeAdj = new Array(nNodes);
    for (let n = 0; n < nNodes; n++) nodeAdj[n] = [];
    for (let e = 0; e < nEdges; e++) {
        const n0 = edges[2*e], n1 = edges[2*e+1];
        nodeAdj[n0].push({ edge: e, neighbor: n1 });
        nodeAdj[n1].push({ edge: e, neighbor: n0 });
    }

    let startNode = -1, bestDist = Infinity;
    for (let n = 0; n < nNodes; n++) {
        if (Math.abs(nodes[2*n+1]) > 1e-12) continue;
        const d = Math.abs(nodes[2*n] - xc);
        if (d < bestDist) { bestDist = d; startNode = n; }
    }
    let endNode = -1; bestDist = Infinity;
    for (let n = 0; n < nNodes; n++) {
        if (!fm.isCondNode[n]) continue;
        const d = Math.abs(nodes[2*n] - xc) + Math.abs(nodes[2*n+1] - condRect.ymin) * 10;
        if (d < bestDist) { bestDist = d; endNode = n; }
    }

    const prev = new Int32Array(nNodes).fill(-1);
    const prevEdge = new Int32Array(nNodes).fill(-1);
    const visited = new Uint8Array(nNodes);
    visited[startNode] = 1;
    const queue = [startNode];
    let qi = 0;
    while (qi < queue.length) {
        const cur = queue[qi++];
        if (cur === endNode) break;
        for (const { edge, neighbor } of nodeAdj[cur]) {
            if (visited[neighbor]) continue;
            const ny = nodes[2*neighbor+1];
            if (ny < nodes[2*cur+1] - 1e-12) continue;
            if (Math.abs(nodes[2*neighbor] - xc) > TOL2) continue;
            visited[neighbor] = 1;
            prev[neighbor] = cur;
            prevEdge[neighbor] = edge;
            queue.push(neighbor);
        }
    }

    // Build edge-to-triangle map for V integration
    const edgeToTri = new Int32Array(nEdges).fill(-1);
    for (let t = 0; t < nTris; t++) {
        for (let le = 0; le < 3; le++) {
            const eIdx = triEdges[3*t+le];
            if (edgeToTri[eIdx] < 0) edgeToTri[eIdx] = t;
        }
    }

    // V from ne1 DOFs: ∫e_t·dl (accumulate in e_t, divide by γ at end)
    let Vre = 0, Vim = 0;
    let pathLen = 0;
    if (visited[endNode]) {
        let cur = endNode;
        while (cur !== startNode) {
            const e = prevEdge[cur];
            const n0 = edges[2*e], n1 = edges[2*e+1];
            const pathSign = (prev[cur] === n0) ? 1 : -1;
            // ne1: ∫ne1·dl = 1
            const ef = edgeF[2*e];
            if (ef >= 0) {
                Vre -= pathSign * vecRe[ef];
                Vim -= pathSign * vecIm[ef];
            }
            // ne3: ∫ne3·dl = 1/3 (antisymmetric like ne1, so pathSign applies)
            if (edgeF3) {
                const ef3 = edgeF3[e];
                if (ef3 >= 0) {
                    Vre -= pathSign * vecRe[ef3] / 3;
                    Vim -= pathSign * vecIm[ef3] / 3;
                }
            }
            pathLen++;
            cur = prev[cur];
        }
        // V_phys = V_et / γ (single complex division on final scalar)
        const Vr = (Vre*gamma.re + Vim*gamma.im) / gMag2;
        const Vi = (Vim*gamma.re - Vre*gamma.im) / gMag2;
        Vre = Vr; Vim = Vi;
        globalThis.__TRI_DEBUG__ && console.log(`  V path: ${pathLen} edges, V=${Math.sqrt(Vre*Vre+Vim*Vim).toExponential(3)}`);
    } else {
        globalThis.__TRI_DEBUG__ && console.log(`  V path: NOT FOUND (start=${startNode}, end=${endNode})`);
    }
    const Vmag = Math.sqrt(Vre*Vre + Vim*Vim);

    // --- Current: offset rectangle contour integral ---
    let HIre = 0, HIim = 0;

    // Build adjacency: for each triangle, which triangles share an edge
    const triNeighbors = new Array(nTris);
    for (let t = 0; t < nTris; t++) triNeighbors[t] = [];
    const edgeToTris = new Map();
    for (let t = 0; t < nTris; t++) {
        for (let le = 0; le < 3; le++) {
            const eIdx = triEdges[3*t+le];
            if (!edgeToTris.has(eIdx)) edgeToTris.set(eIdx, []);
            edgeToTris.get(eIdx).push(t);
        }
    }

    const sym = condRect.symmetry || 1;
    const offset = (condRect.ymax - condRect.ymin) * 0.1;
    const cxmin = Math.max(condRect.xmin - offset, condRect.xmin_domain || -Infinity);
    const cxmax = condRect.xmax + offset;
    const cymin = condRect.ymin - offset, cymax = condRect.ymax + offset;

    // Grid-based spatial index for O(1) point-in-triangle lookup
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    for (let n = 0; n < nNodes; n++) {
        const x = nodes[2*n], y = nodes[2*n+1];
        if (x < xMin) xMin = x; if (x > xMax) xMax = x;
        if (y < yMin) yMin = y; if (y > yMax) yMax = y;
    }
    const gridNx = Math.max(1, Math.ceil(Math.sqrt(nTris / 2)));
    const gridNy = gridNx;
    const gdx = (xMax - xMin) / gridNx + 1e-15;
    const gdy = (yMax - yMin) / gridNy + 1e-15;
    const grid = new Array(gridNx * gridNy);
    for (let i = 0; i < grid.length; i++) grid[i] = [];
    for (let t = 0; t < nTris; t++) {
        const va = tris[3*t], vb = tris[3*t+1], vc = tris[3*t+2];
        let txmin = Math.min(nodes[2*va], nodes[2*vb], nodes[2*vc]);
        let txmax = Math.max(nodes[2*va], nodes[2*vb], nodes[2*vc]);
        let tymin = Math.min(nodes[2*va+1], nodes[2*vb+1], nodes[2*vc+1]);
        let tymax = Math.max(nodes[2*va+1], nodes[2*vb+1], nodes[2*vc+1]);
        const ix0 = Math.max(0, Math.floor((txmin - xMin) / gdx));
        const ix1 = Math.min(gridNx-1, Math.floor((txmax - xMin) / gdx));
        const iy0 = Math.max(0, Math.floor((tymin - yMin) / gdy));
        const iy1 = Math.min(gridNy-1, Math.floor((tymax - yMin) / gdy));
        for (let iy = iy0; iy <= iy1; iy++)
            for (let ix = ix0; ix <= ix1; ix++)
                grid[iy * gridNx + ix].push(t);
    }
    function findTriangle(px, py) {
        const ix = Math.max(0, Math.min(gridNx-1, Math.floor((px - xMin) / gdx)));
        const iy = Math.max(0, Math.min(gridNy-1, Math.floor((py - yMin) / gdy)));
        const cell = grid[iy * gridNx + ix];
        for (let ci = 0; ci < cell.length; ci++) {
            const t = cell[ci];
            const va = tris[3*t], vb = tris[3*t+1], vc = tris[3*t+2];
            const ax=nodes[2*va],ay=nodes[2*va+1],bx=nodes[2*vb],by=nodes[2*vb+1],cx_=nodes[2*vc],cy_=nodes[2*vc+1];
            const d = (by-cy_)*(ax-cx_)+(cx_-bx)*(ay-cy_);
            const l1 = ((by-cy_)*(px-cx_)+(cx_-bx)*(py-cy_))/d;
            const l2 = ((cy_-ay)*(px-cx_)+(ax-cx_)*(py-cy_))/d;
            const l3 = 1-l1-l2;
            if (l1 >= -1e-10 && l2 >= -1e-10 && l3 >= -1e-10) return t;
        }
        return -1;
    }

    // 4 sides of offset rectangle, 8-point Gauss-Legendre per side
    const GL8p = [0.01986, 0.10167, 0.23723, 0.40828, 0.59172, 0.76277, 0.89833, 0.98014];
    const GL8w = [0.05061, 0.11119, 0.15685, 0.18134, 0.18134, 0.15685, 0.11119, 0.05061];
    // In symmetry mode, skip left segment on PMC plane (Ht≡0 analytically,
    // but field eval on boundary elements adds numerical noise)
    const skipLeft = sym > 1 && Math.abs(cxmin - (condRect.xmin_domain || 0)) < 1e-12;
    const contourSegs = [
        { x0: cxmin, y0: cymin, x1: cxmax, y1: cymin }, // bottom
        { x0: cxmax, y0: cymin, x1: cxmax, y1: cymax }, // right
        { x0: cxmax, y0: cymax, x1: cxmin, y1: cymax }, // top
    ];
    if (!skipLeft) contourSegs.push({ x0: cxmin, y0: cymax, x1: cxmin, y1: cymin }); // left

    let nFound = 0, nMissed = 0;
    for (const seg of contourSegs) {
        const sdx = seg.x1-seg.x0, sdy = seg.y1-seg.y0;
        for (let qi = 0; qi < 8; qi++) {
            const px = seg.x0 + GL8p[qi]*sdx, py = seg.y0 + GL8p[qi]*sdy;
            const t = findTriangle(px, py);
            if (t < 0) { nMissed++; continue; }
            nFound++;

            const F = evalFieldsAtPoint(t, px, py, mesh, fm, vecRe, vecIm);
            // H = -(1/(jωμ₀))(∇ez + ẑ×e_t) — e_t form, no γ multiply
            const nHxr=-(F.dezdyr+F.eyr), nHxi=-(F.dezdyi+F.eyi);
            const nHyr=F.exr+F.dezdxr, nHyi=F.exi+F.dezdxi;
            const hxr=nHxi/omu, hxi=-nHxr/omu;
            const hyr=nHyi/omu, hyi=-nHyr/omu;

            HIre += GL8w[qi]*(hxr*sdx + hyr*sdy);
            HIim += GL8w[qi]*(hxi*sdx + hyi*sdy);
        }
    }
    if (nMissed > 0) globalThis.__TRI_DEBUG__ && console.log(`  I contour: ${nFound} found, ${nMissed} missed`);

    // For half-domain symmetry, scale area integrals (P, Wm, Pdiel) to full domain
    // and line integrals (H·dl current) by sym (half-contour captures half the current)
    HIre *= sym; HIim *= sym;
    const RImag = Math.sqrt(HIre*HIre + HIim*HIim);
    Wm *= sym;
    const Pdiel_scaled = Pdiel * sym;
    const P = Math.abs(Ppoynting) * sym;
    const Zpi = RImag > 1e-30 ? 2 * P / (RImag * RImag) : NaN;
    const Zpv = P > 1e-30 ? (Vmag * Vmag) / (2 * P) : NaN;
    const RId2 = HIre*HIre + HIim*HIim;
    const ZviRe = RId2 > 1e-30 ? (Vre*HIre + Vim*HIim) / RId2 : NaN;
    const ZviIm = RId2 > 1e-30 ? (Vim*HIre - Vre*HIim) / RId2 : NaN;

    // Element-interior contour: evaluate H at each ring-triangle centroid (guaranteed
    // inside the element — no inter-element discontinuities in the evaluation).
    const elemSegs = buildElementContour(mesh, fm, condRect, distFrac);
    let ImeshRe = 0, ImeshIm = 0;
    if (elemSegs) {
        for (const seg of elemSegs) {
            const F = evalFieldsAtPoint(seg.tri, seg.xm, seg.ym, mesh, fm, vecRe, vecIm);
            // H = -(1/(jωμ₀))(∇ez + ẑ×e_t)
            const nHxr=-(F.dezdyr+F.eyr), nHxi=-(F.dezdyi+F.eyi);
            const nHyr=F.exr+F.dezdxr, nHyi=F.exi+F.dezdxi;
            const hxr=nHxi/omu, hxi=-nHxr/omu;
            const hyr=nHyi/omu, hyi=-nHyr/omu;
            ImeshRe += hxr*seg.dx + hyr*seg.dy;
            ImeshIm += hxi*seg.dx + hyi*seg.dy;
        }
    }
    // In half-domain symmetry, element contour captures half the current
    ImeshRe *= sym; ImeshIm *= sym;
    const ImeshMag = Math.sqrt(ImeshRe*ImeshRe + ImeshIm*ImeshIm);
    const ZpiMesh = ImeshMag > 1e-30 ? 2 * P / (ImeshMag * ImeshMag) : NaN;

    // Compute ring average distance for diagnostics
    let ringAvgDist = 0;
    if (elemSegs) {
        for (const seg of elemSegs) {
            const cx = Math.max(condRect.xmin, Math.min(condRect.xmax, seg.xm));
            const cy = Math.max(condRect.ymin, Math.min(condRect.ymax, seg.ym));
            ringAvgDist += Math.sqrt((seg.xm-cx)**2 + (seg.ym-cy)**2);
        }
        ringAvgDist /= elemSegs.length;
    }

    // Domain boundary Z_PI via Ampere's law: I = (jω/γ) · ∮ε(ẑ×Et)·dl
    //
    // From ∇×H = jωεE (Ampere) with e^{-γz}:
    //   γHt = -∇_tHz + jωε(ẑ×Et)
    // For a closed contour, ∮∇Hz·dl = 0 (Hz single-valued), so:
    //   γ · ∮Ht·dl = jω · ∮ε(ẑ×Et)·dl
    //   I = (jω/γ) · ∮ε(ẑ×Et)·dl
    //
    // (ẑ×Et)·dl = -Ey·dx + Ex·dy = normal Et flux through the contour.
    // On PEC boundary edges: only one adjacent triangle, no edge-crossing
    // discontinuity. Uses the Nedelec polynomial's normal Et (well-defined).

    // Find domain boundary edges and their adjacent triangles
    const edgeToTriListB = new Array(nEdges);
    for (let e = 0; e < nEdges; e++) edgeToTriListB[e] = [];
    for (let t = 0; t < nTris; t++) {
        for (let le = 0; le < 3; le++) edgeToTriListB[triEdges[3*t+le]].push(t);
    }

    const GL3p = [0.11270, 0.50000, 0.88730];
    const GL3w = [0.27778, 0.44444, 0.27778];
    const eps0_c = 8.854187817e-12;

    // PMC symmetry plane: skip these edges (Ht·dl ≡ 0 analytically, but FEM
    // doesn't enforce this exactly — spurious flux from PMC edges adds error)
    const pmcX = sym > 1 ? (condRect.xmin_domain || 0) : null;

    let ampFluxRe = 0, ampFluxIm = 0;
    for (let e = 0; e < nEdges; e++) {
        if (edgeToTriListB[e].length !== 1) continue; // only boundary edges
        if (fm.isCondEdge && fm.isCondEdge[e]) continue; // skip conductor edges

        const n0 = edges[2*e], n1 = edges[2*e+1];
        const x0 = nodes[2*n0], y0 = nodes[2*n0+1];
        const x1 = nodes[2*n1], y1 = nodes[2*n1+1];

        // Skip PMC boundary edges (both vertices on symmetry plane)
        if (pmcX !== null && Math.abs(x0 - pmcX) < 1e-12 && Math.abs(x1 - pmcX) < 1e-12) continue;
        const t = edgeToTriListB[e][0];
        const eps_t = epsMap[t];
        const epsRe = eps0_c * eps_t.re, epsIm = eps0_c * eps_t.im;

        // CCW direction: interior triangle on LEFT of edge direction
        const tv0 = tris[3*t], tv1 = tris[3*t+1], tv2 = tris[3*t+2];
        let oppV = -1;
        if (tv0 !== n0 && tv0 !== n1) oppV = tv0;
        else if (tv1 !== n0 && tv1 !== n1) oppV = tv1;
        else oppV = tv2;
        const tx = x1-x0, ty = y1-y0;
        const cross = tx*(nodes[2*oppV+1]-y0) - ty*(nodes[2*oppV]-x0);
        const sign = cross > 0 ? 1 : -1;
        const dlx = sign*tx, dly = sign*ty;

        // 3-point Gauss quadrature: ∮ε(ẑ×Et)·dl = Σ ε·(-Ey·dlx + Ex·dly)
        for (let q = 0; q < 3; q++) {
            const s = GL3p[q];
            const px = x0 + s*tx, py = y0 + s*ty;
            const F = evalFieldsAtPoint(t, px, py, mesh, fm, vecRe, vecIm);

            // (ẑ×Et)·dl = -Ey·dlx + Ex·dly
            const fluxRe = -F.eyr*dlx + F.exr*dly;
            const fluxIm = -F.eyi*dlx + F.exi*dly;

            // ε · flux (complex multiply)
            ampFluxRe += GL3w[q] * (epsRe*fluxRe - epsIm*fluxIm);
            ampFluxIm += GL3w[q] * (epsRe*fluxIm + epsIm*fluxRe);
        }
    }

    // I = (jω/γ²) · ∮ε(ẑ×e_t)·dl  (e_t = γ·Et, so extra 1/γ vs physical formula)
    // jω · ampFlux = -ω·ampFlux_im + j·ω·ampFlux_re
    const jwfRe = -omega * ampFluxIm, jwfIm = omega * ampFluxRe;
    // divide by γ² (= gamma2Re + j·gamma2Im)
    const g2m = gamma2Re*gamma2Re + gamma2Im*gamma2Im;
    const IlineRe = (jwfRe*gamma2Re + jwfIm*gamma2Im) / g2m;
    const IlineIm = (jwfIm*gamma2Re - jwfRe*gamma2Im) / g2m;
    const IlineMagHalf = Math.sqrt(IlineRe**2 + IlineIm**2);
    // In half-domain symmetry, boundary contour only captures half the current
    const IlineMag = IlineMagHalf * sym;
    const ZpiLine = IlineMag > 1e-30 ? 2 * P / (IlineMag * IlineMag) : NaN;

    // I+/I- separation for differential pair modes
    // For common mode: all contributions same sign → I+ = total, I- ≈ 0
    // For differential mode: contributions cancel → net ≈ 0, but max(I+,|I-|) = I_trace
    let ampPlusRe = 0, ampPlusIm = 0, ampMinusRe = 0, ampMinusIm = 0;
    for (let e = 0; e < nEdges; e++) {
        if (edgeToTriListB[e].length !== 1) continue;
        if (fm.isCondEdge && fm.isCondEdge[e]) continue;

        const n0 = edges[2*e], n1 = edges[2*e+1];
        const x0 = nodes[2*n0], y0 = nodes[2*n0+1];
        const x1 = nodes[2*n1], y1 = nodes[2*n1+1];

        // Skip PMC boundary edges
        if (pmcX !== null && Math.abs(x0 - pmcX) < 1e-12 && Math.abs(x1 - pmcX) < 1e-12) continue;
        const t = edgeToTriListB[e][0];
        const eps_t = epsMap[t];
        const epsRe2 = eps0_c * eps_t.re, epsIm2 = eps0_c * eps_t.im;

        const tv0 = tris[3*t], tv1 = tris[3*t+1], tv2 = tris[3*t+2];
        let oppV = -1;
        if (tv0 !== n0 && tv0 !== n1) oppV = tv0;
        else if (tv1 !== n0 && tv1 !== n1) oppV = tv1;
        else oppV = tv2;
        const tx = x1-x0, ty = y1-y0;
        const cross = tx*(nodes[2*oppV+1]-y0) - ty*(nodes[2*oppV]-x0);
        const sign = cross > 0 ? 1 : -1;
        const dlx = sign*tx, dly = sign*ty;

        let efRe = 0, efIm = 0;
        for (let q = 0; q < 3; q++) {
            const s = GL3p[q];
            const px = x0 + s*tx, py = y0 + s*ty;
            const F = evalFieldsAtPoint(t, px, py, mesh, fm, vecRe, vecIm);
            const fluxRe = -F.eyr*dlx + F.exr*dly;
            const fluxIm = -F.eyi*dlx + F.exi*dly;
            efRe += GL3w[q] * (epsRe2*fluxRe - epsIm2*fluxIm);
            efIm += GL3w[q] * (epsRe2*fluxIm + epsIm2*fluxRe);
        }
        if (efRe > 0) { ampPlusRe += efRe; ampPlusIm += efIm; }
        else { ampMinusRe += efRe; ampMinusIm += efIm; }
    }

    // Convert I+ and I- e_t flux to current via jω/γ²
    function fluxToI(fRe, fIm) {
        const jRe = -omega * fIm, jIm = omega * fRe;
        return { re: (jRe*gamma2Re + jIm*gamma2Im)/g2m, im: (jIm*gamma2Re - jRe*gamma2Im)/g2m };
    }
    const Iplus = fluxToI(ampPlusRe, ampPlusIm);
    const Iminus = fluxToI(ampMinusRe, ampMinusIm);
    const IplusMag = Math.sqrt(Iplus.re**2 + Iplus.im**2) * sym;
    const IminusMag = Math.sqrt(Iminus.re**2 + Iminus.im**2) * sym;
    const IavgMag = (IplusMag + IminusMag) / 2;
    const ImaxMag = Math.max(IplusMag, IminusMag);
    const ZpiLineAvg = IavgMag > 1e-30 ? 2 * P / (IavgMag * IavgMag) : NaN;

    // Dielectric attenuation from perturbation integral:
    // α_d = P_diel / (2P) where P_diel = (ωε₀/2)∫ε″|Et|² dA
    // This is independent of eigenvalue/shift artifacts and is zero when tand=0.
    const alpha_d = P > 1e-30 ? Pdiel_scaled / (2 * P) : 0;

    return { P, Vmag, Imag: RImag, Vre, Vim, Zpi, Zpv, ZviRe, ZviIm, gamma,
             ImeshMag, ZpiMesh, ringAvgDist, IlineMag, ZpiLine, Wm,
             IplusMag, IminusMag, ZpiLineAvg, alpha_d };
}

// --- Build per-triangle epsilon map ---

export function buildTriEpsMap(mesh, h_sub, eps_r, tand) {
    const epsMap = new Array(mesh.nTris);
    for (let t = 0; t < mesh.nTris; t++) {
        const v0 = mesh.tris[3*t], v1 = mesh.tris[3*t+1], v2 = mesh.tris[3*t+2];
        const yc = (mesh.nodes[2*v0+1] + mesh.nodes[2*v1+1] + mesh.nodes[2*v2+1]) / 3;
        epsMap[t] = yc < h_sub ? { re: eps_r, im: -eps_r * tand } : { re: 1.0, im: 0.0 };
    }
    return epsMap;
}

// --- Compute per-triangle field arrays for visualization ---

export function computeTriFieldArrays(mesh, fm, vecRe, vecIm, gamma2Re, gamma2Im, f, enrichment) {
    const { nodes, tris, nTris } = mesh;
    const omega = 2 * Math.PI * f;
    const mu0 = 4 * Math.PI * 1e-7;
    const omu = omega * mu0;

    let gamma = csqrt(gamma2Re, gamma2Im);
    if (gamma.im < 0) gamma = { re: -gamma.re, im: -gamma.im };
    const gMag2 = gamma.re*gamma.re + gamma.im*gamma.im;

    const ExRe = new Float64Array(nTris);
    const ExIm = new Float64Array(nTris);
    const EyRe = new Float64Array(nTris);
    const EyIm = new Float64Array(nTris);
    const HxRe = new Float64Array(nTris);
    const HxIm = new Float64Array(nTris);
    const HyRe = new Float64Array(nTris);
    const HyIm = new Float64Array(nTris);
    const SzArr = new Float64Array(nTris);

    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3*t], v1 = tris[3*t+1], v2 = tris[3*t+2];
        const cx = (nodes[2*v0] + nodes[2*v1] + nodes[2*v2]) / 3;
        const cy = (nodes[2*v0+1] + nodes[2*v1+1] + nodes[2*v2+1]) / 3;

        const F = evalFieldsAtPoint(t, cx, cy, mesh, fm, vecRe, vecIm, enrichment);

        // evalFieldsAtPoint returns e_t (= γ·Et), convert to physical Et for display
        const igr = gamma.re/gMag2, igi = -gamma.im/gMag2;
        ExRe[t] = F.exr*igr - F.exi*igi;
        ExIm[t] = F.exr*igi + F.exi*igr;
        EyRe[t] = F.eyr*igr - F.eyi*igi;
        EyIm[t] = F.eyr*igi + F.eyi*igr;

        // H = -(1/(jωμ₀))(∇ez + ẑ×e_t) — e_t form, no γ multiply
        const nHxr = -(F.dezdyr + F.eyr), nHxi = -(F.dezdyi + F.eyi);
        const nHyr = F.exr + F.dezdxr, nHyi = F.exi + F.dezdxi;
        const hxr = nHxi/omu, hxi = -nHxr/omu;
        const hyr = nHyi/omu, hyi = -nHyr/omu;

        HxRe[t] = hxr; HxIm[t] = hxi;
        HyRe[t] = hyr; HyIm[t] = hyi;

        // Sz = 0.5 * Re(Et*Hy* - Et*Hx*) using physical Et
        SzArr[t] = 0.5 * (ExRe[t]*hyr + ExIm[t]*hyi - EyRe[t]*hxr - EyIm[t]*hxi);
    }

    return { ExRe, ExIm, EyRe, EyIm, HxRe, HxIm, HyRe, HyIm, SzArr };
}
