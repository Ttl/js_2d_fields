// Resample a triangular-mesh static FEM solution (P2 scalar potential) onto a
// regular grid, producing the SAME { x, y, V[ny][nx], Ex[ny][nx], Ey[ny][nx] }
// shape the rectilinear FDM solver exposes — so plot.js (heatmap/contour) and
// streamlines.js work unchanged on the triangular backend's output.
//
// We resample the static potential φ and E = −∇φ (real-valued, always available,
// and the same quasi-TEM field the FDM plots), not the complex full-wave
// eigenvector. The triangle mesh itself is also returned for the mesh overlay.

import { triCoefficients, lv, le, lvGrad, leGrad } from './tri_fem.js';
import { evalFieldsAtPoint } from './tri_ms_solver.js';

const EDGE_VERTS = [[0, 1], [1, 2], [2, 0]];

// Build the regular sample grid for `domain`, honoring an optional caller-supplied
// grid (the FDM mesher's graded grid). Shared by resampleStatic / resampleModeField.
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
    function locate(x, y) {
        if (x < xmin || x > xmax || y < ymin || y > ymax) return -1;
        const list = buckets[bx(x) * nb + by(y)];
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

// Evaluate φ and ∇φ of a P2 static solution at point (x,y) in triangle t.
function evalPhi(phi, mesh, coeff, t, x, y) {
    const { tris, triEdges } = mesh;
    const v = [tris[3 * t], tris[3 * t + 1], tris[3 * t + 2]];
    let val = 0, gx = 0, gy = 0;
    for (let k = 0; k < 3; k++) {
        const pv = phi.phiVertex[v[k]];
        val += pv * lv(coeff, k, x, y);
        const g = lvGrad(coeff, k, x, y);
        gx += pv * g[0]; gy += pv * g[1];
    }
    for (let k = 0; k < 3; k++) {
        const pe = phi.phiEdge[triEdges[3 * t + k]];
        const [i, j] = EDGE_VERTS[k];
        val += pe * le(coeff, i, j, x, y);
        const g = leGrad(coeff, i, j, x, y);
        gx += pe * g[0]; gy += pe * g[1];
    }
    return { V: val, Ex: -gx, Ey: -gy };
}

// Resample a static solution onto a regular grid spanning `domain`.
// Returns { x:Float64Array(nx), y:Float64Array(ny), V, Ex, Ey } with V/Ex/Ey as
// [ny][nx] arrays (row index = y, matching plot.js / streamlines.js). Points
// outside the mesh (inside conductors, below ground) are left at 0.
export function resampleStatic(mesh, phi, domain, opts = {}) {
    // Sample on the caller-supplied grid (the FDM mesher's graded grid) when given, so
    // the contour plot resolves thin conductors / the ground surface exactly like the
    // rectilinear backend; otherwise a uniform grid spanning the domain.
    const { x, y } = makeGrid(domain, opts);
    const nx = x.length, ny = y.length;

    // For a half-domain (symmetry) solve, mirror x<0 from the meshed x>0 half.
    // parity: 'even' → V even; 'odd' → V odd. E-field parity follows automatically
    // from differentiating the mirrored V grid below.
    const parity = opts.parity || null;
    const { locate, coeffOf } = buildLocator(mesh);
    // Sample only V (the P2 potential — continuous, and conductor interiors carry
    // their exact Dirichlet potential since they are meshed), then compute
    // E = −∇V with the SAME non-uniform central-difference stencil the FDM
    // backend uses for its plots (compute_fields in field_solver.js). This makes
    // the two backends' field plots share the differentiation and rendering
    // characteristics by construction. Evaluating ∇φ per P2 element instead
    // (tried) leaves the element-boundary gradient discontinuities in the data —
    // faceted/kinked |E| contours on coarse elements, and a hard one-sided jump
    // along material-interface rows where the FDM's centered stencil blends the
    // two sides.
    const V = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ex = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ey = Array.from({ length: ny }, () => new Float64Array(nx));
    // Deterministic side selection for on-edge samples (grid rows sit exactly on
    // mesh lines): nudge the LOCATE query by a tiny NE bias; V is continuous so
    // the evaluated value is side-independent, this only makes degenerate
    // on-vertex lookups deterministic.
    const eps = 1e-9 * Math.hypot(
        (domain.x_max ?? 0) - (domain.x_min ?? 0),
        (domain.y_max ?? 0) - (domain.y_min ?? 0)) || 1e-15;
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
            const r = evalPhi(phi, mesh, coeffOf(t), t, qx, y[j]);
            V[j][i] = sV * r.V;
        }
    }
    // E = −∇V by non-uniform central differences (FDM compute_fields stencil).
    // Boundary rows/columns stay 0, exactly like the FDM plot arrays.
    for (let j = 1; j < ny - 1; j++) {
        const dyd = y[j] - y[j - 1], dyu = y[j + 1] - y[j];
        for (let i = 1; i < nx - 1; i++) {
            const dxl = x[i] - x[i - 1], dxr = x[i + 1] - x[i];
            Ex[j][i] = -(
                (dxl / (dxr * (dxl + dxr))) * V[j][i + 1] +
                ((dxr - dxl) / (dxl * dxr)) * V[j][i] -
                (dxr / (dxl * (dxl + dxr))) * V[j][i - 1]
            );
            Ey[j][i] = -(
                (dyd / (dyu * (dyd + dyu))) * V[j + 1][i] +
                ((dyu - dyd) / (dyd * dyu)) * V[j][i] -
                (dyu / (dyd * (dyd + dyu))) * V[j - 1][i]
            );
        }
    }
    return { x, y, V, Ex, Ey };
}

// Recover a CONTINUOUS nodal transverse field by area-weighted averaging of the
// per-element Nedelec field at each node. The Lee–Jin edge-element field e_t is
// DISCONTINUOUS across triangle edges, so sampling it directly at grid points gives
// faceted/discontinuous artifacts (worst near the trace, where the mesh is fine and the
// field is strong). Averaging the four complex components (Ex, Ey × re, im) to nodes and
// interpolating linearly within each triangle yields a smooth field — the same recovery
// resampleStatic uses for the static E-field.
function recoverNodalModeField(mesh, fm, vRe, vIm) {
    const { nodes, tris, nTris, nNodes } = mesh;
    const exr = new Float64Array(nNodes), exi = new Float64Array(nNodes);
    const eyr = new Float64Array(nNodes), eyi = new Float64Array(nNodes);
    const w = new Float64Array(nNodes);
    for (let t = 0; t < nTris; t++) {
        const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
        const ax = nodes[2*v0], ay = nodes[2*v0+1], bx = nodes[2*v1], by = nodes[2*v1+1], cx = nodes[2*v2], cy = nodes[2*v2+1];
        const area = Math.abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) / 2;
        if (!(area > 0)) continue;
        for (const vk of [v0, v1, v2]) {
            const f = evalFieldsAtPoint(t, nodes[2 * vk], nodes[2 * vk + 1], mesh, fm, vRe, vIm);
            exr[vk] += area * f.exr; exi[vk] += area * f.exi;
            eyr[vk] += area * f.eyr; eyi[vk] += area * f.eyi; w[vk] += area;
        }
    }
    for (let i = 0; i < nNodes; i++) {
        if (w[i] > 0) { exr[i] /= w[i]; exi[i] /= w[i]; eyr[i] /= w[i]; eyi[i] /= w[i]; }
    }
    return { exr, exi, eyr, eyi };
}

// Resample a full-wave eigenmode's transverse E-field onto a regular grid.
// The Lee–Jin eigenvector stores e_t = γ·Et, so the spatial pattern of |e_t| is the
// transverse field pattern (the constant γ only scales the whole map), which is what a
// mode plot shows. Returns { x, y, E[ny][nx], Ex[ny][nx], Ey[ny][nx] }: E is the
// transverse magnitude |Et| (for the heatmap); Ex/Ey are the real-part components (for
// an optional quiver/streamline overlay). Points outside the mesh are left at 0.
export function resampleModeField(mesh, fm, vRe, vIm, domain, opts = {}) {
    const { x, y } = makeGrid(domain, opts);
    const nx = x.length, ny = y.length;
    const { locate, coeffOf } = buildLocator(mesh);
    const nodal = recoverNodalModeField(mesh, fm, vRe, vIm);   // continuous field
    const { tris } = mesh;
    const E = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ex = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ey = Array.from({ length: ny }, () => new Float64Array(nx));
    for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
            const t = locate(x[i], y[j]);
            if (t < 0) continue;
            const coeff = coeffOf(t);
            const l0 = coeff[0][0] + coeff[0][1] * x[i] + coeff[0][2] * y[j];
            const l1 = coeff[1][0] + coeff[1][1] * x[i] + coeff[1][2] * y[j];
            const l2 = coeff[2][0] + coeff[2][1] * x[i] + coeff[2][2] * y[j];
            const v0 = tris[3 * t], v1 = tris[3 * t + 1], v2 = tris[3 * t + 2];
            // Barycentric interpolation of the recovered nodal field → continuous.
            const exr = l0 * nodal.exr[v0] + l1 * nodal.exr[v1] + l2 * nodal.exr[v2];
            const exi = l0 * nodal.exi[v0] + l1 * nodal.exi[v1] + l2 * nodal.exi[v2];
            const eyr = l0 * nodal.eyr[v0] + l1 * nodal.eyr[v1] + l2 * nodal.eyr[v2];
            const eyi = l0 * nodal.eyi[v0] + l1 * nodal.eyi[v1] + l2 * nodal.eyi[v2];
            E[j][i] = Math.hypot(Math.hypot(exr, exi), Math.hypot(eyr, eyi));
            Ex[j][i] = exr; Ey[j][i] = eyr;
        }
    }
    return { x, y, E, Ex, Ey };
}
