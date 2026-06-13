// Resample a triangular-mesh static FEM solution (P2 scalar potential) onto a
// regular grid, producing the SAME { x, y, V[ny][nx], Ex[ny][nx], Ey[ny][nx] }
// shape the rectilinear FDM solver exposes — so plot.js (heatmap/contour) and
// streamlines.js work unchanged on the triangular backend's output.
//
// We resample the static potential φ and E = −∇φ (real-valued, always available,
// and the same quasi-TEM field the FDM plots), not the complex full-wave
// eigenvector. The triangle mesh itself is also returned for the mesh overlay.

import { triCoefficients, lv, le, lvGrad, leGrad } from './tri_fem.js';

const EDGE_VERTS = [[0, 1], [1, 2], [2, 0]];

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
    const { x_min, x_max, y_min, y_max } = domain;
    const aspect = (x_max - x_min) / (y_max - y_min || 1);
    const nLong = opts.resolution || 240;
    const nx = aspect >= 1 ? nLong : Math.max(40, Math.round(nLong * aspect));
    const ny = aspect >= 1 ? Math.max(40, Math.round(nLong / aspect)) : nLong;

    const x = new Float64Array(nx), y = new Float64Array(ny);
    for (let i = 0; i < nx; i++) x[i] = x_min + (x_max - x_min) * i / (nx - 1);
    for (let j = 0; j < ny; j++) y[j] = y_min + (y_max - y_min) * j / (ny - 1);

    // For a half-domain (symmetry) solve, mirror x<0 from the meshed x>0 half.
    // parity: 'even' → V even, Ex odd, Ey even; 'odd' → V odd, Ex even, Ey odd.
    const parity = opts.parity || null;
    const { locate, coeffOf } = buildLocator(mesh);
    const V = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ex = Array.from({ length: ny }, () => new Float64Array(nx));
    const Ey = Array.from({ length: ny }, () => new Float64Array(nx));
    for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
            let qx = x[i], sV = 1, sEx = 1, sEy = 1;
            if (parity && qx < 0) {
                qx = -qx;
                if (parity === 'odd') { sV = -1; sEy = -1; }
                else { sEx = -1; }
            }
            const t = locate(qx, y[j]);
            if (t < 0) continue;
            const r = evalPhi(phi, mesh, coeffOf(t), t, qx, y[j]);
            V[j][i] = sV * r.V; Ex[j][i] = sEx * r.Ex; Ey[j][i] = sEy * r.Ey;
        }
    }
    return { x, y, V, Ex, Ey };
}
