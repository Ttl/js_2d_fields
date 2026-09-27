// Exact DC resistance and inductance of a line of rectangular conductors in free space.
//
// At DC the current density is uniform in the metal of each net, J = sigma E0, and
// the return current divides over the grounds by conductance. Its magnetic energy
// follows from the 2D log kernel, whose rectangle by rectangle integral is
// closed-form, so no field mesh is needed:
//   L_dc = -(mu0 / 2 pi) sum_ij J_i J_j int int ln|r - r'| dA dA'.
// Grounds of unlimited width (the domain walls and the grounds reaching an open
// domain edge, see wall_grounds.js) are ideal returns: they carry the whole return
// current with no resistance, and the other grounds carry none of it. The current
// keeps the distribution K of a perfect conductor, uniform through the ground
// thickness d behind each face, which adds the slab internal inductance
// mu0 d / 3 int K^2 dl over the faces inside the domain: the value both backends
// reach below the skin transition. The perfect-conductor limit L_pec of all conductors comes from a
// Galerkin boundary-element solve of the surface currents, and L_dc - L_pec is the
// DC internal inductance. Both are free-space values, independent of any domain
// truncation.
import { platingCoreOf, platedThrough } from './shapes.js';

const MU0 = 4e-7 * Math.PI;

// F(u, v) with d^4 F / du^2 dv^2 = ln sqrt(u^2 + v^2).
function F4(u, v) {
    const u2 = u * u, v2 = v * v, r2 = u2 + v2;
    let s = -25 / 48 * u2 * v2;
    if (r2 > 0) s += (6 * u2 * v2 - u2 * u2 - v2 * v2) / 48 * Math.log(r2);
    if (u !== 0 && v !== 0) s += (u2 * u * v * Math.atan(v / u) + u * v2 * v * Math.atan(u / v)) / 6;
    return s;
}

// int_A int_B ln|r - r'| over two axis-aligned rects {x0, x1, y0, y1}.
export function rectLogIntegral(a, b) {
    const us = [a.x1 - b.x0, a.x1 - b.x1, a.x0 - b.x0, a.x0 - b.x1];
    const vs = [a.y1 - b.y0, a.y1 - b.y1, a.y0 - b.y0, a.y0 - b.y1];
    const sg = [1, -1, -1, 1];
    let s = 0;
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) s += sg[i] * sg[j] * F4(us[i], vs[j]);
    return s;
}

// G(u, v) with d^2 G / du dv = ln sqrt(u^2 + v^2).
function F2(u, v) {
    const r2 = u * u + v * v;
    let s = 0;
    if (r2 > 0) s += u * v * (Math.log(r2) - 3);
    if (u !== 0 && v !== 0) s += u * u * Math.atan(v / u) + v * v * Math.atan(u / v);
    return s / 2;
}

// int_A ln|p - r'| dA' over an axis-aligned rect {x0, x1, y0, y1}.
function rectPointLogIntegral(px, py, a) {
    return F2(px - a.x0, py - a.y0) - F2(px - a.x0, py - a.y1) - F2(px - a.x1, py - a.y0) + F2(px - a.x1, py - a.y1);
}

// int f dl along panel p by Gauss quadrature, 8 points when the source at centre
// offset (dx, dy) with size `size` is close.
function panelQuadrature(p, dx, dy, size, f) {
    const near = Math.hypot(dx, dy) < 2 * (p.len + size);
    let s = 0;
    for (const [t, w] of GAUSS[near ? 8 : 2]) {
        s += w * f(p.x0 + (p.x1 - p.x0) * (t + 1) / 2, p.y0 + (p.y1 - p.y0) * (t + 1) / 2);
    }
    return s * p.len / 2;
}

// int_A int_p ln|r - r'| dA dl over a rect and a panel, Gauss on the panel.
function rectPanelLogIntegral(a, p) {
    const cx = (a.x0 + a.x1) / 2, cy = (a.y0 + a.y1) / 2, size = Math.max(a.x1 - a.x0, a.y1 - a.y0);
    return panelQuadrature(p, (p.x0 + p.x1) / 2 - cx, (p.y0 + p.y1) / 2 - cy, size, (x, y) => rectPointLogIntegral(x, y, a));
}

// int ln|p - r'| dl' over the segment from q0 to q1.
function segLogIntegral(px, py, x0, y0, x1, y1) {
    const L = Math.hypot(x1 - x0, y1 - y0);
    const tx = (x1 - x0) / L, ty = (y1 - y0) / L;
    const s0 = (px - x0) * tx + (py - y0) * ty;
    const h = Math.abs(-(px - x0) * ty + (py - y0) * tx);
    const G = s => {
        const d = s - s0, r2 = d * d + h * h;
        return (r2 > 0 ? 0.5 * d * Math.log(r2) : 0) - d + (h > 0 ? h * Math.atan(d / h) : 0);
    };
    return G(L) - G(0);
}

const GAUSS = {
    2: [[-0.5773502691896257, 1], [0.5773502691896257, 1]],
    8: [[-0.9602898564975363, 0.1012285362903763], [-0.7966664774136267, 0.2223810344533745],
        [-0.5255324099163290, 0.3137066458778873], [-0.1834346424956498, 0.3626837833783620],
        [0.1834346424956498, 0.3626837833783620], [0.5255324099163290, 0.3137066458778873],
        [0.7966664774136267, 0.2223810344533745], [0.9602898564975363, 0.1012285362903763]],
};

// int_p int_q ln|r - r'| dl dl' for two panels {x0, y0, x1, y1, len}: closed-form
// for a panel with itself, the inner integral closed-form and Gauss quadrature on
// the outer one otherwise (8 points when they are close).
function panelLogIntegral(p, q) {
    if (p === q) return p.len * p.len * (Math.log(p.len) - 1.5);
    return panelQuadrature(p, (p.x0 + p.x1 - q.x0 - q.x1) / 2, (p.y0 + p.y1 - q.y0 - q.y1) / 2, q.len,
        (x, y) => segLogIntegral(x, y, q.x0, q.y0, q.x1, q.y1));
}

// Dense solve of A x = b with partial pivoting (A is overwritten).
function solveDense(A, b, n) {
    for (let k = 0; k < n; k++) {
        let piv = k, best = Math.abs(A[k * n + k]);
        for (let i = k + 1; i < n; i++) { const v = Math.abs(A[i * n + k]); if (v > best) { best = v; piv = i; } }
        if (!(best > 0)) return null;
        if (piv !== k) {
            for (let j = 0; j < n; j++) { const t = A[k * n + j]; A[k * n + j] = A[piv * n + j]; A[piv * n + j] = t; }
            const t = b[k]; b[k] = b[piv]; b[piv] = t;
        }
        const d = A[k * n + k];
        for (let i = k + 1; i < n; i++) {
            const m = A[i * n + k] / d;
            if (m === 0) continue;
            for (let j = k + 1; j < n; j++) A[i * n + j] -= m * A[k * n + j];
            b[i] -= m * b[k];
        }
    }
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
        let s = b[i];
        for (let j = i + 1; j < n; j++) s -= A[i * n + j] * x[j];
        x[i] = s / A[i * n + i];
    }
    return x;
}

// (x, y) strictly inside the cell {x0, x1, y0, y1}.
const inCell = (c, x, y) => x > c.x0 && x < c.x1 && y > c.y0 && y < c.y1;

// Net of a conductor: 0 ground, 1 positive trace, 2 negative trace.
const netOf = c => (!c.is_signal ? 0 : c.polarity < 0 ? 2 : 1);

// Disjoint metal cells of the conductors on the grid of their edge lines, merged into
// maximal runs: [{ x0, x1, y0, y1, net, sigma, wall, d }] with the later conductor
// (and a plating core over its plating layer) winning where they overlap. `wall` marks
// the ideal grounds, `d` their slab thickness: the thin side of the conductor, null
// for an absorbed wall (the thickness of the merged cell, stacked slabs added up).
// null when a conductor is a shape.
function metalCells(conductors, sigmaDefault, ideal, absorbed) {
    const pieces = [];
    for (const [ci, c] of conductors.entries()) {
        if (c.shape) return null;
        const net = netOf(c), bulk = c.sigma > 0 ? c.sigma : sigmaDefault, wall = ideal.has(ci);
        const d = wall && !absorbed.has(ci) ? Math.min(c.x_max - c.x_min, c.y_max - c.y_min) : null;
        const box = { x0: c.x_min, x1: c.x_max, y0: c.y_min, y1: c.y_max, net, wall, d };
        const pl = c.plating;
        if (pl && platedThrough(c)) { pieces.push({ ...box, sigma: pl.sigma }); continue; }
        const core = platingCoreOf(c);
        if (core && core.rects) {
            pieces.push({ ...box, sigma: pl.sigma });
            for (const r of core.rects) pieces.push({ x0: r.xmin, x1: r.xmax, y0: r.ymin, y1: r.ymax, net, wall, d, sigma: bulk });
        } else pieces.push({ ...box, sigma: bulk });
    }
    if (!pieces.length) return { cells: [], X: [0], Y: [0] };
    let span = 0;
    for (const p of pieces) span = Math.max(span, Math.abs(p.x0), Math.abs(p.x1), Math.abs(p.y0), Math.abs(p.y1));
    const tol = 1e-9 * span;
    const lines = vals => {
        const s = [...vals].sort((a, b) => a - b), out = [];
        for (const v of s) if (!out.length || v - out[out.length - 1] > tol) out.push(v);
        return out;
    };
    const X = lines(pieces.flatMap(p => [p.x0, p.x1])), Y = lines(pieces.flatMap(p => [p.y0, p.y1]));
    const at = (x, y) => {
        for (let k = pieces.length - 1; k >= 0; k--) if (inCell(pieces[k], x, y)) return pieces[k];
        return null;
    };
    // Runs along x in each band, then bands with the same runs stacked.
    const cells = [];
    let open = new Map();
    for (let b = 0; b + 1 < Y.length; b++) {
        const y = (Y[b] + Y[b + 1]) / 2, runs = [];
        for (let a = 0; a + 1 < X.length; a++) {
            const p = at((X[a] + X[a + 1]) / 2, y);
            if (!p) continue;
            const last = runs[runs.length - 1];
            if (last && last.x1 === X[a] && last.net === p.net && last.sigma === p.sigma && last.wall === p.wall
                && last.d === p.d) last.x1 = X[a + 1];
            else runs.push({ x0: X[a], x1: X[a + 1], net: p.net, sigma: p.sigma, wall: p.wall, d: p.d });
        }
        const next = new Map();
        for (const r of runs) {
            const key = `${r.x0}|${r.x1}|${r.net}|${r.sigma}|${r.wall}|${r.d}`;
            const prev = open.get(key);
            if (prev) { prev.y1 = Y[b + 1]; next.set(key, prev); }
            else { const c = { ...r, y0: Y[b], y1: Y[b + 1] }; cells.push(c); next.set(key, c); }
        }
        open = next;
    }
    return { cells, X, Y };
}

// Corners of the outline of each net's metal in `cells`, as boundary segments
// { net, dir, at, a, b } merged along their lines.
function outlineSegments(cells) {
    const netAt = (x, y) => {
        for (const c of cells) if (inCell(c, x, y)) return c.net;
        return -1;
    };
    const raw = [];
    for (const c of cells) {
        const e = 1e-7 * Math.min(c.x1 - c.x0, c.y1 - c.y0);
        // Split each side at the edge lines of the other cells so each piece has one neighbour.
        const xs = new Set([c.x0, c.x1]), ys = new Set([c.y0, c.y1]);
        for (const o of cells) {
            if (o.x0 > c.x0 && o.x0 < c.x1) xs.add(o.x0);
            if (o.x1 > c.x0 && o.x1 < c.x1) xs.add(o.x1);
            if (o.y0 > c.y0 && o.y0 < c.y1) ys.add(o.y0);
            if (o.y1 > c.y0 && o.y1 < c.y1) ys.add(o.y1);
        }
        const X = [...xs].sort((a, b) => a - b), Y = [...ys].sort((a, b) => a - b);
        for (let i = 0; i + 1 < X.length; i++) {
            const xm = (X[i] + X[i + 1]) / 2;
            if (netAt(xm, c.y0 - e) !== c.net) raw.push({ net: c.net, dir: 'h', at: c.y0, a: X[i], b: X[i + 1] });
            if (netAt(xm, c.y1 + e) !== c.net) raw.push({ net: c.net, dir: 'h', at: c.y1, a: X[i], b: X[i + 1] });
        }
        for (let i = 0; i + 1 < Y.length; i++) {
            const ym = (Y[i] + Y[i + 1]) / 2;
            if (netAt(c.x0 - e, ym) !== c.net) raw.push({ net: c.net, dir: 'v', at: c.x0, a: Y[i], b: Y[i + 1] });
            if (netAt(c.x1 + e, ym) !== c.net) raw.push({ net: c.net, dir: 'v', at: c.x1, a: Y[i], b: Y[i + 1] });
        }
    }
    raw.sort((p, q) => (p.net - q.net) || (p.dir < q.dir ? -1 : p.dir > q.dir ? 1 : 0) || (p.at - q.at) || (p.a - q.a));
    const segs = [];
    for (const r of raw) {
        const last = segs[segs.length - 1];
        if (last && last.net === r.net && last.dir === r.dir && last.at === r.at && last.b === r.a) last.b = r.b;
        else segs.push({ ...r });
    }
    return segs;
}

// int d/3 K^2 dl over the faces of the ideal grounds inside the domain `box`, in the
// scaled units. K on a face is the tangential field just outside it, from the
// potential S of the panel currents K and the trace current densities J: on a
// perfect conductor dS/dn = 2 pi K. The field is insensitive to the ill-conditioned
// split of the panel currents between the two faces of a thin ground. The current
// behind a face on the domain edge (the back of a wall) belongs to the opposite face,
// as the domain has no field behind the wall.
function idealSlabIntegral(cells, panels, K, carrying, J, box, tol) {
    const potential = (x, y) => {
        let s = 0;
        for (let q = 0; q < panels.length; q++) s += K[q] * segLogIntegral(x, y, panels[q].x0, panels[q].y0, panels[q].x1, panels[q].y1);
        for (let i = 0; i < carrying.length; i++) s += J[i] * rectPointLogIntegral(x, y, carrying[i]);
        return s;
    };
    const cellAt = (x, y) => cells.find(c => inCell(c, x, y));
    const onBox = (x, y) => box && (Math.abs(x - box.x0) < tol || Math.abs(x - box.x1) < tol
        || Math.abs(y - box.y0) < tol || Math.abs(y - box.y1) < tol);
    // Signed surface current at (x, y) on a face with outward normal (nx, ny).
    const surfaceK = (x, y, nx, ny, e) => {
        const s0 = potential(x, y), s1 = potential(x + nx * e, y + ny * e), s2 = potential(x + 2 * nx * e, y + 2 * ny * e);
        return (-3 * s0 + 4 * s1 - s2) / (2 * e) / (2 * Math.PI);
    };
    let sum = 0;
    for (const p of panels) {
        const mx = (p.x0 + p.x1) / 2, my = (p.y0 + p.y1) / 2;
        if (onBox(mx, my)) continue;
        const horiz = p.y0 === p.y1, probe = 1e-6 * p.len;
        let nx = 0, ny = 0;
        if (horiz) ny = cellAt(mx, my + probe) ? -1 : 1; else nx = cellAt(mx + probe, my) ? -1 : 1;
        const c = cellAt(mx - nx * probe, my - ny * probe);
        if (!c) continue;
        const thin = Math.min(c.x1 - c.x0, c.y1 - c.y0), d = c.d ?? thin;
        const e = 0.02 * Math.min(p.len, thin);
        let k = surfaceK(mx, my, nx, ny, e);
        // The opposite face of the same cell on the domain edge.
        const ox = nx ? (nx > 0 ? c.x0 : c.x1) : mx, oy = ny ? (ny > 0 ? c.y0 : c.y1) : my;
        if (onBox(ox, oy)) k += surfaceK(ox, oy, -nx, -ny, e);
        sum += d / 3 * k * k * p.len;
    }
    return sum;
}

const segCorners = segs => segs.flatMap(s => (s.dir === 'h' ? [[s.a, s.at], [s.b, s.at]] : [[s.at, s.a], [s.at, s.b]]));

// Boundary panels on segments, graded towards the corners: a panel is at most
// `alpha` times its distance to the nearest corner.
function boundaryPanels(segs, corners, alpha) {
    let hmin = Infinity;
    for (const s of segs) hmin = Math.min(hmin, s.b - s.a);
    hmin *= 0.05;
    const panels = [];
    for (const s of segs) {
        const pt = t => (s.dir === 'h' ? [t, s.at] : [s.at, t]);
        const dCorner = t => {
            const [x, y] = pt(t);
            let d = Infinity;
            for (const [cx, cy] of corners) d = Math.min(d, Math.hypot(x - cx, y - cy));
            return d;
        };
        let t = s.a;
        const cuts = [t];
        while (t < s.b) {
            let h = Math.max(hmin, alpha * dCorner(t));
            // The panel may not outgrow the distance at its far end either.
            h = Math.max(hmin, Math.min(h, alpha * dCorner(Math.min(s.b, t + h))));
            t = Math.min(s.b, t + h);
            if (s.b - t < 0.5 * hmin) t = s.b;
            cuts.push(t);
        }
        for (let i = 0; i + 1 < cuts.length; i++) {
            const [x0, y0] = pt(cuts[i]), [x1, y1] = pt(cuts[i + 1]);
            panels.push({ net: s.net, x0, y0, x1, y1, len: cuts[i + 1] - cuts[i] });
        }
    }
    return panels;
}

// Panel currents K of perfect conductors with one vector potential per net, the
// nets' currents I and a known source: sum_q G_pq K_q + src_p = Lambda_net(p) len_p,
// sum over a net of K_q len_q = I_net. Returns { K, Lambda } or null.
function solvePanels(panels, nets, I, src) {
    const N = panels.length, M = N + nets.length;
    const A = new Float64Array(M * M), b = new Float64Array(M);
    for (let p = 0; p < N; p++) {
        for (let q = p; q < N; q++) {
            const g = panelLogIntegral(panels[p], panels[q]);
            A[p * M + q] = g; A[q * M + p] = g;
        }
        const k = nets.indexOf(panels[p].net);
        A[p * M + N + k] = -panels[p].len;
        A[(N + k) * M + p] = panels[p].len;
        if (src) b[p] = -src[p];
    }
    nets.forEach((k, i) => { b[N + i] = I[k]; });
    const x = solveDense(A, b, M);
    return x ? { K: x.subarray(0, N), Lambda: Array.from(x.subarray(N)) } : null;
}

// Net currents of a mode (per trace current 1) and the number of traces it drives.
function modeCurrents(mode) {
    if (mode === 'odd') return { I: [0, 1, -1], n: 2 };
    if (mode === 'even') return { I: [-2, 1, 1], n: 2 };
    if (mode === 'single') return { I: [-1, 1, 0], n: 1 };
    return null;
}

// DC line parameters of conductors [{ x_min, x_max, y_min, y_max, is_signal, polarity,
// sigma?, plating?, shape? }] for mode 'single', 'odd' or 'even' (per trace):
// { R, L, Lpec, Lint } in ohm/m and H/m, or null for shaped conductors, a missing net
// or a failed solve. `unlimited` holds the indices of the grounds of unlimited width,
// `walls` those of them absorbed into the domain walls, `box` the domain
// { x_min, x_max, y_min, y_max } (faces on its edge carry no slab term), `alpha` sets
// the boundary panel grading.
export function dcLineParameters(conductors, mode, { sigmaDefault = 5.8e7, unlimited = new Set(), walls = new Set(),
    box = null, alpha = 0.25 } = {}) {
    const mc = modeCurrents(mode);
    if (!mc) return null;
    const grid = metalCells(conductors, sigmaDefault, unlimited, walls);
    if (!grid || !grid.cells.length) return null;
    const { I, n } = mc;
    const hasWall = grid.cells.some(c => c.wall);
    // Net conductances; a driven net must have metal.
    const G = [0, 0, 0];
    for (const c of grid.cells) G[c.net] += c.sigma * (c.x1 - c.x0) * (c.y1 - c.y0);
    for (let k = 0; k < 3; k++) if (I[k] !== 0 && !(G[k] > 0)) return null;
    let R = 0;
    for (let k = 1; k < 3; k++) if (I[k] !== 0) R += I[k] * I[k] / G[k];
    if (!hasWall) R += I[0] * I[0] / G[0];
    R /= n;

    // Log kernel on coordinates scaled to O(1); the zero total current removes the scale.
    const ell = Math.max(grid.X[grid.X.length - 1] - grid.X[0], grid.Y[grid.Y.length - 1] - grid.Y[0]);
    const scCell = c => ({ ...c, x0: c.x0 / ell, x1: c.x1 / ell, y0: c.y0 / ell, y1: c.y1 / ell });
    const scPanel = p => ({ ...p, x0: p.x0 / ell, y0: p.y0 / ell, x1: p.x1 / ell, y1: p.y1 / ell, len: p.len / ell });
    const allSegs = outlineSegments(grid.cells);
    const corners = segCorners(allSegs);

    // Uniform current in the metal that carries it: every conductor without walls,
    // the traces with them.
    const carrying = grid.cells.filter(c => I[c.net] !== 0 && !(hasWall && c.net === 0)).map(scCell);
    const J = carrying.map(c => c.sigma * I[c.net] / G[c.net] * ell * ell);
    let E = 0, Eslab = 0;
    for (let i = 0; i < carrying.length; i++) {
        E += J[i] * J[i] * rectLogIntegral(carrying[i], carrying[i]);
        for (let j = i + 1; j < carrying.length; j++) E += 2 * J[i] * J[j] * rectLogIntegral(carrying[i], carrying[j]);
    }
    if (hasWall) {
        // The walls carry the return current as perfect conductors in the field of
        // the traces: E = J G_ss J + J G_sw K + Lambda I_gnd.
        const wallCells = grid.cells.filter(c => c.wall);
        const panels = boundaryPanels(outlineSegments(wallCells), corners, alpha).map(scPanel);
        const src = panels.map(p => {
            let s = 0;
            for (let i = 0; i < carrying.length; i++) s += J[i] * rectPanelLogIntegral(carrying[i], p);
            return s;
        });
        const sol = solvePanels(panels, [0], I, src);
        if (!sol) return null;
        for (let p = 0; p < panels.length; p++) E += sol.K[p] * src[p];
        E += sol.Lambda[0] * I[0];
        // Slab term: d/3 int K^2 dl is scale free in the scaled units.
        const sBox = box && { x0: box.x_min / ell, x1: box.x_max / ell, y0: box.y_min / ell, y1: box.y_max / ell };
        Eslab = idealSlabIntegral(wallCells.map(c => ({ ...scCell(c), d: c.d === null ? null : c.d / ell })),
            panels, sol.K, carrying, J, sBox, 1e-9);
    }
    const L = (-MU0 / (2 * Math.PI) * E + MU0 * Eslab) / n;

    // Perfect-conductor limit of all conductors.
    const nets = [0, 1, 2].filter(k => G[k] > 0);
    const panels = boundaryPanels(allSegs, corners, alpha).map(scPanel);
    const sol = solvePanels(panels, nets, I, null);
    if (!sol) return null;
    let Ep = 0;
    nets.forEach((k, i) => { Ep += sol.Lambda[i] * I[k]; });
    const Lpec = -MU0 / (2 * Math.PI) * Ep / n;
    if (!Number.isFinite(L) || !Number.isFinite(Lpec)) return null;
    return { R, L, Lpec, Lint: L - Lpec };
}
