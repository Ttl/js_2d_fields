// Ground slabs absorbed into the domain walls. A ground rect spanning the domain
// along one side, with nothing beyond it, is the wall: the full-wave mesh clips it
// off and the wall takes a PEC boundary condition, and both backends treat it as an
// ideal return of unlimited width at DC. Absorbing a slab moves that side of the
// domain inwards, so a slab stacked on it or a side wall standing on it is absorbed
// next.
function rectOf(o) {
    if (o.xmin !== undefined) return { xmin: o.xmin, xmax: o.xmax, ymin: o.ymin, ymax: o.ymax };
    return { xmin: o.x_min, xmax: o.x_max, ymin: o.y_min, ymax: o.y_max };
}

// Clipped domain { X0, X1, Y0, Y1 }, the PEC walls, the metal thickness of each
// (Infinity for a bare 'gnd' boundary; stacked slabs add up) and the indices of the
// absorbed conductors.
export function clipDomainWalls(domain, conductors, boundaries, tol) {
    let { x_min: X0, x_max: X1, y_min: Y0, y_max: Y1 } = domain;
    const b = boundaries || ['open', 'open', 'open', 'gnd'];
    const wallPEC = { left: b[0] === 'gnd', right: b[1] === 'gnd', top: b[2] === 'gnd', bottom: b[3] === 'gnd' };
    const wallThick = { left: Infinity, right: Infinity, top: Infinity, bottom: Infinity };
    const absorbed = new Set();
    const absorb = (side, t, ci) => {
        wallThick[side] = (Number.isFinite(wallThick[side]) ? wallThick[side] : 0) + t;
        wallPEC[side] = true;
        absorbed.add(ci);
    };
    let changed = true;
    while (changed) {
        changed = false;
        for (const [ci, c] of conductors.entries()) {
            if (c.is_signal || absorbed.has(ci)) continue;
            // Shaped grounds (e.g. a coax shield) are not full-span slabs: their
            // bounding box spans the domain but their body does not fill it.
            if (c.shape) continue;
            const r = rectOf(c);
            const touchesL = r.xmin <= X0 + tol, touchesR = r.xmax >= X1 - tol;
            const touchesB = r.ymin <= Y0 + tol, touchesT = r.ymax >= Y1 - tol;
            if (touchesB && touchesL && touchesR && r.ymax < Y1 - tol && r.ymax > Y0 + tol) { absorb('bottom', r.ymax - Y0, ci); Y0 = r.ymax; changed = true; continue; }
            if (touchesT && touchesL && touchesR && r.ymin > Y0 + tol && r.ymin < Y1 - tol) { absorb('top', Y1 - r.ymin, ci); Y1 = r.ymin; changed = true; continue; }
            if (touchesL && touchesB && touchesT && r.xmax < X1 - tol && r.xmax > X0 + tol) { absorb('left', r.xmax - X0, ci); X0 = r.xmax; changed = true; continue; }
            if (touchesR && touchesB && touchesT && r.xmin > X0 + tol && r.xmin < X1 - tol) { absorb('right', X1 - r.xmin, ci); X1 = r.xmin; changed = true; continue; }
        }
    }
    return { X0, X1, Y0, Y1, wallPEC, wallThick, absorbed };
}

// Grounds of unlimited width: the walls, and every ground rect that reaches an open
// edge of the domain, a plane cut off by the field truncation (a pour, a via fence or
// the rest of a plane beside a ground cutout). Both are ideal returns at DC. Returns
// { walls, unlimited } as sets of conductor indices.
export function unlimitedGrounds(domain, conductors, boundaries, tol) {
    const { absorbed } = clipDomainWalls(domain, conductors, boundaries, tol);
    const unlimited = new Set(absorbed);
    const b = boundaries || ['open', 'open', 'open', 'gnd'];
    const open = { left: b[0] !== 'gnd', right: b[1] !== 'gnd', top: b[2] !== 'gnd', bottom: b[3] !== 'gnd' };
    for (const [ci, c] of conductors.entries()) {
        if (c.is_signal || c.shape || unlimited.has(ci)) continue;
        const r = rectOf(c);
        if ((open.left && r.xmin <= domain.x_min + tol) || (open.right && r.xmax >= domain.x_max - tol)
            || (open.bottom && r.ymin <= domain.y_min + tol) || (open.top && r.ymax >= domain.y_max - tol)) unlimited.add(ci);
    }
    return { walls: absorbed, unlimited };
}
