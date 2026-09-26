// Custom cross-section built from axis-aligned rectangles, see custom_geometry_text.js
// for the input format. Both backends work from the conductor and dielectric lists,
// so this solver only has to resolve the domain, validate the rectangles and build
// the lists.
import { FieldSolver2D } from './field_solver.js';
import { Dielectric, Conductor, Mesher } from './mesher.js';
import { halfDomainSymmetry, isXSymmetric } from './geometry_symmetry.js';
import { parseAndEvaluate, formatErrors } from './custom_geometry_text.js';
import { bodyDistance, translateShapeX } from './shapes.js';

const BIG = 1e30;
const WALLS = ['left', 'right', 'top', 'bottom'];

// Gap between two rectangles given as {x0, x1, y0, y1}, 0 when they touch or overlap.
function rectDistance(a, b) {
    const dx = Math.max(0, a.x0 - b.x1, b.x0 - a.x1);
    const dy = Math.max(0, a.y0 - b.y1, b.y0 - a.y1);
    return Math.hypot(dx, dy);
}

// Gap between two resolved rectangles or shapes ({x0, x1, y0, y1, shape}). Shapes are
// measured on their polygons, plain rectangles as before.
function bodyGap(a, b) {
    if (!a.shape && !b.shape) return rectDistance(a, b);
    const body = o => ({ xmin: o.x0, xmax: o.x1, ymin: o.y0, ymax: o.y1, shape: o.shape || null });
    return bodyDistance(body(a), body(b));
}

function overlapArea(a, b) {
    const w = Math.min(a.x1, b.x1) - Math.max(a.x0, b.x0);
    const h = Math.min(a.y1, b.y1) - Math.max(a.y0, b.y0);
    return (w > 0 && h > 0) ? { w, h } : null;
}

const shift0 = (shift, tol) => (Math.abs(shift) > tol ? shift : 0);

// True when the intervals cover [lo, hi] without a gap larger than tol.
function intervalsCover(intervals, lo, hi, tol) {
    let reach = lo;
    for (const [a, b] of intervals.slice().sort((p, q) => p[0] - q[0])) {
        if (a > reach + tol) return false;
        reach = Math.max(reach, b);
        if (reach >= hi - tol) return true;
    }
    return reach >= hi - tol;
}

class CustomGeometrySolver extends FieldSolver2D {
    // options:
    //   text / overrides  - geometry text and parameter overrides, or
    //   geometry          - an evaluateGeometry() result
    //   plating           - plating material for rectangles with plating= when the
    //                       text has no plating statement
    //   wall_thickness    - thickness of the ground slab added behind a gnd wall
    //   sigma_cond, freq, nx, ny, rq, mesh_backend, symmetry as for the other solvers
    constructor(options) {
        super();

        const geo = options.geometry ?? parseAndEvaluate(options.text ?? '', options.overrides ?? {});
        if (geo.errors.length) throw new Error('Geometry errors:\n' + formatErrors(geo.errors));

        this.sigma_cond = options.sigma_cond ?? 5.8e7;
        this.freq = options.freq ?? 1e9;
        this.nx = options.nx ?? 300;
        this.ny = options.ny ?? 300;
        this.rq = options.rq ?? 0;
        this.mesh_backend = options.mesh_backend ?? 'rectilinear';
        this.boundaries = [...geo.bounds];
        this.has_side_gnd = (this.boundaries[0] === 'gnd' || this.boundaries[1] === 'gnd');
        this.wall_t = options.wall_thickness ?? 35e-6;
        this.geometry_params = { ...geo.params };

        const platingMaterial = geo.plating ?? (options.plating
            ? { sigma: options.plating.sigma, thickness: options.plating.thickness,
                rq: options.plating.rq, thick_corners: options.plating.thick_corners }
            : null);
        this.plating = platingMaterial;

        const rects = this._resolve_rects(geo);
        this._validate_rects(rects, platingMaterial);
        this._build_lists(rects, platingMaterial);
        // Source lines of the trapezoids and n-gons, which only the full-wave solver takes.
        this.shaped_lines = [...new Set(rects.filter(r => r.shape && !r.image).map(r => r.line))];

        const signals = this.conductors.filter(c => c.is_signal);
        const grounds = this.conductors.filter(c => !c.is_signal);
        this.is_differential = signals.some(c => c.polarity < 0);
        // Sizing hints read by the loss transition model and the triangular mesher.
        this.t = Math.min(...signals.map(c => Math.abs(c.height)));
        const finiteW = rects.filter(r => r.kind.startsWith('sig') && !r.xInf).map(r => r.x1 - r.x0);
        this.w = finiteW.length ? Math.min(...finiteW) : undefined;
        this.t_gnd = grounds.length
            ? Math.min(...grounds.map(c => Math.min(Math.abs(c.width), Math.abs(c.height)))) : null;

        // Mirror symmetry about x=0 lets the mesher build a symmetric grid, and the
        // half domain is used on top of that when the modes allow it.
        const paintOk = !this._paint_order_breaks_symmetry();
        // The grid mirrors whenever the shapes do. The half domain also needs mirrored
        // conductors to share their surface finish, or it reports one side's loss for both.
        const mirrorShape = paintOk && isXSymmetric(this.conductors, this.dielectrics, this.domain_width, { finish: false });
        const mirror = paintOk && isXSymmetric(this.conductors, this.dielectrics, this.domain_width);
        const symInfo = halfDomainSymmetry(this.conductors, this.dielectrics,
                                           this.domain_width, this.is_differential);
        const symAllowed = options.symmetry !== false && mirror;
        if (!symAllowed) this.tri_symmetry = false;
        this.sym_half = this.mesh_backend !== 'triangular' && symAllowed && symInfo.ok;
        // A signal body drawn as two rectangles that meet on the plane is cut by it
        // just like a single rectangle that crosses it.
        const planeTol = this.domain_width * 1e-6;
        const meetsPlane = signals.some(c => Math.abs(c.x_min) <= planeTol || Math.abs(c.x_max) <= planeTol);
        this._sym_signal_straddles = this.sym_half && (symInfo.straddles || meetsPlane);
        this._proximityWarn = this._broadside_proximity_note(signals);
        // Per-conductor finishes need the centred loss quadrature (see
        // calculate_conductor_loss): the default rule only balances over mirrored pairs
        // of equal finish. Geometries without a finish of their own keep the default
        // rule, so a converted fixed type solves exactly as before.
        this.centred_loss_quadrature = this._own_finish || (mirrorShape && !mirror);

        this.mesher = new Mesher(
            this.domain_width, this.domain_height,
            this.nx, this.ny,
            this.conductors, this.dielectrics,
            mirrorShape,
            -this.domain_width / 2,
            this.domain_width / 2,
            this.domain_y_min,
            this.domain_height,
            this.sym_half
        );

        this.x = null;
        this.y = null;
        this.dx = null;
        this.dy = null;
        this.mesh_generated = false;
    }

    // Resolves the domain, pins inf edges to it, adds the ground slabs behind gnd
    // walls and centres the result on x=0. Returns rectangles as
    // { kind, x, y, w, h, x0, x1, y0, y1, er, tand, faces, line, xInf }, where x, y, w, h
    // are the constructor arguments of Conductor / Dielectric and x0..y1 the bounds.
    _resolve_rects(geo) {
        const clamp = v => Math.max(-BIG, Math.min(BIG, v));
        const rects = geo.rects.map(r => ({
            kind: r.kind, src: r, line: r.line, er: r.er, tand: r.tand, thin: r.thin, faces: r.plating,
            sigma: r.sigma, rq: r.rq, platingOwn: r.platingMaterial, image: r.image, shape: r.shape || null,
            x0: clamp(r.x.min), x1: clamp(r.x.max), y0: clamp(r.y.min), y1: clamp(r.y.max),
            xInf: !Number.isFinite(r.x.min) || !Number.isFinite(r.x.max),
        }));
        const conds = rects.filter(r => r.kind !== 'diel');
        const signals = conds.filter(r => r.kind !== 'gnd');
        if (signals.length === 0) throw new Error('The geometry needs at least one signal conductor (sig+).');

        const finite = (list, lo, hi) => {
            const v = [];
            for (const r of list) {
                if (Math.abs(r[lo]) < BIG) v.push(r[lo]);
                if (Math.abs(r[hi]) < BIG) v.push(r[hi]);
            }
            return v;
        };
        const cx = finite(conds, 'x0', 'x1'), ay = finite(rects, 'y0', 'y1');
        if (cx.length === 0 || ay.length === 0) {
            throw new Error('The geometry needs at least one finite conductor edge in x and in y.');
        }
        const cLo = Math.min(...cx), cHi = Math.max(...cx);
        const yLo = Math.min(...ay), yHi = Math.max(...ay);
        const b = this.boundaries, d = geo.domain;

        // Reference length of the fringing field: distance from the signals to the
        // nearest ground rectangle or ground wall.
        const wallRects = [];
        // An auto gnd wall sits on the edge of the stack when a dielectric ends there,
        // as the ground plane of a substrate. Past a conductor it is an enclosure wall
        // at the open-boundary distance.
        const eps = (yHi - yLo) * 1e-9;
        const diels = rects.filter(r => r.kind === 'diel');
        const snapBot = b[3] === 'gnd' && diels.some(r => r.y0 <= yLo + eps);
        const snapTop = b[2] === 'gnd' && diels.some(r => r.y1 >= yHi - eps);
        const yBot = d.y_min ?? (snapBot ? yLo : null);
        const yTop = d.y_max ?? (snapTop ? yHi : null);
        if (b[3] === 'gnd' && yBot !== null) wallRects.push({ x0: -BIG, x1: BIG, y0: -BIG, y1: yBot });
        if (b[2] === 'gnd' && yTop !== null) wallRects.push({ x0: -BIG, x1: BIG, y0: yTop, y1: BIG });
        if (b[0] === 'gnd' && d.x_min !== null) wallRects.push({ x0: -BIG, x1: d.x_min, y0: -BIG, y1: BIG });
        if (b[1] === 'gnd' && d.x_max !== null) wallRects.push({ x0: d.x_max, x1: BIG, y0: -BIG, y1: BIG });
        const gnds = [...conds.filter(r => r.kind === 'gnd'), ...wallRects];
        if (gnds.length === 0) {
            throw new Error('The geometry needs a ground: a gnd rectangle or a gnd boundary.');
        }
        let href = Infinity;
        for (const s of signals) for (const g of gnds) href = Math.min(href, bodyGap(s, g));
        const sx = finite(signals, 'x0', 'x1');
        const span = (sx.length ? Math.max(...sx) - Math.min(...sx) : 0) + 2 * href;
        const sub = finite(rects.filter(r => r.kind === 'diel' && r.er > 1.001), 'y0', 'y1');
        const hsub = sub.length ? Math.max(...sub) - Math.min(...sub) : 0;
        const margin = Math.max(4 * span, 15 * href, 8 * hsub);
        if (!(margin > 0)) throw new Error('Cannot size the domain: a signal touches a ground.');

        const X0 = d.x_min ?? cLo - margin, X1 = d.x_max ?? cHi + margin;
        const Y0 = d.y_min ?? (snapBot ? yLo : yLo - margin);
        const Y1 = d.y_max ?? (snapTop ? yHi : yHi + margin);
        if (!(X1 > X0) || !(Y1 > Y0)) throw new Error('The domain is empty.');
        // A user-sized domain is a physical boundary, like an enclosure.
        this.enclosure_width = (d.x_min !== null && d.x_max !== null) ? X1 - X0 : null;
        this.enclosure_height = (d.y_min !== null && d.y_max !== null) ? Y1 - Y0 : null;

        const tol = Math.max(X1 - X0, Y1 - Y0) * 1e-9;
        const out = [];
        for (const r of rects) {
            const s = r.src;
            const o = { ...r };
            delete o.src;
            // Position and size as written when the axis has no pinned edge.
            const axis = (ax, D0, D1) => {
                if (Number.isFinite(ax.min) && Number.isFinite(ax.max)) return [ax.pos, ax.size];
                const lo = Number.isFinite(ax.min) ? ax.min : D0;
                const hi = Number.isFinite(ax.max) ? ax.max : D1;
                return [lo, hi - lo];
            };
            [o.x, o.w] = axis(s.x, X0, X1);
            [o.y, o.h] = axis(s.y, Y0, Y1);
            if (!(o.w > 0) || o.h === 0 || (o.h < 0 && !Number.isFinite(o.h))) {
                throw new Error(`line ${o.line}: rectangle is empty inside the domain.`);
            }
            o.x0 = o.x; o.x1 = o.x + o.w;
            o.y0 = o.h >= 0 ? o.y : o.y + o.h;
            o.y1 = o.h >= 0 ? o.y + o.h : o.y;
            if (o.shape && (o.x0 < X0 - tol || o.x1 > X1 + tol || o.y0 < Y0 - tol || o.y1 > Y1 + tol)) {
                throw new Error(`line ${o.line}: the ${o.shape.prim === 'trap' ? 'trapezoid' : 'n-gon'} lies outside the domain.`);
            }
            if (o.kind === 'diel') {
                // Dielectrics are clipped to the domain, conductors must fit.
                if (o.x0 < X0 - tol || o.x1 > X1 + tol || o.y0 < Y0 - tol || o.y1 > Y1 + tol) {
                    const nx0 = Math.max(o.x0, X0), nx1 = Math.min(o.x1, X1);
                    const ny0 = Math.max(o.y0, Y0), ny1 = Math.min(o.y1, Y1);
                    if (!(nx1 - nx0 > tol) || !(ny1 - ny0 > tol)) continue;
                    Object.assign(o, { x: nx0, w: nx1 - nx0, y: ny0, h: ny1 - ny0, x0: nx0, x1: nx1, y0: ny0, y1: ny1 });
                }
            } else if (o.x0 < X0 - tol || o.x1 > X1 + tol || o.y0 < Y0 - tol || o.y1 > Y1 + tol) {
                throw new Error(`line ${o.line}: conductor lies outside the domain.`);
            }
            out.push(o);
        }

        // Ground slabs behind gnd walls, outside the box given by the user. Both
        // backends then see the same wall: the FDM has no wall condition of its own,
        // the triangular backend absorbs a slab on a gnd wall into the wall.
        const covered = wall => {
            const horiz = wall === 'top' || wall === 'bottom';
            const at = { left: X0, right: X1, top: Y1, bottom: Y0 }[wall];
            const spans = out.filter(r => r.kind === 'gnd').filter(r => {
                const edge = { left: r.x0, right: r.x1, top: r.y1, bottom: r.y0 }[wall];
                return Math.abs(edge - at) <= tol;
            }).map(r => (horiz ? [r.x0, r.x1] : [r.y0, r.y1]));
            return horiz ? intervalsCover(spans, X0, X1, tol) : intervalsCover(spans, Y0, Y1, tol);
        };
        const need = {};
        WALLS.forEach((wall, i) => { need[wall] = b[i] === 'gnd' && !covered(wall); });
        const t = this.wall_t;
        const TX0 = need.left ? X0 - t : X0, TX1 = need.right ? X1 + t : X1;
        const TY0 = need.bottom ? Y0 - t : Y0, TY1 = need.top ? Y1 + t : Y1;
        const slab = (x0, y0, x1, y1) => out.push({ kind: 'gnd', line: 0, er: 1, tand: 0, faces: null,
            x: x0, y: y0, w: x1 - x0, h: y1 - y0, x0, x1, y0, y1, xInf: true, wall: true });
        if (need.bottom) slab(TX0, TY0, TX1, Y0);
        if (need.top) slab(TX0, Y1, TX1, TY1);
        if (need.left) slab(TX0, Y0, X0, Y1);
        if (need.right) slab(X1, Y0, TX1, Y1);

        // Centre the domain on x=0, which is where both backends put it.
        const shift = (TX0 + TX1) / 2;
        const half = (TX1 - TX0) / 2;
        if (Math.abs(shift) > tol) {
            for (const r of out) {
                r.x -= shift;
                // A shape's bounds are its polygon's, which may start left of x.
                if (r.shape) { r.shape = translateShapeX(r.shape, -shift); r.x0 -= shift; r.x1 -= shift; }
                else { r.x0 = r.x; r.x1 = r.x + r.w; }
            }
        }
        this.x_shift = Math.abs(shift) > tol ? shift : 0;
        this.domain_width = 2 * half;
        this.domain_y_min = TY0;
        this.domain_height = TY1;
        // The box given by the user (without the ground slabs), for the preview.
        this.user_domain = { x_min: X0 - shift0(shift, tol), x_max: X1 - shift0(shift, tol), y_min: Y0, y_max: Y1,
            auto: [d.x_min === null, d.x_max === null, d.y_min === null, d.y_max === null] };

        return out;
    }

    _validate_rects(rects, platingMaterial) {
        const where = r => (r.wall ? 'a gnd boundary' : `line ${r.line}`);
        const conds = rects.filter(r => r.kind !== 'diel');
        if (conds.some(r => r.kind === 'sig-') && !conds.some(r => r.kind === 'sig+')) {
            throw new Error('A sig- conductor needs a sig+ conductor.');
        }
        const tol = this.domain_width * 1e-9;
        for (let i = 0; i < conds.length; i++) {
            for (let j = i + 1; j < conds.length; j++) {
                const a = conds[i], c = conds[j];
                if (a.kind === c.kind) continue;
                // Touching counts as a short too: the grid puts both on the same nodes.
                const dx = Math.max(a.x0 - c.x1, c.x0 - a.x1), dy = Math.max(a.y0 - c.y1, c.y0 - a.y1);
                if (dx <= tol && dy <= tol && (!(a.shape || c.shape) || bodyGap(a, c) <= tol)) {
                    throw new Error(`${where(a)} (${a.kind}) and ${where(c)} (${c.kind}) touch or overlap: the conductors are shorted.`);
                }
            }
        }
        for (const r of conds) {
            if (!r.faces) continue;
            const pm = { ...(platingMaterial || {}), ...(r.platingOwn || {}) };
            if (!(pm.sigma > 0) || !(pm.thickness > 0)) {
                throw new Error(`line ${r.line}: plating= needs a plating material: plating_sigma= and plating_t= on the line, or a plating statement.`);
            }
            const joined = conds.some(o => o !== r && o.kind === r.kind && bodyGap(r, o) <= tol);
            if (joined) {
                throw new Error(`line ${r.line}: plating is not supported on a conductor that touches another ${r.kind} rectangle.`);
            }
        }
    }

    _build_lists(rects, platingMaterial) {
        this.dielectrics = [];
        this.conductors = [];
        const own = v => v !== null && v !== undefined;
        this._own_finish = rects.some(r => r.kind !== 'diel' && (own(r.sigma) || own(r.rq) || r.platingOwn));
        for (const r of rects) {
            // A shape's rectangle is its bounding box.
            const [bx, by, bw, bh] = r.shape ? [r.x0, r.y0, r.x1 - r.x0, r.y1 - r.y0] : [r.x, r.y, r.w, r.h];
            if (r.kind === 'diel') {
                const d = new Dielectric(bx, by, bw, bh, r.er, r.tand, r.shape);
                if (r.thin) d.thin_sheet = true;
                d.src_line = r.line;
                if (r.image) d.src_image = true;
                this.dielectrics.push(d);
                continue;
            }
            // Same fields whether the material comes from the statement, the line or the options.
            const plating = r.faces
                ? { rq: 0, thick_corners: false, ...platingMaterial, ...(r.platingOwn || {}), ...r.faces } : null;
            if (plating) { plating.rq = plating.rq ?? 0; plating.thick_corners = !!plating.thick_corners; }
            const polarity = r.kind === 'sig+' ? 1 : (r.kind === 'sig-' ? -1 : 0);
            const c = new Conductor(bx, by, bw, bh, polarity !== 0, polarity, plating, r.shape);
            c.src_line = r.line;   // source line, the editor highlights the rectangle from it
            if (r.image) c.src_image = true;   // generated by mirror=1
            if (own(r.sigma)) c.sigma = r.sigma;   // own conductivity
            if (own(r.rq)) c.rq = r.rq;            // own surface roughness
            this.conductors.push(c);
        }
    }

    // The quasi-static solver meshes rectangles only.
    ensure_mesh() {
        if (this.mesh_backend !== 'triangular' && this.shaped_lines.length) {
            const n = this.shaped_lines.length;
            throw new Error(`The geometry has non-rectangular shapes (trapezoid or n-gon, line${n > 1 ? 's' : ''} ` +
                `${this.shaped_lines.join(', ')}), which the quasi-static solver does not support. ` +
                'Use the full-wave solver.');
        }
        return super.ensure_mesh();
    }

    // The mirror test compares rectangle sets and cannot see paint order. Where two
    // dielectrics of different material overlap, their mirror images must overlap in
    // the same order, or the painted result is not symmetric.
    _paint_order_breaks_symmetry() {
        const tol = this.domain_width * 1e-6;
        const ds = this.dielectrics;
        const same = (a, c) => a.epsilon_r === c.epsilon_r && a.tan_delta === c.tan_delta;
        const partner = ds.map(d => ds.findIndex(o => same(d, o)
            && Math.abs(o.x_min + d.x_max) <= tol && Math.abs(o.x_max + d.x_min) <= tol
            && Math.abs(o.y_min - d.y_min) <= tol && Math.abs(o.y_max - d.y_max) <= tol));
        for (let i = 0; i < ds.length; i++) {
            for (let j = i + 1; j < ds.length; j++) {
                const a = ds[i], c = ds[j];
                if (same(a, c)) continue;
                const o = overlapArea({ x0: a.x_min, x1: a.x_max, y0: a.y_min, y1: a.y_max },
                                      { x0: c.x_min, x1: c.x_max, y0: c.y_min, y1: c.y_max });
                if (!o || o.w <= tol || o.h <= tol) continue;
                if (partner[i] < 0 || partner[j] < 0 || partner[i] > partner[j]) return true;
            }
        }
        return false;
    }

    // Signal conductors that reach a domain boundary. A ground boundary connects the
    // conductor to ground. An open boundary cuts it off where the domain ends: either a
    // mistake or a slotline-like half plane, which has no quasi-static limit (its
    // capacitance grows with the logarithm of the domain size). Grounds are not
    // reported: a ground run to a boundary is a normal way to draw a reference plane.
    boundaryContactWarnings() {
        const out = [];
        const tol = this.domain_width * 1e-9;
        const b = this.boundaries;
        const X0 = -this.domain_width / 2, X1 = this.domain_width / 2;
        const Y0 = this.domain_y_min, Y1 = this.domain_height;
        for (const c of this.conductors) {
            if (!c.is_signal) continue;
            const touches = [c.x_min <= X0 + tol, c.x_max >= X1 - tol, c.y_max >= Y1 - tol, c.y_min <= Y0 + tol];
            const where = c.src_line ? ` (line ${c.src_line})` : '';
            WALLS.forEach((name, k) => {
                if (!touches[k]) return;
                out.push(b[k] === 'open'
                    ? `Signal conductor${where} reaches the open ${name} boundary and is cut off there. ` +
                      'Give it a finite size or enlarge the domain, unless a slotline-like half plane is intended. ' +
                      'Such a line has no quasi-static limit: the characteristic impedance depends on the domain ' +
                      'size on both solvers, and the quasi-static effective permittivity does too. Use the ' +
                      'full-wave solver for the effective permittivity.'
                    : `Signal conductor${where} touches the ${name} ground boundary and is connected to it.`);
            });
        }
        return out;
    }

    // A signal drawn as rectangles that do not touch is still one net: every sig+
    // rectangle is at the same potential, as if connected outside the cross-section.
    signalBodyWarnings() {
        const out = [];
        const tol = this.domain_width * 1e-9;
        const box = c => ({ x0: c.x_min, x1: c.x_max, y0: c.y_min, y1: c.y_max, shape: c.shape });
        for (const [polarity, kind] of [[1, 'sig+'], [-1, 'sig-']]) {
            const list = this.conductors.filter(c => c.is_signal && c.polarity === polarity);
            // Union of touching rectangles into bodies.
            const body = list.map((_, i) => i);
            const find = i => (body[i] === i ? i : (body[i] = find(body[i])));
            for (let i = 0; i < list.length; i++) {
                for (let j = i + 1; j < list.length; j++) {
                    if (bodyGap(box(list[i]), box(list[j])) <= tol) body[find(i)] = find(j);
                }
            }
            const bodies = new Set(list.map((_, i) => find(i))).size;
            if (bodies < 2) continue;
            const lines = [...new Set(list.map(c => c.src_line).filter(l => l > 0))];
            out.push(`The ${kind} conductor is ${bodies} separate bodies` +
                (lines.length ? ` (lines ${lines.join(', ')})` : '') + '. They are solved as one net at the same ' +
                'potential, as if connected outside the cross-section. Make them touch to draw one conductor.');
        }
        return out;
    }

    openBoundaryWarnings(opts) {
        return [...super.openBoundaryWarnings(opts), ...this.boundaryContactWarnings(), ...this.signalBodyWarnings()];
    }

    // Accuracy note for a strongly coupled broadside pair, see BroadsideStriplineSolver.
    _broadside_proximity_note(signals) {
        if (signals.length !== 2 || !this.is_differential) return null;
        const [a, c] = signals;
        const overlap = Math.min(a.x_max, c.x_max) - Math.max(a.x_min, c.x_min);
        const facing_gap = Math.max(a.y_min - c.y_max, c.y_min - a.y_max);
        const w = Math.min(a.width, c.width);
        if (!(facing_gap > 0) || !(overlap > -w) || w / facing_gap < 1.5) return null;
        return { type: 'accuracy', reason: 'broadside-proximity', mode: 'all', message:
            `Strongly coupled broadside pair (trace width ${(w * 1e6).toFixed(0)} µm vs ` +
            `${(facing_gap * 1e6).toFixed(0)} µm facing gap): conductor loss accuracy is reduced. ` +
            `R typically reads up to 50% high in this regime. The full-wave solver models the ` +
            `broadside proximity effect accurately.` };
    }
}

export { CustomGeometrySolver };
