// Geometry-construction unit tests — no solver, no mesh, milliseconds.
//
// Every downstream result trusts the geometry lists (conductors, dielectrics,
// domain box) that the solver constructors build. Bugs here are cheap to make and
// expensive to notice: the has_side_gnd ordering bug (side walls eating 2·t_gnd out
// of the specified enclosure width) survived until an analytic mode test measured
// the cavity. These tests pin the documented layout directly:
//   y: bottom ground −t_gnd..0, substrate 0..h, trace h..h+t
//   enclosure_width = INNER wall-to-wall width (side slabs added outside it)
//   differential trace_spacing = edge-to-edge gap, mirror-symmetric about x=0
// plus the _clipDomain wall-absorption rule that turns full-span boundary slabs
// into PEC walls (this defines the cavity that box_modes_test computes truth from).
//
// Run: node tests/test_geometry.js
import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';
import { CoaxSolver } from '../src/coax.js';
import { RectWaveguideSolver } from '../src/rect_waveguide.js';
import { _clipDomain } from '../src/tri_solver/occ_to_mesh.js';
import { polyRadiusForArea, circlePolygon, shapePoly, shapeArea, shapeContains,
         shapeSegments, shapeSignedDist, isComplement } from '../src/shapes.js';

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}
const near = (x, y, tol = 1e-12) => Math.abs(x - y) < tol;
const signals = (s) => s.conductors.filter(c => c.is_signal);
const grounds = (s) => s.conductors.filter(c => !c.is_signal);

// Base parameters shared by the microstrip-family cases
const BASE = {
    trace_width: 0.35e-3, substrate_height: 0.21e-3, trace_thickness: 35e-6,
    gnd_thickness: 35e-6, epsilon_r: 4.4, tan_delta: 0.02, sigma_cond: 5.8e7, freq: 1e9,
};
const { trace_width: W, substrate_height: H, trace_thickness: T, gnd_thickness: TG } = BASE;

// ---------- open microstrip ----------
{
    const s = new MicrostripSolver({ ...BASE, boundaries: ['open', 'open', 'open', 'gnd'] });
    const sig = signals(s);
    check('microstrip: exactly one signal conductor', sig.length === 1);
    const tr = sig[0];
    check('microstrip: trace centered, width w', near(tr.x_min, -W / 2) && near(tr.x_max, W / 2));
    check('microstrip: trace sits on the substrate (y: h..h+t)', near(tr.y_min, H) && near(tr.y_max, H + T));
    const gnd = grounds(s).find(c => near(c.y_min, -TG) && near(c.y_max, 0));
    check('microstrip: bottom ground slab spans −t_gnd..0, full width',
        !!gnd && near(gnd.x_min, -s.domain_width / 2) && near(gnd.x_max, s.domain_width / 2));
    const sub = s.dielectrics.find(d => d.epsilon_r === BASE.epsilon_r);
    check('microstrip: substrate dielectric spans 0..h with εr/tanδ',
        !!sub && near(sub.y_min, 0) && near(sub.y_max, H) && sub.tan_delta === BASE.tan_delta);
    check('microstrip: open boundaries → no side/lid slabs',
        !grounds(s).some(c => c !== gnd && (c.x_min <= -s.domain_width / 2 + 1e-12 || c.y_min > H)));
    const inDomain = s.conductors.every(c =>
        c.x_min >= -s.domain_width / 2 - 1e-12 && c.x_max <= s.domain_width / 2 + 1e-12 &&
        c.y_min >= -TG - 1e-12 && c.y_max <= s.domain_height + 1e-12);
    check('microstrip: all conductors inside the domain box', inDomain);
}

// ---------- enclosed microstrip (has_side_gnd ordering regression) ----------
{
    const ENC_W = 4e-3, ENC_H = 4e-3;
    const s = new MicrostripSolver({ ...BASE, enclosure_width: ENC_W, enclosure_height: ENC_H,
        boundaries: ['gnd', 'gnd', 'gnd', 'gnd'] });

    // the fixed bug: side walls are ADDED OUTSIDE the specified inner width
    check('enclosed: domain_width = enclosure_width + 2·t_gnd (has_side_gnd ordering bug)',
        near(s.domain_width, ENC_W + 2 * TG), `domain=${(s.domain_width * 1e3).toFixed(3)} mm`);

    const lidY = H + ENC_H;   // lid inner face: substrate top + enclosure height
    check('enclosed: domain_height = h + enclosure_height + t_gnd', near(s.domain_height, lidY + TG));

    const left = grounds(s).find(c => near(c.x_min, -s.domain_width / 2) && near(c.x_max, -s.domain_width / 2 + TG));
    const right = grounds(s).find(c => near(c.x_max, s.domain_width / 2) && near(c.x_min, s.domain_width / 2 - TG));
    check('enclosed: side ground slabs at both edges, t_gnd thick', !!left && !!right);
    check('enclosed: INNER wall-to-wall width == enclosure_width',
        left && right && near(right.x_min - left.x_max, ENC_W),
        left && right ? `inner=${((right.x_min - left.x_max) * 1e3).toFixed(3)} mm` : 'slabs missing');
    const lid = grounds(s).find(c => near(c.y_min, lidY) && near(c.y_max, lidY + TG));
    check('enclosed: lid slab at h+enclosure_height, full width',
        !!lid && near(lid.x_min, -s.domain_width / 2) && near(lid.x_max, s.domain_width / 2));

    // _clipDomain integration: the meshed cavity is exactly the specified enclosure
    const dom = { x_min: -s.domain_width / 2, x_max: s.domain_width / 2, y_min: -TG, y_max: s.domain_height };
    const clip = _clipDomain(dom, s.conductors, s.boundaries, 1e-9);
    check('enclosed: _clipDomain cavity = enclosure_width × (h + enclosure_height)',
        near(clip.X1 - clip.X0, ENC_W) && near(clip.Y0, 0) && near(clip.Y1, lidY));
    check('enclosed: _clipDomain marks all four walls PEC',
        clip.wallPEC.left && clip.wallPEC.right && clip.wallPEC.top && clip.wallPEC.bottom);
}

// ---------- differential microstrip ----------
{
    const SP = 0.2e-3;
    const s = new MicrostripSolver({ ...BASE, trace_spacing: SP, boundaries: ['open', 'open', 'open', 'gnd'] });
    const sig = signals(s).sort((a, b) => a.x_min - b.x_min);
    check('differential: two signal traces', sig.length === 2);
    check('differential: edge-to-edge gap == trace_spacing',
        sig.length === 2 && near(sig[1].x_min - sig[0].x_max, SP));
    check('differential: mirror-symmetric about x=0',
        sig.length === 2 && near(sig[0].x_min, -sig[1].x_max) && near(sig[0].x_max, -sig[1].x_min));
    check('differential: both traces at y = h..h+t',
        sig.every(c => near(c.y_min, H) && near(c.y_max, H + T)));
    check('differential: opposite polarity', sig.length === 2 && sig[0].polarity * sig[1].polarity === -1);
}

// ---------- stripline ----------
{
    const H_TOP = 0.3e-3;
    const s = new MicrostripSolver({ ...BASE, epsilon_r_top: BASE.epsilon_r, tan_delta_top: BASE.tan_delta,
        enclosure_height: H_TOP, boundaries: ['open', 'open', 'gnd', 'gnd'] });
    const lidY = H + H_TOP;
    const lid = grounds(s).find(c => near(c.y_min, lidY) && near(c.y_max, lidY + TG));
    check('stripline: lid slab at h + stripline_top_h', !!lid);
    const top = s.dielectrics.find(d => near(d.y_min, H) && near(d.y_max, lidY));
    check('stripline: top dielectric fills trace level to the lid with er_top/tand_top',
        !!top && top.epsilon_r === BASE.epsilon_r && top.tan_delta === BASE.tan_delta);
    check('stripline: no side slabs for open side boundaries',
        !grounds(s).some(c => c.y_max > 0 && c.y_min < lidY && (c.x_max - c.x_min) < 2 * TG + 1e-12));
}

// ---------- GCPW ----------
{
    const GAP = 0.2e-3, VIA_GAP = 0.3e-3;
    const s = new MicrostripSolver({ ...BASE, use_coplanar_gnd: true, gap: GAP, via_gap: VIA_GAP,
        use_vias: true, boundaries: ['open', 'open', 'open', 'gnd'] });
    const sig = signals(s);
    check('gcpw: one signal trace', sig.length === 1);
    // coplanar grounds: on the substrate top, inner edges at trace edge + gap
    const cop = grounds(s).filter(c => near(c.y_min, H) && near(c.y_max, H + T));
    const rightCop = cop.find(c => c.x_min > 0), leftCop = cop.find(c => c.x_max < 0);
    check('gcpw: coplanar grounds at trace level on both sides', !!rightCop && !!leftCop);
    check('gcpw: signal-to-ground gap == gap',
        rightCop && leftCop && near(rightCop.x_min, W / 2 + GAP) && near(leftCop.x_max, -W / 2 - GAP));
    // vias: ground conductors from the bottom-ground face up through the substrate
    // to the top of the coplanar grounds (y: 0..h+t), one on each side, inner edge
    // via_gap outside the coplanar-ground inner edge
    const vias = grounds(s).filter(c => near(c.y_min, 0) && near(c.y_max, H + T));
    const rightVia = vias.find(c => c.x_min > 0), leftVia = vias.find(c => c.x_max < 0);
    check('gcpw: vias connect coplanar grounds through the substrate', !!rightVia && !!leftVia);
    check('gcpw: via inner edge at gap + via_gap from the trace edge',
        rightVia && near(rightVia.x_min, W / 2 + GAP + VIA_GAP));
}

// ---------- broadside coupled stripline ----------
{
    const opt = { trace_width: 0.2e-3, trace_thickness: 35e-6, x_offset: 0.05e-3, sigma_cond: 5.8e7,
        h_bottom: 0.2e-3, er_bottom: 4.4, tand_bottom: 0.02,
        h_middle: 0.25e-3, er_middle: 3.5, tand_middle: 0.01,
        h_top: 0.2e-3, er_top: 4.4, tand_top: 0.02,
        freq: 1e9, boundaries: ['open', 'open', 'gnd', 'gnd'] };
    const s = new BroadsideStriplineSolver(opt);
    const sig = signals(s).sort((a, b) => a.y_min - b.y_min);
    check('broadside: two signal traces', sig.length === 2);
    check('broadside: traces vertically stacked (disjoint y, lower below upper)',
        sig.length === 2 && sig[0].y_max <= sig[1].y_min + 1e-12);
    check('broadside: horizontal center offset == x_offset', sig.length === 2 &&
        near(Math.abs((sig[1].x_min + sig[1].x_max) / 2 - (sig[0].x_min + sig[0].x_max) / 2), opt.x_offset));
    check('broadside: opposite polarity', sig.length === 2 && sig[0].polarity * sig[1].polarity === -1);
    // ground-to-ground spacing = h_bottom + h_middle + h_top (trace thickness excluded)
    const lid = grounds(s).find(c => c.y_min > 0);
    check('broadside: ground-to-ground spacing = h_bottom + h_middle + h_top',
        !!lid && near(lid.y_min, opt.h_bottom + opt.h_middle + opt.h_top));
    const ers = s.dielectrics.map(d => d.epsilon_r);
    check('broadside: three dielectric layers with the specified εr stack',
        ers.includes(4.4) && ers.includes(3.5));

    // Auto domain width: x_offset only TRANSLATES the upper trace, so the trace-to-wall
    // margin must not grow with it (it used to scale 8·(w + |x_offset|), which turned a
    // 1 mm offset on a 0.2 mm trace into a 19 mm domain). Clearance from the outermost
    // trace edge to the nearer wall stays exactly what the un-offset pair gets.
    const clearance = (sol) => {
        const c = signals(sol), half = sol.domain_width / 2;
        return Math.min(Math.min(...c.map(k => k.x_min)) + half, half - Math.max(...c.map(k => k.x_max)));
    };
    const base = clearance(new BroadsideStriplineSolver({ ...opt, x_offset: 0 }));
    for (const off of [0.05e-3, 1e-3, -1e-3]) {
        check(`broadside: auto width keeps the un-offset trace-to-wall margin at x_offset=${off * 1e3} mm`,
            near(clearance(new BroadsideStriplineSolver({ ...opt, x_offset: off })), base, 1e-15),
            `${(base * 1e3).toFixed(3)} mm`);
    }
}

// ---------- _clipDomain unit cases (synthetic rects) ----------
{
    const dom = { x_min: -2, x_max: 2, y_min: -0.1, y_max: 4 };
    const slab = (x0, y0, x1, y1, is_signal = false) => ({ x_min: x0, y_min: y0, x_max: x1, y_max: y1, is_signal });

    // open boundaries, no slabs → unchanged, no PEC walls
    let c = _clipDomain(dom, [slab(-0.2, 1, 0.2, 1.1, true)], ['open', 'open', 'open', 'open'], 1e-9);
    check('_clipDomain: open + no slabs → domain unchanged, no PEC',
        c.X0 === -2 && c.X1 === 2 && c.Y0 === -0.1 && c.Y1 === 4 &&
        !c.wallPEC.left && !c.wallPEC.right && !c.wallPEC.top && !c.wallPEC.bottom);

    // full-width bottom slab → absorbed: floor moves to slab top, bottom PEC
    c = _clipDomain(dom, [slab(-2, -0.1, 2, 0)], ['open', 'open', 'open', 'gnd'], 1e-9);
    check('_clipDomain: bottom slab absorbed into the floor', c.Y0 === 0 && c.wallPEC.bottom);

    // a SIGNAL conductor spanning the bottom must NOT be absorbed
    c = _clipDomain(dom, [slab(-2, -0.1, 2, 0, true)], ['open', 'open', 'open', 'open'], 1e-9);
    check('_clipDomain: signal conductor is never absorbed', c.Y0 === -0.1 && !c.wallPEC.bottom);

    // a slab NOT touching the domain edge is an interior conductor, not a wall
    c = _clipDomain(dom, [slab(-2, 3.5, 2, 3.6)], ['open', 'open', 'open', 'open'], 1e-9);
    check('_clipDomain: floating slab (gap above) is not absorbed', c.Y1 === 4 && !c.wallPEC.top);

    // full box: floor + lid at the domain top + side slabs spanning floor..lid
    // (like microstrip builds them) — needs the iterative pass: the sides only
    // become full-span after the floor/lid clips. Open boundaries so wallPEC
    // comes from the absorption itself, not the boundary spec.
    const box = [
        slab(-2, -0.1, 2, 0),        // floor
        slab(-2, 3.5, 2, 4),         // lid (touches domain top)
        slab(-2, -0.1, -1.9, 3.5),   // left wall
        slab(1.9, -0.1, 2, 3.5),     // right wall
    ];
    c = _clipDomain({ ...dom }, box, ['open', 'open', 'open', 'open'], 1e-9);
    check('_clipDomain: full box → cavity between all four inner faces',
        c.X0 === -1.9 && c.X1 === 1.9 && c.Y0 === 0 && c.Y1 === 3.5);
    check('_clipDomain: full box → absorption marks all four walls PEC',
        c.wallPEC.left && c.wallPEC.right && c.wallPEC.top && c.wallPEC.bottom);

    // slab NOT reaching a wall (coplanar ground) → not absorbed
    c = _clipDomain(dom, [slab(-1, -0.1, 2, 0)], ['open', 'open', 'open', 'gnd'], 1e-9);
    check('_clipDomain: partial-span slab is not absorbed', c.Y0 === -0.1);
}

// ---------- shape primitives (shapes.js) ----------
// Circles are materialized as convex polygons and the POLYGON is the geometry, not an
// approximation of it: the mesher, the material tagging, the freedom map and the loss
// integrals all test against these same vertices. Two properties carry that design and
// are pinned here.
{
    const r = 1.475e-3;

    // 1. Area matching. The n-gon must enclose EXACTLY pi*r^2, which is what keeps the
    //    capacitance right (and R_dc, which reads shapeArea directly).
    for (const n of [32, 64, 128]) {
        const disk = { shape: { type: 'circle', cx: 0, cy: 0, r, n, phase: 0 } };
        check(`shapes: n=${n} polygon encloses exactly pi*r^2`,
            near(shapeArea(disk) / (Math.PI * r * r), 1, 1e-12));
        check(`shapes: n=${n} half polygon is exactly half`,
            near(shapeArea(disk, { half: true }) / (shapeArea(disk) / 2), 1, 1e-12));
    }
    check('shapes: area-matched radius sits between in- and circumradius of the circle',
        polyRadiusForArea(r, 64) > r && polyRadiusForArea(r, 64) / Math.cos(Math.PI / 64) > r);

    // 2. Boundary ownership, both polarities. Vertices AND side midpoints must classify
    //    as inside. Midpoints are the load-bearing case: buildTriFreedomMap decides a PEC
    //    edge by its midpoint, so a complement shape (the coax shield) whose side
    //    midpoints read as "outside the metal" would leak the whole outer boundary.
    const n = 64, tol = r * 1e-9;
    const disk = { shape: { type: 'circle', cx: 0, cy: 0, r, n, phase: 0 } };
    const shell = { shape: { type: 'outside_circle', cx: 0, cy: 0, r, n, phase: 0 } };
    const poly = shapePoly(disk.shape);
    let vertexOk = true, midOk = true;
    for (let k = 0; k < n; k++) {
        const j = (k + 1) % n;
        const vx = poly[2 * k], vy = poly[2 * k + 1];
        const mx = (vx + poly[2 * j]) / 2, my = (vy + poly[2 * j + 1]) / 2;
        if (!shapeContains(disk, vx, vy, tol) || !shapeContains(shell, vx, vy, tol)) vertexOk = false;
        if (!shapeContains(disk, mx, my, tol) || !shapeContains(shell, mx, my, tol)) midOk = false;
    }
    check('shapes: both polarities own their polygon vertices', vertexOk);
    check('shapes: both polarities own their side MIDPOINTS', midOk);
    check('shapes: disk and shell partition the plane away from the boundary',
        shapeContains(disk, 0, 0, 0) && !shapeContains(shell, 0, 0, 0) &&
        !shapeContains(disk, 2 * r, 0, 0) && shapeContains(shell, 2 * r, 0, 0));
    check('shapes: a shell has no finite area', shapeArea(shell) === 0 && isComplement(shell.shape));
    check('shapes: signed distance is negative inside, positive outside',
        shapeSignedDist(disk.shape, 0, 0) < 0 && shapeSignedDist(disk.shape, 2 * r, 0) > 0);
    check('shapes: full boundary has n segments, half drops the symmetry chord',
        shapeSegments(disk).length === n && shapeSegments(disk, { half: true }).length === n / 2);
    check('shapes: circlePolygon is deterministic',
        circlePolygon(0, 0, r, n, 0).every((v, i) => v === circlePolygon(0, 0, r, n, 0)[i]));

    // A shapeless object must take the LEGACY bbox path unchanged — this is the
    // invariant that keeps every rectangular medium bit-identical.
    const rect = { x_min: -1, x_max: 2, y_min: -3, y_max: 4, width: 3, height: 7 };
    check('shapes: shapeless objects fall back to the bbox test',
        shapeContains(rect, 0, 0, 0) && shapeContains(rect, 2, 4, 0) &&
        !shapeContains(rect, 2.1, 0, 0) && shapeArea(rect) === 21 &&
        shapeSegments(rect).length === 0);
}

// ---------- coaxial ----------
{
    const d = 0.92e-3, D = 2.95e-3, a = d / 2, b = D / 2;
    const s = new CoaxSolver({ inner_diameter: d, dielectric_diameter: D,
        epsilon_r: 2.1, tan_delta: 2e-4, sigma_cond: 5.8e7 });
    const sig = signals(s), gnd = grounds(s);
    check('coax: one signal conductor and one ground', sig.length === 1 && gnd.length === 1);
    check('coax: centre conductor is a disk of radius d/2',
        sig[0].shape.type === 'circle' && near(sig[0].shape.r, a));
    // The shield is the COMPLEMENT of the dielectric disk, not a ring: it owns the outer
    // boundary (making it PEC and giving it loss edges) without meshing any metal or
    // leaving dead air cavities in the corners of a bounding box.
    check('coax: shield is the complement of the dielectric disk',
        gnd[0].shape.type === 'outside_circle' && near(gnd[0].shape.r, b));
    check('coax: the dielectric and the meshed domain are the SAME polygon object',
        s.dielectrics[0].shape === s.domain_shape);
    check('coax: fully enclosed — every wall is gnd', s.boundaries.every(v => v === 'gnd'));
    check('coax: single-ended', s.is_differential === false);
    check('coax: DC-resistance area is pi*a^2, not the bbox 4a^2',
        near(shapeArea(sig[0]) / (Math.PI * a * a), 1, 1e-12));
    // n % 4 == 0 with phase 0 puts vertices exactly on both axes, which is what makes
    // the x >= 0 half an EXACT half — required by the half-domain symmetry solve.
    check('coax: segment counts are multiples of 4',
        s.n_inner % 4 === 0 && s.n_outer % 4 === 0, `${s.n_inner}, ${s.n_outer}`);
    check('coax: shapes declare mirror symmetry',
        [...s.conductors, ...s.dielectrics].every(o => o.shape.xSymmetric === true));
    check('coax: domain box strictly encloses the meshed disk',
        s.domain_width / 2 > b && s.domain_height > b && s.t_gnd > b);

    // Plating on a circle has no faces to select between; what it selects is WHICH
    // CONDUCTOR carries the layer. A block naming neither (the pre-inner/outer shape, or
    // one filled from a rectangular line's face checkboxes) means the centre conductor.
    const platingOf = (sel) => new CoaxSolver({ inner_diameter: d, dielectric_diameter: D,
        epsilon_r: 2.1, plating: { sigma: 6.3e7, thickness: 4e-6, rq: 0, ...sel } })
        .conductors.map(c => c.plating);
    const plated = platingOf({ top: true, sides: false, bottom: false });
    check('coax: plating is normalized to all-around on the centre conductor',
        plated[0].all === true && plated[0].top && plated[0].sides && plated[0].bottom &&
        plated[0].thick_corners === false);
    check('coax: a plating block naming no conductor plates the centre one only',
        plated[1] === null);
    const outerOnly = platingOf({ inner: false, outer: true });
    check('coax: outer-only plating plates the shield and leaves the centre bare',
        outerOnly[0] === null && outerOnly[1] !== null && outerOnly[1].all === true);
    const bothPlated = platingOf({ inner: true, outer: true });
    check('coax: both conductors can be plated at once',
        bothPlated[0] !== null && bothPlated[1] !== null);
    check('coax: selecting neither conductor is bare metal',
        platingOf({ inner: false, outer: false }).every(p => p === null));

    // A shaped ground must never be absorbed into a wall — its bbox spans the domain but
    // its body does not fill it.
    const cl = _clipDomain({ x_min: -s.domain_width / 2, x_max: s.domain_width / 2,
                             y_min: -s.t_gnd, y_max: s.domain_height },
                           s.conductors, s.boundaries, s.domain_width * 1e-9);
    check('coax: _clipDomain leaves the domain untouched (no full-span slab)',
        near(cl.X0, -s.domain_width / 2) && near(cl.X1, s.domain_width / 2) &&
        near(cl.Y0, -s.t_gnd) && near(cl.Y1, s.domain_height));

    const rejects = (o, label) => {
        try { new CoaxSolver({ inner_diameter: 1e-3, dielectric_diameter: 3e-3, epsilon_r: 2.1, ...o }); }
        catch { check(`coax: rejects ${label}`, true); return; }
        check(`coax: rejects ${label}`, false);
    };
    rejects({ dielectric_diameter: 1e-3 }, 'D == d');
    rejects({ dielectric_diameter: 1.01e-3 }, 'a dielectric annulus too thin to mesh');
    rejects({ inner_diameter: -1e-3 }, 'a negative diameter');
    rejects({ epsilon_r: 0.5 }, 'epsilon_r below 1');
    rejects({ mesh_backend: 'rectilinear' }, 'the quasi-static backend');
}

// ---------- rectangular waveguide ----------
// The first SOURCE-FREE medium: the walls are boundary conditions, not meshed bodies, so
// the conductor list is empty. That is exactly what validateTriMesh's enclosed-domain
// exemption keys on, so the _clipDomain assertion below is load-bearing.
{
    const A = 22.86e-3, B = 10.16e-3;                 // WR-90
    const s = new RectWaveguideSolver({ width: A, height: B, sigma_cond: 5.8e7 });

    check('waveguide: no conductors — the walls are boundary conditions',
        s.conductors.length === 0);
    check('waveguide: one homogeneous dielectric filling the guide',
        s.dielectrics.length === 1 && near(s.dielectrics[0].x_min, -A / 2) &&
        near(s.dielectrics[0].x_max, A / 2) && near(s.dielectrics[0].y_min, 0) &&
        near(s.dielectrics[0].y_max, B));
    check('waveguide: fully enclosed — every wall is gnd', s.boundaries.every(v => v === 'gnd'));
    check('waveguide: single-ended', s.is_differential === false);
    check('waveguide: domain box IS the guide interior',
        near(s.domain_width, A) && near(s.domain_height, B) && s.t_gnd === 0);
    check('waveguide: backend dispatch flags',
        s.mode_type === 'waveguide' && s.tri_symmetry === false &&
        s.has_potential === false && s.allow_dc === false);

    // kc is purely geometric — no er — which is what lets beta(f) be propagated
    // analytically from a single eigensolve.
    check('waveguide: kc = pi/max(a,b)', near(s.kc_analytic, Math.PI / A));
    check('waveguide: WR-90 cutoff is 6.557 GHz', near(s.fc / 1e9, 6.5571, 1e-3), (s.fc / 1e9).toFixed(4));
    check('waveguide: WR-90 second cutoff is 13.114 GHz', near(s.fc2 / 1e9, 13.1143, 1e-3),
        (s.fc2 / 1e9).toFixed(4));
    check('waveguide: the published 8.2-12.4 GHz band is single-mode',
        s.fc < 8.2e9 && s.fc2 > 12.4e9);

    // Dielectric fill lowers both cutoffs by sqrt(er) and leaves kc alone.
    const filled = new RectWaveguideSolver({ width: A, height: B, epsilon_r: 2.1 });
    check('waveguide: er scales fc by 1/sqrt(er), kc unchanged',
        near(filled.kc_analytic, s.kc_analytic) &&
        near(filled.fc * Math.sqrt(2.1), s.fc, 1e-3));

    // b > a: the fundamental is TE01, so max(a,b) must pick the HEIGHT.
    const tall = new RectWaveguideSolver({ width: B, height: A });
    check('waveguide: b > a puts the fundamental at pi/b (TE01)',
        near(tall.kc_analytic, Math.PI / A) && near(tall.fc, s.fc));

    // One continuous wall surface: no faces to select between.
    const plated = new RectWaveguideSolver({ width: A, height: B,
        plating: { sigma: 6.3e7, thickness: 4e-6, rq: 0, top: true, sides: false, bottom: false } });
    check('waveguide: plating is normalized to all-around',
        plated.plating.all === true && plated.plating.top && plated.plating.sides &&
        plated.plating.bottom && plated.plating.thick_corners === false);
    check('waveguide: no plating when none requested', s.plating === null);

    // The enclosed-domain exemption in validateTriMesh requires all four wallPEC true and
    // an untouched domain — with no conductors there is no full-span slab to absorb.
    const cl = _clipDomain({ x_min: -A / 2, x_max: A / 2, y_min: 0, y_max: B },
                           s.conductors, s.boundaries, A * 1e-9);
    check('waveguide: _clipDomain sets all four wallPEC and leaves the domain untouched',
        cl.wallPEC.left && cl.wallPEC.right && cl.wallPEC.top && cl.wallPEC.bottom &&
        near(cl.X0, -A / 2) && near(cl.X1, A / 2) && near(cl.Y0, 0) && near(cl.Y1, B));

    // Warnings, not throws: an out-of-band frequency is legal, just annotated.
    check('waveguide: no warnings in band', s.waveguideWarnings(10e9).length === 0);
    check('waveguide: warns below cutoff', s.waveguideWarnings(5e9).length === 1);
    check('waveguide: warns when over-moded', s.waveguideWarnings(15e9).length === 1);

    const rejects = (o, label) => {
        try { new RectWaveguideSolver({ width: A, height: B, ...o }); }
        catch { check(`waveguide: rejects ${label}`, true); return; }
        check(`waveguide: rejects ${label}`, false);
    };
    rejects({ width: -1e-3 }, 'a negative width');
    rejects({ height: 0 }, 'a zero height');
    rejects({ epsilon_r: 0.5 }, 'epsilon_r below 1');
    rejects({ tan_delta: -1e-3 }, 'a negative loss tangent');
    rejects({ height: A / 200 }, 'an aspect ratio too extreme to mesh');
    rejects({ mesh_backend: 'rectilinear' }, 'the quasi-static backend');
}

console.log(failures === 0 ? '\nGEOMETRY OK' : `\nGEOMETRY: ${failures} FAILURE(S)`);
process.exit(failures === 0 ? 0 : 1);
