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
import { _clipDomain } from '../src/tri_solver/occ_to_mesh.js';

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

console.log(failures === 0 ? '\nGEOMETRY OK' : `\nGEOMETRY: ${failures} FAILURE(S)`);
process.exit(failures === 0 ? 0 : 1);
