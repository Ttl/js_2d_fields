// Custom geometry trapezoids and n-gons without solving: the polygons the text makes,
// mirroring, unit conversion, validation, the quasi-static refusal, the coax
// conversion and the polygon helpers in shapes.js. Runs in well under a second.
import { parseAndEvaluate, parseGeometryText, serializeGeometry, changeUnitsInText, solverToGeometryText }
    from '../src/custom_geometry_text.js';
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { buildSolverFromParams } from '../src/solver_factory.js';
import { isXSymmetric, conductorSwapSymmetric } from '../src/geometry_symmetry.js';
import { bodyDistance, shapeLoops, shapeArea, shapeSegments, shapeContains, shapeFaceAt } from '../src/shapes.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
const near = (a, b, tol = 1e-12) => Math.abs(a - b) <= tol;
const polyOf = (text, i = 0) => parseAndEvaluate(text).rects[i].shape;
const errorsOf = text => parseAndEvaluate(text).errors.map(e => e.message).join('; ');
const U = 'units mm\n';

// --- Trapezoid ---
{
    const s = polyOf(U + 'sig+ trap x=0 y=1 w=2 h=0.5 angle=45');
    const p = Array.from(s.poly).map(v => v * 1e3);
    check('trap: base at y, top narrowed by h*tan(angle) on both sides',
        [0, 1, 2, 1, 1.5, 1.5, 0.5, 1.5].every((v, i) => near(p[i], v, 1e-9)), p.join(' '));
    check('trap: edges named bottom, sides, top, sides', s.faces.join() === 'bottom,sides,top,sides');
    const r = polyOf(U + 'sig+ trap x=0 y=1 w=2 h=0.5 angle=45 angle2=0');
    check('trap: angle2 sets the right side', near(r.poly[4] * 1e3, 2, 1e-9) && near(r.poly[6] * 1e3, 0.5, 1e-9));
    const o = polyOf(U + 'sig+ trap x=0 y=1 w=2 h=0.5 angle2=45');
    check('trap: angle2 alone applies to both sides', near(o.poly[6] * 1e3, 0.5, 1e-9) && near(o.poly[4] * 1e3, 1.5, 1e-9));
    const n = polyOf(U + 'sig+ trap x=0 y=1 w=2 h=-0.5 angle=45');
    check('trap: negative h puts the narrow face below y', near(n.poly[1] * 1e3, 0.5, 1e-9) && near(n.poly[5] * 1e3, 1, 1e-9)
        && near(n.poly[0] * 1e3, 0.5, 1e-9) && n.faces[0] === 'bottom');
    const g = parseAndEvaluate(U + 'sig+ trap x=0 y=1 w=2 h=0.5 angle=-45').rects[0];
    check('trap: negative angle widens the far face, bounds follow', near(g.x.min * 1e3, -0.5, 1e-9) && near(g.x.max * 1e3, 2.5, 1e-9));
    const z = parseAndEvaluate(U + 'sig+ trap x=0 y=1 w=2 h=0.5 angle=0').rects[0];
    check('trap: zero angles give a plain rectangle', z.shape === null && near(z.x.size, 2e-3) && near(z.y.size, 0.5e-3));
    check('trap: angles that close the top are rejected', /no face opposite/.test(errorsOf(U + 'sig+ trap x=0 y=0 w=1 h=1 angle=30')));
    check('trap: needs finite edges', /finite/.test(errorsOf(U + 'sig+ trap x=-inf y=0 w=inf h=1 angle=10')));
    check('trap: ngon keys rejected', /r= does not apply to a trapezoid/.test(errorsOf(U + 'sig+ trap x=0 y=0 w=1 h=1 r=2')));
    check('rect: angle rejected', /angle= does not apply to a rectangle/.test(errorsOf(U + 'sig+ x=0 y=0 w=1 h=1 angle=2')));
}

// --- N-gon ---
{
    const s = polyOf(U + 'gnd ngon x=1 y=2 r=0.5 n=6');
    const p = s.poly;
    check('ngon: n vertices, first on top', p.length === 12 && near(p[0], 1e-3) && near(p[1], 2.5e-3, 1e-15));
    let mirror = true;
    for (let k = 1; k < 6; k++) mirror = mirror && p[2 * (6 - k)] - 1e-3 === -(p[2 * k] - 1e-3) && p[2 * (6 - k) + 1] === p[2 * k + 1];
    check('ngon: vertices mirror exactly about the centre line', mirror);
    const odd = polyOf(U + 'gnd ngon x=0 y=0 r=1 n=5');
    check('ngon: odd n is mirror symmetric too', odd.poly[0] === 0 && odd.poly[2] === -odd.poly[8] && odd.poly[3] === odd.poly[9]);
    const rot = polyOf(U + 'gnd ngon x=0 y=0 r=1 n=4 rot=45');
    check('ngon: rot turns it (a square)', near(rot.poly[0], -Math.SQRT1_2 * 1e-3, 1e-15) && near(rot.poly[1], Math.SQRT1_2 * 1e-3, 1e-15));
    const ring = polyOf(U + 'gnd ngon x=0 y=0 r=1 r_in=0.8 n=16');
    check('ngon: r_in makes a ring', ring.type === 'ring' && ring.hole.length === 32);
    check('ngon: area of a ring', near(shapeArea({ shape: ring }), 0.5 * 16 * Math.sin(2 * Math.PI / 16) * (1 - 0.64) * 1e-6, 1e-15));
    check('ngon: bad n rejected', /whole number/.test(errorsOf(U + 'gnd ngon x=0 y=0 r=1 n=2.5')));
    check('ngon: r_in >= r rejected', /r_in must be/.test(errorsOf(U + 'gnd ngon x=0 y=0 r=1 r_in=1 n=8')));
    check('ngon: plating must be all', /plating=all/.test(errorsOf(U + 'gnd ngon x=0 y=0 r=1 n=8 plating=top')));
    const pl = parseAndEvaluate(U + 'gnd ngon x=0 y=0 r=1 n=8 plating=all').rects[0].plating;
    check('ngon: plating=all plates the whole surface', pl.all && pl.top && pl.sides && pl.bottom);
}

// --- Ellipse ---
{
    const e = polyOf(U + 'diel ellipse x=0 y=0 rx=2 ry=1 n=64 er=2');
    const sides = [];
    for (let i = 0; i < 64; i++) { const j = (i + 1) % 64; sides.push(Math.hypot(e.poly[2 * j] - e.poly[2 * i], e.poly[2 * j + 1] - e.poly[2 * i + 1])); }
    check('ellipse: vertices spaced evenly along the outline', Math.max(...sides) / Math.min(...sides) < 1.02,
        `${(Math.max(...sides) / Math.min(...sides)).toFixed(4)}`);
    check('ellipse: first vertex on top, exact mirror pairs', e.poly[0] === 0 && e.poly[1] === 1e-3
        && e.poly[2 * 63] === -e.poly[2] && e.poly[2 * 63 + 1] === e.poly[3]);
    check('ellipse: area near pi*rx*ry', Math.abs(shapeArea({ shape: e }) / (Math.PI * 2e-6) - 1) < 2e-3);
    const c = polyOf(U + 'gnd ellipse x=0.3 y=0.1 rx=1 ry=1 n=24'), g = polyOf(U + 'gnd ngon x=0.3 y=0.1 r=1 n=24');
    check('ellipse: equal semi-axes give the n-gon exactly', c.poly.every((v, i) => v === g.poly[i]));
    const ring = polyOf(U + 'gnd ellipse x=0 y=0 rx=2 ry=1 rx_in=1.8 ry_in=0.8 n=64');
    check('ellipse: rx_in, ry_in make a ring', ring.type === 'ring' && ring.hole.length === 128);
    check('ellipse: a ring needs both inner semi-axes', /needs rx_in and ry_in/.test(errorsOf(U + 'gnd ellipse x=0 y=0 rx=2 ry=1 rx_in=1 n=64')));
    check('ellipse: plating must be all', /plating=all/.test(errorsOf(U + 'gnd ellipse x=0 y=0 rx=2 ry=1 n=16 plating=top')));
    check('ellipse: r= rejected', /r= does not apply to an ellipse/.test(errorsOf(U + 'gnd ellipse x=0 y=0 r=2 n=16')));
}

// --- Rounded corners and walls ---
{
    const area = t => shapeArea({ shape: polyOf(U + t) }) * 1e6;
    const st = polyOf(U + 'sig+ x=0 y=0 w=2 h=1 radius=0.5');
    // The arcs are inscribed polygons: area short by the segment areas.
    check('radius: half the height gives a stadium', st.type === 'polygon' && Math.abs(area('sig+ x=0 y=0 w=2 h=1 radius=0.5') / (1 + Math.PI / 4) - 1) < 5e-3
        && Math.abs(Math.min(...Array.from(st.poly).filter((_, i) => i % 2 === 1))) < 1e-15);
    const rr = polyOf(U + 'sig+ x=0 y=0 w=2 h=1 radius=0.2 radius_bottom=0');
    const has = (p, x, y) => [...Array(p.length / 2).keys()].some(k => Math.abs(p[2 * k] - x) < 1e-15 && Math.abs(p[2 * k + 1] - y) < 1e-15);
    check('radius_bottom=0: sharp lower corners, rounded upper ones', has(rr.poly, 0, 0) && has(rr.poly, 2e-3, 0) && !has(rr.poly, 2e-3, 1e-3));
    check('radius: arc halves take the faces they run into', rr.faces.filter(f => f === 'top').length > 1 && rr.faces.filter(f => f === 'sides').length > 2
        && shapeFaceAt(rr, 1e-3, 1e-3).face === 'top' && shapeFaceAt(rr, 2e-3, 0.5e-3).face === 'sides');
    const tr = polyOf(U + 'sig+ trap x=0 y=0 w=2 h=0.5 angle=30 radius=0.1');
    check('radius: rounds a trapezoid inside its outline', area('sig+ trap x=0 y=0 w=2 h=0.5 angle=30 radius=0.1')
        < area('sig+ trap x=0 y=0 w=2 h=0.5 angle=30') && tr.faces.includes('top'));
    check('radius: too large rejected', /larger than the sides allow/.test(errorsOf(U + 'sig+ x=0 y=0 w=2 h=1 radius=0.6')));
    const w = polyOf(U + 'gnd x=0 y=0 w=2 h=1 wall=0.1');
    check('wall: a hollow rectangle', w.type === 'ring' && Math.abs(area('gnd x=0 y=0 w=2 h=1 wall=0.1') - (2 - 1.8 * 0.8)) < 1e-9);
    const ws = polyOf(U + 'gnd x=0 y=0 w=2 h=1 radius=0.5 wall=0.1');
    const minY = Math.min(...Array.from(ws.hole).filter((_, i) => i % 2 === 1));
    check('wall: a stadium shell keeps the corner centres', ws.type === 'ring' && Math.abs(minY - 0.1e-3) < 1e-15);
    check('wall: too thick rejected', /too thick/.test(errorsOf(U + 'gnd x=0 y=0 w=2 h=1 wall=0.5')));
    check('wall: plating on a shell must be all', /plating=all/.test(errorsOf(U + 'gnd x=0 y=0 w=2 h=1 wall=0.1 plating=top')));
    check('radius: needs finite edges', /finite/.test(errorsOf(U + 'gnd x=-inf y=0 w=inf h=1 radius=0.1')));
    const um = changeUnitsInText(U + 'rc = 0.1\nsig+ x=0 y=0 w=2 h=1 radius=rc radius_bottom=0.05 wall=0.02\n', 'um');
    check('units: radius and wall are lengths', /rc = 100\b/.test(um) && /radius_bottom=50 wall=20/.test(um), um);
    let qs = '';
    try {
        new CustomGeometrySolver({ text: U + 'bounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4\nsig+ x=-0.1 w=0.2 y=0.2 h=0.035 radius=5um\n' }).ensure_mesh();
    } catch (e) { qs = e.message; }
    check('radius: the quasi-static solver refuses rounded corners', /quasi-static solver does not support/.test(qs), qs);
}

// --- Thick plating: one solver option for every conductor ---
{
    const text = U + 'bounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4\n' +
        'sig- x=-0.4 w=0.3 y=0.2 h=0.035 plating=top plating_sigma=1e7 plating_t=4um\n' +
        'sig+ x=0.1 w=0.3 y=0.2 h=0.035 plating=top,sides plating_sigma=2e7 plating_t=2um\n';
    const th = opt => new CustomGeometrySolver({ text, thick_plating: opt }).conductors.filter(c => c.plating).map(c => c.plating.thick_corners);
    check('thick plating: the solver option sets every conductor, off by default',
        th(undefined).every(v => v === false) && th(true).every(v => v === true));
    check('thick plating: thick_corners= in the text points to the option',
        /Model Thick Plating option/.test(errorsOf(U + 'plating sigma=1e7 t=4um thick_corners=1\n'))
        && /unknown rectangle key 'plating_thick_corners'/.test(errorsOf(U + 'sig+ x=0 y=0 w=1 h=1 plating_thick_corners=1\n')));
}

// --- Overlapping conductors of different metal ---
{
    const T = extra => U + 'bounds open open gnd gnd\ndomain -1 1 0 0.4\ndiel x=-inf w=inf y=0 h=0.4 er=3.5\n' +
        'sig+ x=-0.15 w=0.3 y=0.18 h=0.035 sigma=5.8e7\n' + extra;
    const refuse = t => { try { new CustomGeometrySolver({ text: t }).ensure_mesh(); return ''; } catch (e) { return e.message; } };
    const diff = new CustomGeometrySolver({ text: T('sig+ x=-0.14 w=0.28 y=0.19 h=0.015 sigma=5.8e5\n') });
    check('overlap: different metals found', JSON.stringify(diff.metal_overlaps) === '[[5,6]]', JSON.stringify(diff.metal_overlaps));
    check('overlap: warned, with the rule', diff.openBoundaryWarnings().some(w => /later line's metal fills the overlap/.test(w)));
    check('overlap: the quasi-static solver refuses it', /different metal overlap \(lines 5 and 6\).*quasi-static/.test(refuse(T('sig+ x=-0.14 w=0.28 y=0.19 h=0.015 sigma=5.8e5\n'))));
    check('overlap: the same metal is fine', refuse(T('sig+ x=-0.14 w=0.28 y=0.19 h=0.015\n'.replace('\n', ' sigma=5.8e7\n'))) === '');
    check('overlap: touching blocks of different metal are fine', refuse(T('sig+ x=0.15 w=0.05 y=0.18 h=0.035 sigma=5.8e5\n')) === '');
    check('overlap: different roughness counts as different metal',
        new CustomGeometrySolver({ text: T('sig+ x=-0.14 w=0.28 y=0.19 h=0.015 sigma=5.8e7 rq=1um\n') }).metal_overlaps.length === 1);
}

// --- Mirror ---
{
    const g = parseAndEvaluate(U + 'sig+ trap x=0.1 y=0 w=0.3 h=0.035 angle=30 angle2=10 mirror=1');
    const [a, b] = g.rects;
    check('mirror: a trapezoid gets a sig- image', g.rects.length === 2 && b.kind === 'sig-' && b.image);
    const ok = [0, 1, 2, 3].every(i => {
        const x = a.shape.poly[2 * i], y = a.shape.poly[2 * i + 1];
        return [0, 1, 2, 3].some(k => b.shape.poly[2 * k] === -x && b.shape.poly[2 * k + 1] === y);
    });
    check('mirror: image vertices mirror exactly', ok);
    const face = (s, x, y) => shapeFaceAt(s, x, y).face;
    check('mirror: image faces follow the geometry', face(b.shape, -0.25e-3, 0) === 'bottom'
        && face(b.shape, -0.25e-3, 0.035e-3) === 'top' && face(b.shape, -0.4e-3, 0.01e-3) === 'sides');
    check('mirror: symmetric n-gon on x=0 stays one', parseAndEvaluate(U + 'gnd ngon x=0 y=0 r=1 n=8 mirror=1').rects.length === 1);
    check('mirror: asymmetric trapezoid on x=0 rejected',
        /symmetric about it/.test(errorsOf(U + 'sig+ trap x=-0.1 y=0 w=0.3 h=0.035 angle=30 mirror=1')));
}

// --- Round trip and units ---
{
    const text = U + 'a = 30\nsig+ trap x=-0.1 w=0.2 y=0 h=0.035 angle=a angle2=10\ngnd ngon x=0 y=-1 r=0.3 r_in=0.2 n=12 rot=15\n';
    check('serialize keeps the shape word', serializeGeometry(parseGeometryText(text)) === text);
    const um = changeUnitsInText(text, 'um');
    const g0 = parseAndEvaluate(text), g1 = um && parseAndEvaluate(um);
    const same = g1 && g0.rects.every((r, i) => Array.from(r.shape.poly).every((v, k) => near(v, g1.rects[i].shape.poly[k], 1e-15)));
    check('units: mm -> um keeps angles, n and rot', same && /angle=a/.test(um) && /n=12 rot=15/.test(um) && /r=300/.test(um), um);
}

// --- Solver: validation and the quasi-static refusal ---
{
    const coax = r => U + `bounds open open open open\ndomain -2 2 -2 2\ndiel ngon x=0 y=0 r=1.475 n=64 er=2.1\n` +
        `sig+ ngon x=0 y=0 r=${r} n=32\ngnd ngon x=0 y=0 r=1.6 r_in=1.475 n=64\n`;
    let s = null;
    try { s = new CustomGeometrySolver({ text: coax(0.46) }); } catch (e) { check('solver: conductor inside a ring hole is not a short', false, e.message); }
    if (s) {
        check('solver: conductor inside a ring hole is not a short', true);
        check('solver: centred coax is mirror symmetric', isXSymmetric(s.conductors, s.dielectrics, s.domain_width));
        let msg = '';
        try { s.ensure_mesh(); } catch (e) { msg = e.message; }
        check('solver: quasi-static refuses shapes', /quasi-static solver does not support/.test(msg) && /lines 4, 5, 6/.test(msg), msg);
        s.mesh_backend = 'triangular';
        let ok = true;
        try { s.ensure_mesh(); } catch { ok = false; }
        check('solver: full-wave accepts shapes', ok);
    }
    let short = '';
    try { new CustomGeometrySolver({ text: coax(1.5) }); } catch (e) { short = e.message; }
    check('solver: conductor reaching the ring is a short', /shorted/.test(short), short);
    let outside = '';
    try { new CustomGeometrySolver({ text: U + 'bounds open open open gnd\ndomain -1 1 0 1\nsig+ ngon x=0 y=0.5 r=0.6 n=8\n' }); } catch (e) { outside = e.message; }
    check('solver: shape outside the domain rejected', /outside the domain/.test(outside), outside);
    const rect = new CustomGeometrySolver({ text: U + 'bounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4\nsig+ trap x=-0.1 w=0.2 y=0.2 h=0.035 angle=0\n' });
    let ok = true;
    try { rect.ensure_mesh(); } catch { ok = false; }
    check('solver: a trapezoid with zero angles solves on the quasi-static solver', ok && rect.shaped_lines.length === 0);

    const pair = mirror => new CustomGeometrySolver({ text: U + 'bounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4\n' +
        (mirror ? 'sig+ trap x=0.1 w=0.3 y=0.2 h=0.035 angle=30 angle2=5 mirror=1\n'
            : 'sig+ trap x=0.1 w=0.3 y=0.2 h=0.035 angle=30 angle2=5\nsig- trap x=-0.4 w=0.3 y=0.2 h=0.035 angle=30 angle2=5\n') });
    const m = pair(true), c = pair(false);
    check('symmetry: mirrored trapezoid pair is symmetric', isXSymmetric(m.conductors, m.dielectrics, m.domain_width)
        && conductorSwapSymmetric(m.conductors, m.dielectrics) === true);
    check('symmetry: copied (unmirrored) trapezoid pair is not', !isXSymmetric(c.conductors, c.dielectrics, c.domain_width)
        && conductorSwapSymmetric(c.conductors, c.dielectrics) === false);
}

// --- Polygon helpers ---
{
    const rect = (x0, x1, y0, y1) => ({ x_min: x0, x_max: x1, y_min: y0, y_max: y1 });
    check('bodyDistance: crossing bars overlap', bodyDistance(rect(-2, 2, -0.1, 0.1), rect(-0.1, 0.1, -2, 2)) === 0);
    check('bodyDistance: rectangles apart', near(bodyDistance(rect(0, 1, 0, 1), rect(2, 3, 0, 1)), 1));
    const ring = parseAndEvaluate(U + 'gnd ngon x=0 y=0 r=2 r_in=1.5 n=4').rects[0].shape;
    const disk = parseAndEvaluate(U + 'gnd ngon x=0 y=0 r=0.5 n=4').rects[0].shape;
    check('bodyDistance: disk in a ring hole is apart by the gap',
        near(bodyDistance({ shape: disk }, { shape: ring }), Math.SQRT1_2 * 1e-3, 1e-15));
    check('shapeContains: ring hole is not metal', !shapeContains({ shape: ring }, 0, 0) && shapeContains({ shape: ring }, 1.75e-3, 0));
    // Half of a centred ring is one C-shaped loop with half the area, no segment on the plane.
    const r16 = parseAndEvaluate(U + 'gnd ngon x=0 y=0 r=2 r_in=1.5 n=16').rects[0].shape;
    const loops = shapeLoops(r16, { half: true });
    check('half ring: one loop', loops.length === 1);
    check('half ring: half the area', near(shapeArea({ shape: r16 }, { half: true }), shapeArea({ shape: r16 }) / 2, 1e-18));
    check('half ring: plane is not a surface', shapeSegments({ shape: r16 }, { half: true }).every(sg => !(sg.x0 === 0 && sg.x1 === 0)));
    const left = parseAndEvaluate(U + 'gnd ngon x=-3 y=0 r=1 n=8').rects[0].shape;
    check('half: a shape left of the plane has no loops', shapeLoops(left, { half: true }).length === 0);
}

// --- Coax conversion ---
{
    const p = { tl_type: 'coax', coax_d: 0.92e-3, coax_D: 2.95e-3, coax_er: 2.1, coax_tand: 2e-4, coax_sigma: 5.8e7,
        rq: 0, freq: 1e9, mesh_backend: 'fullwave_mqs', use_plating: true, plating_sigma: 4.1e7, plating_t: 2e-6,
        plating_rq: 0, coax_plating_inner: 1, coax_plating_outer: 0 };
    const native = buildSolverFromParams(p, () => {});
    const text = solverToGeometryText(native, { units: 'mm' });
    let custom = null;
    try { custom = new CustomGeometrySolver({ text }); } catch (e) { check('coax: converted text builds', false, e.message); }
    if (custom) {
        check('coax: converted text builds', true);
        const inner = custom.conductors.find(c => c.is_signal), shield = custom.conductors.find(c => !c.is_signal);
        const area = o => shapeArea(o);
        const relNear = (a, b) => Math.abs(a / b - 1) < 1e-5;
        check('coax: n-gons carry the circle areas', relNear(area(inner), Math.PI * 0.46e-3 ** 2)
            && relNear(shapeArea({ shape: { type: 'polygon', poly: shield.shape.hole } }), Math.PI * 1.475e-3 ** 2));
        check('coax: the shield keeps the open-boundary warning away', custom.openBoundaryWarnings().length === 0,
            custom.openBoundaryWarnings().join(' '));
        check('coax: plating on the centre conductor only', !!(inner.plating && inner.plating.all) && !shield.plating);
        check('coax: dielectric material', custom.dielectrics[0].epsilon_r === 2.1 && custom.dielectrics[0].tan_delta === 2e-4);
    }
}

console.log(failures ? `\n${failures} check(s) failed` : '\nAll checks passed');
process.exit(failures ? 1 : 0);
