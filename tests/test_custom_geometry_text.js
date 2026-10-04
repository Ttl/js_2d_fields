// Custom geometry text without solving: parser, expressions, units, pinned edges,
// round trip, mirror=1, the form's expression arithmetic, trapezoids and n-gons (the
// polygons the text makes, validation, the quasi-static refusal, the coax conversion and
// the polygon helpers in shapes.js), and the rejection cases of CustomGeometrySolver.
// Runs in well under a second.
import { parseGeometryText, evaluateGeometry, parseAndEvaluate, serializeGeometry,
    evaluateExpression, setParamInText, renameParamInText, setStatementInText, replaceStatementInText,
    insertLineInText, moveRectInText, rectStatementText, plausibilityWarnings, formatLength,
    solverToGeometryText, changeUnitsInText, isReservedName, isPlainNumber, isLengthLiteral,
    axisEdges, addExpr } from '../src/custom_geometry_text.js';
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { buildSolverFromParams } from '../src/solver_factory.js';
import { isXSymmetric, conductorSwapSymmetric, halfDomainSymmetry } from '../src/geometry_symmetry.js';
import { bodyDistance, shapeLoops, shapeArea, shapeSegments, shapeContains, shapeFaceAt, platingArea, platedThrough, insideRingHole } from '../src/shapes.js';
import { check, done } from './helpers.js';

const close = (a, b, rel = 1e-12) => Math.abs(a - b) <= rel * Math.max(Math.abs(a), Math.abs(b), 1e-300);

// --- Expressions ---
check('precedence and parentheses', evaluateExpression('1+2*3-(4-6)/2') === 8);
check('unary minus', evaluateExpression('-2*-3') === 6 && evaluateExpression('-(1+1)') === -2);
check('parameters', evaluateExpression('s/2+w', { s: 3, w: 0.5 }) === 2);
check('functions', evaluateExpression('max(1,min(5,3))+abs(-2)+sqrt(16)') === 9);
check('exponent notation', evaluateExpression('2.5e-3*2') === 0.005);
check('unit suffix converts to the declared unit', close(evaluateExpression('35um+1', {}, 1e-3), 1.035));
check('spaced micro-sign unit converts', close(evaluateExpression('17.5 µm', {}, 1e-3), 0.0175));
check('inf', evaluateExpression('-inf') === -Infinity);
for (const bad of ['1+', '2*(3', 'foo', '1 2', '3$', 'nofn(1)', '5furlong', 'inf-inf']) {
    let threw = false;
    try { evaluateExpression(bad); } catch { threw = true; }
    check(`rejects '${bad}'`, threw);
}

// --- A complete geometry ---
const TEXT = `
units um
w = 200; s = 150; t = 35   # trace
h1 = 100; h2 = 4*h1; er2 = 3.0
bounds open open open gnd
domain auto

diel  x=-inf    y=0      w=inf  h=h1   er=4.3 tand=0.02
diel  x=-inf    y=h1     w=inf  h=h2   er=er2   # second layer
sig-  x=-s/2-w  y=h1+h2  w=w    h=t
sig+  x=s/2     y=h1+h2  w=w    h=t
`;
const model = parseGeometryText(TEXT);
check('parses without errors', model.errors.length === 0, JSON.stringify(model.errors));
const geo = evaluateGeometry(model);
check('evaluates without errors', geo.errors.length === 0, JSON.stringify(geo.errors));
check('parameters evaluated in order', geo.params.h2 === 400 && geo.params.er2 === 3);
check('lengths in metres', close(geo.rects[2].x.pos, -275e-6) && close(geo.rects[2].y.pos, 500e-6)
    && close(geo.rects[2].x.size, 200e-6));
check('er is not scaled by the unit', geo.rects[1].er === 3);
check('second trace', close(geo.rects[3].x.min, 75e-6) && close(geo.rects[3].x.max, 275e-6));
check('pinned edges', geo.rects[0].x.min === -Infinity && geo.rects[0].x.max === Infinity);
check('bounds and auto domain', geo.bounds.join() === 'open,open,open,gnd' && geo.domain.x_min === null);

const swept = evaluateGeometry(model, { s: 300 });
check('override replaces a parameter', close(swept.rects[3].x.min, 150e-6) && swept.params.s === 300);
check('override of an unknown parameter is an error', evaluateGeometry(model, { nope: 1 }).errors.length === 1);

const text2 = serializeGeometry(model);
const geo2 = parseAndEvaluate(text2);
check('round trip evaluates identically',
    JSON.stringify(geo2.rects.map(r => [r.kind, r.x, r.y, r.er])) === JSON.stringify(geo.rects.map(r => [r.kind, r.x, r.y, r.er])));
check('round trip keeps comments', text2.includes('# second layer') && text2.includes('# trace'));
check('round trip is stable', serializeGeometry(parseGeometryText(text2)) === text2);

// In-place edits keep the rest of the text as written.
{
    const edited = setParamInText(TEXT, 's', '300');
    check('setParamInText rewrites one statement', edited.includes('w = 200; s = 300; t = 35   # trace')
        && parseAndEvaluate(edited).params.s === 300 && edited.split('\n').length === TEXT.split('\n').length);
    check('setParamInText leaves other parameters alone', setParamInText(TEXT, 'h', '5') === TEXT
        && parseAndEvaluate(setParamInText(TEXT, 'h1', '50')).params.h2 === 200);
    const b = setStatementInText(TEXT, 'bounds', 'bounds gnd gnd gnd gnd');
    check('setStatementInText replaces a statement', b.includes('bounds gnd gnd gnd gnd') && !b.includes('open')
        && parseAndEvaluate(b).bounds.join() === 'gnd,gnd,gnd,gnd');
    const ins = setStatementInText('units mm\nsig+ x=0 y=1 w=1 h=1', 'bounds', 'bounds open open open gnd');
    check('setStatementInText inserts after units', ins.split('\n')[1] === 'bounds open open open gnd');
}

// Statement-level edits used by the form editor.
{
    const m = parseGeometryText(TEXT);
    const sigN = m.statements.find(s => s.type === 'rect' && s.kind === 'sig-');
    const r1 = replaceStatementInText(TEXT, sigN, rectStatementText('sig-', { ...sigN.fields, w: '2*w' }));
    check('replaceStatementInText rewrites one rectangle', close(parseAndEvaluate(r1).rects[2].x.size, 400e-6)
        && r1.split('\n').length === TEXT.split('\n').length);
    const pS = m.statements.find(s => s.type === 'param' && s.name === 's');
    const r2 = replaceStatementInText(TEXT, pS, 's = 99');
    check('replaceStatementInText edits one of several statements on a line',
        /w = 200; s = 99; t = 35   # trace/.test(r2));
    const r3 = replaceStatementInText(TEXT, sigN, null);
    check('removing a statement drops its line', parseAndEvaluate(r3).rects.length === 3
        && r3.split('\n').length === TEXT.split('\n').length - 1);
    const r4 = replaceStatementInText(TEXT, pS, null);
    check('removing one of several statements keeps the others', /w = 200; t = 35   # trace/.test(r4), r4.split('\n')[2]);
    const r5 = insertLineInText(TEXT, 1e9, 'gnd  x=0 y=-10 w=5 h=5');
    check('insertLineInText appends before the final newline', r5.endsWith('gnd  x=0 y=-10 w=5 h=5\n')
        && parseAndEvaluate(r5).rects.length === 5);
    const d2 = m.statements.filter(s => s.type === 'rect')[1];
    const r6 = moveRectInText(TEXT, m, d2, -1);
    const e6 = parseAndEvaluate(r6);
    check('moveRectInText swaps two rectangles', e6.errors.length === 0 && e6.rects[0].er === 3 && e6.rects[1].er === 4.3
        && r6.includes('# second layer'));
}

// Shared edges written as the same expression are the same double.
check('shared edges are exactly equal', geo.rects[0].y.max === geo.rects[1].y.min);

// --- Parse and evaluation errors carry the line ---
const ERR = [
    ['units parsec', 1, 'units'],
    ['w = ', 1, null],
    ['bounds open open gnd', 1, 'bounds'],
    ['sig+ x=0 y=0 w=1', 1, 'h'],
    ['sig+ x=0 x1=0 w=1 y=0 h=1', 1, 'unknown rectangle key'],
    ['sig+ x=0 y=0 w=1 h=0', 1, 'nonzero'],
    ['sig+ x=0 y=0 w=1 h=1 er=2', 1, 'diel only'],
    ['diel x=0 y=0 w=1 h=1', 1, 'er'],
    ['diel x=0 y=0 w=1 h=1 er=0.5', 1, 'er'],
    ['sig+ x=-inf y=0 w=1 h=1', 1, 'w=inf'],
    ['sig+ x=0 y=0 w=1 h=1 colour=red', 1, 'unknown'],
    ['sig+ x=0 y=0 w=1 h=1 color=red', 1, '#rrggbb'],
    ['a = b\nb = 1', 1, 'defined below'],
    ['a = zz', 1, "unknown parameter 'zz'"],
    ['w = w + 1', 1, 'refers to itself'],
    ['a = 1\na = 2', 2, 'twice'],
    ['inf = 3', 1, 'reserved'],
    ['units mm\nunits um', 2, 'more than one'],
    ['\n\nwibble 3', 3, 'unknown statement'],
    ['sig+ x=0 y=0 w=1 h=1 plating=left', 1, 'faces'],
    ['sig+ x=0 y=0 w=1 h=1 plating=top plating_sigma=inf plating_t=0.01', 1, 'finite'],
    ['sig+ x=0 y=0 w=1 h=1 plating=top plating_sigma=1e7 plating_t=inf', 1, 'finite'],
    ['sig+ x=0 y=0 w=1 h=1 plating=top plating_sigma=1e7 plating_t=0.01 plating_rq=inf', 1, 'finite'],
    ['sig+ x=0 y=0 w=1 h=1 plating=top plating_sigma=1e7 plating_t=0.01 plating_rq_iface=inf', 1, 'finite'],
    ['plating sigma=inf t=0.01', 1, 'finite'],
    ['plating sigma=1e7 t=0.01 rq=inf', 1, 'finite'],
];
// A rectangle using a parameter that failed names that parameter, not an unknown one.
{
    const e = parseAndEvaluate('q = 1/0\nsig+ x=0 y=0 w=q h=1').errors;
    check('error: rectangle using a failed parameter', e.length === 2 && e[1].line === 2
        && e[1].message.includes("parameter 'q' has an error"), JSON.stringify(e));
}
for (const [text, line, word] of ERR) {
    const g = parseAndEvaluate(text);
    const e = g.errors[0];
    check(`error: ${JSON.stringify(text)}`,
        !!e && e.line === line && (word === null || e.message.includes(word)), e ? `line ${e.line}: ${e.message}` : 'no error');
}

// --- Plausibility of material values ---
{
    const text = `units mm
plating sigma=4e7 t=5 rq=0
sig+ x=0 y=0 w=0.3 h=0.035 rq=1 sigma=5.8e4
sig- x=1 y=0 w=0.3 h=0.035 rq=1um plating=top plating_t=0.05
gnd x=-1 y=-1 w=3 h=0.5 rq=0.0005 plating=top plating_sigma=5.8e7
`;
    const model = parseGeometryText(text);
    const w = plausibilityWarnings(evaluateGeometry(model), model);
    const on = line => w.filter(o => o.line === line).map(o => o.message);
    check('plausibility: plating statement t=5 under mm is flagged', on(2).length === 1 && on(2)[0].includes('plating t = 5 mm'), JSON.stringify(on(2)));
    check('plausibility: bare rq=1 under mm and a low sigma are flagged', on(3).length === 2
        && on(3).some(m => m.includes('rq = 1 mm')) && on(3).some(m => m.includes('S/m')), JSON.stringify(on(3)));
    check('plausibility: rq=1um and plating thicker than the trace are fine', on(4).length === 0, JSON.stringify(on(4)));
    check('plausibility: typical values are fine', on(5).length === 0, JSON.stringify(on(5)));
    check('formatLength picks a unit', [1, 1e-3, 35e-6, 5e-8].map(formatLength).join() === '1 m,1 mm,35 µm,50 nm');
}

// --- Solver construction and validation ---
const build = (text, extra = {}) => new CustomGeometrySolver({ text, nx: 20, ny: 20, ...extra });
const rejects = (name, text, word, extra) => {
    let msg = null;
    try { build(text, extra); } catch (e) { msg = e.message; }
    check(`rejects: ${name}`, msg !== null && msg.includes(word), msg ?? 'accepted');
};

const s1 = build(TEXT);
check('differential pair recognised', s1.is_differential === true && s1.conductors.filter(c => c.is_signal).length === 2);
check('gnd wall becomes a ground slab below the substrate',
    s1.conductors.some(c => !c.is_signal && close(c.y_max, 0) && close(c.y_min, -35e-6))
    && close(s1.domain_y_min, -35e-6));
check('pinned dielectric spans the domain',
    close(s1.dielectrics[0].x_min, -s1.domain_width / 2) && close(s1.dielectrics[0].x_max, s1.domain_width / 2));
check('symmetric geometry uses the half domain', s1.sym_half === true);
check('sizing hints', close(s1.t, 35e-6) && close(s1.w, 200e-6));

// Off-centre input is centred on x=0 and then solves on the half domain.
const s2 = build(TEXT.replace('x=-s/2-w', 'x=1000-s/2-w').replace('x=s/2 ', 'x=1000+s/2 '));
check('off-centre geometry is recentred', close(s2.x_shift, 1000e-6, 1e-9) && s2.sym_half === true
    && close(s2.conductors.find(c => c.polarity > 0).x_min, 75e-6, 1e-9));

const s3 = build(`units mm\nbounds open open open open\ndiel x=-inf w=inf y=0 h=0.5 er=4\n` +
    `gnd x=-2 w=4 y=-0.035 h=0.035\nsig+ x=-0.2 w=0.4 y=0.5 h=0.035`);
check('all-open box with a finite ground', s3.domain_y_min < -0.035e-3 && s3.conductors.length === 2
    && s3.dielectrics[0].y_min === 0);

const s4 = build(`units mm\nbounds open open open open\ndomain -5 5 -3 3\n` +
    `sig+ x=-0.1 w=-inf y=0 h=0.035\ngnd x=0.1 w=inf y=0 h=0.035\ndiel x=-inf w=inf y=-0.6 h=0.6 er=9.8`);
check('slotline: half planes pinned to the walls', close(s4.conductors[0].x_min, -5e-3) && close(s4.conductors[1].x_max, 5e-3)
    && s4.w === undefined && s4.is_differential === false);

const sEnc = build(`units mm\nbounds gnd gnd gnd gnd\ndiel x=-2 w=4 y=-0.2 h=0.2 er=4.4\n` +
    `sig+ x=-1.15 w=1 y=0 h=0.05\ngnd x=0.15 w=1 y=0 h=0.05`);
check('auto gnd wall clears a conductor on the stack edge', sEnc.user_domain.y_max > 1e-3
    && close(sEnc.user_domain.y_min, -0.2e-3));

const RECT = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.5 er=4\n';
rejects('no signal', RECT + 'gnd x=0 y=1 w=1 h=0.1', 'signal');
rejects('no ground', 'units mm\nsig+ x=0 y=0 w=1 h=0.1', 'ground');
rejects('sig- without sig+', RECT + 'sig- x=0 y=0.5 w=1 h=0.1', 'sig+');
rejects('signal overlapping ground', RECT + 'sig+ x=0 y=0.5 w=1 h=0.1\ngnd x=0.5 y=0.55 w=1 h=0.1', 'shorted');
rejects('signal touching ground', RECT + 'sig+ x=0 y=0.5 w=1 h=0.1\ngnd x=1 y=0.5 w=1 h=0.1', 'shorted');
rejects('signal on a gnd wall', RECT + 'sig+ x=0 y=0 w=1 h=0.1', 'shorted');
rejects('sig+ touching sig-', RECT + 'sig+ x=0 y=0.5 w=1 h=0.1\nsig- x=1 y=0.5 w=1 h=0.1', 'shorted');
rejects('conductor outside the domain', RECT + 'domain -1 1 0 2\nsig+ x=0.5 y=0.5 w=1 h=0.1', 'outside');
rejects('plating without a material', RECT + 'sig+ x=0 y=0.5 w=1 h=0.1 plating=top', 'plating');
rejects('plating on a joined conductor', RECT + 'plating sigma=1e7 t=0.004\n' +
    'sig+ x=0 y=0.5 w=1 h=0.1 plating=top\nsig+ x=1 y=0.5 w=1 h=0.1', 'touches');
rejects('parse errors are reported with lines', RECT + 'sig+ x=0 y=0.5 w=oops h=0.1', 'line 4');

const s5 = build(RECT + 'sig+ x=0 y=0.5 w=1 h=0.1 plating=top', { plating: { sigma: 1e7, thickness: 4e-6, rq: 0 } });
check('plating material from the options', s5.conductors.find(c => c.is_signal).plating.top === true
    && s5.conductors.find(c => c.is_signal).plating.sigma === 1e7);

// Paint order and symmetry: mirrored overlapping dielectrics are only symmetric when
// their order mirrors too.
const OV = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.5 er=4\n';
const sig = 'sig+ x=-0.2 w=0.4 y=0.5 h=0.035\n';
const sOk = build(OV + 'diel x=-0.5 w=1 y=0.1 h=0.2 er=9\n' + sig);
check('symmetric inset keeps the half domain', sOk.sym_half === true);
const sBad = build(OV + 'diel x=-1 w=1.5 y=0.1 h=0.2 er=9\ndiel x=-0.5 w=1.5 y=0.1 h=0.2 er=2\n' + sig);
check('mirrored rectangles painted asymmetrically use the full domain', sBad.sym_half === false && sBad.tri_symmetry === false);

// --- Factory: the path the app and the worker use ---
{
    const { buildSolverFromParams } = await import('../src/solver_factory.js');
    const p = { tl_type: 'custom', custom_geom: TEXT, custom_overrides: { s: 300 }, sigma: 4e7, freq: 2e9,
        nx: 20, ny: 20, rq: 1e-7, use_plating: false, use_causal_materials: true, mesh_backend: 'fullwave_mqs' };
    const s = buildSolverFromParams(p);
    check('factory builds a custom solver', !!s && s.sigma_cond === 4e7 && s.freq === 2e9 && s.rq === 1e-7
        && s.mesh_backend === 'triangular' && s.geometry_params.s === 300);
    let msg = null;
    const bad = buildSolverFromParams({ ...p, custom_geom: 'sig+ x=0' }, m => { msg = m; });
    check('factory reports geometry errors', bad === null && /line 1/.test(msg ?? ''), msg ?? '');
}

// --- A space between a number and its unit ---
{
    const g = parseAndEvaluate('units mm\nr = 1 um\nsig+ x=0 y=0 w=1 h=35 um rq=r plating=top plating_sigma=1e7 plating_t=2 um');
    const c = g.rects[0];
    check('a unit may follow its number after a space', g.errors.length === 0 && Math.abs(c.rq - 1e-6) < 1e-18
        && Math.abs(c.y.size - 35e-6) < 1e-15 && Math.abs(c.platingMaterial.thickness - 2e-6) < 1e-18,
        JSON.stringify(g.errors));
    check('a stray word is still an error', parseAndEvaluate('sig+ x=0 y=0 w=1 h=1 um2').errors.length === 1
        && parseAndEvaluate('sig+ x=0 y=0 w=w1 um h=1').errors.length >= 1);
}

// --- Renaming a parameter ---
{
    const src = 'units mm\nw = 0.3; w2 = 2*w   # w is the width\nsig+ x=-w/2 w=w y=0 h=0.035\ngnd x=-w2 w = w2 y=-1 h=0.5w';
    const out = renameParamInText(src, 'w', 'wt');
    check('renameParamInText renames the definition and its uses, not keys, comments or other names',
        out === 'units mm\nwt = 0.3; w2 = 2*wt   # w is the width\nsig+ x=-wt/2 w=wt y=0 h=0.035\ngnd x=-w2 w = w2 y=-1 h=0.5w', out);
    const a = parseAndEvaluate(src), b = parseAndEvaluate(out);
    check('the renamed geometry evaluates the same', b.errors.length === a.errors.length
        && JSON.stringify(b.rects.map(r => [r.x, r.y])) === JSON.stringify(a.rects.map(r => [r.x, r.y])));
}
// A parameter may share its name with a statement word, a plating face or a shape:
// renaming it changes only expressions.
{
    const src = 'units mm\ntop = 0.3; gnd = 1; domain = 3; ngon = 0.1\nbounds open open open gnd\n'
        + 'domain -domain domain auto auto\nplating sigma=1e7 t=top/100\n'
        + 'diel x=-inf w=inf y=0 h=gnd er=4\nsig+ x=-top/2 w=top y=gnd h=0.035 plating=top,sides color=#abc\n'
        + 'gnd ngon x=0 y=-2 r=ngon n=8\ngnd x=-gnd w=2*gnd y=-0.035 h=0.035';
    const ok = (from, to) => {
        const out = renameParamInText(src, from, to);
        const a = parseAndEvaluate(src), b = parseAndEvaluate(out);
        return { out, same: a.errors.length === 0 && b.errors.length === 0
            && JSON.stringify(b.rects.map(r => [r.x, r.y, r.kind])) === JSON.stringify(a.rects.map(r => [r.x, r.y, r.kind])) };
    };
    const top = ok('top', 'wt'), gnd = ok('gnd', 'hs'), dom = ok('domain', 'dw'), ng = ok('ngon', 'rr');
    check('renaming a parameter named top keeps plating=top,sides', top.same && /plating=top,sides/.test(top.out)
        && /t=wt\/100/.test(top.out), top.out);
    check('renaming a parameter named gnd keeps the gnd statements and bounds', gnd.same
        && /bounds open open open gnd/.test(gnd.out) && /^gnd x=-hs w=2\*hs/m.test(gnd.out), gnd.out);
    check('renaming a parameter named domain renames the domain values, not the statement', dom.same
        && /^domain -dw dw auto auto/m.test(dom.out), dom.out);
    check('renaming a parameter named ngon keeps the shape word', ng.same && /gnd ngon x=0 y=-2 r=rr/.test(ng.out), ng.out);
    const d = setStatementInText(src, 'domain', 'domain -5 5 auto auto');
    check('setStatementInText leaves a parameter named domain alone', /domain = 3/.test(d)
        && /^domain -5 5 auto auto/m.test(d) && parseAndEvaluate(d).errors.length === 0, d);
}

// Conversion writes roughness and plating thickness in micrometres under any declared
// unit, and the text rebuilds the same values.
{
    const src = new CustomGeometrySolver({ nx: 20, ny: 20, text: `units mm
plating sigma=4.1e7 t=0.005 rq=0.0002
bounds open open open gnd
diel x=-inf w=inf y=0 h=0.2 er=4
sig+ x=-0.15 w=0.3 y=0.2 h=0.035 rq=0.0015 plating=top,sides
` });
    const text = solverToGeometryText(src, { units: 'mm', pinWalls: true });
    const back = new CustomGeometrySolver({ nx: 20, ny: 20, text });
    const sig = c => c.find(o => o.is_signal);
    check('conversion writes rq and plating t in um', /\brq=1\.5um\b/.test(text) && /\bt=5um\b/.test(text) && /\brq=0\.2um\b/.test(text), text);
    check('conversion in um rebuilds the same rq and plating', close(sig(back.conductors).rq, 1.5e-6, 1e-9)
        && close(back.plating.thickness, 5e-6, 1e-9) && close(back.plating.rq, 0.2e-6, 1e-9));
}

// A plating switched on with zero thickness has no effect: conversion writes no
// plating instead of failing on the missing plating statement.
{
    const { MicrostripSolver } = await import('../src/microstrip.js');
    const { CoaxSolver } = await import('../src/coax.js');
    const plating = { sigma: 1e7, thickness: 0, rq: 0, top: true, sides: true, bottom: false, inner: true };
    let text = null, err = null;
    try {
        text = solverToGeometryText(new MicrostripSolver({ substrate_height: 0.2e-3, trace_width: 0.3e-3,
            trace_thickness: 35e-6, epsilon_r: 4.3, plating }), { units: 'mm', pinWalls: true });
    } catch (e) { err = e.message; }
    check('conversion with a zero-thickness plating writes no plating', text !== null && !/plating/.test(text)
        && parseAndEvaluate(text).errors.length === 0, err || text);
    let coax = null;
    try { coax = new CoaxSolver({ inner_diameter: 0.6e-3, dielectric_diameter: 2e-3, epsilon_r: 2.1, plating }); } catch (e) { err = e.message; }
    check('a coax with a zero-thickness plating builds unplated', coax && coax.conductors.every(c => !c.plating), err);
}

// --- Change of units ---
{
    const text = `units mm
w = 0.3; s = 0.2; t = 35 um   # widths
n = 2; h1 = 0.1; h2 = h1*n + 0.05
plating sigma=4e7 t=0.005 rq=0.0002
bounds open open open gnd
domain -3 3 auto 2
diel  x=-inf    y=0      w=inf  h=h1  er=4.3  tand=0.02
diel  x=-inf    y=h1     w=inf  h=h2  er=2.2  tand=0.001
sig-  x=-s/2-w  y=h1+h2  w=w    h=t  rq=1 um
sig+  x=s/2     y=h1+h2  w=w    h=t  rq=0.001; gnd x=-1 w=max(2*w,0.5)+0.1 y=-0.2 h=0.1
`;
    const um = changeUnitsInText(text, 'um');
    const a = parseAndEvaluate(text), b = um && parseAndEvaluate(um);
    const same = (p, q) => p.rects.length === q.rects.length && p.rects.every((r, i) =>
        ['min', 'max'].every(k => [[r.x[k], q.rects[i].x[k]], [r.y[k], q.rects[i].y[k]]].every(([u, v]) => u === v || close(u, v, 1e-9)))
        && (r.rq === null ? q.rects[i].rq === null : close(r.rq, q.rects[i].rq, 1e-9)));
    check('unit change mm -> um keeps the geometry', !!um && b.units === 'um' && same(a, b), um ?? 'null');
    check('unit change scales lengths, not factors or suffixed numbers',
        !!um && /^w = 300; s = 200; t = 35 um {3}# widths$/m.test(um) && /^n = 2; h1 = 100; h2 = h1\*n \+ 50$/m.test(um)
        && /x=-s\/2-w/.test(um) && /rq=1 um/.test(um) && /w=max\(2\*w,500\)\+100/.test(um)
        && /^domain -3000 3000 auto 2000$/m.test(um) && /^plating sigma=4e7 t=5 rq=0\.2$/m.test(um), um ?? 'null');
    const mil = changeUnitsInText(text, 'mil');
    check('unit change mm -> mil keeps the geometry', !!mil && same(a, parseAndEvaluate(mil)));
    check('unit change without a units line adds one', /^units um$/m.test(changeUnitsInText('w = 0.3\nsig+ x=0 w=w y=0 h=0.035', 'um') ?? ''));
    check('unit change of a text with errors gives null', changeUnitsInText('units mm\nw = oops', 'um') === null);
    {
        // k is a length in x=k and a factor in w=k*p0: no consistent rewrite, and with 14
        // parameters and 300 shapes the full search would take tens of seconds.
        const params = Array.from({ length: 13 }, (_, i) => `p${i} = 0.${i + 1}`).join('\n');
        const shapes = Array.from({ length: 300 }, (_, i) => `diel x=k+${i} y=p${i % 13} w=k*p0 h=p1 er=4`).join('\n');
        const t0 = Date.now();
        const r = changeUnitsInText(`units mm\nk = 2\n${params}\n${shapes}`, 'um');
        const ms = Date.now() - t0;
        check('unit change without a consistent rewrite gives up within its time budget', r === null && ms < 2000, `${ms} ms`);
    }
}

// --- Names and literals ---
check('reserved names', isReservedName('mm') && isReservedName('inf') && isReservedName('sqrt') && !isReservedName('w'));
check('plain numbers and length literals', isPlainNumber('-1.5e-3') && !isPlainNumber('35um')
    && isLengthLiteral('35um') && isLengthLiteral('35 um') && isLengthLiteral('0.2') && !isLengthLiteral('2*w') && !isLengthLiteral('3 kg'));

// --- Signal drawn as separate bodies ---
{
    const sep = new CustomGeometrySolver({ nx: 10, ny: 10, text: `units mm
bounds open open open gnd
diel x=-inf w=inf y=0 h=0.2 er=4
sig+ x=-1 w=0.3 y=0.2 h=0.035
sig+ x=1 w=0.3 y=0.2 h=0.035
sig+ x=1.3 w=0.3 y=0.2 h=0.035
sig- x=0 w=0.2 y=0.2 h=0.035
` }).openBoundaryWarnings();
    check('separate sig+ bodies are reported once, sig- alone is not', sep.length === 1
        && sep[0].includes('sig+ conductor is 2 separate bodies'), JSON.stringify(sep));
}

// --- mirror=1, pinned edges and expression arithmetic ---
{
    const HEAD = `units mm
    w = 0.2; s = 0.15; g = 0.1; t = 0.035; h = 0.2; wv = 0.3
    bounds open open open gnd
    `;
    const MIRRORED = HEAD + `
    diel x=0 w=inf y=0 h=h er=4.4 tand=0.02 mirror=1
    sig+ x=s/2 w=w y=h h=t mirror=1
    gnd  x=s/2+w+g w=inf y=h h=t mirror=1
    gnd  x=s/2+w+2*g w=wv y=0 h=h mirror=1
    `;
    const HAND = HEAD + `
    diel x=-inf w=inf y=0 h=h er=4.4 tand=0.02
    sig+ x=s/2 w=w y=h h=t
    sig- x=-s/2-w w=w y=h h=t
    gnd  x=s/2+w+g w=inf y=h h=t
    gnd  x=-s/2-w-g w=-inf y=h h=t
    gnd  x=s/2+w+2*g w=wv y=0 h=h
    gnd  x=-s/2-w-2*g-wv w=wv y=0 h=h
    `;
    const geoM = parseAndEvaluate(MIRRORED), geoH = parseAndEvaluate(HAND);
    check('mirrored text evaluates', geoM.errors.length === 0, JSON.stringify(geoM.errors));
    check('mirror emits the images', geoM.rects.length === 7 && geoM.rects.filter(r => r.image).length === 3);
    check('image follows its rectangle, sig+ imaged as sig-',
        geoM.rects[1].kind === 'sig+' && geoM.rects[2].kind === 'sig-' && geoM.rects[2].image && geoM.rects[2].line === geoM.rects[1].line);
    check('touching x=0 gives one full-width rectangle',
        geoM.rects[0].x.min === -Infinity && geoM.rects[0].x.max === Infinity && !geoM.rects[0].image);
    check('pinned edge mirrors to the other wall', geoM.rects[4].x.min === -Infinity && geoM.rects[4].x.max === -geoM.rects[3].x.min);

    const key = o => [o.is_signal ? o.polarity : 'x', o.x_min, o.x_max, o.y_min, o.y_max, o.epsilon_r ?? ''].join();
    const sM = new CustomGeometrySolver({ geometry: geoM, nx: 10, ny: 10 });
    const sH = new CustomGeometrySolver({ geometry: geoH, nx: 10, ny: 10 });
    const same = (a, b) => JSON.stringify(a.map(key).sort()) === JSON.stringify(b.map(key).sort());
    check('same conductors as the hand-written geometry, bit for bit', same(sM.conductors, sH.conductors));
    check('same dielectrics as the hand-written geometry, bit for bit', same(sM.dielectrics, sH.dielectrics));
    check('same domain', sM.domain_width === sH.domain_width && sM.domain_height === sH.domain_height);
    check('half domain used on both', sM.sym_half && sH.sym_half);
    check('images are marked for the preview',
        sM.conductors.filter(c => c.src_image).length === 3 && sH.conductors.every(c => !c.src_image));

    // --- Other mirror cases ---
    const one = body => parseAndEvaluate('units mm\nbounds open open open gnd\n' + body);
    let r = one('sig+ x=0 w=0.1 y=0 h=0.1 mirror=1').rects;
    check('signal touching x=0 stays one sig+', r.length === 1 && r[0].kind === 'sig+'
        && Math.abs(r[0].x.min + 1e-4) < 1e-18 && Math.abs(r[0].x.max - 1e-4) < 1e-18);
    r = one('gnd x=-0.1 w=0.3 y=0 h=0.1 mirror=1').rects;
    check('crossing x=0 spans the wider side', r.length === 1 && Math.abs(r[0].x.min + 2e-4) < 1e-18 && Math.abs(r[0].x.max - 2e-4) < 1e-18);
    r = one('sig- x=-0.3 w=0.1 y=0 h=0.1 mirror=1').rects;
    check('sig- on the left imaged as sig+ on the right', r.length === 2 && r[1].kind === 'sig+'
        && r[1].x.min === -r[0].x.max && r[1].x.max === -r[0].x.min);
    r = one('gnd x=-0.5 w=-inf y=0 h=0.1 mirror=1').rects;
    check('left half-plane mirrors to a right half-plane', r[1].x.min === 5e-4 && r[1].x.max === Infinity);
    check('mirror=0 adds nothing', one('sig+ x=0.1 w=0.1 y=0 h=0.1 mirror=0').rects.length === 1);
    const sP = new CustomGeometrySolver({ geometry: one('diel x=0.1 w=0.2 y=0 h=0.1 er=3 mirror=1\ndiel x=0.2 w=0.2 y=0 h=0.1 er=5 mirror=1\nsig+ x=-0.05 w=0.1 y=0.1 h=0.01'), nx: 10, ny: 10 });
    check('overlapping mirrored dielectrics keep a symmetric paint order', sP.sym_half);
    const rt = parseAndEvaluate(serializeGeometry(parseGeometryText(MIRRORED)));
    check('mirror survives serialization', rt.rects.length === 7);

    // --- Expression arithmetic ---
    check('add', addExpr('s/2', 'w') === 's/2+w' && addExpr('0', 'w') === 'w' && addExpr('a', '-b') === 'a-b'
        && addExpr('0.1', '0.2') === '0.3');
    const vars = { a: 1.5, b: 0.25, c: 3, d: 7, s: 0.15, w: 0.2 };
    for (const [e1, e2] of [['a', 'b-c'], ['a*b', '-c/d'], ['-(a+b)', 'c'], ['a', 'b*-c']]) {
        check(`addExpr(${e1}, ${e2}) evaluates right`, Math.abs(evaluateExpression(addExpr(e1, e2), vars)
            - (evaluateExpression(e1, vars) + evaluateExpression(e2, vars))) < 1e-12);
    }

    // --- Pinned edges ---
    const edges = (f, axis, neg) => JSON.stringify(axisEdges(f, axis, neg));
    check('edges of x,w', edges({ x: 's/2', w: 'w' }, 'x') === '{"lo":"s/2","hi":"s/2+w"}');
    check('edges of a negative h', edges({ y: 'h1+h2', h: '-t' }, 'y', true) === '{"lo":"h1+h2-t","hi":"h1+h2"}');
    check('edges of a pinned axis', edges({ x: '-inf', w: 'inf' }, 'x') === '{"lo":"-inf","hi":"inf"}');
    check('edges of w=inf', edges({ x: '0.1', w: 'inf' }, 'x') === '{"lo":"0.1","hi":"inf"}');
    check('edges of w=-inf', edges({ x: '-0.1', w: '-inf' }, 'x') === '{"lo":"-inf","hi":"-0.1"}');
    r = one('gnd x=-0.5 w=-inf y=0 h=0.1\ndiel x=-inf w=inf y=0 h=-inf er=4').rects;
    check('w=-inf runs to the left wall', r[0].x.min === -Infinity && r[0].x.max === -5e-4);
    check('h=-inf runs to the bottom wall', r[1].y.min === -Infinity && r[1].y.max === 0);
    r = one('sig- x=-0.1 w=-0.2 y=0 h=0.1 mirror=1\ngnd x=0.5 w=0.1 y=0.1 h=-0.1').rects;
    check('negative w puts x at the right edge', Math.abs(r[0].x.min + 3e-4) < 1e-18 && r[0].x.max === -1e-4
        && r[0].x.size > 0 && r[0].x.flipped);
    check('negative w mirrors', r[1].kind === 'sig+' && r[1].x.min === 1e-4 && Math.abs(r[1].x.max - 3e-4) < 1e-18);
    const sigOf = body => new CustomGeometrySolver({ geometry: one('diel x=-inf w=inf y=0 h=0.2 er=4\n' + body), nx: 10, ny: 10 })
        .conductors.find(c => c.is_signal);
    const [cf, cp] = [sigOf('sig+ x=0.3 w=-0.2 y=0.2 h=0.035'), sigOf('sig+ x=0.1 w=0.2 y=0.2 h=0.035')];
    check('flipped x solves like the same rectangle written with positive w',
        cf.width > 0 && Math.abs(cf.x_min - cp.x_min) < 1e-15 && Math.abs(cf.x_max - cp.x_max) < 1e-15);
    check('edges of a negative w', edges({ x: 'a', w: '-w' }, 'x', true) === '{"lo":"a-w","hi":"a"}');
    for (const [body, msg] of [['sig+ x=inf w=1 y=0 h=1', 'cannot be inf'], ['sig+ x=0 y=0 h=1', 'must both be given']]) {
        const errs = one(body).errors;
        check(`rejects '${body}'`, errs.length === 1 && errs[0].message.includes(msg), errs[0]?.message);
    }
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
    // A long ellipse keeps its ends: vertices on x = +-rx, area close to pi rx ry.
    const long = polyOf(U + 'diel ellipse x=0 y=0 rx=10 ry=1 n=32 er=2');
    const xs = Array.from(long.poly).filter((_, i) => i % 2 === 0);
    check('ellipse: a long ellipse keeps its ends', Math.max(...xs) === 10e-3 && Math.min(...xs) === -10e-3
        && Math.abs(shapeArea({ shape: long }) / (Math.PI * 10e-6) - 1) < 0.01, `area ${(shapeArea({ shape: long }) / (Math.PI * 10e-6)).toFixed(4)} of the ellipse`);
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

// --- Review fixes ---
{
    const convex = p => {
        const n = p.length >> 1;
        for (let i = 0; i < n; i++) {
            const a = i, b = (i + 1) % n, c = (i + 2) % n;
            const cr = (p[2 * b] - p[2 * a]) * (p[2 * c + 1] - p[2 * b + 1]) - (p[2 * b + 1] - p[2 * a + 1]) * (p[2 * c] - p[2 * b]);
            if (cr < 0) return false;
        }
        return true;
    };
    check('radius a hair over half a side stays convex', convex(polyOf(U + 'sig+ x=-2 w=4 y=0 h=1 radius=0.5000000004').poly));
    const tr = polyOf(U + 'sig+ trap x=-2 w=4 y=0 h=1 angle=30 radius=0.3');
    const bottom = tr.faces.map((f, i) => [f, i]).filter(([f]) => f === 'bottom');
    const lenOf = i => { const n = tr.faces.length, j = (i + 1) % n; return Math.hypot(tr.poly[2 * j] - tr.poly[2 * i], tr.poly[2 * j + 1] - tr.poly[2 * i + 1]); };
    const leftB = bottom.filter(([, i]) => tr.poly[2 * i] < 0).reduce((a, [, i]) => a + lenOf(i), 0);
    const rightB = bottom.filter(([, i]) => tr.poly[2 * ((i + 1) % tr.faces.length)] > 0).reduce((a, [, i]) => a + lenOf(i), 0);
    check('rounded corner faces split symmetrically', Math.abs(leftB - rightB) < 1e-12, `${leftB} vs ${rightB}`);
    check('radius=0, radius_bottom=0, wall=0 on a boundary-spanning rectangle is a plain rectangle',
        ['radius=0', 'radius_bottom=0', 'wall=0'].every(k => { const g = parseAndEvaluate(U + `gnd x=-inf w=inf y=-1 h=1 ${k}`); return !g.errors.length && g.rects[0].shape === null; }));
    const rm = parseAndEvaluate(U + 'sig+ x=0 w=0.3 y=0 h=0.035 radius=0.01 mirror=1');
    check('mirror=1 widens a rounded rectangle on x=0 to one symmetric rectangle', !rm.errors.length && rm.rects.length === 1
        && rm.rects[0].x.min === -0.3e-3 && rm.rects[0].x.max === 0.3e-3, rm.errors.map(e => e.message).join());
    // A rotated n-gon cut by the symmetry plane: no zero-length side.
    const rot = polyOf(U + 'sig+ ngon x=0 y=0.4 r=0.1 n=4 rot=90');
    const half = shapeLoops(rot, { half: true })[0];
    let minSide = Infinity;
    for (let i = 0; i < half.length / 2; i++) { const j = (i + 1) % (half.length / 2); minSide = Math.min(minSide, Math.hypot(half[2 * j] - half[2 * i], half[2 * j + 1] - half[2 * i + 1])); }
    check('half of a rotated n-gon has no zero-length side', minSide > 1e-9, `shortest ${minSide}`);
    // Plating of a trapezoid's top only: its area is the top face times the thickness.
    const pt = { shape: polyOf(U + 'sig+ trap x=-0.1 w=0.2 y=0 h=0.035 angle=0.0001'), x_min: -0.1e-3, x_max: 0.1e-3, y_min: 0, y_max: 0.035e-3 };
    const top = { sigma: 1e7, thickness: 5e-6, top: true, sides: false, bottom: false };
    check('plating area of a partly plated shape counts the plated faces', Math.abs(platingArea(pt, top) / (0.2e-3 * 5e-6) - 1) < 1e-3,
        `${platingArea(pt, top)} vs ${0.2e-3 * 5e-6}`);
    const narrow = { x_min: -4e-6, x_max: 4e-6, y_min: 0, y_max: 35e-6, width: 8e-6, height: 35e-6 };
    check('plating that fills the width is plating through', platedThrough(narrow, { sigma: 1e7, thickness: 5e-6, top: true, sides: true, bottom: false })
        && !platedThrough(narrow, { sigma: 1e7, thickness: 3e-6, top: true, sides: true, bottom: false }));
    const cx = parseAndEvaluate(U + 'gnd ngon x=0 y=0 r=0.7 r_in=0.625 n=64\nsig+ ngon x=0 y=0 r=0.46 n=64').rects;
    check('a round wire in a tight shield is inside its hole (b/a < sqrt 2)', insideRingHole(cx[0].shape, { shape: cx[1].shape }));
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
    const touch = parseAndEvaluate(U + 'diel ngon x=0.6 y=0 r=0.6 n=64 er=2 mirror=1');
    check('mirror: an n-gon touching x=0 from one side gets its image', !touch.errors.length && touch.rects.length === 2,
        touch.errors.map(e => e.message).join());
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

// --- Half domain of a differential pair: mirror partners must have opposite polarity ---
{
    const U = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4.3\n';
    const half = text => {
        const s = new CustomGeometrySolver({ text: U + text });
        return { s, ok: halfDomainSymmetry(s.conductors, s.dielectrics, s.domain_width, s.is_differential).ok };
    };
    // + - - +: each net maps onto itself, the odd mode is x-symmetric
    const inter = half('sig+ x=-1.1 w=0.2 y=0.2 h=0.035\nsig- x=-0.5 w=0.2 y=0.2 h=0.035\n' +
        'sig- x=0.3 w=0.2 y=0.2 h=0.035\nsig+ x=0.9 w=0.2 y=0.2 h=0.035\n');
    check('symmetry: interleaved + - - + pair is mirror symmetric', isXSymmetric(inter.s.conductors, inter.s.dielectrics, inter.s.domain_width));
    check('symmetry: interleaved + - - + pair gets no half domain', !inter.ok && !inter.s.sym_half);
    const swapped = half('sig+ x=-0.4 w=0.3 y=0.2 h=0.035\nsig- x=0.1 w=0.3 y=0.2 h=0.035\n');
    check('symmetry: pair with sig+ on the left keeps the half domain', swapped.ok && swapped.s.sym_half);
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

// --- Fill colors ---
{
    const text = 'diel x=-inf w=inf y=0 h=1 er=4 color=#1E3A78  # SiO2 # note\n'
        + 'gnd x=-inf w=inf y=0 h=-0.1 color=#ccc\nsig+ x=1 w=1 y=1 h=0.1 mirror=1 color=#c80';
    const model = parseGeometryText(text);
    const rects = model.statements.filter(s => s.type === 'rect');
    check('color: a # after = is the value, the next # the comment', model.errors.length === 0
        && rects[0].fields.color === '#1E3A78' && rects[0].comment === 'SiO2 # note', JSON.stringify(model.errors));
    const geo = evaluateGeometry(model);
    check('color: normalized to #rrggbb, mirror images keep it',
        JSON.stringify(geo.rects.map(r => r.color)) === JSON.stringify(['#1e3a78', '#cccccc', '#cc8800', '#cc8800']));
    const solver = new CustomGeometrySolver({ text });
    check('color: carried to the solver bodies', solver.dielectrics[0].color === '#1e3a78'
        && solver.conductors.map(c => c.color).join() === '#cccccc,#cc8800,#cc8800');
    const edited = replaceStatementInText(text, rects[1], rectStatementText('gnd', { ...rects[1].fields, color: '#abcdef' }));
    check('color: an edit keeps the colors and comments of the other lines', edited.split('\n')[0] === text.split('\n')[0]
        && /color=#abcdef/.test(edited.split('\n')[1]));
}

done();
