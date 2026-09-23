// Custom geometry text: parser, expressions, units, pinned edges, round trip, and the
// rejection cases of CustomGeometrySolver. No solves, runs in well under a second.
import { parseGeometryText, evaluateGeometry, parseAndEvaluate, serializeGeometry,
    evaluateExpression, setParamInText, renameParamInText, setStatementInText, replaceStatementInText,
    insertLineInText, moveRectInText, rectStatementText, plausibilityWarnings, formatLength,
    solverToGeometryText } from '../src/custom_geometry_text.js';
import { CustomGeometrySolver } from '../src/custom_geometry.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
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
    ['sig+ x=0 y=0 w=1 h=1 color=red', 1, 'unknown'],
    ['a = b\nb = 1', 1, 'defined below'],
    ['a = zz', 1, "unknown parameter 'zz'"],
    ['w = w + 1', 1, 'refers to itself'],
    ['a = 1\na = 2', 2, 'twice'],
    ['inf = 3', 1, 'reserved'],
    ['units mm\nunits um', 2, 'more than one'],
    ['\n\nwibble 3', 3, 'unknown statement'],
    ['sig+ x=0 y=0 w=1 h=1 plating=left', 1, 'faces'],
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

console.log(failures === 0 ? '\nALL CUSTOM GEOMETRY TEXT TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
