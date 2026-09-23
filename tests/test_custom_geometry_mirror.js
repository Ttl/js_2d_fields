// Custom geometry mirror=1 and the form's axis edits (x,w <-> x1,x2, expression
// arithmetic). No solves, runs in well under a second.
import { parseGeometryText, parseAndEvaluate, serializeGeometry, evaluateExpression,
    axisForm, axisEdges, toggleAxisForm, addExpr, subExpr, negateExpr } from '../src/custom_geometry_text.js';
import { CustomGeometrySolver } from '../src/custom_geometry.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}

// --- Mirror against the same geometry written out by hand ---
const HEAD = `units mm
w = 0.2; s = 0.15; g = 0.1; t = 0.035; h = 0.2; wv = 0.3
bounds open open open gnd
`;
const MIRRORED = HEAD + `
diel x1=0 x2=inf y=0 h=h er=4.4 tand=0.02 mirror=1
sig+ x=s/2 w=w y=h h=t mirror=1
gnd  x1=s/2+w+g x2=inf y=h h=t mirror=1
gnd  x=s/2+w+2*g w=wv y=0 h=h mirror=1
`;
const HAND = HEAD + `
diel x=-inf w=inf y=0 h=h er=4.4 tand=0.02
sig+ x=s/2 w=w y=h h=t
sig- x=-s/2-w w=w y=h h=t
gnd  x1=s/2+w+g x2=inf y=h h=t
gnd  x1=-inf x2=-s/2-w-g y=h h=t
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
let r = one('sig+ x1=0 x2=0.1 y=0 h=0.1 mirror=1').rects;
check('signal touching x=0 stays one sig+', r.length === 1 && r[0].kind === 'sig+'
    && Math.abs(r[0].x.min + 1e-4) < 1e-18 && Math.abs(r[0].x.max - 1e-4) < 1e-18);
r = one('gnd x=-0.1 w=0.3 y=0 h=0.1 mirror=1').rects;
check('crossing x=0 spans the wider side', r.length === 1 && Math.abs(r[0].x.min + 2e-4) < 1e-18 && Math.abs(r[0].x.max - 2e-4) < 1e-18);
r = one('sig- x=-0.3 w=0.1 y=0 h=0.1 mirror=1').rects;
check('sig- on the left imaged as sig+ on the right', r.length === 2 && r[1].kind === 'sig+'
    && r[1].x.min === -r[0].x.max && r[1].x.max === -r[0].x.min);
r = one('gnd x1=-inf x2=-0.5 y=0 h=0.1 mirror=1').rects;
check('left half-plane mirrors to a right half-plane', r[1].x.min === 5e-4 && r[1].x.max === Infinity);
check('mirror=0 adds nothing', one('sig+ x=0.1 w=0.1 y=0 h=0.1 mirror=0').rects.length === 1);
const sP = new CustomGeometrySolver({ geometry: one('diel x=0.1 w=0.2 y=0 h=0.1 er=3 mirror=1\ndiel x=0.2 w=0.2 y=0 h=0.1 er=5 mirror=1\nsig+ x=-0.05 w=0.1 y=0.1 h=0.01'), nx: 10, ny: 10 });
check('overlapping mirrored dielectrics keep a symmetric paint order', sP.sym_half);
const rt = parseAndEvaluate(serializeGeometry(parseGeometryText(MIRRORED)));
check('mirror survives serialization', rt.rects.length === 7);

// --- Expression arithmetic ---
check('negate flips every term', negateExpr('-s/2-w') === 's/2+w' && negateExpr('a-b*(c+d)') === '-a+b*(c+d)');
check('negate a number', negateExpr('0.25') === '-0.25' && negateExpr('1e-3') === '-0.001');
check('add', addExpr('s/2', 'w') === 's/2+w' && addExpr('0', 'w') === 'w' && addExpr('a', '-b') === 'a-b'
    && addExpr('0.1', '0.2') === '0.3');
check('subtract', subExpr('s/2+w', 's/2') === 'w' && subExpr('a', 'b-c') === 'a-(b-c)' && subExpr('x', 'x') === '0'
    && subExpr('0', 'a+b') === '-a-b');
const vars = { a: 1.5, b: 0.25, c: 3, d: 7, s: 0.15, w: 0.2 };
for (const [e1, e2] of [['a', 'b-c'], ['a*b', '-c/d'], ['-(a+b)', 'c'], ['a', 'b*-c']]) {
    check(`subExpr(${e1}, ${e2}) evaluates right`, Math.abs(evaluateExpression(subExpr(e1, e2), vars)
        - (evaluateExpression(e1, vars) - evaluateExpression(e2, vars))) < 1e-12);
    check(`addExpr(${e1}, ${e2}) evaluates right`, Math.abs(evaluateExpression(addExpr(e1, e2), vars)
        - (evaluateExpression(e1, vars) + evaluateExpression(e2, vars))) < 1e-12);
    check(`negateExpr(${e2}) evaluates right`, Math.abs(evaluateExpression(negateExpr(e2), vars)
        + evaluateExpression(e2, vars)) < 1e-12);
}

// --- Axis form switch and pinned edges ---
let f = { x: 's/2', w: 'w' };
check('x,w to x1,x2', toggleAxisForm(f, 'x') && f.x1 === 's/2' && f.x2 === 's/2+w' && f.x === undefined);
check('and back', toggleAxisForm(f, 'x') && f.x === 's/2' && f.w === 'w' && f.x1 === undefined);
f = { y: 'h1+h2', h: '-t' };
check('negative h to y1,y2', toggleAxisForm(f, 'y', true) && f.y1 === 'h1+h2-t' && f.y2 === 'h1+h2');
f = { x1: '-inf', x2: '2' };
check('-inf with a finite high edge has no x,w form', !toggleAxisForm(f, 'x') && axisForm(f, 'x') === 'bounds');
f = { x: '0.1', w: 'inf' };
check('x,w with w=inf to x1,x2', toggleAxisForm(f, 'x') && f.x1 === '0.1' && f.x2 === 'inf');
check('and back', toggleAxisForm(f, 'x') && f.x === '0.1' && f.w === 'inf');
check('edges of a pinned axis', JSON.stringify(axisEdges({ x: '-inf', w: 'inf' }, 'x')) === '{"lo":"-inf","hi":"inf"}');

console.log(failures === 0 ? '\nALL CUSTOM GEOMETRY MIRROR TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
