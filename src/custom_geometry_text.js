// Text format for custom cross-sections built from axis-aligned rectangles.
//
//   units um
//   w = 200; s = 150; t = 35; h1 = 100
//   bounds open open open gnd          # left right top bottom
//   domain auto                        # or: domain x1 x2 y1 y2, each may be auto
//   plating sigma=4.1e7 t=5 rq=0.1     # default plating material for rects with plating=
//                                      # rq_iface= plating/bulk roughness, default rq
//
//   diel  x=-inf    y=0   w=inf  h=h1  er=4.3 tand=0.02
//   sig-  x=-s/2-w  y=h1  w=w    h=t
//   sig+  x=s/2     y=h1  w=w    h=t   plating=top,sides
//
// Statements are separated by newlines or ';'. '#' starts a comment, except right after
// '=' where it starts a color value. A rectangle is
// given per axis as position and size, x,w and y,h. A negative size flips the rectangle
// to the other side of its position: x is then the right edge, y the top face. x=-inf
// w=inf spans the domain, w=inf runs from x to the right wall and w=-inf from x to the
// left wall (h the same way on y).
// Values in a rectangle, domain or plating statement are expressions without
// whitespace. A parameter definition takes the rest of the statement.
//
// Lengths are in the declared units (default mm). er, tand and sigma are plain numbers.
// A number may carry its own unit (35um or 35 um), which converts it to the declared units.
// A dielectric may conduct (sigma=, S/m, a doped silicon substrate): it is solved with the
// complex permittivity er*(1 - j*tand) - j*sigma/(omega*eps0), see conductive_dielectric.js.
// A conductor may carry its own conductivity (sigma=, S/m), surface roughness (rq=) and
// plating material (plating_sigma=, plating_t=, plating_rq=, plating_rq_iface=), which override the
// solver-wide values and the plating statement for that conductor. Whether plating is
// thick (a layer of metal the full-wave solver meshes) or thin (a layered surface
// impedance) is the solver option Model Thick Plating, the same for every conductor.
// thin=1 on a dielectric marks a thin sheet (a solder mask) whose faces the FDM mesher
// brackets with grid lines the way it does conductor faces.
// mirror=1 adds the mirror image about x=0 right after the rectangle, sig+ imaged as sig-
// and sig- as sig+. A rectangle that touches or crosses x=0 becomes one rectangle
// symmetric about x=0 instead.
// Dielectrics are painted in order, a later one overrides an earlier one where they
// overlap. Conductors override dielectrics.
// color=#rrggbb (or #rgb) on any rectangle sets its fill color in the plots.
//
// A shape word after the kind draws something other than a rectangle. These need the
// full-wave solver.
//   sig+  trap  x=-w/2 y=h w=w h=t angle=30 angle2=20
//     trapezoid on the base x..x+w at y (a PCB trace with an etch angle). The sides lean
//     in by angle (left side) and angle2 (right side, default angle) degrees from the
//     vertical, so a positive angle makes the face at y+h narrower and a negative one
//     wider. A negative h puts that face below y. Zero angles give a plain rectangle.
//   gnd   ngon  x=0 y=0 r=1.5 n=64 r_in=1.2 rot=0
//     regular n-gon centred on (x, y) with vertex radius r, one vertex on top (rotated
//     by rot degrees counterclockwise). r_in cuts a concentric n-gon hole: a ring, such
//     as a coax shield.
//
// parseGeometryText keeps the expressions, evaluateGeometry turns them into metres.
// Parameter overrides go to evaluateGeometry, which is what parameter sweeps use.

import { isMirrorShape, mirrorShapeX, polyRadiusForArea, shapeBBox } from './shapes.js';
import { hexToRGB, rgbToHex } from './body_colors.js';

export const LENGTH_UNITS = {
    m: 1, cm: 1e-2, mm: 1e-3, um: 1e-6, 'µm': 1e-6, nm: 1e-9, mil: 25.4e-6, in: 25.4e-3,
};

const RECT_KINDS = ['diel', 'gnd', 'sig+', 'sig-'];
const BOUND_VALUES = ['open', 'gnd'];
// Order of the bounds values.
export const WALLS = ['left', 'right', 'top', 'bottom'];
const RESERVED = new Set(['inf', 'auto', 'min', 'max', 'abs', 'sqrt']);
// No prototype: a name like constructor or toString is not a function of the format.
const FUNCTIONS = Object.assign(Object.create(null), {
    min: Math.min, max: Math.max, abs: Math.abs, sqrt: Math.sqrt,
});
// Shape words after the kind. A rectangle has none ('rect' may be written).
const SHAPES = ['rect', 'trap', 'ngon', 'ellipse'];
const RECT_KEY_ORDER = ['x', 'w', 'y', 'h', 'r', 'r_in', 'rx', 'ry', 'rx_in', 'ry_in', 'n', 'rot', 'angle', 'angle2',
    'radius', 'radius_bottom', 'wall', 'er', 'tand', 'thin',
    'sigma', 'rq', 'plating', 'plating_sigma', 'plating_t', 'plating_rq', 'plating_rq_iface', 'color', 'mirror'];
const RECT_KEYS = new Set(RECT_KEY_ORDER);
// Keys only a conductor takes.
export const CONDUCTOR_KEYS = ['plating', 'sigma', 'rq', 'plating_sigma', 'plating_t', 'plating_rq', 'plating_rq_iface'];
// Shapes described by a centre and radii, and the keys of each.
export const ROUND_KEYS = { ngon: ['x', 'y', 'r', 'r_in', 'n', 'rot'], ellipse: ['x', 'y', 'rx', 'ry', 'rx_in', 'ry_in', 'n', 'rot'] };
// Geometry keys each shape takes.
const SHAPE_KEYS = {
    rect: ['x', 'y', 'w', 'h', 'radius', 'radius_bottom', 'wall'],
    trap: ['x', 'y', 'w', 'h', 'angle', 'angle2', 'radius', 'radius_bottom', 'wall'],
    ...ROUND_KEYS,
};
const GEOMETRY_KEYS = new Set(Object.values(SHAPE_KEYS).flat());
export const SHAPE_NAMES = { rect: 'rectangle', trap: 'trapezoid', ngon: 'n-gon', ellipse: 'ellipse' };
// Arc segments per quarter turn of a rounded corner.
const CORNER_SEGMENTS = 8;
const MIRROR_KIND = { 'sig+': 'sig-', 'sig-': 'sig+', gnd: 'gnd', diel: 'diel' };
const PLATING_KEYS = new Set(['sigma', 't', 'rq', 'rq_iface']);
export const PLATING_FACES = ['top', 'sides', 'bottom'];

// --- Expressions ------------------------------------------------------------------

function tokenize(src) {
    const tokens = [];
    const re = /\s*(?:(\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)([A-Za-zµ]+)?|([A-Za-z_][A-Za-z_0-9]*|µm)|([-+*/(),]))/y;
    let pos = 0;
    while (pos < src.length) {
        if (/^\s*$/.test(src.slice(pos))) break;
        re.lastIndex = pos;
        const m = re.exec(src);
        if (!m) throw new Error(`unexpected character '${src.slice(pos).trim()[0]}'`);
        if (m[1] !== undefined) {
            // start / end: the digits in src, for rewriting the number in place.
            const start = pos + m[0].length - m[0].trimStart().length;
            tokens.push({ type: 'num', value: parseFloat(m[1]), unit: m[2], start, end: start + m[1].length });
        }
        else if (m[3] !== undefined) {
            // A unit after a space belongs to the number before it ("1 um"). Units are
            // reserved names, so this cannot be a parameter.
            const prev = tokens[tokens.length - 1];
            if (prev && prev.type === 'num' && prev.unit === undefined && LENGTH_UNITS[m[3]] !== undefined) prev.unit = m[3];
            else tokens.push({ type: 'id', value: m[3] });
        }
        else tokens.push({ type: 'op', value: m[4] });
        pos = re.lastIndex;
    }
    return tokens;
}

// Recursive-descent parser. `build` makes the result bottom-up as the parse goes:
// num(tk), id(name), fn(name, args), neg(node), prod(factors) and sum(terms), the
// last two as [{ op, node }] with op '*' or '+' on the first.
function parseExpression(src, build) {
    const tokens = tokenize(src);
    let i = 0;
    const isOp = v => tokens[i] && tokens[i].type === 'op' && tokens[i].value === v;
    function primary() {
        const tk = tokens[i++];
        if (!tk) throw new Error('unexpected end of expression');
        if (tk.type === 'num') {
            if (tk.unit !== undefined && LENGTH_UNITS[tk.unit] === undefined) throw new Error(`unknown unit '${tk.unit}'`);
            return build.num(tk);
        }
        if (tk.type === 'id') {
            if (!isOp('(')) return build.id(tk.value);
            if (!FUNCTIONS[tk.value]) throw new Error(`unknown function '${tk.value}'`);
            i++;
            const args = [sum()];
            while (isOp(',')) { i++; args.push(sum()); }
            if (!isOp(')')) throw new Error("expected ')'");
            i++;
            return build.fn(tk.value, args);
        }
        if (tk.value === '(') {
            const v = sum();
            if (!isOp(')')) throw new Error("expected ')'");
            i++;
            return v;
        }
        throw new Error(`unexpected '${tk.value}'`);
    }
    function unary() {
        if (isOp('-')) { i++; return build.neg(unary()); }
        if (isOp('+')) { i++; return unary(); }
        return primary();
    }
    function chain(next, ops, first, make) {
        const items = [{ op: first, node: next() }];
        while (ops.some(isOp)) { const op = tokens[i++].value; items.push({ op, node: next() }); }
        return items.length === 1 ? items[0].node : make(items);
    }
    const product = () => chain(unary, ['*', '/'], '*', build.prod);
    const sum = () => chain(product, ['+', '-'], '+', build.sum);
    const out = sum();
    if (i < tokens.length) throw new Error(`unexpected '${tokens[i].value}'`);
    return out;
}

// Evaluates `src`. `vars` maps parameter names to numbers, `unitScale` is the declared
// unit in metres (a suffixed number is converted to it).
export function evaluateExpression(src, vars = {}, unitScale = 1) {
    const fold = items => items.slice(1).reduce((v, { op, node }) => (
        op === '*' ? v * node : op === '/' ? v / node : op === '+' ? v + node : v - node), items[0].node);
    const v = parseExpression(src, {
        num: tk => (tk.unit === undefined ? tk.value : tk.value * LENGTH_UNITS[tk.unit] / unitScale),
        id: name => {
            if (name === 'inf') return Infinity;
            if (!Object.prototype.hasOwnProperty.call(vars, name)) throw new Error(`unknown parameter '${name}'`);
            return vars[name];
        },
        fn: (name, args) => FUNCTIONS[name](...args),
        neg: x => -x,
        prod: fold,
        sum: fold,
    });
    if (Number.isNaN(v)) throw new Error('expression is not a number');
    return v;
}

// A name that cannot be a parameter: functions, inf, auto and the length units.
export const isReservedName = name => RESERVED.has(name) || LENGTH_UNITS[name] !== undefined;

// A plain number, no unit: "0.2", "-1e-3".
export const isPlainNumber = expr => PLAIN_NUMBER.test(String(expr).trim());

// A number with an optional length unit: "0.2", "35um", "35 um".
export function isLengthLiteral(expr) {
    const m = /^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?\s*([A-Za-zµ]*)$/.exec(String(expr).trim());
    return !!m && (!m[1] || LENGTH_UNITS[m[1]] !== undefined);
}

// --- Parsing ----------------------------------------------------------------------

// Index of the '#' that starts the comment of a line, -1 for none. A '#' right after '='
// starts a value (color=#c0c0c0).
export function commentStart(line) {
    for (let i = line.indexOf('#'); i >= 0; i = line.indexOf('#', i + 1)) {
        if (line[i - 1] !== '=') return i;
    }
    return -1;
}

function parseFields(words, allowed, what) {
    const fields = {};
    // A bare unit after a value that ends in a number belongs to it ("rq=1 um").
    const merged = [];
    for (const word of words) {
        if (merged.length && LENGTH_UNITS[word] !== undefined && /=.*[\d.]$/.test(merged[merged.length - 1])) {
            merged[merged.length - 1] += word;
        } else merged.push(word);
    }
    for (const word of merged) {
        const m = /^([A-Za-z_][A-Za-z_0-9]*)=(.+)$/.exec(word);
        if (!m) throw new Error(`expected key=value, got '${word}'`);
        if (!allowed.has(m[1])) throw new Error(`unknown ${what} key '${m[1]}'`);
        if (m[1] in fields) throw new Error(`duplicate key '${m[1]}'`);
        fields[m[1]] = m[2];
    }
    return fields;
}

function parseStatement(src) {
    // A parameter statement has nothing but the name before '='.
    const param = /^([A-Za-z_][A-Za-z_0-9]*)\s*=\s*(.+)$/.exec(src);
    if (param) {
        if (isReservedName(param[1])) {
            throw new Error(`'${param[1]}' is a reserved name`);
        }
        return { type: 'param', name: param[1], expr: param[2].trim() };
    }
    const words = src.split(/\s+/);
    const head = words[0];
    if (head === 'units') {
        if (words.length !== 2 || LENGTH_UNITS[words[1]] === undefined) {
            throw new Error(`units must be one of ${Object.keys(LENGTH_UNITS).join(', ')}`);
        }
        return { type: 'units', value: words[1] };
    }
    if (head === 'bounds') {
        const values = words.slice(1);
        if (values.length !== 4 || !values.every(v => BOUND_VALUES.includes(v))) {
            throw new Error('bounds takes four values (left right top bottom), each open or gnd');
        }
        return { type: 'bounds', values };
    }
    if (head === 'domain') {
        let values = words.slice(1);
        if (values.length === 1 && values[0] === 'auto') values = ['auto', 'auto', 'auto', 'auto'];
        if (values.length !== 4) throw new Error('domain takes auto or four values (x1 x2 y1 y2)');
        return { type: 'domain', values };
    }
    if (head === 'plating') {
        if (words.some(w => /^thick_corners=/.test(w))) {
            throw new Error('thick_corners is the Model Thick Plating option in Advanced Options, for every conductor');
        }
        return { type: 'plating', fields: parseFields(words.slice(1), PLATING_KEYS, 'plating') };
    }
    if (RECT_KINDS.includes(head)) {
        let shape;
        let rest = words.slice(1);
        if (rest.length && SHAPES.includes(rest[0])) {
            if (rest[0] !== 'rect') shape = rest[0];
            rest = rest.slice(1);
        }
        const fields = parseFields(rest, RECT_KEYS, SHAPE_NAMES[shape ?? 'rect']);
        for (const k of Object.keys(fields)) {
            if (GEOMETRY_KEYS.has(k) && !SHAPE_KEYS[shape ?? 'rect'].includes(k)) {
                throw new Error(`${k}= does not apply to ${shape === 'ellipse' ? 'an' : 'a'} ${SHAPE_NAMES[shape ?? 'rect']}`);
            }
        }
        return shape ? { type: 'rect', kind: head, shape, fields } : { type: 'rect', kind: head, fields };
    }
    throw new Error(`unknown statement '${head}'`);
}

// Returns { statements, errors }. Every statement carries its 1-based source line.
// Comment-only and blank lines are kept so serializeGeometry can write them back.
export function parseGeometryText(text) {
    const statements = [], errors = [];
    const lines = String(text ?? '').split(/\r?\n/);
    for (let li = 0; li < lines.length; li++) {
        const line = li + 1;
        const hash = commentStart(lines[li]);
        const code = hash >= 0 ? lines[li].slice(0, hash) : lines[li];
        const comment = hash >= 0 ? lines[li].slice(hash + 1).trim() : null;
        const parts = code.split(';').map(s => s.trim()).filter(s => s.length > 0);
        if (parts.length === 0) {
            statements.push({ type: comment !== null ? 'comment' : 'blank', line, comment });
            continue;
        }
        parts.forEach((part, pi) => {
            try {
                const st = parseStatement(part);
                st.line = line;
                st.part = pi;   // position among the ';' separated statements of the line
                st.comment = pi === parts.length - 1 ? comment : null;
                statements.push(st);
            } catch (e) {
                errors.push({ line, message: e.message });
            }
        });
    }
    for (const type of ['units', 'bounds', 'domain', 'plating']) {
        const dup = statements.filter(s => s.type === type);
        for (const s of dup.slice(1)) errors.push({ line: s.line, message: `more than one ${type} statement` });
    }
    const seen = new Set();
    for (const s of statements) {
        if (s.type !== 'param') continue;
        if (seen.has(s.name)) errors.push({ line: s.line, message: `parameter '${s.name}' defined twice` });
        seen.add(s.name);
    }
    return { statements, errors };
}

function fieldsToText(fields, order) {
    return order.filter(k => fields[k] !== undefined).map(k => `${k}=${fields[k]}`).join(' ');
}

export function serializeGeometry(model) {
    const out = [];
    for (const s of model.statements) {
        let code = '';
        if (s.type === 'units') code = `units ${s.value}`;
        else if (s.type === 'param') code = `${s.name} = ${s.expr}`;
        else if (s.type === 'bounds') code = `bounds ${s.values.join(' ')}`;
        else if (s.type === 'domain') {
            code = s.values.every(v => v === 'auto') ? 'domain auto' : `domain ${s.values.join(' ')}`;
        } else if (s.type === 'plating') {
            code = `plating ${fieldsToText(s.fields, ['sigma', 't', 'rq', 'rq_iface'])}`;
        } else if (s.type === 'rect') {
            code = `${s.kind}${s.shape ? ' ' + s.shape : ''} ${fieldsToText(s.fields, RECT_KEY_ORDER)}`;
        }
        const comment = s.comment !== null && s.comment !== undefined ? `# ${s.comment}` : '';
        out.push([code, comment].filter(p => p.length > 0).join('  '));
    }
    while (out.length && out[out.length - 1] === '') out.pop();
    return out.join('\n') + '\n';
}

// --- In-place edits ---------------------------------------------------------------
// The editor keeps the text as written by the user, so a change made from a sidebar
// control rewrites one statement and leaves the rest of the text alone.

function splitComment(line) {
    const hash = commentStart(line);
    return hash >= 0 ? [line.slice(0, hash), line.slice(hash)] : [line, ''];
}

// Rewrites the text line by line: edit(parts, comment, i) gets the ';' separated
// statement sources of line i (its comment split off) and returns true when it changed
// them in place. With `first` the walk stops at the first change. Returns the new
// text, null when nothing changed.
function editStatements(text, edit, first = false) {
    const lines = String(text ?? '').split(/\r?\n/);
    let changed = false;
    for (let i = 0; i < lines.length; i++) {
        const [code, comment] = splitComment(lines[i]);
        const parts = code.split(';');
        if (!edit(parts, comment, i)) continue;
        lines[i] = parts.join(';') + comment;
        changed = true;
        if (first) break;
    }
    return changed ? lines.join('\n') : null;
}

// Index in `parts` of the statement the parser numbered `part`: it counts only the
// non-empty ones. -1 when there is none.
function partIndex(parts, part) {
    let seen = -1;
    return parts.findIndex(p => p.trim().length > 0 && ++seen === part);
}

// Replaces the expression of parameter `name`. Returns the text unchanged when the
// parameter is not defined.
export function setParamInText(text, name, expr) {
    const re = new RegExp(`^(\\s*${name}\\s*=\\s*)(.*?)(\\s*)$`);
    return editStatements(text, parts => parts.some((p, k) => {
        const m = re.exec(p);
        if (m) parts[k] = m[1] + expr + m[3];
        return !!m;
    }), true) ?? text;
}

// A parameter definition: nothing but a name before the first '='.
const PARAM_DEF = /^(\s*)([A-Za-z_][A-Za-z_0-9]*)(\s*=)(.*)$/s;

// Renames parameter `from` to `to`: its definition and every expression that uses it.
// Only expressions change: statement words (gnd, domain, ...), field keys (the w of
// w=...), the face list of plating=, colors and comments are left alone, so a parameter
// may share its name with any of them. Values hold no whitespace, so a statement splits
// into words on it.
export function renameParamInText(text, from, to) {
    const ident = new RegExp(`(?<![A-Za-z_0-9.])${from}(?![A-Za-z_0-9])`, 'g');
    const expr = e => e.replace(ident, to);
    // Field values: key=value words, except the non-expression keys.
    const fieldValue = w => w.replace(/^([A-Za-z_][A-Za-z_0-9]*=)(.*)$/s,
        (m, key, v) => (key === 'plating=' || key === 'color=' ? m : key + expr(v)));
    return editStatements(text, parts => {
        parts.forEach((part, k) => {
            const def = PARAM_DEF.exec(part);
            if (def) {
                parts[k] = def[1] + (def[2] === from ? to : def[2]) + def[3] + expr(def[4]);
                return;
            }
            const words = part.split(/(\s+)/);
            const head = words.find(w => w.trim());
            if (head === 'units' || head === 'bounds') return;
            let seenHead = false;
            parts[k] = words.map(w => {
                if (!w.trim()) return w;
                if (!seenHead) { seenHead = true; return w; }
                // The domain values are bare expressions, everything else key=value.
                return head === 'domain' ? expr(w) : fieldValue(w);
            }).join('');
        });
        return true;
    });
}

// Statement text for a rectangle or shape, the inverse of the parser for one statement.
export function rectStatementText(kind, fields, shape = null) {
    return `${kind}${shape && shape !== 'rect' ? '  ' + shape : ''}  ${fieldsToText(fields, RECT_KEY_ORDER)}`;
}

// Replaces the statement `st` (from parseGeometryText of the same text) by `code`, or
// removes it when code is null. A line left empty by a removal is dropped.
export function replaceStatementInText(text, st, code) {
    const lines = String(text ?? '').split(/\r?\n/);
    const i = st.line - 1;
    if (i < 0 || i >= lines.length) return text;
    const [src, comment] = splitComment(lines[i]);
    const parts = src.split(';');
    const k = partIndex(parts, st.part);
    if (k < 0) return text;
    if (code === null) parts.splice(k, 1);
    else parts[k] = (k > 0 ? ' ' : '') + code;
    const rest = parts.join(';');
    if (code === null && rest.trim().length === 0) lines.splice(i, 1);
    else lines[i] = rest + (comment && code !== null && !rest.endsWith(' ') ? '  ' : '') + comment;
    return lines.join('\n');
}

// Inserts a line after 1-based line `after` (0 puts it first, a value past the end
// appends).
export function insertLineInText(text, after, code) {
    const lines = String(text ?? '').split(/\r?\n/);
    // Append before the trailing empty line of a text that ends with a newline.
    let at = Math.min(Math.max(after, 0), lines.length);
    if (at === lines.length && lines.length && lines[lines.length - 1] === '') at = lines.length - 1;
    lines.splice(at, 0, code);
    return lines.join('\n');
}

// Swaps the statement `st` with the rectangle statement before (dir = -1) or after it
// (dir = 1). Both must be alone on their lines. Order matters between dielectrics.
export function moveRectInText(text, model, st, dir) {
    const rects = model.statements.filter(s => s.type === 'rect');
    const other = rects[rects.indexOf(st) + dir];
    const alone = s => model.statements.filter(o => o.line === s.line && o.type !== 'comment' && o.type !== 'blank').length === 1;
    if (!other || !alone(st) || !alone(other)) return text;
    const lines = String(text ?? '').split(/\r?\n/);
    [lines[st.line - 1], lines[other.line - 1]] = [lines[other.line - 1], lines[st.line - 1]];
    return lines.join('\n');
}

// Replaces the first statement starting with `keyword` (bounds, domain, units, plating)
// by `statement`, or inserts it after the units line (at the top without one).
export function setStatementInText(text, keyword, statement) {
    // A parameter of the same name (domain = 3) is not the statement.
    const startsWith = (s, word) => !PARAM_DEF.test(s) && (s.trim() === word || s.trim().startsWith(word + ' '));
    let unitsLine = -1;
    const out = editStatements(text, (parts, comment, i) => {
        const k = parts.findIndex(s => startsWith(s, keyword));
        if (k < 0) {
            if (unitsLine < 0 && parts.some(s => startsWith(s, 'units'))) unitsLine = i;
            return false;
        }
        parts[k] = /^\s*/.exec(parts[k])[0] + statement + (k < parts.length - 1 || !comment ? '' : '  ');
        return true;
    }, true);
    if (out !== null) return out;
    const lines = String(text ?? '').split(/\r?\n/);
    lines.splice(unitsLine + 1, 0, statement);
    return lines.join('\n');
}

// --- Expression arithmetic for the form's edits --------------------------------------
// Builds new expressions from the ones written, as short as they reasonably get, so a
// rectangle placed from another one keeps the parameters instead of freezing numbers.

const PLAIN_NUMBER = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;
const isNeg = e => e.trim() === '-inf';
const isPosInf = e => /^\+?inf$/.test(e.trim());
const fmtNumber = v => String(parseFloat(v.toPrecision(12)));

// Positions of the + and - signs at the top level of `e` that join terms (not a
// leading sign, not an exponent sign).
function topLevelSigns(e) {
    const out = [];
    let depth = 0;
    for (let i = 0; i < e.length; i++) {
        const c = e[i];
        if (c === '(') depth++;
        else if (c === ')') depth--;
        else if ((c === '+' || c === '-') && depth === 0 && i > 0) {
            if (/[eE]/.test(e[i - 1]) && /[\d.]/.test(e[i - 2] ?? '')) continue;
            if (/[-+*/(]/.test(e[i - 1])) continue;   // unary sign after an operator
            out.push(i);
        }
    }
    return out;
}
const hasTopLevelSum = e => topLevelSigns(e).length > 0;

export function addExpr(a, b) {
    a = a.trim(); b = b.trim();
    if (PLAIN_NUMBER.test(a) && PLAIN_NUMBER.test(b)) return fmtNumber(parseFloat(a) + parseFloat(b));
    if (PLAIN_NUMBER.test(a) && parseFloat(a) === 0) return b;
    if (PLAIN_NUMBER.test(b) && parseFloat(b) === 0) return a;
    if (b.startsWith('-') && !hasTopLevelSum(b)) return `${a}-${b.slice(1)}`;
    return `${a}+${b}`;
}

const AXIS_KEYS = { x: ['x', 'w'], y: ['y', 'h'] };

// Low and high edge of an axis as expressions. `negative` says the size evaluates
// negative (the `flipped` flag of the evaluated axis).
export function axisEdges(fields, axis, negative = false) {
    const [p, s] = AXIS_KEYS[axis];
    const pos = fields[p] ?? '0', size = fields[s] ?? '0';
    if (isNeg(pos)) return { lo: '-inf', hi: 'inf' };
    if (isPosInf(size)) return { lo: pos, hi: 'inf' };
    if (isNeg(size)) return { lo: '-inf', hi: pos };
    return negative ? { lo: addExpr(pos, size), hi: pos } : { lo: pos, hi: addExpr(pos, size) };
}

// --- Evaluation -------------------------------------------------------------------

// One axis of a rectangle in metres: { pos, size, min, max, flipped }. pos and size are
// the constructor arguments of Conductor / Dielectric, kept as written so a geometry
// converted from a solver reproduces its doubles exactly. Those take a negative height
// but not a negative width, so a negative w is turned around to x+w, -w. min / max are
// the bounds, -Infinity / Infinity for a pinned edge (pos and size are then resolved
// against the domain by the solver). flipped says the size was written negative.
function evalAxis(fields, pos, size, ev, keepNegative) {
    if (fields[pos] === undefined || fields[size] === undefined) {
        throw new Error(`${pos} and ${size} must both be given`);
    }
    const p = ev(fields[pos]), s = ev(fields[size]);
    if (p === Infinity) throw new Error(`${pos} cannot be inf`);
    if (p === -Infinity && s !== Infinity) throw new Error(`${pos}=-inf needs ${size}=inf`);
    if (s === 0) throw new Error(`${size} must be nonzero`);
    if (s === -Infinity) return { pos: -Infinity, size: Infinity, min: -Infinity, max: p, flipped: true };
    if (s < 0) return { pos: keepNegative ? p : p + s, size: keepNegative ? s : -s, min: p + s, max: p, flipped: true };
    return { pos: p, size: s, min: p, max: s === Infinity ? Infinity : p + s, flipped: false };
}

// Trapezoid on the base [x0, x1] at y = yb, height h (negative: the other face is
// below), side angles aL, aR in degrees from the vertical, as a convex CCW polygon with
// the face name of each edge. Vertices 0 and 1 are the lower corners.
function trapezoidPolygon(x0, x1, yb, h, aL, aR) {
    for (const a of [aL, aR]) {
        if (!(Math.abs(a) < 89)) throw new Error('angle must be between -89 and 89 degrees');
    }
    const H = Math.abs(h), yt = yb + h;
    const iL = H * Math.tan(aL * Math.PI / 180), iR = H * Math.tan(aR * Math.PI / 180);
    const tl = x0 + iL, tr = x1 - iR;
    if (!(tr - tl > (x1 - x0) * 1e-6)) throw new Error('the side angles leave the trapezoid no face opposite its base');
    const poly = h > 0
        ? [x0, yb, x1, yb, tr, yt, tl, yt]
        : [tl, yt, tr, yt, x1, yb, x0, yb];
    return { poly: new Float64Array(poly), faces: ['bottom', 'sides', 'top', 'sides'] };
}

// The convex CCW polygon moved inward by t along every edge normal: the inside of a
// wall of thickness t. Throws when the wall closes the inside.
function insetPolygon(poly, t) {
    const n = poly.length >> 1;
    const lines = [];
    for (let i = 0; i < n; i++) {
        const j = (i + 1) % n;
        const ex = poly[2 * j] - poly[2 * i], ey = poly[2 * j + 1] - poly[2 * i + 1];
        const l = Math.hypot(ex, ey);
        // Inward normal of a CCW edge.
        lines.push({ px: poly[2 * i] - ey / l * t, py: poly[2 * i + 1] + ex / l * t, dx: ex / l, dy: ey / l });
    }
    const out = new Float64Array(2 * n);
    for (let i = 0; i < n; i++) {
        const a = lines[(i + n - 1) % n], b = lines[i];
        const den = a.dx * b.dy - a.dy * b.dx;
        const u = ((b.px - a.px) * b.dy - (b.py - a.py) * b.dx) / den;
        out[2 * i] = a.px + u * a.dx; out[2 * i + 1] = a.py + u * a.dy;
    }
    for (let i = 0; i < n; i++) {
        const j = (i + 1) % n;
        const dot = (out[2 * j] - out[2 * i]) * lines[i].dx + (out[2 * j + 1] - out[2 * i + 1]) * lines[i].dy;
        if (!(dot > 0)) throw new Error('wall is too thick: it leaves no inside');
    }
    return out;
}

// Rounds the corners of a convex CCW polygon: corner i becomes an arc of radius radii[i]
// tangent to both of its edges (0 keeps the corner). Each half of an arc takes the face
// name of the edge it runs into. Returns { poly, faces }.
function roundPolygon(poly, faces, radii) {
    const n = poly.length >> 1;
    const P = i => [poly[2 * ((i + n) % n)], poly[2 * ((i + n) % n) + 1]];
    const unit = (x, y) => { const l = Math.hypot(x, y); return [x / l, y / l]; };
    // Unit vectors along the two edges of each rounded corner and the angle between them.
    const corner = radii.map((R, i) => {
        if (!(R > 0)) return null;
        const [px, py] = P(i), [ax, ay] = P(i - 1), [bx, by] = P(i + 1);
        const u1 = unit(ax - px, ay - py), u2 = unit(bx - px, by - py);
        return { u1, u2, th: Math.acos(Math.max(-1, Math.min(1, u1[0] * u2[0] + u1[1] * u2[1]))) };
    });
    // Distance from each corner to its tangent points.
    const cut = corner.map((c, i) => (c ? radii[i] / Math.tan(c.th / 2) : 0));
    const scale = Math.max(...Array.from(poly).map(Math.abs));
    for (let i = 0; i < n; i++) {
        const [ax, ay] = P(i), [bx, by] = P(i + 1);
        const side = Math.hypot(bx - ax, by - ay), j = (i + 1) % n;
        if (cut[i] + cut[j] > side * (1 + 1e-9)) throw new Error('corner radius is larger than the sides allow');
        // Two arcs using up a side within rounding meet exactly (a stadium): an
        // overshoot would fold the outline back on itself.
        if (cut[i] + cut[j] > side) { const k = side / (cut[i] + cut[j]); cut[i] *= k; cut[j] *= k; }
    }
    const pts = [];   // { x, y, face of the edge that starts here }
    for (let i = 0; i < n; i++) {
        const [px, py] = P(i);
        const fin = faces[(i + n - 1) % n], fout = faces[i];
        if (!(radii[i] > 0)) { pts.push({ x: px, y: py, face: fout }); continue; }
        const R = radii[i];
        const { u1, u2, th } = corner[i];
        const bis = unit(u1[0] + u2[0], u1[1] + u2[1]);
        const cx = px + bis[0] * R / Math.sin(th / 2), cy = py + bis[1] * R / Math.sin(th / 2);
        const t1x = px + u1[0] * cut[i], t1y = py + u1[1] * cut[i];
        const span = Math.PI - th;
        // An even count splits the arc's faces at its middle.
        const m = 2 * Math.max(1, Math.ceil(CORNER_SEGMENTS * span / (Math.PI / 2) / 2));
        const phi = Math.atan2(t1y - cy, t1x - cx);
        for (let k = 0; k <= m; k++) {
            const a = phi + span * k / m;
            // The end points sit exactly on the edges.
            const x = k === 0 ? t1x : k === m ? px + u2[0] * cut[i] : cx + R * Math.cos(a);
            const y = k === 0 ? t1y : k === m ? py + u2[1] * cut[i] : cy + R * Math.sin(a);
            pts.push({ x, y, face: k < m ? (2 * k < m ? fin : fout) : fout });
        }
    }
    // A side the arcs use up entirely leaves two coincident points.
    const tol = scale * 1e-12;
    const kept = pts.filter((p, i) => {
        const q = pts[(i + 1) % pts.length];
        return Math.hypot(q.x - p.x, q.y - p.y) > tol;
    });
    return { poly: new Float64Array(kept.flatMap(p => [p.x, p.y])), faces: kept.map(p => p.face) };
}

// n vertices on the ellipse of semi-axes rx, ry centred on (cx, cy), at equal steps of
// the ellipse parameter from the top, turned by rot degrees counterclockwise, CCW. The
// steps put the vertices closer together at the tightly curved ends, which keeps the
// polygon on the ellipse there (equal steps along the outline cut the ends off a long
// ellipse). A circle gives the regular n-gon. Without rotation the vertices mirror
// exactly about x = cx.
export function ellipsePolygon(cx, cy, rx, ry, n, rot = 0) {
    // Parameter t of each vertex (x = rx cos t, y = ry sin t, t = pi/2 on top).
    const ts = new Float64Array(n);
    for (let k = 0; k < n; k++) ts[k] = Math.PI / 2 + 2 * Math.PI * k / n;
    const c = Math.cos(rot * Math.PI / 180), s = Math.sin(rot * Math.PI / 180);
    const poly = new Float64Array(2 * n);
    for (let k = 0; k < n; k++) {
        if (rot === 0 && k > n / 2) {
            poly[2 * k] = cx - (poly[2 * (n - k)] - cx);
            poly[2 * k + 1] = poly[2 * (n - k) + 1];
            continue;
        }
        const onAxis = rot === 0 && (k === 0 || 2 * k === n);
        const ex = onAxis ? 0 : rx * Math.cos(ts[k]);
        const ey = onAxis ? (k === 0 ? ry : -ry) : ry * Math.sin(ts[k]);
        poly[2 * k] = rot === 0 ? cx + ex : cx + c * ex - s * ey;
        poly[2 * k + 1] = rot === 0 ? cy + ey : cy + s * ex + c * ey;
    }
    return poly;
}

function bboxAxes(poly) {
    const { xmin, xmax, ymin, ymax } = shapeBBox({ type: 'polygon', poly });
    return { x: { pos: xmin, size: xmax - xmin, min: xmin, max: xmax, flipped: false },
             y: { pos: ymin, size: ymax - ymin, min: ymin, max: ymax, flipped: false } };
}

// Shape of a statement in metres: { x, y axes, shape }, shape null for a plain
// rectangle (no angle, radius or wall). The shape is a convex CCW polygon, or a ring
// { poly, hole }, with the fields shapes.js describes.
// xAxis replaces the x axis of a rectangle (mirror=1 widening it to be symmetric).
function evalShape(st, f, len, num, xAxis = null) {
    const kind = st.shape ?? 'rect';
    if (kind === 'rect' || kind === 'trap') {
        const x = xAxis ?? evalAxis(f, 'x', 'w', len, false);
        const y = evalAxis(f, 'y', 'h', len, true);
        const aL = f.angle !== undefined ? num(f.angle) : (f.angle2 !== undefined ? num(f.angle2) : 0);
        const aR = f.angle2 !== undefined ? num(f.angle2) : aL;
        const rTop = f.radius !== undefined ? len(f.radius) : 0;
        const rBot = f.radius_bottom !== undefined ? len(f.radius_bottom) : rTop;
        const wall = f.wall !== undefined ? len(f.wall) : 0;
        if (!(rTop >= 0) || !(rBot >= 0) || !Number.isFinite(rTop) || !Number.isFinite(rBot)) throw new Error('radius must be non-negative');
        if (!(wall >= 0) || !Number.isFinite(wall)) throw new Error('wall must be non-negative');
        if (aL === 0 && aR === 0 && rTop === 0 && rBot === 0 && wall === 0) return { x, y, shape: null };
        if (![x.min, x.max, y.min, y.max].every(Number.isFinite)) {
            throw new Error(`a ${SHAPE_NAMES[kind]} with ${kind === 'trap' ? 'angles' : 'a radius or wall'} needs finite x, w, y and h`);
        }
        // The base is the face at y, the other face is y + h.
        const yb = y.flipped ? y.max : y.min;
        const tp = trapezoidPolygon(x.min, x.max, yb, y.flipped ? -(y.max - y.min) : y.max - y.min, aL, aR);
        // Vertices 0, 1 are the lower corners, 2, 3 the upper ones.
        const radii = [rBot, rBot, rTop, rTop];
        const outer = roundPolygon(tp.poly, tp.faces, radii);
        const ax = bboxAxes(outer.poly);
        let shape = { type: 'polygon', prim: kind, poly: outer.poly, faces: outer.faces,
                      thickness: Math.min(y.max - y.min, ax.x.size) };
        if (wall > 0) {
            // The inside keeps the corner centres: its radii are the outer ones less the wall.
            const hole = roundPolygon(insetPolygon(tp.poly, wall), tp.faces, radii.map(r => Math.max(r - wall, 0)));
            shape = { type: 'ring', prim: kind, poly: outer.poly, hole: hole.poly, faces: outer.faces, thickness: wall };
        }
        return { x: { ...ax.x, flipped: x.flipped }, y, shape };
    }
    const round = kind === 'ngon' ? ['r', 'r'] : ['rx', 'ry'];
    for (const k of ['x', 'y', ...round, 'n']) {
        if (f[k] === undefined) throw new Error(`${kind === 'ngon' ? 'an n-gon' : 'an ellipse'} needs x, y, ${[...new Set(round)].join(', ')} and n`);
    }
    const cx = len(f.x), cy = len(f.y), rx = len(f[round[0]]), ry = len(f[round[1]]);
    const n = num(f.n), rot = f.rot !== undefined ? num(f.rot) : 0;
    if (![cx, cy].every(Number.isFinite)) throw new Error(`${SHAPE_NAMES[kind]} x and y must be finite`);
    if (!(rx > 0) || !(ry > 0) || !Number.isFinite(rx) || !Number.isFinite(ry)) throw new Error(`${[...new Set(round)].join(' and ')} must be positive`);
    if (!Number.isInteger(n) || n < 3 || n > 1024) throw new Error('n must be a whole number from 3 to 1024');
    if (!Number.isFinite(rot)) throw new Error('rot must be a number');
    const poly = ellipsePolygon(cx, cy, rx, ry, n, rot);
    const cosn = Math.cos(Math.PI / n);
    let shape = { type: 'polygon', prim: kind, poly, round: true, thickness: 2 * Math.min(rx, ry) * cosn };
    // A regular n-gon lies between its inscribed and circumscribed circles, which
    // decides containment without the edge loop away from the boundary.
    if (rx === ry) shape.radial = { cx, cy, rIn: rx * cosn, rOut: rx };
    const inner = kind === 'ngon' ? ['r_in', 'r_in'] : ['rx_in', 'ry_in'];
    if (f[inner[0]] !== undefined || f[inner[1]] !== undefined) {
        if (f[inner[0]] === undefined || f[inner[1]] === undefined) throw new Error('an elliptical ring needs rx_in and ry_in');
        const ix = len(f[inner[0]]), iy = len(f[inner[1]]);
        if (!(ix > 0) || !(ix < rx) || !(iy > 0) || !(iy < ry)) {
            throw new Error(kind === 'ngon' ? 'r_in must be positive and smaller than r' : 'rx_in and ry_in must be positive and smaller than rx and ry');
        }
        shape = { type: 'ring', prim: kind, poly, hole: ellipsePolygon(cx, cy, ix, iy, n, rot),
                  thickness: Math.min(rx - ix, ry - iy) * cosn };
        if (rx === ry && ix === iy) shape.radial = { cx, cy, rIn: rx * cosn, rOut: rx, holeIn: ix * cosn, holeOut: ix };
    }
    return { ...bboxAxes(poly), shape };
}

// The x axis symmetric about x=0 that covers `ax`, which touches or crosses x=0.
function symmetricAxis(ax) {
    const m = Math.max(-ax.min, ax.max);
    return { pos: -m, size: 2 * m, min: -m, max: m };
}

// Mirror image of an x axis about x=0.
function mirrorAxis(ax) {
    if (ax.max === Infinity) return { pos: -Infinity, size: Infinity, min: -Infinity, max: -ax.min };
    return { pos: -ax.max, size: ax.min === -Infinity ? Infinity : ax.size, min: -ax.max, max: -ax.min };
}

// Evaluates a parsed model to numbers in metres.
//   overrides - { name: value } replaces parameter values (in the declared units)
// Returns { errors, units, params, bounds, domain, plating, rects }:
//   domain  - { x_min, x_max, y_min, y_max }, null for auto
//   plating - { sigma, thickness, rq, rq_interface? } or null
//   rects   - { kind, x: axis, y: axis, er, tand, thin, plating: faces|null, sigma: S/m|null, rq: m|null,
//               platingMaterial: { sigma?, thickness?, rq?, rq_interface? }|null, color: '#rrggbb'|null,
//               line, image }
//             image is true on a rectangle generated by mirror=1
export function evaluateGeometry(model, overrides = {}) {
    const errors = [...model.errors];
    const unitsSt = model.statements.find(s => s.type === 'units');
    const units = unitsSt ? unitsSt.value : 'mm';
    const scale = LENGTH_UNITS[units];
    const params = {};
    const defined = new Set(model.statements.filter(s => s.type === 'param').map(s => s.name));
    // A name that is defined but has no value yet is not unknown: say why it has none.
    const explain = (s, message) => message.replace(/unknown parameter '([A-Za-z_][A-Za-z_0-9]*)'/, (m, name) => {
        if (!defined.has(name)) return m;
        if (s.type === 'param' && s.name === name) return `parameter '${name}' refers to itself`;
        if (errors.some(e => e.param === name)) return `parameter '${name}' has an error`;
        return `parameter '${name}' is defined below this line`;
    });
    const fail = (s, e) => errors.push({ line: s.line, message: explain(s, e.message ?? String(e)) });
    const num = src => evaluateExpression(src, params, scale);
    const len = src => num(src) * scale;

    for (const name of Object.keys(overrides)) {
        if (!model.statements.some(s => s.type === 'param' && s.name === name)) {
            errors.push({ line: 0, message: `override for unknown parameter '${name}'` });
        }
    }
    for (const s of model.statements) {
        if (s.type !== 'param') continue;
        try {
            const v = Object.prototype.hasOwnProperty.call(overrides, s.name)
                ? Number(overrides[s.name]) : num(s.expr);
            if (!Number.isFinite(v)) throw new Error(`parameter '${s.name}' is not finite`);
            params[s.name] = v;
        } catch (e) { fail(s, e); errors[errors.length - 1].param = s.name; }
    }

    const boundsSt = model.statements.find(s => s.type === 'bounds');
    const bounds = boundsSt ? [...boundsSt.values] : ['open', 'open', 'open', 'open'];

    const domain = { x_min: null, x_max: null, y_min: null, y_max: null };
    const domainSt = model.statements.find(s => s.type === 'domain');
    if (domainSt) {
        try {
            ['x_min', 'x_max', 'y_min', 'y_max'].forEach((k, i) => {
                if (domainSt.values[i] === 'auto') return;
                const v = len(domainSt.values[i]);
                if (!Number.isFinite(v)) throw new Error('domain values must be finite or auto');
                domain[k] = v;
            });
            if (domain.x_min !== null && domain.x_max !== null && !(domain.x_max > domain.x_min)) {
                throw new Error('domain x2 must be greater than x1');
            }
            if (domain.y_min !== null && domain.y_max !== null && !(domain.y_max > domain.y_min)) {
                throw new Error('domain y2 must be greater than y1');
            }
        } catch (e) { fail(domainSt, e); }
    }

    let plating = null;
    const platingSt = model.statements.find(s => s.type === 'plating');
    if (platingSt) {
        try {
            const f = platingSt.fields;
            if (f.sigma === undefined || f.t === undefined) throw new Error('plating needs sigma and t');
            plating = {
                sigma: num(f.sigma), thickness: len(f.t),
                rq: f.rq !== undefined ? len(f.rq) : 0,
            };
            if (f.rq_iface !== undefined) plating.rq_interface = len(f.rq_iface);
            if (!(plating.sigma > 0) || !(plating.thickness > 0) || !(plating.rq >= 0)
                || (plating.rq_interface !== undefined && !(plating.rq_interface >= 0))) {
                throw new Error('plating sigma and t must be positive, rq and rq_iface non-negative');
            }
        } catch (e) { fail(platingSt, e); }
    }

    const rects = [];
    for (const s of model.statements) {
        if (s.type !== 'rect') continue;
        try {
            const f = s.fields;
            const r = { kind: s.kind, line: s.line, part: s.part, shape: null,
                er: 1, tand: 0, thin: false, plating: null, sigma: null, rq: null, platingMaterial: null, color: null };
            if (s.shape || f.radius !== undefined || f.radius_bottom !== undefined || f.wall !== undefined) {
                Object.assign(r, evalShape(s, f, len, num));
            } else { r.x = evalAxis(f, 'x', 'w', len, false); r.y = evalAxis(f, 'y', 'h', len, true); }
            if (s.kind === 'diel') {
                if (f.er === undefined) throw new Error('diel needs er');
                if (CONDUCTOR_KEYS.some(k => k !== 'sigma' && f[k] !== undefined)) {
                    throw new Error('rq and plating apply to conductors only');
                }
                r.er = num(f.er);
                r.tand = f.tand !== undefined ? num(f.tand) : 0;
                r.thin = f.thin !== undefined ? num(f.thin) !== 0 : false;
                if (!(r.er >= 1) || !Number.isFinite(r.er)) throw new Error('er must be a finite number >= 1');
                if (!(r.tand >= 0) || !Number.isFinite(r.tand)) throw new Error('tand must be non-negative');
                if (f.sigma !== undefined) {
                    r.sigma = num(f.sigma);
                    if (!(r.sigma >= 0) || !Number.isFinite(r.sigma)) throw new Error('sigma must be non-negative');
                }
            } else {
                if (f.er !== undefined || f.tand !== undefined || f.thin !== undefined) {
                    throw new Error('er, tand and thin apply to diel only');
                }
                if (f.plating !== undefined && f.plating !== 'none') {
                    const faces = f.plating === 'all' ? PLATING_FACES : f.plating.split(',');
                    for (const face of faces) {
                        if (!PLATING_FACES.includes(face)) {
                            throw new Error(`plating faces are ${PLATING_FACES.join(', ')}, all or none`);
                        }
                    }
                    r.plating = { top: faces.includes('top'), sides: faces.includes('sides'), bottom: faces.includes('bottom') };
                    // A round shape has one surface and no faces to choose from, a ring
                    // has an inside the faces do not name.
                    if (r.shape && (r.shape.round || r.shape.type === 'ring')) {
                        if (!(r.plating.top && r.plating.sides && r.plating.bottom)) {
                            throw new Error(`plating on ${r.shape.type === 'ring' ? 'a ring' : `an ${SHAPE_NAMES[r.shape.prim]}`} covers its whole surface: plating=all`);
                        }
                        r.plating.all = true;
                    }
                }
                if (f.sigma !== undefined) {
                    r.sigma = num(f.sigma);
                    if (!(r.sigma > 0) || !Number.isFinite(r.sigma)) throw new Error('sigma must be positive');
                }
                if (f.rq !== undefined) {
                    r.rq = len(f.rq);
                    if (!(r.rq >= 0) || !Number.isFinite(r.rq)) throw new Error('rq must be non-negative');
                }
                const pm = {};
                if (f.plating_sigma !== undefined) pm.sigma = num(f.plating_sigma);
                if (f.plating_t !== undefined) pm.thickness = len(f.plating_t);
                if (f.plating_rq !== undefined) pm.rq = len(f.plating_rq);
                if (f.plating_rq_iface !== undefined) pm.rq_interface = len(f.plating_rq_iface);
                if (Object.keys(pm).length) {
                    if ((pm.sigma !== undefined && !(pm.sigma > 0)) || (pm.thickness !== undefined && !(pm.thickness > 0))
                        || (pm.rq !== undefined && !(pm.rq >= 0))
                        || (pm.rq_interface !== undefined && !(pm.rq_interface >= 0))) {
                        throw new Error('plating_sigma and plating_t must be positive, plating_rq and plating_rq_iface non-negative');
                    }
                    r.platingMaterial = pm;
                }
            }
            if (f.color !== undefined) {
                const rgb = hexToRGB(f.color);
                if (!rgb) throw new Error('color is #rrggbb or #rgb, such as color=#c0c0c0');
                r.color = rgbToHex(rgb);
            }
            r.image = false;
            const mirror = f.mirror !== undefined && num(f.mirror) !== 0;
            // A rounded or hollow rectangle touching or crossing x=0 widens to one
            // symmetric about it, like a plain rectangle.
            if (mirror && r.shape && !s.shape) {
                const ax = evalAxis(f, 'x', 'w', len, false);
                if (ax.min <= 0 && ax.max >= 0) {
                    Object.assign(r, evalShape(s, f, len, num, { ...symmetricAxis(ax), flipped: false }));
                    rects.push(r);
                    continue;
                }
            }
            if (mirror && r.shape) {
                // A shape crossing x=0 has to be its own image; one touching it from
                // one side (twinax insulations meeting in the middle) gets an image.
                const tol = (r.x.max - r.x.min) * 1e-9;
                if (r.x.min < -tol && r.x.max > tol) {
                    if (!isMirrorShape(r.shape, r.shape, tol)) {
                        throw new Error(`mirror=1 needs the ${SHAPE_NAMES[s.shape ?? 'rect']} on one side of x=0 or symmetric about it`);
                    }
                    rects.push(r);
                } else {
                    rects.push(r, { ...r, kind: MIRROR_KIND[r.kind], x: mirrorAxis(r.x), shape: mirrorShapeX(r.shape), image: true });
                }
                continue;
            }
            if (mirror) {
                if (r.x.min <= 0 && r.x.max >= 0) {
                    r.x = symmetricAxis(r.x);
                    rects.push(r);
                } else {
                    rects.push(r, { ...r, kind: MIRROR_KIND[r.kind], x: mirrorAxis(r.x), image: true });
                }
                continue;
            }
            rects.push(r);
        } catch (e) { fail(s, e); }
    }

    errors.sort((a, b) => a.line - b.line);
    return { errors, units, params, bounds, domain, plating, rects };
}

// A length in metres with a unit that suits its size: 35 µm, 1.2 mm, 50 nm.
export function formatLength(v) {
    if (!Number.isFinite(v)) return v > 0 ? 'inf' : '-inf';
    const a = Math.abs(v);
    const [unit, k] = a >= 1 ? ['m', 1] : a >= 1e-3 ? ['mm', 1e-3] : a >= 1e-6 || a === 0 ? ['µm', 1e-6] : ['nm', 1e-9];
    return `${parseFloat((v / k).toPrecision(4))} ${unit}`;
}

// Material values far outside their usual range, most likely a bare number read in the
// declared unit (rq=1 under units mm is 1 mm). Returns { line, message } notes.
export function plausibilityWarnings(geo, model = null) {
    const out = [];
    const bare = ` A bare number is in the declared unit (${geo.units}): write 1um for micrometres.`;
    const check = (line, what, v, max, thickness) => {
        if (v === null || v === undefined || !Number.isFinite(v)) return;
        if (v > max) out.push({ line, message: `${what} = ${formatLength(v)} is unusually large.${bare}` });
        else if (thickness > 0 && v >= thickness) {
            out.push({ line, message: `${what} = ${formatLength(v)} is not smaller than the conductor thickness ${formatLength(thickness)}.${bare}` });
        }
    };
    const sigma = (line, what, v) => {
        if (v === null || v === undefined || (v >= 1e5 && v <= 1e9)) return;
        out.push({ line, message: `${what} = ${v.toPrecision(3)} S/m is outside the usual range of metals (1e5 to 1e9 S/m).` });
    };
    const RQ_MAX = 20e-6, PLATING_T_MAX = 100e-6;
    const st = model && model.statements.find(o => o.type === 'plating');
    if (geo.plating) {
        const line = st ? st.line : 0;
        sigma(line, 'plating sigma', geo.plating.sigma);
        check(line, 'plating t', geo.plating.thickness, PLATING_T_MAX, 0);
        check(line, 'plating rq', geo.plating.rq, RQ_MAX, 0);
        check(line, 'plating rq_iface', geo.plating.rq_interface, RQ_MAX, 0);
    }
    for (const r of geo.rects) {
        if (r.image || r.kind === 'diel') continue;
        const sizes = [r.x.size, r.y.size].map(Math.abs).filter(Number.isFinite);
        const thickness = sizes.length ? Math.min(...sizes) : 0;
        sigma(r.line, 'sigma', r.sigma);
        check(r.line, 'rq', r.rq, RQ_MAX, thickness);
        const pm = r.platingMaterial || {};
        sigma(r.line, 'plating_sigma', pm.sigma);
        check(r.line, 'plating_t', pm.thickness, PLATING_T_MAX, 0);
        check(r.line, 'plating_rq', pm.rq, RQ_MAX, 0);
        check(r.line, 'plating_rq_iface', pm.rq_interface, RQ_MAX, 0);
    }
    return out;
}

export function parseAndEvaluate(text, overrides = {}) {
    return evaluateGeometry(parseGeometryText(text), overrides);
}

export function formatErrors(errors) {
    return errors.map(e => (e.line > 0 ? `line ${e.line}: ` : '') + e.message).join('\n');
}

// --- Conversion from a solver ------------------------------------------------------

// Writes the rectangles, domain and boundaries of a rectangular solver as geometry text.
//   units    - 'm' (default) prints the doubles exactly, so the text rebuilds the same
//              lists bit for bit. Other units are rounded to 12 significant digits.
//   pinWalls - write edges that lie on a domain wall as -inf / inf so they follow a
//              later change of the domain.
export function solverToGeometryText(solver, { units = 'm', pinWalls = false } = {}) {
    const scale = LENGTH_UNITS[units];
    if (scale === undefined) throw new Error(`unknown unit '${units}'`);
    const all = [...(solver.dielectrics || []), ...(solver.conductors || [])];
    if (solver.is_coax) return coaxToGeometryText(solver, units);
    if (all.some(o => o.shape)) throw new Error('Only rectangular geometries can be converted.');
    const fmt = v => (units === 'm' ? String(v) : String(Number((v / scale).toPrecision(12))));
    // Roughness and plating thickness in micrometres, their usual unit, whatever the
    // declared one. Metres stay exact.
    const fmtSmall = v => (units === 'm' ? String(v) : `${Number((v / 1e-6).toPrecision(12))}um`);
    const X0 = -solver.domain_width / 2, X1 = solver.domain_width / 2;
    const Y0 = solver.domain_y_min, Y1 = solver.domain_height;
    const tol = Math.max(X1 - X0, Y1 - Y0) * 1e-9;

    const axisText = (pos, size, p, s, d0, d1) => {
        const a = s >= 0 ? p : p + s, b = s >= 0 ? p + s : p;
        const loWall = pinWalls && s > 0 && Math.abs(a - d0) <= tol;
        const hiWall = pinWalls && s > 0 && Math.abs(b - d1) <= tol;
        if (!loWall && !hiWall) return `${pos}=${fmt(p)} ${size}=${fmt(s)}`;
        if (loWall && hiWall) return `${pos}=-inf ${size}=inf`;
        if (hiWall) return `${pos}=${fmt(p)} ${size}=inf`;
        return `${pos}=${fmt(b)} ${size}=-inf`;
    };
    const rectText = o => `${axisText('x', 'w', o.x, o.width, X0, X1)} ` +
        `${axisText('y', 'h', o.y, o.height, Y0, Y1)}`;

    const lines = [`units ${units}`, `bounds ${(solver.boundaries || ['open', 'open', 'open', 'gnd']).join(' ')}`,
        `domain ${fmt(X0)} ${fmt(X1)} ${fmt(Y0)} ${fmt(Y1)}`];
    // The first plating becomes the plating statement, a conductor plated with another
    // material carries its own keys.
    const pl = (solver.conductors || []).map(c => c.plating).find(p => p && p.sigma > 0 && p.thickness > 0);
    if (pl) {
        lines.push(`plating sigma=${pl.sigma} t=${fmtSmall(pl.thickness)} rq=${fmtSmall(pl.rq ?? 0)}`
            + (pl.rq_interface !== undefined ? ` rq_iface=${fmtSmall(pl.rq_interface)}` : ''));
    }
    lines.push('');
    for (const d of (solver.dielectrics || [])) {
        lines.push(`diel ${rectText(d)} er=${d.epsilon_r} tand=${d.tan_delta ?? 0}` + (d.sigma > 0 ? ` sigma=${d.sigma}` : '')
            + (d.thin_sheet ? ' thin=1' : ''));
    }
    for (const c of (solver.conductors || [])) {
        const kind = c.is_signal ? (c.polarity < 0 ? 'sig-' : 'sig+') : 'gnd';
        // A plating of no thickness or conductivity has no effect and writes nothing.
        const real = c.plating && c.plating.sigma > 0 && c.plating.thickness > 0;
        const faces = real ? PLATING_FACES.filter(f => c.plating[f]) : [];
        let extra = (c.sigma !== undefined && c.sigma !== null) ? ` sigma=${c.sigma}` : '';
        if (c.rq !== undefined && c.rq !== null) extra += ` rq=${fmtSmall(c.rq)}`;
        if (faces.length) {
            extra += ` plating=${faces.join(',')}`;
            if (c.plating.sigma !== pl.sigma) extra += ` plating_sigma=${c.plating.sigma}`;
            if (c.plating.thickness !== pl.thickness) extra += ` plating_t=${fmtSmall(c.plating.thickness)}`;
            if ((c.plating.rq ?? 0) !== (pl.rq ?? 0)) extra += ` plating_rq=${fmtSmall(c.plating.rq ?? 0)}`;
            if (c.plating.rq_interface !== pl.rq_interface) {
                extra += ` plating_rq_iface=${fmtSmall(c.plating.rq_interface ?? c.plating.rq ?? 0)}`;
            }
        }
        lines.push(`${kind} ${rectText(c)}` + extra);
    }
    return lines.join('\n') + '\n';
}

// A coaxial line as n-gons: the dielectric disk, the centre conductor and the shield
// as a ring of the given wall thickness. Each n-gon has the area of its circle. The
// domain is open around the shield. CoaxSolver builds its geometry from this text.
// solver: { a, b, n_inner, n_outer, shield_thickness, epsilon_r, tan_delta, plating },
// plating { sigma, thickness, rq, inner, outer } or null.
export function coaxToGeometryText(solver, units) {
    const scale = LENGTH_UNITS[units];
    const fmt = v => String(Number((v / scale).toPrecision(12)));
    // Vertex radius of the n-gon with the area of a circle of radius r.
    const R = polyRadiusForArea;
    const { a, b, n_inner: ni, n_outer: no } = solver;
    const c = b + solver.shield_thickness;
    // A plating of no thickness or conductivity has no effect and writes nothing.
    const pl = solver.plating && solver.plating.sigma > 0 && solver.plating.thickness > 0 ? solver.plating : null;
    const plated = which => (pl && pl[which] ? ' plating=all' : '');
    const lines = [
        `# Inner diameter ${fmt(2 * a)}, dielectric diameter ${fmt(2 * b)}, shield thickness ${fmt(c - b)} ${units}.`,
        '# Each circle is an n-gon of the same area: r1, r2 and r3 are the vertex radii of the',
        '# centre conductor, the dielectric and the outside of the shield.',
        `units ${units}`,
        `r1 = ${fmt(R(a, ni))}; r2 = ${fmt(R(b, no))}; r3 = ${fmt(R(c, no))}`,
        'bounds open open open open',
        'domain -1.1*r3 1.1*r3 -1.1*r3 1.1*r3',
    ];
    if (pl) lines.push(`plating sigma=${pl.sigma} t=${Number((pl.thickness / 1e-6).toPrecision(12))}um rq=${Number(((pl.rq ?? 0) / 1e-6).toPrecision(12))}um`
        + (pl.rq_interface !== undefined ? ` rq_iface=${Number((pl.rq_interface / 1e-6).toPrecision(12))}um` : ''));
    lines.push('',
        `diel  ngon  x=0  y=0  r=r2  n=${no}  er=${solver.epsilon_r}  tand=${solver.tan_delta ?? 0}`,
        `sig+  ngon  x=0  y=0  r=r1  n=${ni}${plated('inner')}`,
        `gnd   ngon  x=0  y=0  r=r3  r_in=r2  n=${no}${plated('outer')}`);
    return lines.join('\n') + '\n';
}

// --- Change of units ----------------------------------------------------------------
// Rewrites the text for another declared unit so that the geometry keeps its size. A
// bare number that is a length is scaled, one that is a factor (the 2 of s/2) is not.
// Which is which follows from the expression: in a sum every term is a length, in a
// product the numbers carry what the parameters leave over. Parameters are taken as
// lengths first, then with growing sets of them as plain factors, until the rewritten
// text evaluates to the same geometry. Numbers with their own unit stay as written.

// Expression tree: num (with its token), id, neg, sum, prod (factors with op), fn.
const expressionTree = src => parseExpression(src, {
    num: tk => ({ type: 'num', tk }),
    id: name => ({ type: 'id', name }),
    fn: (name, args) => ({ type: 'fn', name, args }),
    neg: node => ({ type: 'neg', node }),
    prod: factors => ({ type: 'prod', factors }),
    sum: terms => ({ type: 'sum', items: terms.map(t => t.node) }),
});

// Length degree of a node from its parameters, null when only its numbers could say.
function naturalDegree(node, degreeOf) {
    switch (node.type) {
    case 'num': return null;
    case 'id': return node.name === 'inf' ? null : degreeOf(node.name);
    case 'neg': return naturalDegree(node.node, degreeOf);
    case 'sum': case 'fn': {
        const args = node.type === 'sum' ? node.items : node.args;
        for (const a of args) {
            const d = naturalDegree(a, degreeOf);
            if (d !== null) return node.type === 'fn' && node.name === 'sqrt' ? d / 2 : d;
        }
        return null;
    }
    case 'prod': {
        let d = 0;
        for (const f of node.factors) {
            const fd = naturalDegree(f.node, degreeOf);
            if (fd === null) return null;
            d += f.op === '/' ? -fd : fd;
        }
        return d;
    }
    }
    return null;
}

// Collects the numbers of `node` to scale for it to have length degree d.
function collectScaled(node, d, degreeOf, out) {
    switch (node.type) {
    case 'num': if (node.tk.unit === undefined && d !== 0) out.push({ tk: node.tk, d }); return;
    case 'id': return;
    case 'neg': collectScaled(node.node, d, degreeOf, out); return;
    case 'sum': node.items.forEach(n => collectScaled(n, d, degreeOf, out)); return;
    case 'fn': node.args.forEach(n => collectScaled(n, node.name === 'sqrt' ? 2 * d : d, degreeOf, out)); return;
    case 'prod': {
        const known = node.factors.map(f => naturalDegree(f.node, degreeOf));
        let rest = d;
        node.factors.forEach((f, k) => { if (known[k] !== null) rest -= f.op === '/' ? -known[k] : known[k]; });
        let first = true;
        node.factors.forEach((f, k) => {
            let fd = known[k];
            if (fd === null) { fd = first ? (f.op === '/' ? -rest : rest) : 0; first = false; }
            collectScaled(f.node, fd, degreeOf, out);
        });
    }
    }
}

// `src` rewritten for length degree d with lengths scaled by f.
function scaleExpression(src, d, f, degreeOf) {
    const edits = [];
    collectScaled(expressionTree(src), d, degreeOf, edits);
    let out = src;
    for (const { tk, d: nd } of edits.sort((a, b) => b.tk.start - a.tk.start)) {
        out = out.slice(0, tk.start) + fmtNumber(tk.value * f ** nd) + out.slice(tk.end);
    }
    return out;
}

const RECT_LENGTH_KEYS = ['x', 'y', 'w', 'h', 'r', 'r_in', 'rx', 'ry', 'rx_in', 'ry_in', 'radius', 'radius_bottom', 'wall',
    'rq', 'plating_t', 'plating_rq', 'plating_rq_iface'];
const PLATING_LENGTH_KEYS = ['t', 'rq', 'rq_iface'];

// One statement's source rewritten.
function scaleStatement(src, st, f, to, degreeOf) {
    const lead = /^\s*/.exec(src)[0], trail = /\s*$/.exec(src)[0];
    const scaled = e => scaleExpression(e, 1, f, degreeOf);
    // key=value in place, keeping the spacing. A value followed by a separate unit
    // word ("rq=1 um") already has its unit.
    const keys = (code, list) => code.replace(/(^|\s)([A-Za-z_][A-Za-z_0-9]*)=(\S+)(\s+[A-Za-zµ]+(?=\s|$))?/g,
        (m, pre, key, value, unit) => (list.includes(key) && !(unit && LENGTH_UNITS[unit.trim()] !== undefined)
            ? `${pre}${key}=${scaled(value)}${unit ?? ''}` : m));
    switch (st.type) {
    case 'units': return `${lead}units ${to}${trail}`;
    case 'param': {
        const m = /^(\s*[A-Za-z_][A-Za-z_0-9]*\s*=\s*)(.*?)(\s*)$/.exec(src);
        return m[1] + scaleExpression(m[2], degreeOf(st.name), f, degreeOf) + m[3];
    }
    case 'domain':
        if (st.values.every(v => v === 'auto')) return src;
        return `${lead}domain ${st.values.map(v => (v === 'auto' ? v : scaled(v))).join(' ')}${trail}`;
    case 'plating': return keys(src, PLATING_LENGTH_KEYS);
    case 'rect': return keys(src, RECT_LENGTH_KEYS);
    }
    return src;
}

function rewriteUnits(text, model, f, to, degreeOf) {
    let hasUnits = false;
    const out = editStatements(text, (parts, comment, li) => {
        const sts = model.statements.filter(s => s.line === li + 1 && s.part !== undefined);
        for (const st of sts) {
            const k = partIndex(parts, st.part);
            if (k < 0) continue;
            if (st.type === 'units') hasUnits = true;
            parts[k] = scaleStatement(parts[k], st, f, to, degreeOf);
        }
        return sts.length > 0;
    }) ?? String(text ?? '');
    return hasUnits ? out : setStatementInText(out, 'units', `units ${to}`);
}

// Same geometry in metres, to a relative 1e-9. Source lines may differ.
function sameGeometry(a, b) {
    const eq = (x, y) => {
        if (typeof x === 'number' && typeof y === 'number') {
            return x === y || Math.abs(x - y) <= 1e-9 * Math.max(Math.abs(x), Math.abs(y));
        }
        if (x && y && typeof x === 'object' && typeof y === 'object') {
            const keys = new Set([...Object.keys(x), ...Object.keys(y)].filter(k => k !== 'line'));
            return [...keys].every(k => eq(x[k], y[k]));
        }
        return x === y;
    };
    return eq(a.rects, b.rects) && eq(a.domain, b.domain) && eq(a.plating, b.plating);
}

// The text with its declared unit changed to `to` and every length rewritten so the
// geometry keeps its size, or null when the text has errors or no consistent rewrite
// exists (then the user has to convert by hand).
export function changeUnitsInText(text, to) {
    if (LENGTH_UNITS[to] === undefined) throw new Error(`unknown unit '${to}'`);
    const model = parseGeometryText(text);
    const before = evaluateGeometry(model);
    if (before.errors.length) return null;
    const f = LENGTH_UNITS[before.units] / LENGTH_UNITS[to];
    const names = model.statements.filter(s => s.type === 'param').map(s => s.name);
    // Sets of plain-factor parameters, smallest first, at most a few thousand tries and
    // half a second: each try evaluates the whole text, on the UI thread.
    const MAX_TRIES = 4096, MAX_MS = 500;
    const t0 = Date.now();
    let tries = 0;
    const subsets = function* (start, size, chosen) {
        if (chosen.length === size) { yield new Set(chosen); return; }
        for (let k = start; k < names.length; k++) yield* subsets(k + 1, size, [...chosen, names[k]]);
    };
    for (let size = 0; size <= names.length; size++) {
        for (const factors of subsets(0, size, [])) {
            if (++tries > MAX_TRIES || Date.now() - t0 > MAX_MS) return null;
            const degreeOf = name => (names.includes(name) ? (factors.has(name) ? 0 : 1) : null);
            let out;
            try { out = rewriteUnits(text, model, f, to, degreeOf); } catch { continue; }
            const after = evaluateGeometry(parseGeometryText(out));
            if (!after.errors.length && after.units === to && sameGeometry(before, after)) return out;
        }
    }
    return null;
}
