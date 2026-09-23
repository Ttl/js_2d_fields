// Text format for custom cross-sections built from axis-aligned rectangles.
//
//   units um
//   w = 200; s = 150; t = 35; h1 = 100
//   bounds open open open gnd          # left right top bottom
//   domain auto                        # or: domain x1 x2 y1 y2, each may be auto
//   plating sigma=4.1e7 t=5 rq=0.1     # default plating material for rects with plating=
//
//   diel  x=-inf    y=0   w=inf  h=h1  er=4.3 tand=0.02
//   sig-  x=-s/2-w  y=h1  w=w    h=t
//   sig+  x=s/2     y=h1  w=w    h=t   plating=top,sides
//
// Statements are separated by newlines or ';'. '#' starts a comment. A rectangle is
// given per axis as x,w or x1,x2 (y,h or y1,y2). -inf / inf pins an edge to the domain
// wall. A negative h puts y at the top face (embedded trace convention of Conductor).
// Values in a rectangle, domain or plating statement are expressions without
// whitespace. A parameter definition takes the rest of the statement.
//
// Lengths are in the declared units (default mm). er, tand and sigma are plain numbers.
// A number may carry its own unit (35um or 35 um), which converts it to the declared units.
// A conductor may carry its own conductivity (sigma=, S/m), surface roughness (rq=) and
// plating material (plating_sigma=, plating_t=, plating_rq=), which override the
// solver-wide values and the plating statement for that conductor.
// thin=1 on a dielectric marks a thin sheet (a solder mask) whose faces the FDM mesher
// brackets with grid lines the way it does conductor faces.
// mirror=1 adds the mirror image about x=0 right after the rectangle, sig+ imaged as sig-
// and sig- as sig+. A rectangle that touches or crosses x=0 becomes one rectangle
// symmetric about x=0 instead.
// Dielectrics are painted in order, a later one overrides an earlier one where they
// overlap. Conductors override dielectrics.
//
// parseGeometryText keeps the expressions, evaluateGeometry turns them into metres.
// Parameter overrides go to evaluateGeometry, which is what parameter sweeps use.

export const LENGTH_UNITS = {
    m: 1, cm: 1e-2, mm: 1e-3, um: 1e-6, 'µm': 1e-6, nm: 1e-9, mil: 25.4e-6, in: 25.4e-3,
};

export const RECT_KINDS = ['diel', 'gnd', 'sig+', 'sig-'];
const BOUND_VALUES = ['open', 'gnd'];
const RESERVED = new Set(['inf', 'auto', 'min', 'max', 'abs', 'sqrt']);
const FUNCTIONS = {
    min: Math.min, max: Math.max, abs: Math.abs, sqrt: Math.sqrt,
};
const RECT_KEYS = new Set(['x', 'y', 'w', 'h', 'x1', 'x2', 'y1', 'y2', 'er', 'tand', 'thin', 'sigma', 'rq', 'plating',
    'plating_sigma', 'plating_t', 'plating_rq', 'mirror']);
const RECT_KEY_ORDER = ['x', 'x1', 'x2', 'w', 'y', 'y1', 'y2', 'h', 'er', 'tand', 'thin', 'sigma', 'rq', 'plating',
    'plating_sigma', 'plating_t', 'plating_rq', 'mirror'];
const MIRROR_KIND = { 'sig+': 'sig-', 'sig-': 'sig+', gnd: 'gnd', diel: 'diel' };
const PLATING_KEYS = new Set(['sigma', 't', 'rq', 'thick_corners']);
const PLATING_FACES = ['top', 'sides', 'bottom'];

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
        if (m[1] !== undefined) tokens.push({ type: 'num', value: parseFloat(m[1]), unit: m[2] });
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

// Recursive-descent evaluation. `vars` maps parameter names to numbers, `unitScale`
// is the declared unit in metres (a suffixed number is converted to it).
export function evaluateExpression(src, vars = {}, unitScale = 1) {
    const tokens = tokenize(src);
    let i = 0;
    const peek = () => tokens[i];
    const isOp = v => peek() && peek().type === 'op' && peek().value === v;

    function primary() {
        const tk = tokens[i++];
        if (!tk) throw new Error('unexpected end of expression');
        if (tk.type === 'num') {
            if (tk.unit === undefined) return tk.value;
            const scale = LENGTH_UNITS[tk.unit];
            if (scale === undefined) throw new Error(`unknown unit '${tk.unit}'`);
            return tk.value * scale / unitScale;
        }
        if (tk.type === 'id') {
            if (isOp('(')) {
                const fn = FUNCTIONS[tk.value];
                if (!fn) throw new Error(`unknown function '${tk.value}'`);
                i++;
                const args = [expr()];
                while (isOp(',')) { i++; args.push(expr()); }
                if (!isOp(')')) throw new Error("expected ')'");
                i++;
                return fn(...args);
            }
            if (tk.value === 'inf') return Infinity;
            if (!Object.prototype.hasOwnProperty.call(vars, tk.value)) {
                throw new Error(`unknown parameter '${tk.value}'`);
            }
            return vars[tk.value];
        }
        if (tk.value === '(') {
            const v = expr();
            if (!isOp(')')) throw new Error("expected ')'");
            i++;
            return v;
        }
        throw new Error(`unexpected '${tk.value}'`);
    }
    function unary() {
        if (isOp('-')) { i++; return -unary(); }
        if (isOp('+')) { i++; return unary(); }
        return primary();
    }
    function term() {
        let v = unary();
        while (isOp('*') || isOp('/')) {
            const op = tokens[i++].value;
            const r = unary();
            v = op === '*' ? v * r : v / r;
        }
        return v;
    }
    function expr() {
        let v = term();
        while (isOp('+') || isOp('-')) {
            const op = tokens[i++].value;
            const r = term();
            v = op === '+' ? v + r : v - r;
        }
        return v;
    }

    const v = expr();
    if (i < tokens.length) throw new Error(`unexpected '${tokens[i].value}'`);
    if (Number.isNaN(v)) throw new Error('expression is not a number');
    return v;
}

// --- Parsing ----------------------------------------------------------------------

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
        if (RESERVED.has(param[1]) || LENGTH_UNITS[param[1]] !== undefined) {
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
        return { type: 'plating', fields: parseFields(words.slice(1), PLATING_KEYS, 'plating') };
    }
    if (RECT_KINDS.includes(head)) {
        return { type: 'rect', kind: head, fields: parseFields(words.slice(1), RECT_KEYS, 'rectangle') };
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
        const hash = lines[li].indexOf('#');
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
            code = `plating ${fieldsToText(s.fields, ['sigma', 't', 'rq', 'thick_corners'])}`;
        } else if (s.type === 'rect') {
            code = `${s.kind} ${fieldsToText(s.fields, RECT_KEY_ORDER)}`;
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
    const hash = line.indexOf('#');
    return hash >= 0 ? [line.slice(0, hash), line.slice(hash)] : [line, ''];
}

// Replaces the expression of parameter `name`. Returns the text unchanged when the
// parameter is not defined.
export function setParamInText(text, name, expr) {
    const lines = String(text ?? '').split(/\r?\n/);
    const re = new RegExp(`^(\\s*${name}\\s*=\\s*)(.*?)(\\s*)$`);
    for (let i = 0; i < lines.length; i++) {
        const [code, comment] = splitComment(lines[i]);
        const parts = code.split(';');
        for (let k = 0; k < parts.length; k++) {
            const m = re.exec(parts[k]);
            if (!m) continue;
            parts[k] = m[1] + expr + m[3];
            lines[i] = parts.join(';') + comment;
            return lines.join('\n');
        }
    }
    return text;
}

// Renames parameter `from` to `to`: its definition and every expression that uses it.
// Field keys (the w of w=...) and comments are left alone.
export function renameParamInText(text, from, to) {
    const lines = String(text ?? '').split(/\r?\n/);
    const ident = new RegExp(`(?<![A-Za-z_0-9.])${from}(?![A-Za-z_0-9])`, 'g');
    const keyed = new RegExp(`(?<![A-Za-z_0-9.])${from}(?![A-Za-z_0-9])(?!\\s*=)`, 'g');
    for (let i = 0; i < lines.length; i++) {
        const [code, comment] = splitComment(lines[i]);
        const parts = code.split(';').map(part => {
            // A parameter definition has a bare name left of the first '='.
            const isParam = /^\s*[A-Za-z_][A-Za-z_0-9]*\s*=/.test(part);
            return part.replace(isParam ? ident : keyed, to);
        });
        lines[i] = parts.join(';') + comment;
    }
    return lines.join('\n');
}

// Statement text for a rectangle, the inverse of the parser for one statement.
export function rectStatementText(kind, fields) {
    return `${kind}  ${fieldsToText(fields, RECT_KEY_ORDER)}`;
}

// Replaces the statement `st` (from parseGeometryText of the same text) by `code`, or
// removes it when code is null. A line left empty by a removal is dropped.
export function replaceStatementInText(text, st, code) {
    const lines = String(text ?? '').split(/\r?\n/);
    const i = st.line - 1;
    if (i < 0 || i >= lines.length) return text;
    const [src, comment] = splitComment(lines[i]);
    const parts = src.split(';');
    // st.part counts non-empty statements, the split keeps the empty ones.
    let seen = -1, k = -1;
    for (let j = 0; j < parts.length; j++) {
        if (parts[j].trim().length === 0) continue;
        if (++seen === st.part) { k = j; break; }
    }
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
    const lines = String(text ?? '').split(/\r?\n/);
    const starts = s => s.trim() === keyword || s.trim().startsWith(keyword + ' ');
    let unitsLine = -1;
    for (let i = 0; i < lines.length; i++) {
        const [code, comment] = splitComment(lines[i]);
        const parts = code.split(';');
        const k = parts.findIndex(starts);
        if (k >= 0) {
            const lead = /^\s*/.exec(parts[k])[0];
            parts[k] = lead + statement + (k < parts.length - 1 || !comment ? '' : '  ');
            lines[i] = parts.join(';') + comment;
            return lines.join('\n');
        }
        if (unitsLine < 0 && parts.some(s => s.trim().startsWith('units '))) unitsLine = i;
    }
    lines.splice(unitsLine + 1, 0, statement);
    return lines.join('\n');
}

// --- Expression arithmetic for the form's edits --------------------------------------
// Builds new expressions from the ones written, as short as they reasonably get, so a
// switch between x,w and x1,x2 keeps the parameters instead of freezing numbers.

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

// -e with the sign of every top-level term flipped: -(-s/2-w) is s/2+w.
export function negateExpr(e) {
    e = e.trim();
    if (PLAIN_NUMBER.test(e)) return fmtNumber(-parseFloat(e));
    const cuts = [0, ...topLevelSigns(e), e.length];
    let out = '';
    for (let k = 0; k + 1 < cuts.length; k++) {
        let term = e.slice(cuts[k], cuts[k + 1]);
        let sign = '+';
        if (term.startsWith('-') || term.startsWith('+')) { sign = term[0]; term = term.slice(1); }
        const flipped = sign === '-' ? '+' : '-';
        out += (out === '' && flipped === '+') ? term : flipped + term;
    }
    return out;
}

export function addExpr(a, b) {
    a = a.trim(); b = b.trim();
    if (PLAIN_NUMBER.test(a) && PLAIN_NUMBER.test(b)) return fmtNumber(parseFloat(a) + parseFloat(b));
    if (PLAIN_NUMBER.test(a) && parseFloat(a) === 0) return b;
    if (PLAIN_NUMBER.test(b) && parseFloat(b) === 0) return a;
    if (b.startsWith('-') && !hasTopLevelSum(b)) return `${a}-${b.slice(1)}`;
    return `${a}+${b}`;
}

export function subExpr(a, b) {
    a = a.trim(); b = b.trim();
    if (a === b) return '0';
    if (PLAIN_NUMBER.test(a) && PLAIN_NUMBER.test(b)) return fmtNumber(parseFloat(a) - parseFloat(b));
    if (PLAIN_NUMBER.test(b) && parseFloat(b) === 0) return a;
    if (PLAIN_NUMBER.test(a) && parseFloat(a) === 0) return negateExpr(b);
    // (b+c)-b is c.
    if (a.startsWith(b + '+') && !hasTopLevelSum(b)) return a.slice(b.length + 1);
    if (hasTopLevelSum(b)) return `${a}-(${b})`;
    return b.startsWith('-') ? `${a}+${b.slice(1)}` : `${a}-${b}`;
}

const AXIS_KEYS = { x: ['x', 'w', 'x1', 'x2'], y: ['y', 'h', 'y1', 'y2'] };

// Which form an axis of a rectangle's fields is written in: 'size' (x,w) or 'bounds' (x1,x2).
export function axisForm(fields, axis) {
    const [, , lo, hi] = AXIS_KEYS[axis];
    return fields[lo] !== undefined || fields[hi] !== undefined ? 'bounds' : 'size';
}

// Low and high edge of an axis as expressions. `negative` says the size evaluates
// negative (y,h with h < 0, the top-face convention).
export function axisEdges(fields, axis, negative = false) {
    const [p, s, lo, hi] = AXIS_KEYS[axis];
    if (axisForm(fields, axis) === 'bounds') return { lo: fields[lo] ?? '0', hi: fields[hi] ?? '0' };
    const pos = fields[p] ?? '0', size = fields[s] ?? '0';
    if (isNeg(pos)) return { lo: '-inf', hi: 'inf' };
    if (isPosInf(size)) return { lo: pos, hi: 'inf' };
    return negative ? { lo: addExpr(pos, size), hi: pos } : { lo: pos, hi: addExpr(pos, size) };
}

// Writes the edges of an axis into `fields` in the given form. Returns false, leaving
// the fields alone, when the size form cannot hold them (a -inf low edge with a finite
// high edge).
function setAxisEdges(fields, axis, edges, form) {
    const [p, s, lo, hi] = AXIS_KEYS[axis];
    let out;
    if (form === 'bounds') out = { [lo]: edges.lo, [hi]: edges.hi };
    else if (isNeg(edges.lo)) {
        if (!isPosInf(edges.hi)) return false;
        out = { [p]: '-inf', [s]: 'inf' };
    } else if (isPosInf(edges.hi)) out = { [p]: edges.lo, [s]: 'inf' };
    else out = { [p]: edges.lo, [s]: subExpr(edges.hi, edges.lo) };
    for (const k of AXIS_KEYS[axis]) delete fields[k];
    Object.assign(fields, out);
    return true;
}

// Switches an axis between x,w and x1,x2. Returns false when it cannot be written in
// the other form.
export function toggleAxisForm(fields, axis, negative = false) {
    const edges = axisEdges(fields, axis, negative);
    return setAxisEdges(fields, axis, edges, axisForm(fields, axis) === 'bounds' ? 'size' : 'bounds');
}

// --- Evaluation -------------------------------------------------------------------

// One axis of a rectangle in metres: { pos, size, min, max }. pos and size are the
// constructor arguments of Conductor / Dielectric, kept as written when the axis is
// given as position + size so a geometry converted from a solver reproduces its
// doubles exactly. min / max are the bounds, -Infinity / Infinity for a pinned edge
// (pos and size are then resolved against the domain by the solver).
function evalAxis(fields, pos, size, lo, hi, ev, allowNegative) {
    const hasPS = fields[pos] !== undefined || fields[size] !== undefined;
    const hasLH = fields[lo] !== undefined || fields[hi] !== undefined;
    if (hasPS && hasLH) throw new Error(`give either ${pos},${size} or ${lo},${hi}, not both`);
    if (hasPS) {
        if (fields[pos] === undefined || fields[size] === undefined) {
            throw new Error(`${pos} and ${size} must both be given`);
        }
        const p = ev(fields[pos]), s = ev(fields[size]);
        if (p === Infinity) throw new Error(`${pos} cannot be inf`);
        if (s === -Infinity) throw new Error(`${size} cannot be -inf`);
        if (p === -Infinity && s !== Infinity) {
            throw new Error(`${pos}=-inf needs ${size}=inf, use ${lo},${hi} for a wall-pinned edge`);
        }
        if (s === 0) throw new Error(`${size} must be nonzero`);
        if (s < 0 && !allowNegative) throw new Error(`${size} must be positive`);
        if (s < 0 && p === -Infinity) throw new Error(`negative ${size} cannot be combined with inf`);
        if (s < 0) return { pos: p, size: s, min: p + s, max: p };
        return { pos: p, size: s, min: p, max: s === Infinity ? Infinity : p + s };
    }
    if (fields[lo] === undefined || fields[hi] === undefined) {
        throw new Error(`missing ${pos},${size} or ${lo},${hi}`);
    }
    const a = ev(fields[lo]), b = ev(fields[hi]);
    if (a === Infinity || b === -Infinity || !(b > a)) throw new Error(`${hi} must be greater than ${lo}`);
    return { pos: a, size: b - a, min: a, max: b };
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
//   plating - { sigma, thickness, rq, thick_corners } or null
//   rects   - { kind, x: axis, y: axis, er, tand, thin, plating: faces|null, sigma: S/m|null, rq: m|null,
//               platingMaterial: { sigma?, thickness?, rq? }|null, line, image }
//             image is true on a rectangle generated by mirror=1
export function evaluateGeometry(model, overrides = {}) {
    const errors = [...model.errors];
    const unitsSt = model.statements.find(s => s.type === 'units');
    const units = unitsSt ? unitsSt.value : 'mm';
    const scale = LENGTH_UNITS[units];
    const params = {};
    const fail = (s, e) => errors.push({ line: s.line, message: e.message ?? String(e) });
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
                thick_corners: f.thick_corners !== undefined ? num(f.thick_corners) !== 0 : false,
            };
            if (!(plating.sigma > 0) || !(plating.thickness > 0) || !(plating.rq >= 0)) {
                throw new Error('plating sigma and t must be positive, rq non-negative');
            }
        } catch (e) { fail(platingSt, e); }
    }

    const rects = [];
    for (const s of model.statements) {
        if (s.type !== 'rect') continue;
        try {
            const f = s.fields;
            const r = {
                kind: s.kind, line: s.line,
                x: evalAxis(f, 'x', 'w', 'x1', 'x2', len, false),
                y: evalAxis(f, 'y', 'h', 'y1', 'y2', len, true),
                er: 1, tand: 0, thin: false, plating: null, sigma: null, rq: null, platingMaterial: null,
            };
            if (s.kind === 'diel') {
                if (f.er === undefined) throw new Error('diel needs er');
                if (['plating', 'sigma', 'rq', 'plating_sigma', 'plating_t', 'plating_rq'].some(k => f[k] !== undefined)) {
                    throw new Error('sigma, rq and plating apply to conductors only');
                }
                r.er = num(f.er);
                r.tand = f.tand !== undefined ? num(f.tand) : 0;
                r.thin = f.thin !== undefined ? num(f.thin) !== 0 : false;
                if (!(r.er >= 1) || !Number.isFinite(r.er)) throw new Error('er must be a finite number >= 1');
                if (!(r.tand >= 0) || !Number.isFinite(r.tand)) throw new Error('tand must be non-negative');
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
                if (Object.keys(pm).length) {
                    if ((pm.sigma !== undefined && !(pm.sigma > 0)) || (pm.thickness !== undefined && !(pm.thickness > 0))
                        || (pm.rq !== undefined && !(pm.rq >= 0))) {
                        throw new Error('plating_sigma and plating_t must be positive, plating_rq non-negative');
                    }
                    r.platingMaterial = pm;
                }
            }
            r.image = false;
            if (f.mirror !== undefined && num(f.mirror) !== 0) {
                if (r.x.min <= 0 && r.x.max >= 0) {
                    const m = Math.max(-r.x.min, r.x.max);
                    r.x = m === Infinity ? { pos: -Infinity, size: Infinity, min: -Infinity, max: Infinity }
                        : { pos: -m, size: 2 * m, min: -m, max: m };
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
    if (all.some(o => o.shape)) throw new Error('Only rectangular geometries can be converted.');
    const fmt = v => (units === 'm' ? String(v) : String(Number((v / scale).toPrecision(12))));
    const X0 = -solver.domain_width / 2, X1 = solver.domain_width / 2;
    const Y0 = solver.domain_y_min, Y1 = solver.domain_height;
    const tol = Math.max(X1 - X0, Y1 - Y0) * 1e-9;

    const axisText = (pos, size, lo, hi, p, s, d0, d1) => {
        const a = s >= 0 ? p : p + s, b = s >= 0 ? p + s : p;
        const loWall = pinWalls && s > 0 && Math.abs(a - d0) <= tol;
        const hiWall = pinWalls && s > 0 && Math.abs(b - d1) <= tol;
        if (!loWall && !hiWall) return `${pos}=${fmt(p)} ${size}=${fmt(s)}`;
        if (loWall && hiWall) return `${pos}=-inf ${size}=inf`;
        if (hiWall) return `${pos}=${fmt(p)} ${size}=inf`;
        return `${lo}=-inf ${hi}=${fmt(b)}`;
    };
    const rectText = o => `${axisText('x', 'w', 'x1', 'x2', o.x, o.width, X0, X1)} ` +
        `${axisText('y', 'h', 'y1', 'y2', o.y, o.height, Y0, Y1)}`;

    const lines = [`units ${units}`, `bounds ${(solver.boundaries || ['open', 'open', 'open', 'gnd']).join(' ')}`,
        `domain ${fmt(X0)} ${fmt(X1)} ${fmt(Y0)} ${fmt(Y1)}`];
    // The first plating becomes the plating statement, a conductor plated with another
    // material carries its own keys.
    const pl = (solver.conductors || []).map(c => c.plating).find(p => p && p.sigma > 0 && p.thickness > 0);
    if (pl) {
        lines.push(`plating sigma=${pl.sigma} t=${fmt(pl.thickness)} rq=${fmt(pl.rq ?? 0)}` +
            (pl.thick_corners ? ' thick_corners=1' : ''));
    }
    lines.push('');
    for (const d of (solver.dielectrics || [])) {
        lines.push(`diel ${rectText(d)} er=${d.epsilon_r} tand=${d.tan_delta ?? 0}` + (d.thin_sheet ? ' thin=1' : ''));
    }
    for (const c of (solver.conductors || [])) {
        const kind = c.is_signal ? (c.polarity < 0 ? 'sig-' : 'sig+') : 'gnd';
        const faces = c.plating ? PLATING_FACES.filter(f => c.plating[f]) : [];
        let extra = (c.sigma !== undefined && c.sigma !== null) ? ` sigma=${c.sigma}` : '';
        if (c.rq !== undefined && c.rq !== null) extra += ` rq=${fmt(c.rq)}`;
        if (faces.length) {
            extra += ` plating=${faces.join(',')}`;
            if (c.plating.sigma !== pl.sigma) extra += ` plating_sigma=${c.plating.sigma}`;
            if (c.plating.thickness !== pl.thickness) extra += ` plating_t=${fmt(c.plating.thickness)}`;
            if ((c.plating.rq ?? 0) !== (pl.rq ?? 0)) extra += ` plating_rq=${fmt(c.plating.rq ?? 0)}`;
        }
        lines.push(`${kind} ${rectText(c)}` + extra);
    }
    return lines.join('\n') + '\n';
}
