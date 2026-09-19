// Text format for custom cross-sections built from axis-aligned rectangles.
//
//   units um
//   w = 200; s = 150; t = 35; h1 = 100
//   bounds open open open gnd          # left right top bottom
//   domain auto                        # or: domain x1 x2 y1 y2, each may be auto
//   plating sigma=4.1e7 t=5 rq=0.1     # plating material for rects with plating=
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
// A number may carry its own unit (35um), which converts it to the declared units.
// thin=1 on a dielectric marks a thin sheet (a solder mask) whose faces the FDM mesher
// brackets with grid lines the way it does conductor faces.
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
const RECT_KEYS = new Set(['x', 'y', 'w', 'h', 'x1', 'x2', 'y1', 'y2', 'er', 'tand', 'thin', 'plating']);
const PLATING_KEYS = new Set(['sigma', 't', 'rq', 'thick_corners']);
const PLATING_FACES = ['top', 'sides', 'bottom'];

// --- Expressions ------------------------------------------------------------------

function tokenize(src) {
    const tokens = [];
    const re = /\s*(?:(\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)([A-Za-zµ]+)?|([A-Za-z_][A-Za-z_0-9]*)|([-+*/(),]))/y;
    let pos = 0;
    while (pos < src.length) {
        if (/^\s*$/.test(src.slice(pos))) break;
        re.lastIndex = pos;
        const m = re.exec(src);
        if (!m) throw new Error(`unexpected character '${src.slice(pos).trim()[0]}'`);
        if (m[1] !== undefined) tokens.push({ type: 'num', value: parseFloat(m[1]), unit: m[2] });
        else if (m[3] !== undefined) tokens.push({ type: 'id', value: m[3] });
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
    for (const word of words) {
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
            code = `${s.kind} ${fieldsToText(s.fields,
                ['x', 'x1', 'x2', 'w', 'y', 'y1', 'y2', 'h', 'er', 'tand', 'thin', 'plating'])}`;
        }
        const comment = s.comment !== null && s.comment !== undefined ? `# ${s.comment}` : '';
        out.push([code, comment].filter(p => p.length > 0).join('  '));
    }
    while (out.length && out[out.length - 1] === '') out.pop();
    return out.join('\n') + '\n';
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

// Evaluates a parsed model to numbers in metres.
//   overrides - { name: value } replaces parameter values (in the declared units)
// Returns { errors, units, params, bounds, domain, plating, rects }:
//   domain  - { x_min, x_max, y_min, y_max }, null for auto
//   plating - { sigma, thickness, rq, thick_corners } or null
//   rects   - { kind, x: axis, y: axis, er, tand, thin, plating: faces|null, line }
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
        } catch (e) { fail(s, e); }
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
                er: 1, tand: 0, thin: false, plating: null,
            };
            if (s.kind === 'diel') {
                if (f.er === undefined) throw new Error('diel needs er');
                if (f.plating !== undefined) throw new Error('plating applies to conductors only');
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
        lines.push(`${kind} ${rectText(c)}` + (faces.length ? ` plating=${faces.join(',')}` : ''));
    }
    return lines.join('\n') + '\n';
}
