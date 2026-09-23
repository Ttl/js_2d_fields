// Editor for the custom geometry type. The sidebar holds the scalar inputs (units,
// parameters, boundaries, solved region), the panel under the geometry preview holds the
// rectangles: a form view with one row per rectangle and a text view of the whole
// geometry, plus templates and file / clipboard transfer. The text is the single source
// of truth. Every control rewrites one statement in it and leaves the rest as typed.
import { parseGeometryText, evaluateGeometry, setParamInText, renameParamInText, setStatementInText, replaceStatementInText,
         insertLineInText, moveRectInText, rectStatementText, formatErrors, evaluateExpression, LENGTH_UNITS,
         axisEdges, addExpr, formatLength, plausibilityWarnings, isPlainNumber, isLengthLiteral, isReservedName,
         changeUnitsInText } from './custom_geometry_text.js';
import { CustomGeometrySolver } from './custom_geometry.js';

export const CUSTOM_TEMPLATES = {
    'Stacked-dielectric microstrip': `# Microstrip on two dielectric layers. The ground boundary below the
# substrate is the ground plane, the air above is the rest of the domain.
units mm
w = 0.3; t = 0.035
h1 = 0.1; h2 = 0.4
bounds open open open gnd
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-inf  y=0      w=inf  h=h1  er=4.3  tand=0.02
diel  x=-inf  y=h1     w=inf  h=h2  er=2.2  tand=0.001
sig+  x=-w/2  y=h1+h2  w=w    h=t
`,
    'Stacked-dielectric differential microstrip': `# Differential microstrip on two dielectric layers
units mm
w = 0.3; s = 0.2; t = 0.035
h1 = 0.1; h2 = 0.4
bounds open open open gnd
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-inf    y=0      w=inf  h=h1  er=4.3  tand=0.02
diel  x=-inf    y=h1     w=inf  h=h2  er=2.2  tand=0.001
sig-  x=-s/2-w  y=h1+h2  w=w    h=t
sig+  x=s/2     y=h1+h2  w=w    h=t
`,
    'Microstrip on a finite ground': `# Finite substrate and ground, air and open boundaries all around
units mm
w = 0.35; t = 0.035; h = 0.2104
wsub = 3; wgnd = 0.5
bounds open open open open
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-wsub/2  y=0   w=wsub  h=h  er=4.4  tand=0.02
sig+  x=-w/2     y=h   w=w     h=t
gnd   x=-wgnd/2  y=-t  w=wgnd  h=t
`,
    'CPW over air': `# Coplanar waveguide on a finite substrate, air above and below
units mm
w = 0.2; g = 0.1; t = 0.017; h = 0.635
wgnd = 1.5; wsub = 5
bounds open open open open
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-wsub/2      y=-h  w=wsub  h=h  er=9.8  tand=0.001
gnd   x=-w/2-g-wgnd  y=0   w=wgnd  h=t
gnd   x=w/2+g        y=0   w=wgnd  h=t
sig+  x=-w/2         y=0   w=w     h=t
`,
    'Differential CPW over air': `# Differential coplanar waveguide on a finite substrate, air above and below
units mm
w = 0.2; s = 0.15; g = 0.1; t = 0.017; h = 0.635
wgnd = 1.5; wsub = 5
bounds open open open open
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-wsub/2        y=-h  w=wsub  h=h  er=9.8  tand=0.001
gnd   x=-s/2-w-g-wgnd  y=0   w=wgnd  h=t
gnd   x=s/2+w+g        y=0   w=wgnd  h=t
sig-  x=-s/2-w         y=0   w=w     h=t
sig+  x=s/2            y=0   w=w     h=t
`,
    'Coplanar strips / slotline': `# Two strips on a finite substrate, the slot between them is a slotline.
# One strip is the signal, the other the return conductor.
units mm
w = 1; s = 0.3; t = 0.05; h = 0.2104; wsub = 4
bounds open open open open
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-wsub/2  y=-h  w=wsub  h=h  er=4.4  tand=0.02
sig+  x=-s/2-w   y=0   w=w     h=t
gnd   x=s/2      y=0   w=w     h=t
`,
};

const WALLS = ['left', 'right', 'top', 'bottom'];
const KIND_LABELS = { 'sig+': 'Signal (+)', 'sig-': 'Signal (−)', 'gnd': 'Ground', 'diel': 'Dielectric' };
const FACES = ['top', 'sides', 'bottom'];
const $ = id => document.getElementById(id);
const paramInputId = name => `inp_cgp_${name}`;
// Value of a number with an optional unit in the declared units, NaN when it is not one.
// Parameters defined this way are the ones a sweep can vary.
const literalValue = (expr, unitScale) => {
    if (!isLengthLiteral(expr)) return NaN;
    try { return evaluateExpression(expr, {}, unitScale); } catch { return NaN; }
};
const fmt = v => (Number.isFinite(v) ? String(parseFloat(v.toPrecision(6))) : (v > 0 ? 'inf' : '-inf'));

let onChange = () => {};
let debounceTimer = null;
// Highlight changes arrive in pairs (focus leaves one row and enters the next), so the
// redraw they ask for is coalesced to one per frame.
let highlightFrame = 0;
function redrawSoon() {
    if (highlightFrame) return;
    highlightFrame = requestAnimationFrame(() => { highlightFrame = 0; onChange(); });
}
let lastParamSignature = null;

function el(tag, props = {}, ...children) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(props)) {
        if (k === 'class') e.className = v;
        else if (k === 'text') e.textContent = v;
        else if (k.startsWith('on')) e.addEventListener(k.slice(2), v);
        else if (v !== undefined && v !== null) e.setAttribute(k, v);
    }
    e.append(...children);
    return e;
}

// A modal question with a button per choice. Resolves to the chosen value, null when
// the dialog is dismissed (Escape).
function askChoice(title, message, choices) {
    return new Promise(resolve => {
        const dialog = el('dialog', { class: 'custom-choice' });
        const done = value => { dialog.close(); dialog.remove(); resolve(value); };
        dialog.append(el('div', { class: 'help-modal-content' },
            el('h3', { text: title }),
            ...message.split('\n').map(p => el('p', { text: p })),
            el('div', { class: 'custom-choice-buttons' }, ...choices.map(c =>
                el('button', { type: 'button', class: c.primary ? '' : 'secondary-btn', text: c.label,
                    onclick: () => done(c.value) })))));
        dialog.addEventListener('cancel', (e) => { e.preventDefault(); done(null); });
        document.body.append(dialog);
        dialog.showModal();
    });
}

// Replacing the text loses nothing when it is empty or an untouched template.
function replaceable() {
    const t = $('custom_geom_text');
    return !t.value.trim() || t.value === t.dataset.loaded;
}

export function getCustomGeometryText() {
    const t = $('custom_geom_text');
    return t ? t.value : '';
}

// Parsed and evaluated text, kept until the text changes: one edit reads it from the
// sidebar, the form, the messages, the sweep list and the preview. The model is shared,
// so callers do not modify it.
let analysed = { text: null };
function analyse(text = getCustomGeometryText()) {
    if (analysed.text !== text) {
        const model = parseGeometryText(text);
        analysed = { text, model, geo: evaluateGeometry(model), validation: null };
    }
    return analysed;
}

// Parameter values typed into the sidebar that differ from the text. Normally none: an
// edit rewrites the text. A parameter sweep sets the input directly, which is how the
// swept value reaches the solver.
export function getCustomOverrides() {
    const overrides = {};
    const { model } = analyse();
    const unitsSt = model.statements.find(s => s.type === 'units');
    const unitScale = LENGTH_UNITS[unitsSt ? unitsSt.value : 'mm'];
    for (const s of model.statements) {
        if (s.type !== 'param' || !isLengthLiteral(s.expr)) continue;
        const input = $(paramInputId(s.name));
        if (!input) continue;
        const v = literalValue(input.value, unitScale);
        if (Number.isFinite(v) && v !== literalValue(s.expr, unitScale)) overrides[s.name] = v;
    }
    return overrides;
}

// Sweepable geometry parameters: those defined as a number, with or without a unit.
export function customSweepParams() {
    return analyse().model.statements.filter(s => s.type === 'param' && isLengthLiteral(s.expr))
        .map(s => ({ key: `cgp_${s.name}`, label: `${s.name} (geometry parameter)`, inputId: paramInputId(s.name) }));
}

// Checks the text: parse and evaluation errors with their line, then the solver's own
// validation (shorts, missing ground, conductors outside the domain). Returns
// { errors, warnings, solver }, the solver only when the geometry is valid.
export function validateCustomGeometry(text = getCustomGeometryText()) {
    const a = analyse(text);
    if (!a.validation) {
        a.validation = validate(a.geo);
        if (!a.validation.errors.length) {
            a.validation.warnings.push(...plausibilityWarnings(a.geo, a.model)
                .map(w => (w.line > 0 ? `line ${w.line}: ` : '') + w.message));
        }
    }
    return a.validation;
}

function validate(geo) {
    if (geo.errors.length) return { errors: geo.errors, warnings: [], solver: null };
    try {
        const solver = new CustomGeometrySolver({ geometry: geo, nx: 10, ny: 10 });
        return { errors: [], warnings: solver.openBoundaryWarnings(), solver };
    } catch (e) {
        const m = /line (\d+)/.exec(e.message);
        return { errors: [{ line: m ? parseInt(m[1], 10) : 0, message: e.message.replace(/^line \d+: /, '') }],
                 warnings: [], solver: null };
    }
}

// A button per parameter that `message` reports as unknown, which defines it. A
// parameter used by the parameter on line `beforeLine` is defined above that line.
function unknownParamFixes(message, beforeLine = 0) {
    const names = [...new Set([...message.matchAll(/unknown parameter '([A-Za-z_][A-Za-z_0-9]*)'/g)].map(m => m[1]))]
        .filter(name => !isReservedName(name));
    return names.map(name => el('button', { type: 'button', class: 'secondary-btn custom-quick-fix',
        text: `+ parameter ${name}`, title: `Define ${name} in the parameters`,
        onclick: (e) => { e.stopPropagation(); addParameter(name, beforeLine); } }));
}

// Shows `message` under a sidebar parameter or a form row, or removes it.
function setInlineError(row, message, beforeLine = 0) {
    row.classList.toggle('has-error', !!message);
    let note = row.querySelector(':scope > .custom-inline-error');
    if (!message) { if (note) note.remove(); return; }
    if (!note) { note = el('div', { class: 'custom-inline-error' }); row.append(note); }
    note.textContent = message;
    note.append(...unknownParamFixes(message, beforeLine));
}

// Puts the text cursor on a source line and selects it.
function selectTextLine(line) {
    const t = $('custom_geom_text');
    const lines = t.value.split('\n');
    const start = lines.slice(0, line - 1).reduce((n, l) => n + l.length + 1, 0);
    t.focus();
    t.setSelectionRange(start, start + (lines[line - 1] || '').length);
    // Scroll the line to about a third down the visible part.
    const lineHeight = parseFloat(getComputedStyle(t).lineHeight) || 19;
    t.scrollTop = Math.max(0, (line - 1) * lineHeight - t.clientHeight / 3);
    window.customHighlightLine = line;
    redrawSoon();
}

// Line numbers of the text view, error lines marked. Rebuilt only when the line count
// or the error lines change.
let gutterKey = '';
let errorLines = new Set();
function renderLineNumbers() {
    const gutter = $('custom-line-numbers'), t = $('custom_geom_text');
    if (!gutter || !t) return;
    const count = t.value.split('\n').length;
    const key = `${count}|${[...errorLines].join(',')}`;
    if (key !== gutterKey) {
        gutterKey = key;
        gutter.replaceChildren(...Array.from({ length: count }, (_, i) => {
            const n = i + 1;
            return errorLines.has(n) ? el('span', { class: 'err', text: `${n}\n` }) : `${n}\n`;
        }));
    }
    gutter.scrollTop = t.scrollTop;
}

// Errors are shown where they are made: under the sidebar parameter or the form row of
// their line. The box under the editor lists the rest, and in the text view all of them,
// where a click selects the line.
function renderMessages({ errors, warnings }) {
    const attached = new Set();
    // An evaluation error names its parameter. Other errors go by the line, and when
    // several parameters share it the message has to name the parameter.
    const params = analyse().model.statements.filter(s => s.type === 'param');
    document.querySelectorAll('#custom-param-list [data-param]').forEach(row => {
        const st = params.find(p => p.name === row.dataset.param);
        const alone = st && params.filter(p => p.line === st.line).length === 1;
        const mine = st ? errors.filter(e => e.param ? e.param === st.name
            : e.line === st.line && (alone || e.message.includes(`'${st.name}'`))) : [];
        mine.forEach(e => attached.add(e));
        setInlineError(row, mine.map(e => e.message).join('\n'), st ? st.line : 0);
    });
    const inForm = formVisible();
    document.querySelectorAll('#custom-form [data-line]').forEach(row => {
        const mine = errors.filter(e => String(e.line) === row.dataset.line);
        if (inForm) mine.forEach(e => attached.add(e));
        setInlineError(row, mine.map(e => e.message).join('\n'));
    });

    const listed = inForm ? errors.filter(e => !attached.has(e)) : errors;
    const box = $('custom-geom-errors');
    if (box) {
        box.innerHTML = '';
        for (const e of listed) {
            const item = el('div', { class: 'custom-error-item', text: formatErrors([e]) },
                ...unknownParamFixes(e.message, e.param ? e.line : 0));
            if (e.line > 0) {
                item.classList.add('clickable');
                item.title = 'Show the line in the text view';
                item.addEventListener('click', () => { showView('text'); selectTextLine(e.line); });
            }
            box.append(item);
        }
        box.style.display = listed.length ? 'block' : 'none';
    }
    // Error count in the editor bar: the row with the error may be scrolled out of view.
    errorLines = new Set(errors.map(e => e.line).filter(l => l > 0));
    renderLineNumbers();
    const count = $('custom-error-count');
    if (count) {
        count.textContent = errors.length ? `${errors.length} error${errors.length > 1 ? 's' : ''}` : '';
        count.style.display = errors.length ? '' : 'none';
    }
    const warn = $('custom-geom-warnings');
    if (warn) {
        warn.textContent = errors.length ? '' : warnings.map(w => '⚠ ' + w).join('\n');
        warn.style.display = (!errors.length && warnings.length) ? 'block' : 'none';
    }
    // The preview keeps the last valid geometry, which the badge over it says.
    document.body.classList.toggle('custom-invalid', errors.length > 0);
    const btn = $('btn_solve');
    if (btn && $('tl_type').value === 'custom') {
        btn.dataset.customInvalid = errors.length ? '1' : '';
        if (!btn.classList.contains('stop-mode')) btn.disabled = errors.length > 0;
    }
}

// Brings the first error into view: its form row or sidebar parameter, else its text line.
function showFirstError() {
    const mark = document.querySelector('#custom-param-list .has-error')
        || (formVisible() && document.querySelector('#custom-form .has-error'));
    if (mark) {
        mark.scrollIntoView({ block: 'center', behavior: 'smooth' });
        const field = mark.querySelector('input');
        if (field) field.focus({ preventScroll: true });
        return;
    }
    const first = validateCustomGeometry().errors.find(e => e.line > 0);
    if (first) { showView('text'); selectTextLine(first.line); }
}

// --- Sidebar ------------------------------------------------------------------------

const PARAM_NAME = /^[A-Za-z_][A-Za-z_0-9]*$/;
const DOMAIN_KEYS = ['x1', 'x2', 'y1', 'y2'];

// Statement of parameter `name` in the current text.
function currentParam(name) {
    return analyse().model.statements.find(s => s.type === 'param' && s.name === name);
}

function setText(newText, source = 'sidebar-structure') {
    $('custom_geom_text').value = newText;
    refresh(source);
}

// The name of a parameter row: a click turns it into a field, Enter or leaving it renames
// the parameter everywhere it is used.
function paramNameCell(name, names) {
    const label = el('button', { type: 'button', class: 'custom-param-name-btn', title: 'Rename', text: name });
    label.addEventListener('click', () => {
        const input = el('input', { type: 'text', class: 'custom-param-rename', value: name, spellcheck: 'false',
            autocomplete: 'off', autocapitalize: 'off' });
        let done = false;
        const finish = (commit) => {
            if (done) return;
            done = true;
            const to = input.value.trim();
            const problem = !commit || to === name ? null
                : !PARAM_NAME.test(to) ? `'${to}' is not a name: letters, digits and _, not starting with a digit.`
                : isReservedName(to) ? `'${to}' is reserved (a function, a unit, inf or auto).`
                : names.includes(to) ? `'${to}' is already a parameter.` : null;
            if (commit && to !== name && !problem) {
                setText(renameParamInText(getCustomGeometryText(), name, to));
            } else {
                input.replaceWith(label);
                const row = label.closest('.custom-param-row');
                if (row && problem) setInlineError(row, `Not renamed: ${problem}`);
            }
        };
        input.addEventListener('keydown', (e) => {
            if (e.key === 'Enter') finish(true);
            else if (e.key === 'Escape') finish(false);
        });
        input.addEventListener('blur', () => finish(true));
        label.replaceWith(input);
        input.select();
    });
    return label;
}

function sidebarParamRow(st, names) {
    const input = el('input', { type: 'text', id: paramInputId(st.name), inputmode: 'decimal', spellcheck: 'false',
        autocomplete: 'off', autocapitalize: 'off', title: 'A number or an expression of the parameters above it' });
    input.addEventListener('input', () => {
        const expr = input.value.trim();
        if (!expr) return;
        // ';' and '#' would end the statement in the text.
        if (/[;#]/.test(expr)) {
            setInlineError(input.parentElement, 'A value cannot contain ; or #.');
            return;
        }
        const t = $('custom_geom_text');
        t.value = setParamInText(t.value, st.name, expr);
        scheduleChange('sidebar');
    });
    const del = el('button', { type: 'button', class: 'secondary-btn custom-row-btn', title: 'Delete the parameter', text: '✕' });
    del.addEventListener('click', () => {
        const cur = currentParam(st.name);
        if (cur) setText(replaceStatementInText(getCustomGeometryText(), cur, null));
    });
    return el('div', { class: 'custom-param-row', 'data-param': st.name },
        paramNameCell(st.name, names), input, del,
        el('div', { class: 'custom-hint custom-param-value' }));
}

// Defines a parameter, `name` or the first free p1, p2, ..., with the value 1. It goes
// after the last parameter, or above line `beforeLine` when that is given.
function addParameter(name = null, beforeLine = 0) {
    const { model } = analyse();
    const params = model.statements.filter(s => s.type === 'param');
    const names = new Set(params.map(p => p.name));
    if (!name) { let n = 1; while (names.has(`p${n}`)) n++; name = `p${n}`; }
    if (names.has(name)) return;
    const lastLine = params.length ? Math.max(...params.map(p => p.line))
        : Math.max(0, ...model.statements.filter(s => s.type === 'units').map(s => s.line));
    setText(insertLineInText(getCustomGeometryText(), beforeLine > 0 ? beforeLine - 1 : lastLine, `${name} = 1`));
    const input = $(paramInputId(name));
    if (input) { input.focus(); input.select(); }
}

function renderSidebarParams(model, geo) {
    const list = $('custom-param-list');
    if (!list) return;
    const params = model.statements.filter(s => s.type === 'param');
    const names = params.map(s => s.name);
    const signature = names.join('|');
    if (signature !== lastParamSignature) {
        lastParamSignature = signature;
        list.innerHTML = '';
        if (!params.length) {
            list.append(el('div', { class: 'custom-hint', text: 'No parameters defined.' }));
        }
        for (const s of params) list.append(sidebarParamRow(s, names));
    }
    for (const s of params) {
        const input = $(paramInputId(s.name));
        if (!input) continue;
        if (input !== document.activeElement) input.value = s.expr.trim();
        // The value of an expression is shown under its field.
        const v = geo.params[s.name];
        const shown = isPlainNumber(s.expr) || v === undefined ? '' : `= ${fmt(v)}`;
        const hint = input.parentElement.querySelector('.custom-param-value');
        hint.textContent = shown;
        hint.style.display = shown ? 'block' : 'none';
    }
    const unitSel = $('custom-units');
    if (unitSel && unitSel !== document.activeElement) unitSel.value = geo.units;
}

function renderDomain(model) {
    const st = model.statements.find(s => s.type === 'domain');
    DOMAIN_KEYS.forEach((k, i) => {
        const input = $(`custom_domain_${k}`);
        if (!input || input === document.activeElement) return;
        const v = st ? st.values[i] : 'auto';
        input.value = v === 'auto' ? '' : v;
    });
}

function renderBounds(geo) {
    WALLS.forEach((wall, i) => {
        const sel = $(`custom_bound_${wall}`);
        if (sel && sel !== document.activeElement) sel.value = geo.bounds[i];
    });
}

// --- Form view ----------------------------------------------------------------------

function applyText(newText, structural) {
    $('custom_geom_text').value = newText;
    if (structural) { refresh('form-structure'); return; }
    scheduleChange('form');
}

// An expression input. `commit` writes a new expression. `value` says how the evaluated
// value is shown under the field: 'length' and 'number' when the field holds more than a
// plain number, 'small' (roughness, plating thickness) always, in a unit that suits it.
function exprInput(label, value, commit, { placeholder = '', cls = '', title = '', kind = 'length' } = {}) {
    const input = el('input', { type: 'text', class: `custom-expr ${cls}`, value: value ?? '', placeholder,
        title: title || undefined, spellcheck: 'false', autocomplete: 'off', autocapitalize: 'off' });
    input.dataset.value = kind;
    // Values in a statement cannot hold spaces ("1 um", "w + 2"): the text gets the
    // value without them, the field keeps what was typed.
    input.addEventListener('input', () => commit(input.value.replace(/\s+/g, '')));
    return el('label', { class: 'custom-cell' }, el('span', { class: 'custom-cell-label', text: label }),
        el('span', { class: 'custom-cell-field' }, input, el('span', { class: 'custom-cell-value' })));
}

// Evaluated values under the form's fields, for the current parameters.
function updateFormValues(geo) {
    const k = LENGTH_UNITS[geo.units];
    document.querySelectorAll('#custom-form input.custom-expr').forEach(input => {
        const out = input.parentElement.querySelector('.custom-cell-value');
        if (!out) return;
        const expr = input.value.replace(/\s+/g, '');
        const kind = input.dataset.value;
        let text = '';
        if (expr && (kind === 'small' || !isPlainNumber(expr))) {
            let v = NaN;
            try { v = evaluateExpression(expr, geo.params, k); } catch { /* the row shows the error */ }
            if (Number.isFinite(v)) {
                text = kind === 'small' ? `= ${formatLength(v * k)}`
                    : kind === 'length' ? `= ${fmt(v)} ${geo.units}` : `= ${fmt(v)}`;
            }
        }
        out.textContent = text;
    });
}

// Source lines of the conductor rows whose plating panel is open, kept across form
// rebuilds.
const openPlatingRows = new Set();

// Source line of the form row that last had the focus: new rectangles go on top of it.
let selectedLine = 0;
// Control to focus once the form is rebuilt: { line, key }, key null for the first field.
let pendingFocus = null;

const AXIS_NAMES = { x: ['x', 'w'], y: ['y', 'h'] };
const isMirrored = fields => fields.mirror !== undefined && fields.mirror.trim() !== '0';

// A small button of a form row. `key` names it for restoring the focus after a rebuild.
function rowButton(text, title, handler, { disabled = false, key = null, cls = '' } = {}) {
    const b = el('button', { class: `secondary-btn custom-row-btn ${cls}`, title, text, type: 'button' });
    b.disabled = disabled;
    if (key) b.dataset.focus = key;
    b.addEventListener('click', handler);
    return b;
}

// st - the parsed statement, geoRect - its evaluated rectangle (null when it has an error)
function rectRow(model, st, geoRect, index, count) {
    const fields = { ...st.fields };
    let kind = st.kind;
    const isDiel = () => kind === 'diel';
    // This row's statement in the current text: field edits do not rebuild the form.
    const mine = m => m.statements.find(s => s.type === 'rect' && s.line === st.line && s.part === st.part);
    const write = (structural = false, focusKey = null) => {
        const current = mine(analyse().model);
        if (!current) return;
        if (structural) pendingFocus = { line: st.line, key: focusKey };
        applyText(replaceStatementInText(getCustomGeometryText(), current, rectStatementText(kind, fields)), structural);
    };
    const setField = key => v => { if (v === '') delete fields[key]; else fields[key] = v; write(); };

    const kindSel = el('select', { class: 'custom-kind' },
        ...Object.entries(KIND_LABELS).map(([k, label]) => el('option', { value: k, text: label })));
    kindSel.value = kind;
    kindSel.addEventListener('change', () => {
        const wasDiel = isDiel();
        kind = kindSel.value;
        if (isDiel() && !wasDiel) {
            for (const k of ['plating', 'sigma', 'rq', 'plating_sigma', 'plating_t', 'plating_rq']) delete fields[k];
            fields.er = fields.er ?? '4.4'; fields.tand = fields.tand ?? '0.02';
        }
        if (!isDiel() && wasDiel) { delete fields.er; delete fields.tand; delete fields.thin; }
        write(true);
    });

    // One axis: its position and size fields.
    const axisCells = axis => el('span', { class: 'custom-axis' },
        ...AXIS_NAMES[axis].map(k => exprInput(k, fields[k], setField(k))));
    const units = model.statements.find(o => o.type === 'units')?.value ?? 'mm';
    const lengthTip = (what, empty) => `${what}. Empty: ${empty}. A bare number is in ${units}, ` +
        'a suffix gives another unit: 1um, 500nm.';

    let extra, below = null;
    if (isDiel()) {
        extra = [exprInput('er', fields.er, setField('er'), { cls: 'narrow', kind: 'number' }),
                 exprInput('tand', fields.tand, setField('tand'), { placeholder: '0', cls: 'narrow', kind: 'number' })];
    } else {
        const on = new Set(fields.plating === 'all' ? FACES : (fields.plating && fields.plating !== 'none' ? fields.plating.split(',') : []));
        // The plating options sit in a panel that opens from a button on the row. The
        // button names the plated faces, so a collapsed row still shows its plating.
        const platingBtn = el('button', { class: 'secondary-btn custom-plating-toggle', type: 'button' });
        let platingOpen = openPlatingRows.has(st.line);
        const showPlating = () => {
            const list = FACES.filter(f => on.has(f));
            // Initials keep the row on one line: T S B for top, sides, bottom.
            platingBtn.textContent = `${platingOpen ? '▾' : '▸'} plating${list.length ? ' ' + list.map(f => f[0].toUpperCase()).join('') : ''}`;
            platingBtn.classList.toggle('active', list.length > 0);
            platingBtn.title = list.length ? `Plated faces: ${list.join(', ')}. Open for the faces and the plating material`
                : 'No plating. Open to plate faces of this conductor';
            platingPanel.style.display = platingOpen ? '' : 'none';
        };
        platingBtn.addEventListener('click', () => {
            platingOpen = !platingOpen;
            if (platingOpen) openPlatingRows.add(st.line); else openPlatingRows.delete(st.line);
            showPlating();
        });
        const boxes = FACES.map(face => {
            const cb = el('input', { type: 'checkbox' });
            cb.checked = on.has(face);
            cb.addEventListener('change', () => {
                if (cb.checked) on.add(face); else on.delete(face);
                const list = FACES.filter(f => on.has(f));
                if (list.length) fields.plating = list.join(','); else delete fields.plating;
                showPlating();
                // First plated face with no material anywhere: start from a typical one.
                const needsMaterial = list.length && !fields.plating_sigma && !fields.plating_t
                    && !model.statements.some(s2 => s2.type === 'plating');
                if (needsMaterial) { fields.plating_sigma = '1e7'; fields.plating_t = '4um'; }
                write(needsMaterial);
            });
            return el('label', { class: 'custom-face' }, cb, face);
        });
        // Conductivity, roughness and plating material of this conductor. Empty fields
        // fall back to the Conductivity and Surface Roughness options, and to the
        // plating statement.
        const opt = { cls: 'narrow', placeholder: 'default' };
        const platingPanel = el('div', { class: 'custom-plating-panel' },
            el('div', { class: 'custom-cell custom-plating' }, el('span', { class: 'custom-cell-label', text: 'faces' }),
                el('div', { class: 'custom-faces' }, ...boxes)),
            el('span', { class: 'custom-plating-material' },
                exprInput('σ', fields.plating_sigma, setField('plating_sigma'),
                    { ...opt, kind: 'number', title: 'Plating conductivity in S/m. Empty: the plating statement.' }),
                exprInput('t', fields.plating_t, setField('plating_t'),
                    { ...opt, kind: 'small', title: lengthTip('Plating thickness', 'the plating statement') }),
                exprInput('rq', fields.plating_rq, setField('plating_rq'),
                    { ...opt, kind: 'small', title: lengthTip('Plating surface roughness (rms)', 'the plating statement') })));
        showPlating();
        extra = [exprInput('σ', fields.sigma, setField('sigma'),
                { ...opt, kind: 'number', title: 'Conductivity in S/m. Empty: the Conductivity option.' }),
            exprInput('rq', fields.rq, setField('rq'),
                { ...opt, kind: 'small', title: lengthTip('Surface roughness (rms)', 'the Surface Roughness option') }),
            platingBtn];
        below = platingPanel;
    }

    const mirrored = isMirrored(fields);
    const image = { 'sig+': 'a Signal (−)', 'sig-': 'a Signal (+)', gnd: 'a ground', diel: 'a dielectric' }[st.kind];
    const mirrorBtn = rowButton('⇋ mirror',
        (mirrored ? 'Mirrored about x=0. Click to remove the image.\n' : 'Mirror about x=0.\n') +
        `Adds ${image} image on the other side of x=0 (drawn dashed in the preview). ` +
        'A rectangle that touches or crosses x=0 becomes one rectangle symmetric about x=0 instead.',
        () => { if (mirrored) delete fields.mirror; else fields.mirror = '1'; write(true, 'mirror'); },
        { key: 'mirror', cls: 'custom-mirror' });
    mirrorBtn.classList.toggle('active', mirrored);

    const current = () => analyse().model;
    const actions = el('div', { class: 'custom-row-actions' },
        rowButton('↑', 'Move up (a later dielectric covers an earlier one)', () => { const m = current(); applyText(moveRectInText(getCustomGeometryText(), m, mine(m), -1), true); }, { disabled: index === 0 }),
        rowButton('↓', 'Move down', () => { const m = current(); applyText(moveRectInText(getCustomGeometryText(), m, mine(m), 1), true); }, { disabled: index === count - 1 }),
        rowButton('⧉', 'Duplicate', () => {
            pendingFocus = { line: st.line + 1, key: null };
            applyText(insertLineInText(getCustomGeometryText(), st.line, rectStatementText(kind, fields)), true);
        }),
        rowButton('✕', 'Delete', () => { const m = current(); applyText(replaceStatementInText(getCustomGeometryText(), mine(m), null), true); }));

    const row = el('div', { class: `custom-rect-row kind-${st.kind.replace('+', 'p').replace('-', 'n')}`, 'data-line': st.line },
        kindSel, axisCells('x'), axisCells('y'), mirrorBtn, ...extra, actions, ...(below ? [below] : []));
    if (geoRect) {
        row.title = `x ${fmt(geoRect.x.min / geoRect.scale)} … ${fmt(geoRect.x.max / geoRect.scale)}, ` +
                    `y ${fmt(geoRect.y.min / geoRect.scale)} … ${fmt(geoRect.y.max / geoRect.scale)} ${geoRect.units}` +
                    (mirrored ? ', mirrored about x=0' : '');
    }
    // Focusing anything in the row outlines its rectangle in the preview. The row takes
    // the focus itself when clicked outside its fields.
    row.tabIndex = -1;
    row.addEventListener('focusin', () => { selectedLine = st.line; window.customHighlightLine = st.line; redrawSoon(); });
    row.addEventListener('focusout', (e) => {
        if (row.contains(e.relatedTarget)) return;
        window.customHighlightLine = 0; redrawSoon();
    });
    return row;
}

const niceNumber = v => fmt(parseFloat(v.toPrecision(4)));

// Starting fields of a new rectangle, placed from the rectangle whose row had the focus
// last. A dielectric spans the domain on top of it (or of the stack) and takes the
// height of the dielectric below. A conductor goes right of a selected conductor at
// its height. Otherwise it goes on top of the selected dielectric or the stack, right
// of the conductors already at that height, or centred on x=0. The gap is the parameter
// g or gap for a ground, s or gap for a signal, else the new conductor's width.
function newRectFields(kind, model, geo) {
    const k = LENGTH_UNITS[geo.units];
    const drawn = model.statements.filter(s => s.type === 'rect')
        .map(st => ({ st, r: geo.rects.find(r => r.line === st.line && !r.image) }))
        .filter(o => o.r);
    const stacked = drawn.filter(o => Number.isFinite(o.r.y.max));
    const base = stacked.find(o => o.st.line === selectedLine)
        ?? stacked.reduce((a, o) => (!a || o.r.y.max >= a.r.y.max ? o : a), null);
    const y = base ? axisEdges(base.st.fields, 'y', base.r.y.flipped).hi : '0';
    const yV = base ? base.r.y.max / k : 0;
    const ys = stacked.flatMap(o => [o.r.y.min, o.r.y.max]).filter(Number.isFinite);
    const stack = ys.length ? (Math.max(...ys) - Math.min(...ys)) / k : 0;
    const params = new Set(model.statements.filter(s => s.type === 'param').map(s => s.name));
    // The size expression of the first rectangle in `list` that has a finite positive one.
    const sizeOf = (list, key) => {
        const axis = key === 'w' ? 'x' : 'y';
        const o = list.find(q => !q.r[axis].flipped && Number.isFinite(q.r[axis].size)
            && Number.isFinite(q.r[axis].min));
        return o ? o.st.fields[key] : null;
    };
    const value = expr => { try { return evaluateExpression(expr, geo.params, k); } catch { return NaN; } };

    if (kind === 'diel') {
        const diels = drawn.filter(o => o.st.kind === 'diel');
        const h = (base && base.st.kind === 'diel' ? sizeOf([base], 'h') : null) ?? sizeOf(diels, 'h')
            ?? (stack > 0 ? niceNumber(stack / 2) : '1');
        return { x: '-inf', w: 'inf', y, h, er: '4.4', tand: '0.02' };
    }
    const conds = drawn.filter(o => o.st.kind !== 'diel');
    const sameKind = conds.filter(o => o.st.kind === kind);
    const signals = conds.filter(o => o.st.kind !== 'gnd');
    const h = sizeOf(sameKind, 'h') ?? sizeOf(conds, 'h') ?? (params.has('t') ? 't' : '35um');
    const w = sizeOf(kind === 'gnd' ? sameKind : signals, 'w') ?? sizeOf(conds, 'w')
        ?? (params.has('w') ? 'w' : niceNumber(stack > 0 ? 2 * stack : 1));
    const gap = (kind === 'gnd' ? ['g', 'gap'] : ['s', 'gap']).find(n => params.has(n)) ?? w;
    const selected = drawn.find(o => o.st.line === selectedLine && o.st.kind !== 'diel' && Number.isFinite(o.r.x.max)
        && Number.isFinite(o.r.y.min) && Number.isFinite(o.r.y.max));
    if (selected) {
        const neg = selected.r.y.flipped;
        const edges = axisEdges(selected.st.fields, 'y', neg);
        const hSel = neg ? null : selected.st.fields.h;
        return { x: addExpr(axisEdges(selected.st.fields, 'x', selected.r.x.flipped).hi, gap), w,
                 y: edges.lo, h: hSel ?? niceNumber((selected.r.y.max - selected.r.y.min) / k) };
    }
    const hV = value(h), wV = value(w);
    // Conductors at the new one's height, the rightmost of them with a finite right edge.
    const beside = conds.filter(o => o.r.y.min < (yV + hV) * k && o.r.y.max > yV * k && Number.isFinite(o.r.x.max))
        .reduce((a, o) => (!a || o.r.x.max > a.r.x.max ? o : a), null);
    let x;
    if (beside) x = addExpr(axisEdges(beside.st.fields, 'x', beside.r.x.flipped).hi, gap);
    else if (Number.isFinite(wV) && /^[\d.eE+-]+$/.test(w)) x = niceNumber(-wV / 2);
    else x = /^[A-Za-z_][A-Za-z_0-9]*$/.test(w) ? `-${w}/2` : `-(${w})/2`;
    return { x, w, y, h };
}

function renderForm(model, geo) {
    const form = $('custom-form');
    if (!form) return;
    form.innerHTML = '';
    const text = () => getCustomGeometryText();

    // Rectangles.
    const rects = model.statements.filter(s => s.type === 'rect');
    const byLine = new Map(geo.rects.filter(r => !r.image)
        .map(r => [r.line, { ...r, scale: LENGTH_UNITS[geo.units], units: geo.units }]));
    const adders = Object.entries(KIND_LABELS).map(([kind, label]) => {
        const b = el('button', { class: 'secondary-btn', text: `+ ${label}`,
            title: kind === 'diel' ? 'Adds a layer on top of the rectangle selected last, or on top of the stack'
                : 'Adds a conductor beside the conductor selected last, or on top of the selected dielectric or the stack' });
        b.addEventListener('click', () => {
            // Field edits do not rebuild the form, so the model is read afresh.
            const { model: m, geo: g } = analyse();
            const drawn = m.statements.filter(s => s.type === 'rect');
            const after = drawn.length ? Math.max(...drawn.map(r => r.line)) : 1e9;
            const newText = insertLineInText(text(), after, rectStatementText(kind, newRectFields(kind, m, g)));
            const added = analyse(newText).model.statements.filter(s => s.type === 'rect');
            pendingFocus = { line: Math.max(...added.map(s => s.line)), key: null };
            applyText(newText, true);
        });
        return b;
    });
    form.append(el('div', { class: 'custom-form-section' },
        el('div', { class: 'custom-form-title', text: 'Rectangles' },
            el('span', { class: 'custom-hint', text: '  Fields take expressions of the sidebar parameters: w/2, h1+h2, 35um. -inf / inf runs an edge to the boundary, a negative size flips the rectangle to the other side of its position. ⇋ mirror adds the mirror image about x=0. A later dielectric covers an earlier one.' })),
        el('div', { class: 'custom-rect-list' }, ...rects.map((r, i) => rectRow(model, r, byLine.get(r.line), i, rects.length))),
        el('div', { class: 'custom-adders' }, ...adders)));

    if (pendingFocus) {
        const { line, key } = pendingFocus;
        pendingFocus = null;
        const row = form.querySelector(`.custom-rect-row[data-line="${line}"]`);
        const target = row && ((key && row.querySelector(`[data-focus="${key}"]`)) || row.querySelector('input.custom-expr'));
        if (target) { target.focus(); if (target.select) target.select(); }
    }
}

function renderResolvedDomain(solver, geo) {
    const out = $('custom-domain-resolved');
    if (!out) return;
    if (!solver) { out.textContent = ''; return; }
    const u = solver.user_domain, k = LENGTH_UNITS[geo.units], x0 = solver.x_shift || 0;
    out.textContent = `Solved: x ${fmt((u.x_min + x0) / k)} … ${fmt((u.x_max + x0) / k)}, y ${fmt(u.y_min / k)} … ${fmt(u.y_max / k)} ${geo.units}`;
}

// --- Refresh ------------------------------------------------------------------------

// Source line under the text cursor, for highlighting its rectangle in the preview.
function cursorLine() {
    const t = $('custom_geom_text');
    if (!t || document.activeElement !== t) return 0;
    return t.value.slice(0, t.selectionStart).split('\n').length;
}

const formVisible = () => $('custom-form') && $('custom-form').style.display !== 'none';

// source: 'text', 'sidebar' (a field edit), 'sidebar-structure' (a parameter added, removed
// or renamed), 'form' (a field edit, the form keeps its inputs and focus),
// 'form-structure' (rows added, removed, moved or retyped), 'history' (undo, redo) or 'load'.
// notify = false leaves the redraw to the caller (the type-change handler draws once
// itself, with a zoom reset).
// The geometry text survives a page reload in the session storage (browsers do not
// all restore a hidden textarea).
const SESSION_KEY = 'custom_geom_text';
function saveSession(text) { try { sessionStorage.setItem(SESSION_KEY, text); } catch { /* unavailable */ } }
function loadSession() { try { return sessionStorage.getItem(SESSION_KEY) || ''; } catch { return ''; } }

// Undo history of the geometry text. Every control writes the text, so one stack covers
// both views and the sidebar. Edits of a field that follow each other within a second
// are one step.
const HISTORY_MAX = 200;
const FIELD_SOURCES = new Set(['text', 'form', 'sidebar']);
const history = { undo: [], redo: [], last: null, lastSource: '', lastTime: 0 };

function recordHistory(text, source) {
    if (history.last === null || text === history.last) { history.last = text; return; }
    if (source !== 'history') {
        const now = Date.now();
        const merge = FIELD_SOURCES.has(source) && source === history.lastSource && now - history.lastTime < 1000;
        if (!merge) {
            history.undo.push(history.last);
            if (history.undo.length > HISTORY_MAX) history.undo.shift();
        }
        history.redo.length = 0;
        history.lastSource = source;
        history.lastTime = now;
    }
    history.last = text;
}

function updateHistoryButtons() {
    if ($('btn-custom-undo')) $('btn-custom-undo').disabled = !history.undo.length;
    if ($('btn-custom-redo')) $('btn-custom-redo').disabled = !history.redo.length;
}

function stepHistory(from, to) {
    // A pending field edit becomes the state that is stepped away from.
    if (debounceTimer) { clearTimeout(debounceTimer); debounceTimer = null; refresh('form-structure'); }
    if (!from.length) return;
    to.push(getCustomGeometryText());
    // A focused field keeps what was typed in it, so let go of it first.
    // The focus moves to the panel it was in, so the next shortcut still arrives there.
    const active = document.activeElement;
    if (active && active !== $('custom_geom_text')) {
        const panel = active.closest && active.closest('#custom-editor, #custom-params');
        if (panel) panel.focus({ preventScroll: true }); else if (active.blur) active.blur();
    }
    $('custom_geom_text').value = from.pop();
    history.lastSource = '';
    lastParamSignature = null;
    refresh('history');
}

// Maps the line numbers of `before` to those of `after` for an edit that changed one
// block of lines: lines above and below it keep their text, a block of equal length
// keeps its line numbers (an edit in place) or swaps its two end lines (a row moved
// past comments), and lines removed with the block map to 0.
function lineMap(before, after) {
    const a = before.split('\n'), b = after.split('\n');
    let top = 0;
    while (top < a.length && top < b.length && a[top] === b[top]) top++;
    let bottom = 0;
    while (bottom < a.length - top && bottom < b.length - top
           && a[a.length - 1 - bottom] === b[b.length - 1 - bottom]) bottom++;
    const shift = b.length - a.length, endA = a.length - bottom;
    const first = top + 1, last = endA;
    const swapped = shift === 0 && last > first && a[first - 1] === b[last - 1] && a[last - 1] === b[first - 1];
    return line => {
        if (line <= top) return line;
        if (line > endA) return line + shift;
        if (shift !== 0) return 0;
        if (swapped && line === first) return last;
        if (swapped && line === last) return first;
        return line;
    };
}

// Row state kept by source line follows the text through an edit.
function remapRowState(before, after, source) {
    if (source === 'load') { openPlatingRows.clear(); selectedLine = 0; return; }
    if (before === null || before === after) return;
    const map = lineMap(before, after);
    const open = [...openPlatingRows].map(map).filter(l => l > 0);
    openPlatingRows.clear();
    open.forEach(l => openPlatingRows.add(l));
    selectedLine = map(selectedLine);
}

function refresh(source = 'text', notify = true) {
    const text = getCustomGeometryText();
    remapRowState(history.last, text, source);
    saveSession(text);
    recordHistory(text, source);
    updateHistoryButtons();
    const { model, geo } = analyse(text);
    const result = validateCustomGeometry(text);
    renderSidebarParams(model, geo);
    if (source !== 'sidebar') { renderBounds(geo); renderDomain(model); }
    if (formVisible() && source !== 'form') renderForm(model, geo);
    renderResolvedDomain(result.solver, geo);
    if (formVisible()) updateFormValues(geo);
    renderMessages(result);
    if (source === 'text') window.customHighlightLine = cursorLine();
    // A newly loaded geometry (template, file) gets a fresh view, an edit keeps the zoom.
    if (notify) onChange(source === 'load');
}

function scheduleChange(source = 'text') {
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => { debounceTimer = null; refresh(source); }, 250);
}

function showView(view) {
    const form = view === 'form';
    $('custom-form').style.display = form ? 'block' : 'none';
    $('custom_geom_text').style.display = form ? 'none' : 'block';
    $('custom-text-pane').style.display = form ? 'none' : 'flex';
    if (!form) renderLineNumbers();
    $('btn-custom-format').style.display = form ? 'none' : '';
    $('btn-custom-view-form').classList.toggle('active', form);
    $('btn-custom-view-text').classList.toggle('active', !form);
    // The error box lists different errors in the two views.
    refresh(form ? 'form-structure' : 'text', false);
}

// notify = false when the caller redraws itself (settings restore, conversion).
export function setCustomGeometryText(text, notify = false) {
    const t = $('custom_geom_text');
    if (!t) return;
    t.value = text;
    lastParamSignature = null;
    if ($('tl_type').value === 'custom') refresh('load', notify);
}

// Called when the custom type is selected: fills an empty editor with the first
// template and brings the controls in line with the text.
export function activateCustomGeometry() {
    const t = $('custom_geom_text');
    if (!t) return;
    // dataset.loaded marks an untouched template, which another template may replace
    // without asking.
    if (!t.value.trim()) t.value = loadSession();
    if (!t.value.trim()) { t.value = Object.values(CUSTOM_TEMPLATES)[0]; t.dataset.loaded = t.value; }
    lastParamSignature = null;
    refresh('load', false);
}

// Changes the declared unit. Without any length written in the text there is nothing
// to decide. Otherwise the numbers are either converted, so the geometry keeps its
// size, or kept and read in the new unit.
async function changeUnits(to) {
    const text = getCustomGeometryText();
    const { model, geo } = analyse(text);
    const from = geo.units;
    if (to === from) return;
    const reinterpret = () => setText(setStatementInText(text, 'units', `units ${to}`));
    const hasLengths = model.statements.some(s => s.type === 'param' || s.type === 'rect' || s.type === 'plating'
        || (s.type === 'domain' && s.values.some(v => v !== 'auto')));
    if (!hasLengths) { reinterpret(); return; }
    const choice = await askChoice('Change the length unit',
        `Bare numbers in the geometry are in ${from}. They can be converted to ${to}, so the geometry keeps ` +
        `its size, or kept as they are and read in ${to}, which scales the geometry.\n` +
        'Numbers with their own unit (35um) stay as written.',
        [{ label: `Convert to ${to}`, value: 'convert', primary: true },
         { label: `Keep the numbers`, value: 'reinterpret' },
         { label: 'Cancel', value: null }]);
    // The text may have changed while the dialog was open.
    if (choice === null || getCustomGeometryText() !== text) { $('custom-units').value = analyse().geo.units; return; }
    if (choice === 'reinterpret') { reinterpret(); return; }
    const converted = changeUnitsInText(text, to);
    if (converted === null) {
        window.alert(geo.errors.length
            ? 'The geometry has errors. Fix them before converting the unit.'
            : 'Some expressions could not be converted consistently. The unit was not changed.');
        $('custom-units').value = from;
        return;
    }
    setText(converted);
}

export function initCustomGeometryEditor({ onGeometryChange, log }) {
    onChange = onGeometryChange;
    const text = $('custom_geom_text');
    if (!text) return;

    text.addEventListener('input', () => { renderLineNumbers(); scheduleChange('text'); });
    text.addEventListener('scroll', () => { $('custom-line-numbers').scrollTop = text.scrollTop; });
    // Tab inserts spaces instead of leaving the editor.
    text.addEventListener('keydown', (e) => {
        if (e.key !== 'Tab' || e.shiftKey) return;
        e.preventDefault();
        text.setRangeText('  ', text.selectionStart, text.selectionEnd, 'end');
        scheduleChange('text');
    });
    const moveHighlight = () => {
        const line = cursorLine();
        if (line === window.customHighlightLine) return;
        window.customHighlightLine = line;
        redrawSoon();
    };
    for (const ev of ['click', 'keyup', 'focus']) text.addEventListener(ev, moveHighlight);
    text.addEventListener('blur', () => { window.customHighlightLine = 0; redrawSoon(); });

    $('btn-custom-undo').addEventListener('click', () => stepHistory(history.undo, history.redo));
    $('btn-custom-redo').addEventListener('click', () => stepHistory(history.redo, history.undo));
    // Ctrl+Z / Ctrl+Shift+Z / Ctrl+Y in the editor and the sidebar block. The browser's
    // own undo of a field does not see the edits the controls make to the text.
    for (const id of ['custom-editor', 'custom-params']) {
        $(id).tabIndex = -1;
        $(id).addEventListener('keydown', (e) => {
            if (!(e.ctrlKey || e.metaKey) || e.altKey) return;
            const key = e.key.toLowerCase();
            if (key !== 'z' && key !== 'y') return;
            e.preventDefault();
            if (key === 'y' || e.shiftKey) stepHistory(history.redo, history.undo);
            else stepHistory(history.undo, history.redo);
        });
    }
    $('custom-error-count').addEventListener('click', showFirstError);
    const FORMAT_KEY = 'custom_format_open';
    const setFormat = (open) => {
        $('custom-format-ref').style.display = open ? 'block' : 'none';
        $('btn-custom-format').classList.toggle('active', open);
    };
    let formatOpen = false;
    try { formatOpen = localStorage.getItem(FORMAT_KEY) === '1'; } catch { /* unavailable */ }
    setFormat(formatOpen);
    $('btn-custom-format').addEventListener('click', () => {
        formatOpen = !formatOpen;
        setFormat(formatOpen);
        try { localStorage.setItem(FORMAT_KEY, formatOpen ? '1' : '0'); } catch { /* unavailable */ }
    });

    $('btn-custom-view-form').addEventListener('click', () => showView('form'));
    $('btn-custom-view-text').addEventListener('click', () => showView('text'));

    const sel = $('custom-template');
    for (const name of Object.keys(CUSTOM_TEMPLATES)) sel.append(el('option', { value: name, text: name }));
    sel.addEventListener('change', () => {
        if (!sel.value) return;
        const replace = replaceable() || window.confirm('Replace the current geometry with the template?');
        if (replace) { setCustomGeometryText(CUSTOM_TEMPLATES[sel.value], true); text.dataset.loaded = text.value; }
        sel.value = '';
    });

    const unitSel = $('custom-units');
    for (const u of Object.keys(LENGTH_UNITS).filter(u => u !== 'µm')) unitSel.append(el('option', { value: u, text: u }));
    unitSel.addEventListener('change', () => changeUnits(unitSel.value));
    $('btn-custom-add-param').addEventListener('click', () => addParameter());
    DOMAIN_KEYS.forEach(key => {
        $(`custom_domain_${key}`).addEventListener('input', () => {
            // Values in a statement cannot hold spaces, an empty field is automatic.
            const values = DOMAIN_KEYS.map(k => $(`custom_domain_${k}`).value.replace(/\s+/g, '') || 'auto');
            const stmt = values.every(v => v === 'auto') ? 'domain auto' : `domain ${values.join(' ')}`;
            text.value = setStatementInText(text.value, 'domain', stmt);
            scheduleChange('sidebar');
        });
    });

    WALLS.forEach(wall => {
        $(`custom_bound_${wall}`).addEventListener('change', () => {
            const values = WALLS.map(w => $(`custom_bound_${w}`).value);
            text.value = setStatementInText(text.value, 'bounds', `bounds ${values.join(' ')}`);
            scheduleChange('sidebar');
        });
    });

    $('btn-custom-copy').addEventListener('click', () => {
        navigator.clipboard.writeText(text.value).then(() => {
            const b = $('btn-custom-copy'), old = b.textContent;
            b.textContent = 'Copied!';
            setTimeout(() => { b.textContent = old; }, 1500);
        }).catch(() => log('Could not copy to the clipboard.'));
    });
    $('btn-custom-save').addEventListener('click', () => {
        const url = URL.createObjectURL(new Blob([text.value], { type: 'text/plain' }));
        const a = el('a', { href: url, download: 'geometry.txt' });
        a.click();
        URL.revokeObjectURL(url);
    });
    $('btn-custom-load').addEventListener('click', () => $('custom-file-input').click());
    $('custom-file-input').addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (!file) return;
        if (!replaceable() && !window.confirm(`Replace the current geometry with ${file.name}?`)) return;
        file.text().then(t => { setCustomGeometryText(t, true); log(`Loaded geometry from ${file.name}`); });
        e.target.value = '';
    });
    $('btn-custom-collapse').addEventListener('click', () => {
        const collapsed = $('custom-editor').classList.toggle('collapsed');
        $('btn-custom-collapse').textContent = collapsed ? 'Show editor' : 'Hide editor';
        window.dispatchEvent(new Event('resize'));
    });
}
