// Editor for the custom geometry type. The sidebar holds the scalar inputs (units,
// parameters, boundaries, solved region), the panel under the geometry preview holds the
// rectangles: a form view with one row per rectangle and a text view of the whole
// geometry, plus templates and file / clipboard transfer. The text is the single source
// of truth. Every control rewrites one statement in it and leaves the rest as typed.
import { parseGeometryText, evaluateGeometry, setParamInText, renameParamInText, setStatementInText, replaceStatementInText,
         insertLineInText, moveRectInText, rectStatementText, formatErrors, evaluateExpression, LENGTH_UNITS,
         axisEdges, addExpr, formatLength, plausibilityWarnings, isPlainNumber, isLengthLiteral, isReservedName,
         changeUnitsInText, WALLS, ROUND_KEYS, PLATING_FACES as FACES, SHAPE_NAMES, CONDUCTOR_KEYS } from './custom_geometry_text.js';
import { CustomGeometrySolver } from './custom_geometry.js';
import { CONDUCTOR_COLOR, dielectricRGB, overAir, rgbToHex, hexToRGB } from './body_colors.js';

const CUSTOM_TEMPLATES = {
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
    'Integrated circuit microstrip': `# On-chip microstrip: copper top metal on SiO2 over an aluminium ground plane on
# a conductive silicon substrate. Over it a conformal SiO2 layer (tox on top, tox_s on the
# trace sides) and a conformal SiN passivation tn thick.
# hox is the trace bottom height from the bottom of the ground plane.
units um
w = 11.5; t = 3; wgnd = 46; tg = 0.49
hox = 7.66; hbox = 2.6
tox = 1; tox_s = 0.6; tn = 0.4
bounds open open open open
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-inf           y=-hbox    w=inf             h=-inf          er=11.9  sigma=2  color=#5e5c64   # Si, 50 ohm*cm
diel  x=-inf           y=-hbox    w=inf             h=hbox+hox+tox  er=4.1  tand=0.001  color=#62a0ea   # SiO2
diel  x=-inf           y=hox+tox  w=inf             h=tn            er=6.6  tand=0.001  thin=1  color=#99c1f1   # SiN
diel  x=-w/2-tox_s-tn  y=hox+tox  w=w+2*(tox_s+tn)  h=t+tn          er=6.6  tand=0.001  color=#99c1f1
diel  x=-w/2-tox_s     y=hox      w=w+2*tox_s       h=t+tox         er=4.1  tand=0.001  color=#62a0ea
gnd   x=-wgnd/2        y=0        w=wgnd            h=tg            sigma=3.5e7  color=#c0bfbc   # Al
sig+  x=-w/2           y=hox      w=w               h=t
`,
    'CPW over silicon': `# On-chip coplanar waveguide: copper metal embedded in SiO2 on a conductive silicon
# substrate, air above and below. hb is the SiO2 under the metal, the gap sets about
# 50 ohm above 10 GHz. Rounded metal corners, full-wave solver only.
units um
w = 5; g = 2.8; t = 1; wgnd = 20; rc = 0.5
hsi = 300; hox = 10; hb = 3
bounds open open open open
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-inf          y=-hsi  w=inf   h=hsi  er=11.9  sigma=2  color=#5e5c64   # Si, 50 ohm*cm
diel  x=-inf          y=0     w=inf   h=hox  er=4.1  tand=0.001  color=#62a0ea   # SiO2
gnd   x=-w/2-g-wgnd   y=hb    w=wgnd  h=t    radius=rc
gnd   x=w/2+g         y=hb    w=wgnd  h=t    radius=rc
sig+  x=-w/2          y=hb    w=w     h=t    radius=rc
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
    'Microstrip with etched trace': `# Microstrip whose trace is a trapezoid: the sides lean in by the etch angle
# (degrees from the vertical). Full-wave solver only.
units mm
w = 0.3; t = 0.035; h = 0.2; etch = 30
bounds open open open gnd
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-inf  y=0  w=inf  h=h  er=4.3  tand=0.02
sig+  trap  x=-w/2  y=h  w=w  h=t  angle=etch
`,
    'Differential microstrip with etched traces': `# Differential microstrip, trapezoidal traces. Full-wave solver only.
units mm
w = 0.3; s = 0.2; t = 0.035; h = 0.2; etch = 30
bounds open open open gnd
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-inf  y=0  w=inf  h=h  er=4.3  tand=0.02
sig+  trap  x=s/2  y=h  w=w  h=t  angle=etch  mirror=1
`,
    'Coaxial line': `# Coaxial line from n-gons: dielectric, centre conductor and a shield ring.
# r is the vertex radius. Full-wave solver only.
units mm
d = 0.92; D = 2.95; t_sh = 0.15
bounds open open open open
domain -1.1*(D/2+t_sh) 1.1*(D/2+t_sh) -1.1*(D/2+t_sh) 1.1*(D/2+t_sh)

diel  ngon  x=0  y=0  r=D/2  n=128  er=2.1  tand=0.0002
sig+  ngon  x=0  y=0  r=d/2  n=64
gnd   ngon  x=0  y=0  r=D/2+t_sh  r_in=D/2  n=128
`,
    'Microstrip with rounded trace corners': `# Microstrip whose trace has rounded top corners (radius) and sharp bottom
# corners where it sits on the substrate (radius_bottom=0). Full-wave solver only.
units mm
w = 0.3; t = 0.035; h = 0.2; rc = 10um
bounds open open open gnd
domain auto   # sized from the conductors, or: domain x1 x2 y1 y2

diel  x=-inf  y=0  w=inf  h=h  er=4.3  tand=0.02
sig+  x=-w/2  y=h  w=w  h=t  radius=rc  radius_bottom=0
`,
    'Elliptical coaxial line': `# Round centre conductor in an elliptical shield. a and b are the semi-axes of
# the dielectric, the shield is an elliptical ring around it. Full-wave solver only.
units mm
d = 0.9; a = 2.2; b = 1.4; t_sh = 0.15
bounds open open open open
domain -1.1*(a+t_sh) 1.1*(a+t_sh) -1.1*(b+t_sh) 1.1*(b+t_sh)

diel  ellipse  x=0  y=0  rx=a  ry=b  n=128  er=2.1  tand=0.0002
sig+  ngon  x=0  y=0  r=d/2  n=64
gnd   ellipse  x=0  y=0  rx=a+t_sh  ry=b+t_sh  rx_in=a  ry_in=b  n=128
`,
    'Twinax cable': `# Twinax: two insulated wires side by side, wrapped in a stadium-shaped shield
# that touches the insulation (a rectangle with fully rounded ends, radius = half its
# height, and a wall). Full-wave solver only.
units mm
d = 0.4; di = 1.2; t_sh = 0.03
s = di   # wire spacing: the insulations touch
hw = s/2+di/2; hh = di/2   # half width and half height inside the shield
bounds open open open open
domain -1.2*(hw+t_sh) 1.2*(hw+t_sh) -1.5*(hh+t_sh) 1.5*(hh+t_sh)

diel  ngon  x=s/2  y=0  r=di/2  n=64  er=2.1  tand=0.0003  mirror=1
sig+  ngon  x=s/2  y=0  r=d/2  n=32  mirror=1
gnd   x=-hw-t_sh  w=2*(hw+t_sh)  y=-hh-t_sh  h=2*(hh+t_sh)  radius=hh+t_sh  wall=t_sh
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

const KIND_LABELS = { 'sig+': 'Signal (+)', 'sig-': 'Signal (−)', 'gnd': 'Ground', 'diel': 'Dielectric' };
const CORNER_KEYS = ['radius', 'radius_bottom', 'wall'];
// Side angle a rectangle turned into a trapezoid starts with, degrees from the vertical.
const DEFAULT_ANGLE = '20';
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
let onHighlight = () => {};
let debounceTimer = null, pendingSource = 'text';
// Highlight changes arrive in pairs (focus leaves one row and enters the next), so the
// redraw they ask for is coalesced to one per frame. Only the plot is redrawn, the
// geometry and any solution stay.
let highlightFrame = 0;
function redrawSoon() {
    if (highlightFrame) return;
    highlightFrame = requestAnimationFrame(() => { highlightFrame = 0; onHighlight(); });
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
// Flags a field whose value holds ';' or '#', which the text cannot take in a value.
function markStatementBreak(input, bad) {
    input.classList.toggle('invalid', bad);
    input.title = bad ? 'A value cannot contain ; or #.' : (input.dataset.title || '');
}

function exprInput(label, value, commit, { placeholder = '', cls = '', title = '', kind = 'length' } = {}) {
    const input = el('input', { type: 'text', class: `custom-expr ${cls}`, value: value ?? '', placeholder,
        title: title || undefined, spellcheck: 'false', autocomplete: 'off', autocapitalize: 'off' });
    input.dataset.value = kind;
    input.dataset.title = title;
    // Values in a statement cannot hold spaces ("1 um", "w + 2"): the text gets the
    // value without them, the field keeps what was typed. ';' and '#' would end the
    // statement there, so such a value stays in the field only.
    input.addEventListener('input', () => {
        const bad = /[;#]/.test(input.value);
        markStatementBreak(input, bad);
        if (!bad) commit(input.value.replace(/\s+/g, ''));
    });
    return el('label', { class: 'custom-cell' }, el('span', { class: 'custom-cell-label', text: label }),
        el('span', { class: 'custom-cell-field' }, input, el('span', { class: 'custom-cell-value' })));
}

// Color band on the left of a row: the row's fill color in the plots, or the kind's color
// while it has none. Clicking it opens a popup with the color picker, a hex field and a
// button back to the default. `fallback` is the plot's default color for the row, where
// the picker starts. commit gets '#rrggbb', or null for the default.
function colorBand(value, fallback, commit) {
    // The color set on the row, '#rrggbb' or null. A pick does not rebuild the row, so
    // this follows it.
    let current = hexToRGB(value) ? rgbToHex(hexToRGB(value)) : null;
    const band = el('button', { type: 'button', class: 'custom-color-band' });
    const show = () => {
        band.style.background = current ?? '';
        band.title = current ? `Fill color in the plots, ${current}. Click to change.` : 'Fill color in the plots: the default. Click to change.';
    };
    show();
    band.addEventListener('click', () => {
        const old = band.parentElement?.querySelector('.custom-color-popup');
        if (old) { old.close(); return; }
        const start = current ?? fallback;
        const picker = el('input', { type: 'color', class: 'custom-color', value: start, title: 'Pick a color' });
        const hex = el('input', { type: 'text', class: 'custom-color-hex', value: start, spellcheck: 'false',
            autocomplete: 'off', title: '#rrggbb or #rgb' });
        const set = rgb => { current = rgbToHex(rgb); show(); commit(current); };
        picker.addEventListener('input', () => { hex.value = picker.value; set(hexToRGB(picker.value)); });
        hex.addEventListener('input', () => {
            const rgb = hexToRGB(hex.value);
            hex.classList.toggle('invalid', !rgb);
            if (rgb) { picker.value = rgbToHex(rgb); set(rgb); }
        });
        const popup = el('div', { class: 'custom-color-popup' }, picker, hex,
            rowButton('Default', 'Back to the default color: conductors orange, dielectrics shaded by er',
                () => { close(); current = null; show(); commit(null); }));
        const r = band.getBoundingClientRect();
        popup.style.left = `${r.right + 4}px`;
        popup.style.top = `${r.top}px`;
        const onDown = e => { if (!popup.isConnected || (!popup.contains(e.target) && e.target !== band)) close(true); };
        const onKey = e => { if (e.key === 'Escape') close(); };
        // `away`: a click outside, which deselects the row. A focused element that is
        // removed gives no focusout in every browser, and the row's focusout is what
        // clears its highlight in the plot, so the focus leaves before the popup goes.
        function close(away = false) {
            if (popup.contains(document.activeElement)) {
                if (away) document.activeElement.blur(); else band.focus();
            }
            popup.remove();
            document.removeEventListener('mousedown', onDown, true);
            document.removeEventListener('keydown', onKey, true);
        }
        popup.close = close;
        document.addEventListener('mousedown', onDown, true);
        document.addEventListener('keydown', onKey, true);
        band.after(popup);
        // Keep the popup inside the window.
        const p = popup.getBoundingClientRect();
        if (p.bottom > window.innerHeight) popup.style.top = `${Math.max(0, window.innerHeight - p.height - 4)}px`;
        hex.focus();
        hex.select();
    });
    return band;
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
// Source lines of the rows whose corner panel (radius, wall) is open.
const openCornerRows = new Set();

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

// Button that opens and closes a panel of a form row. The open rows are kept by source
// line in `openSet` across form rebuilds. label() gives the button's { text, title,
// active }. Returns { btn, show }, show() redraws the button and the panel.
function panelToggle(cls, openSet, line, panel, label) {
    const btn = el('button', { class: `secondary-btn ${cls}`, type: 'button' });
    let open = openSet.has(line);
    const show = () => {
        const { text, title, active } = label();
        btn.textContent = `${open ? '▾' : '▸'} ${text}`;
        btn.classList.toggle('active', active);
        btn.title = title;
        panel.style.display = open ? '' : 'none';
    };
    btn.addEventListener('click', () => {
        open = !open;
        if (open) openSet.add(line); else openSet.delete(line);
        show();
    });
    return { btn, show };
}

// st - the parsed statement, geoRect - its evaluated rectangle (null when it has an error)
// Fields of a row whose shape changes from `from` to `to`. Position and size carry over
// through the evaluated bounding box when the two shapes describe them differently.
function shapeFields(fields, from, to, geoRect) {
    const f = { ...fields };
    const box = geoRect && [geoRect.x.min, geoRect.x.max, geoRect.y.min, geoRect.y.max].every(Number.isFinite)
        ? [geoRect.x.min, geoRect.x.max, geoRect.y.min, geoRect.y.max].map(v => v / geoRect.scale) : null;
    const num = v => fmt(parseFloat(v.toPrecision(6)));
    const [x0, x1, y0, y1] = box ?? [-0.5, 0.5, -0.5, 0.5];
    const round = k => !!ROUND_KEYS[k];
    if (round(to) && round(from)) {
        // n-gon <-> ellipse: the radius becomes both semi-axes and back.
        if (to === 'ellipse') {
            if (f.r !== undefined) Object.assign(f, { rx: f.r, ry: f.r });
            if (f.r_in !== undefined) Object.assign(f, { rx_in: f.r_in, ry_in: f.r_in });
            delete f.r; delete f.r_in;
        } else {
            // Equal semi-axes carry over, an elongated ellipse becomes the n-gon that
            // fits inside it.
            f.r = f.rx !== undefined && f.rx === f.ry ? f.rx : num(Math.min(x1 - x0, y1 - y0) / 2);
            if (f.rx_in !== undefined && f.rx_in === f.ry_in) f.r_in = f.rx_in;
            for (const k of ['rx', 'ry', 'rx_in', 'ry_in']) delete f[k];
        }
    } else if (round(to)) {
        for (const k of ['x', 'w', 'y', 'h', 'angle', 'angle2', ...CORNER_KEYS]) delete f[k];
        Object.assign(f, { x: num((x0 + x1) / 2), y: num((y0 + y1) / 2) },
            to === 'ngon' ? { r: num(Math.min(x1 - x0, y1 - y0) / 2) } : { rx: num((x1 - x0) / 2), ry: num((y1 - y0) / 2) },
            { n: '32' });
    } else if (round(from)) {
        for (const k of ROUND_KEYS[from]) delete f[k];
        Object.assign(f, { x: num(x0), w: num(x1 - x0), y: num(y0), h: num(y1 - y0) });
    }
    // A round shape is plated all around or not at all.
    if (round(to) && f.plating && f.plating !== 'none') f.plating = 'all';
    if (to === 'trap' && f.angle === undefined && f.angle2 === undefined) f.angle = DEFAULT_ANGLE;
    if (to !== 'trap') { delete f.angle; delete f.angle2; }
    return f;
}

function rectRow(model, st, geoRect, index, count) {
    let fields = { ...st.fields };
    let kind = st.kind;
    let shape = st.shape || 'rect';
    const isDiel = () => kind === 'diel';
    // This row's statement in the current text: field edits do not rebuild the form.
    const mine = m => m.statements.find(s => s.type === 'rect' && s.line === st.line && s.part === st.part);
    const write = (structural = false, focusKey = null) => {
        const current = mine(analyse().model);
        if (!current) return;
        if (structural) pendingFocus = { line: st.line, key: focusKey };
        applyText(replaceStatementInText(getCustomGeometryText(), current, rectStatementText(kind, fields, shape)), structural);
    };
    const setField = key => v => { if (v === '') delete fields[key]; else fields[key] = v; write(); };

    const kindSel = el('select', { class: 'custom-kind' },
        ...Object.entries(KIND_LABELS).map(([k, label]) => el('option', { value: k, text: label })));
    kindSel.value = kind;
    kindSel.addEventListener('change', () => {
        const wasDiel = isDiel();
        kind = kindSel.value;
        if (isDiel() && !wasDiel) {
            for (const k of CONDUCTOR_KEYS) delete fields[k];
            fields.er = fields.er ?? '4.4'; fields.tand = fields.tand ?? '0.02';
        }
        if (!isDiel() && wasDiel) { delete fields.er; delete fields.tand; delete fields.thin; delete fields.sigma; }
        write(true);
    });

    const shapeSel = el('select', { class: 'custom-shape', title: 'Shape. Anything but a plain rectangle needs the full-wave solver.' },
        ...Object.entries(SHAPE_NAMES).map(([k, name]) => el('option', { value: k, text: name[0].toUpperCase() + name.slice(1) })));
    shapeSel.value = shape;
    shapeSel.addEventListener('change', () => {
        fields = shapeFields(fields, shape, shapeSel.value, geoRect);
        shape = shapeSel.value;
        write(true, 'shape');
    });
    shapeSel.dataset.focus = 'shape';

    // One axis: its position and size fields.
    const axisCells = axis => el('span', { class: 'custom-axis' },
        ...AXIS_NAMES[axis].map(k => exprInput(k, fields[k], setField(k))));
    const units = model.statements.find(o => o.type === 'units')?.value ?? 'mm';
    const lengthTip = (what, empty) => `${what}. Empty: ${empty}. A bare number is in ${units}, ` +
        'a suffix gives another unit: 1um, 500nm.';
    let geometryCells, cornerPanel = null;
    if (shape === 'ngon') {
        geometryCells = [
            el('span', { class: 'custom-axis' },
                exprInput('x', fields.x, setField('x'), { title: 'Centre x' }),
                exprInput('y', fields.y, setField('y'), { title: 'Centre y' })),
            el('span', { class: 'custom-axis' },
                exprInput('r', fields.r, setField('r'), { title: 'Vertex radius: the distance from the centre to each vertex' }),
                exprInput('n', fields.n, setField('n'), { cls: 'narrow', kind: 'number', title: 'Number of vertices, 3 to 1024. One vertex is on top.' }),
                exprInput('r_in', fields.r_in, setField('r_in'), { cls: 'narrow', placeholder: 'solid',
                    title: 'Vertex radius of a concentric hole, which makes a ring (a coax shield). Empty: solid.' })),
        ];
    } else if (shape === 'ellipse') {
        geometryCells = [
            el('span', { class: 'custom-axis' },
                exprInput('x', fields.x, setField('x'), { title: 'Centre x' }),
                exprInput('y', fields.y, setField('y'), { title: 'Centre y' })),
            el('span', { class: 'custom-axis' },
                exprInput('rx', fields.rx, setField('rx'), { title: 'Horizontal semi-axis' }),
                exprInput('ry', fields.ry, setField('ry'), { title: 'Vertical semi-axis' }),
                exprInput('n', fields.n, setField('n'), { cls: 'narrow', kind: 'number',
                    title: 'Number of vertices, 3 to 1024, equally spaced along the outline. One vertex is on top.' })),
            el('span', { class: 'custom-axis' },
                exprInput('rx_in', fields.rx_in, setField('rx_in'), { cls: 'narrow', placeholder: 'solid',
                    title: 'Horizontal semi-axis of a concentric elliptical hole, which makes a ring. Empty: solid.' }),
                exprInput('ry_in', fields.ry_in, setField('ry_in'), { cls: 'narrow', placeholder: 'solid',
                    title: 'Vertical semi-axis of the hole. Empty: solid.' })),
        ];
    } else {
        geometryCells = [axisCells('x'), axisCells('y')];
        if (shape === 'trap') {
            geometryCells.push(el('span', { class: 'custom-axis' },
                exprInput('∠L', fields.angle, setField('angle'), { cls: 'narrow', kind: 'number',
                    title: 'Left side angle in degrees from the vertical. Positive: the face at y+h is narrower than the base at y, negative: wider.' }),
                exprInput('∠R', fields.angle2, setField('angle2'), { cls: 'narrow', kind: 'number', placeholder: '= ∠L',
                    title: 'Right side angle in degrees from the vertical. Empty: the same as the left side.' })));
        }
        // Corner radius and wall sit in a panel that opens from a button, like plating.
        const setCorner = key => v => {
            // A shell has one outside: its plating covers it all.
            if (key === 'wall' && v !== '' && !/^0*\.?0*$/.test(v) && fields.plating && fields.plating !== 'none') fields.plating = 'all';
            setField(key)(v); showCorners();
        };
        cornerPanel = el('div', { class: 'custom-corner-panel' },
            exprInput('radius', fields.radius, setCorner('radius'), { cls: 'narrow', placeholder: '0',
                title: lengthTip('Corner radius, every corner', 'sharp corners') }),
            exprInput('bottom', fields.radius_bottom, setCorner('radius_bottom'), { cls: 'narrow', placeholder: '= radius',
                title: lengthTip('Radius of the two lower corners', 'the radius above. 0 keeps them sharp, as on a trace etched from a foil') }),
            exprInput('wall', fields.wall, setCorner('wall'), { cls: 'narrow', placeholder: 'solid',
                title: lengthTip('Wall thickness: the shape becomes a shell of this thickness, such as a cable shield', 'solid') }));
        const { btn: cornersBtn, show: showCorners } = panelToggle('custom-corners-toggle', openCornerRows, st.line, cornerPanel, () => {
            const set = CORNER_KEYS.filter(k => fields[k] !== undefined && fields[k] !== '');
            return { text: `corners${set.length ? ' ' + set.map(k => ({ radius: 'r', radius_bottom: 'r↓', wall: 'wall' })[k]).join(' ') : ''}`,
                active: set.length > 0,
                title: 'Rounded corners and a wall (a hollow shape). Anything set here needs the full-wave solver.' };
        });
        geometryCells.push(cornersBtn);
        showCorners();
    }

    let extra, below = null;
    if (isDiel()) {
        extra = [exprInput('er', fields.er, setField('er'), { cls: 'narrow', kind: 'number' }),
                 exprInput('tand', fields.tand, setField('tand'), { placeholder: '0', cls: 'narrow', kind: 'number' }),
                 exprInput('σ', fields.sigma, setField('sigma'), { placeholder: '0', cls: 'narrow', kind: 'number',
                     title: 'Conductivity in S/m, for a semiconductor such as a silicon substrate (1/(resistivity in ohm*m)).' })];
    } else {
        const on = new Set(fields.plating === 'all' ? FACES : (fields.plating && fields.plating !== 'none' ? fields.plating.split(',') : []));
        // A round shape or a shell has one outside: its plating is all or nothing. The
        // evaluated shape decides (a wall of 0 is no shell), the typed keys when the row
        // has an error.
        const round = !!ROUND_KEYS[shape] || (geoRect ? !!(geoRect.shape && geoRect.shape.type === 'ring')
            : !!fields.wall || fields.r_in !== undefined || fields.rx_in !== undefined);
        const boxes = (round ? ['all'] : FACES).map(face => {
            const cb = el('input', { type: 'checkbox' });
            cb.checked = round ? on.size === FACES.length : on.has(face);
            cb.addEventListener('change', () => {
                const set = round ? FACES : [face];
                for (const f of set) { if (cb.checked) on.add(f); else on.delete(f); }
                const list = FACES.filter(f => on.has(f));
                if (list.length) fields.plating = round ? 'all' : list.join(','); else delete fields.plating;
                showPlating();
                // First plated face with no material anywhere: start from a typical one.
                const needsMaterial = list.length && !fields.plating_sigma && !fields.plating_t
                    && !model.statements.some(s2 => s2.type === 'plating');
                if (needsMaterial) { fields.plating_sigma = '1e7'; fields.plating_t = '4um'; }
                write(needsMaterial);
            });
            return el('label', { class: 'custom-face' }, cb, round ? 'whole surface' : face);
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
                    { ...opt, kind: 'small', title: lengthTip('Plating surface roughness (rms)', 'the plating statement') }),
                exprInput('rq iface', fields.plating_rq_iface, setField('plating_rq_iface'),
                    { ...opt, kind: 'small', title: lengthTip('Roughness (rms) of the interface between plating and bulk metal',
                        'the plating statement, else the plating roughness') })));
        // The plating options sit in a panel that opens from a button on the row. The
        // button names the plated faces, so a collapsed row still shows its plating.
        const { btn: platingBtn, show: showPlating } = panelToggle('custom-plating-toggle', openPlatingRows, st.line, platingPanel, () => {
            const list = FACES.filter(f => on.has(f));
            // Initials keep the row on one line: T S B for top, sides, bottom.
            return { text: `plating${list.length ? ' ' + (round ? 'all' : list.map(f => f[0].toUpperCase()).join('')) : ''}`,
                active: list.length > 0,
                title: list.length ? `Plated faces: ${list.join(', ')}. Open for the faces and the plating material`
                    : 'No plating. Open to plate faces of this conductor' };
        });
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
        (shape === 'rect' && !CORNER_KEYS.some(k => fields[k]) ? 'A rectangle that touches or crosses x=0 becomes one rectangle symmetric about x=0 instead.'
            : 'A shape on x=0 has to be symmetric about it.'),
        () => { if (mirrored) delete fields.mirror; else fields.mirror = '1'; write(true, 'mirror'); },
        { key: 'mirror', cls: 'custom-mirror' });
    mirrorBtn.classList.toggle('active', mirrored);

    // A structural edit of this row's statement in the current text.
    const edit = fn => () => { const m = analyse().model; applyText(fn(getCustomGeometryText(), m, mine(m)), true); };
    const actions = el('div', { class: 'custom-row-actions' },
        rowButton('↑', 'Move up (a later dielectric covers an earlier one)', edit((t, m, s) => moveRectInText(t, m, s, -1)), { disabled: index === 0 }),
        rowButton('↓', 'Move down', edit((t, m, s) => moveRectInText(t, m, s, 1)), { disabled: index === count - 1 }),
        rowButton('⧉', 'Duplicate', () => {
            pendingFocus = { line: st.line + 1, key: null };
            applyText(insertLineInText(getCustomGeometryText(), st.line, rectStatementText(kind, fields, shape)), true);
        }),
        rowButton('✕', 'Delete', edit((t, m, s) => replaceStatementInText(t, s, null))));

    const band = colorBand(fields.color, isDiel() ? rgbToHex(overAir(dielectricRGB(geoRect ? geoRect.er : parseFloat(fields.er)), 0.8))
        : CONDUCTOR_COLOR, v => { if (v === null) delete fields.color; else fields.color = v; write(v === null); });
    const row = el('div', { class: `custom-rect-row kind-${st.kind.replace('+', 'p').replace('-', 'n')}`, 'data-line': st.line },
        band, kindSel, shapeSel, ...geometryCells, mirrorBtn, ...extra, actions, ...(cornerPanel ? [cornerPanel] : []), ...(below ? [below] : []));
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
    // An n-gon or ellipse has no position and size fields to place a rectangle from.
    const drawn = model.statements.filter(s => s.type === 'rect' && !ROUND_KEYS[s.shape])
        .map(st => ({ st, r: geo.rects.find(r => r.line === st.line && r.part === st.part && !r.image) }))
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
    // Keyed by line and position on it: a line may hold several statements.
    const byStatement = new Map(geo.rects.filter(r => !r.image)
        .map(r => [`${r.line}:${r.part}`, { ...r, scale: LENGTH_UNITS[geo.units], units: geo.units }]));
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
        el('div', { class: 'custom-form-title', text: 'Shapes' },
            el('span', { class: 'custom-hint', text: '  Fields take expressions of the sidebar parameters: w/2, h1+h2, 35um. -inf / inf runs an edge to the boundary, a negative size flips the rectangle to the other side of its position. ⇋ mirror adds the mirror image about x=0. A later dielectric covers an earlier one. Trapezoids, n-gons, ellipses and rounded or hollow shapes need the full-wave solver.' })),
        el('div', { class: 'custom-rect-list' }, ...rects.map((r, i) => rectRow(model, r, byStatement.get(`${r.line}:${r.part}`), i, rects.length))),
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
    if (source === 'load') { openPlatingRows.clear(); openCornerRows.clear(); selectedLine = 0; return; }
    if (before === null || before === after) return;
    const map = lineMap(before, after);
    for (const set of [openPlatingRows, openCornerRows]) {
        const open = [...set].map(map).filter(l => l > 0);
        set.clear();
        open.forEach(l => set.add(l));
    }
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
    pendingSource = source;
    debounceTimer = setTimeout(() => { debounceTimer = null; refresh(source); }, 250);
}

// Applies an edit still waiting on the debounce, so a solve started right after
// typing sees it (and its solver is not replaced by the late refresh).
export function flushCustomGeometryEdits() {
    if (!debounceTimer) return;
    clearTimeout(debounceTimer);
    debounceTimer = null;
    refresh(pendingSource);
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

export function initCustomGeometryEditor({ onGeometryChange, onHighlightChange, log }) {
    onChange = onGeometryChange;
    onHighlight = onHighlightChange;
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
            // Values in a statement cannot hold spaces, an empty field is automatic, and
            // ';' or '#' would end the statement.
            const inputs = DOMAIN_KEYS.map(k => $(`custom_domain_${k}`));
            inputs.forEach(i => markStatementBreak(i, /[;#]/.test(i.value)));
            if (inputs.some(i => /[;#]/.test(i.value))) return;
            const values = inputs.map(i => i.value.replace(/\s+/g, '') || 'auto');
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
