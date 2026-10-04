// Browser end-to-end: the custom geometry type through the real UI. Covers what the node
// tests cannot reach: the editor under the geometry preview, the generated parameter
// inputs, live validation, conversion of a fixed type, the compressed link, a solve on
// both solvers and a parameter sweep over a geometry parameter.
//
// Needs the dev server on localhost:8731 (see tests/run.mjs e2e tier).
import { launch, URL } from './e2e_helpers.mjs';
import { check, done } from './helpers.js';

const { browser, page, errors } = await launch({ viewport: { width: 1400, height: 900 }, clipboard: '__copied' });
page.on('dialog', d => d.accept());
const geomText = () => page.evaluate(() => document.getElementById('custom_geom_text').value);
// The input of the form cell labelled `label` inside `row`.
const cellInput = (row, label, p = page) =>
    row.locator('.custom-cell', { has: p.locator('.custom-cell-label', { hasText: label }) }).first().locator('input');
const logText = () => page.evaluate(() => document.getElementById('console_out').textContent);
const shot = process.env.SHOT_DIR;
async function solveAndRead() {
    await page.click('#btn_solve');
    await page.waitForFunction(() => document.getElementById('btn_solve').textContent === 'Solve'
        && !document.getElementById('btn_solve').disabled, null, { timeout: 180000 });
    return page.evaluate(() => {
        const rows = [...document.querySelectorAll('#results-table tr, .results-table tr')].map(r => r.textContent);
        return { rows, log: document.getElementById('console_out').textContent };
    });
}
// Last "Z0: <value> Ohm" of the log.
const z0FromLog = (text) => {
    const m = [...text.matchAll(/Z0:\s*([\d.]+)\s*Ohm/gi)];
    return m.length ? parseFloat(m[m.length - 1][1]) : NaN;
};

await page.goto(URL, { waitUntil: 'networkidle' });

// ---- native microstrip result, for the conversion check ----
const nativeRun = await solveAndRead();
const zNative = z0FromLog(nativeRun.log);
check('native microstrip solves', zNative > 10, `Z0 ${zNative}`);

// ---- convert to custom ----
await page.click('#btn-convert-custom');
await page.waitForTimeout(600);
const conv = await page.evaluate(() => ({
    type: document.getElementById('tl_type').value,
    text: document.getElementById('custom_geom_text').value,
    editor: getComputedStyle(document.getElementById('custom-editor')).display,
    sidebar: document.getElementById('custom-params').style.display,
    ms: document.getElementById('microstrip-params').style.display,
    errors: document.getElementById('custom-geom-errors').textContent,
    convertBtn: document.getElementById('btn-convert-custom').style.display,
}));
check('convert switches to the custom type', conv.type === 'custom' && conv.editor === 'flex'
    && conv.sidebar === 'block' && conv.ms === 'none' && conv.convertBtn === 'none');
check('converted text has the rectangles', /sig\+/.test(conv.text) && /gnd/.test(conv.text) && /diel/.test(conv.text)
    && conv.errors === '', conv.errors);
const convRun = await solveAndRead();
const zConv = z0FromLog(convRun.log);
check('converted geometry solves to the native impedance', Math.abs(zConv - zNative) / zNative < 2e-3,
    `${zConv} vs ${zNative}`);
if (shot) await page.screenshot({ path: `${shot}/custom_converted.png` });

// ---- template, parameters, live validation ----
await page.selectOption('#custom-template', 'Microstrip on a finite ground');
await page.waitForTimeout(600);
const tpl = await page.evaluate(() => ({
    params: [...document.querySelectorAll('#custom-param-list input')].map(i => [i.id, i.value, i.disabled]),
    bounds: ['left', 'right', 'top', 'bottom'].map(w => document.getElementById('custom_bound_' + w).value),
    sweep: [...document.getElementById('sweep-x-selector').options].map(o => o.value),
    errors: document.getElementById('custom-geom-errors').textContent,
}));
check('template fills the parameter inputs', tpl.params.length === 5
    && tpl.params.some(p => p[0] === 'inp_cgp_w' && p[1] === '0.35'), JSON.stringify(tpl.params));
check('boundary selects follow the text', tpl.bounds.join() === 'open,open,open,open');
check('sweep list offers the geometry parameters', tpl.sweep.includes('cgp_w') && tpl.sweep.includes('cgp_wgnd')
    && !tpl.sweep.includes('w') && tpl.sweep.includes('custom_sigma'), tpl.sweep.join());

// ---- Solve right after a text edit solves the edited text ----
// The edit is still waiting on the editor's debounce when Solve is clicked; the stale
// sidebar input used to override it and the late refresh dropped the field plot.
{
    const zBefore = z0FromLog((await solveAndRead()).log);
    const orig = await geomText();
    // The text view may be hidden behind the form view: write it as a typed edit would.
    const typeText = (v) => page.evaluate((v) => {
        const t = document.getElementById('custom_geom_text');
        t.value = v;
        t.dispatchEvent(new Event('input'));
    }, v);
    await typeText(orig.replace(/^(\s*w\s*=\s*)0\.35/m, '$10.7'));
    const zAfter = z0FromLog((await solveAndRead()).log);
    check('Solve clicked during the edit debounce uses the edited width', zAfter < 0.9 * zBefore,
        `Z0 ${zBefore} -> ${zAfter}`);
    check('and the sidebar input follows the text', await page.inputValue('#inp_cgp_w') === '0.7');
    await typeText(orig);
    await page.waitForTimeout(600);
}

// ---- form view: one row per rectangle, editing the same text ----
const form = await page.evaluate(() => ({
    visible: getComputedStyle(document.getElementById('custom-form')).display !== 'none'
        && getComputedStyle(document.getElementById('custom_geom_text')).display === 'none',
    kinds: [...document.querySelectorAll('#custom-form .custom-rect-row select.custom-kind')].map(s => s.value),
    params: [...document.querySelectorAll('#custom-param-list .custom-param-name-btn')].map(b => b.textContent),
    formSections: document.querySelectorAll('#custom-form .custom-form-section').length,
    domain: document.getElementById('custom-domain-resolved').textContent,
    air: (document.getElementById('sim_canvas').layout?.shapes || []).filter(s => s.line && s.line.dash === 'dot').length,
}));
check('the form view is the default and lists the rectangles', form.visible && form.kinds.join() === 'diel,sig+,gnd'
    && form.params.join() === 'w,t,h,wsub,wgnd' && form.formSections === 1, JSON.stringify(form));
check('the automatic solved region is shown and drawn as air', /x -[\d.]+ … [\d.]+, y -[\d.]+ … [\d.]+ mm/.test(form.domain) && form.air === 1, form.domain);

// ---- sidebar: parameters, units and solved region are edited there only ----
{
    await page.click('#btn-custom-add-param');
    await page.waitForTimeout(400);
    check('+ Add parameter appends a definition and focuses its field', /^p1 = 1$/m.test(await geomText())
        && await page.evaluate(() => document.activeElement?.id === 'inp_cgp_p1'));
    await page.fill('#inp_cgp_p1', 'w + wgnd');
    await page.waitForTimeout(600);
    const expr = await page.evaluate(() => ({
        hint: document.querySelector('#custom-param-list [data-param="p1"] .custom-param-value').textContent,
        sweep: [...document.getElementById('sweep-x-selector').options].map(o => o.value),
    }));
    check('a parameter takes an expression, shows its value and leaves the sweep list',
        /^p1 = w \+ wgnd$/m.test(await geomText()) && expr.hint === '= 0.85' && !expr.sweep.includes('cgp_p1'), JSON.stringify(expr));
    await page.click('#custom-param-list [data-param="wsub"] .custom-param-name-btn');
    await page.fill('.custom-param-rename', 'ws');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(600);
    const renamed = await geomText();
    check('renaming a parameter renames its uses', /wsub/.test(renamed) === false && /ws = 3/.test(renamed)
        && /x=-ws\/2\s+y=0\s+w=ws/.test(renamed)
        && await page.evaluate(() => document.getElementById('custom-geom-errors').textContent === ''), renamed);
    await page.click('#custom-param-list [data-param="ws"] .custom-param-name-btn');
    await page.fill('.custom-param-rename', 'wsub');
    await page.keyboard.press('Enter');
    await page.waitForTimeout(400);
    await page.fill('#inp_cgp_p1', 'nope');
    await page.waitForTimeout(600);
    check('a bad expression marks its sidebar row', await page.evaluate(() =>
        document.querySelector('#custom-param-list [data-param="p1"]').classList.contains('has-error')
        && !document.querySelector('#custom-param-list [data-param="w"]').classList.contains('has-error')));
    await page.click('#custom-param-list [data-param="p1"] .custom-row-btn');
    await page.waitForTimeout(600);
    check('deleting a parameter removes its definition', !/p1/.test(await geomText())
        && await page.evaluate(() => !document.getElementById('inp_cgp_p1')
            && document.getElementById('custom-geom-errors').textContent === ''));

    await page.fill('#custom_domain_y2', '3');
    await page.waitForTimeout(600);
    const dom = await page.evaluate(() => document.getElementById('custom-domain-resolved').textContent);
    check('a solved region field writes the domain statement', /^domain auto auto auto 3/m.test(await geomText()) && /y -[\d.]+ … 3 mm/.test(dom), dom);
    await page.fill('#custom_domain_y2', '');
    await page.waitForTimeout(600);
    check('clearing it returns to domain auto', /^domain auto(\s|$)/m.test(await geomText()));
    const before = await geomText();
    // A unit change asks whether to convert the numbers; converting keeps the size.
    await page.selectOption('#custom-units', 'um');
    await page.click('dialog.custom-choice button:has-text("Convert to um")');
    await page.waitForTimeout(600);
    check('the unit select converts the numbers to the new unit', /^units um$/m.test(await geomText())
        && await page.inputValue('#inp_cgp_w') === '350', await page.inputValue('#inp_cgp_w'));
    await page.selectOption('#custom-units', 'mm');
    await page.click('dialog.custom-choice button:has-text("Convert to mm")');
    await page.waitForTimeout(600);
    check('the sidebar edits leave the text as it was', (await geomText()) === before);
}

// ---- errors on their row, undo and redo ----
{
    const start = await geomText();
    const hCell = cellInput(page.locator('#custom-form .custom-rect-row.kind-sigp'), /^h$/);
    await hCell.fill('tt');
    await page.waitForTimeout(600);
    const inline = await page.evaluate(() => ({
        row: document.querySelector('#custom-form .custom-rect-row.kind-sigp .custom-inline-error')?.textContent ?? '',
        others: document.querySelectorAll('#custom-form .custom-rect-row:not(.kind-sigp) .custom-inline-error').length,
        box: getComputedStyle(document.getElementById('custom-geom-errors')).display,
        count: document.getElementById('custom-error-count').textContent,
        solve: document.getElementById('btn_solve').disabled,
    }));
    check('a form error is shown under its row only, with a count in the bar', /unknown parameter 'tt'/.test(inline.row)
        && !/line \d/.test(inline.row) && inline.others === 0 && inline.box === 'none' && inline.count === '1 error' && inline.solve,
        JSON.stringify(inline));
    await page.click('#btn-custom-undo');
    await page.waitForTimeout(500);
    check('undo takes the edit back and clears the error', (await geomText()) === start && await page.evaluate(() =>
        !document.querySelector('#custom-form .custom-inline-error')
        && getComputedStyle(document.getElementById('custom-error-count')).display === 'none'
        && !document.getElementById('btn_solve').disabled));

    const rows = () => page.evaluate(() => document.querySelectorAll('#custom-form .custom-rect-row').length);
    const n0 = await rows();
    await page.click('#custom-form .custom-rect-row.kind-gnd .custom-row-btn[title="Delete"]');
    await page.waitForTimeout(500);
    const n1 = await rows();
    await page.locator('#custom-form .custom-rect-row').first().click({ position: { x: 3, y: 3 } });
    await page.keyboard.press('Control+z');
    await page.waitForTimeout(500);
    const n2 = await rows();
    await page.keyboard.press('Control+Shift+z');
    await page.waitForTimeout(500);
    const n3 = await rows();
    await page.click('#btn-custom-undo');
    await page.waitForTimeout(500);
    check('Ctrl+Z restores a deleted rectangle, redo deletes it again', n1 === n0 - 1 && n2 === n0 && n3 === n0 - 1
        && (await geomText()) === start, `${n0} ${n1} ${n2} ${n3}`);

    await page.fill('#inp_cgp_w', '0.5');
    await page.waitForTimeout(500);
    await page.click('#btn-custom-undo');
    await page.waitForTimeout(500);
    check('undo covers a sidebar edit and refills its field', (await geomText()) === start
        && await page.evaluate(() => document.getElementById('inp_cgp_w').value === '0.35'));
}

// Edit the trace width expression in its row, the text follows and the preview highlights the row.
const wCell = cellInput(page.locator('#custom-form .custom-rect-row.kind-sigp'), /^w$/);
await wCell.fill('2*w');
await page.waitForTimeout(600);
const afterCell = await page.evaluate(() => ({
    text: document.getElementById('custom_geom_text').value,
    focusKept: document.activeElement?.classList.contains('custom-expr'),
    hl: (document.getElementById('sim_canvas').layout?.shapes || []).filter(s => s.line && /56, 189, 248/.test(s.line.color)).length,
    width: (document.getElementById('sim_canvas').layout?.shapes || []).filter(s => s.fillcolor === '#d97706').map(s => +(s.x1 - s.x0).toFixed(3)),
}));
check('a form field edit rewrites its statement and keeps the focus', /sig\+\s+x=-w\/2 w=2\*w y=h h=t/.test(afterCell.text) && afterCell.focusKept,
    afterCell.text.split('\n').find(l => l.startsWith('sig+')));
check('the edited row is outlined and the preview follows', afterCell.hl === 1 && afterCell.width.includes(0.7), afterCell.width.join());
await wCell.fill('w');
await page.waitForTimeout(600);
// ';' and '#' would end the statement: the field is flagged and the text keeps its value.
{
    const before = await geomText();
    await wCell.fill('w#2');
    await page.waitForTimeout(600);
    check('a form value with # is flagged and not written', await geomText() === before
        && await wCell.evaluate(i => i.classList.contains('invalid')));
    await wCell.fill('w');
    await page.waitForTimeout(600);
    check('and the flag clears with a valid value', await geomText() === before
        && !(await wCell.evaluate(i => i.classList.contains('invalid'))));
}

// A click on the row itself, outside its fields, outlines the rectangle too.
await page.locator('#custom-form .custom-rect-row.kind-gnd').click({ position: { x: 12, y: 2 } });
await page.waitForTimeout(400);
check('clicking a row outside its fields outlines its rectangle', await page.evaluate(() =>
    (document.getElementById('sim_canvas').layout?.shapes || []).filter(s => s.line && /56, 189, 248/.test(s.line.color)).length === 1
    && document.activeElement.classList.contains('custom-rect-row')));

// The color band opens a popup that writes color=, and shows the new color when reopened.
{
    const band = '#custom-form .custom-rect-row.kind-gnd .custom-color-band';
    const shown = () => page.evaluate(() => [document.querySelector('.custom-color-popup input[type=color]')?.value,
        document.querySelector('.custom-color-hex')?.value]);
    await page.click(band);
    await page.fill('.custom-color-hex', '#123456');
    await page.waitForTimeout(600);
    const written = /^gnd .*color=#123456/m.test(await geomText());
    await page.keyboard.press('Escape');
    await page.click(band);
    const reopened = await shown();
    await page.locator('.custom-color-popup button', { hasText: 'Default' }).click();
    await page.waitForTimeout(600);
    const cleared = !/^gnd .*color=/m.test(await geomText());
    await page.click(band);
    const afterDefault = await shown();
    await page.keyboard.press('Escape');
    check('the color popup writes color=, shows it again on reopening and clears it with Default',
        written && reopened.join() === '#123456,#123456' && cleared && afterDefault[0] !== '#123456',
        `${written} ${reopened} ${cleared} ${afterDefault}`);
}

// Per-conductor metal and finish: own conductivity, roughness and plating material on
// the trace row.
const sigRow = page.locator('#custom-form .custom-rect-row.kind-sigp');
await cellInput(sigRow, /^σ$/).fill('4.1e7');
// A unit typed with a space goes into the text without it.
await cellInput(sigRow, /^rq$/).fill('1 um');
await page.waitForTimeout(500);
check('a field value with a space before the unit parses', await page.evaluate(() =>
    /rq=1um/.test(document.getElementById('custom_geom_text').value)
    && document.getElementById('custom-geom-errors').textContent === ''));
await cellInput(sigRow, /^rq$/).fill('0.001');
check('the sidebar has no surface plating block in the custom type, only Model Thick Plating', await page.evaluate(() =>
    getComputedStyle(document.querySelector('.control-group:has(#chk_plating)')).display === 'none'
    && [...document.querySelectorAll('#plating-params > .control-group')]
        .every(g => (getComputedStyle(g).display === 'none') !== !!g.querySelector('#chk_plating_thick_corners'))));
const platingToggle = sigRow.locator('.custom-plating-toggle');
check('the plating options are collapsed on an unplated conductor',
    !(await sigRow.locator('.custom-plating-panel').isVisible()) && /▸ plating$/.test(await platingToggle.textContent()));
await platingToggle.click();
await sigRow.locator('.custom-face', { hasText: 'top' }).locator('input').check();
await page.waitForTimeout(500);
check('the first plated face starts from a typical plating material', await page.evaluate(() =>
    /plating=top plating_sigma=1e7 plating_t=4um/.test(document.getElementById('custom_geom_text').value)
    && document.getElementById('custom-geom-errors').textContent === ''));
await cellInput(sigRow.locator('.custom-plating-material'), /^σ$/).fill('1e7');
await cellInput(sigRow.locator('.custom-plating-material'), /^t$/).fill('4um');
await page.waitForTimeout(600);
const finish = await page.evaluate(() => ({
    line: document.getElementById('custom_geom_text').value.split('\n').find(l => l.startsWith('sig+')),
    errors: document.getElementById('custom-geom-errors').textContent,
    gold: (document.getElementById('sim_canvas').layout?.shapes || []).filter(s => s.type === 'line' && /255, 215, 0/.test(s.line.color)).length,
}));
await platingToggle.click();
check('the collapsed plating button names the plated faces',
    !(await sigRow.locator('.custom-plating-panel').isVisible()) && /▸ plating T$/.test(await platingToggle.textContent()));
await platingToggle.click();
check('conductivity, roughness and plating fields of a row write their keys', /sigma=4\.1e7 rq=0\.001 plating=top plating_sigma=1e7 plating_t=4um/.test(finish.line)
    && finish.errors === '' && finish.gold === 1, finish.line + ' | ' + finish.errors);
for (const label of [/^σ$/, /^t$/]) await cellInput(sigRow.locator('.custom-plating-material'), label).fill('');
await sigRow.locator('.custom-face', { hasText: 'top' }).locator('input').uncheck();
await cellInput(sigRow, /^rq$/).fill('');
await cellInput(sigRow, /^σ$/).fill('');
await page.waitForTimeout(600);
check('clearing the fields removes the keys', await page.evaluate(() =>
    /^sig\+\s+x=-w\/2 w=w y=h h=t$/m.test(document.getElementById('custom_geom_text').value)));

// Add a dielectric, make it a cover layer, then delete it again.
await page.click('#custom-form .custom-adders button:has-text("+ Dielectric")');
await page.waitForTimeout(400);
const added = await page.evaluate(() => [...document.querySelectorAll('#custom-form .custom-rect-row')].length);
await cellInput(page.locator('#custom-form .custom-rect-row').last(), /^y$/).fill('h');
await cellInput(page.locator('#custom-form .custom-rect-row').last(), /^h$/).fill('0.1');
await page.waitForTimeout(600);
const cover = await geomText();
check('a rectangle added in the form appears in the text', added === 4 && /diel\s+x=-inf w=inf y=h h=0\.1 er=4\.4 tand=0\.02/.test(cover));
await page.locator('#custom-form .custom-rect-row').last().locator('button[title="Delete"]').click();
await page.waitForTimeout(400);
check('deleting the row removes its line', await page.evaluate(() =>
    document.querySelectorAll('#custom-form .custom-rect-row').length === 3 && !/y=h h=0\.1/.test(document.getElementById('custom_geom_text').value)));

// A ground run to an open boundary is a normal reference plane and says nothing.
const gndX = cellInput(page.locator('#custom-form .custom-rect-row.kind-gnd'), /^w$/);
await gndX.fill('inf');
await page.waitForTimeout(600);
check('a ground reaching an open boundary does not warn',
    await page.evaluate(() => document.getElementById('custom-geom-warnings').style.display === 'none'));
await gndX.fill('wgnd');
await page.waitForTimeout(600);

// A signal run to an open boundary is reported as a warning.
const sigW = cellInput(sigRow, /^w$/);
await sigW.fill('inf');
await page.waitForTimeout(600);
const warn = await page.evaluate(() => ({ text: document.getElementById('custom-geom-warnings').textContent,
    solve: document.getElementById('btn_solve').disabled }));
check('a signal reaching an open boundary warns without blocking Solve',
    /Signal conductor \(line \d+\) reaches the open right boundary/.test(warn.text) && warn.solve === false, warn.text.slice(0, 120));
await sigW.fill('w');
await page.waitForTimeout(600);
check('the warning clears', await page.evaluate(() => document.getElementById('custom-geom-warnings').style.display === 'none'));
if (shot) await page.screenshot({ path: `${shot}/custom_form.png` });

await page.fill('#inp_cgp_wgnd', '2');
await page.waitForTimeout(600);
const edited = await geomText();
check('a sidebar parameter edit rewrites the text', /wgnd = 2\b/.test(edited) && /wsub = 3;/.test(edited));

await page.selectOption('#custom_bound_bottom', 'gnd');
await page.waitForTimeout(600);
const bnd = await geomText();
check('a boundary select rewrites the bounds statement', /bounds open open open gnd/.test(bnd));
await page.selectOption('#custom_bound_bottom', 'open');
await page.waitForTimeout(600);

await page.click('#btn-custom-view-text');
const good = await geomText();
await page.fill('#custom_geom_text', good.replace('w=wgnd', 'w=oops'));
await page.waitForTimeout(600);
const bad = await page.evaluate(() => ({
    errors: document.getElementById('custom-geom-errors').textContent,
    shown: document.getElementById('custom-geom-errors').style.display,
    solveDisabled: document.getElementById('btn_solve').disabled,
    sweepDisabled: document.getElementById('btn-run-sweep').disabled,
    log: document.getElementById('console_out').textContent,
    badge: getComputedStyle(document.getElementById('custom-stale-badge')).display,
}));
check('the editor error stays out of the log and the preview is marked stale',
    !/oops/.test(bad.log) && bad.badge === 'block', bad.log.slice(-200));
check('an invalid text lists the error with its line and blocks Solve',
    /line \d+: unknown parameter 'oops'/.test(bad.errors) && bad.shown === 'block' && bad.solveDisabled === true, bad.errors);
check('an invalid text blocks Run Sweep too', bad.sweepDisabled === true);
await page.fill('#custom_geom_text', good.replace(/^(gnd .*)y=-t/m, '$1y=h-t/2'));
await page.waitForTimeout(600);
const shorted = await page.evaluate(() => document.getElementById('custom-geom-errors').textContent);
check('a short between signal and ground is reported', /shorted/.test(shorted), shorted);
await page.fill('#custom_geom_text', good);
await page.waitForTimeout(600);
// A click on a listed error selects its line, and the format reference opens beside the text.
{
    const good2 = await geomText();
    await page.fill('#custom_geom_text', good2.replace('w=wgnd', 'w=oops'));
    await page.waitForTimeout(600);
    await page.click('#custom-geom-errors .custom-error-item');
    const sel = await page.evaluate(() => {
        const t = document.getElementById('custom_geom_text');
        return { focus: document.activeElement === t, picked: t.value.slice(t.selectionStart, t.selectionEnd) };
    });
    check('a click on a listed error selects its line', sel.focus && /w=oops/.test(sel.picked), sel.picked);
    await page.click('#btn-custom-format');
    const ref = await page.evaluate(() => ({
        shown: getComputedStyle(document.getElementById('custom-format-ref')).display,
        text: document.getElementById('custom-format-ref').textContent,
        area: document.getElementById('custom_geom_text').getBoundingClientRect().width,
    }));
    check('the format reference opens beside the text', ref.shown === 'block' && /bounds left right top bottom/.test(ref.text) && ref.area > 300,
        `${ref.shown} ${ref.area}`);
    await page.click('#btn-custom-format');
    await page.fill('#custom_geom_text', good2);
    await page.waitForTimeout(600);
}
check('a valid text clears the errors and enables Solve', await page.evaluate(() =>
    document.getElementById('custom-geom-errors').style.display === 'none' && !document.getElementById('btn_solve').disabled
    && !document.getElementById('btn-run-sweep').disabled));
// A text made invalid while a solve runs keeps Solve off once the solve ends.
{
    const good3 = await geomText();
    await page.click('#btn_solve');
    await page.fill('#custom_geom_text', good3.replace('w=wgnd', 'w=oops'));
    await page.waitForFunction(() => document.getElementById('btn_solve').textContent === 'Solve', null, { timeout: 180000 });
    check('a text made invalid during a solve keeps Solve off after it', await page.evaluate(() =>
        document.getElementById('btn_solve').disabled));
    await page.fill('#custom_geom_text', good3);
    await page.waitForTimeout(600);
    check('fixing that text enables Solve again', await page.evaluate(() => !document.getElementById('btn_solve').disabled));
}
check('the stale badge goes with the errors', await page.evaluate(() =>
    getComputedStyle(document.getElementById('custom-stale-badge')).display === 'none'));

// ---- panes: log collapse, log and editor resize ----
{
    const drag = async (sel, dy) => {
        const b = await page.locator(sel).boundingBox();
        await page.mouse.move(b.x + b.width / 2, b.y + b.height / 2);
        await page.mouse.down();
        await page.mouse.move(b.x + b.width / 2, b.y + b.height / 2 + dy, { steps: 4 });
        await page.mouse.up();
    };
    const h = (id) => page.evaluate((i) => document.getElementById(i).getBoundingClientRect().height, id);
    const e0 = await h('custom-editor');
    await drag('#custom-splitter', -80);
    const e1 = await h('custom-editor');
    check('dragging the splitter resizes the editor', Math.abs(e1 - e0 - 80) < 6, `${e0} -> ${e1}`);
    // A solve has opened the log by now. Collapse it, then open it again.
    await page.click('#log-bar');
    check('a click on the log bar collapses the log and shows the last line', await page.evaluate(() =>
        getComputedStyle(document.getElementById('console_out')).display === 'none'
        && getComputedStyle(document.getElementById('log-last')).visibility === 'visible'
        && document.getElementById('log-last').textContent.length > 0));
    await page.click('#log-bar');
    const l0 = await h('console_out');
    await drag('#log-splitter', -60);
    const l1 = await h('console_out');
    check('the log opens on a click and resizes by its handle', l0 > 0 && Math.abs(l1 - l0 - 60) < 6, `${l0} -> ${l1}`);
}

// ---- highlight of the rectangle under the cursor ----
await page.evaluate(() => {
    const el = document.getElementById('custom_geom_text');
    const pos = el.value.indexOf('sig+');
    el.focus(); el.setSelectionRange(pos, pos);
    el.dispatchEvent(new Event('keyup'));
});
await page.waitForTimeout(300);
const hl = await page.evaluate(() => (document.getElementById('sim_canvas').layout?.shapes || [])
    .filter(s => s.line && /56, 189, 248/.test(s.line.color)).length);
check('the rectangle on the cursor line is outlined in the preview', hl === 1, `${hl} highlight shapes`);

// Back in the form view, the rows are rebuilt from the edited text.
await page.click('#btn-custom-view-form');
await page.waitForTimeout(300);
check('switching back to the form rebuilds it from the text', await page.evaluate(() =>
    document.querySelectorAll('#custom-form .custom-rect-row').length === 3));

// ---- solve on both solvers ----
await page.fill('#inp_cgp_wgnd', '0.5');
await page.waitForTimeout(600);
const qs = z0FromLog((await solveAndRead()).log);
check('finite-ground microstrip solves quasi-static', Math.abs(qs - 59.7) < 1.5, `Z0 ${qs}`);
await page.selectOption('#mesh_backend', 'fullwave_mqs');
await page.waitForTimeout(300);
const fwRun = await solveAndRead();
const fw = z0FromLog(fwRun.log);
// The triangular backend logs triangle counts, the FDM node counts.
const fwPasses = (fwRun.log.split('Starting simulation').pop().match(/Tris=\d+/g) || []).length;
check('finite-ground microstrip solves full-wave', Math.abs(fw - 59.8) < 1.5 && fwPasses > 0, `Z0 ${fw}, ${fwPasses} triangular passes`);
if (shot) await page.screenshot({ path: `${shot}/custom_solved.png` });
await page.selectOption('#mesh_backend', 'rectilinear');

// ---- link round trip ----
await page.click('[data-tab="results"]');
await page.click('#copy-link-btn');
await page.waitForFunction(() => !!window.__copied, null, { timeout: 5000 });
const link = await page.evaluate(() => window.__copied);
check('custom geometry link uses the compressed fragment', /#params=z\./.test(link) && link.length < 1500, `${link.length} chars`);
const page2 = await browser.newPage();
await page2.goto(link, { waitUntil: 'networkidle' });
await page2.waitForTimeout(800);
const restored = await page2.evaluate(() => ({
    type: document.getElementById('tl_type').value,
    text: document.getElementById('custom_geom_text').value,
    wgnd: document.getElementById('inp_cgp_wgnd')?.value,
}));
check('the link restores the geometry text', restored.type === 'custom' && restored.text === good.replace('wgnd = 2', 'wgnd = 0.5')
    && restored.wgnd === '0.5', JSON.stringify(restored).slice(0, 200));
check('the applied link leaves the address bar, so a reload keeps later edits', !/#params=/.test(page2.url()), page2.url());
await page2.close();

// A link that inflates far past any geometry is refused, not decompressed into memory.
{
    const { deflateRawSync } = await import('zlib');
    const bomb = 'z.' + deflateRawSync(Buffer.alloc(64 << 20, 32)).toString('base64')
        .replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
    const pb = await browser.newPage();
    await pb.goto(`${URL}#params=${bomb}`, { waitUntil: 'networkidle' });
    await pb.waitForTimeout(1500);
    const st = await pb.evaluate(() => ({ type: document.getElementById('tl_type').value,
        log: document.getElementById('console_out').textContent }));
    check('an oversized link is refused', st.type === 'microstrip' && /Failed to parse URL parameters: the link decompresses to more than/.test(st.log),
        st.log.split('\n').filter(l => /URL/.test(l)).join(' | '));
    await pb.close();
}

// The same link applied to a running page (pasted into its address bar: only the hash
// changes) after a solve of the native type. Two geometry draws in one tick used to
// replot each other's Plotly config forever and hang the tab.
{
    const p4 = await browser.newPage({ viewport: { width: 1400, height: 900 } });
    await p4.goto(URL, { waitUntil: 'networkidle' });
    p4.evaluate(() => {
        const s = document.getElementById('tl_type');
        for (const v of ['stripline', 'microstrip']) { s.value = v; s.dispatchEvent(new Event('change', { bubbles: true })); }
    }).catch(() => {});
    await p4.waitForTimeout(1500);
    const alive = () => Promise.race([p4.evaluate(() => document.getElementById('tl_type').value),
        new Promise(r => setTimeout(() => r('HUNG'), 5000))]);
    check('two line type changes in one tick leave the page responsive', await alive() === 'microstrip');
    p4.evaluate(h => { location.hash = h; }, link.split('#')[1]).catch(() => {});
    await p4.waitForTimeout(1500);
    const t = await alive();
    check('a custom link applied to a running page restores it', t === 'custom', t);
    await p4.close();
}

// An old-style link of a fixed type still loads.
const page3 = await browser.newPage();
await page3.goto(`${URL}?params=${Buffer.from(encodeURIComponent(JSON.stringify({ tl_type: 'gcpw', w: 0.5 }))).toString('base64')}`,
    { waitUntil: 'networkidle' });
const old = await page3.evaluate(() => [document.getElementById('tl_type').value, document.getElementById('inp_w').value]);
check('a ?params= link of a fixed type still loads', old[0] === 'gcpw' && /^0\.5/.test(old[1]), old.join());
await page3.close();

// ---- parameter sweep over a geometry parameter ----
await page.click('[data-tab="sweep"]');
await page.selectOption('#sweep-x-selector', 'cgp_wgnd');
await page.fill('#sweep-x-min', '0.4');
await page.fill('#sweep-x-max', '2');
await page.fill('#sweep-points', '3');
await page.click('#btn-run-sweep');
await page.waitForFunction(() => /Sweep complete: 3 points/.test(document.getElementById('console_out').textContent),
    null, { timeout: 240000 });
const sweep = await page.evaluate(() => {
    const d = document.getElementById('sweep-plot').data || [];
    return { y: d[0] ? Array.from(d[0].y) : [], text: document.getElementById('custom_geom_text').value,
        input: document.getElementById('inp_cgp_wgnd').value };
});
check('sweeping the ground width lowers the impedance', sweep.y.length >= 3 && sweep.y[0] > sweep.y[sweep.y.length - 1] + 5,
    `${sweep.y[0].toFixed(2)} -> ${sweep.y[sweep.y.length - 1].toFixed(2)} over ${sweep.y.length} plot points`);
check('the sweep leaves the text and the input as they were', /wgnd = 0\.5\b/.test(sweep.text) && sweep.input === '0.5');
// Points the solver rejects are skipped, the sweep goes on: an n-gon is valid text but
// full-wave only, so every quasi-static point fails at the mesh.
{
    const plain = sweep.text;
    const setText = t => page.evaluate(v => {
        const el = document.getElementById('custom_geom_text');
        el.value = v;
        el.dispatchEvent(new Event('input', { bubbles: true }));
    }, t);
    await setText(plain + '\ndiel ngon x=0 y=3 r=0.2 n=6 er=2\n');
    await page.waitForTimeout(600);
    const from = (await logText()).length;
    await page.click('#btn-run-sweep');
    await page.waitForFunction(n => /Sweep complete: \d+ points|Sweep error/.test(
        document.getElementById('console_out').textContent.slice(n)), from, { timeout: 240000 });
    const tail = (await logText()).slice(from);
    check('a sweep skips the points the solver rejects', /Point 1 skipped: /.test(tail) && /Point 3 skipped: /.test(tail)
        && /Sweep complete: 0 points/.test(tail), tail.slice(-400));
    await setText(plain);
    await page.waitForTimeout(600);
}
if (shot) await page.screenshot({ path: `${shot}/custom_sweep.png` });

// ---- Modes tab: the field plot covers the air above and below the conductors ----
await page.click('[data-tab="geometry"]');
await page.selectOption('#custom-template', 'Differential CPW over air');
await page.waitForTimeout(600);
await page.click('.tab-button[data-tab="modes"]');
await page.waitForTimeout(500);
await page.evaluate(() => { document.getElementById('modes-freq').value = '10 GHz'; document.getElementById('modes-nev').value = '3'; });
await page.click('#btn-solve-modes');
await page.waitForFunction(() => document.querySelectorAll('#modes-list tr[data-idx]').length > 0, null, { timeout: 180000 });
await page.waitForTimeout(1500);
const modes = await page.evaluate(() => {
    // The field is an image of the mode mesh (the hover markers follow the view, not the mesh).
    const m = document.getElementById('modes-plot')._modesMesh;
    const T = m ? m.blocks[0].tris : [];
    const y = [];
    for (let k = 1; k < T.length; k += 2) y.push(T[k] * 1000);
    return { yMin: Math.min(...y), yMax: Math.max(...y), n: T.length / 6 };
});
// Conductors span y = 0 ... 0.017 mm, the substrate reaches down to -0.635 mm.
check('the mode field extends into the air above and below the structure', modes.n > 10 && modes.yMax > 1 && modes.yMin < -1,
    `y ${modes.yMin.toFixed(2)} … ${modes.yMax.toFixed(2)} mm`);
if (shot) await page.screenshot({ path: `${shot}/custom_modes.png` });

// ---- back to a fixed type ----
await page.selectOption('#tl_type', 'microstrip');
await page.waitForTimeout(300);
check('leaving the custom type restores the fixed sidebar', await page.evaluate(() =>
    getComputedStyle(document.getElementById('custom-editor')).display === 'none'
    && document.getElementById('microstrip-params').style.display === 'block'
    && !document.getElementById('btn_solve').disabled));

// ---- page reload with the custom type selected ----
// Browsers restore the type select on reload, and not all of them the hidden textarea.
// The init script stands in for the restored select; the text comes from the session.
{
    const p2 = await browser.newPage({ viewport: { width: 1400, height: 900 } });
    await p2.addInitScript(() => document.addEventListener('DOMContentLoaded', () => {
        if (sessionStorage.getItem('e2e_reloaded')) document.getElementById('tl_type').value = 'custom';
    }, true));
    await p2.goto(URL);
    await p2.waitForTimeout(1500);
    await p2.selectOption('#tl_type', 'custom');
    await p2.waitForTimeout(800);
    await p2.evaluate(() => {
        const t = document.getElementById('custom_geom_text');
        t.value = t.value.replace(/^w = [^;\n]*/m, 'w = 0.41');
        t.dispatchEvent(new Event('input'));
        sessionStorage.setItem('e2e_reloaded', '1');
    });
    await p2.waitForTimeout(1200);
    await p2.reload();
    await p2.waitForTimeout(2500);
    const re = await p2.evaluate(() => ({
        kept: /^w = 0\.41/m.test(document.getElementById('custom_geom_text').value),
        rows: document.querySelectorAll('#custom-form .custom-rect-row').length,
        shapes: (document.getElementById('sim_canvas').layout?.shapes || []).length,
        err: document.getElementById('console_out').textContent.match(/ERROR[^\n]*/)?.[0] ?? '',
    }));
    check('a reload in the custom type keeps the geometry and draws it', re.kept && re.rows > 0 && re.shapes > 0 && re.err === '',
        JSON.stringify(re));
    await p2.close();
}

// ---- S-parameters of a symmetric geometry with unequal trace metals ----
{
    const p3 = await browser.newPage({ viewport: { width: 1400, height: 900 } });
    p3.on('dialog', d => d.accept());
    await p3.goto(URL, { waitUntil: 'networkidle' });
    await p3.selectOption('#tl_type', 'custom');
    await p3.waitForTimeout(600);
    const names = async (sigma) => {
        await p3.click('[data-tab="geometry"]');
        await p3.evaluate((sg) => {
            const t = document.getElementById('custom_geom_text');
            t.value = 'units mm\nbounds open open open gnd\ndiel x=-inf w=inf y=0 h=0.2 er=4.3 tand=0.02\n'
                + 'sig+ x=-0.4 w=0.3 y=0.2 h=0.035\nsig- x=0.1 w=0.3 y=0.2 h=0.035' + sg + '\n';
            t.dispatchEvent(new Event('input'));
        }, sigma);
        await p3.waitForTimeout(800);
        await p3.click('#btn_solve');
        await p3.waitForFunction(() => document.getElementById('btn_solve').textContent === 'Solve'
            && !document.getElementById('btn_solve').disabled, null, { timeout: 180000 });
        await p3.click('[data-tab="sparams"]');
        await p3.waitForTimeout(800);
        const single = await p3.evaluate(() => (document.getElementById('sparam-plot').data || []).filter(t => t.showlegend !== false).map(t => t.name.split(' ')[0]));
        await p3.check('#sparam-diff');
        await p3.waitForTimeout(600);
        const mixed = await p3.evaluate(() => (document.getElementById('sparam-plot').data || []).filter(t => t.showlegend !== false).map(t => t.name.split(' ')[0]));
        await p3.uncheck('#sparam-diff');
        return { single, mixed };
    };
    const same = await names('');
    const diff = await names(' sigma=1e7');
    check('equal traces plot the first column only', same.single.includes('S41') && !same.single.includes('S22') && !same.mixed.includes('SCD21'),
        same.single.join() + ' | ' + same.mixed.join());
    check('unequal trace metals add S22, S32, S42 and the mode conversion', ['S22', 'S32', 'S42'].every(n => diff.single.includes(n))
        && ['SDC11', 'SCD11', 'SCD21'].every(n => diff.mixed.includes(n)), diff.single.join() + ' | ' + diff.mixed.join());
    await p3.close();
}

// ---- sidebar: the Advanced section with every option enabled is not clipped ----
{
    const p4 = await browser.newPage({ viewport: { width: 1400, height: 800 } });
    await p4.goto(URL, { waitUntil: 'networkidle' });
    await p4.evaluate(() => document.querySelector('#advanced-params .section-header').click());
    await p4.waitForTimeout(400);
    for (const id of ['chk_solder_mask', 'chk_top_diel', 'chk_gnd_cut', 'chk_enclosure', 'chk_plating']) {
        await p4.evaluate((i) => { const c = document.getElementById(i); c.scrollIntoView(); if (!c.checked) c.click(); }, id);
        await p4.waitForTimeout(400);
    }
    const adv = await p4.evaluate(() => {
        const e = document.querySelector('#advanced-params .section-content');
        return { content: e.scrollHeight, shown: e.clientHeight };
    });
    check('the Advanced section shows all of its content with every option enabled',
        adv.content > 1000 && adv.shown >= adv.content, `${adv.shown} of ${adv.content} px`);
    await p4.close();
}

// ---- trapezoids and n-gons: form, preview, the quasi-static refusal, coax conversion ----
// Own page: the refused solve reports its error on the console.
{
    const p5 = await browser.newPage({ viewport: { width: 1400, height: 900 } });
    p5.on('dialog', d => d.accept());
    await p5.goto(URL, { waitUntil: 'networkidle' });
    await p5.selectOption('#tl_type', 'custom');
    await p5.waitForTimeout(400);
    await p5.selectOption('#custom-template', 'Differential microstrip with etched traces');
    await p5.waitForTimeout(600);
    const row = p5.locator('#custom-form .custom-rect-row.kind-sigp');
    const form = await p5.evaluate(() => {
        const r = document.querySelector('#custom-form .custom-rect-row.kind-sigp');
        return { shape: r.querySelector('select.custom-shape').value,
                 labels: [...r.querySelectorAll('.custom-cell-label')].map(l => l.textContent),
                 paths: (document.getElementById('sim_canvas').layout.shapes || []).filter(s => s.type === 'path').length };
    });
    check('a trapezoid row has the angle fields and draws as a polygon', form.shape === 'trap'
        && form.labels.includes('∠L') && form.labels.includes('∠R') && form.paths >= 2, JSON.stringify(form));
    await cellInput(row, /^∠R$/, p5).fill('10');
    await p5.waitForTimeout(500);
    check('the right angle field writes angle2', /angle=etch angle2=10/.test(await p5.inputValue('#custom_geom_text')));

    await p5.selectOption('#mesh_backend', 'rectilinear');
    await p5.click('#btn_solve');
    await p5.waitForFunction(() => /does not support/.test(document.getElementById('console_out').textContent), null, { timeout: 30000 });
    const refused = await p5.evaluate(() => document.getElementById('console_out').textContent);
    check('the quasi-static solver refuses the trapezoids', /other than plain rectangles.*line 8.*quasi-static solver does not support/.test(refused));

    await row.locator('select.custom-shape').selectOption('ngon');
    await p5.waitForTimeout(500);
    const ngonText = await p5.inputValue('#custom_geom_text');
    check('switching the shape to n-gon rewrites the line', /sig\+ {2}ngon {2}x=\S+ y=\S+ r=\S+ n=32 mirror=1/.test(ngonText),
        ngonText.split('\n').find(l => l.startsWith('sig+')));
    await p5.locator('#custom-form .custom-rect-row.kind-sigp select.custom-shape').selectOption('rect');
    await p5.waitForTimeout(500);
    check('and back to a rectangle', /sig\+ {2}x=\S+ w=\S+ y=\S+ h=\S+ mirror=1/.test(await p5.inputValue('#custom_geom_text')));

    await p5.locator('#custom-form .custom-rect-row.kind-sigp select.custom-shape').selectOption('ellipse');
    await p5.waitForTimeout(500);
    check('switching to an ellipse writes its semi-axes', /sig\+ {2}ellipse {2}x=\S+ y=\S+ rx=\S+ ry=\S+ n=32/.test(await p5.inputValue('#custom_geom_text')));

    // Corner radius and wall from the corners panel of a rectangle row.
    await p5.selectOption('#custom-template', 'Microstrip with rounded trace corners');
    await p5.waitForTimeout(600);
    const sig = p5.locator('#custom-form .custom-rect-row.kind-sigp');
    const btnText = await sig.locator('.custom-corners-toggle').textContent();
    check('the corners button names what is set', /corners r r↓/.test(btnText), btnText);
    await sig.locator('.custom-corners-toggle').click();
    await cellInput(sig.locator('.custom-corner-panel'), /^radius$/, p5).fill('5um');
    await p5.waitForTimeout(500);
    check('the radius field rewrites the line', /sig\+ {2}x=-w\/2 w=w y=h h=t radius=5um radius_bottom=0/.test(await p5.inputValue('#custom_geom_text')));
    await p5.selectOption('#custom-template', 'Twinax cable');
    await p5.waitForTimeout(600);
    const tw = await p5.evaluate(() => ({ shapes: [...document.querySelectorAll('#custom-form select.custom-shape')].map(s => s.value).join(),
        paths: (document.getElementById('sim_canvas').layout.shapes || []).filter(s => s.type === 'path').length,
        errors: document.getElementById('custom-geom-errors').textContent }));
    check('the twinax template: n-gons and a stadium shell', tw.shapes === 'ngon,ngon,rect' && tw.paths >= 5 && !tw.errors, JSON.stringify(tw));

    // Model Thick Plating applies to every conductor: the only plating control in the sidebar.
    const side = await p5.evaluate(() => {
        const vis = id => { const e = document.getElementById(id); return !!e && e.getClientRects().length > 0 && getComputedStyle(e).visibility !== 'hidden'; };
        return { thick: vis('chk_plating_thick_corners'), faces: vis('chk_plating_top'), sigma: vis('inp_plating_sigma') };
    });
    check('the sidebar keeps Model Thick Plating, not the other plating controls', side.thick && !side.faces && !side.sigma, JSON.stringify(side));

    // Coax: converts to n-gons and solves on the full-wave solver.
    await p5.selectOption('#tl_type', 'coax');
    await p5.waitForTimeout(500);
    await p5.click('#btn-convert-custom');
    await p5.waitForTimeout(800);
    const coax = await p5.evaluate(() => ({ type: document.getElementById('tl_type').value,
        text: document.getElementById('custom_geom_text').value,
        rows: [...document.querySelectorAll('#custom-form select.custom-shape')].map(s => s.value).join(),
        warnings: document.getElementById('custom-geom-warnings').textContent }));
    check('a coax converts to n-gons', coax.type === 'custom' && coax.rows === 'ngon,ngon,ngon' && /r_in=r2/.test(coax.text), coax.rows);
    check('the shielded coax has no open-boundary warning', !/open boundary/.test(coax.warnings), coax.warnings.slice(0, 80));
    await p5.click('#btn_solve');
    await p5.waitForFunction(() => document.getElementById('btn_solve').textContent === 'Solve'
        && !document.getElementById('btn_solve').disabled, null, { timeout: 300000 });
    const z = z0FromLog(await p5.evaluate(() => document.getElementById('console_out').textContent));
    const zRef = 376.730313668 / Math.sqrt(2.1) / (2 * Math.PI) * Math.log(2.95 / 0.92);
    check('the converted coax solves to the closed-form impedance', Math.abs(z - zRef) < 0.05, `Z0 ${z} vs ${zRef.toFixed(2)}`);
    await p5.close();
}

// First selection of the custom type in a fresh page: the template fills the text and
// the sweep tab lists its parameters right away.
{
    const p6 = await browser.newPage({ viewport: { width: 1400, height: 900 } });
    await p6.goto(URL, { waitUntil: 'networkidle' });
    await p6.selectOption('#tl_type', 'custom');
    await p6.waitForTimeout(500);
    const opts = await p6.evaluate(() => [...document.getElementById('sweep-x-selector').options].map(o => o.value));
    check('selecting the custom type lists the template parameters in the sweep tab',
        opts.filter(v => v.startsWith('cgp_')).length >= 3 && opts[0].startsWith('cgp_'), opts.join());
    // A shape change after a field edit of the same row starts from the edited box:
    // the trace x=-w/2 w=2 (w = 0.3) spans -0.15 ... 1.85, an n-gon centred on it.
    const sigRow = p6.locator('#custom-form .custom-rect-row.kind-sigp');
    await cellInput(sigRow, /^w$/, p6).fill('2');
    await p6.waitForTimeout(600);
    await sigRow.locator('select.custom-shape').selectOption('ngon');
    await p6.waitForTimeout(600);
    const ngon = (await p6.evaluate(() => document.getElementById('custom_geom_text').value))
        .split('\n').find(l => l.startsWith('sig+')) || '';
    check('a shape change after a field edit starts from the edited box', /\bngon\b/.test(ngon) && /\bx=0\.85\b/.test(ngon), ngon);
    await p6.close();
}

check('no console errors', errors.length === 0, errors.slice(0, 3).join(' | '));
await browser.close();
done();
