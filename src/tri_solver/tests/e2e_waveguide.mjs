// Browser end-to-end: the rectangular waveguide media type through the real UI.
//
// This type touches more UI wiring than any other addition — the Solver dropdown lock, the
// interpolating-sweep lock, the S-parameter reference swap, the plating-face hiding, the
// missing "Potential" view — and none of it is reachable from the node tests. It also
// exercises the two behaviours that only show up in a sweep across cutoff: NaN gaps in the
// Results traces, and below-cutoff points dropped from the S-parameters.
//
// Needs the dev server on localhost:8731 (see tests/run.mjs e2e tier).
import { chromium } from 'playwright-core';

const URL = 'http://localhost:8731/field_solver.html';
const browser = await chromium.launch({ executablePath: '/snap/bin/chromium', args: ['--no-sandbox'] });
const page = await browser.newPage();
// The copy-link round-trip below needs the URL the button builds. Capture it by stubbing
// clipboard.writeText rather than granting clipboard permissions: headless Chromium
// rejects the real write, which sends copySettingsLink() into its prompt() fallback — and
// a modal prompt blocks the page forever. This still exercises the real settingsToURL path,
// which is the part under test.
await page.addInitScript(() => {
    Object.defineProperty(navigator, 'clipboard', {
        configurable: true,
        value: { writeText: (t) => { window.__copiedURL = t; return Promise.resolve(); } },
    });
});
const errors = [];
page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource/i.test(m.text())) errors.push(m.text());
});
page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));
page.on('response', r => { if (r.status() === 404 && !/favicon/.test(r.url())) errors.push('404 ' + r.url()); });

let failures = 0;
function check(name, cond, detail = '') {
    console.log(`${cond ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!cond) failures++;
}

await page.goto(URL, { waitUntil: 'networkidle' });
await page.selectOption('#tl_type', 'rect_waveguide');
await page.waitForTimeout(200);

// ---- the type switch reconfigures the sidebar ----
const ui = await page.evaluate(() => ({
    backend: document.getElementById('mesh_backend').value,
    rectDisabled: document.querySelector('#mesh_backend option[value="rectilinear"]').disabled,
    wgVisible: document.getElementById('wg-params').style.display,
    msVisible: document.getElementById('microstrip-params').style.display,
    interpChecked: document.getElementById('chk_interp_sweep').checked,
    interpDisabled: document.getElementById('chk_interp_sweep').disabled,
    tolGroup: document.getElementById('interp-tolerance-group').style.display,
    zRefBox: document.getElementById('sparam-z-ref-label').style.display,
    zRefModal: document.getElementById('sparam-z-ref-modal').style.display,
    platingFaces: document.querySelector('#plating-params .control-group:has(#chk_plating_top)').style.display,
    sweepOpts: [...document.getElementById('sweep-x-selector').options].map(o => o.value),
}));
check('solver locks to full-wave', ui.backend === 'fullwave_mqs' && ui.rectDisabled === true);
check('waveguide params replace the microstrip block',
    ui.wgVisible === 'block' && ui.msVisible === 'none');
check('interpolating sweep is off and locked',
    ui.interpChecked === false && ui.interpDisabled === true && ui.tolGroup === 'none');
check('S-parameter reference shows the modal label, not the box',
    ui.zRefBox === 'none' && ui.zRefModal !== 'none');
check('per-face plating selection is hidden (one continuous wall)', ui.platingFaces === 'none');
check('sweep parameters are the waveguide set + roughness',
    ui.sweepOpts.join(',') === 'rq,wg_a,wg_b,wg_er,wg_tand,wg_sigma', ui.sweepOpts.join(','));

// ---- solve a sweep that straddles the 6.557 GHz WR-90 cutoff ----
await page.evaluate(() => {
    document.getElementById('freq-start').value = '5 GHz';
    document.getElementById('freq-stop').value = '13 GHz';
    document.getElementById('freq-points').value = '17';
});
await page.click('#btn_solve');
try {
    await page.waitForFunction(
        () => document.getElementById('btn_solve').textContent.trim() === 'Solve',
        { timeout: 180000 });
} catch {
    check('solve completed', false, 'timed out');
}
await page.waitForTimeout(1500);

const out = await page.evaluate(() => document.getElementById('console_out').textContent);
check('logs the fundamental-mode-only limitation and the single-mode band',
    /fundamental mode only \(cutoff 6\.55\d GHz, single-mode up to 13\.11\d GHz\)/.test(out));
check('warns about the below-cutoff points', /below the 6\.55\d GHz cutoff/.test(out));
check('summary quotes the first PROPAGATING point, not NaN',
    /Z0: \d+\.\d+ Ohm \(at 7\.00 GHz/.test(out) && !/Z0: NaN/.test(out), out.match(/Z0:.*/)?.[0]);

// ---- Results: the trace breaks at cutoff ----
await page.click('[data-tab="results"]');
await page.waitForTimeout(800);
const res = await page.evaluate(() => {
    const d = document.getElementById('results-plot').data;
    const y = d[0].y;
    return { points: y.length, nan: y.filter(v => v == null || Number.isNaN(v)).length,
             last: y[y.length - 1] };
});
// 5.0, 5.5, 6.0, 6.5 GHz are below the 6.557 GHz cutoff; 7.0 GHz up is propagating.
check('Results has all 17 points with the 4 below-cutoff ones blank',
    res.points === 17 && res.nan === 4, `${res.nan}/${res.points} blank`);
check('impedance falls toward 2(b/a)*eta0 well above cutoff',
    res.last > 330 && res.last < 450, res.last && res.last.toFixed(1));

// ---- S-parameters: self-referenced, below-cutoff rows dropped ----
await page.click('[data-tab="sparams"]');
await page.waitForTimeout(800);
const sp = await page.evaluate(() => {
    const d = document.getElementById('sparam-plot').data;
    return { points: d[0].y.length, s11max: Math.max(...d[0].y), s21max: Math.max(...d[1].y) };
});
check('S-parameters drop the below-cutoff points', sp.points === 13, `${sp.points} points`);
// Referenced to its own (complex) modal impedance a uniform guide is matched by
// definition, so this must be the -300 dB log-of-zero floor at EVERY point — not merely
// small. A real-only reference would leave ~-80 dB of spurious reflection here.
check('S11 is zero at every point (self-referenced)', sp.s11max === -300,
    `max ${sp.s11max.toFixed(1)} dB`);
check('S21 is a small insertion loss', sp.s21max < 0 && sp.s21max > -1, `max ${sp.s21max.toFixed(4)} dB`);

// ---- Geometry view has no Potential button (there is no static potential) ----
await page.click('[data-tab="geometry"]');
await page.waitForTimeout(800);
const geo = await page.evaluate(() => {
    const l = document.getElementById('sim_canvas').layout;
    const shapes = l.shapes || [];
    return {
        buttons: l.updatemenus ? l.updatemenus[0].buttons.map(b => b.label) : null,
        x: l.xaxis.range, y: l.yaxis.range,
        rings: shapes.filter(s => s.type === 'path' && s.fillrule === 'evenodd').length,
    };
});
check('geometry view offers Geometry and |E| Field only',
    JSON.stringify(geo.buttons) === '["Geometry","|E| Field"]', JSON.stringify(geo.buttons));
// The walls are not in solver.conductors (they live in the boundary conditions), so they
// are drawn from the cosmetic enclosure_walls descriptor as a ring. Without it the view is
// a bare rectangle of dielectric with nothing marking the metal.
check('the enclosing walls are drawn', geo.rings === 1, `${geo.rings} ring shape(s)`);
// WR-90 spans x = ±11.43 mm and y = 0..10.16 mm, so a centred view is symmetric about
// x = 0 and about y = 5.08, and reaches past the walls on every side.
check('the guide is centred in the view',
    Math.abs(geo.x[0] + geo.x[1]) < 0.05 &&
    Math.abs((geo.y[0] + geo.y[1]) / 2 - 5.08) < 1.5,
    `x ${geo.x.map(v => v.toFixed(2))}, y ${geo.y.map(v => v.toFixed(2))}`);
check('the view frames the whole cross-section including the walls',
    geo.x[1] > 11.43 && geo.y[1] > 10.16 && geo.y[0] < 0);

// ---- copy-link round trip ----
// The wg_ prefix has to be in TYPE_ONLY_KEYS and EXCLUDED_BY_TYPE for this to survive, and
// tl_type is restored before mesh_backend, so the backend lock has to re-run afterwards.
{
    // The Copy link button lives in the Results tab's bottom controls.
    await page.click('[data-tab="results"]');
    await page.waitForTimeout(300);
    // One field per evaluate, with a beat between. Each 'input' triggers a geometry
    // rebuild and redraw, and firing four of them inside ONE synchronous evaluate wedges
    // the page — verified to happen for microstrip too, so it is a property of driving the
    // app synchronously from a test, not of this media type. A real user cannot produce
    // that, and this loop matches what they can.
    for (const [id, v] of [['inp_wg_a', '19.05 mm'],   // WR-75-ish, non-default
                           ['inp_wg_b', '9.525 mm'],
                           ['inp_wg_er', '2.2'],
                           ['inp_rq', '0.8 um']]) {
        await page.evaluate(([id, v]) => {
            const el = document.getElementById(id);
            el.value = v;
            el.dispatchEvent(new Event('input', { bubbles: true }));
        }, [id, v]);
        await page.waitForTimeout(150);
    }
    await page.click('#copy-link-btn');
    await page.waitForTimeout(300);
    const url = await page.evaluate(() => window.__copiedURL || '');
    check('copy link produced a params URL', /\?params=/.test(url), url.slice(0, 60));
    if (!/\?params=/.test(url)) { await browser.close(); process.exit(1); }

    await page.goto(url, { waitUntil: 'networkidle' });
    await page.waitForTimeout(500);
    const restored = await page.evaluate(() => ({
        type: document.getElementById('tl_type').value,
        backend: document.getElementById('mesh_backend').value,
        a: document.getElementById('inp_wg_a').value,
        b: document.getElementById('inp_wg_b').value,
        er: document.getElementById('inp_wg_er').value,
        rq: document.getElementById('inp_rq').value,
        interpDisabled: document.getElementById('chk_interp_sweep').disabled,
    }));
    check('link restores the waveguide type and stays full-wave',
        restored.type === 'rect_waveguide' && restored.backend === 'fullwave_mqs'
        && restored.interpDisabled === true);
    check('link restores every wg_ parameter',
        parseFloat(restored.a) === 19.05 && parseFloat(restored.b) === 9.525
        && parseFloat(restored.er) === 2.2, JSON.stringify(restored));
    check('link restores the shared surface roughness', parseFloat(restored.rq) === 0.8, restored.rq);
}

check('no console errors', errors.length === 0, errors.join(' | '));
await browser.close();
console.log(failures === 0 ? '\nWAVEGUIDE E2E OK' : `\nWAVEGUIDE E2E: ${failures} FAILURE(S)`);
process.exit(failures === 0 ? 0 : 1);
