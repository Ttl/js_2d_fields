// Browser end-to-end: the UI must stay responsive while a solve runs.
//
// Before the worker migration every solve ran on the main thread, and a single
// uninterruptible WASM call blocked it for seconds at a time (measured: 7.3 s for one FDM
// certification LU, 7.5 s for one modes eigensolve) — a frozen tab and Firefox's "this
// tab is slowing down the browser" dialog. This test pins the fix: it drives a real solve
// and samples requestAnimationFrame throughout, then asserts the worst frame gap stays
// far below the ~5 s at which browsers start complaining.
//
// Needs the dev server on localhost:8731 (see tests/run.mjs e2e tier) and chromium.
import { chromium } from 'playwright-core';

const URL = 'http://localhost:8731/field_solver.html';
const browser = await chromium.launch({ executablePath: '/snap/bin/chromium', args: ['--no-sandbox'] });
const page = await browser.newPage();
const errors = [];
page.on('console', m => {
    if (m.type() === 'error' && !/Failed to load resource/i.test(m.text())) errors.push(m.text());
});
page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));

let failures = 0;
const check = (name, ok, detail = '') => {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
};

await page.goto(URL, { waitUntil: 'networkidle' });

// A dense-ish quasi-static solve: enough refinement that the pre-worker build would
// certify on a large bisected grid and block for seconds.
await page.evaluate(() => {
    document.getElementById('mesh_backend').value = 'rectilinear';
    document.getElementById('mesh_backend').dispatchEvent(new Event('change', { bubbles: true }));
    const cb = document.getElementById('chk_interp_sweep'); if (cb) cb.checked = false;
    document.getElementById('freq-points').value = '3';
    document.getElementById('freq-start').value = '1 GHz';
    document.getElementById('freq-stop').value = '20 GHz';
    document.getElementById('inp_tolerance').value = '0.1';   // percent — forces real refinement
    document.getElementById('inp_max_nodes').value = '60';
});

// Frame-gap probe. rAF only fires when the main thread is free, so the largest gap
// between consecutive callbacks IS the longest UI freeze the user would feel. Installed
// once and restarted before each phase — a probe that is not running records nothing, and
// a check against an empty sample would pass vacuously.
await page.evaluate(() => {
    window.__probeStart = () => {
        if (window.__probe) cancelAnimationFrame(window.__probe);
        window.__gaps = [];
        let last = performance.now();
        const tick = () => {
            const now = performance.now();
            window.__gaps.push(now - last);
            last = now;
            window.__probe = requestAnimationFrame(tick);
        };
        window.__probe = requestAnimationFrame(tick);
    };
    window.__probeStop = () => {
        cancelAnimationFrame(window.__probe);
        window.__probe = null;
        const g = window.__gaps.slice().sort((a, b) => b - a);
        return { maxGap: g[0] || 0, top: g.slice(0, 5), frames: window.__gaps.length };
    };
});
await page.evaluate(() => window.__probeStart());

const t0 = Date.now();
await page.click('#btn_solve');
// Wait for the solve to finish (button flips back from "Stop").
await page.waitForFunction(() => document.getElementById('btn_solve').textContent === 'Solve',
    null, { timeout: 300000 });
const wall = (Date.now() - t0) / 1000;

const { maxGap, top: gaps, frames } = await page.evaluate(() => window.__probeStop());

console.log(`solve wall time ${wall.toFixed(1)}s, ${frames} animation frames observed`);
console.log(`worst frame gaps: ${gaps.map(v => v.toFixed(0) + 'ms').join(', ')}`);

// Generous bound: browsers start warning around 5 s, and Plotly redraws on the main
// thread legitimately cost a few hundred ms. Anything under a second means the solve
// itself is no longer blocking.
check('main thread never blocks for ~1s during a solve', maxGap < 1000, `worst gap ${maxGap.toFixed(0)}ms`);
// A frozen main thread produces almost no frames; a live one produces many.
check('UI kept animating throughout the solve', frames > wall * 5,
    `${frames} frames over ${wall.toFixed(1)}s`);

const logText = await page.textContent('#console_out');
check('solve produced results', /RESULTS:/.test(logText), logText.slice(-120).replace(/\n/g, ' | '));

// The geometry tab must be able to plot the returned fields: they crossed a thread
// boundary as structured-clone copies and are grafted back onto the main-thread solver.
const fieldsOK = await page.evaluate(() => {
    const el = document.getElementById('sim_canvas');
    return !!(el && el.data && el.data.length);
});
check('field plot rendered from worker-returned arrays', fieldsOK);

// ---- Modes tab: the single worst blocking call in the old build ----
await page.click('[data-tab="modes"]');
await page.evaluate(() => {
    document.getElementById('modes-nev').value = '4';
    document.getElementById('modes-mesh-density').value = '8';
});
await page.evaluate(() => window.__probeStart());
await page.click('#btn-solve-modes');
await page.waitForFunction(() => !document.getElementById('btn-solve-modes').disabled,
    null, { timeout: 300000 });
const modesProbe = await page.evaluate(() => window.__probeStop());
const modesGap = modesProbe.maxGap;
const modesStatus = await page.textContent('#modes-status');
check('modes solve completed', /converged|No modes|failed/i.test(modesStatus), modesStatus);
check('modes solve kept the UI responsive', modesProbe.frames > 20 && modesGap < 1000,
    `worst gap ${modesGap.toFixed(0)}ms over ${modesProbe.frames} frames`);
// Selecting a mode fetches its field grid from the worker asynchronously.
const modePlotOK = await page.evaluate(async () => {
    const rows = document.querySelectorAll('#modes-list tr[data-idx]');
    if (!rows.length) return 'no-rows';
    rows[0].click();
    for (let i = 0; i < 100; i++) {
        await new Promise(r => setTimeout(r, 100));
        const el = document.getElementById('modes-plot');
        if (el && el.data && el.data.length) return 'ok';
    }
    return 'timeout';
});
check('mode field grid fetched from worker and plotted', modePlotOK === 'ok', modePlotOK);

// ---- Parameter sweep ----
await page.click('[data-tab="sweep"]');
await page.evaluate(() => {
    document.getElementById('sweep-points').value = '3';
});
await page.evaluate(() => window.__probeStart());
await page.click('#btn-run-sweep');
await page.waitForFunction(() => document.getElementById('btn-run-sweep').style.display !== 'none',
    null, { timeout: 300000 });
const sweepProbe = await page.evaluate(() => window.__probeStop());
const sweepGap = sweepProbe.maxGap;
const sweepLog = await page.textContent('#console_out');
check('parameter sweep completed', /Sweep complete: \d+ points/.test(sweepLog),
    (sweepLog.match(/Sweep complete[^\n]*/) || ['(no completion line)'])[0]);
check('parameter sweep kept the UI responsive', sweepProbe.frames > 20 && sweepGap < 1000,
    `worst gap ${sweepGap.toFixed(0)}ms over ${sweepProbe.frames} frames`);

check('no console/page errors', errors.length === 0, errors.slice(0, 2).join(' | '));

await browser.close();
console.log(failures === 0 ? '\nALL RESPONSIVENESS TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
