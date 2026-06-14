// Browser end-to-end: load the app, select the triangular backend, solve, and
// confirm the WASM loads, the solve completes, and results render — no console errors.
import { chromium } from 'playwright-core';

const URL = 'http://localhost:8731/field_solver.html';
const browser = await chromium.launch({ executablePath: '/snap/bin/chromium', args: ['--no-sandbox'] });
const page = await browser.newPage();
const errors = [];
const logs = [];
page.on('console', m => {
    logs.push(m.text());
    // Ignore favicon / generic resource-load 404s (covered by the response handler).
    if (m.type() === 'error' && !/Failed to load resource/i.test(m.text())) errors.push(m.text());
});
page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));
page.on('response', r => { if (r.status() === 404 && !/favicon/.test(r.url())) errors.push('404 ' + r.url()); });

await page.goto(URL, { waitUntil: 'networkidle' });
console.log('page loaded');

// Select the full-wave (triangular) backend
await page.selectOption('#mesh_backend', 'fullwave_mqs');
console.log('selected full-wave backend:', await page.inputValue('#mesh_backend'));

// Isolate: discrete sweep, few points, to measure per-frequency solve cost.
await page.evaluate(() => {
    const cb = document.getElementById('chk_interp_sweep'); if (cb) cb.checked = false;
    document.getElementById('freq-points').value = '3';
    document.getElementById('freq-start').value = '1 GHz';
    document.getElementById('freq-stop').value = '10 GHz';
});

// Find and click the solve button
const solveBtn = await page.$('#btn_solve') || await page.$('button:has-text("Solve")');
if (!solveBtn) { console.log('SOLVE BUTTON NOT FOUND'); await browser.close(); process.exit(1); }
await solveBtn.click();
console.log('clicked solve — waiting for results...');

// Wait for a result value to appear (Z0 / impedance text), up to 90s (gmsh wasm + solve)
let ok = false;
try {
    await page.waitForFunction(() => {
        const t = document.body.innerText;
        return /Z0|Z₀|Impedance|Ω|ohm/i.test(t) && /\d/.test(t);
    }, { timeout: 90000 });
    ok = true;
} catch { ok = false; }

// Wait for the frequency sweep to finish (full-wave eigensolve per point is slow).
try {
    await page.waitForFunction(() => {
        const t = document.getElementById('console_out')?.textContent || '';
        return /exact solves|frequency sweep \(/i.test(t) || /ERROR/i.test(t);
    }, { timeout: 150000 });
} catch {}
await page.waitForTimeout(3000);

const hasPlot = await page.$$eval('.js-plotly-plot, .plotly', els => els.length);
const consoleOut = await page.$eval('#console_out', el => el.textContent).catch(() => '');
const noticeVisible = await page.$eval('#results-notice', el => el.offsetParent !== null).catch(() => false);

console.log('plotly plots:', hasPlot);
console.log('--- app console_out (tail) ---');
console.log(consoleOut.split('\n').slice(-14).map(l => '   ' + l).join('\n'));
console.log('results "stale/no-data" notice visible:', noticeVisible);
console.log('--- console errors:', errors.length, '---');
errors.slice(0, 10).forEach(e => console.log('  ERR:', e.slice(0, 200)));

// Success: a numeric Z0 result was rendered and no console errors occurred.
const m = consoleOut.match(/Z0:\s*([\d.]+)\s*Ohm/i);
const solved = !!m && parseFloat(m[1]) > 0 && !/ERROR:/.test(consoleOut);
console.log('parsed Z0:', m ? m[1] : '(none)');
const pass = errors.length === 0 && solved && hasPlot >= 2;
await browser.close();
console.log(pass ? '\nE2E PASS' : '\nE2E NEEDS REVIEW');
process.exit(pass ? 0 : 1);
