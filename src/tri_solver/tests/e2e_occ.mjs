// Browser end-to-end for the OCC mesher option: load the app, select the
// "Full-wave (OCC mesh)" backend, solve, confirm WASM loads + results render + mesh
// overlay plots + no console errors.
import { chromium } from 'playwright-core';

const URL = 'http://localhost:8731/field_solver.html';
const browser = await chromium.launch({ executablePath: '/snap/bin/chromium', args: ['--no-sandbox'] });
const page = await browser.newPage();
const errors = [];
page.on('console', m => { if (m.type() === 'error' && !/Failed to load resource/i.test(m.text())) errors.push(m.text()); });
page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));
page.on('response', r => { if (r.status() === 404 && !/favicon/.test(r.url())) errors.push('404 ' + r.url()); });

await page.goto(URL, { waitUntil: 'networkidle' });
console.log('page loaded');
await page.selectOption('#mesh_backend', 'fullwave_mqs');
console.log('selected backend:', await page.inputValue('#mesh_backend'));

// enable solder mask (exercise the OCC conformity win) if the control exists
await page.evaluate(() => {
    const cb = document.getElementById('chk_interp_sweep'); if (cb) cb.checked = false;
    document.getElementById('freq-points').value = '3';
    document.getElementById('freq-start').value = '1 GHz';
    document.getElementById('freq-stop').value = '10 GHz';
});

const solveBtn = await page.$('#btn_solve') || await page.$('button:has-text("Solve")');
if (!solveBtn) { console.log('SOLVE BUTTON NOT FOUND'); await browser.close(); process.exit(1); }
await solveBtn.click();
console.log('clicked solve — waiting for results...');

// Wait for the solve to actually FINISH, signalled by the button flipping back out of
// "Stop" mode. The previous condition scanned the log for /Z0:/ or /ERROR/i — and
// /ERROR/i matches the "Energy error=..." in every refinement progress line, so it
// returned on the first pass and then read a half-finished log. That only ever passed
// because the solve froze the main thread and Playwright's poll could not run until it
// was over; once the solve moved to a worker the race became visible. Do not match log
// prose here: completion is a UI state, not a string.
try {
    await page.waitForFunction(
        () => document.getElementById('btn_solve')?.textContent === 'Solve',
        { timeout: 150000 });
} catch {}
await page.waitForTimeout(1000);

const hasPlot = await page.$$eval('.js-plotly-plot, .plotly', els => els.length);
const consoleOut = await page.$eval('#console_out', el => el.textContent).catch(() => '');
console.log('plotly plots:', hasPlot);
console.log('--- app console_out (tail) ---');
console.log(consoleOut.split('\n').slice(-16).map(l => '   ' + l).join('\n'));
console.log('--- console errors:', errors.length, '---');
errors.slice(0, 10).forEach(e => console.log('  ERR:', e.slice(0, 200)));

const m = consoleOut.match(/Z0:\s*([\d.]+)\s*Ohm/i);
const solved = !!m && parseFloat(m[1]) > 0 && !/ERROR:/.test(consoleOut);
console.log('parsed Z0:', m ? m[1] : '(none)');
const pass = errors.length === 0 && solved && hasPlot >= 2;
await browser.close();
console.log(pass ? '\nOCC E2E PASS' : '\nOCC E2E NEEDS REVIEW');
process.exit(pass ? 0 : 1);
