// Browser end-to-end: load the app, select the full-wave backend, solve, and confirm the
// WASM loads, the solve completes, the results and mesh plots render, and there are no
// console errors.
import { launch, URL, printErrors, finish } from './e2e_helpers.mjs';

const { browser, page, errors } = await launch();
await page.goto(URL, { waitUntil: 'networkidle' });
await page.selectOption('#mesh_backend', 'fullwave_mqs');
console.log('selected backend:', await page.inputValue('#mesh_backend'));

// Discrete sweep, few points.
await page.evaluate(() => {
    document.getElementById('chk_interp_sweep').checked = false;
    document.getElementById('freq-points').value = '3';
    document.getElementById('freq-start').value = '1 GHz';
    document.getElementById('freq-stop').value = '10 GHz';
});
await page.click('#btn_solve');
console.log('clicked solve, waiting for results...');

// Completion is the Solve button leaving "Stop" mode. Log text is no signal: /ERROR/i
// matches the "Energy error=..." of every refinement pass.
try {
    await page.waitForFunction(() => document.getElementById('btn_solve')?.textContent === 'Solve', null, { timeout: 150000 });
} catch {}
await page.waitForTimeout(1000);

const hasPlot = await page.$$eval('.js-plotly-plot, .plotly', els => els.length);
const consoleOut = await page.$eval('#console_out', el => el.textContent).catch(() => '');
console.log('plotly plots:', hasPlot);
console.log('--- app console_out (tail) ---');
console.log(consoleOut.split('\n').slice(-14).map(l => '   ' + l).join('\n'));
printErrors(errors);

const m = consoleOut.match(/Zc:\s*([\d.]+)/i);
const solved = !!m && parseFloat(m[1]) > 0 && !/ERROR:/.test(consoleOut);
console.log('parsed Zc:', m ? m[1] : '(none)');
await finish(browser, errors.length === 0 && solved && hasPlot >= 2);
