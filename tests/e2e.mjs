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

// Stop right after Solve, with the solve queued behind a plot job: a new plot frequency
// starts a plot job on the kept solve and Solve (after an input change) queues behind it.
// Stop targets the solve job by id, so it holds whenever it reaches the worker.
const stopStart = Date.now();
await page.evaluate(() => {
    const pf = document.getElementById('plot-freq');
    pf.value = '5 GHz';
    pf.dispatchEvent(new Event('change'));
    document.getElementById('inp_w').value = '0.36 mm';
    const btn = document.getElementById('btn_solve');
    btn.click();
    btn.click();
});
try {
    await page.waitForFunction(() => document.getElementById('btn_solve')?.textContent === 'Solve', null, { timeout: 150000 });
} catch {}
const stopOut = await page.$eval('#console_out', el => el.textContent).catch(() => '');
const stopOk = /Simulation stopped by user/.test(stopOut.slice(stopOut.lastIndexOf('Stop requested')));
console.log(`queued solve stopped: ${stopOk} (${Date.now() - stopStart} ms)`);
await page.evaluate(() => {
    document.getElementById('plot-freq').value = '';
    document.getElementById('inp_w').value = '0.35 mm';
});

// A DC-only solve: the summary reports the DC row (Zc infinite), not "below cutoff".
await page.evaluate(() => {
    document.getElementById('freq-points').value = '1';
    document.getElementById('freq-start').value = '0';
    document.getElementById('freq-stop').value = '0';
});
await page.click('#btn_solve');
try {
    await page.waitForFunction(() => document.getElementById('btn_solve')?.textContent === 'Solve', null, { timeout: 150000 });
} catch {}
await page.waitForTimeout(500);
const dcOut = await page.$eval('#console_out', el => el.textContent).catch(() => '');
const dcSummary = dcOut.slice(dcOut.lastIndexOf('RESULTS:'));
const dcOk = /Frequency: 0 Hz/.test(dcSummary) && /Z0: \d+\.\d+ Ohm/.test(dcSummary) && /Zc: ∞ Ohm/.test(dcSummary)
    && /eps_eff: \d+\.\d+/.test(dcSummary) && !/Below cutoff/.test(dcSummary);
console.log('DC-only summary:', dcOk ? 'ok' : dcSummary.slice(0, 200));

// An interpolating sweep from DC: the DC row is solved on its own and its ground-spreading
// note (the microstrip ground return spreads below about 1 MHz) reaches the log.
await page.evaluate(() => {
    document.getElementById('chk_interp_sweep').checked = true;
    document.getElementById('freq-points').value = '5';
    document.getElementById('freq-start').value = '0';
    document.getElementById('freq-stop').value = '10 GHz';
    document.getElementById('inp_w').value = '0.351 mm';
});
const spreadFrom = (await page.$eval('#console_out', el => el.textContent)).length;
await page.click('#btn_solve');
try {
    await page.waitForFunction(() => document.getElementById('btn_solve')?.textContent === 'Solve', null, { timeout: 150000 });
} catch {}
const spreadOk = /spreads sideways/.test((await page.$eval('#console_out', el => el.textContent)).slice(spreadFrom));
console.log('interpolating sweep from DC logs the ground-spreading note:', spreadOk);

await finish(browser, errors.length === 0 && solved && hasPlot >= 2 && stopOk && dcOk && spreadOk);
