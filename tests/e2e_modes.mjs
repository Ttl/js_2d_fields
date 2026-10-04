// Browser end-to-end for the Modes tab: load the app, open the Modes tab, solve the
// eigenmodes, confirm the mode list + field plot render, and clicking a mode replots —
// all with no console errors.
import { launch, URL, printErrors, finish } from './e2e_helpers.mjs';

const { browser, page, errors } = await launch();

await page.goto(URL, { waitUntil: 'networkidle' });
console.log('page loaded');

// Open the Modes tab
await page.click('.tab-button[data-tab="modes"]');
const tabVisible = await page.$eval('#tab-modes', el => el.classList.contains('active'));
console.log('modes tab active:', tabVisible);

await page.evaluate(() => {
    document.getElementById('modes-freq').value = '20 GHz';
    document.getElementById('modes-nev').value = '5';
});

await page.click('#btn-solve-modes');
console.log('clicked Solve Modes — waiting for the mode list...');

let ok = false;
try {
    await page.waitForFunction(() => document.querySelectorAll('#modes-list tr[data-idx]').length > 0,
        { timeout: 120000 });
    ok = true;
} catch { ok = false; }

await page.waitForTimeout(2000);

const nRows = await page.$$eval('#modes-list tr[data-idx]', els => els.length);
const status = await page.$eval('#modes-status', el => el.textContent).catch(() => '');
const hasPlot = await page.$$eval('#modes-plot.js-plotly-plot, #modes-plot .plotly', els => els.length);
const nSelected = await page.$$eval('#modes-list tr.selected', els => els.length);
console.log('mode rows:', nRows);
console.log('status:', status);
console.log('plotly plot in modes-plot:', hasPlot);
console.log('auto-selected rows:', nSelected);

// Click the last row and confirm the selection moves (replot path runs).
let clickOk = false;
if (nRows > 1) {
    await page.click('#modes-list tr[data-idx="' + (nRows - 1) + '"]');
    await page.waitForTimeout(1000);
    clickOk = await page.$eval('#modes-list tr[data-idx="' + (nRows - 1) + '"]',
        el => el.classList.contains('selected'));
}
console.log('clicking a mode selects it:', clickOk);

// A row clicked right before a new modes solve: its field fetch still runs against the
// old solve and must not end up in the new solve's cache. The new solve is of a wider
// trace, so a stale field shows as a mesh different from the other rows'.
let cacheOk = nRows < 3;
if (nRows >= 3) {
    await page.evaluate(() => {
        document.getElementById('inp_w').value = '0.6 mm';
        document.querySelector('#modes-list tr[data-idx="1"]').click();
        document.getElementById('btn-solve-modes').click();
    });
    await page.waitForFunction(() => !document.getElementById('btn-solve-modes').disabled
        && document.querySelectorAll('#modes-list tr[data-idx]').length >= 3, null, { timeout: 120000 });
    const meshOf = async idx => {
        await page.click(`#modes-list tr[data-idx="${idx}"]`);
        await page.waitForTimeout(1000);
        return page.evaluate(() => {
            const m = document.getElementById('modes-plot')._modesMesh;
            return m ? m.blocks[0].tris.length / 6 : -1;
        });
    };
    const n1 = await meshOf(1), n2 = await meshOf(2);
    cacheOk = n1 > 0 && n1 === n2;
    console.log(`row clicked before a new solve shows the new solve's field: ${cacheOk} (${n1} vs ${n2} triangles)`);
}

printErrors(errors);
await finish(browser, ok && nRows > 0 && hasPlot >= 1 && nSelected === 1 && (nRows <= 1 || clickOk) && cacheOk && errors.length === 0);
