// Browser end-to-end for the Eye tab: open it, simulate a pulse (the tab runs its own
// channel sweep), then a clock steady state reusing that sweep; the plot, metrics and
// the clock eye render with no console errors. The Eye settings round-trip through a
// share link.
import { launch, URL, printErrors, finish } from './e2e_helpers.mjs';

const { browser, page, errors } = await launch({ clipboard: '__copiedURL' });
await page.goto(URL, { waitUntil: 'networkidle' });

await page.click('.tab-button[data-tab="eye"]');
const active = await page.$eval('#tab-eye', el => el.classList.contains('active'));
console.log('eye tab active:', active);

const runAndWait = async () => {
    await page.click('#btn-eye-run');
    await page.waitForFunction(() => /Done|failed/.test(document.getElementById('eye-status').textContent)
        || document.getElementById('eye-warnings').style.display === 'block' && !document.getElementById('eye-status').textContent,
        { timeout: 180000 });
    return page.$eval('#eye-status', el => el.textContent);
};

const s1 = await runAndWait();
console.log('pulse status:', s1);
const metrics1 = await page.$$eval("#eye-metrics .eye-metric-grid span", els => els.map(e => e.textContent));
console.log('pulse metrics:', metrics1.join(' | '));
const traces1 = await page.$eval('#eye-plot', el => (el.data || []).length);
console.log('pulse traces:', traces1);
const info1 = await page.$eval('#eye-metrics .eye-info', el => el.textContent).catch(() => '');
console.log('info:', info1);

await page.selectOption('#eye-stim', 'clock');
const s2 = await runAndWait();
console.log('clock status:', s2);
const info2 = await page.$eval('#eye-metrics .eye-info', el => el.textContent).catch(() => '');
console.log('clock info:', info2);
await page.selectOption('#eye-view', 'eye');
await page.waitForTimeout(500);
const eyeTraces = await page.$eval('#eye-plot', el => (el.data || []).length);
console.log('eye traces:', eyeTraces);
await page.selectOption('#eye-view', 'spectrum');
await page.waitForTimeout(500);
const specTraces = await page.$eval('#eye-plot', el => (el.data || []).length);
console.log('spectrum traces:', specTraces);
const warn = await page.$eval('#eye-warnings', el => el.style.display === 'block' ? el.textContent : '');
console.log('warnings:', warn);

// Bit pattern: PRBS7, the density eye and its metrics.
await page.selectOption('#eye-stim', 'pattern');
const s3 = await runAndWait();
await page.selectOption('#eye-view', 'eye');
await page.waitForTimeout(500);
const heat = await page.$eval('#eye-plot', el => (el.data || []).some(t => t.type === 'heatmap'));
const eyeMetrics = await page.$$eval('#eye-metrics .eye-metric-grid span', els => els.map(e => e.textContent));
const patternOk = /Done/.test(s3) && heat && eyeMetrics.includes('Eye height (centre)') && eyeMetrics.includes('Eye width (at threshold)');
console.log('PRBS7 eye:', patternOk, s3, eyeMetrics.slice(0, 12).join(' | '));
await page.selectOption('#eye-view', 'wave');
await page.selectOption('#eye-stim', 'clock');

// Share link: non-default Eye settings come back.
await page.fill('#eye-tr', '20 ps');
await page.selectOption('#eye-term', 'vtt');
await page.fill('#eye-vtt', '0.8 V');
await page.check('#eye-rx-ac');
await page.evaluate(() => document.getElementById('copy-link-btn').click());
await page.waitForTimeout(300);
const linkPage = await browser.newPage();
await linkPage.goto(await page.evaluate(() => window.__copiedURL), { waitUntil: 'networkidle' });
const restored = await linkPage.evaluate(() => ({
    tr: document.getElementById('eye-tr').value, term: document.getElementById('eye-term').value,
    vtt: document.getElementById('eye-vtt').value, rxAc: document.getElementById('eye-rx-ac').checked,
    stim: document.getElementById('eye-stim').value,
}));
const linkOk = restored.tr === '20 ps' && restored.term === 'vtt' && restored.vtt === '0.8 V' && restored.rxAc && restored.stim === 'clock';
console.log('share link restores the Eye settings:', linkOk, restored);

// Link contents: SI numbers, no inactive settings, nothing at all for the defaults.
const linkJson = async (p) => {
    await p.evaluate(() => { window.__copiedURL = null; document.getElementById('copy-link-btn').click(); });
    await p.waitForTimeout(300);
    const u = new globalThis.URL(await p.evaluate(() => window.__copiedURL));
    const enc = u.searchParams.get('params');
    return enc ? JSON.parse(decodeURIComponent(atob(enc))) : {};
};
const j1 = await linkJson(page);
await page.selectOption('#eye-term', 'r');
const j2 = await linkJson(page);
const fresh = await browser.newPage();
await fresh.addInitScript(() => Object.defineProperty(navigator, 'clipboard', { configurable: true,
    value: { writeText: (t) => { window.__copiedURL = t; return Promise.resolve(); } } }));
await fresh.goto(URL, { waitUntil: 'networkidle' });
const j0 = await linkJson(fresh);
const contentOk = j1.eye_tr === 2e-11 && j1.eye_vtt === 0.8 && !('eye_width' in j1) && !('eye_cycles' in j1)
    && !('eye_vtt' in j2) && Object.keys(j0).length === 0;
console.log('link contents:', contentOk, JSON.stringify(j1), '| after term=r:', JSON.stringify(j2), '| defaults:', JSON.stringify(j0));

// Pattern settings in a link: the bit rate as a number, the custom bits only when used.
await page.selectOption('#eye-stim', 'pattern');
await page.fill('#eye-bitrate', '16 Gb/s');
await page.selectOption('#eye-pattern', 'custom');
await page.fill('#eye-bits', '0101 1100');
const j3 = await linkJson(page);
await page.selectOption('#eye-pattern', 'prbs9');
const j4 = await linkJson(page);
const fresh2 = await browser.newPage();
await fresh2.goto(URL + '?params=' + btoa(encodeURIComponent(JSON.stringify(j3))), { waitUntil: 'networkidle' });
const back = await fresh2.evaluate(() => [document.getElementById('eye-bitrate').value, document.getElementById('eye-pattern').value, document.getElementById('eye-bits').value]);
const patLinkOk = j3.eye_bitrate === 16e9 && j3.eye_bits === '0101 1100' && !('eye_bits' in j4) && j4.eye_pattern === 'prbs9'
    && back[0] === '16 Gb/s' && back[1] === 'custom' && back[2] === '0101 1100';
console.log('pattern link:', patLinkOk, JSON.stringify(j3), '| prbs9:', JSON.stringify(j4), '| restored', back);

// A link from before the numeric form (unit strings) still loads.
const legacy = btoa(encodeURIComponent(JSON.stringify({ eye_length: '500 mm', eye_tr: '5 ps', eye_stim: 'step' })));
await fresh.goto(URL + '?params=' + legacy, { waitUntil: 'networkidle' });
const leg = await fresh.evaluate(() => [document.getElementById('eye-length').value, document.getElementById('eye-tr').value, document.getElementById('eye-stim').value]);
const legacyOk = leg[0] === '500 mm' && leg[1] === '5 ps' && leg[2] === 'step';
console.log('legacy string link loads:', legacyOk, leg);

// A waveguide has no time-domain simulation.
await page.selectOption('#tl_type', 'rect_waveguide');
await page.click('.tab-button[data-tab="eye"]');
const wgOff = await page.evaluate(() => document.getElementById('btn-eye-run').disabled
    && document.getElementById('eye-notice').style.display === 'block');
console.log('waveguide disables the Eye tab:', wgOff);

printErrors(errors);
const pass = linkOk && contentOk && legacyOk && wgOff && patternOk && patLinkOk && active && /Done/.test(s1) && /Done/.test(s2) && traces1 > 0 && metrics1.length > 4
    && eyeTraces > 0 && specTraces > 0 && /channel sweep/.test(info1) && /channel sweep/.test(info2) && errors.length === 0;
await finish(browser, pass);
