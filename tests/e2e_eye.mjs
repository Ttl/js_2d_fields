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

// Interface presets on a differential line: LVDS (current mode) and PCIe Gen3 (TX preset,
// CTLE and DFE chosen automatically), the eye after equalization and the pulse response.
await page.click('.tab-button[data-tab="geometry"]');
await page.selectOption('#tl_type', 'diff_microstrip');
await page.click('.tab-button[data-tab="eye"]');
const presetRun = async (preset) => {
    await page.selectOption('#eye-preset', preset);
    const st = await runAndWait();
    const m = await page.$$eval('#eye-metrics .eye-metric-grid > div', els => Object.fromEntries(els.map(e => [e.children[0].textContent, e.children[1].textContent])));
    return { st, m };
};
const lvds = await presetRun('lvds');
// Current drive: the channel response is relative to its low-frequency level (not V/A).
await page.selectOption('#eye-view', 'spectrum');
await page.waitForTimeout(300);
const lvdsH = await page.$eval('#eye-plot', el => { const t = el.data[0]; return { name: t.name, y0: t.y[0], max: Math.max(...t.y) }; });
const lvdsHOk = /relative to low frequency/.test(lvdsH.name) && Math.abs(lvdsH.y0) < 0.5 && lvdsH.max < 3;
console.log('LVDS channel response near 0 dB:', lvdsHOk, lvdsH);
await page.selectOption('#eye-view', 'wave');
const lvdsPkg = await page.evaluate(() => ['eye-pkgtx', 'eye-pkgrx', 'eye-padrx'].map(id => document.getElementById(id).value));
const lvdsOk = /Done/.test(lvds.st) && /mV/.test(lvds.m['Eye height (centre)'] || '') && JSON.stringify(lvdsPkg) === JSON.stringify(['rl', 'rl', '4.5 pF']);
console.log('LVDS preset:', lvdsOk, lvds.st, lvds.m['Eye height (centre)'], lvds.m['Eye amplitude']);
const pcie = await presetRun('pcie3');
await page.selectOption('#eye-view', 'pulse');
await page.waitForTimeout(500);
const pulseTraces = await page.$eval('#eye-plot', el => (el.data || []).map(t => t.name));
const pcieOk = /Done/.test(pcie.st) && /^P\d+: (no TX equalization|[\d.]+ dB boost)/.test(pcie.m['TX FIR (chosen)'] || '') && /dB/.test(pcie.m['CTLE DC gain (chosen)'] || '')
    && 'DFE taps' in pcie.m && pulseTraces.includes('Pulse response') && pulseTraces.includes('DFE taps');
console.log('PCIe Gen3 preset:', pcieOk, pcie.st, pcie.m['TX FIR (chosen)'], pcie.m['CTLE DC gain (chosen)'], pcie.m['Eye height (centre)'], pulseTraces);
await page.selectOption('#eye-view', 'wave');
const gen3Pkg = await page.evaluate(() => ['eye-pkgtx', 'eye-pkgrx', 'eye-padtx', 'eye-padrx'].map(id => document.getElementById(id).value));
const gen3PkgOk = JSON.stringify(gen3Pkg) === JSON.stringify(['line', 'line', '1 pF', '0.8 pF']);
console.log('PCIe Gen3 packages (8 GT/s):', gen3PkgOk, gen3Pkg);

// PCIe Gen5: 32 GT/s root port package at TX.
const gen5 = await presetRun('pcie5');
const gen5Pkg = await page.evaluate(() => ['eye-pkgtx-z0', 'eye-padtx', 'eye-pkgrx-z0'].map(id => document.getElementById(id).value));
const gen5Ok = /Done/.test(gen5.st) && JSON.stringify(gen5Pkg) === JSON.stringify(['38.3', '0.44 pF', '41.3']);
console.log('PCIe Gen5 preset with packages:', gen5Ok, gen5.st, gen5.m['Eye height (centre)'], gen5Pkg);

// PCIe Gen4: the behavioral packages are filled in, simulated and kept in a link.
const gen4 = await presetRun('pcie4');
const jPkg = await linkJson(page);
const pkgPage = await browser.newPage();
await pkgPage.goto(URL + '?params=' + btoa(encodeURIComponent(JSON.stringify(jPkg))), { waitUntil: 'networkidle' });
const pkgBack = await pkgPage.evaluate(() => ['eye-pkgtx', 'eye-pkgtx-z0', 'eye-pkgrx-delay', 'eye-padrx'].map(id => document.getElementById(id).value));
const gen4Ok = /Done/.test(gen4.st) && /mV/.test(gen4.m['Eye height (centre)'] || '') && jPkg.eye_pkgtx === 'line' && jPkg.eye_pkgrx === 'line'
    && !('eye_pkgtx_l' in jPkg) && pkgBack[0] === 'line' && pkgBack[1] === '35.8' && pkgBack[2] === '135 ps' && pkgBack[3] === '0.5 pF';
// Package S-parameter view: drawn from the inputs, redrawn while one is edited.
await page.selectOption('#eye-view', 'package');
await page.waitForTimeout(300);
const pkgTraces = await page.$eval('#eye-plot', el => el.data.map(t => t.name));
const s21At = () => page.$eval('#eye-plot', el => { const t = el.data.find(d => d.name === 'TX SDD21'); return t.y[t.y.length >> 1]; });
const before = await s21At();
await page.fill('#eye-pkgtx-lossc', '3');
await page.waitForTimeout(300);
const after = await s21At();
await page.screenshot({ path: '/tmp/eye_package_view.png' });
const pkgViewOk = ['TX SDD21', 'TX SDD11 (pin)', 'RX SDD21', 'RX SDD11 (pin)'].every(n => pkgTraces.includes(n)) && after < before - 1;
console.log('package view, live update:', pkgViewOk, pkgTraces, `${before.toFixed(2)} -> ${after.toFixed(2)} dB`);
// R + L template: the IBIS values and C_comp as the pad capacitance, the plot follows.
await page.selectOption('#eye-pkgrx-fill', 'fpga18');
await page.waitForTimeout(300);
const fpga = await page.evaluate(() => ['eye-pkgrx', 'eye-pkgrx-cpin', 'eye-pkgrx-r', 'eye-pkgrx-l', 'eye-padrx', 'eye-pkgrx-fill']
    .map(id => document.getElementById(id).value));
const fpgaRows = await page.evaluate(() => ['l', 'z0'].map(k => document.getElementById(`eye-row-pkgrx-${k}`).style.display !== 'none'));
const fpgaOk = JSON.stringify(fpga) === JSON.stringify(['rl', '1 pF', '0.3', '5 nH', '4.5 pF', '']) && fpgaRows[0] && !fpgaRows[1];
console.log('1.8 V FPGA I/O template:', fpgaOk, fpga, fpgaRows);
await page.selectOption('#eye-view', 'wave');
// Table 8-14 transmitter: jitter, Gaussian 10-90 % edges and 804 mV swing, in the link too.
const gen4JitOk = /ps rms \+ .* ps pp/.test(gen4.m['TX jitter per edge'] || '') && 'TX jitter of this run' in gen4.m
    && Math.abs(jPkg.eye_rj - 0.71e-12) < 1e-18 && Math.abs(jPkg.eye_pwdj - 5e-12) < 1e-18 && jPkg.eye_edge_shape === 'gauss'
    && Math.abs(jPkg.eye_tr - 21.875e-12) < 1e-18 && Math.abs(jPkg.eye_vhigh - 0.804) < 1e-12;
console.log('PCIe Gen4 TX jitter and edges:', gen4JitOk, gen4.m['TX jitter per edge'], gen4.m['TX jitter of this run'], gen4.m['Crossing jitter pk-pk / rms']);
console.log('PCIe Gen4 preset with packages:', gen4Ok, gen4.st, gen4.m['Eye height (centre)'], JSON.stringify(jPkg), pkgBack);

// LVCMOS presets: the 1.8 V FPGA I/O package at both ends, on a single-ended line.
await page.click('.tab-button[data-tab="geometry"]');
await page.selectOption('#tl_type', 'microstrip');
await page.click('.tab-button[data-tab="eye"]');
// LVPECL: terminated on the board at the RX pin; the location row only with an RX package.
const pecl = await presetRun('lvpecl');
const vis = () => page.$eval('#eye-row-termloc', el => el.style.display !== 'none');
const peclLoc = await page.$eval('#eye-termloc', el => el.value), peclVis = await vis();
const jPecl = await linkJson(page);
await page.selectOption('#eye-pkgrx', 'none');
const hidNoPkg = !(await vis());
const jNoPkg = await linkJson(page);
const jitShown = await page.$eval('#eye-jitter-set', el => el.style.display !== 'none');
const peclOk = jitShown && !('eye_rj' in jPecl) && /Done/.test(pecl.st) && peclLoc === 'board' && peclVis && hidNoPkg && jPecl.eye_termloc === 'board' && !('eye_termloc' in jNoPkg);
console.log('LVPECL board termination:', peclOk, pecl.st, peclLoc, peclVis, hidNoPkg, jPecl.eye_termloc, 'eye_termloc' in jNoPkg);
// A steady clock with jitter: TIE and period jitter; a burst hides the jitter inputs.
await page.fill('#eye-rj', '5 ps');
const peclJ = await runAndWait();
const peclJm = await page.$$eval('#eye-metrics .eye-metric-grid > div', els => Object.fromEntries(els.map(e => [e.children[0].textContent, e.children[1].textContent])));
await page.selectOption('#eye-clock-mode', 'burst');
const burstHidden = await page.$eval('#eye-jitter-set', el => el.style.display === 'none');
const jBurst = await linkJson(page);
await page.selectOption('#eye-clock-mode', 'steady');
const clockJitOk = /Done/.test(peclJ) && /ps \/ .*ps/.test(peclJm['Rising-edge TIE pk-pk / rms'] || '') && 'Period jitter pk-pk / rms' in peclJm && burstHidden && !('eye_rj' in jBurst);
console.log('LVPECL clock with jitter:', clockJitOk, peclJ, peclJm['Rising-edge TIE pk-pk / rms'], peclJm['Period jitter pk-pk / rms'], burstHidden);
const cmos = await presetRun('lvcmos33');
const cmosPkg = await page.evaluate(() => ['eye-pkgtx', 'eye-pkgrx', 'eye-pkgtx-l', 'eye-padtx', 'eye-padrx'].map(id => document.getElementById(id).value));
const cmosOk = /Done/.test(cmos.st) && JSON.stringify(cmosPkg) === JSON.stringify(['rl', 'rl', '5 nH', '4.5 pF', '4.5 pF']);
console.log('LVCMOS 3.3 V preset with packages:', cmosOk, cmos.st, cmosPkg);

// A waveguide has no time-domain simulation.
await page.selectOption('#tl_type', 'rect_waveguide');
await page.click('.tab-button[data-tab="eye"]');
const wgOff = await page.evaluate(() => document.getElementById('btn-eye-run').disabled
    && document.getElementById('eye-notice').style.display === 'block');
console.log('waveguide disables the Eye tab:', wgOff);

printErrors(errors);
const pass = clockJitOk && gen4JitOk && peclOk && lvdsHOk && gen3PkgOk && gen5Ok && cmosOk && gen4Ok && pkgViewOk && fpgaOk && linkOk && contentOk && legacyOk && wgOff && patternOk && patLinkOk && lvdsOk && pcieOk && active && /Done/.test(s1) && /Done/.test(s2) && traces1 > 0 && metrics1.length > 4
    && eyeTraces > 0 && specTraces > 0 && /channel sweep/.test(info1) && /channel sweep/.test(info2) && errors.length === 0;
await finish(browser, pass);
