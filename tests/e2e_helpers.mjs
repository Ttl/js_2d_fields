// Shared browser setup for the e2e tests. Needs the dev server on localhost:8731 (see
// tests/run.mjs e2e tier) and chromium.
import { chromium } from 'playwright-core';

export const URL = 'http://localhost:8731/field_solver.html';

// Launch chromium and open a page that collects console errors, page errors and 404s.
// clipboard: name of a window property that receives navigator.clipboard.writeText
// text. Headless Chromium rejects the real write, which sends copySettingsLink() into
// its prompt() fallback, and a modal prompt blocks the page forever.
export async function launch({ viewport, clipboard, onConsole } = {}) {
    const browser = await chromium.launch({ executablePath: '/snap/bin/chromium', args: ['--no-sandbox'] });
    const page = await browser.newPage(viewport ? { viewport } : {});
    if (clipboard) {
        await page.addInitScript(name => {
            Object.defineProperty(navigator, 'clipboard', {
                configurable: true,
                value: { writeText: (t) => { window[name] = t; return Promise.resolve(); } },
            });
        }, clipboard);
    }
    const errors = [];
    page.on('console', m => {
        onConsole?.(m);
        if (m.type() === 'error' && !/Failed to load resource/i.test(m.text())) errors.push(m.text());
    });
    page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));
    page.on('response', r => { if (r.status() === 404 && !/favicon/.test(r.url())) errors.push('404 ' + r.url()); });
    return { browser, page, errors };
}

// Print the collected browser errors.
export function printErrors(errors) {
    console.log('--- console errors:', errors.length, '---');
    errors.slice(0, 10).forEach(e => console.log('  ERR:', e.slice(0, 200)));
}

// Close the browser, print the verdict and exit.
export async function finish(browser, pass) {
    await browser.close();
    console.log(pass ? '\n✓ PASS' : '\n✗ FAIL');
    process.exit(pass ? 0 : 1);
}
