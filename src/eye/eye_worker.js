// Eye worker: runs the time-domain link simulation off the UI thread.
//
// main to worker: { id, type: 'run', table, cfg }
// worker to main: { id, type: 'done', result } or { id, type: 'error', message }
import { loadFFT } from './fft.js';
import { simulateLink } from './transient.js';

let fftPromise = null;

// Every typed array of the result, for a zero-copy transfer.
function buffers(obj, out = new Set()) {
    if (!obj || typeof obj !== 'object') return out;
    if (ArrayBuffer.isView(obj)) { out.add(obj.buffer); return out; }
    for (const v of Object.values(obj)) buffers(v, out);
    return out;
}

self.onmessage = async (e) => {
    const { id, type, table, cfg } = e.data;
    if (type !== 'run') { self.postMessage({ id, type: 'error', message: `unknown job "${type}"` }); return; }
    try {
        if (!fftPromise) fftPromise = loadFFT();
        const fft = await fftPromise;
        const t0 = Date.now();
        const result = simulateLink(table, cfg, fft);
        result.elapsed = Date.now() - t0;
        // Only what the page shows: the full-resolution arrays stay here.
        delete result.t; delete result.traces;
        self.postMessage({ id, type: 'done', result }, [...buffers(result)]);
    } catch (err) {
        self.postMessage({ id, type: 'error', message: (err && err.message) || String(err) });
    }
};
