// Real FFTs from the WASM module (Eigen FFT, kissfft backend).
//
// Half spectra are Float64Arrays of n/2 + 1 interleaved (re, im) pairs. n must be even;
// sizes with only the factors 2, 3 and 5 are fast (see smoothSize).
import createFFTModule from '../wasm_solver/fft.js';

let modulePromise = null;

function getModule() {
    if (!modulePromise) modulePromise = createFFTModule().catch((err) => { modulePromise = null; throw err; });
    return modulePromise;
}

// Copy into a fresh WASM allocation. HEAPF64 is re-read after every malloc because
// memory growth replaces the buffer.
function withBuffers(M, sizes, fn) {
    const ptrs = sizes.map(n => M._malloc(8 * n));
    try {
        return fn(ptrs);
    } finally {
        for (const p of ptrs) M._free(p);
    }
}

async function loadFFT() {
    const M = await getModule();
    return {
        // nbatch inverse transforms, spectra back to back. Returns nbatch * n samples,
        // with the 1/n factor applied.
        irfftBatch(n, nbatch, spectra) {
            const nh = n / 2 + 1;
            if (spectra.length !== 2 * nh * nbatch) throw new Error('irfftBatch: spectrum length mismatch');
            return withBuffers(M, [2 * nh * nbatch, n * nbatch], ([ps, po]) => {
                M.HEAPF64.set(spectra, ps / 8);
                if (M._irfft_batch(n, nbatch, ps, po) !== 0) throw new Error(`irfft: bad size ${n}`);
                return M.HEAPF64.slice(po / 8, po / 8 + n * nbatch);
            });
        },
        // Forward transform of n real samples, unscaled.
        rfft(x) {
            const n = x.length, nh = n / 2 + 1;
            return withBuffers(M, [n, 2 * nh], ([pi, ps]) => {
                M.HEAPF64.set(x, pi / 8);
                if (M._rfft(n, pi, ps) !== 0) throw new Error(`rfft: bad size ${n}`);
                return M.HEAPF64.slice(ps / 8, ps / 8 + 2 * nh);
            });
        },
    };
}

// Smallest even n >= m whose only prime factors are 2, 3 and 5, and a multiple of `step`
// (even). A step of 16 keeps radix-2 stages in large transforms: 2 * 3^12 runs at twice the
// time per point of its smooth neighbours.
function smoothSize(m, step = 2) {
    let n = Math.max(step, Math.ceil(m / step) * step);
    for (;; n += step) {
        let k = n;
        for (const p of [2, 3, 5]) while (k % p === 0) k /= p;
        if (k === 1) return n;
    }
}

export { loadFFT, smoothSize };
