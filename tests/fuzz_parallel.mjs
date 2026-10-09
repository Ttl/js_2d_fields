// Parallel QS vs full-wave fuzzer: the cases of tests/fuzz_qs_vs_fullwave.js spread over
// worker threads. Same arguments, same cases (drawn once, here), same report: the
// per-case lines print in case order and the summary is printSummary's, so the output
// matches a serial run of the same N and seed. ONLY= works as there.
//
//   node tests/fuzz_parallel.mjs [-j jobs] [N] [seed] [threshold%] [R-threshold%]
//
// Cases go to the workers one at a time, the next to whichever finishes first, since
// the solve time of a case varies by orders of magnitude. Default jobs: cpus - 2, at
// most one per 2 GB of memory (a worker holds its own WASM solvers). A worker that dies
// mid-case (WASM abort, out of memory) is replaced and the case reported as a
// one-sided failure. Progress goes to stderr when it is a terminal.
import { Worker, isMainThread, parentPort, workerData } from 'worker_threads';
import os from 'os';
import { parseArgs, printHeader, drawCases, evaluateCase, printSummary } from './fuzz_qs_vs_fullwave.js';

if (isMainThread) {
    await runMain();
} else {
    parentPort.on('message', async ({ i, spec }) => {
        parentPort.postMessage(await evaluateCase(i, spec, workerData.cfg));
    });
    parentPort.postMessage(null);   // ready
}

async function runMain() {
    const argv = process.argv.slice(2);
    let jobs = Math.max(1, Math.min(os.cpus().length - 2, Math.floor(os.totalmem() / 2 ** 31)));
    const jIdx = argv.findIndex(a => a === '-j' || a.startsWith('-j'));
    if (jIdx >= 0) {
        const v = argv[jIdx] === '-j' ? argv.splice(jIdx, 2)[1] : argv.splice(jIdx, 1)[0].slice(2);
        jobs = Math.max(1, parseInt(v) || 1);
    }
    const cfg = parseArgs(argv);
    printHeader(cfg);
    const cases = drawCases(cfg);
    jobs = Math.min(jobs, cases.length) || 1;
    console.error(`${cases.length} cases on ${jobs} worker threads`);

    const results = new Array(cases.length);
    const pos = new Map(cases.map((c, k) => [c.i, k]));
    let next = 0, printed = 0, done = 0;
    const t0 = Date.now();
    const tty = process.stderr.isTTY;
    const progress = () => {
        if (!tty) return;
        const el = (Date.now() - t0) / 1000;
        const eta = done ? el / done * (cases.length - done) : 0;
        const fmt = s => s >= 3600 ? `${(s / 3600).toFixed(1)}h` : s >= 60 ? `${(s / 60).toFixed(0)}m` : `${s.toFixed(0)}s`;
        process.stderr.write(`\r\x1b[K${done}/${cases.length} done, ${fmt(el)} elapsed, ETA ${fmt(eta)}`);
    };
    // Lines print in case order: everything up to the first unfinished case.
    const flush = () => {
        let out = '';
        while (printed < cases.length && results[printed]) {
            for (const line of results[printed].lines) out += line + '\n';
            printed++;
        }
        if (out) {
            if (tty) process.stderr.write('\r\x1b[K');
            process.stdout.write(out);
        }
        progress();
    };
    const record = r => { results[pos.get(r.i)] = r; done++; flush(); };

    await new Promise(resolve => {
        let live = 0;
        const spawn = () => {
            const w = new Worker(new URL(import.meta.url), { workerData: { cfg } });
            let current = null, ready = false;
            live++;
            const feed = () => {
                if (next < cases.length) { current = cases[next++]; w.postMessage(current); }
                else { current = null; w.terminate(); }
            };
            w.on('message', r => { ready = true; if (r) record(r); feed(); });
            w.on('error', e => {
                if (!current) return;
                const { i, spec } = current;
                current = null;
                const who = `fuzz worker died: ${e && e.message || e}`;
                record({ i, tl: spec.tl, desc: '', lines: [`[${i}] ✗ ONE-SIDED FAILURE — ${who}`],
                    status: 'crash', ok: false, badMesh: false, discrepancy: null, loss: null, rlgc: null,
                    floating: null, crash: who });
            });
            w.on('exit', code => {
                live--;
                // A worker that never came up would fail again: stop instead of respawning.
                if (!ready) {
                    console.error(`fuzz worker failed to start (exit ${code})`);
                    process.exit(2);
                }
                if (next < cases.length) spawn();
                else if (live === 0) resolve();
            });
        };
        for (let k = 0; k < jobs; k++) spawn();
    });
    if (tty) process.stderr.write('\r\x1b[K');
    console.error(`finished in ${((Date.now() - t0) / 1000).toFixed(0)} s`);
    printSummary(results.filter(Boolean), cfg);
}
