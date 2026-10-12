// Time-domain link simulation (src/eye) against ngspice.
//
// The driver EMF goes to ngspice as a finely sampled PWL of the same edges, the
// receiver and TX pad voltages are compared over the whole run.
//
//   1. Lossy line (LTRA, constant R, L, C), mismatched source and load, pulse.
//   2. AC coupling (1 nF) at both ends, pulse with its droop over 40 ns (the slow-mode
//      split).
//   3. Terminations: V_TT, R + C to ground, open; step.
//   4. Clock burst on a mismatched lossy line.
//   5. Packages: R + L form at TX, line form (DC resistance only, an LTRA) at RX, with pin
//      and pad capacitances.
//
// These cover the time-domain machinery (transform, DC point, AC coupling split,
// terminations, settling window). The line transfer itself is checked exactly in
// test_eye_transient.mjs: pairs against the MTL S-parameters, frequency-dependent loss
// against the analytic skin-effect step. ngspice has no exact reference for those
// (its CPL misses DC levels by ~2.6 %, ladders need thousands of sections).
//
// Needs ngspice 47 or newer: in 45 the LTRA drifts off after a fixed number of internal
// time steps (about 8 ns at a 1 ps step: a discharged line jumps to -14 mV and stays
// near -6 mV) and long LTRA runs are very slow. reltol 1e-6 stalls ngspice on a clock
// burst, so the runs use 1e-5. Skipped when ngspice is missing or older.
import { execFileSync } from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';
import { check, done } from './helpers.js';
import { loadFFT } from '../src/eye/fft.js';
import { simulateLink } from '../src/eye/transient.js';
import { buildStimulus, legValue } from '../src/eye/stimulus.js';

let ngVersion = 0;
try {
    const m = execFileSync('ngspice', ['--version'], { encoding: 'utf8' }).match(/ngspice-(\d+)/);
    ngVersion = m ? parseInt(m[1], 10) : 0;
} catch {
    console.log('ngspice not found: skipped');
    process.exit(0);
}
if (ngVersion < 47) {
    console.log(`ngspice ${ngVersion || '(unknown version)'} is older than 47 (LTRA drifts on long runs): skipped`);
    process.exit(0);
}

const fft = await loadFFT();
const t0 = Date.now();
const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'eye-ngspice-'));

// Run a netlist body; returns { t, v: { name: Float64Array } } for the node names.
function spice(body, nodes, tStop, tMax) {
    const file = path.join(tmp, 'c.cir'), out = path.join(tmp, 'out.txt');
    const net = `eye comparison\n${body}\n.options reltol=1e-5 abstol=1e-12 vntol=1e-9\n` +
        `.tran ${tMax / 4} ${tStop} 0 ${tMax}\n.control\nset wr_singlescale\nrun\n` +
        `wrdata ${out} ${nodes.map(n => `v(${n})`).join(' ')}\n.endc\n.end\n`;
    fs.writeFileSync(file, net);
    execFileSync('ngspice', ['-b', file], { stdio: 'ignore' });
    const rows = fs.readFileSync(out, 'utf8').trim().split('\n').map(l => l.trim().split(/\s+/).map(Number));
    const t = Float64Array.from(rows, r => r[0]);
    const v = Object.fromEntries(nodes.map((n, k) => [n, Float64Array.from(rows, r => r[k + 1])]));
    return { t, v };
}

// PWL source of a driver leg: dense samples across every edge, straight between.
function pwl(stim, leg, tStop) {
    const pts = [[0, legValue(stim, leg, 0)]];
    const reps = isFinite(stim.count) ? stim.count : 1;
    const edges = [];
    for (let k = 0; k < reps; k++) for (const e of leg.edges) edges.push({ t: e.t + k * (stim.period || 0), h: e.edge.half });
    edges.sort((a, b) => a.t - b.t);
    for (const e of edges) {
        for (let i = 0; i <= 80; i++) {
            const t = e.t - e.h + 2 * e.h * i / 80;
            if (t > pts[pts.length - 1][0]) pts.push([t, legValue(stim, leg, t)]);
        }
    }
    pts.push([tStop, legValue(stim, leg, tStop)]);
    return `PWL(${pts.map(([t, v]) => `${t.toExponential(9)} ${v.toPrecision(12)}`).join(' ')})`;
}

// Max difference of an Eye trace and an ngspice node over the Eye time grid up to the
// end of the ngspice run, as a fraction of the trace's swing.
function compare(res, name, sp, node) {
    const v = res.traces[name], t = res.t, st = sp.t, sv = sp.v[node];
    let lo = Infinity, hi = -Infinity, e = 0, j = 1;
    for (let i = 0; i < t.length; i++) {
        if (t[i] < 0 || t[i] > st[st.length - 1]) continue;
        while (j < st.length - 1 && st[j] < t[i]) j++;
        const a = (t[i] - st[j - 1]) / (st[j] - st[j - 1]);
        const s = sv[j - 1] + a * (sv[j] - sv[j - 1]);
        e = Math.max(e, Math.abs(v[i] - s));
        lo = Math.min(lo, v[i]); hi = Math.max(hi, v[i]);
    }
    return e / Math.max(hi - lo, 1e-12);
}

const FREQS = [0, 1e6, 1e8, 1e10, 1e12];
const seTable = (fn, freqs = FREQS) => {
    const data = { R: [], L: [], G: [], C: [] };
    for (const f of freqs) { const v = fn(f); for (const k in data) data[k].push(v[k]); }
    return { differential: false, freqs: Float64Array.from(freqs), keys: ['R', 'L', 'G', 'C'],
        data: Object.fromEntries(Object.entries(data).map(([k, v]) => [k, Float64Array.from(v)])),
        lext: { L: fn(freqs[freqs.length - 1]).L } };
};
const STIM = { type: 'pulse', width: 1e-9, vLow: 0, vHigh: 1, tr: 30e-12, tf: 30e-12, edgeDef: '10-90', shape: 'rcos', delay: 0 };
const run = (table, over) => {
    const cfg = { length: 0.1, samplesPerEdge: 16, ...over,
        network: { rs: 50, tx: { ac: false }, rx: { ac: false }, term: { type: 'r', rt: 50 }, ...(over.network || {}) },
        stimulus: { ...STIM, ...(over.stimulus || {}) } };
    const res = simulateLink(table, cfg, fft);
    const stim = buildStimulus({ ...cfg.stimulus, differential: table.differential });
    return { res, stim, cfg };
};
const tEnd = (res) => res.t[res.t.length - 1];
const TOL = 2e-3;   // of the swing

// Single-ended netlist around an LTRA line from the Eye config.
function seNetlist(stim, cfg, line, tStop) {
    const n = cfg.network, l = [];
    l.push(`V1 src 0 ${pwl(stim, stim.legs[0], tStop)}`);
    l.push(`Rs src a ${n.rs}`);
    if (n.tx.ac) l.push(`Ctx a b ${n.tx.c}`); else l.push('Vtx a b 0');
    l.push(line);
    if (n.rx.ac) l.push(`Crx c d ${n.rx.c}`); else l.push('Vrx c d 0');
    const t = n.term;
    if (t.type === 'r') l.push(`Rt d 0 ${t.rt}`);
    if (t.type === 'vtt') l.push(`Rt d tt ${t.rt}`, `Vtt tt 0 ${t.vtt}`);
    if (t.type === 'rc') l.push(`Rt d m ${t.rt}`, `Ct m 0 ${t.ct}`);
    // The Eye model leaks floating nodes to ground with 1 GOhm.
    if (n.tx.ac && n.rx.ac) l.push('Rleak b 0 1e9');
    if (t.type === 'open' || (n.rx.ac && t.type === 'diff')) l.push('Rleakd d 0 1e9');
    return l.join('\n');
}
const ltra = (R, L, C, len) => `O1 b 0 c 0 lmod\n.model lmod ltra R=${R} L=${L} G=0 C=${C} LEN=${len}`;

const lossy = seTable(() => ({ R: 20, L: 250e-9, G: 0, C: 100e-12 }));

// 1. Lossy, mismatched.
{
    const { res, stim, cfg } = run(lossy, { network: { rs: 25, term: { type: 'r', rt: 100 } } });
    const sp = spice(seNetlist(stim, cfg, ltra(20, 250e-9, 100e-12, 0.1), tEnd(res)), ['d', 'a'], tEnd(res), 1e-12);
    const eR = compare(res, 'rx', sp, 'd'), eT = compare(res, 'tx', sp, 'a');
    check('LTRA lossy mismatched line: RX and TX pad', eR < TOL && eT < TOL, `RX ${(100 * eR).toFixed(3)} %, TX ${(100 * eT).toFixed(3)} % of swing`);
}

// 2. AC coupling at both ends.
{
    const { res, stim, cfg } = run(lossy, { span: 40e-9, network: { tx: { ac: true, c: 1e-9 }, rx: { ac: true, c: 1e-9 } } });
    const sp = spice(seNetlist(stim, cfg, ltra(20, 250e-9, 100e-12, 0.1), tEnd(res)), ['d', 'a'], tEnd(res), 1e-12);
    const eR = compare(res, 'rx', sp, 'd'), eT = compare(res, 'tx', sp, 'a');
    check('AC coupling, 1 nF, 40 ns: RX and TX pad', eR < TOL && eT < TOL, `RX ${(100 * eR).toFixed(3)} %, TX ${(100 * eT).toFixed(3)} %, ${res.info.slowModes} slow modes`);
}

// 3. Terminations, step.
for (const term of [{ type: 'vtt', rt: 50, vtt: 1.2 }, { type: 'rc', rt: 50, ct: 20e-12 }, { type: 'open' }]) {
    const { res, stim, cfg } = run(lossy, { stimulus: { type: 'step' }, span: 10e-9, network: { term } });
    const sp = spice(seNetlist(stim, cfg, ltra(20, 250e-9, 100e-12, 0.1), tEnd(res)), ['d', 'a'], tEnd(res), 1e-12);
    const eR = compare(res, 'rx', sp, 'd'), eT = compare(res, 'tx', sp, 'a');
    check(`termination ${term.type}, step: RX and TX pad`, eR < TOL && eT < TOL, `RX ${(100 * eR).toFixed(3)} %, TX ${(100 * eT).toFixed(3)} %`);
}

// 4. Clock burst.
{
    const { res, stim, cfg } = run(lossy, { stimulus: { type: 'clock', clockFreq: 2.5e9, duty: 0.5, cycles: 8 }, network: { rs: 30, term: { type: 'r', rt: 80 } } });
    const sp = spice(seNetlist(stim, cfg, ltra(20, 250e-9, 100e-12, 0.1), tEnd(res)), ['d', 'a'], tEnd(res), 1e-12);
    const eR = compare(res, 'rx', sp, 'd');
    check('clock burst, 8 cycles at 2.5 GHz: RX', eR < TOL, `${(100 * eR).toFixed(3)} %`);
}

// 5. Packages at both ends.
{
    const pkgTx = { form: 'rl', cpin: 0.25e-12, r: 0.5, l: 1e-9 };
    const pkgRx = { form: 'line', cpin: 0.25e-12, z0: 40, delay: 140e-12, r: 0.3 };
    const { res, stim, cfg } = run(lossy, { network: { rs: 50, padTx: 1e-12, padRx: 0.5e-12, pkgTx, pkgRx, term: { type: 'r', rt: 50 } } });
    const body = [
        `V1 src 0 ${pwl(stim, stim.legs[0], tEnd(res))}`, 'Rs src a 50', 'Cpadt a 0 1e-12',
        'Rpt a x 0.5', 'Lpt x b 1e-9', 'Cpint b 0 0.25e-12',
        ltra(20, 250e-9, 100e-12, 0.1),
        'Cpinr c 0 0.25e-12', `O2 c 0 d 0 pmod\n.model pmod ltra R=0.3 L=${40 * 140e-12} G=0 C=${140e-12 / 40} LEN=1`,
        'Cpadr d 0 0.5e-12', 'Rt d 0 50',
    ].join('\n');
    const sp = spice(body, ['d', 'a'], tEnd(res), 1e-12);
    const eR = compare(res, 'rx', sp, 'd'), eT = compare(res, 'tx', sp, 'a');
    check('packages (R + L at TX, line at RX): RX and TX die pads', eR < TOL && eT < TOL, `RX ${(100 * eR).toFixed(3)} %, TX ${(100 * eT).toFixed(3)} %`);
}

fs.rmSync(tmp, { recursive: true, force: true });
console.log(`(${((Date.now() - t0) / 1000).toFixed(0)} s)`);
done();
