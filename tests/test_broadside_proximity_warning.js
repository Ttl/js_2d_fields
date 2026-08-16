// Broadside strong-coupling accuracy warning (QS backend only).
//
// The rectilinear vacuum-field loss integrand over-weights the facing surfaces
// of a strongly coupled broadside pair (loss ∝ charge² vs the true current
// distribution): measured vs the tri multi-drive MQS reference, w/facing-gap
// ≥ 1.75 reads R +10…+45% high. The QS surfaces a type:'accuracy' warning there
// (recommending the full-wave backend); weakly coupled pairs and the triangular
// backend must stay clean. See docs/backend_agreement_handover.md issue 1.
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';

let failures = 0;
function check(name, ok, detail = '') {
    console.log(`${ok ? '✓ PASS' : '✗ FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
    if (!ok) failures++;
}
const quiet = async (fn) => {
    const log = console.log, warn = console.warn;
    console.log = () => {}; console.warn = () => {};
    try { return await fn(); } finally { console.log = log; console.warn = warn; }
};

const base = {
    trace_thickness: 17e-6, gnd_thickness: 35e-6,
    h_bottom: 0.2e-3, h_top: 0.2e-3,
    er_bottom: 4.1, er_middle: 4.1, er_top: 4.1,
    tand_bottom: 0.02, tand_middle: 0.02, tand_top: 0.02,
    sigma_cond: 5.8e7, freq: 2e9, rq: 0, nx: 30, ny: 30,
};
// Coupled: w/(h_middle − 2t) = 0.3/0.066 ≈ 4.5 (well past the 1.75 threshold).
// Weak: 0.1/0.366 ≈ 0.27.
const coupled = { ...base, trace_width: 0.3e-3, h_middle: 0.1e-3 };
const weak = { ...base, trace_width: 0.1e-3, h_middle: 0.4e-3 };

const proximityWarns = r => (r.warnings || []).filter(w =>
    w.type === 'accuracy' && /broadside/i.test(w.message));

async function solveQS(o) {
    const s = new BroadsideStriplineSolver({ ...o, mesh_backend: 'rectilinear' });
    return quiet(() => s.solve_adaptive({
        max_iters: 6, energy_tol: 0.01, param_tol: 0.05, max_nodes: 8000, min_converged_passes: 1,
    }));
}

{
    const r = await solveQS(coupled);
    check('QS coupled broadside (w/gap ≈ 4.5) surfaces the proximity accuracy warning',
        proximityWarns(r).length === 1,
        (r.warnings || []).map(w => w.type).join(',') || 'no warnings');
}
{
    const r = await solveQS(weak);
    check('QS weakly coupled broadside (w/gap ≈ 0.27) stays clean',
        proximityWarns(r).length === 0,
        (r.warnings || []).map(w => w.type).join(',') || 'no warnings');
}
{
    // Offset past the trace width (zero vertical overlap) must still warn — the
    // worst measured case (s29#0, +45%) had no overlap; proximity acts diagonally.
    const r = await solveQS({ ...coupled, x_offset: 0.35e-3 });
    check('QS coupled broadside with zero-overlap offset still warns',
        proximityWarns(r).length === 1,
        (r.warnings || []).map(w => w.type).join(',') || 'no warnings');
}
{
    // Triangular backend on the same coupled geometry: multi-drive MQS models the
    // proximity accurately — no proximity warning may leak into its results.
    const s = new BroadsideStriplineSolver({ ...coupled, mesh_backend: 'triangular' });
    const r = await quiet(() => s.solve_adaptive({
        max_iters: 6, energy_tol: 0.01, param_tol: 0.05, max_nodes: 8000, min_converged_passes: 1,
    }));
    check('triangular backend result carries no proximity warning',
        proximityWarns(r).length === 0,
        (r.warnings || []).map(w => w.type).join(',') || 'no warnings');
}

console.log(failures === 0 ? '\nALL BROADSIDE PROXIMITY WARNING TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
