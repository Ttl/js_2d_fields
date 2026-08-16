// Validates openBoundaryWarnings(): each 'open' wall must be at least 3 substrate
// heights from the nearest conductor for the open-boundary approximation to hold.
//  - auto-sized geometries (large margins) produce no warnings,
//  - a cramped enclosure warns for the walls that are actually close,
//  - grounded walls never warn (the check only applies to 'open'),
//  - conductors touching a wall (coplanar ground pours, full-width ground planes)
//    don't trip the check for that wall.
import { MicrostripSolver } from '../src/microstrip.js';
import { BroadsideStriplineSolver } from '../src/broadside_stripline.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
    console.log(`  ${cond ? 'PASS' : 'FAIL'}: ${name}${detail && !cond ? ` — ${detail}` : ''}`);
    cond ? pass++ : fail++;
}

const h = 0.5e-3;
function ms(opts = {}) {
    return new MicrostripSolver({
        substrate_height: h, trace_width: 0.5e-3, trace_thickness: 35e-6,
        gnd_thickness: 35e-6, epsilon_r: 4.3, freq: 1e9,
        boundaries: ["open", "open", "open", "gnd"], ...opts,
    });
}

// Auto-sized microstrip: ample clearance, no warnings.
check('auto-sized microstrip: no warnings', ms().openBoundaryWarnings().length === 0);

// Cramped enclosure: trace fills most of the width and the top sits within 3·h.
// Trace is 0.5 mm wide; a 2 mm enclosure leaves 0.75 mm ≈ 1.5·h to each side wall;
// enclosure_height (air above the substrate) 1.2 mm leaves ~1.17 mm ≈ 2.3·h above
// the trace top.
const cramped = ms({ enclosure_width: 2e-3, enclosure_height: 1.2e-3 });
const wc = cramped.openBoundaryWarnings();
check('cramped microstrip: ONE warning naming "sides and top"',
    wc.length === 1 && /on the sides and top/.test(wc[0]), `got ${wc.length}: ${wc.join(' | ')}`);
check('warning states the rule', /3×/.test(wc[0]) && /open/i.test(wc[0]), wc[0]);

// Only ONE side open and cramped: keep naming that side, no merge.
const rightOnly = ms({ enclosure_width: 2e-3, enclosure_height: 1.2e-3,
    boundaries: ["gnd", "open", "gnd", "gnd"] });
const wr = rightOnly.openBoundaryWarnings();
check('single open side keeps its name', wr.length === 1 && /on the right/.test(wr[0]),
    wr.join(' | '));

// Asymmetric mix: left grounded, right + top open and cramped → "right and top".
const asymMix = ms({ enclosure_width: 2e-3, enclosure_height: 1.2e-3,
    boundaries: ["gnd", "open", "open", "gnd"] });
const wa = asymMix.openBoundaryWarnings();
check('asymmetric mix reads "right and top"', wa.length === 1 && /on the right and top/.test(wa[0]),
    wa.join(' | '));

// Same cramped enclosure with grounded walls: nothing to warn about.
const grounded = ms({ enclosure_width: 2e-3, enclosure_height: 1.2e-3,
    boundaries: ["gnd", "gnd", "gnd", "gnd"] });
check('grounded walls: no warnings', grounded.openBoundaryWarnings().length === 0);

// Only the top left open: exactly one warning.
const topOpen = ms({ enclosure_width: 2e-3, enclosure_height: 1.2e-3,
    boundaries: ["gnd", "gnd", "open", "gnd"] });
const wt = topOpen.openBoundaryWarnings();
check('only open walls are checked', wt.length === 1 && /on the top/.test(wt[0]),
    wt.join(' | '));

// Wide-but-low enclosure: sides are fine, only the top warns.
const lowTop = ms({ enclosure_width: 8e-3, enclosure_height: 1.2e-3 });
const wl = lowTop.openBoundaryWarnings();
check('sides clear, low top: warns only for top', wl.length === 1 && /on the top/.test(wl[0]) && !/side/.test(wl[0]),
    wl.join(' | '));

// The full-width ground plane touches the side walls; it must not trip the side
// check by itself (distance 0). With a huge enclosure nothing else is close.
const wide = ms({ enclosure_width: 30e-3, enclosure_height: 20e-3 });
check('wall-touching ground plane ignored', wide.openBoundaryWarnings().length === 0);

// GCPW shielding: the coplanar pour touches the side walls AND lies between the
// trace and the wall at trace level, so the open side boundary sees no field —
// no side warning even when the raw trace-to-wall clearance is under 3·h. (The
// real-world trigger: solder mask inflates the dielectric-stack h to the trace
// top + mask, and the UI-default GCPW's clearance — measured THROUGH the pour —
// fell under the threshold.) The same clearance without the pour must still warn.
const smOpts = { use_sm: true, sm_t_sub: 20e-6, sm_t_trace: 20e-6, sm_t_side: 20e-6, sm_er: 3.5, sm_tand: 0.02 };
const gcpwGeom = {
    substrate_height: 0.21e-3, trace_width: 0.35e-3, trace_thickness: 35e-6,
    gnd_thickness: 35e-6, epsilon_r: 4.4, freq: 1e9,
    boundaries: ["open", "open", "open", "gnd"], ...smOpts,
    enclosure_width: 1.75e-3,   // trace 0.70 mm from each side wall < 3·h = 0.795 mm
};
const gcpwShielded = new MicrostripSolver({ ...gcpwGeom,
    use_coplanar_gnd: true, gap: 0.1e-3, via_gap: 0.1e-3, use_vias: true });
const wg = gcpwShielded.openBoundaryWarnings();
check('gcpw: pour shields the open sides (no warning)', wg.length === 0, wg.join(' | '));
const msSameClearance = new MicrostripSolver(gcpwGeom);
const wm = msSameClearance.openBoundaryWarnings();
check('same clearance without the pour still warns for the sides',
    wm.length === 1 && /on the sides/.test(wm[0]), wm.join(' | '));

// Broadside pair (different line type, same rule): auto-sized is quiet, a narrow
// enclosure warns for the open sides. Substrate stack = h_bottom+h_middle+h_top = 0.9 mm.
function bs(opts = {}) {
    return new BroadsideStriplineSolver({
        trace_width: 0.2e-3, trace_thickness: 35e-6, sigma_cond: 5.8e7,
        h_bottom: 0.3e-3, er_bottom: 4.4, tand_bottom: 0.02,
        h_middle: 0.3e-3, er_middle: 3.0, tand_middle: 0.01,
        h_top: 0.3e-3, er_top: 4.4, tand_top: 0.02,
        x_offset: 0, freq: 1e9, boundaries: ["open", "open", "gnd", "gnd"], ...opts,
    });
}
check('auto-sized broadside: no warnings', bs().openBoundaryWarnings().length === 0);
const bsNarrow = bs({ enclosure_width: 1.5e-3 });   // ~0.65 mm ≈ 0.7 stack heights to sides
const wb = bsNarrow.openBoundaryWarnings();
check('narrow broadside: one merged "sides" warning',
    wb.length === 1 && /on the sides/.test(wb[0]), wb.join(' | '));

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
