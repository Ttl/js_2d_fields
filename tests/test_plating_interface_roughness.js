// Roughness of the plating/bulk interface (plating.rq_interface, sidebar Interface
// Roughness RMS, custom text rq_iface= / plating_rq_iface=).
//
//   1 the layered model: left out it is the plating roughness, bit-identical to passing
//     it; a smooth interface changes Zs where the current reaches it (t ~ delta) only.
//   2 the three plating models see the same interface: on a wide microstrip the
//     quasi-static and full-wave layered surface impedance and the full-wave meshed
//     plating (whose mesh holds the interface smooth, with the 1D increment added on
//     the plated faces) move R by the same amount between a rough and a smooth
//     interface.
//   3 the custom geometry text carries the key both ways.
//
// Run: node tests/test_plating_interface_roughness.js
import { calculate_Zrough_layered } from '../src/surface_roughness.js';
import { MicrostripSolver } from '../src/microstrip.js';
import { CustomGeometrySolver } from '../src/custom_geometry.js';
import { parseAndEvaluate, solverToGeometryText } from '../src/custom_geometry_text.js';
import { check, quiet, APP, done } from './helpers.js';

// --- 1 ---
{
    console.log('1 layered model');
    let same = true;
    for (const f of [1e8, 1e9, 1e10, 1e11]) {
        const a = calculate_Zrough_layered(f, 5.8e7, 1e-6, 8.7e6, 3e-6);
        const b = calculate_Zrough_layered(f, 5.8e7, 1e-6, 8.7e6, 3e-6, 1e-6);
        if (a.re !== b.re || a.im !== b.im) same = false;
    }
    check('interface roughness left out equals the plating roughness, bit-identical', same);
    const ratio = f => calculate_Zrough_layered(f, 5.8e7, 1e-6, 8.7e6, 3e-6).re
        / calculate_Zrough_layered(f, 5.8e7, 1e-6, 8.7e6, 3e-6, 0).re;
    check('3 um Sn on 1 um roughness: a rough interface changes Re(Zs) by > 3% at 3 GHz',
        Math.abs(ratio(3e9) - 1) > 0.03, ratio(3e9).toFixed(4));
    check('and by < 0.2% at 30 GHz (the current no longer reaches it)',
        Math.abs(ratio(30e9) - 1) < 0.002, ratio(30e9).toFixed(4));
    const thick = calculate_Zrough_layered(1e9, 5.8e7, 1e-6, 4e7, 30e-6, 0).re
        / calculate_Zrough_layered(1e9, 5.8e7, 1e-6, 4e7, 30e-6).re;
    check('30 um plating at 1 GHz: interface roughness has no effect', Math.abs(thick - 1) < 1e-6, thick.toFixed(8));
    const z = calculate_Zrough_layered(3e9, 5.8e7, 0, 8.7e6, 1e-6, 2e-6);
    check('interface rougher than the plating is thick stays finite', Number.isFinite(z.re) && z.re > 0, z.re);
}

// --- 2 ---
{
    console.log('2 quasi-static layered, full-wave layered and full-wave meshed plating');
    const MS = { trace_width: 1.5e-3, substrate_height: 0.254e-3, trace_thickness: 35e-6, gnd_thickness: 35e-6,
        epsilon_r: 4, tan_delta: 0.01, sigma_cond: 5.8e7, freq: 3e9, rq: 0, nx: 30, ny: 30,
        boundaries: ['open', 'open', 'open', 'gnd'] };
    const pl = (thick, rqi) => ({ sigma: 8.7e6, thickness: 3e-6, rq: 1e-6, rq_interface: rqi,
        top: true, sides: true, bottom: true, thick_corners: thick });
    const R = async (backend, thick, rqi) => {
        const s = new MicrostripSolver({ ...MS, plating: pl(thick, rqi), mesh_backend: backend });
        s.use_causal_materials = false;
        const m = (await quiet(() => s.solve_adaptive({ ...APP }))).modes[0];
        return { R: m.RLGC.R, via: m.lossVia };
    };
    const shift = {};
    for (const [name, backend, thick] of [['QS layered', 'rectilinear', false],
        ['FW layered', 'triangular', false], ['FW meshed', 'triangular', true]]) {
        const rough = await R(backend, thick, undefined), smooth = await R(backend, thick, 0);
        shift[name] = rough.R / smooth.R - 1;
        if (backend === 'triangular') check(`${name}: MQS loss`, rough.via === 'mqs', rough.via);
    }
    const txt = Object.entries(shift).map(([n, v]) => `${n} ${(v * 100).toFixed(2)}%`).join(', ');
    check('rough interface lowers R at 3 GHz on every model', Object.values(shift).every(v => v < -0.05), txt);
    check('full-wave meshed matches full-wave layered within 1 pp', Math.abs(shift['FW meshed'] - shift['FW layered']) < 0.01, txt);
    check('full-wave layered matches quasi-static within 0.5 pp', Math.abs(shift['FW layered'] - shift['QS layered']) < 0.005, txt);
}

// --- 3 ---
{
    console.log('3 custom geometry text');
    const g = parseAndEvaluate('units um\nplating sigma=1e7 t=4 rq=1 rq_iface=0.5\n'
        + 'sig+ x=0 y=0 w=100 h=35 plating=top plating_rq_iface=0.2\nsig+ x=200 y=0 w=100 h=35 plating=top');
    check('rq_iface and plating_rq_iface evaluate to metres', g.errors.length === 0
        && Math.abs(g.plating.rq_interface - 0.5e-6) < 1e-15
        && Math.abs(g.rects[0].platingMaterial.rq_interface - 0.2e-6) < 1e-15, JSON.stringify(g.errors));
    check('negative rq_iface is an error',
        parseAndEvaluate('units um\nplating sigma=1e7 t=4 rq_iface=-1').errors.length === 1);
    const native = rqi => new MicrostripSolver({ trace_width: 0.3e-3, substrate_height: 0.2e-3, trace_thickness: 35e-6,
        gnd_thickness: 35e-6, epsilon_r: 4, tan_delta: 0.01, sigma_cond: 5.8e7, freq: 1e9,
        boundaries: ['open', 'open', 'open', 'gnd'],
        plating: { sigma: 1e7, thickness: 4e-6, rq: 1e-6, rq_interface: rqi, top: true, sides: true, bottom: false } });
    const text = solverToGeometryText(native(0.3e-6));
    check('native geometry text carries rq_iface', /rq_iface=/.test(text), text.split('\n').find(l => l.startsWith('plating')));
    const back = new CustomGeometrySolver({ text, nx: 30, ny: 30, freq: 1e9 });
    const c = back.conductors.find(o => o.is_signal);
    check('and it builds the same interface roughness', c && Math.abs(c.plating.rq_interface - 0.3e-6) < 1e-15,
        c && c.plating.rq_interface);
    const plain = solverToGeometryText(native(undefined));
    check('without an interface roughness the text has no rq_iface', !/rq_iface/.test(plain));
}

done();
