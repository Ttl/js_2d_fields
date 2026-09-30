// Accumulator for surface current plot segments { x0, y0, x1, y1, K }.
export function segmentBuffer() {
    const seg = { x0: [], y0: [], x1: [], y1: [], K: [] };
    return {
        push(xa, ya, xb, yb, K) { seg.x0.push(xa); seg.y0.push(ya); seg.x1.push(xb); seg.y1.push(yb); seg.K.push(K); },
        // Typed arrays, K multiplied by `scale`.
        out(scale = 1) {
            const o = {};
            for (const k of ['x0', 'y0', 'x1', 'y1']) o[k] = Float64Array.from(seg[k]);
            o.K = scale === 1 ? Float64Array.from(seg.K) : Float64Array.from(seg.K, v => v * scale);
            return o;
        },
    };
}
