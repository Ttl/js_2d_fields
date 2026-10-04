// Fill colors of conductors and dielectrics in the plots, shared with the custom geometry
// editor so its color swatches show what the plot draws.

export const CONDUCTOR_COLOR = '#d97706';

// Plot background and air as the geometry view shows it: white at 0.8 over #1a1a1a.
const PLOT_BG = 0x1a;
export const AIR_GREY = 0.8 * 255 + 0.2 * PLOT_BG;

// Default color of a dielectric as [r, g, b]: white for air, green shades getting lighter
// with er. Drawn with transparency.
// A conducting dielectric (sigma in S/m) is tinted towards brown, more the more it
// conducts: 1 S/m a little, 1e4 S/m fully.
export function dielectricRGB(er, sigma = 0) {
    const s = sigma > 0 ? Math.min(1, 0.25 + Math.log10(1 + sigma) / 5.3) : 0;
    if (er <= 1.01 && !s) return [255, 255, 255];
    const g = er <= 1.01 ? 160 : Math.min(255, 100 + (er - 1) * 30);
    return [100 + 60 * s, g - (g - 90) * 0.5 * s, 100 - 40 * s];
}

// [r, g, b] drawn at `alpha` over air, as an opaque color.
export const overAir = (rgb, alpha) => rgb.map(c => Math.round(alpha * c + (1 - alpha) * AIR_GREY));

const hex2 = v => Math.round(v).toString(16).padStart(2, '0');
export const rgbToHex = rgb => '#' + rgb.map(hex2).join('');

// [r, g, b] of a '#rgb' or '#rrggbb' color, null for anything else.
export function hexToRGB(s) {
    const m = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(String(s ?? '').trim());
    if (!m) return null;
    const h = m[1].length === 3 ? m[1].replace(/./g, c => c + c) : m[1];
    return [0, 2, 4].map(i => parseInt(h.slice(i, i + 2), 16));
}
