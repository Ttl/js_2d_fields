# RLGC Matrix for Differential Lines

For differential transmission lines, RLGC parameters can be represented in three ways:

1. Modal (odd/even modes)
2. Physical (per-trace coupling)
3. Differential (diff/common impedances)

## Modal Representation

Solver computes odd and even mode parameters separately.

### Odd Mode (Differential)
- Traces driven with opposite polarity: V1 = +V, V2 = -V
- Return current flows between traces
- Higher capacitance (opposite charges attract)
- Lower inductance (opposite currents cancel flux)

### Even Mode (Common)
- Traces driven with same polarity: V1 = V2 = V
- Return current flows through ground plane
- Lower capacitance (same charges repel)
- Higher inductance (same-direction currents add flux)

### Output in solve_adaptive()
```javascript
{
  modes: [
    { mode: 'odd', Z0, eps_eff, RLGC, ... },
    { mode: 'even', Z0, eps_eff, RLGC, ... }
  ],
  Z_diff: 2 * Z_odd,
  Z_common: Z_even / 2
}
```

Results plot shows modal parameters. Modal parameters diagonalize the RLGC
system, so coupling appears only when transforming back to the physical
(per-trace) representation.

## Physical Representation

2x2 matrices describing per-trace behavior with mutual coupling.

### Transformation from Modal to Physical

```
Self terms (diagonal):     X11 = X22 = (X_odd + X_even) / 2
Mutual terms (off-diag):   X12 = X21 = (X_even - X_odd) / 2
```

where X represents R, L, G, or C and the transmission line is symmetrical.

Function `_modal_to_physical_rlgc()` in field_solver.js performs this conversion.

### Telegrapher's Equations

```
V1' = -(R11*I1 + R12*I2 + jω(L11*I1 + L12*I2))
V2' = -(R21*I1 + R22*I2 + jω(L21*I1 + L22*I2))

I1' = -(G11*V1 + G12*V2 + jω(C11*V1 + C12*V2))
I2' = -(G21*V1 + G22*V2 + jω(C21*V1 + C22*V2))
```

where ' denotes derivate in length direction. 

### Matrix Elements

Diagonal (11, 22) - self parameters:
- R11, R22: series resistance per unit length (Ohm/m)
- L11, L22: self-inductance per unit length (H/m)
- G11, G22: shunt conductance per unit length (S/m)
- C11, C22: self-capacitance per unit length (F/m)

Off-diagonal (12, 21) - mutual coupling:
- R12, R21: mutual resistance (Ohm/m)
- L12, L21: mutual inductance (H/m)
- G12, G21: mutual conductance (S/m)
- C12, C21: mutual capacitance (F/m)

### Sign Convention

Mutual capacitance C12 and C21 are negative:
- C_odd > C_even (differential excitation increases capacitance)
- C12 = (C_even - C_odd)/2 < 0

Mutual inductance L12 and L21 are positive (current defined in same direction):
- L_even > L_odd (common-mode currents increase inductance)
- L12 = (L_even - L_odd)/2 > 0

### Output in solve_adaptive()
```javascript
{
  RLGC_matrix: {
    R: [[R11, R12], [R21, R22]],
    L: [[L11, L12], [L21, L22]],
    G: [[G11, G12], [G21, G22]],
    C: [[C11, C12], [C21, C22]]
  }
}
```

### Output in solve_sweep()
```javascript
{
  frequencies: [...],
  RLGC_matrix: {
    R: { R11: [...], R12: [...], R21: [...], R22: [...] },
    L: { L11: [...], L12: [...], L21: [...], L22: [...] },
    G: { G11: [...], G12: [...], G21: [...], G22: [...] },
    C: { C11: [...], C12: [...], C21: [...], C22: [...] }
  }
}
```

## Differential RLGC parameters

```
X_dd = 2*X_odd
X_cc = 0.5*X_even
X_dc = X_dc = 0
```

Assumes symmetrical transmission line.
