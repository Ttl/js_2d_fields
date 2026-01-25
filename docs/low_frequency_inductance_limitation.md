# Low-Frequency Inductance Limitation

Reference EM simulator results at 1 kHz with different ground plane widths:

| Ground Width | L @ 1 kHz | Change |
|--------------|----------------|--------|
| 30 mm | 447 nH/m | baseline |
| 60 mm | 553 nH/m | +24% |
| 100 mm | 641 nH/m | +43% |

Inductance increases with ground plane width at low frequency.

## Physical Explanation

### At Low Frequency (< 1 MHz)

At DC and low frequencies:

1. Return current spreads out across the ground plane
2. Wider ground -> current can spread further
3. Larger current loop -> more magnetic flux -> higher inductance

Why this happens:
- At DC, current minimizes DC resistance
- Spreading reduces resistance: R proportional to 1/width
- Return current uses all available ground width
- More spreading -> larger magnetic flux loop -> higher L

### At High Frequency (> 100 MHz)

At RF/microwave:

1. Skin effect confines return current to surface
2. Current flows in narrow path directly under trace
3. Ground width far away from the trace has minimal effect
4. Inductance is constant for wide-enough ground

## Static Solver Limitation

```javascript
L_external = Z0^2 * C
```

Where:
- `Z0 = sqrt(L/C), v = 1/sqrt(L*C)` For TEM: `v = c/sqrt(eps_eff) -> Z0 = 1/(c * sqrt(C * C0))`.
- `C` from electrostatic field (Laplace equation)
- `C0` vacuum capacitance

### Why It's Wrong at Low Frequency

Capacitance is an electrostatic quantity:
- Determined by Laplace equation
- Fields reach steady state when ground is very wide
- Doesn't change much with wider ground

TEM approximation assumes L and C are related by:
- Z0 = sqrt(L/C) -> L = Z0^2 * C
- Valid for propagating waves (high frequency)
- Solver captures the electrostatic limit, but low-frequency inductance is governed by the magnetostatic limit

### Test Results

This solver's L decreases with wider ground (opposite of physics) because:
- Wider ground -> slightly higher C
- Higher C -> lower `Z0 = 1/(c*sqrt(C*C0))`
- Lower Z0 -> lower `L = Z0^2 * C`

## Frequency Dependence

### High Frequency (1 GHz)

At 1 GHz, ground width has minimal effect (both solvers agree):

| Ground Width | L @ 1 GHz | Reference L @ 1 GHz |
|--------------|---------------|----------------|
| 30 mm | 276 nH/m | 280 nH/m |
| 100 mm | 280 nH/m | 282 nH/m |

### This Solver (Quasi-Static Electric)

```
∇²V = 0  ->  E-field  ->  C  ->  Z0 = 1/(c*sqrt(C*C0)) ->  L = Z0^2 * C
```
- Solves for electric potential only
- Assumes TEM relationship between E and H fields
- Valid for propagating waves (high frequency)

### What's Needed (Magnetostatic)

```
∇xH = J  ->  H-field  ->  Magnetic flux (Phi)  ->  L = Phi/I
```

- Need to solve for current distribution in conductors
- Calculate magnetic field from currents
- Integrate magnetic flux through loop
- Account for finite ground boundaries

Needs big changes in solver.

## Why We Can't Easily Fix This

To calculate low-frequency inductance correctly, we would need to:

1. Solve for current distribution in both signal and ground
   - Not just surface current (skin effect)
   - Full current density in conductors
   - Depends on conductor geometry and resistivity

2. Calculate magnetic field from current distribution
   - Integration over conductor volumes

3. Integrate magnetic flux through the loop
   - Depends on actual current paths

4. Get inductance: L = Phi/I

## Practical Implications

### When Solver is Accurate

- RF and Microwave (> 100 MHz)
- Skin effect confines return current
- TEM approximation valid
- Ground width effect minimal

### When Solver Has Limitations

- Low Frequency (< 1 MHz), but DC is special cased
- Return current spreads
- Ground width strongly affects L
- TEM approximation invalid
- This is only limitation when analyzing transmission line behaviour at these
  frequencies, loss is roughly correct so for typical PCB lines this should have
  minimal importance.
