# 2D Transmission Line Field Solver

A browser-based quasi-static 2D field solver for transmission line analysis. Computes characteristic impedance, effective permittivity, RLGC parameters, losses, and S-parameters.

## Features

- Transmission Line Types: Microstrip, stripline, GCPW (single-ended and differential)
- Electric Field Solving: Calculates characteristic impedance
- Full RLGC Extraction: Resistance, inductance, capacitance, conductance per unit length
- Loss Modeling: Conductor losses (skin effect, surface roughness) and dielectric losses
- S-Parameter Export: Touchstone .s2p and .s4p file generation
- Visualization: 2D potential plots, E-field streamlines, frequency-dependent plots
- Adaptive Meshing: Automatic mesh refinement for accurate field solutions

## Quick Start

Build WASM sparse matrix solver with "make" in "src/wasm_solver".
Open `src/field_solver.html` in a browser.

## Solution Flow

1. Geometry Setup: Define conductors and dielectric regions. Only rectangles are supported currently.
2. Mesh Generation: Create non-uniform coarse grid
3. Laplace Solve: Solve ∇²V = 0. Refine mesh until solution converges
4. Parameter Extraction:
   - Capacitance from field energy
   - Losses from perturbation methods
   - RLGC extraction
5. Frequency Sweep: Repeat for multiple frequencies
6. S-Parameters: Convert RLGC to S-parameters via ABCD matrix

### Validity

- Optimized for PCB Transmission Lines: Designed for microstrip, stripline, and coplanar waveguide structures commonly used in PCB RF and high-speed digital designs.
- TEM/Quasi-TEM Regime Accuracy: Provides reliable results when the transmission line supports TEM or quasi-TEM propagation, which covers most practical PCB geometries below the onset of higher-order modes.
- Validated Against Full-Wave Solvers: Results have been checked against full-wave EM solver and actual measurement data with different geometries and transmission line types, showing close agreement with small error in typical use cases.
- Wide Frequency Applicability: Accurate from RF through microwave and high-speed digital frequencies where return currents are confined and skin effect is significant.
- Sufficient for Most Practical Designs: Suitable for impedance control, loss estimation, and S-parameter generation in the vast majority of PCB transmission line applications.

### Limitations

- Mode Limitations (TEM/Quasi-TEM Only): Higher-order modes, dispersion, and cutoff behavior are not modeled. Structures that support non-TEM modes (e.g., waveguides) are outside the solver’s validity range.
- Conductor Current Distribution: Current is modeled at the surface for AC and DC resistance is blended smoothly at low frequencies. Full 2D/3D current density inside conductors is not solved, which can reduce accuracy at frequencies where skin depth is comparable to conductor thickness.
- Low-Frequency Inductance Accuracy: At DC and low frequencies (~<1 MHz), return current spreads over the ground plane. Since the solver infers inductance from capacitance, partial inductance and finite ground width effects may be inaccurate.
- 2D Cross-Section Only: Longitudinal variations, bends, vias, tapers, connectors, and transitions are not supported. Results apply to uniform, infinitely long transmission lines.
- Radiation is not modeled

## Common Tasks

### Run the Web App

```bash
# Open in browser
open src/field_solver.html
```

### Run Tests

```bash
node tests/test_vs_hfss.js
```

### Build WASM Solver

```bash
cd src/wasm_solver
make
```
