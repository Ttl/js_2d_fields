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

1. Open `src/field_solver.html` in a browser

## Solution Flow

1. Geometry Setup: Define conductors and dielectric regions. Only rectangles are supported currently.
2. Mesh Generation: Create non-uniform coarse grid
3. Laplace Solve: Solve ∇²V = 0. Refine mesh until solution converges
4. Parameter Extraction:
   - Capacitance from field energy
   - Losses from perturbation method
   - RLGC extraction
5. Frequency Sweep: Repeat for multiple frequencies
6. S-Parameters: Convert RLGC to S-parameters via ABCD matrix

### Validity

- Designed for microstrip, stripline, and coplanar waveguide structures commonly used in PCB RF and high-speed digital designs.
- Provides good results when the transmission line supports TEM or quasi-TEM propagation, which covers most practical PCB geometries below the onset of higher-order modes.
- Results have been checked against EM solver and actual measurement data with different geometries and transmission line types, showing close agreement with small error in typical use cases (`tests` folder).
- Accurate from RF through microwave and high-speed digital frequencies where return currents are confined and skin effect is significant.
- Sufficient for Most Practical Designs. Suitable for impedance control, loss estimation, and S-parameter generation in the vast majority of PCB transmission line applications.

### Limitations

- Higher-order modes, dispersion, and cutoff behavior are not modeled. Structures that support non-TEM modes (e.g., waveguides) are outside the solver’s validity range.
- Current is modeled at the surface for AC and DC resistance is blended smoothly at low frequencies. Full 2D/3D current density inside conductors is not solved, which can reduce accuracy at frequencies where skin depth is comparable to conductor thickness.
- At DC and low frequencies (~<1 MHz), return current spreads over the ground plane. Since the solver infers inductance from capacitance, partial inductance and finite ground width effects may be inaccurate.
- Results apply to uniform, infinitely long transmission lines. Bends can often be approximated to behave similarly to straight line if they curve smoothly.
- Radiation is not modeled.

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

WebAssembly (WASM) module is used for high-performance sparse matrix solving. This module is built using Emscripten and the Eigen C++ library.
Compiled js and wasm is already included. Compiling is only needed if changes to
WASM module are necessary.

#### Prerequisites

**Emscripten Compiler (emcc)**

Install the Emscripten SDK:

```bash
# Clone the emsdk repository
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk

# Install and activate the latest SDK
./emsdk install latest
./emsdk activate latest

# Add to PATH (add this to your .bashrc or .zshrc for permanent use)
source ./emsdk_env.sh
```

Verify installation:
```bash
emcc --version
```

**Eigen Library**

The Eigen library is included as a git submodule. Initialize it:

```bash
git submodule update --init --recursive
```

#### Build Steps

```bash
cd src/wasm_solver
make
```

This compiles `solver.cpp` using Emscripten and generates:
- `solver.js` - JavaScript wrapper for the WASM module
- `solver.wasm` - Compiled WebAssembly binary
