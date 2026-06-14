# 2D Transmission Line Field Solver

![Header](https://github.com/Ttl/js_2d_fields/blob/master/docs/header.png?raw=true)

A browser-based quasi-static and full-wave 2D field solver for transmission line analysis. Computes characteristic impedance, effective permittivity, RLGC parameters, losses, and S-parameters.

Try it online: https://hforsten.com/field_solver.html

## Features

- Transmission Line Types: Microstrip, stripline, GCPW (single-ended and differential)
- Electric Field Solving: Calculates characteristic impedance
- Full RLGC Extraction: Resistance, inductance, capacitance, conductance per unit length
- Loss Modeling: Conductor losses (skin effect, surface roughness) and dielectric losses
- S-Parameter Export: Touchstone .s2p and .s4p file generation
- Visualization: 2D potential plots, E-field streamlines, frequency-dependent plots
- Adaptive Meshing: Automatic mesh refinement for accurate field solutions

## Quick Start

1. Host src folder on a web server (for example `python -m http.server 8000`). Open `src/field_solver.html` in a browser.

### Validity

- Designed for microstrip, stripline, and coplanar waveguide structures commonly used in PCB RF and high-speed digital designs.
- Provides good results when the transmission line supports TEM or quasi-TEM propagation, which covers most practical PCB geometries below the onset of higher-order modes.
- Results have been checked against EM solver and actual measurement data with different geometries and transmission line types, showing close agreement with small error in typical use cases (`tests` folder).
- Accurate from RF through microwave and high-speed digital frequencies where return currents are confined and skin effect is significant.
- Sufficient for Most Practical Designs. Suitable for impedance control, loss estimation, and S-parameter generation in the vast majority of PCB transmission line applications.
- Dispersion and higher-order modes are support with full-wave solver.

## Common Tasks

### Run Tests

```bash
node tests/test_vs_ref.js
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

**gmsh**

Needed for full-wave solver mesh generation.

**Spectra**

Needed for full-wave solver mesh generation.

#### Build Steps

```bash
cd src/wasm_solver
make
```

This compiles `solver.cpp` using Emscripten and generates:
- `solver.js` - JavaScript wrapper for the WASM module
- `solver.wasm` - Compiled WebAssembly binary
