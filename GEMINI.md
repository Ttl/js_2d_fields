# Javascript 2D transmission line field solver

Solves for transmission line electric potential with Laplace solver, then
calculates electric field from potential. Effective permittivity is calculated
by solving the mesh once with dielectrics and second time without dielectrics,
and dividing the calculated capacitance: er_eff = C/C_0. Z0 is also calculated
from capacitances.

Dielectric loss is calculated by integrating loss in dielectrics. Conductor loss
is solved with perturbation method (TODO).

The reference implementation is a Python program "field_solver_ref.py". The
Javascript solver output should match it.

The goal is to make browser based 2D transmission line solver that solves Z0,
eps_eff, losses and RLGC parameters.

Web app is in app.js. Front end is "field_solver.html".

# Run command line tests

$ node test_microstrip.js

# WASM solver

wasm_solver directory has WASM matrix solver that uses Eigen. Use make to
compile it. C ++ source is in solver.cpp file.
