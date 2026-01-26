#!/usr/bin/env bash
# Optional build script to minify files
set -e

rm -rf dist
mkdir dist

#make -C src/wasm_solver/

# JS
npx esbuild src/app_solver.js \
  --bundle \
  --minify \
  --format=esm \
  --sourcemap \
  --platform=node \
  --outfile=dist/app_solver.js

# CSS
npx postcss src/solver-style.css \
  --use cssnano \
  -o dist/solver-style.css

# HTML
cp src/field_solver.html dist/
cp src/plotly-3.3.0.min.js dist/
cp src/wasm_solver/solver.wasm dist/solver.wasm
