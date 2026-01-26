#!/usr/bin/env bash
# Optional build script to minify files
set -e

rm -rf dist
mkdir dist
mkdir dist/wasm_solver

#make -C src/wasm_solver/

# JS
npx esbuild src/app.js \
  --bundle \
  --minify \
  --format=esm \
  --sourcemap \
  --platform=node \
  --outfile=dist/app.js

# CSS
npx postcss src/style.css \
  --use cssnano \
  -o dist/style.css

# HTML
cp src/field_solver.html dist/
cp src/plotly-3.3.0.min.js dist/
cp src/wasm_solver/solver.js dist/wasm_solver/solver.js
cp src/wasm_solver/solver.wasm dist/solver.wasm
