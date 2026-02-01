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

# HTML - swap to CDN first, then minify (order matters!)
# Copy and swap CDN URL first
cp src/field_solver.html dist/field_solver.html
sed -i "s|const PRIMARY_PLOTLY_SRC = 'plotly-3.3.0.min.js';|const PRIMARY_PLOTLY_SRC = 'https://cdn.plot.ly/plotly-3.3.0.min.js';|g" dist/field_solver.html

# Then minify (which will inline the CDN URL)
npx html-minifier-terser \
  --collapse-whitespace \
  --remove-comments \
  --minify-js true \
  --minify-css true \
  dist/field_solver.html \
  -o dist/field_solver.html

# Copy plotly as fallback in case CDN is unreachable
cp src/plotly-3.3.0.min.js dist/

cp src/wasm_solver/solver.wasm dist/solver.wasm
