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
cp src/.htaccess dist/.htaccess

# Cache-busting: append content hashes as query strings so browsers fetch
# new versions when files change. The HTML must be served with no-cache headers
# (e.g. Cache-Control: no-cache) so users always get the latest asset URLs.
BUILD_DATE=$(date -u '+%Y-%m-%d')
GIT_HASH=$(git rev-parse --short HEAD)
sed -i "s|BUILD_VERSION|${BUILD_DATE} (${GIT_HASH})|g" dist/field_solver.html

JS_HASH=$(sha256sum dist/app_solver.js | cut -c1-8)
CSS_HASH=$(sha256sum dist/solver-style.css | cut -c1-8)
WASM_HASH=$(sha256sum dist/solver.wasm | cut -c1-8)

# Inject WASM hash into bundled JS (solver.wasm is fetched at runtime)
sed -i "s|\"solver\.wasm\"|\"solver.wasm?v=${WASM_HASH}\"|g" dist/app_solver.js

# Inject JS and CSS hashes into HTML
sed -i \
  -e "s|\"app_solver\.js\"|\"app_solver.js?v=${JS_HASH}\"|g" \
  -e "s|\"solver-style\.css\"|\"solver-style.css?v=${CSS_HASH}\"|g" \
  dist/field_solver.html
