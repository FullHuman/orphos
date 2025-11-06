#!/bin/bash
set -e

echo "Building Orphos WebAssembly module..."

# Check if wasm-pack is installed
if ! command -v wasm-pack &> /dev/null; then
    echo "Error: wasm-pack is not installed."
    echo "Install it with: cargo install wasm-pack"
    exit 1
fi

# Build the WASM module
echo "Compiling Rust to WebAssembly..."
wasm-pack build --target web --out-dir www/pkg

echo "✓ Build complete!"
echo ""
echo "To start the development server, run:"
echo "  cd prodigal-wasm && npm run serve"
echo ""
echo "Then open http://localhost:8080 in your browser"
