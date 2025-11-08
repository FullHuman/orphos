#!/bin/bash
set -e

# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source $HOME/.cargo/env

# Install wasm-pack
curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh

# Build the WASM module
echo "Compiling Rust to WebAssembly..."
wasm-pack build --target web --out-dir www/pkg
npm run build:release
