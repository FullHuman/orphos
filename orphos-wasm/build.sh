#!/bin/bash
set -e

# Install Rust only if it is not already installed
if ! command -v rustc >/dev/null 2>&1; then
	echo "Rust not found. Installing Rust..."
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
	source "$HOME/.cargo/env"
else
	echo "Rust is already installed."
fi

# Install wasm-pack only if it is not already installed
if ! command -v wasm-pack >/dev/null 2>&1; then
	echo "wasm-pack not found. Installing wasm-pack..."
	curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
else
	echo "wasm-pack is already installed."
fi

# Build the WASM module
echo "Compiling Rust to WebAssembly..."
wasm-pack build --target web --out-dir www/pkg
npm run build:release
