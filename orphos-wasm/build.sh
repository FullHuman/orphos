#!/bin/bash
set -e

# Prefer rustup-managed toolchain binaries when available.
export PATH="$HOME/.cargo/bin:$PATH"

# Ensure rustup is available (required to add the wasm target).
if ! command -v rustup >/dev/null 2>&1; then
	echo "rustup not found. Installing rustup..."
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
	source "$HOME/.cargo/env"
else
	echo "rustup is already installed."
fi

# Ensure Rust compiler is available.
if ! command -v rustc >/dev/null 2>&1; then
	echo "rustc not found. Installing stable Rust toolchain..."
	rustup toolchain install stable --profile minimal
	rustup default stable
else
	echo "Rust is already installed."
fi

# Ensure the WebAssembly target is available.
if ! rustup target list --installed | grep -q '^wasm32-unknown-unknown$'; then
	echo "Installing Rust target wasm32-unknown-unknown..."
	rustup target add wasm32-unknown-unknown
else
	echo "Rust target wasm32-unknown-unknown is already installed."
fi

# Install wasm-pack only if it is not already installed
if ! command -v wasm-pack >/dev/null 2>&1; then
	echo "wasm-pack not found. Installing wasm-pack..."
	cargo install wasm-pack --locked
else
	echo "wasm-pack is already installed."
fi

# Build the WASM module
echo "Compiling Rust to WebAssembly..."
wasm-pack build --target web --out-dir www/pkg
npm run build:release
