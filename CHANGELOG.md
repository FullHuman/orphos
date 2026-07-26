# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-07-26

### Added

- Parallel processing of multi-record FASTA files using Rayon.
- BED output support in the WebAssembly interface.
- RefSeq CDS validation fixtures for *Escherichia coli*, *Helicobacter pylori*,
  *Bacillus subtilis*, and *Salmonella enterica*.
- A single-threaded CLI benchmark and a contributor guide with a maintainer
  release checklist.

### Changed

- `OrphosAnalyzer` analysis methods now use shared references, so callers no
  longer need mutable analyzer access.
- Optimized nucleotide encoding, reverse complements, GC and k-mer
  calculations, motif scoring, and display-score lookup.
- Redesigned the web interface with improved accessibility, responsive layout,
  focus states, animations, and maintainable CSS custom properties.
- Updated the Rust dependency graph, including `bio` 4, `pyo3` 0.29, and the
  latest compatible releases of the remaining dependencies.
- Declared Rust 1.89 as the minimum supported version and updated CI, release,
  crates.io, PyPI, and npm workflows to current GitHub Actions releases.

### Fixed

- Corrected encoded-sequence operations, producing exact Prodigal output
  matches for two additional reference genomes.
- Made Rayon global thread-pool initialization safe across repeated analyzer
  construction.
- Corrected the Node.js package name to `@fullhuman/orphos` and included the
  GPL-3.0-or-later license in the WebAssembly package.

## [0.2.0] - 2026-04-11

### Added

- Support for circular sequences (https://github.com/FullHuman/orphos/issues/5)
- Support for BED output (https://github.com/FullHuman/orphos/issues/4)

[0.3.0]: https://github.com/FullHuman/orphos/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/FullHuman/orphos/releases/tag/v0.2.0
