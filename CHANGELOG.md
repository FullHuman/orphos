# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Gene prediction on a single genome reuses the nodes training created and
  sorted, with their scores reset, instead of creating and sorting them a
  second time; nodes are sorted in place instead of being gathered into a
  second array; the node vector is sized from the sequence length; and
  sequences are encoded with the vector encoder when runs of `N` are masked
  too. The nodes of the two strands are created on two threads while the
  GC frame plot is computed on a third, and the three passes of the raw
  coding score, the GC bias record, the ORF GC content and the scans of the
  Shine-Dalgarno start training run per strand on two threads, each strand
  through its own view of the shared nodes so that nothing is copied or
  written back. On a pool of one thread the same passes run in place on the
  calling thread, as fast as before or faster. Output is unchanged; E. coli
  K-12 is predicted in about 0.6 s instead of 0.9 s on four cores.

## [0.4.0] - 2026-10-11

### Added

- `TrainedOrphos::training()` to access the trained model.
- Reference fixtures for *Mycoplasma genitalium*, *Rickettsia prowazekii* and
  *Nanoarchaeum equitans*, and an `ORPHOS_PRODIGAL_BIN` override for the
  equivalence tests.

### Changed

- Predictions are now byte-identical to Prodigal (GitHub master) in single and
  metagenomic mode, verified on 41 genomes. Single-genome output had differed
  on a few genes per genome: the GC frame bias counted the last bases of a
  sequence differently, the coding length factor used `powi` instead of
  `pow`, ORF GC content ignored a trailing partial codon, start tweaking
  ranked alternatives differently from `tweak_final_starts`, and a start
  codon failing the length check skipped the edge-gene test. Expressions the
  reference Prodigal build evaluates with fused multiply-adds now use
  `mul_add`, so near-ties resolve identically. Metagenomic output follows
  Prodigal's master branch and pyrodigal; the 2.6.3 release applies its
  short-gene penalty to every gene because of a typo fixed upstream.
- Single-genome prediction of *E. coli* K-12 runs about 6.5x faster
  single-threaded and metagenomic mode about 8.5x faster. The
  dynamic-programming loop visits only compatible source nodes, resolves
  constant-score connection types with sliding-window maxima and scores the
  rest on a compact node view; node creation, node sorting, Shine-Dalgarno
  scoring, the coding, ORF GC, GC-frame and background k-mer passes and start
  training were rewritten around lookup tables and rolling k-mers.
- Metagenomic mode computes the ORF GC content once per node set instead of
  once per preset.

### Fixed

- `-g`/`--translation-table` (and `translation_table` in the library and
  Python configuration) was parsed but never used; training always ran with
  table 11. The requested table is now used, matching `prodigal -g`.
- `-c`, `-r`, `-m` and `-q` are plain flags again, as documented; they no
  longer require a value.
- Non-SD training and metagenomic mode no longer overflow the 1 MiB
  main-thread stack on Windows: half-megabyte motif tables are no longer
  built on the stack, and the CLI runs its analysis on a thread with a
  64 MiB stack.

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

[0.4.0]: https://github.com/FullHuman/orphos/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/FullHuman/orphos/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/FullHuman/orphos/releases/tag/v0.2.0
