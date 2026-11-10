# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `OrphosAnalyzer::analyze_sequences` analyzes the records of one input in
  order, numbering them; `analyze_fasta_file` and the Python and
  WebAssembly bindings go through it.
- `SequenceInfo::number`, the number of a sequence in its input, counted
  from 1. Code that builds a `SequenceInfo` has to set it.

### Fixed

- The sequences of a multi-sequence input are numbered as Prodigal numbers
  them. Every sequence used to be number 1, so the genes of every sequence
  were `1_1`, `1_2`, ... and every GFF sequence header said `seqnum=1`; the
  first gene of the second sequence is now `2_1` and its header `seqnum=2`.
  This shows in GFF output and in the name column of BED output. Output for
  a single sequence is unchanged.
- GFF output of several sequences has one `##gff-version` line, at the top,
  instead of one per sequence.
- With `-m`, a run of `N` that lasts to the end of a sequence is not masked,
  as in Prodigal, which records a run when it reads the base after it;
  genes reaching into such a run were dropped.
- A start codon of the reverse strand within twenty bases of the end of the
  sequence is given the Shine-Dalgarno motifs Prodigal finds for it;
  Prodigal also searches from positions past the end there, and such a
  search can find a motif that no other does.
- In metagenomic mode a sequence in which no preset model finds a gene is
  reported under the first model tried for its GC content, as in Prodigal,
  instead of under model 0.
- Tests compare the GFF output with files written by a build of Prodigal's
  master branch, line by line, for inputs of many short sequences, runs of
  `N`, closed ends, masking and translation table 4
  (`tests/prodigal_reference_outputs.rs`).

### Changed

- **Breaking:** `OrphosResults::training_used` is an `Arc<Training>`
  instead of a `Training`. A model is half a megabyte and every result held
  its own copy, so 400 contigs held 250 MB and a large metagenome could not
  fit in memory; in metagenomic mode the results of a preset now share its
  model (24 MB for the same 400 contigs). Reading the model is unchanged;
  code that builds an `OrphosResults` wraps the model in an `Arc`, and code
  that changes it goes through `Arc::make_mut`.
  `metagenomic::get_preset_training_shared` returns a preset's shared model.
- Metagenomic mode no longer rebuilds the Shine-Dalgarno motif table for
  every preset model of every sequence; the table depends on the model
  alone and is now built once per model, from the 64 patterns of matching
  bases instead of the 4096 6-mers. The same holds for the gene length
  factors of a model and the codon table of a genetic code, and the motifs
  upstream of the start codons are read once per sequence rather than once
  per model. A metagenome of 400 contigs of 1 to 8 kb is predicted in 0.29 s
  instead of 4.9 s on one thread, and in 0.05 s instead of 0.63 s on ten.
- In metagenomic mode the preset models of a long sequence are scored on
  several threads, each on its own copy of the nodes, with at most 256 MB of
  copies at a time over all sequences: *E. coli* K-12 takes 0.55 s instead
  of 2.2 s on ten threads, and 530 MB instead of 245 MB. Sequences of fewer
  than 4096 nodes, and runs on one thread, are scored as before.
- Single-genome prediction is about twice as fast, on one thread and on
  many (*E. coli* K-12: 0.22 s instead of 0.39 s on one thread, 0.15 s
  instead of 0.29 s on ten), and uses less memory. The dynamic programming
  visits only the sources that can connect to a target: the latest ones for
  the connections that reach back a bounded distance, the reverse stops
  whose open reading frame is still open, and for the 3'-to-3' connections
  across strands the few sources that can overlap a gene, the rest being
  resolved by a sliding-window maximum. Node positions are created and
  sorted on their own, and every node is then written once; the strand
  scans visit only start and stop codons; the raw coding scores computed
  for training are kept for the prediction instead of being computed
  twice; the GC frame plot is kept as running totals per codon, so the GC
  bias of a node is a subtraction instead of a walk over its codons; the
  ORF GC content counts bases 32 at a time; and the coding score reads one
  codon per step and remembers the length factor of each gene length.
- Start training without Shine-Dalgarno motifs is about three times as
  fast, so that *M. genitalium* is predicted in 0.047 s instead of 0.11 s:
  it searched the motifs of every start node twice per iteration, and now
  once, on compact copies of the nodes and of the sequence windows the
  searches read.
- None of the above changes the output.
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
