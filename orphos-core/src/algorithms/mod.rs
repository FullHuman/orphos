//! Core gene-finding algorithms.
//!
//! This module contains the fundamental algorithms used for prokaryotic gene prediction:
//!
//! ## Modules
//!
//! - [`dynamic_programming`]: Optimal gene path selection using dynamic programming
//! - [`gene_finding`]: Gene construction and annotation from prediction nodes
//!
//! ## Algorithm Overview
//!
//! Gene prediction follows a multi-stage pipeline:
//!
//! 1. **Node Generation**: Identify all potential start/stop codons
//! 2. **Node Scoring**: Calculate quality scores for each potential gene
//!    - Coding potential (dicodon frequencies)
//!    - Start codon preference (ATG vs GTG vs TTG)
//!    - Ribosome binding site strength
//!    - Upstream sequence composition
//! 3. **Dynamic Programming**: Select optimal non-overlapping gene set
//! 4. **Gene Refinement**: Fine-tune start positions and annotations
//!
//! ## Dynamic Programming
//!
//! The core algorithm uses dynamic programming to find the highest-scoring set
//! of non-overlapping genes. For each position, it calculates:
//!
//! ```text
//! score[i] = max(
//!     score[i-1],                    // Skip this gene
//!     score[j] + gene_score(j, i)    // Take this gene starting at j
//! )
//! ```
//!
//! This ensures optimal gene selection while respecting biological constraints.

pub mod dynamic_programming;
pub mod gene_finding;
