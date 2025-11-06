//! # Orphos Gene Finder - Rust Implementation
//!
//! A high-performance Rust implementation of the Orphos prokaryotic gene finding algorithm.
//! This library provides both single genome and metagenomic gene prediction capabilities with
//! support for parallel processing.
//!
//! ## Overview
//!
//! Orphos (Prokaryotic Dynamic Programming Gene-finding Algorithm) is an unsupervised machine
//! learning method for finding genes in prokaryotic genomes. This Rust implementation maintains
//! compatibility with the original C version while offering improved performance and safety.
//!
//! ## Features
//!
//! - **Single Genome Mode**: Train on a complete genome for optimal gene prediction
//! - **Metagenomic Mode**: Predict genes in fragmented or mixed sequences
//! - **Multiple Output Formats**: Support for GenBank, GFF, GCA, and SCO formats
//! - **Parallel Processing**: Multi-threaded execution using Rayon
//! - **Type Safety**: Compile-time guarantees for training states
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use orphos_core::{OrphosAnalyzer, config::OrphosConfig};
//!
//! // Create analyzer with default configuration
//! let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
//!
//! // Analyze a genome sequence
//! let results = analyzer.analyze_sequence(
//!     "ATGCGATCGATCG...",
//!     Some("MyGenome".to_string())
//! )?;
//!
//! println!("Found {} genes", results.genes.len());
//! # Ok::<(), orphos_core::types::OrphosError>(())
//! ```
//!
//! ## Architecture
//!
//! The library uses a type-state pattern to ensure training is performed before gene prediction:
//!
//! ```rust,no_run
//! use orphos_core::engine::{UntrainedOrphos, Orphos, Untrained};
//! use orphos_core::config::OrphosConfig;
//! use orphos_core::sequence::encoded::EncodedSequence;
//!
//! // Create an untrained analyzer
//! let mut untrained = UntrainedOrphos::with_config(OrphosConfig::default())?;
//!
//! // Encode the sequence
//! let encoded = EncodedSequence::without_masking(b"ATGCGATCGATCG...");
//!
//! // Train on the genome (type changes to TrainedOrphos)
//! let trained = untrained.train_single_genome(&encoded)?;
//!
//! // Use the higher-level API to find genes
//! use orphos_core::OrphosAnalyzer;
//! let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
//! let results = analyzer.analyze_sequence("ATGCGATCGATCG...", None)?;
//! println!("Found {} genes", results.genes.len());
//! # Ok::<(), orphos_core::types::OrphosError>(())
//! ```
//!
//! ## Module Organization
//!
//! - [`config`]: Configuration options for analysis
//! - [`engine`]: Main analysis engine and training logic
//! - [`types`]: Core data types and structures
//! - [`results`]: Gene prediction results
//! - [`sequence`]: Sequence encoding and manipulation
//! - [`training`]: Training algorithms for gene models
//! - [`node`]: Gene node management and scoring
//! - [`algorithms`]: Core gene-finding algorithms
//! - [`output`]: Output formatting for various file types
//! - [`bitmap`]: Efficient sequence encoding utilities
//!
//! ## Output Formats
//!
//! The library supports multiple output formats configured via [`config::OutputFormat`]:
//!
//! - **GenBank (GBK)**: Rich feature annotation format
//! - **GFF3**: General Feature Format version 3
//! - **GCA**: Gene coordinate annotation
//! - **SCO**: Simple coordinate output
//!
//! ## Error Handling
//!
//! All fallible operations return [`Result<T, OrphosError>`](types::OrphosError),
//! providing detailed error information for:
//!
//! - Invalid sequences (too short, invalid characters)
//! - I/O errors during file operations
//! - Training failures
//! - Configuration errors

pub mod algorithms;
pub mod bitmap;
pub mod config;
pub mod constants;
pub mod engine;
pub mod metagenomic;
pub mod node;
pub mod output;
pub mod results;
pub mod sequence;
pub mod training;
pub mod types;

pub use engine::OrphosAnalyzer;

use crate::{node::rbs_score, types::OrphosError};
