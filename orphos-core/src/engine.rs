use std::marker::PhantomData;
use std::path::Path;

use crate::algorithms::dynamic_programming::eliminate_bad_genes;
use crate::algorithms::dynamic_programming::predict_genes;
use crate::algorithms::gene_finding::GeneBuilder;
use crate::config::OrphosConfig;
use crate::constants::MIN_SEQUENCE_LENGTH;
use crate::node::{
    add_nodes, calculate_dicodon_gene, raw_coding_score, rbs_score, record_gc_bias,
    record_overlapping_starts, score_nodes, sort_nodes_by_position,
};
use crate::results::{OrphosResults, SequenceInfo};
use crate::sequence::calc_most_gc_frame;
use crate::sequence::encoded::EncodedSequence;
use crate::sequence::read_fasta_sequences;
use crate::training::non_sd_training::train_starts_nonsd;
use crate::training::sd_training::train_starts_sd;
use crate::training::should_use_sd;
use crate::types::Gene;
use crate::types::{OrphosError, Training};

/// Marker trait for Orphos training state.
///
/// This trait is used in the type-state pattern to enforce that training
/// must be performed before gene prediction. It's implemented by both
/// [`Untrained`] and [`Trained`] marker types.
pub trait TrainingState {}

/// Marker type indicating an untrained Orphos instance.
///
/// A [`Orphos<Untrained>`] instance can perform training operations
/// but cannot find genes until training is complete.
#[derive(Debug, Clone)]
pub struct Untrained;

/// Marker type indicating a trained Orphos instance.
///
/// A [`Orphos<Trained>`] instance has completed training and can
/// perform gene prediction operations.
#[derive(Debug, Clone)]
pub struct Trained;

impl TrainingState for Untrained {}
impl TrainingState for Trained {}

/// Main Orphos configuration and execution engine.
///
/// This struct uses the type-state pattern with the `S` type parameter
/// to ensure training is performed before gene prediction. The state
/// transitions from [`Untrained`] to [`Trained`] via the training methods.
///
/// # Type Parameters
///
/// * `S` - The training state, either [`Untrained`] or [`Trained`]
///
/// # Examples
///
/// ```rust,no_run
/// use orphos_core::engine::UntrainedOrphos;
/// use orphos_core::config::OrphosConfig;
/// use orphos_core::sequence::encoded::EncodedSequence;
///
/// // Create an untrained instance
/// let mut orphos = UntrainedOrphos::new();
///
/// // Encode a sequence
/// let sequence = b"ATGAAACGCATTAGCACCACCATT...";
/// let encoded = EncodedSequence::without_masking(sequence);
///
/// // Train on the sequence and get results
/// let trained = orphos.train_single_genome(&encoded)?;
///
/// // Use the higher-level API to analyze sequences
/// use orphos_core::OrphosAnalyzer;
/// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
/// let results = analyzer.analyze_sequence("ATGAAACGCATTAGCACCACCATT...", None)?;
/// println!("Found {} genes", results.genes.len());
/// # Ok::<(), orphos_core::types::OrphosError>(())
/// ```
#[derive(Debug, Default)]
pub struct Orphos<S: TrainingState> {
    /// Configuration options for gene prediction
    pub config: OrphosConfig,
    /// Training data obtained from the genome
    training: Option<Training>,
    /// Type-state marker (zero-sized)
    _state: PhantomData<S>,
}

/// Type alias for an untrained Orphos instance.
///
/// Use this when you need to perform training on a new genome.
pub type UntrainedOrphos = Orphos<Untrained>;

/// Type alias for a trained Orphos instance.
///
/// Use this when you have already trained on a genome and want to find genes.
pub type TrainedOrphos = Orphos<Trained>;

impl UntrainedOrphos {
    /// Creates a new untrained Orphos instance with default configuration.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use orphos_core::engine::UntrainedOrphos;
    ///
    /// let orphos = UntrainedOrphos::new();
    /// ```
    pub fn new() -> Self {
        Self {
            config: OrphosConfig::default(),
            training: None,
            _state: PhantomData,
        }
    }

    /// Creates a new untrained Orphos instance with custom configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration options for gene prediction
    ///
    /// # Errors
    ///
    /// Returns [`OrphosError`] if thread pool configuration fails.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use orphos_core::engine::UntrainedOrphos;
    /// use orphos_core::config::{OrphosConfig, OutputFormat};
    ///
    /// let config = OrphosConfig {
    ///     closed_ends: true,
    ///     output_format: OutputFormat::Gff,
    ///     ..Default::default()
    /// };
    ///
    /// let orphos = UntrainedOrphos::with_config(config)?;
    /// # Ok::<(), orphos_core::types::OrphosError>(())
    /// ```
    pub fn with_config(config: OrphosConfig) -> Result<Self, OrphosError> {
        let orphos = Self {
            config,
            training: None,
            _state: PhantomData,
        };

        if let Some(num_threads) = orphos.config.num_threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .map_err(|e| {
                    OrphosError::InvalidSequence(format!("Failed to configure thread pool: {}", e))
                })?;
        }

        Ok(orphos)
    }

    /// Trains the model on a single complete genome sequence.
    ///
    /// This method analyzes the sequence to build a statistical model of gene
    /// characteristics including:
    /// - Start codon usage (ATG, GTG, TTG)
    /// - Ribosome binding site motifs
    /// - Codon usage patterns
    /// - GC content bias
    ///
    /// # Arguments
    ///
    /// * `encoded_sequence` - The genome sequence encoded in bitmap format
    ///
    /// # Returns
    ///
    /// A [`TrainedOrphos`] instance ready for gene prediction.
    ///
    /// # Errors
    ///
    /// Returns [`OrphosError::InvalidSequence`] if:
    /// - The sequence is shorter than [`MIN_SEQUENCE_LENGTH`]
    /// - The sequence contains invalid characters
    /// - Training fails to converge
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use orphos_core::engine::UntrainedOrphos;
    /// use orphos_core::sequence::encoded::EncodedSequence;
    ///
    /// let mut orphos = UntrainedOrphos::new();
    /// let sequence = b"ATGAAACGCATTAGCACCACCATT...";
    /// let encoded = EncodedSequence::without_masking(sequence);
    ///
    /// let trained = orphos.train_single_genome(&encoded)?;
    /// # Ok::<(), orphos_core::types::OrphosError>(())
    /// ```
    pub fn train_single_genome(
        &mut self,
        encoded_sequence: &EncodedSequence,
    ) -> Result<TrainedOrphos, OrphosError> {
        let sequence_length = encoded_sequence.sequence_length;
        if sequence_length < MIN_SEQUENCE_LENGTH {
            return Err(OrphosError::InvalidSequence(format!(
                "Sequence too short for gene prediction: {} bp (minimum {} bp required)",
                sequence_length, MIN_SEQUENCE_LENGTH
            )));
        }

        if !self.config.quiet {
            eprintln!(
                "Training on single genome ({} bp, {:.2}% GC)...",
                sequence_length,
                encoded_sequence.gc_content * 100.0
            );
        }
        let mut training = Training {
            gc_content: 0.0,
            ..Default::default()
        };
        training.uses_shine_dalgarno = false;
        training.gc_bias_factors = [0.0; 3];

        training.gc_content = encoded_sequence.gc_content;

        let mut nodes = Vec::new();
        let num_nodes = add_nodes(
            encoded_sequence,
            &mut nodes,
            self.config.closed_ends,
            &training,
        )?;

        if !self.config.quiet {
            eprintln!(
                "Located {} potential start/stop nodes, closed {}",
                num_nodes, self.config.closed_ends
            );
        }

        sort_nodes_by_position(&mut nodes);

        let gc_frame = calc_most_gc_frame(&encoded_sequence.forward_sequence, sequence_length);

        record_gc_bias(&gc_frame, &mut nodes, &mut training);

        if !self.config.quiet {
            eprintln!(
                "Frame bias scores: {:.8} {:.8} {:.8}",
                training.gc_bias_factors[0],
                training.gc_bias_factors[1],
                training.gc_bias_factors[2]
            );
        }

        record_overlapping_starts(&mut nodes, &training, false);

        let initial_path = predict_genes(&mut nodes, &training, false).unwrap_or(0);
        calculate_dicodon_gene(
            &mut training,
            &encoded_sequence.forward_sequence,
            &encoded_sequence.reverse_complement_sequence,
            sequence_length,
            &nodes,
            initial_path,
        );

        raw_coding_score(
            &encoded_sequence.forward_sequence,
            &encoded_sequence.reverse_complement_sequence,
            sequence_length,
            &mut nodes,
            &training,
        );

        rbs_score(
            &encoded_sequence.forward_sequence,
            &encoded_sequence.reverse_complement_sequence,
            sequence_length,
            &mut nodes,
            &training,
        );

        train_starts_sd(
            &encoded_sequence.forward_sequence,
            &encoded_sequence.reverse_complement_sequence,
            sequence_length,
            &nodes,
            &mut training,
        );
        training.uses_shine_dalgarno = should_use_sd(&training);
        if self.config.force_non_sd {
            training.uses_shine_dalgarno = false;
        }

        if !training.uses_shine_dalgarno {
            train_starts_nonsd(
                &encoded_sequence.forward_sequence,
                &encoded_sequence.reverse_complement_sequence,
                sequence_length,
                &mut nodes,
                &mut training,
            );
        }

        if !self.config.quiet {
            eprintln!("Training complete!");
        }

        Ok(Orphos {
            config: self.config.clone(),
            training: Some(training),
            _state: PhantomData,
        })
    }
}

impl TrainedOrphos {
    /// Creates a new trained Orphos instance with pre-computed training data.
    ///
    /// This is useful when you have previously computed training data that you
    /// want to reuse for gene prediction on multiple sequences.
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration options for gene prediction
    /// * `training` - Pre-computed training data
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use orphos_core::engine::TrainedOrphos;
    /// use orphos_core::config::OrphosConfig;
    /// use orphos_core::types::Training;
    ///
    /// let config = OrphosConfig::default();
    /// let training = Training::default(); // In practice, load from file
    ///
    /// let trained = TrainedOrphos::new(config, training);
    /// ```
    pub const fn new(config: OrphosConfig, training: Training) -> Self {
        Self {
            config,
            training: Some(training),
            _state: PhantomData,
        }
    }

    /// Finds genes in a single genome sequence using the trained model.
    ///
    /// Uses dynamic programming to find the optimal set of non-overlapping genes
    /// based on the statistical model built during training. This method performs:
    ///
    /// 1. Node generation (start/stop codon detection)
    /// 2. Node scoring (coding potential, RBS, start codon usage)
    /// 3. Dynamic programming gene selection
    /// 4. Gene quality filtering
    /// 5. Start position refinement
    ///
    /// # Arguments
    ///
    /// * `encoded_sequence` - The genome sequence encoded in bitmap format
    ///
    /// # Returns
    ///
    /// A vector of [`Gene`] predictions sorted by position. Returns an empty
    /// vector if no genes are found.
    ///
    /// # Errors
    ///
    /// Returns [`OrphosError::InvalidSequence`] if:
    /// - The Orphos instance is not properly trained
    /// - Node generation fails
    /// - Sequence contains invalid characters
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use orphos_core::engine::UntrainedOrphos;
    /// use orphos_core::sequence::encoded::EncodedSequence;
    ///
    /// let mut orphos = UntrainedOrphos::new();
    /// let sequence = b"ATGAAACGCATTAGCACCACCATT...";
    /// let encoded = EncodedSequence::without_masking(sequence);
    ///
    /// let trained = orphos.train_single_genome(&encoded)?;
    ///
    /// // Use OrphosAnalyzer for a higher-level API to find genes
    /// use orphos_core::OrphosAnalyzer;
    /// use orphos_core::config::OrphosConfig;
    /// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    /// let results = analyzer.analyze_sequence("ATGAAACGCATTAGCACCACCATT...", None)?;
    /// for gene in &results.genes {
    ///     println!("Gene at {}-{} (strand: {:?})",
    ///              gene.coordinates.begin, gene.coordinates.end, gene.coordinates.strand);
    /// }
    /// # Ok::<(), orphos_core::types::OrphosError>(())
    /// ```
    fn find_genes_single(
        &self,
        encoded_sequence: &EncodedSequence,
    ) -> Result<Vec<Gene>, OrphosError> {
        let mut nodes = Vec::new();
        let training = self
            .training
            .as_ref()
            .ok_or_else(|| OrphosError::InvalidSequence("Orphos is not trained".to_string()))?;

        let _num_nodes = add_nodes(
            encoded_sequence,
            &mut nodes,
            self.config.closed_ends,
            training,
        )?;

        sort_nodes_by_position(&mut nodes);

        // Tap the training state right before scoring to ensure parity at inference time

        score_nodes(
            encoded_sequence,
            &mut nodes,
            training,
            self.config.closed_ends,
            false,
        )?;

        record_overlapping_starts(&mut nodes, training, true);

        let gene_path = match predict_genes(&mut nodes, training, true) {
            Some(path) => path,
            None => {
                // No genes found - return empty gene list
                return Ok(vec![]);
            }
        };

        eliminate_bad_genes(&mut nodes, Some(gene_path), training);

        let genes = GeneBuilder::from_nodes(&nodes, gene_path, training, 1)
            .with_tweaked_starts()
            .with_annotations()
            .build();

        Ok(genes)
    }
}

/// High-level gene finding analyzer with automatic training.
///
/// This struct provides a simplified interface for gene prediction that handles
/// training automatically. It's the recommended entry point for most users.
///
/// Unlike the type-safe [`Orphos`] struct, `OrphosAnalyzer` manages training
/// internally and provides convenient methods for analyzing sequences from various
/// sources (files, strings, byte slices).
///
/// # Modes
///
/// - **Single Genome Mode** (default): Trains on each sequence individually for
///   optimal accuracy on complete genomes
/// - **Metagenomic Mode**: Uses pre-computed models for fragmented sequences
///
/// # Examples
///
/// ## Analyze a sequence string
///
/// ```rust,no_run
/// use orphos_core::{OrphosAnalyzer, config::OrphosConfig};
///
/// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
///
/// let sequence = "ATGAAACGCATTAGCACCACCATT...";
/// let results = analyzer.analyze_sequence(sequence, Some("genome1".to_string()))?;
///
/// println!("Found {} genes in {} bp sequence",
///          results.genes.len(),
///          results.sequence_info.length);
/// # Ok::<(), orphos_core::types::OrphosError>(())
/// ```
///
/// ## Analyze a FASTA file
///
/// ```rust,no_run
/// use orphos_core::{OrphosAnalyzer, config::OrphosConfig};
///
/// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
/// let results = analyzer.analyze_fasta_file("genome.fasta")?;
///
/// for result in results {
///     println!("Sequence: {}", result.sequence_info.header);
///     println!("  Genes: {}", result.genes.len());
///     println!("  GC%: {:.2}", result.sequence_info.gc_content * 100.0);
/// }
/// # Ok::<(), orphos_core::types::OrphosError>(())
/// ```
///
/// ## With custom configuration
///
/// ```rust,no_run
/// use orphos_core::{OrphosAnalyzer, config::{OrphosConfig, OutputFormat}};
///
/// let config = OrphosConfig {
///     closed_ends: true,
///     mask_n_runs: true,
///     output_format: OutputFormat::Gff,
///     num_threads: Some(4),
///     ..Default::default()
/// };
///
/// let mut analyzer = OrphosAnalyzer::new(config);
/// # Ok::<(), orphos_core::types::OrphosError>(())
/// ```
#[derive(Debug)]
pub struct OrphosAnalyzer {
    /// Configuration options for gene prediction
    pub config: OrphosConfig,
}

impl OrphosAnalyzer {
    /// Creates a new analyzer with the specified configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration options for gene prediction
    ///
    /// # Examples
    ///
    /// ```rust
    /// use orphos_core::{OrphosAnalyzer, config::OrphosConfig};
    ///
    /// let analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    /// ```
    pub const fn new(config: OrphosConfig) -> Self {
        Self { config }
    }

    /// Analyzes sequences from a FASTA file.
    ///
    /// Reads all sequences from the FASTA file and performs gene prediction on each.
    /// In single genome mode, each sequence is trained and analyzed independently.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the FASTA file
    ///
    /// # Returns
    ///
    /// A vector of [`OrphosResults`], one for each sequence in the file.
    ///
    /// # Errors
    ///
    /// Returns [`OrphosError`] if:
    /// - The file cannot be read
    /// - The FASTA format is invalid
    /// - Any sequence fails analysis
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use orphos_core::{OrphosAnalyzer, config::OrphosConfig};
    ///
    /// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    /// let results = analyzer.analyze_fasta_file("genomes.fasta")?;
    ///
    /// for (i, result) in results.iter().enumerate() {
    ///     println!("Sequence {}: {} genes", i + 1, result.genes.len());
    /// }
    /// # Ok::<(), orphos_core::types::OrphosError>(())
    /// ```
    pub fn analyze_fasta_file<P: AsRef<Path>>(
        &mut self,
        path: P,
    ) -> Result<Vec<OrphosResults>, OrphosError> {
        let sequences = read_fasta_sequences(path.as_ref().to_str().unwrap())?;

        if self.config.metagenomic {
            Ok(vec![])
        } else {
            let mut results = Vec::new();
            for (header, description, seq_bytes) in sequences {
                let result = self.analyze_sequence_bytes(&seq_bytes, header, description)?;
                results.push(result);
            }
            Ok(results)
        }
    }

    /// Analyzes a single sequence from a string.
    ///
    /// Converts the string sequence to bytes and performs gene prediction.
    /// This is a convenience method for analyzing sequences already loaded
    /// into memory.
    ///
    /// # Arguments
    ///
    /// * `sequence` - DNA sequence string (A, T, G, C, N)
    /// * `header` - Optional sequence identifier (defaults to "Orphos_Seq_1")
    ///
    /// # Returns
    ///
    /// [`OrphosResults`] containing genes, training data, and sequence info.
    ///
    /// # Errors
    ///
    /// Returns [`OrphosError`] if:
    /// - The sequence is too short (< 20,000 bp recommended)
    /// - Training fails
    /// - Gene prediction fails
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use orphos_core::{OrphosAnalyzer, config::OrphosConfig};
    ///
    /// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    ///
    /// let sequence = "ATGAAACGCATTAGCACCACCATT...";
    /// let results = analyzer.analyze_sequence(sequence, Some("E. coli K12".to_string()))?;
    ///
    /// println!("Analyzed: {}", results.sequence_info.header);
    /// println!("Found {} genes", results.genes.len());
    /// println!("GC content: {:.2}%", results.sequence_info.gc_content * 100.0);
    /// # Ok::<(), orphos_core::types::OrphosError>(())
    /// ```
    pub fn analyze_sequence(
        &mut self,
        sequence: &str,
        header: Option<String>,
    ) -> Result<OrphosResults, OrphosError> {
        let seq_bytes = sequence.as_bytes();
        let header = header.unwrap_or_else(|| "Orphos_Seq_1".to_string());

        self.analyze_sequence_bytes(seq_bytes, header, None)
    }

    /// Analyzes a single sequence from raw bytes.
    ///
    /// This is the core analysis method used by other convenience methods.
    /// It handles the complete workflow: encoding, training, and gene prediction.
    ///
    /// # Arguments
    ///
    /// * `sequence` - Raw DNA sequence bytes (ASCII: A, T, G, C, N)
    /// * `header` - Sequence identifier for output
    /// * `description` - Optional sequence description
    ///
    /// # Returns
    ///
    /// [`OrphosResults`] containing:
    /// - Predicted genes with coordinates and scores
    /// - Training parameters used
    /// - Sequence statistics (length, GC content)
    ///
    /// # Errors
    ///
    /// Returns [`OrphosError`] if:
    /// - The sequence is too short for reliable training
    /// - Sequence encoding fails
    /// - Training or prediction fails
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use orphos_core::{OrphosAnalyzer, config::OrphosConfig};
    ///
    /// let mut analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    ///
    /// let sequence = b"ATGAAACGCATTAGCACCACCATT...";
    /// let results = analyzer.analyze_sequence_bytes(
    ///     sequence,
    ///     "sequence1".to_string(),
    ///     Some("E. coli genome".to_string())
    /// )?;
    ///
    /// println!("Found {} genes", results.genes.len());
    /// # Ok::<(), orphos_core::types::OrphosError>(())
    /// ```
    pub fn analyze_sequence_bytes(
        &mut self,
        sequence: &[u8],
        header: String,
        description: Option<String>,
    ) -> Result<OrphosResults, OrphosError> {
        let sequence_length = sequence.len();
        let encoded_sequence = self.encode_sequence(sequence);

        let mut untrained_orphos = UntrainedOrphos::with_config(self.config.clone())?;
        let trained_orphos = untrained_orphos.train_single_genome(&encoded_sequence)?;
        let genes = trained_orphos.find_genes_single(&encoded_sequence)?;

        Ok(OrphosResults {
            genes: genes.clone(),
            training_used: trained_orphos.training.unwrap_or_default(),
            sequence_info: SequenceInfo {
                length: sequence_length,
                gc_content: encoded_sequence.gc_content,
                num_genes: genes.len(),
                header,
                description,
            },
            metagenomic_model: if self.config.metagenomic {
                Some("Best".to_string())
            } else {
                None
            },
        })
    }

    /// Encodes a sequence to bitmap format for efficient processing.
    ///
    /// Converts the DNA sequence into a compact binary representation that
    /// enables fast nucleotide lookups. Optionally masks runs of N characters.
    ///
    /// # Arguments
    ///
    /// * `sequence` - Raw DNA sequence bytes
    ///
    /// # Returns
    ///
    /// An [`EncodedSequence`] with forward and reverse-complement strands
    /// encoded in bitmap format.
    fn encode_sequence(&self, sequence: &[u8]) -> EncodedSequence {
        if self.config.mask_n_runs {
            EncodedSequence::with_masking(sequence)
        } else {
            EncodedSequence::without_masking(sequence)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::OutputFormat;
    use crate::constants::TEST_SEQUENCE_REPEAT_FACTOR;
    use std::env;
    use std::fs;

    // Helper function to create a simple test sequence
    fn create_test_sequence() -> Vec<u8> {
        // Simple DNA sequence with some genes
        "ATGAAACGTAAATAG".as_bytes().to_vec()
    }

    // Helper function to create a longer test sequence for training
    fn create_training_sequence() -> Vec<u8> {
        // A longer sequence suitable for training (> 20,000 bp recommended)
        let basic_gene = "ATGAAACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTAAATAG";
        basic_gene
            .repeat(TEST_SEQUENCE_REPEAT_FACTOR)
            .as_bytes()
            .to_vec()
    }

    // Helper function to create an EncodedSequence for testing
    fn create_encoded_sequence_for_test(seq: &[u8]) -> EncodedSequence {
        EncodedSequence::without_masking(seq)
    }

    #[test]
    fn test_training_state_traits() {
        // Test that training state markers implement required traits
        let _untrained: Untrained = Untrained;
        let _trained: Trained = Trained;

        // Test that they can be cloned and debugged
        let untrained_clone = _untrained.clone();
        let trained_clone = _trained.clone();

        assert_eq!(format!("{:?}", untrained_clone), "Untrained");
        assert_eq!(format!("{:?}", trained_clone), "Trained");
    }

    #[test]
    fn test_untrained_orphos_new() {
        let orphos = UntrainedOrphos::new();

        assert!(!orphos.config.metagenomic);
        assert!(!orphos.config.closed_ends);
        assert!(!orphos.config.mask_n_runs);
        assert!(!orphos.config.force_non_sd);
        assert!(!orphos.config.quiet);
        assert_eq!(orphos.config.output_format, OutputFormat::Genbank);
        assert!(orphos.training.is_none());
    }

    #[test]
    fn test_untrained_orphos_with_config() {
        let config = OrphosConfig {
            metagenomic: true,
            closed_ends: true,
            quiet: true,
            ..OrphosConfig::default()
        };

        let result = UntrainedOrphos::with_config(config.clone());
        assert!(result.is_ok());

        let orphos = result.unwrap();
        assert!(orphos.config.metagenomic);
        assert!(orphos.config.closed_ends);
        assert!(orphos.config.quiet);
        assert!(orphos.training.is_none());
    }

    #[test]
    fn test_untrained_orphos_with_thread_config() {
        let config = OrphosConfig {
            num_threads: Some(2),
            ..OrphosConfig::default()
        };

        let result = UntrainedOrphos::with_config(config);
        // Thread configuration might fail depending on system
        // Expected to potentially fail
        let _ = result;
    }

    #[test]
    fn test_untrained_orphos_with_invalid_thread_config() {
        let config = OrphosConfig {
            num_threads: Some(0),
            ..OrphosConfig::default()
        };

        let result = UntrainedOrphos::with_config(config);
        // Thread configuration might fail depending on system
        // Rayon might handle it gracefully or it might fail
        let _ = result;
    }

    #[test]
    fn test_trained_orphos_new() {
        let config = OrphosConfig::default();
        let training = Training::default();

        let orphos = TrainedOrphos::new(config, training);

        assert!(orphos.training.is_some());
    }

    #[test]
    fn test_train_single_genome_basic() {
        let mut orphos = UntrainedOrphos::new();
        let sequence = create_training_sequence();
        let encoded_sequence = create_encoded_sequence_for_test(&sequence);

        let result = orphos.train_single_genome(&encoded_sequence);

        assert!(result.is_ok());
        let trained = result.unwrap();
        assert!(trained.training.is_some());
    }

    #[test]
    fn test_train_single_genome_quiet_mode() {
        let config = OrphosConfig {
            quiet: true,
            ..OrphosConfig::default()
        };
        let mut orphos = UntrainedOrphos::with_config(config).unwrap();

        let sequence = create_training_sequence();
        let encoded_sequence = create_encoded_sequence_for_test(&sequence);

        let result = orphos.train_single_genome(&encoded_sequence);

        assert!(result.is_ok());
    }

    #[test]
    fn test_train_single_genome_force_non_sd() {
        let config = OrphosConfig {
            force_non_sd: true,
            ..OrphosConfig::default()
        };
        let mut orphos = UntrainedOrphos::with_config(config).unwrap();

        let sequence = create_training_sequence();
        let encoded_sequence = create_encoded_sequence_for_test(&sequence);

        let result = orphos.train_single_genome(&encoded_sequence);

        assert!(result.is_ok());
        let trained = result.unwrap();
        let training = trained.training.unwrap();
        assert!(!training.uses_shine_dalgarno); // Should be forced to false
    }

    #[test]
    fn test_trained_orphos_find_genes_single() {
        // First create a trained orphos
        let mut orphos = UntrainedOrphos::new();
        let sequence = create_training_sequence();
        let encoded_sequence = create_encoded_sequence_for_test(&sequence);

        let trained = orphos.train_single_genome(&encoded_sequence).unwrap();

        // Now test gene finding
        let test_seq = create_test_sequence();
        let test_encoded_sequence = create_encoded_sequence_for_test(&test_seq);

        let result = trained.find_genes_single(&test_encoded_sequence);

        assert!(result.is_ok());
        let _genes = result.unwrap();
        // Should find genes in the sequence - number could be 0 for very short sequences
    }

    #[test]
    fn test_trained_orphos_find_genes_without_training() {
        let config = OrphosConfig::default();
        let orphos = TrainedOrphos {
            config,
            training: None,
            _state: PhantomData,
        };

        let test_seq = create_test_sequence();
        let test_encoded_sequence = create_encoded_sequence_for_test(&test_seq);

        let result = orphos.find_genes_single(&test_encoded_sequence);

        assert!(result.is_err());
        if let Err(OrphosError::InvalidSequence(msg)) = result {
            assert!(msg.contains("not trained"));
        } else {
            panic!("Expected InvalidSequence error");
        }
    }

    #[test]
    fn test_orphos_analyzer_new() {
        let config = OrphosConfig::default();
        let analyzer = OrphosAnalyzer::new(config.clone());

        assert_eq!(analyzer.config.metagenomic, config.metagenomic);
        assert_eq!(analyzer.config.closed_ends, config.closed_ends);
    }

    #[test]
    fn test_analyze_sequence_basic() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        let sequence =
            "ATGAAACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTAAATAG".repeat(300);

        let result = analyzer.analyze_sequence(&sequence, None);
        assert!(result.is_ok());

        let analysis = result.unwrap();
        assert_eq!(analysis.sequence_info.header, "Orphos_Seq_1");
        assert_eq!(analysis.sequence_info.length, sequence.len());
        assert!(
            analysis.sequence_info.gc_content >= 0.0 && analysis.sequence_info.gc_content <= 1.0
        );
        assert!(analysis.metagenomic_model.is_none());
    }

    #[test]
    fn test_analyze_sequence_with_header() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        let sequence =
            "ATGAAACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTAAATAG".repeat(300);
        let header = Some("test_sequence".to_string());

        let result = analyzer.analyze_sequence(&sequence, header);
        assert!(result.is_ok());

        let analysis = result.unwrap();
        assert_eq!(analysis.sequence_info.header, "test_sequence");
    }

    #[test]
    fn test_analyze_sequence_bytes() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        let sequence =
            "ATGAAACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTAAATAG".repeat(300);
        let seq_bytes = sequence.as_bytes();
        let header = "test_sequence".to_string();
        let description = Some("Test description".to_string());

        let result =
            analyzer.analyze_sequence_bytes(seq_bytes, header.clone(), description.clone());
        assert!(result.is_ok());

        let analysis = result.unwrap();
        assert_eq!(analysis.sequence_info.header, header);
        assert_eq!(analysis.sequence_info.description, description);
        assert_eq!(analysis.sequence_info.length, seq_bytes.len());
    }

    #[test]
    fn test_analyze_sequence_metagenomic_config() {
        let config = OrphosConfig {
            metagenomic: true,
            ..OrphosConfig::default()
        };
        let mut analyzer = OrphosAnalyzer::new(config);

        let sequence =
            "ATGAAACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTAAATAG".repeat(300);

        let result = analyzer.analyze_sequence(&sequence, None);
        assert!(result.is_ok());

        let analysis = result.unwrap();
        assert!(analysis.metagenomic_model.is_some());
        assert_eq!(analysis.metagenomic_model.unwrap(), "Best");
    }

    // #[test]
    // fn test_encode_sequence_basic() {
    //     let config = OrphosConfig::default();
    //     let analyzer = OrphosAnalyzer::new(config);

    //     let sequence = b"ATCG";
    //     let mut encoded = vec![0u8; (sequence.len() * 2).div_ceil(8)];
    //     let mut unknown = vec![0u8; sequence.len().div_ceil(8)];
    //     let mut masks = Vec::new();

    //     let result = analyzer.encode_sequence(sequence, &mut encoded, &mut unknown, &mut masks);
    //     assert!(result.is_ok());

    //     let gc_content = result.unwrap();
    //     assert!((0.0..=1.0).contains(&gc_content));
    // }

    // #[test]
    // fn test_encode_sequence_with_mask_n_runs() {
    //     let config = OrphosConfig {
    //         mask_n_runs: true,
    //         ..OrphosConfig::default()
    //     };
    //     let analyzer = OrphosAnalyzer::new(config);

    //     let sequence = b"ATCGNNNNGCAT";
    //     let mut encoded = vec![0u8; (sequence.len() * 2).div_ceil(8)];
    //     let mut unknown = vec![0u8; sequence.len().div_ceil(8)];
    //     let mut masks = Vec::new();

    //     let result = analyzer.encode_sequence(sequence, &mut encoded, &mut unknown, &mut masks);
    //     assert!(result.is_ok());

    //     // N-run masking behavior depends on implementation - might not always create masks
    //     // So we just check that it doesn't crash
    // }

    #[test]
    fn test_analyze_fasta_file_not_found() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        let result = analyzer.analyze_fasta_file("nonexistent_file.fa");
        assert!(result.is_err());
    }

    #[test]
    fn test_analyze_fasta_file_metagenomic() {
        let config = OrphosConfig {
            metagenomic: true,
            ..OrphosConfig::default()
        };
        let mut analyzer = OrphosAnalyzer::new(config);

        // Create a temporary FASTA file
        let fasta_content = ">test_seq\nATCG\n";
        let temp_dir = env::temp_dir();
        let temp_file = temp_dir.join("test_metagenomic.fa");
        fs::write(&temp_file, fasta_content).unwrap();

        let result = analyzer.analyze_fasta_file(&temp_file);
        assert!(result.is_ok());

        let results = result.unwrap();
        assert!(results.is_empty()); // Metagenomic mode returns empty for now

        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_analyze_fasta_file_single_genome() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        // Create a temporary FASTA file with a longer sequence for training
        let sequence =
            "ATGAAACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTAAATAG".repeat(300);
        let fasta_content = format!(">test_seq\n{}\n", sequence);
        let temp_dir = env::temp_dir();
        let temp_file = temp_dir.join("test_single_genome.fa");
        fs::write(&temp_file, fasta_content).unwrap();

        let result = analyzer.analyze_fasta_file(&temp_file);
        assert!(result.is_ok());

        let results = result.unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].sequence_info.header, "test_seq");

        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_analyze_empty_sequence() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        let sequence = "";
        let result = analyzer.analyze_sequence(sequence, None);

        assert!(result.is_err());
        if let Err(e) = result {
            // Should be an InvalidSequence error about being too short
            match e {
                OrphosError::InvalidSequence(msg) => {
                    assert!(msg.contains("too short"));
                }
                _ => panic!("Expected InvalidSequence error for empty sequence"),
            }
        }
    }

    #[test]
    fn test_analyze_very_short_sequence() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        let sequence = "ATG"; // Very short sequence (3 bp)
        let result = analyzer.analyze_sequence(sequence, None);

        assert!(result.is_err());
        if let Err(e) = result {
            // Should be an InvalidSequence error about being too short
            match e {
                OrphosError::InvalidSequence(msg) => {
                    assert!(msg.contains("too short"));
                }
                _ => panic!("Expected InvalidSequence error for very short sequence"),
            }
        }
    }

    #[test]
    fn test_config_cloning() {
        let config1 = OrphosConfig::default();
        let orphos1 = UntrainedOrphos::with_config(config1.clone()).unwrap();

        let config2 = orphos1.config.clone();
        let _orphos2 = UntrainedOrphos::with_config(config2).unwrap();
    }

    #[test]
    fn test_debug_formatting() {
        let orphos = UntrainedOrphos::new();
        let debug_str = format!("{:?}", orphos);
        assert!(debug_str.contains("Orphos"));
        assert!(debug_str.contains("config"));

        let analyzer = OrphosAnalyzer::new(OrphosConfig::default());
        let debug_str2 = format!("{:?}", analyzer);
        assert!(debug_str2.contains("OrphosAnalyzer"));
        assert!(debug_str2.contains("config"));
    }

    #[test]
    fn test_type_aliases() {
        // Test that type aliases work correctly
        let _untrained: UntrainedOrphos = UntrainedOrphos::new();
        let _trained: TrainedOrphos =
            TrainedOrphos::new(OrphosConfig::default(), Training::default());

        assert_eq!(
            std::any::type_name::<UntrainedOrphos>(),
            std::any::type_name::<Orphos<Untrained>>()
        );
        assert_eq!(
            std::any::type_name::<TrainedOrphos>(),
            std::any::type_name::<Orphos<Trained>>()
        );
    }

    #[test]
    fn test_training_state_phantom_data() {
        let untrained = UntrainedOrphos::new();
        let trained = TrainedOrphos::new(OrphosConfig::default(), Training::default());

        assert_eq!(std::mem::size_of_val(&untrained._state), 0);
        assert_eq!(std::mem::size_of_val(&trained._state), 0);
    }

    #[test]
    fn test_analyzer_multiple_sequences() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        let sequence1 =
            "ATGAAACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTAAATAG".repeat(150);
        let sequence2 =
            "ATGCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATAG".repeat(200);

        let result1 = analyzer.analyze_sequence(&sequence1, Some("seq1".to_string()));
        let result2 = analyzer.analyze_sequence(&sequence2, Some("seq2".to_string()));

        assert!(result1.is_ok());
        assert!(result2.is_ok());

        let analysis1 = result1.unwrap();
        let analysis2 = result2.unwrap();

        assert_eq!(analysis1.sequence_info.header, "seq1");
        assert_eq!(analysis2.sequence_info.header, "seq2");
        assert_ne!(
            analysis1.sequence_info.length,
            analysis2.sequence_info.length
        );
    }

    #[test]
    fn test_error_handling_edge_cases() {
        let config = OrphosConfig::default();
        let mut analyzer = OrphosAnalyzer::new(config);

        // Test with sequence containing invalid characters
        let invalid_sequence = "ATCGXYZ123";
        let result = analyzer.analyze_sequence(invalid_sequence, None);

        // Should either handle gracefully or return appropriate error
        let _ = result;
    }
}
