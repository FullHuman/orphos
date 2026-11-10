use std::marker::PhantomData;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use crate::algorithms::dynamic_programming::eliminate_bad_genes;
use crate::algorithms::dynamic_programming::predict_genes;
use crate::algorithms::dynamic_programming::{DpWorkspace, predict_genes_in};
use crate::algorithms::gene_finding::GeneBuilder;
use crate::config::OrphosConfig;
use crate::constants::DEFAULT_TRANSLATION_TABLE;
use crate::constants::MIN_SEQUENCE_LENGTH;
use crate::metagenomic::bins;
use crate::metagenomic::{get_preset_training_ref, get_preset_training_shared};
use crate::metagenomic::{preset_length_factors, preset_rbs_lookup};
use crate::node::{
    Given, MotifSite, add_sorted_nodes, calc_orf_gc, calculate_dicodon_gene, motif_sites,
    raw_coding_score, rbs_score, record_gc_bias_of, record_overlapping_starts, reset_node_scores,
    reset_nodes_keeping_coding_scores, score_nodes_given,
};
use crate::parallel;
use crate::results::{OrphosResults, SequenceInfo};
use crate::sequence::GcFramePlot;
use crate::sequence::encoded::EncodedSequence;
use crate::sequence::{FastaRecord, read_fasta_sequences};
use crate::training::non_sd_training::train_starts_nonsd;
use crate::training::sd_training::train_starts_sd;
use crate::training::should_use_sd;
use crate::types::Gene;
use crate::types::{CodonType, Node, OrphosError, Training};
use bio::bio_types::strand::Strand;
use rayon::prelude::*;

static RAYON_POOL_INIT: OnceLock<()> = OnceLock::new();

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
        if config.circular && config.closed_ends {
            return Err(OrphosError::InvalidSequence(
                "Invalid configuration: circular and closed_ends cannot both be true".to_string(),
            ));
        }

        let orphos = Self {
            config,
            training: None,
            _state: PhantomData,
        };

        if let Some(num_threads) = orphos.config.num_threads {
            RAYON_POOL_INIT.get_or_init(|| {
                // Ignore error if pool was already initialized by a previous call
                let _ = rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build_global();
            });
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
        self.train_single_genome_keeping_nodes(encoded_sequence)
            .map(|(trained, _)| trained)
    }

    /// [`Self::train_single_genome`], also returning the nodes it trained
    /// on, sorted by position. Reset with `reset_nodes_to_pristine`, they are
    /// the nodes `add_nodes` builds for the same sequence with the trained
    /// model, so the gene prediction that follows can start from them. They
    /// hold the raw coding scores of the trained model, which the prediction
    /// can keep as well (`reset_nodes_keeping_coding_scores`).
    fn train_single_genome_keeping_nodes(
        &mut self,
        encoded_sequence: &EncodedSequence,
    ) -> Result<(TrainedOrphos, Vec<Node>), OrphosError> {
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
        // Prodigal trains with translation table 11 unless the user requested
        // another one (`-g`); the table decides which codons are starts and
        // stops throughout training and prediction.
        let mut training = Training {
            gc_content: encoded_sequence.gc_content,
            translation_table: self
                .config
                .translation_table
                .map_or(DEFAULT_TRANSLATION_TABLE, i32::from),
            uses_shine_dalgarno: false,
            gc_bias_factors: [0.0; 3],
            ..Default::default()
        };

        // The GC frame plot depends on the sequence alone, so it is computed
        // while the nodes are created and sorted.
        let mut nodes = Vec::new();
        let (created, gc_frame) = parallel::join(
            || {
                add_sorted_nodes(
                    encoded_sequence,
                    &mut nodes,
                    self.config.closed_ends,
                    self.config.circular,
                    &training,
                )
            },
            || GcFramePlot::new(&encoded_sequence.forward_sequence, sequence_length),
        );
        let num_nodes = created?;

        if !self.config.quiet {
            eprintln!(
                "Located {} potential start/stop nodes, closed {}",
                num_nodes, self.config.closed_ends
            );
        }

        record_gc_bias_of(&gc_frame, &mut nodes, &mut training);

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

        Ok((
            Orphos {
                config: self.config.clone(),
                training: Some(training),
                _state: PhantomData,
            },
            nodes,
        ))
    }

    pub fn train_meta_genome(
        &mut self,
        _encoded_sequence: &EncodedSequence,
    ) -> Result<TrainedOrphos, OrphosError> {
        if !self.config.quiet {
            eprintln!("Request: Metagenomic, Phase: Training");
            eprintln!("Initializing training files...");
        }
        // initialize metagenomic bins;
        if !self.config.quiet {
            eprintln!("Metagenomic training initialized.");
        }
        Ok(Orphos {
            config: self.config.clone(),
            training: None,
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

    /// The training data used for gene finding, if any.
    #[must_use]
    pub const fn training(&self) -> Option<&Training> {
        self.training.as_ref()
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
        sequence_number: usize,
    ) -> Result<Vec<Gene>, OrphosError> {
        let mut nodes = Vec::new();
        let training = self
            .training
            .as_ref()
            .ok_or_else(|| OrphosError::InvalidSequence("Orphos is not trained".to_string()))?;

        add_sorted_nodes(
            encoded_sequence,
            &mut nodes,
            self.config.closed_ends,
            self.config.circular,
            training,
        )?;

        self.find_genes_single_in(encoded_sequence, &mut nodes, false, sequence_number)
    }

    /// [`Self::find_genes_single`] on `nodes`, the nodes of `encoded_sequence`
    /// as `add_nodes` creates them with this model, sorted by position.
    /// `has_coding_scores` tells that they hold the raw coding scores of this
    /// model already. The genes are identified as those of sequence number
    /// `sequence_number`.
    fn find_genes_single_in(
        &self,
        encoded_sequence: &EncodedSequence,
        nodes: &mut [Node],
        has_coding_scores: bool,
        sequence_number: usize,
    ) -> Result<Vec<Gene>, OrphosError> {
        let training = self
            .training
            .as_ref()
            .ok_or_else(|| OrphosError::InvalidSequence("Orphos is not trained".to_string()))?;

        // Tap the training state right before scoring to ensure parity at inference time

        score_nodes_given(
            encoded_sequence,
            nodes,
            training,
            self.config.closed_ends,
            false,
            Given {
                coding_scores: has_coding_scores,
                ..Given::default()
            },
        )?;

        record_overlapping_starts(nodes, training, true);

        let gene_path = match predict_genes(nodes, training, true) {
            Some(path) => path,
            None => {
                // No genes found - return empty gene list
                return Ok(vec![]);
            }
        };

        eliminate_bad_genes(nodes, Some(gene_path), training);

        let genes = GeneBuilder::from_nodes(nodes, gene_path, training, sequence_number)
            .with_tweaked_starts()
            .with_annotations()
            .build();

        Ok(genes)
    }

    /// Scores `nodes` with the model of `preset` and predicts its genes.
    ///
    /// Returns the score of the prediction, if there is one, and the genes
    /// unless a prediction with a higher score is known already.
    fn score_preset(
        &self,
        sequence: &PresetSequence,
        nodes: &mut [Node],
        workspace: &mut DpWorkspace,
        preset: usize,
    ) -> Result<Option<PresetCandidate>, OrphosError> {
        let PresetSequence {
            encoded_sequence,
            node_motif_sites,
            best_so_far,
            sequence_number,
        } = *sequence;
        let preset_training = get_preset_training_ref(preset).unwrap();
        reset_node_scores(nodes);
        score_nodes_given(
            encoded_sequence,
            nodes,
            preset_training,
            self.config.closed_ends,
            true,
            Given {
                orf_gc: true,
                rbs_lookup: preset_rbs_lookup(preset),
                length_factors: preset_length_factors(preset),
                motif_sites: node_motif_sites,
                ..Given::default()
            },
        )?;
        record_overlapping_starts(nodes, preset_training, true);
        let gene_path = predict_genes_in(nodes, preset_training, true, workspace);
        let Some((path, score)) =
            gene_path.and_then(|path| Some((path, nodes.get(path)?.scores.total_score)))
        else {
            // Without a path Prodigal reads the score of the node before
            // the first one, which is zero in practice: the preset counts
            // as a prediction of no genes with a score of zero, so that a
            // sequence no preset finds a gene in is reported under the
            // first preset tried, not under preset 0.
            return Ok(Some(PresetCandidate {
                preset,
                score: 0.0,
                genes: Some(Vec::new()),
            }));
        };

        let may_be_best = {
            let mut best_so_far = best_so_far.lock().unwrap();
            if score > *best_so_far {
                *best_so_far = score;
            }
            score >= *best_so_far
        };
        let genes = may_be_best.then(|| {
            eliminate_bad_genes(nodes, gene_path, preset_training);
            GeneBuilder::from_nodes(nodes, path, preset_training, sequence_number)
                .with_tweaked_starts()
                .with_annotations()
                .build()
        });
        Ok(Some(PresetCandidate {
            preset,
            score,
            genes,
        }))
    }

    fn find_genes_meta(
        &self,
        encoded_sequence: &EncodedSequence,
        sequence_number: usize,
    ) -> Result<(Vec<Gene>, usize), OrphosError> {
        if !self.config.quiet {
            eprintln!("Request: Metagenomic, Phase: Gene Finding");
        }
        let mut low = 0.88495 * encoded_sequence.gc_content - 0.0102337;
        if low > 0.65 {
            low = 0.65;
        }
        let mut high = 0.86596 * encoded_sequence.gc_content + 0.1131991;
        if high < 0.35 {
            high = 0.35;
        }
        let mut nodes = Vec::new();
        let mut workspace = DpWorkspace::default();
        // The best score any preset has reached so far. A preset whose
        // score is below it cannot be the best, so its genes are not built.
        let best_so_far = Mutex::new(-100.0_f64);
        let mut candidates: Vec<PresetCandidate> = Vec::new();

        // Consecutive presets with the same translation table score the
        // same nodes.
        let preset_count = bins().len();
        let mut group_start = 0;
        while group_start < preset_count {
            let translation_table = get_preset_training_ref(group_start)
                .unwrap()
                .translation_table;
            let group_end = (group_start + 1..preset_count)
                .find(|&preset| {
                    get_preset_training_ref(preset).unwrap().translation_table != translation_table
                })
                .unwrap_or(preset_count);
            let in_gc_range: Vec<usize> = (group_start..group_end)
                .filter(|&preset| {
                    let gc_content = get_preset_training_ref(preset).unwrap().gc_content;
                    !(gc_content < low || gc_content > high)
                })
                .collect();

            group_start = group_end;
            // Only the presets close to the sequence in GC content are
            // tried; without one, the group's nodes are not needed.
            let Some((&first, rest)) = in_gc_range.split_first() else {
                continue;
            };
            add_sorted_nodes(
                encoded_sequence,
                &mut nodes,
                self.config.closed_ends,
                self.config.circular,
                get_preset_training_ref(first).unwrap(),
            )?;

            // The ORF GC content depends only on the nodes and the
            // sequence, so compute it once per node set rather than once
            // per preset. So does what the upstream motif searches of the
            // presets without Shine-Dalgarno motifs read of the sequence.
            calc_orf_gc(
                &encoded_sequence.forward_sequence,
                encoded_sequence.sequence_length,
                &mut nodes,
            );
            let node_motif_sites = in_gc_range
                .iter()
                .any(|&preset| !get_preset_training_ref(preset).unwrap().uses_shine_dalgarno)
                .then(|| {
                    motif_sites(
                        &encoded_sequence.forward_sequence,
                        &encoded_sequence.reverse_complement_sequence,
                        encoded_sequence.sequence_length,
                        &nodes,
                    )
                });
            let sequence = PresetSequence {
                encoded_sequence,
                node_motif_sites: node_motif_sites.as_deref(),
                best_so_far: &best_so_far,
                sequence_number,
            };

            // Scoring turns the start nodes at the ends of the sequence
            // into edge nodes, for good, and some nodes are scored by what
            // the nodes after them are. So the first preset scores the
            // nodes as they were created, and every other preset the nodes
            // as the first one left them, in whatever order they run.
            candidates.extend(self.score_preset(&sequence, &mut nodes, &mut workspace, first)?);

            let copies = NodeCopies::reserve(rest.len(), &nodes);
            if copies.count < 2 {
                for &preset in rest {
                    candidates.extend(self.score_preset(
                        &sequence,
                        &mut nodes,
                        &mut workspace,
                        preset,
                    )?);
                }
            } else {
                // Each worker scores its share of the presets on its own
                // copy of the nodes.
                let nodes = &nodes;
                let scored: Vec<Vec<PresetCandidate>> = (0..copies.count)
                    .into_par_iter()
                    .map(|worker| {
                        let mut nodes = nodes.clone();
                        let mut workspace = DpWorkspace::default();
                        let mut scored = Vec::new();
                        for &preset in rest.iter().skip(worker).step_by(copies.count) {
                            scored.extend(self.score_preset(
                                &sequence,
                                &mut nodes,
                                &mut workspace,
                                preset,
                            )?);
                        }
                        Ok(scored)
                    })
                    .collect::<Result<_, OrphosError>>()?;
                candidates.extend(scored.into_iter().flatten());
            }
        }

        // The best preset is the first, in preset order, with the highest
        // score.
        candidates.sort_by_key(|candidate| candidate.preset);
        let mut max_score = -100.0;
        let mut best: Option<PresetCandidate> = None;
        for candidate in candidates {
            if candidate.score > max_score {
                max_score = candidate.score;
                best = Some(candidate);
            }
        }
        let (max_phase, mut genes) = match best {
            Some(best) => (
                best.preset,
                best.genes
                    .expect("the genes of a preset are built unless a better one is known"),
            ),
            None => (0, Vec::new()),
        };

        nodes.clear();
        let best_training = get_preset_training_ref(max_phase).unwrap();
        add_sorted_nodes(
            encoded_sequence,
            &mut nodes,
            self.config.closed_ends,
            self.config.circular,
            best_training,
        )?;
        score_nodes_given(
            encoded_sequence,
            &mut nodes,
            best_training,
            self.config.closed_ends,
            true,
            Given {
                rbs_lookup: preset_rbs_lookup(max_phase),
                length_factors: preset_length_factors(max_phase),
                ..Given::default()
            },
        )?;
        // Update display_score from re-scored nodes for GFF column 6
        // This matches Prodigal's behavior where column 6 uses re-scored values
        // but the attribute scores stay from the original scoring during the loop
        update_display_scores(&mut genes, &nodes);
        Ok((genes, max_phase))
    }
}

/// A sequence whose nodes preset models score, and what the models share.
#[derive(Clone, Copy)]
struct PresetSequence<'a> {
    encoded_sequence: &'a EncodedSequence,
    /// What the upstream motif searches read of the sequence, if a preset
    /// without Shine-Dalgarno motifs scores the nodes.
    node_motif_sites: Option<&'a [MotifSite]>,
    /// The best score any preset has reached so far. A preset whose score
    /// is below it cannot be the best, so its genes are not built.
    best_so_far: &'a Mutex<f64>,
    /// Number of the sequence in its input, for the identifiers of its genes.
    sequence_number: usize,
}

/// The prediction of one preset model for a sequence.
struct PresetCandidate {
    preset: usize,
    /// Score of the best path through the nodes.
    score: f64,
    /// The genes on that path, unless a better prediction was known when
    /// this one was made.
    genes: Option<Vec<Gene>>,
}

/// Most bytes of node copies in use at a time, over all the sequences
/// being analyzed, for scoring the presets of a sequence in parallel.
const NODE_COPY_BUDGET: usize = 256 << 20;

/// Fewest nodes worth scoring in parallel: with fewer, the presets are
/// scored faster than the copies are made, and sequences this short come by
/// the thousand and already keep every thread busy.
const MIN_NODES_FOR_PARALLEL_PRESETS: usize = 4096;

/// Bytes of node copies in use.
static NODE_COPY_BYTES: AtomicUsize = AtomicUsize::new(0);

/// The right to hold `count` copies of the nodes of a sequence at once,
/// each scored by one thread; given back when dropped.
///
/// The presets of a sequence are independent of each other once the first
/// has been scored, but each scoring writes to every node, so presets
/// scored at the same time need a copy of the nodes each. The copies of a
/// genome's nodes are tens of megabytes each, hence the budget.
struct NodeCopies {
    count: usize,
    bytes: usize,
}

impl NodeCopies {
    /// Reserves up to one copy of `nodes` per thread and per preset in
    /// `preset_count`, as the budget allows; none if `nodes` are few.
    fn reserve(preset_count: usize, nodes: &[Node]) -> Self {
        let none = Self { count: 0, bytes: 0 };
        let wanted = preset_count.min(rayon::current_num_threads());
        if wanted < 2 || nodes.len() < MIN_NODES_FOR_PARALLEL_PRESETS {
            return none;
        }
        let bytes_per_copy = std::mem::size_of_val(nodes);
        let mut in_use = NODE_COPY_BYTES.load(Ordering::Relaxed);
        loop {
            let count = wanted.min(NODE_COPY_BUDGET.saturating_sub(in_use) / bytes_per_copy);
            if count < 2 {
                return none;
            }
            let bytes = count * bytes_per_copy;
            match NODE_COPY_BYTES.compare_exchange_weak(
                in_use,
                in_use + bytes,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => return Self { count, bytes },
                Err(current) => in_use = current,
            }
        }
    }
}

impl Drop for NodeCopies {
    fn drop(&mut self) {
        NODE_COPY_BYTES.fetch_sub(self.bytes, Ordering::Relaxed);
    }
}

/// Update gene display_score from re-scored nodes for GFF column 6 output.
/// Matches Prodigal's behavior where column 6 uses re-scored values.
fn update_display_scores(genes: &mut [Gene], nodes: &[Node]) {
    for gene in genes.iter_mut() {
        // Find the start node by matching position and strand
        let target_pos = if gene.coordinates.strand == Strand::Forward {
            gene.coordinates.begin.saturating_sub(1)
        } else {
            gene.coordinates.end.saturating_sub(1)
        };

        let lo = nodes.partition_point(|n| n.position.index < target_pos);
        if let Some(start_node) = nodes[lo..]
            .iter()
            .take_while(|n| n.position.index == target_pos)
            .find(|n| {
                n.position.strand == gene.coordinates.strand
                    && n.position.codon_type != CodonType::Stop
            })
        {
            gene.display_score =
                Some(start_node.scores.coding_score + start_node.scores.start_score);
        }
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
/// let analyzer = OrphosAnalyzer::new(OrphosConfig::default());
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
/// let analyzer = OrphosAnalyzer::new(config);
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
    /// let analyzer = OrphosAnalyzer::new(OrphosConfig::default());
    /// let results = analyzer.analyze_fasta_file("genomes.fasta")?;
    ///
    /// for (i, result) in results.iter().enumerate() {
    ///     println!("Sequence {}: {} genes", i + 1, result.genes.len());
    /// }
    /// # Ok::<(), orphos_core::types::OrphosError>(())
    /// ```
    pub fn analyze_fasta_file<P: AsRef<Path>>(
        &self,
        path: P,
    ) -> Result<Vec<OrphosResults>, OrphosError> {
        let sequences = read_fasta_sequences(path.as_ref().to_str().unwrap())?;

        self.analyze_sequences(sequences)
    }

    /// Analyzes the sequences of one input, in order.
    ///
    /// The sequences are numbered from 1 in the order given, as Prodigal
    /// numbers the sequences of a file: the number is the `seqnum` of the
    /// sequence in GFF output and the first part of the identifier of its
    /// genes (`2_1` is the first gene of the second sequence). Each record is
    /// the identifier, the description and the bases of a sequence.
    ///
    /// # Errors
    ///
    /// Returns [`OrphosError`] if any sequence fails analysis.
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use orphos_core::{OrphosAnalyzer, config::OrphosConfig};
    ///
    /// let analyzer = OrphosAnalyzer::new(OrphosConfig {
    ///     metagenomic: true,
    ///     ..OrphosConfig::default()
    /// });
    /// let results = analyzer.analyze_sequences(vec![
    ///     ("contig_1".to_string(), None, b"ATGAAACGCATTAGC...".to_vec()),
    ///     ("contig_2".to_string(), None, b"ATGGCTGAAACCGGT...".to_vec()),
    /// ])?;
    ///
    /// assert_eq!(results[1].sequence_info.number, 2);
    /// # Ok::<(), orphos_core::types::OrphosError>(())
    /// ```
    pub fn analyze_sequences(
        &self,
        sequences: Vec<FastaRecord>,
    ) -> Result<Vec<OrphosResults>, OrphosError> {
        sequences
            .into_par_iter()
            .enumerate()
            .map(|(index, (header, description, seq_bytes))| {
                self.analyze_numbered_sequence(&seq_bytes, header, description, index + 1)
            })
            .collect()
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
    /// let analyzer = OrphosAnalyzer::new(OrphosConfig::default());
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
        &self,
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
    /// let analyzer = OrphosAnalyzer::new(OrphosConfig::default());
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
        &self,
        sequence: &[u8],
        header: String,
        description: Option<String>,
    ) -> Result<OrphosResults, OrphosError> {
        self.analyze_numbered_sequence(sequence, header, description, 1)
    }

    /// [`Self::analyze_sequence_bytes`] for the sequence number
    /// `sequence_number` of its input, counted from 1.
    fn analyze_numbered_sequence(
        &self,
        sequence: &[u8],
        header: String,
        description: Option<String>,
        sequence_number: usize,
    ) -> Result<OrphosResults, OrphosError> {
        let sequence_length = sequence.len();
        let encoded_sequence = self.encode_sequence(sequence);

        let mut untrained_orphos = UntrainedOrphos::with_config(self.config.clone())?;
        // Training on a single genome creates and sorts its nodes and
        // computes their raw coding scores; the gene prediction reuses both
        // once the other scores are reset, instead of creating, sorting and
        // scoring the same nodes again.
        let (mut trained_orphos, mut trained_nodes) = if self.config.metagenomic {
            (untrained_orphos.train_meta_genome(&encoded_sequence)?, None)
        } else {
            let (trained, nodes) =
                untrained_orphos.train_single_genome_keeping_nodes(&encoded_sequence)?;
            (trained, Some(nodes))
        };
        let (genes, metagenomic_model, best_training) = if self.config.metagenomic {
            let (genes, best_phase) =
                trained_orphos.find_genes_meta(&encoded_sequence, sequence_number)?;
            let bin = &bins()[best_phase];
            let training = get_preset_training_ref(best_phase).unwrap();
            let model_desc = format!(
                "{}|{}|{}|{:.1}|{}|{}",
                bin.id,
                bin.name,
                bin.domain,
                bin.gc_percent,
                training.translation_table,
                if training.uses_shine_dalgarno { 1 } else { 0 }
            );
            let shared = get_preset_training_shared(best_phase).unwrap();
            (genes, Some(model_desc), Arc::clone(shared))
        } else {
            let genes = match trained_nodes.as_mut() {
                Some(nodes) => {
                    reset_nodes_keeping_coding_scores(nodes);
                    trained_orphos.find_genes_single_in(
                        &encoded_sequence,
                        nodes,
                        true,
                        sequence_number,
                    )?
                }
                None => trained_orphos.find_genes_single(&encoded_sequence, sequence_number)?,
            };
            // The model trained on this sequence is of no further use here.
            let training = trained_orphos.training.take().unwrap_or_default();
            (genes, None, Arc::new(training))
        };

        let num_genes = genes.len();
        Ok(OrphosResults {
            genes,
            training_used: best_training,
            sequence_info: SequenceInfo {
                number: sequence_number,
                length: sequence_length,
                gc_content: encoded_sequence.gc_content,
                num_genes,
                header,
                description,
            },
            metagenomic_model,
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
    fn test_untrained_orphos_with_conflicting_topology_config() {
        let config = OrphosConfig {
            circular: true,
            closed_ends: true,
            ..OrphosConfig::default()
        };

        let result = UntrainedOrphos::with_config(config);
        assert!(result.is_err());
        if let Err(OrphosError::InvalidSequence(msg)) = result {
            assert!(msg.contains("circular"));
            assert!(msg.contains("closed_ends"));
        } else {
            panic!("Expected InvalidSequence error");
        }
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

        let result = trained.find_genes_single(&test_encoded_sequence, 1);

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

        let result = orphos.find_genes_single(&test_encoded_sequence, 1);

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
        let analyzer = OrphosAnalyzer::new(config);

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
        let analyzer = OrphosAnalyzer::new(config);

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
        let analyzer = OrphosAnalyzer::new(config);

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
        let analyzer = OrphosAnalyzer::new(config);

        let sequence =
            "ATGAAACGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTAAATAG".repeat(300);

        let result = analyzer.analyze_sequence(&sequence, None);
        assert!(result.is_ok());

        let analysis = result.unwrap();
        assert!(analysis.metagenomic_model.is_some());
        // Metagenomic model string format: "id|name|domain|gc_percent|transl_table|uses_sd"
        let model = analysis.metagenomic_model.unwrap();
        assert!(
            model.contains("|"),
            "Expected metagenomic model description with | separator"
        );
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
        let analyzer = OrphosAnalyzer::new(config);

        let result = analyzer.analyze_fasta_file("nonexistent_file.fa");
        assert!(result.is_err());
    }

    #[test]
    fn test_analyze_fasta_file_metagenomic() {
        let config = OrphosConfig {
            metagenomic: true,
            ..OrphosConfig::default()
        };
        let analyzer = OrphosAnalyzer::new(config);

        // Create a temporary FASTA file
        let fasta_content = ">test_seq\nATCG\n";
        let temp_dir = env::temp_dir();
        let temp_file = temp_dir.join("test_metagenomic.fa");
        fs::write(&temp_file, fasta_content).unwrap();

        let result = analyzer.analyze_fasta_file(&temp_file);
        assert!(result.is_ok());

        let results = result.unwrap();
        assert_eq!(results.len(), 1); // One result per sequence in the file
        assert!(results[0].genes.is_empty()); // But no genes found in 4bp sequence

        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_train_single_genome_uses_requested_translation_table() {
        let seq = create_training_sequence();
        let encoded = create_encoded_sequence_for_test(&seq);

        let mut default_orphos = UntrainedOrphos::new();
        let trained = default_orphos.train_single_genome(&encoded).unwrap();
        assert_eq!(trained.training().unwrap().translation_table, 11);

        let config = OrphosConfig {
            translation_table: Some(4),
            quiet: true,
            ..Default::default()
        };
        let mut orphos = UntrainedOrphos::with_config(config).unwrap();
        let trained = orphos.train_single_genome(&encoded).unwrap();
        assert_eq!(trained.training().unwrap().translation_table, 4);
    }

    /// A sequence of open reading frames of varied codons in both
    /// orientations, far enough apart to give training something to choose
    /// between.
    fn create_gene_rich_sequence(seed: u64, length: usize) -> Vec<u8> {
        let mut state = seed;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let complement = |base: u8| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        };
        let mut sequence: Vec<u8> = Vec::new();
        while sequence.len() < length {
            let mut gene = b"AGGAGGTAACAATG".to_vec();
            for _ in 0..100 + next() % 300 {
                let codon = [
                    b"GCT", b"GAA", b"AAA", b"CTG", b"GGC", b"ACC", b"GTT", b"CGT",
                ][(next() % 8) as usize];
                gene.extend_from_slice(codon);
            }
            gene.extend_from_slice(b"TAA");
            if next() % 2 == 0 {
                gene.reverse();
                gene.iter_mut().for_each(|base| *base = complement(*base));
            }
            sequence.extend_from_slice(&gene);
            for _ in 0..20 + next() % 80 {
                sequence.push(b"ACGT"[(next() % 4) as usize]);
            }
        }
        sequence
    }

    /// Gene prediction on the nodes training leaves behind, with their
    /// raw coding scores, must find the genes that prediction from scratch
    /// finds with the trained model.
    #[test]
    fn test_prediction_on_training_nodes_matches_prediction_from_scratch() {
        let sequence = create_gene_rich_sequence(0x6A09_E667_F3BC_C908, 60_000);

        for circular in [false, true] {
            let config = OrphosConfig {
                circular,
                quiet: true,
                ..OrphosConfig::default()
            };
            let encoded = EncodedSequence::without_masking(&sequence);
            let (trained, _) = UntrainedOrphos::with_config(config.clone())
                .unwrap()
                .train_single_genome_keeping_nodes(&encoded)
                .unwrap();
            let from_scratch = trained.find_genes_single(&encoded, 1).unwrap();

            let results = OrphosAnalyzer::new(config)
                .analyze_sequence_bytes(&sequence, "genes".to_string(), None)
                .unwrap();

            assert!(from_scratch.len() > 50, "{} genes", from_scratch.len());
            assert_eq!(
                format!("{from_scratch:?}"),
                format!("{:?}", results.genes),
                "circular {circular}"
            );
        }
    }

    /// The preset models of a long sequence are scored on as many threads
    /// as there are; the prediction must not depend on how many.
    #[test]
    fn test_metagenomic_prediction_does_not_depend_on_the_thread_count() {
        let sequence = create_gene_rich_sequence(0xBB67_AE85_84CA_A73B, 150_000);
        let config = OrphosConfig {
            metagenomic: true,
            quiet: true,
            ..OrphosConfig::default()
        };

        // Long enough for its presets to be scored in parallel.
        let encoded = EncodedSequence::without_masking(&sequence);
        let mut nodes = Vec::new();
        add_sorted_nodes(
            &encoded,
            &mut nodes,
            false,
            false,
            get_preset_training_ref(10).unwrap(),
        )
        .unwrap();
        assert!(
            nodes.len() >= MIN_NODES_FOR_PARALLEL_PRESETS,
            "{} nodes",
            nodes.len()
        );

        let predict_on = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            pool.install(|| {
                OrphosAnalyzer::new(config.clone())
                    .analyze_sequence_bytes(&sequence, "genes".to_string(), None)
                    .unwrap()
            })
        };
        let on_one_thread = predict_on(1);
        assert!(
            on_one_thread.genes.len() > 100,
            "{} genes",
            on_one_thread.genes.len()
        );
        for threads in [2, 3, 8] {
            let on_many_threads = predict_on(threads);
            assert_eq!(
                on_one_thread.metagenomic_model, on_many_threads.metagenomic_model,
                "{threads} threads"
            );
            assert_eq!(
                format!("{:?}", on_one_thread.genes),
                format!("{:?}", on_many_threads.genes),
                "{threads} threads"
            );
        }
    }

    /// The sequences of an input are numbered in order, and the number
    /// shows in the identifiers of their genes and in the GFF header.
    #[test]
    fn test_sequences_are_numbered_in_input_order() {
        use crate::output::write_results;

        let config = OrphosConfig {
            metagenomic: true,
            quiet: true,
            output_format: OutputFormat::Gff,
            ..OrphosConfig::default()
        };
        let analyzer = OrphosAnalyzer::new(config);
        // The third sequence is too short to have a gene; it is counted.
        let records: Vec<FastaRecord> = vec![
            (
                "first".to_string(),
                None,
                create_gene_rich_sequence(1, 4_000),
            ),
            (
                "second".to_string(),
                Some("of four".to_string()),
                create_gene_rich_sequence(2, 9_000),
            ),
            ("third".to_string(), None, b"ACGTACGTAC".to_vec()),
            (
                "fourth".to_string(),
                None,
                create_gene_rich_sequence(4, 5_000),
            ),
        ];

        let results = analyzer.analyze_sequences(records.clone()).unwrap();

        assert_eq!(results.len(), 4);
        for (index, result) in results.iter().enumerate() {
            let number = index + 1;
            assert_eq!(result.sequence_info.header, records[index].0);
            assert_eq!(result.sequence_info.number, number);
            for (gene_index, gene) in result.genes.iter().enumerate() {
                assert_eq!(
                    gene.annotation.identifier,
                    format!("{number}_{}", gene_index + 1)
                );
            }

            let mut gff = Vec::new();
            write_results(&mut gff, result, OutputFormat::Gff).unwrap();
            let gff = String::from_utf8(gff).unwrap();
            assert!(
                gff.contains(&format!("# Sequence Data: seqnum={number};seqlen=")),
                "{gff}"
            );
        }
        assert!(results[2].genes.is_empty());
        for index in [0, 1, 3] {
            assert!(results[index].genes.len() > 3, "sequence {}", index + 1);
        }

        // Sequences predicted with the same preset share its model rather
        // than each holding a copy of it.
        let same_model =
            |a: usize, b: usize| results[a].metagenomic_model == results[b].metagenomic_model;
        let (a, b) = [(0, 1), (0, 3), (1, 3)]
            .into_iter()
            .find(|&(a, b)| same_model(a, b))
            .expect("two of three similar sequences are predicted with the same preset");
        assert!(Arc::ptr_eq(
            &results[a].training_used,
            &results[b].training_used
        ));

        // A sequence analyzed on its own is the first of its input.
        let (header, description, sequence) = records[1].clone();
        let alone = analyzer
            .analyze_sequence_bytes(&sequence, header, description)
            .unwrap();
        assert_eq!(alone.sequence_info.number, 1);
        assert_eq!(alone.genes[0].annotation.identifier, "1_1");
        assert_eq!(alone.genes.len(), results[1].genes.len());
    }

    #[test]
    fn test_analyze_fasta_file_single_genome() {
        let config = OrphosConfig::default();
        let analyzer = OrphosAnalyzer::new(config);

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
        let analyzer = OrphosAnalyzer::new(config);

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
        let analyzer = OrphosAnalyzer::new(config);

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
        let analyzer = OrphosAnalyzer::new(config);

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
        let analyzer = OrphosAnalyzer::new(config);

        // Test with sequence containing invalid characters
        let invalid_sequence = "ATCGXYZ123";
        let result = analyzer.analyze_sequence(invalid_sequence, None);

        // Should either handle gracefully or return appropriate error
        let _ = result;
    }
}
