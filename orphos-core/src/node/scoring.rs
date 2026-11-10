use crate::{
    constants::{
        CODING_SCORE_SENTINEL, CODING_SCORE_THRESHOLD, CODON_LENGTH, DICODON_SIZE, EDGE_BONUS,
        EDGE_POSITION_OFFSET, EDGE_POSITION_THRESHOLD, EDGE_UPSTREAM, GENE_SIZE_SCALING_FACTOR,
        LENGTH_FACTOR_MULTIPLIER, LENGTH_FACTOR_THRESHOLD, MAX_GENE_SIZE_CODONS,
        METAGENOMIC_LENGTH_THRESHOLD, METAGENOMIC_MIN_LENGTH_FACTOR, METAGENOMIC_PENALTY,
        METAGENOMIC_PENALTY_DIVISOR, MIN_GENE_SIZE_CODONS, MIN_META_GENE_LENGTH,
        NEGATIVE_SCORE_PENALTY, NO_MOTIF_THRESHOLD, NODE_SEARCH_WINDOW, SHORT_GENE_THRESHOLD,
        UPSTREAM_COMPOSITION_WEIGHT, UPSTREAM_SCAN_RANGE, UPSTREAM_SKIP_END, UPSTREAM_SKIP_START,
    },
    node::motifs::{MotifSite, RbsLookup, rbs_score_with},
    node::{calc_orf_gc, find_best_upstream_motif},
    sequence::encoded::EncodedSequence,
    sequence::{calculate_kmer_index, codon_index, is_stop},
    types::{CodonType, Node, OrphosError, Training},
};

use std::sync::atomic::{AtomicU64, Ordering};

use bio::bio_types::strand::Strand;

use crate::node::strands::{StrandNodes, has_start_before_stop, split_by_strand};
use crate::parallel;

/// Context struct containing sequence data and metadata
#[derive(Debug)]
pub struct ScoringContext<'a> {
    /// Encoded sequence containing forward and reverse strands
    pub encoded_sequence: &'a EncodedSequence,
    /// Training data for scoring
    pub training: &'a Training,
    /// Whether the sequence is in closed reading frame mode
    pub closed: bool,
    /// Whether this is metagenomic mode
    pub is_meta: bool,
    /// Cached upstream composition scan indices (for performance)
    pub upstream_scan_indices: Vec<usize>,
}

impl<'a> ScoringContext<'a> {
    /// Create a new scoring context
    pub fn new(
        encoded_sequence: &'a EncodedSequence,
        training: &'a Training,
        closed: bool,
        is_meta: bool,
    ) -> Self {
        // Pre-calculate upstream scan indices to avoid repeated computation
        let upstream_scan_indices: Vec<usize> = (1..UPSTREAM_SCAN_RANGE)
            .filter(|&i| i <= UPSTREAM_SKIP_START || i >= UPSTREAM_SKIP_END)
            .collect();

        Self {
            encoded_sequence,
            training,
            closed,
            is_meta,
            upstream_scan_indices,
        }
    }
}

/// Context struct containing scoring parameters and factors
#[derive(Debug, Clone)]
pub struct ScoringParams {
    /// Edge bonus factor
    pub edge_bonus_factor: f64,
    /// Cached start weight factor
    pub start_weight_factor: f64,
    /// Cached no-stop probability for length calculations
    pub no_stop_probability: f64,
}

impl ScoringParams {
    /// Create scoring parameters from training data
    pub fn from_training(training: &Training) -> Self {
        let start_weight_factor = training.start_weight_factor;
        let edge_bonus_factor = EDGE_BONUS * start_weight_factor;
        let no_stop_probability =
            calculate_no_stop_probability(training.gc_content, training.translation_table);

        Self {
            edge_bonus_factor,
            start_weight_factor,
            no_stop_probability,
        }
    }
}

/// Cached calculations for a node to avoid repetitive computations
#[derive(Debug, Clone)]
pub struct NodeCache {
    /// Gene length in base pairs
    pub gene_length: usize,
    /// Whether the node is at sequence boundaries
    pub is_boundary_node: bool,
    /// Calculated short gene factors (negative, positive)
    pub short_gene_factors: Option<(f64, f64)>,
    /// Whether metagenomic penalties apply
    pub needs_metagenomic_penalty: bool,
}

impl NodeCache {
    /// Create a new cache for a node
    pub fn new(node: &Node, context: &ScoringContext) -> Self {
        let gene_length = calculate_gene_length(node);

        let is_boundary_node = is_node_at_boundary(node, context);

        let short_gene_factors = if gene_length < SHORT_GENE_THRESHOLD {
            Some(calculate_short_gene_factors(gene_length))
        } else {
            None
        };

        let needs_metagenomic_penalty = context.is_meta
            && context.encoded_sequence.sequence_length < METAGENOMIC_LENGTH_THRESHOLD
            && (node.scores.coding_score < CODING_SCORE_THRESHOLD
                || gene_length < MIN_META_GENE_LENGTH);

        Self {
            gene_length,
            is_boundary_node,
            short_gene_factors,
            needs_metagenomic_penalty,
        }
    }
}

/// Check if a node is at sequence boundaries
fn is_node_at_boundary(node: &Node, context: &ScoringContext) -> bool {
    (node.position.index <= EDGE_POSITION_THRESHOLD && node.position.strand == Strand::Forward)
        || (node.position.index >= context.encoded_sequence.sequence_length - EDGE_POSITION_OFFSET
            && node.position.strand == Strand::Reverse)
}

/// Calculate edge gene status for a node
///
/// Returns the number of edge gene indicators:
/// - 1 if the node is already marked as edge
/// - +1 if the stop codon is invalid (suggesting edge condition)
fn calculate_edge_gene_status(node: &Node, context: &ScoringContext) -> usize {
    let mut edge_gene = 0;

    if node.position.is_edge {
        edge_gene += 1;
    }

    // Check if stop codon is valid
    if (node.position.strand == Strand::Forward
        && (node.position.stop_value < 0
            || !is_stop(
                &context.encoded_sequence.forward_sequence,
                node.position.stop_value as usize,
                context.training,
            )))
        || (node.position.strand == Strand::Reverse
            && (node.position.stop_value < 0
                || !is_stop(
                    &context.encoded_sequence.reverse_complement_sequence,
                    context.encoded_sequence.sequence_length
                        - 1
                        - node.position.stop_value as usize,
                    context.training,
                )))
    {
        edge_gene += 1;
    }

    edge_gene
}

/// Calculate type scores for a node
///
/// For edge genes: Uses edge bonus factor divided by edge gene count
/// For regular genes: Uses training weights based on codon type
fn calculate_type_scores(
    node: &mut Node,
    edge_gene: usize,
    params: &ScoringParams,
    training: &Training,
) {
    if node.position.is_edge {
        node.scores.type_score = params.edge_bonus_factor / edge_gene as f64;
        node.scores.upstream_score = 0.0;
        node.scores.ribosome_binding_score = 0.0;
    } else {
        // Type score using cached weight factor
        node.scores.type_score = training.start_type_weights[node.position.codon_type.to_index()]
            * params.start_weight_factor;
    }
}

/// Calculate RBS (ribosome binding site) scores for a node
///
/// Uses either Shine-Dalgarno detection or motif scoring based on training parameters
fn calculate_rbs_scores(node: &mut Node, training: &Training, params: &ScoringParams) {
    if node.position.is_edge {
        return; // Already handled in calculate_type_scores
    }

    // RBS score using cached weight factor
    let rbs1 = training.rbs_weights[node.motif_info.ribosome_binding_sites[0]];
    let rbs2 = training.rbs_weights[node.motif_info.ribosome_binding_sites[1]];
    let sd_score = rbs1.max(rbs2) * params.start_weight_factor;

    if training.uses_shine_dalgarno {
        node.scores.ribosome_binding_score = sd_score;
    } else {
        node.scores.ribosome_binding_score =
            params.start_weight_factor * node.motif_info.best_motif.score;
        if node.scores.ribosome_binding_score < sd_score
            && training.no_motif_weight > NO_MOTIF_THRESHOLD
        {
            node.scores.ribosome_binding_score = sd_score;
        }
    }
}

/// Calculate upstream scores for a node
///
/// Analyzes nucleotide composition upstream of the start codon and applies
/// edge bonuses for nodes near sequence boundaries
fn calculate_upstream_scores(
    node: &mut Node,
    context: &ScoringContext,
    params: &ScoringParams,
    cache: &NodeCache,
) {
    if node.position.is_edge {
        return; // Already handled in calculate_type_scores
    }

    if node.position.strand == Strand::Forward {
        node.scores.upstream_score =
            score_upstream_composition(&context.encoded_sequence.forward_sequence, node, context);
    } else {
        node.scores.upstream_score = score_upstream_composition(
            &context.encoded_sequence.reverse_complement_sequence,
            node,
            context,
        );
    }

    // Apply edge upstream bonus using cached boundary check
    if !context.closed && cache.is_boundary_node {
        node.scores.upstream_score =
            EDGE_UPSTREAM.mul_add(params.start_weight_factor, node.scores.upstream_score);
    }
}

/// The stop values of the edge nodes near the two ends of the node array,
/// where a start node looks for an edge node of its own open reading frame.
///
/// A forward start among the first [`NODE_SEARCH_WINDOW`] nodes looks at the
/// nodes before it, a reverse start among the last ones at the nodes after
/// it. There are only a handful of edge nodes, so their stop values are
/// kept aside instead of every such start walking over its neighbours.
struct EdgeStops {
    /// Stop values of the edge nodes among the nodes scored so far, while
    /// those are among the first [`NODE_SEARCH_WINDOW`] nodes.
    scored: Vec<isize>,
    /// Index and stop value of the edge nodes among the last
    /// [`NODE_SEARCH_WINDOW`] nodes, as they are before any is scored.
    last: Vec<(usize, isize)>,
}

impl EdgeStops {
    fn new(nodes: &[Node]) -> Self {
        let first_of_last = nodes.len().saturating_sub(NODE_SEARCH_WINDOW);
        Self {
            scored: Vec::new(),
            last: nodes[first_of_last..]
                .iter()
                .enumerate()
                .filter(|(_, node)| node.position.is_edge)
                .map(|(offset, node)| (first_of_last + offset, node.position.stop_value))
                .collect(),
        }
    }

    /// Records node `index` once it is scored, or passed over.
    #[inline]
    fn scored(&mut self, index: usize, node: &Node) {
        if index < NODE_SEARCH_WINDOW && node.position.is_edge {
            self.scored.push(node.position.stop_value);
        }
    }
}

/// Apply edge upstream bonus based on neighboring nodes
///
/// Searches for edge nodes with the same stop position within a limited window
/// and applies bonus to the current node if found
fn apply_edge_upstream_bonus_for_neighbors(
    node: &mut Node,
    current_index: usize,
    node_count: usize,
    params: &ScoringParams,
    boundary_bonus_applied: bool,
    edge_stops: &EdgeStops,
) {
    // Prodigal applies the boundary bonus and the neighbour search as one
    // if/else-if chain, so a node that received the boundary bonus is never
    // searched.
    if node.position.is_edge || boundary_bonus_applied {
        return;
    }

    let stop_value = node.position.stop_value;
    let has_edge_neighbour =
        if current_index < NODE_SEARCH_WINDOW && node.position.strand == Strand::Forward {
            // Forward strand: an edge node among the previous nodes
            edge_stops.scored.contains(&stop_value)
        } else if current_index >= node_count.saturating_sub(NODE_SEARCH_WINDOW)
            && node.position.strand == Strand::Reverse
        {
            // Reverse strand: an edge node among the following nodes
            edge_stops
                .last
                .iter()
                .any(|&(index, edge_stop)| index > current_index && edge_stop == stop_value)
        } else {
            false
        };
    if has_edge_neighbour {
        node.scores.upstream_score =
            EDGE_UPSTREAM.mul_add(params.start_weight_factor, node.scores.upstream_score);
    }
}

/// Convert nodes to edge genes if they are at sequence boundaries
///
/// In open reading frame mode (!closed), nodes near sequence start/end
/// are converted to edge genes with special scoring
fn convert_boundary_nodes_to_edge_genes(
    node: &mut Node,
    context: &ScoringContext,
    edge_gene: &mut usize,
    params: &ScoringParams,
    cache: &NodeCache,
) {
    // Convert starts at base 1 and sequence_length to edge genes if closed = 0
    if cache.is_boundary_node && !node.position.is_edge && !context.closed {
        *edge_gene += 1;
        node.position.is_edge = true;
        node.scores.type_score = 0.0;
        let new_upstream_score = params.edge_bonus_factor / *edge_gene as f64;
        node.scores.upstream_score = new_upstream_score;
        node.scores.ribosome_binding_score = 0.0;
    }
}

/// Apply penalties and bonuses to node scores
///
/// Handles various scoring adjustments including:
/// - Penalties for genes without proper stop codons
/// - Short gene penalties/bonuses
/// - Metagenomic-specific penalties
/// - Negative coding score penalties
fn apply_penalties_and_bonuses(
    node: &mut Node,
    edge_gene: usize,
    context: &ScoringContext,
    params: &ScoringParams,
    cache: &NodeCache,
) {
    // Apply penalty for starts with no stop codon - this happens AFTER boundary conversion
    // (`uscore -= 0.5*EDGE_BONUS*st_wt`, one fused multiply-add in the reference build)
    if !node.position.is_edge && edge_gene == 1 {
        node.scores.upstream_score = (-(NEGATIVE_SCORE_PENALTY * EDGE_BONUS))
            .mul_add(params.start_weight_factor, node.scores.upstream_score);
    }

    // Apply short gene penalties/bonuses using cached factors
    if edge_gene == 0
        && let Some((negative_factor, positive_factor)) = cache.short_gene_factors
    {
        if node.scores.ribosome_binding_score < 0.0 {
            node.scores.ribosome_binding_score *= negative_factor;
        }
        if node.scores.upstream_score < 0.0 {
            node.scores.upstream_score *= negative_factor;
        }
        if node.scores.type_score < 0.0 {
            node.scores.type_score *= negative_factor;
        }

        if node.scores.ribosome_binding_score > 0.0 {
            node.scores.ribosome_binding_score *= positive_factor;
        }
        if node.scores.upstream_score > 0.0 {
            node.scores.upstream_score *= positive_factor;
        }
        if node.scores.type_score > 0.0 {
            node.scores.type_score *= positive_factor;
        }
    }

    // Apply metagenomic penalties using cached calculation
    // (`cscore -= META_PEN*dmax(0, (3000-slen)/2700.0)`)
    if cache.needs_metagenomic_penalty && edge_gene == 0 {
        let shortfall = (METAGENOMIC_LENGTH_THRESHOLD as f64
            - context.encoded_sequence.sequence_length as f64)
            / METAGENOMIC_PENALTY_DIVISOR;
        node.scores.coding_score =
            (-METAGENOMIC_PENALTY).mul_add(shortfall.max(0.0), node.scores.coding_score);
    }

    // Base start score calculation - this happens in C code at line 498
    node.scores.start_score =
        node.scores.type_score + node.scores.ribosome_binding_score + node.scores.upstream_score;

    // Penalize negative coding scores
    if node.scores.coding_score < 0.0 {
        if edge_gene > 0 && !node.position.is_edge {
            if !context.is_meta
                || context.encoded_sequence.sequence_length > METAGENOMIC_MIN_LENGTH_FACTOR
            {
                node.scores.start_score -= params.start_weight_factor;
            } else {
                // C: nod[i].sscore -= (10.31 - 0.004*slen);
                let penalty =
                    (-0.004_f64).mul_add(context.encoded_sequence.sequence_length as f64, 10.31);
                node.scores.start_score -= penalty;
            }
        } else if context.is_meta
            && context.encoded_sequence.sequence_length < METAGENOMIC_LENGTH_THRESHOLD
            && node.position.is_edge
        {
            let min_meta_length =
                calculate_min_metagenomic_length(context.encoded_sequence.sequence_length);
            if cache.gene_length as f64 >= min_meta_length {
                if node.scores.coding_score >= 0.0 {
                    node.scores.coding_score = -1.0;
                }
                node.scores.start_score = 0.0;
                node.scores.upstream_score = 0.0;
            }
        } else {
            node.scores.start_score -= NEGATIVE_SCORE_PENALTY;
        }
    } else if node.scores.coding_score < CODING_SCORE_THRESHOLD
        && context.is_meta
        && cache.gene_length < MIN_META_GENE_LENGTH
        && node.scores.start_score < 0.0
    {
        node.scores.start_score -= params.start_weight_factor;
    }
}

/// Calculate final combined scores for a node
///
/// Computes the total score from coding + start scores
/// Start score is already calculated in apply_penalties_and_bonuses
fn calculate_final_scores(node: &mut Node) {
    node.scores.total_score = node.scores.coding_score + node.scores.start_score;
}

/// What the caller of [`score_nodes_given`] has computed already.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Given<'a> {
    /// The nodes hold their ORF GC content ([`calc_orf_gc`]).
    pub orf_gc: bool,
    /// The nodes hold the raw coding scores ([`raw_coding_score`]) of the
    /// training. Training on a single genome computes them to train the
    /// start model, and they depend on nothing start training changes.
    pub coding_scores: bool,
    /// The Shine-Dalgarno motif searches for the weights of the training,
    /// tabulated. The table depends on the weights alone, so it is built
    /// once for a model that scores many sequences.
    pub rbs_lookup: Option<&'a RbsLookup>,
    /// The length factors of the training ([`LengthFactors::for_training`]),
    /// likewise shared by the sequences a model scores.
    pub length_factors: Option<&'a LengthFactors>,
    /// What the upstream motif searches read of the sequence, for every
    /// node ([`motif_sites`]). It depends on the nodes and the sequence
    /// alone, so it is read once for nodes that several models score.
    ///
    /// [`motif_sites`]: crate::node::motifs::motif_sites
    pub motif_sites: Option<&'a [MotifSite]>,
}

/// Score nodes for gene finding
pub fn score_nodes(
    encoded_sequence: &EncodedSequence,
    nodes: &mut [Node],
    training: &Training,
    closed: bool,
    is_meta: bool,
) -> Result<(), OrphosError> {
    score_nodes_given(
        encoded_sequence,
        nodes,
        training,
        closed,
        is_meta,
        Given::default(),
    )
}

/// [`score_nodes`] for nodes whose ORF GC content ([`calc_orf_gc`]) has
/// already been computed and is still valid.
///
/// The GC content only depends on the node positions and the sequence, so
/// when the same nodes are scored with several training sets (metagenomic
/// mode) it needs to be computed once.
pub fn score_nodes_with_orf_gc(
    encoded_sequence: &EncodedSequence,
    nodes: &mut [Node],
    training: &Training,
    closed: bool,
    is_meta: bool,
) -> Result<(), OrphosError> {
    score_nodes_given(
        encoded_sequence,
        nodes,
        training,
        closed,
        is_meta,
        Given {
            orf_gc: true,
            ..Given::default()
        },
    )
}

/// [`score_nodes`], skipping what `given` says the caller has computed.
pub(crate) fn score_nodes_given(
    encoded_sequence: &EncodedSequence,
    nodes: &mut [Node],
    training: &Training,
    closed: bool,
    is_meta: bool,
    given: Given,
) -> Result<(), OrphosError> {
    let node_count = nodes.len();

    // Early exit for empty nodes
    if node_count == 0 {
        return Ok(());
    }

    // Create context structs to reduce parameter passing
    let context = ScoringContext::new(encoded_sequence, training, closed, is_meta);
    let params = ScoringParams::from_training(training);

    if !given.orf_gc {
        calc_orf_gc(
            &encoded_sequence.forward_sequence,
            encoded_sequence.sequence_length,
            nodes,
        );
    }
    if !given.coding_scores {
        let own_length_factors;
        let length_factors = match given.length_factors {
            Some(length_factors) => length_factors,
            None => {
                own_length_factors = LengthFactors::new(params.no_stop_probability);
                &own_length_factors
            }
        };
        score_strands(
            &encoded_sequence.forward_sequence,
            &encoded_sequence.reverse_complement_sequence,
            encoded_sequence.sequence_length,
            nodes,
            training,
            length_factors,
        );
    }

    // Calculate RBS scores using optimized approach
    if training.uses_shine_dalgarno {
        let tabulated;
        let rbs_lookup = match given.rbs_lookup {
            Some(rbs_lookup) => rbs_lookup,
            None => {
                tabulated = RbsLookup::new(&*training.rbs_weights);
                &tabulated
            }
        };
        rbs_score_with(
            &encoded_sequence.forward_sequence,
            &encoded_sequence.reverse_complement_sequence,
            encoded_sequence.sequence_length,
            nodes,
            rbs_lookup,
        );
    } else if let Some(motif_sites) = given.motif_sites {
        for (node, site) in nodes.iter_mut().zip(motif_sites) {
            if node.position.codon_type != CodonType::Stop && !node.position.is_edge {
                node.motif_info.best_motif = site.best_motif(training, 2);
            }
        }
    } else {
        // For non-SD mode, process motif finding with bounds checking
        for node in nodes.iter_mut() {
            if node.position.codon_type != CodonType::Stop && !node.position.is_edge {
                find_best_upstream_motif(
                    training,
                    &encoded_sequence.forward_sequence,
                    &encoded_sequence.reverse_complement_sequence,
                    encoded_sequence.sequence_length,
                    node,
                    2,
                );
            }
        }
    }

    // Score each start node using extracted helper functions with cached calculations
    let mut edge_stops = EdgeStops::new(nodes);
    for (i, node) in nodes.iter_mut().enumerate() {
        if node.position.codon_type == CodonType::Stop {
            edge_stops.scored(i, node);
            continue;
        }

        // Create cache for repetitive calculations
        let cache = NodeCache::new(node, &context);

        // Calculate edge gene status
        let mut edge_gene = calculate_edge_gene_status(node, &context);

        // Calculate various score components
        calculate_type_scores(node, edge_gene, &params, training);

        calculate_rbs_scores(node, training, &params);

        calculate_upstream_scores(node, &context, &params, &cache);

        // Apply edge upstream bonus based on neighboring nodes
        let boundary_bonus_applied = !context.closed && cache.is_boundary_node;
        apply_edge_upstream_bonus_for_neighbors(
            node,
            i,
            node_count,
            &params,
            boundary_bonus_applied,
            &edge_stops,
        );

        // Convert boundary nodes to edge genes if necessary
        convert_boundary_nodes_to_edge_genes(node, &context, &mut edge_gene, &params, &cache);

        // Calculate base start score before applying penalties
        calculate_final_scores(node);

        // Apply penalties and bonuses (includes recalculating start score)
        apply_penalties_and_bonuses(node, edge_gene, &context, &params, &cache);

        // Calculate final total score
        node.scores.total_score = node.scores.coding_score + node.scores.start_score;
        edge_stops.scored(i, node);
    }

    Ok(())
}

/// Score the upstream composition for a given node with optimized scanning
fn score_upstream_composition(seq: &[u8], node: &Node, context: &ScoringContext) -> f64 {
    let start = if node.position.strand == Strand::Forward {
        node.position.index
    } else {
        context.encoded_sequence.sequence_length - 1 - node.position.index
    };

    let mut upstream_score = 0.0;
    let mut count = 0;

    // Use pre-calculated scan indices to avoid repeated filtering
    for &i in &context.upstream_scan_indices {
        if start < i {
            continue;
        }

        let pos = start - i;

        // Get the nucleotide index at position
        let mer_idx = calculate_kmer_index(1, seq, pos);

        // Add to upstream score (`uscore += 0.4*st_wt*ups_comp[..]`, fused)
        upstream_score = (UPSTREAM_COMPOSITION_WEIGHT * context.training.start_weight_factor)
            .mul_add(
                context.training.upstream_composition[count][mer_idx],
                upstream_score,
            );
        count += 1;
    }
    upstream_score
}

/// Calculate the length of a gene in base pairs
///
/// Returns the absolute difference between start and stop positions
const fn calculate_gene_length(node: &Node) -> usize {
    if node.position.stop_value < 0 {
        (node.position.index as isize - node.position.stop_value) as usize
    } else {
        (node.position.stop_value as usize).abs_diff(node.position.index)
    }
}

/// Calculate the length of a gene in codons
///
/// Converts gene length in base pairs to codons, adding 3 to account for partial codons
fn calculate_gene_size_codons(node: &Node) -> f64 {
    ((calculate_gene_length(node) + 3) as f64) / 3.0
}

/// Calculate short gene penalty factors
///
/// Returns (negative_factor, positive_factor) for penalizing/rewarding short genes
fn calculate_short_gene_factors(gene_length: usize) -> (f64, f64) {
    let negative_factor = SHORT_GENE_THRESHOLD as f64 / gene_length as f64;
    let positive_factor = gene_length as f64 / SHORT_GENE_THRESHOLD as f64;
    (negative_factor, positive_factor)
}

/// Calculate minimum metagenomic gene length threshold
///
/// Uses sequence length to determine minimum acceptable gene length for metagenomic mode
fn calculate_min_metagenomic_length(sequence_length: usize) -> f64 {
    (sequence_length as f64).sqrt() * 5.0
}

/// Calculate the probability that a random codon is not a stop codon
///
/// This depends on the genetic code and GC content of the organism.
/// For the standard genetic code (11), there are 3 stop codons out of 64.
/// For other genetic codes, the number varies.
fn calculate_no_stop_probability(gc_content: f64, translation_table: i32) -> f64 {
    let at_content = 1.0 - gc_content;
    let at_squared = at_content * at_content;

    if translation_table != 11 {
        // For non-standard genetic codes, calculate stop codon probability

        // TAA and TAG (AT-rich stop codons)
        let at_stop_prob = (at_squared * gc_content) / 8.0;
        // TGA (mixed stop codon)
        let mixed_stop_prob = (at_squared * at_content) / 8.0;

        1.0 - (at_stop_prob + mixed_stop_prob)
    } else {
        // Standard genetic code has 3 stop codons: TAA, TAG, TGA

        // TAA and TAG probability
        let tag_taa_prob = (at_squared * gc_content) / 4.0;
        // TGA probability
        let tga_prob = (at_squared * at_content) / 8.0;

        1.0 - (tag_taa_prob + tga_prob)
    }
}

/// The three passes of the raw coding score over the nodes of one strand:
/// cumulative dicodon scores from stop to start, normalisation to the
/// frame's maximum, length factors. The first pass starts from, and leaves
/// behind, the position of the last node it saw in each frame.
#[allow(clippy::too_many_arguments)]
fn strand_coding_scores(
    nodes: &mut StrandNodes,
    encoded_sequence: &[u8],
    sequence_length: usize,
    strand: Strand,
    training: &Training,
    length_factors: &LengthFactors,
    last_positions: &mut [usize; 3],
) {
    let mut frame_scores = [0.0; 3];
    process_strand_first_pass(
        nodes,
        encoded_sequence,
        sequence_length,
        strand,
        training,
        &mut frame_scores,
        last_positions,
    );
    process_strand_second_pass(nodes, strand, &mut frame_scores);
    process_strand_third_pass(nodes, strand, length_factors, &mut frame_scores);
}

/// Score the two strands with [`strand_coding_scores`], each on its own
/// thread: the passes of a strand touch that strand's nodes only.
///
/// Prodigal walks the reverse strand from the per-frame positions the
/// forward strand's walk left behind. Only a reverse start node that no
/// reverse stop of its frame precedes sees them (a gene wrapping around a
/// circular sequence); when there is one, the strands are scored one after
/// the other, as before.
fn score_strands(
    encoded_sequence: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &mut [Node],
    training: &Training,
    length_factors: &LengthFactors,
) {
    let (mut forward, mut reverse) = split_by_strand(nodes);
    let mut last_positions = [0usize; 3];
    if has_start_before_stop(&reverse) {
        strand_coding_scores(
            &mut forward,
            encoded_sequence,
            sequence_length,
            Strand::Forward,
            training,
            length_factors,
            &mut last_positions,
        );
        strand_coding_scores(
            &mut reverse,
            reverse_complement_encoded_sequence,
            sequence_length,
            Strand::Reverse,
            training,
            length_factors,
            &mut last_positions,
        );
    } else {
        parallel::join(
            || {
                strand_coding_scores(
                    &mut forward,
                    encoded_sequence,
                    sequence_length,
                    Strand::Forward,
                    training,
                    length_factors,
                    &mut last_positions,
                )
            },
            || {
                strand_coding_scores(
                    &mut reverse,
                    reverse_complement_encoded_sequence,
                    sequence_length,
                    Strand::Reverse,
                    training,
                    length_factors,
                    &mut [0usize; 3],
                )
            },
        );
    }
}

/// Process nodes for a specific strand in the first pass (dicodon scoring)
///
/// This calculates cumulative dicodon scores from stop codons to start codons
fn process_strand_first_pass(
    nodes: &mut StrandNodes,
    encoded_sequence: &[u8],
    sequence_length: usize,
    strand: Strand,
    training: &Training,
    frame_scores: &mut [f64; 3],
    last_positions: &mut [usize; 3],
) {
    frame_scores.fill(0.0);
    let dicodon_table: &[f64; 1 << (2 * DICODON_SIZE)] = &training.gene_dicodon_table;

    // Walking an ORF codon by codon, every dicodon shares its second codon
    // with the dicodon read three bases downstream, so only the new codon is
    // read from the sequence. The running score lives in a local while an
    // ORF is walked; the additions happen in the same order as ever.
    if strand == Strand::Forward {
        // Process forward strand in reverse order (from end to start)
        for node in nodes.iter_mut().rev() {
            let index = node.position.index;
            let frame_index = index % 3;

            if node.position.codon_type == CodonType::Stop {
                last_positions[frame_index] = index;
                frame_scores[frame_index] = 0.0;
            } else {
                // Calculate dicodon scores from last stop to current start
                let last = last_positions[frame_index];
                if last >= index + 3 {
                    let mut j = last - 3;
                    let mut mer_idx = calculate_kmer_index(DICODON_SIZE, encoded_sequence, j);
                    let mut score = frame_scores[frame_index];
                    loop {
                        score += dicodon_table[mer_idx & DICODON_INDEX_MASK];
                        if j < index + 3 {
                            break;
                        }
                        j -= 3;
                        mer_idx = codon_index(encoded_sequence, j)
                            | ((mer_idx & CODON_INDEX_MASK) << (2 * CODON_LENGTH));
                    }
                    frame_scores[frame_index] = score;
                }
                node.scores.coding_score = frame_scores[frame_index];
                last_positions[frame_index] = index;
            }
        }
    } else {
        // Process reverse strand in forward order
        for node in nodes.iter_mut() {
            let index = node.position.index;
            let frame_index = index % 3;

            if node.position.codon_type == CodonType::Stop {
                last_positions[frame_index] = index;
                frame_scores[frame_index] = 0.0;
            } else {
                // Calculate dicodon scores from last stop to current start
                let mut j = last_positions[frame_index] + 3;
                if j <= index {
                    let mut reverse_seq_pos = sequence_length - j - 1;
                    let mut mer_idx =
                        calculate_kmer_index(DICODON_SIZE, encoded_sequence, reverse_seq_pos);
                    let mut score = frame_scores[frame_index];
                    loop {
                        score += dicodon_table[mer_idx & DICODON_INDEX_MASK];
                        j += 3;
                        if j > index {
                            break;
                        }
                        reverse_seq_pos -= 3;
                        mer_idx = codon_index(encoded_sequence, reverse_seq_pos)
                            | ((mer_idx & CODON_INDEX_MASK) << (2 * CODON_LENGTH));
                    }
                    frame_scores[frame_index] = score;
                }
                node.scores.coding_score = frame_scores[frame_index];
                last_positions[frame_index] = index;
            }
        }
    }
}

const CODON_INDEX_MASK: usize = (1 << (2 * CODON_LENGTH)) - 1;
const DICODON_INDEX_MASK: usize = (1 << (2 * DICODON_SIZE)) - 1;

/// Process nodes for a specific strand in the second pass (score normalization)
///
/// This normalizes scores relative to the maximum score in each frame
fn process_strand_second_pass(
    nodes: &mut StrandNodes,
    strand: Strand,
    frame_scores: &mut [f64; 3],
) {
    frame_scores.fill(CODING_SCORE_SENTINEL);

    let process_node = |node: &mut Node| {
        let frame_index = node.position.index % 3;

        if node.position.codon_type == CodonType::Stop {
            frame_scores[frame_index] = CODING_SCORE_SENTINEL;
        } else if node.scores.coding_score > frame_scores[frame_index] {
            frame_scores[frame_index] = node.scores.coding_score;
        } else {
            node.scores.coding_score -= frame_scores[frame_index] - node.scores.coding_score;
        }
    };

    if strand == Strand::Forward {
        // Process forward strand in forward order
        nodes.iter_mut().for_each(process_node);
    } else {
        // Process reverse strand in reverse order
        nodes.iter_mut().rev().for_each(process_node);
    }
}

/// Calculate the length-based scoring factor for a gene
///
/// This implements the biological principle that very short or very long genes
/// are less likely to be real genes
///
/// Powers use `powf` (libm `pow`) rather than `powi`: Prodigal calls `pow`
/// for the integer exponents too, and the two round differently.
fn calculate_length_factor(gene_size_codons: f64, no_stop_probability: f64) -> f64 {
    let max_codons = MAX_GENE_SIZE_CODONS as f64;
    let min_codons = MIN_GENE_SIZE_CODONS as f64;
    if gene_size_codons > max_codons {
        // For very long genes, use a scaled length factor
        let mut factor = ((1.0 - no_stop_probability.powf(max_codons))
            / no_stop_probability.powf(max_codons))
        .ln();
        factor -= ((1.0 - no_stop_probability.powf(min_codons))
            / no_stop_probability.powf(min_codons))
        .ln();
        factor *= (gene_size_codons - min_codons) / GENE_SIZE_SCALING_FACTOR;
        factor
    } else {
        // For normal-sized genes, use direct calculation
        let mut factor = ((1.0 - no_stop_probability.powf(gene_size_codons))
            / no_stop_probability.powf(gene_size_codons))
        .ln();
        factor -= ((1.0 - no_stop_probability.powf(min_codons))
            / no_stop_probability.powf(min_codons))
        .ln();
        factor
    }
}

/// [`calculate_length_factor`] for one probability of a codon not being a
/// stop, remembering the factor of every gene length it is asked for.
///
/// The factor is a function of the gene length alone, and costs two `pow`
/// and two `ln` calls, while the start nodes of a sequence share a few
/// thousand distinct gene lengths. The factors are filled in as they are
/// asked for, by whichever thread asks first, so the two strands of a
/// sequence, and every sequence a preset model scores, share one table.
#[derive(Debug)]
pub(crate) struct LengthFactors {
    no_stop_probability: f64,
    /// Factor per gene length in bases, as the bits of an `f64`;
    /// [`Self::UNKNOWN`] until computed.
    by_gene_length: Box<[AtomicU64]>,
}

impl LengthFactors {
    /// Gene lengths below this are remembered: the lengths of genes up to
    /// [`MAX_GENE_SIZE_CODONS`], past which the factor is cheap to compute.
    const REMEMBERED_LENGTHS: usize = 3 * MAX_GENE_SIZE_CODONS as usize;
    /// Not the bits of any factor: a factor that is not a number has the
    /// bits the arithmetic gives it, which are not all ones.
    const UNKNOWN: u64 = u64::MAX;

    /// The length factors for the model of `training`.
    pub(crate) fn for_training(training: &Training) -> Self {
        Self::new(calculate_no_stop_probability(
            training.gc_content,
            training.translation_table,
        ))
    }

    fn new(no_stop_probability: f64) -> Self {
        Self {
            no_stop_probability,
            by_gene_length: (0..Self::REMEMBERED_LENGTHS)
                .map(|_| AtomicU64::new(Self::UNKNOWN))
                .collect(),
        }
    }

    #[inline]
    fn factor(&self, node: &Node) -> f64 {
        let gene_length = calculate_gene_length(node);
        let Some(remembered) = self.by_gene_length.get(gene_length) else {
            return calculate_length_factor(
                calculate_gene_size_codons(node),
                self.no_stop_probability,
            );
        };
        // Threads that find the factor unknown compute the same value, so
        // it does not matter whose is stored.
        let bits = remembered.load(Ordering::Relaxed);
        if bits != Self::UNKNOWN {
            return f64::from_bits(bits);
        }
        let factor =
            calculate_length_factor(calculate_gene_size_codons(node), self.no_stop_probability);
        remembered.store(factor.to_bits(), Ordering::Relaxed);
        factor
    }
}

/// Apply length-based adjustments to coding scores for a specific strand
///
/// This adds length factors and applies score thresholds based on gene size
fn process_strand_third_pass(
    nodes: &mut StrandNodes,
    strand: Strand,
    length_factors: &LengthFactors,
    frame_scores: &mut [f64; 3],
) {
    frame_scores.fill(CODING_SCORE_SENTINEL);

    let process_node = |node: &mut Node| {
        let frame_index = node.position.index % 3;

        if node.position.codon_type == CodonType::Stop {
            frame_scores[frame_index] = CODING_SCORE_SENTINEL;
        } else {
            let length_factor = length_factors.factor(node);
            let mut adjusted_length_factor = length_factor;

            // Apply frame-based adjustments
            if adjusted_length_factor > frame_scores[frame_index] {
                frame_scores[frame_index] = adjusted_length_factor;
            } else {
                let temp = (frame_scores[frame_index] - adjusted_length_factor)
                    .min(adjusted_length_factor)
                    .max(0.0);
                adjusted_length_factor -= temp;
            }

            // Apply length factor threshold
            if adjusted_length_factor > LENGTH_FACTOR_THRESHOLD
                && node.scores.coding_score < LENGTH_FACTOR_MULTIPLIER * adjusted_length_factor
            {
                node.scores.coding_score = LENGTH_FACTOR_MULTIPLIER * adjusted_length_factor;
            }
            node.scores.coding_score += adjusted_length_factor;
        }
    };

    if strand == Strand::Forward {
        // Process forward strand in forward order
        nodes.iter_mut().for_each(process_node);
    } else {
        // Process reverse strand in reverse order
        nodes.iter_mut().rev().for_each(process_node);
    }
}

/// Calculate raw coding scores
///
/// This is a three-pass algorithm:
/// 1. First pass: Calculate cumulative dicodon scores from stop to start
/// 2. Second pass: Normalize scores relative to maximum in each frame
/// 3. Third pass: Add length-based factors and apply thresholds
pub fn raw_coding_score(
    encoded_sequence: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &mut [Node],
    training: &Training,
) {
    let no_stop_probability =
        calculate_no_stop_probability(training.gc_content, training.translation_table);
    score_strands(
        encoded_sequence,
        reverse_complement_encoded_sequence,
        sequence_length,
        nodes,
        training,
        &LengthFactors::new(no_stop_probability),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{constants::EXPECTED_NO_STOP_PROB, types::*};
    use bio::bio_types::strand::Strand;

    fn create_test_node(
        strand: Strand,
        index: usize,
        codon_type: CodonType,
        stop_value: isize,
    ) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value,
                is_edge: false,
            },
            scores: NodeScores::default(),
            state: NodeState::default(),
            motif_info: NodeMotifInfo {
                ribosome_binding_sites: [0; 2],
                best_motif: Motif::default(),
            },
        }
    }

    fn create_test_training() -> Training {
        Training {
            gc_content: 0.5,
            translation_table: 11,
            uses_shine_dalgarno: true,
            start_type_weights: [2.0, 1.5, 1.0],
            rbs_weights: Box::new([1.0; 28]),
            upstream_composition: Box::new([[0.25; 4]; 32]),
            motif_weights: Box::new([[[1.0; 4096]; 4]; 4]),
            no_motif_weight: 0.5,
            start_weight_factor: 4.35,
            gc_bias_factors: [1.0; 3],
            gene_dicodon_table: Box::new([1.0; 4096]),
            total_dicodons: 0,
        }
    }

    #[test]
    fn test_score_nodes_basic() {
        let seq = vec![0; 60];
        let reverse_seq = vec![3; 60];
        let training = create_test_training();

        let mut nodes = vec![create_test_node(Strand::Forward, 30, CodonType::Atg, 45)];

        let encoded_sequence = EncodedSequence {
            forward_sequence: seq,
            reverse_complement_sequence: reverse_seq,
            unknown_sequence: vec![0; 60],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: 60,
        };

        let result = score_nodes(&encoded_sequence, &mut nodes, &training, false, false);

        assert!(result.is_ok());
        assert!(nodes[0].scores.total_score.is_finite());
    }

    #[test]
    fn test_score_nodes_edge_gene() {
        let seq = vec![0; 60];
        let reverse_seq = vec![3; 60];
        let training = create_test_training();

        // Create a node that will become an edge gene
        let mut nodes = vec![Node {
            position: NodePosition {
                index: 2,
                stop_value: 30,
                strand: Strand::Forward,
                codon_type: CodonType::Atg,
                is_edge: true, // Mark as edge to skip RBS scoring
            },
            scores: NodeScores::default(),
            state: NodeState::default(),
            motif_info: NodeMotifInfo {
                ribosome_binding_sites: [0; 2],
                best_motif: Motif::default(),
            },
        }];

        let encoded_sequence = EncodedSequence {
            forward_sequence: seq,
            reverse_complement_sequence: reverse_seq,
            unknown_sequence: vec![0; 60],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: 60,
        };

        let result = score_nodes(&encoded_sequence, &mut nodes, &training, false, false);

        assert!(result.is_ok());
        assert!(nodes[0].position.is_edge);
    }

    #[test]
    fn test_score_nodes_closed_mode() {
        let seq = vec![0; 40];
        let reverse_seq = vec![3; 40];
        let training = create_test_training();

        let mut nodes = vec![create_test_node(Strand::Forward, 25, CodonType::Atg, 35)];

        let encoded_sequence = EncodedSequence {
            forward_sequence: seq,
            reverse_complement_sequence: reverse_seq,
            unknown_sequence: vec![0; 40],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: 40,
        };

        let result = score_nodes(&encoded_sequence, &mut nodes, &training, true, false);

        assert!(result.is_ok());
        assert!(!nodes[0].position.is_edge);
    }

    #[test]
    fn test_score_nodes_metagenomic() {
        let seq = vec![0; 100]; // Short sequence for metagenomic
        let reverse_seq = vec![3; 100];
        let training = create_test_training();

        let mut nodes = vec![create_test_node(Strand::Forward, 40, CodonType::Atg, 60)];

        let encoded_sequence = EncodedSequence {
            forward_sequence: seq,
            reverse_complement_sequence: reverse_seq,
            unknown_sequence: vec![0; 100],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: 100,
        };

        let result = score_nodes(&encoded_sequence, &mut nodes, &training, false, true);

        assert!(result.is_ok());
        assert!(nodes[0].scores.total_score.is_finite());
    }

    #[test]
    fn test_score_nodes_non_shine_dalgarno() {
        let seq = vec![0; 60];
        let reverse_seq = vec![3; 60];
        let mut training = create_test_training();
        training.uses_shine_dalgarno = false;

        let mut nodes = vec![create_test_node(Strand::Forward, 30, CodonType::Atg, 45)];

        let encoded_sequence = EncodedSequence {
            forward_sequence: seq,
            reverse_complement_sequence: reverse_seq,
            unknown_sequence: vec![0; 60],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: 60,
        };

        let result = score_nodes(&encoded_sequence, &mut nodes, &training, false, false);

        assert!(result.is_ok());
        assert!(nodes[0].scores.ribosome_binding_score.is_finite());
    }

    #[test]
    fn test_score_nodes_stop_codon_skip() {
        let seq = vec![0; 40];
        let reverse_seq = vec![3; 40];
        let training = create_test_training();

        let mut nodes = vec![create_test_node(Strand::Forward, 30, CodonType::Stop, 30)];

        let encoded_sequence = EncodedSequence {
            forward_sequence: seq,
            reverse_complement_sequence: reverse_seq,
            unknown_sequence: vec![0; 40],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: 40,
        };

        let result = score_nodes(&encoded_sequence, &mut nodes, &training, false, false);

        assert!(result.is_ok());
        assert_eq!(nodes[0].scores.total_score, 0.0);
    }

    #[test]
    fn test_score_upstream_composition() {
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0,
            1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ];
        let training = create_test_training();
        let encoded_sequence = EncodedSequence {
            forward_sequence: seq.clone(),
            reverse_complement_sequence: seq.clone(),
            unknown_sequence: vec![0; seq.len()],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: seq.len(),
        };
        let context = ScoringContext::new(&encoded_sequence, &training, false, false);
        let node = create_test_node(Strand::Forward, 45, CodonType::Atg, 60);

        let score = score_upstream_composition(&seq, &node, &context);

        assert!(score.is_finite());
    }

    #[test]
    fn test_score_upstream_composition_reverse() {
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0,
            1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ];
        let training = create_test_training();
        let encoded_sequence = EncodedSequence {
            forward_sequence: seq.clone(),
            reverse_complement_sequence: seq.clone(),
            unknown_sequence: vec![0; seq.len()],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: seq.len(),
        };
        let context = ScoringContext::new(&encoded_sequence, &training, false, false);
        let node = create_test_node(Strand::Reverse, 45, CodonType::Atg, 30);

        let score = score_upstream_composition(&seq, &node, &context);

        assert!(score.is_finite());
    }

    #[test]
    fn test_score_upstream_composition_near_start() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let training = create_test_training();
        let encoded_sequence = EncodedSequence {
            forward_sequence: seq.clone(),
            reverse_complement_sequence: seq.clone(),
            unknown_sequence: vec![0; seq.len()],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: seq.len(),
        };
        let context = ScoringContext::new(&encoded_sequence, &training, false, false);
        let node = create_test_node(Strand::Forward, 3, CodonType::Atg, 6);

        let score = score_upstream_composition(&seq, &node, &context);

        // Should handle nodes near sequence start
        assert!(score.is_finite());
    }

    #[test]
    fn test_raw_coding_score_forward() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 12),
            create_test_node(Strand::Forward, 12, CodonType::Stop, 12),
        ];

        raw_coding_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert!(nodes[0].scores.coding_score.is_finite());
    }

    #[test]
    fn test_raw_coding_score_reverse() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Reverse, 12, CodonType::Atg, 3),
            create_test_node(Strand::Reverse, 3, CodonType::Stop, 3),
        ];

        raw_coding_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert!(nodes[0].scores.coding_score.is_finite());
    }

    #[test]
    fn test_raw_coding_score_mixed_strands() {
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ];
        let reverse_seq = vec![
            3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0,
        ];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 15),
            create_test_node(Strand::Reverse, 18, CodonType::Gtg, 6),
            create_test_node(Strand::Forward, 15, CodonType::Stop, 15),
            create_test_node(Strand::Reverse, 6, CodonType::Stop, 6),
        ];

        raw_coding_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert!(nodes[0].scores.coding_score.is_finite());
        assert!(nodes[1].scores.coding_score.is_finite());
    }

    #[test]
    fn test_raw_coding_score_different_translation_table() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let mut training = create_test_training();
        training.translation_table = 4; // Non-standard table

        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 9),
            create_test_node(Strand::Forward, 9, CodonType::Stop, 9),
        ];

        raw_coding_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        // Should handle different translation tables
        assert!(nodes[0].scores.coding_score.is_finite());
    }

    #[test]
    fn test_raw_coding_score_long_gene() {
        let seq = vec![0; 3000]; // Long sequence for large gene
        let reverse_seq = vec![3; 3000];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Forward, 100, CodonType::Atg, 2900),
            create_test_node(Strand::Forward, 2900, CodonType::Stop, 2900),
        ];

        raw_coding_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert!(nodes[0].scores.coding_score.is_finite());
    }

    #[test]
    fn test_raw_coding_score_empty_nodes() {
        let seq = vec![0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0];
        let training = create_test_training();

        let mut nodes = vec![];

        raw_coding_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        // Should handle empty nodes gracefully
        assert_eq!(nodes.len(), 0);
    }

    #[test]
    fn test_calculate_gene_length() {
        let node = create_test_node(Strand::Forward, 100, CodonType::Atg, 400);
        let length = calculate_gene_length(&node);
        assert_eq!(length, 300);

        let node_reverse = create_test_node(Strand::Reverse, 400, CodonType::Atg, 100);
        let length_reverse = calculate_gene_length(&node_reverse);
        assert_eq!(length_reverse, 300);
    }

    #[test]
    fn test_calculate_gene_size_codons() {
        let node = create_test_node(Strand::Forward, 100, CodonType::Atg, 403); // 303 bp = 102 codons
        let size_codons = calculate_gene_size_codons(&node);
        assert_eq!(size_codons, 102.0); // (303 + 3) / 3 = 102

        // Test with partial codon
        let node_partial = create_test_node(Strand::Forward, 100, CodonType::Atg, 401); // 301 bp
        let size_codons_partial = calculate_gene_size_codons(&node_partial);
        assert_eq!(size_codons_partial, 101.33333333333333); // (301 + 3) / 3 = 101.33
    }

    #[test]
    fn test_calculate_short_gene_factors() {
        let gene_length = 100;
        let (negative_factor, positive_factor) = calculate_short_gene_factors(gene_length);

        // For 100 bp gene with SHORT_GENE_THRESHOLD = 250
        assert_eq!(negative_factor, 2.5); // 250 / 100
        assert_eq!(positive_factor, 0.4); // 100 / 250

        let (neg_edge, pos_edge) = calculate_short_gene_factors(SHORT_GENE_THRESHOLD);
        assert_eq!(neg_edge, 1.0);
        assert_eq!(pos_edge, 1.0);
    }

    #[test]
    fn test_calculate_min_metagenomic_length() {
        let sequence_length = 1000;
        let min_length = calculate_min_metagenomic_length(sequence_length);
        assert_eq!(min_length, (1000.0_f64).sqrt() * 5.0); // sqrt(1000) * 5 ≈ 158.11

        // Test different sequence length
        let min_length_small = calculate_min_metagenomic_length(100);
        assert_eq!(min_length_small, 50.0); // sqrt(100) * 5 = 50
    }

    #[test]
    fn test_score_nodes_short_gene_penalty() {
        let seq = vec![0; 60];
        let reverse_seq = vec![3; 60];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Forward, 30, CodonType::Atg, 35), // Short gene (length 5)
        ];

        let encoded_sequence = EncodedSequence {
            forward_sequence: seq,
            reverse_complement_sequence: reverse_seq,
            unknown_sequence: vec![0; 60],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: 60,
        };

        let result = score_nodes(&encoded_sequence, &mut nodes, &training, false, false);

        assert!(result.is_ok());
        assert!(nodes[0].scores.total_score.is_finite());
    }

    #[test]
    fn test_score_nodes_negative_coding_score() {
        let seq = vec![0; 40];
        let reverse_seq = vec![3; 40];
        let mut training = create_test_training();

        // Set negative dicodon scores to force negative coding scores
        for score in training.gene_dicodon_table.iter_mut() {
            *score = -1.0;
        }

        let mut nodes = vec![create_test_node(Strand::Forward, 25, CodonType::Atg, 35)];

        let encoded_sequence = EncodedSequence {
            forward_sequence: seq,
            reverse_complement_sequence: reverse_seq,
            unknown_sequence: vec![0; 40],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: 40,
        };

        let result = score_nodes(&encoded_sequence, &mut nodes, &training, false, false);

        assert!(result.is_ok());
        assert!(nodes[0].scores.total_score.is_finite());
    }

    #[test]
    fn test_calculate_no_stop_probability_standard_genetic_code() {
        let gc_content = 0.5;
        let translation_table = 11;

        let prob = calculate_no_stop_probability(gc_content, translation_table);

        // Should be a valid probability between 0 and 1
        assert!(prob > 0.0 && prob < 1.0);
        assert!((prob - EXPECTED_NO_STOP_PROB).abs() < 0.01);
    }

    #[test]
    fn test_calculate_no_stop_probability_non_standard_genetic_code() {
        let gc_content = 0.3;
        let translation_table = 4; // Non-standard

        let prob = calculate_no_stop_probability(gc_content, translation_table);

        // Should be a valid probability between 0 and 1
        assert!(prob > 0.0 && prob < 1.0);
    }

    #[test]
    fn test_calculate_length_factor_normal_gene() {
        let gene_size_codons = 200.0;
        let no_stop_probability = 0.95;

        let factor = calculate_length_factor(gene_size_codons, no_stop_probability);

        // Should return a finite value
        assert!(factor.is_finite());
    }

    #[test]
    fn test_calculate_length_factor_long_gene() {
        let gene_size_codons = 1500.0; // Longer than MAX_GENE_SIZE_CODONS
        let no_stop_probability = 0.95;

        let factor = calculate_length_factor(gene_size_codons, no_stop_probability);

        // Should return a finite value with scaling applied
        assert!(factor.is_finite());
    }

    #[test]
    fn test_process_strand_first_pass_forward() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let training = create_test_training();
        let mut frame_scores = [0.0; 3];
        let mut last_positions = [0usize; 3];

        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 12),
            create_test_node(Strand::Forward, 12, CodonType::Stop, 12),
        ];

        let (mut forward, _) = split_by_strand(&mut nodes);
        process_strand_first_pass(
            &mut forward,
            &seq,
            seq.len(),
            Strand::Forward,
            &training,
            &mut frame_scores,
            &mut last_positions,
        );

        // Should have calculated coding scores
        assert!(nodes[0].scores.coding_score.is_finite());
    }

    #[test]
    fn test_process_strand_second_pass_forward() {
        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 12),
            create_test_node(Strand::Forward, 6, CodonType::Gtg, 15),
        ];

        // Set initial coding scores
        nodes[0].scores.coding_score = 10.0;
        nodes[1].scores.coding_score = 5.0;

        let mut frame_scores = [0.0; 3];

        let (mut forward, _) = split_by_strand(&mut nodes);
        process_strand_second_pass(&mut forward, Strand::Forward, &mut frame_scores);

        // Should have adjusted scores based on maximums
        assert!(nodes[0].scores.coding_score.is_finite());
        assert!(nodes[1].scores.coding_score.is_finite());
    }

    #[test]
    fn test_process_strand_third_pass_forward() {
        let mut nodes = vec![create_test_node(Strand::Forward, 3, CodonType::Atg, 300)];

        // Set initial coding score
        nodes[0].scores.coding_score = 5.0;

        let no_stop_probability = 0.95;
        let mut frame_scores = [0.0; 3];

        let (mut forward, _) = split_by_strand(&mut nodes);
        process_strand_third_pass(
            &mut forward,
            Strand::Forward,
            &LengthFactors::new(no_stop_probability),
            &mut frame_scores,
        );

        // Should have added length factor to coding score
        assert!(nodes[0].scores.coding_score > 5.0);
        assert!(nodes[0].scores.coding_score.is_finite());
    }

    #[test]
    fn length_factors_remember_what_calculate_length_factor_gives() {
        let no_stop_probability = calculate_no_stop_probability(0.43, 11);
        let length_factors = LengthFactors::new(no_stop_probability);
        // Lengths on both sides of the longest remembered one, each asked
        // for twice: computed, then remembered.
        for gene_length in [
            0usize, 3, 87, 90, 1500, 2996, 2997, 2999, 3000, 3003, 90_000,
        ] {
            let node = create_test_node(
                Strand::Forward,
                30,
                CodonType::Atg,
                30 + gene_length as isize,
            );
            let expected =
                calculate_length_factor(((gene_length + 3) as f64) / 3.0, no_stop_probability);
            for _ in 0..2 {
                assert_eq!(
                    length_factors.factor(&node).to_bits(),
                    expected.to_bits(),
                    "gene length {gene_length}"
                );
            }
        }
    }

    /// The neighbour search as a walk over the nodes: a forward start among
    /// the first nodes looks back for an edge node with its stop value, a
    /// reverse start among the last nodes looks ahead.
    fn has_edge_neighbour_walking(nodes: &[Node], current_index: usize) -> bool {
        let current = &nodes[current_index];
        let is_edge_of_frame = |node: &Node| {
            node.position.is_edge && node.position.stop_value == current.position.stop_value
        };
        if current_index < NODE_SEARCH_WINDOW && current.position.strand == Strand::Forward {
            nodes[..current_index].iter().any(is_edge_of_frame)
        } else if current_index >= nodes.len().saturating_sub(NODE_SEARCH_WINDOW)
            && current.position.strand == Strand::Reverse
        {
            nodes[current_index + 1..].iter().any(is_edge_of_frame)
        } else {
            false
        }
    }

    #[test]
    fn edge_stops_find_the_edge_neighbours_a_walk_finds() {
        let params = ScoringParams::from_training(&create_test_training());
        let mut random = crate::test_random::XorShift(0xA076_1D64_78BD_642F_u64);
        let mut next = || random.next();

        // Fewer nodes than the search window, exactly as many, and enough
        // for the middle nodes to search neither end.
        let mut bonuses = 0;
        for node_count in [1, 2, 37, NODE_SEARCH_WINDOW, 3 * NODE_SEARCH_WINDOW + 5] {
            let mut nodes: Vec<Node> = (0..node_count)
                .map(|index| {
                    let strand = if next() % 2 == 0 {
                        Strand::Forward
                    } else {
                        Strand::Reverse
                    };
                    let codon_type = if next() % 4 == 0 {
                        CodonType::Stop
                    } else {
                        CodonType::Atg
                    };
                    // Few distinct stop values, so that nodes share them.
                    let mut node =
                        create_test_node(strand, 3 * index, codon_type, (next() % 12) as isize);
                    node.position.is_edge = next() % 9 == 0;
                    node
                })
                .collect();

            let mut edge_stops = EdgeStops::new(&nodes);
            for index in 0..node_count {
                if nodes[index].position.codon_type == CodonType::Stop {
                    edge_stops.scored(index, &nodes[index]);
                    continue;
                }
                let boundary_bonus_applied = next() % 11 == 0;
                let expected = !nodes[index].position.is_edge
                    && !boundary_bonus_applied
                    && has_edge_neighbour_walking(&nodes, index);

                nodes[index].scores.upstream_score = 1.0;
                apply_edge_upstream_bonus_for_neighbors(
                    &mut nodes[index],
                    index,
                    node_count,
                    &params,
                    boundary_bonus_applied,
                    &edge_stops,
                );
                let bonus_applied = nodes[index].scores.upstream_score != 1.0;
                assert_eq!(bonus_applied, expected, "{node_count} nodes, node {index}");
                bonuses += usize::from(bonus_applied);

                // Scoring turns some start nodes into edge nodes, which the
                // nodes after them then see.
                if next() % 7 == 0 {
                    nodes[index].position.is_edge = true;
                }
                edge_stops.scored(index, &nodes[index]);
            }
        }
        assert!(bonuses > 50, "only {bonuses} bonuses");
    }
}
