use crate::{
    constants::{
        CODING_SCORE_SENTINEL, CODING_SCORE_THRESHOLD, DICODON_SIZE, EDGE_BONUS,
        EDGE_POSITION_OFFSET, EDGE_POSITION_THRESHOLD, EDGE_UPSTREAM, GENE_SIZE_SCALING_FACTOR,
        LENGTH_FACTOR_MULTIPLIER, LENGTH_FACTOR_THRESHOLD, MAX_GENE_SIZE_CODONS,
        METAGENOMIC_LENGTH_THRESHOLD, METAGENOMIC_MIN_LENGTH_FACTOR, METAGENOMIC_PENALTY,
        METAGENOMIC_PENALTY_DIVISOR, MIN_GENE_SIZE_CODONS, MIN_META_GENE_LENGTH,
        NEGATIVE_SCORE_PENALTY, NO_MOTIF_THRESHOLD, NODE_SEARCH_WINDOW, SHORT_GENE_THRESHOLD,
        UPSTREAM_COMPOSITION_WEIGHT, UPSTREAM_SCAN_RANGE, UPSTREAM_SKIP_END, UPSTREAM_SKIP_START,
    },
    node::{calc_orf_gc, find_best_upstream_motif},
    rbs_score,
    sequence::encoded::EncodedSequence,
    sequence::{calculate_kmer_index, is_stop},
    types::{CodonType, Node, OrphosError, Training},
};

use bio::bio_types::strand::Strand;

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
    /// Edge upstream bonus factor
    pub edge_upstream_bonus: f64,
    /// Edge bonus factor
    pub edge_bonus_factor: f64,
    /// Penalty factor for various score adjustments
    pub penalty_factor: f64,
    /// Cached start weight factor
    pub start_weight_factor: f64,
    /// Cached metagenomic penalty coefficient
    pub metagenomic_penalty_coeff: f64,
    /// Cached no-stop probability for length calculations
    pub no_stop_probability: f64,
}

impl ScoringParams {
    /// Create scoring parameters from training data
    pub fn from_training(training: &Training) -> Self {
        let start_weight_factor = training.start_weight_factor;
        let edge_upstream_bonus = EDGE_UPSTREAM * start_weight_factor;
        let edge_bonus_factor = EDGE_BONUS * start_weight_factor;
        let penalty_factor = NEGATIVE_SCORE_PENALTY * edge_bonus_factor;
        let metagenomic_penalty_coeff = METAGENOMIC_PENALTY / METAGENOMIC_PENALTY_DIVISOR;
        let no_stop_probability =
            calculate_no_stop_probability(training.gc_content, training.translation_table);

        Self {
            edge_upstream_bonus,
            edge_bonus_factor,
            penalty_factor,
            start_weight_factor,
            metagenomic_penalty_coeff,
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
        node.scores.upstream_score += params.edge_upstream_bonus;
    }
}

/// Apply edge upstream bonus based on neighboring nodes
///
/// Searches for edge nodes with the same stop position within a limited window
/// and applies bonus to the current node if found
fn apply_edge_upstream_bonus_for_neighbors(
    nodes: &mut [Node],
    current_index: usize,
    node_count: usize,
    params: &ScoringParams,
) {
    let current_node = &nodes[current_index];

    if current_node.position.is_edge {
        return;
    }

    if current_index < NODE_SEARCH_WINDOW && current_node.position.strand == Strand::Forward {
        // Forward strand: check previous nodes
        for j in (0..current_index).rev() {
            if nodes[j].position.is_edge
                && current_node.position.stop_value == nodes[j].position.stop_value
            {
                nodes[current_index].scores.upstream_score += params.edge_upstream_bonus;
                break;
            }
        }
    } else if current_index >= node_count.saturating_sub(NODE_SEARCH_WINDOW)
        && current_node.position.strand == Strand::Reverse
    {
        // Reverse strand: check following nodes
        for j in (current_index + 1)..node_count {
            if nodes[j].position.is_edge
                && current_node.position.stop_value == nodes[j].position.stop_value
            {
                nodes[current_index].scores.upstream_score += params.edge_upstream_bonus;
                break;
            }
        }
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
    if !node.position.is_edge && edge_gene == 1 {
        node.scores.upstream_score -= params.penalty_factor;
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
    if cache.needs_metagenomic_penalty && edge_gene == 0 {
        let penalty = params.metagenomic_penalty_coeff
            * (METAGENOMIC_LENGTH_THRESHOLD - context.encoded_sequence.sequence_length) as f64;
        node.scores.coding_score -= penalty;
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
                let penalty =
                    10.31f64.mul_add(-(context.encoded_sequence.sequence_length as f64), 0.004);
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

/// Score nodes for gene finding
pub fn score_nodes(
    encoded_sequence: &EncodedSequence,
    nodes: &mut [Node],
    training: &Training,
    closed: bool,
    is_meta: bool,
) -> Result<(), OrphosError> {
    let node_count = nodes.len();

    // Early exit for empty nodes
    if node_count == 0 {
        return Ok(());
    }

    // Create context structs to reduce parameter passing
    let context = ScoringContext::new(encoded_sequence, training, closed, is_meta);
    let params = ScoringParams::from_training(training);

    calc_orf_gc(
        &encoded_sequence.forward_sequence,
        encoded_sequence.sequence_length,
        nodes,
    );
    raw_coding_score_optimized(
        &encoded_sequence.forward_sequence,
        &encoded_sequence.reverse_complement_sequence,
        encoded_sequence.sequence_length,
        nodes,
        training,
        params.no_stop_probability,
    );

    // Calculate RBS scores using optimized approach
    if training.uses_shine_dalgarno {
        rbs_score(
            &encoded_sequence.forward_sequence,
            &encoded_sequence.reverse_complement_sequence,
            encoded_sequence.sequence_length,
            nodes,
            training,
        );
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
    for i in 0..node_count {
        if nodes[i].position.codon_type == CodonType::Stop {
            continue;
        }

        // Create cache for repetitive calculations
        let cache = NodeCache::new(&nodes[i], &context);

        // Calculate edge gene status
        let mut edge_gene = calculate_edge_gene_status(&nodes[i], &context);

        // Calculate various score components
        calculate_type_scores(&mut nodes[i], edge_gene, &params, training);

        calculate_rbs_scores(&mut nodes[i], training, &params);

        calculate_upstream_scores(&mut nodes[i], &context, &params, &cache);

        // Apply edge upstream bonus based on neighboring nodes
        apply_edge_upstream_bonus_for_neighbors(nodes, i, node_count, &params);

        // Convert boundary nodes to edge genes if necessary
        convert_boundary_nodes_to_edge_genes(
            &mut nodes[i],
            &context,
            &mut edge_gene,
            &params,
            &cache,
        );

        // Calculate base start score before applying penalties
        calculate_final_scores(&mut nodes[i]);

        // Apply penalties and bonuses (includes recalculating start score)
        apply_penalties_and_bonuses(&mut nodes[i], edge_gene, &context, &params, &cache);

        // Calculate final total score
        nodes[i].scores.total_score = nodes[i].scores.coding_score + nodes[i].scores.start_score;
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

        // Add to upstream score
        upstream_score += UPSTREAM_COMPOSITION_WEIGHT
            * context.training.start_weight_factor
            * context.training.upstream_composition[count][mer_idx];
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

/// Process nodes for a specific strand in the first pass (dicodon scoring)
///
/// This calculates cumulative dicodon scores from stop codons to start codons
fn process_strand_first_pass(
    nodes: &mut [Node],
    encoded_sequence: &[u8],
    sequence_length: usize,
    strand: Strand,
    training: &Training,
    frame_scores: &mut [f64; 3],
    last_positions: &mut [usize; 3],
) {
    frame_scores.fill(0.0);

    if strand == Strand::Forward {
        // Process forward strand in reverse order (from end to start)
        for i in (0..nodes.len()).rev() {
            if nodes[i].position.strand != Strand::Forward {
                continue;
            }

            let frame_index = nodes[i].position.index % 3;

            if nodes[i].position.codon_type == CodonType::Stop {
                last_positions[frame_index] = nodes[i].position.index;
                frame_scores[frame_index] = 0.0;
            } else {
                // Calculate dicodon scores from last stop to current start
                if last_positions[frame_index] >= 3 {
                    let mut j = last_positions[frame_index] - 3;
                    loop {
                        if j < nodes[i].position.index {
                            break;
                        }
                        let mer_idx = calculate_kmer_index(DICODON_SIZE, encoded_sequence, j);
                        let dicodon_score = training.gene_dicodon_table[mer_idx];
                        frame_scores[frame_index] += dicodon_score;

                        if j < 3 {
                            break;
                        }
                        j -= 3;
                    }
                }
                nodes[i].scores.coding_score = frame_scores[frame_index];
                last_positions[frame_index] = nodes[i].position.index;
            }
        }
    } else {
        // Process reverse strand in forward order
        for node in nodes.iter_mut() {
            if node.position.strand != Strand::Reverse {
                continue;
            }

            let frame_index = node.position.index % 3;

            if node.position.codon_type == CodonType::Stop {
                last_positions[frame_index] = node.position.index;
                frame_scores[frame_index] = 0.0;
            } else {
                // Calculate dicodon scores from last stop to current start
                let mut j = last_positions[frame_index] + 3;
                while j <= node.position.index {
                    let reverse_seq_pos = sequence_length - j - 1;
                    let mer_idx =
                        calculate_kmer_index(DICODON_SIZE, encoded_sequence, reverse_seq_pos);
                    let dicodon_score = training.gene_dicodon_table[mer_idx];
                    frame_scores[frame_index] += dicodon_score;
                    j += 3;
                }
                node.scores.coding_score = frame_scores[frame_index];
                last_positions[frame_index] = node.position.index;
            }
        }
    }
}

/// Process nodes for a specific strand in the second pass (score normalization)
///
/// This normalizes scores relative to the maximum score in each frame
fn process_strand_second_pass(nodes: &mut [Node], strand: Strand, frame_scores: &mut [f64; 3]) {
    frame_scores.fill(CODING_SCORE_SENTINEL);

    if strand == Strand::Forward {
        // Process forward strand in forward order
        for node in nodes.iter_mut() {
            if node.position.strand != Strand::Forward {
                continue;
            }

            let frame_index = node.position.index % 3;

            if node.position.codon_type == CodonType::Stop {
                frame_scores[frame_index] = CODING_SCORE_SENTINEL;
            } else if node.scores.coding_score > frame_scores[frame_index] {
                frame_scores[frame_index] = node.scores.coding_score;
            } else {
                node.scores.coding_score -= frame_scores[frame_index] - node.scores.coding_score;
            }
        }
    } else {
        // Process reverse strand in reverse order
        for i in (0..nodes.len()).rev() {
            if nodes[i].position.strand != Strand::Reverse {
                continue;
            }

            let frame_index = nodes[i].position.index % 3;

            if nodes[i].position.codon_type == CodonType::Stop {
                frame_scores[frame_index] = CODING_SCORE_SENTINEL;
            } else if nodes[i].scores.coding_score > frame_scores[frame_index] {
                frame_scores[frame_index] = nodes[i].scores.coding_score;
            } else {
                nodes[i].scores.coding_score -=
                    frame_scores[frame_index] - nodes[i].scores.coding_score;
            }
        }
    }
}

/// Calculate the length-based scoring factor for a gene
///
/// This implements the biological principle that very short or very long genes
/// are less likely to be real genes
fn calculate_length_factor(gene_size_codons: f64, no_stop_probability: f64) -> f64 {
    if gene_size_codons > MAX_GENE_SIZE_CODONS as f64 {
        // For very long genes, use a scaled length factor
        let mut factor = ((1.0 - no_stop_probability.powi(MAX_GENE_SIZE_CODONS))
            / no_stop_probability.powi(MAX_GENE_SIZE_CODONS))
        .ln();
        factor -= ((1.0 - no_stop_probability.powi(MIN_GENE_SIZE_CODONS))
            / no_stop_probability.powi(MIN_GENE_SIZE_CODONS))
        .ln();
        factor *= (gene_size_codons - MIN_GENE_SIZE_CODONS as f64) / GENE_SIZE_SCALING_FACTOR;
        factor
    } else {
        // For normal-sized genes, use direct calculation
        let mut factor = ((1.0 - no_stop_probability.powf(gene_size_codons))
            / no_stop_probability.powf(gene_size_codons))
        .ln();
        factor -= ((1.0 - no_stop_probability.powi(MIN_GENE_SIZE_CODONS))
            / no_stop_probability.powi(MIN_GENE_SIZE_CODONS))
        .ln();
        factor
    }
}

/// Apply length-based adjustments to coding scores for a specific strand
///
/// This adds length factors and applies score thresholds based on gene size
fn process_strand_third_pass(
    nodes: &mut [Node],
    strand: Strand,
    no_stop_probability: f64,
    frame_scores: &mut [f64; 3],
) {
    frame_scores.fill(CODING_SCORE_SENTINEL);

    let process_node = |node: &mut Node, frame_scores: &mut [f64; 3]| {
        if node.position.strand != strand {
            return;
        }

        let frame_index = node.position.index % 3;

        if node.position.codon_type == CodonType::Stop {
            frame_scores[frame_index] = CODING_SCORE_SENTINEL;
        } else {
            let gene_size_codons = calculate_gene_size_codons(node);
            let length_factor = calculate_length_factor(gene_size_codons, no_stop_probability);
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
        for node in nodes.iter_mut() {
            process_node(node, frame_scores);
        }
    } else {
        // Process reverse strand in reverse order
        for i in (0..nodes.len()).rev() {
            process_node(&mut nodes[i], frame_scores);
        }
    }
}

/// Calculate raw coding scores with cached no_stop_probability
///
/// This is a three-pass algorithm:
/// 1. First pass: Calculate cumulative dicodon scores from stop to start
/// 2. Second pass: Normalize scores relative to maximum in each frame
/// 3. Third pass: Add length-based factors and apply thresholds
pub fn raw_coding_score_optimized(
    encoded_sequence: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &mut [Node],
    training: &Training,
    no_stop_probability: f64,
) {
    let mut frame_scores = [0.0; 3];
    let mut last_positions = [0usize; 3];

    // FIRST PASS: Score coding potential (start->stop) for both strands
    process_strand_first_pass(
        nodes,
        encoded_sequence,
        sequence_length,
        Strand::Forward,
        training,
        &mut frame_scores,
        &mut last_positions,
    );

    process_strand_first_pass(
        nodes,
        reverse_complement_encoded_sequence,
        sequence_length,
        Strand::Reverse,
        training,
        &mut frame_scores,
        &mut last_positions,
    );

    // SECOND PASS: Normalize scores relative to maximum in each frame
    process_strand_second_pass(nodes, Strand::Forward, &mut frame_scores);
    process_strand_second_pass(nodes, Strand::Reverse, &mut frame_scores);

    // THIRD PASS: Add length-based factors and apply thresholds
    process_strand_third_pass(
        nodes,
        Strand::Forward,
        no_stop_probability,
        &mut frame_scores,
    );
    process_strand_third_pass(
        nodes,
        Strand::Reverse,
        no_stop_probability,
        &mut frame_scores,
    );
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
    let mut frame_scores = [0.0; 3];
    let mut last_positions = [0usize; 3];

    // FIRST PASS: Score coding potential (start->stop) for both strands
    process_strand_first_pass(
        nodes,
        encoded_sequence,
        sequence_length,
        Strand::Forward,
        training,
        &mut frame_scores,
        &mut last_positions,
    );

    process_strand_first_pass(
        nodes,
        reverse_complement_encoded_sequence,
        sequence_length,
        Strand::Reverse,
        training,
        &mut frame_scores,
        &mut last_positions,
    );

    // SECOND PASS: Normalize scores relative to maximum in each frame
    process_strand_second_pass(nodes, Strand::Forward, &mut frame_scores);
    process_strand_second_pass(nodes, Strand::Reverse, &mut frame_scores);

    // THIRD PASS: Add length-based factors and apply thresholds
    process_strand_third_pass(
        nodes,
        Strand::Forward,
        no_stop_probability,
        &mut frame_scores,
    );
    process_strand_third_pass(
        nodes,
        Strand::Reverse,
        no_stop_probability,
        &mut frame_scores,
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

        process_strand_first_pass(
            &mut nodes,
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

        process_strand_second_pass(&mut nodes, Strand::Forward, &mut frame_scores);

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

        process_strand_third_pass(
            &mut nodes,
            Strand::Forward,
            no_stop_probability,
            &mut frame_scores,
        );

        // Should have added length factor to coding score
        assert!(nodes[0].scores.coding_score > 5.0);
        assert!(nodes[0].scores.coding_score.is_finite());
    }
}
