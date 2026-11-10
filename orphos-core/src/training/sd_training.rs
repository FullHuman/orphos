use bio::bio_types::strand::Strand;

use crate::{
    parallel,
    training::{
        GC_HIGH_AT_FREQ, GC_HIGH_GC_FREQ, GC_LOW_AT_FREQ, GC_LOW_GC_FREQ, GENE_RATIO_THRESHOLD,
        INITIAL_SCORE_THRESHOLD, MAX_GC_CONTENT, MAX_TRAINING_ITERATIONS_SD, MIN_GC_CONTENT,
        NUM_BASES, NUM_CODON_TYPES, NUM_RBS_WEIGHTS, RBS_WEIGHT_THRESHOLD_HIGH, THRESHOLD_DIVISOR,
        UPSTREAM_POSITIONS, WEIGHT_CLAMP_MAX, WEIGHT_CLAMP_MIN, count_upstream_composition_into,
    },
    types::{CodonType, Node, Training},
};

/// The node fields read by [`train_starts_sd`], packed so that the nodes of
/// a strand are scanned in a few megabytes: every training iteration scans
/// them all.
#[derive(Debug, Clone, Copy)]
struct SdNode {
    coding_score: f64,
    /// The two ribosome binding site motifs of the node, as one index into
    /// [`RbsChoices`].
    rbs_pair: u16,
    codon_type_index: u8,
    frame: u8,
    is_stop: bool,
    is_edge: bool,
}

impl From<&Node> for SdNode {
    fn from(node: &Node) -> Self {
        let [exact, mismatch] = node.motif_info.ribosome_binding_sites;
        assert!(
            exact < NUM_RBS_WEIGHTS && mismatch < NUM_RBS_WEIGHTS,
            "ribosome binding site motif out of range"
        );
        Self {
            coding_score: node.scores.coding_score,
            rbs_pair: (exact * NUM_RBS_WEIGHTS + mismatch) as u16,
            codon_type_index: node.position.codon_type.to_index() as u8,
            frame: (node.position.index % NUM_CODON_TYPES) as u8,
            is_stop: node.position.codon_type == CodonType::Stop,
            is_edge: node.position.is_edge,
        }
    }
}

/// The nodes of one strand, in node order.
struct SdStrand {
    nodes: Vec<SdNode>,
    /// Sequence position of each node.
    positions: Vec<usize>,
}

impl SdStrand {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(capacity),
            positions: Vec::with_capacity(capacity),
        }
    }

    fn push(&mut self, node: &Node) {
        self.nodes.push(SdNode::from(node));
        self.positions.push(node.position.index);
    }
}

/// For every pair of ribosome binding site motifs a node can have (exact
/// match, match with one mismatch), the motif an iteration's weights choose
/// and its weight.
struct RbsChoices {
    site: [u8; NUM_RBS_WEIGHTS * NUM_RBS_WEIGHTS],
    weight: [f64; NUM_RBS_WEIGHTS * NUM_RBS_WEIGHTS],
}

impl RbsChoices {
    fn new(rbs_weights: &[f64; NUM_RBS_WEIGHTS]) -> Self {
        let mut choices = Self {
            site: [0; NUM_RBS_WEIGHTS * NUM_RBS_WEIGHTS],
            weight: [0.0; NUM_RBS_WEIGHTS * NUM_RBS_WEIGHTS],
        };
        for exact in 0..NUM_RBS_WEIGHTS {
            for mismatch in 0..NUM_RBS_WEIGHTS {
                let site = if rbs_weights[exact] > rbs_weights[mismatch] + RBS_WEIGHT_THRESHOLD_HIGH
                    || mismatch == 0
                {
                    exact
                } else if rbs_weights[exact] < rbs_weights[mismatch] - RBS_WEIGHT_THRESHOLD_HIGH
                    || exact == 0
                {
                    mismatch
                } else {
                    exact.max(mismatch)
                };
                let pair = exact * NUM_RBS_WEIGHTS + mismatch;
                choices.site[pair] = site as u8;
                choices.weight[pair] = rbs_weights[site];
            }
        }
        choices
    }
}

/// Train start codon recognition using Shine-Dalgarno motifs
pub fn train_starts_sd(
    sequence: &[u8],
    reverse_sequence: &[u8],
    sequence_length: usize,
    nodes: &[Node],
    training: &mut Training,
) {
    // The iterations below only read a few fields of every node; compact
    // copies keep the repeated scans cache friendly, one per strand since
    // every scan is of one strand. (`add_nodes` gives every node a strand.)
    let node_count = nodes.len();
    let mut forward_nodes = SdStrand::with_capacity(node_count / 2 + 1);
    let mut reverse_nodes = SdStrand::with_capacity(node_count / 2 + 1);
    for node in nodes {
        match node.position.strand {
            Strand::Forward => forward_nodes.push(node),
            Strand::Reverse => reverse_nodes.push(node),
            Strand::Unknown => {}
        }
    }
    let all_nodes = || forward_nodes.nodes.iter().chain(&reverse_nodes.nodes);
    let weight_factor = training.start_weight_factor;
    let mut type_background = [0.0; NUM_CODON_TYPES];
    let mut ribosome_binding_site_background = [0.0; NUM_RBS_WEIGHTS];

    training.start_type_weights = [0.0; NUM_CODON_TYPES];
    *training.rbs_weights = [0.0; NUM_RBS_WEIGHTS];
    *training.upstream_composition = [[0.0; NUM_BASES]; UPSTREAM_POSITIONS];

    let mut total_sum = 0.0;
    for node in all_nodes() {
        if node.is_stop {
            continue;
        }
        type_background[usize::from(node.codon_type_index)] += 1.0;
        total_sum += 1.0;
    }

    if total_sum > 0.0 {
        for weight in type_background.iter_mut().take(NUM_CODON_TYPES) {
            *weight /= total_sum;
        }
    }

    let mut score_threshold = INITIAL_SCORE_THRESHOLD;

    for iteration in 0..MAX_TRAINING_ITERATIONS_SD {
        // The motif a node counts for depends on the weights of the
        // previous iteration and on its pair of motifs only.
        let rbs_choices = RbsChoices::new(&training.rbs_weights);

        // The two strands are scanned on two threads. Every count is a whole
        // number, so adding the strands' counts gives exactly what one scan
        // over both strands gave.
        let count_upstream = iteration == MAX_TRAINING_ITERATIONS_SD - 1;
        let (forward, reverse) = parallel::join(
            || {
                scan_strand(
                    &forward_nodes,
                    Strand::Forward,
                    sequence,
                    sequence_length,
                    training,
                    &rbs_choices,
                    weight_factor,
                    score_threshold,
                    count_upstream,
                )
            },
            || {
                scan_strand(
                    &reverse_nodes,
                    Strand::Reverse,
                    reverse_sequence,
                    sequence_length,
                    training,
                    &rbs_choices,
                    weight_factor,
                    score_threshold,
                    count_upstream,
                )
            },
        );

        // Background: the motif of every start node that is not an edge.
        for (total, (forward_count, reverse_count)) in ribosome_binding_site_background
            .iter_mut()
            .zip(forward.rbs_background.iter().zip(&reverse.rbs_background))
        {
            *total = forward_count + reverse_count;
        }
        let ribosome_binding_site_background_sum =
            forward.rbs_background_sum + reverse.rbs_background_sum;
        if ribosome_binding_site_background_sum > 0.0 {
            for weight in ribosome_binding_site_background
                .iter_mut()
                .take(NUM_RBS_WEIGHTS)
            {
                *weight /= ribosome_binding_site_background_sum;
            }
        }

        let mut ribosome_binding_site_real = [0.0; NUM_RBS_WEIGHTS];
        for (total, (forward_count, reverse_count)) in ribosome_binding_site_real
            .iter_mut()
            .zip(forward.rbs_real.iter().zip(&reverse.rbs_real))
        {
            *total = forward_count + reverse_count;
        }
        let mut type_real = [0.0; NUM_CODON_TYPES];
        for (total, (forward_count, reverse_count)) in type_real
            .iter_mut()
            .zip(forward.type_real.iter().zip(&reverse.type_real))
        {
            *total = forward_count + reverse_count;
        }
        if count_upstream {
            for (position, (forward_counts, reverse_counts)) in
                training.upstream_composition.iter_mut().zip(
                    forward
                        .upstream_composition
                        .iter()
                        .zip(&reverse.upstream_composition),
                )
            {
                for (total, (forward_count, reverse_count)) in position
                    .iter_mut()
                    .zip(forward_counts.iter().zip(reverse_counts))
                {
                    *total += forward_count + reverse_count;
                }
            }
        }

        let total_real_ribosome_binding_sites: f64 = ribosome_binding_site_real.iter().sum();
        if total_real_ribosome_binding_sites == 0.0 {
            training
                .rbs_weights
                .iter_mut()
                .for_each(|weight| *weight = 0.0);
        } else {
            ribosome_binding_site_real
                .iter_mut()
                .for_each(|val| *val /= total_real_ribosome_binding_sites);

            training
                .rbs_weights
                .iter_mut()
                .enumerate()
                .for_each(|(weight_index, weight)| {
                    *weight = if ribosome_binding_site_background[weight_index] != 0.0 {
                        (ribosome_binding_site_real[weight_index]
                            / ribosome_binding_site_background[weight_index])
                            .ln()
                    } else {
                        WEIGHT_CLAMP_MIN
                    };
                    *weight = weight.clamp(WEIGHT_CLAMP_MIN, WEIGHT_CLAMP_MAX);
                });
        }

        let total_real_types: f64 = type_real.iter().sum();

        if total_real_types == 0.0 {
            training
                .start_type_weights
                .iter_mut()
                .for_each(|weight| *weight = 0.0);
        } else {
            type_real
                .iter_mut()
                .for_each(|val| *val /= total_real_types);

            training
                .start_type_weights
                .iter_mut()
                .enumerate()
                .for_each(|(type_index, weight)| {
                    *weight = if type_background[type_index] != 0.0 {
                        (type_real[type_index] / type_background[type_index]).ln()
                    } else {
                        WEIGHT_CLAMP_MIN
                    };
                    *weight = weight.clamp(WEIGHT_CLAMP_MIN, WEIGHT_CLAMP_MAX);
                });
        }

        if total_real_types <= node_count as f64 / GENE_RATIO_THRESHOLD {
            score_threshold /= THRESHOLD_DIVISOR;
        }
    }

    training
        .upstream_composition
        .iter_mut()
        .enumerate()
        .take(UPSTREAM_POSITIONS)
        .for_each(|(_, position_data)| {
            let position_sum: f64 = position_data.iter().sum();

            if position_sum == 0.0 {
                position_data.fill(0.0);
            } else {
                for (base_index, base_value) in position_data.iter_mut().enumerate().take(NUM_BASES)
                {
                    let original_frequency = *base_value / position_sum;
                    *base_value = original_frequency;

                    let (background_denominator, _base_name) = if training.gc_content
                        > MIN_GC_CONTENT
                        && training.gc_content < MAX_GC_CONTENT
                    {
                        if base_index == 0 || base_index == 3 {
                            (
                                1.0 - training.gc_content,
                                if base_index == 0 { "A" } else { "T" },
                            )
                        } else {
                            (training.gc_content, if base_index == 1 { "G" } else { "C" })
                        }
                    } else if training.gc_content <= MIN_GC_CONTENT {
                        if base_index == 0 || base_index == 3 {
                            (GC_LOW_AT_FREQ, if base_index == 0 { "A" } else { "T" })
                        } else {
                            (GC_LOW_GC_FREQ, if base_index == 1 { "G" } else { "C" })
                        }
                    } else if base_index == 0 || base_index == 3 {
                        (GC_HIGH_AT_FREQ, if base_index == 0 { "A" } else { "T" })
                    } else {
                        (GC_HIGH_GC_FREQ, if base_index == 1 { "G" } else { "C" })
                    };

                    let log_ratio = (*base_value * 2.0 / background_denominator).ln();
                    *base_value = log_ratio;
                    *base_value = base_value.clamp(WEIGHT_CLAMP_MIN, WEIGHT_CLAMP_MAX);
                }
            }
        });
}

/// What one scan of a strand counts in an iteration of [`train_starts_sd`].
struct StrandCounts {
    /// The motif of every start node that is not an edge.
    rbs_background: [f64; NUM_RBS_WEIGHTS],
    /// The number of those nodes.
    rbs_background_sum: f64,
    rbs_real: [f64; NUM_RBS_WEIGHTS],
    type_real: [f64; NUM_CODON_TYPES],
    upstream_composition: [[f64; NUM_BASES]; UPSTREAM_POSITIONS],
}

/// One iteration's scan of `strand_nodes`, the nodes of `strand`: for every
/// frame, the best start in the ORF ending at each stop is counted, with its
/// RBS motif and start codon type, when its score reaches `score_threshold`;
/// with `count_upstream`, the bases upstream of those starts are counted
/// too. The motif of every start that is not an edge is counted as
/// background. The forward strand is scanned from the first node, the
/// reverse strand from the last.
#[allow(clippy::too_many_arguments)]
fn scan_strand(
    strand_nodes: &SdStrand,
    strand: Strand,
    sequence: &[u8],
    sequence_length: usize,
    training: &Training,
    rbs_choices: &RbsChoices,
    weight_factor: f64,
    score_threshold: f64,
    count_upstream: bool,
) -> StrandCounts {
    let nodes = strand_nodes.nodes.as_slice();
    let mut counts = StrandCounts {
        rbs_background: [0.0; NUM_RBS_WEIGHTS],
        rbs_background_sum: 0.0,
        rbs_real: [0.0; NUM_RBS_WEIGHTS],
        type_real: [0.0; NUM_CODON_TYPES],
        upstream_composition: [[0.0; NUM_BASES]; UPSTREAM_POSITIONS],
    };
    let mut best_scores = [0.0; NUM_CODON_TYPES];
    let mut best_node_indices: [Option<usize>; NUM_CODON_TYPES] = [None; NUM_CODON_TYPES];
    let mut ribosome_binding_sites = [0; NUM_CODON_TYPES];
    let mut codon_types = [0; NUM_CODON_TYPES];

    let mut visit = |node_index: usize, node: &SdNode| {
        if !node.is_stop && node.is_edge {
            return;
        }

        let frame = usize::from(node.frame);

        if node.is_stop {
            if best_scores[frame] >= score_threshold && best_node_indices[frame].is_some() {
                let best_index = best_node_indices[frame].unwrap();
                if best_index < nodes.len() && usize::from(nodes[best_index].frame) == frame {
                    counts.rbs_real[ribosome_binding_sites[frame]] += 1.0;
                    counts.type_real[codon_types[frame]] += 1.0;

                    if count_upstream {
                        count_upstream_composition_into(
                            sequence,
                            sequence_length,
                            strand,
                            strand_nodes.positions[best_index],
                            &mut counts.upstream_composition,
                        );
                    }
                }
            }

            best_scores[frame] = 0.0;
            best_node_indices[frame] = None;
            ribosome_binding_sites[frame] = 0;
            codon_types[frame] = 0;
        } else {
            let rbs_pair = usize::from(node.rbs_pair);
            let max_ribosome_binding_site = usize::from(rbs_choices.site[rbs_pair]);
            let codon_type_index = usize::from(node.codon_type_index);
            counts.rbs_background[max_ribosome_binding_site] += 1.0;
            counts.rbs_background_sum += 1.0;

            // `cscore + wt*rbs_wt + wt*type_wt` with the two fused
            // multiply-adds of the reference build.
            let score = weight_factor.mul_add(
                training.start_type_weights[codon_type_index],
                weight_factor.mul_add(rbs_choices.weight[rbs_pair], node.coding_score),
            );

            if score >= best_scores[frame] {
                best_scores[frame] = score;
                best_node_indices[frame] = Some(node_index);
                codon_types[frame] = codon_type_index;
                ribosome_binding_sites[frame] = max_ribosome_binding_site;
            }
        }
    };

    match strand {
        Strand::Forward => {
            for (node_index, node) in nodes.iter().enumerate() {
                visit(node_index, node);
            }
        }
        _ => {
            for (node_index, node) in nodes.iter().enumerate().rev() {
                visit(node_index, node);
            }
        }
    }
    counts
}
