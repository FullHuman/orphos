use bio::bio_types::strand::Strand;
use rayon::prelude::*;

use crate::{
    training::{
        GC_HIGH_AT_FREQ, GC_HIGH_GC_FREQ, GC_LOW_AT_FREQ, GC_LOW_GC_FREQ, GENE_RATIO_THRESHOLD,
        INITIAL_SCORE_THRESHOLD, MAX_GC_CONTENT, MAX_TRAINING_ITERATIONS_SD, MIN_GC_CONTENT,
        NUM_BASES, NUM_CODON_TYPES, NUM_RBS_WEIGHTS, RBS_WEIGHT_THRESHOLD_HIGH, THRESHOLD_DIVISOR,
        UPSTREAM_POSITIONS, WEIGHT_CLAMP_MAX, WEIGHT_CLAMP_MIN, count_upstream_composition_into,
    },
    types::{CodonType, Node, Training},
};

/// Train start codon recognition using Shine-Dalgarno motifs
/// The node fields read by [`train_starts_sd`].
#[derive(Debug, Clone, Copy)]
struct SdNode {
    index: usize,
    strand: Strand,
    coding_score: f64,
    ribosome_binding_sites: [usize; 2],
    codon_type_index: usize,
    is_stop: bool,
    is_edge: bool,
}

impl From<&Node> for SdNode {
    fn from(node: &Node) -> Self {
        Self {
            index: node.position.index,
            strand: node.position.strand,
            coding_score: node.scores.coding_score,
            ribosome_binding_sites: node.motif_info.ribosome_binding_sites,
            codon_type_index: node.position.codon_type.to_index(),
            is_stop: node.position.codon_type == CodonType::Stop,
            is_edge: node.position.is_edge,
        }
    }
}

pub fn train_starts_sd(
    sequence: &[u8],
    reverse_sequence: &[u8],
    sequence_length: usize,
    nodes: &[Node],
    training: &mut Training,
) {
    // The iterations below only read a few fields of every node; a compact
    // copy keeps the repeated scans cache friendly.
    let compact: Vec<SdNode> = nodes.iter().map(SdNode::from).collect();
    let nodes = &compact;
    let weight_factor = training.start_weight_factor;
    let mut type_background = [0.0; NUM_CODON_TYPES];
    let mut ribosome_binding_site_background = [0.0; NUM_RBS_WEIGHTS];

    training.start_type_weights = [0.0; NUM_CODON_TYPES];
    *training.rbs_weights = [0.0; NUM_RBS_WEIGHTS];
    *training.upstream_composition = [[0.0; NUM_BASES]; UPSTREAM_POSITIONS];

    let mut total_sum = 0.0;
    for node in nodes {
        if node.is_stop {
            continue;
        }
        type_background[node.codon_type_index] += 1.0;
        total_sum += 1.0;
    }

    if total_sum > 0.0 {
        for weight in type_background.iter_mut().take(NUM_CODON_TYPES) {
            *weight /= total_sum;
        }
    }

    let mut score_threshold = INITIAL_SCORE_THRESHOLD;

    for iteration in 0..MAX_TRAINING_ITERATIONS_SD {
        ribosome_binding_site_background.fill(0.0);
        let mut ribosome_binding_site_background_sum = 0.0;

        for node in nodes {
            if node.is_stop || node.is_edge {
                continue;
            }

            let max_ribosome_binding_site = if training.rbs_weights[node.ribosome_binding_sites[0]]
                > training.rbs_weights[node.ribosome_binding_sites[1]] + RBS_WEIGHT_THRESHOLD_HIGH
                || node.ribosome_binding_sites[1] == 0
            {
                node.ribosome_binding_sites[0]
            } else if training.rbs_weights[node.ribosome_binding_sites[0]]
                < training.rbs_weights[node.ribosome_binding_sites[1]] - RBS_WEIGHT_THRESHOLD_HIGH
                || node.ribosome_binding_sites[0] == 0
            {
                node.ribosome_binding_sites[1]
            } else {
                node.ribosome_binding_sites[0].max(node.ribosome_binding_sites[1])
            };

            ribosome_binding_site_background[max_ribosome_binding_site] += 1.0;
            ribosome_binding_site_background_sum += 1.0;
        }

        if ribosome_binding_site_background_sum > 0.0 {
            for weight in ribosome_binding_site_background
                .iter_mut()
                .take(NUM_RBS_WEIGHTS)
            {
                *weight /= ribosome_binding_site_background_sum;
            }
        }

        // The two strands are scanned on two threads. Every count is a whole
        // number, so adding the strands' counts gives exactly what one scan
        // over both strands gave.
        let count_upstream = iteration == MAX_TRAINING_ITERATIONS_SD - 1;
        let (forward, reverse) = rayon::join(
            || {
                scan_strand(
                    nodes,
                    Strand::Forward,
                    sequence,
                    sequence_length,
                    training,
                    weight_factor,
                    score_threshold,
                    count_upstream,
                )
            },
            || {
                scan_strand(
                    nodes,
                    Strand::Reverse,
                    reverse_sequence,
                    sequence_length,
                    training,
                    weight_factor,
                    score_threshold,
                    count_upstream,
                )
            },
        );
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
                .par_iter_mut()
                .for_each(|weight| *weight = 0.0);
        } else {
            ribosome_binding_site_real
                .par_iter_mut()
                .for_each(|val| *val /= total_real_ribosome_binding_sites);

            training
                .rbs_weights
                .par_iter_mut()
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
                .par_iter_mut()
                .for_each(|weight| *weight = 0.0);
        } else {
            type_real
                .par_iter_mut()
                .for_each(|val| *val /= total_real_types);

            training
                .start_type_weights
                .par_iter_mut()
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

        if total_real_types <= nodes.len() as f64 / GENE_RATIO_THRESHOLD {
            score_threshold /= THRESHOLD_DIVISOR;
        }
    }

    training
        .upstream_composition
        .par_iter_mut()
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
    rbs_real: [f64; NUM_RBS_WEIGHTS],
    type_real: [f64; NUM_CODON_TYPES],
    upstream_composition: [[f64; NUM_BASES]; UPSTREAM_POSITIONS],
}

/// One iteration's scan of `strand`: for every frame, the best start in the
/// ORF ending at each stop is counted, with its RBS motif and start codon
/// type, when its score reaches `score_threshold`; with `count_upstream`,
/// the bases upstream of those starts are counted too. The forward strand
/// is scanned from the first node, the reverse strand from the last.
#[allow(clippy::too_many_arguments)]
fn scan_strand(
    nodes: &[SdNode],
    strand: Strand,
    sequence: &[u8],
    sequence_length: usize,
    training: &Training,
    weight_factor: f64,
    score_threshold: f64,
    count_upstream: bool,
) -> StrandCounts {
    let mut counts = StrandCounts {
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

        let frame = node.index % NUM_CODON_TYPES;

        if node.is_stop && node.strand == strand {
            if best_scores[frame] >= score_threshold && best_node_indices[frame].is_some() {
                let best_index = best_node_indices[frame].unwrap();
                if best_index < nodes.len() && nodes[best_index].index % NUM_CODON_TYPES == frame {
                    counts.rbs_real[ribosome_binding_sites[frame]] += 1.0;
                    counts.type_real[codon_types[frame]] += 1.0;

                    if count_upstream {
                        count_upstream_composition_into(
                            sequence,
                            sequence_length,
                            strand,
                            nodes[best_index].index,
                            &mut counts.upstream_composition,
                        );
                    }
                }
            }

            best_scores[frame] = 0.0;
            best_node_indices[frame] = None;
            ribosome_binding_sites[frame] = 0;
            codon_types[frame] = 0;
        } else if node.strand == strand {
            let max_ribosome_binding_site = if training.rbs_weights[node.ribosome_binding_sites[0]]
                > training.rbs_weights[node.ribosome_binding_sites[1]] + RBS_WEIGHT_THRESHOLD_HIGH
                || node.ribosome_binding_sites[1] == 0
            {
                node.ribosome_binding_sites[0]
            } else if training.rbs_weights[node.ribosome_binding_sites[0]]
                < training.rbs_weights[node.ribosome_binding_sites[1]] - RBS_WEIGHT_THRESHOLD_HIGH
                || node.ribosome_binding_sites[0] == 0
            {
                node.ribosome_binding_sites[1]
            } else {
                node.ribosome_binding_sites[0].max(node.ribosome_binding_sites[1])
            };

            // `cscore + wt*rbs_wt + wt*type_wt` with the two fused
            // multiply-adds of the reference build.
            let score = weight_factor.mul_add(
                training.start_type_weights[node.codon_type_index],
                weight_factor.mul_add(
                    training.rbs_weights[max_ribosome_binding_site],
                    node.coding_score,
                ),
            );

            if score >= best_scores[frame] {
                best_scores[frame] = score;
                best_node_indices[frame] = Some(node_index);
                codon_types[frame] = node.codon_type_index;
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
