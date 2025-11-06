use bio::bio_types::strand::Strand;
use rayon::prelude::*;

use crate::{
    training::{
        GC_HIGH_AT_FREQ, GC_HIGH_GC_FREQ, GC_LOW_AT_FREQ, GC_LOW_GC_FREQ, GENE_RATIO_THRESHOLD,
        INITIAL_SCORE_THRESHOLD, MAX_GC_CONTENT, MAX_TRAINING_ITERATIONS_SD, MIN_GC_CONTENT,
        NUM_BASES, NUM_CODON_TYPES, NUM_RBS_WEIGHTS, RBS_WEIGHT_THRESHOLD_HIGH, THRESHOLD_DIVISOR,
        UPSTREAM_POSITIONS, WEIGHT_CLAMP_MAX, WEIGHT_CLAMP_MIN, count_upstream_composition,
    },
    types::{CodonType, Node, Training},
};

/// Train start codon recognition using Shine-Dalgarno motifs
pub fn train_starts_sd(
    sequence: &[u8],
    reverse_sequence: &[u8],
    sequence_length: usize,
    nodes: &[Node],
    training: &mut Training,
) {
    let weight_factor = training.start_weight_factor;
    let mut type_background = [0.0; NUM_CODON_TYPES];
    let mut ribosome_binding_site_background = [0.0; NUM_RBS_WEIGHTS];

    training.start_type_weights = [0.0; NUM_CODON_TYPES];
    training.rbs_weights = Box::new([0.0; NUM_RBS_WEIGHTS]);
    training.upstream_composition = Box::new([[0.0; NUM_BASES]; UPSTREAM_POSITIONS]);

    let mut total_sum = 0.0;
    for node in nodes {
        if node.position.codon_type == CodonType::Stop {
            continue;
        }
        type_background[node.position.codon_type.to_index()] += 1.0;
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
            if node.position.codon_type == CodonType::Stop || node.position.is_edge {
                continue;
            }

            let max_ribosome_binding_site = if training.rbs_weights
                [node.motif_info.ribosome_binding_sites[0]]
                > training.rbs_weights[node.motif_info.ribosome_binding_sites[1]]
                    + RBS_WEIGHT_THRESHOLD_HIGH
                || node.motif_info.ribosome_binding_sites[1] == 0
            {
                node.motif_info.ribosome_binding_sites[0]
            } else if training.rbs_weights[node.motif_info.ribosome_binding_sites[0]]
                < training.rbs_weights[node.motif_info.ribosome_binding_sites[1]]
                    - RBS_WEIGHT_THRESHOLD_HIGH
                || node.motif_info.ribosome_binding_sites[0] == 0
            {
                node.motif_info.ribosome_binding_sites[1]
            } else {
                node.motif_info.ribosome_binding_sites[0]
                    .max(node.motif_info.ribosome_binding_sites[1])
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

        let mut ribosome_binding_site_real = [0.0; NUM_RBS_WEIGHTS];
        let mut type_real = [0.0; NUM_CODON_TYPES];

        let mut best_scores = [0.0; NUM_CODON_TYPES];
        let mut best_node_indices: [Option<usize>; NUM_CODON_TYPES] = [None; NUM_CODON_TYPES];
        let mut ribosome_binding_sites = [0; NUM_CODON_TYPES];
        let mut codon_types = [0; NUM_CODON_TYPES];

        for (node_index, node) in nodes.iter().enumerate() {
            if node.position.codon_type != CodonType::Stop && node.position.is_edge {
                continue;
            }

            let frame = node.position.index % NUM_CODON_TYPES;

            if node.position.codon_type == CodonType::Stop
                && node.position.strand == Strand::Forward
            {
                if best_scores[frame] >= score_threshold && best_node_indices[frame].is_some() {
                    let best_index = best_node_indices[frame].unwrap();
                    if best_index < nodes.len()
                        && nodes[best_index].position.index % NUM_CODON_TYPES == frame
                    {
                        ribosome_binding_site_real[ribosome_binding_sites[frame]] += 1.0;
                        type_real[codon_types[frame]] += 1.0;

                        if iteration == MAX_TRAINING_ITERATIONS_SD - 1 {
                            count_upstream_composition(
                                sequence,
                                sequence_length,
                                Strand::Forward,
                                nodes[best_index].position.index,
                                training,
                            );
                        }
                    }
                }

                best_scores[frame] = 0.0;
                best_node_indices[frame] = None;
                ribosome_binding_sites[frame] = 0;
                codon_types[frame] = 0;
            } else if node.position.strand == Strand::Forward {
                let max_ribosome_binding_site = if training.rbs_weights
                    [node.motif_info.ribosome_binding_sites[0]]
                    > training.rbs_weights[node.motif_info.ribosome_binding_sites[1]]
                        + RBS_WEIGHT_THRESHOLD_HIGH
                    || node.motif_info.ribosome_binding_sites[1] == 0
                {
                    node.motif_info.ribosome_binding_sites[0]
                } else if training.rbs_weights[node.motif_info.ribosome_binding_sites[0]]
                    < training.rbs_weights[node.motif_info.ribosome_binding_sites[1]]
                        - RBS_WEIGHT_THRESHOLD_HIGH
                    || node.motif_info.ribosome_binding_sites[0] == 0
                {
                    node.motif_info.ribosome_binding_sites[1]
                } else {
                    node.motif_info.ribosome_binding_sites[0]
                        .max(node.motif_info.ribosome_binding_sites[1])
                };

                let score = node.scores.coding_score
                    + weight_factor * training.rbs_weights[max_ribosome_binding_site]
                    + weight_factor
                        * training.start_type_weights[node.position.codon_type.to_index()];

                if score >= best_scores[frame] {
                    best_scores[frame] = score;
                    best_node_indices[frame] = Some(node_index);
                    codon_types[frame] = node.position.codon_type.to_index();
                    ribosome_binding_sites[frame] = max_ribosome_binding_site;
                }
            }
        }

        best_scores.fill(0.0);
        best_node_indices.fill(None);
        ribosome_binding_sites.fill(0);
        codon_types.fill(0);

        for (node_index, node) in nodes.iter().enumerate().rev() {
            if node.position.codon_type != CodonType::Stop && node.position.is_edge {
                continue;
            }

            let frame = node.position.index % NUM_CODON_TYPES;

            if node.position.codon_type == CodonType::Stop
                && node.position.strand == Strand::Reverse
            {
                if best_scores[frame] >= score_threshold && best_node_indices[frame].is_some() {
                    let best_index = best_node_indices[frame].unwrap();
                    if best_index < nodes.len()
                        && nodes[best_index].position.index % NUM_CODON_TYPES == frame
                    {
                        ribosome_binding_site_real[ribosome_binding_sites[frame]] += 1.0;
                        type_real[codon_types[frame]] += 1.0;

                        if iteration == MAX_TRAINING_ITERATIONS_SD - 1 {
                            count_upstream_composition(
                                reverse_sequence,
                                sequence_length,
                                Strand::Reverse,
                                nodes[best_index].position.index,
                                training,
                            );
                        }
                    }
                }

                best_scores[frame] = 0.0;
                best_node_indices[frame] = None;
                ribosome_binding_sites[frame] = 0;
                codon_types[frame] = 0;
            } else if node.position.strand == Strand::Reverse {
                let max_ribosome_binding_site = if training.rbs_weights
                    [node.motif_info.ribosome_binding_sites[0]]
                    > training.rbs_weights[node.motif_info.ribosome_binding_sites[1]]
                        + RBS_WEIGHT_THRESHOLD_HIGH
                    || node.motif_info.ribosome_binding_sites[1] == 0
                {
                    node.motif_info.ribosome_binding_sites[0]
                } else if training.rbs_weights[node.motif_info.ribosome_binding_sites[0]]
                    < training.rbs_weights[node.motif_info.ribosome_binding_sites[1]]
                        - RBS_WEIGHT_THRESHOLD_HIGH
                    || node.motif_info.ribosome_binding_sites[0] == 0
                {
                    node.motif_info.ribosome_binding_sites[1]
                } else {
                    node.motif_info.ribosome_binding_sites[0]
                        .max(node.motif_info.ribosome_binding_sites[1])
                };

                let score = node.scores.coding_score
                    + weight_factor * training.rbs_weights[max_ribosome_binding_site]
                    + weight_factor
                        * training.start_type_weights[node.position.codon_type.to_index()];

                if score >= best_scores[frame] {
                    best_scores[frame] = score;
                    best_node_indices[frame] = Some(node_index);
                    codon_types[frame] = node.position.codon_type.to_index();
                    ribosome_binding_sites[frame] = max_ribosome_binding_site;
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
