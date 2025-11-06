use crate::constants::{GENE_RATIO_THRESHOLD, INITIAL_SCORE_THRESHOLD, THRESHOLD_DIVISOR};
use crate::training::count_upstream_composition;
use crate::{
    node::find_best_upstream_motif,
    training::{
        MAX_MOTIF_INDEX, MAX_TRAINING_ITERATIONS_NONSD, NUM_BASES, NUM_CODON_TYPES,
        NUM_MOTIF_SIZES, UPSTREAM_POSITIONS, WEIGHT_CLAMP_MAX, WEIGHT_CLAMP_MIN,
        build_coverage_map, common::normalize_and_log_transform,
        convert_upstream_composition_to_log_scores, update_motif_counts, update_motif_weights,
    },
    types::{CodonType, Node, Training},
};
use bio::bio_types::strand::Strand;

/// Train start codon recognition without Shine-Dalgarno motifs
pub fn train_starts_nonsd(
    sequence: &[u8],
    reverse_sequence: &[u8],
    sequence_length: usize,
    nodes: &mut [Node],
    training: &mut Training,
) {
    let mut type_background = [0.0; NUM_CODON_TYPES];
    let mut score_threshold = INITIAL_SCORE_THRESHOLD;
    let weight_factor = training.start_weight_factor;

    training.start_type_weights = [0.0; NUM_CODON_TYPES];
    training.upstream_composition = Box::new([[0.0; NUM_BASES]; UPSTREAM_POSITIONS]);
    training.motif_weights = Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
    training.no_motif_weight = 0.0;

    let mut total_starts = 0.0;
    for node in nodes.iter() {
        if node.position.codon_type != CodonType::Stop {
            type_background[node.position.codon_type.to_index()] += 1.0;
            total_starts += 1.0;
        }
    }

    if total_starts > 0.0 {
        for type_background_item in type_background.iter_mut().take(NUM_CODON_TYPES) {
            *type_background_item /= total_starts;
        }
    }

    let mut motif_good = Box::new([[[0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);

    for iteration in 0..MAX_TRAINING_ITERATIONS_NONSD {
        let stage = if iteration < 4 {
            0
        } else if iteration < 12 {
            1
        } else {
            2
        };

        let mut motif_background =
            Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let mut zero_motif_background = 0.0;
        for node in nodes.iter_mut() {
            if node.position.codon_type == CodonType::Stop || node.position.is_edge {
                continue;
            }
            find_best_upstream_motif(
                training,
                sequence,
                reverse_sequence,
                sequence_length,
                node,
                stage,
            );
            update_motif_counts(
                &mut motif_background,
                &mut zero_motif_background,
                sequence,
                reverse_sequence,
                sequence_length,
                node,
                stage,
            );
        }

        let mut bg_sum = zero_motif_background;
        for motif_length_index in 0..NUM_MOTIF_SIZES {
            for spacer_index in 0..NUM_MOTIF_SIZES {
                for motif_index in 0..MAX_MOTIF_INDEX {
                    bg_sum += motif_background[motif_length_index][spacer_index][motif_index];
                }
            }
        }
        if bg_sum > 0.0 {
            for motif_length_index in 0..NUM_MOTIF_SIZES {
                for spacer_index in 0..NUM_MOTIF_SIZES {
                    for motif_index in 0..MAX_MOTIF_INDEX {
                        motif_background[motif_length_index][spacer_index][motif_index] /= bg_sum;
                    }
                }
            }
            zero_motif_background /= bg_sum;
        }

        for node in nodes.iter_mut() {
            if node.position.codon_type != CodonType::Stop && !node.position.is_edge {
                find_best_upstream_motif(
                    training,
                    sequence,
                    reverse_sequence,
                    sequence_length,
                    node,
                    stage,
                );
            }
        }

        let mut motif_real = Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let mut zero_motif_real = 0.0;
        let mut type_real = [0.0; NUM_CODON_TYPES];
        let mut number_of_genes = 0.0;
        let mut selected_nodes: Vec<usize> = Vec::new();

        // Diagnostics accumulators for score summaries
        let mut _seg_count: usize = 0;
        let mut _seg_best_sum: f64 = 0.0;
        let mut seg_best_max: f64 = f64::NEG_INFINITY;
        let mut _selected_count: usize = 0;
        let mut _selected_best_sum: f64 = 0.0;
        let mut cscore_max: f64 = f64::NEG_INFINITY;

        // Forward strand pass: choose best starts per frame between STOPs
        let mut best_scores = [0.0f64; NUM_CODON_TYPES];
        let mut best_indices: [Option<usize>; NUM_CODON_TYPES] = [None; NUM_CODON_TYPES];
        for (j, node) in nodes.iter().enumerate() {
            if node.position.codon_type != CodonType::Stop && node.position.is_edge {
                continue;
            }
            let frame = node.position.index % NUM_CODON_TYPES;
            if node.position.codon_type == CodonType::Stop
                && node.position.strand == Strand::Forward
            {
                // TAP forward pass diagnostics removed
                if best_scores[frame] >= score_threshold
                    && let Some(bi) = best_indices[frame]
                {
                    number_of_genes += 1.0;
                    type_real[nodes[bi].position.codon_type.to_index()] += 1.0;
                    update_motif_counts(
                        &mut motif_real,
                        &mut zero_motif_real,
                        sequence,
                        reverse_sequence,
                        sequence_length,
                        &nodes[bi],
                        stage,
                    );
                    selected_nodes.push(bi);
                    _selected_count += 1;
                    _selected_best_sum += best_scores[frame];
                    if iteration == MAX_TRAINING_ITERATIONS_NONSD - 1 {
                        count_upstream_composition(
                            sequence,
                            sequence_length,
                            Strand::Forward,
                            nodes[bi].position.index,
                            training,
                        );
                    }
                }

                _seg_count += 1;
                _seg_best_sum += best_scores[frame];
                if best_scores[frame] > seg_best_max {
                    seg_best_max = best_scores[frame];
                }
                best_scores[frame] = 0.0;
                best_indices[frame] = None;
            } else if node.position.strand == Strand::Forward {
                let score = node.scores.coding_score
                    + weight_factor * node.motif_info.best_motif.score
                    + weight_factor
                        * training.start_type_weights[node.position.codon_type.to_index()];
                if score >= best_scores[frame] {
                    best_scores[frame] = score;
                    best_indices[frame] = Some(j);
                }
                if node.scores.coding_score > cscore_max {
                    cscore_max = node.scores.coding_score;
                }
            }
        }

        best_scores.fill(0.0);
        best_indices.fill(None);
        for (rev_j, node) in nodes.iter().enumerate().rev() {
            if node.position.codon_type != CodonType::Stop && node.position.is_edge {
                continue;
            }
            let frame = node.position.index % NUM_CODON_TYPES;
            if node.position.codon_type == CodonType::Stop
                && node.position.strand == Strand::Reverse
            {
                // TAP reverse pass diagnostics removed
                if best_scores[frame] >= score_threshold
                    && let Some(bi) = best_indices[frame]
                {
                    number_of_genes += 1.0;
                    type_real[nodes[bi].position.codon_type.to_index()] += 1.0;
                    update_motif_counts(
                        &mut motif_real,
                        &mut zero_motif_real,
                        sequence,
                        reverse_sequence,
                        sequence_length,
                        &nodes[bi],
                        stage,
                    );
                    selected_nodes.push(bi);
                    _selected_count += 1;
                    _selected_best_sum += best_scores[frame];
                    if iteration == MAX_TRAINING_ITERATIONS_NONSD - 1 {
                        count_upstream_composition(
                            reverse_sequence,
                            sequence_length,
                            Strand::Reverse,
                            nodes[bi].position.index,
                            training,
                        );
                    }
                }

                _seg_count += 1;
                _seg_best_sum += best_scores[frame];
                if best_scores[frame] > seg_best_max {
                    seg_best_max = best_scores[frame];
                }
                best_scores[frame] = 0.0;
                best_indices[frame] = None;
            } else if node.position.strand == Strand::Reverse {
                let score = node.scores.coding_score
                    + weight_factor * node.motif_info.best_motif.score
                    + weight_factor
                        * training.start_type_weights[node.position.codon_type.to_index()];
                if score >= best_scores[frame] {
                    best_scores[frame] = score;
                    best_indices[frame] = Some(rev_j);
                }
                if node.scores.coding_score > cscore_max {
                    cscore_max = node.scores.coding_score;
                }
            }
        }

        // Iteration window count summary removed (debug only)

        if stage < 2 {
            build_coverage_map(&motif_real, &mut motif_good, number_of_genes, stage);
            {
                // Summarize coverage map and how many real counts fall into good vs bad motifs (debug only)
                let mut _good_cnt = 0usize;
                let mut _good1_cnt = 0usize;
                let mut _good2_cnt = 0usize;
                let mut _bad_real_sum = 0.0f64;
                let mut _good_real_sum = 0.0f64;
                for l in 0..NUM_MOTIF_SIZES {
                    for s in 0..NUM_MOTIF_SIZES {
                        for idx in 0..MAX_MOTIF_INDEX {
                            match motif_good[l][s][idx] {
                                0 => {
                                    _bad_real_sum += motif_real[l][s][idx];
                                }
                                1 => {
                                    _good_cnt += 1;
                                    _good1_cnt += 1;
                                    _good_real_sum += motif_real[l][s][idx];
                                }
                                2 => {
                                    _good_cnt += 1;
                                    _good2_cnt += 1;
                                    _good_real_sum += motif_real[l][s][idx];
                                }
                                _ => {}
                            }
                        }
                    }
                }
                // coverage_map_summary logging removed
            }
        }

        update_motif_weights(
            &motif_real,
            &motif_background,
            zero_motif_real,
            zero_motif_background,
            &motif_good,
            stage,
            training,
        );

        // Update type weights using common function
        normalize_and_log_transform(
            &mut type_real,
            &type_background,
            WEIGHT_CLAMP_MIN,
            WEIGHT_CLAMP_MAX,
        );
        training.start_type_weights.copy_from_slice(&type_real);

        // Lower threshold if not enough genes found
        if number_of_genes <= nodes.len() as f64 / GENE_RATIO_THRESHOLD {
            score_threshold /= THRESHOLD_DIVISOR;
        }
    }

    convert_upstream_composition_to_log_scores(training);
}
