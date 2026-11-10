use crate::constants::{GENE_RATIO_THRESHOLD, INITIAL_SCORE_THRESHOLD, THRESHOLD_DIVISOR};
use crate::training::count_upstream_composition;
use crate::{
    node::{MotifSite, motif_site},
    training::{
        MAX_MOTIF_INDEX, MAX_TRAINING_ITERATIONS_NONSD, NUM_BASES, NUM_CODON_TYPES,
        NUM_MOTIF_SIZES, UPSTREAM_POSITIONS, WEIGHT_CLAMP_MAX, WEIGHT_CLAMP_MIN,
        build_coverage_map, common::normalize_and_log_transform,
        convert_upstream_composition_to_log_scores, count_motif, update_motif_weights,
    },
    types::{CodonType, Motif, MotifWeights, Node, Training},
};
use bio::bio_types::strand::Strand;

/// The node fields the iterations of [`train_starts_nonsd`] read, packed so
/// that every iteration scans the nodes in a fraction of their size.
#[derive(Debug, Clone, Copy)]
struct ScanNode {
    coding_score: f64,
    /// Score of the node's best upstream motif in the current iteration.
    motif_score: f64,
    /// Number of the node among the start nodes whose motifs are searched
    /// (`starts`); meaningless for the other nodes.
    start: u32,
    frame: u8,
    codon_type_index: u8,
    is_stop: bool,
    is_edge: bool,
    strand: Strand,
}

/// A start node whose upstream motifs are searched: every start node that
/// is not an edge.
struct SearchedStart {
    /// Index of the node.
    node: usize,
    /// What the motif searches read of the sequence, the same in every
    /// iteration.
    site: MotifSite,
    /// Position of the start codon on its own strand.
    strand_position: usize,
    is_forward: bool,
    /// The best upstream motif in the current iteration.
    motif: Motif,
}

/// The counts one iteration takes of the starts it selects as real.
struct RealCounts {
    motifs: Box<MotifWeights>,
    zero_motif: f64,
    types: [f64; NUM_CODON_TYPES],
    genes: f64,
}

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
    *training.upstream_composition = [[0.0; NUM_BASES]; UPSTREAM_POSITIONS];
    for motif_length in training.motif_weights.iter_mut() {
        for spacer in motif_length.iter_mut() {
            spacer.fill(0.0);
        }
    }
    training.no_motif_weight = 0.0;

    // The iterations read a few fields of every node and, of the sequence,
    // the same windows upstream of every start node. Both are gathered
    // once; the iterations then run on compact copies.
    let mut scan_nodes: Vec<ScanNode> = Vec::with_capacity(nodes.len());
    let mut starts: Vec<SearchedStart> = Vec::new();
    let mut total_starts = 0.0;
    for (node_index, node) in nodes.iter().enumerate() {
        let is_stop = node.position.codon_type == CodonType::Stop;
        if !is_stop {
            type_background[node.position.codon_type.to_index()] += 1.0;
            total_starts += 1.0;
        }
        let is_searched = !is_stop && !node.position.is_edge;
        scan_nodes.push(ScanNode {
            coding_score: node.scores.coding_score,
            motif_score: node.motif_info.best_motif.score,
            start: starts.len() as u32,
            frame: (node.position.index % NUM_CODON_TYPES) as u8,
            codon_type_index: node.position.codon_type.to_index() as u8,
            is_stop,
            is_edge: node.position.is_edge,
            strand: node.position.strand,
        });
        if is_searched {
            let is_forward = match node.position.strand {
                Strand::Forward => true,
                Strand::Reverse => false,
                Strand::Unknown => unreachable!(),
            };
            starts.push(SearchedStart {
                node: node_index,
                site: motif_site(sequence, reverse_sequence, sequence_length, node),
                strand_position: if is_forward {
                    node.position.index
                } else {
                    sequence_length - 1 - node.position.index
                },
                is_forward,
                motif: Motif::default(),
            });
        }
    }

    if total_starts > 0.0 {
        for type_background_item in type_background.iter_mut().take(NUM_CODON_TYPES) {
            *type_background_item /= total_starts;
        }
    }

    let strand_sequence = |is_forward: bool| {
        if is_forward {
            sequence
        } else {
            reverse_sequence
        }
    };

    let mut motif_good = Box::new([[[0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
    let mut motif_background =
        Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
    let mut zero_motif_background = 0.0;

    for iteration in 0..MAX_TRAINING_ITERATIONS_NONSD {
        let stage = if iteration < 4 {
            0
        } else if iteration < 12 {
            1
        } else {
            2
        };

        // The best motif of every start node under the current weights.
        // In the first stage the background counts every motif upstream of
        // every start node, whatever its best motif is, so the background
        // of the first iteration is that of the whole stage.
        let count_background = stage != 0 || iteration == 0;
        if count_background {
            for motif_length in motif_background.iter_mut() {
                for spacer in motif_length.iter_mut() {
                    spacer.fill(0.0);
                }
            }
            zero_motif_background = 0.0;
        }
        for start in starts.iter_mut() {
            start.motif = start.site.best_motif(training, stage);
            scan_nodes[start.node].motif_score = start.motif.score;
            if count_background {
                if stage == 0 && start.motif.length != 0 {
                    start.site.count_windows(&mut motif_background);
                } else {
                    count_motif(
                        &mut motif_background,
                        &mut zero_motif_background,
                        &start.motif,
                        strand_sequence(start.is_forward),
                        start.strand_position,
                        stage,
                    );
                }
            }
        }

        if count_background {
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
                            motif_background[motif_length_index][spacer_index][motif_index] /=
                                bg_sum;
                        }
                    }
                }
                zero_motif_background /= bg_sum;
            }
        }

        let mut real = RealCounts {
            motifs: Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]),
            zero_motif: 0.0,
            types: [0.0; NUM_CODON_TYPES],
            genes: 0.0,
        };
        let count_upstream = iteration == MAX_TRAINING_ITERATIONS_NONSD - 1;

        // Counts the best start of an open reading frame as a real start.
        let mut count_real = |best_node: usize, training: &mut Training| {
            let node = &scan_nodes[best_node];
            let start = &starts[node.start as usize];
            real.genes += 1.0;
            real.types[usize::from(node.codon_type_index)] += 1.0;
            count_motif(
                &mut real.motifs,
                &mut real.zero_motif,
                &start.motif,
                strand_sequence(start.is_forward),
                start.strand_position,
                stage,
            );
            if count_upstream {
                count_upstream_composition(
                    strand_sequence(start.is_forward),
                    sequence_length,
                    if start.is_forward {
                        Strand::Forward
                    } else {
                        Strand::Reverse
                    },
                    nodes[best_node].position.index,
                    training,
                );
            }
        };

        // Choose the best start per frame between stops: the forward strand
        // from the first node, then the reverse strand from the last.
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut best_scores = [0.0f64; NUM_CODON_TYPES];
            let mut best_indices: [Option<usize>; NUM_CODON_TYPES] = [None; NUM_CODON_TYPES];
            let mut visit = |node_index: usize, node: &ScanNode| {
                if !node.is_stop && node.is_edge {
                    return;
                }
                let frame = usize::from(node.frame);
                if node.is_stop && node.strand == strand {
                    if best_scores[frame] >= score_threshold
                        && let Some(best_node) = best_indices[frame]
                    {
                        count_real(best_node, training);
                    }
                    best_scores[frame] = 0.0;
                    best_indices[frame] = None;
                } else if node.strand == strand {
                    // `cscore + wt*mot.score + wt*type_wt` with the two fused
                    // multiply-adds of the reference build.
                    let score = weight_factor.mul_add(
                        training.start_type_weights[usize::from(node.codon_type_index)],
                        weight_factor.mul_add(node.motif_score, node.coding_score),
                    );
                    if score >= best_scores[frame] {
                        best_scores[frame] = score;
                        best_indices[frame] = Some(node_index);
                    }
                }
            };
            if strand == Strand::Forward {
                for (node_index, node) in scan_nodes.iter().enumerate() {
                    visit(node_index, node);
                }
            } else {
                for (node_index, node) in scan_nodes.iter().enumerate().rev() {
                    visit(node_index, node);
                }
            }
        }

        if stage < 2 {
            build_coverage_map(&real.motifs, &mut motif_good, real.genes, stage);
        }

        update_motif_weights(
            &real.motifs,
            &motif_background,
            real.zero_motif,
            zero_motif_background,
            &motif_good,
            stage,
            training,
        );

        // Update type weights using common function
        normalize_and_log_transform(
            &mut real.types,
            &type_background,
            WEIGHT_CLAMP_MIN,
            WEIGHT_CLAMP_MAX,
        );
        training.start_type_weights.copy_from_slice(&real.types);

        // Lower threshold if not enough genes found
        if real.genes <= nodes.len() as f64 / GENE_RATIO_THRESHOLD {
            score_threshold /= THRESHOLD_DIVISOR;
        }
    }

    // The nodes keep the motifs of the last iteration.
    for start in starts {
        nodes[start.node].motif_info.best_motif = start.motif;
    }

    convert_upstream_composition_to_log_scores(training);
}
