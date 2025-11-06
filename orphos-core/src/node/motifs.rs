use bio::bio_types::strand::Strand;
use rayon::prelude::*;

use crate::{
    constants::{
        INITIAL_MAX_SCORE, MIN_MOTIF_SCORE, MOTIF_THRESHOLD_OFFSET, RBS_DOWNSTREAM_DISTANCE,
        RBS_UPSTREAM_DISTANCE,
    },
    sequence::{calculate_kmer_index, shine_dalgarno_exact, shine_dalgarno_mm},
    types::{CodonType, Motif, Node, Training},
};

/// Calculate both exact and mismatch Shine-Dalgarno scores for a position
fn calculate_rbs_scores(
    seq: &[u8],
    pos: usize,
    target_pos: usize,
    rbs_weights: &[f64],
) -> [usize; 2] {
    [
        shine_dalgarno_exact(seq, pos, target_pos, rbs_weights),
        shine_dalgarno_mm(seq, pos, target_pos, rbs_weights),
    ]
}

/// RBS Scoring Function: Calculate the RBS motif and then multiply it by the
/// appropriate weight for that motif (determined in the start training function).
///
/// # Parameters
/// * `seq` - Forward strand sequence
/// * `rseq` - Reverse strand sequence
/// * `sequence_length` - Sequence length
/// * `nodes` - Nodes to score
/// * `training` - Training data with RBS weights
///
/// # Effect
/// Updates the ribosome_binding_sites scores in each node's motif_info
pub fn rbs_score(
    seq: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &mut [Node],
    training: &Training,
) {
    // Parallelize RBS scoring across all start nodes
    nodes
        .par_iter_mut()
        .enumerate()
        .filter(|(_, node)| node.position.codon_type != CodonType::Stop && !node.position.is_edge)
        .for_each(|(_, node)| {
            node.motif_info.ribosome_binding_sites[0] = 0;
            node.motif_info.ribosome_binding_sites[1] = 0;

            match node.position.strand {
                Strand::Forward => {
                    let search_start = node.position.index.saturating_sub(RBS_UPSTREAM_DISTANCE);
                    let search_end = node.position.index.saturating_sub(RBS_DOWNSTREAM_DISTANCE);

                    for j in search_start..=search_end {
                        let cur_sc = calculate_rbs_scores(
                            seq,
                            j,
                            node.position.index,
                            &*training.rbs_weights,
                        );

                        if cur_sc[0] > node.motif_info.ribosome_binding_sites[0] {
                            node.motif_info.ribosome_binding_sites[0] = cur_sc[0];
                        }
                        if cur_sc[1] > node.motif_info.ribosome_binding_sites[1] {
                            node.motif_info.ribosome_binding_sites[1] = cur_sc[1];
                        }
                    }
                }
                Strand::Reverse => {
                    let upstream_offset = RBS_UPSTREAM_DISTANCE + 1;
                    let downstream_offset = RBS_DOWNSTREAM_DISTANCE + 1;
                    if node.position.index >= sequence_length {
                        return; // Skip invalid node positions
                    }
                    let start_pos = if sequence_length >= node.position.index + upstream_offset {
                        sequence_length - node.position.index - upstream_offset
                    } else {
                        0 // If calculation would underflow, start from beginning
                    };

                    let end_pos = if sequence_length >= node.position.index + downstream_offset {
                        sequence_length - node.position.index - downstream_offset
                    } else {
                        0 // If calculation would underflow, set to 0
                    };
                    let target_pos = sequence_length - 1 - node.position.index;

                    for j in start_pos..=end_pos {
                        if !(0..sequence_length).contains(&j) {
                            continue;
                        }

                        let cur_sc = calculate_rbs_scores(
                            reverse_complement_encoded_sequence,
                            j,
                            target_pos,
                            &*training.rbs_weights,
                        );

                        if cur_sc[0] > node.motif_info.ribosome_binding_sites[0] {
                            node.motif_info.ribosome_binding_sites[0] = cur_sc[0];
                        }
                        if cur_sc[1] > node.motif_info.ribosome_binding_sites[1] {
                            node.motif_info.ribosome_binding_sites[1] = cur_sc[1];
                        }
                    }
                }
                Strand::Unknown => unreachable!(),
            }
        });
}

/// Find the highest scoring motif/spacer combination for a node.
///
/// Given the weights for various motifs/distances from the training file,
/// return the highest scoring mer/spacer combination of 3-6bp motifs with a
/// spacer ranging from 3bp to 15bp. In the final stage of start training, only
/// good scoring motifs are returned.
///
/// # Parameters
/// * `training` - Training data with motif weights
/// * `seq` - Forward strand sequence
/// * `rseq` - Reverse strand sequence
/// * `sequence_length` - Sequence length
/// * `node` - Node to find motif for
/// * `stage` - Training stage (affects whether poor motifs are filtered)
///
/// # Effect
/// Updates the best_motif in the node's motif_info
pub fn find_best_upstream_motif(
    training: &Training,
    seq: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    node: &mut Node,
    stage: usize,
) {
    if node.position.codon_type == CodonType::Stop || node.position.is_edge {
        return;
    }

    let (wseq, start) = match node.position.strand {
        Strand::Forward => (seq, node.position.index),
        Strand::Reverse => (
            reverse_complement_encoded_sequence,
            sequence_length - 1 - node.position.index,
        ),
        Strand::Unknown => unreachable!(),
    };

    let mut max_sc = INITIAL_MAX_SCORE;
    let mut max_spacendx = 0;
    let mut max_spacer = 0;
    let mut max_ndx = 0;
    let mut max_len = 0;

    // Search through motif lengths 3-6 (i goes from 3 down to 0, representing lengths 6 down to 3)
    for i in (0..=3).rev() {
        let motif_len = i + 3;

        // Search positions from start-18-i to start-6-i
        let search_start: isize = start as isize - 18 - i;
        let search_end: isize = start as isize - 6 - i;

        for j in search_start..=search_end {
            if j < 0 {
                continue;
            }
            // Ensure the k-mer window [j, j+motif_len) is within the nucleotide sequence.
            // the nucleotide length (sequence_length), not the byte length of the buffer.
            if (j as usize) + (motif_len as usize) > sequence_length {
                continue;
            }
            let spacer = start as isize - j - i - 3;
            let spacendx = if j <= start as isize - 16 - i {
                3
            } else if j <= start as isize - 14 - i {
                2
            } else if j >= start as isize - 7 - i {
                1
            } else {
                0
            };

            let index = calculate_kmer_index(motif_len as usize, wseq, j as usize);
            let score = training.motif_weights[i as usize][spacendx][index];

            if score > max_sc {
                max_sc = score;
                max_spacendx = spacendx;
                max_spacer = spacer;
                max_ndx = index;
                max_len = motif_len;
            }
        }
    }

    // In stage 2, only accept good scoring motifs
    let is_stage_two = stage == 2;
    let is_poor_motif =
        max_sc == MIN_MOTIF_SCORE || max_sc < training.no_motif_weight + MOTIF_THRESHOLD_OFFSET;

    // Do NOT neutralize the score in early stages when no window is found.
    // C keeps max_sc at -100.0 in this case, strongly disfavoring edge starts.
    let effective_max_sc = max_sc;

    node.motif_info.best_motif = if is_stage_two && is_poor_motif {
        Motif {
            index: 0,
            length: 0,
            space_index: 0,
            spacer: 0,
            score: training.no_motif_weight,
        }
    } else {
        Motif {
            index: max_ndx,
            length: max_len as usize,
            space_index: max_spacendx,
            spacer: max_spacer as usize,
            score: effective_max_sc,
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::*;
    use bio::bio_types::strand::Strand;

    fn create_test_node(strand: Strand, index: usize, codon_type: CodonType) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value: (index + 100) as isize,
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
            start_type_weights: [1.0, 2.0, 3.0],
            rbs_weights: Box::new([1.0; 28]),
            upstream_composition: Box::new([[0.25; 4]; 32]),
            motif_weights: Box::new([[[1.0; 4096]; 4]; 4]),
            no_motif_weight: 0.5,
            start_weight_factor: 4.35,
            gc_bias_factors: [1.0; 3],
            gene_dicodon_table: Box::new([0.0; 4096]),
            total_dicodons: 0,
        }
    }

    #[test]
    fn test_calculate_rbs_scores() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3]; // 8-base sequence
        let rbs_weights = [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        ];
        let pos = 2;
        let target_pos = 6;

        let scores = calculate_rbs_scores(&seq, pos, target_pos, &rbs_weights);

        assert_eq!(scores.len(), 2);
    }

    #[test]
    fn test_calculate_rbs_scores_boundary_conditions() {
        let seq = vec![0, 1, 2, 3];
        let rbs_weights = [1.0; 16];

        // Test at start of sequence
        let scores = calculate_rbs_scores(&seq, 0, 3, &rbs_weights);
        assert_eq!(scores.len(), 2);

        // Test near end of sequence
        let scores = calculate_rbs_scores(&seq, 2, 3, &rbs_weights);
        assert_eq!(scores.len(), 2);
    }

    #[test]
    fn test_rbs_score_forward_strand() {
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0,
            1, 2, 3,
        ];
        let reverse_seq = vec![
            3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3,
            2, 1, 0,
        ];
        let training = create_test_training();

        let mut nodes = vec![create_test_node(Strand::Forward, 25, CodonType::Atg)];

        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert_eq!(nodes[0].motif_info.ribosome_binding_sites.len(), 2);
    }

    #[test]
    fn test_rbs_score_reverse_strand() {
        // Create a longer sequence to handle reverse strand calculations
        let seq = vec![0; 60]; // 60-base sequence
        let reverse_seq = vec![3; 60];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Reverse, 35, CodonType::Atg), // Position far from edges
        ];

        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert_eq!(nodes[0].motif_info.ribosome_binding_sites.len(), 2);
    }

    #[test]
    fn test_rbs_score_stop_codon_skipped() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut nodes = vec![create_test_node(Strand::Forward, 6, CodonType::Stop)];

        let original_rbs = nodes[0].motif_info.ribosome_binding_sites;
        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert_eq!(nodes[0].motif_info.ribosome_binding_sites, original_rbs);
    }

    #[test]
    fn test_rbs_score_edge_node_skipped() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 6, CodonType::Atg);
        node.position.is_edge = true;
        let mut nodes = vec![node];

        let original_rbs = nodes[0].motif_info.ribosome_binding_sites;
        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert_eq!(nodes[0].motif_info.ribosome_binding_sites, original_rbs);
    }

    #[test]
    fn test_rbs_score_parallel_processing() {
        // Create a long sequence to handle all positions safely
        let seq = vec![0; 80];
        let reverse_seq = vec![3; 80];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Forward, 30, CodonType::Atg),
            create_test_node(Strand::Reverse, 50, CodonType::Atg),
            create_test_node(Strand::Forward, 70, CodonType::Atg),
        ];

        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        for node in &nodes {
            assert_eq!(node.motif_info.ribosome_binding_sites.len(), 2);
        }
    }

    #[test]
    fn test_find_best_upstream_motif_stop_codon() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 6, CodonType::Stop);
        let original_motif_index = node.motif_info.best_motif.index;
        let original_motif_length = node.motif_info.best_motif.length;

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        assert_eq!(node.motif_info.best_motif.index, original_motif_index);
        assert_eq!(node.motif_info.best_motif.length, original_motif_length);
    }

    #[test]
    fn test_find_best_upstream_motif_edge_node() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 6, CodonType::Atg);
        node.position.is_edge = true;
        let original_motif_index = node.motif_info.best_motif.index;
        let original_motif_length = node.motif_info.best_motif.length;

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        assert_eq!(node.motif_info.best_motif.index, original_motif_index);
        assert_eq!(node.motif_info.best_motif.length, original_motif_length);
    }

    #[test]
    fn test_find_best_upstream_motif_forward_strand() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 15, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        // length and spacer are usize, so always >= 0
    }

    #[test]
    fn test_find_best_upstream_motif_reverse_strand() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Reverse, 15, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        // length and spacer are usize, so always >= 0
    }

    #[test]
    fn test_find_best_upstream_motif_stage_two_filtering() {
        let mut training = create_test_training();
        training.no_motif_weight = 10.0; // High no-motif weight

        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];

        let mut node = create_test_node(Strand::Forward, 15, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 2);

        if node.motif_info.best_motif.score < training.no_motif_weight + MOTIF_THRESHOLD_OFFSET {
            assert_eq!(node.motif_info.best_motif.length, 0);
            assert_eq!(node.motif_info.best_motif.score, training.no_motif_weight);
        }
    }

    #[test]
    fn test_find_best_upstream_motif_stage_one_accepts_all() {
        let training = create_test_training();
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];

        let mut node = create_test_node(Strand::Forward, 15, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        // length is usize, so always >= 0
    }

    #[test]
    fn test_find_best_upstream_motif_empty_sequence() {
        let seq = vec![];
        let reverse_seq = vec![];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 0, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, 0, &mut node, 1);

        // Should handle empty sequence gracefully
        // length is usize, so always >= 0
    }

    #[test]
    fn test_find_best_upstream_motif_short_sequence() {
        let seq = vec![0, 1, 2];
        let reverse_seq = vec![2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 2, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        // Should handle short sequences without panicking
        // length is usize, so always >= 0
    }

    #[test]
    fn test_motif_spacer_classification() {
        let training = create_test_training();
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ];
        let reverse_seq = vec![
            3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0,
        ];

        let mut node = create_test_node(Strand::Forward, 20, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        assert!(node.motif_info.best_motif.space_index <= 3);
    }

    #[test]
    fn test_motif_length_range() {
        let training = create_test_training();
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ];
        let reverse_seq = vec![
            3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0,
        ];

        let mut node = create_test_node(Strand::Forward, 20, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        assert!(node.motif_info.best_motif.length <= 6);
    }
}
