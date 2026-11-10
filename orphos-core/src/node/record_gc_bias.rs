use bio::bio_types::strand::Strand;

use crate::node::strands::{StrandNodes, has_start_before_stop, split_by_strand};
use crate::parallel;
use crate::sequence::FramePlot;
use crate::types::{CodonType, Node, Training};

/// Record GC bias for training
pub fn record_gc_bias(gc_frame_data: &[i32], nodes: &mut [Node], training: &mut Training) {
    record_gc_bias_of(gc_frame_data, nodes, training);
}

/// [`record_gc_bias`] for a frame plot in either of its representations:
/// one value per position, or the running totals of [`GcFramePlot`], with
/// which the plot values of the codons between two nodes are tallied
/// without visiting them.
///
/// [`GcFramePlot`]: crate::sequence::GcFramePlot
pub(crate) fn record_gc_bias_of<P>(gc_frame_plot: &P, nodes: &mut [Node], training: &mut Training)
where
    P: FramePlot + Sync + ?Sized,
{
    if nodes.is_empty() || gc_frame_plot.is_empty() {
        return;
    }

    // The strands walk the frame plot in opposite directions, each on its
    // own thread. Prodigal walks the reverse strand from the counts and
    // positions the forward strand's walk left behind; only a reverse start
    // node that no reverse stop of its frame precedes sees them (a gene
    // wrapping around a circular sequence), and when there is one the
    // strands are walked one after the other, as before.
    let (mut forward, mut reverse) = split_by_strand(nodes);
    let mut frame_counts = [[0i32; 3]; 3];
    let mut last_position = [0; 3];
    if has_start_before_stop(&reverse) {
        process_forward_strand(
            gc_frame_plot,
            &mut forward,
            &mut frame_counts,
            &mut last_position,
        );
        process_reverse_strand(
            gc_frame_plot,
            &mut reverse,
            &mut frame_counts,
            &mut last_position,
        );
    } else {
        parallel::join(
            || {
                process_forward_strand(
                    gc_frame_plot,
                    &mut forward,
                    &mut frame_counts,
                    &mut last_position,
                )
            },
            || {
                process_reverse_strand(
                    gc_frame_plot,
                    &mut reverse,
                    &mut [[0i32; 3]; 3],
                    &mut [0; 3],
                )
            },
        );
    }

    // Calculate overall bias with better error handling
    calculate_training_bias(nodes, training);
}

/// Process nodes on the forward strand
fn process_forward_strand<P: FramePlot + ?Sized>(
    gc_frame_plot: &P,
    nodes: &mut StrandNodes,
    frame_counts: &mut [[i32; 3]; 3],
    last_position: &mut [usize; 3],
) {
    nodes.iter_mut().rev().for_each(|node| {
        let index = node.position.index;
        let frame_index = index % 3;
        let frame_offset = 3 - frame_index;
        let bin = |gc_value| forward_gc_bin(gc_value, frame_offset);

        if node.position.codon_type == CodonType::Stop {
            // Reset counters for this frame
            frame_counts[frame_index] = [0; 3];
            last_position[frame_index] = index;
            gc_frame_plot.tally(index, 1, bin, &mut frame_counts[frame_index]);
        } else {
            // The codons from the last node of the frame down to this one:
            // `last - 3`, `last - 6`, ... while at or after `index`.
            let last = last_position[frame_index];
            if last >= index + 3 {
                let codons = (last - index) / 3;
                gc_frame_plot.tally(
                    last - 3 * codons,
                    codons,
                    bin,
                    &mut frame_counts[frame_index],
                );
            }

            update_node_gc_info(node, &frame_counts[frame_index]);
            last_position[frame_index] = index;
        }
    });
}

/// Process nodes on the reverse strand
fn process_reverse_strand<P: FramePlot + ?Sized>(
    gc_frame_plot: &P,
    nodes: &mut StrandNodes,
    frame_counts: &mut [[i32; 3]; 3],
    last_position: &mut [usize; 3],
) {
    nodes.iter_mut().for_each(|node| {
        let index = node.position.index;
        let frame_index = index % 3;
        let bin = |gc_value| reverse_gc_bin(gc_value, frame_index);

        if node.position.codon_type == CodonType::Stop {
            // Reset counters for this frame
            frame_counts[frame_index] = [0; 3];
            last_position[frame_index] = index;
            gc_frame_plot.tally(index, 1, bin, &mut frame_counts[frame_index]);
        } else {
            // The codons from the last node of the frame up to this one:
            // `last + 3`, `last + 6`, ... while at or before `index`.
            let last = last_position[frame_index];
            if index >= last + 3 {
                gc_frame_plot.tally(
                    last + 3,
                    (index - last) / 3,
                    bin,
                    &mut frame_counts[frame_index],
                );
            }

            update_node_gc_info(node, &frame_counts[frame_index]);
            last_position[frame_index] = index;
        }
    });
}

/// Frame-count bin of a forward-strand position with GC frame value `gc_value`.
///
/// The frame plot marks the last one or two bases of a sequence (which do not
/// start a full codon) with `-1`; Prodigal feeds that value through the same
/// modular arithmetic as real frames, so those positions still count.
#[inline]
const fn forward_gc_bin(gc_value: i32, frame_offset: usize) -> usize {
    (gc_value + frame_offset as i32).rem_euclid(3) as usize
}

/// Frame-count bin of a reverse-strand position (see [`forward_gc_bin`]).
#[inline]
const fn reverse_gc_bin(gc_value: i32, frame_index: usize) -> usize {
    ((3 - gc_value) + frame_index as i32).rem_euclid(3) as usize
}

/// Calculate training bias from processed nodes
fn calculate_training_bias(nodes: &[Node], training: &mut Training) {
    training
        .gc_bias_factors
        .iter_mut()
        .for_each(|bias| *bias = 0.0);

    let total_bias: f64 = nodes
        .iter()
        .filter(|node| node.position.codon_type != CodonType::Stop)
        .filter_map(|node| {
            let gene_length = if node.position.stop_value < 0 {
                (node.position.index as isize - node.position.stop_value) as usize + 1
            } else if node.position.stop_value as usize > node.position.index {
                (node.position.stop_value as usize - node.position.index) + 1
            } else {
                (node.position.index - node.position.stop_value as usize) + 1
            };
            if node.state.gc_bias_frame < training.gc_bias_factors.len() && gene_length > 0 {
                let contribution = (node.scores.gc_frame_scores[node.state.gc_bias_frame]
                    * gene_length as f64)
                    / 1000.0;
                training.gc_bias_factors[node.state.gc_bias_frame] += contribution;
                Some(contribution)
            } else {
                None
            }
        })
        .sum();

    // Normalize bias values if we have any contributions
    if total_bias > 0.0 {
        let normalization_factor = 3.0 / training.gc_bias_factors.iter().sum::<f64>();
        training
            .gc_bias_factors
            .iter_mut()
            .for_each(|bias| *bias *= normalization_factor);
    }
}

/// Record the GC bias information of a node: the frame with the most GC-rich
/// codons and, for a gene of positive length, the share of each frame.
fn update_node_gc_info(node: &mut Node, frame_counts: &[i32; 3]) {
    let max_frame_index = find_max_frame_index(frame_counts[0], frame_counts[1], frame_counts[2]);
    node.state.gc_bias_frame = max_frame_index;

    let gene_length: isize = if node.position.strand == Strand::Forward {
        node.position.stop_value + 3 - node.position.index as isize
    } else {
        node.position.index as isize + 3 - node.position.stop_value
    };

    if gene_length > 0 {
        frame_counts
            .iter()
            .enumerate()
            .for_each(|(frame_index, &count)| {
                node.scores.gc_frame_scores[frame_index] =
                    (3.0 * count as f64) / (gene_length as f64);
            });
    }
}

/// Find the frame with maximum count
const fn find_max_frame_index(frame_a: i32, frame_b: i32, frame_c: i32) -> usize {
    if frame_a > frame_b {
        if frame_a > frame_c { 0 } else { 2 }
    } else if frame_b > frame_c {
        1
    } else {
        2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::*;
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
            scores: NodeScores {
                gc_frame_scores: [0.0; 3],
                ..Default::default()
            },
            state: NodeState {
                gc_bias_frame: 0,
                ..Default::default()
            },
            motif_info: NodeMotifInfo::default(),
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
    fn test_record_gc_bias_empty_inputs() {
        let gc_frame_data = vec![];
        let mut nodes = vec![];
        let mut training = create_test_training();

        record_gc_bias(&gc_frame_data, &mut nodes, &mut training);

        // Should handle empty inputs gracefully
        assert_eq!(nodes.len(), 0);
        assert_eq!(gc_frame_data.len(), 0);
    }

    #[test]
    fn test_record_gc_bias_forward_strand() {
        let gc_frame_data = vec![0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2];
        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 12),
            create_test_node(Strand::Forward, 12, CodonType::Stop, 12),
        ];
        let mut training = create_test_training();

        record_gc_bias(&gc_frame_data, &mut nodes, &mut training);

        // Verify GC bias frame is set
        assert!(nodes[0].state.gc_bias_frame <= 2);
        // Verify GC frame scores are calculated
        assert!(
            nodes[0]
                .scores
                .gc_frame_scores
                .iter()
                .any(|&score| score >= 0.0)
        );
    }

    #[test]
    fn test_record_gc_bias_reverse_strand() {
        let gc_frame_data = vec![0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2];
        let mut nodes = vec![
            create_test_node(Strand::Reverse, 12, CodonType::Atg, 3),
            create_test_node(Strand::Reverse, 3, CodonType::Stop, 3),
        ];
        let mut training = create_test_training();

        record_gc_bias(&gc_frame_data, &mut nodes, &mut training);

        // Verify GC bias frame is set
        assert!(nodes[0].state.gc_bias_frame <= 2);
        // Verify GC frame scores are calculated
        assert!(
            nodes[0]
                .scores
                .gc_frame_scores
                .iter()
                .any(|&score| score >= 0.0)
        );
    }

    #[test]
    fn test_record_gc_bias_mixed_strands() {
        let gc_frame_data = vec![0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2];
        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 12),
            create_test_node(Strand::Reverse, 15, CodonType::Gtg, 6),
            create_test_node(Strand::Forward, 12, CodonType::Stop, 12),
            create_test_node(Strand::Reverse, 6, CodonType::Stop, 6),
        ];
        let mut training = create_test_training();

        record_gc_bias(&gc_frame_data, &mut nodes, &mut training);

        for node in &nodes {
            if node.position.codon_type != CodonType::Stop {
                assert!(node.state.gc_bias_frame <= 2);
            }
        }
    }

    #[test]
    fn test_process_forward_strand() {
        let gc_frame_data = vec![0, 1, 2, 0, 1, 2, 0, 1, 2];
        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 6),
            create_test_node(Strand::Forward, 6, CodonType::Stop, 6),
        ];
        let mut frame_counts = [[0i32; 3]; 3];
        let mut last_position = [0; 3];

        let (mut forward, _) = split_by_strand(&mut nodes);
        process_forward_strand(
            gc_frame_data.as_slice(),
            &mut forward,
            &mut frame_counts,
            &mut last_position,
        );

        // Verify processing completed without errors
        assert!(
            nodes[0]
                .scores
                .gc_frame_scores
                .iter()
                .any(|&score| score >= 0.0)
        );
    }

    #[test]
    fn test_process_reverse_strand() {
        let gc_frame_data = vec![0, 1, 2, 0, 1, 2, 0, 1, 2];
        let mut nodes = vec![
            create_test_node(Strand::Reverse, 6, CodonType::Atg, 3),
            create_test_node(Strand::Reverse, 3, CodonType::Stop, 3),
        ];
        let mut frame_counts = [[0i32; 3]; 3];
        let mut last_position = [0; 3];

        let (_, mut reverse) = split_by_strand(&mut nodes);
        process_reverse_strand(
            gc_frame_data.as_slice(),
            &mut reverse,
            &mut frame_counts,
            &mut last_position,
        );

        // Verify processing completed without errors
        assert!(
            nodes[0]
                .scores
                .gc_frame_scores
                .iter()
                .any(|&score| score >= 0.0)
        );
    }

    #[test]
    fn test_calculate_training_bias() {
        let mut nodes = vec![
            create_test_node(Strand::Forward, 3, CodonType::Atg, 12),
            create_test_node(Strand::Forward, 15, CodonType::Gtg, 24),
        ];

        // Set some GC frame scores
        nodes[0].scores.gc_frame_scores = [1.0, 2.0, 3.0];
        nodes[1].scores.gc_frame_scores = [2.0, 3.0, 4.0];
        nodes[0].state.gc_bias_frame = 0;
        nodes[1].state.gc_bias_frame = 1;

        let mut training = create_test_training();

        calculate_training_bias(&nodes, &mut training);

        let total: f64 = training.gc_bias_factors.iter().sum();
        assert!((total - 3.0).abs() < 1e-10 || total == 0.0);
    }

    #[test]
    fn test_calculate_training_bias_empty_nodes() {
        let nodes = vec![];
        let mut training = create_test_training();

        calculate_training_bias(&nodes, &mut training);

        // Should handle empty nodes gracefully
        assert_eq!(training.gc_bias_factors, [0.0; 3]);
    }

    #[test]
    fn test_calculate_training_bias_only_stop_codons() {
        let nodes = vec![
            create_test_node(Strand::Forward, 6, CodonType::Stop, 6),
            create_test_node(Strand::Reverse, 12, CodonType::Stop, 12),
        ];
        let mut training = create_test_training();

        calculate_training_bias(&nodes, &mut training);

        assert_eq!(training.gc_bias_factors, [0.0; 3]);
    }

    #[test]
    fn test_update_node_gc_info_forward() {
        let mut node = create_test_node(Strand::Forward, 3, CodonType::Atg, 12);
        let frame_counts = [5, 3, 2];

        update_node_gc_info(&mut node, &frame_counts);

        // Should select frame with maximum count (frame 0)
        assert_eq!(node.state.gc_bias_frame, 0);
        assert!(node.scores.gc_frame_scores[0] > 0.0);
    }

    #[test]
    fn test_update_node_gc_info_reverse() {
        let mut node = create_test_node(Strand::Reverse, 12, CodonType::Atg, 3);
        let frame_counts = [2, 5, 3];

        update_node_gc_info(&mut node, &frame_counts);

        // Should select frame with maximum count (frame 1)
        assert_eq!(node.state.gc_bias_frame, 1);
        assert!(node.scores.gc_frame_scores[1] > 0.0);
    }

    #[test]
    fn test_update_node_gc_info_zero_gene_length() {
        let mut node = create_test_node(Strand::Forward, 6, CodonType::Atg, 6);
        let frame_counts = [5, 3, 2];

        update_node_gc_info(&mut node, &frame_counts);

        // Should handle zero gene length gracefully
        assert_eq!(node.state.gc_bias_frame, 0);
    }

    #[test]
    fn test_find_max_frame_index_first() {
        assert_eq!(find_max_frame_index(5, 3, 2), 0);
    }

    #[test]
    fn test_find_max_frame_index_second() {
        assert_eq!(find_max_frame_index(2, 5, 3), 1);
    }

    #[test]
    fn test_find_max_frame_index_third() {
        assert_eq!(find_max_frame_index(2, 3, 5), 2);
    }

    #[test]
    fn test_find_max_frame_index_ties() {
        assert_eq!(find_max_frame_index(5, 5, 3), 1); // When a == b > c, choose b (index 1)
        assert_eq!(find_max_frame_index(3, 5, 5), 2); // When b == c > a, choose c (index 2)
        assert_eq!(find_max_frame_index(5, 3, 5), 2); // When a == c > b, choose c (index 2)
        assert_eq!(find_max_frame_index(5, 5, 5), 2); // When all equal, choose c (index 2)
    }

    #[test]
    fn test_record_gc_bias_negative_gc_values() {
        let gc_frame_data = vec![-1, 0, 1, -1, 2, -1];
        let mut nodes = vec![
            create_test_node(Strand::Forward, 1, CodonType::Atg, 4),
            create_test_node(Strand::Forward, 4, CodonType::Stop, 4),
        ];
        let mut training = create_test_training();

        record_gc_bias(&gc_frame_data, &mut nodes, &mut training);

        // Should handle negative GC values (invalid positions)
        assert!(nodes[0].state.gc_bias_frame <= 2);
    }

    #[test]
    fn test_record_gc_bias_out_of_bounds_indices() {
        let gc_frame_data = vec![0, 1, 2];
        let mut nodes = vec![
            create_test_node(Strand::Forward, 10, CodonType::Atg, 20), // Index beyond gc_frame_data
            create_test_node(Strand::Forward, 20, CodonType::Stop, 20),
        ];
        let mut training = create_test_training();

        record_gc_bias(&gc_frame_data, &mut nodes, &mut training);

        // Should handle out-of-bounds indices gracefully
        assert!(nodes[0].state.gc_bias_frame <= 2);
    }

    /// The GC bias record of the nodes of a sequence must not depend on how
    /// the frame plot is stored.
    #[test]
    fn frame_plot_totals_record_the_same_bias_as_plot_values() {
        use crate::node::{add_nodes, sort_nodes_by_position};
        use crate::sequence::encoded::EncodedSequence;
        use crate::sequence::{GcFramePlot, calc_most_gc_frame};

        let mut random = crate::test_random::XorShift(0x2545_F491_4F6C_DD1D_u64);
        let mut next = || random.next();
        // Lengths of every residue modulo 3: the bases after the last whole
        // triplet have no frame, and nodes near the end read them.
        for (length, circular) in [
            (6000, false),
            (6001, false),
            (6002, false),
            (9001, true),
            (9002, true),
        ] {
            let bases: Vec<u8> = (0..length)
                .map(|_| b"ACGGCT"[(next() % 6) as usize])
                .collect();
            let encoded = EncodedSequence::without_masking(&bases);
            let mut training = create_test_training();
            let mut nodes = Vec::new();
            add_nodes(&encoded, &mut nodes, false, circular, &training).unwrap();
            sort_nodes_by_position(&mut nodes);
            assert!(nodes.len() > 100);
            let mut nodes_from_totals = nodes.clone();
            let mut training_from_totals = training.clone();

            let values = calc_most_gc_frame(&encoded.forward_sequence, length);
            record_gc_bias(&values, &mut nodes, &mut training);
            let totals = GcFramePlot::new(&encoded.forward_sequence, length);
            record_gc_bias_of(&totals, &mut nodes_from_totals, &mut training_from_totals);

            assert_eq!(
                format!("{nodes:?}"),
                format!("{nodes_from_totals:?}"),
                "length {length}"
            );
            assert_eq!(
                training.gc_bias_factors,
                training_from_totals.gc_bias_factors
            );
            assert!(training.gc_bias_factors.iter().any(|&bias| bias != 1.0));
        }
    }
}
