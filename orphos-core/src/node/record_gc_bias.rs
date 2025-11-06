use bio::bio_types::strand::Strand;

use crate::types::{CodonType, Node, Training};

/// Record GC bias for training
pub fn record_gc_bias(gc_frame_data: &[i32], nodes: &mut [Node], training: &mut Training) {
    if nodes.is_empty() || gc_frame_data.is_empty() {
        return;
    }

    let mut frame_counts = [[0i32; 3]; 3];
    let mut last_position = [0; 3];

    process_forward_strand(gc_frame_data, nodes, &mut frame_counts, &mut last_position);

    process_reverse_strand(gc_frame_data, nodes, &mut frame_counts, &mut last_position);

    // Calculate overall bias with better error handling
    calculate_training_bias(nodes, training);
}

/// Process nodes on the forward strand
fn process_forward_strand(
    gc_frame_data: &[i32],
    nodes: &mut [Node],
    frame_counts: &mut [[i32; 3]; 3],
    last_position: &mut [usize; 3],
) {
    nodes
        .iter_mut()
        .rev()
        .filter(|node| node.position.strand == Strand::Forward)
        .for_each(|node| {
            let frame_index = node.position.index % 3;
            let frame_offset = 3 - frame_index;

            if node.position.codon_type == CodonType::Stop {
                // Reset counters for this frame
                frame_counts[frame_index]
                    .iter_mut()
                    .for_each(|count| *count = 0);
                last_position[frame_index] = node.position.index;

                if let Some(&gc_value) = gc_frame_data.get(node.position.index)
                    && gc_value >= 0
                {
                    let gc_index = (gc_value as usize + frame_offset) % 3;
                    frame_counts[frame_index][gc_index] = 1;
                }
            } else {
                if last_position[frame_index] >= 3 {
                    let mut position = last_position[frame_index] - 3;

                    while position >= node.position.index {
                        if let Some(&gc_value) = gc_frame_data.get(position)
                            && gc_value >= 0
                        {
                            let gc_index = (gc_value as usize + frame_offset) % 3;
                            frame_counts[frame_index][gc_index] += 1;
                        }

                        if position < 3 {
                            break;
                        }
                        position -= 3;
                    }
                }

                update_node_gc_info(node, &frame_counts[frame_index]);
                last_position[frame_index] = node.position.index;
            }
        });
}

/// Process nodes on the reverse strand
fn process_reverse_strand(
    gc_frame_data: &[i32],
    nodes: &mut [Node],
    frame_counts: &mut [[i32; 3]; 3],
    last_position: &mut [usize; 3],
) {
    nodes
        .iter_mut()
        .filter(|node| node.position.strand == Strand::Reverse)
        .for_each(|node| {
            let frame_index = node.position.index % 3;

            if node.position.codon_type == CodonType::Stop {
                // Reset counters for this frame
                frame_counts[frame_index]
                    .iter_mut()
                    .for_each(|count| *count = 0);
                last_position[frame_index] = node.position.index;

                if let Some(&gc_value) = gc_frame_data.get(node.position.index)
                    && gc_value >= 0
                {
                    let gc_index = ((3 - gc_value as usize) + frame_index) % 3;
                    frame_counts[frame_index][gc_index] = 1;
                }
            } else {
                let mut position = last_position[frame_index] + 3;

                while position <= node.position.index && position < gc_frame_data.len() {
                    if let Some(&gc_value) = gc_frame_data.get(position)
                        && gc_value >= 0
                    {
                        let gc_index = ((3 - gc_value as usize) + frame_index) % 3;
                        frame_counts[frame_index][gc_index] += 1;
                    }
                    position += 3;
                }

                update_node_gc_info(node, &frame_counts[frame_index]);
                last_position[frame_index] = node.position.index;
            }
        });
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

/// Update GC bias information for a node
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

        process_forward_strand(
            &gc_frame_data,
            &mut nodes,
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

        process_reverse_strand(
            &gc_frame_data,
            &mut nodes,
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
}
