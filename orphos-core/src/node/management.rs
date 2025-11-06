use bio::bio_types::strand::Strand;

use crate::sequence::is_gc;
use crate::types::{CodonType, Node};

pub fn calc_orf_gc(encoded_sequence: &[u8], sequence_length: usize, nodes: &mut [Node]) {
    // Early termination for empty nodes or sequences
    if nodes.is_empty() || sequence_length == 0 {
        return;
    }

    process_strand(encoded_sequence, sequence_length, nodes, Strand::Forward);
    process_strand(encoded_sequence, sequence_length, nodes, Strand::Reverse);
}

fn process_strand(
    encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &mut [Node],
    strand: Strand,
) {
    let mut gc = [0.0f64; 3];
    let mut last = [0usize; 3];

    // Process nodes in the appropriate order without heap allocation
    match strand {
        Strand::Forward => {
            for node in nodes.iter_mut().rev() {
                if node.position.strand != strand {
                    continue;
                }

                let frame = node.position.index % 3;

                if node.position.codon_type == CodonType::Stop {
                    update_stop_codon_gc(
                        encoded_sequence,
                        sequence_length,
                        node,
                        &mut gc,
                        &mut last,
                        frame,
                        strand,
                    );
                } else {
                    update_start_codon_gc(
                        encoded_sequence,
                        sequence_length,
                        node,
                        &mut gc,
                        &mut last,
                        frame,
                        strand,
                    );
                }
            }
        }
        Strand::Reverse => {
            for node in nodes.iter_mut() {
                if node.position.strand != strand {
                    continue;
                }

                let frame = node.position.index % 3;

                if node.position.codon_type == CodonType::Stop {
                    update_stop_codon_gc(
                        encoded_sequence,
                        sequence_length,
                        node,
                        &mut gc,
                        &mut last,
                        frame,
                        strand,
                    );
                } else {
                    update_start_codon_gc(
                        encoded_sequence,
                        sequence_length,
                        node,
                        &mut gc,
                        &mut last,
                        frame,
                        strand,
                    );
                }
            }
        }
        Strand::Unknown => (), // Handle gracefully instead of panic
    }
}

fn update_stop_codon_gc(
    seq: &[u8],
    sequence_length: usize,
    node: &Node,
    gc: &mut [f64; 3],
    last: &mut [usize; 3],
    frame: usize,
    strand: Strand,
) {
    last[frame] = node.position.index;
    gc[frame] = calculate_codon_gc(seq, sequence_length, node.position.index, strand);
}

fn update_start_codon_gc(
    seq: &[u8],
    sequence_length: usize,
    node: &mut Node,
    gc: &mut [f64; 3],
    last: &mut [usize; 3],
    frame: usize,
    strand: Strand,
) {
    match strand {
        Strand::Forward => {
            if last[frame] >= 3 {
                gc[frame] +=
                    accumulate_gc_forward(seq, sequence_length, node.position.index, last[frame]);
            }
        }
        Strand::Reverse => {
            gc[frame] +=
                accumulate_gc_reverse(seq, sequence_length, node.position.index, last[frame]);
        }
        Strand::Unknown => todo!(),
    }

    // Calculate GC content for this ORF
    let orf_length = calculate_orf_length(node.position.stop_value, node.position.index);
    node.scores.gc_content = if orf_length > 0.0 {
        gc[frame] / orf_length
    } else {
        0.0
    };

    last[frame] = node.position.index;
}

fn calculate_codon_gc(
    encoded_sequence: &[u8],
    sequence_length: usize,
    index: usize,
    strand: Strand,
) -> f64 {
    let positions = match strand {
        Strand::Forward => [index, index + 1, index + 2],
        Strand::Reverse => [index, index.saturating_sub(1), index.saturating_sub(2)],
        Strand::Unknown => return 0.0,
    };

    // Optimized GC counting with early bounds check
    if strand == Strand::Forward && index + 2 >= sequence_length {
        // For forward strand, check if we're near the end
        return positions
            .iter()
            .take_while(|&&pos| pos < sequence_length)
            .map(|&pos| {
                if is_gc(encoded_sequence, pos) {
                    1.0
                } else {
                    0.0
                }
            })
            .sum();
    }

    if strand == Strand::Reverse && index < 2 {
        // For reverse strand, check if we're near the beginning
        return positions
            .iter()
            .filter(|&&pos| pos <= index)
            .map(|&pos| {
                if is_gc(encoded_sequence, pos) {
                    1.0
                } else {
                    0.0
                }
            })
            .sum();
    }

    // Normal case - all positions are valid
    positions
        .iter()
        .map(|&pos| {
            if is_gc(encoded_sequence, pos) {
                1.0
            } else {
                0.0
            }
        })
        .sum()
}

fn accumulate_gc_forward(
    encoded_sequence: &[u8],
    sequence_length: usize,
    current_index: usize,
    last_index: usize,
) -> f64 {
    if last_index < 3 || last_index <= current_index {
        return 0.0;
    }

    let mut gc_count = 0.0;
    let end = last_index.saturating_sub(3);

    // Process in steps of 3, going backwards
    let mut j = end;
    while j >= current_index {
        gc_count += calculate_triplet_gc(encoded_sequence, sequence_length, j);
        if j < 3 {
            break;
        }
        j = j.saturating_sub(3);
    }

    gc_count
}

fn accumulate_gc_reverse(
    encoded_sequence: &[u8],
    sequence_length: usize,
    current_index: usize,
    last_index: usize,
) -> f64 {
    if last_index + 3 > current_index {
        return 0.0;
    }

    let mut gc_count = 0.0;
    let start = last_index + 3;
    let end = current_index.min(sequence_length.saturating_sub(3));

    // Process in steps of 3, going forwards
    let mut j = start;
    while j <= end && j < sequence_length {
        gc_count += calculate_triplet_gc(encoded_sequence, sequence_length, j);
        j += 3;
    }

    gc_count
}

fn calculate_triplet_gc(encoded_sequence: &[u8], sequence_length: usize, start_pos: usize) -> f64 {
    // Early bounds check to avoid unnecessary computation
    if start_pos >= sequence_length {
        return 0.0;
    }

    let end_pos = (start_pos + 3).min(sequence_length);
    let mut gc_count = 0.0;

    // Unrolled loop for better performance when we have all 3 positions
    if start_pos + 2 < sequence_length {
        // Fast path: all 3 positions are valid
        if is_gc(encoded_sequence, start_pos) {
            gc_count += 1.0;
        }
        if is_gc(encoded_sequence, start_pos + 1) {
            gc_count += 1.0;
        }
        if is_gc(encoded_sequence, start_pos + 2) {
            gc_count += 1.0;
        }
    } else {
        // Slow path: handle boundary cases
        for pos in start_pos..end_pos {
            if is_gc(encoded_sequence, pos) {
                gc_count += 1.0;
            }
        }
    }

    gc_count
}

const fn calculate_orf_length(stop_value: isize, index: usize) -> f64 {
    // stop_value may be negative (edge sentinel). Treat distance as (index - stop_value) when negative.
    let length = if stop_value < 0 {
        (index as isize - stop_value) as usize
    } else {
        (stop_value as usize).abs_diff(index)
    };
    (length + 3) as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::encode_sequence;
    use crate::types::{NodeMotifInfo, NodePosition, NodeScores, NodeState};

    fn get_encoded_sequence(input: &[u8]) -> Vec<u8> {
        let sequence_length = input.len();
        let mut seq = vec![0u8; (sequence_length * 2).div_ceil(8)];
        let mut unknown_sequence = vec![0u8; sequence_length.div_ceil(8)];
        let mut masks = Vec::new();
        let _ = encode_sequence(input, &mut seq, &mut unknown_sequence, &mut masks, false).unwrap();
        seq
    }

    fn create_test_node(
        index: usize,
        strand: Strand,
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
            motif_info: NodeMotifInfo::default(),
        }
    }

    #[test]
    fn test_calc_orf_gc_empty_nodes() {
        let sequence = b"ATGCGCGCGCGC";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes: Vec<Node> = vec![];

        // Should not panic with empty nodes
        calc_orf_gc(&encoded_seq, sequence.len(), &mut nodes);
        assert!(nodes.is_empty());
    }

    #[test]
    fn test_calc_orf_gc_forward_nodes() {
        let sequence = b"ATGCGCGCGCGC";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes = vec![
            create_test_node(0, Strand::Forward, CodonType::Atg, 9),
            create_test_node(9, Strand::Forward, CodonType::Stop, 0),
        ];

        calc_orf_gc(&encoded_seq, sequence.len(), &mut nodes);

        // Check that GC content was calculated for the start node
        let start_node = &nodes[0];
        assert!(start_node.scores.gc_content >= 0.0);
        assert!(start_node.scores.gc_content <= 1.0);
    }

    #[test]
    fn test_calc_orf_gc_reverse_nodes() {
        let sequence = b"ATGCGCGCGCGC";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes = vec![
            create_test_node(9, Strand::Reverse, CodonType::Stop, 0),
            create_test_node(0, Strand::Reverse, CodonType::Atg, 9),
        ];

        calc_orf_gc(&encoded_seq, sequence.len(), &mut nodes);

        // Check that GC content was calculated for the start node
        let start_node = &nodes[1];
        assert!(start_node.scores.gc_content >= 0.0);
        assert!(start_node.scores.gc_content <= 1.0);
    }

    #[test]
    fn test_calc_orf_gc_mixed_strands() {
        let sequence = b"ATGCGCGCGCGC";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes = vec![
            create_test_node(0, Strand::Forward, CodonType::Atg, 6),
            create_test_node(6, Strand::Forward, CodonType::Stop, 0),
            create_test_node(9, Strand::Reverse, CodonType::Stop, 3),
            create_test_node(3, Strand::Reverse, CodonType::Atg, 9),
        ];

        calc_orf_gc(&encoded_seq, sequence.len(), &mut nodes);

        assert!(nodes[0].scores.gc_content >= 0.0);
        assert!(nodes[3].scores.gc_content >= 0.0);
    }

    #[test]
    fn test_calculate_codon_gc_forward() {
        let sequence = b"GGG"; // All GC
        let encoded_seq = get_encoded_sequence(sequence);

        let gc_content = calculate_codon_gc(&encoded_seq, sequence.len(), 0, Strand::Forward);
        assert_eq!(gc_content, 3.0);
    }

    #[test]
    fn test_calculate_codon_gc_reverse() {
        let sequence = b"GGG"; // All GC
        let encoded_seq = get_encoded_sequence(sequence);

        let gc_content = calculate_codon_gc(&encoded_seq, sequence.len(), 2, Strand::Reverse);
        assert_eq!(gc_content, 3.0);
    }

    #[test]
    fn test_calculate_codon_gc_no_gc() {
        let sequence = b"AAA"; // No GC
        let encoded_seq = get_encoded_sequence(sequence);

        let gc_content = calculate_codon_gc(&encoded_seq, sequence.len(), 0, Strand::Forward);
        assert_eq!(gc_content, 0.0);
    }

    #[test]
    fn test_calculate_triplet_gc() {
        let sequence = b"GCA"; // 2 GC out of 3
        let encoded_seq = get_encoded_sequence(sequence);

        let gc_content = calculate_triplet_gc(&encoded_seq, sequence.len(), 0);
        assert_eq!(gc_content, 2.0);
    }

    #[test]
    fn test_calculate_triplet_gc_edge_case() {
        let sequence = b"GC"; // Only 2 nucleotides
        let encoded_seq = get_encoded_sequence(sequence);

        let gc_content = calculate_triplet_gc(&encoded_seq, sequence.len(), 0);
        assert_eq!(gc_content, 2.0);
    }

    #[test]
    fn test_calculate_orf_length() {
        assert_eq!(calculate_orf_length(9, 0), 12.0); // 9 - 0 + 3 = 12
        assert_eq!(calculate_orf_length(0, 9), 12.0); // 9 - 0 + 3 = 12
        assert_eq!(calculate_orf_length(5, 5), 3.0); // Same position + 3
        assert_eq!(calculate_orf_length(-6, 0), 9.0); // edge case negative stop
    }

    #[test]
    fn test_accumulate_gc_forward() {
        let sequence = b"GGGAAAGGG"; // GC at start and end
        let encoded_seq = get_encoded_sequence(sequence);

        let gc_count = accumulate_gc_forward(&encoded_seq, sequence.len(), 3, 9);
        assert!(gc_count >= 0.0);
    }

    #[test]
    fn test_accumulate_gc_reverse() {
        let sequence = b"GGGAAAGGG"; // GC at start and end
        let encoded_seq = get_encoded_sequence(sequence);

        let gc_count = accumulate_gc_reverse(&encoded_seq, sequence.len(), 6, 0);
        assert!(gc_count >= 0.0);
    }

    #[test]
    fn test_process_strand_forward_only() {
        let sequence = b"ATGCGCGCGCGC";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes = vec![
            create_test_node(0, Strand::Forward, CodonType::Atg, 9),
            create_test_node(9, Strand::Forward, CodonType::Stop, 0),
        ];

        process_strand(&encoded_seq, sequence.len(), &mut nodes, Strand::Forward);

        // Should process forward strand nodes
        assert!(nodes[0].scores.gc_content >= 0.0);
    }

    #[test]
    fn test_process_strand_reverse_only() {
        let sequence = b"ATGCGCGCGCGC";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes = vec![
            create_test_node(9, Strand::Reverse, CodonType::Stop, 0),
            create_test_node(0, Strand::Reverse, CodonType::Atg, 9),
        ];

        process_strand(&encoded_seq, sequence.len(), &mut nodes, Strand::Reverse);

        // Should process reverse strand nodes
        assert!(nodes[1].scores.gc_content >= 0.0);
    }

    #[test]
    fn test_update_stop_codon_gc() {
        let sequence = b"GGGAAATAG";
        let encoded_seq = get_encoded_sequence(sequence);
        let node = create_test_node(6, Strand::Forward, CodonType::Stop, 0);
        let mut gc = [0.0; 3];
        let mut last = [0; 3];

        update_stop_codon_gc(
            &encoded_seq,
            sequence.len(),
            &node,
            &mut gc,
            &mut last,
            0,
            Strand::Forward,
        );

        assert_eq!(last[0], 6);
        assert!(gc[0] >= 0.0);
    }

    #[test]
    fn test_update_start_codon_gc() {
        let sequence = b"GGGAAATAG";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut node = create_test_node(0, Strand::Forward, CodonType::Atg, 6);
        let mut gc = [3.0; 3]; // Some initial GC count
        let mut last = [9; 3];

        update_start_codon_gc(
            &encoded_seq,
            sequence.len(),
            &mut node,
            &mut gc,
            &mut last,
            0,
            Strand::Forward,
        );

        assert_eq!(last[0], 0);
        assert!(node.scores.gc_content >= 0.0);
        assert!(node.scores.gc_content <= 1.0);
    }
}
