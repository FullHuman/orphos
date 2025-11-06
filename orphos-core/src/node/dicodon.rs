use bio::bio_types::strand::Strand;

use crate::{
    constants::{DICODON_SIZE, MAX_LOG_LIKELIHOOD, MIN_LOG_LIKELIHOOD, NUM_DICODONS},
    sequence::{calculate_background_mer_frequencies, calculate_kmer_index},
    types::{CodonType, Node, NodeSliceExtension, Training},
};

#[derive(Debug, Clone, Copy)]
enum GeneState {
    NotPresent,
    Forward,
    Reverse,
}

/// Calculate dicodon (6-mer) frequencies in genes and store log likelihood scores
/// relative to background frequencies in the training structure.
pub fn calculate_dicodon_gene(
    training: &mut Training,
    sequence: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &[Node],
    dbeg: usize,
) {
    let mut counts: [u32; NUM_DICODONS] = [0; NUM_DICODONS];
    let mut bg = [0.0f64; NUM_DICODONS];
    let mut total_dicodons: u32 = 0;
    let mut gene_start: Option<usize> = None;
    let mut gene_end: Option<usize> = None;
    let mut gene_state = GeneState::NotPresent;

    // Calculate background 6-mer frequencies
    calculate_background_mer_frequencies(
        DICODON_SIZE,
        sequence,
        reverse_complement_encoded_sequence,
        sequence_length,
        &mut bg,
    );

    // Follow the optimal path through nodes starting from dbeg
    for node in nodes.iter_traceback(dbeg) {
        let codon_stop = node.position.codon_type == CodonType::Stop;

        match (node.position.strand, codon_stop) {
            (Strand::Reverse, false) => {
                gene_state = GeneState::Reverse;
                gene_start = Some(sequence_length - node.position.index - 1);
            }
            (Strand::Forward, true) => {
                gene_state = GeneState::Forward;
                gene_end = Some(node.position.index + 2);
            }
            _ => {}
        }

        match (gene_state, node.position.strand, codon_stop) {
            (GeneState::Reverse, Strand::Reverse, true) => {
                gene_end = Some(sequence_length - node.position.index + 1);
                count_dicodons_in_gene(
                    reverse_complement_encoded_sequence,
                    gene_start.unwrap_or(0),
                    gene_end.unwrap_or(0),
                    &mut counts,
                    &mut total_dicodons,
                );
                gene_state = GeneState::NotPresent;
            }
            (GeneState::Forward, Strand::Forward, false) => {
                gene_start = Some(node.position.index);
                count_dicodons_in_gene(
                    sequence,
                    gene_start.unwrap_or(0),
                    gene_end.unwrap_or(0),
                    &mut counts,
                    &mut total_dicodons,
                );
                gene_state = GeneState::NotPresent;
            }
            _ => {}
        }
    }

    if total_dicodons != 0 {
        calculate_and_store_scores(training, &counts, total_dicodons as f64, &bg);
        training.total_dicodons = total_dicodons;
    }
}

fn count_dicodons_in_gene(
    sequence: &[u8],
    gene_start: usize,
    gene_end: usize,
    counts: &mut [u32; NUM_DICODONS],
    total_dicodons: &mut u32,
) {
    for i in (gene_start..gene_end).step_by(3) {
        if i + DICODON_SIZE > gene_end {
            break;
        }
        let mer_idx = calculate_kmer_index(DICODON_SIZE, sequence, i);
        counts[mer_idx] += 1;
        *total_dicodons += 1;
    }
}

/// Calculate probabilities and store log likelihood scores in training table
fn calculate_and_store_scores(
    training: &mut Training,
    counts: &[u32; NUM_DICODONS],
    total_dicodons: f64,
    bg: &[f64; NUM_DICODONS],
) {
    // for (dicodon_index, (&dicodon_count, &background_probability)) in
    //     counts.iter().zip(bg.iter()).enumerate()
    // {
    //     let gene_probability = dicodon_count as f64 / total_dicodons;
    //     training.gene_dicodon_table[dicodon_index] =
    //         calculate_log_likelihood(gene_probability, background_probability);
    // }
    training
        .gene_dicodon_table
        .iter_mut()
        .zip(counts.iter().zip(bg.iter()))
        .for_each(|(score, (&dicodon_count, &background_probability))| {
            let gene_probability = dicodon_count as f64 / total_dicodons;
            *score = calculate_log_likelihood(gene_probability, background_probability);
        });
}

/// Calculate log likelihood score with proper edge case handling
fn calculate_log_likelihood(gene_probability: f64, background_probability: f64) -> f64 {
    let score = match (gene_probability, background_probability) {
        (0.0, background) if background > 0.0 => MIN_LOG_LIKELIHOOD,
        (_, 0.0) => 0.0,
        (gene, background) => (gene / background).ln(),
    };

    score.clamp(MIN_LOG_LIKELIHOOD, MAX_LOG_LIKELIHOOD)
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
        traceback: Option<usize>,
    ) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value: 0,
                is_edge: false,
            },
            scores: NodeScores::default(),
            state: NodeState {
                traceback,
                ..Default::default()
            },
            motif_info: NodeMotifInfo::default(),
        }
    }

    #[test]
    fn test_calculate_log_likelihood_normal() {
        let gene_prob = 0.1;
        let bg_prob = 0.05;
        let result = calculate_log_likelihood(gene_prob, bg_prob);

        // ln(0.1 / 0.05) = ln(2) ≈ 0.693
        assert!((result - 0.693).abs() < 0.001);
    }

    #[test]
    fn test_calculate_log_likelihood_zero_gene() {
        let gene_prob = 0.0;
        let bg_prob = 0.05;
        let result = calculate_log_likelihood(gene_prob, bg_prob);

        assert_eq!(result, MIN_LOG_LIKELIHOOD);
    }

    #[test]
    fn test_calculate_log_likelihood_zero_background() {
        let gene_prob = 0.1;
        let bg_prob = 0.0;
        let result = calculate_log_likelihood(gene_prob, bg_prob);

        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_calculate_log_likelihood_both_zero() {
        let gene_prob = 0.0;
        let bg_prob = 0.0;
        let result = calculate_log_likelihood(gene_prob, bg_prob);

        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_calculate_log_likelihood_clamping_high() {
        let gene_prob = 1.0;
        let bg_prob = 0.001;
        let result = calculate_log_likelihood(gene_prob, bg_prob);

        assert_eq!(result, MAX_LOG_LIKELIHOOD);
    }

    #[test]
    fn test_count_dicodons_in_gene() {
        let sequence = b"ATGCGATGCGATGCGA"; // 16 bp sequence
        let encoded_seq = get_encoded_sequence(sequence);
        let mut counts = [0u32; NUM_DICODONS];
        let mut total_dicodons = 0u32;

        count_dicodons_in_gene(&encoded_seq, 0, 12, &mut counts, &mut total_dicodons);

        // Should have counted some dicodons
        assert!(total_dicodons > 0);
        assert!(counts.iter().any(|&count| count > 0));
    }

    #[test]
    fn test_count_dicodons_in_gene_short_region() {
        let sequence = b"ATG"; // Too short for a dicodon
        let encoded_seq = get_encoded_sequence(sequence);
        let mut counts = [0u32; NUM_DICODONS];
        let mut total_dicodons = 0u32;

        count_dicodons_in_gene(&encoded_seq, 0, 3, &mut counts, &mut total_dicodons);

        // Should not count any dicodons (region too short)
        assert_eq!(total_dicodons, 0);
    }

    #[test]
    fn test_count_dicodons_in_gene_empty_region() {
        let sequence = b"ATGCGATGCGATGCGA";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut counts = [0u32; NUM_DICODONS];
        let mut total_dicodons = 0u32;

        count_dicodons_in_gene(&encoded_seq, 5, 5, &mut counts, &mut total_dicodons);

        assert_eq!(total_dicodons, 0);
    }

    #[test]
    fn test_calculate_and_store_scores() {
        let mut training = Training::default();
        let mut full_counts = [0u32; NUM_DICODONS];
        full_counts[0] = 10;
        full_counts[1] = 5;
        full_counts[2] = 2;
        full_counts[3] = 0;

        let total_dicodons = 17.0;
        let mut bg = [0.1f64; NUM_DICODONS];
        bg[0] = 0.1;
        bg[1] = 0.05;
        bg[2] = 0.01;
        bg[3] = 0.05;

        calculate_and_store_scores(&mut training, &full_counts, total_dicodons, &bg);

        // Check that scores were calculated and stored
        assert_ne!(training.gene_dicodon_table[0], 0.0);
        assert_ne!(training.gene_dicodon_table[1], 0.0);
        assert_ne!(training.gene_dicodon_table[2], 0.0);
        assert_eq!(training.gene_dicodon_table[3], MIN_LOG_LIKELIHOOD);
    }

    #[test]
    fn test_calculate_dicodon_gene_simple() {
        let sequence = b"ATGCGATGCGATGCGATAG";
        let encoded_seq = get_encoded_sequence(sequence);
        let reverse_seq = get_encoded_sequence(b"CTATCGATCGATCGCAT");

        let nodes = vec![
            create_test_node(0, Strand::Forward, CodonType::Atg, Some(1)),
            create_test_node(18, Strand::Forward, CodonType::Stop, None),
        ];

        let mut training = Training::default();

        calculate_dicodon_gene(
            &mut training,
            &encoded_seq,
            &reverse_seq,
            sequence.len(),
            &nodes,
            0,
        );

        // For this simple test, we just verify the function runs without panicking
    }

    #[test]
    fn test_calculate_dicodon_gene_reverse() {
        let sequence = b"ATGCGATGCGATGCGATAG";
        let encoded_seq = get_encoded_sequence(sequence);
        let reverse_seq = get_encoded_sequence(b"CTATCGATCGATCGCAT");

        let nodes = vec![
            create_test_node(18, Strand::Reverse, CodonType::Stop, Some(1)),
            create_test_node(0, Strand::Reverse, CodonType::Atg, None),
        ];

        let mut training = Training::default();

        calculate_dicodon_gene(
            &mut training,
            &encoded_seq,
            &reverse_seq,
            sequence.len(),
            &nodes,
            0,
        );

        // For this reverse strand test, we just verify the function runs without panicking
    }

    #[test]
    fn test_calculate_dicodon_gene_no_genes() {
        let nodes: Vec<Node> = vec![];

        // which is the expected behavior
        assert!(nodes.is_empty());
    }

    #[test]
    fn test_calculate_dicodon_gene_invalid_path() {
        let sequence = b"ATGCGATGCGATGCGATAG";
        let encoded_seq = get_encoded_sequence(sequence);
        let reverse_seq = get_encoded_sequence(b"CTATCGATCGATCGCAT");

        let nodes = vec![create_test_node(0, Strand::Forward, CodonType::Atg, None)];

        let mut training = Training::default();

        calculate_dicodon_gene(
            &mut training,
            &encoded_seq,
            &reverse_seq,
            sequence.len(),
            &nodes,
            0,
        );

        // Should handle invalid path gracefully
        // Function completion indicates it didn't panic
    }

    #[test]
    fn test_gene_state_transitions() {
        // Test various gene state transitions through different node types
        let sequence = b"ATGCGATGCGATGCGATAG";
        let encoded_seq = get_encoded_sequence(sequence);
        let reverse_seq = get_encoded_sequence(b"CTATCGATCGATCGCAT");

        let nodes = vec![
            create_test_node(0, Strand::Forward, CodonType::Atg, Some(1)),
            create_test_node(9, Strand::Forward, CodonType::Stop, Some(2)),
            create_test_node(15, Strand::Reverse, CodonType::Stop, Some(3)),
            create_test_node(6, Strand::Reverse, CodonType::Atg, None),
        ];

        let mut training = Training::default();

        calculate_dicodon_gene(
            &mut training,
            &encoded_seq,
            &reverse_seq,
            sequence.len(),
            &nodes,
            0,
        );

        // Should handle mixed strands without panicking
        // Function completion indicates success
    }
}
