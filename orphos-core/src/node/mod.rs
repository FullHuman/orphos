//! Node management and scoring for gene prediction
//!
//! This module contains all functionality related to creating, scoring, and managing
//! gene prediction nodes. Nodes represent potential start and stop codons in the
//! sequence and are used in dynamic programming to find optimal gene predictions.

mod creation;
mod dicodon;
mod management;
mod motifs;
mod overlaps;
mod record_gc_bias;
mod scoring;
mod utilities;
mod validation;

pub use creation::add_nodes;
pub use dicodon::calculate_dicodon_gene;
pub use management::*;
pub use motifs::{find_best_upstream_motif, rbs_score};
pub use overlaps::{intergenic_mod, record_overlapping_starts};
pub use record_gc_bias::record_gc_bias;
pub use scoring::{raw_coding_score, score_nodes};
pub use utilities::{reset_node_scores, sort_nodes_by_position};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::*;
    use bio::bio_types::strand::Strand;

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

    fn create_test_node() -> Node {
        Node {
            position: NodePosition {
                index: 25,
                stop_value: 35,
                strand: Strand::Forward,
                codon_type: CodonType::Atg,
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

    #[test]
    fn test_module_exports_exist() {
        // Test that all the main module exports are accessible

        use crate::sequence::encoded::EncodedSequence;

        let seq = vec![0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 3, 0, 0];
        let reverse_seq = vec![3; seq.len()];
        let training = create_test_training();
        let mut nodes = vec![create_test_node()];

        // Create encoded sequence for the new API
        let encoded_sequence = EncodedSequence {
            forward_sequence: seq.clone(),
            reverse_complement_sequence: reverse_seq.clone(),
            unknown_sequence: vec![0; seq.len()],
            masks: vec![],
            gc_content: 0.5,
            sequence_length: seq.len(),
        };

        // Test that we can call the main exported functions
        // (We won't actually run them as they have complex signatures)

        // score_nodes is accessible
        let _ = score_nodes(&encoded_sequence, &mut nodes, &training, false, false);

        // raw_coding_score is accessible
        raw_coding_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        // calc_orf_gc is accessible
        calc_orf_gc(&seq, seq.len(), &mut nodes);

        // reset_node_scores is accessible
        reset_node_scores(&mut nodes);

        // sort_nodes_by_position is accessible
        sort_nodes_by_position(&mut nodes);

        // Function completion indicates success
    }

    #[test]
    fn test_reset_node_scores_functionality() {
        let mut nodes = vec![Node {
            position: NodePosition {
                index: 0,
                stop_value: 12,
                strand: Strand::Forward,
                codon_type: CodonType::Atg,
                is_edge: false,
            },
            scores: NodeScores {
                coding_score: 10.0,
                upstream_score: 5.0,
                ribosome_binding_score: 3.0,
                type_score: 2.0,
                total_score: 20.0,
                gc_content: 0.5,
                start_score: 1.0,
                gc_frame_scores: [0.4, 0.5, 0.6],
            },
            state: NodeState::default(),
            motif_info: NodeMotifInfo {
                ribosome_binding_sites: [0; 2],
                best_motif: Motif::default(),
            },
        }];

        reset_node_scores(&mut nodes);

        assert_eq!(nodes[0].scores.coding_score, 0.0);
        assert_eq!(nodes[0].scores.upstream_score, 0.0);
        assert_eq!(nodes[0].scores.ribosome_binding_score, 0.0);
        assert_eq!(nodes[0].scores.type_score, 0.0);
        assert_eq!(nodes[0].scores.total_score, 0.0);
    }

    #[test]
    fn test_sort_nodes_by_position_functionality() {
        let mut nodes = vec![
            Node {
                position: NodePosition {
                    index: 30,
                    stop_value: 40,
                    strand: Strand::Forward,
                    codon_type: CodonType::Atg,
                    is_edge: false,
                },
                scores: NodeScores::default(),
                state: NodeState::default(),
                motif_info: NodeMotifInfo {
                    ribosome_binding_sites: [0; 2],
                    best_motif: Motif::default(),
                },
            },
            Node {
                position: NodePosition {
                    index: 10,
                    stop_value: 20,
                    strand: Strand::Forward,
                    codon_type: CodonType::Atg,
                    is_edge: false,
                },
                scores: NodeScores::default(),
                state: NodeState::default(),
                motif_info: NodeMotifInfo {
                    ribosome_binding_sites: [0; 2],
                    best_motif: Motif::default(),
                },
            },
        ];

        sort_nodes_by_position(&mut nodes);

        assert!(nodes[0].position.index <= nodes[1].position.index);
        assert_eq!(nodes[0].position.index, 10);
        assert_eq!(nodes[1].position.index, 30);
    }

    #[test]
    fn test_node_module_integration() {
        // Test that the node module can perform basic operations
        let seq = vec![
            0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 3, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ];
        let reverse_seq = vec![3; seq.len()];
        let training = create_test_training();

        let mut nodes = vec![create_test_node()];

        // Test that we can perform basic scoring operations
        raw_coding_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);
        assert!(nodes[0].scores.coding_score.is_finite());

        // Test GC calculation
        calc_orf_gc(&seq, seq.len(), &mut nodes);

        reset_node_scores(&mut nodes);
        assert_eq!(nodes[0].scores.total_score, 0.0);

        sort_nodes_by_position(&mut nodes);

        // If we get here, all basic node operations work
        // Function completion indicates success
    }
}
