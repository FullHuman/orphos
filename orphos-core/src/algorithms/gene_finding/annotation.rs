use bio::bio_types::strand::Strand;

use crate::{
    constants::{MAX_CONFIDENCE_SCORE, RBS_DESCRIPTIONS},
    sequence::mer_text,
    types::{Gene, GeneAnnotation, GeneScore, Node, StartType, Training},
};

/// Record gene annotation and score.
///
/// This function creates a `GeneAnnotation` and `GeneScore` for a given gene based on the
/// nodes and training data. It calculates the start type, partial flags, RBS motif,
/// and confidence score. The function returns a tuple containing the annotation and score.
pub fn record_gene_annotation_and_score(
    gene: &Gene,
    nodes: &[Node],
    training: &Training,
    sequence_number: usize,
    index: usize,
) -> (GeneAnnotation, GeneScore) {
    let start_node = &nodes[gene.coordinates.start_index];
    let stop_node = &nodes[gene.coordinates.stop_index];

    let id = format!("{}_{}", sequence_number, index + 1);
    let partial_left = is_partial_left(start_node, stop_node);
    let partial_right = is_partial_right(start_node, stop_node);
    let start_type = if start_node.position.is_edge {
        StartType::Edge
    } else {
        StartType::from(start_node.position.codon_type)
    };
    let mut annotation = GeneAnnotation::new(
        id,
        partial_left,
        partial_right,
        start_type,
        start_node.scores.gc_content,
    );

    let rbs1_score = training.rbs_weights[start_node.motif_info.ribosome_binding_sites[0]]
        * training.start_weight_factor;
    let rbs2_score = training.rbs_weights[start_node.motif_info.ribosome_binding_sites[1]]
        * training.start_weight_factor;
    if training.uses_shine_dalgarno {
        let (motif, spacer) = rbs_motif_and_spacer(start_node, rbs1_score, rbs2_score);
        annotation = annotation.with_rbs(motif.into(), spacer.into());
    } else {
        // Non-SD mode - use motif information
        if training.no_motif_weight > -0.5
            && ((rbs1_score > rbs2_score
                && rbs1_score
                    > start_node.motif_info.best_motif.score * training.start_weight_factor)
                || (rbs2_score >= rbs1_score
                    && rbs2_score
                        > start_node.motif_info.best_motif.score * training.start_weight_factor))
        {
            let (motif, spacer) = rbs_motif_and_spacer(start_node, rbs1_score, rbs2_score);
            annotation = annotation.with_rbs(motif.into(), spacer.into());
        } else if start_node.motif_info.best_motif.length == 0 {
        } else {
            let motif_text = mer_text(
                start_node.motif_info.best_motif.length,
                start_node.motif_info.best_motif.index,
            );
            annotation = annotation.with_rbs(
                motif_text,
                format!("{}bp", start_node.motif_info.best_motif.spacer),
            );
        }
    }

    // Calculate confidence and build score data
    let total_score = start_node.scores.coding_score + start_node.scores.start_score;
    let confidence = calculate_confidence(total_score, training.start_weight_factor);
    let gene_score = GeneScore {
        confidence,
        total_score,
        coding_score: start_node.scores.coding_score,
        start_score: start_node.scores.start_score,
        ribosome_binding_score: start_node.scores.ribosome_binding_score,
        upstream_score: start_node.scores.upstream_score,
        type_score: start_node.scores.type_score,
    };

    (annotation, gene_score)
}

/// Get the rbs motif and the rbs spacer
/// This function retrieves the RBS motif and spacer length based on the RBS index.
/// It uses a predefined list of RBS descriptions to find the corresponding motif and spacer.
fn rbs_motif_and_spacer(
    start_node: &Node,
    rbs1_score: f64,
    rbs2_score: f64,
) -> (&'static str, &'static str) {
    let rbs_idx = if rbs1_score > rbs2_score {
        start_node.motif_info.ribosome_binding_sites[0]
    } else {
        start_node.motif_info.ribosome_binding_sites[1]
    };

    RBS_DESCRIPTIONS[rbs_idx.min(RBS_DESCRIPTIONS.len() - 1)]
}

fn is_partial_left(start_node: &Node, stop_node: &Node) -> bool {
    (start_node.position.is_edge && start_node.position.strand == Strand::Forward)
        || (stop_node.position.is_edge && start_node.position.strand == Strand::Reverse)
}

fn is_partial_right(start_node: &Node, stop_node: &Node) -> bool {
    (stop_node.position.is_edge && start_node.position.strand == Strand::Forward)
        || (start_node.position.is_edge && start_node.position.strand == Strand::Reverse)
}

/// Convert score to percent confidence
fn calculate_confidence(score: f64, start_weight: f64) -> f64 {
    let normalized_score = score / start_weight;

    if normalized_score < 41.0 {
        let exp_score = normalized_score.exp();
        let confidence = (exp_score / (exp_score + 1.0)) * 100.0;
        confidence.max(50.0)
    } else {
        MAX_CONFIDENCE_SCORE
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{
        CodonType, GeneCoordinates, Motif, NodeMotifInfo, NodePosition, NodeScores, NodeState,
    };

    /// Helper function to create a test node
    #[allow(clippy::too_many_arguments)]
    fn create_test_node(
        index: usize,
        strand: Strand,
        codon_type: CodonType,
        is_edge: bool,
        gc_content: f64,
        coding_score: f64,
        start_score: f64,
        ribosome_binding_score: f64,
        upstream_score: f64,
        type_score: f64,
        rbs_sites: [usize; 2],
        best_motif: Motif,
    ) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value: 0,
                is_edge,
            },
            scores: NodeScores {
                gc_content,
                coding_score,
                start_score,
                ribosome_binding_score,
                type_score,
                upstream_score,
                total_score: coding_score + start_score,
                ..Default::default()
            },
            state: NodeState::default(),
            motif_info: NodeMotifInfo {
                ribosome_binding_sites: rbs_sites,
                best_motif,
            },
        }
    }

    /// Helper function to create a test gene
    fn create_test_gene(start_index: usize, stop_index: usize, strand: Strand) -> Gene {
        Gene {
            coordinates: GeneCoordinates {
                begin: 100,
                end: 300,
                strand,
                start_index,
                stop_index,
            },
            score: GeneScore::default(),
            annotation: GeneAnnotation::default(),
        }
    }

    /// Helper function to create test training data
    fn create_test_training(
        uses_shine_dalgarno: bool,
        start_weight_factor: f64,
        no_motif_weight: f64,
    ) -> Training {
        let mut training = Training {
            uses_shine_dalgarno,
            start_weight_factor,
            no_motif_weight,
            ..Training::default()
        };

        // Set some non-zero RBS weights for testing
        training.rbs_weights[1] = 2.0;
        training.rbs_weights[5] = 1.5;
        training.rbs_weights[10] = 3.0;

        training
    }

    #[test]
    fn test_calculate_confidence_low_score() {
        let confidence = calculate_confidence(10.0, 4.35);
        assert!(confidence >= 50.0);
        assert!(confidence < 99.99);
    }

    #[test]
    fn test_calculate_confidence_high_score() {
        let confidence = calculate_confidence(200.0, 4.35);
        assert_eq!(confidence, 99.99);
    }

    #[test]
    fn test_calculate_confidence_boundary_score() {
        // Test at boundary (41.0 * start_weight)
        let confidence = calculate_confidence(41.0 * 4.35, 4.35);
        assert_eq!(confidence, 99.99);
    }

    #[test]
    fn test_calculate_confidence_zero_score() {
        let confidence = calculate_confidence(0.0, 4.35);
        assert_eq!(confidence, 50.0); // Should be clamped to minimum 50.0
    }

    #[test]
    fn test_calculate_confidence_negative_score() {
        let confidence = calculate_confidence(-10.0, 4.35);
        assert_eq!(confidence, 50.0); // Should be clamped to minimum 50.0
    }

    #[test]
    fn test_is_partial_left_forward_edge_start() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            true,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );

        assert!(is_partial_left(&start_node, &stop_node));
    }

    #[test]
    fn test_is_partial_left_reverse_edge_stop() {
        let start_node = create_test_node(
            300,
            Strand::Reverse,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            100,
            Strand::Reverse,
            CodonType::Stop,
            true,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );

        assert!(is_partial_left(&start_node, &stop_node));
    }

    #[test]
    fn test_is_partial_left_false() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );

        assert!(!is_partial_left(&start_node, &stop_node));
    }

    #[test]
    fn test_is_partial_right_forward_edge_stop() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            true,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );

        assert!(is_partial_right(&start_node, &stop_node));
    }

    #[test]
    fn test_is_partial_right_reverse_edge_start() {
        let start_node = create_test_node(
            300,
            Strand::Reverse,
            CodonType::Atg,
            true,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            100,
            Strand::Reverse,
            CodonType::Stop,
            false,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );

        assert!(is_partial_right(&start_node, &stop_node));
    }

    #[test]
    fn test_is_partial_right_false() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );

        assert!(!is_partial_right(&start_node, &stop_node));
    }

    #[test]
    fn test_rbs_motif_and_spacer_first_higher() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );

        let (motif, spacer) = rbs_motif_and_spacer(&start_node, 3.0, 2.0);
        assert_eq!(motif, "GGA/GAG/AGG");
        assert_eq!(spacer, "3-4bp");
    }

    #[test]
    fn test_rbs_motif_and_spacer_second_higher() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );

        let (motif, spacer) = rbs_motif_and_spacer(&start_node, 1.0, 2.0);
        assert_eq!(motif, "AGxAG");
        assert_eq!(spacer, "3-4bp");
    }

    #[test]
    fn test_rbs_motif_and_spacer_equal_scores() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );

        let (motif, spacer) = rbs_motif_and_spacer(&start_node, 2.0, 2.0);
        assert_eq!(motif, "AGxAG"); // Second one wins in case of tie
        assert_eq!(spacer, "3-4bp");
    }

    #[test]
    fn test_rbs_motif_and_spacer_boundary_index() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [27, 100],
            Motif::default(), // Test with boundary indices
        );

        let (motif, spacer) = rbs_motif_and_spacer(&start_node, 3.0, 2.0);
        assert_eq!(motif, "AGGAGG"); // Index 27 (last valid)
        assert_eq!(spacer, "5-10bp");
    }

    #[test]
    fn test_record_gene_annotation_and_score_basic_forward() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.45,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.45,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );
        let nodes = vec![start_node, stop_node];
        let gene = create_test_gene(0, 1, Strand::Forward);
        let training = create_test_training(true, 4.35, 0.0);

        let (annotation, score) = record_gene_annotation_and_score(&gene, &nodes, &training, 1, 0);

        assert_eq!(annotation.identifier, "1_1");
        assert!(!annotation.is_partial_left);
        assert!(!annotation.is_partial_right);
        assert!(matches!(annotation.start_type, StartType::Atg));
        assert_eq!(annotation.gc_content, 0.45);
        assert!(annotation.ribosome_binding_motif.is_some());
        assert!(annotation.ribosome_binding_spacer.is_some());

        assert_eq!(score.coding_score, 10.0);
        assert_eq!(score.start_score, 5.0);
        assert_eq!(score.total_score, 15.0);
        assert!(score.confidence >= 50.0);
    }

    #[test]
    fn test_record_gene_annotation_and_score_edge_gene() {
        let start_node = create_test_node(
            0,
            Strand::Forward,
            CodonType::Atg,
            true,
            0.6,
            12.0,
            6.0,
            3.0,
            2.0,
            4.0,
            [10, 15],
            Motif::default(),
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.6,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );
        let nodes = vec![start_node, stop_node];
        let gene = create_test_gene(0, 1, Strand::Forward);
        let training = create_test_training(true, 4.35, 0.0);

        let (annotation, _score) = record_gene_annotation_and_score(&gene, &nodes, &training, 2, 5);

        assert_eq!(annotation.identifier, "2_6");
        assert!(annotation.is_partial_left); // Edge start for forward gene
        assert!(!annotation.is_partial_right);
        assert!(matches!(annotation.start_type, StartType::Edge));
        assert_eq!(annotation.gc_content, 0.6);
    }

    #[test]
    fn test_record_gene_annotation_and_score_reverse_gene() {
        let start_node = create_test_node(
            300,
            Strand::Reverse,
            CodonType::Gtg,
            false,
            0.55,
            9.0,
            4.0,
            1.5,
            0.8,
            2.5,
            [2, 7],
            Motif::default(),
        );
        let stop_node = create_test_node(
            100,
            Strand::Reverse,
            CodonType::Stop,
            false,
            0.55,
            7.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );
        let nodes = vec![start_node, stop_node];
        let gene = create_test_gene(0, 1, Strand::Reverse);
        let training = create_test_training(true, 4.35, 0.0);

        let (annotation, _score) = record_gene_annotation_and_score(&gene, &nodes, &training, 1, 2);

        assert_eq!(annotation.identifier, "1_3");
        assert!(!annotation.is_partial_left);
        assert!(!annotation.is_partial_right);
        assert!(matches!(annotation.start_type, StartType::Gtg));
        assert_eq!(annotation.gc_content, 0.55);
    }

    #[test]
    fn test_record_gene_annotation_and_score_non_sd_mode() {
        let motif = Motif {
            index: 123,
            length: 4,
            space_index: 0,
            spacer: 8,
            score: 2.5,
        };
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Ttg,
            false,
            0.5,
            8.0,
            3.0,
            1.0,
            0.5,
            1.5,
            [1, 5],
            motif,
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.5,
            6.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );
        let nodes = vec![start_node, stop_node];
        let gene = create_test_gene(0, 1, Strand::Forward);
        let training = create_test_training(false, 4.35, 1.0); // Non-SD mode

        let (annotation, _score) = record_gene_annotation_and_score(&gene, &nodes, &training, 1, 0);

        assert_eq!(annotation.identifier, "1_1");
        assert!(matches!(annotation.start_type, StartType::Ttg));
        // In non-SD mode, might use motif information differently
        assert_eq!(annotation.gc_content, 0.5);
    }

    #[test]
    fn test_record_gene_annotation_and_score_non_sd_with_rbs() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );
        let nodes = vec![start_node, stop_node];
        let gene = create_test_gene(0, 1, Strand::Forward);
        let mut training = create_test_training(false, 4.35, 1.0); // Non-SD mode
        training.rbs_weights[1] = 5.0; // Make RBS score high enough

        let (annotation, _score) = record_gene_annotation_and_score(&gene, &nodes, &training, 1, 0);

        // Should use RBS since RBS score is higher than motif score
        assert!(annotation.ribosome_binding_motif.is_some());
        assert!(annotation.ribosome_binding_spacer.is_some());
    }

    #[test]
    fn test_record_gene_annotation_and_score_non_sd_with_motif() {
        let motif = Motif {
            index: 123,
            length: 4,
            space_index: 0,
            spacer: 8,
            score: 10.0,
        };
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            motif,
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );
        let nodes = vec![start_node, stop_node];
        let gene = create_test_gene(0, 1, Strand::Forward);
        let training = create_test_training(false, 4.35, 1.0); // Non-SD mode

        let (annotation, _score) = record_gene_annotation_and_score(&gene, &nodes, &training, 1, 0);

        // Should use motif since motif score is higher than RBS
        assert!(annotation.ribosome_binding_motif.is_some());
        assert_eq!(annotation.ribosome_binding_spacer, Some("8bp".to_string()));
    }

    #[test]
    fn test_record_gene_annotation_and_score_non_sd_no_motif() {
        let motif = Motif {
            index: 0,
            length: 0,
            space_index: 0,
            spacer: 0,
            score: 0.0,
        };
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            10.0,
            5.0,
            2.0,
            1.0,
            3.0,
            [0, 0],
            motif, // Low RBS indices
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );
        let nodes = vec![start_node, stop_node];
        let gene = create_test_gene(0, 1, Strand::Forward);
        let training = create_test_training(false, 4.35, 1.0); // Non-SD mode

        let (annotation, _score) = record_gene_annotation_and_score(&gene, &nodes, &training, 1, 0);

        // Should not set RBS/motif information
        assert!(annotation.ribosome_binding_motif.is_none());
        assert!(annotation.ribosome_binding_spacer.is_none());
    }

    #[test]
    fn test_record_gene_annotation_and_score_high_confidence() {
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            false,
            0.5,
            100.0,
            80.0,
            2.0,
            1.0,
            3.0,
            [1, 5],
            Motif::default(),
        );
        let stop_node = create_test_node(
            300,
            Strand::Forward,
            CodonType::Stop,
            false,
            0.5,
            8.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0, 0],
            Motif::default(),
        );
        let nodes = vec![start_node, stop_node];
        let gene = create_test_gene(0, 1, Strand::Forward);
        let training = create_test_training(true, 4.35, 0.0);

        let (_annotation, score) = record_gene_annotation_and_score(&gene, &nodes, &training, 1, 0);

        assert_eq!(score.confidence, 99.99); // High score should give max confidence
        assert_eq!(score.total_score, 180.0); // coding + start
        assert_eq!(score.coding_score, 100.0);
        assert_eq!(score.start_score, 80.0);
    }

    #[test]
    fn test_record_gene_annotation_and_score_all_partial_combinations() {
        // Test all combinations of partial flags
        let test_cases = vec![
            (true, false, false, false),
            (false, true, false, false),
            (false, false, true, false),
            (false, false, false, true),
        ];

        for (start_edge, stop_edge, is_reverse_start, is_reverse_stop) in test_cases {
            let (start_strand, stop_strand) = if is_reverse_start || is_reverse_stop {
                (Strand::Reverse, Strand::Reverse)
            } else {
                (Strand::Forward, Strand::Forward)
            };

            let start_node = create_test_node(
                if is_reverse_stop { 300 } else { 100 },
                start_strand,
                CodonType::Atg,
                start_edge,
                0.5,
                10.0,
                5.0,
                2.0,
                1.0,
                3.0,
                [1, 5],
                Motif::default(),
            );
            let stop_node = create_test_node(
                if is_reverse_stop { 100 } else { 300 },
                stop_strand,
                CodonType::Stop,
                stop_edge,
                0.5,
                8.0,
                0.0,
                0.0,
                0.0,
                0.0,
                [0, 0],
                Motif::default(),
            );
            let nodes = vec![start_node, stop_node];
            let gene = create_test_gene(0, 1, start_strand);
            let training = create_test_training(true, 4.35, 0.0);

            let (annotation, _score) =
                record_gene_annotation_and_score(&gene, &nodes, &training, 1, 0);

            let expected_partial_left = (start_edge && start_strand == Strand::Forward)
                || (stop_edge && start_strand == Strand::Reverse);
            let expected_partial_right = (stop_edge && start_strand == Strand::Forward)
                || (start_edge && start_strand == Strand::Reverse);

            assert_eq!(annotation.is_partial_left, expected_partial_left);
            assert_eq!(annotation.is_partial_right, expected_partial_right);
        }
    }
}
