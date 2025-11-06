use bio::bio_types::strand::Strand;

use crate::{
    node::intergenic_mod,
    types::{CodonType, Node, Training},
};

/// Sometimes bad genes creep into the model due to the node distance constraint
/// in the dynamic programming routine. This routine does a sweep through
/// the genes and eliminates ones with negative scores.
pub fn eliminate_bad_genes(
    nodes: &mut [Node],
    dynamic_programming_start: Option<usize>,
    training: &Training,
) {
    let Some(dynamic_programming_start) = dynamic_programming_start else {
        return;
    };

    // Find the beginning of the traceback path
    let path_start = find_path_start(nodes, dynamic_programming_start);

    // First pass: Add intergenic modifiers to start scores
    add_intergenic_bonuses(nodes, path_start, training);

    // Second pass: Mark genes with negative scores for elimination
    mark_bad_genes_for_elimination(nodes, path_start);
}

/// Find the start of the traceback path by following traceback pointers
fn find_path_start(nodes: &[Node], start_idx: usize) -> usize {
    let mut path = start_idx;
    while let Some(next_path) = nodes[path].state.traceback {
        path = next_path;
    }
    path
}

/// Add intergenic bonuses to start scores along the forward path
fn add_intergenic_bonuses(nodes: &mut [Node], path_start: usize, training: &Training) {
    let mut path = path_start;

    while let Some(next_path) = nodes[path].state.trace_forward {
        // Check bounds to prevent index out of bounds
        if next_path >= nodes.len() {
            break;
        }

        let is_path_stop = nodes[path].position.codon_type == CodonType::Stop;

        // Apply intergenic bonuses based on connection type
        match (nodes[path].position.strand, is_path_stop) {
            (Strand::Forward, true) => {
                let intergenic_bonus = intergenic_mod(&nodes[path], &nodes[next_path], training);
                nodes[next_path].scores.start_score += intergenic_bonus;
            }
            (Strand::Reverse, false) => {
                let intergenic_bonus = intergenic_mod(&nodes[path], &nodes[next_path], training);
                nodes[path].scores.start_score += intergenic_bonus;
            }
            _ => {} // Other connection types don't need bonuses here
        }

        path = next_path;
    }
}

/// Mark genes with negative total scores for elimination
fn mark_bad_genes_for_elimination(nodes: &mut [Node], path_start: usize) {
    let mut path = path_start;

    while let Some(next_path) = nodes[path].state.trace_forward {
        if next_path >= nodes.len() {
            break;
        }

        // Check and eliminate bad genes based on their type and scores
        let is_path_stop = nodes[path].position.codon_type == CodonType::Stop;
        match (nodes[path].position.strand, is_path_stop) {
            (Strand::Forward, false) => {
                if should_eliminate_gene(&nodes[path]) {
                    eliminate_gene_pair(nodes, path, next_path);
                }
            }
            (Strand::Reverse, true) => {
                if should_eliminate_gene(&nodes[next_path]) {
                    eliminate_gene_pair(nodes, path, next_path);
                }
            }
            _ => {}
        }

        path = next_path;
    }
}

/// Check if a gene should be eliminated based on its scores
#[inline]
fn should_eliminate_gene(node: &Node) -> bool {
    node.scores.coding_score + node.scores.start_score < 0.0
}

/// Mark both nodes of a gene pair for elimination
#[inline]
fn eliminate_gene_pair(nodes: &mut [Node], idx1: usize, idx2: usize) {
    nodes[idx1].state.is_eliminated = true;
    nodes[idx2].state.is_eliminated = true;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{NodeMotifInfo, NodePosition, NodeScores, NodeState};

    /// Helper function to create a test node
    fn create_test_node(
        index: usize,
        strand: Strand,
        codon_type: CodonType,
        stop_value: isize,
        coding_score: f64,
        start_score: f64,
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
                coding_score,
                start_score,
                ..Default::default()
            },
            state: NodeState::default(),
            motif_info: NodeMotifInfo::default(),
        }
    }

    /// Helper function to create test training data
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
    fn test_eliminate_bad_genes_none_input() {
        let mut nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            10.0,
            5.0,
        )];
        let training = create_test_training();

        eliminate_bad_genes(&mut nodes, None, &training);

        // Should return early and not modify nodes
        assert!(!nodes[0].state.is_eliminated);
    }

    #[test]
    fn test_eliminate_bad_genes_single_node() {
        let mut nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Stop,
            50,
            10.0,
            5.0,
        )];
        let training = create_test_training();

        eliminate_bad_genes(&mut nodes, Some(0), &training);

        assert!(!nodes[0].state.is_eliminated);
    }

    #[test]
    fn test_eliminate_bad_genes_positive_scores() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, 10.0, 5.0),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, 8.0, 3.0),
        ];
        let training = create_test_training();

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);

        eliminate_bad_genes(&mut nodes, Some(0), &training);

        assert!(!nodes[0].state.is_eliminated);
        assert!(!nodes[1].state.is_eliminated);
    }

    #[test]
    fn test_eliminate_bad_genes_negative_forward_gene() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, -10.0, -5.0),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, 8.0, 3.0),
        ];
        let training = create_test_training();

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);

        eliminate_bad_genes(&mut nodes, Some(0), &training);

        assert!(nodes[0].state.is_eliminated);
        assert!(nodes[1].state.is_eliminated);
    }

    #[test]
    fn test_eliminate_bad_genes_negative_reverse_gene() {
        let mut nodes = vec![
            create_test_node(200, Strand::Reverse, CodonType::Stop, 100, 8.0, 3.0),
            create_test_node(100, Strand::Reverse, CodonType::Atg, 200, -10.0, -5.0),
        ];
        let training = create_test_training();

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);

        eliminate_bad_genes(&mut nodes, Some(0), &training);

        assert!(nodes[0].state.is_eliminated);
        assert!(nodes[1].state.is_eliminated);
    }

    #[test]
    fn test_find_path_start_single_node() {
        let nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            10.0,
            5.0,
        )];

        let start = find_path_start(&nodes, 0);
        assert_eq!(start, 0);
    }

    #[test]
    fn test_find_path_start_chain() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, 10.0, 5.0),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, 8.0, 3.0),
            create_test_node(300, Strand::Forward, CodonType::Atg, 0, 7.0, 4.0),
        ];

        // Set up traceback chain: 2 -> 1 -> 0
        nodes[2].state.traceback = Some(1);
        nodes[1].state.traceback = Some(0);

        let start = find_path_start(&nodes, 2);
        assert_eq!(start, 0); // Should find the beginning of the chain
    }

    #[test]
    fn test_add_intergenic_bonuses_forward_gene() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Stop, 50, 8.0, 3.0),
            create_test_node(130, Strand::Forward, CodonType::Atg, 0, 10.0, 5.0), // Distance 30 < OPERON_DISTANCE (60)
        ];
        let training = create_test_training();

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);
        let original_start_score = nodes[1].scores.start_score;

        add_intergenic_bonuses(&mut nodes, 0, &training);

        assert_ne!(nodes[1].scores.start_score, original_start_score);
    }

    #[test]
    fn test_add_intergenic_bonuses_reverse_gene() {
        let mut nodes = vec![
            create_test_node(130, Strand::Reverse, CodonType::Atg, 100, 10.0, 5.0),
            create_test_node(170, Strand::Reverse, CodonType::Stop, 200, 8.0, 3.0), // Distance 40, no overlap
        ];
        let training = create_test_training();

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);
        let original_start_score = nodes[0].scores.start_score;

        add_intergenic_bonuses(&mut nodes, 0, &training);

        assert_ne!(nodes[0].scores.start_score, original_start_score);
    }

    #[test]
    fn test_mark_bad_genes_for_elimination_forward_gene() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, -10.0, -5.0), // Negative
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, 8.0, 3.0),
        ];

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);

        mark_bad_genes_for_elimination(&mut nodes, 0);

        assert!(nodes[0].state.is_eliminated);
        assert!(nodes[1].state.is_eliminated);
    }

    #[test]
    fn test_mark_bad_genes_for_elimination_reverse_gene() {
        let mut nodes = vec![
            create_test_node(200, Strand::Reverse, CodonType::Stop, 100, 8.0, 3.0),
            create_test_node(100, Strand::Reverse, CodonType::Atg, 200, -10.0, -5.0), // Negative
        ];

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);

        mark_bad_genes_for_elimination(&mut nodes, 0);

        assert!(nodes[0].state.is_eliminated);
        assert!(nodes[1].state.is_eliminated);
    }

    #[test]
    fn test_mark_bad_genes_no_elimination_needed() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, 10.0, 5.0), // Positive
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, 8.0, 3.0),
        ];

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);

        mark_bad_genes_for_elimination(&mut nodes, 0);

        assert!(!nodes[0].state.is_eliminated);
        assert!(!nodes[1].state.is_eliminated);
    }

    #[test]
    fn test_should_eliminate_gene_positive() {
        let node = create_test_node(100, Strand::Forward, CodonType::Atg, 0, 10.0, 5.0);
        assert!(!should_eliminate_gene(&node));
    }

    #[test]
    fn test_should_eliminate_gene_negative() {
        let node = create_test_node(100, Strand::Forward, CodonType::Atg, 0, -10.0, -5.0);
        assert!(should_eliminate_gene(&node));
    }

    #[test]
    fn test_should_eliminate_gene_zero() {
        let node = create_test_node(100, Strand::Forward, CodonType::Atg, 0, 0.0, 0.0);
        assert!(!should_eliminate_gene(&node)); // Zero is not negative
    }

    #[test]
    fn test_should_eliminate_gene_mixed_scores() {
        let node1 = create_test_node(100, Strand::Forward, CodonType::Atg, 0, 10.0, -15.0);
        let node2 = create_test_node(200, Strand::Forward, CodonType::Atg, 0, -5.0, 10.0);

        assert!(should_eliminate_gene(&node1));
        assert!(!should_eliminate_gene(&node2));
    }

    #[test]
    fn test_eliminate_gene_pair() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, 10.0, 5.0),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, 8.0, 3.0),
        ];

        eliminate_gene_pair(&mut nodes, 0, 1);

        assert!(nodes[0].state.is_eliminated);
        assert!(nodes[1].state.is_eliminated);
    }

    #[test]
    fn test_eliminate_bad_genes_complex_path() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, 10.0, 5.0),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, 8.0, 3.0),
            create_test_node(300, Strand::Forward, CodonType::Atg, 0, -10.0, -5.0),
            create_test_node(400, Strand::Forward, CodonType::Stop, 300, 8.0, 3.0),
            create_test_node(500, Strand::Forward, CodonType::Atg, 0, 12.0, 6.0),
            create_test_node(600, Strand::Forward, CodonType::Stop, 500, 9.0, 4.0),
        ];
        let training = create_test_training();

        // Set up forward path
        nodes[0].state.trace_forward = Some(1);
        nodes[1].state.trace_forward = Some(2);
        nodes[2].state.trace_forward = Some(3);
        nodes[3].state.trace_forward = Some(4);
        nodes[4].state.trace_forward = Some(5);

        eliminate_bad_genes(&mut nodes, Some(0), &training);

        assert!(!nodes[0].state.is_eliminated);
        assert!(!nodes[1].state.is_eliminated);
        assert!(nodes[2].state.is_eliminated);
        assert!(nodes[3].state.is_eliminated);
        assert!(!nodes[4].state.is_eliminated);
        assert!(!nodes[5].state.is_eliminated);
    }

    #[test]
    fn test_eliminate_bad_genes_out_of_bounds_protection() {
        let mut nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            10.0,
            5.0,
        )];
        let training = create_test_training();

        // Set up invalid forward path (out of bounds)
        nodes[0].state.trace_forward = Some(10); // Index 10 doesn't exist

        eliminate_bad_genes(&mut nodes, Some(0), &training);

        assert!(!nodes[0].state.is_eliminated);
    }
}
