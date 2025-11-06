use bio::bio_types::strand::Strand;

use crate::{
    algorithms::dynamic_programming::{
        overlap_resolution::{resolve_simple_overlaps, resolve_triple_overlaps},
        scoring::handlers::score_connection,
    },
    constants::MAX_NODE_DIST,
    types::{CodonType, Node, Training},
};

/// Basic dynamic programming routine for predicting genes.
///
/// The `flag` variable is set to false for the initial dynamic programming routine
/// based solely on GC frame plot (used to construct a training set). If the flag
/// is set to true, the routine does the final dynamic programming based on coding,
/// RBS scores, etc.
///
/// * `is_final_pass = false` – first pass, GC‑frame only (training set)
/// * `is_final_pass = true`  – second pass, full scoring (final prediction)
pub fn predict_genes(
    nodes: &mut [Node],
    training: &Training,
    is_final_pass: bool,
) -> Option<usize> {
    let node_count = nodes.len();
    if node_count == 0 {
        return None;
    }

    initialize_nodes(nodes);
    compute_dynamic_programming_scores(nodes, training, is_final_pass);
    let max_index = find_maximum_terminal_node(nodes)?;

    resolve_triple_overlaps(nodes, max_index);
    resolve_simple_overlaps(nodes, max_index);
    setup_forward_pointers(nodes, max_index);

    Some(max_index)
}

/// Initialize all nodes for dynamic programming
fn initialize_nodes(nodes: &mut [Node]) {
    for node in nodes.iter_mut() {
        node.scores.total_score = 0.0;
        node.state.traceback = None;
        node.state.trace_forward = None;
    }
}

/// Perform the main dynamic programming computation
fn compute_dynamic_programming_scores(
    nodes: &mut [Node],
    training: &Training,
    is_final_pass: bool,
) {
    for i in 0..nodes.len() {
        let min = calculate_distance_constraints(nodes, i);
        if min >= i {
            continue;
        }

        for j in min..i {
            score_connection(nodes, j, i, training, is_final_pass);
        }
    }
}

#[inline]
fn calculate_distance_constraints(nodes: &[Node], i: usize) -> usize {
    let mut min = i.saturating_sub(MAX_NODE_DIST);

    let node = &nodes[i];

    if node.position.strand == Strand::Reverse
        && node.position.codon_type != CodonType::Stop
        && nodes[min].position.index as isize >= node.position.stop_value
    {
        while min > 0 && node.position.index as isize != node.position.stop_value {
            min -= 1;
        }
    }
    if node.position.strand == Strand::Forward
        && node.position.codon_type == CodonType::Stop
        && nodes[min].position.index as isize >= node.position.stop_value
    {
        while min > 0 && node.position.index as isize != node.position.stop_value {
            min -= 1;
        }
    }

    min.saturating_sub(MAX_NODE_DIST)
}

fn is_terminal(node: &Node) -> bool {
    (node.position.strand == Strand::Forward && node.position.codon_type == CodonType::Stop)
        || (node.position.strand == Strand::Reverse && node.position.codon_type != CodonType::Stop)
}

fn find_maximum_terminal_node(nodes: &[Node]) -> Option<usize> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| is_terminal(n))
        .max_by(|(_, a), (_, b)| {
            a.scores
                .total_score
                .partial_cmp(&b.scores.total_score)
                .unwrap()
        })
        .map(|(i, _)| i)
}

fn setup_forward_pointers(nodes: &mut [Node], start_index: usize) {
    let mut path = start_index;

    while let Some(traceback_idx) = nodes[path].state.traceback {
        nodes[traceback_idx].state.trace_forward = Some(path);
        path = traceback_idx;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{NodeMotifInfo, NodePosition, NodeScores, NodeState};

    fn create_test_node(
        index: usize,
        strand: Strand,
        codon_type: CodonType,
        stop_value: isize,
        is_edge: bool,
    ) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value,
                is_edge,
            },
            scores: NodeScores::default(),
            state: NodeState::default(),
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
    fn test_predict_genes_empty_nodes() {
        let mut nodes = vec![];
        let training = create_test_training();

        let result = predict_genes(&mut nodes, &training, true);
        assert!(result.is_none());
    }

    #[test]
    fn test_predict_genes_single_node() {
        let mut nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            false,
        )];
        let training = create_test_training();

        let result = predict_genes(&mut nodes, &training, true);
        assert!(result.is_none());
    }

    #[test]
    fn test_predict_genes_simple_forward_gene() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(400, Strand::Forward, CodonType::Stop, 403, false),
        ];
        let training = create_test_training();

        nodes[0].scores.total_score = 5.0;
        nodes[1].scores.total_score = 10.0;

        let result = predict_genes(&mut nodes, &training, false);
        assert!(result.is_some());
        if let Some(max_idx) = result {
            assert_eq!(max_idx, 1); // Should be the stop node
        }
    }

    #[test]
    fn test_predict_genes_reverse_gene() {
        let mut nodes = vec![
            create_test_node(400, Strand::Reverse, CodonType::Stop, 0, false),
            create_test_node(100, Strand::Reverse, CodonType::Atg, 403, false),
        ];
        let training = create_test_training();

        // Set up scores to ensure connection
        nodes[0].scores.total_score = 5.0;
        nodes[1].scores.total_score = 10.0;

        let result = predict_genes(&mut nodes, &training, false);
        assert!(result.is_some());
        if let Some(max_idx) = result {
            assert_eq!(max_idx, 1); // Should be the reverse start node
        }
    }

    #[test]
    fn test_initialize_nodes() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, false),
        ];

        // Set some initial values
        nodes[0].scores.total_score = 5.0;
        nodes[0].state.traceback = Some(1);
        nodes[1].state.trace_forward = Some(0);

        initialize_nodes(&mut nodes);

        // Check that initialization cleared the values
        assert_eq!(nodes[0].scores.total_score, 0.0);
        assert!(nodes[0].state.traceback.is_none());
        assert!(nodes[0].state.trace_forward.is_none());
        assert_eq!(nodes[1].scores.total_score, 0.0);
        assert!(nodes[1].state.traceback.is_none());
        assert!(nodes[1].state.trace_forward.is_none());
    }

    #[test]
    fn test_calculate_distance_constraints_small_index() {
        let nodes = vec![
            create_test_node(10, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(20, Strand::Forward, CodonType::Stop, 10, false),
        ];

        let min = calculate_distance_constraints(&nodes, 1);
        assert_eq!(min, 0); // Small index should return 0
    }

    #[test]
    fn test_calculate_distance_constraints_large_index() {
        let mut nodes = Vec::new();
        for i in 0..1000 {
            nodes.push(create_test_node(
                i * 10,
                Strand::Forward,
                CodonType::Atg,
                0,
                false,
            ));
        }

        let min = calculate_distance_constraints(&nodes, 600); // Use index > MAX_NODE_DIST
        assert!(min <= 600); // Should not exceed the index
        // Test that it returns a reasonable value (since usize is always >= 0)
        assert!(min < nodes.len());
    }

    #[test]
    fn test_find_maximum_terminal_node_forward_stops() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, false),
            create_test_node(300, Strand::Forward, CodonType::Stop, 150, false),
        ];

        // Set different scores for stop nodes
        nodes[1].scores.total_score = 5.0;
        nodes[2].scores.total_score = 10.0;

        let max_idx = find_maximum_terminal_node(&nodes);
        assert_eq!(max_idx, Some(2)); // Should pick the highest scoring terminal
    }

    #[test]
    fn test_find_maximum_terminal_node_reverse_starts() {
        let mut nodes = vec![
            create_test_node(300, Strand::Reverse, CodonType::Stop, 100, false),
            create_test_node(200, Strand::Reverse, CodonType::Atg, 300, false),
            create_test_node(100, Strand::Reverse, CodonType::Gtg, 300, false),
        ];

        // Set different scores for start nodes
        nodes[1].scores.total_score = 8.0;
        nodes[2].scores.total_score = 12.0;

        let max_idx = find_maximum_terminal_node(&nodes);
        assert_eq!(max_idx, Some(2)); // Should pick the highest scoring reverse start
    }

    #[test]
    fn test_find_maximum_terminal_node_no_terminals() {
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(200, Strand::Reverse, CodonType::Stop, 100, false),
        ];

        let max_idx = find_maximum_terminal_node(&nodes);
        assert!(max_idx.is_none());
    }

    #[test]
    fn test_find_maximum_terminal_node_negative_scores() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Stop, 50, false),
            create_test_node(200, Strand::Forward, CodonType::Stop, 150, false),
        ];

        nodes[0].scores.total_score = -5.0;
        nodes[1].scores.total_score = -10.0;

        let max_idx = find_maximum_terminal_node(&nodes);
        assert_eq!(max_idx, Some(0)); // Should pick the least negative
    }

    #[test]
    fn test_setup_forward_pointers_simple_chain() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, false),
            create_test_node(300, Strand::Forward, CodonType::Atg, 0, false),
        ];

        // Set up traceback chain: 2 -> 1 -> 0
        nodes[2].state.traceback = Some(1);
        nodes[1].state.traceback = Some(0);

        setup_forward_pointers(&mut nodes, 2);

        // Check forward pointers were set correctly
        assert_eq!(nodes[0].state.trace_forward, Some(1));
        assert_eq!(nodes[1].state.trace_forward, Some(2));
        assert!(nodes[2].state.trace_forward.is_none());
    }

    #[test]
    fn test_setup_forward_pointers_single_node() {
        let mut nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Stop,
            50,
            false,
        )];

        setup_forward_pointers(&mut nodes, 0);

        assert!(nodes[0].state.trace_forward.is_none());
    }

    #[test]
    fn test_setup_forward_pointers_complex_chain() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, false),
            create_test_node(300, Strand::Reverse, CodonType::Stop, 150, false),
            create_test_node(150, Strand::Reverse, CodonType::Atg, 300, false),
            create_test_node(400, Strand::Forward, CodonType::Stop, 350, false),
        ];

        // Set up traceback chain: 4 -> 3 -> 2 -> 1 -> 0
        nodes[4].state.traceback = Some(3);
        nodes[3].state.traceback = Some(2);
        nodes[2].state.traceback = Some(1);
        nodes[1].state.traceback = Some(0);

        setup_forward_pointers(&mut nodes, 4);

        // Check all forward pointers
        assert_eq!(nodes[0].state.trace_forward, Some(1));
        assert_eq!(nodes[1].state.trace_forward, Some(2));
        assert_eq!(nodes[2].state.trace_forward, Some(3));
        assert_eq!(nodes[3].state.trace_forward, Some(4));
        assert!(nodes[4].state.trace_forward.is_none());
    }

    #[test]
    fn test_compute_dynamic_programming_scores_empty() {
        let mut nodes = vec![];
        let training = create_test_training();

        compute_dynamic_programming_scores(&mut nodes, &training, true);
    }

    #[test]
    fn test_compute_dynamic_programming_scores_single_node() {
        let mut nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            false,
        )];
        let training = create_test_training();

        compute_dynamic_programming_scores(&mut nodes, &training, true);

        assert_eq!(nodes[0].scores.total_score, 0.0);
        assert!(nodes[0].state.traceback.is_none());
    }

    #[test]
    fn test_predict_genes_with_forward_pointers() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(200, Strand::Forward, CodonType::Stop, 103, false),
            create_test_node(300, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(400, Strand::Forward, CodonType::Stop, 303, false),
        ];
        let training = create_test_training();

        // Set up scores to make the last stop node the best terminal
        nodes[1].scores.total_score = 5.0;
        nodes[3].scores.total_score = 15.0;

        let result = predict_genes(&mut nodes, &training, false);
        assert!(result.is_some());

        // Just verify we get a result (forward pointers are complex to test in isolation)
        if let Some(max_idx) = result {
            assert_eq!(max_idx, 3); // Should pick the highest scoring terminal
        }
    }
}
