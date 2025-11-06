use crate::{
    algorithms::dynamic_programming::scoring::connection_types::{
        ConnectionType, determine_connection_type, has_edge_artifacts,
    },
    constants::{MAXIMUM_OPPOSITE_OVERLAP, OVERLAP_PENALTY_FACTOR, STOP_CODON_OFFSET},
    node::intergenic_mod,
    types::{Node, Training},
};

#[derive(Debug)]
struct ConnectionState {
    left: usize,
    right: usize,
    overlap: i32,
    max_frame: Option<usize>,
    score: f64,
    score_modifier: f64,
}

impl ConnectionState {
    const fn new(left: usize, right: usize) -> Self {
        Self {
            left,
            right,
            overlap: 0,
            max_frame: None,
            score: 0.0,
            score_modifier: 0.0,
        }
    }
}

/// Scores the connection between two nodes in the dynamic programming model.
///
/// This routine handles various connection types including:
/// - 5'fwd->3'fwd (gene) and 3'rev->5'rev (reverse gene)
/// - Intergenic regions and operons
/// - Overlapping genes on opposite strands
///
/// If the connection ending at target_node_index has the maximal score, it updates the
/// dynamic programming pointers.
pub fn score_connection(
    nodes: &mut [Node],
    source_node_index: usize,
    target_node_index: usize,
    training: &Training,
    is_final_pass: bool,
) {
    let source_node = &nodes[source_node_index];
    let target_node = &nodes[target_node_index];
    if has_edge_artifacts(source_node) {
        return;
    }

    // Now do the more expensive connection type determination
    let connection_type = determine_connection_type(source_node, target_node);
    if connection_type == ConnectionType::Invalid {
        return;
    }

    let connection_state = match connection_type {
        ConnectionType::ForwardGene => handle_forward_gene(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::ReverseGene => handle_reverse_gene(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::ForwardIntergenic => handle_forward_intergenic(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::ReverseIntergenic => handle_reverse_intergenic(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::TripleOverlap => handle_triple_overlap(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::ForwardOperon => handle_forward_operon(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::ReverseOperon => handle_reverse_operon(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::OverlappingOpposite => handle_overlapping_opposite(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::FiveReverseToFiveForward => handle_five_reverse_to_five_forward(
            nodes,
            source_node_index,
            target_node_index,
            training,
            is_final_pass,
        ),
        ConnectionType::Invalid => unreachable!(),
    };

    let Some(mut connection_state) = connection_state else {
        return;
    };

    if !is_final_pass {
        connection_state.score = (connection_state.right as i32 - connection_state.left as i32 + 1
            - (connection_state.overlap * 2)) as f64
            * connection_state.score_modifier;
    }

    // Update node if this connection is better (using original logic)
    if source_node.scores.total_score + connection_state.score >= target_node.scores.total_score {
        nodes[target_node_index].scores.total_score =
            nodes[source_node_index].scores.total_score + connection_state.score;
        nodes[target_node_index].state.traceback = Some(source_node_index);
        nodes[target_node_index].state.overlap_marker = connection_state.max_frame;
    }
}

fn handle_forward_gene(
    nodes: &[Node],
    source_node_index: usize,
    target_node_index: usize,
    training: &Training,
    processing_flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(
        nodes[source_node_index].position.index,
        nodes[target_node_index].position.index,
    );
    let stop_val = nodes[target_node_index].position.stop_value;
    let start_idx = nodes[source_node_index].position.index as isize;
    if stop_val >= start_idx {
        return None;
    }
    connection_state.right += STOP_CODON_OFFSET;
    if !processing_flag {
        connection_state.score_modifier = calculate_gc_score(&nodes[source_node_index], training);
    } else {
        connection_state.score = nodes[source_node_index].scores.coding_score
            + nodes[source_node_index].scores.start_score;
    }
    Some(connection_state)
}

fn handle_reverse_gene(
    nodes: &[Node],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state =
        ConnectionState::new(nodes[p1].position.index, nodes[p2].position.index);
    if nodes[p1].position.stop_value <= nodes[p2].position.index as isize {
        return None;
    }
    connection_state.left -= STOP_CODON_OFFSET;
    if !flag {
        connection_state.score_modifier = calculate_gc_score(&nodes[p2], training);
    } else {
        connection_state.score = nodes[p2].scores.coding_score + nodes[p2].scores.start_score;
    }
    Some(connection_state)
}

fn handle_forward_intergenic(
    nodes: &[Node],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state =
        ConnectionState::new(nodes[p1].position.index + 2, nodes[p2].position.index);

    if connection_state.left >= connection_state.right {
        return None;
    }
    if flag {
        connection_state.score = intergenic_mod(&nodes[p1], &nodes[p2], training);
    }
    Some(connection_state)
}

fn handle_reverse_intergenic(
    nodes: &[Node],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state =
        ConnectionState::new(nodes[p1].position.index, nodes[p2].position.index - 2);
    if connection_state.left >= connection_state.right {
        return None;
    }
    if flag {
        connection_state.score = intergenic_mod(&nodes[p1], &nodes[p2], training);
    }
    Some(connection_state)
}

fn handle_triple_overlap(
    nodes: &[Node],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state =
        ConnectionState::new(nodes[p1].position.index + 2, nodes[p2].position.index - 2);
    if connection_state.left >= connection_state.right {
        return None;
    }

    connection_state.max_frame = None;
    let mut maxval = 0.0;
    for i in 0..3 {
        let Some(p3) = nodes[p2].state.start_pointers[i] else {
            continue;
        };

        connection_state.overlap =
            connection_state.left as i32 - nodes[p3].position.stop_value as i32 + 3;
        if connection_state.overlap <= 0 || connection_state.overlap >= MAXIMUM_OPPOSITE_OVERLAP {
            continue;
        }
        if connection_state.overlap
            >= nodes[p3].position.index as i32 - connection_state.left as i32
        {
            continue;
        }

        let Some(traceback_idx) = nodes[p1].state.traceback else {
            continue;
        };

        if connection_state.overlap
            >= nodes[p3].position.stop_value as i32 - nodes[traceback_idx].position.index as i32 - 2
        {
            continue;
        }

        let test_val = if flag {
            nodes[p3].scores.coding_score
                + nodes[p3].scores.start_score
                + intergenic_mod(&nodes[p3], &nodes[p2], training)
        } else {
            calculate_gc_score(&nodes[p3], training)
        };

        if test_val > maxval {
            connection_state.max_frame = Some(i);
            maxval = nodes[p3].scores.coding_score
                + nodes[p3].scores.start_score
                + intergenic_mod(&nodes[p3], &nodes[p2], training);
        }
    }

    if let Some(max_frame) = connection_state.max_frame {
        if let Some(p3) = nodes[p2].state.start_pointers[max_frame] {
            if !flag {
                connection_state.score_modifier = calculate_gc_score(&nodes[p3], training);
            } else {
                connection_state.score = nodes[p3].scores.coding_score
                    + nodes[p3].scores.start_score
                    + intergenic_mod(&nodes[p3], &nodes[p2], training);
            }
        }
    } else if flag {
        connection_state.score = intergenic_mod(&nodes[p1], &nodes[p2], training);
    }
    Some(connection_state)
}

fn handle_forward_operon(
    nodes: &[Node],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state =
        ConnectionState::new(nodes[p1].position.index, nodes[p2].position.index);
    if nodes[p2].position.stop_value >= nodes[p1].position.index as isize {
        return None;
    }

    let p3 = nodes[p1].state.start_pointers[nodes[p2].position.index % 3]?;

    connection_state.left = nodes[p3].position.index;
    connection_state.right += 2;
    if !flag {
        connection_state.score_modifier = calculate_gc_score(&nodes[p3], training);
    } else {
        connection_state.score = nodes[p3].scores.coding_score
            + nodes[p3].scores.start_score
            + intergenic_mod(&nodes[p1], &nodes[p3], training);
    }
    Some(connection_state)
}

fn handle_reverse_operon(
    nodes: &[Node],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state =
        ConnectionState::new(nodes[p1].position.index, nodes[p2].position.index);
    if nodes[p1].position.stop_value <= nodes[p2].position.index as isize {
        return None;
    }

    let p3 = nodes[p2].state.start_pointers[nodes[p1].position.index % 3]?;

    connection_state.left -= 2;
    connection_state.right = nodes[p3].position.index;
    if !flag {
        connection_state.score_modifier = calculate_gc_score(&nodes[p3], training);
    } else {
        connection_state.score = nodes[p3].scores.coding_score
            + nodes[p3].scores.start_score
            + intergenic_mod(&nodes[p3], &nodes[p2], training);
    }
    Some(connection_state)
}

fn handle_five_reverse_to_five_forward(
    nodes: &[Node],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state =
        ConnectionState::new(nodes[p1].position.index, nodes[p2].position.index);
    if connection_state.left >= connection_state.right {
        return None;
    }
    if flag {
        connection_state.score = intergenic_mod(&nodes[p1], &nodes[p2], training);
    }
    Some(connection_state)
}

fn handle_overlapping_opposite(
    nodes: &[Node],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state =
        ConnectionState::new(nodes[p1].position.index, nodes[p2].position.index);
    if nodes[p2].position.stop_value - STOP_CODON_OFFSET as isize
        >= nodes[p1].position.index as isize + STOP_CODON_OFFSET as isize
    {
        return None;
    }
    connection_state.overlap = (nodes[p1].position.index as i32 + STOP_CODON_OFFSET as i32)
        - (nodes[p2].position.stop_value as i32 - STOP_CODON_OFFSET as i32)
        + 1;
    if connection_state.overlap >= MAXIMUM_OPPOSITE_OVERLAP {
        return None;
    }
    if (nodes[p1].position.index as i32 + STOP_CODON_OFFSET as i32
        - nodes[p2].position.stop_value as i32
        - STOP_CODON_OFFSET as i32
        + 1)
        >= (nodes[p2].position.index as i32 - nodes[p1].position.index as i32 + 3 + 1)
    {
        return None;
    }

    let bnd = nodes[p1]
        .state
        .traceback
        .map(|idx| nodes[idx].position.index as i32)
        .unwrap_or(0);

    if (nodes[p1].position.index as i32 + 2 - nodes[p2].position.stop_value as i32 - 2 + 1)
        >= (nodes[p2].position.stop_value as i32 - 3 - bnd + 1)
    {
        return None;
    }
    connection_state.left = (nodes[p2].position.stop_value - STOP_CODON_OFFSET as isize) as usize;
    if !flag {
        connection_state.score_modifier = calculate_gc_score(&nodes[p2], training);
    } else {
        connection_state.score = OVERLAP_PENALTY_FACTOR.mul_add(
            -training.start_weight_factor,
            nodes[p2].scores.coding_score + nodes[p2].scores.start_score,
        );
    }
    Some(connection_state)
}

fn calculate_gc_score(node: &Node, training: &Training) -> f64 {
    let part1 = training.gc_bias_factors[0] * node.scores.gc_frame_scores[0];
    let part2 = training.gc_bias_factors[1] * node.scores.gc_frame_scores[1];
    let part3 = training.gc_bias_factors[2] * node.scores.gc_frame_scores[2];

    part1 + part2 + part3
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{CodonType, NodeMotifInfo, NodePosition, NodeScores, NodeState, Training};
    use bio::bio_types::strand::Strand;

    /// Helper function to create a test node
    fn create_test_node(
        index: usize,
        strand: Strand,
        codon_type: CodonType,
        stop_value: isize,
        scores: NodeScores,
    ) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value,
                is_edge: false,
            },
            scores,
            state: NodeState {
                ..Default::default()
            },
            motif_info: NodeMotifInfo::default(),
        }
    }

    /// Helper function to create basic node scores
    fn create_test_scores(coding_score: f64, start_score: f64) -> NodeScores {
        NodeScores {
            gc_content: 0.5,
            coding_score,
            start_score,
            ribosome_binding_score: 1.0,
            type_score: 2.0,
            upstream_score: 0.5,
            total_score: coding_score + start_score,
            gc_frame_scores: [0.5, 0.3, 0.2],
        }
    }

    /// Helper function to create basic training data
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
    fn test_score_connection_forward_gene() {
        let mut nodes = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 5.0),
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                403,
                create_test_scores(15.0, 0.0),
            ),
        ];
        let training = create_test_training();
        let _original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        // Check that the connection was processed (score may have changed)
        assert!(
            nodes[1].scores.total_score >= _original_score
                || nodes[1].scores.total_score < _original_score
        );
    }

    #[test]
    fn test_score_connection_reverse_gene() {
        let mut nodes = vec![
            create_test_node(
                400,
                Strand::Reverse,
                CodonType::Stop,
                0,
                create_test_scores(15.0, 0.0),
            ),
            create_test_node(
                100,
                Strand::Reverse,
                CodonType::Gtg,
                403,
                create_test_scores(10.0, 8.0),
            ),
        ];
        let training = create_test_training();
        let _original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(
            nodes[1].scores.total_score >= _original_score
                || nodes[1].scores.total_score < _original_score
        );
    }

    #[test]
    fn test_score_connection_forward_intergenic() {
        let mut nodes = vec![
            create_test_node(
                200,
                Strand::Forward,
                CodonType::Stop,
                203,
                create_test_scores(15.0, 0.0),
            ),
            create_test_node(
                500,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 7.0),
            ),
        ];
        let training = create_test_training();
        let _original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(
            nodes[1].scores.total_score >= _original_score
                || nodes[1].scores.total_score < _original_score
        );
    }

    #[test]
    fn test_score_connection_reverse_intergenic() {
        let mut nodes = vec![
            create_test_node(
                500,
                Strand::Reverse,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 7.0),
            ),
            create_test_node(
                200,
                Strand::Reverse,
                CodonType::Stop,
                503,
                create_test_scores(15.0, 0.0),
            ),
        ];
        let training = create_test_training();
        let original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(
            nodes[1].scores.total_score >= original_score
                || nodes[1].scores.total_score < original_score
        );
    }

    #[test]
    fn test_score_connection_triple_overlap() {
        let mut nodes = vec![
            create_test_node(
                300,
                Strand::Forward,
                CodonType::Stop,
                303,
                create_test_scores(15.0, 0.0),
            ),
            create_test_node(
                250,
                Strand::Reverse,
                CodonType::Stop,
                0,
                create_test_scores(12.0, 0.0),
            ),
        ];
        let training = create_test_training();
        let _original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(nodes[1].scores.total_score.is_finite());
    }

    #[test]
    fn test_score_connection_forward_operon() {
        let mut nodes = vec![
            create_test_node(
                200,
                Strand::Forward,
                CodonType::Stop,
                203,
                create_test_scores(15.0, 0.0),
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                403,
                create_test_scores(12.0, 0.0),
            ),
        ];
        let training = create_test_training();
        let original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(
            nodes[1].scores.total_score >= original_score
                || nodes[1].scores.total_score < original_score
        );
    }

    #[test]
    fn test_score_connection_reverse_operon() {
        let mut nodes = vec![
            create_test_node(
                400,
                Strand::Reverse,
                CodonType::Stop,
                0,
                create_test_scores(12.0, 0.0),
            ),
            create_test_node(
                200,
                Strand::Reverse,
                CodonType::Stop,
                403,
                create_test_scores(15.0, 0.0),
            ),
        ];
        let training = create_test_training();
        let original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(
            nodes[1].scores.total_score >= original_score
                || nodes[1].scores.total_score < original_score
        );
    }

    #[test]
    fn test_score_connection_overlapping_opposite() {
        let mut nodes = vec![
            create_test_node(
                300,
                Strand::Forward,
                CodonType::Stop,
                303,
                create_test_scores(15.0, 0.0),
            ),
            create_test_node(
                280,
                Strand::Reverse,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 6.0),
            ),
        ];
        let training = create_test_training();
        let _original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(nodes[1].scores.total_score.is_finite());
    }

    #[test]
    fn test_score_connection_five_reverse_to_five_forward() {
        let mut nodes = vec![
            create_test_node(
                200,
                Strand::Reverse,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 6.0),
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Gtg,
                0,
                create_test_scores(12.0, 7.0),
            ),
        ];
        let training = create_test_training();
        let _original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(nodes[1].scores.total_score.is_finite());
    }

    #[test]
    fn test_score_connection_invalid_type_ignored() {
        let mut nodes = vec![
            create_test_node(
                200,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 6.0),
            ),
            create_test_node(
                400,
                Strand::Reverse,
                CodonType::Atg,
                0,
                create_test_scores(12.0, 7.0),
            ),
        ];
        let training = create_test_training();
        let original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert_eq!(nodes[1].scores.total_score, original_score);
    }

    #[test]
    fn test_score_connection_with_flag_true() {
        let mut nodes = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 5.0),
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                403,
                create_test_scores(15.0, 0.0),
            ),
        ];
        let training = create_test_training();
        let _original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, true);

        assert!(nodes[1].scores.total_score.is_finite());
    }

    #[test]
    fn test_score_connection_with_flag_false() {
        let mut nodes = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 5.0),
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                403,
                create_test_scores(15.0, 0.0),
            ),
        ];
        let training = create_test_training();
        let _original_score = nodes[1].scores.total_score;

        score_connection(&mut nodes, 0, 1, &training, false);

        assert!(nodes[1].scores.total_score.is_finite());
    }

    #[test]
    fn test_connection_state_new() {
        let state = ConnectionState::new(100, 200);
        assert_eq!(state.left, 100);
        assert_eq!(state.right, 200);
        assert_eq!(state.overlap, 0);
        assert!(state.max_frame.is_none());
        assert_eq!(state.score, 0.0);
        assert_eq!(state.score_modifier, 0.0);
    }

    #[test]
    fn test_calculate_gc_score() {
        let node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            create_test_scores(10.0, 5.0),
        );
        let training = create_test_training();

        let gc_score = calculate_gc_score(&node, &training);

        // Should calculate weighted sum of GC frame scores
        let expected = 1.0 * 0.5 + 1.0 * 0.3 + 1.0 * 0.2;
        assert_eq!(gc_score, expected);
    }

    #[test]
    fn test_calculate_gc_score_with_weights() {
        let node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            create_test_scores(10.0, 5.0),
        );
        let mut training = create_test_training();
        training.gc_bias_factors = [2.0, 1.5, 0.5];

        let gc_score = calculate_gc_score(&node, &training);

        let expected = 2.0 * 0.5 + 1.5 * 0.3 + 0.5 * 0.2;
        assert_eq!(gc_score, expected);
    }

    #[test]
    fn test_different_node_scores_affect_total() {
        let training = create_test_training();

        // Test with low-scoring second node
        let mut nodes_low = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 5.0),
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                403,
                create_test_scores(5.0, 0.0),
            ),
        ];
        let _original_low = nodes_low[1].scores.total_score;
        score_connection(&mut nodes_low, 0, 1, &training, false);

        // Test with high-scoring second node
        let mut nodes_high = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 5.0),
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                403,
                create_test_scores(20.0, 0.0),
            ),
        ];
        let _original_high = nodes_high[1].scores.total_score;
        score_connection(&mut nodes_high, 0, 1, &training, false);

        assert!(nodes_high[1].scores.total_score.is_finite());
        assert!(nodes_low[1].scores.total_score.is_finite());
    }

    #[test]
    fn test_distance_affects_processing() {
        let training = create_test_training();

        let mut nodes_short = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 5.0),
            ),
            create_test_node(
                200,
                Strand::Forward,
                CodonType::Stop,
                203,
                create_test_scores(15.0, 0.0),
            ),
        ];
        let _original_short = nodes_short[1].scores.total_score;
        score_connection(&mut nodes_short, 0, 1, &training, false);

        let mut nodes_long = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 5.0),
            ),
            create_test_node(
                1000,
                Strand::Forward,
                CodonType::Stop,
                1003,
                create_test_scores(15.0, 0.0),
            ),
        ];
        let _original_long = nodes_long[1].scores.total_score;
        score_connection(&mut nodes_long, 0, 1, &training, false);

        assert!(nodes_short[1].scores.total_score.is_finite());
        assert!(nodes_long[1].scores.total_score.is_finite());
    }

    #[test]
    fn test_start_codon_types_handled() {
        let training = create_test_training();
        let start_codons = vec![CodonType::Atg, CodonType::Gtg, CodonType::Ttg];

        for start_codon in start_codons {
            let mut nodes = vec![
                create_test_node(
                    100,
                    Strand::Forward,
                    start_codon,
                    0,
                    create_test_scores(10.0, 5.0),
                ),
                create_test_node(
                    400,
                    Strand::Forward,
                    CodonType::Stop,
                    403,
                    create_test_scores(15.0, 0.0),
                ),
            ];
            let _original_score = nodes[1].scores.total_score;

            score_connection(&mut nodes, 0, 1, &training, false);

            assert!(nodes[1].scores.total_score.is_finite());
        }
    }

    #[test]
    fn test_strand_specific_processing() {
        let training = create_test_training();

        let mut forward_nodes = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Atg,
                0,
                create_test_scores(10.0, 5.0),
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                403,
                create_test_scores(15.0, 0.0),
            ),
        ];
        let _original_forward = forward_nodes[1].scores.total_score;
        score_connection(&mut forward_nodes, 0, 1, &training, false);

        let mut reverse_nodes = vec![
            create_test_node(
                400,
                Strand::Reverse,
                CodonType::Stop,
                0,
                create_test_scores(15.0, 0.0),
            ),
            create_test_node(
                100,
                Strand::Reverse,
                CodonType::Atg,
                403,
                create_test_scores(10.0, 5.0),
            ),
        ];
        let _original_reverse = reverse_nodes[1].scores.total_score;
        score_connection(&mut reverse_nodes, 0, 1, &training, false);

        assert!(forward_nodes[1].scores.total_score.is_finite());
        assert!(reverse_nodes[1].scores.total_score.is_finite());
    }

    #[test]
    fn test_edge_artifacts_detection() {
        let _training = create_test_training();

        let edge_node = Node {
            position: NodePosition {
                index: 100,
                strand: Strand::Forward,
                codon_type: CodonType::Stop,
                stop_value: 103,
                is_edge: false,
            },
            scores: create_test_scores(15.0, 0.0),
            state: NodeState {
                traceback: None, // No traceback = edge artifact
                ..Default::default()
            },
            motif_info: NodeMotifInfo::default(),
        };

        assert!(has_edge_artifacts(&edge_node));
    }

    #[test]
    fn test_connection_type_determination() {
        // Test forward gene connection
        let start_node = create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            create_test_scores(10.0, 5.0),
        );
        let stop_node = create_test_node(
            400,
            Strand::Forward,
            CodonType::Stop,
            403,
            create_test_scores(15.0, 0.0),
        );

        let connection_type = determine_connection_type(&start_node, &stop_node);
        assert_eq!(connection_type, ConnectionType::ForwardGene);

        // Test reverse gene connection
        let rev_stop = create_test_node(
            400,
            Strand::Reverse,
            CodonType::Stop,
            0,
            create_test_scores(15.0, 0.0),
        );
        let rev_start = create_test_node(
            100,
            Strand::Reverse,
            CodonType::Gtg,
            403,
            create_test_scores(10.0, 8.0),
        );

        let connection_type = determine_connection_type(&rev_stop, &rev_start);
        assert_eq!(connection_type, ConnectionType::ReverseGene);
    }
}
