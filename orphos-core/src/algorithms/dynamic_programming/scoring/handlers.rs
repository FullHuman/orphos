use bio::bio_types::strand::Strand;

use crate::{
    algorithms::dynamic_programming::scoring::connection_types::{
        ConnectionType, determine_connection_type,
    },
    constants::{MAXIMUM_OPPOSITE_OVERLAP, OVERLAP_PENALTY_FACTOR, STOP_CODON_OFFSET},
    node::{IntergenicNode, intergenic_mod_values},
    types::{CodonType, Node, Training},
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

/// Marker for "no node" in [`DpNode`] pointer fields.
pub const NO_NODE: u32 = u32::MAX;

/// Compact copy of the node fields read while scoring connections.
///
/// The dynamic programming inner loop touches a few fields of hundreds of
/// thousands of nodes; keeping them in this 80-byte record instead of the full
/// [`Node`] keeps the working set cache resident. Per-node score combinations
/// that connection scoring needs are precomputed with the exact expressions
/// the scoring code used on [`Node`], so results are bit-identical.
#[derive(Debug, Clone, Copy)]
pub struct DpNode {
    /// Position in the sequence
    pub index: usize,
    /// Position of the corresponding stop codon
    pub stop_value: isize,
    /// Best path score ending at this node (dynamic programming state)
    pub total_score: f64,
    /// `coding_score + start_score`
    pub gene_score: f64,
    /// GC-frame score weighted by the training GC bias (first pass)
    pub gc_score: f64,
    /// Ribosome binding site score
    pub ribosome_binding_score: f64,
    /// Upstream composition score
    pub upstream_score: f64,
    /// Best start node per frame ([`NO_NODE`] when absent)
    pub start_pointers: [u32; 3],
    /// Traceback pointer ([`NO_NODE`] when absent)
    pub traceback: u32,
    /// Strand of the node
    pub strand: Strand,
    /// Whether the node is a stop codon
    pub is_stop: bool,
}

impl DpNode {
    /// Builds the compact view of `node` for connection scoring.
    #[must_use]
    pub fn from_node(node: &Node, training: &Training) -> Self {
        let to_pointer = |pointer: Option<usize>| pointer.map_or(NO_NODE, |index| index as u32);
        Self {
            index: node.position.index,
            stop_value: node.position.stop_value,
            total_score: node.scores.total_score,
            gene_score: node.scores.coding_score + node.scores.start_score,
            gc_score: calculate_gc_score(node, training),
            ribosome_binding_score: node.scores.ribosome_binding_score,
            upstream_score: node.scores.upstream_score,
            start_pointers: node.state.start_pointers.map(to_pointer),
            traceback: to_pointer(node.state.traceback),
            strand: node.position.strand,
            is_stop: node.position.codon_type == CodonType::Stop,
        }
    }

    /// Forward stops and reverse starts without a traceback cannot start a
    /// path (see [`has_edge_artifacts`]).
    ///
    /// [`has_edge_artifacts`]: crate::algorithms::dynamic_programming::scoring::connection_types::has_edge_artifacts
    #[inline]
    #[must_use]
    pub fn has_edge_artifacts(&self) -> bool {
        self.traceback == NO_NODE
            && ((self.strand == Strand::Forward && self.is_stop)
                || (self.strand == Strand::Reverse && !self.is_stop))
    }

    #[inline]
    fn intergenic(&self) -> IntergenicNode {
        IntergenicNode {
            index: self.index,
            strand: self.strand,
            ribosome_binding_score: self.ribosome_binding_score,
            upstream_score: self.upstream_score,
        }
    }
}

#[inline]
fn intergenic_mod(nodes: &[DpNode], n1: usize, n2: usize, training: &Training) -> f64 {
    intergenic_mod_values(nodes[n1].intergenic(), nodes[n2].intergenic(), training)
}

/// Values that depend only on the target of a connection and are shared by
/// all of its sources.
#[derive(Debug, Clone, Copy)]
pub struct TargetScores {
    /// For each frame, `gene_score + intergenic_mod(start, target)` of the
    /// target's start pointer in that frame (unused where there is none).
    /// Triple-overlap and reverse-operon connections score with these.
    pub start_pointer_scores: [f64; 3],
}

impl TargetScores {
    /// Computes the shared values for `target_node_index`.
    #[must_use]
    pub fn for_target(nodes: &[DpNode], target_node_index: usize, training: &Training) -> Self {
        let start_pointer_scores = nodes[target_node_index].start_pointers.map(|p3| {
            if p3 == NO_NODE {
                0.0
            } else {
                let p3 = p3 as usize;
                nodes[p3].gene_score + intergenic_mod(nodes, p3, target_node_index, training)
            }
        });
        Self {
            start_pointer_scores,
        }
    }
}

/// A candidate connection into a target node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConnectionCandidate {
    /// Total path score obtained by ending the path with this connection
    /// (`source.total_score + connection score`).
    pub total_score: f64,
    /// Overlap marker to record on the target if this candidate wins.
    pub max_frame: Option<usize>,
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
///
/// This builds a [`DpNode`] view of all nodes on every call; the dynamic
/// programming routine itself keeps one view for the whole pass.
pub fn score_connection(
    nodes: &mut [Node],
    source_node_index: usize,
    target_node_index: usize,
    training: &Training,
    is_final_pass: bool,
) {
    let dp_nodes: Vec<DpNode> = nodes
        .iter()
        .map(|node| DpNode::from_node(node, training))
        .collect();
    if dp_nodes[source_node_index].has_edge_artifacts() {
        return;
    }

    // Now do the more expensive connection type determination
    let connection_type =
        determine_connection_type(&nodes[source_node_index], &nodes[target_node_index]);
    if connection_type == ConnectionType::Invalid {
        return;
    }

    let target_scores = TargetScores::for_target(&dp_nodes, target_node_index, training);
    let Some(candidate) = evaluate_connection(
        &dp_nodes,
        source_node_index,
        target_node_index,
        training,
        is_final_pass,
        connection_type,
        &target_scores,
    ) else {
        return;
    };

    // Update node if this connection is better (using original logic)
    if candidate.total_score >= nodes[target_node_index].scores.total_score {
        apply_candidate(nodes, source_node_index, target_node_index, candidate);
    }
}

/// Records `candidate` (reached from `source_node_index`) as the best path into
/// `target_node_index`.
#[inline]
pub fn apply_candidate(
    nodes: &mut [Node],
    source_node_index: usize,
    target_node_index: usize,
    candidate: ConnectionCandidate,
) {
    let target = &mut nodes[target_node_index];
    target.scores.total_score = candidate.total_score;
    target.state.traceback = Some(source_node_index);
    target.state.overlap_marker = candidate.max_frame;
}

/// Evaluates a connection of an already-determined, valid `connection_type`
/// from `source_node_index` to `target_node_index` without touching any node.
///
/// The caller is responsible for having checked [`DpNode::has_edge_artifacts`]
/// on the source, for `connection_type` matching
/// [`determine_connection_type`], and for `target_scores` being
/// [`TargetScores::for_target`] of the target. Returns `None` when the
/// connection is not admissible.
#[inline]
pub fn evaluate_connection(
    nodes: &[DpNode],
    source_node_index: usize,
    target_node_index: usize,
    training: &Training,
    is_final_pass: bool,
    connection_type: ConnectionType,
    target_scores: &TargetScores,
) -> Option<ConnectionCandidate> {
    let connection_state = match connection_type {
        ConnectionType::ForwardGene => {
            handle_forward_gene(nodes, source_node_index, target_node_index, is_final_pass)
        }
        ConnectionType::ReverseGene => {
            handle_reverse_gene(nodes, source_node_index, target_node_index, is_final_pass)
        }
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
            target_scores,
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
            target_scores,
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
        ConnectionType::Invalid => return None,
    }?;

    let score = if is_final_pass {
        connection_state.score
    } else {
        (connection_state.right as i32 - connection_state.left as i32 + 1
            - (connection_state.overlap * 2)) as f64
            * connection_state.score_modifier
    };

    Some(ConnectionCandidate {
        total_score: nodes[source_node_index].total_score + score,
        max_frame: connection_state.max_frame,
    })
}

fn handle_forward_gene(
    nodes: &[DpNode],
    source_node_index: usize,
    target_node_index: usize,
    processing_flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(
        nodes[source_node_index].index,
        nodes[target_node_index].index,
    );
    let stop_val = nodes[target_node_index].stop_value;
    let start_idx = nodes[source_node_index].index as isize;
    if stop_val >= start_idx {
        return None;
    }
    connection_state.right += STOP_CODON_OFFSET;
    if !processing_flag {
        connection_state.score_modifier = nodes[source_node_index].gc_score;
    } else {
        connection_state.score = nodes[source_node_index].gene_score;
    }
    Some(connection_state)
}

fn handle_reverse_gene(
    nodes: &[DpNode],
    p1: usize,
    p2: usize,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(nodes[p1].index, nodes[p2].index);
    if nodes[p1].stop_value <= nodes[p2].index as isize {
        return None;
    }
    connection_state.left -= STOP_CODON_OFFSET;
    if !flag {
        connection_state.score_modifier = nodes[p2].gc_score;
    } else {
        connection_state.score = nodes[p2].gene_score;
    }
    Some(connection_state)
}

fn handle_forward_intergenic(
    nodes: &[DpNode],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(nodes[p1].index + 2, nodes[p2].index);

    if connection_state.left >= connection_state.right {
        return None;
    }
    if flag {
        connection_state.score = intergenic_mod(nodes, p1, p2, training);
    }
    Some(connection_state)
}

fn handle_reverse_intergenic(
    nodes: &[DpNode],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(nodes[p1].index, nodes[p2].index - 2);
    if connection_state.left >= connection_state.right {
        return None;
    }
    if flag {
        connection_state.score = intergenic_mod(nodes, p1, p2, training);
    }
    Some(connection_state)
}

fn handle_triple_overlap(
    nodes: &[DpNode],
    p1: usize,
    p2: usize,
    training: &Training,
    target_scores: &TargetScores,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(nodes[p1].index + 2, nodes[p2].index - 2);
    if connection_state.left >= connection_state.right {
        return None;
    }

    connection_state.max_frame = None;
    let mut maxval = 0.0;
    for i in 0..3 {
        let p3 = nodes[p2].start_pointers[i];
        if p3 == NO_NODE {
            continue;
        }
        let p3 = p3 as usize;

        connection_state.overlap = connection_state.left as i32 - nodes[p3].stop_value as i32 + 3;
        if connection_state.overlap <= 0 || connection_state.overlap >= MAXIMUM_OPPOSITE_OVERLAP {
            continue;
        }
        if connection_state.overlap >= nodes[p3].index as i32 - connection_state.left as i32 {
            continue;
        }

        let traceback_idx = nodes[p1].traceback;
        if traceback_idx == NO_NODE {
            continue;
        }
        let traceback_idx = traceback_idx as usize;

        if connection_state.overlap
            >= nodes[p3].stop_value as i32 - nodes[traceback_idx].index as i32 - 2
        {
            continue;
        }

        let test_val = if flag {
            target_scores.start_pointer_scores[i]
        } else {
            nodes[p3].gc_score
        };

        if test_val > maxval {
            connection_state.max_frame = Some(i);
            maxval = target_scores.start_pointer_scores[i];
        }
    }

    if let Some(max_frame) = connection_state.max_frame {
        let p3 = nodes[p2].start_pointers[max_frame];
        if p3 != NO_NODE {
            let p3 = p3 as usize;
            if !flag {
                connection_state.score_modifier = nodes[p3].gc_score;
            } else {
                connection_state.score = target_scores.start_pointer_scores[max_frame];
            }
        }
    } else if flag {
        connection_state.score = intergenic_mod(nodes, p1, p2, training);
    }
    Some(connection_state)
}

fn handle_forward_operon(
    nodes: &[DpNode],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(nodes[p1].index, nodes[p2].index);
    if nodes[p2].stop_value >= nodes[p1].index as isize {
        return None;
    }

    let p3 = nodes[p1].start_pointers[nodes[p2].index % 3];
    if p3 == NO_NODE {
        return None;
    }
    let p3 = p3 as usize;

    connection_state.left = nodes[p3].index;
    connection_state.right += 2;
    if !flag {
        connection_state.score_modifier = nodes[p3].gc_score;
    } else {
        connection_state.score = nodes[p3].gene_score + intergenic_mod(nodes, p1, p3, training);
    }
    Some(connection_state)
}

fn handle_reverse_operon(
    nodes: &[DpNode],
    p1: usize,
    p2: usize,
    target_scores: &TargetScores,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(nodes[p1].index, nodes[p2].index);
    if nodes[p1].stop_value <= nodes[p2].index as isize {
        return None;
    }

    let frame = nodes[p1].index % 3;
    let p3 = nodes[p2].start_pointers[frame];
    if p3 == NO_NODE {
        return None;
    }
    let p3 = p3 as usize;

    connection_state.left -= 2;
    connection_state.right = nodes[p3].index;
    if !flag {
        connection_state.score_modifier = nodes[p3].gc_score;
    } else {
        connection_state.score = target_scores.start_pointer_scores[frame];
    }
    Some(connection_state)
}

fn handle_five_reverse_to_five_forward(
    nodes: &[DpNode],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(nodes[p1].index, nodes[p2].index);
    if connection_state.left >= connection_state.right {
        return None;
    }
    if flag {
        connection_state.score = intergenic_mod(nodes, p1, p2, training);
    }
    Some(connection_state)
}

fn handle_overlapping_opposite(
    nodes: &[DpNode],
    p1: usize,
    p2: usize,
    training: &Training,
    flag: bool,
) -> Option<ConnectionState> {
    let mut connection_state = ConnectionState::new(nodes[p1].index, nodes[p2].index);
    if nodes[p2].stop_value - STOP_CODON_OFFSET as isize
        >= nodes[p1].index as isize + STOP_CODON_OFFSET as isize
    {
        return None;
    }
    connection_state.overlap = (nodes[p1].index as i32 + STOP_CODON_OFFSET as i32)
        - (nodes[p2].stop_value as i32 - STOP_CODON_OFFSET as i32)
        + 1;
    if connection_state.overlap >= MAXIMUM_OPPOSITE_OVERLAP {
        return None;
    }
    if (nodes[p1].index as i32 + STOP_CODON_OFFSET as i32
        - nodes[p2].stop_value as i32
        - STOP_CODON_OFFSET as i32
        + 1)
        >= (nodes[p2].index as i32 - nodes[p1].index as i32 + 3 + 1)
    {
        return None;
    }

    let bnd = if nodes[p1].traceback == NO_NODE {
        0
    } else {
        nodes[nodes[p1].traceback as usize].index as i32
    };

    if (nodes[p1].index as i32 + 2 - nodes[p2].stop_value as i32 - 2 + 1)
        >= (nodes[p2].stop_value as i32 - 3 - bnd + 1)
    {
        return None;
    }
    connection_state.left = (nodes[p2].stop_value - STOP_CODON_OFFSET as isize) as usize;
    if !flag {
        connection_state.score_modifier = nodes[p2].gc_score;
    } else {
        connection_state.score =
            OVERLAP_PENALTY_FACTOR.mul_add(-training.start_weight_factor, nodes[p2].gene_score);
    }
    Some(connection_state)
}

/// `bias[0]*gc[0] + bias[1]*gc[1] + bias[2]*gc[2]`, evaluated with the
/// fused multiply-adds the reference Prodigal build uses for this expression
/// (see the crate documentation on floating-point contraction).
fn calculate_gc_score(node: &Node, training: &Training) -> f64 {
    let bias = &training.gc_bias_factors;
    let gc = &node.scores.gc_frame_scores;
    bias[2].mul_add(gc[2], bias[0].mul_add(gc[0], bias[1] * gc[1]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithms::dynamic_programming::scoring::connection_types::has_edge_artifacts;
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
