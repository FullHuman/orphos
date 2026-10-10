use std::collections::VecDeque;

use bio::bio_types::strand::Strand;

use crate::{
    algorithms::dynamic_programming::{
        overlap_resolution::{resolve_simple_overlaps, resolve_triple_overlaps},
        scoring::{
            connection_types::ConnectionType,
            handlers::{ConnectionCandidate, DpNode, TargetScores, evaluate_connection},
        },
    },
    constants::{MAX_NODE_DIST, MAXIMUM_OPPOSITE_OVERLAP, OPERON_DISTANCE, STOP_CODON_OFFSET},
    node::{IntergenicNode, intergenic_mod_values},
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

/// Perform the main dynamic programming computation.
///
/// Every node `i` considers connections from the nodes `j` in the window
/// `[min, i)` chosen by [`calculate_distance_constraints`]. Only a few
/// strand/codon-type combinations form valid connections, so instead of
/// scanning the whole window and rejecting most pairs, sources are kept in
/// per-kind (strand × stop/start) and per-frame lists and only compatible
/// lists are visited for each target. Sources that can never form a
/// connection (see [`DpNode::has_edge_artifacts`]) are never entered into the
/// lists.
///
/// Connections whose score does not depend on the pair (5'rev->5'fwd, and
/// intergenic connections farther apart than the operon distance, which all
/// receive the same penalty, as well as every intergenic connection in the
/// first pass, which scores zero) are resolved with sliding-window maxima
/// over the source totals instead of being evaluated one by one. Connection
/// types whose handlers reject sources outside an index range are only
/// evaluated inside that range.
///
/// The result is identical to the exhaustive scan: a target takes the
/// candidate with the highest `source.total_score + connection score`, and on
/// exact ties the source with the largest index wins, which matches the
/// `>=` update rule of the sequential scan. A target is only updated if the
/// best candidate is at least `0.0`, the initial total score.
fn compute_dynamic_programming_scores(
    nodes: &mut [Node],
    training: &Training,
    is_final_pass: bool,
) {
    let mut dp_nodes: Vec<DpNode> = nodes
        .iter()
        .map(|node| DpNode::from_node(node, training))
        .collect();
    let mut sources = SourceLists::new();
    // Score of a connection that does not depend on the pair: zero in the
    // first pass (no score modifier), the distant-pair intergenic modifier in
    // the final pass.
    let constant_score = if is_final_pass {
        distant_intergenic_mod(training)
    } else {
        0.0
    };
    let mut previous_min = 0;

    for i in 0..dp_nodes.len() {
        let min = calculate_distance_constraints(&dp_nodes, i);
        if min < i {
            if min < previous_min {
                sources.reset_windows();
            }
            sources.advance_to(min, previous_min);
            previous_min = min;

            let target = dp_nodes[i];
            let target_kind = NodeKind::of(&target);
            let target_frame = target.index % 3;
            let mut best = BestCandidate::default();
            // Only reverse stop targets have connection types that use the
            // start-pointer scores.
            let target_scores = if matches!(target_kind, NodeKind::ReverseStop) {
                TargetScores::for_target(&dp_nodes, i, training)
            } else {
                TargetScores {
                    start_pointer_scores: [0.0; 3],
                }
            };

            {
                let dp_nodes = &dp_nodes;
                let mut consider = |source_list: &[usize], connection_type: ConnectionType| {
                    for &j in source_list {
                        if let Some(candidate) = evaluate_connection(
                            dp_nodes,
                            j,
                            i,
                            training,
                            is_final_pass,
                            connection_type,
                            &target_scores,
                        ) {
                            best.offer(j, candidate);
                        }
                    }
                };

                let is_constant =
                    |source: &DpNode| intergenic_is_constant(source, &target, is_final_pass);

                match target_kind {
                    NodeKind::ForwardStart => {
                        // 3'fwd->5'fwd: distant sources share one score; near
                        // sources are scored individually in the final pass.
                        let forward_stops = sources.window(NodeKind::ForwardStop);
                        let split = intergenic_split(forward_stops, dp_nodes, is_constant);
                        consider(&forward_stops[split..], ConnectionType::ForwardIntergenic);
                        sources.admit(FORWARD_INTERGENIC_WINDOW, dp_nodes, is_constant);

                        // 5'rev->5'fwd always scores the constant.
                        sources.admit(FIVE_REVERSE_FIVE_FORWARD_WINDOW, dp_nodes, |source| {
                            source.index < target.index
                        });
                    }
                    NodeKind::ForwardStop => {
                        // 5'fwd->3'fwd rejects starts at or before the previous stop.
                        let starts = sources.window_in_frame(NodeKind::ForwardStart, target_frame);
                        let first = starts
                            .partition_point(|&j| dp_nodes[j].index as isize <= target.stop_value);
                        consider(&starts[first..], ConnectionType::ForwardGene);

                        // 3'fwd->3'fwd rejects stops at or before the previous stop.
                        let stops = sources.window(NodeKind::ForwardStop);
                        let first = stops
                            .partition_point(|&j| dp_nodes[j].index as isize <= target.stop_value);
                        consider(&stops[first..], ConnectionType::ForwardOperon);
                    }
                    NodeKind::ReverseStop => {
                        // 5'rev->3'rev: as for the forward intergenic case.
                        let reverse_starts = sources.window(NodeKind::ReverseStart);
                        let split = intergenic_split(reverse_starts, dp_nodes, is_constant);
                        consider(&reverse_starts[split..], ConnectionType::ReverseIntergenic);
                        sources.admit(REVERSE_INTERGENIC_WINDOW, dp_nodes, is_constant);

                        consider(
                            sources.window(NodeKind::ForwardStop),
                            ConnectionType::TripleOverlap,
                        );
                        consider(
                            sources.window(NodeKind::ReverseStop),
                            ConnectionType::ReverseOperon,
                        );
                    }
                    NodeKind::ReverseStart => {
                        consider(
                            sources.window_in_frame(NodeKind::ReverseStop, target_frame),
                            ConnectionType::ReverseGene,
                        );

                        // 3'fwd->5'rev rejects stops outside
                        // [stop_value - 3, stop_value + MAXIMUM_OPPOSITE_OVERLAP - 6].
                        let stops = sources.window(NodeKind::ForwardStop);
                        let lowest = target.stop_value - 3;
                        let highest = target.stop_value + MAXIMUM_OPPOSITE_OVERLAP as isize - 6;
                        let first =
                            stops.partition_point(|&j| (dp_nodes[j].index as isize) < lowest);
                        let end = first
                            + stops[first..]
                                .partition_point(|&j| dp_nodes[j].index as isize <= highest);
                        consider(&stops[first..end], ConnectionType::OverlappingOpposite);
                    }
                }

                let windows: &[usize] = match target_kind {
                    NodeKind::ForwardStart => {
                        &[FORWARD_INTERGENIC_WINDOW, FIVE_REVERSE_FIVE_FORWARD_WINDOW]
                    }
                    NodeKind::ReverseStop => &[REVERSE_INTERGENIC_WINDOW],
                    _ => &[],
                };
                for &window in windows {
                    if let Some((j, total_score)) = sources.window_best(window) {
                        best.offer(
                            j,
                            ConnectionCandidate {
                                total_score,
                                max_frame: None,
                            },
                        );
                    }
                }
            }

            if let Some(source) = best.source {
                dp_nodes[i].total_score = best.candidate.total_score;
                dp_nodes[i].traceback = source as u32;
                let node = &mut nodes[i];
                node.scores.total_score = best.candidate.total_score;
                node.state.traceback = Some(source);
                node.state.overlap_marker = best.candidate.max_frame;
            }
        }

        // Node `i` is final now; it can act as a source for later targets.
        let node = dp_nodes[i];
        if !node.has_edge_artifacts() {
            sources.push(
                NodeKind::of(&node),
                node.index % 3,
                i,
                node.total_score + constant_score,
            );
        }
    }
}

/// Intergenic modifier for a pair of nodes that are on opposite strands or
/// farther apart than three times the operon distance. The modifier is the
/// same for every such pair.
fn distant_intergenic_mod(training: &Training) -> f64 {
    let distant = |index, strand| IntergenicNode {
        index,
        strand,
        ribosome_binding_score: 0.0,
        upstream_score: 0.0,
    };
    intergenic_mod_values(
        distant(0, Strand::Reverse),
        distant(usize::MAX / 2, Strand::Forward),
        training,
    )
}

/// Distance beyond which `intergenic_mod` only applies the distant-pair
/// penalty (`dist > 3.0 * OPERON_DISTANCE`).
const DISTANT_INTERGENIC_DISTANCE: usize = 3 * OPERON_DISTANCE as usize;

/// Whether an intergenic connection from `source` to `target` scores the
/// pass constant: in the first pass every admissible connection does, in the
/// final pass those farther apart than three times the operon distance.
#[inline]
fn intergenic_is_constant(source: &DpNode, target: &DpNode, is_final_pass: bool) -> bool {
    if is_final_pass {
        source.index + DISTANT_INTERGENIC_DISTANCE < target.index
    } else {
        source.index + STOP_CODON_OFFSET < target.index
    }
}

/// Position in `sources` (sorted by index) of the first source whose
/// intergenic connection must be scored individually; the ones before it
/// score the pass constant (or are rejected by the handler in the first pass,
/// where nothing needs individual scoring).
#[inline]
fn intergenic_split(
    sources: &[usize],
    dp_nodes: &[DpNode],
    is_constant: impl Fn(&DpNode) -> bool,
) -> usize {
    sources.partition_point(|&j| is_constant(&dp_nodes[j]))
}

/// Strand × start/stop classification of a node for connection dispatch.
#[derive(Clone, Copy)]
enum NodeKind {
    ForwardStart = 0,
    ForwardStop = 1,
    ReverseStop = 2,
    ReverseStart = 3,
}

impl NodeKind {
    #[inline]
    fn of(node: &DpNode) -> Self {
        match (node.strand, node.is_stop) {
            (Strand::Reverse, true) => Self::ReverseStop,
            (Strand::Reverse, false) => Self::ReverseStart,
            (_, true) => Self::ForwardStop,
            (_, false) => Self::ForwardStart,
        }
    }
}

const FORWARD_INTERGENIC_WINDOW: usize = 0;
const REVERSE_INTERGENIC_WINDOW: usize = 1;
const FIVE_REVERSE_FIVE_FORWARD_WINDOW: usize = 2;

/// Candidate sources, grouped by kind (and additionally by reading frame),
/// in increasing node order, together with a cursor per list marking the
/// start of the current distance window, and the sliding-window maxima for
/// the constant-score connection types.
struct SourceLists {
    by_kind: [Vec<usize>; 4],
    by_kind_and_frame: [[Vec<usize>; 3]; 4],
    cursors: [usize; 4],
    frame_cursors: [[usize; 3]; 4],
    windows: [ConstantWindow; 3],
}

impl SourceLists {
    fn new() -> Self {
        Self {
            by_kind: Default::default(),
            by_kind_and_frame: Default::default(),
            cursors: [0; 4],
            frame_cursors: [[0; 3]; 4],
            windows: [
                ConstantWindow::new(NodeKind::ForwardStop),
                ConstantWindow::new(NodeKind::ReverseStart),
                ConstantWindow::new(NodeKind::ReverseStart),
            ],
        }
    }

    #[inline]
    fn push(&mut self, kind: NodeKind, frame: usize, index: usize, constant_total: f64) {
        self.by_kind[kind as usize].push(index);
        self.by_kind_and_frame[kind as usize][frame].push(index);
        for window in &mut self.windows {
            if window.kind == kind as usize {
                window.keys.push(constant_total);
            }
        }
    }

    /// Moves every cursor to the first source with index `>= min` and drops
    /// sources before `min` from the sliding windows.
    #[inline]
    fn advance_to(&mut self, min: usize, previous_min: usize) {
        for (list, cursor) in self.by_kind.iter().zip(self.cursors.iter_mut()) {
            advance_cursor(list, cursor, min, previous_min);
        }
        for (lists, cursors) in self
            .by_kind_and_frame
            .iter()
            .zip(self.frame_cursors.iter_mut())
        {
            for (list, cursor) in lists.iter().zip(cursors.iter_mut()) {
                advance_cursor(list, cursor, min, previous_min);
            }
        }
        for window in &mut self.windows {
            window.evict_before(&self.by_kind[window.kind], min);
        }
    }

    /// Forgets the sliding-window state after the window lower bound moved
    /// backwards (the Prodigal giant-ORF exception); it is rebuilt on the
    /// next `admit`.
    fn reset_windows(&mut self) {
        for window in &mut self.windows {
            window.deque.clear();
            window.admitted = 0;
        }
    }

    /// Admits the sources of the window's kind that satisfy `eligible` (a
    /// predicate that, once false for a source, is false for every later
    /// source) into the sliding-window maximum.
    #[inline]
    fn admit(&mut self, window: usize, dp_nodes: &[DpNode], eligible: impl Fn(&DpNode) -> bool) {
        let window = &mut self.windows[window];
        let list = &self.by_kind[window.kind];
        window.admit(list, self.cursors[window.kind], dp_nodes, eligible);
    }

    /// Best `(source index, total score)` in the sliding window.
    #[inline]
    fn window_best(&self, window: usize) -> Option<(usize, f64)> {
        let window = &self.windows[window];
        window.best(&self.by_kind[window.kind])
    }

    #[inline]
    fn window(&self, kind: NodeKind) -> &[usize] {
        &self.by_kind[kind as usize][self.cursors[kind as usize]..]
    }

    #[inline]
    fn window_in_frame(&self, kind: NodeKind, frame: usize) -> &[usize] {
        &self.by_kind_and_frame[kind as usize][frame][self.frame_cursors[kind as usize][frame]..]
    }
}

#[inline]
fn advance_cursor(list: &[usize], cursor: &mut usize, min: usize, previous_min: usize) {
    if min >= previous_min {
        while *cursor < list.len() && list[*cursor] < min {
            *cursor += 1;
        }
    } else {
        *cursor = list.partition_point(|&j| j < min);
    }
}

/// Sliding-window maximum of `source.total_score + constant` over the sources
/// of one kind, for connection types whose score is that constant.
///
/// `keys[p]` is the key of the `p`-th source of the kind. The deque holds
/// list positions with strictly decreasing keys; a new source evicts every
/// earlier source with a key `<=` its own, so that among equal keys the
/// latest source wins, matching the `>=` rule of the sequential scan.
struct ConstantWindow {
    kind: usize,
    keys: Vec<f64>,
    deque: VecDeque<usize>,
    admitted: usize,
}

impl ConstantWindow {
    fn new(kind: NodeKind) -> Self {
        Self {
            kind: kind as usize,
            keys: Vec::new(),
            deque: VecDeque::new(),
            admitted: 0,
        }
    }

    #[inline]
    fn admit(
        &mut self,
        list: &[usize],
        cursor: usize,
        dp_nodes: &[DpNode],
        eligible: impl Fn(&DpNode) -> bool,
    ) {
        if self.admitted < cursor {
            self.admitted = cursor;
        }
        while self.admitted < list.len() {
            let position = self.admitted;
            if !eligible(&dp_nodes[list[position]]) {
                break;
            }
            self.admitted += 1;
            let key = self.keys[position];
            // A NaN total can never be selected by the sequential scan.
            if key.is_nan() {
                continue;
            }
            while let Some(&back) = self.deque.back() {
                if self.keys[back] <= key {
                    self.deque.pop_back();
                } else {
                    break;
                }
            }
            self.deque.push_back(position);
        }
    }

    #[inline]
    fn evict_before(&mut self, list: &[usize], min: usize) {
        while let Some(&front) = self.deque.front() {
            if list[front] < min {
                self.deque.pop_front();
            } else {
                break;
            }
        }
    }

    #[inline]
    fn best(&self, list: &[usize]) -> Option<(usize, f64)> {
        self.deque
            .front()
            .map(|&position| (list[position], self.keys[position]))
    }
}

/// Tracks the best candidate seen for a target, replicating the `>=` update
/// rule of a sequential scan over increasing source indices.
struct BestCandidate {
    candidate: ConnectionCandidate,
    source: Option<usize>,
}

impl Default for BestCandidate {
    fn default() -> Self {
        Self {
            candidate: ConnectionCandidate {
                total_score: 0.0,
                max_frame: None,
            },
            source: None,
        }
    }
}

impl BestCandidate {
    #[inline]
    fn offer(&mut self, source: usize, candidate: ConnectionCandidate) {
        let score = candidate.total_score;
        let best = self.candidate.total_score;
        let better =
            score > best || (score == best && self.source.is_none_or(|current| source > current));
        if better {
            self.candidate = candidate;
            self.source = Some(source);
        }
    }
}

/// Lower bound of the source window for target `i`.
///
/// Prodigal widens the window to the whole prefix for reverse starts and
/// forward stops whose stop codon lies before the regular window; its loop
/// compares the node with itself and therefore always runs down to zero, which
/// is reproduced here in constant time.
#[inline]
fn calculate_distance_constraints(nodes: &[DpNode], i: usize) -> usize {
    let mut min = i.saturating_sub(MAX_NODE_DIST);

    let node = &nodes[i];

    let lookback = match (node.strand, node.is_stop) {
        (Strand::Reverse, false) | (Strand::Forward, true) => {
            nodes[min].index as isize >= node.stop_value
        }
        _ => false,
    };
    if lookback && node.index as isize != node.stop_value {
        min = 0;
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
        let nodes = [
            create_test_node(10, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(20, Strand::Forward, CodonType::Stop, 10, false),
        ];

        let dp_nodes: Vec<DpNode> = nodes
            .iter()
            .map(|node| DpNode::from_node(node, &create_test_training()))
            .collect();
        let min = calculate_distance_constraints(&dp_nodes, 1);
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

        let dp_nodes: Vec<DpNode> = nodes
            .iter()
            .map(|node| DpNode::from_node(node, &create_test_training()))
            .collect();
        let min = calculate_distance_constraints(&dp_nodes, 600); // Use index > MAX_NODE_DIST
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
