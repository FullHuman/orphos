use std::collections::VecDeque;

use bio::bio_types::strand::Strand;

use crate::{
    algorithms::dynamic_programming::{
        overlap_resolution::{resolve_simple_overlaps, resolve_triple_overlaps},
        scoring::{
            connection_types::ConnectionType,
            handlers::{ConnectionCandidate, DpNode, NO_NODE, TargetScores, evaluate_connection},
        },
    },
    constants::{MAX_NODE_DIST, MAXIMUM_OPPOSITE_OVERLAP, OPERON_DISTANCE, STOP_CODON_OFFSET},
    node::{IntergenicNode, intergenic_mod_values},
    types::{Node, Training},
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
///
/// `nodes` must be sorted by position, as `sort_nodes_by_position` leaves
/// them.
pub fn predict_genes(
    nodes: &mut [Node],
    training: &Training,
    is_final_pass: bool,
) -> Option<usize> {
    predict_genes_in(nodes, training, is_final_pass, &mut DpWorkspace::default())
}

/// [`predict_genes`] in `workspace`, whose buffers a caller that predicts
/// genes many times (once per preset model in metagenomic mode) keeps from
/// one prediction to the next.
pub(crate) fn predict_genes_in(
    nodes: &mut [Node],
    training: &Training,
    is_final_pass: bool,
    workspace: &mut DpWorkspace,
) -> Option<usize> {
    let node_count = nodes.len();
    if node_count == 0 {
        return None;
    }

    initialize_nodes(nodes);
    compute_dynamic_programming_scores(nodes, training, is_final_pass, workspace);
    let max_index = find_maximum_terminal_node(&workspace.dp_nodes)?;

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

/// The buffers of the dynamic programming: the compact view of the nodes
/// and the lists of candidate sources.
#[derive(Default)]
pub(crate) struct DpWorkspace {
    dp_nodes: Vec<DpNode>,
    sources: SourceLists,
}

/// Perform the main dynamic programming computation, leaving the compact
/// view of the nodes it ran on, with their final total scores, in
/// `workspace`.
///
/// Every node `i` considers connections from the nodes `j` in the window
/// `[min, i)` chosen by [`calculate_distance_constraints`]. Only a few
/// strand/codon-type combinations form valid connections, and most of those
/// only for sources in a narrow range of positions, so instead of scanning
/// the whole window and rejecting most pairs, sources are kept in per-kind
/// (strand × stop/start) lists in node order and each target visits only the
/// sources that can connect to it:
///
/// * the lists are sorted by position, so sources a handler accepts only
///   inside a position range are found from the end of the list or with a
///   binary search;
/// * a reverse stop connects to a later reverse node only while that node
///   lies before the stop's own stop value, so reverse stops are dropped from
///   the lists the reverse gene and reverse operon connections read once a
///   target has passed it ([`SourceLists::live_reverse_stops`]);
/// * connections whose score does not depend on the pair are resolved with
///   sliding-window maxima over the source totals instead of being evaluated
///   one by one: 5'rev->5'fwd; intergenic connections farther apart than the
///   operon distance, which all receive the same penalty, and every
///   intergenic connection in the first pass, which scores zero; and
///   3'fwd->3'rev connections whose source lies outside the short ranges in
///   which it can overlap one of the target's genes.
///
/// Sources that can never form a connection (see
/// [`DpNode::has_edge_artifacts`]) are never entered into the lists.
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
    workspace: &mut DpWorkspace,
) {
    let DpWorkspace { dp_nodes, sources } = workspace;
    dp_nodes.clear();
    dp_nodes.extend(nodes.iter().map(|node| DpNode::from_node(node, training)));
    sources.clear();
    // Score of a connection that does not depend on the pair: zero in the
    // first pass (no score modifier), the distant-pair intergenic modifier in
    // the final pass.
    let constant_score = if is_final_pass {
        distant_intergenic_mod(training)
    } else {
        0.0
    };
    // An intergenic connection scores the pass constant when its nodes are
    // farther apart than this: in the first pass every admissible connection
    // does, in the final pass those farther apart than three times the
    // operon distance.
    let intergenic_reach = if is_final_pass {
        DISTANT_INTERGENIC_DISTANCE
    } else {
        STOP_CODON_OFFSET
    };
    let mut previous_min = 0;

    for i in 0..dp_nodes.len() {
        let min = calculate_distance_constraints(dp_nodes, i);
        if min < i {
            if min < previous_min {
                sources.reset_windows();
            }
            previous_min = min;

            let target = dp_nodes[i];
            let target_kind = NodeKind::of(&target);
            let mut best = BestCandidate::default();
            // Only reverse stop targets have connection types that use the
            // start-pointer scores.
            let target_scores = if matches!(target_kind, NodeKind::ReverseStop) {
                TargetScores::for_target(dp_nodes, i, training)
            } else {
                TargetScores {
                    start_pointer_scores: [0.0; 3],
                }
            };

            // Scores the connection from source `$source` and offers it.
            macro_rules! connect {
                ($source:expr, $connection_type:expr) => {{
                    let source: usize = $source;
                    if let Some(candidate) = evaluate_connection(
                        dp_nodes,
                        source,
                        i,
                        training,
                        is_final_pass,
                        $connection_type,
                        &target_scores,
                    ) {
                        best.offer(source, candidate);
                    }
                }};
            }
            // Scores the connections from the sources of `$kind` at the end
            // of its list, back to the first one in the window for which
            // `$rejected` holds: a condition on the source position that,
            // once true going backwards, holds for every earlier source.
            macro_rules! connect_latest {
                ($kind:expr, $connection_type:expr, |$position:ident| $rejected:expr) => {{
                    let list = &sources.kinds[$kind as usize];
                    for (&id, &$position) in list.ids.iter().zip(&list.positions).rev() {
                        if (id as usize) < min || $rejected {
                            break;
                        }
                        connect!(id as usize, $connection_type);
                    }
                }};
            }
            // Offers the best source of a constant-score window.
            macro_rules! offer_window_best {
                ($window:expr) => {{
                    if let Some(position) = sources.window_best($window, min) {
                        let list = sources.window_sources($window);
                        best.offer(
                            list.ids[position] as usize,
                            ConnectionCandidate {
                                total_score: list.keys[position],
                                max_frame: None,
                            },
                        );
                    }
                }};
            }

            match target_kind {
                NodeKind::ForwardStart => {
                    // 3'fwd->5'fwd: distant sources share one score; near
                    // sources are scored individually (and all rejected in
                    // the first pass).
                    connect_latest!(
                        NodeKind::ForwardStop,
                        ConnectionType::ForwardIntergenic,
                        |position| position + intergenic_reach < target.index
                    );
                    sources.admit(FORWARD_INTERGENIC_WINDOW, min, |position| {
                        position + intergenic_reach < target.index
                    });
                    offer_window_best!(FORWARD_INTERGENIC_WINDOW);

                    // 5'rev->5'fwd always scores the constant.
                    sources.admit(FIVE_REVERSE_FIVE_FORWARD_WINDOW, min, |position| {
                        position < target.index
                    });
                    offer_window_best!(FIVE_REVERSE_FIVE_FORWARD_WINDOW);
                }
                NodeKind::ForwardStop => {
                    // 5'fwd->3'fwd rejects starts at or before the previous
                    // stop, and starts of the other frames.
                    for &id in sources.forward_starts_by_frame[target.index % 3]
                        .iter()
                        .rev()
                    {
                        let id = id as usize;
                        if id < min || dp_nodes[id].index as isize <= target.stop_value {
                            break;
                        }
                        connect!(id, ConnectionType::ForwardGene);
                    }

                    // 3'fwd->3'fwd rejects stops at or before the previous stop.
                    connect_latest!(
                        NodeKind::ForwardStop,
                        ConnectionType::ForwardOperon,
                        |position| position as isize <= target.stop_value
                    );
                }
                NodeKind::ReverseStop => {
                    // 5'rev->3'rev: as for the forward intergenic case.
                    connect_latest!(
                        NodeKind::ReverseStart,
                        ConnectionType::ReverseIntergenic,
                        |position| position + intergenic_reach < target.index
                    );
                    sources.admit(REVERSE_INTERGENIC_WINDOW, min, |position| {
                        position + intergenic_reach < target.index
                    });
                    offer_window_best!(REVERSE_INTERGENIC_WINDOW);

                    // 3'fwd->3'rev: a source overlaps the gene of one of the
                    // target's start pointers only if it lies in a short
                    // range around that gene's stop; such sources are scored
                    // individually. Every other source the handler accepts
                    // (those at least five bases before the target) scores
                    // the pass constant.
                    let mut overlap_ranges = [(0isize, -1isize); 3];
                    for (range, &start_pointer) in
                        overlap_ranges.iter_mut().zip(&target.start_pointers)
                    {
                        if start_pointer != NO_NODE {
                            let stop_value = dp_nodes[start_pointer as usize].stop_value;
                            *range = (
                                stop_value - 4,
                                stop_value + MAXIMUM_OPPOSITE_OVERLAP as isize - 6,
                            );
                        }
                    }
                    let in_overlap_range = |position: usize| {
                        overlap_ranges.iter().any(|&(low, high)| {
                            low <= position as isize && position as isize <= high
                        })
                    };
                    let window_start = sources.window_start(NodeKind::ForwardStop, min);
                    for &(low, high) in &overlap_ranges {
                        if low > high {
                            continue;
                        }
                        let list = &sources.kinds[NodeKind::ForwardStop as usize];
                        let first = window_start
                            + list.positions[window_start..]
                                .partition_point(|&position| (position as isize) < low);
                        for (&id, &position) in
                            list.ids[first..].iter().zip(&list.positions[first..])
                        {
                            if position as isize > high {
                                break;
                            }
                            connect!(id as usize, ConnectionType::TripleOverlap);
                        }
                    }
                    sources.admit(TRIPLE_OVERLAP_WINDOW, min, |position| {
                        position + 2 * STOP_CODON_OFFSET < target.index
                    });
                    if let Some(position) = sources.window_best(TRIPLE_OVERLAP_WINDOW, min) {
                        let list = &sources.kinds[NodeKind::ForwardStop as usize];
                        // The window covers the sources in the ranges too.
                        // If its best is one of them, the best of the others
                        // is found by looking at each.
                        let position = if in_overlap_range(list.positions[position]) {
                            let admitted = sources.windows[TRIPLE_OVERLAP_WINDOW].admitted;
                            let mut best_outside: Option<usize> = None;
                            for candidate in window_start..admitted {
                                if in_overlap_range(list.positions[candidate]) {
                                    continue;
                                }
                                let key = list.keys[candidate];
                                let better = match best_outside {
                                    Some(current) => key >= list.keys[current],
                                    None => !key.is_nan(),
                                };
                                if better {
                                    best_outside = Some(candidate);
                                }
                            }
                            best_outside
                        } else {
                            Some(position)
                        };
                        if let Some(position) = position {
                            best.offer(
                                list.ids[position] as usize,
                                ConnectionCandidate {
                                    total_score: list.keys[position],
                                    max_frame: None,
                                },
                            );
                        }
                    }

                    // 3'rev->3'rev rejects sources whose stop value is at or
                    // before the target.
                    for live in &mut sources.live_reverse_stops {
                        live.retain(|&id| {
                            let id = id as usize;
                            if dp_nodes[id].stop_value <= target.index as isize {
                                return false;
                            }
                            if id >= min {
                                connect!(id, ConnectionType::ReverseOperon);
                            }
                            true
                        });
                    }
                }
                NodeKind::ReverseStart => {
                    // 3'rev->5'rev rejects sources whose stop value is at or
                    // before the target, and stops of the other frames.
                    sources.live_reverse_stops[target.index % 3].retain(|&id| {
                        let id = id as usize;
                        if dp_nodes[id].stop_value <= target.index as isize {
                            return false;
                        }
                        if id >= min {
                            connect!(id, ConnectionType::ReverseGene);
                        }
                        true
                    });

                    // 3'fwd->5'rev rejects stops outside
                    // [stop_value - 3, stop_value + MAXIMUM_OPPOSITE_OVERLAP - 6].
                    let lowest = target.stop_value - 3;
                    let highest = target.stop_value + MAXIMUM_OPPOSITE_OVERLAP as isize - 6;
                    let window_start = sources.window_start(NodeKind::ForwardStop, min);
                    let list = &sources.kinds[NodeKind::ForwardStop as usize];
                    let first = window_start
                        + list.positions[window_start..]
                            .partition_point(|&position| (position as isize) < lowest);
                    for (&id, &position) in list.ids[first..].iter().zip(&list.positions[first..]) {
                        if position as isize > highest {
                            break;
                        }
                        connect!(id as usize, ConnectionType::OverlappingOpposite);
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
                node.index,
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
const TRIPLE_OVERLAP_WINDOW: usize = 3;

/// The candidate sources of one kind, in increasing node order: for each,
/// its node index, its position in the sequence and its key in the
/// constant-score windows (`total_score + constant`).
#[derive(Default)]
struct KindSources {
    ids: Vec<u32>,
    positions: Vec<usize>,
    keys: Vec<f64>,
    /// List position of the first source inside the distance window it was
    /// last asked for ([`SourceLists::window_start`]).
    window_start: usize,
}

/// Candidate sources, grouped by kind, together with the views of them that
/// the connection types read and the sliding-window maxima for the
/// constant-score connection types.
struct SourceLists {
    kinds: [KindSources; 4],
    /// Forward starts by reading frame, in increasing node order.
    forward_starts_by_frame: [Vec<u32>; 3],
    /// Reverse stops by reading frame, in increasing node order, without
    /// those whose stop value a target has reached: targets come in position
    /// order, so no later target can connect to them as a reverse gene or a
    /// reverse operon.
    live_reverse_stops: [Vec<u32>; 3],
    windows: [ConstantWindow; 4],
}

impl Default for SourceLists {
    fn default() -> Self {
        Self {
            kinds: Default::default(),
            forward_starts_by_frame: Default::default(),
            live_reverse_stops: Default::default(),
            windows: [
                ConstantWindow::new(NodeKind::ForwardStop),
                ConstantWindow::new(NodeKind::ReverseStart),
                ConstantWindow::new(NodeKind::ReverseStart),
                ConstantWindow::new(NodeKind::ForwardStop),
            ],
        }
    }
}

impl SourceLists {
    /// Empties the lists, keeping their buffers.
    fn clear(&mut self) {
        for list in &mut self.kinds {
            list.ids.clear();
            list.positions.clear();
            list.keys.clear();
            list.window_start = 0;
        }
        for list in self
            .forward_starts_by_frame
            .iter_mut()
            .chain(&mut self.live_reverse_stops)
        {
            list.clear();
        }
        self.reset_windows();
    }

    #[inline]
    fn push(&mut self, kind: NodeKind, position: usize, index: usize, constant_total: f64) {
        let list = &mut self.kinds[kind as usize];
        list.ids.push(index as u32);
        list.positions.push(position);
        list.keys.push(constant_total);
        match kind {
            NodeKind::ForwardStart => self.forward_starts_by_frame[position % 3].push(index as u32),
            NodeKind::ReverseStop => self.live_reverse_stops[position % 3].push(index as u32),
            NodeKind::ForwardStop | NodeKind::ReverseStart => {}
        }
    }

    /// List position of the first source of `kind` with node index `>= min`.
    #[inline]
    fn window_start(&mut self, kind: NodeKind, min: usize) -> usize {
        let list = &mut self.kinds[kind as usize];
        let ids = &list.ids;
        let mut start = list.window_start;
        if start > 0 && ids[start - 1] as usize >= min {
            // The window moved backwards (the Prodigal giant-ORF exception).
            start = ids[..start].partition_point(|&id| (id as usize) < min);
        } else {
            // It usually moves forward by a few sources at most.
            let mut steps = 0;
            while start < ids.len() && (ids[start] as usize) < min {
                start += 1;
                steps += 1;
                if steps == 8 {
                    start += ids[start..].partition_point(|&id| (id as usize) < min);
                    break;
                }
            }
        }
        list.window_start = start;
        start
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

    /// Admits the sources of the window's kind whose position satisfies
    /// `eligible` (a predicate that, once false for a source, is false for
    /// every later source) into the sliding-window maximum.
    #[inline]
    fn admit(&mut self, window: usize, min: usize, eligible: impl Fn(usize) -> bool) {
        let kind = self.windows[window].kind;
        let start = self.window_start(kind, min);
        self.windows[window].admit(&self.kinds[kind as usize], start, eligible);
    }

    /// List position of the best source in the sliding window, among the
    /// admitted sources with node index `>= min`.
    #[inline]
    fn window_best(&mut self, window: usize, min: usize) -> Option<usize> {
        let window = &mut self.windows[window];
        window.best(&self.kinds[window.kind as usize], min)
    }

    /// The sources a window slides over.
    #[inline]
    fn window_sources(&self, window: usize) -> &KindSources {
        &self.kinds[self.windows[window].kind as usize]
    }
}

/// Sliding-window maximum of `source.total_score + constant` over the sources
/// of one kind, for connection types whose score is that constant.
///
/// The deque holds list positions with strictly decreasing keys; a new
/// source evicts every earlier source with a key `<=` its own, so that among
/// equal keys the latest source wins, matching the `>=` rule of the
/// sequential scan. Sources that have left the distance window are dropped
/// from the front when the best is asked for.
struct ConstantWindow {
    kind: NodeKind,
    deque: VecDeque<u32>,
    admitted: usize,
}

impl ConstantWindow {
    fn new(kind: NodeKind) -> Self {
        Self {
            kind,
            deque: VecDeque::new(),
            admitted: 0,
        }
    }

    /// Admits the eligible sources of `list` that have not been admitted
    /// yet, starting no earlier than list position `start`, the first source
    /// inside the distance window.
    #[inline]
    fn admit(&mut self, list: &KindSources, start: usize, eligible: impl Fn(usize) -> bool) {
        if self.admitted < start {
            self.admitted = start;
        }
        while self.admitted < list.ids.len() {
            let position = self.admitted;
            if !eligible(list.positions[position]) {
                break;
            }
            self.admitted += 1;
            let key = list.keys[position];
            // A NaN total can never be selected by the sequential scan.
            if key.is_nan() {
                continue;
            }
            while let Some(&back) = self.deque.back() {
                if list.keys[back as usize] <= key {
                    self.deque.pop_back();
                } else {
                    break;
                }
            }
            self.deque.push_back(position as u32);
        }
    }

    #[inline]
    fn best(&mut self, list: &KindSources, min: usize) -> Option<usize> {
        while let Some(&front) = self.deque.front() {
            if (list.ids[front as usize] as usize) < min {
                self.deque.pop_front();
            } else {
                return Some(front as usize);
            }
        }
        None
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

fn is_terminal(node: &DpNode) -> bool {
    (node.strand == Strand::Forward && node.is_stop)
        || (node.strand == Strand::Reverse && !node.is_stop)
}

/// The terminal node (forward stop or reverse start) with the highest total
/// score; the last one among equals.
fn find_maximum_terminal_node(nodes: &[DpNode]) -> Option<usize> {
    nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| is_terminal(n))
        .max_by(|(_, a), (_, b)| a.total_score.partial_cmp(&b.total_score).unwrap())
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
    use crate::types::{CodonType, NodeMotifInfo, NodePosition, NodeScores, NodeState};

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

    fn dp_nodes_of(nodes: &[Node]) -> Vec<DpNode> {
        let training = create_test_training();
        nodes
            .iter()
            .map(|node| DpNode::from_node(node, &training))
            .collect()
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

    /// Deterministic pseudo-random numbers for the reference comparison.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }

        /// A value in `[low, high)`.
        fn float(&mut self, low: f64, high: f64) -> f64 {
            low + (high - low) * (self.next() % 1_000_000) as f64 / 1_000_000.0
        }
    }

    /// A sequence whose stretches differ in GC content: stop codons are rare
    /// in the GC-rich ones, which gives long, overlapping open reading
    /// frames on both strands, some of them longer than the node window.
    fn mixed_gc_sequence(seed: u64, length: usize) -> Vec<u8> {
        let mut random = Lcg(seed);
        let mut sequence = Vec::with_capacity(length);
        while sequence.len() < length {
            let gc_percent = [45, 60, 75, 90][(random.next() % 4) as usize];
            let stretch = 500 + (random.next() % 6000) as usize;
            for _ in 0..stretch.min(length - sequence.len()) {
                let strong = random.next() % 100 < gc_percent;
                let first = random.next().is_multiple_of(2);
                sequence.push(match (strong, first) {
                    (true, true) => b'G',
                    (true, false) => b'C',
                    (false, true) => b'A',
                    (false, false) => b'T',
                });
            }
        }

        // One open reading frame per strand without a stop codon for 15 kb.
        for codon in sequence[5_000..20_000].chunks_exact_mut(3) {
            if matches!(&*codon, b"TAA" | b"TAG" | b"TGA") {
                codon[0] = b'C';
            }
        }
        for codon in sequence[30_000..45_000].chunks_exact_mut(3) {
            if matches!(&*codon, b"TTA" | b"CTA" | b"TCA") {
                codon[2] = b'G';
            }
        }
        sequence
    }

    /// Sorted nodes of a mixed sequence with pseudo-random scores and the
    /// start pointers of `record_overlapping_starts`.
    fn scored_nodes(
        seed: u64,
        circular: bool,
        training: &Training,
        is_final_pass: bool,
    ) -> Vec<Node> {
        use crate::node::{add_nodes, record_overlapping_starts, sort_nodes_by_position};
        use crate::sequence::encoded::EncodedSequence;

        let encoded = EncodedSequence::without_masking(&mixed_gc_sequence(seed, 60_000));
        let mut nodes = Vec::new();
        add_nodes(&encoded, &mut nodes, false, circular, training).unwrap();
        sort_nodes_by_position(&mut nodes);

        let mut random = Lcg(seed ^ 0x9e37_79b9);
        for node in nodes.iter_mut() {
            if node.position.codon_type == CodonType::Stop {
                continue;
            }
            node.scores.coding_score = random.float(-30.0, 60.0);
            node.scores.ribosome_binding_score = random.float(-6.0, 6.0);
            node.scores.upstream_score = random.float(-3.0, 3.0);
            node.scores.type_score = random.float(-4.0, 4.0);
            node.scores.start_score = node.scores.ribosome_binding_score
                + node.scores.upstream_score
                + node.scores.type_score;
            node.scores.gc_frame_scores = [
                random.float(0.0, 1.5),
                random.float(0.0, 1.5),
                random.float(0.0, 1.5),
            ];
        }
        record_overlapping_starts(&mut nodes, training, is_final_pass);
        nodes
    }

    /// The dynamic programming as Prodigal runs it: every target looks at
    /// every node of its window, in order, and takes a connection that is at
    /// least as good as the best so far.
    fn exhaustive_dynamic_programming(
        nodes: &mut [Node],
        training: &Training,
        is_final_pass: bool,
    ) {
        use crate::algorithms::dynamic_programming::scoring::connection_types::determine_connection_type;

        let mut dp_nodes = dp_nodes_with(nodes, training);
        for target in 0..dp_nodes.len() {
            let min = calculate_distance_constraints(&dp_nodes, target);
            let target_scores = TargetScores::for_target(&dp_nodes, target, training);
            for source in min..target {
                if dp_nodes[source].has_edge_artifacts() {
                    continue;
                }
                let connection_type = determine_connection_type(&nodes[source], &nodes[target]);
                let Some(candidate) = evaluate_connection(
                    &dp_nodes,
                    source,
                    target,
                    training,
                    is_final_pass,
                    connection_type,
                    &target_scores,
                ) else {
                    continue;
                };
                if candidate.total_score >= dp_nodes[target].total_score {
                    dp_nodes[target].total_score = candidate.total_score;
                    dp_nodes[target].traceback = source as u32;
                    nodes[target].scores.total_score = candidate.total_score;
                    nodes[target].state.traceback = Some(source);
                    nodes[target].state.overlap_marker = candidate.max_frame;
                }
            }
        }
    }

    fn dp_nodes_with(nodes: &[Node], training: &Training) -> Vec<DpNode> {
        nodes
            .iter()
            .map(|node| DpNode::from_node(node, training))
            .collect()
    }

    #[test]
    fn dynamic_programming_matches_the_exhaustive_scan() {
        let training = Training {
            gc_bias_factors: [1.3, 0.4, 1.1],
            ..create_test_training()
        };

        let mut workspace = DpWorkspace::default();
        for (seed, circular) in [(1, false), (2, false), (3, true), (4, false), (5, true)] {
            for is_final_pass in [false, true] {
                let mut expected = scored_nodes(seed, circular, &training, is_final_pass);
                let mut actual = expected.clone();
                assert!(expected.len() > 2000, "too few nodes: {}", expected.len());

                initialize_nodes(&mut expected);
                exhaustive_dynamic_programming(&mut expected, &training, is_final_pass);
                initialize_nodes(&mut actual);
                // One workspace for every run, as the metagenomic mode keeps it.
                compute_dynamic_programming_scores(
                    &mut actual,
                    &training,
                    is_final_pass,
                    &mut workspace,
                );

                let mut connected = 0;
                let mut overlapping = 0;
                let mut whole_prefix = 0;
                let dp_nodes = dp_nodes_with(&actual, &training);
                for (index, (expected, actual)) in expected.iter().zip(&actual).enumerate() {
                    let context = format!(
                        "seed {seed}, circular {circular}, final pass {is_final_pass}, node {index}"
                    );
                    assert_eq!(
                        expected.scores.total_score.to_bits(),
                        actual.scores.total_score.to_bits(),
                        "total score, {context}"
                    );
                    assert_eq!(
                        expected.state.traceback, actual.state.traceback,
                        "traceback, {context}"
                    );
                    assert_eq!(
                        expected.state.overlap_marker, actual.state.overlap_marker,
                        "overlap marker, {context}"
                    );
                    connected += usize::from(actual.state.traceback.is_some());
                    overlapping += usize::from(actual.state.overlap_marker.is_some());
                    whole_prefix += usize::from(
                        index > 2 * MAX_NODE_DIST
                            && calculate_distance_constraints(&dp_nodes, index) == 0,
                    );
                }
                // The comparison is only worth something if the nodes
                // exercise the connections it is about.
                assert!(
                    connected > expected.len() / 2,
                    "{connected} connected nodes"
                );
                assert!(overlapping > 0, "no overlapping genes");
                assert!(whole_prefix > 0, "no node looks back to the first node");
            }
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

        let max_idx = find_maximum_terminal_node(&dp_nodes_of(&nodes));
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

        let max_idx = find_maximum_terminal_node(&dp_nodes_of(&nodes));
        assert_eq!(max_idx, Some(2)); // Should pick the highest scoring reverse start
    }

    #[test]
    fn test_find_maximum_terminal_node_no_terminals() {
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false),
            create_test_node(200, Strand::Reverse, CodonType::Stop, 100, false),
        ];

        let max_idx = find_maximum_terminal_node(&dp_nodes_of(&nodes));
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

        let max_idx = find_maximum_terminal_node(&dp_nodes_of(&nodes));
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

        compute_dynamic_programming_scores(
            &mut nodes,
            &training,
            true,
            &mut DpWorkspace::default(),
        );
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

        compute_dynamic_programming_scores(
            &mut nodes,
            &training,
            true,
            &mut DpWorkspace::default(),
        );

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
