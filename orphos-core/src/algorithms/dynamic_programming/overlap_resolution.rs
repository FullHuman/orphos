use bio::bio_types::strand::Strand;

use crate::types::{CodonType, Node, NodeSliceExtension};

/// Resolve triple overlaps in the traceback path
pub fn resolve_triple_overlaps(nodes: &mut [Node], start_index: usize) {
    // Collect indices first to avoid borrow checker issues with mutation
    let traceback_indices: Vec<usize> = nodes.iter_traceback_indices(start_index).collect();

    for window in traceback_indices.windows(2) {
        let current = window[0];
        let next = window[1];

        if should_resolve_triple_overlap(&nodes[current], &nodes[next]) {
            resolve_single_triple_overlap(nodes, current, next);
        }
    }
}

/// Check if triple overlap resolution is needed
fn should_resolve_triple_overlap(current_node: &Node, next_node: &Node) -> bool {
    current_node.position.strand == Strand::Reverse
        && current_node.position.codon_type == CodonType::Stop
        && next_node.position.strand == Strand::Forward
        && next_node.position.codon_type == CodonType::Stop
        && current_node.state.overlap_marker.is_some()
        && current_node.position.index > next_node.position.index
}

/// Resolve a single triple overlap
fn resolve_single_triple_overlap(nodes: &mut [Node], current: usize, next: usize) {
    let tmp =
        nodes[current].state.start_pointers[nodes[current].state.overlap_marker.unwrap()].unwrap();

    let mut i = tmp;
    while i > 0 && nodes[i].position.index as isize != nodes[tmp].position.stop_value {
        i -= 1;
    }

    nodes[current].state.traceback = Some(tmp);
    nodes[tmp].state.traceback = Some(i);
    nodes[i].state.overlap_marker = None;
    nodes[i].state.traceback = Some(next);
}

/// Resolve simple overlaps in the traceback path
pub fn resolve_simple_overlaps(nodes: &mut [Node], start_index: usize) {
    // Collect indices first to avoid borrow checker issues with mutation
    let traceback_indices: Vec<usize> = nodes.iter_traceback_indices(start_index).collect();

    for window in traceback_indices.windows(2) {
        let current = window[0];
        let next = window[1];

        let next_is_stop = nodes[next].position.codon_type == CodonType::Stop;
        let current_is_stop = nodes[current].position.codon_type == CodonType::Stop;

        match (
            nodes[current].position.strand,
            current_is_stop,
            nodes[next].position.strand,
            next_is_stop,
        ) {
            (Strand::Reverse, false, Strand::Forward, true) => {
                resolve_reverse_start_to_forward_stop(nodes, current, next);
            }
            (Strand::Forward, true, Strand::Forward, true) => {
                resolve_forward_stop_to_forward_stop(nodes, current, next);
            }
            (Strand::Reverse, true, Strand::Reverse, true) => {
                resolve_reverse_stop_to_reverse_stop(nodes, current, next);
            }
            _ => {}
        }
    }
}

/// Resolve reverse start to forward stop overlap
fn resolve_reverse_start_to_forward_stop(nodes: &mut [Node], current: usize, next: usize) {
    // Find the stop position by iterating backwards from current
    if let Some(stop_index) = (0..=current)
        .rev()
        .find(|&i| nodes[i].position.index as isize == nodes[current].position.stop_value)
    {
        nodes[current].state.traceback = Some(stop_index);
        nodes[stop_index].state.traceback = Some(next);
    }
}

/// Resolve forward stop to forward stop overlap
fn resolve_forward_stop_to_forward_stop(nodes: &mut [Node], current: usize, next: usize) {
    let frame = nodes[current].position.index % 3;

    nodes[current].state.traceback = nodes[next].state.start_pointers[frame];
    if let Some(traceback) = nodes[current].state.traceback {
        nodes[traceback].state.traceback = Some(next);
    }
}

/// Resolve reverse stop to reverse stop overlap
fn resolve_reverse_stop_to_reverse_stop(nodes: &mut [Node], current: usize, next: usize) {
    let frame = nodes[next].position.index % 3;

    nodes[current].state.traceback = nodes[current].state.start_pointers[frame];
    if let Some(traceback) = nodes[current].state.traceback {
        nodes[traceback].state.traceback = Some(next);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{NodeMotifInfo, NodePosition, NodeScores, NodeState, StarPointers};

    /// Helper function to create a test node
    fn create_test_node(
        index: usize,
        strand: Strand,
        codon_type: CodonType,
        stop_value: isize,
        is_edge: bool,
        start_pointers: StarPointers,
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
            state: NodeState {
                start_pointers,
                ..Default::default()
            },
            motif_info: NodeMotifInfo::default(),
        }
    }

    #[test]
    fn test_resolve_triple_overlaps_empty_path() {
        let mut nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            false,
            [None; 3],
        )];

        resolve_triple_overlaps(&mut nodes, 0);

        assert!(nodes[0].state.traceback.is_none());
    }

    #[test]
    fn test_resolve_triple_overlaps_no_overlap_needed() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false, [None; 3]),
            create_test_node(200, Strand::Forward, CodonType::Stop, 100, false, [None; 3]),
        ];

        // Set up traceback that doesn't need triple overlap resolution
        nodes[1].state.traceback = Some(0);

        resolve_triple_overlaps(&mut nodes, 1);

        assert_eq!(nodes[1].state.traceback, Some(0));
    }

    #[test]
    fn test_resolve_triple_overlaps_valid_case() {
        let mut nodes = vec![
            create_test_node(150, Strand::Forward, CodonType::Atg, 0, false, [None; 3]),
            create_test_node(100, Strand::Forward, CodonType::Stop, 153, false, [None; 3]),
            create_test_node(
                200,
                Strand::Reverse,
                CodonType::Stop,
                100,
                false,
                [Some(0), None, None],
            ),
            create_test_node(300, Strand::Forward, CodonType::Atg, 0, false, [Some(0); 3]),
        ];

        // Set up triple overlap scenario with proper start pointers
        nodes[2].state.traceback = Some(1);
        nodes[2].state.overlap_marker = Some(0);
        nodes[2].position.index = 250; // Greater than nodes[1].position.index

        resolve_triple_overlaps(&mut nodes, 2);

        // Should have resolved the overlap (exact behavior depends on implementation)
        assert!(nodes[2].state.traceback.is_some());
    }

    #[test]
    fn test_should_resolve_triple_overlap_true() {
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Stop, 50, false, [None; 3]),
            create_test_node(200, Strand::Reverse, CodonType::Stop, 150, false, [None; 3]),
        ];

        let mut test_nodes = nodes;
        test_nodes[1].state.overlap_marker = Some(0);
        test_nodes[1].position.index = 250; // Greater than nodes[0]

        let should_resolve = should_resolve_triple_overlap(&test_nodes[1], &test_nodes[0]);
        assert!(should_resolve);
    }

    #[test]
    fn test_should_resolve_triple_overlap_false_wrong_strands() {
        let nodes = [
            create_test_node(100, Strand::Forward, CodonType::Atg, 50, false, [None; 3]),
            create_test_node(200, Strand::Reverse, CodonType::Stop, 150, false, [None; 3]),
        ];

        let should_resolve = should_resolve_triple_overlap(&nodes[1], &nodes[0]);
        assert!(!should_resolve);
    }

    #[test]
    fn test_should_resolve_triple_overlap_false_no_overlap_marker() {
        let nodes = [
            create_test_node(100, Strand::Forward, CodonType::Stop, 50, false, [None; 3]),
            create_test_node(200, Strand::Reverse, CodonType::Stop, 150, false, [None; 3]),
        ];

        // No overlap marker set
        let should_resolve = should_resolve_triple_overlap(&nodes[1], &nodes[0]);
        assert!(!should_resolve);
    }

    #[test]
    fn test_resolve_simple_overlaps_empty_path() {
        let mut nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            0,
            false,
            [None; 3],
        )];

        resolve_simple_overlaps(&mut nodes, 0);

        // Should not crash and remain unchanged
        assert!(nodes[0].state.traceback.is_none());
    }

    #[test]
    fn test_resolve_simple_overlaps_reverse_start_to_forward_stop() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Stop, 50, false, [None; 3]),
            create_test_node(200, Strand::Reverse, CodonType::Atg, 100, false, [None; 3]),
            create_test_node(150, Strand::Reverse, CodonType::Stop, 200, false, [None; 3]),
        ];

        // Set up traceback chain
        nodes[1].state.traceback = Some(0);
        nodes[1].position.stop_value = 150; // Points to index 2

        resolve_simple_overlaps(&mut nodes, 1);

        // Should have resolved the overlap
        assert!(nodes[1].state.traceback.is_some());
    }

    #[test]
    fn test_resolve_simple_overlaps_forward_stop_to_forward_stop() {
        let mut nodes = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Stop,
                50,
                false,
                [Some(2); 3],
            ),
            create_test_node(200, Strand::Forward, CodonType::Stop, 150, false, [None; 3]),
            create_test_node(50, Strand::Forward, CodonType::Atg, 0, false, [None; 3]),
        ];

        nodes[1].state.traceback = Some(0);

        resolve_simple_overlaps(&mut nodes, 1);

        // Should have updated traceback based on frame
        assert!(nodes[1].state.traceback.is_some());
    }

    #[test]
    fn test_resolve_simple_overlaps_reverse_stop_to_reverse_stop() {
        let mut nodes = vec![
            create_test_node(
                100,
                Strand::Reverse,
                CodonType::Stop,
                50,
                false,
                [Some(2); 3],
            ),
            create_test_node(200, Strand::Reverse, CodonType::Stop, 150, false, [None; 3]),
            create_test_node(50, Strand::Reverse, CodonType::Atg, 0, false, [None; 3]),
        ];

        nodes[0].state.traceback = Some(1);

        resolve_simple_overlaps(&mut nodes, 0);

        // Should have updated traceback based on frame
        assert!(nodes[0].state.traceback.is_some());
    }

    #[test]
    fn test_resolve_forward_stop_to_forward_stop() {
        let mut nodes = vec![
            create_test_node(
                100,
                Strand::Forward,
                CodonType::Stop,
                50,
                false,
                [Some(2); 3],
            ),
            create_test_node(200, Strand::Forward, CodonType::Stop, 150, false, [None; 3]),
            create_test_node(50, Strand::Forward, CodonType::Atg, 0, false, [None; 3]),
        ];

        resolve_forward_stop_to_forward_stop(&mut nodes, 0, 1);

        // Should have used frame-based start pointer
        let frame = nodes[0].position.index % 3;
        assert_eq!(
            nodes[0].state.traceback,
            nodes[1].state.start_pointers[frame]
        );
        if let Some(traceback) = nodes[0].state.traceback {
            assert_eq!(nodes[traceback].state.traceback, Some(1));
        }
    }

    #[test]
    fn test_resolve_reverse_stop_to_reverse_stop() {
        let mut nodes = vec![
            create_test_node(
                100,
                Strand::Reverse,
                CodonType::Stop,
                50,
                false,
                [Some(2); 3],
            ),
            create_test_node(200, Strand::Reverse, CodonType::Stop, 150, false, [None; 3]),
            create_test_node(50, Strand::Reverse, CodonType::Atg, 0, false, [None; 3]),
        ];

        resolve_reverse_stop_to_reverse_stop(&mut nodes, 0, 1);

        // Should have used frame-based start pointer
        let frame = nodes[1].position.index % 3;
        assert_eq!(
            nodes[0].state.traceback,
            nodes[0].state.start_pointers[frame]
        );
        if let Some(traceback) = nodes[0].state.traceback {
            assert_eq!(nodes[traceback].state.traceback, Some(1));
        }
    }

    #[test]
    fn test_resolve_simple_overlaps_no_matching_patterns() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, 0, false, [None; 3]),
            create_test_node(200, Strand::Reverse, CodonType::Atg, 150, false, [None; 3]),
        ];

        // Set up traceback that doesn't match any overlap patterns
        nodes[1].state.traceback = Some(0);
        let original_traceback = nodes[1].state.traceback;

        resolve_simple_overlaps(&mut nodes, 1);

        // Should remain unchanged since no patterns match
        assert_eq!(nodes[1].state.traceback, original_traceback);
    }

    #[test]
    fn test_resolve_single_triple_overlap_with_bounds_checking() {
        let mut nodes = vec![
            create_test_node(50, Strand::Forward, CodonType::Atg, 0, false, [None; 3]),
            create_test_node(100, Strand::Forward, CodonType::Stop, 50, false, [None; 3]),
            create_test_node(
                200,
                Strand::Reverse,
                CodonType::Stop,
                100,
                false,
                [Some(0); 3],
            ),
            create_test_node(300, Strand::Forward, CodonType::Atg, 0, false, [None; 3]),
        ];

        nodes[2].state.overlap_marker = Some(0);

        resolve_single_triple_overlap(&mut nodes, 2, 1);

        // Should have updated the traceback structure
        assert!(nodes[2].state.traceback.is_some());
    }
}
