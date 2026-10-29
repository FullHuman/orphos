use bio::bio_types::strand::Strand;

use crate::types::{Motif, Node, NodeState};

/// Reset all scoring data for nodes to default values
///
/// Clears all calculated scores, motif information, and dynamic programming
/// state to prepare nodes for re-scoring with different parameters.
pub fn reset_node_scores(nodes: &mut [Node]) {
    for node in nodes {
        node.scores.gc_frame_scores.fill(0.0);
        node.motif_info.ribosome_binding_sites.fill(0);

        node.scores.total_score = 0.0;
        node.scores.coding_score = 0.0;
        node.scores.start_score = 0.0;
        node.scores.ribosome_binding_score = 0.0;
        node.scores.type_score = 0.0;
        node.scores.upstream_score = 0.0;

        node.motif_info.best_motif = Motif::default();
        node.state = NodeState::default();
    }
}

/// Sort nodes by genomic position with strand ordering
///
/// Sorts nodes first by sequence position, then by strand (forward before reverse).
/// This ordering is critical for dynamic programming algorithms that process
/// nodes in genomic order.
pub fn sort_nodes_by_position(nodes: &mut [Node]) {
    // Nodes are large, so sort a compact key array and then move each node
    // once, in place, instead of moving them while sorting or gathering them
    // into a second array as large as the first. Keys are unique for nodes
    // created by `add_nodes`, which makes the order fully determined.
    let mut keys: Vec<(u64, u32)> = nodes
        .iter()
        .enumerate()
        .map(|(original, node)| {
            let strand_rank = match node.position.strand {
                Strand::Forward => 0,
                Strand::Unknown => 1,
                Strand::Reverse => 2,
            };
            (
                (node.position.index as u64) << 2 | strand_rank,
                original as u32,
            )
        })
        .collect();
    keys.sort_unstable();

    // `source[i]` is the original index of the node that belongs at `i`.
    // Each cycle of the permutation is followed from its first position,
    // moving every node of the cycle once; `source[i] = i` marks a position
    // as done.
    let mut source: Vec<u32> = keys.into_iter().map(|(_, original)| original).collect();
    for start in 0..source.len() {
        if source[start] as usize == start {
            continue;
        }
        let displaced = std::mem::take(&mut nodes[start]);
        let mut position = start;
        loop {
            let from = source[position] as usize;
            source[position] = position as u32;
            if from == start {
                nodes[position] = displaced;
                break;
            }
            nodes[position] = std::mem::take(&mut nodes[from]);
            position = from;
        }
    }
}

/// Return nodes to the state `add_nodes` created them in: their positions
/// stay, every score, motif and dynamic programming field goes back to its
/// default. Nodes sorted once can then be scored again, as for the gene
/// prediction that follows training on the same sequence, without being
/// created and sorted a second time.
pub fn reset_nodes_to_pristine(nodes: &mut [Node]) {
    for node in nodes {
        node.scores = Default::default();
        node.state = Default::default();
        node.motif_info = Default::default();
    }
}
#[cfg(test)]
mod tests {
    use crate::types::{NodeMotifInfo, NodePosition, NodeScores};

    use super::*;

    fn create_test_node(index: usize, strand: Strand) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type: crate::types::CodonType::Atg,
                stop_value: 2,
                is_edge: true,
            },
            scores: NodeScores {
                gc_frame_scores: [1.0, 2.0, 3.0],
                total_score: 10.0,
                coding_score: 5.0,
                start_score: 2.0,
                ribosome_binding_score: 1.0,
                type_score: 1.5,
                upstream_score: 0.5,
                gc_content: 0.4,
            },
            motif_info: NodeMotifInfo {
                ribosome_binding_sites: [1, 5],
                best_motif: Motif {
                    index: 10,
                    length: 20,
                    space_index: 0,
                    spacer: 23,
                    score: 0.8,
                },
            },
            state: NodeState {
                start_pointers: [None, None, None],
                traceback: None,
                trace_forward: None,
                overlap_marker: None,
                is_eliminated: true,
                gc_bias_frame: 10,
            },
        }
    }

    #[test]
    fn test_reset_node_scores_single_node() {
        let mut nodes = vec![create_test_node(100, Strand::Forward)];

        reset_node_scores(&mut nodes);

        let node = &nodes[0];
        assert_eq!(node.scores.gc_frame_scores, [0.0, 0.0, 0.0]);
        assert_eq!(node.motif_info.ribosome_binding_sites, [0, 0]);
        assert_eq!(node.scores.total_score, 0.0);
        assert_eq!(node.scores.coding_score, 0.0);
        assert_eq!(node.scores.start_score, 0.0);
        assert_eq!(node.scores.ribosome_binding_score, 0.0);
        assert_eq!(node.scores.type_score, 0.0);
        assert_eq!(node.scores.upstream_score, 0.0);
        assert_eq!(node.motif_info.best_motif, Motif::default());
        assert_eq!(node.state, NodeState::default());
    }

    #[test]
    fn test_reset_node_scores_multiple_nodes() {
        let mut nodes = vec![
            create_test_node(50, Strand::Forward),
            create_test_node(100, Strand::Reverse),
            create_test_node(75, Strand::Forward),
        ];

        reset_node_scores(&mut nodes);

        for node in &nodes {
            assert_eq!(node.scores.total_score, 0.0);
            assert_eq!(node.scores.gc_frame_scores, [0.0, 0.0, 0.0]);
            assert_eq!(node.motif_info.ribosome_binding_sites, [0, 0]);
        }
    }

    #[test]
    fn test_reset_node_scores_empty_slice() {
        let mut nodes: Vec<Node> = vec![];
        reset_node_scores(&mut nodes);
        assert!(nodes.is_empty());
    }

    #[test]
    fn pristine_nodes_match_freshly_created_ones() {
        use crate::node::add_nodes;
        use crate::sequence::encoded::EncodedSequence;
        use crate::types::Training;

        let sequence: Vec<u8> =
            b"ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGATAG"
                .iter()
                .cycle()
                .take(4000)
                .copied()
                .collect();
        let encoded = EncodedSequence::with_masking(&sequence);
        let training = Training {
            translation_table: 11,
            ..Default::default()
        };
        let mut nodes = Vec::new();
        add_nodes(&encoded, &mut nodes, true, false, &training).unwrap();
        sort_nodes_by_position(&mut nodes);
        let fresh = format!("{nodes:?}");

        for (index, node) in nodes.iter_mut().enumerate() {
            node.scores.coding_score = index as f64;
            node.scores.gc_frame_scores = [1.0, 2.0, 3.0];
            node.state.traceback = Some(index);
            node.state.gc_bias_frame = 2;
            node.motif_info.ribosome_binding_sites = [3, 4];
            node.motif_info.best_motif.score = 5.0;
        }
        reset_nodes_to_pristine(&mut nodes);

        assert_eq!(format!("{nodes:?}"), fresh);
    }

    #[test]
    fn test_sort_nodes_by_position_basic() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward),
            create_test_node(50, Strand::Forward),
            create_test_node(75, Strand::Forward),
        ];

        sort_nodes_by_position(&mut nodes);

        assert_eq!(nodes[0].position.index, 50);
        assert_eq!(nodes[1].position.index, 75);
        assert_eq!(nodes[2].position.index, 100);
    }

    #[test]
    fn test_sort_nodes_strand_priority() {
        let mut nodes = vec![
            create_test_node(100, Strand::Reverse),
            create_test_node(100, Strand::Forward),
        ];

        sort_nodes_by_position(&mut nodes);

        assert_eq!(nodes[0].position.strand, Strand::Forward);
        assert_eq!(nodes[1].position.strand, Strand::Reverse);
    }

    #[test]
    fn test_sort_nodes_mixed_positions_and_strands() {
        let mut nodes = vec![
            create_test_node(200, Strand::Reverse),
            create_test_node(100, Strand::Reverse),
            create_test_node(100, Strand::Forward),
            create_test_node(150, Strand::Forward),
        ];

        sort_nodes_by_position(&mut nodes);

        assert_eq!(nodes[0].position.index, 100);
        assert_eq!(nodes[0].position.strand, Strand::Forward);
        assert_eq!(nodes[1].position.index, 100);
        assert_eq!(nodes[1].position.strand, Strand::Reverse);
        assert_eq!(nodes[2].position.index, 150);
        assert_eq!(nodes[3].position.index, 200);
    }

    #[test]
    fn test_sort_nodes_empty_slice() {
        let mut nodes: Vec<Node> = vec![];
        sort_nodes_by_position(&mut nodes);
        assert!(nodes.is_empty());
    }

    #[test]
    fn test_sort_nodes_single_node() {
        let mut nodes = vec![create_test_node(42, Strand::Forward)];
        sort_nodes_by_position(&mut nodes);
        assert_eq!(nodes.len(), 1);
        assert_eq!(nodes[0].position.index, 42);
    }

    #[test]
    fn test_sort_nodes_same_strand_same_position() {
        let mut nodes = vec![
            create_test_node(100, Strand::Forward),
            create_test_node(100, Strand::Forward),
        ];

        sort_nodes_by_position(&mut nodes);

        assert_eq!(nodes[0].position.index, 100);
        assert_eq!(nodes[1].position.index, 100);
        assert_eq!(nodes[0].position.strand, Strand::Forward);
        assert_eq!(nodes[1].position.strand, Strand::Forward);
    }
}
