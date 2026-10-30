//! The nodes of each strand as two disjoint views of one node slice.
//!
//! Several passes over the nodes walk one strand at a time and read and
//! write that strand's nodes only, so the two strands can be processed on
//! two threads. The nodes of the two strands are interleaved in the slice,
//! so [`split_by_strand`] gathers the references to each strand's nodes in a
//! vector of its own; the two vectors refer to disjoint nodes, so each can
//! be used on its own thread.

use bio::bio_types::strand::Strand;

use crate::types::{CodonType, Node};

/// The nodes of one strand within a slice of nodes, in slice order.
///
/// Created by [`split_by_strand`], which hands out one view per strand.
/// The views refer to disjoint nodes, so each can be used on its own thread
/// while the other is in use.
pub struct StrandNodes<'a>(Vec<&'a mut Node>);

impl StrandNodes<'_> {
    /// The nodes of the strand, in slice order.
    pub fn iter(&self) -> impl DoubleEndedIterator<Item = &Node> + ExactSizeIterator {
        self.0.iter().map(|node| &**node)
    }

    /// The nodes of the strand, mutably, in slice order.
    pub fn iter_mut(&mut self) -> impl DoubleEndedIterator<Item = &mut Node> + ExactSizeIterator {
        self.0.iter_mut().map(|node| &mut **node)
    }
}

/// Split `nodes` into the nodes of the forward strand and those of the
/// reverse strand, each in slice order. Nodes of an unknown strand are in
/// neither view.
pub fn split_by_strand(nodes: &mut [Node]) -> (StrandNodes<'_>, StrandNodes<'_>) {
    let mut forward = Vec::with_capacity(nodes.len() / 2 + 1);
    let mut reverse = Vec::with_capacity(nodes.len() / 2 + 1);
    for node in nodes.iter_mut() {
        match node.position.strand {
            Strand::Forward => forward.push(node),
            Strand::Reverse => reverse.push(node),
            Strand::Unknown => {}
        }
    }
    (StrandNodes(forward), StrandNodes(reverse))
}

/// Whether a start node of `nodes` comes, in slice order, before every stop
/// node of its frame.
///
/// The walks over the nodes that reset their per-frame state at each stop
/// node start the reverse strand from the state the forward strand's walk
/// left behind, as Prodigal does; such a start node is the only node that
/// sees it. `add_nodes` creates one only for a gene that wraps around a
/// circular sequence.
pub fn has_start_before_stop(nodes: &StrandNodes) -> bool {
    let mut stop_seen = [false; 3];
    for node in nodes.iter() {
        let frame = node.position.index % 3;
        if node.position.codon_type == CodonType::Stop {
            stop_seen[frame] = true;
            if stop_seen.iter().all(|&seen| seen) {
                return false;
            }
        } else if !stop_seen[frame] {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::NodePosition;

    fn node(index: usize, strand: Strand) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type: CodonType::Atg,
                stop_value: 0,
                is_edge: false,
            },
            ..Default::default()
        }
    }

    fn mixed_nodes() -> Vec<Node> {
        vec![
            node(0, Strand::Forward),
            node(1, Strand::Reverse),
            node(2, Strand::Reverse),
            node(3, Strand::Unknown),
            node(4, Strand::Forward),
            node(5, Strand::Reverse),
        ]
    }

    #[test]
    fn views_hold_each_strand_in_order() {
        let mut nodes = mixed_nodes();
        let (mut forward, mut reverse) = split_by_strand(&mut nodes);
        let positions: Vec<usize> = forward.iter_mut().map(|node| node.position.index).collect();
        assert_eq!(positions, [0, 4]);
        let positions: Vec<usize> = reverse.iter_mut().map(|node| node.position.index).collect();
        assert_eq!(positions, [1, 2, 5]);
        let reversed: Vec<usize> = reverse
            .iter_mut()
            .rev()
            .map(|node| node.position.index)
            .collect();
        assert_eq!(reversed, [5, 2, 1]);
        assert_eq!(forward.iter_mut().len(), 2);
    }

    #[test]
    fn writes_through_the_views_reach_the_slice() {
        let mut nodes = mixed_nodes();
        let (mut forward, mut reverse) = split_by_strand(&mut nodes);
        for node in forward.iter_mut() {
            node.scores.coding_score = 1.0;
        }
        for node in reverse.iter_mut().rev() {
            node.scores.coding_score = -1.0;
        }
        let scores: Vec<f64> = nodes.iter().map(|node| node.scores.coding_score).collect();
        assert_eq!(scores, [1.0, -1.0, -1.0, 0.0, 1.0, -1.0]);
    }

    #[test]
    fn views_can_be_written_on_two_threads() {
        let mut nodes: Vec<Node> = (0..10_000)
            .map(|index| {
                node(
                    index,
                    if index % 3 == 0 {
                        Strand::Forward
                    } else {
                        Strand::Reverse
                    },
                )
            })
            .collect();
        let (mut forward, mut reverse) = split_by_strand(&mut nodes);
        std::thread::scope(|scope| {
            scope.spawn(|| {
                for node in forward.iter_mut() {
                    node.scores.coding_score = node.position.index as f64;
                }
            });
            scope.spawn(|| {
                for node in reverse.iter_mut() {
                    node.scores.coding_score = -(node.position.index as f64);
                }
            });
        });
        for (index, node) in nodes.iter().enumerate() {
            let expected = if index % 3 == 0 {
                index as f64
            } else {
                -(index as f64)
            };
            assert_eq!(node.scores.coding_score, expected);
        }
    }

    #[test]
    fn start_before_stop_is_found_per_frame() {
        let stop = |index: usize| {
            let mut node = node(index, Strand::Reverse);
            node.position.codon_type = CodonType::Stop;
            node
        };
        // Stops open every frame before any start: nothing to find.
        let mut nodes = vec![stop(0), stop(1), stop(2), node(3, Strand::Reverse)];
        let (_, reverse) = split_by_strand(&mut nodes);
        assert!(!has_start_before_stop(&reverse));

        // A start of frame 2 comes before the stop of frame 2.
        let mut nodes = vec![stop(0), stop(1), node(5, Strand::Reverse), stop(8)];
        let (_, reverse) = split_by_strand(&mut nodes);
        assert!(has_start_before_stop(&reverse));

        // A frame without any stop at all counts too, and the forward
        // strand's nodes are not looked at.
        let mut nodes = vec![stop(0), node(1, Strand::Forward), node(4, Strand::Reverse)];
        let (forward, reverse) = split_by_strand(&mut nodes);
        assert!(has_start_before_stop(&reverse));
        assert!(has_start_before_stop(&forward));
        let mut nodes: Vec<Node> = Vec::new();
        let (_, reverse) = split_by_strand(&mut nodes);
        assert!(!has_start_before_stop(&reverse));
    }

    #[test]
    fn empty_slice_gives_empty_views() {
        let mut nodes: Vec<Node> = Vec::new();
        let (mut forward, mut reverse) = split_by_strand(&mut nodes);
        assert_eq!(forward.iter_mut().count(), 0);
        assert_eq!(reverse.iter_mut().count(), 0);
    }
}
