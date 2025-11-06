use bio::bio_types::strand::Strand;

use crate::types::{CodonType, Node};

/// Represents the type of connection between two nodes in the gene prediction graph
#[derive(PartialEq, Eq, Debug)]
pub enum ConnectionType {
    /// 5' forward start -> 3' forward stop (normal gene)
    ForwardGene,
    /// 3' reverse stop -> 5' reverse start (reverse gene)
    ReverseGene,
    /// 3' forward stop -> 5' forward start (intergenic region)
    ForwardIntergenic,
    /// 5' reverse start -> 3' reverse stop (reverse intergenic region)
    ReverseIntergenic,
    /// 3' forward stop -> 3' reverse stop (triple overlap case)
    TripleOverlap,
    /// 3' forward stop -> 3' forward stop (forward operon)
    ForwardOperon,
    /// 3' reverse stop -> 3' reverse stop (reverse operon)
    ReverseOperon,
    /// 3' forward stop -> 5' reverse start (overlapping opposite strands)
    OverlappingOpposite,
    /// 5' reverse start -> 5' forward start
    FiveReverseToFiveForward,
    /// Invalid connection type
    Invalid,
}

#[must_use]
pub fn has_edge_artifacts(node: &Node) -> bool {
    if node.state.traceback.is_none()
        && node.position.strand == Strand::Forward
        && node.position.codon_type == CodonType::Stop
    {
        return true;
    }
    if node.state.traceback.is_none()
        && node.position.strand == Strand::Reverse
        && node.position.codon_type != CodonType::Stop
    {
        return true;
    }
    false
}

/// Determines the type of connection between two nodes based on their strand and codon type
pub fn determine_connection_type(node1: &Node, node2: &Node) -> ConnectionType {
    use CodonType::*;
    use Strand::*;

    let n1_strand = node1.position.strand;
    let n1_is_stop = node1.position.codon_type == Stop;
    let n1_frame = node1.position.index % 3;
    let n2_strand = node2.position.strand;
    let n2_is_stop = node2.position.codon_type == Stop;
    let n2_frame = node2.position.index % 3;

    if (n1_strand == n2_strand) && (n1_frame != n2_frame) {
        if (n1_strand == Forward) && !n1_is_stop && n2_is_stop {
            return ConnectionType::Invalid;
        }
        if (n1_strand == Reverse) && n1_is_stop && !n2_is_stop {
            return ConnectionType::Invalid;
        }
    }

    match (n1_strand, n1_is_stop, n2_strand, n2_is_stop) {
        (Forward, false, Forward, true) => ConnectionType::ForwardGene,
        (Reverse, true, Reverse, false) => ConnectionType::ReverseGene,
        (Forward, true, Forward, false) => ConnectionType::ForwardIntergenic,
        (Reverse, false, Reverse, true) => ConnectionType::ReverseIntergenic,
        (Forward, true, Reverse, true) => ConnectionType::TripleOverlap,
        (Forward, true, Forward, true) => ConnectionType::ForwardOperon,
        (Reverse, true, Reverse, true) => ConnectionType::ReverseOperon,
        (Forward, true, Reverse, false) => ConnectionType::OverlappingOpposite,
        (Reverse, false, Forward, false) => ConnectionType::FiveReverseToFiveForward,
        _ => ConnectionType::Invalid,
    }
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
        traceback: Option<usize>,
    ) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value: 0,
                is_edge: false,
            },
            scores: NodeScores::default(),
            state: NodeState {
                traceback,
                ..Default::default()
            },
            motif_info: NodeMotifInfo::default(),
        }
    }

    #[test]
    fn test_has_edge_artifacts_forward_stop_no_traceback() {
        let node = create_test_node(100, Strand::Forward, CodonType::Stop, None);
        assert!(has_edge_artifacts(&node));
    }

    #[test]
    fn test_has_edge_artifacts_reverse_start_no_traceback() {
        let node = create_test_node(100, Strand::Reverse, CodonType::Atg, None);
        assert!(has_edge_artifacts(&node));
    }

    #[test]
    fn test_has_edge_artifacts_forward_stop_with_traceback() {
        let node = create_test_node(100, Strand::Forward, CodonType::Stop, Some(0));
        assert!(!has_edge_artifacts(&node));
    }

    #[test]
    fn test_has_edge_artifacts_reverse_start_with_traceback() {
        let node = create_test_node(100, Strand::Reverse, CodonType::Atg, Some(0));
        assert!(!has_edge_artifacts(&node));
    }

    #[test]
    fn test_has_edge_artifacts_forward_start() {
        let node = create_test_node(100, Strand::Forward, CodonType::Atg, None);
        assert!(!has_edge_artifacts(&node));
    }

    #[test]
    fn test_has_edge_artifacts_reverse_stop() {
        let node = create_test_node(100, Strand::Reverse, CodonType::Stop, None);
        assert!(!has_edge_artifacts(&node));
    }

    #[test]
    fn test_determine_connection_type_forward_intergenic() {
        let stop_node = create_test_node(200, Strand::Forward, CodonType::Stop, None);
        let start_node = create_test_node(300, Strand::Forward, CodonType::Atg, None);

        let connection_type = determine_connection_type(&stop_node, &start_node);
        assert_eq!(connection_type, ConnectionType::ForwardIntergenic);
    }

    #[test]
    fn test_determine_connection_type_reverse_intergenic() {
        let start_node = create_test_node(300, Strand::Reverse, CodonType::Atg, None);
        let stop_node = create_test_node(200, Strand::Reverse, CodonType::Stop, None);

        let connection_type = determine_connection_type(&start_node, &stop_node);
        assert_eq!(connection_type, ConnectionType::ReverseIntergenic);
    }

    #[test]
    fn test_determine_connection_type_triple_overlap() {
        let forward_stop = create_test_node(200, Strand::Forward, CodonType::Stop, None);
        let reverse_stop = create_test_node(100, Strand::Reverse, CodonType::Stop, None);

        let connection_type = determine_connection_type(&forward_stop, &reverse_stop);
        assert_eq!(connection_type, ConnectionType::TripleOverlap);
    }

    #[test]
    fn test_determine_connection_type_forward_operon() {
        let stop_node1 = create_test_node(200, Strand::Forward, CodonType::Stop, None);
        let stop_node2 = create_test_node(300, Strand::Forward, CodonType::Stop, None);

        let connection_type = determine_connection_type(&stop_node1, &stop_node2);
        assert_eq!(connection_type, ConnectionType::ForwardOperon);
    }

    #[test]
    fn test_determine_connection_type_reverse_operon() {
        let stop_node1 = create_test_node(300, Strand::Reverse, CodonType::Stop, None);
        let stop_node2 = create_test_node(200, Strand::Reverse, CodonType::Stop, None);

        let connection_type = determine_connection_type(&stop_node1, &stop_node2);
        assert_eq!(connection_type, ConnectionType::ReverseOperon);
    }

    #[test]
    fn test_determine_connection_type_overlapping_opposite() {
        let forward_stop = create_test_node(200, Strand::Forward, CodonType::Stop, None);
        let reverse_start = create_test_node(150, Strand::Reverse, CodonType::Atg, None);

        let connection_type = determine_connection_type(&forward_stop, &reverse_start);
        assert_eq!(connection_type, ConnectionType::OverlappingOpposite);
    }

    #[test]
    fn test_determine_connection_type_five_reverse_to_five_forward() {
        let reverse_start = create_test_node(200, Strand::Reverse, CodonType::Atg, None);
        let forward_start = create_test_node(300, Strand::Forward, CodonType::Gtg, None);

        let connection_type = determine_connection_type(&reverse_start, &forward_start);
        assert_eq!(connection_type, ConnectionType::FiveReverseToFiveForward);
    }

    #[test]
    fn test_determine_connection_type_invalid_patterns() {
        // Test various invalid connection patterns
        let test_cases = vec![
            // Reverse stop to forward stop (invalid)
            (
                create_test_node(200, Strand::Reverse, CodonType::Stop, None),
                create_test_node(300, Strand::Forward, CodonType::Stop, None),
            ),
            // Forward start to reverse start (invalid)
            (
                create_test_node(200, Strand::Forward, CodonType::Atg, None),
                create_test_node(300, Strand::Reverse, CodonType::Gtg, None),
            ),
            // Forward start to forward start (invalid)
            (
                create_test_node(200, Strand::Forward, CodonType::Atg, None),
                create_test_node(300, Strand::Forward, CodonType::Gtg, None),
            ),
            // Reverse start to reverse start (invalid)
            (
                create_test_node(300, Strand::Reverse, CodonType::Atg, None),
                create_test_node(200, Strand::Reverse, CodonType::Gtg, None),
            ),
        ];

        for (node1, node2) in test_cases {
            let connection_type = determine_connection_type(&node1, &node2);
            assert_eq!(connection_type, ConnectionType::Invalid);
        }
    }

    #[test]
    fn test_connection_type_debug_impl() {
        // Test that all connection types can be formatted for debugging
        let connection_types = vec![
            ConnectionType::ForwardGene,
            ConnectionType::ReverseGene,
            ConnectionType::ForwardIntergenic,
            ConnectionType::ReverseIntergenic,
            ConnectionType::TripleOverlap,
            ConnectionType::ForwardOperon,
            ConnectionType::ReverseOperon,
            ConnectionType::OverlappingOpposite,
            ConnectionType::FiveReverseToFiveForward,
            ConnectionType::Invalid,
        ];

        for connection_type in connection_types {
            let debug_str = format!("{:?}", connection_type);
            assert!(!debug_str.is_empty());
        }
    }

    #[test]
    fn test_connection_type_equality() {
        assert_eq!(ConnectionType::ForwardGene, ConnectionType::ForwardGene);
        assert_ne!(ConnectionType::ForwardGene, ConnectionType::ReverseGene);
        assert_ne!(ConnectionType::Invalid, ConnectionType::ForwardGene);
    }
}
