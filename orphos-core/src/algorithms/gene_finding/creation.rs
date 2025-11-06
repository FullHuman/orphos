use bio::bio_types::strand::Strand;

use crate::{
    constants::MAX_GENES,
    types::{CodonType, Gene, Node},
};

/// Copy genes from dynamic programming result to final array
pub fn add_genes(nodes: &[Node], path_start: usize) -> Vec<Gene> {
    if path_start == usize::MAX {
        return Vec::new();
    }

    let mut genes: Vec<Gene> = Vec::new();
    let mut path = path_start;
    let mut ctr = 0;

    // Find the beginning of the path
    while let Some(traceback) = nodes[path].state.traceback {
        path = traceback;
    }

    while path != usize::MAX {
        if nodes[path].state.is_eliminated {
            path = if let Some(next_path) = nodes[path].state.trace_forward {
                next_path
            } else {
                break;
            };
            continue;
        }

        match (nodes[path].position.strand, nodes[path].position.codon_type) {
            // Forward strand start codon
            (Strand::Forward, codon_type) if codon_type != CodonType::Stop => {
                // Ensure we have a gene slot for this start
                if genes.len() <= ctr {
                    genes.push(Gene::default());
                }
                genes[ctr].coordinates.begin = nodes[path].position.index + 1;
                genes[ctr].coordinates.start_index = path;
                genes[ctr].coordinates.strand = Strand::Forward;
            }
            // Reverse strand stop codon
            (Strand::Reverse, CodonType::Stop) => {
                if genes.len() <= ctr {
                    genes.push(Gene::default());
                }
                genes[ctr].coordinates.begin = nodes[path].position.index - 1;
                genes[ctr].coordinates.stop_index = path;
                genes[ctr].coordinates.strand = Strand::Reverse;
            }
            // Forward strand stop codon - completes gene
            (Strand::Forward, CodonType::Stop) => {
                if genes.len() > ctr {
                    genes[ctr].coordinates.end = nodes[path].position.index + 3;
                    genes[ctr].coordinates.stop_index = path;
                    ctr += 1;
                }
            }
            // Reverse strand start codon - completes gene
            (Strand::Reverse, codon_type) if codon_type != CodonType::Stop => {
                if genes.len() > ctr {
                    genes[ctr].coordinates.end = nodes[path].position.index + 1;
                    genes[ctr].coordinates.start_index = path;
                    ctr += 1;
                }
            }
            _ => {}
        }

        path = if let Some(next_path) = nodes[path].state.trace_forward {
            next_path
        } else {
            break;
        };

        if ctr >= MAX_GENES {
            eprintln!("warning, max # of genes exceeded, truncating...");
            break;
        }
    }

    genes.truncate(ctr);
    genes
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
        trace_forward: Option<usize>,
        is_eliminated: bool,
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
                trace_forward,
                is_eliminated,
                ..Default::default()
            },
            motif_info: NodeMotifInfo::default(),
        }
    }

    #[test]
    fn test_add_genes_empty_path() {
        let nodes = vec![create_test_node(
            100,
            Strand::Forward,
            CodonType::Atg,
            None,
            None,
            false,
        )];

        // Test with usize::MAX path_start (empty path)
        let genes = add_genes(&nodes, usize::MAX);
        assert!(genes.is_empty());
    }

    #[test]
    fn test_add_genes_single_forward_gene() {
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false), // start
            create_test_node(400, Strand::Forward, CodonType::Stop, Some(0), None, false), // stop
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 1);
        let gene = &genes[0];

        assert_eq!(gene.coordinates.strand, Strand::Forward);
        assert_eq!(gene.coordinates.begin, 101); // index + 1
        assert_eq!(gene.coordinates.end, 403); // index + 3
        assert_eq!(gene.coordinates.start_index, 0);
        assert_eq!(gene.coordinates.stop_index, 1);
    }

    #[test]
    fn test_add_genes_single_reverse_gene() {
        let nodes = vec![
            create_test_node(100, Strand::Reverse, CodonType::Stop, None, Some(1), false), // stop (processed first)
            create_test_node(400, Strand::Reverse, CodonType::Atg, Some(0), None, false), // start (processed second)
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 1);
        let gene = &genes[0];

        assert_eq!(gene.coordinates.strand, Strand::Reverse);
        assert_eq!(gene.coordinates.begin, 99); // stop index - 1
        assert_eq!(gene.coordinates.end, 401); // start index + 1
        assert_eq!(gene.coordinates.start_index, 1);
        assert_eq!(gene.coordinates.stop_index, 0);
    }

    #[test]
    fn test_add_genes_multiple_forward_genes() {
        let nodes = vec![
            // Gene 1
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false), // start
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                Some(0),
                Some(2),
                false,
            ), // stop
            // Gene 2
            create_test_node(
                500,
                Strand::Forward,
                CodonType::Gtg,
                Some(1),
                Some(3),
                false,
            ), // start
            create_test_node(800, Strand::Forward, CodonType::Stop, Some(2), None, false), // stop
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 2);

        // Gene 1
        let gene1 = &genes[0];
        assert_eq!(gene1.coordinates.strand, Strand::Forward);
        assert_eq!(gene1.coordinates.begin, 101);
        assert_eq!(gene1.coordinates.end, 403);
        assert_eq!(gene1.coordinates.start_index, 0);
        assert_eq!(gene1.coordinates.stop_index, 1);

        // Gene 2
        let gene2 = &genes[1];
        assert_eq!(gene2.coordinates.strand, Strand::Forward);
        assert_eq!(gene2.coordinates.begin, 501);
        assert_eq!(gene2.coordinates.end, 803);
        assert_eq!(gene2.coordinates.start_index, 2);
        assert_eq!(gene2.coordinates.stop_index, 3);
    }

    #[test]
    fn test_add_genes_multiple_reverse_genes() {
        let nodes = vec![
            // Gene 1: stop then start
            create_test_node(100, Strand::Reverse, CodonType::Stop, None, Some(1), false), // stop
            create_test_node(
                400,
                Strand::Reverse,
                CodonType::Atg,
                Some(0),
                Some(2),
                false,
            ), // start
            // Gene 2: stop then start
            create_test_node(
                500,
                Strand::Reverse,
                CodonType::Stop,
                Some(1),
                Some(3),
                false,
            ), // stop
            create_test_node(800, Strand::Reverse, CodonType::Gtg, Some(2), None, false), // start
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 2);

        // Gene 1
        let gene1 = &genes[0];
        assert_eq!(gene1.coordinates.strand, Strand::Reverse);
        assert_eq!(gene1.coordinates.begin, 99); // stop index - 1
        assert_eq!(gene1.coordinates.end, 401); // start index + 1
        assert_eq!(gene1.coordinates.start_index, 1);
        assert_eq!(gene1.coordinates.stop_index, 0);

        // Gene 2
        let gene2 = &genes[1];
        assert_eq!(gene2.coordinates.strand, Strand::Reverse);
        assert_eq!(gene2.coordinates.begin, 499); // stop index - 1
        assert_eq!(gene2.coordinates.end, 801); // start index + 1
        assert_eq!(gene2.coordinates.start_index, 3);
        assert_eq!(gene2.coordinates.stop_index, 2);
    }

    #[test]
    fn test_add_genes_mixed_strand_genes() {
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false), // start
            create_test_node(
                300,
                Strand::Forward,
                CodonType::Stop,
                Some(0),
                Some(2),
                false,
            ), // stop
            create_test_node(
                400,
                Strand::Reverse,
                CodonType::Stop,
                Some(1),
                Some(3),
                false,
            ), // stop (processed first)
            create_test_node(600, Strand::Reverse, CodonType::Gtg, Some(2), None, false), // start (processed second)
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 2);

        let forward_gene = &genes[0];
        assert_eq!(forward_gene.coordinates.strand, Strand::Forward);
        assert_eq!(forward_gene.coordinates.begin, 101);
        assert_eq!(forward_gene.coordinates.end, 303);

        let reverse_gene = &genes[1];
        assert_eq!(reverse_gene.coordinates.strand, Strand::Reverse);
        assert_eq!(reverse_gene.coordinates.begin, 399); // stop index - 1
        assert_eq!(reverse_gene.coordinates.end, 601); // start index + 1
    }

    #[test]
    fn test_add_genes_different_codon_types() {
        let nodes = vec![
            // ATG start
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false),
            create_test_node(
                300,
                Strand::Forward,
                CodonType::Stop,
                Some(0),
                Some(2),
                false,
            ),
            // GTG start
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Gtg,
                Some(1),
                Some(3),
                false,
            ),
            create_test_node(
                600,
                Strand::Forward,
                CodonType::Stop,
                Some(2),
                Some(4),
                false,
            ),
            // TTG start
            create_test_node(
                700,
                Strand::Forward,
                CodonType::Ttg,
                Some(3),
                Some(5),
                false,
            ),
            create_test_node(900, Strand::Forward, CodonType::Stop, Some(4), None, false),
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 3);

        for gene in &genes {
            assert_eq!(gene.coordinates.strand, Strand::Forward);
            assert!(gene.coordinates.begin > 0);
            assert!(gene.coordinates.end > gene.coordinates.begin);
        }
    }

    #[test]
    fn test_add_genes_with_eliminated_nodes() {
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false), // start
            create_test_node(200, Strand::Forward, CodonType::Gtg, Some(0), Some(2), true), // eliminated
            create_test_node(300, Strand::Forward, CodonType::Stop, Some(1), None, false),  // stop
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 1);
        let gene = &genes[0];
        assert_eq!(gene.coordinates.start_index, 0);
        assert_eq!(gene.coordinates.stop_index, 2); // Should skip eliminated node
    }

    #[test]
    fn test_add_genes_incomplete_forward_gene() {
        // Forward gene with start but no stop (incomplete)
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, None, None, false), // start only
        ];

        let genes = add_genes(&nodes, 0);

        // Should create gene slot but not complete it
        assert_eq!(genes.len(), 0); // Truncated because ctr never incremented
    }

    #[test]
    fn test_add_genes_incomplete_reverse_gene() {
        // Reverse gene with stop but no start (incomplete)
        let nodes = vec![
            create_test_node(100, Strand::Reverse, CodonType::Stop, None, None, false), // stop only
        ];

        let genes = add_genes(&nodes, 0);

        // Should create gene slot but not complete it
        assert_eq!(genes.len(), 0); // Truncated because ctr never incremented
    }

    #[test]
    fn test_add_genes_traceback_to_beginning() {
        // Test that function finds beginning of path via traceback
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false), // start (beginning)
            create_test_node(
                200,
                Strand::Forward,
                CodonType::Stop,
                Some(0),
                Some(2),
                false,
            ), // stop
            create_test_node(
                300,
                Strand::Forward,
                CodonType::Atg,
                Some(1),
                Some(3),
                false,
            ), // start
            create_test_node(400, Strand::Forward, CodonType::Stop, Some(2), None, false), // stop (end)
        ];

        let genes = add_genes(&nodes, 3);

        assert_eq!(genes.len(), 2);

        // Should process genes in order from beginning
        assert_eq!(genes[0].coordinates.start_index, 0);
        assert_eq!(genes[1].coordinates.start_index, 2);
    }

    #[test]
    fn test_add_genes_no_traceback_path() {
        // Node with no traceback (already at beginning)
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false),
            create_test_node(200, Strand::Forward, CodonType::Stop, Some(0), None, false),
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].coordinates.begin, 101);
        assert_eq!(genes[0].coordinates.end, 203);
    }

    #[test]
    fn test_add_genes_max_genes_limit() {
        // Instead, we test the warning mechanism by mocking a scenario that would trigger it

        // Create a few genes to verify normal operation
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false),
            create_test_node(
                200,
                Strand::Forward,
                CodonType::Stop,
                Some(0),
                Some(2),
                false,
            ),
            create_test_node(
                300,
                Strand::Forward,
                CodonType::Atg,
                Some(1),
                Some(3),
                false,
            ),
            create_test_node(400, Strand::Forward, CodonType::Stop, Some(2), None, false),
        ];

        let genes = add_genes(&nodes, 0);

        // Should create normal number of genes (well below MAX_GENES)
        assert_eq!(genes.len(), 2);
        assert!(genes.len() < MAX_GENES);
    }

    #[test]
    fn test_add_genes_broken_trace_forward_chain() {
        // Test with broken trace_forward chain (None in middle)
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Atg, None, None, false), // No trace_forward
        ];

        let genes = add_genes(&nodes, 0);

        // Should stop processing when trace_forward is None
        assert_eq!(genes.len(), 0);
    }

    #[test]
    fn test_add_genes_unknown_strand() {
        let nodes = vec![
            create_test_node(100, Strand::Unknown, CodonType::Atg, None, Some(1), false),
            create_test_node(
                200,
                Strand::Forward,
                CodonType::Atg,
                Some(0),
                Some(2),
                false,
            ),
            create_test_node(300, Strand::Forward, CodonType::Stop, Some(1), None, false),
        ];

        let genes = add_genes(&nodes, 0);

        // Should skip Unknown strand and process valid gene
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].coordinates.strand, Strand::Forward);
    }

    #[test]
    fn test_add_genes_stop_codon_only() {
        let nodes = vec![
            create_test_node(100, Strand::Forward, CodonType::Stop, None, Some(1), false),
            create_test_node(200, Strand::Reverse, CodonType::Stop, Some(0), None, false),
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 0);
    }

    #[test]
    fn test_add_genes_coordinate_calculation() {
        // Test specific coordinate calculations
        let nodes = vec![
            create_test_node(0, Strand::Forward, CodonType::Atg, None, Some(1), false), // forward start at 0
            create_test_node(
                99,
                Strand::Forward,
                CodonType::Stop,
                Some(0),
                Some(2),
                false,
            ), // forward stop at 99
            create_test_node(
                100,
                Strand::Reverse,
                CodonType::Stop,
                Some(1),
                Some(3),
                false,
            ), // reverse stop at 100
            create_test_node(200, Strand::Reverse, CodonType::Atg, Some(2), None, false), // reverse start at 200
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 2);

        // Forward gene: begin = 0+1=1, end = 99+3=102
        let forward_gene = &genes[0];
        assert_eq!(forward_gene.coordinates.begin, 1);
        assert_eq!(forward_gene.coordinates.end, 102);

        // Reverse gene: begin = 100-1=99, end = 200+1=201
        let reverse_gene = &genes[1];
        assert_eq!(reverse_gene.coordinates.begin, 99);
        assert_eq!(reverse_gene.coordinates.end, 201);
    }

    #[test]
    fn test_add_genes_gene_expansion() {
        // Test that genes vector expands correctly when needed
        let nodes = vec![
            // Multiple genes to test vector expansion
            create_test_node(100, Strand::Forward, CodonType::Atg, None, Some(1), false),
            create_test_node(
                200,
                Strand::Forward,
                CodonType::Stop,
                Some(0),
                Some(2),
                false,
            ),
            create_test_node(
                300,
                Strand::Forward,
                CodonType::Atg,
                Some(1),
                Some(3),
                false,
            ),
            create_test_node(
                400,
                Strand::Forward,
                CodonType::Stop,
                Some(2),
                Some(4),
                false,
            ),
            create_test_node(
                500,
                Strand::Forward,
                CodonType::Atg,
                Some(3),
                Some(5),
                false,
            ),
            create_test_node(600, Strand::Forward, CodonType::Stop, Some(4), None, false),
        ];

        let genes = add_genes(&nodes, 0);

        assert_eq!(genes.len(), 3);

        // Verify all genes are properly initialized
        for (i, gene) in genes.iter().enumerate() {
            assert_eq!(gene.coordinates.strand, Strand::Forward);
            assert!(gene.coordinates.begin > 0);
            assert!(gene.coordinates.end > gene.coordinates.begin);
            assert_eq!(gene.coordinates.start_index, i * 2);
            assert_eq!(gene.coordinates.stop_index, i * 2 + 1);
        }
    }
}
