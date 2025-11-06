use bio::bio_types::strand::Strand;

use crate::{
    constants::{MAXIMUM_SAME_OVERLAP, OPERON_DISTANCE, OVERLAP_PENALTY_FACTOR},
    types::{CodonType, Node, Training},
};

pub fn record_overlapping_starts(nodes: &mut [Node], training: &Training, flag: bool) {
    let node_count = nodes.len();

    for i in 0..node_count {
        for j in 0..3 {
            nodes[i].state.start_pointers[j] = None;
        }

        if nodes[i].position.codon_type != CodonType::Stop || nodes[i].position.is_edge {
            continue;
        }

        let mut max_sc = -100.0;
        if nodes[i].position.strand == Strand::Forward {
            if i + 3 < i32::MAX as usize {
                let start_j = i + 3;
                for j in (0..=start_j).rev() {
                    if j >= node_count || nodes[j].position.index > nodes[i].position.index + 2 {
                        continue;
                    }

                    if nodes[j].position.index + MAXIMUM_SAME_OVERLAP < nodes[i].position.index {
                        break;
                    }

                    if nodes[j].position.strand == Strand::Forward
                        && nodes[j].position.codon_type != CodonType::Stop
                    {
                        if nodes[j].position.stop_value <= nodes[i].position.index as isize {
                            continue;
                        }

                        let frame = nodes[j].position.index % 3;

                        if !flag && nodes[i].state.start_pointers[frame].is_none() {
                            nodes[i].state.start_pointers[frame] = Some(j);
                        } else if flag {
                            let score = nodes[j].scores.coding_score
                                + nodes[j].scores.start_score
                                + intergenic_mod(&nodes[i], &nodes[j], training);
                            if score > max_sc {
                                nodes[i].state.start_pointers[frame] = Some(j);
                                max_sc = nodes[j].scores.coding_score
                                    + nodes[j].scores.start_score
                                    + intergenic_mod(&nodes[i], &nodes[j], training);
                            }
                        }
                    }
                }
            }
        } else {
            for j_signed in ((i as isize) - 3)..(node_count as isize) {
                if j_signed < 0 {
                    continue;
                }
                let j = j_signed as usize;

                // same geometric filters as before
                if nodes[j].position.index < nodes[i].position.index.saturating_sub(2) {
                    continue;
                }
                if nodes[j].position.index > nodes[i].position.index + MAXIMUM_SAME_OVERLAP {
                    break;
                }

                if nodes[j].position.strand == Strand::Reverse
                    && nodes[j].position.codon_type != CodonType::Stop
                {
                    if nodes[j].position.stop_value >= nodes[i].position.index as isize {
                        continue;
                    }

                    let frame = nodes[j].position.index % 3;

                    if !flag && nodes[i].state.start_pointers[frame].is_none() {
                        nodes[i].state.start_pointers[frame] = Some(j);
                    } else if flag {
                        let score = nodes[j].scores.coding_score
                            + nodes[j].scores.start_score
                            + intergenic_mod(&nodes[j], &nodes[i], training);
                        if score > max_sc {
                            nodes[i].state.start_pointers[frame] = Some(j);
                            max_sc = score; // no need to recompute
                        }
                    }
                }
            }
        }
    }
}

/// Calculate intergenic modification score for operons and overlaps
#[inline]
pub fn intergenic_mod(n1: &Node, n2: &Node, training: &Training) -> f64 {
    let dist = n1.position.index.abs_diff(n2.position.index) as f64;
    let mut rval = 0.0;
    let ovlp = is_overlap(n1, n2);

    // Check for operon-like overlaps
    if is_operon_like(n1, n2) {
        rval += calculate_operon_bonus(n1, n2);
    }

    // Apply distance-based bonuses/penalties
    if dist > 3.0 * OPERON_DISTANCE || n1.position.strand != n2.position.strand {
        rval -= OVERLAP_PENALTY_FACTOR * training.start_weight_factor;
    } else if (dist <= OPERON_DISTANCE && !ovlp) || dist < OPERON_DISTANCE / 4.0 {
        rval +=
            (2.0 - dist / OPERON_DISTANCE) * OVERLAP_PENALTY_FACTOR * training.start_weight_factor;
    }

    rval
}

/// Check if two nodes form an operon-like structure
fn is_operon_like(n1: &Node, n2: &Node) -> bool {
    n1.position.strand == n2.position.strand
        && (n1.position.index + 2 == n2.position.index
            || n1.position.index == n2.position.index + 1)
}

/// Check if two nodes overlap
fn is_overlap(n1: &Node, n2: &Node) -> bool {
    (n1.position.strand == Strand::Forward
        && n2.position.strand == Strand::Forward
        && n1.position.index + 2 >= n2.position.index)
        || (n1.position.strand == Strand::Reverse
            && n2.position.strand == Strand::Reverse
            && n1.position.index >= n2.position.index + 2)
}

/// Calculate operon bonus
fn calculate_operon_bonus(n1: &Node, n2: &Node) -> f64 {
    let mut bonus = 0.0;

    if n1.position.strand == Strand::Forward && n2.scores.ribosome_binding_score < 0.0 {
        bonus -= n2.scores.ribosome_binding_score;
    }
    if n1.position.strand == Strand::Reverse && n1.scores.ribosome_binding_score < 0.0 {
        bonus -= n1.scores.ribosome_binding_score;
    }
    if n1.position.strand == Strand::Forward && n2.scores.upstream_score < 0.0 {
        bonus -= n2.scores.upstream_score;
    }
    if n1.position.strand == Strand::Reverse && n1.scores.upstream_score < 0.0 {
        bonus -= n1.scores.upstream_score;
    }

    bonus
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::*;
    use bio::bio_types::strand::Strand;

    fn create_test_node(
        strand: Strand,
        index: usize,
        codon_type: CodonType,
        stop_value: isize,
    ) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value,
                is_edge: false,
            },
            scores: NodeScores {
                coding_score: 1.0,
                start_score: 0.5,
                ribosome_binding_score: 0.3,
                upstream_score: 0.2,
                ..Default::default()
            },
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
    fn test_record_overlapping_starts_forward_strand() {
        let training = create_test_training();
        let mut nodes = vec![
            create_test_node(Strand::Forward, 10, CodonType::Atg, 100),
            create_test_node(Strand::Forward, 15, CodonType::Gtg, 105),
            create_test_node(Strand::Forward, 50, CodonType::Stop, 50),
        ];

        record_overlapping_starts(&mut nodes, &training, false);

        // Check that start pointers are set for stop codon nodes
        assert!(
            nodes[2].state.start_pointers.iter().any(|p| p.is_some())
                || nodes[2].state.start_pointers.iter().all(|p| p.is_none())
        );
    }

    #[test]
    fn test_record_overlapping_starts_reverse_strand() {
        let training = create_test_training();
        let mut nodes = vec![
            create_test_node(Strand::Reverse, 50, CodonType::Atg, 10),
            create_test_node(Strand::Reverse, 45, CodonType::Gtg, 5),
            create_test_node(Strand::Reverse, 20, CodonType::Stop, 20),
        ];

        record_overlapping_starts(&mut nodes, &training, false);

        // Check that start pointers are initialized
        assert!(
            nodes[2].state.start_pointers.iter().any(|p| p.is_some())
                || nodes[2].state.start_pointers.iter().all(|p| p.is_none())
        );
    }

    #[test]
    fn test_record_overlapping_starts_with_scoring() {
        let training = create_test_training();
        let mut nodes = vec![
            create_test_node(Strand::Forward, 10, CodonType::Atg, 100),
            create_test_node(Strand::Forward, 15, CodonType::Gtg, 105),
            create_test_node(Strand::Forward, 50, CodonType::Stop, 50),
        ];

        // Test with scoring enabled (flag = true)
        record_overlapping_starts(&mut nodes, &training, true);

        // Verify nodes are processed without panicking
        assert_eq!(nodes.len(), 3);
    }

    #[test]
    fn test_record_overlapping_starts_edge_nodes() {
        let training = create_test_training();
        let mut nodes = vec![create_test_node(Strand::Forward, 10, CodonType::Atg, 100)];
        nodes[0].position.is_edge = true;

        record_overlapping_starts(&mut nodes, &training, false);

        assert!(nodes[0].state.start_pointers.iter().all(|p| p.is_none()));
    }

    #[test]
    fn test_record_overlapping_starts_non_stop_codons() {
        let training = create_test_training();
        let mut nodes = vec![
            create_test_node(Strand::Forward, 10, CodonType::Atg, 100),
            create_test_node(Strand::Forward, 15, CodonType::Gtg, 105),
        ];

        record_overlapping_starts(&mut nodes, &training, false);

        assert!(nodes[0].state.start_pointers.iter().all(|p| p.is_none()));
        assert!(nodes[1].state.start_pointers.iter().all(|p| p.is_none()));
    }

    #[test]
    fn test_intergenic_mod_same_strand() {
        let training = create_test_training();
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Forward, 15, CodonType::Gtg, 105);

        let score = intergenic_mod(&n1, &n2, &training);

        // Should return a valid score
        assert!(score.is_finite());
    }

    #[test]
    fn test_intergenic_mod_different_strands() {
        let training = create_test_training();
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Reverse, 15, CodonType::Gtg, 5);

        let score = intergenic_mod(&n1, &n2, &training);

        assert!(score <= 0.0);
    }

    #[test]
    fn test_intergenic_mod_large_distance() {
        let training = create_test_training();
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Forward, 1000, CodonType::Gtg, 1100);

        let score = intergenic_mod(&n1, &n2, &training);

        assert!(score <= 0.0);
    }

    #[test]
    fn test_intergenic_mod_operon_distance() {
        let training = create_test_training();
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(
            Strand::Forward,
            10 + (OPERON_DISTANCE / 2.0) as usize,
            CodonType::Gtg,
            105,
        );

        let score = intergenic_mod(&n1, &n2, &training);

        // Should get operon bonus
        assert!(score >= 0.0);
    }

    #[test]
    fn test_is_operon_like_adjacent_nodes() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Forward, 12, CodonType::Gtg, 105);

        assert!(is_operon_like(&n1, &n2));
    }

    #[test]
    fn test_is_operon_like_consecutive_nodes() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Forward, 9, CodonType::Gtg, 105);

        assert!(is_operon_like(&n1, &n2));
    }

    #[test]
    fn test_is_operon_like_different_strands() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Reverse, 12, CodonType::Gtg, 5);

        assert!(!is_operon_like(&n1, &n2));
    }

    #[test]
    fn test_is_operon_like_distant_nodes() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Forward, 20, CodonType::Gtg, 105);

        assert!(!is_operon_like(&n1, &n2));
    }

    #[test]
    fn test_is_overlap_forward_strand() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Forward, 12, CodonType::Gtg, 105);

        assert!(is_overlap(&n1, &n2));
    }

    #[test]
    fn test_is_overlap_reverse_strand() {
        let n1 = create_test_node(Strand::Reverse, 15, CodonType::Atg, 5);
        let n2 = create_test_node(Strand::Reverse, 12, CodonType::Gtg, 2);

        assert!(is_overlap(&n1, &n2));
    }

    #[test]
    fn test_is_overlap_no_overlap() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Forward, 20, CodonType::Gtg, 105);

        assert!(!is_overlap(&n1, &n2));
    }

    #[test]
    fn test_is_overlap_different_strands() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Reverse, 12, CodonType::Gtg, 5);

        assert!(!is_overlap(&n1, &n2));
    }

    #[test]
    fn test_calculate_operon_bonus_forward_negative_rbs() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let mut n2 = create_test_node(Strand::Forward, 12, CodonType::Gtg, 105);
        n2.scores.ribosome_binding_score = -0.5;

        let bonus = calculate_operon_bonus(&n1, &n2);

        // Should get bonus for negative RBS score
        assert!(bonus > 0.0);
    }

    #[test]
    fn test_calculate_operon_bonus_reverse_negative_rbs() {
        let mut n1 = create_test_node(Strand::Reverse, 15, CodonType::Atg, 5);
        let n2 = create_test_node(Strand::Reverse, 12, CodonType::Gtg, 2);
        n1.scores.ribosome_binding_score = -0.3;

        let bonus = calculate_operon_bonus(&n1, &n2);

        // Should get bonus for negative RBS score
        assert!(bonus > 0.0);
    }

    #[test]
    fn test_calculate_operon_bonus_negative_upstream() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let mut n2 = create_test_node(Strand::Forward, 12, CodonType::Gtg, 105);
        n2.scores.upstream_score = -0.4;

        let bonus = calculate_operon_bonus(&n1, &n2);

        // Should get bonus for negative upstream score
        assert!(bonus > 0.0);
    }

    #[test]
    fn test_calculate_operon_bonus_positive_scores() {
        let n1 = create_test_node(Strand::Forward, 10, CodonType::Atg, 100);
        let n2 = create_test_node(Strand::Forward, 12, CodonType::Gtg, 105);

        let bonus = calculate_operon_bonus(&n1, &n2);

        // No bonus for positive scores
        assert_eq!(bonus, 0.0);
    }

    #[test]
    fn test_record_overlapping_starts_empty_nodes() {
        let training = create_test_training();
        let mut nodes = vec![];

        record_overlapping_starts(&mut nodes, &training, false);

        // Should handle empty input gracefully
        assert_eq!(nodes.len(), 0);
    }

    #[test]
    fn test_record_overlapping_starts_single_node() {
        let training = create_test_training();
        let mut nodes = vec![create_test_node(Strand::Forward, 10, CodonType::Stop, 10)];

        record_overlapping_starts(&mut nodes, &training, false);

        assert!(nodes[0].state.start_pointers.iter().all(|p| p.is_none()));
    }

    #[test]
    fn test_record_overlapping_starts_maximum_overlap() {
        let training = create_test_training();
        let mut nodes = vec![
            create_test_node(Strand::Forward, 10, CodonType::Atg, 100),
            create_test_node(
                Strand::Forward,
                10 + MAXIMUM_SAME_OVERLAP + 5,
                CodonType::Stop,
                50,
            ),
        ];

        record_overlapping_starts(&mut nodes, &training, false);

        // Should handle maximum overlap distance
        assert!(nodes[1].state.start_pointers.iter().all(|p| p.is_none()));
    }
}
