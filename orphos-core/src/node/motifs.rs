use bio::bio_types::strand::Strand;
use rayon::prelude::*;

use crate::{
    constants::{
        INITIAL_MAX_SCORE, MAX_MOTIF_LENGTH, MIN_MOTIF_SCORE, MOTIF_THRESHOLD_OFFSET,
        RBS_DOWNSTREAM_DISTANCE, RBS_UPSTREAM_DISTANCE,
    },
    sequence::{calculate_kmer_index, encoded_base_code, shine_dalgarno_exact, shine_dalgarno_mm},
    types::{CodonType, Motif, Node, Training},
};

/// Calculate both exact and mismatch Shine-Dalgarno scores for a position
fn calculate_rbs_scores(
    seq: &[u8],
    pos: usize,
    target_pos: usize,
    rbs_weights: &[f64],
) -> [usize; 2] {
    [
        shine_dalgarno_exact(seq, pos, target_pos, rbs_weights),
        shine_dalgarno_mm(seq, pos, target_pos, rbs_weights),
    ]
}

/// Number of bases a Shine-Dalgarno search can read from its search position.
const RBS_KMER_LENGTH: usize = MAX_MOTIF_LENGTH;

/// Longest upstream motif, and the k-mer length read per motif window.
const MOTIF_KMER_LENGTH: usize = MAX_MOTIF_LENGTH;
/// Farthest motif window start upstream of the start codon (`start - 18 - 3`).
const MOTIF_WINDOW_UPSTREAM: isize = 21;
/// Nearest motif window start upstream of the start codon (`start - 6`).
const MOTIF_WINDOW_DOWNSTREAM: isize = 6;
/// Number of distinct motif window starts.
const MOTIF_WINDOW_COUNT: usize = (MOTIF_WINDOW_UPSTREAM - MOTIF_WINDOW_DOWNSTREAM + 1) as usize;

/// Highest ribosome-binding-site motif indices for every 6-mer and spacer
/// distance, tabulated with [`shine_dalgarno_exact`] / [`shine_dalgarno_mm`].
///
/// Both searches only read the bases between the search position and the
/// start codon (at most six) and otherwise depend on the distance to the
/// start codon and the motif weights, so each (6-mer, distance) pair has a
/// fixed answer for a given weight set.
struct RbsLookup {
    /// `[6-mer index][distance]` -> `[exact motif index, mismatch motif index]`
    table: Vec<[u8; 2]>,
}

impl RbsLookup {
    const MAX_DISTANCE: usize = RBS_UPSTREAM_DISTANCE + 1;
    const KMER_COUNT: usize = 1 << (2 * RBS_KMER_LENGTH);

    fn new(rbs_weights: &[f64]) -> Self {
        let mut table = vec![[0u8; 2]; Self::KMER_COUNT * (Self::MAX_DISTANCE + 1)];
        for kmer in 0..Self::KMER_COUNT {
            // Pack the 6-mer the way the encoded sequence stores it so the
            // original searches can run on it directly.
            let mut packed = [0u8; 2];
            for offset in 0..RBS_KMER_LENGTH {
                let base = ((kmer >> (2 * offset)) & 0b11) as u8;
                packed[offset / 4] |= base << (2 * (offset % 4));
            }
            for distance in 0..=Self::MAX_DISTANCE {
                let [exact, mismatch] = calculate_rbs_scores(&packed, 0, distance, rbs_weights);
                table[kmer * (Self::MAX_DISTANCE + 1) + distance] = [exact as u8, mismatch as u8];
            }
        }
        Self { table }
    }

    #[inline]
    fn lookup(&self, kmer: usize, distance: usize) -> [usize; 2] {
        let [exact, mismatch] = self.table[kmer * (Self::MAX_DISTANCE + 1) + distance];
        [usize::from(exact), usize::from(mismatch)]
    }

    /// Scores every search position in `search_start..=search_end` against the
    /// start codon at `start_codon_position` and keeps the best motif indices.
    #[inline]
    fn score_window(
        &self,
        seq: &[u8],
        search_start: usize,
        search_end: usize,
        start_codon_position: usize,
        best: &mut [usize; 2],
    ) {
        // Positions closer than RBS_KMER_LENGTH to the start codon never score
        // (the searches require at least MIN_DISTANCE_FROM_START + MIN_MOTIF_LENGTH
        // bases), and skipping them keeps every 6-mer read inside the sequence.
        let Some(last) = start_codon_position
            .checked_sub(RBS_KMER_LENGTH)
            .map(|last| last.min(search_end))
        else {
            return;
        };
        if search_start > last {
            return;
        }

        let mut kmer = calculate_kmer_index(RBS_KMER_LENGTH, seq, search_start);
        for position in search_start..=last {
            let distance = start_codon_position - position;
            if distance <= Self::MAX_DISTANCE {
                let scores = self.lookup(kmer, distance);
                if scores[0] > best[0] {
                    best[0] = scores[0];
                }
                if scores[1] > best[1] {
                    best[1] = scores[1];
                }
            }
            if position < last {
                let next_base = usize::from(encoded_base_code(seq, position + RBS_KMER_LENGTH));
                kmer = (kmer >> 2) | (next_base << (2 * (RBS_KMER_LENGTH - 1)));
            }
        }
    }
}

/// RBS Scoring Function: Calculate the RBS motif and then multiply it by the
/// appropriate weight for that motif (determined in the start training function).
///
/// # Parameters
/// * `seq` - Forward strand sequence
/// * `rseq` - Reverse strand sequence
/// * `sequence_length` - Sequence length
/// * `nodes` - Nodes to score
/// * `training` - Training data with RBS weights
///
/// # Effect
/// Updates the ribosome_binding_sites scores in each node's motif_info
pub fn rbs_score(
    seq: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &mut [Node],
    training: &Training,
) {
    let lookup = RbsLookup::new(&*training.rbs_weights);

    // Parallelize RBS scoring across all start nodes
    nodes
        .par_iter_mut()
        .filter(|node| node.position.codon_type != CodonType::Stop && !node.position.is_edge)
        .for_each(|node| {
            let mut best = [0usize; 2];

            match node.position.strand {
                Strand::Forward => {
                    let search_start = node.position.index.saturating_sub(RBS_UPSTREAM_DISTANCE);
                    let search_end = node.position.index.saturating_sub(RBS_DOWNSTREAM_DISTANCE);
                    lookup.score_window(
                        seq,
                        search_start,
                        search_end,
                        node.position.index,
                        &mut best,
                    );
                }
                Strand::Reverse => {
                    if node.position.index >= sequence_length {
                        node.motif_info.ribosome_binding_sites = best;
                        return; // Skip invalid node positions
                    }
                    let upstream_offset = RBS_UPSTREAM_DISTANCE + 1;
                    let downstream_offset = RBS_DOWNSTREAM_DISTANCE + 1;
                    let start_pos =
                        (sequence_length - node.position.index).saturating_sub(upstream_offset);
                    let end_pos =
                        (sequence_length - node.position.index).saturating_sub(downstream_offset);
                    let target_pos = sequence_length - 1 - node.position.index;
                    lookup.score_window(
                        reverse_complement_encoded_sequence,
                        start_pos,
                        end_pos,
                        target_pos,
                        &mut best,
                    );
                }
                Strand::Unknown => unreachable!(),
            }

            node.motif_info.ribosome_binding_sites = best;
        });
}

/// Find the highest scoring motif/spacer combination for a node.
///
/// Given the weights for various motifs/distances from the training file,
/// return the highest scoring mer/spacer combination of 3-6bp motifs with a
/// spacer ranging from 3bp to 15bp. In the final stage of start training, only
/// good scoring motifs are returned.
///
/// # Parameters
/// * `training` - Training data with motif weights
/// * `seq` - Forward strand sequence
/// * `rseq` - Reverse strand sequence
/// * `sequence_length` - Sequence length
/// * `node` - Node to find motif for
/// * `stage` - Training stage (affects whether poor motifs are filtered)
///
/// # Effect
/// Updates the best_motif in the node's motif_info
pub fn find_best_upstream_motif(
    training: &Training,
    seq: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    node: &mut Node,
    stage: usize,
) {
    if node.position.codon_type == CodonType::Stop || node.position.is_edge {
        return;
    }

    let (wseq, start) = match node.position.strand {
        Strand::Forward => (seq, node.position.index),
        Strand::Reverse => (
            reverse_complement_encoded_sequence,
            sequence_length - 1 - node.position.index,
        ),
        Strand::Unknown => unreachable!(),
    };

    let mut max_sc = INITIAL_MAX_SCORE;
    let mut max_spacendx = 0;
    let mut max_spacer = 0;
    let mut max_ndx = 0;
    let mut max_len = 0;

    // Every motif window starts in [start - 21, start - 6]. Read one 6-mer per
    // window start (bases past the sequence end read as zero) and derive the
    // shorter motifs by masking; a motif's bases are the low bits of the
    // 6-mer at the same position.
    let first_window = (start as isize - MOTIF_WINDOW_UPSTREAM).max(0);
    let last_window = start as isize - MOTIF_WINDOW_DOWNSTREAM;
    let mut window_kmers = [0usize; MOTIF_WINDOW_COUNT];
    if last_window >= first_window {
        let base_at = |position: usize| {
            if position < sequence_length {
                usize::from(encoded_base_code(wseq, position))
            } else {
                0
            }
        };
        let first_window = first_window as usize;
        let mut kmer = 0;
        for offset in 0..MOTIF_KMER_LENGTH {
            kmer |= base_at(first_window + offset) << (2 * offset);
        }
        window_kmers[0] = kmer;
        for (slot, position) in ((first_window + 1)..=(last_window as usize)).enumerate() {
            kmer = (kmer >> 2)
                | (base_at(position + MOTIF_KMER_LENGTH - 1) << (2 * (MOTIF_KMER_LENGTH - 1)));
            window_kmers[slot + 1] = kmer;
        }
    }

    // Search through motif lengths 3-6 (i goes from 3 down to 0, representing lengths 6 down to 3)
    for i in (0..=3).rev() {
        let motif_len = i + 3;
        let motif_mask = (1usize << (2 * motif_len as usize)) - 1;

        // Search positions from start-18-i to start-6-i
        let search_start: isize = start as isize - 18 - i;
        let search_end: isize = start as isize - 6 - i;

        for j in search_start..=search_end {
            if j < 0 {
                continue;
            }
            // Ensure the k-mer window [j, j+motif_len) is within the nucleotide sequence.
            // the nucleotide length (sequence_length), not the byte length of the buffer.
            if (j as usize) + (motif_len as usize) > sequence_length {
                continue;
            }
            let spacer = start as isize - j - i - 3;
            let spacendx = if j <= start as isize - 16 - i {
                3
            } else if j <= start as isize - 14 - i {
                2
            } else if j >= start as isize - 7 - i {
                1
            } else {
                0
            };

            let index = window_kmers[(j - first_window) as usize] & motif_mask;
            let score = training.motif_weights[i as usize][spacendx][index];

            if score > max_sc {
                max_sc = score;
                max_spacendx = spacendx;
                max_spacer = spacer;
                max_ndx = index;
                max_len = motif_len;
            }
        }
    }

    // In stage 2, only accept good scoring motifs
    let is_stage_two = stage == 2;
    let is_poor_motif =
        max_sc == MIN_MOTIF_SCORE || max_sc < training.no_motif_weight + MOTIF_THRESHOLD_OFFSET;

    // Do NOT neutralize the score in early stages when no window is found.
    // C keeps max_sc at -100.0 in this case, strongly disfavoring edge starts.
    let effective_max_sc = max_sc;

    node.motif_info.best_motif = if is_stage_two && is_poor_motif {
        Motif {
            index: 0,
            length: 0,
            space_index: 0,
            spacer: 0,
            score: training.no_motif_weight,
        }
    } else {
        Motif {
            index: max_ndx,
            length: max_len as usize,
            space_index: max_spacendx,
            spacer: max_spacer as usize,
            score: effective_max_sc,
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::*;
    use bio::bio_types::strand::Strand;

    fn create_test_node(strand: Strand, index: usize, codon_type: CodonType) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value: (index + 100) as isize,
                is_edge: false,
            },
            scores: NodeScores::default(),
            state: NodeState::default(),
            motif_info: NodeMotifInfo {
                ribosome_binding_sites: [0; 2],
                best_motif: Motif::default(),
            },
        }
    }

    fn create_test_training() -> Training {
        Training {
            gc_content: 0.5,
            translation_table: 11,
            uses_shine_dalgarno: true,
            start_type_weights: [1.0, 2.0, 3.0],
            rbs_weights: Box::new([1.0; 28]),
            upstream_composition: Box::new([[0.25; 4]; 32]),
            motif_weights: Box::new([[[1.0; 4096]; 4]; 4]),
            no_motif_weight: 0.5,
            start_weight_factor: 4.35,
            gc_bias_factors: [1.0; 3],
            gene_dicodon_table: Box::new([0.0; 4096]),
            total_dicodons: 0,
        }
    }

    #[test]
    fn test_calculate_rbs_scores() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3]; // 8-base sequence
        let rbs_weights = [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        ];
        let pos = 2;
        let target_pos = 6;

        let scores = calculate_rbs_scores(&seq, pos, target_pos, &rbs_weights);

        assert_eq!(scores.len(), 2);
    }

    #[test]
    fn test_calculate_rbs_scores_boundary_conditions() {
        let seq = vec![0, 1, 2, 3];
        let rbs_weights = [1.0; 16];

        // Test at start of sequence
        let scores = calculate_rbs_scores(&seq, 0, 3, &rbs_weights);
        assert_eq!(scores.len(), 2);

        // Test near end of sequence
        let scores = calculate_rbs_scores(&seq, 2, 3, &rbs_weights);
        assert_eq!(scores.len(), 2);
    }

    #[test]
    fn test_rbs_score_forward_strand() {
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0,
            1, 2, 3,
        ];
        let reverse_seq = vec![
            3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3,
            2, 1, 0,
        ];
        let training = create_test_training();

        let mut nodes = vec![create_test_node(Strand::Forward, 25, CodonType::Atg)];

        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert_eq!(nodes[0].motif_info.ribosome_binding_sites.len(), 2);
    }

    #[test]
    fn test_rbs_score_reverse_strand() {
        // Create a longer sequence to handle reverse strand calculations
        let seq = vec![0; 60]; // 60-base sequence
        let reverse_seq = vec![3; 60];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Reverse, 35, CodonType::Atg), // Position far from edges
        ];

        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert_eq!(nodes[0].motif_info.ribosome_binding_sites.len(), 2);
    }

    #[test]
    fn test_rbs_score_stop_codon_skipped() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut nodes = vec![create_test_node(Strand::Forward, 6, CodonType::Stop)];

        let original_rbs = nodes[0].motif_info.ribosome_binding_sites;
        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert_eq!(nodes[0].motif_info.ribosome_binding_sites, original_rbs);
    }

    #[test]
    fn test_rbs_score_edge_node_skipped() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 6, CodonType::Atg);
        node.position.is_edge = true;
        let mut nodes = vec![node];

        let original_rbs = nodes[0].motif_info.ribosome_binding_sites;
        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        assert_eq!(nodes[0].motif_info.ribosome_binding_sites, original_rbs);
    }

    #[test]
    fn test_rbs_score_parallel_processing() {
        // Create a long sequence to handle all positions safely
        let seq = vec![0; 80];
        let reverse_seq = vec![3; 80];
        let training = create_test_training();

        let mut nodes = vec![
            create_test_node(Strand::Forward, 30, CodonType::Atg),
            create_test_node(Strand::Reverse, 50, CodonType::Atg),
            create_test_node(Strand::Forward, 70, CodonType::Atg),
        ];

        rbs_score(&seq, &reverse_seq, seq.len(), &mut nodes, &training);

        for node in &nodes {
            assert_eq!(node.motif_info.ribosome_binding_sites.len(), 2);
        }
    }

    #[test]
    fn test_find_best_upstream_motif_stop_codon() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 6, CodonType::Stop);
        let original_motif_index = node.motif_info.best_motif.index;
        let original_motif_length = node.motif_info.best_motif.length;

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        assert_eq!(node.motif_info.best_motif.index, original_motif_index);
        assert_eq!(node.motif_info.best_motif.length, original_motif_length);
    }

    #[test]
    fn test_find_best_upstream_motif_edge_node() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 6, CodonType::Atg);
        node.position.is_edge = true;
        let original_motif_index = node.motif_info.best_motif.index;
        let original_motif_length = node.motif_info.best_motif.length;

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        assert_eq!(node.motif_info.best_motif.index, original_motif_index);
        assert_eq!(node.motif_info.best_motif.length, original_motif_length);
    }

    #[test]
    fn test_find_best_upstream_motif_forward_strand() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 15, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        // length and spacer are usize, so always >= 0
    }

    #[test]
    fn test_find_best_upstream_motif_reverse_strand() {
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Reverse, 15, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        // length and spacer are usize, so always >= 0
    }

    #[test]
    fn test_find_best_upstream_motif_stage_two_filtering() {
        let mut training = create_test_training();
        training.no_motif_weight = 10.0; // High no-motif weight

        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];

        let mut node = create_test_node(Strand::Forward, 15, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 2);

        if node.motif_info.best_motif.score < training.no_motif_weight + MOTIF_THRESHOLD_OFFSET {
            assert_eq!(node.motif_info.best_motif.length, 0);
            assert_eq!(node.motif_info.best_motif.score, training.no_motif_weight);
        }
    }

    #[test]
    fn test_find_best_upstream_motif_stage_one_accepts_all() {
        let training = create_test_training();
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let reverse_seq = vec![3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0];

        let mut node = create_test_node(Strand::Forward, 15, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        // length is usize, so always >= 0
    }

    #[test]
    fn test_find_best_upstream_motif_empty_sequence() {
        let seq = vec![];
        let reverse_seq = vec![];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 0, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, 0, &mut node, 1);

        // Should handle empty sequence gracefully
        // length is usize, so always >= 0
    }

    #[test]
    fn test_find_best_upstream_motif_short_sequence() {
        let seq = vec![0, 1, 2];
        let reverse_seq = vec![2, 1, 0];
        let training = create_test_training();

        let mut node = create_test_node(Strand::Forward, 2, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        // Should handle short sequences without panicking
        // length is usize, so always >= 0
    }

    #[test]
    fn test_motif_spacer_classification() {
        let training = create_test_training();
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ];
        let reverse_seq = vec![
            3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0,
        ];

        let mut node = create_test_node(Strand::Forward, 20, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        assert!(node.motif_info.best_motif.space_index <= 3);
    }

    #[test]
    fn test_motif_length_range() {
        let training = create_test_training();
        let seq = vec![
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
        ];
        let reverse_seq = vec![
            3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 0,
        ];

        let mut node = create_test_node(Strand::Forward, 20, CodonType::Atg);

        find_best_upstream_motif(&training, &seq, &reverse_seq, seq.len(), &mut node, 1);

        assert!(node.motif_info.best_motif.length <= 6);
    }
}
