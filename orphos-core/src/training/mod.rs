//! Training algorithms for gene prediction models.
//!
//! This module implements the unsupervised machine learning algorithms that
//! train Orphos's statistical models from genome sequences.
//!
//! ## Overview
//!
//! Training extracts statistical patterns from genes predicted in an initial pass:
//!
//! 1. **Initial gene finding**: Find high-confidence genes using basic models
//! 2. **Codon usage**: Calculate dicodon frequencies in predicted genes
//! 3. **Start codon preference**: Learn ATG/GTG/TTG usage patterns
//! 4. **RBS detection**: Identify ribosome binding site motifs (Shine-Dalgarno)
//! 5. **Upstream composition**: Analyze nucleotide patterns near start codons
//! 6. **GC bias**: Detect reading frame preferences based on GC content
//!
//! ## Training Modes
//!
//! - **Shine-Dalgarno (SD)**: For organisms with canonical RBS motifs
//! - **Non-SD**: For organisms without RBS or with alternative start recognition
//!
//! The mode is auto-detected based on the strength of SD signals in the training data.
//!
//! ## Modules
//!
//! - [`sd_training`]: Shine-Dalgarno motif training
//! - [`non_sd_training`]: Alternative start recognition training
//! - [`common`]: Shared training utilities
//!
//! ## Examples
//!
//! Training is normally performed automatically by the `OrphosAnalyzer`, but
//! can be done manually for advanced use cases:
//!
//! ```rust,no_run
//! use orphos_core::engine::UntrainedOrphos;
//! use orphos_core::config::OrphosConfig;
//! use orphos_core::sequence::encoded::EncodedSequence;
//!
//! let mut orphos = UntrainedOrphos::new();
//! let sequence = b"ATGAAACGCATTAGCACCACCATT...";
//! let encoded = EncodedSequence::without_masking(sequence);
//!
//! // Train on the genome
//! let trained = orphos.train_single_genome(&encoded)?;
//!
//! // Training data is now stored in the TrainedOrphos instance
//! # Ok::<(), orphos_core::types::OrphosError>(())
//! ```

pub mod common;
pub mod non_sd_training;
pub mod sd_training;

use bio::bio_types::strand::Strand;
use rayon::prelude::*;

use crate::{
    constants::{
        GC_HIGH_AT_FREQ, GC_HIGH_GC_FREQ, GC_LOW_AT_FREQ, GC_LOW_GC_FREQ, GENE_RATIO_THRESHOLD,
        HIGH_GC_FREQ, INITIAL_SCORE_THRESHOLD, LOW_GC_FREQ, MAX_GC_CONTENT, MAX_MOTIF_INDEX,
        MAX_TRAINING_ITERATIONS_NONSD, MAX_TRAINING_ITERATIONS_SD, MIN_GC_CONTENT,
        MIN_MOTIF_LENGTH, NUM_BASES, NUM_CODON_TYPES, NUM_MOTIF_SIZES, NUM_RBS_WEIGHTS,
        RBS_WEIGHT_STRONG_THRESHOLD, RBS_WEIGHT_THRESHOLD_HIGH, RBS_WEIGHT_THRESHOLD_LOW,
        THRESHOLD_DIVISOR, UPSTREAM_END_POS, UPSTREAM_MOTIF_COVERAGE_THRESHOLD, UPSTREAM_POSITIONS,
        UPSTREAM_SKIP_END, UPSTREAM_START_POS, WEIGHT_CLAMP_MAX, WEIGHT_CLAMP_MIN,
    },
    sequence::calculate_kmer_index,
    types::{CodonType, MotifWeights, Node, OrphosError, Training},
};

/// Checks if training should use Shine-Dalgarno motifs.
///
/// Analyzes the RBS weights learned during training to determine whether
/// the organism uses canonical Shine-Dalgarno ribosome binding sites.
///
/// # Arguments
///
/// * `training` - Training data with RBS weights
///
/// # Returns
///
/// `true` if SD motifs should be used, `false` for non-SD start recognition.
///
/// # Algorithm
///
/// The decision is based on:
/// - Strength of canonical SD motifs (AGGAGG variants)
/// - Absence of strong SD signals indicates non-SD organism
/// - Threshold-based classification using empirically determined cutoffs
#[must_use]
pub fn should_use_sd(training: &Training) -> bool {
    if training.rbs_weights[0] >= 0.0 {
        return false;
    }
    if training.rbs_weights[16] < RBS_WEIGHT_THRESHOLD_HIGH
        && training.rbs_weights[13] < RBS_WEIGHT_THRESHOLD_HIGH
        && training.rbs_weights[15] < RBS_WEIGHT_THRESHOLD_HIGH
        && (training.rbs_weights[0] >= RBS_WEIGHT_THRESHOLD_LOW
            || (training.rbs_weights[22] < RBS_WEIGHT_STRONG_THRESHOLD
                && training.rbs_weights[24] < RBS_WEIGHT_STRONG_THRESHOLD
                && training.rbs_weights[27] < RBS_WEIGHT_STRONG_THRESHOLD))
    {
        return false;
    }
    true
}

fn count_upstream_composition(
    sequence: &[u8],
    sequence_length: usize,
    strand: Strand,
    position: usize,
    training: &mut Training,
) {
    let start_position = match strand {
        Strand::Forward => position,
        Strand::Reverse => sequence_length - 1 - position,
        Strand::Unknown => unreachable!(),
    };

    let mut count = 0;
    for upstream_index in UPSTREAM_START_POS..UPSTREAM_END_POS {
        if upstream_index > 2 && upstream_index < UPSTREAM_SKIP_END {
            continue;
        }

        if start_position >= upstream_index {
            let base_index = calculate_kmer_index(1, sequence, start_position - upstream_index);
            training.upstream_composition[count][base_index] += 1.0;
        }
        count += 1;
    }
}

/// Convert upstream composition counts to log likelihood scores
fn convert_upstream_composition_to_log_scores(training: &mut Training) {
    training
        .upstream_composition
        .par_iter_mut()
        .take(UPSTREAM_POSITIONS)
        .for_each(|position_data| {
            let sum: f64 = position_data.iter().sum();

            if sum > 0.0 {
                for (j, value) in position_data.iter_mut().enumerate().take(NUM_BASES) {
                    *value /= sum;

                    let bg_freq = if training.gc_content > MIN_GC_CONTENT
                        && training.gc_content < MAX_GC_CONTENT
                    {
                        if j == 0 || j == 3 {
                            (1.0 - training.gc_content) / 2.0
                        } else {
                            training.gc_content / 2.0
                        }
                    } else if training.gc_content <= MIN_GC_CONTENT {
                        if j == 0 || j == 3 {
                            LOW_GC_FREQ
                        } else {
                            HIGH_GC_FREQ
                        }
                    } else if j == 0 || j == 3 {
                        HIGH_GC_FREQ
                    } else {
                        LOW_GC_FREQ
                    };

                    *value = (*value / bg_freq)
                        .ln()
                        .clamp(WEIGHT_CLAMP_MIN, WEIGHT_CLAMP_MAX);
                }
            } else {
                *position_data = [0.0; NUM_BASES];
            }
        });
}

/// Update motif counts during training
fn update_motif_counts(
    motif_counts: &mut MotifWeights,
    zero_motif_count: &mut f64,
    sequence: &[u8],
    reverse_sequence: &[u8],
    sequence_length: usize,
    node: &Node,
    stage: usize,
) {
    if node.position.codon_type == CodonType::Stop || node.position.is_edge {
        return;
    }

    let is_zero = node.motif_info.best_motif.length == 0;
    if is_zero {
        *zero_motif_count += 1.0;
        return;
    }

    let (working_sequence, start_position) = match node.position.strand {
        Strand::Forward => (sequence, node.position.index),
        Strand::Reverse => (reverse_sequence, sequence_length - 1 - node.position.index),
        Strand::Unknown => unreachable!(),
    };

    match stage {
        0 => {
            for (motif_length_index, motif_counts_length) in
                motif_counts.iter_mut().enumerate().take(NUM_MOTIF_SIZES)
            {
                let motif_length = motif_length_index + MIN_MOTIF_LENGTH;
                let j_start = start_position as isize - 18 - motif_length_index as isize;
                let j_end = start_position as isize - 6 - motif_length_index as isize;
                for j in j_start..=j_end {
                    if j < 0 {
                        continue;
                    }
                    let motif_index =
                        calculate_kmer_index(motif_length, working_sequence, j as usize);
                    for spacer_data in motif_counts_length.iter_mut().take(NUM_MOTIF_SIZES) {
                        spacer_data[motif_index] += 1.0;
                    }
                }
            }
        }
        1 => {
            let motif = &node.motif_info.best_motif;
            motif_counts[motif.length - MIN_MOTIF_LENGTH][motif.space_index][motif.index] += 1.0;

            for (submotif_length_index, motif_counts_length) in motif_counts
                .iter_mut()
                .enumerate()
                .take(motif.length - MIN_MOTIF_LENGTH)
            {
                let submotif_length = submotif_length_index + MIN_MOTIF_LENGTH;
                let j_start = start_position as isize - (motif.spacer + motif.length) as isize;
                let j_end = start_position as isize - (motif.spacer + submotif_length) as isize;
                for j in j_start..=j_end {
                    if j < 0 {
                        continue;
                    }
                    let spacer_index =
                        get_spacer_index(j as usize, start_position, submotif_length_index);
                    let motif_index =
                        calculate_kmer_index(submotif_length, working_sequence, j as usize);
                    motif_counts_length[spacer_index][motif_index] += 1.0;
                }
            }
        }
        2 => {
            let motif = &node.motif_info.best_motif;
            motif_counts[motif.length - MIN_MOTIF_LENGTH][motif.space_index][motif.index] += 1.0;
        }
        _ => {}
    }
}

#[allow(clippy::needless_range_loop)]
fn build_coverage_map(
    real_motifs: &[[[f64; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES],
    good_motifs: &mut [[[i32; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES],
    number_of_genes: f64,
    _stage: usize,
) {
    let threshold = UPSTREAM_MOTIF_COVERAGE_THRESHOLD;

    // Initialize all as not good
    *good_motifs = [[[0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES];

    // 3-base motifs: mark as good if above threshold
    for spacer_index in 0..NUM_MOTIF_SIZES {
        for motif_index in 0..64 {
            if real_motifs[0][spacer_index][motif_index] / number_of_genes >= threshold {
                for alternative_spacer_index in 0..NUM_MOTIF_SIZES {
                    good_motifs[0][alternative_spacer_index][motif_index] = 1;
                }
            }
        }
    }

    // 4-base motifs: must contain two valid 3-base motifs
    for spacer_index in 0..NUM_MOTIF_SIZES {
        for motif_index in 0..256 {
            let decomposition_0 = (motif_index & 252) >> 2;
            let decomposition_1 = motif_index & 63;
            if good_motifs[0][spacer_index][decomposition_0] != 0
                && good_motifs[0][spacer_index][decomposition_1] != 0
            {
                good_motifs[1][spacer_index][motif_index] = 1;
            }
        }
    }

    // 5-base motifs: interior mismatch allowed if entire 5-base motif
    // represents 3 valid 3-base motifs (if mismatch converted)
    for spacer_index in 0..NUM_MOTIF_SIZES {
        for motif_index in 0..1024 {
            let decomp0 = (motif_index & 1008) >> 4; // top 3 bases
            let decomp1 = (motif_index & 252) >> 2; // middle 3 bases
            let decomp2 = motif_index & 63; // bottom 3 bases
            if good_motifs[0][spacer_index][decomp0] == 0
                || good_motifs[0][spacer_index][decomp1] == 0
                || good_motifs[0][spacer_index][decomp2] == 0
            {
                continue;
            }
            good_motifs[2][spacer_index][motif_index] = 1;

            let mut tmp = motif_index;
            for k in (0..=16).step_by(16) {
                tmp ^= k;
                for l in (0..=32).step_by(32) {
                    tmp ^= l;
                    if good_motifs[2][spacer_index][tmp] == 0 {
                        good_motifs[2][spacer_index][tmp] = 2; // good with mismatch
                    }
                }
            }
        }
    }

    // 6-base motifs: must contain two valid 5-base motifs
    for spacer_index in 0..NUM_MOTIF_SIZES {
        for motif_index in 0..4096 {
            let decomp0 = (motif_index & 4092) >> 2; // top 5 bases
            let decomp1 = motif_index & 1023; // bottom 5 bases
            if good_motifs[2][spacer_index][decomp0] == 0
                || good_motifs[2][spacer_index][decomp1] == 0
            {
                continue;
            }
            if good_motifs[2][spacer_index][decomp0] == 1
                && good_motifs[2][spacer_index][decomp1] == 1
            {
                good_motifs[3][spacer_index][motif_index] = 1;
            } else {
                good_motifs[3][spacer_index][motif_index] = 2; // good with mismatch
            }
        }
    }
}

/// Update motif weights based on real vs background counts
fn update_motif_weights(
    motif_real: &[[[f64; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES],
    motif_background: &[[[f64; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES],
    zero_motif_real: f64,
    zero_motif_background: f64,
    motif_good: &[[[i32; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES],
    stage: usize,
    training: &mut Training,
) {
    // Stage is unused in the weight update after aligning with the C logic
    // (zbg accumulation no longer depends on stage).
    let _ = stage;
    // 1) sum = sum(mreal) + zreal
    // 2) For bad motifs: zreal += mreal; zbg += mreal; set mreal=0 and mbg=0
    // 3) Normalize mreal by sum and compute mot_wt = log(mreal/bg) or -4 if bg==0
    // 4) zreal /= sum; no_mot = log(zreal/zbg) or -4 if zbg==0; clamp in [-4,4]

    let mut sum_real = zero_motif_real;
    for motif in motif_real.iter() {
        for motif_row in motif.iter() {
            for &value in motif_row.iter().take(MAX_MOTIF_INDEX) {
                sum_real += value;
            }
        }
    }

    if sum_real == 0.0 {
        // If no real counts, zero all motif weights and no_mot as in C
        training
            .motif_weights
            .fill([[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]);
        training.no_motif_weight = 0.0;
        return;
    }

    // Make a mutable copy of background so we can zero out bad motifs (matches C)
    let mut bg = *motif_background;

    // counts of bad motifs to zbg), then zero their slots in the background matrix.
    let mut zreal = zero_motif_real;
    let mut zbg = zero_motif_background; // already normalized earlier

    let mut _good_real_sum = 0.0f64;
    let mut _bad_real_sum = 0.0f64;

    // We'll accumulate a temporary array for normalized real frequencies of good motifs
    let mut real_freqs = [[[0.0f64; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES];

    for l in 0..NUM_MOTIF_SIZES {
        for s in 0..NUM_MOTIF_SIZES {
            for idx in 0..MAX_MOTIF_INDEX {
                let r = motif_real[l][s][idx];
                if motif_good[l][s][idx] == 0 {
                    // Bad motif: fold raw real counts into zreal and zbg
                    // (C code: zreal += mreal; zbg += mreal)
                    zreal += r;
                    zbg += r;
                    _bad_real_sum += r;
                    // Zero out background slot (so bg==0 triggers -4.0)
                    bg[l][s][idx] = 0.0;
                    // real_freqs remains 0 here
                } else {
                    // Good motif: keep and normalize by sum later
                    _good_real_sum += r;
                    real_freqs[l][s][idx] = r / sum_real;
                }
            }
        }
    }

    for l in 0..NUM_MOTIF_SIZES {
        for s in 0..NUM_MOTIF_SIZES {
            for idx in 0..MAX_MOTIF_INDEX {
                let rf = real_freqs[l][s][idx];
                let b = bg[l][s][idx];
                let w = if b != 0.0 {
                    (rf / b).ln()
                } else {
                    WEIGHT_CLAMP_MIN
                };
                training.motif_weights[l][s][idx] = w.clamp(WEIGHT_CLAMP_MIN, WEIGHT_CLAMP_MAX);
            }
        }
    }

    let zreal_freq = zreal / sum_real;
    let mut no_mot = if zbg != 0.0 {
        (zreal_freq / zbg).ln()
    } else {
        WEIGHT_CLAMP_MIN
    };
    no_mot = no_mot.clamp(WEIGHT_CLAMP_MIN, WEIGHT_CLAMP_MAX);
    training.no_motif_weight = no_mot;

    // removed TAP motif_weights_debug
}

/// Get spacer index based on position
const fn get_spacer_index(
    position_index: usize,
    start_position: usize,
    motif_length_index: usize,
) -> usize {
    if position_index + 16 + motif_length_index <= start_position {
        3 // 13-15bp spacer
    } else if position_index + 14 + motif_length_index <= start_position {
        2 // 11-12bp spacer
    } else if position_index + 7 + motif_length_index >= start_position {
        1 // 3-4bp spacer
    } else {
        0 // 5-10bp spacer
    }
}

pub fn load_training_file(_file: &str) -> Training {
    Training::default()
}

pub const fn write_training_file(_file: &str, _training: &Training) -> Result<(), OrphosError> {
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Motif, NodeMotifInfo, NodePosition, NodeScores, NodeState};
    use bio::bio_types::strand::Strand;

    fn create_test_training() -> Training {
        Training {
            gc_content: 0.5,
            translation_table: 11,
            uses_shine_dalgarno: true,
            start_type_weights: [0.0; 3],
            rbs_weights: Box::new([0.0; 28]),
            upstream_composition: Box::new([[0.0; 4]; 32]),
            motif_weights: Box::new([[[0.0; 4096]; 4]; 4]),
            no_motif_weight: 0.0,
            start_weight_factor: 4.35,
            gc_bias_factors: [1.0; 3],
            gene_dicodon_table: Box::new([0.0; 4096]),
            total_dicodons: 0,
        }
    }

    fn create_test_node_with_motif(
        index: usize,
        strand: Strand,
        codon_type: CodonType,
        motif_length: usize,
        motif_index: usize,
    ) -> Node {
        Node {
            position: NodePosition {
                index,
                strand,
                codon_type,
                stop_value: (index + 100) as isize,
                is_edge: false,
            },
            scores: NodeScores {
                gc_content: 0.5,
                coding_score: 5.0,
                start_score: 2.0,
                ribosome_binding_score: 1.0,
                type_score: 1.5,
                upstream_score: 0.5,
                total_score: 10.0,
                gc_frame_scores: [1.0, 2.0, 3.0],
            },
            state: NodeState::default(),
            motif_info: NodeMotifInfo {
                ribosome_binding_sites: [0, 0],
                best_motif: Motif {
                    index: motif_index,
                    length: motif_length,
                    space_index: 0,
                    spacer: 8,
                    score: 2.0,
                },
            },
        }
    }

    #[test]
    fn test_should_use_sd_positive_rbs_weight() {
        let mut training = create_test_training();
        training.rbs_weights[0] = 1.0;

        assert!(!should_use_sd(&training));
    }

    #[test]
    fn test_should_use_sd_low_key_motifs() {
        let mut training = create_test_training();
        training.rbs_weights[0] = -1.0;
        training.rbs_weights[16] = 0.5;
        training.rbs_weights[13] = 0.5;
        training.rbs_weights[15] = 0.5;

        assert!(!should_use_sd(&training));
    }

    #[test]
    fn test_should_use_sd_high_key_motifs() {
        let mut training = create_test_training();
        training.rbs_weights[0] = -1.0;
        training.rbs_weights[16] = 2.0;
        training.rbs_weights[13] = 2.0;
        training.rbs_weights[15] = 2.0;

        assert!(should_use_sd(&training));
    }

    #[test]
    fn test_count_upstream_composition_forward() {
        let sequence = b"ATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCG";
        let sequence_length = sequence.len();
        let mut training = create_test_training();
        let position = 20;

        count_upstream_composition(
            sequence,
            sequence_length,
            Strand::Forward,
            position,
            &mut training,
        );

        // Check that some composition was counted
        let total_counts: f64 = training
            .upstream_composition
            .iter()
            .flat_map(|row| row.iter())
            .sum();

        assert!(total_counts > 0.0);
    }

    #[test]
    fn test_count_upstream_composition_reverse() {
        let sequence = b"ATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCG";
        let sequence_length = sequence.len();
        let mut training = create_test_training();
        let position = 20;

        count_upstream_composition(
            sequence,
            sequence_length,
            Strand::Reverse,
            position,
            &mut training,
        );

        // Check that some composition was counted
        let total_counts: f64 = training
            .upstream_composition
            .iter()
            .flat_map(|row| row.iter())
            .sum();

        assert!(total_counts > 0.0);
    }

    #[test]
    fn test_count_upstream_composition_edge_position() {
        let sequence = b"ATGCGATGC";
        let sequence_length = sequence.len();
        let mut training = create_test_training();
        let position = 2; // Near start of sequence

        count_upstream_composition(
            sequence,
            sequence_length,
            Strand::Forward,
            position,
            &mut training,
        );
    }

    #[test]
    fn test_convert_upstream_composition_to_log_scores() {
        let mut training = create_test_training();

        training.upstream_composition[0] = [10.0, 5.0, 3.0, 2.0];
        training.upstream_composition[1] = [1.0, 1.0, 1.0, 1.0];

        convert_upstream_composition_to_log_scores(&mut training);

        assert_ne!(training.upstream_composition[0], [10.0, 5.0, 3.0, 2.0]);

        for position in &*training.upstream_composition {
            for &value in position {
                assert!(value >= WEIGHT_CLAMP_MIN);
                assert!(value <= WEIGHT_CLAMP_MAX);
            }
        }
    }

    #[test]
    fn test_convert_upstream_composition_zero_sum() {
        let mut training = create_test_training();

        training.upstream_composition[0] = [0.0, 0.0, 0.0, 0.0];

        convert_upstream_composition_to_log_scores(&mut training);

        assert_eq!(training.upstream_composition[0], [0.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_update_motif_counts_stage_0() {
        let sequence = b"ATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCG";
        let reverse_sequence = b"CGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATG";
        let sequence_length = sequence.len();
        let mut motif_counts =
            Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let mut zero_motif_count = 0.0;

        let node = create_test_node_with_motif(20, Strand::Forward, CodonType::Atg, 4, 10);

        update_motif_counts(
            &mut motif_counts,
            &mut zero_motif_count,
            sequence,
            reverse_sequence,
            sequence_length,
            &node,
            0,
        );

        let total_counts: f64 = motif_counts
            .iter()
            .flat_map(|l| l.iter())
            .flat_map(|s| s.iter())
            .sum();

        assert!(total_counts > 0.0);
    }

    #[test]
    fn test_update_motif_counts_zero_motif() {
        let sequence = b"ATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCG";
        let reverse_sequence = b"CGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATG";
        let sequence_length = sequence.len();
        let mut motif_counts =
            Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let mut zero_motif_count = 0.0;

        let mut node = create_test_node_with_motif(20, Strand::Forward, CodonType::Atg, 0, 0);
        node.motif_info.best_motif.length = 0;

        update_motif_counts(
            &mut motif_counts,
            &mut zero_motif_count,
            sequence,
            reverse_sequence,
            sequence_length,
            &node,
            0,
        );

        // Should increment zero motif count
        assert_eq!(zero_motif_count, 1.0);
    }

    #[test]
    fn test_update_motif_counts_stop_codon() {
        let sequence = b"ATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCG";
        let reverse_sequence = b"CGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATGCGATG";
        let sequence_length = sequence.len();
        let mut motif_counts =
            Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let mut zero_motif_count = 0.0;

        let node = create_test_node_with_motif(20, Strand::Forward, CodonType::Stop, 4, 10);

        update_motif_counts(
            &mut motif_counts,
            &mut zero_motif_count,
            sequence,
            reverse_sequence,
            sequence_length,
            &node,
            0,
        );

        // Should not count anything for stop codons
        let total_counts: f64 = motif_counts
            .iter()
            .flat_map(|l| l.iter())
            .flat_map(|s| s.iter())
            .sum();

        assert_eq!(total_counts, 0.0);
        assert_eq!(zero_motif_count, 0.0);
    }

    #[test]
    fn test_build_coverage_map() {
        let real_motifs = Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let mut good_motifs = Box::new([[[0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let number_of_genes = 100.0;

        build_coverage_map(&real_motifs, &mut good_motifs, number_of_genes, 0);

        // Should initialize all motifs as not good (0)
        let has_good_motifs = good_motifs
            .iter()
            .flat_map(|l| l.iter())
            .flat_map(|s| s.iter())
            .any(|&val| val != 0);

        assert!(!has_good_motifs);
    }

    #[test]
    fn test_get_spacer_index() {
        let start_position = 100;
        let motif_length_index = 1;

        // Test different spacer distances
        assert_eq!(get_spacer_index(83, start_position, motif_length_index), 3); // 13-15bp
        assert_eq!(get_spacer_index(85, start_position, motif_length_index), 2); // 11-12bp
        assert_eq!(get_spacer_index(95, start_position, motif_length_index), 1); // 3-4bp
        assert_eq!(get_spacer_index(90, start_position, motif_length_index), 0); // 5-10bp
    }

    #[test]
    fn test_load_training_file() {
        let training = load_training_file("dummy_file.txt");

        // Should return default training
        assert_eq!(training.translation_table, 11);
        assert_eq!(training.gc_content, 0.5);
    }

    #[test]
    fn test_write_training_file() {
        let training = create_test_training();

        // Should not panic (placeholder implementation)
        let result = write_training_file("dummy_output.txt", &training);
        assert!(result.is_ok());
    }

    #[test]
    fn test_update_motif_weights_basic() {
        let motif_real = Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let motif_background =
            Box::new([[[0.0; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let zero_motif_real = 1.0;
        let zero_motif_background = 0.5;
        let motif_good = Box::new([[[1; MAX_MOTIF_INDEX]; NUM_MOTIF_SIZES]; NUM_MOTIF_SIZES]);
        let mut training = create_test_training();

        update_motif_weights(
            &motif_real,
            &motif_background,
            zero_motif_real,
            zero_motif_background,
            &motif_good,
            0,
            &mut training,
        );

        // Should complete without panicking
        // Function completion indicates success
    }
}
