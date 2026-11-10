use std::sync::LazyLock;

use bio::bio_types::strand::Strand;
use rayon::prelude::*;

use crate::{
    constants::{
        INITIAL_MAX_SCORE, MAX_MOTIF_LENGTH, MIN_MOTIF_SCORE, MOTIF_THRESHOLD_OFFSET,
        RBS_DOWNSTREAM_DISTANCE, RBS_UPSTREAM_DISTANCE,
    },
    sequence::{
        calculate_kmer_index, encoded_base_code, is_a, is_g, shine_dalgarno_exact,
        shine_dalgarno_exact_from, shine_dalgarno_mm, shine_dalgarno_mm_from,
    },
    types::{CodonType, Motif, MotifWeights, Node, Training},
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

/// For every 6-mer, the bases that match the Shine-Dalgarno pattern the
/// searches look for (`A` at offsets 0 and 3, `G` at the others), one bit
/// per offset.
static RBS_MATCH_MASKS: LazyLock<[u8; RbsLookup::KMER_COUNT]> = LazyLock::new(|| {
    let mut masks = [0u8; RbsLookup::KMER_COUNT];
    for (kmer, mask) in masks.iter_mut().enumerate() {
        let packed = pack_kmer(kmer);
        for offset in 0..RBS_KMER_LENGTH {
            let matches = if offset % 3 == 0 {
                is_a(&packed, offset)
            } else {
                is_g(&packed, offset)
            };
            *mask |= u8::from(matches) << offset;
        }
    }
    masks
});

/// The 6-mer `kmer`, packed the way the encoded sequence stores it.
fn pack_kmer(kmer: usize) -> [u8; 2] {
    let mut packed = [0u8; 2];
    for offset in 0..RBS_KMER_LENGTH {
        let base = ((kmer >> (2 * offset)) & 0b11) as u8;
        packed[offset / 4] |= base << (2 * (offset % 4));
    }
    packed
}

/// Highest ribosome-binding-site motif indices for every 6-mer and spacer
/// distance, tabulated with [`shine_dalgarno_exact`] / [`shine_dalgarno_mm`].
///
/// Both searches only read the bases between the search position and the
/// start codon (at most six) and otherwise depend on the distance to the
/// start codon and the motif weights, so each (6-mer, distance) pair has a
/// fixed answer for a given weight set. Of the bases they only ask whether
/// each matches the pattern ([`RBS_MATCH_MASKS`]), so the 6-mers with the
/// same matches share one entry per distance.
#[derive(Debug)]
pub(crate) struct RbsLookup {
    /// `[match mask][distance]` -> `[exact motif index, mismatch motif index]`
    table: [[u8; 2]; Self::MASK_COUNT * (Self::MAX_DISTANCE + 1)],
    match_masks: &'static [u8; Self::KMER_COUNT],
    /// The motif weights, for the few searches that are not tabulated.
    rbs_weights: Vec<f64>,
}

impl RbsLookup {
    const MAX_DISTANCE: usize = RBS_UPSTREAM_DISTANCE + 1;
    const KMER_COUNT: usize = 1 << (2 * RBS_KMER_LENGTH);
    const MASK_COUNT: usize = 1 << RBS_KMER_LENGTH;

    /// Tabulates the searches for the motif weights `rbs_weights`.
    pub(crate) fn new(rbs_weights: &[f64]) -> Self {
        let match_masks: &'static [u8; Self::KMER_COUNT] = &RBS_MATCH_MASKS;
        let mut table = [[0u8; 2]; Self::MASK_COUNT * (Self::MAX_DISTANCE + 1)];
        let mut tabulated = [false; Self::MASK_COUNT];
        for (kmer, &mask) in match_masks.iter().enumerate() {
            // The first 6-mer with each set of matches stands for them all.
            let mask = usize::from(mask);
            if std::mem::replace(&mut tabulated[mask], true) {
                continue;
            }
            let packed = pack_kmer(kmer);
            for distance in 0..=Self::MAX_DISTANCE {
                let [exact, mismatch] = calculate_rbs_scores(&packed, 0, distance, rbs_weights);
                table[mask * (Self::MAX_DISTANCE + 1) + distance] = [exact as u8, mismatch as u8];
            }
        }
        Self {
            table,
            match_masks,
            rbs_weights: rbs_weights.to_vec(),
        }
    }

    #[inline]
    fn lookup(&self, kmer: usize, distance: usize) -> [usize; 2] {
        let mask = usize::from(self.match_masks[kmer & (Self::KMER_COUNT - 1)]);
        let [exact, mismatch] = self.table[mask * (Self::MAX_DISTANCE + 1) + distance];
        [usize::from(exact), usize::from(mismatch)]
    }

    /// Scores the search positions before the start of the sequence, down to
    /// [`RBS_UPSTREAM_DISTANCE`] bases upstream of the start codon at
    /// `start_codon_position`, and keeps the best motif indices.
    ///
    /// Prodigal skips those positions on the forward strand but not on the
    /// reverse strand, where a search from one of them can find a motif at
    /// the very start of the strand that no search from inside finds: the
    /// motif pattern is aligned to the search position. Only a start codon
    /// within twenty bases of the start of its strand has any.
    fn score_before_sequence(
        &self,
        seq: &[u8],
        start_codon_position: usize,
        best: &mut [usize; 2],
    ) {
        let first = start_codon_position as isize - RBS_UPSTREAM_DISTANCE as isize;
        for position in first..0 {
            let exact =
                shine_dalgarno_exact_from(seq, position, start_codon_position, &self.rbs_weights);
            let mismatch =
                shine_dalgarno_mm_from(seq, position, start_codon_position, &self.rbs_weights);
            best[0] = best[0].max(exact);
            best[1] = best[1].max(mismatch);
        }
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
    rbs_score_with(
        seq,
        reverse_complement_encoded_sequence,
        sequence_length,
        nodes,
        &RbsLookup::new(&*training.rbs_weights),
    );
}

/// [`rbs_score`] with the motif searches of the training's weights already
/// tabulated in `lookup`.
pub(crate) fn rbs_score_with(
    seq: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &mut [Node],
    lookup: &RbsLookup,
) {
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
                    lookup.score_before_sequence(
                        reverse_complement_encoded_sequence,
                        target_pos,
                        &mut best,
                    );
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

/// One of the motif windows [`find_best_upstream_motif`] scores for a start
/// codon: a motif of `length_index + 3` bases, `spacer` bases upstream of
/// the start codon.
#[derive(Clone, Copy)]
struct MotifCandidate {
    /// Motif length minus three.
    length_index: u8,
    /// Position of the window among the 6-mers of a [`MotifSite`].
    slot: u8,
    /// Spacer class: 0 for 5-10 bases, 1 for 3-4, 2 for 11-12, 3 for 13-15.
    space_index: u8,
    /// Bases between the motif and the start codon.
    spacer: u8,
}

/// Motif lengths: 3 to 6 bases.
const MOTIF_LENGTHS: usize = 4;
/// Spacers between a motif and the start codon: 3 to 15 bases.
const MOTIF_SPACERS: usize = 13;
/// Motif lengths times spacers.
const MOTIF_CANDIDATE_COUNT: usize = MOTIF_LENGTHS * MOTIF_SPACERS;
/// The `inside` bits of a [`MotifSite`] whose windows are all in the sequence.
const ALL_MOTIF_CANDIDATES: u64 = (1 << MOTIF_CANDIDATE_COUNT) - 1;

/// The candidate windows in the order the search visits them: motif lengths
/// from six bases down to three, and for each the windows from the farthest
/// upstream to the nearest. The first of several equally good windows wins.
const MOTIF_CANDIDATES: [MotifCandidate; MOTIF_CANDIDATE_COUNT] = {
    let mut candidates = [MotifCandidate {
        length_index: 0,
        slot: 0,
        space_index: 0,
        spacer: 0,
    }; MOTIF_CANDIDATE_COUNT];
    let mut count = 0;
    let mut length_index = 4;
    while length_index > 0 {
        length_index -= 1;
        // The window starts `distance + length_index` bases upstream of the
        // start codon, `distance` running from 18 down to 6.
        let mut distance = 18;
        while distance >= 6 {
            candidates[count] = MotifCandidate {
                length_index: length_index as u8,
                slot: (MOTIF_WINDOW_UPSTREAM as usize - distance - length_index) as u8,
                space_index: if distance >= 16 {
                    3
                } else if distance >= 14 {
                    2
                } else if distance <= 7 {
                    1
                } else {
                    0
                },
                spacer: (distance - 3) as u8,
            };
            count += 1;
            distance -= 1;
        }
    }
    candidates
};

/// What a motif search reads of the sequence upstream of a start codon: the
/// 6-mer at each window start, of which a shorter motif is the low bits.
///
/// The windows start 21 to 6 bases upstream of the start codon. The sites of
/// a set of nodes do not depend on the motif weights, so start training,
/// which searches every start node's motifs once per iteration, reads them
/// from the sequence once.
#[derive(Debug, Clone, Copy)]
pub(crate) struct MotifSite {
    /// 6-mer at window start `start - 21 + slot` (bases past the sequence
    /// end read as zero); unspecified for windows before the sequence start.
    kmers: [u16; MOTIF_WINDOW_COUNT],
    /// One bit per candidate window that lies inside the sequence.
    inside: u64,
}

impl MotifSite {
    /// A site without a window.
    const EMPTY: Self = Self {
        kmers: [0; MOTIF_WINDOW_COUNT],
        inside: 0,
    };

    /// The site of the start codon at `start` of the strand `strand_sequence`.
    pub(crate) fn new(strand_sequence: &[u8], start: usize, sequence_length: usize) -> Self {
        let mut site = Self {
            kmers: [0; MOTIF_WINDOW_COUNT],
            inside: 0,
        };

        let first_window = (start as isize - MOTIF_WINDOW_UPSTREAM).max(0);
        let last_window = start as isize - MOTIF_WINDOW_DOWNSTREAM;
        if last_window >= first_window {
            let base_at = |position: usize| {
                if position < sequence_length {
                    usize::from(encoded_base_code(strand_sequence, position))
                } else {
                    0
                }
            };
            let first_window = first_window as usize;
            let first_slot =
                (first_window as isize - (start as isize - MOTIF_WINDOW_UPSTREAM)) as usize;
            let mut kmer = 0;
            for offset in 0..MOTIF_KMER_LENGTH {
                kmer |= base_at(first_window + offset) << (2 * offset);
            }
            site.kmers[first_slot] = kmer as u16;
            for (slot, position) in ((first_window + 1)..=(last_window as usize)).enumerate() {
                kmer = (kmer >> 2)
                    | (base_at(position + MOTIF_KMER_LENGTH - 1) << (2 * (MOTIF_KMER_LENGTH - 1)));
                site.kmers[first_slot + slot + 1] = kmer as u16;
            }
        }

        for (number, candidate) in MOTIF_CANDIDATES.iter().enumerate() {
            let window = start as isize - MOTIF_WINDOW_UPSTREAM + isize::from(candidate.slot as i8);
            let motif_length = usize::from(candidate.length_index) + 3;
            if window >= 0 && window as usize + motif_length <= sequence_length {
                site.inside |= 1 << number;
            }
        }
        site
    }

    /// The highest scoring motif of the site under the weights of
    /// `training`; in stage 2 of start training, and for gene prediction,
    /// only if it scores well.
    pub(crate) fn best_motif(&self, training: &Training, stage: usize) -> Motif {
        // Score every window, then find the first of the best. A window
        // outside the sequence scores NaN, which is never the best.
        let mut scores = [f64::NAN; MOTIF_CANDIDATE_COUNT];
        if self.inside == ALL_MOTIF_CANDIDATES {
            // The candidates in the order of `MOTIF_CANDIDATES`, spelled
            // out as loops of fixed length so that their slots, masks and
            // weight rows are constants in the generated code.
            for (group, group_scores) in scores.chunks_exact_mut(MOTIF_SPACERS).enumerate() {
                let length_index = MOTIF_LENGTHS - 1 - group;
                let weight_rows = &training.motif_weights[length_index];
                let motif_mask = (1usize << (2 * (length_index + 3))) - 1;
                for (nearer, score) in group_scores.iter_mut().enumerate() {
                    let candidate = &MOTIF_CANDIDATES[group * MOTIF_SPACERS + nearer];
                    *score = weight_rows[usize::from(candidate.space_index)]
                        [usize::from(self.kmers[usize::from(candidate.slot)]) & motif_mask];
                }
            }
        } else {
            for (number, (score, candidate)) in scores.iter_mut().zip(&MOTIF_CANDIDATES).enumerate()
            {
                if self.inside & (1 << number) != 0 {
                    *score = training.motif_weights[usize::from(candidate.length_index)]
                        [usize::from(candidate.space_index)][self.motif_index(candidate)];
                }
            }
        }

        // The scores do not depend on each other, so the best is looked for
        // in four interleaved lanes, each keeping the first of its best,
        // without one comparison waiting for the previous one.
        const LANES: usize = 4;
        let mut lane_score = [INITIAL_MAX_SCORE; LANES];
        let mut lane_best = [usize::MAX; LANES];
        for (group, group_scores) in scores.chunks_exact(LANES).enumerate() {
            for lane in 0..LANES {
                let better = group_scores[lane] > lane_score[lane];
                lane_score[lane] = if better {
                    group_scores[lane]
                } else {
                    lane_score[lane]
                };
                lane_best[lane] = if better {
                    group * LANES + lane
                } else {
                    lane_best[lane]
                };
            }
        }
        let mut max_sc = INITIAL_MAX_SCORE;
        let mut best = usize::MAX;
        for lane in 0..LANES {
            if lane_best[lane] != usize::MAX
                && (lane_score[lane] > max_sc
                    || (lane_score[lane] == max_sc && lane_best[lane] < best))
            {
                max_sc = lane_score[lane];
                best = lane_best[lane];
            }
        }

        // In stage 2, only accept good scoring motifs
        let is_stage_two = stage == 2;
        let is_poor_motif =
            max_sc == MIN_MOTIF_SCORE || max_sc < training.no_motif_weight + MOTIF_THRESHOLD_OFFSET;

        // Do NOT neutralize the score in early stages when no window is found.
        // C keeps max_sc at -100.0 in this case, strongly disfavoring edge starts.
        if is_stage_two && is_poor_motif {
            Motif {
                index: 0,
                length: 0,
                space_index: 0,
                spacer: 0,
                score: training.no_motif_weight,
            }
        } else if let Some(candidate) = MOTIF_CANDIDATES.get(best) {
            Motif {
                index: self.motif_index(candidate),
                length: usize::from(candidate.length_index) + 3,
                space_index: usize::from(candidate.space_index),
                spacer: usize::from(candidate.spacer),
                score: max_sc,
            }
        } else {
            Motif {
                index: 0,
                length: 0,
                space_index: 0,
                spacer: 0,
                score: max_sc,
            }
        }
    }

    /// Adds one to the count of the motif in every window of the site, for
    /// every spacer class: what the first stage of start training counts for
    /// a start node.
    pub(crate) fn count_windows(&self, motif_counts: &mut MotifWeights) {
        for (number, candidate) in MOTIF_CANDIDATES.iter().enumerate() {
            if self.inside & (1 << number) != 0 {
                let motif_index = self.motif_index(candidate);
                for spacer_counts in &mut motif_counts[usize::from(candidate.length_index)] {
                    spacer_counts[motif_index] += 1.0;
                }
            }
        }
    }

    /// Index of the motif in `candidate`'s window: the low bits of the
    /// 6-mer that starts there.
    #[inline]
    fn motif_index(&self, candidate: &MotifCandidate) -> usize {
        let motif_mask = (1usize << (2 * (usize::from(candidate.length_index) + 3))) - 1;
        usize::from(self.kmers[usize::from(candidate.slot)]) & motif_mask
    }
}

/// The motif sites of `nodes`, one per node: empty for the stop nodes, whose
/// motifs are never searched.
pub(crate) fn motif_sites(
    seq: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    nodes: &[Node],
) -> Vec<MotifSite> {
    nodes
        .iter()
        .map(|node| {
            if node.position.codon_type == CodonType::Stop {
                MotifSite::EMPTY
            } else {
                motif_site(
                    seq,
                    reverse_complement_encoded_sequence,
                    sequence_length,
                    node,
                )
            }
        })
        .collect()
}

/// The motif site of `node`, a start node.
pub(crate) fn motif_site(
    seq: &[u8],
    reverse_complement_encoded_sequence: &[u8],
    sequence_length: usize,
    node: &Node,
) -> MotifSite {
    let (wseq, start) = match node.position.strand {
        Strand::Forward => (seq, node.position.index),
        Strand::Reverse => (
            reverse_complement_encoded_sequence,
            sequence_length - 1 - node.position.index,
        ),
        Strand::Unknown => unreachable!(),
    };
    MotifSite::new(wseq, start, sequence_length)
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

    let site = motif_site(
        seq,
        reverse_complement_encoded_sequence,
        sequence_length,
        node,
    );
    node.motif_info.best_motif = site.best_motif(training, stage);
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

    /// The motif search as Prodigal writes it: every motif length from six
    /// bases down to three, every window from the farthest to the nearest.
    fn best_motif_reference(
        training: &Training,
        strand_sequence: &[u8],
        start: usize,
        sequence_length: usize,
        stage: usize,
    ) -> Motif {
        let mut best = Motif {
            score: INITIAL_MAX_SCORE,
            ..Motif::default()
        };
        for length_index in (0..=3isize).rev() {
            let motif_length = (length_index + 3) as usize;
            for window in (start as isize - 18 - length_index)..=(start as isize - 6 - length_index)
            {
                if window < 0 || window as usize + motif_length > sequence_length {
                    continue;
                }
                let spacer = start as isize - window - length_index - 3;
                let space_index = if window <= start as isize - 16 - length_index {
                    3
                } else if window <= start as isize - 14 - length_index {
                    2
                } else if window >= start as isize - 7 - length_index {
                    1
                } else {
                    0
                };
                let index = calculate_kmer_index(motif_length, strand_sequence, window as usize);
                let score = training.motif_weights[length_index as usize][space_index][index];
                if score > best.score {
                    best = Motif {
                        index,
                        length: motif_length,
                        space_index,
                        spacer: spacer as usize,
                        score,
                    };
                }
            }
        }
        if stage == 2
            && (best.score == MIN_MOTIF_SCORE
                || best.score < training.no_motif_weight + MOTIF_THRESHOLD_OFFSET)
        {
            Motif {
                score: training.no_motif_weight,
                ..Motif::default()
            }
        } else {
            best
        }
    }

    #[test]
    fn motif_site_search_matches_the_window_by_window_search() {
        use crate::sequence::encoded::EncodedSequence;

        let mut random = crate::test_random::XorShift(0x94D0_49BB_1331_11EB_u64);
        let mut next = || random.next();
        let length = 257;
        let bases: Vec<u8> = (0..length)
            .map(|_| b"ACGT"[(next() % 4) as usize])
            .collect();
        let encoded = EncodedSequence::without_masking(&bases);

        // Weights with few distinct values, so that many windows tie and
        // the first of the best has to be the one found; one set in which
        // every motif is poor.
        let mut trainings = Vec::new();
        for distinct in [1u64, 2, 3, 17, 1 << 20] {
            let mut training = create_test_training();
            for weight in training.motif_weights.iter_mut().flatten().flatten() {
                *weight = (next() % distinct) as f64 / distinct as f64 * 8.0 - 4.0;
            }
            training.no_motif_weight = (next() % 5) as f64 - 2.0;
            trainings.push(training);
        }

        for training in &trainings {
            for strand_sequence in [
                &encoded.forward_sequence,
                &encoded.reverse_complement_sequence,
            ] {
                // Every start position, including those too close to the
                // sequence start for some or all of the windows.
                for start in 0..length {
                    let site = MotifSite::new(strand_sequence, start, length);
                    for stage in 0..3 {
                        assert_eq!(
                            site.best_motif(training, stage),
                            best_motif_reference(training, strand_sequence, start, length, stage),
                            "start {start}, stage {stage}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn rbs_lookup_matches_the_searches_for_every_kmer() {
        let mut random = crate::test_random::XorShift(0x853C_49E6_748F_EA9B_u64);
        let mut next = || random.next();
        let mut weight_sets = vec![[0.0f64; 28], [1.0f64; 28]];
        for _ in 0..3 {
            let mut weights = [0.0f64; 28];
            for weight in weights.iter_mut() {
                // A few distinct values, so that ties between motifs occur.
                *weight = (next() % 9) as f64 - 4.0;
            }
            weight_sets.push(weights);
        }

        for weights in &weight_sets {
            let lookup = RbsLookup::new(weights);
            for kmer in 0..RbsLookup::KMER_COUNT {
                let packed = pack_kmer(kmer);
                for distance in 0..=RbsLookup::MAX_DISTANCE {
                    assert_eq!(
                        lookup.lookup(kmer, distance),
                        calculate_rbs_scores(&packed, 0, distance, weights),
                        "kmer {kmer}, distance {distance}"
                    );
                }
            }
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
