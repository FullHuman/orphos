use super::*;
use crate::constants::{
    MAX_MOTIF_LENGTH, MAX_RIBOSOME_DISTANCE, MIN_CUMULATIVE_SCORE, MIN_DISTANCE_FROM_START,
    MIN_MOTIF_LENGTH, READING_FRAMES, SLIDING_WINDOW_SIZE,
};

/// Calculate most GC-rich frame for each position
///
/// For every position `p` this counts the G/C bases in the same reading
/// frame within half a sliding window on either side (`|q - p| <= 57` for the
/// 120 bp window), and assigns each triplet the frame with the highest count.
pub fn calc_most_gc_frame(sequence: &[u8], sequence_length: usize) -> Vec<i32> {
    if sequence_length < READING_FRAMES {
        return vec![-1; sequence_length];
    }

    let total_gc_counts = calculate_window_gc_counts(sequence, sequence_length);

    assign_gc_rich_frames(&total_gc_counts, sequence_length)
}

/// Per position, the number of G/C bases in the same frame at most
/// `SLIDING_WINDOW_SIZE / 2 - READING_FRAMES` bases away.
fn calculate_window_gc_counts(sequence: &[u8], sequence_length: usize) -> Vec<i32> {
    let half_window = SLIDING_WINDOW_SIZE / 2;
    let reach = half_window - READING_FRAMES;
    let is_gc_at: Vec<i32> = (0..sequence_length)
        .map(|position| i32::from(is_gc(sequence, position)))
        .collect();

    let mut counts = vec![0; sequence_length];
    for position in 0..sequence_length {
        let mut count = if position < READING_FRAMES {
            // Start of the frame: count the bases at or after `position`.
            (position..=(position + reach).min(sequence_length - 1))
                .step_by(READING_FRAMES)
                .map(|q| is_gc_at[q])
                .sum()
        } else {
            counts[position - READING_FRAMES]
        };
        if position >= READING_FRAMES {
            // Slide the window by one codon within the frame.
            if position + reach < sequence_length {
                count += is_gc_at[position + reach];
            }
            if position >= half_window {
                count -= is_gc_at[position - half_window];
            }
        }
        counts[position] = count;
    }

    counts
}

fn assign_gc_rich_frames(total_gc_counts: &[i32], sequence_length: usize) -> Vec<i32> {
    let mut gc_rich_frames = vec![-1; sequence_length];

    for triplet_start in (0..sequence_length.saturating_sub(2)).step_by(READING_FRAMES) {
        let counts = [
            total_gc_counts[triplet_start],
            total_gc_counts.get(triplet_start + 1).copied().unwrap_or(0),
            total_gc_counts.get(triplet_start + 2).copied().unwrap_or(0),
        ];

        let max_gc_frame = find_max_reading_frame(counts[0], counts[1], counts[2]) as i32;

        for frame_offset in 0..READING_FRAMES.min(sequence_length - triplet_start) {
            gc_rich_frames[triplet_start + frame_offset] = max_gc_frame;
        }
    }

    gc_rich_frames
}

/// Find exact Shine-Dalgarno motifs
#[must_use]
pub fn shine_dalgarno_exact(
    sequence: &[u8],
    search_position: usize,
    start_codon_position: usize,
    ribosome_weights: &[f64],
) -> usize {
    if start_codon_position <= search_position + MIN_DISTANCE_FROM_START {
        return 0;
    }

    let search_limit =
        MAX_MOTIF_LENGTH.min(start_codon_position - MIN_DISTANCE_FROM_START - search_position);
    let base_scores = calculate_exact_base_scores(sequence, search_position, search_limit);

    find_best_exact_motif(
        &base_scores,
        search_limit,
        search_position,
        start_codon_position,
        ribosome_weights,
    )
}

fn calculate_exact_base_scores(
    sequence: &[u8],
    search_position: usize,
    search_limit: usize,
) -> [f64; MAX_MOTIF_LENGTH] {
    let mut base_scores = [0.0; MAX_MOTIF_LENGTH];

    for (pattern_index, score) in base_scores.iter_mut().take(search_limit).enumerate() {
        let sequence_position = search_position + pattern_index;
        *score = match pattern_index % 3 {
            0 if is_a(sequence, sequence_position) => 2.0,
            1 | 2 if is_g(sequence, sequence_position) => 3.0,
            _ => -10.0,
        };
    }

    base_scores
}

fn find_best_exact_motif(
    base_scores: &[f64],
    search_limit: usize,
    search_position: usize,
    start_codon_position: usize,
    ribosome_weights: &[f64],
) -> usize {
    let mut score_prefix = [0.0; MAX_MOTIF_LENGTH + 1];
    let mut mismatch_prefix = [0_u8; MAX_MOTIF_LENGTH + 1];

    for position in 0..search_limit {
        let score = base_scores[position];
        score_prefix[position + 1] = score_prefix[position] + score;
        mismatch_prefix[position + 1] = mismatch_prefix[position] + u8::from(score < 0.0);
    }

    let mut best_motif_index = 0;
    let start_distance = start_codon_position - search_position;

    for motif_length in (MIN_MOTIF_LENGTH..=search_limit).rev() {
        for motif_start_offset in 0..=(search_limit - motif_length) {
            let end = motif_start_offset + motif_length;

            if mismatch_prefix[end] != mismatch_prefix[motif_start_offset] {
                continue;
            }

            let cumulative_score = score_prefix[end] - score_prefix[motif_start_offset] - 2.0;
            let ribosome_distance = start_distance - end;

            if ribosome_distance > MAX_RIBOSOME_DISTANCE || cumulative_score < MIN_CUMULATIVE_SCORE
            {
                continue;
            }

            let motif_index = map_score_to_motif_index(
                cumulative_score as i32,
                categorize_distance(ribosome_distance, motif_length),
            );

            if is_better_motif(motif_index, best_motif_index, ribosome_weights) {
                best_motif_index = motif_index;
            }
        }
    }

    best_motif_index
}

/// Find Shine-Dalgarno motifs with single mismatch
#[must_use]
pub fn shine_dalgarno_mm(
    sequence: &[u8],
    search_position: usize,
    start_codon_position: usize,
    ribosome_weights: &[f64],
) -> usize {
    if start_codon_position <= search_position + MIN_DISTANCE_FROM_START {
        return 0;
    }

    let search_limit =
        MAX_MOTIF_LENGTH.min(start_codon_position - MIN_DISTANCE_FROM_START - search_position);
    let base_scores = calculate_mismatch_base_scores(sequence, search_position, search_limit);

    find_best_mismatch_motif(
        &base_scores,
        search_limit,
        search_position,
        start_codon_position,
        ribosome_weights,
    )
}

fn calculate_mismatch_base_scores(
    sequence: &[u8],
    search_position: usize,
    search_limit: usize,
) -> [f64; MAX_MOTIF_LENGTH] {
    let mut base_scores = [0.0; MAX_MOTIF_LENGTH];

    for (pattern_index, score) in base_scores.iter_mut().take(search_limit).enumerate() {
        let sequence_position = search_position + pattern_index;
        *score = match pattern_index % 3 {
            0 if is_a(sequence, sequence_position) => 2.0,
            0 => -3.0,
            _ if is_g(sequence, sequence_position) => 3.0,
            _ => -2.0,
        };
    }

    base_scores
}

fn find_best_mismatch_motif(
    base_scores: &[f64],
    search_limit: usize,
    search_position: usize,
    start_codon_position: usize,
    ribosome_weights: &[f64],
) -> usize {
    let mut score_prefix = [0.0; MAX_MOTIF_LENGTH + 1];
    let mut mismatch_prefix = [0_u8; MAX_MOTIF_LENGTH + 1];

    for position in 0..search_limit {
        let score = base_scores[position];
        score_prefix[position + 1] = score_prefix[position] + score;
        mismatch_prefix[position + 1] = mismatch_prefix[position] + u8::from(score < 0.0);
    }

    let mut best_motif_index = 0;
    let start_distance = start_codon_position - search_position;

    for motif_length in (5..=search_limit).rev() {
        for motif_start_offset in 0..=(search_limit - motif_length) {
            let end = motif_start_offset + motif_length;

            if mismatch_prefix[end] - mismatch_prefix[motif_start_offset] != 1 {
                continue;
            }

            let mut cumulative_score = score_prefix[end] - score_prefix[motif_start_offset] - 2.0;

            if base_scores[motif_start_offset] < 0.0
                || base_scores[motif_start_offset + 1] < 0.0
                || base_scores[end - 2] < 0.0
                || base_scores[end - 1] < 0.0
            {
                cumulative_score -= 10.0;
            }

            let ribosome_distance = start_distance - end;
            if ribosome_distance > MAX_RIBOSOME_DISTANCE || cumulative_score < MIN_CUMULATIVE_SCORE
            {
                continue;
            }

            let motif_index = map_mismatch_score_to_motif_index(
                cumulative_score as i32,
                categorize_mismatch_distance(ribosome_distance),
            );

            if is_better_motif(motif_index, best_motif_index, ribosome_weights) {
                best_motif_index = motif_index;
            }
        }
    }

    best_motif_index
}

const fn categorize_distance(ribosome_distance: usize, motif_length: usize) -> usize {
    match ribosome_distance {
        0..=4 => {
            if motif_length < 5 {
                2
            } else {
                1
            }
        }
        5..=10 => 0,
        11..=12 => {
            if motif_length < 5 {
                1
            } else {
                2
            }
        }
        _ => 3,
    }
}

const fn categorize_mismatch_distance(ribosome_distance: usize) -> usize {
    match ribosome_distance {
        0..=4 => 1,
        5..=10 => 0,
        11..=12 => 2,
        _ => 3,
    }
}

const fn map_score_to_motif_index(score: i32, distance_category: usize) -> usize {
    match (score, distance_category) {
        (6, 2) => 1,
        (6, 3) => 2,
        (8 | 9, 3) => 3,
        (6, 1) => 6,
        (11 | 12 | 14, 3) => 10,
        (8 | 9, 2) => 11,
        (8 | 9, 1) => 12,
        (6, 0) => 13,
        (8, 0) => 15,
        (9, 0) => 16,
        (11 | 12, 2) => 20,
        (11, 1) => 21,
        (11, 0) => 22,
        (12, 1) => 23,
        (12, 0) => 24,
        (14, 2) => 25,
        (14, 1) => 26,
        (14, 0) => 27,
        _ => 0,
    }
}

const fn map_mismatch_score_to_motif_index(score: i32, distance_category: usize) -> usize {
    match (score, distance_category) {
        (6 | 7, 3) => 2,
        (9, 3) => 3,
        (6, 2) => 4,
        (6, 1) => 5,
        (6, 0) => 9,
        (7, 2) => 7,
        (7, 1) => 8,
        (7, 0) => 14,
        (9, 2) => 17,
        (9, 1) => 18,
        (9, 0) => 19,
        _ => 0,
    }
}

fn is_better_motif(current_index: usize, best_index: usize, ribosome_weights: &[f64]) -> bool {
    ribosome_weights[current_index] > ribosome_weights[best_index]
        || (ribosome_weights[current_index] == ribosome_weights[best_index]
            && current_index > best_index)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::encoded::EncodedSequence;

    /// Reference implementation of the GC frame plot with explicit forward,
    /// backward and windowed counts, as previously implemented.
    fn calc_most_gc_frame_reference(sequence: &[u8], sequence_length: usize) -> Vec<i32> {
        if sequence_length < READING_FRAMES {
            return vec![-1; sequence_length];
        }
        let mut forward = vec![0; sequence_length];
        let mut backward = vec![0; sequence_length];
        for frame in 0..READING_FRAMES {
            for position in (frame..sequence_length).step_by(READING_FRAMES) {
                forward[position] = if position < READING_FRAMES {
                    i32::from(is_gc(sequence, position))
                } else {
                    forward[position - READING_FRAMES] + i32::from(is_gc(sequence, position))
                };
                let reverse_position = sequence_length - position - 1;
                backward[reverse_position] = if position < READING_FRAMES {
                    i32::from(is_gc(sequence, reverse_position))
                } else {
                    backward[reverse_position + READING_FRAMES]
                        + i32::from(is_gc(sequence, reverse_position))
                };
            }
        }
        let total: Vec<i32> = (0..sequence_length)
            .map(|position| {
                let mut total =
                    forward[position] + backward[position] - i32::from(is_gc(sequence, position));
                if position >= SLIDING_WINDOW_SIZE / 2 {
                    total -= forward[position - SLIDING_WINDOW_SIZE / 2];
                }
                if position + SLIDING_WINDOW_SIZE / 2 < sequence_length {
                    total -= backward[position + SLIDING_WINDOW_SIZE / 2];
                }
                total
            })
            .collect();
        assign_gc_rich_frames(&total, sequence_length)
    }

    #[test]
    fn test_calc_most_gc_frame_matches_reference() {
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for length in [
            3usize, 4, 5, 59, 60, 61, 119, 120, 121, 122, 123, 500, 1001, 4096,
        ] {
            for _ in 0..4 {
                let bases: Vec<u8> = (0..length)
                    .map(|_| b"ACGT"[(next() % 4) as usize])
                    .collect();
                let encoded = EncodedSequence::without_masking(&bases);
                assert_eq!(
                    calc_most_gc_frame(&encoded.forward_sequence, length),
                    calc_most_gc_frame_reference(&encoded.forward_sequence, length),
                    "length {length}"
                );
            }
        }
    }

    #[test]
    fn test_calc_most_gc_frame_basic() {
        let sequence = b"ATCGGCGCGCTAATCGGCGC";
        let result = calc_most_gc_frame(sequence, sequence.len());
        assert_eq!(result.len(), sequence.len());
        for &frame in &result {
            assert!((-1..=2).contains(&frame));
        }
    }

    #[test]
    fn test_calc_most_gc_frame_encoded_sequence() {
        let encoded = EncodedSequence::without_masking(b"ATCGGCGCGCTAATCGGCGC");
        let result = calc_most_gc_frame(&encoded.forward_sequence, encoded.sequence_length);

        assert_eq!(result.len(), encoded.sequence_length);
        for &frame in &result {
            assert!((-1..=2).contains(&frame));
        }
    }

    #[test]
    fn test_calc_most_gc_frame_empty() {
        let sequence = b"";
        let result = calc_most_gc_frame(sequence, 0);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_shine_dalgarno_exact_basic() {
        let ribosome_weights = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let sequence = b"AGGAGGTATGATCGGC";

        let result = shine_dalgarno_exact(sequence, 5, 12, &ribosome_weights);
        assert!(result < ribosome_weights.len());
    }

    #[test]
    fn test_shine_dalgarno_exact_encoded_sequence() {
        let ribosome_weights = vec![1.0; 28];
        let encoded = EncodedSequence::without_masking(b"AGGAGGTATGATCGGC");

        let result = shine_dalgarno_exact(&encoded.forward_sequence, 0, 12, &ribosome_weights);
        assert!(result < ribosome_weights.len());
    }

    #[test]
    fn test_shine_dalgarno_mm_basic() {
        let ribosome_weights = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let sequence = b"AGGAGGTATGATCGGC";

        let result = shine_dalgarno_mm(sequence, 0, 8, &ribosome_weights);
        assert!(result < ribosome_weights.len());
    }

    #[test]
    fn test_shine_dalgarno_mm_encoded_sequence() {
        let ribosome_weights = vec![1.0; 28];
        let encoded = EncodedSequence::without_masking(b"AGGAGGTATGATCGGC");

        let result = shine_dalgarno_mm(&encoded.forward_sequence, 0, 12, &ribosome_weights);
        assert!(result < ribosome_weights.len());
    }

    #[test]
    fn test_assign_gc_rich_frames_basic() {
        let total_gc_counts = vec![5, 3, 8, 2, 7, 1, 9, 4];
        let result = assign_gc_rich_frames(&total_gc_counts, 8);
        assert_eq!(result.len(), 8);

        for &frame in &result {
            assert!((-1..=2).contains(&frame));
        }
    }
}
