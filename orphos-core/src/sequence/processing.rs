use super::*;
use crate::constants::{
    MAX_MOTIF_LENGTH, MAX_RIBOSOME_DISTANCE, MIN_CUMULATIVE_SCORE, MIN_DISTANCE_FROM_START,
    MIN_MOTIF_LENGTH, READING_FRAMES, SLIDING_WINDOW_SIZE,
};

/// Calculate most GC-rich frame for each position
pub fn calc_most_gc_frame(sequence: &[u8], sequence_length: usize) -> Vec<i32> {
    if sequence_length < READING_FRAMES {
        return vec![-1; sequence_length];
    }

    let forward_gc_counts = calculate_forward_gc_counts(sequence, sequence_length);
    let backward_gc_counts = calculate_backward_gc_counts(sequence, sequence_length);
    let total_gc_counts = calculate_total_gc_counts(
        &forward_gc_counts,
        &backward_gc_counts,
        sequence,
        sequence_length,
    );

    assign_gc_rich_frames(&total_gc_counts, sequence_length)
}

fn calculate_forward_gc_counts(sequence: &[u8], sequence_length: usize) -> Vec<i32> {
    let mut counts = vec![0; sequence_length];

    for reading_frame in 0..READING_FRAMES {
        for position in (reading_frame..sequence_length).step_by(READING_FRAMES) {
            counts[position] = if position < READING_FRAMES {
                i32::from(is_gc(sequence, position))
            } else {
                counts[position - READING_FRAMES] + i32::from(is_gc(sequence, position))
            };
        }
    }

    counts
}

fn calculate_backward_gc_counts(sequence: &[u8], sequence_length: usize) -> Vec<i32> {
    let mut counts = vec![0; sequence_length];

    for reading_frame in 0..READING_FRAMES {
        for position in (reading_frame..sequence_length).step_by(READING_FRAMES) {
            let reverse_position = sequence_length - position - 1;
            counts[reverse_position] = if position < READING_FRAMES {
                i32::from(is_gc(sequence, reverse_position))
            } else {
                counts[reverse_position + READING_FRAMES]
                    + i32::from(is_gc(sequence, reverse_position))
            };
        }
    }

    counts
}

fn calculate_total_gc_counts(
    forward_counts: &[i32],
    backward_counts: &[i32],
    sequence: &[u8],
    sequence_length: usize,
) -> Vec<i32> {
    (0..sequence_length)
        .map(|position| {
            let mut total = forward_counts[position] + backward_counts[position]
                - i32::from(is_gc(sequence, position));

            if position >= SLIDING_WINDOW_SIZE / 2 {
                total -= forward_counts[position - SLIDING_WINDOW_SIZE / 2];
            }

            if position + SLIDING_WINDOW_SIZE / 2 < sequence_length {
                total -= backward_counts[position + SLIDING_WINDOW_SIZE / 2];
            }

            total
        })
        .collect()
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
) -> Vec<f64> {
    (0..search_limit)
        .map(|pattern_index| {
            let sequence_position = search_position + pattern_index;
            match pattern_index % 3 {
                0 if is_a(sequence, sequence_position) => 2.0,
                1 | 2 if is_g(sequence, sequence_position) => 3.0,
                _ => -10.0,
            }
        })
        .collect()
}

fn find_best_exact_motif(
    base_scores: &[f64],
    search_limit: usize,
    search_position: usize,
    start_codon_position: usize,
    ribosome_weights: &[f64],
) -> usize {
    let mut best_motif_index = 0;

    for motif_length in (MIN_MOTIF_LENGTH..=search_limit).rev() {
        for motif_start_offset in 0..=(search_limit - motif_length) {
            if let Some(motif_index) = evaluate_exact_motif(
                base_scores,
                motif_start_offset,
                motif_length,
                search_position,
                start_codon_position,
            ) && is_better_motif(motif_index, best_motif_index, ribosome_weights)
            {
                best_motif_index = motif_index;
            }
        }
    }

    best_motif_index
}

fn evaluate_exact_motif(
    base_scores: &[f64],
    motif_start_offset: usize,
    motif_length: usize,
    search_position: usize,
    start_codon_position: usize,
) -> Option<usize> {
    let start = motif_start_offset;
    let end = motif_start_offset + motif_length;
    let window = &base_scores[start..end];

    if window.iter().any(|&score| score < 0.0) {
        return None;
    }

    let cumulative_score: f64 = window.iter().copied().sum::<f64>() - 2.0;
    let ribosome_distance =
        start_codon_position - (search_position + motif_start_offset + motif_length);

    if ribosome_distance > MAX_RIBOSOME_DISTANCE || cumulative_score < MIN_CUMULATIVE_SCORE {
        return None;
    }

    let distance_category = categorize_distance(ribosome_distance, motif_length);
    Some(map_score_to_motif_index(
        cumulative_score as i32,
        distance_category,
    ))
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
) -> Vec<f64> {
    (0..search_limit)
        .map(|pattern_index| {
            let sequence_position = search_position + pattern_index;
            match pattern_index % 3 {
                0 => {
                    if is_a(sequence, sequence_position) {
                        2.0
                    } else {
                        -3.0
                    }
                }
                _ => {
                    if is_g(sequence, sequence_position) {
                        3.0
                    } else {
                        -2.0
                    }
                }
            }
        })
        .collect()
}

fn find_best_mismatch_motif(
    base_scores: &[f64],
    search_limit: usize,
    search_position: usize,
    start_codon_position: usize,
    ribosome_weights: &[f64],
) -> usize {
    let mut best_motif_index = 0;

    for motif_length in (5..=search_limit).rev() {
        for motif_start_offset in 0..=(search_limit - motif_length) {
            if let Some(motif_index) = evaluate_mismatch_motif(
                base_scores,
                motif_start_offset,
                motif_length,
                search_position,
                start_codon_position,
            ) && is_better_motif(motif_index, best_motif_index, ribosome_weights)
            {
                best_motif_index = motif_index;
            }
        }
    }

    best_motif_index
}

fn evaluate_mismatch_motif(
    base_scores: &[f64],
    motif_start_offset: usize,
    motif_length: usize,
    search_position: usize,
    start_codon_position: usize,
) -> Option<usize> {
    let mut cumulative_score = -2.0;
    let mut mismatch_count = 0;

    let start = motif_start_offset;
    let end = motif_start_offset + motif_length;
    for (pos_in_motif, &score) in base_scores[start..end].iter().enumerate() {
        cumulative_score += score;
        if score < 0.0 {
            mismatch_count += 1;
            if pos_in_motif <= 1 || pos_in_motif >= motif_length - 2 {
                cumulative_score -= 10.0;
            }
        }
    }

    if mismatch_count != 1 {
        return None;
    }

    let ribosome_distance =
        start_codon_position - (search_position + motif_start_offset + motif_length);

    if ribosome_distance > MAX_RIBOSOME_DISTANCE || cumulative_score < MIN_CUMULATIVE_SCORE {
        return None;
    }

    let distance_category = categorize_mismatch_distance(ribosome_distance);
    Some(map_mismatch_score_to_motif_index(
        cumulative_score as i32,
        distance_category,
    ))
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
    fn test_shine_dalgarno_mm_basic() {
        let ribosome_weights = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let sequence = b"AGGAGGTATGATCGGC";

        let result = shine_dalgarno_mm(sequence, 0, 8, &ribosome_weights);
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
