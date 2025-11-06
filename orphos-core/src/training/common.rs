/// Normalize array and convert to log ratios with clamping
pub fn normalize_and_log_transform(
    real_counts: &mut [f64],
    background: &[f64],
    weight_clamp_min: f64,
    weight_clamp_max: f64,
) {
    let total_real: f64 = real_counts.iter().sum();

    if total_real == 0.0 {
        real_counts.fill(0.0);
        return;
    }

    for count in real_counts.iter_mut() {
        *count /= total_real;
    }

    for (i, count) in real_counts.iter_mut().enumerate() {
        *count = if background[i] > 0.0 {
            (*count / background[i])
                .ln()
                .clamp(weight_clamp_min, weight_clamp_max)
        } else {
            weight_clamp_min
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_and_log_transform_basic() {
        let mut real_counts = vec![10.0, 20.0, 30.0];
        let background = vec![0.1, 0.2, 0.3];
        let weight_clamp_min = -10.0;
        let weight_clamp_max = 10.0;

        normalize_and_log_transform(
            &mut real_counts,
            &background,
            weight_clamp_min,
            weight_clamp_max,
        );

        // Check that values are within clamp range
        for &value in &real_counts {
            assert!(value >= weight_clamp_min);
            assert!(value <= weight_clamp_max);
        }

        assert_ne!(real_counts[0], 10.0);
        assert_ne!(real_counts[1], 20.0);
        assert_ne!(real_counts[2], 30.0);
    }

    #[test]
    fn test_normalize_and_log_transform_zero_total() {
        let mut real_counts = vec![0.0, 0.0, 0.0];
        let background = vec![0.1, 0.2, 0.3];
        let weight_clamp_min = -10.0;
        let weight_clamp_max = 10.0;

        normalize_and_log_transform(
            &mut real_counts,
            &background,
            weight_clamp_min,
            weight_clamp_max,
        );

        assert_eq!(real_counts, vec![0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_normalize_and_log_transform_zero_background() {
        let mut real_counts = vec![10.0, 20.0];
        let background = vec![0.0, 0.2]; // First background is zero
        let weight_clamp_min = -5.0;
        let weight_clamp_max = 5.0;

        normalize_and_log_transform(
            &mut real_counts,
            &background,
            weight_clamp_min,
            weight_clamp_max,
        );

        assert_eq!(real_counts[0], weight_clamp_min);

        assert!(real_counts[1] >= weight_clamp_min);
        assert!(real_counts[1] <= weight_clamp_max);
    }

    #[test]
    fn test_normalize_and_log_transform_clamping() {
        let mut real_counts = vec![1000.0, 0.001]; // Values that will produce extreme log ratios
        let background = vec![0.001, 1000.0];
        let weight_clamp_min = -2.0;
        let weight_clamp_max = 2.0;

        normalize_and_log_transform(
            &mut real_counts,
            &background,
            weight_clamp_min,
            weight_clamp_max,
        );

        for &value in &real_counts {
            assert!(value >= weight_clamp_min);
            assert!(value <= weight_clamp_max);
        }
    }

    #[test]
    fn test_normalize_and_log_transform_empty() {
        let mut real_counts: Vec<f64> = vec![];
        let background: Vec<f64> = vec![];
        let weight_clamp_min = -10.0;
        let weight_clamp_max = 10.0;

        // Should not panic with empty arrays
        normalize_and_log_transform(
            &mut real_counts,
            &background,
            weight_clamp_min,
            weight_clamp_max,
        );

        assert!(real_counts.is_empty());
    }

    #[test]
    fn test_normalize_and_log_transform_single_element() {
        let mut real_counts = vec![42.0];
        let background = vec![0.5];
        let weight_clamp_min = -10.0;
        let weight_clamp_max = 10.0;

        normalize_and_log_transform(
            &mut real_counts,
            &background,
            weight_clamp_min,
            weight_clamp_max,
        );

        // log(1.0 / 0.5) = log(2.0) ≈ 0.693
        assert!((real_counts[0] - 0.693).abs() < 0.001);
    }
}
