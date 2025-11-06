use criterion::Criterion;
use std::time::Duration;

pub fn configure_criterion() -> Criterion {
    Criterion::default()
        .measurement_time(Duration::from_secs(300))
        .warm_up_time(Duration::from_secs(3))
        .sample_size(10)
        .significance_level(0.05)
        .noise_threshold(0.02)
}
