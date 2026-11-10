//! Deterministic pseudo-random numbers for tests.

/// An xorshift generator: the same numbers for the same seed, everywhere.
pub(crate) struct XorShift(pub(crate) u64);

impl XorShift {
    pub(crate) fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    /// A value in `[low, high)`.
    pub(crate) fn float(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * (self.next() % 1_000_000) as f64 / 1_000_000.0
    }
}
