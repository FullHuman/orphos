//! Work that runs on two threads when the thread pool has two.

/// [`rayon::join`], except that on a pool of one thread both closures run
/// on the calling thread, `a` first, without a hand-off to the pool.
pub(crate) fn join<A, B, RA, RB>(a: A, b: B) -> (RA, RB)
where
    A: FnOnce() -> RA + Send,
    B: FnOnce() -> RB + Send,
    RA: Send,
    RB: Send,
{
    if rayon::current_num_threads() > 1 {
        rayon::join(a, b)
    } else {
        (a(), b())
    }
}

#[cfg(test)]
mod tests {
    use super::join;

    #[test]
    fn join_returns_both_results() {
        let (a, b) = join(|| 1 + 1, || "two");
        assert_eq!(a, 2);
        assert_eq!(b, "two");
    }

    #[test]
    fn join_runs_both_closures_on_a_pool_of_one_thread() {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let (a, b) = pool.install(|| join(rayon::current_num_threads, || 3));
        assert_eq!((a, b), (1, 3));
    }
}
