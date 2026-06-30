//! Minimal data-parallel `map` over a slice, using only `std` (no rayon dependency).
//!
//! Genome digestion, sketch construction and pairwise clustering are embarrassingly
//! parallel; `par_map` spreads a slice across scoped threads and returns results in input
//! order. Thread count = `STRAIN2BSCAN_THREADS` if set (used to benchmark single- vs
//! multi-threaded scaling), else `std::thread::available_parallelism()`.

use std::thread;

/// Resolve the worker count (env override → available parallelism → 1).
pub fn num_threads() -> usize {
    if let Ok(v) = std::env::var("STRAIN2BSCAN_THREADS") {
        if let Ok(n) = v.parse::<usize>() {
            if n >= 1 {
                return n;
            }
        }
    }
    thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

/// Parallel map preserving input order. `f` runs on `&T` and must be `Sync`.
pub fn par_map<T, R, F>(items: &[T], f: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> R + Sync,
{
    let n = items.len();
    if n == 0 {
        return Vec::new();
    }
    let nt = num_threads().min(n);
    if nt <= 1 {
        return items.iter().map(&f).collect();
    }
    let chunk = n.div_ceil(nt);
    let fr = &f;
    thread::scope(|s| {
        let handles: Vec<_> = items
            .chunks(chunk)
            .map(|c| s.spawn(move || c.iter().map(fr).collect::<Vec<R>>()))
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn par_map_preserves_order_and_values() {
        let v: Vec<u64> = (0..1000).collect();
        let got = par_map(&v, |&x| x * x);
        let want: Vec<u64> = (0..1000).map(|x| x * x).collect();
        assert_eq!(got, want);
    }

    #[test]
    fn par_map_empty() {
        let v: Vec<u64> = Vec::new();
        let got: Vec<u64> = par_map(&v, |&x| x);
        assert!(got.is_empty());
    }
}
