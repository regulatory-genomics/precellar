//! Rayon thread pool shared by all parallel stages of one pipeline run.
//!
//! [`with_threads`] installs a pool of exactly `num_threads` threads for the
//! calling thread. While it runs, [`run`] executes parallel work in that pool,
//! so every stage (barcode annotation, alignment, transcriptome annotation,
//! sorting, counting) is capped at the thread count the user asked for instead
//! of using rayon's global pool, which spans all cores of the machine.
//!
//! A helper thread, such as the alignment producer thread, takes the pool with
//! [`current`] on the spawning thread and installs it with [`with_pool`].
//! Stages that build their own pool, such as the external sorter, size it with
//! [`current_num_threads`].

use std::cell::RefCell;
use std::sync::Arc;

use anyhow::{Context, Result};

thread_local! {
    static POOL: RefCell<Option<Arc<rayon::ThreadPool>>> = const { RefCell::new(None) };
}

/// Resets the calling thread's pool when dropped, including on unwind.
struct Restore(Option<Arc<rayon::ThreadPool>>);

impl Drop for Restore {
    fn drop(&mut self) {
        let previous = self.0.take();
        POOL.with(|pool| *pool.borrow_mut() = previous);
    }
}

/// Builds a pool of `num_threads` threads (at least 1).
pub fn build_pool(num_threads: usize) -> Result<rayon::ThreadPool> {
    let num_threads = num_threads.max(1);
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .thread_name(|i| format!("precellar-{i}"))
        .build()
        .with_context(|| format!("failed to build a pool of {num_threads} threads"))
}

/// Runs `f` with a pool of `num_threads` threads installed for the calling
/// thread. The previous pool, if any, is restored afterwards.
pub fn with_threads<R>(num_threads: usize, f: impl FnOnce() -> R) -> Result<R> {
    let pool = Arc::new(build_pool(num_threads)?);
    Ok(with_pool(Some(pool), f))
}

/// Runs `f` with `pool` installed for the calling thread. The previous pool,
/// if any, is restored afterwards. A helper thread passes the [`current`] pool
/// of the thread that spawned it, so both threads share one pool.
pub fn with_pool<R>(pool: Option<Arc<rayon::ThreadPool>>, f: impl FnOnce() -> R) -> R {
    let previous = POOL.with(|slot| std::mem::replace(&mut *slot.borrow_mut(), pool));
    let _restore = Restore(previous);
    f()
}

/// Returns the pool installed by the innermost [`with_threads`] on the calling
/// thread, if any.
pub fn current() -> Option<Arc<rayon::ThreadPool>> {
    POOL.with(|pool| pool.borrow().clone())
}

/// Runs `op` in the installed pool, or in rayon's current pool when none is
/// installed. Rayon parallel iterators inside `op` use the same pool.
pub fn run<R: Send>(op: impl FnOnce() -> R + Send) -> R {
    match current() {
        Some(pool) => pool.install(op),
        None => op(),
    }
}

/// Number of threads that [`run`] executes on.
pub fn current_num_threads() -> usize {
    current().map_or_else(rayon::current_num_threads, |pool| {
        pool.current_num_threads()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rayon::prelude::*;
    use std::collections::HashSet;
    use std::sync::Mutex;

    #[test]
    fn run_uses_the_installed_pool_and_restores_it() {
        assert!(current().is_none());
        let used = with_threads(3, || {
            assert_eq!(current_num_threads(), 3);
            let threads = Mutex::new(HashSet::new());
            run(|| {
                (0..10_000).into_par_iter().for_each(|_| {
                    threads.lock().unwrap().insert(std::thread::current().id());
                })
            });
            threads.into_inner().unwrap().len()
        })
        .unwrap();
        assert!(used <= 3, "{used} threads for a pool of 3");
        assert!(current().is_none());
    }
}
