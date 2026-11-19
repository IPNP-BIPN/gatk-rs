//! Work spread over threads, its results taken back in the order of the items (#1467): what lets
//! HaplotypeCaller prepare assembly regions on every core while everything observable (the
//! generator's draws, the genotyping engine's state, the writers) still happens region by region.

use std::collections::BTreeMap;
use std::panic::{catch_unwind, resume_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{mpsc, Mutex};

/// `RAYON_NUM_THREADS` when it is a positive count, else the cores: the one knob
/// `tools/coverage/thread_independence.py` turns, as for the Spark tools' BWA.
pub fn threads() -> usize {
    std::env::var("RAYON_NUM_THREADS")
        .ok()
        .and_then(|text| text.parse().ok())
        .filter(|threads| *threads > 0)
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |count| count.get()))
}

/// `finish(work(item))` for each item, in the items' order, as a plain loop would.
///
/// With more than one thread, `work` runs on `threads` threads, at most `2 * threads` items ahead
/// of the one being finished, which bounds the memory. The items are pulled, and `finish` runs, on
/// the calling thread only. The first error, of a pull or of `finish`, is returned once every
/// earlier item is finished and no later one has been: what the loop would have done, provided
/// `work` has no effect beyond its result. A panic in `work` is resumed here, at its item's turn.
pub fn map_in_order<T: Send, R: Send, E>(
    mut items: impl Iterator<Item = Result<T, E>>,
    threads: usize,
    work: impl Fn(T) -> R + Sync,
    mut finish: impl FnMut(R) -> Result<(), E>,
) -> Result<(), E> {
    if threads <= 1 {
        for item in items {
            finish(work(item?))?;
        }
        return Ok(());
    }
    let window = 2 * threads;
    let (jobs, queue) = mpsc::channel::<(usize, T)>();
    let queue = Mutex::new(queue);
    let (results, done) = mpsc::channel::<(usize, std::thread::Result<R>)>();
    // Set once the caller stops early: the jobs still queued are dropped, not worked.
    let stop = AtomicBool::new(false);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            let (queue, work, stop, results) = (&queue, &work, &stop, results.clone());
            scope.spawn(move || loop {
                // The guard dies with this statement: only the wait for a job is serialised.
                let job = queue.lock().expect("the job queue").recv();
                let Ok((index, item)) = job else {
                    break;
                };
                if stop.load(Ordering::Relaxed) {
                    continue;
                }
                let result = catch_unwind(AssertUnwindSafe(|| work(item)));
                if results.send((index, result)).is_err() {
                    break;
                }
            });
        }
        let mut pending = BTreeMap::new();
        let (mut pulled, mut next) = (0usize, 0usize);
        let mut more = true;
        let mut pull_error = None;
        let outcome = loop {
            while more && pull_error.is_none() && pulled - next < window {
                match items.next() {
                    Some(Ok(item)) => {
                        jobs.send((pulled, item))
                            .expect("the workers outlive the jobs");
                        pulled += 1;
                    }
                    Some(Err(error)) => pull_error = Some(error),
                    None => more = false,
                }
            }
            if next == pulled {
                break pull_error.map_or(Ok(()), Err);
            }
            let result = loop {
                if let Some(result) = pending.remove(&next) {
                    break result;
                }
                let (index, result) = done.recv().expect("a worker answers every job");
                pending.insert(index, result);
            };
            next += 1;
            let finished = match result {
                Ok(value) => finish(value),
                Err(panic) => {
                    stop.store(true, Ordering::Relaxed);
                    drop(jobs);
                    resume_unwind(panic)
                }
            };
            if let Err(error) = finished {
                break Err(error);
            }
        };
        stop.store(true, Ordering::Relaxed);
        drop(jobs);
        outcome
    })
}

#[cfg(test)]
mod tests {
    use super::map_in_order;

    /// Uneven work finishes in item order, an error stops at its item with every earlier one
    /// finished, at one thread and at several.
    #[test]
    fn in_order_and_stops_where_the_loop_would() {
        for threads in [1, 2, 7] {
            let items = (0..200u64).map(|i| if i == 150 { Err(i) } else { Ok(i) });
            let mut seen = Vec::new();
            let outcome = map_in_order(
                items,
                threads,
                |i| {
                    std::thread::sleep(std::time::Duration::from_micros((i * 7919) % 300));
                    i * 2
                },
                |r| {
                    seen.push(r);
                    if r == 2 * 120 {
                        Err(1000)
                    } else {
                        Ok(())
                    }
                },
            );
            assert_eq!(outcome, Err(1000));
            assert_eq!(seen, (0..=120).map(|i| i * 2).collect::<Vec<_>>());

            let mut seen = Vec::new();
            let items = (0..200u64).map(|i| if i == 150 { Err(i) } else { Ok(i) });
            let outcome = map_in_order(
                items,
                threads,
                |i| i,
                |r| {
                    seen.push(r);
                    Ok::<_, u64>(())
                },
            );
            assert_eq!(outcome, Err(150));
            assert_eq!(seen, (0..150).collect::<Vec<_>>());
        }
    }

    #[test]
    fn a_panic_in_work_reaches_the_caller() {
        let caught = std::panic::catch_unwind(|| {
            map_in_order(
                (0..50).map(Ok::<_, ()>),
                4,
                |i| assert_ne!(i, 30),
                |_| Ok(()),
            )
        });
        assert!(caught.is_err());
    }
}
