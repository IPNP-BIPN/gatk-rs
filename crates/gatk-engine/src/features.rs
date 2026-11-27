//! Ported from `org.broadinstitute.hellbender.engine.FeatureCache`,
//! `org.broadinstitute.hellbender.engine.FeatureDataSource` (its query path) and
//! `org.broadinstitute.hellbender.engine.FeatureContext` (GATK 4.7.0.0; all three files are
//! unchanged since 4.6.2.0).
//!
//! What a tool sees when it asks "which variants are here" is not what the file contains at that
//! position: it is what a lookahead cache decided to keep. The cache is not an optimisation that a
//! port may skip, because two of its behaviours are observable in the answers:
//!
//!  * **a hit returns a subset of a wider prefetch, not a fresh query.** `refillQueryCache` extends
//!    the query end by `queryLookaheadBases` and caches everything it finds, so a later query
//!    inside that window is answered from memory. Since `fill` keeps whatever the reader returned,
//!    the answer to a narrow query is the *prefetched* set trimmed by start position, and
//!    `getCachedFeaturesUpToStopPosition` filters only on `start > stopPosition`. Neither step
//!    re-tests overlap with the query: the answer equals a fresh overlap query only because the
//!    trim drops every feature that starts before the new start and ends before it too, and the
//!    cache is sorted by start. `tests/data/feature_cache.txt.gz` shows the two agree: its ten
//!    queries return the same features at lookahead 0 and at lookahead 100;
//!  * **the trim preserves the reader's order rather than sorting.** `trimToNewStartPosition` pops
//!    features that start before the new start, keeps the ones that still overlap it, and pushes
//!    them back **in reverse**, so the relative order the file had survives. A port that re-sorted
//!    would agree on every set and differ on every list.
//!
//! The window arithmetic in `FeatureContext` is the third: the leading expansion is clamped at 1
//! and the trailing one is not clamped to the contig at all, because a Feature query is allowed to
//! run off the end.
//!
//! # Coordinates
//!
//! Every position here is 1-based and inclusive, Tribble's convention after the codec has
//! converted (a BED line `chr1 9 20` arrives as `chr1:10-20`). Nothing in this file converts.
//!
//! # Reading this file without reading Rust
//!
//! - `VecDeque` is Java's `ArrayDeque`: `front`, `pop_front`, `push_front` are `getFirst`,
//!   `removeFirst`, `addFirst`;
//! - `Option<T>` is "a `T` or nothing", Java's nullable reference with the null check enforced;
//!   `Result<T, String>` is "a `T` or an error message", how a reference exception is modelled;
//! - `?` after a `Result` returns its error to the caller at once;
//! - `let Some(x) = y else { ... }` unwraps `y` or takes the early exit in the braces;
//! - `a.checked_add(b).expect(..)` is `Math.addExact(a, b)`: `None` on overflow, and `expect`
//!   turns that into a panic, the Rust counterpart of the `ArithmeticException`;
//! - a `trait` is a Java interface, and `<R: FeatureReader>` is a generic bounded by it.

use crate::interval::SimpleInterval;
use std::collections::VecDeque;

/// `htsjdk.tribble.Feature`, reduced to what the cache and the context need: a location.
///
/// The payload is deliberately opaque. Decoding a VCF, a BED or an interval list is htsjdk's job
/// and belongs in `htsjdk-rs`; what belongs here is what GATK does with whatever comes back.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub contig: String,
    /// 1-based, inclusive, as Tribble reports it after its own codec has converted.
    pub start: i32,
    pub end: i32,
    /// An identifier carried through so a test can name which feature came back.
    pub name: String,
}

impl Feature {
    pub fn new(contig: &str, start: i32, end: i32, name: &str) -> Self {
        Feature {
            contig: contig.to_string(),
            start,
            end,
            name: name.to_string(),
        }
    }
}

/// `FeatureCache`: the features overlapping one interval, in reader order, plus that interval.
///
/// The hit and miss counters are the reference's `numCacheHits` / `numCacheMisses`, which it
/// reports through `printCacheStatistics`; they never change an answer.
#[derive(Debug, Default)]
pub struct FeatureCache {
    cache: VecDeque<Feature>,
    /// The interval the cache currently holds every overlapping feature for.
    cached_interval: Option<SimpleInterval>,
    hits: usize,
    misses: usize,
}

impl FeatureCache {
    pub fn new() -> Self {
        FeatureCache::default()
    }

    pub fn hits(&self) -> usize {
        self.hits
    }

    pub fn misses(&self) -> usize {
        self.misses
    }

    /// `FeatureCache.fill`: replace the contents wholesale, and record what they cover.
    ///
    /// Note what it does *not* do: it does not filter what the reader returned, so the cache holds
    /// whatever the query produced, including features that only overlap the lookahead tail.
    pub fn fill(&mut self, features: Vec<Feature>, interval: SimpleInterval) {
        // `.into()` converts the `Vec` into the deque, keeping the order; Java's `cache.clear()`
        // then `cache.add` per feature.
        self.cache = features.into();
        self.cached_interval = Some(interval);
    }

    /// `FeatureCache.cacheHit`: containment, not overlap. A query one base past the cached end is
    /// a miss even though almost everything it needs is already in memory.
    ///
    /// ```java
    /// final boolean cacheHit = cachedInterval != null && cachedInterval.contains(interval);
    /// ```
    ///
    /// `contains` also requires the same contig, which is the first test below. A query that
    /// starts before the cached start is a miss too, so the cache never answers backwards.
    pub fn cache_hit(&mut self, interval: &SimpleInterval) -> bool {
        // `match` on the `Option`: `Some(cached)` binds the interval, `None` is Java's null.
        let hit = match &self.cached_interval {
            Some(cached) => {
                cached.contig == interval.contig
                    && cached.start <= interval.start
                    && interval.end <= cached.end
            }
            None => false,
        };
        if hit {
            self.hits += 1;
        } else {
            self.misses += 1;
        }
        hit
    }

    /// `FeatureCache.trimToNewStartPosition`.
    ///
    /// Returns `Err` where the reference throws `GATKException`: trimming past the cached end is
    /// called a bug upstream rather than handled.
    ///
    /// ```java
    /// if ( featureBeforeNewStart.getEnd() >= newStart ) {
    /// ```
    ///
    /// is the keep test, and the survivors go back with `cache.addFirst` from the last to the
    /// first, so their relative order is the reader's.
    pub fn trim_to_new_start_position(&mut self, new_start: i32) -> Result<(), String> {
        // The reference dereferences `cachedInterval` without a null check, which would be a
        // `NullPointerException`; it is unreachable there because the only caller trims after a
        // hit. `ok_or_else(..)?` turns the `None` into an early `Err` instead.
        let cached = self
            .cached_interval
            .clone()
            .ok_or_else(|| "BUG: trimming an empty cache".to_string())?;
        if new_start > cached.end {
            return Err(format!(
                "BUG: attempted to trim Feature cache to an improper new start position ({new_start}). Cache stop = {}",
                cached.end
            ));
        }

        // Pop everything that starts before the new start, keeping those that still overlap it.
        // The loop stops at the first feature starting on or after the new start, because the
        // features are assumed sorted by start position: an unsorted file is not detected here.
        let mut kept: Vec<Feature> = Vec::new();
        // `while let Some(front) = ...` loops while the deque has a first element, binding it;
        // Java's `while ( ! cache.isEmpty() && ... )` with the `getFirst` folded in.
        while let Some(front) = self.cache.front() {
            if front.start >= new_start {
                break;
            }
            let feature = self.cache.pop_front().expect("checked above");
            if feature.end >= new_start {
                kept.push(feature);
            }
        }
        // Pushed back in reverse, which restores the order they were popped in.
        for feature in kept.into_iter().rev() {
            self.cache.push_front(feature);
        }

        self.cached_interval = Some(SimpleInterval {
            contig: cached.contig,
            start: new_start,
            end: cached.end,
        });
        Ok(())
    }

    /// `FeatureCache.getCachedFeaturesUpToStopPosition`.
    ///
    /// The only test is `start > stopPosition`, and the loop **breaks** rather than continuing, so
    /// one feature starting past the stop hides every feature behind it. Nothing here tests the
    /// end: excluding a feature that ends before the query started is left to
    /// [`FeatureCache::trim_to_new_start_position`], which drops each one it meets while the
    /// cache is sorted by start.
    ///
    /// ```java
    /// if ( candidateFeature.getStart() > stopPosition ) {
    ///     break; // No more possible matches among the remaining cached Features, so stop looking
    /// }
    /// ```
    pub fn cached_features_up_to_stop_position(&self, stop_position: i32) -> Vec<Feature> {
        let mut matching = Vec::new();
        for feature in &self.cache {
            if feature.start > stop_position {
                break;
            }
            matching.push(feature.clone());
        }
        matching
    }

    /// The cache's own contents, in order. Not part of the reference's API; the conformance suite
    /// compares it so that a divergence lands on the cache rather than on a later query.
    pub fn contents(&self) -> Vec<Feature> {
        self.cache.iter().cloned().collect()
    }

    pub fn cached_interval(&self) -> Option<&SimpleInterval> {
        self.cached_interval.as_ref()
    }
}

/// Where the features come from. A real one reads a file through a codec; the suite uses an
/// in-memory one, because what is measured here is the cache and not the parsing.
pub trait FeatureReader {
    /// Every feature overlapping the interval, in file order.
    fn query(&self, interval: &SimpleInterval) -> Vec<Feature>;
}

/// A `FeatureDataSource`'s query path: the cache, the lookahead, and the reader under both.
pub struct FeatureDataSource<R: FeatureReader> {
    reader: R,
    cache: FeatureCache,
    /// The constructor's `queryLookaheadBases`. GATK's default is
    /// `DEFAULT_QUERY_LOOKAHEAD_BASES = 1000`, and `IntervalWalker.initializeFeatures` sets it to
    /// **0** on purpose, because its query intervals are guaranteed not to overlap. There is no
    /// engine-wide command-line argument for it; a tool that exposes one names its own
    /// (`AnnotateIntervals`' `--feature-query-lookahead`, Funcotator's `--lookahead-cache-bp`).
    lookahead_bases: i32,
}

impl<R: FeatureReader> FeatureDataSource<R> {
    pub fn new(reader: R, lookahead_bases: i32) -> Self {
        FeatureDataSource {
            reader,
            cache: FeatureCache::new(),
            lookahead_bases,
        }
    }

    pub fn cache(&self) -> &FeatureCache {
        &self.cache
    }

    /// `FeatureDataSource.queryAndPrefetch`.
    ///
    /// Returns `Err` only where the reference throws: an improper trim. A query past the end of a
    /// contig is fine, because the reader is not aware of contig boundaries and the reference says
    /// so in a comment. The reference's first check, a `UserException` for a source that does not
    /// support random access, has no counterpart: a [`FeatureReader`] can always be queried.
    ///
    /// ```java
    /// if (queryCache.cacheHit(interval)) {
    ///     queryCache.trimToNewStartPosition(interval.getStart());
    /// }
    /// ```
    pub fn query_and_prefetch(
        &mut self,
        interval: &SimpleInterval,
    ) -> Result<Vec<Feature>, String> {
        if self.cache.cache_hit(interval) {
            self.cache.trim_to_new_start_position(interval.start)?;
        } else {
            self.refill_query_cache(interval);
        }
        Ok(self.cache.cached_features_up_to_stop_position(interval.end))
    }

    /// `FeatureDataSource.refillQueryCache`: query the reader over the interval extended by the
    /// lookahead, and keep everything it returns.
    fn refill_query_cache(&mut self, interval: &SimpleInterval) {
        let query_interval = SimpleInterval {
            contig: interval.contig.clone(),
            start: interval.start,
            // `Math.addExact`, which the reference chose so an overflow blows up rather than
            // travelling downstream as a negative end.
            end: interval
                .end
                .checked_add(self.lookahead_bases)
                .expect("interval end overflowed, which addExact refuses upstream"),
        };
        let features = self.reader.query(&query_interval);
        self.cache.fill(features, query_interval);
    }
}

/// `FeatureContext`: the interval a walker is at, and the window arithmetic around it.
pub struct FeatureContext {
    interval: Option<SimpleInterval>,
}

impl FeatureContext {
    pub fn new(interval: Option<SimpleInterval>) -> Self {
        FeatureContext { interval }
    }

    pub fn interval(&self) -> Option<&SimpleInterval> {
        self.interval.as_ref()
    }

    /// `FeatureContext.getQueryInterval`.
    ///
    /// Returns `Err` where the reference's `Utils.validateArg` throws: a negative expansion. The
    /// asymmetry is deliberate upstream: the leading edge is clamped at 1, and the trailing edge is
    /// **not** clamped to the contig, because a Feature query past the end of a contig is legal.
    ///
    /// ```java
    /// return Math.max(locus.getStart() - windowLeadingBases, 1);
    /// ```
    ///
    /// ```java
    /// //Note: queries past the end of contig are handled but overflow
    /// // blows up too late so we change it here to blow up early.
    /// return Math.addExact(locus.getEnd(), windowTrailingBases);//blow up on overflow
    /// ```
    ///
    /// The two argument checks come before the null-interval test, as in the reference, so a
    /// negative window is refused even when there is no interval to widen.
    pub fn query_interval(
        &self,
        leading: i32,
        trailing: i32,
    ) -> Result<Option<SimpleInterval>, String> {
        if leading < 0 {
            return Err("Window starts after the current interval".to_string());
        }
        if trailing < 0 {
            return Err("Window ends before the current interval".to_string());
        }
        let Some(interval) = &self.interval else {
            return Ok(None);
        };
        if leading == 0 && trailing == 0 {
            return Ok(Some(interval.clone()));
        }
        Ok(Some(SimpleInterval {
            contig: interval.contig.clone(),
            start: (interval.start - leading).max(1),
            end: interval
                .end
                .checked_add(trailing)
                .expect("window end overflowed, which addExact refuses upstream"),
        }))
    }

    /// `FeatureContext.getValues`: an empty list when there is no interval or no source, which is
    /// how a walker running without the optional Feature input gets an answer rather than an error.
    ///
    /// ```java
    /// if (featureManager == null || queryInterval == null || featureDescriptor == null) {
    ///     return Collections.emptyList();
    /// }
    /// ```
    ///
    /// `source` is `None` for the first and third of those; a `None` window is the second. The
    /// window is computed only once a source exists, so with no source a negative window is NOT
    /// refused here, while the reference computes `getQueryInterval` first and would refuse it.
    pub fn values<R: FeatureReader>(
        &self,
        source: Option<&mut FeatureDataSource<R>>,
        leading: i32,
        trailing: i32,
    ) -> Result<Vec<Feature>, String> {
        let Some(source) = source else {
            return Ok(Vec::new());
        };
        match self.query_interval(leading, trailing)? {
            None => Ok(Vec::new()),
            Some(interval) => source.query_and_prefetch(&interval),
        }
    }
}
