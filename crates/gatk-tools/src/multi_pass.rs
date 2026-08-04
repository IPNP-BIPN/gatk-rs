//! The multi-pass walkers, ported from `org.broadinstitute.hellbender.engine`
//! (GATK 4.6.2.0): `MultiplePassVariantWalker`, `TwoPassVariantWalker` and
//! `MultiplePassReadWalker`.
//!
//! Every other walker makes one traversal. These make several, and what they do *between* passes
//! is where they disagree with each other. All three disagreements are observable, and none is
//! visible from a single-pass run.
//!
//! # One filter for the whole run, or one per pass
//!
//! ```java
//! // MultiplePassVariantWalker.traverse
//! final CountingVariantFilter countingVariantFilter = makeVariantFilter();
//! final CountingReadFilter readFilter = makeReadFilter();
//! for (int n = 0; n < numberOfPasses(); n++) { ... }
//! logger.info(countingVariantFilter.getSummaryLine());
//! ```
//!
//! ```java
//! // MultiplePassReadWalker.forEachRead
//! if (passCount > 1) {
//!     countedFilter = makeReadFilter();
//!     resetReadsDataSource();
//! }
//! ```
//!
//! The variant walker builds both filters once, before the loop, and the counts therefore
//! **accumulate**: the same file filtered by the same predicate reports one drop after a one-pass
//! run, two after a two-pass run and three after a three-pass run, because it is one counter that
//! nobody reset. The read walker builds a **new** filter at the top of every pass after the first,
//! so each pass reports only its own drops. The two classes are one directory apart.
//!
//! The port makes the difference structural rather than remembered: [`traverse_multiple_pass`]
//! takes one `&mut` counter for the whole run, and [`MultiplePassReadWalker`] owns the vector of
//! counters it built.
//!
//! # `afterNthPass` runs after **every** pass, including the last
//!
//! It is inside the loop, after the traversal, not after the loop. A two-pass run calls it twice.
//! Its name and its javadoc ("Process the data collected during the first pass. This method is
//! called between the two traversals") both suggest otherwise.
//!
//! # `TwoPassVariantWalker` has no `afterSecondPass`, and does not say so
//!
//! ```java
//! final protected void afterNthPass(final int n) {
//!     if (n == 0) { afterFirstPass(); }
//!     else if (n > 1) { throw new GATKException.ShouldNeverReachHereException(...); }
//! }
//! ```
//!
//! `n == 1` matches neither branch, so the call after the second pass falls through and does
//! nothing at all. The guard reads as "anything but zero is an error", and the one value it
//! actually receives besides zero is the one it silently ignores. [`two_pass_after_route`] answers
//! all three cases so the gap is a value rather than an omission.
//!
//! # Zero passes is a legal traversal
//!
//! `numberOfPasses()` is consulted by `n < numberOfPasses()`, so a walker returning zero traverses
//! nothing, calls `afterNthPass` never, and still logs both filter summaries. Nothing rejects it.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_readfilter::counting::Counting;
use htsjdk_bam::record::BamRecord;

/// `CountingVariantFilter`, as far as a multi-pass traversal observes it: a predicate and the
/// number of records it rejected.
///
/// The read side already has [`gatk_readfilter::counting::Counting`], whose composition and summary
/// text are ported there. This is the variant side's counterpart, and only the counting is here:
/// the filter **library** for variants is its own slice, and the default one
/// (`VariantFilterLibrary.ALLOW_ALL_VARIANTS`) never rejects anything, so a port that stopped at
/// the default could not tell a reused counter from a rebuilt one.
pub struct CountingVariantFilter<V> {
    predicate: Box<dyn Fn(&V) -> bool>,
    filtered_count: u64,
}

impl<V> CountingVariantFilter<V> {
    /// `new CountingVariantFilter(filter)`.
    pub fn new(predicate: impl Fn(&V) -> bool + 'static) -> Self {
        Self {
            predicate: Box::new(predicate),
            filtered_count: 0,
        }
    }

    /// `getFilteredCount()`.
    pub fn filtered_count(&self) -> u64 {
        self.filtered_count
    }

    /// `resetFilteredCount()`. Nothing in the multi-pass traversal calls it, which is the point.
    pub fn reset_filtered_count(&mut self) {
        self.filtered_count = 0;
    }

    /// `test(variant)`: the count goes up on **rejection**, so it is a count of what did not reach
    /// `apply`.
    pub fn test(&mut self, variant: &V) -> bool {
        let keep = (self.predicate)(variant);
        if !keep {
            self.filtered_count += 1;
        }
        keep
    }
}

/// `MultiplePassVariantWalker.traverse`.
///
/// `filter` is taken by reference for the whole run because the reference builds it once before the
/// loop. Passing a fresh one per pass would be the read walker's behaviour, and the counts are how
/// the two are told apart.
///
/// `after_nth_pass` is called after every pass, the last one included.
pub fn traverse_multiple_pass<V>(
    passes: usize,
    variants: &[V],
    filter: &mut CountingVariantFilter<V>,
    nth_pass_apply: &mut dyn FnMut(&V, usize),
    after_nth_pass: &mut dyn FnMut(usize),
) {
    for pass in 0..passes {
        for variant in variants {
            // `.filter(variantFilter).forEach(...)`: the filter runs first, and a rejected variant
            // reaches neither `apply` nor the progress meter.
            if filter.test(variant) {
                nth_pass_apply(variant, pass);
            }
        }
        after_nth_pass(pass);
    }
}

/// Where `TwoPassVariantWalker.nthPassApply` sends pass `n`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PassApply {
    /// `n == 0`.
    FirstPassApply,
    /// `n == 1`.
    SecondPassApply,
    /// Anything else: `GATKException.ShouldNeverReachHereException`.
    Refused,
}

/// Where `TwoPassVariantWalker.afterNthPass` sends pass `n`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AfterPass {
    /// `n == 0`.
    AfterFirstPass,
    /// `n == 1`: **nothing happens**. Not an error, not a callback. See the module note.
    Nothing,
    /// `n > 1`: `GATKException.ShouldNeverReachHereException`.
    Refused,
}

/// `TwoPassVariantWalker.nthPassApply`, which is `final` and routes by `n`.
pub fn two_pass_apply_route(n: usize) -> PassApply {
    match n {
        0 => PassApply::FirstPassApply,
        1 => PassApply::SecondPassApply,
        _ => PassApply::Refused,
    }
}

/// `TwoPassVariantWalker.afterNthPass`, whose middle case is a hole rather than a branch.
pub fn two_pass_after_route(n: usize) -> AfterPass {
    if n == 0 {
        AfterPass::AfterFirstPass
    } else if n > 1 {
        AfterPass::Refused
    } else {
        AfterPass::Nothing
    }
}

/// `MultiplePassReadWalker`.
///
/// The state that matters is `passCount`, which starts at **one** and is incremented at the end of
/// `forEachRead`, so the test `passCount > 1` is false exactly on the first call. The filter for
/// that first call is the one `traverse()` built.
pub struct MultiplePassReadWalker {
    /// `passCount`, one-based as in the reference.
    pass_count: u32,
    /// Every filter this run built, in order: the one from `traverse()` and then one per pass after
    /// the first. Their counts are the observable that separates this class from
    /// `MultiplePassVariantWalker`.
    filters: Vec<Counting>,
}

impl MultiplePassReadWalker {
    /// `traverse()`, which builds the first filter and then hands control to `traverseReads()`.
    ///
    /// The filter is built here even by a tool whose `traverseReads` calls `forEachRead` zero
    /// times, so a zero-pass run still reports one filter built and no reads.
    pub fn new(make_read_filter: impl FnOnce() -> Counting) -> Self {
        Self {
            pass_count: 1,
            filters: vec![make_read_filter()],
        }
    }

    /// `forEachRead(readHandler)`: one traversal of the input reads.
    ///
    /// ```java
    /// if (passCount > 1) {
    ///     countedFilter = makeReadFilter();
    ///     resetReadsDataSource();
    /// }
    /// ```
    ///
    /// Both halves of that branch are here, but only one of them is observable in this port.
    /// `resetReadsDataSource` exists because htsjdk's reader is consumed by iterating it, so a
    /// second pass over the same source would see nothing without it; the ported
    /// [`crate::read_walker::traverse_with_bounds_mut`]
    /// re-queries the source per call and is already repeatable, so the reset has nothing to undo.
    /// What is observable is the **new filter**, and that is reproduced exactly: a fresh counter
    /// per pass after the first, each reporting only its own drops.
    ///
    /// The bounds are passed in again rather than remembered, which matches
    /// `resetReadsDataSource` re-initialising the source **with** its traversal parameters: pass
    /// two of a `-L`-bounded run sees the same reads as pass one, not the whole file.
    pub fn for_each_read(
        &mut self,
        source: &ReadsDataSource,
        intervals: &[SimpleInterval],
        traverse_unmapped: bool,
        make_read_filter: impl FnOnce() -> Counting,
        handler: &mut dyn FnMut(&BamRecord),
    ) -> Result<(), ReadsError> {
        if self.pass_count > 1 {
            self.filters.push(make_read_filter());
        }
        // The filter in force is the last one built, which on the first pass is `traverse()`'s.
        let filter_index = self.filters.len() - 1;

        let records = {
            let filter = &mut self.filters[filter_index];
            crate::read_walker::traverse_with_bounds_mut(
                source,
                intervals,
                traverse_unmapped,
                &mut |read| filter.test(read),
            )?
        };
        for read in &records {
            handler(read);
        }

        self.pass_count += 1;
        Ok(())
    }

    /// `passCount`, after however many passes have run.
    pub fn pass_count(&self) -> u32 {
        self.pass_count
    }

    /// The filters this run built, in order.
    pub fn filters(&self) -> &[Counting] {
        &self.filters
    }
}
