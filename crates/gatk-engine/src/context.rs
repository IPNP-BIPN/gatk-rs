//! Ported from `org.broadinstitute.hellbender.engine.ReferenceContext` and `ReadsContext`
//! (GATK 4.6.2.0).
//!
//! These are what a walker is handed at each position: the bases around it and the reads over it.
//! Every annotation and every caller reads its inputs through them, so a window off by one base is
//! an off-by-one in every number downstream, and it is invisible in the tool's own output.
//!
//! # The window is not the interval
//!
//! A `ReferenceContext` carries two spans. The **interval** is where the walker is; the **window**
//! is what it can see, the interval expanded by a leading and a trailing count and then cropped to
//! the contig. Three consequences, all measured in the conformance golden:
//!
//!  * near a contig edge the window is **smaller than asked for**, silently, and
//!    `numWindowLeadingBases` reports what was actually obtained rather than what was requested;
//!  * `getBases(leading, trailing)` expands from the **window**, not from the interval. After
//!    `setWindow(10, 10)`, asking for five more bases each side gives thirty bases, not twenty:
//!    the expansions compose rather than replace;
//!  * `getBase()` indexes the window by `interval.start - window.start`, so it is the base at the
//!    interval's start whatever the window did, including when cropping shortened the lead.
//!
//! And what comes back is not the FASTA's bytes: every query goes through the upper-casing and
//! IUPAC flattening documented in [`crate::reference`].

use crate::interval::SimpleInterval;
use crate::reference::{ReferenceError, ReferenceFileSource};

/// What a context refuses to answer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContextError {
    /// `setWindow`: a negative lead or trail. A `GATKException`, not a clamp to zero.
    NegativeWindow,
    /// `trimToContigLength`: the reference has no such contig. A `UserException`.
    ContigNotInReference(String),
    /// The window does not contain the interval, or a padding is negative.
    InvalidArgument,
    /// The query itself failed; see [`ReferenceError`].
    Reference(ReferenceError),
}

impl From<ReferenceError> for ContextError {
    fn from(error: ReferenceError) -> Self {
        ContextError::Reference(error)
    }
}

/// `ReferenceContext`.
///
/// The data source is passed in per call rather than held, because a query needs it mutably and a
/// walker holds one source for many contexts. The cached sequence, and the fact that changing the
/// window drops it, is the reference's and is kept.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceContext {
    interval: Option<SimpleInterval>,
    window: Option<SimpleInterval>,
    cached: Option<Vec<u8>>,
}

impl ReferenceContext {
    /// `new ReferenceContext(dataSource, interval, leading, trailing)`.
    ///
    /// `interval == None` is the "no known location" case, which answers every query with an
    /// empty array rather than an error.
    pub fn new(
        source: &ReferenceFileSource,
        interval: Option<SimpleInterval>,
        leading: i32,
        trailing: i32,
    ) -> Result<ReferenceContext, ContextError> {
        let mut context = ReferenceContext {
            interval,
            window: None,
            cached: None,
        };
        context.set_window(source, leading, trailing)?;
        Ok(context)
    }

    /// `new ReferenceContext(thatContext, interval)`: the *sizes* of the other context's window
    /// are carried over, not its coordinates.
    ///
    /// Since those sizes are what the other context actually obtained, a context built beside a
    /// contig edge propagates the cropped window rather than the requested one.
    pub fn with_interval(
        &self,
        source: &ReferenceFileSource,
        interval: Option<SimpleInterval>,
    ) -> Result<ReferenceContext, ContextError> {
        let leading = self.num_window_leading_bases();
        let trailing = self.num_window_trailing_bases();
        ReferenceContext::new(source, interval, leading, trailing)
    }

    /// `new ReferenceContext(dataSource, interval, window)`: an explicit window, still cropped.
    pub fn with_window(
        source: &ReferenceFileSource,
        interval: Option<SimpleInterval>,
        window: Option<SimpleInterval>,
    ) -> Result<ReferenceContext, ContextError> {
        if interval.is_none() && window.is_some() {
            return Err(ContextError::InvalidArgument);
        }
        let window = match (&interval, &window) {
            (Some(interval), Some(window)) => {
                // `window.contains(interval)`, which is same contig and covering span.
                if window.contig != interval.contig
                    || window.start > interval.start
                    || window.end < interval.end
                {
                    return Err(ContextError::InvalidArgument);
                }
                Some(SimpleInterval {
                    contig: interval.contig.clone(),
                    start: trim_to_contig_start(window.start),
                    end: trim_to_contig_length(source, &interval.contig, window.end)?,
                })
            }
            // The windowless case: the window *is* the interval.
            _ => interval.clone(),
        };
        Ok(ReferenceContext {
            interval,
            window,
            cached: None,
        })
    }

    pub fn interval(&self) -> Option<&SimpleInterval> {
        self.interval.as_ref()
    }

    pub fn window(&self) -> Option<&SimpleInterval> {
        self.window.as_ref()
    }

    /// `setWindow`: negative counts throw, and a `(0, 0)` window is the interval itself.
    pub fn set_window(
        &mut self,
        source: &ReferenceFileSource,
        leading: i32,
        trailing: i32,
    ) -> Result<(), ContextError> {
        if leading < 0 || trailing < 0 {
            return Err(ContextError::NegativeWindow);
        }
        self.window = match &self.interval {
            None => None,
            Some(interval) if leading == 0 && trailing == 0 => Some(interval.clone()),
            Some(interval) => Some(SimpleInterval {
                contig: interval.contig.clone(),
                start: trim_to_contig_start(interval.start - leading),
                end: trim_to_contig_length(source, &interval.contig, interval.end + trailing)?,
            }),
        };
        // Changing the window invalidates the cached query.
        self.cached = None;
        Ok(())
    }

    /// `numWindowLeadingBases`: what was obtained, which near a contig start is less than what was
    /// asked for.
    pub fn num_window_leading_bases(&self) -> i32 {
        match (&self.interval, &self.window) {
            (Some(interval), Some(window)) => interval.start - window.start,
            _ => 0,
        }
    }

    /// `numWindowTrailingBases`.
    pub fn num_window_trailing_bases(&self) -> i32 {
        match (&self.interval, &self.window) {
            (Some(interval), Some(window)) => window.end - interval.end,
            _ => 0,
        }
    }

    /// `getBases()`: the whole window, cached.
    pub fn bases(&mut self, source: &mut ReferenceFileSource) -> Result<Vec<u8>, ContextError> {
        let Some(window) = self.window.clone() else {
            return Ok(Vec::new());
        };
        if self.cached.is_none() {
            self.cached = Some(source.query(&window.contig, window.start, window.end)?);
        }
        Ok(self.cached.clone().unwrap_or_default())
    }

    /// `getBases(window)`: an arbitrary window, trimmed, and **not** cached.
    pub fn bases_of(
        &self,
        source: &mut ReferenceFileSource,
        window: &SimpleInterval,
    ) -> Result<Vec<u8>, ContextError> {
        let start = trim_to_contig_start(window.start);
        let end = trim_to_contig_length(source, &window.contig, window.end)?;
        Ok(source.query(&window.contig, start, end)?)
    }

    /// `getBases(leading, trailing)`: expands the **window**, not the interval.
    ///
    /// This is the one that surprises. On a context already built with `setWindow(10, 10)`, asking
    /// for five more each side returns the interval plus fifteen on each side, because the counts
    /// are added to the window's bounds rather than the interval's.
    pub fn bases_expanded(
        &self,
        source: &mut ReferenceFileSource,
        leading: i32,
        trailing: i32,
    ) -> Result<Vec<u8>, ContextError> {
        let Some(window) = &self.window else {
            return Ok(Vec::new());
        };
        let start = trim_to_contig_start(window.start - leading);
        let end = trim_to_contig_length(source, &window.contig, window.end + trailing)?;
        Ok(source.query(&window.contig, start, end)?)
    }

    /// `getForwardBases()`: from the interval's start to the end of the window.
    pub fn forward_bases(
        &mut self,
        source: &mut ReferenceFileSource,
    ) -> Result<Vec<u8>, ContextError> {
        let bases = self.bases(source)?;
        let (Some(interval), Some(window)) = (&self.interval, &self.window) else {
            return Ok(Vec::new());
        };
        let mid = (interval.start - window.start) as usize;
        if mid > bases.len() {
            // `String.substring` past the end throws rather than returning empty.
            return Err(ContextError::InvalidArgument);
        }
        Ok(bases[mid..].to_vec())
    }

    /// `getBase()`: the base at the interval's start, indexed inside the window.
    pub fn base(&mut self, source: &mut ReferenceFileSource) -> Result<u8, ContextError> {
        let bases = self.bases(source)?;
        let (Some(interval), Some(window)) = (&self.interval, &self.window) else {
            return Err(ContextError::InvalidArgument);
        };
        let index = (interval.start - window.start) as usize;
        bases
            .get(index)
            .copied()
            .ok_or(ContextError::InvalidArgument)
    }

    /// `getKmerAround(center, numBasesOnEachSide)`, or `None` where the reference returns null.
    ///
    /// The null case is a real answer rather than a failure: at a contig edge the window cannot
    /// expand to the requested size, and the reference declines to return a shorter kmer. A port
    /// that padded or truncated would give every caller a kmer the reference never produced.
    pub fn kmer_around(
        &self,
        source: &mut ReferenceFileSource,
        center: i32,
        bases_each_side: i32,
    ) -> Result<Option<Vec<u8>>, ContextError> {
        let Some(window) = &self.window else {
            return Err(ContextError::InvalidArgument);
        };
        if center < 1 || center < window.start || center > window.end || bases_each_side < 0 {
            return Err(ContextError::InvalidArgument);
        }
        let contig_length = contig_length(source, &window.contig)? as i32;
        // `expandWithinContig` then `trimIntervalToContig`, which returns null off the contig.
        let start = (center - bases_each_side).max(1);
        let end = (center + bases_each_side).min(contig_length);
        if start > contig_length || end < 1 {
            return Err(ContextError::InvalidArgument);
        }
        if end - start < 2 * bases_each_side {
            return Ok(None);
        }
        Ok(Some(source.query(&window.contig, start, end)?))
    }
}

/// `trimToContigStart`.
fn trim_to_contig_start(start: i32) -> i32 {
    start.max(1)
}

fn contig_length(source: &ReferenceFileSource, contig: &str) -> Result<usize, ContextError> {
    source
        .sequence_length(contig)
        .ok_or_else(|| ContextError::ContigNotInReference(contig.to_string()))
}

/// `trimToContigLength`, which throws for a contig the reference does not have rather than
/// returning the untrimmed end.
fn trim_to_contig_length(
    source: &ReferenceFileSource,
    contig: &str,
    end: i32,
) -> Result<i32, ContextError> {
    Ok(end.min(contig_length(source, contig)? as i32))
}

/// `ReadsContext`: the reads over an interval, optionally filtered.
///
/// Thin by design, and the thinness is the point: it is a query plus a filter, so the decisions
/// live in [`crate::reads`] and in the filter, not here. Both a missing source and a missing
/// interval answer with no reads rather than with every read.
pub struct ReadsContext<'a> {
    source: Option<&'a crate::reads::ReadsDataSource>,
    interval: Option<SimpleInterval>,
    filter: Option<ReadFilterFn<'a>>,
}

/// The `ReadFilter` a context may carry, which is `null` for most walkers.
pub type ReadFilterFn<'a> = Box<dyn Fn(&htsjdk_bam::record::BamRecord) -> bool + 'a>;

impl<'a> ReadsContext<'a> {
    pub fn new(
        source: Option<&'a crate::reads::ReadsDataSource>,
        interval: Option<SimpleInterval>,
        filter: Option<ReadFilterFn<'a>>,
    ) -> ReadsContext<'a> {
        ReadsContext {
            source,
            interval,
            filter,
        }
    }

    pub fn has_backing_data_source(&self) -> bool {
        self.source.is_some()
    }

    pub fn interval(&self) -> Option<&SimpleInterval> {
        self.interval.as_ref()
    }

    /// `iterator()`, and `iterator(interval)` when given one.
    pub fn reads(
        &self,
        interval: Option<&SimpleInterval>,
    ) -> Result<Vec<htsjdk_bam::record::BamRecord>, crate::reads::ReadsError> {
        let interval = interval.or(self.interval.as_ref());
        let (Some(source), Some(interval)) = (self.source, interval) else {
            return Ok(Vec::new());
        };
        let records = source.query(std::slice::from_ref(interval))?;
        Ok(match &self.filter {
            None => records,
            Some(filter) => records.into_iter().filter(|r| filter(r)).collect(),
        })
    }
}
