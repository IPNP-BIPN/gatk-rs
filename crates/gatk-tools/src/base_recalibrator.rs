//! Ported from `org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator` (GATK 4.6.2.0).
//!
//! The tool that writes the table [`crate::apply_bqsr`] reads, and the one that closes the BQSR
//! cycle. Its body is three calls: count every read that survives the filters, finalise the tables,
//! and write the report.
//!
//! # Its default read filters are seven, not one
//!
//! ```java
//! public static List<ReadFilter> getBQSRSpecificReadFilterList() {
//!     filters.add(ReadFilterLibrary.MAPPING_QUALITY_NOT_ZERO);
//!     filters.add(ReadFilterLibrary.MAPPING_QUALITY_AVAILABLE);
//!     filters.add(ReadFilterLibrary.MAPPED);
//!     filters.add(ReadFilterLibrary.NOT_SECONDARY_ALIGNMENT);
//!     filters.add(ReadFilterLibrary.NOT_DUPLICATE);
//!     filters.add(ReadFilterLibrary.PASSES_VENDOR_QUALITY_CHECK);
//! }
//! ```
//!
//! plus `WellformedReadFilter`. That is a fifth pattern of `getDefaultReadFilters` across the ported
//! tools, after taking the engine's default, replacing it with `ALLOW_ALL_READS`, extending it with
//! four more, and replacing it with a single non-Wellformed filter. It is also what decides which
//! reads are counted at all, so it is part of the table's contents and not a detail of the
//! traversal.
//!
//! # And the known sites must be indexed
//!
//! `--known-sites` is required and takes any feature file, VCF or BED. Without an index the
//! reference refuses before reading a single read, because the traversal queries it by interval.

use gatk_engine::base_recalibration_engine::{
    BaseRecalibrationEngine, EngineArguments, EngineError, ReferenceWindow,
};
use gatk_engine::interval::SimpleInterval;
use gatk_engine::qual_quantizer::{QualQuantizer, QuantizationInfo, MIN_USABLE_Q_SCORE};
use gatk_engine::reads::{ReadsDataSource, ReadsError};
use gatk_engine::recal_utils::output_recalibration_report;
use htsjdk_bam::record::BamRecord;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK BaseRecalibrator";

/// `RecalibrationArgumentCollection.QUANTIZING_LEVELS`.
pub const QUANTIZING_LEVELS: i32 = 16;

/// The reference as the counting pass asks for it: the bases `[start, end)` of the named contig,
/// zero-based and cut at its end, and the contig's length. A window, not the whole contig, is what
/// a run without BAQ holds: chr1's 249 million bases were held whole for reads that need a few
/// hundred at a time.
pub type ContigBases<'a> =
    dyn FnMut(&str, usize, usize) -> Result<(Vec<u8>, usize), BaseRecalibratorError> + 'a;

/// What the tool stops with.
#[derive(Debug, Clone, PartialEq)]
pub enum BaseRecalibratorError {
    Engine(EngineError),
    Reads(ReadsError),
}

impl BaseRecalibratorError {
    pub fn message(&self) -> String {
        match self {
            BaseRecalibratorError::Engine(error) => error.message(),
            BaseRecalibratorError::Reads(error) => format!("{error:?}"),
        }
    }
}

/// `BaseRecalibrator`: the recalibration table this run produced, as text.
///
/// `contig_bases` answers some of a reference contig by name (see [`ContigBases`]): a window around
/// the reads, or with BAQ the whole contig, whose model reads past the read. Each read is compared
/// against ITS contig: taking the reads' first one compared a chr2 read against chr1's bases.
///
/// The quantization is computed **after** `finalizeData`, so it depends on every read the traversal
/// kept, and the report is written from the finalised tables.
/// `contig_bases` is a CLOSURE and not a slice, because the order matters: the reference is opened
/// at startup and queried only once the traversal has produced reads, so a file whose reads cannot
/// be read fails on the reads and not on the reference. A row that hands a BAM the wrong index and
/// a reference the reads' contig is missing from answers the read failure, which is what the
/// reference answers.
#[allow(clippy::too_many_arguments)]
pub fn base_recalibrator(
    source: &ReadsDataSource,
    contig_bases: &mut ContigBases<'_>,
    known_sites: &KnownSites,
    arguments: &EngineArguments,
    quantizing_levels: i32,
    filter: &dyn Fn(&BamRecord) -> bool,
    intervals: &[SimpleInterval],
) -> Result<String, BaseRecalibratorError> {
    base_recalibrator_on(
        source,
        contig_bases,
        known_sites,
        arguments,
        quantizing_levels,
        filter,
        intervals,
        default_workers(),
        BATCH,
        WINDOW,
    )
}

/// [`base_recalibrator`] on `workers` counting threads in batches of `batch` reads, against
/// reference windows of `window` bases, which is how the test holds the threaded and windowed path
/// to the one-thread path on the golden's reads. Fewer than two workers, or BAQ on, counts on the
/// traversal's thread.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn base_recalibrator_on(
    source: &ReadsDataSource,
    contig_bases: &mut ContigBases<'_>,
    known_sites: &KnownSites,
    arguments: &EngineArguments,
    quantizing_levels: i32,
    filter: &dyn Fn(&BamRecord) -> bool,
    intervals: &[SimpleInterval],
    workers: usize,
    batch: usize,
    window: usize,
) -> Result<String, BaseRecalibratorError> {
    let header = source.header().clone();
    // The traversal is bounded like any other walker's. Hardcoding no intervals here counted the
    // reads a command line had excluded: `--intervals chr1:1-6000 --exclude-intervals chr1:1-500`
    // leaves seven of the corpus's eight reads, and the eighth's nine countable bases were nine
    // observations the recalibration table carried too many.
    // Counted as it is decoded, at most a few batches in flight: the BAM is never held (#1453
    // found 8 GiB at chr20 before the step stopped). A failure surfaces at the earliest read that
    // fails, which is the reference's own order.
    let records = crate::read_walker::traverse_iter(source, intervals, false, filter)
        .map_err(BaseRecalibratorError::Reads)?;
    base_recalibrator_counted(
        &header,
        records.map(|record| record.map_err(BaseRecalibratorError::Reads)),
        contig_bases,
        known_sites,
        arguments,
        quantizing_levels,
        workers,
        batch,
        window,
    )
}

/// Counting threads besides the traversal's: one per spare core, at most [`MAX_WORKERS`], and none
/// on a single core.
fn default_workers() -> usize {
    std::thread::available_parallelism()
        .map(|threads| threads.get().saturating_sub(1).min(MAX_WORKERS))
        .unwrap_or(0)
}

/// Past this many counting threads, decoding the BAM on one thread is the limit.
const MAX_WORKERS: usize = 6;
/// Reads per batch handed to a counting thread.
const BATCH: usize = 1024;

/// The table over `records`, counted as [`base_recalibrator_on`] says.
#[allow(clippy::too_many_arguments)]
fn base_recalibrator_counted<R: std::borrow::Borrow<BamRecord> + Send>(
    header: &htsjdk_bam::header::SamHeader,
    records: impl Iterator<Item = Result<R, BaseRecalibratorError>>,
    contig_bases: &mut ContigBases<'_>,
    known_sites: &KnownSites,
    arguments: &EngineArguments,
    quantizing_levels: i32,
    workers: usize,
    batch: usize,
    window: usize,
) -> Result<String, BaseRecalibratorError> {
    let fetched = Bases::new(arguments, window);
    let engine = if workers < 2 || arguments.enable_baq {
        count_here(
            header,
            records,
            contig_bases,
            fetched,
            known_sites,
            arguments,
        )?
    } else {
        count_on_threads(
            header,
            records,
            contig_bases,
            fetched,
            known_sites,
            arguments,
            workers,
            batch,
        )?
    };
    report(engine, arguments, quantizing_levels)
}

/// The table over reads something else has already chosen, which is what a pipeline hands
/// `BaseRecalibratorSparkFn`.
pub fn table_of_reads(
    header: &htsjdk_bam::header::SamHeader,
    records: &[BamRecord],
    contig_bases: &mut ContigBases<'_>,
    known_sites: &KnownSites,
    arguments: &EngineArguments,
    quantizing_levels: i32,
) -> Result<String, BaseRecalibratorError> {
    let engine = count_here(
        header,
        records.iter().map(Ok),
        contig_bases,
        Bases::new(arguments, WINDOW),
        known_sites,
        arguments,
    )?;
    report(engine, arguments, quantizing_levels)
}

/// The counting pass on this thread, read by read in the traversal's order.
fn count_here<R: std::borrow::Borrow<BamRecord>>(
    header: &htsjdk_bam::header::SamHeader,
    records: impl Iterator<Item = Result<R, BaseRecalibratorError>>,
    contig_bases: &mut ContigBases<'_>,
    mut fetched: Bases,
    known_sites: &KnownSites,
    arguments: &EngineArguments,
) -> Result<BaseRecalibrationEngine, BaseRecalibratorError> {
    let mut engine = BaseRecalibrationEngine::new(arguments.clone(), header)
        .map_err(BaseRecalibratorError::Engine)?;
    for record in records {
        let record = record?;
        let record = record.borrow();
        let bases = fetched.of(record, header, contig_bases)?;
        let overlapping = known_sites.overlapping(record, header);
        engine
            .process_read_in(record, header, &bases.view(), &overlapping)
            .map_err(BaseRecalibratorError::Engine)?;
    }
    Ok(engine)
}

/// The counting pass on `workers` threads, with the same tables as [`count_here`].
///
/// This thread decodes and hands batches of reads round robin to the workers; each worker counts
/// its reads into tables of its own, and their counts are added together at the end. That is
/// exact, so the split cannot show: the tables are indexed by key, so arrival order is nothing,
/// and without BAQ every increment is a whole error or none, whose sums do not round (the caller
/// keeps BAQ on one thread; `RecalDatum::absorb` says why). GATK counts on one thread and so
/// fails at the first read that fails; so does this: every failure carries its read's position,
/// and the earliest one is returned, whether a worker, the traversal or the reference raised it.
#[allow(clippy::too_many_arguments)]
fn count_on_threads<R: std::borrow::Borrow<BamRecord> + Send>(
    header: &htsjdk_bam::header::SamHeader,
    mut records: impl Iterator<Item = Result<R, BaseRecalibratorError>>,
    contig_bases: &mut ContigBases<'_>,
    mut fetched: Bases,
    known_sites: &KnownSites,
    arguments: &EngineArguments,
    workers: usize,
    batch: usize,
) -> Result<BaseRecalibrationEngine, BaseRecalibratorError> {
    use std::sync::{mpsc, Arc};
    // A read, its position in the traversal, and its contig's bases.
    type Item<R> = (usize, R, Arc<Window>);
    let mut engine = BaseRecalibrationEngine::new(arguments.clone(), header)
        .map_err(BaseRecalibratorError::Engine)?;

    let (counted, here) = std::thread::scope(|scope| {
        let lanes: Vec<_> = (0..workers)
            .map(|_| {
                let (to_worker, batches) = mpsc::sync_channel::<Vec<Item<R>>>(2);
                let worker = scope.spawn(move || {
                    // Built from the arguments and header the engine above was built from, so it
                    // cannot fail where that one did not.
                    let mut engine = BaseRecalibrationEngine::new(arguments.clone(), header)
                        .expect("the engine above was built from the same arguments");
                    let mut failed = None;
                    for items in batches {
                        // After a failure the lane is drained, not closed, so the traversal never
                        // blocks on it.
                        if failed.is_some() {
                            continue;
                        }
                        for (position, record, bases) in items {
                            let record = record.borrow();
                            let overlapping = known_sites.overlapping(record, header);
                            if let Err(error) =
                                engine.process_read_in(record, header, &bases.view(), &overlapping)
                            {
                                failed = Some((position, BaseRecalibratorError::Engine(error)));
                                break;
                            }
                        }
                    }
                    (engine.counted(), failed)
                });
                (to_worker, worker)
            })
            .collect();

        // The traversal: decode, fetch the reference the reads need, batch. A failure here ends it at
        // its position.
        let mut here = None;
        let mut items = Vec::with_capacity(batch);
        let mut lane = 0;
        for (position, record) in records.by_ref().enumerate() {
            let record = match record {
                Ok(record) => record,
                Err(error) => {
                    here = Some((position, error));
                    break;
                }
            };
            let bases = match fetched.of(record.borrow(), header, contig_bases) {
                Ok(bases) => bases,
                Err(error) => {
                    here = Some((position, error));
                    break;
                }
            };
            items.push((position, record, bases));
            if items.len() == batch {
                // A send fails only if the worker panicked, which the join below reports.
                let _ = lanes[lane]
                    .0
                    .send(std::mem::replace(&mut items, Vec::with_capacity(batch)));
                lane = (lane + 1) % workers;
            }
        }
        if !items.is_empty() {
            let _ = lanes[lane].0.send(items);
        }
        let counted: Vec<_> = lanes
            .into_iter()
            .map(|(to_worker, worker)| {
                drop(to_worker);
                worker.join().expect("a counting thread panicked")
            })
            .collect();
        (counted, here)
    });

    let mut earliest = here;
    let mut shares = Vec::with_capacity(counted.len());
    for (share, failed) in counted {
        if let Some((position, error)) = failed {
            if earliest.as_ref().is_none_or(|(first, _)| position < *first) {
                earliest = Some((position, error));
            }
        }
        shares.push(share);
    }
    if let Some((_, error)) = earliest {
        return Err(error);
    }
    for share in shares {
        engine
            .absorb(share)
            .map_err(BaseRecalibratorError::Engine)?;
    }
    Ok(engine)
}

/// Bases of the contig the reads are on, fetched once a read needs them, for the reason in
/// [`base_recalibrator`]'s doc comment: the whole contig once with BAQ on, otherwise a window of
/// [`WINDOW`] bases from the first read the last window does not hold, which a coordinate-sorted
/// traversal then fills read after read.
struct Bases {
    whole: bool,
    window: usize,
    current: Option<(String, std::sync::Arc<Window>)>,
}

/// Bases of one contig from the zero-based `start` on.
#[derive(Default)]
struct Window {
    start: usize,
    bases: Vec<u8>,
    contig_length: usize,
}

impl Window {
    fn view(&self) -> ReferenceWindow<'_> {
        ReferenceWindow {
            bases: &self.bases,
            start: self.start,
            contig_length: self.contig_length,
        }
    }

    /// Whether the window holds `[from, to)` once `to` is cut at the contig's end, which is where
    /// the engine cuts a read's span.
    fn holds(&self, from: usize, to: usize) -> bool {
        from >= self.start && to.min(self.contig_length) <= self.start + self.bases.len()
    }
}

/// How many bases a window fetches from the read that needs it.
const WINDOW: usize = 1 << 20;

impl Bases {
    fn new(arguments: &EngineArguments, window: usize) -> Self {
        Bases {
            whole: arguments.enable_baq,
            window,
            current: None,
        }
    }

    /// Bases that hold `record`'s aligned span, and so the span of the read the engine makes of
    /// it, since its transforms only clip.
    fn of(
        &mut self,
        record: &BamRecord,
        header: &htsjdk_bam::header::SamHeader,
        contig_bases: &mut ContigBases<'_>,
    ) -> Result<std::sync::Arc<Window>, BaseRecalibratorError> {
        let contig = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| header.sequences.get(index))
            .map(|sequence| sequence.name.as_str());
        let Some(contig) = contig else {
            return Ok(std::sync::Arc::default());
        };
        let from = gatk_engine::read_utils::start(record).max(1) as usize - 1;
        let to = usize::try_from(gatk_engine::read_utils::end(record))
            .unwrap_or(0)
            .max(from);
        if let Some((name, window)) = &self.current {
            if name == contig && (self.whole || window.holds(from, to)) {
                return Ok(window.clone());
            }
        }
        let (start, end) = if self.whole {
            (0, usize::MAX)
        } else {
            (from, to.max(from + self.window))
        };
        let (bases, contig_length) = contig_bases(contig, start, end)?;
        let window = std::sync::Arc::new(Window {
            start: start.min(contig_length),
            bases,
            contig_length,
        });
        self.current = Some((contig.to_string(), window.clone()));
        Ok(window)
    }
}

/// `finalizeData`, then the report. The quantization is computed **after** `finalizeData`, so it
/// depends on every read the traversal kept.
fn report(
    mut engine: BaseRecalibrationEngine,
    arguments: &EngineArguments,
    quantizing_levels: i32,
) -> Result<String, BaseRecalibratorError> {
    engine
        .finalize_data()
        .map_err(BaseRecalibratorError::Engine)?;
    let quantization = quantization_info(&engine, quantizing_levels);
    Ok(output_recalibration_report(
        &arguments.covariates,
        quantizing_levels,
        &quantization,
        &engine.tables,
        &engine.covariates,
    ))
}

/// The known sites, sorted by contig and start so that one read's are found by binary search.
/// Scanning them all for every read was 1.7 million dbSNP sites times 17.7 million reads on chr20.
///
/// Each site is kept as its two coordinates under its contig's name, not as a `SimpleInterval`:
/// a contig name per site was 1.4 million strings for dbSNP at chr20. Sites are added with
/// [`KnownSites::push`] and searched once [`KnownSites::sort`] has run.
#[derive(Debug, Default)]
pub struct KnownSites {
    by_contig: std::collections::HashMap<String, Contig>,
    sorted: bool,
}

#[derive(Debug)]
struct Contig {
    /// `(start, end)`, sorted by start, ties in the order given.
    sites: Vec<(i32, i32)>,
    /// The longest site, which bounds how far before a read's start an overlapping site can begin.
    longest: i64,
}

impl KnownSites {
    /// The sites of `known_sites`, sorted.
    pub fn new(known_sites: &[SimpleInterval]) -> Self {
        let mut sites = KnownSites::default();
        for site in known_sites {
            sites.push(&site.contig, site.start, site.end);
        }
        sites.sort();
        sites
    }

    /// One more site, in the order the files give them.
    pub fn push(&mut self, contig: &str, start: i32, end: i32) {
        let contig = match self.by_contig.get_mut(contig) {
            Some(sites) => sites,
            None => self.by_contig.entry(contig.to_string()).or_insert(Contig {
                sites: Vec::new(),
                longest: 0,
            }),
        };
        contig.sites.push((start, end));
        contig.longest = contig.longest.max(i64::from(end) - i64::from(start));
        self.sorted = false;
    }

    /// Sorts every contig's sites by start, stably, which the search needs.
    pub fn sort(&mut self) {
        for contig in self.by_contig.values_mut() {
            contig.sites.sort_by_key(|&(start, _)| start);
            contig.sites.shrink_to_fit();
        }
        self.sorted = true;
    }

    /// The known sites the engine is handed for one read: those on its contig overlapping
    /// `new SimpleInterval(read)`, its aligned span, which is the walker's `FeatureContext` and
    /// `JoinReadsWithVariants`' query alike. A site in a soft clip alone is not one of them, and a
    /// read whose span is no valid interval is handed none. The engine ORs them into one mask, so
    /// their order does not matter.
    fn overlapping(
        &self,
        record: &BamRecord,
        header: &htsjdk_bam::header::SamHeader,
    ) -> Vec<SimpleInterval> {
        let Some(contig) = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| header.sequences.get(index))
        else {
            return Vec::new();
        };
        let start = record.alignment_start;
        let end = record.alignment_end();
        if start < 1 || end < start {
            return Vec::new();
        }
        self.in_span(&contig.name, start, end)
    }

    fn in_span(&self, contig: &str, start: i32, end: i32) -> Vec<SimpleInterval> {
        assert!(
            self.sorted,
            "the known sites are searched before they are sorted"
        );
        let Some(sites) = self.by_contig.get(contig) else {
            return Vec::new();
        };
        let earliest = i64::from(start) - sites.longest;
        let first = sites
            .sites
            .partition_point(|&(site_start, _)| i64::from(site_start) < earliest);
        sites.sites[first..]
            .iter()
            .take_while(|&&(site_start, _)| site_start <= end)
            .filter(|&&(_, site_end)| site_end >= start)
            .map(|&(start, end)| SimpleInterval {
                contig: contig.to_string(),
                start,
                end,
            })
            .collect()
    }
}

/// `new QuantizationInfo(finalTables, QUANTIZING_LEVELS)`: the empirical quality histogram of the
/// quality score table, quantized.
///
/// It is built from the **finalised** tables, so the empirical qualities it bins are the ones the
/// rounding left behind.
fn quantization_info(engine: &BaseRecalibrationEngine, levels: i32) -> QuantizationInfo {
    let mut histogram = vec![0i64; 94];
    for (_, datum) in engine.tables.quality_score_table().all_leaves() {
        let empirical =
            gatk_engine::math_utils::fast_round(datum.borrow_mut().empirical_quality()) as usize;
        if empirical < histogram.len() {
            histogram[empirical] += datum.borrow().num_observations();
        }
    }
    match QualQuantizer::new(&histogram, levels, MIN_USABLE_Q_SCORE) {
        Ok(quantizer) => QuantizationInfo::new(quantizer.original_to_quantized_map, histogram),
        // `nLevels = 0` is the reference's own null dereference; the tool never passes it, because
        // its argument's default is sixteen.
        Err(_) => QuantizationInfo::new(vec![0; 94], histogram),
    }
}

#[cfg(test)]
mod known_sites_tests {
    use super::*;

    #[test]
    fn the_index_finds_what_a_scan_of_every_site_finds() {
        // A small LCG, so the test needs no dependency and fails the same way every time.
        let mut state = 7u64;
        let mut next = |bound: i32| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) % bound as u64) as i32
        };
        let sites: Vec<SimpleInterval> = (0..400)
            .map(|_| {
                let start = 1 + next(2000);
                SimpleInterval {
                    contig: if next(4) == 0 { "b" } else { "a" }.to_string(),
                    start,
                    end: start + if next(10) == 0 { next(300) } else { next(3) },
                }
            })
            .collect();
        let index = KnownSites::new(&sites);
        for _ in 0..2000 {
            let contig = if next(3) == 0 { "b" } else { "a" };
            let start = 1 + next(2100);
            let end = start + next(150);
            let mut scanned: Vec<_> = sites
                .iter()
                .filter(|site| site.contig == contig && site.start <= end && site.end >= start)
                .cloned()
                .collect();
            let mut found = index.in_span(contig, start, end);
            scanned.sort_by_key(|site| (site.start, site.end));
            found.sort_by_key(|site| (site.start, site.end));
            assert_eq!(found, scanned, "{contig}:{start}-{end}");
        }
        assert!(index.in_span("c", 1, 10).is_empty());
    }
}
