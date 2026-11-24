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
    BaseRecalibrationEngine, EngineArguments, EngineError,
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
/// `contig_bases` answers a whole reference contig by name, because the counting pass needs the
/// read's span for its comparison and a wider window for BAQ. Each read is compared against ITS
/// contig: taking the reads' first one compared a chr2 read against chr1's bases.
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
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, BaseRecalibratorError>,
    known_sites: &[SimpleInterval],
    arguments: &EngineArguments,
    quantizing_levels: i32,
    filter: &dyn Fn(&BamRecord) -> bool,
    intervals: &[SimpleInterval],
) -> Result<String, BaseRecalibratorError> {
    let header = source.header().clone();
    // The traversal is bounded like any other walker's. Hardcoding no intervals here counted the
    // reads a command line had excluded: `--intervals chr1:1-6000 --exclude-intervals chr1:1-500`
    // leaves seven of the corpus's eight reads, and the eighth's nine countable bases were nine
    // observations the recalibration table carried too many.
    // One read at a time, counted as it is decoded: the BAM is never held (#1453 found 8 GiB at
    // chr20 before the step stopped). A decode failure now surfaces at the read that fails, after
    // the reads before it were counted, which is the reference's own order.
    let records = crate::read_walker::traverse_iter(source, intervals, false, filter)
        .map_err(BaseRecalibratorError::Reads)?;
    table_of_stream(
        &header,
        records.map(|record| record.map_err(BaseRecalibratorError::Reads)),
        contig_bases,
        known_sites,
        arguments,
        quantizing_levels,
    )
}

/// The table over reads something else has already chosen, which is what a pipeline hands
/// `BaseRecalibratorSparkFn`.
pub fn table_of_reads(
    header: &htsjdk_bam::header::SamHeader,
    records: &[BamRecord],
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, BaseRecalibratorError>,
    known_sites: &[SimpleInterval],
    arguments: &EngineArguments,
    quantizing_levels: i32,
) -> Result<String, BaseRecalibratorError> {
    table_of_stream(
        header,
        records.iter().map(Ok),
        contig_bases,
        known_sites,
        arguments,
        quantizing_levels,
    )
}

fn table_of_stream<R: std::borrow::Borrow<BamRecord>>(
    header: &htsjdk_bam::header::SamHeader,
    records: impl Iterator<Item = Result<R, BaseRecalibratorError>>,
    contig_bases: &mut dyn FnMut(&str) -> Result<Vec<u8>, BaseRecalibratorError>,
    known_sites: &[SimpleInterval],
    arguments: &EngineArguments,
    quantizing_levels: i32,
) -> Result<String, BaseRecalibratorError> {
    let sites = KnownSites::new(known_sites);
    // The bases come once the traversal has produced a read on the contig, for the reason in the
    // doc comment above, and once per contig.
    let mut fetched: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();

    let mut engine = BaseRecalibrationEngine::new(arguments.clone(), header)
        .map_err(BaseRecalibratorError::Engine)?;
    for record in records {
        let record = record?;
        let record = record.borrow();
        let contig = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| header.sequences.get(index))
            .map(|sequence| sequence.name.clone());
        let bases: &[u8] = match contig {
            Some(contig) => {
                if !fetched.contains_key(&contig) {
                    let bases = contig_bases(&contig)?;
                    fetched.insert(contig.clone(), bases);
                }
                &fetched[&contig]
            }
            None => &[],
        };
        let overlapping = sites.overlapping(record, header);
        engine
            .process_read(record, header, bases, &overlapping)
            .map_err(BaseRecalibratorError::Engine)?;
    }
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
struct KnownSites<'a> {
    by_contig: std::collections::HashMap<&'a str, Contig<'a>>,
}

struct Contig<'a> {
    /// Sorted by start, ties in the order given.
    sites: Vec<&'a SimpleInterval>,
    /// The longest site, which bounds how far before a read's start an overlapping site can begin.
    longest: i64,
}

impl<'a> KnownSites<'a> {
    fn new(known_sites: &'a [SimpleInterval]) -> Self {
        let mut by_contig: std::collections::HashMap<&str, Contig> =
            std::collections::HashMap::new();
        for site in known_sites {
            let contig = by_contig.entry(site.contig.as_str()).or_insert(Contig {
                sites: Vec::new(),
                longest: 0,
            });
            contig.sites.push(site);
            contig.longest = contig
                .longest
                .max(i64::from(site.end) - i64::from(site.start));
        }
        for contig in by_contig.values_mut() {
            contig.sites.sort_by_key(|site| site.start);
        }
        KnownSites { by_contig }
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
        let Some(sites) = self.by_contig.get(contig) else {
            return Vec::new();
        };
        let earliest = i64::from(start) - sites.longest;
        let first = sites
            .sites
            .partition_point(|site| i64::from(site.start) < earliest);
        sites.sites[first..]
            .iter()
            .take_while(|site| site.start <= end)
            .filter(|site| site.end >= start)
            .map(|site| (*site).clone())
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
