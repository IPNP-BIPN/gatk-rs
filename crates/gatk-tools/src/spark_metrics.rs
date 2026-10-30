//! The metrics the Spark pipelines compute themselves.
//!
//! `MeanQualityByCycleSpark`, `QualityScoreDistributionSpark`,
//! `CollectBaseDistributionByCycleSpark` and `CollectQualityYieldMetricsSpark` are not their Picard
//! namesakes on another engine: each aggregates its own counts in `tools.spark.pipelines.metrics`
//! or `hellbender.metrics`, under a metric class of GATK's, so they are ported from there and read
//! one partition's reads in order.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram, MetricBean, MetricsFile, Value};

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const REVERSE: u16 = 0x10;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const FAILS_VENDOR_QUALITY: u16 = 0x200;
const SUPPLEMENTARY: u16 = 0x800;

/// `MetricsReadFilter`, whose secondary and supplementary checks no argument turns off.
pub fn metrics_read_filter(record: &BamRecord, pf_only: bool, aligned_only: bool) -> bool {
    (!pf_only || record.flags & FAILS_VENDOR_QUALITY == 0)
        && (!aligned_only || record.flags & UNMAPPED == 0)
        && record.flags & (SECONDARY | SUPPLEMENTARY) == 0
}

/// `ReadUtils.getOriginalBaseQualities`: an absent or empty `OQ` is none.
fn original_qualities(record: &BamRecord) -> Option<Vec<u8>> {
    match record.tags.get(Tag::new(b"OQ")) {
        Some(TagValue::Str(text)) if !text.is_empty() => {
            Some(htsjdk_bam::fastq::fastq_to_phred(text))
        }
        _ => None,
    }
}

fn second_of_pair(record: &BamRecord) -> bool {
    record.flags & PAIRED != 0 && record.flags & SECOND_OF_PAIR != 0
}

fn with_headers(headers: &[String]) -> MetricsFile {
    let mut file = MetricsFile::new();
    for header in headers {
        file.add_header(header);
    }
    file
}

/// `MeanQualityByCycleSpark.HistogramGenerator`.
#[derive(Default)]
struct MeanQualities {
    first_totals: Vec<f64>,
    first_counts: Vec<u64>,
    second_totals: Vec<f64>,
    second_counts: Vec<u64>,
    longest: usize,
}

impl MeanQualities {
    fn add(&mut self, qualities: &[u8], reverse: bool, second: bool) {
        let length = qualities.len();
        if length + 1 > self.longest {
            self.longest = length + 1;
            self.first_totals.resize(self.longest, 0.0);
            self.first_counts.resize(self.longest, 0);
            self.second_totals.resize(self.longest, 0.0);
            self.second_counts.resize(self.longest, 0);
        }
        for (index, quality) in qualities.iter().enumerate() {
            let cycle = if reverse { length - index } else { index + 1 };
            if second {
                self.second_totals[cycle] += f64::from(*quality);
                self.second_counts[cycle] += 1;
            } else {
                self.first_totals[cycle] += f64::from(*quality);
                self.first_counts[cycle] += 1;
            }
        }
    }

    /// `getMeanQualityHistogram`. A first-read cycle counts only if its total is positive, which
    /// is not the test the second read's cycles get, and a second read's cycle is added to
    /// whatever the first read's left at that key.
    fn histogram(&self, label: &str) -> Histogram {
        let mut bins: std::collections::BTreeMap<usize, f64> = std::collections::BTreeMap::new();
        let mut first_length = 0;
        for cycle in 0..self.first_totals.len() {
            if self.first_totals[cycle] > 0.0 {
                *bins.entry(cycle).or_default() +=
                    self.first_totals[cycle] / self.first_counts[cycle] as f64;
                first_length = cycle;
            }
        }
        for index in 0..self.second_totals.len() {
            if self.second_counts[index] > 0 {
                *bins.entry(first_length + index).or_default() +=
                    self.second_totals[index] / self.second_counts[index] as f64;
            }
        }
        Histogram {
            bin_label: "CYCLE".to_string(),
            value_label: label.to_string(),
            key_class: "java.lang.Integer".to_string(),
            bins: bins
                .into_iter()
                .map(|(cycle, mean)| (cycle.to_string(), mean))
                .collect(),
        }
    }
}

/// `MeanQualityByCycleSpark.calculateMeanQualityByCycle`, over reads already filtered.
pub fn mean_quality_by_cycle<'a>(
    records: impl Iterator<Item = &'a BamRecord>,
    headers: &[String],
) -> String {
    let (mut qualities, mut originals) = (MeanQualities::default(), MeanQualities::default());
    for record in records {
        let reverse = record.flags & REVERSE != 0;
        let second = second_of_pair(record);
        qualities.add(&record.base_qualities, reverse, second);
        if let Some(original) = original_qualities(record) {
            originals.add(&original, reverse, second);
        }
    }
    let mut file = with_headers(headers);
    file.histograms.push(qualities.histogram("MEAN_QUALITY"));
    if originals.longest != 0 {
        file.histograms
            .push(originals.histogram("MEAN_ORIGINAL_QUALITY"));
    }
    file.write()
}

/// `QualityScoreDistributionSpark`, over reads already filtered.
///
/// `SequenceUtil.isNoCall` is `N`, `n` or `.`, and `OQ` is indexed by the qualities' length.
pub fn quality_score_distribution<'a>(
    records: impl Iterator<Item = &'a BamRecord>,
    include_no_calls: bool,
    headers: &[String],
) -> String {
    let mut counts = [0u64; 128];
    let mut original_counts = [0u64; 128];
    for record in records {
        let original = original_qualities(record);
        for (index, quality) in record.base_qualities.iter().enumerate() {
            let base = record.read_bases.get(index).copied().unwrap_or(b'N');
            if include_no_calls || !matches!(base, b'N' | b'n' | b'.') {
                counts[usize::from(*quality)] += 1;
                if let Some(original) = &original {
                    original_counts[usize::from(original[index])] += 1;
                }
            }
        }
    }
    let histogram = |label: &str, counts: &[u64; 128]| Histogram {
        bin_label: "QUALITY".to_string(),
        value_label: label.to_string(),
        key_class: "java.lang.Byte".to_string(),
        bins: counts
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .map(|(quality, count)| (quality.to_string(), *count as f64))
            .collect(),
    };
    let mut file = with_headers(headers);
    file.histograms.push(histogram("COUNT_OF_Q", &counts));
    let originals = histogram("COUNT_OF_OQ", &original_counts);
    if !originals.bins.is_empty() {
        file.histograms.push(originals);
    }
    file.write()
}

/// `BaseDistributionByCycleMetrics`.
struct BaseDistributionMetric {
    read_end: i64,
    cycle: i64,
    percents: [f64; 5],
}

impl MetricBean for BaseDistributionMetric {
    fn class_name(&self) -> &str {
        "org.broadinstitute.hellbender.metrics.analysis.BaseDistributionByCycleMetrics"
    }

    fn columns(&self) -> &[&'static str] {
        &[
            "READ_END", "CYCLE", "PCT_A", "PCT_C", "PCT_G", "PCT_T", "PCT_N",
        ]
    }

    fn values(&self) -> Vec<Value> {
        let mut values = vec![Value::Long(self.read_end), Value::Long(self.cycle)];
        values.extend(self.percents.iter().map(|percent| Value::Double(*percent)));
        values
    }
}

/// `CollectBaseDistributionByCycleSpark.HistogramGenerator`, over reads already filtered.
pub fn base_distribution_by_cycle<'a>(
    records: impl Iterator<Item = &'a BamRecord>,
    headers: &[String],
) -> String {
    let (mut first, mut second): (Vec<[u64; 5]>, Vec<[u64; 5]>) = (Vec::new(), Vec::new());
    let mut seen_second = false;
    for record in records {
        let length = record.read_bases.len();
        if length + 1 > first.len() {
            first.resize(length + 1, [0; 5]);
            second.resize(length + 1, [0; 5]);
        }
        let reverse = record.flags & REVERSE != 0;
        let is_second = second_of_pair(record);
        seen_second |= is_second;
        let totals = if is_second { &mut second } else { &mut first };
        for (index, base) in record.read_bases.iter().enumerate() {
            let cycle = if reverse { length - index } else { index + 1 };
            let bin = match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            };
            totals[cycle][bin] += 1;
        }
    }
    let metric = |read_end: i64, cycle: usize, offset: usize, totals: &[u64; 5]| {
        let count: u64 = totals.iter().sum();
        BaseDistributionMetric {
            read_end,
            cycle: (cycle + offset) as i64,
            percents: totals.map(|total| 100.0 * total as f64 / count as f64),
        }
    };
    let mut file = with_headers(headers);
    let mut last_first_cycle = 0;
    for (cycle, totals) in first.iter().enumerate() {
        if totals.iter().sum::<u64>() != 0 {
            file.add_metric(&metric(1, cycle, 0, totals));
            last_first_cycle = cycle;
        }
    }
    if seen_second {
        for (cycle, totals) in second.iter().enumerate() {
            if totals.iter().sum::<u64>() != 0 {
                file.add_metric(&metric(2, cycle, last_first_cycle, totals));
            }
        }
    }
    file.write()
}

/// `hellbender.metrics.QualityYieldMetrics`, finished.
#[derive(Default)]
struct QualityYield {
    total_reads: i64,
    pf_reads: i64,
    read_length: i64,
    total_bases: i64,
    pf_bases: i64,
    q20_bases: i64,
    pf_q20_bases: i64,
    q30_bases: i64,
    pf_q30_bases: i64,
    q20_equivalent_yield: i64,
    pf_q20_equivalent_yield: i64,
}

impl MetricBean for QualityYield {
    fn class_name(&self) -> &str {
        "org.broadinstitute.hellbender.metrics.QualityYieldMetrics"
    }

    fn columns(&self) -> &[&'static str] {
        &[
            "TOTAL_READS",
            "PF_READS",
            "READ_LENGTH",
            "TOTAL_BASES",
            "PF_BASES",
            "Q20_BASES",
            "PF_Q20_BASES",
            "Q30_BASES",
            "PF_Q30_BASES",
            "Q20_EQUIVALENT_YIELD",
            "PF_Q20_EQUIVALENT_YIELD",
        ]
    }

    fn values(&self) -> Vec<Value> {
        [
            self.total_reads,
            self.pf_reads,
            self.read_length,
            self.total_bases,
            self.pf_bases,
            self.q20_bases,
            self.pf_q20_bases,
            self.q30_bases,
            self.pf_q30_bases,
            self.q20_equivalent_yield,
            self.pf_q20_equivalent_yield,
        ]
        .into_iter()
        .map(Value::Long)
        .collect()
    }
}

/// `QualityYieldMetricsCollectorSpark`, with no filter of its own beyond the tool's.
pub fn quality_yield<'a>(
    records: impl Iterator<Item = &'a BamRecord>,
    use_original: bool,
    headers: &[String],
) -> String {
    let mut yields = QualityYield::default();
    for record in records {
        yields.total_reads += 1;
        let length = record.read_bases.len() as i64;
        let passes = record.flags & FAILS_VENDOR_QUALITY == 0;
        if passes {
            yields.pf_reads += 1;
            yields.pf_bases += length;
        }
        yields.total_bases += length;
        let original = use_original.then(|| original_qualities(record)).flatten();
        for quality in original.as_deref().unwrap_or(&record.base_qualities) {
            let value = i64::from(*quality);
            yields.q20_equivalent_yield += value;
            yields.q20_bases += i64::from(value >= 20);
            yields.q30_bases += i64::from(value >= 30);
            if passes {
                yields.pf_q20_equivalent_yield += value;
                yields.pf_q20_bases += i64::from(value >= 20);
                yields.pf_q30_bases += i64::from(value >= 30);
            }
        }
    }
    if yields.total_reads != 0 {
        yields.read_length = yields.total_bases / yields.total_reads;
    }
    yields.q20_equivalent_yield /= 20;
    yields.pf_q20_equivalent_yield /= 20;
    let mut file = with_headers(headers);
    file.add_metric(&yields);
    file.write()
}
