//! `GroundTruthScorer`: every read scored against the reference it aligns to, and the report that
//! summarises them.
//!
//! The flow-based scoring is not ported. What is ported is the report the tool builds out of the
//! scores: the accumulators, the four table shapes, the phred each row carries, and the bins the
//! deviation and the base are folded into.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.groundtruth.GroundTruthScorer` in
//! GATK 4.6.2.0.

use crate::series_stats::SeriesStats;

/// The bounds the report is allocated with.
pub const QUAL_VALUE_MAX: usize = 60;
pub const HMER_VALUE_MAX: usize = 100;
/// `FlowBasedRead.DEFAULT_FLOW_ORDER.length() - 1`, the flow order being `TGCA`.
pub const BASE_VALUE_MAX: usize = 3;
pub const DEFAULT_FLOW_ORDER: &str = "TGCA";

/// `NORMALIZED_SCORE_THRESHOLD_DEFAULT`, which is NEGATIVE: the scores it bounds are.
pub const NORMALIZED_SCORE_THRESHOLD_DEFAULT: f64 = -0.1;
pub const DEFAULT_RATIO_THRESHOLD: f64 = 0.003;

/// The percentile columns the report carries when none are asked for.
pub const DEFAULT_QUALITY_PERCENTILES: &str = "10,25,50,75,90";

/// One cell of the report: how many observations were true and how many false.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Accumulator {
    pub false_count: u64,
    pub true_count: u64,
}

impl Accumulator {
    pub fn add(&mut self, value: bool) {
        if value {
            self.true_count += 1;
        } else {
            self.false_count += 1;
        }
    }

    pub fn count(&self) -> u64 {
        self.false_count + self.true_count
    }

    /// The FALSE rate, which is zero for an empty cell rather than a NaN.
    pub fn false_rate(&self) -> f64 {
        if self.count() == 0 {
            0.0
        } else {
            self.false_count as f64 / self.count() as f64
        }
    }
}

/// `deviationToBin`: `0,-1,1,-2,2...` become `0,1,2,3,4...`.
pub fn deviation_to_bin(deviation: i32) -> usize {
    if deviation >= 0 {
        (deviation * 2) as usize
    } else {
        ((-deviation * 2) - 1) as usize
    }
}

/// `binToDeviation`, which is NOT the inverse of `deviation_to_bin` as a number: it is a STRING,
/// and a positive deviation is written with a leading `+` while zero and the negatives are not.
pub fn bin_to_deviation(bin: usize) -> String {
    if bin == 0 {
        "0".to_string()
    } else if bin.is_multiple_of(2) {
        format!("+{}", bin / 2)
    } else {
        format!("{}", -((bin as i64 + 1) / 2))
    }
}

/// `binToBase`, which indexes the flow order.
pub fn bin_to_base(bin: usize) -> char {
    DEFAULT_FLOW_ORDER.chars().nth(bin).expect("a base")
}

/// The phred a rate becomes in the one-level table.
///
/// A rate of zero over a NON-EMPTY cell is replaced by the probability threshold before the
/// logarithm, so a cell that saw observations and no errors reports the threshold's phred rather
/// than zero. A rate of zero over an EMPTY cell is left alone and reports zero.
pub fn phred(rate: f64, count: u64, probability_threshold: f64) -> i64 {
    let effective = if rate == 0.0 && count != 0 && probability_threshold != 0.0 {
        probability_threshold
    } else {
        rate
    };
    if effective != 0.0 {
        (-10.0 * jmath::math::log10(effective)).ceil() as i64
    } else {
        0
    }
}

/// The four table names the report carries, built from the column names.
pub fn one_level_table_name(name: &str) -> String {
    format!("{name}Report")
}

pub fn two_level_table_name(first: &str, second: &str) -> String {
    format!("{first}_{second}Report")
}

pub fn four_level_table_name(first: &str, second: &str, third: &str, fourth: &str) -> String {
    format!("{first}_{second}_{third}_{fourth}_Report")
}

/// Whether a row is written, given where it sits and what it holds.
///
/// The ORIGIN is always written: `omit_zeros` skips an empty row only when at least one of its
/// indices is non-zero, so the first row of every table survives however empty it is.
pub fn keeps_row(omit_zeros: bool, indices: &[usize], count: u64) -> bool {
    if !omit_zeros {
        return true;
    }
    let at_origin = indices.iter().all(|index| *index == 0);
    at_origin || count != 0
}

/// How many rows the four-level table has when its zeros are kept.
///
/// It is the product of the four allocated dimensions and not of what was observed, which is why
/// the option that omits the zeros is not really optional.
pub fn four_level_row_count() -> usize {
    (QUAL_VALUE_MAX + 1)
        * (HMER_VALUE_MAX + 1)
        * deviation_to_bin(HMER_VALUE_MAX as i32 + 1)
        * (BASE_VALUE_MAX + 1)
}

/// The percentile table's columns, which are fixed except for the percentiles themselves.
pub const PERCENTILE_TABLE_NAME: &str = "PhredBinAccumulator";
pub const PERCENTILE_FIXED_COLUMNS: [&str; 7] =
    ["flow", "count", "min", "max", "mean", "median", "std"];

/// The percentile table's column names for a given `--quality-percentiles`.
pub fn percentile_columns(quality_percentiles: &str) -> Vec<String> {
    let mut columns: Vec<String> = PERCENTILE_FIXED_COLUMNS
        .iter()
        .map(|name| name.to_string())
        .collect();
    for percentile in quality_percentiles.split(',') {
        columns.push(format!("p{percentile}"));
    }
    columns
}

/// One flow's percentile row, which is a `SeriesStats` fed in PHRED space.
///
/// `addProb` takes a probability and stores `-10 * log10(p)`, so the series holds phreds and the
/// percentiles are over those rather than over the probabilities they came from.
#[derive(Debug, Clone, Default)]
pub struct PercentileReport {
    pub stats: SeriesStats,
}

impl PercentileReport {
    pub fn add_probability(&mut self, probability: f64) {
        self.stats.add(-10.0 * jmath::math::log10(probability));
    }

    /// One row of the percentile table, in the column order above.
    pub fn row(&self, index: usize, quality_percentiles: &str) -> Vec<f64> {
        let mut values = vec![
            index as f64,
            self.stats.count() as f64,
            self.stats.min(),
            self.stats.max(),
            self.stats.mean(),
            self.stats.median(),
            self.stats.std(),
        ];
        for percentile in quality_percentiles.split(',') {
            values.push(
                self.stats
                    .percentile(percentile.parse().expect("a percentile")),
            );
        }
        values
    }
}

/// The arguments that decide what is scored and what is written.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub use_softclipped_bases: bool,
    pub normalized_score_threshold: f64,
    pub add_mean_call: bool,
    pub no_output: bool,
    pub omit_zeros_from_report: bool,
    pub quality_percentiles: String,
    pub exclude_zero_flows: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            use_softclipped_bases: false,
            normalized_score_threshold: NORMALIZED_SCORE_THRESHOLD_DEFAULT,
            add_mean_call: false,
            no_output: false,
            omit_zeros_from_report: false,
            quality_percentiles: DEFAULT_QUALITY_PERCENTILES.to_string(),
            exclude_zero_flows: false,
        }
    }
}

/// Whether a read's normalized score keeps it.
///
/// The comparison is a strict less-than against a NEGATIVE default, so a read is dropped when its
/// score falls below the threshold rather than above it.
pub fn keeps_read(normalized_score: f64, arguments: &Arguments) -> bool {
    normalized_score >= arguments.normalized_score_threshold
}

/// The two columns `--add-mean-call` appends to the CSV, in the order it appends them.
///
/// The probabilities come FIRST and the mean call second, which is the other way round from the
/// argument's own name.
pub const MEAN_CALL_COLUMNS: [&str; 2] = ["ReadProbs", "ReadMeanCall"];

// ================================================================================================
// The walk: every read's flow score, error probabilities and report entries.
// ================================================================================================

use crate::flow_based_read::FlowRead;
use gatk_engine::gatk_report::{Report, Sorting, Table, Value};
use htsjdk_bam::cigar::{CigarElement as BamCigarElement, Op};

/// `isSoftClipped`: a soft clip at EXACTLY one end. A read clipped at both ends is scored with its
/// clips in place, as if it were not clipped at all.
pub fn is_soft_clipped(unmapped: bool, cigar: &[BamCigarElement]) -> bool {
    if unmapped {
        return false;
    }
    let (Some(first), Some(last)) = (cigar.first(), cigar.last()) else {
        return false;
    };
    (first.op == Op::S) != (last.op == Op::S)
}

/// `GenomePriorDB`: one row per base, in file order, a later row for the same base replacing the
/// earlier one's counts.
///
/// The first COUNT is read twice: `prior[0]` and `prior[1]` are both column one, and every later
/// `prior[i]` is column `i`, so column `HMER_VALUE_MAX + 1` is never read at all.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GenomePrior {
    rows: Vec<(u8, Vec<i64>)>,
}

impl GenomePrior {
    /// `new GenomePriorDB(path)` over the file's text, which is read by opencsv: a line is split
    /// on commas, and a number that does not parse is `NumberFormatException` with its text.
    pub fn parse(text: &str) -> Result<GenomePrior, (String, String)> {
        let mut rows: Vec<(u8, Vec<i64>)> = Vec::new();
        for line in text.lines() {
            let columns: Vec<&str> = line.split(',').collect();
            let base = columns[0].as_bytes().first().copied().ok_or_else(|| {
                (
                    "java.lang.ArrayIndexOutOfBoundsException".to_string(),
                    "Index 0 out of bounds for length 0".to_string(),
                )
            })?;
            let mut prior = vec![0i64; HMER_VALUE_MAX + 1];
            for (i, slot) in prior.iter_mut().enumerate() {
                let at = if i == 0 { 1 } else { i };
                let cell = columns.get(at).ok_or_else(|| {
                    (
                        "java.lang.ArrayIndexOutOfBoundsException".to_string(),
                        format!("Index {at} out of bounds for length {}", columns.len()),
                    )
                })?;
                *slot = cell.parse::<i64>().map_err(|_| {
                    (
                        "java.lang.NumberFormatException".to_string(),
                        format!("For input string: \"{cell}\""),
                    )
                })?;
            }
            match rows.iter_mut().find(|(known, _)| *known == base) {
                Some(row) => row.1 = prior,
                None => rows.push((base, prior)),
            }
        }
        Ok(GenomePrior { rows })
    }

    fn for_base(&self, base: u8) -> Option<&[i64]> {
        self.rows
            .iter()
            .find(|(known, _)| *known == base)
            .map(|(_, prior)| prior.as_slice())
    }
}

/// `computeErrorProb`: per flow, one less the normalized probability of the called length, scaled
/// by the genome prior of the flow's base when there is one.
///
/// Each value is also added to its flow's percentile series, unless the flow called zero and zero
/// calls are excluded.
pub fn error_probabilities(
    read: &FlowRead,
    prior: Option<&GenomePrior>,
    percentiles: Option<&mut Vec<PercentileReport>>,
    exclude_zero_flows: bool,
) -> Vec<f64> {
    let max_hmer = read.max_hmer;
    let mut column = vec![0.0f64; max_hmer as usize + 1];
    let mut result = vec![0.0f64; read.key.len()];
    let mut percentiles = percentiles;
    for i in 0..read.key.len() {
        let mut sum = 0.0;
        for (j, cell) in column.iter_mut().enumerate() {
            *cell = read.prob(i, j as i32);
            sum += *cell;
        }
        if sum != 0.0 {
            for cell in column.iter_mut() {
                *cell /= sum;
            }
        }
        let called = read.key[i].min(max_hmer) as usize;
        match prior {
            Some(prior) => {
                let mut sum = 0.0;
                if let Some(counts) = prior.for_base(read.flow_order[i]) {
                    for (j, cell) in column.iter_mut().enumerate() {
                        *cell *= counts[j] as f64;
                        sum += *cell;
                    }
                }
                result[i] = if sum != 0.0 {
                    1.0 - column[(column.len() - 1).min(called)] / sum
                } else {
                    1.0 - column[called]
                };
            }
            None => result[i] = 1.0 - column[called],
        }
        if let Some(reports) = percentiles.as_deref_mut() {
            if read.key[i] != 0 || !exclude_zero_flows {
                while reports.len() < i + 1 {
                    reports.push(PercentileReport::default());
                }
                reports[i].add_probability(result[i]);
            }
        }
    }
    result
}

/// `computeLowestQBaseTP`: per flow, the `tp` of the lowest-quality base in the first half of its
/// hmer, the first such base winning a tie. The `tp` array is the record's tag as it came, which a
/// hard clip does not shorten, while the qualities are the clipped read's.
pub fn lowest_q_base_tp(key: &[i32], tp: &[i8], quals: &[u8]) -> Result<Vec<i8>, String> {
    let out_of_bounds =
        |index: usize, length: usize| format!("Index {index} out of bounds for length {length}");
    let mut result = vec![0i8; key.len()];
    let mut seq = 0usize;
    for (i, &hmer) in key.iter().enumerate() {
        if hmer == 0 {
            result[i] = 0;
            continue;
        }
        result[i] = *tp.get(seq).ok_or_else(|| out_of_bounds(seq, tp.len()))?;
        let mut lowest = *quals
            .get(seq)
            .ok_or_else(|| out_of_bounds(seq, quals.len()))? as i8;
        let scan = ((hmer + 1) / 2) as usize;
        for j in 1..scan {
            let q = *quals
                .get(seq + j)
                .ok_or_else(|| out_of_bounds(seq + j, quals.len()))? as i8;
            if q < lowest {
                result[i] = *tp
                    .get(seq + j)
                    .ok_or_else(|| out_of_bounds(seq + j, tp.len()))?;
                lowest = q;
            }
        }
        seq += hmer as usize;
    }
    Ok(result)
}

/// `FlowBasedReadUtils.CycleSkipStatus`, whose priority is its order here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CycleSkipStatus {
    NS,
    PCS,
    CS,
}

impl CycleSkipStatus {
    pub fn name(self) -> &'static str {
        match self {
            CycleSkipStatus::NS => "NS",
            CycleSkipStatus::PCS => "PCS",
            CycleSkipStatus::CS => "CS",
        }
    }
}

/// `FlowBasedReadUtils.getCycleSkipStatus`: every `M` element whose bases differ from the
/// reference's, compared in flow space. Only `M` is looked at, and the subarrays are clamped at
/// either array's end as `ArrayUtils.subarray` clamps them.
pub fn cycle_skip_status(
    bases: &[u8],
    cigar: &[BamCigarElement],
    reference: &[u8],
    flow_order: &str,
) -> CycleSkipStatus {
    let sub = |array: &[u8], from: usize, length: usize| -> Vec<u8> {
        let start = from.min(array.len());
        let end = (from + length).min(array.len());
        array[start..end].to_vec()
    };
    let mut status = CycleSkipStatus::NS;
    let mut read_offset = 0usize;
    let mut ref_offset = 0usize;
    for element in cigar {
        let length = element.length as usize;
        if element.op == Op::M {
            let read_bases = sub(bases, read_offset, length);
            let ref_bases = sub(reference, ref_offset, length);
            if read_bases != ref_bases {
                let alt = crate::flow_based_read::base_array_to_key(&read_bases, flow_order);
                let refk = crate::flow_based_read::base_array_to_key(&ref_bases, flow_order);
                let value = match (alt, refk) {
                    (Some(alt), Some(refk)) => {
                        if alt.len() != refk.len() {
                            CycleSkipStatus::CS
                        } else if refk.iter().zip(&alt).any(|(r, a)| (*r == 0) != (*a == 0)) {
                            CycleSkipStatus::PCS
                        } else {
                            CycleSkipStatus::NS
                        }
                    }
                    _ => CycleSkipStatus::CS,
                };
                if value > status {
                    status = value;
                }
            }
        }
        if element.op.consumes_read_bases() {
            read_offset += length;
        }
        if element.op.consumes_reference_bases() {
            ref_offset += length;
        }
        if status == CycleSkipStatus::CS {
            break;
        }
    }
    status
}

/// `collectReadProbs`: every cell of the matrix flow by flow, and each flow's mean call, which
/// divides by the SUM OF THE LENGTHS rather than of the probabilities.
pub fn read_probs_and_mean_call(read: &FlowRead) -> (Vec<f64>, Vec<f64>) {
    let mut probs = Vec::with_capacity(read.key.len() * (read.max_hmer as usize + 1));
    let mut mean = Vec::with_capacity(read.key.len());
    for flow in 0..read.key.len() {
        let mut mc = 0.0;
        let mut mc_sum = 0i32;
        for hmer in 0..=read.max_hmer {
            let p = read.prob(flow, hmer);
            probs.push(p);
            mc += p * hmer as f64;
            mc_sum += hmer;
        }
        mean.push(mc / mc_sum as f64);
    }
    (probs, mean)
}

/// `new DecimalFormat("0.0#####")`: at least one fraction digit and at most six, half to even.
pub fn error_format(value: f64) -> String {
    let text = gatk_annotation::decimal_format::DecimalFormat::new(6).format(value);
    if value.is_finite() && !text.contains('.') {
        format!("{text}.0")
    } else {
        text
    }
}

/// `Precision.round(value, scale)`: the value's `Double.toString`, rounded half up.
fn precision_round(value: f64, scale: usize) -> f64 {
    if !value.is_finite() {
        return value;
    }
    gatk_engine::java_format::format_decimals(value, scale)
        .parse()
        .unwrap_or(value)
}

/// `BASE_VALUE_MAX + 1`, the base bins.
const BASE_BINS: usize = BASE_VALUE_MAX + 1;
/// `deviationToBin(HMER_VALUE_MAX + 1)`, the deviation bins.
const DEVIATION_BINS: usize = 2 * (HMER_VALUE_MAX + 1);

/// The four-level `BooleanAccumulator` tree, kept sparse: the level that is never reported (qual,
/// hmer, deviation) is not kept at all.
#[derive(Debug, Clone, Default)]
pub struct QualReport {
    qual: Vec<Accumulator>,
    qual_hmer: std::collections::BTreeMap<(usize, usize), Accumulator>,
    four: std::collections::BTreeMap<(usize, usize, usize, usize), Accumulator>,
}

/// A bin clipped into range, as `BooleanAccumulator.add` clips one it warns about.
fn clip_bin(bin: i64, bins: usize) -> usize {
    bin.max(0).min(bins as i64 - 1) as usize
}

impl QualReport {
    pub fn new() -> QualReport {
        QualReport {
            qual: vec![Accumulator::default(); QUAL_VALUE_MAX + 1],
            ..Default::default()
        }
    }

    fn add(&mut self, same: bool, qual: usize, hmer: i64, deviation: i64, base: i64) {
        self.qual[qual].add(same);
        let hmer = clip_bin(hmer, HMER_VALUE_MAX + 1);
        self.qual_hmer.entry((qual, hmer)).or_default().add(same);
        let deviation = clip_bin(deviation, DEVIATION_BINS);
        let base = clip_bin(base, BASE_BINS);
        self.four
            .entry((qual, hmer, deviation, base))
            .or_default()
            .add(same);
    }

    /// `addToQualReport`: the read's key against the reference's, flow by flow, when the two have
    /// the same length. `flow_order` is the read's four-base cycle.
    pub fn add_read(
        &mut self,
        read: &FlowRead,
        reference: &[u8],
        flow_order: &str,
        reverse: bool,
        error: &[f64],
    ) {
        let Some(haplotype) = crate::flow_pairhmm_align_reads_to_haplotypes::FlowHaplotype::new(
            reference, flow_order,
        ) else {
            return;
        };
        if read.key.len() != haplotype.key.len() {
            return;
        }
        let order = flow_order.as_bytes();
        for flow in 0..read.key.len() {
            let prob = precision_round(error[flow], QUAL_VALUE_MAX / 10 + 1);
            let qual = (-10.0 * jmath::math::log10(prob)).ceil() as i32;
            let deviation = read.key[flow] - haplotype.key[flow];
            if qual >= 0 && (qual as usize) < self.qual.len() {
                let base = order[flow % order.len()];
                let base = if reverse { complement(base) } else { base };
                let bin = DEFAULT_FLOW_ORDER
                    .bytes()
                    .position(|known| known == base)
                    .map(|at| at as i64)
                    .unwrap_or(-1);
                self.add(
                    deviation == 0,
                    qual as usize,
                    read.key[flow] as i64,
                    deviation_to_bin(deviation) as i64,
                    bin,
                );
            }
        }
    }

    fn qual_table(&self, omit_zeros: bool) -> Table {
        let mut table = Table::new("qualReport", "error rate per qual", Sorting::DoNotSort);
        for (name, format) in [
            ("qual", "%d"),
            ("count", "%d"),
            ("error", "%f"),
            ("phred", "%d"),
        ] {
            table.add_column(name, format);
        }
        let mut row = 0;
        for (i, cell) in self.qual.iter().enumerate() {
            if omit_zeros && i != 0 && cell.count() == 0 {
                continue;
            }
            let rate = cell.false_rate();
            let key = row.to_string();
            table.set(&key, "qual", Value::Int(i as i64));
            table.set(&key, "count", Value::Int(cell.count() as i64));
            table.set(&key, "error", Value::Double(rate));
            table.set(
                &key,
                "phred",
                Value::Int(phred(rate, cell.count(), DEFAULT_RATIO_THRESHOLD)),
            );
            row += 1;
        }
        table
    }

    fn qual_hmer_table(&self, omit_zeros: bool) -> Table {
        let mut table = Table::new(
            "qual_hmerReport",
            "error rate per qual by hmer",
            Sorting::DoNotSort,
        );
        for (name, format) in [
            ("qual", "%d"),
            ("hmer", "%d"),
            ("count", "%d"),
            ("error", "%f"),
        ] {
            table.add_column(name, format);
        }
        let empty = Accumulator::default();
        let mut row = 0;
        for i in 0..self.qual.len() {
            for j in 0..=HMER_VALUE_MAX {
                let cell = self.qual_hmer.get(&(i, j)).unwrap_or(&empty);
                if omit_zeros && (i != 0 || j != 0) && cell.count() == 0 {
                    continue;
                }
                let key = row.to_string();
                table.set(&key, "qual", Value::Int(i as i64));
                table.set(&key, "hmer", Value::Int(j as i64));
                table.set(&key, "count", Value::Int(cell.count() as i64));
                table.set(&key, "error", Value::Double(cell.false_rate()));
                row += 1;
            }
        }
        table
    }

    fn four_level_table(&self, omit_zeros: bool) -> Table {
        let mut table = Table::new(
            "qual_hmer_deviation_base_Report",
            "error rate per qual by hmer and deviation",
            Sorting::DoNotSort,
        );
        for (name, format) in [
            ("qual", "%d"),
            ("hmer", "%d"),
            ("deviation", "%s"),
            ("base", "%s"),
            ("count", "%d"),
        ] {
            table.add_column(name, format);
        }
        let mut row = 0;
        let mut emit =
            |table: &mut Table, (i, j, k, m): (usize, usize, usize, usize), count: u64| {
                let key = row.to_string();
                table.set(&key, "qual", Value::Int(i as i64));
                table.set(&key, "hmer", Value::Int(j as i64));
                table.set(&key, "deviation", Value::Str(bin_to_deviation(k)));
                table.set(&key, "base", Value::Str(bin_to_base(m).to_string()));
                table.set(&key, "count", Value::Int(count as i64));
                row += 1;
            };
        let count_at = |at: (usize, usize, usize, usize)| {
            self.four.get(&at).map(Accumulator::count).unwrap_or(0)
        };
        if omit_zeros {
            // The origin always, then every cell that saw something, in the nested loops' order.
            let origin = (0, 0, 0, 0);
            emit(&mut table, origin, count_at(origin));
            for (&at, cell) in &self.four {
                if at != origin && cell.count() != 0 {
                    emit(&mut table, at, cell.count());
                }
            }
        } else {
            for i in 0..self.qual.len() {
                for j in 0..=HMER_VALUE_MAX {
                    for k in 0..DEVIATION_BINS {
                        for m in 0..BASE_BINS {
                            emit(&mut table, (i, j, k, m), count_at((i, j, k, m)));
                        }
                    }
                }
            }
        }
        table
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

/// `PercentileReport.newReportTable`.
fn percentile_table(reports: &[PercentileReport], quality_percentiles: &str) -> Table {
    let mut table = Table::new(
        PERCENTILE_TABLE_NAME,
        PERCENTILE_TABLE_NAME,
        Sorting::DoNotSort,
    );
    let columns = percentile_columns(quality_percentiles);
    for (index, name) in columns.iter().enumerate() {
        table.add_column(name, if index < 2 { "%d" } else { "%f" });
    }
    for (index, report) in reports.iter().enumerate() {
        let key = index.to_string();
        let values = report.row(index, quality_percentiles);
        for (column, (name, value)) in columns.iter().zip(values).enumerate() {
            let cell = if column < 2 {
                Value::Int(value as i64)
            } else {
                Value::Double(value)
            };
            table.set(&key, name, cell);
        }
    }
    table
}

/// The report file: its four tables in the order `GATKReport`'s `TreeMap` writes them, which is
/// the order of their names.
pub fn report_text(
    qual: &QualReport,
    percentiles: &[PercentileReport],
    omit_zeros: bool,
    quality_percentiles: &str,
) -> String {
    let mut tables = vec![
        qual.qual_table(omit_zeros),
        qual.qual_hmer_table(omit_zeros),
        qual.four_level_table(omit_zeros),
        percentile_table(percentiles, quality_percentiles),
    ];
    tables.sort_by(|a, b| a.name.cmp(&b.name));
    let mut report = Report::new();
    for table in tables {
        report.add_table(table);
    }
    report.write()
}
