//! `GroundTruthScorer`: every read scored against the reference it aligns to, and the report that
//! summarises them.
//!
//! The report the tool builds out of the scores (the accumulators, the four table shapes, the
//! phred each row carries, and the bins the deviation and the base are folded into), and what the
//! scoring needs besides `computeLikelihoodLocal`, which is `FlowFeatureMapper`'s: the cycle-skip
//! status, the error probabilities with and without a genome prior, `LowestQBaseTP`, the read
//! probabilities, and the two number formats the CSV is written in. The walk itself is the
//! `gatk-cli` runner's.
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

/// The CSV's columns, in the order `onTraversalStart` fixes them.
pub const CSV_FIELD_ORDER_BASIC: [&str; 15] = [
    "ReadName",
    "ReadKey",
    "ReadIsReversed",
    "ReadMQ",
    "ReadRQ",
    "GroundTruthKey",
    "ReadSequence",
    "Score",
    "NormalizedScore",
    "ErrorProbability",
    "ReadKeyLength",
    "GroundTruthKeyLength",
    "CycleSkipStatus",
    "Cigar",
    "LowestQBaseTP",
];

/// `FlowBasedReadUtils.CycleSkipStatus`, in the order of its priority.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CycleSkipStatus {
    /// Not a cycle skip.
    NS,
    /// Possibly one: the keys are as long, and a flow is zero in one and not the other.
    PCS,
    /// One: the keys differ in length.
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

/// `baseArrayToKey`'s period guard, which a base outside the flow order trips.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PeriodGuard {
    pub bases: Vec<u8>,
}

/// `FlowBasedReadUtils.getCycleSkipStatus`: each `M` element whose bases differ from the reference
/// under it, its two keys compared in the read's flow order, the worst answer kept.
///
/// `ArrayUtils.subarray` clamps rather than throws, so an element running past either array is
/// compared over what is there.
pub fn cycle_skip_status(
    bases: &[u8],
    cigar: &[(char, usize)],
    reference: &[u8],
    flow_order: &str,
) -> Result<CycleSkipStatus, PeriodGuard> {
    let sub = |array: &[u8], from: usize, length: usize| -> Vec<u8> {
        let start = from.min(array.len());
        let end = (from + length).min(array.len());
        array[start..end].to_vec()
    };
    let mut status = CycleSkipStatus::NS;
    let (mut read_offset, mut ref_offset) = (0usize, 0usize);
    for &(op, length) in cigar {
        if op == 'M' {
            let element = sub(bases, read_offset, length);
            let under = sub(reference, ref_offset, length);
            if element != under {
                let key = |bases: &[u8]| {
                    crate::flow_based_read::base_array_to_key(bases, flow_order).ok_or(
                        PeriodGuard {
                            bases: bases.to_vec(),
                        },
                    )
                };
                let alt_key = key(&element)?;
                let ref_key = key(&under)?;
                let mut value = if ref_key.len() != alt_key.len() {
                    CycleSkipStatus::CS
                } else {
                    CycleSkipStatus::NS
                };
                if value == CycleSkipStatus::NS
                    && ref_key
                        .iter()
                        .zip(&alt_key)
                        .any(|(r, a)| (*r == 0) ^ (*a == 0))
                {
                    value = CycleSkipStatus::PCS;
                }
                status = status.max(value);
            }
        }
        if matches!(op, 'M' | 'I' | 'S' | '=' | 'X') {
            read_offset += length;
        }
        if matches!(op, 'M' | 'D' | 'N' | '=' | 'X') {
            ref_offset += length;
        }
        if status == CycleSkipStatus::CS {
            break;
        }
    }
    Ok(status)
}

/// `GenomePriorDB`: per base, a hundred and one hmer frequencies.
///
/// The row's first frequency is read TWICE, into the first two slots, and every later slot takes
/// the column of its own index, so the row's last column is never read.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GenomePrior {
    pub rows: Vec<(u8, Vec<i64>)>,
}

/// What reading the prior refused, as the JVM raises it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PriorError {
    pub class: &'static str,
    pub message: String,
}

impl GenomePrior {
    /// The CSV, one row per base. A quoted field is refused as the port's: opencsv would unquote it.
    pub fn parse(text: &str) -> Result<GenomePrior, PriorError> {
        let mut rows: Vec<(u8, Vec<i64>)> = Vec::new();
        for line in text.lines() {
            if line.contains('"') {
                return Err(PriorError {
                    class: "",
                    message: "a quoted genome-prior field".to_string(),
                });
            }
            let fields: Vec<&str> = line.split(',').collect();
            let base = *fields[0].as_bytes().first().ok_or_else(|| PriorError {
                class: "java.lang.ArrayIndexOutOfBoundsException",
                message: "Index 0 out of bounds for length 0".to_string(),
            })?;
            let mut prior = vec![0i64; HMER_VALUE_MAX + 1];
            for (i, slot) in prior.iter_mut().enumerate() {
                let column = if i == 0 { 1 } else { i };
                let text = fields.get(column).ok_or_else(|| PriorError {
                    class: "java.lang.ArrayIndexOutOfBoundsException",
                    message: format!("Index {column} out of bounds for length {}", fields.len()),
                })?;
                *slot = text.parse::<i64>().map_err(|_| PriorError {
                    class: "java.lang.NumberFormatException",
                    message: format!("For input string: \"{text}\""),
                })?;
            }
            // `LinkedHashMap.put`: a repeated base keeps its slot and takes the new row.
            match rows.iter_mut().find(|(other, _)| *other == base) {
                Some(row) => row.1 = prior,
                None => rows.push((base, prior)),
            }
        }
        Ok(GenomePrior { rows })
    }

    fn prior_for_base(&self, base: u8) -> Option<&[i64]> {
        self.rows
            .iter()
            .find(|(other, _)| *other == base)
            .map(|(_, prior)| prior.as_slice())
    }
}

/// `computeErrorProb`: per flow, one less the normalized probability of the flow's own call, with
/// the column rescaled by the genome prior of the flow's base when there is one.
///
/// Each value also feeds that flow's percentile row, when the report is kept, unless the call is
/// zero and `--exclude-zero-flows` leaves those out.
pub fn error_probabilities(
    read: &crate::flow_based_read::FlowRead,
    prior: Option<&GenomePrior>,
    mut percentiles: Option<&mut Vec<PercentileReport>>,
    exclude_zero_flows: bool,
) -> Vec<f64> {
    let max_hmer = read.max_hmer;
    let mut column = vec![0.0f64; max_hmer as usize + 1];
    let mut result = vec![0.0f64; read.key.len()];
    for (i, slot) in result.iter_mut().enumerate() {
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
        let call = read.key[i].min(max_hmer) as usize;
        *slot = match prior {
            Some(prior) => {
                let mut sum = 0.0;
                if let Some(frequencies) = prior.prior_for_base(read.flow_order[i]) {
                    for (j, cell) in column.iter_mut().enumerate() {
                        *cell *= frequencies[j] as f64;
                        sum += *cell;
                    }
                }
                if sum != 0.0 {
                    1.0 - column[(column.len() - 1).min(call)] / sum
                } else {
                    1.0 - column[call]
                }
            }
            None => 1.0 - column[call],
        };
        if let Some(reports) = percentiles.as_deref_mut() {
            if read.key[i] != 0 || !exclude_zero_flows {
                while reports.len() < i + 1 {
                    reports.push(PercentileReport::default());
                }
                reports[i].add_probability(*slot);
            }
        }
    }
    result
}

/// `computeLowestQBaseTP`: per flow, the `tp` of the lowest-quality base in the first half of the
/// hmer, the key walked over the bases in the order both are stored in.
///
/// `None` where the reference indexes past an array, which it does not catch.
pub fn lowest_q_base_tp(key: &[i32], tp: &[i8], qualities: &[u8]) -> Result<Vec<i8>, usize> {
    let mut result = vec![0i8; key.len()];
    let mut sequence = 0usize;
    for (i, &hmer) in key.iter().enumerate() {
        if hmer == 0 {
            continue;
        }
        let at = |array_len: usize, index: usize| {
            if index < array_len {
                Ok(index)
            } else {
                Err(index)
            }
        };
        result[i] = tp[at(tp.len(), sequence)?];
        let mut lowest = qualities[at(qualities.len(), sequence)?];
        let scan = (hmer as usize).div_ceil(2);
        for j in 1..scan {
            let q = qualities[at(qualities.len(), sequence + j)?];
            if q < lowest {
                result[i] = tp[at(tp.len(), sequence + j)?];
                lowest = q;
            }
        }
        sequence += hmer as usize;
    }
    Ok(result)
}

/// `collectReadProbs`: the matrix flattened flow by flow, and each flow's mean call, which divides
/// by the sum of the hmer lengths rather than of the probabilities.
pub fn read_probs(read: &crate::flow_based_read::FlowRead) -> (Vec<f64>, Vec<f64>) {
    let max_hmer = read.max_hmer;
    let mut probs = Vec::with_capacity(read.key.len() * (max_hmer as usize + 1));
    let mut mean_call = Vec::with_capacity(read.key.len());
    for flow in 0..read.key.len() {
        let mut call = 0.0;
        let mut total = 0i32;
        for hmer in 0..=max_hmer {
            let p = read.prob(flow, hmer);
            probs.push(p);
            call += p * f64::from(hmer);
            total += hmer;
        }
        mean_call.push(call / f64::from(total));
    }
    (probs, mean_call)
}

/// `new DecimalFormat("0.0#####")`: at least one integer and one fraction digit, at most six
/// fraction digits, half to even on the shortest decimal.
pub fn error_format(value: f64) -> String {
    let body = gatk_annotation::decimal_format::DecimalFormat::new(6).format(value);
    if !value.is_finite() {
        return body;
    }
    let (sign, digits) = match body.strip_prefix('-') {
        Some(rest) => ("-", rest),
        None => ("", body.as_str()),
    };
    let (whole, fraction) = digits.split_once('.').unwrap_or((digits, ""));
    let whole = if whole.is_empty() { "0" } else { whole };
    let fraction = if fraction.is_empty() { "0" } else { fraction };
    format!("{sign}{whole}.{fraction}")
}

/// `Precision.round(x, scale)`: `Double.toString(x)` as a `BigDecimal`, set to `scale` places
/// HALF_UP, back to a double, with a zero taking the sign of `x`.
pub fn precision_round(x: f64, scale: usize) -> f64 {
    if x.is_infinite() {
        return x;
    }
    if x.is_nan() {
        return f64::NAN;
    }
    let text = gatk_engine::tsv_table::java_double_to_string(x);
    let (negative, text) = match text.strip_prefix('-') {
        Some(rest) => (true, rest.to_string()),
        None => (false, text),
    };
    let (mantissa, exponent) = match text.split_once('E') {
        Some((m, e)) => (m.to_string(), e.parse::<i64>().unwrap_or(0)),
        None => (text.clone(), 0),
    };
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((&mantissa, ""));
    let digits: Vec<u8> = format!("{whole}{fraction}")
        .bytes()
        .map(|b| b - b'0')
        .collect();
    // value = digits * 10^-(places)
    let places = fraction.len() as i64 - exponent;
    let rounded = if places <= scale as i64 {
        x.abs()
    } else {
        let drop = (places - scale as i64) as usize;
        let keep = digits.len().saturating_sub(drop);
        let mut kept: Vec<u8> = digits[..keep].to_vec();
        let first_dropped = if drop <= digits.len() {
            digits[keep]
        } else {
            0
        };
        if first_dropped >= 5 {
            let mut carry = true;
            for digit in kept.iter_mut().rev() {
                if !carry {
                    break;
                }
                if *digit == 9 {
                    *digit = 0;
                } else {
                    *digit += 1;
                    carry = false;
                }
            }
            if carry {
                kept.insert(0, 1);
            }
        }
        let integer: String = if kept.is_empty() {
            "0".to_string()
        } else {
            kept.iter().map(|d| (b'0' + d) as char).collect()
        };
        format!("{integer}e-{scale}").parse::<f64>().unwrap_or(0.0)
    };
    let signed = if negative { -rounded } else { rounded };
    if signed == 0.0 {
        0.0 * x
    } else {
        signed
    }
}

/// The qual report: `BooleanAccumulator.newReport(61, 101, 202, 4)`, a tree four levels deep kept
/// here as one flat array of cells per level.
pub struct QualReport {
    top: Vec<Accumulator>,
    hmer: Vec<Accumulator>,
    deviation: Vec<Accumulator>,
    base: Vec<Accumulator>,
}

/// The dimensions of the four levels.
const QUALS: usize = QUAL_VALUE_MAX + 1;
const HMERS: usize = HMER_VALUE_MAX + 1;
const DEVIATIONS: usize = (HMER_VALUE_MAX + 1) * 2;
const BASES: usize = BASE_VALUE_MAX + 1;

impl Default for QualReport {
    fn default() -> Self {
        QualReport {
            top: vec![Accumulator::default(); QUALS],
            hmer: vec![Accumulator::default(); QUALS * HMERS],
            deviation: vec![Accumulator::default(); QUALS * HMERS * DEVIATIONS],
            base: vec![Accumulator::default(); QUALS * HMERS * DEVIATIONS * BASES],
        }
    }
}

/// `BooleanAccumulator.add`'s clipping: a bin out of range is logged and clamped.
fn clamp(bin: i64, length: usize) -> usize {
    bin.max(0).min(length as i64 - 1) as usize
}

impl QualReport {
    /// `qualReport[qual].add(same, hmer, deviationBin, baseBin)`.
    pub fn add(&mut self, same: bool, qual: usize, hmer: i64, deviation: i64, base: i64) {
        self.top[qual].add(same);
        let h = qual * HMERS + clamp(hmer, HMERS);
        self.hmer[h].add(same);
        let d = h * DEVIATIONS + clamp(deviation, DEVIATIONS);
        self.deviation[d].add(same);
        let b = d * BASES + clamp(base, BASES);
        self.base[b].add(same);
    }
}

/// `addToQualReport`: each flow's error probability rounded to seven places and made a phred, and
/// whether the read's call there is the reference's, filed by the call, the deviation and the
/// flow's base on the forward strand. A read whose key is not as long as the reference's files
/// nothing.
pub fn add_to_qual_report(
    report: &mut QualReport,
    read_key: &[i32],
    reference_key: &[i32],
    flow_order: &[u8],
    reverse: bool,
    error_probabilities: &[f64],
) {
    if read_key.len() != reference_key.len() {
        return;
    }
    for flow in 0..read_key.len() {
        let probability = precision_round(error_probabilities[flow], QUAL_VALUE_MAX / 10 + 1);
        let qual = crate::series_stats::java_double_to_int(
            (-10.0 * jmath::math::log10(probability)).ceil(),
        );
        let deviation = read_key[flow] - reference_key[flow];
        if qual >= 0 && (qual as usize) < QUALS {
            let base = flow_order[flow % flow_order.len()];
            let base = if reverse {
                match base {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'C' => b'G',
                    b'G' => b'C',
                    other => other,
                }
            } else {
                base
            };
            let base_bin = DEFAULT_FLOW_ORDER
                .bytes()
                .position(|b| b == base)
                .map_or(-1, |at| at as i64);
            report.add(
                deviation == 0,
                qual as usize,
                i64::from(read_key[flow]),
                deviation_to_bin(deviation) as i64,
                base_bin,
            );
        }
    }
}

/// One table of the report, written as `GATKReportTable.write` writes it, from rows produced
/// twice: once to size the columns and once to print them. The four-level table has five million
/// rows when its zeros are kept, which a table held in memory cell by cell cannot afford.
fn write_table<F>(
    out: &mut String,
    name: &str,
    description: &str,
    columns: &[(&str, &str)],
    rows: F,
) where
    F: Fn(&mut dyn FnMut(Vec<String>)),
{
    use gatk_engine::gatk_report::is_right_align;
    use std::fmt::Write as _;
    let mut widths: Vec<usize> = columns.iter().map(|(name, _)| name.len()).collect();
    let mut right: Vec<bool> = vec![true; columns.len()];
    let mut count = 0usize;
    rows(&mut |row: Vec<String>| {
        count += 1;
        for (i, value) in row.iter().enumerate() {
            widths[i] = widths[i].max(value.chars().count());
            if !is_right_align(value) {
                right[i] = false;
            }
        }
    });
    let _ = write!(out, "#:GATKTable:{}:{}", columns.len(), count);
    for (_, format) in columns {
        let _ = write!(out, ":{format}");
    }
    out.push_str(":;\n");
    let _ = writeln!(out, "#:GATKTable:{name}:{description}");
    for (i, (column, _)) in columns.iter().enumerate() {
        if i > 0 {
            out.push_str("  ");
        }
        let _ = write!(out, "{:<width$}", column, width = widths[i]);
    }
    out.push('\n');
    rows(&mut |row: Vec<String>| {
        for (i, value) in row.iter().enumerate() {
            if i > 0 {
                out.push_str("  ");
            }
            if right[i] {
                let _ = write!(out, "{:>width$}", value, width = widths[i]);
            } else {
                let _ = write!(out, "{:<width$}", value, width = widths[i]);
            }
        }
        out.push('\n');
    });
    out.push('\n');
}

/// `%f`: six places, HALF_UP.
fn six(value: f64) -> String {
    if value.is_finite() {
        gatk_engine::java_format::format_decimals(value, 6)
    } else {
        gatk_engine::tsv_table::java_double_to_string(value)
    }
}

/// `closeTool`'s report: a `GATKReport` over four tables, which it keeps in a `TreeMap` by name,
/// so the percentile table comes first whatever order they were handed in.
pub fn report_text(
    report: &QualReport,
    percentiles: &[PercentileReport],
    quality_percentiles: &str,
    omit_zeros: bool,
) -> String {
    let mut out = String::from("#:GATKReport.v1.1:4\n");

    // PhredBinAccumulator
    let names = percentile_columns(quality_percentiles);
    let mut columns: Vec<(&str, &str)> = Vec::new();
    for (i, name) in names.iter().enumerate() {
        columns.push((name.as_str(), if i < 2 { "%d" } else { "%f" }));
    }
    write_table(
        &mut out,
        PERCENTILE_TABLE_NAME,
        PERCENTILE_TABLE_NAME,
        &columns,
        |emit| {
            for (index, row) in percentiles.iter().enumerate() {
                let values = row.row(index, quality_percentiles);
                let mut cells = vec![index.to_string(), row.stats.count().to_string()];
                cells.extend(values[2..].iter().map(|value| six(*value)));
                emit(cells);
            }
        },
    );

    // qualReport
    write_table(
        &mut out,
        &one_level_table_name("qual"),
        "error rate per qual",
        &[
            ("qual", "%d"),
            ("count", "%d"),
            ("error", "%f"),
            ("phred", "%d"),
        ],
        |emit| {
            for (i, cell) in report.top.iter().enumerate() {
                if omit_zeros && i != 0 && cell.count() == 0 {
                    continue;
                }
                let rate = cell.false_rate();
                emit(vec![
                    i.to_string(),
                    cell.count().to_string(),
                    six(rate),
                    phred(rate, cell.count(), DEFAULT_RATIO_THRESHOLD).to_string(),
                ]);
            }
        },
    );

    // qual_hmerReport
    write_table(
        &mut out,
        &two_level_table_name("qual", "hmer"),
        "error rate per qual by hmer",
        &[
            ("qual", "%d"),
            ("hmer", "%d"),
            ("count", "%d"),
            ("error", "%f"),
        ],
        |emit| {
            for i in 0..QUALS {
                for j in 0..HMERS {
                    let cell = &report.hmer[i * HMERS + j];
                    if omit_zeros && (i != 0 || j != 0) && cell.count() == 0 {
                        continue;
                    }
                    emit(vec![
                        i.to_string(),
                        j.to_string(),
                        cell.count().to_string(),
                        six(cell.false_rate()),
                    ]);
                }
            }
        },
    );

    // qual_hmer_deviation_base_Report
    write_table(
        &mut out,
        &four_level_table_name("qual", "hmer", "deviation", "base"),
        "error rate per qual by hmer and deviation",
        &[
            ("qual", "%d"),
            ("hmer", "%d"),
            ("deviation", "%s"),
            ("base", "%s"),
            ("count", "%d"),
        ],
        |emit| {
            for i in 0..QUALS {
                for j in 0..HMERS {
                    for k in 0..DEVIATIONS {
                        for m in 0..BASES {
                            let cell = &report.base[((i * HMERS + j) * DEVIATIONS + k) * BASES + m];
                            if omit_zeros
                                && (i != 0 || j != 0 || k != 0 || m != 0)
                                && cell.count() == 0
                            {
                                continue;
                            }
                            emit(vec![
                                i.to_string(),
                                j.to_string(),
                                bin_to_deviation(k),
                                bin_to_base(m).to_string(),
                                cell.count().to_string(),
                            ]);
                        }
                    }
                }
            }
        },
    );
    out
}
