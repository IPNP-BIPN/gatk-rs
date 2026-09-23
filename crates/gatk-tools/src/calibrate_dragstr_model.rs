//! `CalibrateDragstrModel`: the DRAGstr parameter table estimated from the reads over a
//! reference's repeats.
//!
//! What is ported is the estimation and the table it writes: the precomputed error probabilities,
//! the search that picks a GP and an API per period, the grouping of repeat lengths that have too
//! little data, and the file's own layout. Scanning the reference and piling up the reads are not
//! ported.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.dragstr.DragstrParametersEstimator`,
//! `org.broadinstitute.hellbender.tools.dragstr.DragstrHyperParameters` and
//! `org.broadinstitute.hellbender.utils.dragstr.DragstrParams` in GATK 4.6.2.0.

use gatk_engine::math_utils::log10_sum_log10;

/// `MathUtils.log10OneMinusPow10`, which is `log10(1 - 10^a)`.
///
/// A positive argument is a NaN and a zero one is negative infinity, both being answered before
/// any arithmetic, so the two are exact rather than the limits of a computation.
pub fn log10_one_minus_pow10(a: f64) -> f64 {
    if a > 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return f64::NEG_INFINITY;
    }
    log1mexp(a * std::f64::consts::LN_10) / std::f64::consts::LN_10
}

/// `NaturalLogUtils.log1mexp`, whose branch is at `log(0.5)`.
///
/// Below the threshold `log1p(-exp(a))` keeps its precision and above it `log(-expm1(a))` does,
/// which is why the function is written as two formulae rather than one.
fn log1mexp(a: f64) -> f64 {
    if a > 0.0 {
        return f64::NAN;
    }
    if a == 0.0 {
        return f64::NEG_INFINITY;
    }
    if a < 0.5f64.ln() {
        (-a.exp()).ln_1p()
    } else {
        (-a.exp_m1()).ln()
    }
}

/// `MathUtils.LOG10_ONE_HALF`, which is `Math.log10(0.5)` and so the negation of `log10(2)`.
pub const LOG10_ONE_HALF: f64 = -std::f64::consts::LOG10_2;

/// The hyper-parameters, which decide the table's shape as well as the search's range.
#[derive(Debug, Clone, PartialEq)]
pub struct HyperParameters {
    /// `--gp-values`, in Phred scale.
    pub phred_gp_values: Vec<f64>,
    /// `--api-values`, in Phred scale.
    pub phred_api_values: Vec<f64>,
    /// `--gop-values`, which are NOT searched for: they are the GOP column as written.
    pub phred_gop_values: Vec<f64>,
    pub het_to_hom_ratio: f64,
    pub min_loci_count: usize,
    pub api_mono_threshold: f64,
    pub max_period: usize,
    pub max_repeat_length: usize,
}

/// `<start>:<step>:<end>`, which is how the three value arguments are written.
///
/// The end is INCLUDED when the step divides the range, so `10:1.0:50` is forty-one values and
/// not forty.
pub fn value_range(start: f64, step: f64, end: f64) -> Vec<f64> {
    let mut values = Vec::new();
    let mut i = 0;
    loop {
        let value = start + step * i as f64;
        if value > end + 1e-9 {
            break;
        }
        values.push(value);
        i += 1;
    }
    values
}

impl Default for HyperParameters {
    /// The defaults the tool ships: `10:1.0:50`, `0:1.0:40`, `10:.25:50`, and a table eight
    /// periods by twenty repeat lengths.
    fn default() -> Self {
        HyperParameters {
            phred_gp_values: value_range(10.0, 1.0, 50.0),
            phred_api_values: value_range(0.0, 1.0, 40.0),
            phred_gop_values: value_range(10.0, 0.25, 50.0),
            het_to_hom_ratio: 2.0,
            min_loci_count: 50,
            api_mono_threshold: 3.0,
            max_period: 8,
            max_repeat_length: 20,
        }
    }
}

/// The tables the estimator computes once, before any data is looked at.
#[derive(Debug, Clone, PartialEq)]
pub struct Precomputed {
    /// `[gp index][period - 1][repeats - 1]`.
    pub log10_p_error: Vec<Vec<Vec<f64>>>,
    pub log10_p_correct: Vec<Vec<Vec<f64>>>,
    /// The first GP index each period's search may start at.
    pub min_gp_index_by_period: Vec<usize>,
}

/// The precomputation, whose two tables depend on the GP values and the table's shape alone.
///
/// The per-position correct probability is `log10(1 - 10^(-log10(0.5) + log10Gp))`, and a repeat's
/// is that raised to its length IN BASES, which is the repeat count times the period. The error
/// probability is one less that, so the two always sum to one.
pub fn precompute(parameters: &HyperParameters) -> Precomputed {
    let log10_gp: Vec<f64> = parameters
        .phred_gp_values
        .iter()
        .map(|phred| -0.1 * phred)
        .collect();
    let mut log10_p_error =
        vec![vec![vec![0.0; parameters.max_repeat_length]; parameters.max_period]; log10_gp.len()];
    let mut log10_p_correct = log10_p_error.clone();
    for (i, log10_gp) in log10_gp.iter().enumerate() {
        for k in 0..parameters.max_period {
            let period = k + 1;
            let per_position = log10_one_minus_pow10(-LOG10_ONE_HALF + log10_gp);
            for j in 0..parameters.max_repeat_length {
                let bases = ((j + 1) * period) as f64;
                log10_p_correct[i][k][j] = bases * per_position;
                log10_p_error[i][k][j] = log10_one_minus_pow10(log10_p_correct[i][k][j]);
            }
        }
    }
    Precomputed {
        log10_p_error,
        log10_p_correct,
        min_gp_index_by_period: (0..parameters.max_period)
            .map(|i| min_gp_index(parameters, i + 1))
            .collect(),
    }
}

/// The first GP index a period's search may start at.
///
/// The formula is a transcription of Illumina's own script, and the `20.0` it was written with is
/// the DEFAULT maximum repeat length rather than a constant, so changing `--max-repeats` moves
/// every period's floor. The search for the value then falls back on the insertion point, and the
/// tolerance it applies is read from the API values rather than the GP ones.
pub fn min_gp_index(parameters: &HyperParameters, period: usize) -> usize {
    let gp_min = (-10.0
        * (1.0 - 0.5f64.powf((1.0 / (parameters.max_repeat_length * period) as f64) / 2.0))
            .log10())
    .ceil();
    match binary_search(&parameters.phred_gp_values, gp_min) {
        Ok(index) => index,
        Err(insertion) => {
            if insertion >= 2 && (gp_min - parameters.phred_api_values[insertion - 2]).abs() < 0.001
            {
                insertion - 2
            } else {
                insertion
            }
        }
    }
}

/// `Arrays.binarySearch` over a sorted array of doubles, as its two outcomes.
fn binary_search(values: &[f64], target: f64) -> Result<usize, usize> {
    values.binary_search_by(|value| value.partial_cmp(&target).expect("no NaN"))
}

/// One site the estimator was given: how many reads covered it and how many carried an indel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Case {
    pub depth: i32,
    pub indels: i32,
}

/// The cases, stratified by period and repeat length, which is how the estimator reads them.
#[derive(Debug, Clone, PartialEq)]
pub struct Cases {
    pub max_period: usize,
    pub max_repeat_length: usize,
    /// `[period - 1][repeats - 1]`.
    pub by_shape: Vec<Vec<Vec<Case>>>,
}

impl Cases {
    pub fn empty(max_period: usize, max_repeat_length: usize) -> Cases {
        Cases {
            max_period,
            max_repeat_length,
            by_shape: vec![vec![Vec::new(); max_repeat_length]; max_period],
        }
    }

    pub fn get(&self, period: usize, repeats: usize) -> &[Case] {
        &self.by_shape[period - 1][repeats - 1]
    }

    pub fn add(&mut self, period: usize, repeats: usize, case: Case) {
        self.by_shape[period - 1][repeats - 1].push(case);
    }
}

/// `log10ProbFunc`: the likelihood of one site under one GP and one API.
///
/// The three terms are the three genotypes. The homozygous-variant one is only included when
/// EVERY read carried the indel, the reference's own comment allowing that this is not quite
/// right: an error that reverts to the reference is possible and unaccounted for.
pub fn log10_prob(
    depth: i32,
    indels: i32,
    log10_p_error: f64,
    log10_p_correct: f64,
    log10_p_hom_ref: f64,
    log10_p_het: f64,
    log10_p_hom_var: f64,
) -> f64 {
    log10_sum_log10(&[
        log10_p_hom_ref + indels as f64 * log10_p_error + (depth - indels) as f64 * log10_p_correct,
        log10_p_het + depth as f64 * LOG10_ONE_HALF,
        if depth == indels {
            log10_p_hom_var
        } else {
            f64::NEG_INFINITY
        },
    ])
}

/// What one repeat-length group's search settles on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Estimate {
    pub gp: f64,
    pub gcp: f64,
    pub api: f64,
}

/// `estimatePeriodRepeatInterval`: the GP and API that maximise the likelihood over a group.
///
/// The inner loop ABORTS a GP as soon as its running total falls below the best found so far,
/// which is safe only because every term is negative: the accumulator can never climb back. The
/// abort is what makes the search's answer depend on the ORDER the cases are visited in, though
/// not on which answer it reaches.
///
/// GCP is not searched for at all: it is ten over the period, always.
pub fn estimate_interval(
    period: usize,
    repeats: std::ops::RangeInclusive<usize>,
    parameters: &HyperParameters,
    precomputed: &Precomputed,
    cases: &Cases,
) -> Estimate {
    let mut max_api_index = 0usize;
    let mut max_gp_index = 0usize;
    let mut max_log10_prob = f64::NEG_INFINITY;
    let period_index = period - 1;
    let log10_het_over_hom_var = parameters.het_to_hom_ratio.log10();
    let max_log10_p_het = log10_het_over_hom_var - (1.0 + parameters.het_to_hom_ratio).log10();
    for i in 0..parameters.phred_api_values.len() {
        let log10_api = -0.1 * parameters.phred_api_values[i];
        let log10_p_het = log10_api.min(max_log10_p_het);
        let log10_p_hom_var = log10_p_het - log10_het_over_hom_var;
        let log10_p_hom_ref =
            log10_one_minus_pow10(log10_sum_log10(&[log10_p_het, log10_p_hom_var]));
        'gp: for j in
            precomputed.min_gp_index_by_period[period_index]..parameters.phred_gp_values.len()
        {
            let mut accumulator = 0.0;
            for r in repeats.clone() {
                let log10_p_error = precomputed.log10_p_error[j][period_index][r - 1];
                let log10_p_correct = precomputed.log10_p_correct[j][period_index][r - 1];
                for case in cases.get(period, r) {
                    accumulator += log10_prob(
                        case.depth,
                        case.indels,
                        log10_p_error,
                        log10_p_correct,
                        log10_p_hom_ref,
                        log10_p_het,
                        log10_p_hom_var,
                    );
                    if accumulator < max_log10_prob {
                        continue 'gp;
                    }
                }
            }
            if accumulator > max_log10_prob {
                max_api_index = i;
                max_gp_index = j;
                max_log10_prob = accumulator;
            }
        }
    }
    Estimate {
        gp: parameters.phred_gp_values[max_gp_index],
        gcp: 10.0 / period as f64,
        api: parameters.phred_api_values[max_api_index],
    }
}

/// The two flanks of repeat lengths that hold too little data to be estimated on their own.
///
/// Each is found by accumulating sizes inwards until `--min-loci-count` is reached, so a period
/// with no data at all leaves the left flank at the maximum and the right at one.
pub fn flanks(period: usize, parameters: &HyperParameters, cases: &Cases) -> (usize, usize) {
    let mut accumulated = 0usize;
    let mut left = 0usize;
    while left < parameters.max_repeat_length {
        left += 1;
        accumulated += cases.get(period, left).len();
        if accumulated >= parameters.min_loci_count {
            break;
        }
    }
    let mut accumulated = 0usize;
    let mut right = parameters.max_repeat_length;
    while right > 1 {
        right -= 1;
        accumulated += cases.get(period, right).len();
        if accumulated >= parameters.min_loci_count {
            break;
        }
    }
    (left, right)
}

/// The repeat-length groups one period's estimation starts from.
///
/// When the two flanks have not crossed the groups are `[1..left]`, then every single repeat
/// length up to the right flank, then `[right+1..max]`. When they HAVE crossed there is too
/// little data to split at all and the whole range is one group.
pub fn initial_groups(
    parameters: &HyperParameters,
    left: usize,
    right: usize,
) -> Vec<std::ops::RangeInclusive<usize>> {
    if right < left {
        return vec![1..=parameters.max_repeat_length];
    }
    let mut groups = vec![1..=left];
    for single in (left + 1)..=right {
        groups.push(single..=single);
    }
    groups.push((right + 1)..=parameters.max_repeat_length);
    groups
}

/// One period's estimation: the groups, and the merging that happens when a group's estimate
/// fails to decrease.
///
/// A group whose GP is above the previous group's, or whose API is above it by more than
/// `--api-mono-threshold`, is MERGED BACK into the previous group and re-estimated, which is what
/// makes both columns monotone across the row.
pub fn estimate_period(
    period: usize,
    parameters: &HyperParameters,
    precomputed: &Precomputed,
    cases: &Cases,
) -> Vec<(std::ops::RangeInclusive<usize>, Estimate)> {
    let (left, right) = flanks(period, parameters, cases);
    let mut pending: std::collections::VecDeque<std::ops::RangeInclusive<usize>> =
        initial_groups(parameters, left, right).into();
    let mut done: Vec<(std::ops::RangeInclusive<usize>, Estimate)> = Vec::new();
    while let Some(next) = pending.pop_front() {
        let estimate = estimate_interval(period, next.clone(), parameters, precomputed, cases);
        let accept = match done.last() {
            None => true,
            Some((_, last)) => {
                last.gp >= estimate.gp && last.api + parameters.api_mono_threshold >= estimate.api
            }
        };
        if accept {
            done.push((next, estimate));
        } else {
            let (last, _) = done.pop().expect("a previous group");
            pending.push_front(*last.start()..=*next.end());
        }
    }
    done
}

// ================================================================================================
// The file.
// ================================================================================================

/// The three blocks the parameter file carries, in this order.
pub const BLOCKS: [&str; 3] = ["GOP", "GCP", "API"];

/// The width each number is written into, and the separator between them.
///
/// Five is exactly what `10.00` takes, so a value of ten or more starts flush against the
/// separator while a single-digit one carries a leading space of its own.
pub const COLUMN_WIDTH: usize = 5;
pub const COLUMN_SEPARATOR: &str = "  ";

/// One row of one block.
pub fn row(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| format!("{value:>COLUMN_WIDTH$.2}"))
        .collect::<Vec<_>>()
        .join(COLUMN_SEPARATOR)
}

/// The column header: the repeat lengths, right-aligned in the same width.
pub fn column_header(max_repeat_length: usize) -> String {
    (1..=max_repeat_length)
        .map(|repeats| format!("{repeats:>COLUMN_WIDTH$}"))
        .collect::<Vec<_>>()
        .join(COLUMN_SEPARATOR)
}

/// The whole table: the shape is the hyper-parameters', so a period with no data still has a row.
///
/// A period the estimation never reached keeps the DEFAULTS, and the file gives no sign of which
/// rows those are: an estimated row and a default one are written the same way.
pub fn table(parameters: &HyperParameters, rows: &[(Vec<f64>, Vec<f64>, Vec<f64>)]) -> String {
    let mut text = column_header(parameters.max_repeat_length);
    text.push('\n');
    for (index, name) in BLOCKS.iter().enumerate() {
        text.push_str(name);
        text.push_str(":\n");
        for row_values in rows {
            let values = match index {
                0 => &row_values.0,
                1 => &row_values.1,
                _ => &row_values.2,
            };
            text.push_str(&row(values));
            text.push('\n');
        }
    }
    text
}

/// The GCP row of one period, which is ten over the period repeated across the row.
pub fn gcp_row(period: usize, max_repeat_length: usize) -> Vec<f64> {
    vec![10.0 / period as f64; max_repeat_length]
}

// ================================================================================================
// The traversal: the sites, the reads over them, the downsampling and the file.
// ================================================================================================

/// `DragstrParams.DEFAULT_GOP`, Illumina's table, eight periods by twenty repeat lengths.
pub const DEFAULT_GOP: [[f64; 20]; 8] = [
    [
        45.00, 45.00, 45.00, 45.00, 45.00, 45.00, 40.50, 33.50, 28.00, 24.00, 21.75, 21.75, 21.75,
        21.75, 21.75, 21.75, 21.75, 21.75, 21.75, 21.75,
    ],
    [
        39.50, 39.50, 39.50, 39.50, 36.00, 30.00, 27.25, 25.00, 24.25, 24.75, 26.25, 26.25, 26.25,
        26.25, 26.25, 26.25, 26.25, 26.25, 26.25, 26.75,
    ],
    [
        38.50, 41.00, 41.00, 41.00, 41.00, 37.50, 35.25, 34.75, 34.75, 33.25, 33.25, 33.25, 32.50,
        30.75, 28.50, 29.00, 29.00, 29.00, 29.00, 29.00,
    ],
    [
        37.50, 39.00, 39.00, 37.75, 34.00, 34.00, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25, 30.25,
        30.25, 30.25, 31.75, 31.75, 31.75, 31.75, 31.75,
    ],
    [
        37.00, 40.00, 40.00, 40.00, 36.00, 35.00, 24.50, 24.50, 24.50, 24.50, 22.50, 22.50, 22.50,
        23.50, 23.50, 23.50, 23.50, 23.50, 23.50, 23.50,
    ],
    [
        36.25, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00,
        40.00, 40.00, 40.00, 40.00, 40.00, 40.00, 40.00,
    ],
    [
        36.00, 40.50, 40.50, 40.50, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75,
        20.75, 20.75, 20.75, 20.75, 20.75, 20.75, 20.75,
    ],
    [
        36.25, 39.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75,
        32.75, 32.75, 32.75, 32.75, 32.75, 32.75, 32.75,
    ],
];

/// `DragstrParams.DEFAULT_API`.
pub const DEFAULT_API: [[f64; 20]; 8] = [
    [
        39.00, 39.00, 37.00, 35.00, 32.00, 26.00, 20.00, 16.00, 12.00, 10.00, 8.00, 7.00, 7.00,
        6.00, 6.00, 5.00, 5.00, 4.00, 4.00, 4.00,
    ],
    [
        30.00, 30.00, 29.00, 22.00, 17.00, 14.00, 11.00, 8.00, 6.00, 5.00, 4.00, 4.00, 3.00, 3.00,
        3.00, 3.00, 3.00, 3.00, 2.00, 2.00,
    ],
    [
        27.00, 27.00, 25.00, 18.00, 14.00, 12.00, 9.00, 7.00, 5.00, 4.00, 3.00, 3.00, 3.00, 3.00,
        2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    ],
    [
        27.00, 27.00, 18.00, 9.00, 9.00, 9.00, 9.00, 3.00, 3.00, 3.00, 3.00, 3.00, 2.00, 2.00,
        2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    ],
    [
        29.00, 29.00, 18.00, 8.00, 8.00, 8.00, 4.00, 3.00, 3.00, 3.00, 2.00, 2.00, 2.00, 2.00,
        2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    ],
    [
        25.00, 25.00, 10.00, 10.00, 10.00, 4.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00,
        3.00, 3.00, 3.00, 3.00, 3.00, 3.00,
    ],
    [
        21.00, 21.00, 11.00, 11.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00,
        5.00, 5.00, 5.00, 5.00, 5.00, 5.00,
    ],
    [
        18.00, 18.00, 10.00, 6.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00,
        4.00, 4.00, 4.00, 4.00, 4.00, 4.00,
    ],
];

/// `DragstrParams.DEFAULT`: the three blocks, eight periods by twenty, whatever the hyper-parameters
/// said. GCP is `Math.round(1000.0 / period) / 100.0`.
pub fn default_rows() -> Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    (0..8)
        .map(|index| {
            let period = (index + 1) as f64;
            let gcp = (1000.0 / period + 0.5).floor() / 100.0;
            (
                DEFAULT_GOP[index].to_vec(),
                vec![gcp; 20],
                DEFAULT_API[index].to_vec(),
            )
        })
        .collect()
}

/// `MathUtils.doubles(start, limit, step)`, which a `DoubleSequence` is built from.
pub fn java_doubles(start: f64, limit: f64, step: f64) -> Result<Vec<f64>, String> {
    if !start.is_finite() || !limit.is_finite() || !step.is_finite() {
        return Err("the start, limit and step must be finite".to_string());
    }
    let tolerance = 10f64.powf((0f64.min(step.abs().log10() - 3.0)).floor());
    let diff = limit - start;
    if diff.abs() < tolerance {
        return Ok(vec![start]);
    }
    if diff * step <= 0.0 {
        return Err(
            "the difference between start and end must have the same sign as the step".to_string(),
        );
    }
    let length = (1.0 + tolerance + diff / step).floor() as usize;
    let mut result: Vec<f64> = (0..length).map(|i| start + step * i as f64).collect();
    if let Some(last) = result.last_mut() {
        if (*last - limit).abs() <= tolerance {
            *last = limit;
        }
    }
    Ok(result)
}

/// A `DoubleSequence`'s text, `start:step:limit`, parsed.
pub fn double_sequence(text: &str) -> Result<Vec<f64>, String> {
    let parts: Vec<&str> = text.split(':').collect();
    let number = |part: &str| -> Option<f64> {
        let valid = !part.is_empty()
            && part
                .trim_start_matches(['+', '-'])
                .chars()
                .all(|c| c.is_ascii_digit() || matches!(c, '.' | 'e' | 'E' | '+' | '-'));
        if valid {
            part.parse::<f64>().ok()
        } else {
            None
        }
    };
    match parts.as_slice() {
        [start, step, limit] => match (number(start), number(step), number(limit)) {
            (Some(start), Some(step), Some(limit)) => java_doubles(start, limit, step),
            _ => Err(format!("invalid double sequence specificatior: {text}")),
        },
        _ => Err(format!("invalid double sequence specificatior: {text}")),
    }
}

/// One site of the STR table, as `DragstrLocus` holds it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Site {
    pub contig: i32,
    pub start: i64,
    pub period: i32,
    /// The span in bases, a short on disk.
    pub length: i32,
    pub mask: i64,
}

impl Site {
    pub fn end(&self) -> i64 {
        self.start + i64::from(self.length) - 1
    }

    /// `getRepeats`: whole units in the span.
    pub fn repeats(&self) -> i32 {
        if self.period == 0 {
            0
        } else {
            self.length / self.period
        }
    }
}

/// `sites.bin`, 23 big-endian bytes a site.
pub fn read_sites(bytes: &[u8]) -> Vec<Site> {
    bytes
        .chunks_exact(23)
        .map(|chunk| Site {
            contig: i32::from_be_bytes(chunk[0..4].try_into().expect("four bytes")),
            start: i64::from_be_bytes(chunk[4..12].try_into().expect("eight bytes")),
            period: i32::from(chunk[12] as i8),
            length: i32::from(i16::from_be_bytes(
                chunk[13..15].try_into().expect("two bytes"),
            )),
            mask: i64::from_be_bytes(chunk[15..23].try_into().expect("eight bytes")),
        })
        .collect()
}

/// One read as the case collector reads it, after the extended-MQ transformer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PileRead {
    pub start: i64,
    pub end: i64,
    pub mapping_quality: i32,
    pub supplementary: bool,
    /// The CIGAR as operator letters and lengths.
    pub cigar: Vec<(u8, i64)>,
}

/// `DragstrLocusCase`: what the reads over one site said.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiteCase {
    pub site: Site,
    pub depth: i32,
    pub indels: i32,
    pub min_mq: i32,
    pub n_sup: i32,
}

impl SiteCase {
    /// `qualifies`.
    pub fn qualifies(&self, min_depth: i32, min_mq: i32, max_sup: i32) -> bool {
        self.depth >= min_depth && self.n_sup <= max_sup && self.min_mq >= min_mq
    }

    /// `outputSiteDetails`' line.
    pub fn line(&self, fate: &str) -> String {
        format!(
            "{}:{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.site.contig,
            self.site.start - 1,
            self.site.period,
            self.site.repeats(),
            self.depth,
            self.indels,
            self.min_mq,
            self.n_sup,
            fate
        )
    }
}

/// `DragstrLocusCaseCollector` over the reads overlapping one site.
///
/// Only a read spanning the site PADDED on both sides counts. Every insertion that starts inside
/// the repeat or right after it, and every deletion that touches it, adds one to the indel count,
/// so a read with two such events counts twice.
pub fn collect(site: &Site, reads: &[&PileRead], padding: i64, contig_length: i64) -> SiteCase {
    let str_start = site.start;
    let str_end = site.end();
    let str_end_plus_one = str_end + 1;
    let padded_start = (str_start - padding).max(1);
    let padded_end = (str_end + padding).min(contig_length);
    let (mut n, mut k, mut n_sup, mut min_mq) = (0, 0, 0, 255);
    for read in reads {
        if read.start <= padded_start && read.end >= padded_end {
            if read.supplementary {
                n_sup += 1;
            }
            min_mq = min_mq.min(read.mapping_quality);
            let mut ref_pos = read.start;
            for (op, length) in &read.cigar {
                let insertion = *op == b'I' && ref_pos >= str_start && ref_pos <= str_end_plus_one;
                let deletion = *op == b'D' && ref_pos + length > str_start && ref_pos <= str_end;
                if insertion || deletion {
                    k += 1;
                }
                if matches!(op, b'M' | b'D' | b'N' | b'=' | b'X') {
                    ref_pos += length;
                }
                if ref_pos > str_end_plus_one {
                    break;
                }
            }
            n += 1;
        }
    }
    SiteCase {
        site: *site,
        depth: n,
        indels: k,
        min_mq,
        n_sup,
    }
}

/// `StratifiedDragstrLocusCases`: `[period - 1][repeats - 1]`, the repeats capped at the maximum.
#[derive(Debug, Clone, PartialEq)]
pub struct Stratified {
    pub cells: Vec<Vec<Vec<SiteCase>>>,
}

impl Stratified {
    pub fn new(max_period: usize, max_repeats: usize) -> Self {
        Stratified {
            cells: vec![vec![Vec::new(); max_repeats]; max_period],
        }
    }

    /// `add`, which indexes the period without a check: a site whose period is past the maximum is
    /// the reference's `ArrayIndexOutOfBoundsException`, reported here as the index and the length.
    pub fn add(&mut self, case: SiteCase) -> Result<(), (i64, usize)> {
        let period_index = i64::from(case.site.period) - 1;
        if period_index < 0 || period_index as usize >= self.cells.len() {
            return Err((period_index, self.cells.len()));
        }
        let row = &mut self.cells[period_index as usize];
        let repeat_index = ((case.site.repeats() - 1).max(0) as usize).min(row.len() - 1);
        row[repeat_index].push(case);
        Ok(())
    }

    /// `get`, whose repeat index is capped.
    pub fn get(&self, period: usize, repeats: usize) -> &[SiteCase] {
        let row = &self.cells[period - 1];
        &row[(repeats - 1).min(row.len() - 1)]
    }

    /// `qualifyingOnly`.
    pub fn qualifying(&self, min_depth: i32, min_mq: i32, max_sup: i32) -> Stratified {
        Stratified {
            cells: self
                .cells
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|cell| {
                            cell.iter()
                                .filter(|case| case.qualifies(min_depth, min_mq, max_sup))
                                .copied()
                                .collect()
                        })
                        .collect()
                })
                .collect(),
        }
    }

    /// The cases as the estimator reads them.
    pub fn as_cases(&self) -> Cases {
        let max_period = self.cells.len();
        let max_repeat_length = self.cells.first().map_or(0, Vec::len);
        let mut cases = Cases::empty(max_period, max_repeat_length);
        for (p, row) in self.cells.iter().enumerate() {
            for (r, cell) in row.iter().enumerate() {
                for case in cell {
                    cases.add(
                        p + 1,
                        r + 1,
                        Case {
                            depth: case.depth,
                            indels: case.indels,
                        },
                    );
                }
            }
        }
        cases
    }
}

/// `DECIMATION_MASKS_BY_BIT[j]`, which the static initializer leaves as `~(1 << j)`.
fn decimation_mask_by_bit(j: usize) -> i64 {
    !(1i64 << j)
}

/// `downSample` over one cell. Returns the cases kept and appends a `downsampled-out` line for each
/// one dropped. `Err` is the reference's array overrun, index and length.
pub fn downsample_cell(
    cell: &[SiteCase],
    min_decimation_bit: usize,
    downsample_size: usize,
    out: &mut Vec<String>,
) -> Result<Vec<SiteCase>, (usize, usize)> {
    if cell.len() <= downsample_size {
        return Ok(cell.to_vec());
    }
    let length = 64usize.saturating_sub(min_decimation_bit);
    let mut count_by_first_bit = vec![0usize; length];
    let mut zero_depth = 0usize;
    for case in cell {
        if case.depth <= 0 {
            zero_depth += 1;
            continue;
        }
        let mask = case.site.mask;
        let mut j = min_decimation_bit;
        while mask != 0 && j < 64 {
            let new_mask = mask & decimation_mask_by_bit(j);
            if new_mask != mask {
                if j >= length {
                    return Err((j, length));
                }
                count_by_first_bit[j] += 1;
                break;
            }
            j += 1;
        }
    }
    let mut final_size = cell.len() - zero_depth;
    let mut filter_mask = 0i64;
    let mut j = min_decimation_bit;
    while final_size > downsample_size && j < 64 {
        if j >= length {
            return Err((j, length));
        }
        final_size = final_size.saturating_sub(count_by_first_bit[j]);
        filter_mask |= !decimation_mask_by_bit(j);
        j += 1;
    }
    let mut kept = Vec::new();
    for case in cell {
        if (case.site.mask & filter_mask) == 0 && case.depth > 0 {
            kept.push(*case);
        } else {
            out.push(case.line("downsampled-out"));
        }
    }
    Ok(kept)
}

/// `MINIMUM_CASES_BY_PERIOD_AND_LENGTH`.
const MINIMUM_CASES: [&[usize]; 9] = [
    &[],
    &[0, 200, 200, 200, 200, 200, 200, 200, 200, 200, 0],
    &[0, 0, 200, 200, 200, 200, 0, 0, 0, 0, 0],
    &[0, 0, 200, 200, 200, 0, 0, 0, 0, 0, 0],
    &[0, 0, 200, 200, 0, 0, 0, 0, 0, 0, 0],
    &[0, 0, 200, 0, 0, 0, 0, 0, 0, 0, 0],
    &[0, 0, 200, 0, 0, 0, 0, 0, 0, 0, 0],
    &[0, 0, 200, 0, 0, 0, 0, 0, 0, 0, 0],
    &[0, 0, 200, 0, 0, 0, 0, 0, 0, 0, 0],
];

/// `isThereEnoughCases`, before `--force-estimation` is considered.
pub fn enough_cases(sites: &Stratified, max_period: usize, max_repeat_length: usize) -> bool {
    let max_p = max_period.min(MINIMUM_CASES.len() - 1);
    for (i, minimums) in MINIMUM_CASES.iter().enumerate().take(max_p + 1).skip(1) {
        let max_l = max_repeat_length.min(minimums.len() - 1);
        for (j, minimum) in minimums.iter().enumerate().take(max_l + 1).skip(1) {
            if sites.get(i, j).len() < *minimum {
                return false;
            }
        }
    }
    true
}

/// `Math.pow`, at run time and never folded, like `gatk_engine::math_utils::pow10`.
fn java_pow(base: f64, exponent: f64) -> f64 {
    std::hint::black_box(base).powf(exponent)
}

/// `DragstrParamsBuilder.gopCalculation`, a grid search from zero whose step accumulates.
fn gop_calculation(gp: f64, gcp: f64, period: usize, max_gop: f64, step: f64) -> f64 {
    let gp_prob = java_pow(10.0, -0.1 * gp);
    let mut best_gop = 0.0;
    let c = java_pow(10.0, -0.1 * gcp);
    let mut best_cost = f64::INFINITY;
    let mut gop = 0.0;
    while max_gop - gop > -0.001 {
        let g = java_pow(10.0, -0.1 * gop);
        let pr_gap = g * java_pow(c, (period - 1) as f64) * (1.0 - c);
        let pr_no_gap = java_pow(1.0 - 2.0 * g, (period + 1) as f64);
        let cost = (pr_gap / pr_no_gap - gp_prob).abs();
        if cost < best_cost {
            best_gop = gop;
            best_cost = cost;
        }
        gop += step;
    }
    best_gop
}

/// `DragstrParametersEstimator.estimate`, then `DragstrParamsBuilder.make(gopValues)`: every
/// period's groups, and GOP derived from each cell's GP by the grid over the GOP values.
pub fn estimate(
    parameters: &HyperParameters,
    cases: &Cases,
) -> Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let precomputed = precompute(parameters);
    let values = &parameters.phred_gop_values;
    let (min_gop, max_gop) = match values.len() {
        0 => (f64::NAN, f64::NAN),
        1 | 2 => (values[0], values[0]),
        n => (values[0].min(values[n - 1]), values[0].max(values[n - 1])),
    };
    let step = if values.len() < 2 {
        f64::NAN
    } else {
        values[1] - values[0]
    };
    (1..=parameters.max_period)
        .map(|period| {
            let mut gp = vec![0.0; parameters.max_repeat_length];
            let mut gcp = vec![0.0; parameters.max_repeat_length];
            let mut api = vec![0.0; parameters.max_repeat_length];
            for (range, estimate) in estimate_period(period, parameters, &precomputed, cases) {
                for repeats in range {
                    gp[repeats - 1] = estimate.gp;
                    gcp[repeats - 1] = estimate.gcp;
                    api[repeats - 1] = estimate.api;
                }
            }
            let gop = gp
                .iter()
                .zip(&gcp)
                .map(|(gp, gcp)| min_gop.max(gop_calculation(*gp, *gcp, period, max_gop, step)))
                .collect();
            (gop, gcp, api)
        })
        .collect()
}

/// `DragstrParamUtils.print`: a banner with the annotations, then the table.
pub fn params_file(
    annotations: &[(&str, String)],
    max_repeat_length: usize,
    rows: &[(Vec<f64>, Vec<f64>, Vec<f64>)],
) -> String {
    let rule = "############################################################################################\n";
    let mut out = String::from(rule);
    out.push_str("# DragstrParams\n# -------------------------\n");
    for (name, value) in annotations {
        out.push_str(&format!("# {name} = {value}\n"));
    }
    out.push_str(rule);
    let parameters = HyperParameters {
        max_period: rows.len(),
        max_repeat_length,
        ..HyperParameters::default()
    };
    out.push_str(&table(&parameters, rows));
    out
}
