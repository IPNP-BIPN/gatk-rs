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
