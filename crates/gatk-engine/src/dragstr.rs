//! Ported from `org.broadinstitute.hellbender.utils.dragstr.DragstrParams` and `DragstrParamUtils`
//! (the reader), `org.broadinstitute.hellbender.utils.pairhmm.DragstrReadSTRAnalyzer`, and the
//! haplotypecaller's `DragstrPairHMMInputScoreImputator` and
//! `NonSymmetricalPairHMMInputScoreImputator` (GATK 4.6.2.0): DRAGstr in the read likelihoods.
//!
//! With DRAGstr parameters, the PairHMM's gap penalties come from the read's own sequence rather
//! than from its indel qualities: at every base but the last, the period that repeats most there
//! (the smallest on a tie) and its repeat count look a gap-open penalty (capped at 40) and a gap
//! continuation penalty up in the parameter table, both rounded half up; the last base gets 45 and
//! 10. The table clamps a period or a repeat count past its edge to the last row or column.
//!
//! What a reader would not guess about the STR analyzer:
//!
//!  * a position inside a run is credited with the **whole** run's repeat count, and a position
//!    between two runs with the larger of the two (the period-one pass propagates a longer run back
//!    onto the base before it; the longer periods slide a window of `period + 1` run lengths);
//!  * the window's maximum is kept in a heap whose entries are searched for linearly by value, and
//!    the first run length of the read is lowered by one when it is above one, so a read that starts
//!    inside a repeat does not count the partial unit it starts with.
//!
//! What a reader would not guess about the table file:
//!
//!  * the header is the first line not starting with `#`, and must read `1 2 3 ...`;
//!  * a table ends where the next line ending in `:` starts, so an empty line is an index error;
//!  * a row is split on whitespace **including the empty string a leading space makes**, and that
//!    empty string counts toward the "wrong number of columns" check, so a row with a leading space
//!    and one value too few is accepted and padded with zeros;
//!  * a negative, infinite or NaN value is refused with a `NullPointerException` the reader throws
//!    on purpose.

use std::collections::HashMap;

/// `DragstrHyperParameters.DEFAULT_MAX_PERIOD`.
pub const DEFAULT_MAX_PERIOD: usize = 8;
/// `DragstrHyperParameters.DEFAULT_MAX_REPEAT_LENGTH`.
pub const DEFAULT_MAX_REPEAT_LENGTH: usize = 20;
/// `DragstrPairHMMInputScoreImputator.GOP_AT_THE_END_OF_READ`.
const GOP_AT_THE_END_OF_READ: u8 = 45;
/// `DragstrPairHMMInputScoreImputator.GCP_AT_THE_END_OF_READ`.
const GCP_AT_THE_END_OF_READ: u8 = 10;
/// `DragstrPairHMMInputScoreImputator.MAX_GOP_IN_READ`.
const MAX_GOP_IN_READ: f64 = 40.0;

const DEFAULT_GOP: [[f64; 20]; 8] = [
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

const DEFAULT_API: [[f64; 20]; 8] = [
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

/// What the parameters and their reader refuse, with the exception the reference throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DragstrError {
    /// `UserException.BadInput`; the message without its "Bad input: " prefix.
    BadInput(String),
    /// `IllegalArgumentException` from `ParamUtils.isPositive` and the analyzer's bounds checks.
    IllegalArgument(String),
    /// `IllegalStateException` from `Utils.validate`.
    IllegalState(String),
    /// `NullPointerException`, which the reader throws for a negative, infinite or NaN value.
    NullPointer,
    /// `StringIndexOutOfBoundsException`: an empty line in a table.
    StringIndexOutOfBounds,
    /// `ArrayIndexOutOfBoundsException`: a row shorter than the repeat count, or an empty read.
    IndexOutOfBounds,
}

impl DragstrError {
    /// The simple name of the exception.
    pub fn class(&self) -> &'static str {
        match self {
            DragstrError::BadInput(_) => "BadInput",
            DragstrError::IllegalArgument(_) => "IllegalArgumentException",
            DragstrError::IllegalState(_) => "IllegalStateException",
            DragstrError::NullPointer => "NullPointerException",
            DragstrError::StringIndexOutOfBounds => "StringIndexOutOfBoundsException",
            DragstrError::IndexOutOfBounds => "ArrayIndexOutOfBoundsException",
        }
    }
}

/// `DragstrParams`: GOP, GCP and API by period and repeat count.
#[derive(Debug, Clone, PartialEq)]
pub struct DragstrParams {
    pub name: String,
    max_period: usize,
    max_repeats: usize,
    gop: Vec<Vec<f64>>,
    gcp: Vec<Vec<f64>>,
    api: Vec<Vec<f64>>,
}

/// `Math.round(double)`: the nearest integer, a half rounding up.
fn java_math_round(x: f64) -> i64 {
    if x.is_nan() {
        return 0;
    }
    let floor = x.floor();
    (if x - floor >= 0.5 { floor + 1.0 } else { floor }) as i64
}

impl DragstrParams {
    /// `DragstrParams.DEFAULT`. The GCP is `Math.round(1000.0 / period) / 100.0` for every repeat
    /// count, which keeps two decimals.
    pub fn default_params() -> DragstrParams {
        let gcp = (1..=DEFAULT_MAX_PERIOD)
            .map(|period| {
                let value = java_math_round(1000.0 / period as f64) as f64 / 100.0;
                vec![value; DEFAULT_MAX_REPEAT_LENGTH]
            })
            .collect();
        DragstrParams {
            name: "<default>".to_string(),
            max_period: DEFAULT_MAX_PERIOD,
            max_repeats: DEFAULT_MAX_REPEAT_LENGTH,
            gop: DEFAULT_GOP.iter().map(|r| r.to_vec()).collect(),
            gcp,
            api: DEFAULT_API.iter().map(|r| r.to_vec()).collect(),
        }
    }

    /// `DragstrParams.of(maxPeriod, maxRepeats, gop, gcp, api, name)`.
    pub fn of(
        max_period: usize,
        max_repeats: usize,
        gop: Vec<Vec<f64>>,
        gcp: Vec<Vec<f64>>,
        api: Vec<Vec<f64>>,
        name: &str,
    ) -> Result<DragstrParams, DragstrError> {
        if max_period == 0 {
            return Err(DragstrError::IllegalArgument(
                "max period must be a positive".to_string(),
            ));
        }
        if max_repeats == 0 {
            return Err(DragstrError::IllegalArgument(
                "max repeats must be a positive".to_string(),
            ));
        }
        for (matrix, label) in [(&gop, "gop"), (&gcp, "gcp"), (&api, "api")] {
            if matrix.len() != max_period {
                return Err(DragstrError::IllegalState(format!(
                    "input {label} length must match maxPeriod"
                )));
            }
        }
        for i in 0..max_period {
            for j in 0..max_repeats {
                for (matrix, label) in [(&gop, "gop"), (&gcp, "gcp"), (&api, "api")] {
                    let value = *matrix[i].get(j).ok_or(DragstrError::IndexOutOfBounds)?;
                    if !(value >= 0.0 && value.is_finite()) {
                        return Err(DragstrError::IllegalState(format!(
                            "bad {label} value: {value}"
                        )));
                    }
                }
            }
        }
        for (matrix, label) in [(&gop, "GOP"), (&gcp, "GCP"), (&api, "API")] {
            if matrix.iter().any(|row| row.len() != max_repeats) {
                return Err(DragstrError::BadInput(format!(
                    "the {label} matrix contains rows with length that does not match the max repeat length"
                )));
            }
        }
        Ok(DragstrParams {
            name: name.to_string(),
            max_period,
            max_repeats,
            gop,
            gcp,
            api,
        })
    }

    /// `lookup`: both coordinates must be positive, and are clamped to the table.
    fn lookup(&self, matrix: &[Vec<f64>], period: i32, repeats: i32) -> Result<f64, DragstrError> {
        if period <= 0 {
            return Err(DragstrError::IllegalArgument("period".to_string()));
        }
        if repeats <= 0 {
            return Err(DragstrError::IllegalArgument(
                "repeat length in units".to_string(),
            ));
        }
        let period_index = (period as usize).min(self.max_period) - 1;
        let repeat_index = (repeats as usize).min(self.max_repeats) - 1;
        Ok(matrix[period_index][repeat_index])
    }

    pub fn gop(&self, period: i32, repeats: i32) -> Result<f64, DragstrError> {
        self.lookup(&self.gop, period, repeats)
    }

    pub fn gcp(&self, period: i32, repeats: i32) -> Result<f64, DragstrError> {
        self.lookup(&self.gcp, period, repeats)
    }

    pub fn api(&self, period: i32, repeats: i32) -> Result<f64, DragstrError> {
        self.lookup(&self.api, period, repeats)
    }

    pub fn maximum_period(&self) -> usize {
        self.max_period
    }

    pub fn maximum_repeats(&self) -> usize {
        self.max_repeats
    }
}

/// `String.split("\\s+")`: a leading separator gives a leading empty string, trailing empty strings
/// are dropped.
fn java_split_whitespace(line: &str) -> Vec<&str> {
    let mut parts: Vec<&str> = line
        .split(|c: char| c.is_ascii_whitespace() || c == '\u{b}' || c == '\u{c}')
        .collect();
    // `split` makes an empty string between consecutive separators; `\s+` swallows them, keeping
    // only a leading one.
    let leading = parts.first().is_some_and(|p| p.is_empty()) && !line.is_empty();
    parts.retain(|p| !p.is_empty());
    // A line of separators only is all trailing empties, which Java drops too.
    if leading && !parts.is_empty() {
        parts.insert(0, "");
    }
    parts
}

/// `Double.parseDouble`, including the `d`/`f` suffixes Java accepts.
fn java_parse_double(text: &str) -> Option<f64> {
    let trimmed = text.strip_suffix(['d', 'D', 'f', 'F']).unwrap_or(text);
    match trimmed {
        "NaN" | "+NaN" | "-NaN" => Some(f64::NAN),
        "Infinity" | "+Infinity" => Some(f64::INFINITY),
        "-Infinity" => Some(f64::NEG_INFINITY),
        t if t
            .chars()
            .all(|c| c.is_ascii_digit() || matches!(c, '.' | '-' | '+' | 'e' | 'E')) =>
        {
            t.parse().ok()
        }
        _ => None,
    }
}

/// `linesToMatrix`.
fn lines_to_matrix(lines: &[&str], columns: usize) -> Result<Vec<Vec<f64>>, DragstrError> {
    let mut result = vec![vec![0.0; columns]; lines.len()];
    for (i, line) in lines.iter().enumerate() {
        let parts = java_split_whitespace(line);
        if parts.len() < columns {
            return Err(DragstrError::BadInput(
                "line has the wrong number of columns".to_string(),
            ));
        }
        let mut k = 0;
        for (j, part) in parts.iter().enumerate() {
            if part.is_empty() {
                continue;
            }
            if k >= columns {
                return Err(DragstrError::BadInput(
                    "line has the wrong number of columns".to_string(),
                ));
            }
            let Some(value) = java_parse_double(part) else {
                return Err(DragstrError::BadInput(format!(
                    "score is not a valid Phred value ({},{}) == {part}",
                    i + 1,
                    j + 1
                )));
            };
            if value.is_nan() || value.is_infinite() || value < 0.0 {
                return Err(DragstrError::NullPointer);
            }
            result[i][k] = value;
            k += 1;
        }
    }
    Ok(result)
}

/// `DragstrParamUtils.parse`: the table file `CalibrateDragstrModel` writes, named `name`.
pub fn parse(text: &str, name: &str) -> Result<DragstrParams, DragstrError> {
    let mut lines = text.lines();
    let header = loop {
        match lines.next() {
            None => {
                return Err(DragstrError::BadInput(format!(
                    "there is no content in the dragstr-params file {name}"
                )))
            }
            Some(line) if line.starts_with('#') => continue,
            Some(line) => break line,
        }
    };
    let mut repeats = Vec::new();
    for part in java_split_whitespace(header)
        .into_iter()
        .filter(|p| !p.is_empty())
    {
        repeats.push(
            part.parse::<i64>()
                .map_err(|_| DragstrError::BadInput("bad format for an integer".to_string()))?,
        );
    }
    if repeats.iter().enumerate().any(|(i, &r)| r != i as i64 + 1) {
        let shown: Vec<String> = repeats.iter().map(i64::to_string).collect();
        return Err(DragstrError::BadInput(format!(
            "the DRAGstr parameter file header line must contain integers starting at 1 [{}]",
            shown.join(", ")
        )));
    }
    let max_repeats = repeats.len();
    let Some(first) = lines.next() else {
        return Err(DragstrError::BadInput(
            "end of table list before expected".to_string(),
        ));
    };
    let mut tables: HashMap<String, Vec<Vec<f64>>> = HashMap::new();
    let mut table_name = first.strip_suffix(':').unwrap_or(first).to_string();
    let mut table_lines: Vec<&str> = Vec::new();
    for line in lines {
        match line.chars().last() {
            None => return Err(DragstrError::StringIndexOutOfBounds),
            Some(':') => {
                tables.insert(table_name, lines_to_matrix(&table_lines, max_repeats)?);
                table_name = line.strip_suffix(':').unwrap_or(line).to_string();
                table_lines.clear();
            }
            Some(_) => table_lines.push(line),
        }
    }
    if table_name.is_empty() {
        return Err(DragstrError::BadInput("table with no name".to_string()));
    }
    tables.insert(table_name, lines_to_matrix(&table_lines, max_repeats)?);
    let mut take = |key: &str| {
        tables
            .remove(key)
            .ok_or_else(|| DragstrError::BadInput(format!("missing matrix {key}")))
    };
    let gop = take("GOP")?;
    let gcp = take("GCP")?;
    let api = take("API")?;
    let max_period = gop.len();
    DragstrParams::of(max_period, max_repeats, gop, gcp, api, name)
}

/// `DragstrReadSTRAnalyzer`: per position, the repeats of each period up to the maximum and the
/// period with the most.
#[derive(Debug, Clone)]
pub struct DragstrReadStrAnalyzer {
    repeats_by_period_and_position: Vec<Vec<i32>>,
    period_with_most_repeats: Vec<i32>,
    max_period: usize,
    seq_length: usize,
}

impl DragstrReadStrAnalyzer {
    /// `of(bases, maxPeriod)`.
    pub fn of(bases: &[u8], max_period: usize) -> Result<Self, DragstrError> {
        if max_period == 0 {
            return Err(DragstrError::IllegalArgument(
                "the input max period must be 1 or greater".to_string(),
            ));
        }
        if bases.is_empty() {
            // `calculateRepeatsForPeriodOne` reads `bases[-1]`.
            return Err(DragstrError::IndexOutOfBounds);
        }
        let n = bases.len();
        let mut repeats = vec![vec![0i32; n]; max_period];
        repeats_for_period_one(bases, &mut repeats[0]);
        let mut run_length_buffer = vec![0i32; n + 1];
        let mut heap = vec![0i32; max_period + 1];
        for period in 2..=max_period {
            repeats_for_period_two_and_above(
                bases,
                period,
                &mut run_length_buffer,
                &mut heap,
                &mut repeats[period - 1],
            );
        }
        let mut period_with_most = vec![1i32; n];
        let mut most = repeats[0].clone();
        for (period_index, values) in repeats.iter().enumerate().skip(1) {
            for position in 0..n {
                let value = values[position];
                if value > most[position] {
                    most[position] = value;
                    period_with_most[position] = period_index as i32 + 1;
                }
            }
        }
        Ok(DragstrReadStrAnalyzer {
            repeats_by_period_and_position: repeats,
            period_with_most_repeats: period_with_most,
            max_period,
            seq_length: n,
        })
    }

    /// `numberOfRepeats(position, period)`: 0 for a period outside `1..=maxPeriod`.
    pub fn number_of_repeats(&self, position: usize, period: i32) -> Result<i32, DragstrError> {
        if period <= 0 || period as usize > self.max_period {
            return Ok(0);
        }
        if position >= self.seq_length {
            return Err(DragstrError::IllegalArgument(
                "cannot query outside requested boundaries".to_string(),
            ));
        }
        Ok(self.repeats_by_period_and_position[period as usize - 1][position])
    }

    /// `mostRepeatedPeriod(position)`.
    pub fn most_repeated_period(&self, position: usize) -> Result<i32, DragstrError> {
        self.period_with_most_repeats
            .get(position)
            .copied()
            .ok_or_else(|| {
                DragstrError::IllegalArgument(format!(
                    "cannot query outside requested boundaries [0, {}): {position}",
                    self.seq_length
                ))
            })
    }

    /// `numberOfMostRepeats(position)`.
    pub fn number_of_most_repeats(&self, position: usize) -> Result<i32, DragstrError> {
        let period = self.most_repeated_period(position).map_err(|_| {
            DragstrError::IllegalArgument("cannot query outside requested boundaries".to_string())
        })?;
        Ok(self.repeats_by_period_and_position[period as usize - 1][position])
    }
}

/// `calculateRepeatsForPeriodOne`: each position gets its homopolymer run's length, and the base
/// before a run gets the run's length if that is longer than its own.
fn repeats_for_period_one(bases: &[u8], output: &mut [i32]) {
    let right_margin = bases.len() - 1;
    let mut last = bases[right_margin];
    let mut carry_back = 1;
    output[right_margin] = 1;
    for position in (0..right_margin).rev() {
        let next = bases[position];
        if next == last {
            carry_back += 1;
        } else {
            carry_back = 1;
        }
        output[position] = carry_back;
        last = next;
    }
    let mut prev_run_length = output[0];
    for position in 1..=right_margin {
        let next = bases[position];
        if next == last {
            output[position] = prev_run_length;
        } else {
            let this_run_length = output[position];
            if prev_run_length < this_run_length {
                output[position - 1] = this_run_length;
            }
            last = next;
            prev_run_length = this_run_length;
        }
    }
}

/// `calculateRepeatsForPeriodTwoAndAbove`.
fn repeats_for_period_two_and_above(
    bases: &[u8],
    period: usize,
    run_length_buffer: &mut [i32],
    heap: &mut [i32],
    output: &mut [i32],
) {
    let seq_length = bases.len();
    if seq_length < period {
        output[..seq_length].fill(0);
        return;
    }
    // The last `period - 1` positions cannot start a repeat.
    let mut position = seq_length as i64 - 1;
    let mut cycle_index = period;
    while cycle_index > 1 {
        run_length_buffer[position as usize] = 0;
        position -= 1;
        cycle_index -= 1;
    }
    let mut carry_back = 1;
    run_length_buffer[position as usize] = 1;
    position -= 1;
    let mut matched_cycles = 0;
    while position >= 0 {
        let p = position as usize;
        if bases[p] == bases[p + period] {
            matched_cycles += 1;
            if matched_cycles == period {
                carry_back += 1;
                run_length_buffer[p] = carry_back;
                matched_cycles = 0;
            } else {
                run_length_buffer[p] = carry_back;
            }
        } else {
            carry_back = 1;
            run_length_buffer[p] = 1;
            matched_cycles = 0;
        }
        position -= 1;
    }
    // Every position of a run gets the run's total.
    for cycle in 0..period {
        let mut p = cycle;
        while p < seq_length {
            let total = run_length_buffer[p];
            for _ in 1..total {
                p += period;
                run_length_buffer[p] = total;
            }
            p += period;
        }
    }
    if run_length_buffer[0] > 1 {
        run_length_buffer[0] -= 1;
    }
    let heap_size = period + 1;
    heap[..heap_size].fill(0);
    heap[0] = run_length_buffer[0];
    let mut current_max = heap[0];
    let stop0 = period.min(seq_length - 1);
    let mut position = 0usize;
    while position < stop0 {
        // `output[position] = currentMax = Math.max(currentMax, heap[++position] = runLengthBuffer[position])`:
        // the output index is read before the increment, the buffer index after it.
        let out = position;
        position += 1;
        heap[position] = run_length_buffer[position];
        current_max = current_max.max(heap[position]);
        output[out] = current_max;
        fix_heap(heap, position, heap_size);
    }
    run_length_buffer[seq_length] = 1;
    let mut out_position = 0usize;
    while position < seq_length {
        let value_out = run_length_buffer[out_position];
        out_position += 1;
        position += 1;
        let value_in = run_length_buffer[position];
        if value_in != value_out {
            // `ArrayUtils.indexOf` over the whole array.
            let index = heap
                .iter()
                .position(|&v| v == value_out)
                .expect("the value leaving is in the heap");
            heap[index] = value_in;
            fix_heap(heap, index, heap_size);
            current_max = heap[0];
        }
        output[position - 1] = current_max;
    }
}

fn fix_heap(heap: &mut [i32], index: usize, heap_size: usize) {
    if index == 0 || heap[(index - 1) >> 1] > heap[index] {
        fix_heap_down(heap, index, heap_size);
    } else {
        fix_heap_up(heap, index);
    }
}

fn fix_heap_up(heap: &mut [i32], mut index: usize) {
    let value = heap[index];
    loop {
        let up = (index - 1) >> 1;
        let up_value = heap[up];
        if up_value >= value {
            break;
        }
        heap[index] = up_value;
        index = up;
        if index == 0 {
            break;
        }
    }
    heap[index] = value;
}

fn fix_heap_down(heap: &mut [i32], mut index: usize, heap_size: usize) {
    let value = heap[index];
    loop {
        let right = (index + 1) << 1;
        let right_value = if right < heap_size { heap[right] } else { -1 };
        let left = right - 1;
        let left_value = if left < heap_size { heap[left] } else { -1 };
        if right_value > value {
            if right_value > left_value {
                heap[index] = right_value;
                index = right;
            } else {
                heap[index] = left_value;
                index = left;
            }
        } else if left_value > value {
            heap[index] = left_value;
            index = left;
        } else {
            break;
        }
    }
    heap[index] = value;
}

/// `DragstrPairHMMInputScoreImputator.impute`: the gap-open penalties (the same for insertions and
/// deletions) and the gap continuation penalties of a read.
pub fn dragstr_impute(
    params: &DragstrParams,
    bases: &[u8],
) -> Result<(Vec<u8>, Vec<u8>), DragstrError> {
    let length = bases.len();
    let analyzer = DragstrReadStrAnalyzer::of(bases, params.maximum_period())?;
    let mut gop = vec![0u8; length];
    let mut gcp = vec![0u8; length];
    for i in 0..length - 1 {
        let period = analyzer.most_repeated_period(i)?;
        let repeats = analyzer.number_of_most_repeats(i)?;
        gop[i] = java_math_round(MAX_GOP_IN_READ.min(params.gop(period, repeats)?)) as u8;
        gcp[i] = java_math_round(params.gcp(period, repeats)?) as u8;
    }
    gop[length - 1] = GOP_AT_THE_END_OF_READ;
    gcp[length - 1] = GCP_AT_THE_END_OF_READ;
    Ok((gop, gcp))
}

/// `NonSymmetricalPairHMMInputScoreImputator.impute`: flat deletion and insertion penalties over the
/// read's base qualities, and a flat continuation penalty over its bases, as
/// `(deletion, insertion, continuation)`.
pub fn non_symmetrical_impute(
    constant_gcp: u8,
    flat_insertion_qual: u8,
    flat_deletion_qual: u8,
    base_quality_count: usize,
    read_length: usize,
) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    (
        vec![flat_deletion_qual; base_quality_count],
        vec![flat_insertion_qual; base_quality_count],
        vec![constant_gcp; read_length],
    )
}
