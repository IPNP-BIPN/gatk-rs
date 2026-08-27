//! `DepthOfCoverage`: how deep the reads lie over a reference.
//!
//! The tool writes a family of files rather than one, and which of them appear is decided by four
//! `omit` arguments rather than by the data. The partition is the SAMPLE, so read groups collapse
//! into one column each.
//!
//! Reading the BAM and the intervals is not ported. Which bases are counted, how they are counted
//! and which files a set of arguments produces are.

/// One read, reduced to what the counter reads off it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    /// The read group's SAMPLE, which is what the counter partitions by.
    pub sample: String,
    pub contig: String,
    pub start: i32,
    pub bases: Vec<u8>,
    pub base_qualities: Vec<i32>,
}

impl Read {
    /// The last reference position the read covers. Only `M` cigars are ported.
    pub fn end(&self) -> i32 {
        self.start + self.bases.len() as i32 - 1
    }

    /// The base and its quality at one reference position, if the read covers it.
    pub fn at(&self, position: i32) -> Option<(u8, i32)> {
        if position < self.start || position > self.end() {
            return None;
        }
        let offset = (position - self.start) as usize;
        Some((self.bases[offset], self.base_qualities[offset]))
    }
}

/// `MIN_BASE_QUALITY` and `MAX_BASE_QUALITY`, and the byte range the argument itself enforces.
pub const DEFAULT_MIN_BASE_QUALITY: i32 = 0;
pub const DEFAULT_MAX_BASE_QUALITY: i32 = 127;
pub const BASE_QUALITY_MIN: i32 = 0;
pub const BASE_QUALITY_MAX: i32 = 127;

/// A base counts when its quality is within BOTH bounds: the ceiling is a filter no other coverage
/// tool has.
pub fn base_counts(quality: i32, minimum: i32, maximum: i32) -> bool {
    quality >= minimum && quality <= maximum
}

/// The five columns `--print-base-counts` breaks a depth into, in the order it writes them.
pub const BASE_COUNT_ORDER: &[u8] = b"ACGTN";

/// One locus of the per-base table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub contig: String,
    pub position: i32,
    /// One entry per sample, in the order the samples are listed.
    pub depths: Vec<i32>,
    /// One `A: C: G: T: N:` breakdown per sample, when `--print-base-counts` was given.
    pub base_counts: Vec<[i32; 5]>,
}

impl Locus {
    pub fn total_depth(&self) -> i32 {
        self.depths.iter().sum()
    }

    /// `Average_Depth_sample`, which divides by the number of SAMPLES rather than by the number
    /// that carried anything.
    pub fn average_depth(&self) -> f64 {
        if self.depths.is_empty() {
            return 0.0;
        }
        self.total_depth() as f64 / self.depths.len() as f64
    }
}

/// The per-locus table over an interval.
///
/// EVERY base of the interval is a row, including those no read reaches, so the table's length is
/// the interval's and not the reads'.
pub fn per_locus(
    reads: &[Read],
    samples: &[String],
    contig: &str,
    start: i32,
    end: i32,
    minimum_base_quality: i32,
    maximum_base_quality: i32,
) -> Vec<Locus> {
    (start..=end)
        .map(|position| {
            let mut depths = vec![0; samples.len()];
            let mut counts = vec![[0; 5]; samples.len()];
            for read in reads {
                if read.contig != contig {
                    continue;
                }
                let Some(index) = samples.iter().position(|sample| *sample == read.sample) else {
                    continue;
                };
                let Some((base, quality)) = read.at(position) else {
                    continue;
                };
                if !base_counts(quality, minimum_base_quality, maximum_base_quality) {
                    continue;
                }
                depths[index] += 1;
                let at = BASE_COUNT_ORDER
                    .iter()
                    .position(|expected| *expected == base.to_ascii_uppercase())
                    // Anything that is not one of the four bases is counted as N.
                    .unwrap_or(4);
                counts[index][at] += 1;
            }
            Locus {
                contig: contig.to_string(),
                position,
                depths,
                base_counts: counts,
            }
        })
        .collect()
}

/// The header of the per-locus table, whose sample columns depend on `--print-base-counts`.
pub fn per_locus_header(samples: &[String], print_base_counts: bool) -> String {
    let mut out = String::from("Locus,Total_Depth,Average_Depth_sample");
    for sample in samples {
        out.push_str(&format!(",Depth_for_{sample}"));
        if print_base_counts {
            out.push_str(&format!(",{sample}_base_counts"));
        }
    }
    out
}

/// One row of the per-locus table, rendered as the writer renders it.
pub fn per_locus_row(locus: &Locus, print_base_counts: bool) -> String {
    let mut out = format!(
        "{}:{},{},{:.2}",
        locus.contig,
        locus.position,
        locus.total_depth(),
        locus.average_depth()
    );
    for (index, depth) in locus.depths.iter().enumerate() {
        out.push_str(&format!(",{depth}"));
        if print_base_counts {
            let counts = locus.base_counts[index];
            let mut text = String::new();
            for (at, base) in BASE_COUNT_ORDER.iter().enumerate() {
                // Each pair is followed by a space, so the field ends with one.
                text.push_str(&format!("{}:{} ", *base as char, counts[at]));
            }
            out.push_str(&format!(",{text}"));
        }
    }
    out
}

/// The four `omit` arguments, each of which removes its own files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Omissions {
    pub locus_table: bool,
    pub depth_output_at_each_base: bool,
    pub per_sample_statistics: bool,
    pub interval_statistics: bool,
}

/// The suffixes a run writes, sorted the way the dump sorts them.
///
/// Which files appear is a function of the ARGUMENTS alone: no omission depends on whether the
/// data would have filled the file.
pub fn written_suffixes(omissions: &Omissions) -> Vec<&'static str> {
    let mut out: Vec<&'static str> = Vec::new();
    if !omissions.depth_output_at_each_base {
        // The base file itself, whose suffix is empty.
        out.push("");
    }
    if !omissions.locus_table {
        out.push(".sample_cumulative_coverage_counts");
        out.push(".sample_cumulative_coverage_proportions");
    }
    if !omissions.interval_statistics {
        out.push(".sample_interval_statistics");
        out.push(".sample_interval_summary");
    }
    if !omissions.per_sample_statistics {
        out.push(".sample_statistics");
        out.push(".sample_summary");
    }
    out.sort();
    out
}

/// What the tool and its argument parser refuse.
#[derive(Debug, Clone, PartialEq)]
pub enum CoverageError {
    /// The two interval arguments, which are declared mutually exclusive.
    MutuallyExclusive { argument: String, other: String },
    /// A value that will not fit in a `Byte` at all.
    BadByte { argument: String, value: String },
    /// A value that fits but is outside the declared range.
    OutOfRange {
        argument: String,
        value: String,
        minimum: f64,
        maximum: f64,
    },
}

impl CoverageError {
    pub fn message(&self) -> String {
        match self {
            CoverageError::MutuallyExclusive { argument, other } => format!(
                "Argument '{argument}' cannot be used in conjunction with argument(s) {other}"
            ),
            CoverageError::BadByte { argument, value } => format!(
                "Argument {argument} has a bad value: {value}. Failure constructing 'Byte' from \
                 the string '{value}'."
            ),
            CoverageError::OutOfRange {
                argument,
                value,
                minimum,
                maximum,
            } => format!(
                "Argument {argument} has a bad value: {value}. allowed range [{minimum:.1}, \
                 {maximum:.1}]."
            ),
        }
    }
}

/// The argument parser's own check, which runs before a read is seen. A value outside the byte
/// range fails one of TWO ways depending on which side it is out on.
pub fn check_base_quality(argument: &str, value: i64) -> Result<(), CoverageError> {
    if value > i8::MAX as i64 || value < i8::MIN as i64 {
        return Err(CoverageError::BadByte {
            argument: argument.to_string(),
            value: value.to_string(),
        });
    }
    if !(BASE_QUALITY_MIN as i64..=BASE_QUALITY_MAX as i64).contains(&value) {
        return Err(CoverageError::OutOfRange {
            argument: argument.to_string(),
            value: value.to_string(),
            minimum: BASE_QUALITY_MIN as f64,
            maximum: BASE_QUALITY_MAX as f64,
        });
    }
    Ok(())
}
