//! `DepthOfCoverage`: how deep the reads lie over a reference.
//!
//! The tool writes a family of files rather than one, and which of them appear is decided by four
//! `omit` arguments rather than by the data. The partition is the SAMPLE unless `--partition-type`
//! says otherwise, so by default read groups collapse
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
    /// The last reference position the read covers, for an all-`M` read. The runner counts off the
    /// engine's pileup instead, which reads the cigar.
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

/// The columns `--print-base-counts` breaks a depth into, in the order it writes them:
/// `BaseUtils.BASES_EXTENDED`, whose `D` is written only under `--include-deletions`.
pub const BASE_COUNT_ORDER: &[u8] = b"ACGTND";

/// One locus of the per-base table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub contig: String,
    pub position: i32,
    /// One entry per sample, in the order the samples are listed.
    pub depths: Vec<i32>,
    /// One `A: C: G: T: N: D:` breakdown per sample.
    pub base_counts: Vec<[i32; 6]>,
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
            let mut counts = vec![[0; 6]; samples.len()];
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
                counts[index][base_index(base)] += 1;
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

/// Where `CoverageUtils.updateCounts` counts a read base: anything that is not one of the four
/// bases is counted as N.
fn base_index(base: u8) -> usize {
    BASE_COUNT_ORDER[..4]
        .iter()
        .position(|expected| *expected == base.to_ascii_uppercase())
        .unwrap_or(4)
}

/// `DoCOutputType.Partition`, in its declaration order, which is the order an `EnumSet` of them
/// iterates in and so the order of the partitions' columns.
pub const PARTITIONS: [&str; 8] = [
    "readgroup",
    "sample",
    "library",
    "platform",
    "center",
    "sample_by_platform",
    "sample_by_center",
    "sample_by_platform_by_center",
];

/// `CoverageUtils.getTypeID`: the column a read group counts into under one partition. A field
/// the group lacks is Java's `null`: the word, inside a concatenation, and `None` on its own.
pub fn type_id(
    partition: &str,
    id: &str,
    field: &dyn Fn(&str) -> Option<String>,
) -> Option<String> {
    let get = |key: &str| field(key).unwrap_or_else(|| "null".to_string());
    Some(match partition {
        "readgroup" => format!("{}_rg_{id}", get("SM")),
        "library" => return field("LB"),
        "platform" => return field("PL"),
        "center" => return field("CN"),
        "sample_by_platform" => format!("{}_pl_{}", get("SM"), get("PL")),
        "sample_by_center" => format!("{}_cn_{}", get("SM"), get("CN")),
        "sample_by_platform_by_center" => {
            format!("{}_pl_{}_cn_{}", get("SM"), get("PL"), get("CN"))
        }
        _ => return field("SM"),
    })
}

/// `makeGlobalIdentifierMap` for one partition: the distinct ids, sorted by `String::compareTo`.
///
/// A `null` id is first in its `HashSet`, so with any other id beside it the sort hands it to
/// `compareTo` as the argument and throws; alone it is never compared, and its column is named
/// `null`.
pub fn sorted_ids(mut ids: Vec<Option<String>>) -> Result<Vec<String>, &'static str> {
    ids.sort();
    ids.dedup();
    if ids.len() > 1 && ids[0].is_none() {
        return Err("Cannot read field \"value\" because \"anotherString\" is null");
    }
    Ok(ids
        .into_iter()
        .map(|id| id.unwrap_or_else(|| "null".to_string()))
        .collect())
}

/// One pileup element as the counter reads it: the columns its read group counts into, one per
/// partition, its base, quality and whether it is a deletion.
pub struct Element {
    pub columns: Vec<usize>,
    pub base: u8,
    pub quality: i32,
    pub deletion: bool,
}

/// `CoverageUtils.getBaseCountsByPartition` over one locus's pileup, into `width` columns.
///
/// A deletion is in the pileup only under `--include-deletions`, and then it counts whatever its
/// quality, into `D` and into the depth.
pub fn count_pileup(
    contig: &str,
    position: i32,
    width: usize,
    elements: impl IntoIterator<Item = Element>,
    minimum_base_quality: i32,
    maximum_base_quality: i32,
) -> Locus {
    let mut depths = vec![0; width];
    let mut counts = vec![[0; 6]; width];
    for element in elements {
        if !element.deletion
            && !base_counts(element.quality, minimum_base_quality, maximum_base_quality)
        {
            continue;
        }
        let at = if element.deletion {
            5
        } else {
            base_index(element.base)
        };
        for column in element.columns {
            depths[column] += 1;
            counts[column][at] += 1;
        }
    }
    Locus {
        contig: contig.to_string(),
        position,
        depths,
        base_counts: counts,
    }
}

/// The header of the per-locus table over partitions, each with its sorted identifiers.
pub fn per_locus_header_for(
    partitions: &[(&str, Vec<String>)],
    print_base_counts: bool,
    separator: char,
) -> String {
    let mut out = format!("Locus{separator}Total_Depth");
    for (partition, _) in partitions {
        out.push_str(&format!("{separator}Average_Depth_{partition}"));
    }
    for (_, ids) in partitions {
        for id in ids {
            out.push_str(&format!("{separator}Depth_for_{id}"));
            if print_base_counts {
                out.push_str(&format!("{separator}{id}_base_counts"));
            }
        }
    }
    out
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
    per_locus_row_for(locus, &[locus.depths.len()], print_base_counts, false, ',')
}

/// `writePerLocusDepthSummary` over partitions of these widths, whose columns `locus` holds one
/// after the other. The total is the first partition's, and every partition's average divides it.
pub fn per_locus_row_for(
    locus: &Locus,
    widths: &[usize],
    print_base_counts: bool,
    deletions: bool,
    separator: char,
) -> String {
    let total: i32 = locus.depths[..widths.first().copied().unwrap_or(0)]
        .iter()
        .sum();
    let mut out = format!("{}:{}{separator}{total}", locus.contig, locus.position);
    for width in widths {
        // `DOUBLE_FORMAT_2PLACES`, HALF_UP.
        let average = if *width == 0 {
            0.0
        } else {
            f64::from(total) / *width as f64
        };
        out.push_str(&format!(
            "{separator}{}",
            gatk_engine::java_format::format_decimals(average, 2)
        ));
    }
    for (index, depth) in locus.depths.iter().enumerate() {
        out.push_str(&format!("{separator}{depth}"));
        if print_base_counts {
            let counts = locus.base_counts[index];
            // `getBaseCountsString`: a space after every pair but the sixth, so without the `D`
            // the field ends with one.
            let mut text = String::new();
            for (at, base) in BASE_COUNT_ORDER.iter().enumerate() {
                if at == 5 && !deletions {
                    continue;
                }
                text.push_str(&format!("{}:{}", *base as char, counts[at]));
                if at < 5 {
                    text.push(' ');
                }
            }
            out.push_str(&format!("{separator}{text}"));
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
