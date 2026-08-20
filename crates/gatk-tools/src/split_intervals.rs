//! `SplitIntervals`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.SplitIntervals` (GATK 4.6.2.0).
//!
//! The shards a scatter-gather run is given: an interval list divided by one of the five modes
//! [`gatk_engine::interval_list_scatter`] pins, written one file per shard.
//!
//! # The file name's width comes from a logarithm of one less than the count
//!
//! ```java
//! final int maxNumberOfPlaces = Math.max((int)Math.floor(Math.log10(scatterCount-1))+1, numDigits);
//! ```
//!
//! At a scatter count of one that is `log10(0)`, which is negative infinity. The cast floors it to
//! `Integer.MIN_VALUE` and the `+ 1` leaves it enormously negative, so the `max` with
//! `--interval-file-num-digits` is the only thing that saves the name. The default is four, which
//! is why every ordinary run writes `0000-scattered.interval_list` whatever it was asked for, and
//! only a scatter count above ten thousand widens the name on its own.
//!
//! The addition is written here as a saturating one. In Java it overflows to `Integer.MIN_VALUE + 1`
//! and stays negative, which is the only property the `max` needs; in Rust the same addition would
//! panic in a debug build, and a port that panicked where the reference wrote a file would be
//! wrong in the loudest possible way.
//!
//! # The contig filter is only for the whole reference
//!
//! `--min-contig-size` is applied to the contigs the tool makes intervals from when no `-L` was
//! given. With any `-L` at all the argument is read and never used.
//!
//! # The contig split happens after the scatter
//!
//! `--dont-mix-contigs` regroups each shard by contig, so the requested count is a lower bound on
//! the number of files. The sublists come out ordered by the contig's index in the sequence
//! dictionary rather than by the order the intervals arrived in.

use gatk_engine::interval_list_scatter::{scatter, ScatterError, ScatterMode};
use htsjdk_bam::interval::{Interval, IntervalList};

/// `DEFAULT_EXTENSION`, `DEFAULT_PREFIX` and `DEFAULT_NUMBER_OF_DIGITS`.
pub const DEFAULT_EXTENSION: &str = "-scattered.interval_list";
pub const DEFAULT_PREFIX: &str = "";
pub const DEFAULT_NUMBER_OF_DIGITS: i32 = 4;

/// What the tool refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitError {
    /// `ParamUtils.isPositive(scatterCount, ...)`.
    ScatterCountNotPositive,
    /// The parser's bound on the number of digits, which is not the tool's code.
    DigitsOutOfRange(i32),
    /// The scatterer underneath, which this tool's own check makes unreachable.
    Scatter(ScatterError),
}

impl SplitError {
    pub fn java_class(&self) -> &'static str {
        match self {
            SplitError::DigitsOutOfRange(_) => {
                "org.broadinstitute.barclay.argparser.CommandLineException$OutOfRangeArgumentValue"
            }
            _ => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            SplitError::ScatterCountNotPositive => "scatter-count must be > 0.".to_string(),
            SplitError::DigitsOutOfRange(value) => format!(
                "Argument interval-file-num-digits has a bad value: {value}. minimum allowed value 1"
            ),
            SplitError::Scatter(error) => error.message(),
        }
    }
}

/// The tool's arguments, with the reference's defaults.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub scatter_count: i32,
    pub min_contig_size: i32,
    pub subdivision_mode: ScatterMode,
    pub prefix: String,
    pub extension: String,
    pub num_digits: i32,
    pub dont_mix_contigs: bool,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            scatter_count: 1,
            min_contig_size: 0,
            subdivision_mode: ScatterMode::IntervalSubdivision,
            prefix: DEFAULT_PREFIX.to_string(),
            extension: DEFAULT_EXTENSION.to_string(),
            num_digits: DEFAULT_NUMBER_OF_DIGITS,
            dont_mix_contigs: false,
        }
    }
}

/// `maxNumberOfPlaces`: the width of the number in each file's name.
///
/// `(int) Math.floor(Math.log10(n - 1)) + 1`, where a scatter count of one makes the logarithm
/// negative infinity and the cast `Integer.MIN_VALUE`. The `max` is what decides in that case, and
/// in every case below ten thousand with the default four digits.
pub fn max_number_of_places(scatter_count: i32, num_digits: i32) -> i32 {
    let from_count = if scatter_count <= 1 {
        // `(int) Math.floor(-Infinity)` is `Integer.MIN_VALUE`, and `+ 1` in Java wraps no
        // further than `MIN_VALUE + 1`. Saturating here keeps the sign, which is all the max
        // reads, without a debug-build panic.
        i32::MIN.saturating_add(1)
    } else {
        (f64::from(scatter_count - 1).log10().floor() as i32).saturating_add(1)
    };
    from_count.max(num_digits)
}

/// One shard's file name: the prefix, the zero-padded index and the extension, concatenated raw.
pub fn file_name(index: usize, arguments: &Arguments) -> String {
    let places =
        max_number_of_places(arguments.scatter_count, arguments.num_digits).max(0) as usize;
    format!(
        "{}{:0width$}{}",
        arguments.prefix,
        index,
        arguments.extension,
        width = places
    )
}

/// `getAllIntervalsForReference` filtered by `--min-contig-size`, which is the input when no `-L`
/// was given.
pub fn intervals_for_reference(sequences: &[(String, i32)], min_contig_size: i32) -> Vec<Interval> {
    sequences
        .iter()
        .filter(|(_, length)| *length >= min_contig_size)
        .map(|(name, length)| Interval::new(name, 1, *length))
        .collect()
}

/// `--dont-mix-contigs`: each shard regrouped by contig, in dictionary order.
pub fn split_by_contig(shard: &IntervalList) -> Vec<IntervalList> {
    let mut out: Vec<IntervalList> = Vec::new();
    for contig in &shard.dictionary {
        let intervals: Vec<Interval> = shard
            .intervals
            .iter()
            .filter(|interval| &interval.contig == contig)
            .cloned()
            .collect();
        if intervals.is_empty() {
            continue;
        }
        let mut list = IntervalList::new(shard.dictionary.clone());
        list.intervals = intervals;
        out.push(list);
    }
    out
}

/// The whole run: the shards, each with the name it is written under.
///
/// `intervals` is `None` for a run with no `-L`, which is the whole reference filtered by
/// `--min-contig-size`. The intervals are expected already merged by the interval argument
/// collection, whose default rule is ALL: this tool does no merging of its own.
pub fn split(
    intervals: Option<&[Interval]>,
    sequences: &[(String, i32)],
    arguments: &Arguments,
) -> Result<Vec<(String, IntervalList)>, SplitError> {
    if arguments.scatter_count <= 0 {
        return Err(SplitError::ScatterCountNotPositive);
    }
    if arguments.num_digits < 1 {
        return Err(SplitError::DigitsOutOfRange(arguments.num_digits));
    }
    let dictionary: Vec<String> = sequences.iter().map(|(name, _)| name.clone()).collect();
    let mut list = IntervalList::new(dictionary);
    list.intervals = match intervals {
        Some(given) => given.to_vec(),
        None => intervals_for_reference(sequences, arguments.min_contig_size),
    };
    let scattered = scatter(&list, arguments.subdivision_mode, arguments.scatter_count)
        .map_err(SplitError::Scatter)?;
    let final_shards: Vec<IntervalList> = if arguments.dont_mix_contigs {
        scattered.iter().flat_map(split_by_contig).collect()
    } else {
        scattered
    };
    Ok(final_shards
        .into_iter()
        .enumerate()
        .map(|(index, shard)| (file_name(index, arguments), shard))
        .collect())
}
