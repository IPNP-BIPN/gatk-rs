//! `CondenseDepthEvidence`, ported from
//! `org.broadinstitute.hellbender.tools.sv.CondenseDepthEvidence` and `DepthEvidenceCodec`
//! (GATK 4.6.2.0).
//!
//! Adjacent depth-evidence bins merged.
//!
//! # The maximum is not a maximum
//!
//! ```java
//! final int intervalLength = accumulator.getLengthOnReference();
//! if ( !isAdjacent(accumulator, feature) || intervalLength >= maxIntervalLength ) { ... }
//! ```
//!
//! The length tested is the one ALREADY accumulated, before the next bin is added, so the check
//! fires one bin late and the interval that is written is longer than the limit. With hundred-base
//! bins, a maximum of 150 and a maximum of 200 produce the same file: intervals of exactly 200.
//!
//! # The minimum drops records rather than merging them
//!
//! A run shorter than the minimum is not written at all, and the same test is applied again to the
//! last accumulator, so a trailing short interval disappears too. Nothing says so in the output.
//!
//! # The file is zero-based on disk and one-based inside
//!
//! `DepthEvidenceCodec.decode` adds one to the start it reads, and `encode` subtracts it again, so
//! a bin written `0 100` is the closed interval 1..=100 and is a hundred bases long.

/// One bin, as the tool holds it: one-based and closed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DepthEvidence {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub counts: Vec<i32>,
}

impl DepthEvidence {
    /// `getLengthOnReference()`.
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
}

/// The arguments, with the tool's own defaults.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Arguments {
    pub max_interval_length: i32,
    pub min_interval_length: i32,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            max_interval_length: 1000,
            min_interval_length: 0,
        }
    }
}

/// What the run refuses, both before a record is read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CondenseError {
    MinimumAboveMaximum,
    /// No codec answers for the output's name, which `FeatureOutputCodecFinder.find` refuses.
    NoOutputCodec {
        path: String,
    },
    /// The output's extension implies another feature type, which the message names.
    WrongOutputType {
        path: String,
        found: String,
    },
    /// `MathUtils.addToArrayInPlace`, over two adjacent bins with different numbers of counts.
    CountMismatch,
}

impl CondenseError {
    pub fn java_class(&self) -> &'static str {
        match self {
            CondenseError::CountMismatch => "java.lang.IllegalArgumentException",
            _ => "org.broadinstitute.hellbender.exceptions.UserException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            CondenseError::MinimumAboveMaximum => {
                "Minimum interval length exceeds maximum interval length.".to_string()
            }
            CondenseError::NoOutputCodec { path } => {
                crate::sv_feature_codecs::no_output_codec(path)
            }
            CondenseError::WrongOutputType { path, found } => format!(
                "Output file {path} implies Feature subtype {found}, but this tool expects to \
                 write DepthEvidence."
            ),
            CondenseError::CountMismatch => "Arrays must have same length".to_string(),
        }
    }
}

/// `onTraversalStart`'s first check.
pub fn check_lengths(arguments: &Arguments) -> Result<(), CondenseError> {
    if arguments.min_interval_length > arguments.max_interval_length {
        return Err(CondenseError::MinimumAboveMaximum);
    }
    Ok(())
}

/// `onTraversalStart`'s second check, which is the finder's refusal and then the tool's.
pub fn check_output(path: &str) -> Result<(), CondenseError> {
    match crate::sv_feature_codecs::find(path) {
        Some(codec) if codec.feature_type == "DepthEvidence" => Ok(()),
        Some(codec) => Err(CondenseError::WrongOutputType {
            path: path.to_string(),
            found: codec.feature_type.to_string(),
        }),
        None => Err(CondenseError::NoOutputCodec {
            path: path.to_string(),
        }),
    }
}

/// `isAdjacent`.
fn adjacent(left: &DepthEvidence, right: &DepthEvidence) -> bool {
    left.contig == right.contig && left.end + 1 == right.start
}

/// `apply` over every record, then `onTraversalSuccess`: the records the sink was handed.
///
/// The counts are summed as Java ints, so a sum past `i32::MAX` wraps, and two adjacent bins with
/// different numbers of counts are refused where they meet rather than widened.
pub fn condense(
    records: &[DepthEvidence],
    arguments: &Arguments,
) -> Result<Vec<DepthEvidence>, CondenseError> {
    let mut written = Vec::new();
    let mut accumulator: Option<DepthEvidence> = None;
    for record in records {
        let held = match accumulator {
            None => {
                accumulator = Some(record.clone());
                continue;
            }
            Some(ref held) => held.clone(),
        };
        // The length of what is already held, tested BEFORE this record joins it.
        let length = held.length();
        if !adjacent(&held, record) || length >= arguments.max_interval_length {
            if length >= arguments.min_interval_length {
                written.push(held);
            }
            accumulator = Some(record.clone());
            continue;
        }
        if held.counts.len() != record.counts.len() {
            return Err(CondenseError::CountMismatch);
        }
        let counts = held
            .counts
            .iter()
            .zip(&record.counts)
            .map(|(left, right)| left.wrapping_add(*right))
            .collect();
        accumulator = Some(DepthEvidence {
            contig: record.contig.clone(),
            start: held.start,
            end: record.end,
            counts,
        });
    }
    if let Some(held) = accumulator {
        if held.length() >= arguments.min_interval_length {
            written.push(held);
        }
    }
    Ok(written)
}

/// `DepthEvidenceCodec.readActualHeader` and then `decode` over every other line: the header's
/// sample names and the records.
///
/// The header is the FIRST line whatever it holds, its columns after the third being the samples.
/// Only a later line beginning `#Chr` is skipped; any other line is a record, and a count is
/// `Integer.parseUnsignedInt`, so a value up to `2^32 - 1` is read and wraps to a negative int.
///
/// `Err` is a file the codec refuses, carrying what is wrong with it. The refusal the reference
/// prints is Tribble's wrapping of the codec's exception, which names an iterator by its identity
/// hash for a bad record, so it is described here rather than reproduced.
pub fn read(text: &str) -> Result<(Vec<String>, Vec<DepthEvidence>), String> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| "the file has no header line".to_string())?;
    let columns: Vec<&str> = header.split('\t').collect();
    if columns.len() < 3 {
        return Err(format!("the header has {} columns", columns.len()));
    }
    let samples = columns[3..].iter().map(|name| name.to_string()).collect();
    let unsigned = |field: &str| {
        field
            .parse::<u32>()
            .map(|value| value as i32)
            .map_err(|_| format!("{field:?} is not an unsigned int"))
    };
    let mut records = Vec::new();
    for line in lines {
        if line.starts_with("#Chr") {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 3 {
            return Err(format!("a record has {} columns", columns.len()));
        }
        records.push(DepthEvidence {
            contig: columns[0].to_string(),
            // Zero-based on disk, one-based here.
            start: unsigned(columns[1])?.wrapping_add(1),
            end: unsigned(columns[2])?,
            counts: columns[3..]
                .iter()
                .map(|count| unsigned(count))
                .collect::<Result<_, _>>()?,
        });
    }
    Ok((samples, records))
}

/// `DepthEvidenceCodec.encode` plus the header the sink writes first.
pub fn write(samples: &[String], records: &[DepthEvidence]) -> String {
    let mut out = String::from("#Chr\tStart\tEnd");
    for sample in samples {
        out.push('\t');
        out.push_str(sample);
    }
    out.push('\n');
    for record in records {
        out.push_str(&format!(
            "{}\t{}\t{}",
            record.contig,
            record.start - 1,
            record.end
        ));
        for count in &record.counts {
            out.push('\t');
            out.push_str(&count.to_string());
        }
        out.push('\n');
    }
    out
}
