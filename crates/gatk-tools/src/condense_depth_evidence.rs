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
    /// The output's extension implies another feature type, which the message names.
    WrongOutputType {
        path: String,
        found: String,
    },
}

impl CondenseError {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
        match self {
            CondenseError::MinimumAboveMaximum => {
                "Minimum interval length exceeds maximum interval length.".to_string()
            }
            CondenseError::WrongOutputType { path, found } => format!(
                "Output file {path} implies Feature subtype {found}, but this tool expects to \
                 write DepthEvidence."
            ),
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

/// `FeatureOutputCodecFinder.find` for the extensions this tool can meet, which is the name and
/// nothing else.
pub fn output_feature_type(path: &str) -> Option<&'static str> {
    let stripped = path.strip_suffix(".gz").unwrap_or(path);
    if stripped.ends_with(".rd.txt") {
        Some("DepthEvidence")
    } else if stripped.ends_with(".baf.txt") {
        Some("BafEvidence")
    } else if stripped.ends_with(".sr.txt") {
        Some("SplitReadEvidence")
    } else if stripped.ends_with(".pe.txt") {
        Some("DiscordantPairEvidence")
    } else {
        None
    }
}

/// `onTraversalStart`'s second check.
pub fn check_output(path: &str) -> Result<(), CondenseError> {
    match output_feature_type(path) {
        Some("DepthEvidence") => Ok(()),
        Some(found) => Err(CondenseError::WrongOutputType {
            path: path.to_string(),
            found: found.to_string(),
        }),
        None => Ok(()),
    }
}

/// `isAdjacent`.
fn adjacent(left: &DepthEvidence, right: &DepthEvidence) -> bool {
    left.contig == right.contig && left.end + 1 == right.start
}

/// `apply` over every record, then `onTraversalSuccess`: the records the sink was handed.
pub fn condense(records: &[DepthEvidence], arguments: &Arguments) -> Vec<DepthEvidence> {
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
        let mut counts = held.counts.clone();
        for (index, count) in record.counts.iter().enumerate() {
            match counts.get_mut(index) {
                Some(slot) => *slot += count,
                None => counts.push(*count),
            }
        }
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
    written
}

/// `DepthEvidenceCodec.decode` over a whole file: the header's sample names and the records.
pub fn read(text: &str) -> (Vec<String>, Vec<DepthEvidence>) {
    let mut samples = Vec::new();
    let mut records = Vec::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if line.starts_with("#Chr") {
            samples = columns[3..].iter().map(|name| name.to_string()).collect();
            continue;
        }
        if line.starts_with('#') || columns.len() < 3 {
            continue;
        }
        records.push(DepthEvidence {
            contig: columns[0].to_string(),
            // Zero-based on disk, one-based here.
            start: columns[1].parse::<i32>().expect("a start") + 1,
            end: columns[2].parse().expect("an end"),
            counts: columns[3..]
                .iter()
                .map(|count| count.parse().expect("a count"))
                .collect(),
        });
    }
    (samples, records)
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
