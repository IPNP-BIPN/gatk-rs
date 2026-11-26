//! `PrintSVEvidence`, ported from `org.broadinstitute.hellbender.tools.sv.PrintSVEvidence`,
//! `DepthEvidenceSortMerger` and `DepthEvidence.extractSamples` (GATK 4.6.2.0).
//!
//! Several SV evidence files merged into one. The walk is [`gatk_engine::multi_feature_walker`];
//! what is here is which samples the output declares and what the sort merger does when two files
//! speak about the same bin.
//!
//! # Merging is only ever a widening
//!
//! Every record is rewritten against the union of the inputs' sample lists, `extractSamples`
//! filling a column its own file does not have with `MISSING_DATA`, which is `-1`. The merger
//! fills only those:
//!
//! ```java
//! if ( count != MISSING_DATA ) {
//!     if ( mergedCounts[idx] == MISSING_DATA ) {
//!         mergedCounts[idx] = count;
//!     } else {
//!         throw new UserException("Multiple sources for count of sample#" + (idx+1) +
//!                 " at " + evidence.getContig() + ":" + evidence.getStart() + "-" +
//!                 evidence.getEnd());
//!     }
//! }
//! ```
//!
//! So two files that both report a sample at one bin are refused, and the message names the sample
//! by ONE-BASED index and the bin by its ONE-BASED interval while the file writes that same bin
//! zero-based.
//!
//! # The sample list is alphabetical, not file order
//!
//! The walker accumulates the headers' sample names into a `TreeSet`, so three files naming zulu,
//! alpha and mike produce the columns `alpha mike zulu`, and the output's header is the header of
//! no input. `--sample-names` replaces that list wholesale, subsetting and reordering the columns,
//! and a name no file knows becomes a column of `-1` rather than a refusal.
//!
//! # Reading this file without reading Rust
//!
//! - `&[T]` is a borrowed, read-only view of an array (Java's `T[]` or `List<T>` passed in); the
//!   caller keeps ownership. `&mut T` is a writable borrow, of which only one may exist at a time;
//! - `Option<T>` is "a value or nothing" (Java's nullable reference, checked by the compiler) and
//!   `Result<T, E>` is "a value or an error", which is how a Java `throw` is modelled here;
//! - `?` after a call that returns a `Result` hands its error straight back to this function's
//!   caller, which is Java's exception propagating through the call site;
//! - a function's last expression, written without `;`, is its return value.

use std::collections::BTreeSet;

/// `DepthEvidence.MISSING_DATA`: the count of a sample the record's own file does not have.
pub const MISSING_DATA: i32 = -1;

/// `DepthEvidence`: one depth-evidence record, one-based and closed as the codec holds it.
///
/// The file writes `start` zero-based (see [`write`]); the refusal message writes it one-based.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DepthEvidence {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub counts: Vec<i32>,
}

/// One input file: its header's sample names, and its records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceFile {
    pub samples: Vec<String>,
    pub records: Vec<DepthEvidence>,
}

/// What the run refuses, one variant per Java `throw` that the ported paths reach.
///
/// A Rust `enum` variant can carry fields, so each refusal holds exactly what its message needs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrintError {
    /// Two files report the same sample at the same bin: `DepthEvidenceSortMerger.merge`.
    MultipleSources {
        /// One-based, as the message writes it.
        sample: usize,
        contig: String,
        start: i32,
        end: i32,
    },
    /// An input's codec produces another feature type than the output path implies:
    /// `onTraversalStart`'s "Incompatible feature input" `UserException`.
    IncompatibleInput {
        path: String,
        input_type: String,
        output_type: String,
        output_path: String,
    },
    /// The output path names no feature type at all, which fails before the tool's own check.
    NoOutputCodec { path: String },
}

impl PrintError {
    /// The class `Main` names: all three are a plain `UserException`.
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    /// The exception's message, built from the same concatenation as the Java `throw`.
    pub fn message(&self) -> String {
        // `match` must name every variant, so a new refusal cannot be added without a message.
        match self {
            PrintError::MultipleSources {
                sample,
                contig,
                start,
                end,
            } => format!("Multiple sources for count of sample#{sample} at {contig}:{start}-{end}"),
            PrintError::IncompatibleInput {
                path,
                input_type,
                output_type,
                output_path,
            } => format!(
                "Incompatible feature input {path} produces features of type {input_type} rather \
                 than features of type {output_type} as dictated by the output path {output_path}"
            ),
            PrintError::NoOutputCodec { path } => {
                format!("No feature output codec found for {path}")
            }
        }
    }
}

/// `FeatureOutputCodecFinder.find`: the feature type the name selects, or `None` where Java
/// throws because no codec claims the path.
pub fn output_feature_type(path: &str) -> Option<&'static str> {
    crate::sv_feature_codecs::find(path).map(|codec| codec.feature_type)
}

/// `onTraversalStart`'s two checks, in the order they fire. `inputs` are the input paths, all of
/// which here produce `DepthEvidence`.
///
/// Java's check between the two (`SVFeature.class.isAssignableFrom(outputClass)`) is not written
/// out: every `FeatureOutputCodec` in GATK 4.7.0.0 is an SV evidence codec, so it cannot fire.
///
/// The Java comparison is `isAssignableFrom`; with no subclassing among the SV feature types an
/// equality of names is the same test.
pub fn check_types(output_path: &str, inputs: &[String]) -> Result<(), PrintError> {
    let output_type = match output_feature_type(output_path) {
        Some(kind) => kind,
        // `find` throws before the tool's own SVFeature check is ever reached.
        None => {
            return Err(PrintError::NoOutputCodec {
                path: output_path.to_string(),
            })
        }
    };
    for path in inputs {
        // `unwrap_or`: take the value, or this default when there is none. An input path the
        // finder does not recognise is still read as depth evidence by the caller.
        let input_type = output_feature_type(path).unwrap_or("DepthEvidence");
        if input_type != output_type {
            return Err(PrintError::IncompatibleInput {
                path: path.clone(),
                input_type: input_type.to_string(),
                output_type: output_type.to_string(),
                output_path: output_path.to_string(),
            });
        }
    }
    Ok(())
}

/// The sample list the output declares: the argument's, or the union of the headers'.
///
/// The union is a `TreeSet`, so it is alphabetical whatever order the files were given in. An
/// empty result turns sample filtering off entirely rather than dropping every record.
///
/// Java's `--sample-names` field is a `LinkedHashSet`, so a repeated name is one column.
/// `requested` is taken as given: the caller in `gatk-cli` drops the repeats before it gets here.
pub fn sample_names(requested: &[String], files: &[EvidenceFile]) -> Vec<String> {
    if !requested.is_empty() {
        return requested.to_vec();
    }
    // `BTreeSet` is Java's `TreeSet`: sorted by the strings' natural order. `flat_map` walks
    // every file's names as one stream; `collect` gathers it into the set's type.
    let union: BTreeSet<&String> = files.iter().flat_map(|file| &file.samples).collect();
    union.into_iter().cloned().collect()
}

/// `DepthEvidence.extractSamples`: one column per requested name, `-1` where the file has none.
///
/// `SVFeaturesHeader.getSampleIndex` is a `HashMap` filled in header order, so a name the header
/// repeats answers with its LAST column; hence `rposition` (search from the end) rather than
/// `position`.
///
/// ```java
/// final Integer sampleIndex = header.getSampleIndex(sampleName);
/// newCounts[idx++] = sampleIndex == null ? MISSING_DATA : counts[sampleIndex];
/// ```
pub fn extract_samples(
    record: &DepthEvidence,
    file_samples: &[String],
    wanted: &[String],
) -> DepthEvidence {
    DepthEvidence {
        contig: record.contig.clone(),
        start: record.start,
        end: record.end,
        counts: wanted
            .iter()
            .map(
                |name| match file_samples.iter().rposition(|own| own == name) {
                    Some(index) => record.counts[index],
                    None => MISSING_DATA,
                },
            )
            .collect(),
    }
}

/// `DepthEvidenceSortMerger.merge`: only a missing column is filled; see the module header.
///
/// Java first throws a `GATKException` when the two count arrays differ in length. Here every
/// record has been rewritten against one sample list, or (with no sample filtering) the caller
/// in `gatk-cli` checks the lengths, so that check is not repeated.
fn merge(held: &mut DepthEvidence, incoming: &DepthEvidence) -> Result<(), PrintError> {
    for (index, count) in incoming.counts.iter().enumerate() {
        if *count == MISSING_DATA {
            continue;
        }
        if held.counts[index] != MISSING_DATA {
            return Err(PrintError::MultipleSources {
                // One-based, as the message writes it.
                sample: index + 1,
                contig: incoming.contig.clone(),
                start: incoming.start,
                end: incoming.end,
            });
        }
        held.counts[index] = *count;
    }
    Ok(())
}

/// The whole run: the walk's merged stream, rewritten against the sample list and sort-merged.
///
/// `merged` is the order the walker hands records over, which for these fixtures is the loci in
/// order; each entry names the file it came from.
///
/// The `held` record is `DepthEvidenceSortMerger.mergedEvidence`: a record at the same contig,
/// start and end (`compareLocatables == 0`) merges into it, anything else flushes it. Java's
/// "features not presented in dictionary order" `GATKException` is not checked here, because the
/// walk only ever hands records over in order.
pub fn run(
    files: &[EvidenceFile],
    merged: &[(usize, DepthEvidence)],
    requested: &[String],
) -> Result<(Vec<String>, Vec<DepthEvidence>), PrintError> {
    let samples = sample_names(requested, files);
    let no_sample_filtering = samples.is_empty();

    let mut written = Vec::new();
    let mut held: Option<DepthEvidence> = None;
    for (file, record) in merged {
        let rewritten = if no_sample_filtering {
            record.clone()
        } else {
            extract_samples(record, &files[*file].samples, &samples)
        };
        // `as_mut` borrows the held record writably in place, so `merge` fills its columns
        // without a copy, as Java mutates `mergedEvidence.getCounts()`.
        match held.as_mut() {
            None => held = Some(rewritten),
            Some(current) => {
                if current.contig == rewritten.contig
                    && current.start == rewritten.start
                    && current.end == rewritten.end
                {
                    merge(current, &rewritten)?;
                } else {
                    written.push(current.clone());
                    held = Some(rewritten);
                }
            }
        }
    }
    // `close()`: the last held record is written when the stream ends.
    if let Some(current) = held {
        written.push(current);
    }
    Ok((samples, written))
}

/// `DepthEvidenceCodec.encode` plus the header the sink writes first, rewritten from the sample
/// list rather than from any input's own header.
pub fn write(samples: &[String], records: &[DepthEvidence]) -> String {
    let mut out = String::from("#Chr\tStart\tEnd");
    for sample in samples {
        out.push('\t');
        out.push_str(sample);
    }
    out.push('\n');
    for record in records {
        // The record holds a one-based start; the file's column is zero-based.
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
