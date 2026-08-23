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
//!     if ( mergedCounts[idx] == MISSING_DATA ) { mergedCounts[idx] = count; }
//!     else { throw new UserException("Multiple sources for count of sample#" + (idx+1) + ...); }
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

use std::collections::BTreeSet;

/// `DepthEvidence.MISSING_DATA`.
pub const MISSING_DATA: i32 = -1;

/// One depth-evidence record, one-based and closed as the codec holds it.
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

/// What the run refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrintError {
    /// Two files report the same sample at the same bin.
    MultipleSources {
        /// One-based, as the message writes it.
        sample: usize,
        contig: String,
        start: i32,
        end: i32,
    },
    /// An input's codec produces another feature type than the output path implies.
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
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
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

/// `FeatureOutputCodecFinder.find`, for the extensions this tool can meet.
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
    } else if stripped.ends_with(".sd.txt") {
        Some("SiteDepth")
    } else {
        None
    }
}

/// `onTraversalStart`'s two checks, in the order they fire. `inputs` are the input paths, all of
/// which here produce `DepthEvidence`.
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
pub fn sample_names(requested: &[String], files: &[EvidenceFile]) -> Vec<String> {
    if !requested.is_empty() {
        return requested.to_vec();
    }
    let union: BTreeSet<&String> = files.iter().flat_map(|file| &file.samples).collect();
    union.into_iter().cloned().collect()
}

/// `DepthEvidence.extractSamples`: one column per requested name, `-1` where the file has none.
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
                |name| match file_samples.iter().position(|own| own == name) {
                    Some(index) => record.counts[index],
                    None => MISSING_DATA,
                },
            )
            .collect(),
    }
}

/// `DepthEvidenceSortMerger.merge`: only a missing column is filled.
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
