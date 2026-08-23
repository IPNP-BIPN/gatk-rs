//! `GtfToBed`, ported from `org.broadinstitute.hellbender.tools.walkers.conversion.GtfToBed` and
//! `GtfInfo` (GATK 4.6.2.0).
//!
//! A Gencode GTF reduced to one row per gene, or one per transcript. The whole tool is a map keyed
//! by gene or transcript id, a comparator, and four columns.
//!
//! # What it calls a BED is one-based
//!
//! ```java
//! String line = info.getInterval().getContig() + "\t" +
//!         info.getInterval().getStart() + "\t" +
//!         info.getInterval().getEnd() + "\t" +
//!         info.getGeneName();
//! ```
//!
//! `Interval.getStart()` is the GTF's own start, which is one-based and closed, and nothing
//! converts it. Every coordinate in the file is therefore one greater than the same feature would
//! carry in any other BED, and the golden says so on every row.
//!
//! # A gene is as wide as its transcripts
//!
//! Each transcript takes the gene's start down and its end up, so a gene row can be wider than the
//! gene line the file carried. That matters twice, because `--use-basic-transcript` is applied
//! BEFORE the widening: a transcript it drops never widens its gene, so the flag changes the GENE
//! rows as well as which rows are written.
//!
//! # The flag that sounds like a sort is a filter
//!
//! `--sort-by-transcript` decides which of the two kinds of row is written and nothing else. The
//! order is the same either way: the dictionary's contig INDEX, then the start, then the key, so
//! two features beginning at one position are separated by their gene or transcript id compared as
//! STRINGS.

use std::collections::BTreeMap;

/// Which kind of row an entry will be.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntryType {
    Gene,
    Transcript,
}

/// One line of the GTF, reduced to what the tool reads from it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub kind: EntryType,
    pub gene_id: String,
    /// Empty for a gene line.
    pub transcript_id: String,
    pub gene_name: String,
    /// The values of the `tag` attributes, in the order the line carried them.
    pub tags: Vec<String>,
}

/// `GtfInfo`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GtfInfo {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub gene_name: String,
    pub kind: EntryType,
}

/// What the run refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GtfError {
    /// `onTraversalStart`.
    NoDictionary,
    /// The comparator, which is where a contig the dictionary does not know is caught.
    UnknownContig { contig: String },
    /// A transcript whose gene was never seen. The reference has no check for this and dies with a
    /// NullPointerException; no case in the golden reaches it.
    MissingGene { gene_id: String },
}

impl GtfError {
    pub fn java_class(&self) -> &str {
        match self {
            GtfError::NoDictionary => "org.broadinstitute.hellbender.exceptions.UserException",
            GtfError::UnknownContig { .. } => "java.lang.IllegalArgumentException",
            GtfError::MissingGene { .. } => "java.lang.NullPointerException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            GtfError::NoDictionary => {
                "Sequence Dictionary must be specified (sequence-dictionary).".to_string()
            }
            GtfError::UnknownContig { contig } => format!("could not get sequence for {contig}"),
            GtfError::MissingGene { gene_id } => {
                format!("no entry for gene {gene_id}")
            }
        }
    }
}

/// The map the traversal fills: gene and transcript ids to their intervals.
///
/// `use_basic` is `--use-basic-transcript`, which keeps only transcripts carrying a `tag` whose
/// value is `basic`. The reference processes such a transcript once per matching tag, which is
/// idempotent, so a transcript tagged twice lands once.
pub fn collect(
    features: &[Feature],
    use_basic: bool,
) -> Result<BTreeMap<String, GtfInfo>, GtfError> {
    let mut map: BTreeMap<String, GtfInfo> = BTreeMap::new();
    for feature in features {
        match feature.kind {
            EntryType::Gene => {
                map.insert(
                    feature.gene_id.clone(),
                    GtfInfo {
                        contig: feature.contig.clone(),
                        start: feature.start,
                        end: feature.end,
                        gene_name: feature.gene_name.clone(),
                        kind: EntryType::Gene,
                    },
                );
            }
            EntryType::Transcript => {
                let matching = feature.tags.iter().filter(|tag| *tag == "basic").count();
                if use_basic && matching == 0 {
                    continue;
                }
                map.insert(
                    feature.transcript_id.clone(),
                    GtfInfo {
                        contig: feature.contig.clone(),
                        start: feature.start,
                        end: feature.end,
                        gene_name: feature.gene_name.clone(),
                        kind: EntryType::Transcript,
                    },
                );
                let gene = map
                    .get_mut(&feature.gene_id)
                    .ok_or_else(|| GtfError::MissingGene {
                        gene_id: feature.gene_id.clone(),
                    })?;
                if feature.start < gene.start {
                    gene.start = feature.start;
                }
                if feature.end > gene.end {
                    gene.end = feature.end;
                }
            }
        }
    }
    Ok(map)
}

/// `GtfInfoComparator`: contig index, then start, then the key as a string.
pub fn sorted(
    map: &BTreeMap<String, GtfInfo>,
    dictionary: &[String],
) -> Result<Vec<(String, GtfInfo)>, GtfError> {
    let index = |contig: &str| -> Result<usize, GtfError> {
        dictionary
            .iter()
            .position(|name| name == contig)
            .ok_or_else(|| GtfError::UnknownContig {
                contig: contig.to_string(),
            })
    };
    let mut entries: Vec<(usize, i32, String, GtfInfo)> = Vec::new();
    for (key, info) in map {
        entries.push((index(&info.contig)?, info.start, key.clone(), info.clone()));
    }
    entries.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then(left.1.cmp(&right.1))
            .then_with(|| left.2.cmp(&right.2))
    });
    Ok(entries
        .into_iter()
        .map(|(_, _, key, info)| (key, info))
        .collect())
}

/// `formatBedLine`, whose start is the GTF's own.
pub fn bed_line(key: &str, info: &GtfInfo, kind: EntryType) -> String {
    let mut line = format!(
        "{}\t{}\t{}\t{}",
        info.contig, info.start, info.end, info.gene_name
    );
    if kind == EntryType::Transcript {
        line.push(',');
        line.push_str(key);
    }
    line
}

/// The whole tool, from the GTF's features to the file it writes.
///
/// `dictionary` is `--sequence-dictionary`, whose ORDER is what the rows are sorted by. The line
/// separator is the platform's, which on the runner that produced the golden is a newline.
pub fn run(
    features: &[Feature],
    dictionary: Option<&[String]>,
    sort_by_transcript: bool,
    use_basic: bool,
) -> Result<String, GtfError> {
    let Some(dictionary) = dictionary else {
        return Err(GtfError::NoDictionary);
    };
    let map = collect(features, use_basic)?;
    let selected = if sort_by_transcript {
        EntryType::Transcript
    } else {
        EntryType::Gene
    };
    let mut out = String::new();
    for (key, info) in sorted(&map, dictionary)? {
        if info.kind != selected {
            continue;
        }
        out.push_str(&bed_line(&key, &info, selected));
        out.push('\n');
    }
    Ok(out)
}
