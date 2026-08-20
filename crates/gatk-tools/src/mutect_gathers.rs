//! The three gathers a scattered Mutect run ends with, ported from
//! `MergeMutectStats`, `GatherPileupSummaries` and `GatherNormalArtifactData` (GATK 4.6.2.0).
//!
//! The record formats are already ported: [`gatk_engine::pileup_summary`] and
//! [`crate::get_normal_artifact_data`]. What is here is what each tool does with several shards.
//!
//! # The three disagree about everything else
//!
//! `MergeMutectStats` sums, and refuses any statistic its aggregation map does not hold: that map
//! carries `callable` and nothing else, so a shard with a second statistic ends the run.
//!
//! `GatherPileupSummaries` sorts its FILES by their first record against the sequence dictionary,
//! having first dropped the files with no records at all, and does not sort within a file.
//!
//! `GatherNormalArtifactData` concatenates in the order it was given, empty shards included.
//!
//! # The sum is a double
//!
//! `MutectStats` holds a `double` and the writer formats it as one, so shards of `1` and `2` come
//! back as `3.0`.

use gatk_engine::pileup_summary::{self, PileupSummary, PileupSummaryError};
use gatk_engine::tsv_table::java_double_to_string;

/// `Mutect2Engine.CALLABLE_SITES_NAME`, the one statistic the aggregation map holds.
pub const CALLABLE_SITES_NAME: &str = "callable";

/// `MutectStats`'s two columns.
pub const STATS_COLUMNS: [&str; 2] = ["statistic", "value"];

/// What a gather refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum GatherError {
    /// `Utils.validate` on a statistic the aggregation map does not hold.
    UnknownStatistic(String),
    /// The pileup-summary reader or writer refused.
    PileupSummary(PileupSummaryError),
}

impl GatherError {
    pub fn java_class(&self) -> &'static str {
        match self {
            GatherError::UnknownStatistic(_) => "java.lang.IllegalStateException",
            GatherError::PileupSummary(error) => error.java_class(),
        }
    }

    pub fn message(&self) -> String {
        match self {
            GatherError::UnknownStatistic(statistic) => {
                format!("aggregations list missing key {statistic}")
            }
            // `UserException.BadInput.getMessage()` prefixes `Bad input: `, which the pileup
            // summary's own error carries without: its suite compares the message and this one
            // compares what the tool printed.
            GatherError::PileupSummary(error) => {
                if self.java_class().ends_with("$BadInput") {
                    format!("Bad input: {}", error.message())
                } else {
                    error.message()
                }
            }
        }
    }
}

/// One row of a stats file.
#[derive(Debug, Clone, PartialEq)]
pub struct MutectStat {
    pub statistic: String,
    pub value: f64,
}

/// `MutectStats.readFromFile`, which is the two columns and nothing else.
pub fn read_stats(text: &str) -> Vec<MutectStat> {
    text.lines()
        .skip(1)
        .filter(|line| !line.is_empty())
        .filter_map(|line| line.split_once('\t'))
        .map(|(statistic, value)| MutectStat {
            statistic: statistic.to_string(),
            value: value.parse().unwrap_or(f64::NAN),
        })
        .collect()
}

/// `MergeMutectStats.doWork`: every shard's statistics summed, in the order the keys were first
/// seen.
///
/// The reference accumulates into a `HashMap` and writes the entry set, so a file with more than
/// one statistic would have Java's hash order; every file this port has met carries `callable`
/// alone, which is the only key the aggregation map allows anyway, so first-seen order is enough.
pub fn merge_stats(shards: &[&str]) -> Result<String, GatherError> {
    let mut order: Vec<String> = Vec::new();
    let mut totals: Vec<f64> = Vec::new();
    for shard in shards {
        for stat in read_stats(shard) {
            match order.iter().position(|name| *name == stat.statistic) {
                Some(index) => totals[index] += stat.value,
                None => {
                    order.push(stat.statistic);
                    totals.push(stat.value);
                }
            }
        }
    }
    for statistic in &order {
        if statistic != CALLABLE_SITES_NAME {
            return Err(GatherError::UnknownStatistic(statistic.clone()));
        }
    }
    let mut text = STATS_COLUMNS.join("\t");
    text.push('\n');
    for (statistic, total) in order.iter().zip(totals.iter()) {
        text.push_str(&format!("{statistic}\t{}\n", java_double_to_string(*total)));
    }
    Ok(text)
}

/// `GatherPileupSummaries.doWork`: the empty files dropped, the rest sorted by their FIRST record,
/// then concatenated.
///
/// Each input is its text and the name the messages use.
pub fn gather_pileup_summaries(
    inputs: &[(&str, &str)],
    dictionary: &[String],
) -> Result<String, GatherError> {
    let mut kept: Vec<(&str, &str, PileupSummary)> = Vec::new();
    for (text, source) in inputs {
        let (_, records) =
            pileup_summary::read_from_file(text, source).map_err(GatherError::PileupSummary)?;
        // `removeEmptyFiles` runs before the sort, so an empty shard is never a first record.
        if let Some(first) = records.first() {
            kept.push((text, source, first.clone()));
        }
    }
    kept.sort_by(|left, right| pileup_summary::compare(dictionary, &left.2, &right.2));
    let ordered: Vec<(&str, &str)> = kept
        .iter()
        .map(|(text, source, _)| (*text, *source))
        .collect();
    pileup_summary::gather(&ordered).map_err(GatherError::PileupSummary)
}

/// `GatherNormalArtifactData.doWork`: every shard's records, in the order given.
///
/// The reader and the writer are the record's own, so an empty shard contributes its nothing and
/// the header is written once.
pub fn gather_normal_artifact_data(shards: &[&str]) -> String {
    let mut records = Vec::new();
    for shard in shards {
        for line in shard.lines().skip(1).filter(|line| !line.is_empty()) {
            records.push(line.to_string());
        }
    }
    let mut text = crate::get_normal_artifact_data::COLUMNS.join("\t");
    text.push('\n');
    for record in records {
        text.push_str(&record);
        text.push('\n');
    }
    text
}
