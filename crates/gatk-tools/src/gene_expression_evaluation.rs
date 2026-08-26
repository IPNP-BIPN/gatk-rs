//! `GeneExpressionEvaluation`, ported from the tool (GATK 4.6.2.0).
//!
//! RNA-seq fragments counted against gff3 features. A read contributes a WEIGHT rather than a
//! count, and the two ways of computing that weight disagree about what a read half of which lands
//! outside any gene is worth.
//!
//! # A single-end BAM cannot be processed at all
//!
//! ```java
//! boolean ret = !read.mateIsUnmapped() && read.isProperlyPaired() && ...
//! ```
//!
//! `mateIsUnmapped()` is asked first, and it throws `IllegalStateException: Cannot get mate
//! information for an unpaired read`. So every read this tool can see has to be paired; a read
//! that should count alone is a paired read that is not PROPERLY paired, which fails the second
//! test instead. [`in_good_pair`] returns that refusal rather than an answer.
//!
//! # PROPORTIONAL counts the bases the read does not cover
//!
//! ```java
//! summedUnNormalizedWeights += 1.0 - (double)totalCoveredBases/basesOnReference;
//! ```
//!
//! The uncovered fraction of the read is added to the denominator as though it were another
//! feature's share, so a read half of which is intergenic gives its gene about half a count.
//! `EQUAL` divides by the number of features and never looks at a length at all, so the same read
//! gives a whole one.
//!
//! # A good pair is counted once, from read one
//!
//! `apply` returns unless the read is first of pair or not in a good pair, and the intervals it
//! counts are the union of BOTH mates' alignment blocks, taken from the mate's `MC` tag. A pair
//! whose `MQ` tag was never written is a `GATKException` rather than a skipped read.
//!
//! # EQUAL multi-mapping moves the mapping quality filter
//!
//! ```java
//! if(multiMapMethod == MultiMapMethod.EQUAL) { mappingQualityFilter.minMappingQualityScore = 0; }
//! ```
//!
//! Choosing how multi-mapped reads are weighted also changes which reads survive the filter, and
//! which pairs count as good, since `inGoodPair` compares the mate's quality against that same
//! field.
//!
//! # An unstranded feature emits one row, and everything over it is sense
//!
//! `isSense` rewrites `Strand.NONE` as positive before comparing, and `onTraversalSuccess` writes
//! the antisense row only when the strand is not `NONE`.
//!
//! # The integral branch of the writer is dead code
//!
//! ```java
//! final long rounded = Math.round(value);
//! if (rounded == value) { set(index, Long.toString(rounded)); }
//! else { set(index, Double.toString(value)); }
//! return set(index, Double.toString(value));
//! ```
//!
//! The `if` writes `1`, and the line after it overwrites that with `1.0` every time. So a count of
//! exactly one is written `1.0`, and the branch that would have written `1` can never be observed.

use gatk_engine::base_recalibration_engine::round_to_n_decimal_places;
use gatk_engine::tsv_table::{java_double_to_string, quote_if_needed};
use std::collections::BTreeMap;

/// `htsjdk.tribble.annotation.Strand`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Strand {
    Positive,
    Negative,
    None,
}

impl Strand {
    /// `Strand.decode`, over the three spellings a gff3 column carries.
    pub fn decode(text: &str) -> Strand {
        match text {
            "+" => Strand::Positive,
            "-" => Strand::Negative,
            _ => Strand::None,
        }
    }

    /// `Strand.encode`, which is what the count table writes.
    pub fn encode(self) -> &'static str {
        match self {
            Strand::Positive => "+",
            Strand::Negative => "-",
            Strand::None => ".",
        }
    }
}

/// A closed interval on one contig, as `htsjdk.samtools.util.Interval` is used here.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Interval {
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }

    /// `getIntersectionLength`, which is zero for intervals that do not meet.
    pub fn intersection_length(&self, other: &Interval) -> i32 {
        if self.contig != other.contig {
            return 0;
        }
        (self.end.min(other.end) - self.start.max(other.start) + 1).max(0)
    }

    fn overlaps(&self, other: &Interval) -> bool {
        self.intersection_length(other) > 0
    }

    /// `withinDistanceOf(other, 1)`, which is what makes abutting intervals merge.
    fn abuts(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.end + 1 && other.start <= self.end + 1
    }
}

/// `Gff3BaseData`, shrunk as `shrinkBaseData` shrinks it: every attribute but the label is dropped
/// before the feature becomes a map key, so two features that differ only elsewhere collapse.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct BaseData {
    pub contig: String,
    pub source: String,
    pub kind: String,
    pub start: i32,
    pub end: i32,
    pub strand: Strand,
    pub attributes: BTreeMap<String, Vec<String>>,
}

/// `FeatureLabelType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureLabel {
    Name,
    Id,
}

impl FeatureLabel {
    pub fn key(self) -> &'static str {
        match self {
            FeatureLabel::Name => "Name",
            FeatureLabel::Id => "ID",
        }
    }

    pub fn value(self, base: &BaseData) -> Option<&str> {
        base.attributes
            .get(self.key())
            .and_then(|values| values.first())
            .map(String::as_str)
    }
}

/// One grouping feature and the intervals its chosen descendants cover.
#[derive(Debug, Clone, PartialEq)]
pub struct GroupingFeature {
    pub base: BaseData,
    pub overlaps: Vec<Interval>,
}

/// One read, reduced to what the counting reads off it.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub name: String,
    pub contig: String,
    pub start: i32,
    /// Alignment blocks, already decoded from the cigar: the runs of reference the read covers.
    pub blocks: Vec<Interval>,
    pub end: i32,
    pub reverse: bool,
    pub paired: bool,
    pub proper_pair: bool,
    pub first_of_pair: bool,
    pub mate_unmapped: bool,
    pub mate_contig: Option<String>,
    pub mate_start: Option<i32>,
    pub mate_blocks: Option<Vec<Interval>>,
    pub mate_reverse: bool,
    /// The `MQ` tag, absent when it was never written.
    pub mate_quality: Option<i32>,
    pub mapping_quality: i32,
    /// The `NH` tag, whose absence means one.
    pub hits: Option<i32>,
    pub fragment_length: i32,
}

/// What the tool refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CountError {
    /// `mateIsUnmapped()` on a read that is not paired.
    UnpairedRead,
    /// The `MQ` tag missing on a read that reached the good-pair test.
    MissingMateQuality,
    /// A grouping feature carrying no value under the label key.
    NoLabel {
        key: String,
        contig: String,
        start: i32,
        end: i32,
    },
}

impl CountError {
    pub fn java_class(&self) -> &'static str {
        match self {
            CountError::UnpairedRead => "java.lang.IllegalStateException",
            CountError::MissingMateQuality => {
                "org.broadinstitute.hellbender.exceptions.GATKException"
            }
            CountError::NoLabel { .. } => "org.broadinstitute.hellbender.exceptions.UserException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            CountError::UnpairedRead => {
                "Cannot get mate information for an unpaired read".to_string()
            }
            // Two spaces after the full stop, as the reference writes it.
            CountError::MissingMateQuality => {
                "Mate quality must be included.  Consider running FixMateInformation.".to_string()
            }
            CountError::NoLabel {
                key,
                contig,
                start,
                end,
            } => format!("no geneid field {key} found in feature at {contig}:{start}-{end}"),
        }
    }
}

/// `ReadStrands`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadStrands {
    ForwardForward,
    ForwardReverse,
    ReverseForward,
    ReverseReverse,
}

impl ReadStrands {
    fn r1_transcription_strand(self) -> bool {
        matches!(
            self,
            ReadStrands::ForwardForward | ReadStrands::ForwardReverse
        )
    }

    fn r2_transcription_strand(self) -> bool {
        matches!(
            self,
            ReadStrands::ForwardForward | ReadStrands::ReverseForward
        )
    }

    fn expect_reads_on_same_strand(self) -> bool {
        self.r1_transcription_strand() == self.r2_transcription_strand()
    }

    /// `isSense`, which calls an unstranded feature positive rather than skipping it.
    pub fn is_sense(self, read: &Read, feature: &BaseData) -> bool {
        let sense_strand = if feature.strand == Strand::None {
            Strand::Positive
        } else {
            feature.strand
        };
        let is_transcription_strand = if read.first_of_pair {
            self.r1_transcription_strand()
        } else {
            self.r2_transcription_strand()
        };
        let transcription_is_reverse = is_transcription_strand == read.reverse;
        let transcription_strand = if transcription_is_reverse {
            Strand::Negative
        } else {
            Strand::Positive
        };
        transcription_strand == sense_strand
    }
}

/// `MultiOverlapMethod`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiOverlapMethod {
    Equal,
    Proportional,
}

/// `MultiMapMethod`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiMapMethod {
    Ignore,
    Equal,
}

/// `getMergedIntervals`: sorted, then merged, with ABUTTING intervals merged too.
pub fn merged_intervals(intervals: &[Interval]) -> Vec<Interval> {
    let mut sorted = intervals.to_vec();
    sorted.sort();
    let mut merged: Vec<Interval> = Vec::new();
    for interval in sorted {
        match merged.last_mut() {
            Some(current) if current.overlaps(&interval) || current.abuts(&interval) => {
                current.start = current.start.min(interval.start);
                current.end = current.end.max(interval.end);
            }
            _ => merged.push(interval),
        }
    }
    merged
}

/// `inGoodPair`, which asks about the mate before it asks whether there is one.
pub fn in_good_pair(
    read: &Read,
    minimum_mapping_quality: i32,
    read_strands: ReadStrands,
) -> Result<bool, CountError> {
    if !read.paired {
        return Err(CountError::UnpairedRead);
    }
    let mut answer = !read.mate_unmapped
        && read.proper_pair
        && read.mate_contig.as_deref() == Some(read.contig.as_str());
    if answer {
        let reads_on_same_strand = read.reverse == read.mate_reverse;
        answer = read_strands.expect_reads_on_same_strand() == reads_on_same_strand;
    }
    if answer {
        let Some(mate_quality) = read.mate_quality else {
            return Err(CountError::MissingMateQuality);
        };
        answer = mate_quality >= minimum_mapping_quality;
    }
    Ok(answer)
}

/// `getAlignmentIntervals`.
pub fn alignment_intervals(
    read: &Read,
    unspliced: bool,
    minimum_mapping_quality: i32,
    read_strands: ReadStrands,
) -> Result<Vec<Interval>, CountError> {
    if !unspliced {
        let mut intervals = read.blocks.clone();
        if in_good_pair(read, minimum_mapping_quality, read_strands)? {
            let mate_blocks = read
                .mate_blocks
                .as_ref()
                .expect("a good pair carries its mate cigar");
            intervals.extend(mate_blocks.iter().cloned());
        }
        return Ok(merged_intervals(&intervals));
    }
    let good = in_good_pair(read, minimum_mapping_quality, read_strands)?;
    // Unspliced replaces the blocks with one interval, so an intron counts as covered.
    let start = if good {
        read.start
            .min(read.mate_start.expect("a good pair has a mate"))
    } else {
        read.start
    };
    let end = if good {
        start + read.fragment_length.abs() - 1
    } else {
        read.end
    };
    Ok(vec![Interval {
        contig: read.contig.clone(),
        start,
        end,
    }])
}

/// `MultiOverlapMethod.getWeights`, over the features in the order they were declared.
///
/// The weights come back in feature order rather than in the order the reference's `HashSet`
/// happens to give them, which is safe because every consumer of this map is keyed by feature.
pub fn overlap_weights(
    method: MultiOverlapMethod,
    alignment_intervals: &[Interval],
    features: &[GroupingFeature],
) -> Vec<(usize, f64)> {
    match method {
        MultiOverlapMethod::Equal => {
            let overlapping: Vec<usize> = features
                .iter()
                .enumerate()
                .filter(|(_, feature)| {
                    alignment_intervals.iter().any(|interval| {
                        feature
                            .overlaps
                            .iter()
                            .any(|overlap| overlap.overlaps(interval))
                    })
                })
                .map(|(index, _)| index)
                .collect();
            let count = overlapping.len();
            overlapping
                .into_iter()
                .map(|index| (index, 1.0 / count as f64))
                .collect()
        }
        MultiOverlapMethod::Proportional => {
            let merged = merged_intervals(alignment_intervals);
            let bases_on_reference: i32 = merged.iter().map(Interval::length).sum();
            let mut total_covered_bases = 0;
            let mut summed = 0.0;
            let mut weights: Vec<(usize, f64)> = Vec::new();

            for alignment in &merged {
                let mut all_overlapping: Vec<Interval> = Vec::new();
                let mut per_feature: Vec<(usize, Vec<Interval>)> = Vec::new();
                for (index, feature) in features.iter().enumerate() {
                    let hits: Vec<Interval> = feature
                        .overlaps
                        .iter()
                        .filter(|overlap| overlap.overlaps(alignment))
                        .cloned()
                        .collect();
                    if hits.is_empty() {
                        continue;
                    }
                    all_overlapping.extend(hits.iter().cloned());
                    per_feature.push((index, hits));
                }

                let all_merged = merged_intervals(&all_overlapping);
                total_covered_bases += all_merged
                    .iter()
                    .map(|interval| alignment.intersection_length(interval))
                    .sum::<i32>();

                for (index, hits) in per_feature {
                    let hits_merged = merged_intervals(&hits);
                    let covered: i32 = hits_merged
                        .iter()
                        .map(|interval| alignment.intersection_length(interval))
                        .sum();
                    let weight = f64::from(covered) / f64::from(bases_on_reference);
                    match weights.iter_mut().find(|(existing, _)| *existing == index) {
                        Some((_, value)) => *value += weight,
                        None => weights.push((index, weight)),
                    }
                    summed += weight;
                }
            }

            // The share of the read that covered nothing is added as though it were a feature.
            summed += 1.0 - f64::from(total_covered_bases) / f64::from(bases_on_reference);
            let normalization = 1.0 / summed;
            for (_, weight) in weights.iter_mut() {
                *weight *= normalization;
            }
            weights
        }
    }
}

/// `MultiMapMethod.getWeights`.
pub fn multi_map_weights(
    method: MultiMapMethod,
    hits: i32,
    previous: Vec<(usize, f64)>,
) -> Vec<(usize, f64)> {
    if hits == 1 {
        return previous;
    }
    match method {
        MultiMapMethod::Ignore => Vec::new(),
        MultiMapMethod::Equal => previous
            .into_iter()
            .map(|(index, weight)| (index, weight / f64::from(hits)))
            .collect(),
    }
}

/// The two counts a feature accumulates.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Coverage {
    pub sense: f64,
    pub antisense: f64,
}

/// The settings one run counts under.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Settings {
    pub multi_overlap_method: MultiOverlapMethod,
    pub multi_map_method: MultiMapMethod,
    pub read_strands: ReadStrands,
    pub unspliced: bool,
    pub feature_label: FeatureLabel,
    /// `MappingQualityReadFilter.minMappingQualityScore`, before onTraversalStart touches it.
    pub minimum_mapping_quality: i32,
}

impl Settings {
    /// What `onTraversalStart` does to the filter: EQUAL multi-mapping drops it to zero, which
    /// changes which reads are counted as well as how they are weighted.
    pub fn effective_minimum_mapping_quality(&self) -> i32 {
        if self.multi_map_method == MultiMapMethod::Equal {
            0
        } else {
            self.minimum_mapping_quality
        }
    }
}

/// `apply` over every read, then the counts in feature declaration order.
pub fn count(
    features: &[GroupingFeature],
    reads: &[Read],
    settings: &Settings,
) -> Result<Vec<Coverage>, CountError> {
    // The label is read once per feature when the feature is registered, and a missing one is a
    // refusal before any read is seen.
    for feature in features {
        if settings.feature_label.value(&feature.base).is_none() {
            return Err(CountError::NoLabel {
                key: format!("{:?}", settings.feature_label).to_uppercase(),
                contig: feature.base.contig.clone(),
                start: feature.base.start,
                end: feature.base.end,
            });
        }
    }

    let minimum = settings.effective_minimum_mapping_quality();
    let mut coverages = vec![Coverage::default(); features.len()];
    for read in reads {
        if read.mapping_quality < minimum {
            continue;
        }
        if !(read.first_of_pair || !in_good_pair(read, minimum, settings.read_strands)?) {
            continue;
        }
        let intervals =
            alignment_intervals(read, settings.unspliced, minimum, settings.read_strands)?;
        let initial = overlap_weights(settings.multi_overlap_method, &intervals, features);
        let hits = read.hits.unwrap_or(1);
        let final_weights = multi_map_weights(settings.multi_map_method, hits, initial);
        for (index, weight) in final_weights {
            if settings.read_strands.is_sense(read, &features[index].base) {
                coverages[index].sense += weight;
            } else {
                coverages[index].antisense += weight;
            }
        }
    }
    Ok(coverages)
}

/// `FragmentCountWriter`, including the metadata comments `onTraversalSuccess` writes first.
pub fn write_counts(
    features: &[GroupingFeature],
    coverages: &[Coverage],
    sample: &str,
    label: FeatureLabel,
    input_bams: &[String],
    annotation_file: &str,
) -> String {
    let mut out = String::new();
    for (index, bam) in input_bams.iter().enumerate() {
        out.push_str(&format!("#<METADATA>input_bam_{index}={bam}\n"));
    }
    out.push_str(&format!("#<METADATA>annotation_file={annotation_file}\n"));
    out.push_str(&format!(
        "gene_label\tcontig\tstart\tstop\tstrand\tsense_antisense\t{sample}_counts\n"
    ));
    for (feature, coverage) in features.iter().zip(coverages) {
        out.push_str(&row(feature, label, coverage.sense, true));
        if feature.base.strand != Strand::None {
            out.push_str(&row(feature, label, coverage.antisense, false));
        }
    }
    out
}

fn row(feature: &GroupingFeature, label: FeatureLabel, count: f64, sense: bool) -> String {
    let values = [
        label.value(&feature.base).unwrap_or("").to_string(),
        feature.base.contig.clone(),
        feature.base.start.to_string(),
        feature.base.end.to_string(),
        feature.base.strand.encode().to_string(),
        if sense { "sense" } else { "antisense" }.to_string(),
        // Rounded to two places, then written by a `Double.toString` whose integral branch is
        // dead: a count of one is `1.0`, never `1`.
        java_double_to_string(round_to_n_decimal_places(count, 2).expect("two places")),
    ];
    let quoted: Vec<String> = values.iter().map(|value| quote_if_needed(value)).collect();
    format!("{}\n", quoted.join("\t"))
}
