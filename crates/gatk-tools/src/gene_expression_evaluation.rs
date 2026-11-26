//! `GeneExpressionEvaluation`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.rnaseq.GeneExpressionEvaluation`
//! (GATK 4.7.0.0; `GeneExpressionEvaluation.java` is unchanged since 4.6.2.0).
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
//! if(multiMapMethod == MultiMapMethod.EQUAL) {
//!     mappingQualityFilter.minMappingQualityScore = 0;
//! }
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
//! `DataLine.set(int, double)`, which the count writer goes through:
//!
//! ```java
//! final long rounded = Math.round(value);
//! if (rounded == value) {
//!     set(index, Long.toString(rounded));
//! } else {
//!     set(index, Double.toString(value));
//! }
//! return set(index, Double.toString(value));
//! ```
//!
//! The `if` writes `1`, and the line after it overwrites that with `1.0` every time. So a count of
//! exactly one is written `1.0`, and the branch that would have written `1` can never be observed.
//!
//! # Reading this file without reading Rust
//!
//! - `&[T]` is a borrowed, read-only view of a list; `&mut` is a writable one, only one at a time;
//! - `Option<T>` is a value or nothing (`Some(v)` / `None`), Java's nullable with the check
//!   enforced; `Result<T, CountError>` is a value or a refusal, and `?` passes the refusal up;
//! - `let Some(x) = y else { .. }` unwraps or leaves through the block;
//! - `match` covers every case, and `Some(x) if cond =>` is a case with a guard;
//! - `#[derive(PartialOrd, Ord)]` orders a struct field by field in declaration order, which is
//!   how `Interval` sorts by contig name, then start, then end, as htsjdk's `compareTo` does;
//! - `.iter().map(..).collect()` is a Java stream; `|x| ..` is a lambda;
//! - `BTreeMap` is Java's `TreeMap`, keys in sorted order;
//! - `f64::from(i)` widens an `i32` to `double` exactly; `x as f64` does the same for `usize`.

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
    /// `getLengthOnReference()`: both ends are included.
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

    /// Whether the two share at least one base.
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
    /// The gff3 attribute the label is read from.
    pub fn key(self) -> &'static str {
        match self {
            FeatureLabel::Name => "Name",
            FeatureLabel::Id => "ID",
        }
    }

    /// The label's first value on `base`, or `None` when the attribute is missing or empty.
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
    /// The feature, shrunk to its label.
    pub base: BaseData,
    /// Its descendants of the overlap types, as `getDescendents()` filtered on
    /// `overlapType.contains(f.getType())` gives them.
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
    /// A good pair with no `MC` tag, when the read is spliced.
    MissingMateCigar,
    /// A grouping feature carrying no value under the label key.
    NoLabel {
        key: String,
        contig: String,
        start: i32,
        end: i32,
    },
}

impl CountError {
    /// The fully qualified class of the exception the reference throws.
    pub fn java_class(&self) -> &'static str {
        match self {
            CountError::UnpairedRead => "java.lang.IllegalStateException",
            CountError::MissingMateQuality | CountError::MissingMateCigar => {
                "org.broadinstitute.hellbender.exceptions.GATKException"
            }
            CountError::NoLabel { .. } => "org.broadinstitute.hellbender.exceptions.UserException",
        }
    }

    /// The exception's message, as the reference builds it.
    pub fn message(&self) -> String {
        match self {
            CountError::UnpairedRead => {
                "Cannot get mate information for an unpaired read".to_string()
            }
            // Two spaces after the full stop, as the reference writes it.
            CountError::MissingMateQuality => {
                "Mate quality must be included.  Consider running FixMateInformation.".to_string()
            }
            CountError::MissingMateCigar => {
                "Mate cigar must be present if using spliced reads".to_string()
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

/// `ReadStrands`: the strand each mate is read from, relative to the transcript.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadStrands {
    ForwardForward,
    ForwardReverse,
    ReverseForward,
    ReverseReverse,
}

impl ReadStrands {
    /// `r1TranscriptionStrand`: whether read one is on the transcribed strand.
    fn r1_transcription_strand(self) -> bool {
        matches!(
            self,
            ReadStrands::ForwardForward | ReadStrands::ForwardReverse
        )
    }

    /// `r2TranscriptionStrand`, the same for read two.
    fn r2_transcription_strand(self) -> bool {
        matches!(
            self,
            ReadStrands::ForwardForward | ReadStrands::ReverseForward
        )
    }

    /// `this.expectReadsOnSameStrand = r1TranscriptionStrand == r2TranscriptionStrand;`
    fn expect_reads_on_same_strand(self) -> bool {
        self.r1_transcription_strand() == self.r2_transcription_strand()
    }

    /// `isSense`, which calls an unstranded feature positive rather than skipping it.
    ///
    /// The reference's comment on that line: "call forward strand sense for unstranded
    /// features".
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
///
/// The reference is `Collections.sort(intervalsCopy)` and then
/// `new IntervalList.IntervalMergerIterator(intervalsCopy.iterator(), true, true, false)`, whose
/// second argument is the one that combines abutting intervals.
pub fn merged_intervals(intervals: &[Interval]) -> Vec<Interval> {
    let mut sorted = intervals.to_vec();
    sorted.sort();
    let mut merged: Vec<Interval> = Vec::new();
    // `last_mut()` is a writable reference to the last merged interval, if any, so the guard can
    // extend it in place; otherwise the interval starts a new run.
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
///
/// The reference first asks `!read.mateIsUnmapped() && read.isProperlyPaired()` and
/// `read.getContig().equals(read.getMateContig())`, then the strand test, then the `MQ` test.
/// `!read.paired` stands for the throw inside `mateIsUnmapped()` (see the module note).
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
    // `as_deref` views the optional mate contig as an optional `&str`; a missing one is not equal.
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
///
/// Spliced: the read's alignment blocks, plus the mate's (from its `MC` cigar) for a good pair,
/// merged. Unspliced: one interval from the leftmost start of the pair over the absolute fragment
/// length, or the read's own span when not a good pair. The reference's
/// `if (read.isUnmapped())` early return of an empty list in the unspliced branch is not here.
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
                .ok_or(CountError::MissingMateCigar)?;
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
///
/// One sum does see that order: PROPORTIONAL's `summedUnNormalizedWeights` adds each feature's
/// weight in the `HashSet` order of `featureOverlapDetector.getOverlaps(..)`, and `summed` here
/// adds them in feature order. Floating-point addition is not associative, so the normalisation
/// factor can differ in its last bit once more than two weights go into the sum; the two-place
/// rounding of the output hides that except at a rounding boundary.
///
/// EQUAL gives each overlapped feature `1.0d/nOverlappingFeatures`. PROPORTIONAL gives each the
/// share of the read's bases it covers, then divides all of them by their sum plus the share
/// nothing covers (see the module note).
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
            // `reduce(0, Integer::sum)` over the merged lengths.
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
                    // The covered bases as a double, divided by `basesOnReference`.
                    let weight = f64::from(covered) / f64::from(bases_on_reference);
                    // `weights.compute(overlap.getKey(), (k,v) -> v == null? weight : v + weight);`
                    match weights.iter_mut().find(|(existing, _)| *existing == index) {
                        Some((_, value)) => *value += weight,
                        None => weights.push((index, weight)),
                    }
                    summed += weight;
                }
            }

            // The share of the read that covered nothing is added as though it were a feature.
            summed += 1.0 - f64::from(total_covered_bases) / f64::from(bases_on_reference);
            // `1.0d/summedUnNormalizedWeights`, then each weight multiplied by it, not divided.
            let normalization = 1.0 / summed;
            for (_, weight) in weights.iter_mut() {
                *weight *= normalization;
            }
            weights
        }
    }
}

/// `MultiMapMethod.getWeights`: unchanged for a uniquely mapped read; for a multi-mapped one,
/// dropped (IGNORE) or divided by the hit count (EQUAL).
///
/// The reference first refuses `nHits < 1` with `GATKException("nHits = " + nHits + ...)`; that
/// check is not here, so an `NH:i:0` read is dropped or divided by zero instead.
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

/// The two counts a feature accumulates, `Coverage.sense_count` and `antisense_count`.
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
    /// Whether `count` applies the mapping quality filter itself. A caller that already ran the
    /// read filters, the mapping quality one among them, passes false.
    pub filter_mapping_quality: bool,
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
///
/// A read is counted when it is first of pair or not in a good pair, so a good pair counts once:
///
/// The reference tests `read.isFirstOfPair()`, then `inGoodPair` with the filter's
/// `minMappingQualityScore`.
/// Its weights go to sense or antisense by `isSense`. A missing `NH` counts as one hit.
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
        if settings.filter_mapping_quality && read.mapping_quality < minimum {
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

/// One line of `FragmentCountWriter`: label, contig, start, stop, strand, sense or antisense, and
/// the count, each quoted only if it needs it.
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

// ================================================================================================
// The gff3 the features come from.
// ================================================================================================

/// One line of a gff3 as `Gff3Codec` decodes it, with the links it makes to other lines.
#[derive(Debug, Clone, PartialEq)]
pub struct GffFeature {
    pub base: BaseData,
    pub score: String,
    pub phase: String,
    /// Indices of the features whose `Parent` names this one's `ID`, in the order they were linked.
    pub children: Vec<usize>,
}

/// `URLDecoder.decode(text, "UTF-8")`: `+` is a space and `%XX` a byte.
///
/// A `%` with fewer than two characters after it is copied as it is, and a `%XX` that is not hex
/// is copied as its three characters; the bytes are then read as UTF-8, with any invalid
/// sequence replaced by U+FFFD.
fn url_decode(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'+' => {
                out.push(b' ');
                i += 1;
            }
            b'%' if i + 2 < bytes.len() => {
                match u8::from_str_radix(&text[i + 1..i + 3], 16) {
                    Ok(byte) => out.push(byte),
                    Err(_) => out.extend_from_slice(&bytes[i..i + 3]),
                }
                i += 3;
            }
            other => {
                out.push(other);
                i += 1;
            }
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// `Gff3Codec` over a whole file, DEEP: every feature line, in file order, each linked to its
/// parents by `ID` among the lines since the last `###` directive. Every feature is returned, a
/// child as much as a top-level one, which is what the codec flushes.
///
/// A `>` line or `##FASTA` ends the features; `###` forgets every `ID` seen so far, so a `Parent`
/// after it cannot reach a feature before it. Links work in both directions: a child seen before
/// its parent is attached when the parent arrives (`waiting_for_parent`). Attribute keys and
/// values, and the first three columns, are URL-decoded; score and phase are kept as text.
pub fn parse_gff3(text: &str) -> Result<Vec<GffFeature>, String> {
    let mut features: Vec<GffFeature> = Vec::new();
    let mut with_id: std::collections::HashMap<String, Vec<usize>> =
        std::collections::HashMap::new();
    let mut waiting_for_parent: std::collections::HashMap<String, Vec<usize>> =
        std::collections::HashMap::new();
    for (number, line) in text.lines().enumerate() {
        if line.starts_with('>') || line.starts_with("##FASTA") {
            break;
        }
        if line.starts_with("###") {
            with_id.clear();
            waiting_for_parent.clear();
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 {
            return Err(format!(
                "Found an invalid number of columns in the given Gff3 file at line + {} - Given: {} Expected: 9 : {line}",
                number + 1,
                fields.len()
            ));
        }
        // Both columns must parse as integers, or the line is refused.
        let (Ok(start), Ok(end)) = (fields[3].parse::<i32>(), fields[4].parse::<i32>()) else {
            return Err(format!(
                "Cannot read integer value for start/end position from line {}.  Line is: {line}",
                number + 1
            ));
        };
        let mut attributes: BTreeMap<String, Vec<String>> = BTreeMap::new();
        if fields[8] != "." {
            for attribute in fields[8].split(';').filter(|a| !a.is_empty()) {
                let parts: Vec<&str> = attribute.split('=').collect();
                if parts.len() != 2 {
                    return Err(format!("Attribute string {} is invalid", fields[8]));
                }
                let values = parts[1]
                    .trim()
                    .split(',')
                    .map(|value| url_decode(value.trim()))
                    .collect();
                attributes.insert(url_decode(parts[0].trim()), values);
            }
        }
        let this = features.len();
        let id = attributes.get("ID").and_then(|ids| ids.first()).cloned();
        let parents = attributes.get("Parent").cloned().unwrap_or_default();
        features.push(GffFeature {
            base: BaseData {
                contig: url_decode(fields[0]),
                source: url_decode(fields[1]),
                kind: url_decode(fields[2]),
                start,
                end,
                strand: Strand::decode(fields[6]),
                attributes,
            },
            score: fields[5].to_string(),
            phase: fields[7].to_string(),
            children: Vec::new(),
        });
        for parent in &parents {
            if let Some(owners) = with_id.get(parent) {
                for owner in owners.clone() {
                    if !features[owner].children.contains(&this) {
                        features[owner].children.push(this);
                    }
                }
            }
            waiting_for_parent
                .entry(parent.clone())
                .or_default()
                .push(this);
        }
        if let Some(id) = id {
            with_id.entry(id.clone()).or_default().push(this);
            if let Some(children) = waiting_for_parent.get(&id) {
                for child in children.clone() {
                    if child != this && !features[this].children.contains(&child) {
                        features[this].children.push(child);
                    }
                }
            }
        }
    }
    Ok(features)
}

/// `getDescendents`: the children, then each child's own descendants, once each and never back
/// through a feature already in the lineage.
pub fn descendants(features: &[GffFeature], index: usize) -> Vec<usize> {
    // `lineage` is the path from `index` down to the current feature, so a cycle in the
    // `Parent` links stops instead of recursing for ever; `out` keeps first-seen order.
    fn walk(features: &[GffFeature], index: usize, lineage: &mut Vec<usize>, out: &mut Vec<usize>) {
        for child in &features[index].children {
            if lineage.contains(child) {
                continue;
            }
            if !out.contains(child) {
                out.push(*child);
            }
            lineage.push(*child);
            walk(features, *child, lineage, out);
            lineage.pop();
        }
    }
    let mut out = Vec::new();
    walk(features, index, &mut vec![index], &mut out);
    out
}
