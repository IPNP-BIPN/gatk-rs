//! `FlowFeatureMapper`: which bases of a read become features, and what each feature's record
//! carries.
//!
//! Ported: the walk that finds the features, the surround test that keeps or drops each one, the
//! flow score of the read's haplotype against the reference's, the per-read counts the records
//! carry, the bounds that take a record away again, and the two queues that decide when a record
//! is written. The GVCF reference-confidence mode is not.
//!
//! Ported from
//! `org.broadinstitute.hellbender.tools.walkers.featuremapping.SNVMapper`,
//! `org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapper` and
//! `org.broadinstitute.hellbender.tools.walkers.featuremapping.FlowFeatureMapperArgumentCollection`
//! in GATK 4.6.2.0.

/// One cigar element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarElement {
    pub operator: char,
    pub length: i32,
}

impl CigarElement {
    pub fn consumes_read_bases(self) -> bool {
        matches!(self.operator, 'M' | 'I' | 'S' | '=' | 'X')
    }

    pub fn consumes_reference_bases(self) -> bool {
        matches!(self.operator, 'M' | 'D' | 'N' | '=' | 'X')
    }
}

pub fn parse_cigar(text: &str) -> Vec<CigarElement> {
    let mut elements = Vec::new();
    let mut length = 0i32;
    for character in text.chars() {
        if let Some(digit) = character.to_digit(10) {
            length = length * 10 + digit as i32;
        } else {
            elements.push(CigarElement {
                operator: character,
                length,
            });
            length = 0;
        }
    }
    elements
}

/// One read, reduced to what the mapper reads off it.
#[derive(Debug, Clone, PartialEq)]
pub struct Read {
    pub name: String,
    pub contig: String,
    /// The alignment start, which is past any soft clip.
    pub start: i32,
    pub cigar: Vec<CigarElement>,
    /// Every base of the read, soft-clipped ones included.
    pub bases: Vec<u8>,
    pub flags: i32,
    pub mapping_quality: i32,
}

impl Read {
    /// The bases the cigar aligns, which is the read without its soft clips.
    pub fn aligned_bases(&self) -> &[u8] {
        let leading = self
            .cigar
            .first()
            .filter(|element| element.operator == 'S')
            .map_or(0, |element| element.length) as usize;
        let trailing = self
            .cigar
            .last()
            .filter(|element| element.operator == 'S')
            .map_or(0, |element| element.length) as usize;
        &self.bases[leading..self.bases.len() - trailing]
    }

    /// `getUnclippedEnd() - getUnclippedStart() + 1`, which is what `X_LENGTH` carries.
    ///
    /// It counts the soft clips, so a clipped read reports more bases than it aligns.
    pub fn unclipped_length(&self) -> i32 {
        self.cigar
            .iter()
            .filter(|element| element.consumes_reference_bases() || element.operator == 'S')
            .map(|element| element.length)
            .sum()
    }

    pub fn is_duplicate(&self) -> bool {
        self.flags & 1024 != 0
    }
}

/// How many identical bases a feature needs on each side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Surround {
    pub before: i32,
    pub after: i32,
}

impl Surround {
    /// `--snv-identical-bases` and `--snv-identical-bases-after`.
    ///
    /// A zero AFTER means "the same as before", so the two arguments are not symmetric: leaving
    /// the second out is not leaving it at zero.
    pub fn new(before: i32, after: i32) -> Surround {
        Surround {
            before,
            after: if after == 0 { before } else { after },
        }
    }

    /// The shortest cigar element the walk will look inside.
    ///
    /// An element shorter than this is skipped WHOLE, so a mismatch in a three-base match run is
    /// seen at a surround of one and not at a surround of two.
    pub fn minimum_element_length(&self) -> i32 {
        self.before + 1 + self.after
    }
}

impl Default for Surround {
    fn default() -> Self {
        Surround {
            before: 1,
            after: 1,
        }
    }
}

/// One base that became a feature.
#[derive(Debug, Clone, PartialEq)]
pub struct Feature {
    /// The position on the reference.
    pub start: i32,
    pub reference_base: u8,
    pub read_base: u8,
    /// `X_INDEX`: the offset in the WHOLE read, soft clip included.
    pub index: i32,
}

/// Whether a base is surrounded by bases that match the reference.
///
/// An index that falls off either array counts as NOT surrounded, so a mismatch on the first base
/// of the aligned read has nothing before it and never becomes a feature.
pub fn is_surrounded(
    bases: &[u8],
    reference: &[u8],
    read_offset: i32,
    reference_offset: i32,
    surround: Surround,
) -> bool {
    for i in 0..surround.before {
        let base = read_offset - 1 - i;
        let reference_index = reference_offset - 1 - i;
        if base < 0
            || base as usize >= bases.len()
            || reference_index < 0
            || reference_index as usize >= reference.len()
            || bases[base as usize] != reference[reference_index as usize]
        {
            return false;
        }
    }
    for i in 0..surround.after {
        let base = read_offset + 1 + i;
        let reference_index = reference_offset + 1 + i;
        if base < 0
            || base as usize >= bases.len()
            || reference_index < 0
            || reference_index as usize >= reference.len()
            || bases[base as usize] != reference[reference_index as usize]
        {
            return false;
        }
    }
    true
}

/// `nonIdentMBases`: the read's mismatches over its match elements, which is what `X_FC1` carries.
///
/// An `N` in the reference is not a mismatch, so a read over a run of them reports fewer than the
/// bases that differ.
pub fn mismatch_count(read: &Read, reference: &[u8]) -> i32 {
    let mut count = 0;
    let mut read_offset = 0usize;
    let mut reference_offset = 0usize;
    for element in &read.cigar {
        let length = element.length as usize;
        if element.consumes_read_bases() && element.consumes_reference_bases() {
            for offset in 0..length {
                if reference[reference_offset + offset] != b'N'
                    && read.bases[read_offset + offset] != reference[reference_offset + offset]
                {
                    count += 1;
                }
            }
        }
        if element.consumes_read_bases() {
            read_offset += length;
        }
        if element.consumes_reference_bases() {
            reference_offset += length;
        }
    }
    count
}

/// The features one read carries, in read order.
///
/// The walk skips an element shorter than the surround needs, then steps over the surround at
/// each end of the element it does look at, so a mismatch inside the surround of an element's
/// edge is never reached.
pub fn features(read: &Read, reference: &[u8], surround: Surround) -> Vec<Feature> {
    let mut features = Vec::new();
    let mut read_offset = 0i32;
    let mut reference_offset = 0i32;
    for element in &read.cigar {
        let length = element.length;
        if length >= surround.minimum_element_length()
            && element.consumes_read_bases()
            && element.consumes_reference_bases()
        {
            read_offset += surround.before;
            reference_offset += surround.before;
            let mut offset = surround.before;
            while offset < length - surround.after {
                let base = read.bases[read_offset as usize];
                let reference_base = reference[reference_offset as usize];
                if reference_base != b'N'
                    && base != reference_base
                    && is_surrounded(
                        &read.bases,
                        reference,
                        read_offset,
                        reference_offset,
                        surround,
                    )
                {
                    features.push(Feature {
                        start: read.start + reference_offset,
                        reference_base,
                        read_base: base,
                        index: read_offset,
                    });
                }
                offset += 1;
                read_offset += 1;
                reference_offset += 1;
            }
            // The walk stopped `after` short of the element's end, so both offsets catch up.
            read_offset += surround.after;
            reference_offset += surround.after;
        } else {
            if element.consumes_read_bases() {
                read_offset += length;
            }
            if element.consumes_reference_bases() {
                reference_offset += length;
            }
        }
    }
    features
}

/// The Levenshtein distance `X_EDIST` carries, between the read's aligned bases and the reference
/// the walker handed it.
///
/// It is not the mismatch count: an `N` in the reference is a difference here even though it is
/// not a mismatch there.
pub fn edit_distance(a: &[u8], b: &[u8]) -> i32 {
    let mut previous: Vec<i32> = (0..=b.len() as i32).collect();
    let mut current = vec![0i32; b.len() + 1];
    for (i, left) in a.iter().enumerate() {
        current[0] = i as i32 + 1;
        for (j, right) in b.iter().enumerate() {
            let substitution = previous[j] + i32::from(left != right);
            current[j + 1] = substitution.min(previous[j + 1] + 1).min(current[j] + 1);
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[b.len()]
}

/// The arguments that decide whether a record survives.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub surround: Surround,
    pub include_duplicate_reads: bool,
    pub minimum_score: f64,
    pub maximum_score: f64,
    pub exclude_nan_scores: bool,
    /// `--copy-attr`, each `<name>,<type>,<description>`.
    pub copy_attributes: Vec<String>,
    pub copy_attribute_prefix: String,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            surround: Surround::default(),
            include_duplicate_reads: false,
            minimum_score: f64::NEG_INFINITY,
            maximum_score: f64::INFINITY,
            exclude_nan_scores: false,
            copy_attributes: Vec::new(),
            copy_attribute_prefix: String::new(),
        }
    }
}

/// Whether a read is walked at all.
pub fn keeps_read(read: &Read, arguments: &Arguments) -> bool {
    arguments.include_duplicate_reads || !read.is_duplicate()
}

/// Whether a scored feature is written.
///
/// Both bounds are INCLUSIVE at the far side and exclusive at the near one: a score equal to
/// `--max-score` is dropped and a score equal to `--min-score` is kept.
pub fn keeps_score(score: f64, arguments: &Arguments) -> bool {
    if score.is_nan() {
        return !arguments.exclude_nan_scores;
    }
    score <= arguments.maximum_score && score >= arguments.minimum_score
}

/// One `--copy-attr` argument, split.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CopyAttribute {
    pub name: String,
    /// The VCF header type, `String` when the argument does not name one.
    pub kind: String,
    /// The description, which is `copy-attr: <name>` when the argument does not carry one.
    pub description: String,
}

impl CopyAttribute {
    /// `<name>,<type>,<description>`, where the description may itself hold commas: everything
    /// after the second field is joined back together.
    pub fn parse(spec: &str) -> CopyAttribute {
        let parts: Vec<&str> = spec.split(',').collect();
        CopyAttribute {
            name: parts[0].to_string(),
            kind: parts.get(1).map_or("String", |kind| kind).to_string(),
            description: if parts.len() > 2 {
                parts[2..].join(",")
            } else {
                format!("copy-attr: {}", parts[0])
            },
        }
    }

    /// The key the record carries, which is the prefix and the tag's own name.
    pub fn key(&self, prefix: &str) -> String {
        format!("{prefix}{}", self.name)
    }
}

/// The INFO keys every record carries, whatever the arguments.
pub const READ_NAME_KEY: &str = "X_RN";
pub const SCORE_KEY: &str = "X_SCORE";
pub const FLAGS_KEY: &str = "X_FLAGS";
pub const MAPPING_QUALITY_KEY: &str = "X_MAPQ";
pub const CIGAR_KEY: &str = "X_CIGAR";
pub const READ_COUNT_KEY: &str = "X_READ_COUNT";
pub const FILTERED_COUNT_KEY: &str = "X_FILTERED_COUNT";
/// The read's MISMATCH count.
pub const FC1_KEY: &str = "X_FC1";
/// The read's FEATURE count, which is the lower of the two whenever a mismatch failed the
/// surround test.
pub const FC2_KEY: &str = "X_FC2";
pub const LENGTH_KEY: &str = "X_LENGTH";
pub const EDIT_DISTANCE_KEY: &str = "X_EDIST";
pub const INDEX_KEY: &str = "X_INDEX";

/// The INFO column of one record, its keys sorted the way the VCF writer sorts them.
pub fn info_column(pairs: &[(String, String)]) -> String {
    let mut pairs = pairs.to_vec();
    pairs.sort();
    pairs
        .into_iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect::<Vec<_>>()
        .join(";")
}

// ================================================================================================
// The run: SNVMapper over records, the flow score, the two queues and the records they emit.
// ================================================================================================

use crate::flow_based_read::{FlowArguments, FlowRead, FlowReadError};
use htsjdk_bam::cigar::Op;
use htsjdk_bam::record::BamRecord;

/// `LOWEST_PROB`, which stands in for a probability of zero in the log.
const LOWEST_PROB: f64 = 0.0001;
/// `SequenceUtil.VALID_BASES_UPPER`, the alternates `--report-all-alts` scores.
pub const VALID_BASES_UPPER: [u8; 4] = *b"ACGT";

/// Everything the command line decides.
#[derive(Debug, Clone, PartialEq)]
pub struct Settings {
    pub surround_before: i32,
    pub surround_after: i32,
    pub smq_size: Option<i32>,
    pub smq_size_mean: Option<i32>,
    pub report_all_alts: bool,
    pub tag_bases_with_adjacent_ref_diff: bool,
    pub limit_score: f64,
    pub max_score: f64,
    pub min_score: f64,
    pub exclude_nan_scores: bool,
    pub include_dup_reads: bool,
    pub keep_negatives: bool,
    pub keep_supplementary_alignments: bool,
    pub include_qc_failed_reads: bool,
    pub copy_attributes: Vec<CopyAttribute>,
    pub copy_attribute_prefix: String,
    pub flow: FlowArguments,
}

impl Settings {
    /// `SNVMapper`'s spans: zero on both sides when every base is reported or tagged.
    fn ignore_surround(&self) -> bool {
        self.report_all_alts || self.tag_bases_with_adjacent_ref_diff
    }

    fn span_before(&self) -> i32 {
        if self.ignore_surround() {
            0
        } else {
            self.surround_before
        }
    }

    fn span_after(&self) -> i32 {
        if self.ignore_surround() {
            0
        } else {
            self.surround_after
        }
    }

    fn min_cigar_element_length(&self) -> i32 {
        self.span_before() + 1 + self.span_after()
    }

    /// `FlowBasedRead.setMinimalReadLength(1 + 1 + spanAfter)`, set by the mapper's constructor.
    fn minimal_read_length(&self) -> i32 {
        1 + 1 + self.span_after()
    }
}

/// A read with the reference bases its context covers.
pub struct ReadWithReference<'a> {
    pub read: &'a BamRecord,
    pub contig: String,
    /// The read's own span on the reference, which is its context's window.
    pub reference: Vec<u8>,
}

fn consumes_read(op: Op) -> bool {
    matches!(op, Op::M | Op::I | Op::S | Op::Eq | Op::X)
}

fn consumes_reference(op: Op) -> bool {
    matches!(op, Op::M | Op::D | Op::N | Op::Eq | Op::X)
}

fn alignment_end(read: &BamRecord) -> i32 {
    read.alignment_start + read.cigar.reference_length() as i32 - 1
}

fn at(array: &[u8], index: i64) -> Result<u8, FlowReadError> {
    if index < 0 || index >= array.len() as i64 {
        Err(out_of_bounds(index, array.len()))
    } else {
        Ok(array[index as usize])
    }
}

fn out_of_bounds(index: i64, length: usize) -> FlowReadError {
    FlowReadError {
        class: "java.lang.ArrayIndexOutOfBoundsException",
        message: format!("Index {index} out of bounds for length {length}"),
        port_limitation: false,
    }
}

fn gatk_exception(message: &str) -> FlowReadError {
    FlowReadError {
        class: "org.broadinstitute.hellbender.exceptions.GATKException",
        message: message.to_string(),
        port_limitation: false,
    }
}

/// One feature on its way to a record.
#[derive(Debug, Clone)]
pub struct MappedFeature {
    /// The feature's read, as an index into the traversal.
    pub read: usize,
    pub contig: String,
    pub read_base: u8,
    pub ref_base: u8,
    pub read_bases_offset: i32,
    pub start: i32,
    pub offset_delta: i32,
    pub score: f64,
    pub read_count: i32,
    pub filtered_count: i32,
    pub non_ident_m_bases_on_read: i32,
    pub features_on_read: i32,
    pub ref_edit_distance: i32,
    pub index: i32,
    pub smq_left: i32,
    pub smq_right: i32,
    pub smq_left_mean: i32,
    pub smq_right_mean: i32,
    pub score_for_base: Option<[f64; 4]>,
    pub adjacent_ref_diff: bool,
}

/// `surrounded`, the loop both walks share: an index off either array is NOT surrounded.
fn surrounded(
    bases: &[u8],
    reference: &[u8],
    read_offset: i32,
    reference_offset: i32,
    before: i32,
    after: i32,
) -> bool {
    let matches = |b: i32, r: i32| {
        b >= 0
            && (b as usize) < bases.len()
            && r >= 0
            && (r as usize) < reference.len()
            && bases[b as usize] == reference[r as usize]
    };
    (0..before).all(|i| matches(read_offset - 1 - i, reference_offset - 1 - i))
        && (0..after).all(|i| matches(read_offset + 1 + i, reference_offset + 1 + i))
}

/// `calcSmq`: the median or the truncated mean of a clamped window, which `Arrays.copyOfRange`
/// pads with zeros where it runs one past the end.
fn calc_smq(quals: &[u8], from: i32, to: i32, median: bool) -> Result<i32, FlowReadError> {
    let length = quals.len() as i32;
    let from = from.min(length).max(0);
    let to = (to - 1).min(length).max(0);
    if from > to {
        return Err(gatk_exception("invalid qualities range: from > to"));
    }
    let mut range: Vec<i32> = (from..=to)
        .map(|i| quals.get(i as usize).map_or(0, |q| i32::from(*q as i8)))
        .collect();
    if median {
        range.sort();
        let mid = range.len() / 2;
        Ok(if range.len() % 2 == 1 {
            range[mid]
        } else {
            (range[mid - 1] + range[mid]) / 2
        })
    } else {
        Ok(range.iter().sum::<i32>() / range.len() as i32)
    }
}

/// `SNVMapper.forEachOnRead`: the features of one read, each with its per-read counts.
pub fn features_on_read(
    index: usize,
    entry: &ReadWithReference,
    settings: &Settings,
) -> Result<Vec<MappedFeature>, FlowReadError> {
    let read = entry.read;
    let bases = &read.read_bases;
    let reference = &entry.reference;
    let elements = &read.cigar.elements;

    // getSoftStart/getSoftEnd, and the unclipped span `hardLength` measures.
    let leading = |ops: &[Op]| -> i32 {
        let mut total = 0;
        for element in elements {
            if ops.contains(&element.op) {
                total += element.length as i32;
            } else if element.op != Op::H {
                break;
            }
        }
        total
    };
    let trailing = |ops: &[Op]| -> i32 {
        let mut total = 0;
        for element in elements.iter().rev() {
            if ops.contains(&element.op) {
                total += element.length as i32;
            } else if element.op != Op::H {
                break;
            }
        }
        total
    };
    let start_soft_clip = leading(&[Op::S]);
    let end_soft_clip = trailing(&[Op::S]);
    let bases_string: &[u8] = if start_soft_clip == 0 && end_soft_clip == 0 {
        bases
    } else {
        let from = start_soft_clip as usize;
        let to = bases.len().saturating_sub(end_soft_clip as usize);
        if from > to {
            return Err(FlowReadError {
                class: "java.lang.IllegalArgumentException",
                message: format!("{from} > {to}"),
                port_limitation: false,
            });
        }
        &bases[from..to]
    };
    let ref_edit_distance = edit_distance(bases_string, reference);

    let mut non_ident = 0;
    let mut read_offset = 0i64;
    let mut ref_offset = 0i64;
    for element in elements {
        let length = element.length as i64;
        if consumes_read(element.op) && consumes_reference(element.op) {
            for offset in 0..length {
                let r = at(reference, ref_offset + offset)?;
                if r != b'N' && at(bases, read_offset + offset)? != r {
                    non_ident += 1;
                }
            }
        }
        if consumes_read(element.op) {
            read_offset += length;
        }
        if consumes_reference(element.op) {
            ref_offset += length;
        }
    }
    let unclipped_start = read.alignment_start - leading(&[Op::S, Op::H]);
    let unclipped_end = alignment_end(read) + trailing(&[Op::S, Op::H]);
    let hard_length = unclipped_end - unclipped_start + 1;
    let reverse = read.flags & 0x10 != 0;

    let (span_before, span_after) = (settings.span_before(), settings.span_after());
    let mut features = Vec::new();
    let mut read_offset = 0i32;
    let mut ref_offset = 0i32;
    for element in elements {
        let length = element.length as i32;
        if length >= settings.min_cigar_element_length()
            && consumes_read(element.op)
            && consumes_reference(element.op)
        {
            read_offset += span_before;
            ref_offset += span_before;
            let mut offset = span_before;
            while offset < length - span_after {
                let r = at(reference, ref_offset as i64)?;
                if r != b'N'
                    && (settings.report_all_alts || at(bases, read_offset as i64)? != r)
                {
                    let is_surrounded = surrounded(
                        bases,
                        reference,
                        read_offset,
                        ref_offset,
                        settings.surround_before,
                        settings.surround_after,
                    );
                    if settings.ignore_surround() || is_surrounded {
                        let mut feature = MappedFeature {
                            read: index,
                            contig: entry.contig.clone(),
                            read_base: at(bases, read_offset as i64)?,
                            ref_base: r,
                            read_bases_offset: read_offset,
                            start: read.alignment_start + ref_offset,
                            offset_delta: read_offset - ref_offset,
                            score: 0.0,
                            read_count: 0,
                            filtered_count: 0,
                            non_ident_m_bases_on_read: non_ident,
                            features_on_read: 0,
                            ref_edit_distance,
                            index: if reverse {
                                hard_length - read_offset
                            } else {
                                read_offset
                            },
                            smq_left: 0,
                            smq_right: 0,
                            smq_left_mean: 0,
                            smq_right_mean: 0,
                            score_for_base: None,
                            adjacent_ref_diff: settings.ignore_surround() && !is_surrounded,
                        };
                        if settings.smq_size.is_some() || settings.smq_size_mean.is_some() {
                            let mut quals = read.base_qualities.clone();
                            if reverse {
                                quals.reverse();
                            }
                            if let Some(size) = settings.smq_size {
                                feature.smq_left = calc_smq(
                                    &quals,
                                    feature.index - 1 - size,
                                    feature.index - 1,
                                    true,
                                )?;
                                feature.smq_right = calc_smq(
                                    &quals,
                                    feature.index + 1,
                                    feature.index + 1 + size,
                                    true,
                                )?;
                                if reverse {
                                    std::mem::swap(&mut feature.smq_left, &mut feature.smq_right);
                                }
                            }
                            if let Some(size) = settings.smq_size_mean {
                                feature.smq_left_mean = calc_smq(
                                    &quals,
                                    feature.index - 1 - size,
                                    feature.index - 1,
                                    false,
                                )?;
                                feature.smq_right_mean = calc_smq(
                                    &quals,
                                    feature.index + 1,
                                    feature.index + 1 + size,
                                    false,
                                )?;
                                if reverse {
                                    std::mem::swap(
                                        &mut feature.smq_left_mean,
                                        &mut feature.smq_right_mean,
                                    );
                                }
                            }
                        }
                        features.push(feature);
                    }
                }
                offset += 1;
                read_offset += 1;
                ref_offset += 1;
            }
            read_offset += span_after;
            ref_offset += span_after;
        } else {
            if consumes_read(element.op) {
                read_offset += length;
            }
            if consumes_reference(element.op) {
                ref_offset += length;
            }
        }
    }
    let count = features.len() as i32;
    for feature in &mut features {
        feature.features_on_read = count;
    }
    Ok(features)
}

/// `FeatureMapper.FilterStatus`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FilterStatus {
    None,
    Filtered,
    NoFeatureAndFiltered,
}

/// `SNVMapper.noFeatureButFilterAt`, `continue` and all: a position whose surround fails goes on
/// to the next cigar element WITHOUT advancing past the one it was in, so the offsets the rest of
/// the walk uses are the position's own.
fn no_feature_but_filter_at(
    entry: &ReadWithReference,
    start: i32,
    settings: &Settings,
) -> Result<FilterStatus, FlowReadError> {
    let read = entry.read;
    let bases = &read.read_bases;
    let reference = &entry.reference;
    let context_start = read.alignment_start;
    let mut read_offset = 0i32;
    let mut ref_offset = 0i32;
    for element in &read.cigar.elements {
        let length = element.length as i32;
        let includes =
            start >= context_start + ref_offset && start < context_start + ref_offset + length;
        if includes
            && length >= settings.min_cigar_element_length()
            && consumes_read(element.op)
            && consumes_reference(element.op)
        {
            if start < context_start + ref_offset + settings.span_before()
                || start >= context_start + ref_offset + length - settings.span_after()
            {
                return Ok(FilterStatus::Filtered);
            }
            let delta = start - (context_start + ref_offset);
            read_offset += delta;
            ref_offset += delta;
            let no_feature = at(bases, read_offset as i64)? == at(reference, ref_offset as i64)?;
            if !surrounded(
                bases,
                reference,
                read_offset,
                ref_offset,
                settings.surround_before,
                settings.surround_after,
            ) {
                continue;
            }
            return Ok(if no_feature {
                FilterStatus::NoFeatureAndFiltered
            } else {
                FilterStatus::Filtered
            });
        } else {
            if consumes_read(element.op) {
                read_offset += length;
            }
            if consumes_reference(element.op) {
                ref_offset += length;
            }
        }
    }
    Ok(FilterStatus::None)
}

/// `computeLikelihoodLocal`.
fn compute_likelihood_local(
    read: &FlowRead,
    haplotype: &crate::flow_pairhmm_align_reads_to_haplotypes::FlowHaplotype,
    hap_key_length: usize,
) -> Result<f64, FlowReadError> {
    let first = *read.flow_order.first().ok_or_else(|| out_of_bounds(0, 0))?;
    let starting_point = haplotype
        .flow_order
        .iter()
        .position(|base| *base == first)
        .unwrap_or(0);
    let mut result = 0.0;
    for i in 0..read.key.len() {
        let index = i + starting_point;
        if index >= hap_key_length {
            break;
        }
        let location = (haplotype.key[index] & 0xff).min(read.max_hmer + 1);
        let mut prob = read.prob(i, location);
        // Precision.equals(prob, 0.0): within one ulp of zero.
        if prob.abs() <= f64::from_bits(1) {
            prob = LOWEST_PROB;
        }
        result += std::hint::black_box(prob).log10();
    }
    Ok(result)
}

/// `scoreFeature(fr, altBase)`: the read's haplotype against the reference's, from the hmer before
/// the feature to the end of the read. `alt_base` 0 is the reference's own base.
fn score_feature(
    feature: &MappedFeature,
    entry: &ReadWithReference,
    header: &htsjdk_bam::header::SamHeader,
    settings: &Settings,
    alt_base: u8,
) -> Result<f64, FlowReadError> {
    use crate::flow_pairhmm_align_reads_to_haplotypes::FlowHaplotype;
    let read = entry.read;
    if !crate::flow_based_read::has_flow_tags(read) {
        return Err(FlowReadError {
            class: "java.lang.IllegalArgumentException",
            message: "a read without flow tags".to_string(),
            port_limitation: true,
        });
    }
    let info = crate::flow_based_read::read_group_info(read, header)?;

    // buildHaplotypes
    let bases = &read.read_bases;
    let mut offset = feature.read_bases_offset;
    let mut ref_start = feature.start;
    let mut ref_mod_ofs = 0usize;
    if offset > 0 {
        offset -= 1;
        ref_mod_ofs += 1;
        ref_start -= 1;
        let hmer_base = bases[offset as usize];
        while offset > 0 && bases[(offset - 1) as usize] == hmer_base {
            offset -= 1;
            ref_mod_ofs += 1;
            ref_start -= 1;
        }
    }
    let alt_bases = bases[offset as usize..].to_vec();
    let mut ref_bases = alt_bases.clone();
    ref_bases[ref_mod_ofs] = if alt_base != 0 {
        alt_base
    } else {
        feature.ref_base
    };
    let hap_end = ref_start + alt_bases.len() as i32 - 1;
    let key_error = |haplotype: &[u8]| {
        gatk_exception(&format!(
            "baseArrayToKey periodGuard tripped, on {}, flowOrder: {} This probably indicates the presence of a base (value) in the sequence that is not included in the provided flow order",
            String::from_utf8_lossy(haplotype),
            info.flow_order
        ))
    };
    let alt =
        FlowHaplotype::new(&alt_bases, &info.flow_order).ok_or_else(|| key_error(&alt_bases))?;
    let reference =
        FlowHaplotype::new(&ref_bases, &info.flow_order).ok_or_else(|| key_error(&ref_bases))?;

    let mut flow_read = FlowRead::new(read, &info.flow_order, info.max_class, &settings.flow)?;
    let diff_left = ref_start - read.alignment_start + feature.offset_delta;
    let diff_right = alignment_end(read) - hap_end;
    flow_read.apply_base_clipping(
        diff_left.max(0),
        diff_right.max(0),
        false,
        bases.len() as i32,
        settings.minimal_read_length(),
    )?;
    if !flow_read.valid {
        return Ok(-1.0);
    }
    let hap_key_length = alt.key.len().min(reference.key.len());
    let read_score = compute_likelihood_local(&flow_read, &alt, hap_key_length)?;
    let ref_score = compute_likelihood_local(&flow_read, &reference, hap_key_length)?;
    let mut score = read_score - ref_score;
    if !settings.limit_score.is_nan() {
        score = score.min(settings.limit_score);
    }
    if score < 0.0 && !settings.keep_negatives && score != -1.0 {
        score = 0.0;
    }
    Ok(score)
}

/// `filterFeature`.
fn filter_feature(score: f64, settings: &Settings) -> bool {
    if settings.exclude_nan_scores && score.is_nan() {
        false
    } else {
        !(score > settings.max_score || score < settings.min_score)
    }
}

/// `java.util.PriorityQueue` over a comparator: `siftUp` and `siftDown` transcribed, because the
/// order equal elements come out in is the heap's and not the comparator's.
struct JavaHeap<T> {
    queue: Vec<T>,
}

impl<T> JavaHeap<T> {
    fn new() -> Self {
        JavaHeap { queue: Vec::new() }
    }

    fn offer(&mut self, element: T, compare: &dyn Fn(&T, &T) -> std::cmp::Ordering) {
        let mut k = self.queue.len();
        self.queue.push(element);
        while k > 0 {
            let parent = (k - 1) >> 1;
            if compare(&self.queue[k], &self.queue[parent]) != std::cmp::Ordering::Less {
                break;
            }
            self.queue.swap(k, parent);
            k = parent;
        }
    }

    fn peek(&self) -> Option<&T> {
        self.queue.first()
    }

    fn poll(&mut self, compare: &dyn Fn(&T, &T) -> std::cmp::Ordering) -> Option<T> {
        if self.queue.is_empty() {
            return None;
        }
        let result = self.queue.swap_remove(0);
        let n = self.queue.len();
        if n > 0 {
            let mut k = 0;
            let half = n >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let right = child + 1;
                if right < n
                    && compare(&self.queue[child], &self.queue[right])
                        == std::cmp::Ordering::Greater
                {
                    child = right;
                }
                if compare(&self.queue[k], &self.queue[child]) != std::cmp::Ordering::Greater {
                    break;
                }
                self.queue.swap(k, child);
                k = child;
            }
        }
        Some(result)
    }
}

/// The run: every feature that survives, in the order the queue emits it.
pub fn map_features(
    entries: &[ReadWithReference],
    header: &htsjdk_bam::header::SamHeader,
    settings: &Settings,
) -> Result<Vec<MappedFeature>, FlowReadError> {
    use gatk_engine::java_hash::compare_strings;
    let feature_order = |a: &MappedFeature, b: &MappedFeature| {
        compare_strings(&a.contig, &b.contig).then((a.start - b.start).cmp(&0))
    };
    let read_order = |a: &usize, b: &usize| {
        let (x, y) = (&entries[*a], &entries[*b]);
        compare_strings(&x.contig, &y.contig)
            .then(x.read.alignment_start.cmp(&y.read.alignment_start))
            .then(alignment_end(x.read).cmp(&alignment_end(y.read)))
    };
    let mut features: JavaHeap<MappedFeature> = JavaHeap::new();
    let mut reads: JavaHeap<usize> = JavaHeap::new();
    let mut emitted = Vec::new();

    // enrichFeature: every queued read over the feature counts, and the filtered ones twice.
    let enrich = |mut feature: MappedFeature,
                  reads: &JavaHeap<usize>|
     -> Result<MappedFeature, FlowReadError> {
        for index in &reads.queue {
            let entry = &entries[*index];
            if entry.contig == feature.contig
                && entry.read.alignment_start <= feature.start
                && feature.start <= alignment_end(entry.read)
            {
                feature.read_count += 1;
                if no_feature_but_filter_at(entry, feature.start, settings)? != FilterStatus::None {
                    feature.filtered_count += 1;
                }
            }
        }
        Ok(feature)
    };

    for (index, entry) in entries.iter().enumerate() {
        let read = entry.read;
        if read.flags & 0x400 != 0 && !settings.include_dup_reads {
            continue;
        }
        if read.flags & 0x800 != 0 && !settings.keep_supplementary_alignments {
            continue;
        }
        if read.flags & 0x200 != 0 && !settings.include_qc_failed_reads {
            continue;
        }
        // flushQueue(read): the read is queued BEFORE the features ahead of it are written.
        reads.offer(index, &read_order);
        while let Some(head) = features.peek() {
            if head.contig != entry.contig || head.start < read.alignment_start {
                let feature = features.poll(&feature_order).expect("peeked");
                emitted.push(enrich(feature, &reads)?);
            } else {
                break;
            }
        }
        while let Some(head) = reads.peek() {
            let other = &entries[*head];
            if other.contig != entry.contig || alignment_end(other.read) < read.alignment_start {
                reads.poll(&read_order);
            } else {
                break;
            }
        }

        for mut feature in features_on_read(index, entry, settings)? {
            feature.score = score_feature(&feature, entry, header, settings, 0)?;
            if settings.report_all_alts {
                let mut scores = [f64::NAN; 4];
                for (slot, base) in VALID_BASES_UPPER.iter().enumerate() {
                    if *base != feature.read_base {
                        scores[slot] = score_feature(&feature, entry, header, settings, *base)?;
                    }
                }
                feature.score_for_base = Some(scores);
            }
            if filter_feature(feature.score, settings) {
                features.offer(feature, &feature_order);
            }
        }
    }
    while let Some(feature) = features.poll(&feature_order) {
        emitted.push(enrich(feature, &reads)?);
    }
    Ok(emitted)
}
