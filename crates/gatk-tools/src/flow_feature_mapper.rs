//! `FlowFeatureMapper`: which bases of a read become features, and what each feature's record
//! carries.
//!
//! The flow matrix that produces the score is not ported. What is ported is the walk that finds
//! the features, the surround test that keeps or drops each one, the per-read counts the records
//! carry, and the bounds that take a record away again.
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
