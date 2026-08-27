//! `ReblockGVCF`: how a GVCF's reference blocks are coarsened and its weak variants demoted.
//!
//! Reading and writing the GVCF are not ported, nor are the annotations the tool recomputes. The
//! banding, the two thresholds and the annotation arguments are.

/// One record of a GVCF: either a reference block or a variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    /// Empty for a reference block, whose only alternate is `<NON_REF>`.
    pub alternates: Vec<String>,
    pub depth: i32,
    pub minimum_depth: i32,
    pub genotype_quality: i32,
    /// `PL[0]`, the likelihood of the reference genotype, which is what the threshold reads.
    pub reference_likelihood: i32,
    pub info: Vec<(String, String)>,
    pub format: Vec<String>,
}

/// The allele every GVCF record carries and which nothing removes.
pub const NON_REF: &str = "<NON_REF>";

impl Record {
    pub fn is_reference_block(&self) -> bool {
        self.alternates.is_empty()
    }
}

/// `GVCFBlock`'s default bounds.
pub const DEFAULT_GQ_BANDS: &[i32] = &[20, 100];

/// Which band a genotype quality falls in: the index of the first bound above it.
///
/// A SINGLE bound of 60 therefore puts 0 and 25 in the same band, which is why the GQ0 block at the
/// end of the fixture merges with the one before it.
pub fn band_of(quality: i32, bands: &[i32]) -> usize {
    bands
        .iter()
        .position(|bound| quality < *bound)
        .unwrap_or(bands.len())
}

/// The band's lower bound, which `--floor-blocks` writes instead of the observed quality.
pub fn band_floor(quality: i32, bands: &[i32]) -> i32 {
    let index = band_of(quality, bands);
    if index == 0 {
        0
    } else {
        bands[index - 1]
    }
}

/// The arguments that change what survives.
#[derive(Debug, Clone, PartialEq)]
pub struct Arguments {
    pub gq_bands: Vec<i32>,
    /// A variant whose reference likelihood is BELOW this becomes a one-base GQ0 block.
    pub rgq_threshold: i32,
    /// Removes a GQ0 reference BLOCK. On its own it does not touch a weak variant, but it removes
    /// the GQ0 block a demoted one becomes, so the two arguments together lose the record.
    pub drop_low_quals: bool,
    pub floor_blocks: bool,
    pub annotations_to_keep: Vec<String>,
    pub format_annotations_to_remove: Vec<String>,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            gq_bands: DEFAULT_GQ_BANDS.to_vec(),
            rgq_threshold: 0,
            drop_low_quals: false,
            floor_blocks: false,
            annotations_to_keep: Vec::new(),
            format_annotations_to_remove: Vec::new(),
        }
    }
}

/// Whether a variant is demoted to a reference block.
///
/// The comparison is strictly less than, so a likelihood exactly at the threshold survives.
pub fn is_demoted(record: &Record, rgq_threshold: i32) -> bool {
    !record.is_reference_block() && record.reference_likelihood < rgq_threshold
}

/// `demote`: the variant becomes a ONE-BASE block at genotype quality zero. It keeps its own start
/// and end, so it never spans anything it did not already.
pub fn demote(record: &Record) -> Record {
    Record {
        alternates: Vec::new(),
        end: record.start,
        genotype_quality: 0,
        minimum_depth: record.depth,
        reference_likelihood: 0,
        info: Vec::new(),
        format: vec![
            "GT".to_string(),
            "DP".to_string(),
            "GQ".to_string(),
            "MIN_DP".to_string(),
            "PL".to_string(),
        ],
        ..record.clone()
    }
}

/// The whole rewrite, in the order the tool makes it.
///
/// A demoted variant does NOT merge with the blocks either side of it, though it is a block in the
/// same band: the merge only ever joins records that arrived as blocks.
pub fn reblock(records: &[Record], arguments: &Arguments) -> Vec<Record> {
    let mut out: Vec<Record> = Vec::new();
    for record in records {
        if record.is_reference_block() {
            if arguments.drop_low_quals && record.genotype_quality == 0 {
                continue;
            }
            // Merge into the previous record only when that record ALSO arrived as a block and
            // sits in the same band.
            let merges = out.last().is_some_and(|previous| {
                previous.is_reference_block()
                    && previous.demoted_from_variant().is_none()
                    && previous.end + 1 == record.start
                    && band_of(previous.genotype_quality, &arguments.gq_bands)
                        == band_of(record.genotype_quality, &arguments.gq_bands)
            });
            if merges {
                let previous = out.last_mut().expect("a previous record");
                previous.end = record.end;
                // The merged block carries the LOWEST quality and the LOWEST minimum depth.
                previous.genotype_quality = previous.genotype_quality.min(record.genotype_quality);
                previous.minimum_depth = previous.minimum_depth.min(record.minimum_depth);
                continue;
            }
            out.push(record.clone());
        } else if is_demoted(record, arguments.rgq_threshold) {
            // The demoted variant is a GQ0 block, and --drop-low-quals removes a GQ0 block
            // whatever produced it: the two arguments together lose the record entirely.
            if arguments.drop_low_quals {
                continue;
            }
            out.push(demote(record));
        } else {
            out.push(record.clone());
        }
    }
    out
}

impl Record {
    /// A demoted variant is a block whose start and end are the same base. Nothing else in a GVCF
    /// is, which is how the merge tells them apart.
    fn demoted_from_variant(&self) -> Option<()> {
        (self.is_reference_block() && self.start == self.end && self.genotype_quality == 0)
            .then_some(())
    }
}

/// What the annotation arguments refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReblockError {
    /// A key asked for by `--annotations-to-keep` that the input's INFO header does not declare.
    /// A FORMAT key is exactly such a case, which is why the two arguments are not symmetric.
    NotInHeader { key: String },
}

impl ReblockError {
    pub fn message(&self) -> String {
        match self {
            ReblockError::NotInHeader { key } => format!(
                "{key} is not in header of input GVCF but was requested to be kept by \
                 annotations-to-keep argument."
            ),
        }
    }
}

/// The check `--annotations-to-keep` makes against the input's INFO keys, before a record is read.
pub fn check_annotations_to_keep(
    requested: &[String],
    info_keys: &[String],
) -> Result<(), ReblockError> {
    for key in requested {
        if !info_keys.contains(key) {
            return Err(ReblockError::NotInHeader { key: key.clone() });
        }
    }
    Ok(())
}
