//! Ported from `org.broadinstitute.hellbender.tools.walkers.RevertBaseQualityScores` (GATK 4.6.2.0).
//!
//! The third whole tool here, and the other half of G2's calibration gate. Same archetype as
//! [`crate::unmark_duplicates`], same two-line `apply`, and a materially different failure mode.
//!
//! ```java
//! final byte[] originalQuals = ReadUtils.getOriginalBaseQualities(read);
//! if ( originalQuals != null ){
//!     read.setBaseQualities(originalQuals);
//! } else {
//!     throw new UserException("RevertQualityScores can only be applied to SAM/BAM files with "
//!         + "original quality scores, caused by read: " + read.getName());
//! }
//! ```
//!
//! # One read without `OQ` aborts the whole run
//!
//! Not "is skipped", not "is written unchanged": the traversal throws, so the output file is
//! whatever had been flushed before the offending read and the tool exits non-zero. A port that
//! quietly passed such a read through would produce a *larger* and apparently healthier output
//! than the reference, which is the worst shape a divergence can take.
//!
//! The message names the read, so it is part of the observable behaviour and is reproduced. Note
//! the class in it says `RevertQualityScores`, without `Base`, which is not the tool's name — it
//! is transcribed rather than corrected.
//!
//! # Three ways to have no original qualities, and they are not the same test
//!
//! `ReadUtils.getOriginalBaseQualities` returns `null` when the tag is **absent**, and also when it
//! is present and the string is **empty**:
//!
//! ```java
//! if ( ! read.hasAttribute(ORIGINAL_BASE_QUALITIES_TAG) ) { return null; }
//! final String oqString = read.getAttributeAsString(ORIGINAL_BASE_QUALITIES_TAG);
//! return !oqString.isEmpty() ? SAMUtils.fastqToPhred(oqString) : null;
//! ```
//!
//! So an empty `OQ` reaches the tool's own exception rather than producing an empty quality array,
//! even though `fastqToPhred("")` is perfectly happy and returns one. Measured in the oracle
//! rather than inferred.
//!
//! The third way is not a `null` at all: a character outside the printable range makes
//! `fastqToPhred` throw `IllegalArgumentException: Invalid fastq character: <c>` from inside
//! htsjdk, before the tool's own check can run. The valid range is `!` (33) to `~` (126)
//! inclusive, both endpoints measured.

use gatk_engine::reads::{ReadsDataSource, ReadsError};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use crate::sam_output::{header_for_sam_writer, write_records, Options};

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK RevertBaseQualityScores";

/// `SAMTag.OQ`, `ReadUtils.ORIGINAL_BASE_QUALITIES_TAG`.
pub fn original_base_qualities_tag() -> Tag {
    Tag::new(b"OQ")
}

/// Why a read could not be reverted, kept apart because the two come from different code and read
/// differently in a log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RevertError {
    /// The tool's own `UserException`. The message is the reference's, transcribed including the
    /// class name it names, which is not this tool's.
    NoOriginalQualities { read_name: String },
    /// htsjdk's `IllegalArgumentException` out of `SAMUtils.fastqToPhred`, which fires before the
    /// tool's check.
    InvalidFastqCharacter { character: char },
}

impl RevertError {
    /// The message the reference produces, which is observable in the tool's stderr.
    pub fn message(&self) -> String {
        match self {
            RevertError::NoOriginalQualities { read_name } => format!(
                "RevertQualityScores can only be applied to SAM/BAM files with original quality \
                 scores, caused by read: {read_name}"
            ),
            RevertError::InvalidFastqCharacter { character } => {
                format!("Invalid fastq character: {character}")
            }
        }
    }

    /// The Java class the reference throws.
    pub fn class(&self) -> &'static str {
        match self {
            RevertError::NoOriginalQualities { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            RevertError::InvalidFastqCharacter { .. } => "java.lang.IllegalArgumentException",
        }
    }
}

/// `SAMUtils.fastqToPhred`: ASCII minus 33, refusing anything outside the printable range.
pub fn fastq_to_phred(text: &str) -> Result<Vec<u8>, RevertError> {
    let mut phred = Vec::with_capacity(text.len());
    for character in text.chars() {
        let code = character as u32;
        if !(33..=126).contains(&code) {
            return Err(RevertError::InvalidFastqCharacter { character });
        }
        phred.push((code - 33) as u8);
    }
    Ok(phred)
}

/// `ReadUtils.getOriginalBaseQualities`.
///
/// `None` covers both an absent tag and a present-but-empty one, which is the reference's own
/// conflation and is what sends both to the same exception.
pub fn original_base_qualities(read: &BamRecord) -> Result<Option<Vec<u8>>, RevertError> {
    let Some(value) = read.tags.get(original_base_qualities_tag()) else {
        return Ok(None);
    };
    // `getAttributeAsString` on anything other than a string would stringify it; nothing in this
    // archetype's corpus writes `OQ` as another type, and a port that invented a conversion would
    // be guessing.
    let TagValue::Str(text) = value else {
        return Ok(None);
    };
    if text.is_empty() {
        return Ok(None);
    }
    fastq_to_phred(text).map(Some)
}

/// `apply`: revert one read, or fail the run.
pub fn revert(read: &mut BamRecord) -> Result<(), RevertError> {
    match original_base_qualities(read)? {
        Some(qualities) => {
            read.base_qualities = qualities;
            Ok(())
        }
        None => Err(RevertError::NoOriginalQualities {
            read_name: read.read_name.clone(),
        }),
    }
}

/// What a run produces: the written bytes, the tool's own refusal, or a failure to read at all.
///
/// The nesting is the two layers, not an accident. The outer `Result` is the source failing —
/// nothing was ever traversed. The inner one is the tool refusing a read it did traverse, which is
/// a `UserException` in the reference and exits non-zero with a partial file on disk. Collapsing
/// them into one enum would lose which of those happened, and they need different answers.
pub type RunResult = Result<Result<(Vec<u8>, Option<Vec<u8>>), RevertError>, ReadsError>;

/// `RevertBaseQualityScores`, which either reverts every read or fails.
///
/// The reference writes as it goes, so a failure part-way leaves a partial file behind. This
/// returns the error instead of a truncated BAM, because a caller that wants the partial bytes can
/// ask for them and one that does not should not have to notice.
pub fn revert_base_quality_scores(
    source: &ReadsDataSource,
    options: &Options,
    filter: &dyn Fn(&BamRecord) -> bool,
) -> RunResult {
    let mut records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    for record in &mut records {
        if let Err(error) = revert(record) {
            return Ok(Err(error));
        }
    }
    let header = header_for_sam_writer(source.header(), TOOL_NAME, options);
    Ok(Ok(write_records(
        &header,
        &records,
        options.create_output_bam_index,
    )?))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_with(oq: Option<&str>) -> BamRecord {
        let mut read = BamRecord {
            read_name: "r0".to_string(),
            base_qualities: vec![30, 30, 30, 30],
            ..BamRecord::default()
        };
        if let Some(oq) = oq {
            read.tags
                .insert(original_base_qualities_tag(), TagValue::Str(oq.to_string()));
        }
        read
    }

    #[test]
    fn the_tag_is_ascii_minus_thirty_three() {
        assert_eq!(fastq_to_phred("!#5I"), Ok(vec![0, 2, 20, 40]));
        // Both endpoints of the printable range, measured in the oracle.
        assert_eq!(fastq_to_phred("!"), Ok(vec![0]));
        assert_eq!(fastq_to_phred("~"), Ok(vec![93]));
    }

    #[test]
    fn a_character_outside_the_printable_range_is_htsjdks_error_not_the_tools() {
        let space = fastq_to_phred(" ").expect_err("space is one below '!'");
        assert_eq!(space, RevertError::InvalidFastqCharacter { character: ' ' });
        assert_eq!(space.class(), "java.lang.IllegalArgumentException");
        let del = fastq_to_phred("\u{7f}").expect_err("del is one past '~'");
        assert!(matches!(del, RevertError::InvalidFastqCharacter { .. }));
    }

    /// An absent tag and an empty one reach the same exception, which is the reference's own
    /// conflation rather than a simplification here.
    #[test]
    fn absent_and_empty_both_abort_the_run() {
        for oq in [None, Some("")] {
            let mut read = read_with(oq);
            let error = revert(&mut read).expect_err("must abort");
            assert_eq!(
                error,
                RevertError::NoOriginalQualities {
                    read_name: "r0".to_string()
                }
            );
            assert!(error.message().contains("caused by read: r0"));
            // The message names a class that is not this tool's, and that is transcribed.
            assert!(error.message().starts_with("RevertQualityScores can only"));
        }
    }

    #[test]
    fn a_present_tag_replaces_the_qualities() {
        let mut read = read_with(Some("!#5I"));
        revert(&mut read).expect("reverts");
        assert_eq!(read.base_qualities, vec![0, 2, 20, 40]);
    }
}
