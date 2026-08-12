//! The four BQSR covariates, ported from
//! `org.broadinstitute.hellbender.utils.recalibration.covariates` (GATK 4.6.2.0).
//!
//! A recalibration table is indexed by covariate keys: `BaseRecalibrator` counts errors against
//! them and `ApplyBQSR` looks corrections up by them, so these encodings decide what a table means.
//! Four covariates, in this order and no other, because the first two are special in the report and
//! the other two are the "additional" ones:
//!
//! | index | covariate | key |
//! |---|---|---|
//! | 0 | `ReadGroup` | the position of the read group's **platform unit** in the header |
//! | 1 | `QualityScore` | the reported quality itself |
//! | 2 | `Context` | the preceding bases packed two bits each, with the length in the low four |
//! | 3 | `Cycle` | the signed cycle, with the sign in the **low bit** |
//!
//! # The key matrix is reused between reads, and it does not leak
//!
//! The reference's `PerReadCovariateMatrix` takes a `CovariateKeyCache` and, on a hit, hands back
//! the same `int[][][]` a previous read of the same length filled in. Nothing clears it. That
//! reads like stale data waiting to happen, and the golden settles it: the dump runs the whole
//! corpus twice, once with one shared cache as BQSR does and once with a fresh cache per read, and
//! the two matrices agree in every cell. Every covariate writes a key at every position of the
//! read, and `ContextCovariate` zeroes the whole array first on the one path that would leave a
//! gap, the low-quality clip shortening the read.
//!
//! So the cache is an allocation optimisation and nothing else, and this port allocates per read.
//! That is a decision the measurement supports rather than an omission.
//!
//! # A missing read group has three different ends
//!
//! `MISSING_READ_GROUP_KEY` is -1 and is documented as what a read outside the table gets. It is
//! reached from [`ReadGroupCovariate::key_from_value`] with an identifier the covariate's own table
//! does not hold. A read whose `RG` names a group the **header** does not declare never gets there:
//! `ReadUtils.getSAMReadGroupRecord` answers null and `getReadGroupIdentifier` dereferences it, so
//! the reference throws `NullPointerException`. And asking the covariate to *format* an unknown key
//! is a third thing, `IllegalStateException: missing key 99`. All three are in the golden.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use crate::base_utils::{
    base_index_to_simple_base, simple_base_to_base_index, simple_reverse_complement,
};
use crate::clipping::{clip_low_qual_ends, ClipError, ClippingRepresentation};
use crate::read;
use crate::read_group;
use crate::read_utils;
use crate::recal_datum::EventType;

/// `ReadGroupCovariate.MISSING_READ_GROUP_KEY`.
pub const MISSING_READ_GROUP_KEY: i32 = -1;

/// `ContextCovariate.UNKNOWN_OR_ERROR_CONTEXT_CODE`.
pub const UNKNOWN_OR_ERROR_CONTEXT_CODE: i32 = -1;

/// `CycleCovariate.CUSHION_FOR_INDELS`: how many bases at each end have no indel cycle key.
pub const CUSHION_FOR_INDELS: i32 = 4;

/// The four bits of a context key that hold the context's length.
const LENGTH_BITS: i32 = 4;

/// The mask for those four bits.
const LENGTH_MASK: i32 = 15;

/// The largest context the covariate accepts: two bits a base, four bits of length, and the leftmost
/// bit left free so no key is negative.
pub const MAX_DNA_CONTEXT: i32 = 13;

/// `QualityUtils.MAX_SAM_QUAL_SCORE`, which is the quality covariate's maximum key.
pub const MAX_SAM_QUAL_SCORE: i32 = 93;

/// The fixed positions of the four covariates in the standard list.
pub const READ_GROUP_COVARIATE_DEFAULT_INDEX: usize = 0;
pub const BASE_QUALITY_COVARIATE_DEFAULT_INDEX: usize = 1;
pub const CONTEXT_COVARIATE_DEFAULT_INDEX: usize = 2;
pub const CYCLE_COVARIATE_DEFAULT_INDEX: usize = 3;
/// `StandardCovariateList.NUM_REQUIRED_COVARITES`, spelled as the reference spells it.
pub const NUM_REQUIRED_COVARITES: usize = 2;

/// The arguments of `RecalibrationArgumentCollection` the covariates read, with their defaults.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecalibrationArguments {
    /// `--mismatches-context-size`.
    pub mismatches_context_size: i32,
    /// `--indels-context-size`.
    pub indels_context_size: i32,
    /// `--maximum-cycle-value`.
    pub maximum_cycle_value: i32,
    /// `--low-quality-tail`, a byte, which the context covariate overwrites with N.
    pub low_qual_tail: u8,
}

impl Default for RecalibrationArguments {
    fn default() -> Self {
        RecalibrationArguments {
            mismatches_context_size: 2,
            indels_context_size: 3,
            maximum_cycle_value: 500,
            low_qual_tail: 2,
        }
    }
}

/// Everything the covariates refuse, with the words the reference refuses it in.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CovariateError {
    /// `CommandLineException.BadArgumentValue` on a context size above [`MAX_DNA_CONTEXT`].
    ContextSizeTooBig { argument: &'static str, size: i32 },
    /// `CommandLineException` on a context size of zero or less. It names **both** sizes, which is
    /// why it carries both.
    ContextSizeNotPositive { mismatches: i32, indels: i32 },
    /// `Utils.validate` inside `ReadGroupCovariate.formatKey`, an `IllegalStateException`.
    MissingKey(i32),
    /// `GATKException` from `contextFromKey` on a negative key.
    NegativeContextKey,
    /// `UserException` from `keyFromCycle`. The cycle here is the **absolute** value, because that
    /// is what the reference's message reports.
    CycleTooBig { maximum: i32, cycle: i32 },
    /// The reference's `NullPointerException`: the read's `RG` names a group the header does not
    /// declare, so `getSAMReadGroupRecord` answers null and `getReadGroupIdentifier` dereferences
    /// it. See the module note.
    NoReadGroupInHeader,
    /// The reference's `ArrayIndexOutOfBoundsException`: the low-quality clip indexes the quality
    /// array by the read's length, and a read with bases and no qualities runs off the end.
    Clip(ClipError),
}

impl CovariateError {
    /// The exact message, which the golden compares character for character.
    pub fn message(&self) -> String {
        match self {
            CovariateError::ContextSizeTooBig { argument, size } => format!(
                "Argument {argument} has a bad value: context size cannot be bigger than {MAX_DNA_CONTEXT}, but was {size}"
            ),
            CovariateError::ContextSizeNotPositive { mismatches, indels } => format!(
                "Context size must be positive. Mismatches: {mismatches} Indels: {indels}"
            ),
            CovariateError::MissingKey(key) => format!("missing key {key}"),
            CovariateError::NegativeContextKey => {
                "dna conversion cannot handle negative numbers. Possible overflow?".to_string()
            }
            // Two spaces after the full stop, because the reference's string has two.
            CovariateError::CycleTooBig { maximum, cycle } => format!(
                "The maximum allowed value for the cycle is {maximum}, but a larger cycle ({cycle}) \
                 was detected.  Please use the --maximum-cycle-value argument (when creating the \
                 recalibration table in BaseRecalibrator) to increase this value (at the expense of \
                 requiring more memory to run)"
            ),
            CovariateError::NoReadGroupInHeader => {
                "Cannot invoke \"htsjdk.samtools.SAMReadGroupRecord.getPlatformUnit()\" \
                 because \"rg\" is null"
                    .to_string()
            }
            CovariateError::Clip(_) => "Index out of bounds".to_string(),
        }
    }
}

/// `PerReadCovariateMatrix`: the keys of one read, indexed by event type, position and covariate.
///
/// The reference's shape is `int[event][position][covariate]` and so is this one. See the module
/// note for why the key cache is not ported.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PerReadCovariateMatrix {
    /// `[event][position][covariate]`.
    covariates: Vec<Vec<Vec<i32>>>,
    /// The column `add_covariate` writes into, set by the list before each covariate runs.
    current_covariate_index: usize,
}

impl PerReadCovariateMatrix {
    pub fn new(read_length: usize, number_of_covariates: usize) -> PerReadCovariateMatrix {
        PerReadCovariateMatrix {
            covariates: vec![
                vec![vec![0; number_of_covariates]; read_length];
                EventType::VALUES.len()
            ],
            current_covariate_index: 0,
        }
    }

    pub fn set_covariate_index(&mut self, index: usize) {
        self.current_covariate_index = index;
    }

    /// `addCovariate(mismatch, insertion, deletion, readOffset)`.
    ///
    /// The reference performs **no bounds check**, "for performance reasons", and its own comment
    /// says an offset past the end is an `ArrayIndexOutOfBoundsException`. Indexing a `Vec` in Rust
    /// panics for the same reason, so the shape of the failure is the same.
    pub fn add_covariate(
        &mut self,
        mismatch: i32,
        insertion: i32,
        deletion: i32,
        read_offset: usize,
    ) {
        let column = self.current_covariate_index;
        self.covariates[EventType::BaseSubstitution.ordinal()][read_offset][column] = mismatch;
        self.covariates[EventType::BaseInsertion.ordinal()][read_offset][column] = insertion;
        self.covariates[EventType::BaseDeletion.ordinal()][read_offset][column] = deletion;
    }

    /// `getCovariatesAtOffset(readPosition, errorModel)`.
    pub fn covariates_at_offset(&self, read_position: usize, error_model: EventType) -> &[i32] {
        &self.covariates[error_model.ordinal()][read_position]
    }

    /// `getMatrixForErrorModel(errorModel)`.
    pub fn matrix_for_error_model(&self, error_model: EventType) -> &[Vec<i32>] {
        &self.covariates[error_model.ordinal()]
    }
}

/// `ReadGroupCovariate`: the read group's platform unit, as a position in the header's list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGroupCovariate {
    /// Identifier to key, in first-seen order and without duplicates, which is the reference's
    /// `LinkedHashMap` built by skipping identifiers already present.
    lookup: Vec<String>,
}

impl ReadGroupCovariate {
    pub fn new(read_groups: &[String]) -> ReadGroupCovariate {
        let mut lookup: Vec<String> = Vec::new();
        for group in read_groups {
            if !lookup.iter().any(|seen| seen == group) {
                lookup.push(group.clone());
            }
        }
        ReadGroupCovariate { lookup }
    }

    /// `getReadGroupIdentifier(rg)`: **the platform unit**, and the id only when there is none.
    ///
    /// This is the whole point of the covariate: a table keyed by `ID` is a different table from
    /// one keyed by `PU`, and two groups sharing a platform unit share a key.
    pub fn read_group_identifier(group: &htsjdk_bam::header::ReadGroup) -> String {
        match group.attributes.get("PU") {
            Some(platform_unit) => platform_unit.to_string(),
            None => group.id.clone(),
        }
    }

    /// `getReadGroupIDs(header)`, whose name says ID and whose values are platform units.
    pub fn read_group_ids(header: &SamHeader) -> Vec<String> {
        header
            .read_groups
            .iter()
            .map(ReadGroupCovariate::read_group_identifier)
            .collect()
    }

    /// `keyForReadGroup`: the position, or [`MISSING_READ_GROUP_KEY`].
    fn key_for_read_group(&self, identifier: &str) -> i32 {
        match self.lookup.iter().position(|seen| seen == identifier) {
            Some(index) => index as i32,
            None => MISSING_READ_GROUP_KEY,
        }
    }

    /// `recordValues`: the same key at every position of the read, for all three event types.
    pub fn record_values(
        &self,
        record: &BamRecord,
        header: &SamHeader,
        values: &mut PerReadCovariateMatrix,
        _record_indel_values: bool,
    ) -> Result<(), CovariateError> {
        // The reference's NullPointerException. See the module note: this is not the -1 key.
        let group =
            read_group::resolve(record, header).ok_or(CovariateError::NoReadGroupInHeader)?;
        let key = self.key_for_read_group(&ReadGroupCovariate::read_group_identifier(group));
        for offset in 0..record.read_bases.len() {
            values.add_covariate(key, key, key, offset);
        }
        Ok(())
    }

    /// `formatKey`: the identifier back. An unknown key is an error here and -1 in
    /// [`ReadGroupCovariate::key_from_value`], which is the asymmetry the golden carries.
    pub fn format_key(&self, key: i32) -> Result<&str, CovariateError> {
        usize::try_from(key)
            .ok()
            .and_then(|index| self.lookup.get(index))
            .map(String::as_str)
            .ok_or(CovariateError::MissingKey(key))
    }

    pub fn key_from_value(&self, value: &str) -> i32 {
        self.key_for_read_group(value)
    }

    pub fn maximum_key_value(&self) -> i32 {
        self.lookup.len() as i32 - 1
    }
}

/// `QualityScoreCovariate`: the reported quality, unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct QualityScoreCovariate;

impl QualityScoreCovariate {
    /// `recordValues`.
    ///
    /// **The loop runs over the base quality COUNT and not the read length**, so a read whose
    /// qualities are absent writes nothing here. It does not survive to be a problem: the context
    /// covariate's low-quality clip indexes the same empty array by the read's length and throws.
    pub fn record_values(
        &self,
        record: &BamRecord,
        values: &mut PerReadCovariateMatrix,
        record_indel_values: bool,
    ) {
        let base_quality_count = record.base_qualities.len();
        if record_indel_values {
            let insertions = read_utils::base_insertion_qualities(record);
            let deletions = read_utils::base_deletion_qualities(record);
            for i in 0..base_quality_count {
                values.add_covariate(
                    record.base_qualities[i] as i32,
                    insertions[i] as i32,
                    deletions[i] as i32,
                    i,
                );
            }
        } else {
            for i in 0..base_quality_count {
                values.add_covariate(record.base_qualities[i] as i32, 0, 0, i);
            }
        }
    }

    pub fn format_key(&self, key: i32) -> String {
        key.to_string()
    }

    pub fn maximum_key_value(&self) -> i32 {
        MAX_SAM_QUAL_SCORE
    }
}

/// `ContextCovariate`: the read's own preceding bases, not the reference's.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContextCovariate {
    mismatches_context_size: i32,
    indels_context_size: i32,
    mismatches_key_mask: i32,
    indels_key_mask: i32,
    low_qual_tail: u8,
}

impl ContextCovariate {
    /// The constructor, whose two checks run in the reference's order: the size ceiling first, for
    /// mismatches then for indels, and the positivity check afterwards for both at once.
    pub fn new(arguments: &RecalibrationArguments) -> Result<ContextCovariate, CovariateError> {
        let mismatches = arguments.mismatches_context_size;
        let indels = arguments.indels_context_size;
        if mismatches > MAX_DNA_CONTEXT {
            return Err(CovariateError::ContextSizeTooBig {
                argument: "mismatches_context_size",
                size: mismatches,
            });
        }
        if indels > MAX_DNA_CONTEXT {
            return Err(CovariateError::ContextSizeTooBig {
                argument: "indels_context_size",
                size: indels,
            });
        }
        if mismatches <= 0 || indels <= 0 {
            return Err(CovariateError::ContextSizeNotPositive { mismatches, indels });
        }
        Ok(ContextCovariate {
            mismatches_context_size: mismatches,
            indels_context_size: indels,
            mismatches_key_mask: create_mask(mismatches),
            indels_key_mask: create_mask(indels),
            low_qual_tail: arguments.low_qual_tail,
        })
    }

    /// `recordValues`.
    ///
    /// Three things happen here that the name does not suggest:
    ///
    ///  * the bases are **clipped and possibly reverse-complemented** first, so the context of a
    ///    negative-strand read is of the reverse complement;
    ///  * if the clip shortened the read, **the whole matrix column is zeroed** before anything is
    ///    written, because the context cannot reach the positions that were clipped away. This is
    ///    the one path that would otherwise leave a gap, and it is why the key cache does not leak;
    ///  * the position written to is the **stranded** offset, counted from the far end for a
    ///    negative-strand read, while the context itself was read forwards along the reverse
    ///    complement.
    pub fn record_values(
        &self,
        record: &BamRecord,
        header: Option<&SamHeader>,
        values: &mut PerReadCovariateMatrix,
        record_indel_values: bool,
    ) -> Result<(), CovariateError> {
        let original_read_length = record.read_bases.len();
        let stranded_clipped_bases = self.stranded_clipped_bytes(record, header)?;

        let context_at_each_cycle = read_context_at_each_position(
            &stranded_clipped_bases,
            self.mismatches_context_size,
            self.mismatches_key_mask,
        );

        let read_length_after_clipping = stranded_clipped_bases.len();
        if read_length_after_clipping != original_read_length {
            for offset in 0..original_read_length {
                values.add_covariate(0, 0, 0, offset);
            }
        }

        let negative_strand = read::is_reverse_strand(record);
        if record_indel_values {
            let indel_keys = read_context_at_each_position(
                &stranded_clipped_bases,
                self.indels_context_size,
                self.indels_key_mask,
            );
            for (i, context_key) in context_at_each_cycle.iter().enumerate() {
                let read_offset = stranded_offset(negative_strand, i, read_length_after_clipping);
                let indel_key = indel_keys[i];
                values.add_covariate(*context_key, indel_key, indel_key, read_offset);
            }
        } else {
            for (i, context_key) in context_at_each_cycle.iter().enumerate() {
                let read_offset = stranded_offset(negative_strand, i, read_length_after_clipping);
                values.add_covariate(*context_key, 0, 0, read_offset);
            }
        }
        Ok(())
    }

    /// `getStrandedClippedBytes`: N over the low-quality tail, then reverse-complemented if the
    /// read is on the negative strand.
    pub fn stranded_clipped_bytes(
        &self,
        record: &BamRecord,
        header: Option<&SamHeader>,
    ) -> Result<Vec<u8>, CovariateError> {
        let clipped = clip_low_qual_ends(
            record,
            header,
            self.low_qual_tail,
            ClippingRepresentation::WriteNs,
        )
        .map_err(CovariateError::Clip)?;
        if read::is_reverse_strand(record) {
            Ok(simple_reverse_complement(&clipped.read_bases))
        } else {
            Ok(clipped.read_bases)
        }
    }

    /// `formatKey`: `None` for -1, which the reference spells as a Java null and never writes to a
    /// report.
    pub fn format_key(&self, key: i32) -> Result<Option<String>, CovariateError> {
        if key == -1 {
            return Ok(None);
        }
        context_from_key(key).map(Some)
    }

    pub fn key_from_value(&self, value: &str) -> i32 {
        key_from_context(value.as_bytes())
    }

    /// `maximumKeyValue()`: every base a `T`, over the longer of the two context sizes.
    pub fn maximum_key_value(&self) -> i32 {
        let length = self.mismatches_context_size.max(self.indels_context_size);
        let mut key = length;
        let mut bit_offset = LENGTH_BITS;
        for _ in 0..length {
            key |= 3 << bit_offset;
            bit_offset += 2;
        }
        key
    }
}

/// `getStrandedOffset`: counted from the far end for a negative-strand read.
pub fn stranded_offset(is_negative_strand: bool, offset: usize, read_length: usize) -> usize {
    if is_negative_strand {
        read_length - offset - 1
    } else {
        offset
    }
}

/// `createMask(contextSize)`: two bits a base, shifted past the four length bits.
fn create_mask(context_size: i32) -> i32 {
    let mut mask = 0;
    for _ in 0..context_size {
        mask = (mask << 2) | 3;
    }
    mask << LENGTH_BITS
}

/// `keyFromContext(dna, start, end)`: the length in the low four bits, then two bits a base.
///
/// A single non-ACGT base makes the whole key -1, and the length is `end - start` rather than the
/// number of bases actually packed, which is what lets [`context_from_key`] read it back.
pub fn key_from_context(dna: &[u8]) -> i32 {
    let mut key = dna.len() as i32;
    let mut bit_offset = LENGTH_BITS;
    for base in dna {
        let base_index = simple_base_to_base_index(*base);
        if base_index == -1 {
            return -1;
        }
        // Java's `<<` on an `int` uses the shift distance modulo 32 and never faults. Rust's
        // panics once the distance reaches the width, so the wrapping form is the faithful one: a
        // context longer than the covariate allows still produces a key rather than a crash.
        key |= base_index.wrapping_shl(bit_offset as u32);
        bit_offset += 2;
    }
    key
}

/// `contextFromKey(key)`: the bases back out, as many as the length nibble claims.
///
/// **The length is trusted.** A key whose low four bits say fifteen decodes fifteen bases whether or
/// not it holds them, and the ones past the end come out as `.` from
/// [`base_index_to_simple_base`]. The golden carries `contextFromKey(4095) = TTTTAAAAAAAAAAA`,
/// eleven of whose characters are read out of bits that were never written.
pub fn context_from_key(key: i32) -> Result<String, CovariateError> {
    if key < 0 {
        return Err(CovariateError::NegativeContextKey);
    }
    let length = key & LENGTH_MASK;
    let mut mask: i32 = 48;
    let mut offset = LENGTH_BITS;
    let mut dna = String::with_capacity(length as usize);
    for _ in 0..length {
        // Both shifts are Java's, which take the distance modulo 32 and let the value overflow
        // silently. A length of fifteen runs the offset past the width of an `int`, and Rust's own
        // shift would panic there where the reference simply wraps and reads zero.
        let base_index = (key & mask).wrapping_shr(offset as u32);
        dna.push(base_index_to_simple_base(base_index) as char);
        mask = mask.wrapping_shl(2);
        offset += 2;
    }
    Ok(dna)
}

/// `getReadContextAtEachPosition`: the preceding n-base context at every position of the read.
///
/// Two things here are the reference's rather than the obvious implementation's.
///
/// **The list can be shorter than the read.** The first `contextSize - 1` positions have no context
/// and get -1, and if the read is shorter than the context size the function returns there, with
/// one entry per base and no key at all. The callers index it by `i` up to the read's length, which
/// is safe only because those two counts coincide.
///
/// **The recovery after an N walks backwards.** When the very first context contains a non-ACGT
/// base the key is -1, and the loop that rebuilds it starts at `contextSize - 1` and walks *down*
/// while the bases are ACGT, filling the key from the high bits down. It is not a restart from the
/// N: it is the same key, assembled from the other end.
pub fn read_context_at_each_position(bases: &[u8], context_size: i32, mask: i32) -> Vec<i32> {
    let read_length = bases.len() as i32;
    let mut keys: Vec<i32> = Vec::with_capacity(bases.len());

    // The first contextSize-1 bases have nothing in front of them.
    let mut i = 1;
    while i < context_size && i <= read_length {
        keys.push(UNKNOWN_OR_ERROR_CONTEXT_CODE);
        i += 1;
    }

    if read_length < context_size {
        return keys;
    }

    let new_base_offset = 2 * (context_size - 1) + LENGTH_BITS;

    let mut current_key = key_from_context(&bases[0..context_size as usize]);
    keys.push(current_key);

    // A non-ACGT base in the first context: rebuild the key from the penalty position downwards.
    let mut current_n_penalty = 0;
    if current_key == -1 {
        current_key = 0;
        current_n_penalty = context_size - 1;
        let mut offset = new_base_offset;
        loop {
            let base_index = simple_base_to_base_index(bases[current_n_penalty as usize]);
            if base_index == -1 {
                break;
            }
            current_key |= base_index << offset;
            offset -= 2;
            current_n_penalty -= 1;
        }
    }

    for current_index in context_size..read_length {
        let base_index = simple_base_to_base_index(bases[current_index as usize]);
        if base_index == -1 {
            current_n_penalty = context_size;
            current_key = 0;
        } else {
            current_key = (current_key >> 2) & mask;
            current_key |= base_index << new_base_offset;
            current_key |= context_size;
        }

        if current_n_penalty == 0 {
            keys.push(current_key);
        } else {
            current_n_penalty -= 1;
            keys.push(-1);
        }
    }

    keys
}

/// `CycleCovariate`: the position in the read, signed, with the sign in the low bit of the key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CycleCovariate {
    maximum_cycle_value: i32,
}

impl CycleCovariate {
    pub fn new(arguments: &RecalibrationArguments) -> CycleCovariate {
        CycleCovariate {
            maximum_cycle_value: arguments.maximum_cycle_value,
        }
    }

    pub fn record_values(
        &self,
        record: &BamRecord,
        values: &mut PerReadCovariateMatrix,
        record_indel_values: bool,
    ) -> Result<(), CovariateError> {
        let read_length = record.read_bases.len();
        for i in 0..read_length {
            let substitution = cycle_key(i as i32, record, false, self.maximum_cycle_value)?;
            if record_indel_values {
                let indel = cycle_key(i as i32, record, true, self.maximum_cycle_value)?;
                values.add_covariate(substitution, indel, indel, i);
            } else {
                values.add_covariate(substitution, 0, 0, i);
            }
        }
        Ok(())
    }

    pub fn format_key(&self, key: i32) -> String {
        cycle_from_key(key).to_string()
    }

    pub fn key_from_value(&self, cycle: i32) -> Result<i32, CovariateError> {
        key_from_cycle(cycle, self.maximum_cycle_value)
    }

    pub fn maximum_key_value(&self) -> i32 {
        (self.maximum_cycle_value << 1) + 1
    }
}

/// `cycleKey(baseNumber, read, indel, maxCycle)`.
///
/// The cycle counts up from one for a forward first-of-pair read, down from the read's length for a
/// negative-strand one, and the whole thing is negated for a second-of-pair read, which is what
/// makes all four corners different. An indel key is additionally -1 within [`CUSHION_FOR_INDELS`]
/// bases of either end.
pub fn cycle_key(
    base_number: i32,
    record: &BamRecord,
    indel: bool,
    max_cycle: i32,
) -> Result<i32, CovariateError> {
    let is_neg_strand = read::is_reverse_strand(record);
    let is_second_in_pair = read::is_paired(record) && read::is_second_of_pair(record);
    let read_length = record.read_bases.len() as i32;

    let read_order_factor = if is_second_in_pair { -1 } else { 1 };
    let (mut cycle, increment) = if is_neg_strand {
        (read_length * read_order_factor, -read_order_factor)
    } else {
        (read_order_factor, read_order_factor)
    };

    cycle += base_number * increment;

    if !indel {
        return key_from_cycle(cycle, max_cycle);
    }
    let max_cycle_for_indels = read_length - CUSHION_FOR_INDELS - 1;
    if base_number < CUSHION_FOR_INDELS || base_number > max_cycle_for_indels {
        Ok(-1)
    } else {
        key_from_cycle(cycle, max_cycle)
    }
}

/// `cycleFromKey(key)`: shift the sign bit off, then apply it.
///
/// Key 1 is cycle "negative zero", which is zero, so two keys decode to the same cycle and only one
/// of them is ever encoded.
pub fn cycle_from_key(key: i32) -> i32 {
    let mut cycle = key >> 1;
    if key & 1 != 0 {
        cycle *= -1;
    }
    cycle
}

/// `keyFromCycle(cycle, maxCycle)`.
///
/// The check is on the **absolute** value and so is the message, which is why refusing cycle -501
/// reports "a larger cycle (501)".
pub fn key_from_cycle(cycle: i32, max_cycle: i32) -> Result<i32, CovariateError> {
    let mut result = cycle.abs();
    if result > max_cycle {
        return Err(CovariateError::CycleTooBig {
            maximum: max_cycle,
            cycle: result,
        });
    }
    result <<= 1;
    if cycle < 0 {
        result += 1;
    }
    Ok(result)
}

/// One of the four, named the way the report names them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CovariateKind {
    ReadGroup,
    QualityScore,
    Context,
    Cycle,
}

impl CovariateKind {
    /// The Java simple class name, which is what `getStandardCovariateClassNames` returns.
    pub fn class_name(&self) -> &'static str {
        match self {
            CovariateKind::ReadGroup => "ReadGroupCovariate",
            CovariateKind::QualityScore => "QualityScoreCovariate",
            CovariateKind::Context => "ContextCovariate",
            CovariateKind::Cycle => "CycleCovariate",
        }
    }

    /// `parseNameForReport()`: the class name split on "Covariate", first part.
    pub fn parsed_name(&self) -> &'static str {
        match self {
            CovariateKind::ReadGroup => "ReadGroup",
            CovariateKind::QualityScore => "QualityScore",
            CovariateKind::Context => "Context",
            CovariateKind::Cycle => "Cycle",
        }
    }
}

/// `StandardCovariateList`: the four, in the one order a recalibration table is written in.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StandardCovariateList {
    pub read_group: ReadGroupCovariate,
    pub quality_score: QualityScoreCovariate,
    pub context: ContextCovariate,
    pub cycle: CycleCovariate,
}

impl StandardCovariateList {
    pub fn new(
        arguments: &RecalibrationArguments,
        all_read_groups: &[String],
    ) -> Result<StandardCovariateList, CovariateError> {
        Ok(StandardCovariateList {
            read_group: ReadGroupCovariate::new(all_read_groups),
            quality_score: QualityScoreCovariate,
            context: ContextCovariate::new(arguments)?,
            cycle: CycleCovariate::new(arguments),
        })
    }

    /// The other constructor, which takes the read groups from the header.
    pub fn from_header(
        arguments: &RecalibrationArguments,
        header: &SamHeader,
    ) -> Result<StandardCovariateList, CovariateError> {
        StandardCovariateList::new(arguments, &ReadGroupCovariate::read_group_ids(header))
    }

    /// Four.
    pub fn size(&self) -> usize {
        4
    }

    /// Two: the read group and the quality score, which the report writes differently.
    pub fn number_of_special_covariates(&self) -> usize {
        NUM_REQUIRED_COVARITES
    }

    pub fn kinds(&self) -> [CovariateKind; 4] {
        [
            CovariateKind::ReadGroup,
            CovariateKind::QualityScore,
            CovariateKind::Context,
            CovariateKind::Cycle,
        ]
    }

    /// The two that are not special, in order.
    pub fn additional_covariates(&self) -> [CovariateKind; 2] {
        [CovariateKind::Context, CovariateKind::Cycle]
    }

    /// `covariateNames()`: the class names joined by commas, which is what the report's header
    /// carries.
    pub fn covariate_names(&self) -> String {
        self.kinds()
            .iter()
            .map(|kind| kind.class_name())
            .collect::<Vec<_>>()
            .join(",")
    }

    /// `indexByClass`, which is the position in the list, or -1.
    pub fn index_by_class(&self, kind: CovariateKind) -> i32 {
        self.kinds()
            .iter()
            .position(|seen| *seen == kind)
            .map(|index| index as i32)
            .unwrap_or(-1)
    }

    /// `getCovariateByParsedName`, which is how a report's column name finds its covariate.
    pub fn covariate_by_parsed_name(&self, name: &str) -> Option<CovariateKind> {
        self.kinds()
            .into_iter()
            .find(|kind| kind.parsed_name() == name)
    }

    pub fn maximum_key_value(&self, kind: CovariateKind) -> i32 {
        match kind {
            CovariateKind::ReadGroup => self.read_group.maximum_key_value(),
            CovariateKind::QualityScore => self.quality_score.maximum_key_value(),
            CovariateKind::Context => self.context.maximum_key_value(),
            CovariateKind::Cycle => self.cycle.maximum_key_value(),
        }
    }

    /// `populatePerReadCovariateMatrix`: every covariate over every position, in list order.
    ///
    /// The column each covariate writes into is set by the list before the covariate runs, which
    /// the reference's own `TODO` calls a pattern to avoid. It is kept because the column index is
    /// the covariate's position in this list and nothing else defines it.
    pub fn populate_per_read_covariate_matrix(
        &self,
        record: &BamRecord,
        header: &SamHeader,
        values: &mut PerReadCovariateMatrix,
        record_indel_values: bool,
    ) -> Result<(), CovariateError> {
        values.set_covariate_index(READ_GROUP_COVARIATE_DEFAULT_INDEX);
        self.read_group
            .record_values(record, header, values, record_indel_values)?;
        values.set_covariate_index(BASE_QUALITY_COVARIATE_DEFAULT_INDEX);
        self.quality_score
            .record_values(record, values, record_indel_values);
        values.set_covariate_index(CONTEXT_COVARIATE_DEFAULT_INDEX);
        self.context
            .record_values(record, Some(header), values, record_indel_values)?;
        values.set_covariate_index(CYCLE_COVARIATE_DEFAULT_INDEX);
        self.cycle
            .record_values(record, values, record_indel_values)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_context_key_carries_its_length() {
        // A single A is 1: length one and base index zero.
        assert_eq!(key_from_context(b"A"), 1);
        // Two As is 2, the length again, with both base indices zero.
        assert_eq!(key_from_context(b"AA"), 2);
        // AC is 66: length two, then C's index one shifted six bits.
        assert_eq!(key_from_context(b"AC"), 66);
        // And a single non-ACGT base poisons the whole key.
        assert_eq!(key_from_context(b"AN"), -1);
        // Lower case is a base, and comes back upper.
        assert_eq!(key_from_context(b"ac"), 66);
        assert_eq!(context_from_key(66).unwrap(), "AC");
    }

    #[test]
    fn a_length_nibble_past_the_bases_decodes_to_dots_and_then_to_a() {
        // Key 4095 is fifteen ones in the length nibble and ones everywhere above, so the first
        // four bases read as T and the rest out of bits that were never written.
        assert_eq!(context_from_key(4095).unwrap(), "TTTTAAAAAAAAAAA");
        assert_eq!(context_from_key(0).unwrap(), "");
        assert_eq!(
            context_from_key(-1),
            Err(CovariateError::NegativeContextKey)
        );
    }

    #[test]
    fn the_cycle_sign_is_the_low_bit() {
        assert_eq!(key_from_cycle(3, 500).unwrap(), 6);
        assert_eq!(key_from_cycle(-3, 500).unwrap(), 7);
        assert_eq!(cycle_from_key(6), 3);
        assert_eq!(cycle_from_key(7), -3);
        // Key 1 is negative zero, which is zero, so it decodes to the same cycle as key 0.
        assert_eq!(cycle_from_key(1), 0);
        assert_eq!(cycle_from_key(0), 0);
    }

    #[test]
    fn the_cycle_refusal_reports_the_absolute_value() {
        let error = key_from_cycle(-501, 500).unwrap_err();
        assert!(
            error.message().contains("a larger cycle (501)"),
            "{}",
            error.message()
        );
        // Two spaces after the full stop, which is the reference's string.
        assert!(error.message().contains("was detected.  Please use"));
    }

    #[test]
    fn a_context_size_past_thirteen_names_the_argument_that_was_too_big() {
        let arguments = RecalibrationArguments {
            mismatches_context_size: 14,
            ..RecalibrationArguments::default()
        };
        assert_eq!(
            ContextCovariate::new(&arguments).unwrap_err().message(),
            "Argument mismatches_context_size has a bad value: context size cannot be bigger than 13, but was 14"
        );
        // And the positivity check names both sizes, not just the one that failed.
        let zero = RecalibrationArguments {
            mismatches_context_size: 0,
            ..RecalibrationArguments::default()
        };
        assert_eq!(
            ContextCovariate::new(&zero).unwrap_err().message(),
            "Context size must be positive. Mismatches: 0 Indels: 3"
        );
    }

    #[test]
    fn an_unknown_read_group_is_minus_one_going_in_and_an_error_coming_out() {
        let covariate = ReadGroupCovariate::new(&["unit-rg1".to_string(), "unit-rg2".to_string()]);
        assert_eq!(covariate.key_from_value("unit-rg2"), 1);
        assert_eq!(covariate.key_from_value("nonesuch"), MISSING_READ_GROUP_KEY);
        assert_eq!(
            covariate.format_key(99).unwrap_err().message(),
            "missing key 99"
        );
        assert_eq!(covariate.maximum_key_value(), 1);
    }

    #[test]
    fn duplicate_read_groups_get_one_key() {
        // Two groups sharing a platform unit share a key, which is what identifying by PU means.
        let covariate =
            ReadGroupCovariate::new(&["unit".to_string(), "unit".to_string(), "other".to_string()]);
        assert_eq!(covariate.key_from_value("unit"), 0);
        assert_eq!(covariate.key_from_value("other"), 1);
        assert_eq!(covariate.maximum_key_value(), 1);
    }

    #[test]
    fn the_context_list_is_shorter_than_the_read_when_the_read_is_shorter_than_the_context() {
        // Two bases and a context of three: the early return leaves one entry per base and no key.
        let keys = read_context_at_each_position(b"AC", 3, create_mask(3));
        assert_eq!(keys, vec![-1, -1]);
        // And with a context of two the second position gets a real key.
        let keys = read_context_at_each_position(b"AC", 2, create_mask(2));
        assert_eq!(keys, vec![-1, 66]);
    }
}
