//! `BaseRecalibrationEngine`, ported from
//! `org.broadinstitute.hellbender.utils.recalibration` (GATK 4.6.2.0).
//!
//! BQSR's counting pass, and the other half of the cycle [`crate::bqsr_transformer`] closes. For
//! every read: mark which bases disagree with the reference, weigh each disagreement by the BAQ
//! array, skip the known sites, and increment the datums the covariates point at.
//!
//! # BAQ is off by default
//!
//! ```java
//! final byte[] baqArray = (nErrors == 0 || !recalArgs.enableBAQ) ? flatBAQArray(read) : calculateBAQArray(read, refDS);
//! ```
//!
//! `enableBAQ` is **false** by default, so the ordinary run uses a flat array of the constant 64 and
//! the whole hidden Markov model is skipped. It is also skipped when the read has no errors at all,
//! whatever the flag says, "for efficiency reasons": for Illumina data about 85% of reads take that
//! branch.
//!
//! # An indel is marked at a different base on each strand
//!
//! A deletion marks `readPos - 1` on a forward read and `readPos` on a reverse one; an insertion
//! marks `readPos - 1` **before** advancing on a forward read and `readPos` **after** advancing on a
//! reverse one. Both are the base on the far side of the event, and both are clamped away rather
//! than wrapping when the event is at an end: `1D9M` and `1I9M` mark nothing, while `9M1D` and
//! `9M1I` mark the last base.
//!
//! # Only substitutions are counted, by default
//!
//! ```java
//! cachedEventTypes = recalArgs.computeIndelBQSRTables ? EventType.values() : new EventType[]{EventType.BASE_SUBSTITUTION};
//! ```
//!
//! `computeIndelBQSRTables` is a **hidden** argument and its default is false, so the insertion and
//! deletion tables exist for the whole run and stay empty. Every datum of a default run has event
//! index zero, which the golden shows: 46 datums and not one of them an indel. A port that looped
//! over all three event types would fill tables the reference never writes.
//!
//! # And the read group table is not counted, it is collapsed
//!
//! `finalizeData` marginalises the quality score table over the reported quality, which is the only
//! place [`crate::recal_datum::RecalDatum::combine`] runs in BQSR. That is why a read group's
//! reported quality is an **estimate** recomputed from expected errors rather than an average, and
//! why the NaN that `combine` can produce is reachable from here.

use std::cell::RefCell;
use std::rc::Rc;

use htsjdk_bam::cigar::Op;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

use crate::baq::Baq;
use crate::covariates::{
    CovariateError, PerReadCovariateMatrix, RecalibrationArguments, StandardCovariateList,
    BASE_QUALITY_COVARIATE_DEFAULT_INDEX, NUM_REQUIRED_COVARITES,
    READ_GROUP_COVARIATE_DEFAULT_INDEX,
};
use crate::interval::SimpleInterval;
use crate::qual_quantizer::MIN_USABLE_Q_SCORE;
use crate::recal_datum::{EventType, RecalDatum, RecalDatumError};
use crate::recalibration_tables::{NestedArrayError, RecalibrationTables, SharedDatum};

/// `BAQ.NO_BAQ_UNCERTAINTY`, the value a flat BAQ array is filled with.
pub const NO_BAQ_UNCERTAINTY: u8 = 64;

/// `ReadUtils.DEFAULT_INSERTION_DELETION_QUAL`, used when a read carries no indel qualities.
pub const DEFAULT_INSERTION_DELETION_QUAL: u8 = 45;

/// `RecalUtils.NUMBER_ERRORS_DECIMAL_PLACES`.
pub const NUMBER_ERRORS_DECIMAL_PLACES: i32 = 2;

/// `RecalUtils.REPORTED_QUALITY_DECIMAL_PLACES`.
pub const REPORTED_QUALITY_DECIMAL_PLACES: i32 = 4;

/// The arguments this engine reads beyond the covariates'.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineArguments {
    pub covariates: RecalibrationArguments,
    /// `--enable-baq`, whose default is **false**.
    pub enable_baq: bool,
    /// `--compute-indel-bqsr-tables`, a hidden argument whose default is **false**.
    ///
    /// With it off, only `BASE_SUBSTITUTION` is counted, so the insertion and deletion tables exist
    /// and stay empty for the whole run. A port that looped over all three event types would fill
    /// them and produce a table the reference never writes.
    pub compute_indel_bqsr_tables: bool,
    /// `--preserve-qscores-less-than`: a base below this is skipped, not counted.
    pub preserve_qscores_less_than: i32,
    /// `--default-base-qualities`, negative for off.
    pub default_base_qualities: i8,
    /// `--use-original-qualities`.
    pub use_original_base_qualities: bool,
}

impl Default for EngineArguments {
    fn default() -> Self {
        EngineArguments {
            covariates: RecalibrationArguments::default(),
            enable_baq: false,
            compute_indel_bqsr_tables: false,
            preserve_qscores_less_than: MIN_USABLE_Q_SCORE,
            default_base_qualities: -1,
            use_original_base_qualities: false,
        }
    }
}

/// What the counting pass refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum EngineError {
    /// `GATKException("Array length mismatch detected. Malformed read?")`.
    ArrayLengthMismatch,
    /// `IllegalArgumentException` from rounding to zero places.
    RoundToZeroPlaces,
    /// A cigar operator the counting loop does not know.
    UnsupportedCigarOperator(String),
    Covariate(CovariateError),
    Nested(NestedArrayError),
    Datum(RecalDatumError),
}

impl EngineError {
    pub fn message(&self) -> String {
        match self {
            EngineError::ArrayLengthMismatch => {
                "Array length mismatch detected. Malformed read?".to_string()
            }
            EngineError::RoundToZeroPlaces => {
                "must round to at least one decimal place".to_string()
            }
            EngineError::UnsupportedCigarOperator(op) => {
                format!("Unsupported cigar operator: {op}")
            }
            EngineError::Covariate(error) => error.message(),
            EngineError::Nested(error) => error.message(),
            EngineError::Datum(error) => error.message(),
        }
    }
}

/// `MathUtils.roundToNDecimalPlaces(in, n)`.
///
/// **The ulp is added before the rounding**, not after, which is not the same as rounding twice. It
/// exists so that a table kept in memory equals one written to a file and read back, and the golden
/// carries the values on both sides of a half to show it.
pub fn round_to_n_decimal_places(value: f64, places: i32) -> Result<f64, EngineError> {
    if places <= 0 {
        return Err(EngineError::RoundToZeroPlaces);
    }
    let multiplier = crate::math_utils::pow10(places as f64);
    // `Math.ulp(in)`: the distance to the next representable double.
    let ulp = ulp(value);
    Ok(jmath::math::round((value + ulp) * multiplier) as f64 / multiplier)
}

/// `Math.ulp(x)`.
fn ulp(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    let x = x.abs();
    if x == 0.0 {
        return f64::from_bits(1);
    }
    let next = f64::from_bits(x.to_bits() + 1);
    next - x
}

/// `flatBAQArray(read)`: no alignment uncertainty anywhere.
pub fn flat_baq_array(length: usize) -> Vec<u8> {
    vec![NO_BAQ_UNCERTAINTY; length]
}

/// `calculateIsSNPOrIndel(read, ref, snp, isIns, isDel)`.
///
/// Fills the three arrays and answers the number of events. See the module note: an indel's mark
/// lands on a different base on each strand, and it is clamped away at the read's ends.
pub fn calculate_is_snp_or_indel(
    record: &BamRecord,
    reference_bases: &[u8],
    snp: &mut [i32],
    is_ins: &mut [i32],
    is_del: &mut [i32],
) -> Result<i32, EngineError> {
    let mut read_pos = 0usize;
    let mut ref_pos = 0usize;
    let mut events = 0;
    let reverse = crate::read::is_reverse_strand(record);

    for element in &record.cigar.elements {
        let length = element.length as usize;
        match element.op {
            Op::M | Op::Eq | Op::X => {
                for _ in 0..length {
                    // `BaseUtils.basesAreEqual`, which upper-cases both sides.
                    let equal =
                        record.read_bases[read_pos].eq_ignore_ascii_case(&reference_bases[ref_pos]);
                    let value = if equal { 0 } else { 1 };
                    snp[read_pos] = value;
                    events += value;
                    read_pos += 1;
                    ref_pos += 1;
                }
            }
            Op::D => {
                let index = if reverse {
                    read_pos as i64
                } else {
                    read_pos as i64 - 1
                };
                update_indel(is_del, index);
                ref_pos += length;
            }
            Op::N => ref_pos += length,
            Op::I => {
                // Forward marks before advancing, reverse marks after: both name the base on the
                // far side of the insertion.
                if !reverse {
                    update_indel(is_ins, read_pos as i64 - 1);
                }
                read_pos += length;
                if reverse {
                    update_indel(is_ins, read_pos as i64);
                }
            }
            // The reference context does not carry the soft-clipped bases.
            Op::S => read_pos += length,
            Op::H | Op::P => {}
        }
    }

    // Summed afterwards, because two events can mark the same base and it must count once.
    events += is_del.iter().sum::<i32>() + is_ins.iter().sum::<i32>();
    Ok(events)
}

/// `updateIndel(indel, index)`: clamped away at both ends rather than wrapping.
fn update_indel(indel: &mut [i32], index: i64) {
    if index >= 0 && (index as usize) < indel.len() {
        indel[index as usize] = 1;
    }
}

/// `calculateFractionalErrorArray(errorArray, baqArray)`.
///
/// Each error inside a block of BAQ uncertainty is spread evenly over the block, and **the block
/// starts one base before** the first uncertain one: `Math.max(0, blockStartIndex - 1)`. With a flat
/// BAQ array there are no blocks and the fractions are the marks themselves, which is what the
/// default run produces.
pub fn calculate_fractional_error_array(
    error_array: &[i32],
    baq_array: &[u8],
) -> Result<Vec<f64>, EngineError> {
    if error_array.len() != baq_array.len() {
        return Err(EngineError::ArrayLengthMismatch);
    }
    const BLOCK_START_UNSET: i64 = -1;
    let mut fractional = vec![0.0f64; baq_array.len()];
    let mut in_block = false;
    let mut block_start = BLOCK_START_UNSET;
    let mut i = 0usize;
    while i < fractional.len() {
        if baq_array[i] == NO_BAQ_UNCERTAINTY {
            if !in_block {
                fractional[i] = error_array[i] as f64;
            } else {
                store_errors_in_block(i as i64, block_start, error_array, &mut fractional);
                in_block = false;
                block_start = BLOCK_START_UNSET;
            }
        } else {
            in_block = true;
            if block_start == BLOCK_START_UNSET {
                block_start = i as i64;
            }
        }
        i += 1;
    }
    if in_block {
        // `i` has already run past the end, so the reference closes the block at `i - 1`.
        store_errors_in_block(i as i64 - 1, block_start, error_array, &mut fractional);
    }
    Ok(fractional)
}

fn store_errors_in_block(i: i64, block_start: i64, error_array: &[i32], fractional: &mut [f64]) {
    let from = (block_start - 1).max(0);
    let mut total = 0;
    for j in from..=i {
        total += error_array[j as usize];
    }
    let denominator = (i - from + 1) as f64;
    for j in from..=i {
        fractional[j as usize] = total as f64 / denominator;
    }
}

/// `calculateKnownSites(read, knownSites)`: which read offsets a known site covers.
///
/// The conversion from a reference coordinate to a read offset goes through
/// `getReadIndexForReferenceCoordinate`, whose **deletion case steps back one base**, and a start
/// past the read's length collapses the whole range to the end.
pub fn calculate_known_sites(record: &BamRecord, known_sites: &[SimpleInterval]) -> Vec<bool> {
    let read_length = record.read_bases.len();
    let mut covered = vec![false; read_length];
    let soft_start = crate::read_utils::soft_start(record);
    let soft_end = crate::read_utils::soft_end(record);

    for site in known_sites {
        if site.end < soft_start || site.start > soft_end {
            continue;
        }
        let (start_index, start_op) = crate::read_utils::read_index_for_read(record, site.start);
        let mut feature_start = if start_index == -1 { 0 } else { start_index };
        if start_op == Some(Op::D) {
            feature_start -= 1;
        }
        let (end_index, _) = crate::read_utils::read_index_for_read(record, site.end);
        let mut feature_end = if end_index == -1 {
            read_length as i32
        } else {
            end_index
        };
        if feature_start > read_length as i32 {
            feature_start = read_length as i32;
            feature_end = read_length as i32;
        }
        let from = feature_start.max(0) as usize;
        let to = (feature_end + 1).min(read_length as i32).max(0) as usize;
        for slot in covered.iter_mut().take(to).skip(from) {
            *slot = true;
        }
    }
    covered
}

/// `calculateSkipArray(read, knownSites)`: a base is skipped when it is not a regular base, when its
/// quality is below the preserve threshold, or when a known site covers it.
pub fn calculate_skip_array(
    record: &BamRecord,
    known_sites: &[SimpleInterval],
    preserve_qscores_less_than: i32,
) -> Vec<bool> {
    let known = calculate_known_sites(record, known_sites);
    (0..record.read_bases.len())
        .map(|i| {
            let base = record.read_bases[i];
            let regular = matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
            let quality = record.base_qualities.get(i).copied().unwrap_or(0) as i32;
            !regular || quality < preserve_qscores_less_than || known[i]
        })
        .collect()
}

/// `BaseRecalibrationEngine`.
pub struct BaseRecalibrationEngine {
    pub arguments: EngineArguments,
    pub covariates: StandardCovariateList,
    pub tables: RecalibrationTables,
    baq: Baq,
    finalized: bool,
    num_reads_processed: u64,
}

impl BaseRecalibrationEngine {
    pub fn new(
        arguments: EngineArguments,
        header: &SamHeader,
    ) -> Result<BaseRecalibrationEngine, EngineError> {
        let covariates = StandardCovariateList::from_header(&arguments.covariates, header)
            .map_err(EngineError::Covariate)?;
        let tables = RecalibrationTables::new(&covariates).map_err(EngineError::Nested)?;
        Ok(BaseRecalibrationEngine {
            arguments,
            covariates,
            tables,
            baq: Baq::default(),
            finalized: false,
            num_reads_processed: 0,
        })
    }

    pub fn num_reads_processed(&self) -> u64 {
        self.num_reads_processed
    }

    /// `processRead(read, refDS, knownSites)`.
    ///
    /// The read is transformed first: the cigar is consolidated, default qualities are filled in,
    /// original qualities are restored, and then the adaptor **and the soft clips** are hard clipped
    /// away. The read the covariates see is not the read in the file.
    pub fn process_read(
        &mut self,
        original: &BamRecord,
        header: &SamHeader,
        contig_bases: &[u8],
        known_sites: &[SimpleInterval],
    ) -> Result<(), EngineError> {
        let read = self.transform(original, header)?;
        if read.read_bases.is_empty() {
            // The whole read was inside the adaptor.
            return Ok(());
        }

        let length = read.read_bases.len();
        let mut is_snp = vec![0i32; length];
        let mut is_insertion = vec![0i32; length];
        let mut is_deletion = vec![0i32; length];
        // `queryAndPrefetch(read.getContig(), read.getStart(), read.getEnd())`: the read's own span.
        let span_start = crate::read_utils::start(&read).max(1) as usize - 1;
        let span_end = (crate::read_utils::end(&read) as usize).min(contig_bases.len());
        let reference_bases = &contig_bases[span_start..span_end];
        let n_errors = calculate_is_snp_or_indel(
            &read,
            reference_bases,
            &mut is_snp,
            &mut is_insertion,
            &mut is_deletion,
        )?;

        // See the module note: the model runs only when there is something to marginalise over AND
        // the flag is on, and the flag is off by default.
        let baq_array = if n_errors == 0 || !self.arguments.enable_baq {
            flat_baq_array(length)
        } else {
            match self.calculate_baq_array(&read, contig_bases) {
                Some(array) => array,
                // "some reads just can't be BAQ'ed". The read is still counted as processed and
                // then dropped, because the reference's `if (baqArray != null)` guards only the
                // counting and `numReadsProcessed++` is outside it.
                None => {
                    self.num_reads_processed += 1;
                    return Ok(());
                }
            }
        };

        let mut matrix = PerReadCovariateMatrix::new(length, self.covariates.size());
        self.covariates
            .populate_per_read_covariate_matrix(&read, header, &mut matrix, true)
            .map_err(EngineError::Covariate)?;
        let skip = calculate_skip_array(
            &read,
            known_sites,
            self.arguments.preserve_qscores_less_than,
        );
        let snp_errors = calculate_fractional_error_array(&is_snp, &baq_array)?;
        let insertion_errors = calculate_fractional_error_array(&is_insertion, &baq_array)?;
        let deletion_errors = calculate_fractional_error_array(&is_deletion, &baq_array)?;

        self.update_tables(
            &read,
            &matrix,
            &skip,
            &snp_errors,
            &insertion_errors,
            &deletion_errors,
        )?;
        self.num_reads_processed += 1;
        Ok(())
    }

    /// `calculateBAQArray(read, refDS)`, which is **not** the BAQ array.
    ///
    /// ```java
    /// baq.baqRead(read, refDS, BAQ.CalculationMode.RECALCULATE, BAQ.QualityMode.ADD_TAG);
    /// return BAQ.getBAQTag(read);
    /// ```
    ///
    /// `ADD_TAG` writes the **encoded** tag, `(char)(quality - baq + 64)`, and `getBAQTag` reads
    /// those characters straight back out. So what reaches `calculateFractionalErrorArray` is the
    /// tag's encoding and not the model's qualities, and `NO_BAQ_UNCERTAINTY` being 64 is exactly
    /// the character for "the BAQ equals the quality". A port that passed the model's `bq` array
    /// would detect completely different blocks.
    ///
    /// The model runs over the **reference window**, not the read's span, and answers nothing when
    /// that window runs past the contig, which is the "some reads just can't be BAQ'ed" case.
    fn calculate_baq_array(&self, read: &BamRecord, contig_bases: &[u8]) -> Option<Vec<u8>> {
        // `excludeReadFromBAQ`: unmapped, vendor-failed or duplicate reads never get a tag, so the
        // caller sees a null and drops the read.
        if crate::read::is_unmapped(read)
            || crate::read::fails_vendor_quality_check(read)
            || crate::read::is_duplicate(read)
        {
            return None;
        }
        let window = crate::baq::reference_window_for_read(read, "", self.baq.band_width())?;
        if window.end as usize > contig_bases.len() {
            return None;
        }
        let from = window.start.max(1) as usize - 1;
        let to = window.end as usize;
        let offset = window.start - crate::read_utils::start(read);
        let result = self
            .baq
            .calc_baq_from_hmm(read, &contig_bases[from..to], offset)
            .ok()??;
        Some(crate::baq::encode_bq_tag(&read.base_qualities, &result.bq).into_bytes())
    }

    /// `makeReadTransform()`: four transforms in this order.
    fn transform(&self, read: &BamRecord, header: &SamHeader) -> Result<BamRecord, EngineError> {
        let mut out = read.clone();
        // `consolidateCigar`, which collapses zero-length and repeated elements.
        // `new CigarBuilder()`, whose no-argument form keeps deletions at the ends.
        let mut builder = crate::cigar_builder::CigarBuilder::new(false);
        let mut consolidated = Ok(out.cigar.clone());
        for element in &out.cigar.elements {
            if builder.add(*element).is_err() {
                consolidated = Err(());
                break;
            }
        }
        if consolidated.is_ok() {
            if let Ok(cigar) = builder.make(true) {
                out.cigar = cigar;
            }
        }
        if self.arguments.default_base_qualities >= 0
            && out.base_qualities.len() < out.read_bases.len()
        {
            out.base_qualities =
                vec![self.arguments.default_base_qualities as u8; out.read_bases.len()];
        }
        let clipped = crate::clipping::hard_clip_adaptor_sequence(&out, Some(header))
            .and_then(|read| crate::clipping::hard_clip_soft_clipped_bases(&read, Some(header), 0));
        Ok(clipped.unwrap_or(out))
    }

    /// `updateRecalTablesForRead`: only the quality score table and the additional ones are counted.
    /// The read group table is filled by [`BaseRecalibrationEngine::finalize_data`].
    fn update_tables(
        &mut self,
        read: &BamRecord,
        matrix: &PerReadCovariateMatrix,
        skip: &[bool],
        snp_errors: &[f64],
        insertion_errors: &[f64],
        deletion_errors: &[f64],
    ) -> Result<(), EngineError> {
        let insertion_quals = crate::read_utils::base_insertion_qualities(read);
        let deletion_quals = crate::read_utils::base_deletion_qualities(read);

        for offset in 0..read.read_bases.len() {
            if skip[offset] {
                continue;
            }
            // `cachedEventTypes`, which is one element unless the hidden indel flag is set.
            let events: &[EventType] = if self.arguments.compute_indel_bqsr_tables {
                &EventType::VALUES
            } else {
                &[EventType::BaseSubstitution]
            };
            for event in events.iter().copied() {
                let keys = matrix.covariates_at_offset(offset, event);
                let event_index = event.ordinal() as i32;
                let (qual, is_error) = match event {
                    EventType::BaseSubstitution => {
                        (read.base_qualities[offset], snp_errors[offset])
                    }
                    EventType::BaseInsertion => (insertion_quals[offset], insertion_errors[offset]),
                    EventType::BaseDeletion => (deletion_quals[offset], deletion_errors[offset]),
                };
                let read_group = keys[READ_GROUP_COVARIATE_DEFAULT_INDEX];
                let base_quality = keys[BASE_QUALITY_COVARIATE_DEFAULT_INDEX];

                increment_datum(
                    &mut self.tables.all_tables[1],
                    qual,
                    is_error,
                    &[read_group, base_quality, event_index],
                )?;

                // Indexed rather than iterated because the index IS the table this covariate
                // belongs to, which is what `recalTables.getTable(i)` means.
                #[allow(clippy::needless_range_loop)]
                for i in NUM_REQUIRED_COVARITES..self.covariates.size() {
                    let special = keys[i];
                    if special >= 0 {
                        increment_datum(
                            &mut self.tables.all_tables[i],
                            qual,
                            is_error,
                            &[read_group, base_quality, special, event_index],
                        )?;
                    }
                }
            }
        }
        Ok(())
    }

    /// `finalizeData()`: collapse the quality score table into the read group table, then round.
    pub fn finalize_data(&mut self) -> Result<(), EngineError> {
        collapse_quality_score_table_to_read_group_table(&mut self.tables)?;
        round_table_values(&mut self.tables)?;
        self.finalized = true;
        Ok(())
    }

    pub fn is_finalized(&self) -> bool {
        self.finalized
    }
}

/// `RecalUtils.incrementDatum*keys`: increment where there is a datum, create one where there is
/// not. The created datum has **one** observation and the reported quality of the base.
fn increment_datum(
    table: &mut crate::recalibration_tables::NestedIntegerArray,
    qual: u8,
    is_error: f64,
    keys: &[i32],
) -> Result<(), EngineError> {
    let existing = match keys.len() {
        3 => table.get3_keys(keys[0], keys[1], keys[2]),
        4 => table.get4_keys(keys[0], keys[1], keys[2], keys[3]),
        _ => table.get(keys),
    }
    .map_err(EngineError::Nested)?;

    match existing {
        Some(datum) => datum.borrow_mut().increment(1, is_error),
        None => {
            let datum: SharedDatum = Rc::new(RefCell::new(
                RecalDatum::new(1, is_error, qual as i8).map_err(EngineError::Datum)?,
            ));
            table.put(datum, keys).map_err(EngineError::Nested)?;
        }
    }
    Ok(())
}

/// `collapseQualityScoreTableToReadGroupTable`.
///
/// The read group table is never counted into: it is marginalised out of the quality score table
/// here, which is the only place `RecalDatum::combine` runs in BQSR. That is why a read group's
/// reported quality is an estimate recomputed from expected errors rather than an average.
pub fn collapse_quality_score_table_to_read_group_table(
    tables: &mut RecalibrationTables,
) -> Result<(), EngineError> {
    let leaves = tables.quality_score_table().all_leaves();
    for (keys, qual_datum) in leaves {
        let read_group = keys[0];
        let event = keys[2];
        let existing = tables
            .read_group_table()
            .get(&[read_group, event])
            .map_err(EngineError::Nested)?;
        match existing {
            None => {
                // `new RecalDatum(qualDatum)`, a copy: the two tables do not share this one.
                let copy = qual_datum.borrow().clone();
                tables
                    .read_group_table_mut()
                    .put(Rc::new(RefCell::new(copy)), &[read_group, event])
                    .map_err(EngineError::Nested)?;
            }
            Some(datum) => {
                let other = qual_datum.borrow().clone();
                datum
                    .borrow_mut()
                    .combine(&other)
                    .map_err(EngineError::Datum)?;
            }
        }
    }
    Ok(())
}

/// `roundTableValues(rt)`: trim every datum to what a file would hold.
///
/// The empirical quality is an integer and needs no rounding, which the reference's own comment
/// says.
pub fn round_table_values(tables: &mut RecalibrationTables) -> Result<(), EngineError> {
    for table in &tables.all_tables {
        for (_, datum) in table.all_leaves() {
            let mismatches = round_to_n_decimal_places(
                datum.borrow().num_mismatches(),
                NUMBER_ERRORS_DECIMAL_PLACES,
            )?;
            let reported = round_to_n_decimal_places(
                datum.borrow().reported_quality(),
                REPORTED_QUALITY_DECIMAL_PLACES,
            )?;
            datum
                .borrow_mut()
                .set_num_mismatches(mismatches)
                .map_err(EngineError::Datum)?;
            datum
                .borrow_mut()
                .set_reported_quality(reported)
                .map_err(EngineError::Datum)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_indel_mark_is_clamped_away_at_the_ends() {
        let mut is_del = vec![0i32; 4];
        update_indel(&mut is_del, -1);
        assert_eq!(is_del, vec![0, 0, 0, 0]);
        update_indel(&mut is_del, 4);
        assert_eq!(is_del, vec![0, 0, 0, 0]);
        update_indel(&mut is_del, 3);
        assert_eq!(is_del, vec![0, 0, 0, 1]);
    }

    #[test]
    fn a_flat_baq_array_leaves_the_marks_alone() {
        let errors = vec![0, 1, 0, 0];
        let flat = flat_baq_array(4);
        assert_eq!(
            calculate_fractional_error_array(&errors, &flat).unwrap(),
            vec![0.0, 1.0, 0.0, 0.0]
        );
    }

    #[test]
    fn a_block_reaches_one_base_before_it_and_one_after_it() {
        // The uncertain run is at indexes 1 and 2. The block starts one base EARLIER, at 0, and the
        // base that closes it, at 3, is inside it too: the store runs from `blockStartIndex - 1` to
        // the index that ended the block. So the one error is spread over four bases and not three,
        // and nothing outside the read is left at zero.
        let errors = vec![1, 0, 0, 0];
        let baq = vec![64, 60, 60, 64];
        let out = calculate_fractional_error_array(&errors, &baq).unwrap();
        assert_eq!(out, vec![0.25, 0.25, 0.25, 0.25]);
    }

    #[test]
    fn the_lengths_must_agree() {
        assert_eq!(
            calculate_fractional_error_array(&[0; 3], &[64; 4])
                .unwrap_err()
                .message(),
            "Array length mismatch detected. Malformed read?"
        );
    }

    #[test]
    fn rounding_adds_the_ulp_before_it_rounds() {
        assert_eq!(round_to_n_decimal_places(1.005, 2).unwrap(), 1.01);
        assert_eq!(round_to_n_decimal_places(0.0, 2).unwrap(), 0.0);
        assert_eq!(
            round_to_n_decimal_places(1.0, 0).unwrap_err().message(),
            "must round to at least one decimal place"
        );
    }
}
