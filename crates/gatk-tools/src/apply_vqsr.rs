//! `ApplyVQSR`, ported from `org.broadinstitute.hellbender.tools.walkers.vqsr.ApplyVQSR`
//! (GATK 4.6.2.0).
//!
//! This module is `onTraversalStart` alone: the tranches file read by
//! [`gatk_engine::tranches`], cut at `--truth-sensitivity-filter-level`, and turned into FILTER
//! header lines. The per-record filtering is a slice of its own.
//!
//! # The last tranche never becomes a filter, and the first becomes two
//!
//! ```java
//! if( tranches.size() >= 2 ) {
//!     for( int iii = 0; iii < tranches.size() - 1; iii++ ) {
//!         final TruthSensitivityTranche t = tranches.get(iii);
//!         hInfo.add(new VCFFilterHeaderLine(t.name, ... + t.minVQSLod + " <= x < " + tranches.get(iii+1).minVQSLod));
//!     }
//! }
//! if( tranches.size() >= 1 ) {
//!     hInfo.add(new VCFFilterHeaderLine(tranches.get(0).name + "+", ... + tranches.get(0).minVQSLod));
//! }
//! ```
//!
//! The loop stops one short, so the tranche whose variants are all kept has no line at all, and the
//! first one after the reversal gets two: one interval line bounded above by the **next** tranche's
//! `minVQSLod`, and one open-ended line under its own name with a `+` appended. A file of three
//! tranches therefore produces three lines naming two tranches.
//!
//! # The IDs are the file's, not the tool's
//!
//! `t.name` is the `filterName` column, so a tranches file naming its tranches `loose`, `middling`
//! and `tight` produces FILTER IDs `middling`, `tight` and `tight+`. Nothing is synthesized from the
//! model or the sensitivities at this point, unlike the string `VariantRecalibrator` writes into the
//! column in the first place.
//!
//! # One LOD, one band, one filter
//!
//! [`filter_string`] walks the kept tranches backwards and answers `PASS` only for the last of them,
//! the one whose variants are all kept, so the bands are exactly the intervals the header lines
//! describe. [`site_specific_filtering`] is what feeds it: the recal record must agree with the
//! variant on **both** ends, by two different mechanisms, and the LOD is reparsed and written back
//! through htsjdk's own double format, which sends every negative value to scientific notation.
//!
//! # One entry per alternate allele, whatever the mode
//!
//! With `-AS` the three lists are always as long as the alternate allele list:
//! [`allele_specific_filtering`] **pads** an allele of the other mode and a spanning deletion rather
//! than skipping either, and the spanning deletion is padded without its recal record ever being
//! looked for. The class test behind that is a length comparison, so `*` counts as a SNP.
//! [`site_filter_for_a_single_mode`] then leaves a mixed site unfiltered, because the other mode has
//! not run yet.
//!
//! # The order is by sensitivity, then cut, then reversed
//!
//! `readTranches` sorts by `targetTruthSensitivity`; `onTraversalStart` keeps those at or above the
//! level and reverses, so the list runs from the least sensitive tranche to the most. A level above
//! every tranche is a refusal rather than an empty filter set, and the mutual exclusion of the two
//! cutoff arguments is checked **after** the file has been read, so a broken file is reported in
//! preference to a contradictory command line.

use gatk_engine::tranches::TruthSensitivityTranche;
use htsjdk_vcf::variant::format_vcf_double;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK ApplyVQSR";

/// `LOW_VQSLOD_FILTER_NAME`.
pub const LOW_VQSLOD_FILTER_NAME: &str = "LOW_VQSLOD";

/// `DEFAULT_VQSLOD_CUTOFF`.
pub const DEFAULT_VQSLOD_CUTOFF: f64 = 0.0;

/// What `onTraversalStart` refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum ApplyVqsrError {
    /// A level no tranche reaches.
    NoTranches(f64),
    /// Both cutoffs at once, checked only after the tranches file has been read.
    MutuallyExclusiveCutoffs,
}

impl ApplyVqsrError {
    /// A plain `UserException` either way.
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
        match self {
            ApplyVqsrError::NoTranches(level) => format!(
                "No tranches were found in the file or were above the truth sensitivity filter level {}",
                gatk_engine::tsv_table::java_double_to_string(*level)
            ),
            ApplyVqsrError::MutuallyExclusiveCutoffs => "Arguments --truth-sensitivity-filter-level and --lod-score-cutoff are mutually exclusive. Please only specify one option.".to_string(),
        }
    }
}

/// One `##FILTER` line, as its two fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilterHeaderLine {
    pub id: String,
    pub description: String,
}

impl FilterHeaderLine {
    /// The line as htsjdk writes it.
    pub fn to_line(&self) -> String {
        format!(
            "##FILTER=<ID={},Description=\"{}\">",
            self.id, self.description
        )
    }
}

/// The tranches kept by a level: at or above it, then reversed.
pub fn keep(tranches: &[TruthSensitivityTranche], level: f64) -> Vec<&TruthSensitivityTranche> {
    let mut kept: Vec<&TruthSensitivityTranche> = tranches
        .iter()
        .filter(|tranche| tranche.target_truth_sensitivity >= level)
        .collect();
    kept.reverse();
    kept
}

/// The FILTER lines of a run given `--truth-sensitivity-filter-level`.
///
/// `kept` is what [`keep`] returned, which is why the emptiness check is here rather than there:
/// the reference refuses inside `onTraversalStart`, after the header lines would have been built.
pub fn tranche_filter_lines(
    kept: &[&TruthSensitivityTranche],
    level: f64,
) -> Result<Vec<FilterHeaderLine>, ApplyVqsrError> {
    let mut lines = Vec::new();
    if kept.len() >= 2 {
        // One short of the end: the last tranche is the one kept whole and named by nothing.
        for (index, tranche) in kept.iter().enumerate().take(kept.len() - 1) {
            lines.push(FilterHeaderLine {
                id: tranche.name.clone(),
                description: format!(
                    "Truth sensitivity tranche level for {} model at VQS Lod: {} <= x < {}",
                    tranche.model.name(),
                    lod(tranche.min_vqslod),
                    lod(kept[index + 1].min_vqslod)
                ),
            });
        }
    }
    let Some(first) = kept.first() else {
        return Err(ApplyVqsrError::NoTranches(level));
    };
    lines.push(FilterHeaderLine {
        id: format!("{}+", first.name),
        description: format!(
            "Truth sensitivity tranche level for {} model at VQS Lod < {}",
            first.model.name(),
            lod(first.min_vqslod)
        ),
    });
    Ok(lines)
}

/// The one FILTER line of a run given no level at all.
pub fn low_vqslod_filter_line(cutoff: Option<f64>) -> FilterHeaderLine {
    let cutoff = cutoff.unwrap_or(DEFAULT_VQSLOD_CUTOFF);
    FilterHeaderLine {
        id: LOW_VQSLOD_FILTER_NAME.to_string(),
        description: format!("VQSLOD < {}", lod(cutoff)),
    }
}

/// A LOD in a description, which is string concatenation and therefore `Double.toString`.
fn lod(value: f64) -> String {
    gatk_engine::tsv_table::java_double_to_string(value)
}

/// `VQS_LOD_KEY`, `CULPRIT_KEY` and the two training labels.
pub const VQS_LOD_KEY: &str = "VQSLOD";
pub const CULPRIT_KEY: &str = "culprit";
pub const POSITIVE_LABEL_KEY: &str = "POSITIVE_TRAIN_SITE";
pub const NEGATIVE_LABEL_KEY: &str = "NEGATIVE_TRAIN_SITE";

/// What `doSiteSpecificFiltering` refuses, each of them quoting the whole record.
///
/// The reference builds these as `UserException` and the walker rethrows each one as a
/// `GATKException` naming the locus, so the wording below is the **cause**'s. `record` is the
/// `VariantContext.toString()` the message ends with, which is not ported: a caller supplies it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SiteFilteringError {
    /// Nothing in the recal file starts and ends where this record does.
    NoRecalRecord { record: String },
    /// A recal record with no `VQSLOD` attribute at all.
    NoLod { record: String },
    /// A `VQSLOD` that `Double.valueOf` will not take.
    UnreadableLod { record: String },
    /// An alternate allele with no recal record of its own, which is the `-AS` wording.
    NoRecalAllele { record: String },
}

impl SiteFilteringError {
    /// The class of the cause, before the walker wraps it.
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
        match self {
            SiteFilteringError::NoRecalRecord { record } => format!(
                "Encountered input variant which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyVQSR were run on the same set of input variants. First seen at: {record}"
            ),
            SiteFilteringError::NoLod { record } => format!(
                "Encountered a malformed record in the input recal file. There is no lod for the record at: {record}"
            ),
            SiteFilteringError::UnreadableLod { record } => format!(
                "Encountered a malformed record in the input recal file. The lod is unreadable for the record at: {record}"
            ),
            SiteFilteringError::NoRecalAllele { record } => format!(
                "Encountered input allele which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyVQSR were run on the same set of input variants with flag -AS. First seen at: {record}"
            ),
        }
    }

    /// The class the walker rethrows it as, which is not a user error any more.
    pub fn wrapper_class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.GATKException"
    }

    /// The wrapper's message, whose own wording holds the locus and the record rather than the
    /// tool's explanation.
    pub fn wrapper_message(&self, contig: &str, start: i32, record: &str) -> String {
        format!("Exception thrown at {contig}:{start} {record}")
    }
}

/// As much of a recal record as the pairing and the filtering look at.
pub trait RecalRecord {
    fn start(&self) -> i32;
    fn end(&self) -> i32;
    /// `getAttributeAsString(VQS_LOD_KEY, null)`: the text, not a number.
    fn lod_string(&self) -> Option<String>;
    fn culprit(&self) -> Option<String>;
    fn has_positive_label(&self) -> bool;
    fn has_negative_label(&self) -> bool;
}

/// The recal record a variant is paired with, or nothing.
///
/// Two mechanisms, not one: `featureContext.getValues(recal, vc.getStart())` has already dropped
/// every record that does not **start** where the variant does, and this takes the first of what is
/// left that **ends** where it does. Both halves are here so that a caller passing the whole file
/// gets the reference's answer rather than the first record at the locus.
pub fn matching_recal<T: RecalRecord>(start: i32, end: i32, recals: &[T]) -> Option<&T> {
    recals
        .iter()
        .filter(|recal| recal.start() == start)
        .find(|recal| recal.end() == end)
}

/// What `doSiteSpecificFiltering` writes onto the record, beside the filter it returns.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    /// `VQSLOD`, rendered as htsjdk renders a double: scientific for anything below `0.01`,
    /// **including every negative value**, because the branch is on the signed number.
    pub vqslod: String,
    /// `culprit`, which is `.` when the recal record carried none.
    pub culprit: String,
    pub positive_label: bool,
    pub negative_label: bool,
}

/// `doSiteSpecificFiltering`: the pairing, the LOD, the annotations and the filter.
pub fn site_specific_filtering<T: RecalRecord>(
    start: i32,
    end: i32,
    recals: &[T],
    record: &str,
) -> Result<(Annotation, f64), SiteFilteringError> {
    let Some(recal) = matching_recal(start, end, recals) else {
        return Err(SiteFilteringError::NoRecalRecord {
            record: record.to_string(),
        });
    };
    let Some(text) = recal.lod_string() else {
        return Err(SiteFilteringError::NoLod {
            record: record.to_string(),
        });
    };
    let Ok(lod) = text.trim().parse::<f64>() else {
        return Err(SiteFilteringError::UnreadableLod {
            record: record.to_string(),
        });
    };
    Ok((
        Annotation {
            vqslod: format_vcf_double(lod),
            // `builder.attribute(CULPRIT_KEY, recalDatum.getAttribute(CULPRIT_KEY))` with no
            // default, and a null attribute is written as a missing value.
            culprit: recal.culprit().unwrap_or_else(|| ".".to_string()),
            positive_label: recal.has_positive_label(),
            negative_label: recal.has_negative_label(),
        },
        lod,
    ))
}

/// `generateFilterString(lod)` under a truth-sensitivity level.
///
/// The walk is backwards and the last tranche is the one kept whole, so it is the only one that
/// answers `PASS`; below every tranche the answer is the first tranche's name with a `+`.
pub fn filter_string(kept: &[&TruthSensitivityTranche], lod: f64) -> String {
    for (index, tranche) in kept.iter().enumerate().rev() {
        if lod >= tranche.min_vqslod {
            return if index == kept.len() - 1 {
                PASSES_FILTERS.to_string()
            } else {
                tranche.name.clone()
            };
        }
    }
    format!("{}+", kept[0].name)
}

/// `generateFilterString(lod)` under a plain LOD cutoff.
pub fn filter_string_by_cutoff(lod: f64, cutoff: f64) -> String {
    if lod < cutoff {
        LOW_VQSLOD_FILTER_NAME.to_string()
    } else {
        PASSES_FILTERS.to_string()
    }
}

/// `VCFConstants.PASSES_FILTERS_v4`.
pub const PASSES_FILTERS: &str = "PASS";

/// `apply`'s gate: whether the record is recalibrated at all.
///
/// `of_this_mode` is `VariantDataManager.checkVariationClass(vc, MODE)`. The filter half is the
/// reference's own three-way test, whose middle clause is why an unfiltered record always passes.
pub fn recalibrates(
    of_this_mode: bool,
    filters: &[String],
    ignore_all_filters: bool,
    ignored_filters: &[String],
) -> bool {
    let not_filtered = ignore_all_filters
        || filters.is_empty()
        || (!ignored_filters.is_empty()
            && filters
                .iter()
                .all(|filter| ignored_filters.contains(filter)));
    of_this_mode && not_filtered
}

/// `AS_FILTER_STATUS_KEY`, `AS_VQS_LOD_KEY` and `AS_CULPRIT_KEY`.
pub const AS_FILTER_STATUS_KEY: &str = "AS_FilterStatus";
pub const AS_VQS_LOD_KEY: &str = "AS_VQSLOD";
pub const AS_CULPRIT_KEY: &str = "AS_culprit";

/// `emptyStringValue` and `emptyFloatValue`, the two paddings.
pub const EMPTY_STRING_VALUE: &str = "NA";
pub const EMPTY_FLOAT_VALUE: &str = "NaN";

/// `AnnotationUtils.LIST_DELIMITER`.
pub const LIST_DELIMITER: &str = ",";

/// `VariantRecalibratorEngine.MIN_ACCEPTABLE_LOD_SCORE`, which is where `bestLod` starts.
pub const MIN_ACCEPTABLE_LOD_SCORE: f64 = -20000.0;

/// `VCFConstants.UNFILTERED`, the FILTER column of a site left for the other mode.
pub const UNFILTERED: &str = ".";

/// `GATKVCFConstants.isSpanningDeletion`, which knows the deprecated spelling too.
pub fn is_spanning_deletion(allele: &str) -> bool {
    allele == "*" || allele == "<*:DEL>"
}

/// `VariantDataManager.checkVariationClass(vc, allele, mode)`.
///
/// A **length** comparison rather than a type: SNP mode keeps every allele as long as the reference,
/// which is why the reference's own comment says a spanning deletion counts as a SNP. INDEL mode
/// keeps the rest, and anything symbolic.
pub fn allele_is_of_mode(reference: &str, allele: &str, snp_mode: bool) -> bool {
    // `isSymbolic()` is the angle brackets alone: the golden's INDEL-mode run pads `*` on both
    // sides, which it would not if the spanning deletion counted as symbolic here.
    let symbolic = allele.starts_with('<');
    if snp_mode {
        reference.len() == allele.len()
    } else {
        reference.len() != allele.len() || symbolic
    }
}

/// The three lists `doAlleleSpecificFiltering` writes, and the best LOD it saw.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleleAnnotations {
    pub filter_status: Vec<String>,
    pub vqslod: Vec<String>,
    pub culprit: Vec<String>,
    /// Copied from whichever allele's recal record carried them.
    pub positive_label: bool,
    pub negative_label: bool,
}

impl AlleleAnnotations {
    /// The three attributes as they are written, each one a `,`-joined list.
    pub fn attributes(&self) -> Vec<(String, String)> {
        vec![
            (
                AS_FILTER_STATUS_KEY.to_string(),
                self.filter_status.join(LIST_DELIMITER),
            ),
            (AS_VQS_LOD_KEY.to_string(), self.vqslod.join(LIST_DELIMITER)),
            (
                AS_CULPRIT_KEY.to_string(),
                self.culprit.join(LIST_DELIMITER),
            ),
        ]
    }
}

/// What one record brings to [`allele_specific_filtering`]: its locus and its alleles.
pub struct AlleleSpecificSite<'a> {
    pub start: i32,
    pub end: i32,
    pub reference: &'a str,
    pub alternates: &'a [String],
    /// `VariantContext.toString()`, which only the refusal uses.
    pub record: &'a str,
}

/// `doAlleleSpecificFiltering` for a first run, where no previous mode has left anything behind.
///
/// One entry per alternate allele, whatever mode is running: an allele of the other mode and a
/// spanning deletion are **padded** rather than skipped, and the spanning deletion is padded without
/// its recal record ever being looked for. The LOD is written with `%.4f`, and nothing site-level is
/// written at all: in this mode there is no `VQSLOD` and no `culprit`.
pub fn allele_specific_filtering<T: RecalRecord + AllelicRecalRecord>(
    site: &AlleleSpecificSite,
    recals: &[T],
    snp_mode: bool,
    kept: &[&TruthSensitivityTranche],
) -> Result<(AlleleAnnotations, f64), SiteFilteringError> {
    let AlleleSpecificSite {
        start,
        end,
        reference,
        alternates,
        record,
    } = *site;
    let mut best_lod = MIN_ACCEPTABLE_LOD_SCORE;
    let mut annotations = AlleleAnnotations {
        filter_status: Vec::new(),
        vqslod: Vec::new(),
        culprit: Vec::new(),
        positive_label: false,
        negative_label: false,
    };

    for allele in alternates {
        if !allele_is_of_mode(reference, allele, snp_mode) {
            // The first run has no previous lists to copy from, so the padding is all there is.
            annotations
                .filter_status
                .push(EMPTY_STRING_VALUE.to_string());
            annotations.vqslod.push(EMPTY_FLOAT_VALUE.to_string());
            annotations.culprit.push(EMPTY_STRING_VALUE.to_string());
            continue;
        }
        if is_spanning_deletion(allele) {
            // Kept by the class test and then padded anyway, with no lookup at all.
            annotations
                .filter_status
                .push(EMPTY_STRING_VALUE.to_string());
            annotations.vqslod.push(EMPTY_FLOAT_VALUE.to_string());
            annotations.culprit.push(EMPTY_STRING_VALUE.to_string());
            continue;
        }
        let Some(recal) = matching_allele_recal(start, end, allele, recals) else {
            return Err(SiteFilteringError::NoRecalAllele {
                record: record.to_string(),
            });
        };
        // `getAttributeAsDouble(VQS_LOD_KEY, MIN_ACCEPTABLE_LOD_SCORE)`: a missing key is the floor
        // here rather than the refusal the site-level path makes of it.
        let lod = recal
            .lod_string()
            .and_then(|text| text.trim().parse::<f64>().ok())
            .unwrap_or(MIN_ACCEPTABLE_LOD_SCORE);
        if lod > best_lod {
            best_lod = lod;
        }
        annotations.filter_status.push(filter_string(kept, lod));
        annotations.vqslod.push(format!("{lod:.4}"));
        annotations
            .culprit
            .push(recal.culprit().unwrap_or_else(|| ".".to_string()));
        annotations.positive_label |= recal.has_positive_label();
        annotations.negative_label |= recal.has_negative_label();
    }
    Ok((annotations, best_lod))
}

/// `getMatchingRecalVC` in allele-specific mode, where the allele is compared too.
pub fn matching_allele_recal<'a, T: RecalRecord + AllelicRecalRecord>(
    start: i32,
    end: i32,
    allele: &str,
    recals: &'a [T],
) -> Option<&'a T> {
    recals
        .iter()
        .filter(|recal| recal.start() == start)
        .find(|recal| recal.end() == end && recal.first_alternate() == allele)
}

/// The one thing a recal record needs beyond [`RecalRecord`] for allele-specific matching.
pub trait AllelicRecalRecord {
    /// `getAlternateAllele(0)`, which is the allele the record was written for.
    fn first_alternate(&self) -> String;
}

/// `generateFilterStringFromAlleles` for a first run, where neither mode has left a filter behind.
///
/// `both_modes_were_run` is false by construction here, so the site is left **unfiltered** unless
/// `onlyOneModeNeeded` holds: the record must not be mixed and must be of the running mode. A mixed
/// site therefore comes out with a FILTER of `.` while its `AS_FilterStatus` already names a tranche.
/// The second run, which reads the first run's filter names back out of the header, is a slice of
/// its own.
pub fn site_filter_for_a_single_mode(
    is_mixed: bool,
    record_is_of_mode: bool,
    best_lod: f64,
    kept: &[&TruthSensitivityTranche],
) -> String {
    let only_one_mode_needed = !is_mixed && record_is_of_mode;
    if !only_one_mode_needed {
        return UNFILTERED.to_string();
    }
    filter_string(kept, best_lod)
}

/// Whether a record reaches the output.
///
/// `--exclude-filtered` is consulted **only** for a record the tool recalibrated: one emitted
/// untouched by [`recalibrates`] is written whatever its filters say.
pub fn writes_out(recalibrated: bool, filter: &str, exclude_filtered: bool) -> bool {
    if !recalibrated {
        return true;
    }
    !exclude_filtered || filter == PASSES_FILTERS
}

#[cfg(test)]
mod tests {
    use super::*;
    use gatk_engine::tranches::{read_tranches, Mode};

    const COLUMNS: &str = "targetTruthSensitivity,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity";

    fn three() -> Vec<TruthSensitivityTranche> {
        let text = format!(
            "# Variant quality score tranches file\n{COLUMNS}\n\
             100.00,30,15,1.9000,1.7000,-0.5000,VQSRTrancheSNP99.00to100.00,SNP,100,100,1.0000\n\
             90.00,10,5,2.1000,1.9000,3.5000,VQSRTrancheSNP0.00to90.00,SNP,100,90,0.9000\n\
             99.00,20,9,2.0000,1.8000,1.5000,VQSRTrancheSNP90.00to99.00,SNP,100,99,0.9900\n"
        );
        read_tranches("tranches", &text).expect("a good file")
    }

    #[test]
    fn the_last_tranche_never_becomes_a_filter_and_the_first_becomes_two() {
        let tranches = three();
        let kept = keep(&tranches, 0.0);
        // Sorted by sensitivity, then reversed: the most sensitive first.
        assert_eq!(
            kept.iter()
                .map(|tranche| tranche.target_truth_sensitivity)
                .collect::<Vec<_>>(),
            vec![100.0, 99.0, 90.0]
        );
        let lines = tranche_filter_lines(&kept, 0.0).expect("three tranches");
        assert_eq!(
            lines.iter().map(|line| line.id.clone()).collect::<Vec<_>>(),
            vec![
                "VQSRTrancheSNP99.00to100.00",
                "VQSRTrancheSNP90.00to99.00",
                "VQSRTrancheSNP99.00to100.00+"
            ]
        );
        // The 90.00 tranche, kept whole, is named by nothing.
        assert!(!lines
            .iter()
            .any(|line| line.id.starts_with("VQSRTrancheSNP0.00")));
        assert_eq!(
            lines[0].to_line(),
            "##FILTER=<ID=VQSRTrancheSNP99.00to100.00,Description=\"Truth sensitivity tranche level for SNP model at VQS Lod: -0.5 <= x < 1.5\">"
        );
        assert_eq!(
            lines[2].to_line(),
            "##FILTER=<ID=VQSRTrancheSNP99.00to100.00+,Description=\"Truth sensitivity tranche level for SNP model at VQS Lod < -0.5\">"
        );
    }

    #[test]
    fn a_level_keeps_the_tranches_at_or_above_it() {
        let tranches = three();
        let kept = keep(&tranches, 99.0);
        assert_eq!(kept.len(), 2);
        let lines = tranche_filter_lines(&kept, 99.0).expect("two tranches");
        // One interval line and the open-ended one, both naming the same tranche.
        assert_eq!(
            lines.iter().map(|line| line.id.clone()).collect::<Vec<_>>(),
            vec![
                "VQSRTrancheSNP99.00to100.00",
                "VQSRTrancheSNP99.00to100.00+"
            ]
        );
    }

    #[test]
    fn a_level_above_every_tranche_is_a_refusal() {
        let tranches = three();
        let kept = keep(&tranches, 100.1);
        assert!(kept.is_empty());
        let error = tranche_filter_lines(&kept, 100.1).expect_err("nothing survived");
        assert_eq!(
            error.message(),
            "No tranches were found in the file or were above the truth sensitivity filter level 100.1"
        );
        assert_eq!(
            error.class(),
            "org.broadinstitute.hellbender.exceptions.UserException"
        );
    }

    #[test]
    fn the_ids_are_the_files_own_names() {
        let text = format!(
            "# Variant quality score tranches file\n{COLUMNS}\n\
             99.00,20,9,2.0000,1.8000,1.5000,middling,SNP,100,99,0.9900\n\
             100.00,30,15,1.9000,1.7000,-0.5000,tight,SNP,100,100,1.0000\n"
        );
        let tranches = read_tranches("tranches", &text).expect("a good file");
        let kept = keep(&tranches, 0.0);
        let lines = tranche_filter_lines(&kept, 0.0).expect("two tranches");
        assert_eq!(
            lines.iter().map(|line| line.id.clone()).collect::<Vec<_>>(),
            vec!["tight", "tight+"]
        );
        assert_eq!(tranches[0].model, Mode::Snp);
    }

    struct Recal {
        start: i32,
        end: i32,
        /// The record's first alternate allele, which allele-specific matching compares.
        alternate: &'static str,
        lod: Option<&'static str>,
        culprit: Option<&'static str>,
        positive: bool,
        negative: bool,
    }

    impl RecalRecord for Recal {
        fn start(&self) -> i32 {
            self.start
        }
        fn end(&self) -> i32 {
            self.end
        }
        fn lod_string(&self) -> Option<String> {
            self.lod.map(|text| text.to_string())
        }
        fn culprit(&self) -> Option<String> {
            self.culprit.map(|text| text.to_string())
        }
        fn has_positive_label(&self) -> bool {
            self.positive
        }
        fn has_negative_label(&self) -> bool {
            self.negative
        }
    }

    fn recal(
        start: i32,
        end: i32,
        lod: Option<&'static str>,
        culprit: Option<&'static str>,
    ) -> Recal {
        Recal {
            start,
            end,
            alternate: "C",
            lod,
            culprit,
            positive: false,
            negative: false,
        }
    }

    #[test]
    fn the_bands_are_the_intervals_the_header_lines_describe() {
        let tranches = three();
        let kept = keep(&tranches, 0.0);
        // Above the last tranche's minVQSLod, which is the one kept whole.
        assert_eq!(filter_string(&kept, 5.0), "PASS");
        // The boundary belongs to the tranche below it.
        assert_eq!(filter_string(&kept, 3.5), "PASS");
        assert_eq!(filter_string(&kept, 2.0), "VQSRTrancheSNP90.00to99.00");
        assert_eq!(filter_string(&kept, 1.5), "VQSRTrancheSNP90.00to99.00");
        assert_eq!(filter_string(&kept, 0.0), "VQSRTrancheSNP99.00to100.00");
        assert_eq!(filter_string(&kept, -0.5), "VQSRTrancheSNP99.00to100.00");
        // Below every tranche: the first one's name, with a plus.
        assert_eq!(filter_string(&kept, -3.0), "VQSRTrancheSNP99.00to100.00+");
    }

    #[test]
    fn the_recal_record_must_agree_on_both_ends() {
        let recals = [
            // The same start, another end: skipped.
            recal(800, 805, Some("9.0"), Some("SKIPPED")),
            recal(800, 800, Some("5.0"), Some("TAKEN")),
            recal(850, 855, Some("9.0"), Some("NEVER")),
        ];
        let (annotation, lod) =
            site_specific_filtering(800, 800, &recals, "[VC]").expect("the second one agrees");
        assert_eq!(annotation.culprit, "TAKEN");
        assert_eq!(lod, 5.0);
        // Nothing at 850 ends where the input record does.
        assert_eq!(
            site_specific_filtering(850, 850, &recals, "[VC]").expect_err("no partner"),
            SiteFilteringError::NoRecalRecord {
                record: "[VC]".to_string()
            }
        );
    }

    #[test]
    fn every_negative_lod_is_written_in_scientific_notation() {
        let recals = [
            recal(400, 400, Some("-3.0"), None),
            recal(100, 100, Some("5.0"), Some("QD")),
            recal(300, 300, Some("0.0"), Some("FS")),
        ];
        let (negative, _) = site_specific_filtering(400, 400, &recals, "[VC]").expect("a lod");
        // The branch is on the signed value, so anything below 0.01 goes scientific.
        assert_eq!(negative.vqslod, "-3.000e+00");
        // And a recal record with no culprit writes a missing value rather than nothing.
        assert_eq!(negative.culprit, ".");
        let (positive, _) = site_specific_filtering(100, 100, &recals, "[VC]").expect("a lod");
        assert_eq!(positive.vqslod, "5.00");
        let (zero, _) = site_specific_filtering(300, 300, &recals, "[VC]").expect("a lod");
        assert_eq!(zero.vqslod, "0.00");
    }

    #[test]
    fn the_two_malformed_recal_records_are_told_apart() {
        let recals = [recal(900, 900, None, Some("QD"))];
        assert_eq!(
            site_specific_filtering(900, 900, &recals, "[VC]").expect_err("no lod at all"),
            SiteFilteringError::NoLod {
                record: "[VC]".to_string()
            }
        );
        let recals = [recal(900, 900, Some("nonsense"), Some("QD"))];
        let error = site_specific_filtering(900, 900, &recals, "[VC]")
            .expect_err("a lod that will not parse");
        assert!(error.message().starts_with(
            "Encountered a malformed record in the input recal file. The lod is unreadable"
        ));
        // A user error on the way out, an internal one by the time it is reported.
        assert_eq!(
            error.class(),
            "org.broadinstitute.hellbender.exceptions.UserException"
        );
        assert_eq!(
            error.wrapper_class(),
            "org.broadinstitute.hellbender.exceptions.GATKException"
        );
        assert_eq!(
            error.wrapper_message("chr1", 900, "[VC]"),
            "Exception thrown at chr1:900 [VC]"
        );
    }

    #[test]
    fn a_filtered_record_is_untouched_and_exclude_filtered_does_not_drop_it() {
        let weak = vec!["weak".to_string()];
        assert!(!recalibrates(true, &weak, false, &[]));
        assert!(recalibrates(true, &weak, true, &[]));
        assert!(recalibrates(true, &weak, false, &["weak".to_string()]));
        // An unfiltered record of the wrong class is still emitted untouched.
        assert!(!recalibrates(false, &[], false, &[]));
        assert!(recalibrates(true, &[], false, &[]));

        // The exclusion only reaches records the tool filtered itself.
        assert!(writes_out(false, "weak", true));
        assert!(!writes_out(true, "VQSRTrancheSNP99.00to100.00+", true));
        assert!(writes_out(true, "PASS", true));
    }

    impl AllelicRecalRecord for Recal {
        fn first_alternate(&self) -> String {
            self.alternate.to_string()
        }
    }

    fn allelic(
        start: i32,
        end: i32,
        alternate: &'static str,
        lod: &'static str,
        culprit: &'static str,
    ) -> Recal {
        Recal {
            start,
            end,
            alternate,
            lod: Some(lod),
            culprit: Some(culprit),
            positive: false,
            negative: false,
        }
    }

    fn alternates(list: &[&str]) -> Vec<String> {
        list.iter().map(|allele| allele.to_string()).collect()
    }

    #[test]
    fn a_spanning_deletion_counts_as_a_snp_and_is_padded_anyway() {
        // The class test is a length comparison, so `*` is kept by SNP mode...
        assert!(allele_is_of_mode("A", "*", true));
        assert!(!allele_is_of_mode("A", "*", false));
        assert!(is_spanning_deletion("*"));

        // ...and then padded without its recal record being looked for.
        let tranches = three();
        let kept = keep(&tranches, 0.0);
        let recals = [allelic(300, 300, "C", "2.0", "MQ")];
        let alleles = alternates(&["*", "C"]);
        let (annotations, best) = allele_specific_filtering(
            &AlleleSpecificSite {
                start: 300,
                end: 300,
                reference: "A",
                alternates: &alleles,
                record: "[VC]",
            },
            &recals,
            true,
            &kept,
        )
        .expect("the deletion needs nothing");
        assert_eq!(
            annotations.filter_status,
            vec!["NA", "VQSRTrancheSNP90.00to99.00"]
        );
        assert_eq!(annotations.vqslod, vec!["NaN", "2.0000"]);
        assert_eq!(annotations.culprit, vec!["NA", "MQ"]);
        assert_eq!(best, 2.0);
    }

    #[test]
    fn an_allele_of_the_other_mode_is_padded_rather_than_skipped() {
        let tranches = three();
        let kept = keep(&tranches, 0.0);
        let recals = [
            allelic(200, 200, "C", "-3.0", "FS"),
            allelic(200, 200, "ACC", "4.0", "SOR"),
        ];
        // SNP mode: the insertion is padded, the SNP is scored.
        let alleles = alternates(&["C", "ACC"]);
        let site = AlleleSpecificSite {
            start: 200,
            end: 200,
            reference: "A",
            alternates: &alleles,
            record: "[VC]",
        };
        let (snp, _) =
            allele_specific_filtering(&site, &recals, true, &kept).expect("one of the two");
        assert_eq!(
            snp.filter_status,
            vec!["VQSRTrancheSNP99.00to100.00+", "NA"]
        );
        assert_eq!(snp.vqslod, vec!["-3.0000", "NaN"]);
        // INDEL mode: the padding swaps sides, and the lists keep one entry per allele.
        let (indel, best) =
            allele_specific_filtering(&site, &recals, false, &kept).expect("the other one");
        assert_eq!(indel.filter_status, vec!["NA", "PASS"]);
        assert_eq!(indel.vqslod, vec!["NaN", "4.0000"]);
        assert_eq!(best, 4.0);
    }

    #[test]
    fn a_mixed_site_is_left_unfiltered_until_both_modes_have_run() {
        let tranches = three();
        let kept = keep(&tranches, 0.0);
        // Mixed: neither `bothModesWereRun` nor `onlyOneModeNeeded`, so the FILTER column is `.`
        // even though the allele's own status names a tranche.
        assert_eq!(site_filter_for_a_single_mode(true, true, -3.0, &kept), ".");
        // A pure indel under a SNP-mode run: not of this mode, so also unfiltered.
        assert_eq!(site_filter_for_a_single_mode(false, false, 4.0, &kept), ".");
        // A pure SNP under a SNP-mode run: the best allele's band.
        assert_eq!(
            site_filter_for_a_single_mode(false, true, 2.0, &kept),
            "VQSRTrancheSNP90.00to99.00"
        );
    }

    #[test]
    fn an_allele_with_no_recal_record_has_a_refusal_of_its_own() {
        let tranches = three();
        let kept = keep(&tranches, 0.0);
        // The right locus, the wrong allele.
        let recals = [allelic(500, 500, "G", "5.0", "QD")];
        let alleles = alternates(&["C"]);
        let error = allele_specific_filtering(
            &AlleleSpecificSite {
                start: 500,
                end: 500,
                reference: "A",
                alternates: &alleles,
                record: "[VC]",
            },
            &recals,
            true,
            &kept,
        )
        .expect_err("no record for C");
        assert!(error
            .message()
            .starts_with("Encountered input allele which isn't found in the input recal file."));
        assert!(error.message().contains("with flag -AS."));
    }

    #[test]
    fn a_plain_cutoff_has_two_answers() {
        assert_eq!(filter_string_by_cutoff(2.0, 1.0), "PASS");
        assert_eq!(filter_string_by_cutoff(0.0, 1.0), "LOW_VQSLOD");
        // The boundary passes, the test being `lod < cutoff`.
        assert_eq!(filter_string_by_cutoff(1.0, 1.0), "PASS");
    }

    #[test]
    fn with_no_level_there_is_one_line_and_a_default_cutoff() {
        assert_eq!(
            low_vqslod_filter_line(None).to_line(),
            "##FILTER=<ID=LOW_VQSLOD,Description=\"VQSLOD < 0.0\">"
        );
        assert_eq!(
            low_vqslod_filter_line(Some(0.5)).description,
            "VQSLOD < 0.5"
        );
    }
}
