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
//! # The order is by sensitivity, then cut, then reversed
//!
//! `readTranches` sorts by `targetTruthSensitivity`; `onTraversalStart` keeps those at or above the
//! level and reverses, so the list runs from the least sensitive tranche to the most. A level above
//! every tranche is a refusal rather than an empty filter set, and the mutual exclusion of the two
//! cutoff arguments is checked **after** the file has been read, so a broken file is reported in
//! preference to a contradictory command line.

use gatk_engine::tranches::TruthSensitivityTranche;

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
