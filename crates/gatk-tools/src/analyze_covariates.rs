//! `AnalyzeCovariates`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.bqsr.AnalyzeCovariates` (GATK 4.6.2.0) with
//! `RecalUtils.generateCsv`.
//!
//! The tool draws plots and the plots are R's, so what the tool itself produces is the intermediate
//! csv the R script reads. That csv is the whole observable, and this is it: ten columns per row,
//! one row per key of a table built by folding every covariate table into one.
//!
//! # The quality score table is filed under an index that is not a covariate's
//!
//! ```java
//! newCovs[1] = covariates.size(); // replace the covariate name with an arbitrary (unused)
//!                                  // index for QualityScore. This is a HACK.
//! ```
//!
//! (The trailing comment is one line in the reference, broken in two here.)
//!
//! One past the last real index, and the comment is the reference's own. Reading it back:
//!
//! ```java
//! final Covariate covariate = (covariateIndex == covariates.size())
//!         ? covariates.getQualityScoreCovariate() : covariates.get(covariateIndex);
//! ```
//!
//! (One line in the reference, broken in two here.)
//!
//! so that out-of-range index is what names the row `QualityScore`. A port that filed the quality
//! score under its own index, 1, would collide with the covariate that really lives there.
//!
//! # The optional covariates lose the quality score from their key
//!
//! ```java
//! covs[2] = leaf.keys[2];
//! ```
//!
//! `keys[1]` is the reported quality and it is skipped, so every context or cycle row is **summed
//! over the reported quality**: two data at the same context and two different qualities become
//! one row whose observations are the sum and whose average reported quality is the weighted mean.
//! The read group table is not folded in at all and never reaches the csv.
//!
//! # Reading this file without reading Rust
//!
//! - `&[T]` is a borrowed, read-only view of an array (Java's `T[]` that the callee cannot
//!   modify); `&mut T` is a writable borrow, of which only one may exist at a time;
//! - `Result<T, E>` is "a value or an error", which is how a reference exception is modelled:
//!   `Err(AnalyzeCovariatesError::NoReport)` stands for the `UserException` thrown, and `?` after
//!   a call hands an `Err` straight back to this function's caller, Java's `throws` at one site;
//! - `RoleReport<'a>` carries a lifetime, `'a`: it holds borrowed references that must not outlive
//!   the reports they point to. It has no Java meaning beyond "this does not own the report";
//! - `match` must cover every case, so adding an error variant fails to compile until every
//!   `match` over it says what to do;
//! - `let _ = writeln!(out, ...)` appends a formatted line to a `String`. Writing to a `String`
//!   cannot fail, and `let _ =` discards the `Result` the macro returns;
//! - a function's last expression is its return value; `return` is written only for early exits.

use gatk_engine::covariates::{CovariateKind, StandardCovariateList};
use gatk_engine::java_format::format_decimals;
use gatk_engine::recal_datum::RecalDatum;
use gatk_engine::recalibration_report::RecalibrationReport;
use gatk_engine::recalibration_tables::RecalibrationTables;
use std::fmt::Write as _;

/// The csv header, in the order `printHeader` writes it.
pub const HEADER: [&str; 10] = [
    "ReadGroup",
    "CovariateValue",
    "CovariateName",
    "EventType",
    "Observations",
    "Errors",
    "EmpiricalQuality",
    "AverageReportedQuality",
    "Accuracy",
    "Recalibration",
];

/// The roles, in the order `buildReportFileMap` puts them into its `LinkedHashMap`.
///
/// This is the order the csv comes out in, whatever order the arguments were given in: `-after`
/// before `-before` on the command line still prints every `Before` row first.
pub const ROLES: [&str; 3] = ["BQSR", "Before", "After"];

/// What the tool refuses before it reads anything.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnalyzeCovariatesError {
    /// No report at all. The message's parenthesis is never closed, which is the reference's.
    NoReport,
    /// Neither a csv nor a plot was asked for.
    NoOutput,
    /// Two reports whose relevant arguments differ, each difference named.
    IncompatibleParameters(Vec<String>),
}

impl AnalyzeCovariatesError {
    /// The message the reference carries.
    pub fn message(&self) -> String {
        match self {
            // The `(` is opened and never closed. It reaches the user exactly like this.
            AnalyzeCovariatesError::NoReport => {
                "you must provide at least one recalibration report \
                 file (arguments -bqsr, -before or -after"
                    .to_string()
            }
            AnalyzeCovariatesError::NoOutput => "you need to request at least one output: the \
                 intermediate csv file (-csv FILE) or the final plot file (-plots FILE)."
                .to_string(),
            AnalyzeCovariatesError::IncompatibleParameters(differences) => format!(
                "There are differences in relevant arguments of two or more input recalibration \
                 reports. Please make sure they have been created using the same recalibration \
                 parameters. {}",
                differences.join("// ")
            ),
        }
    }

    /// The class the reference throws, which the runner prints as the exception's name:
    /// `UserException.IncompatibleRecalibrationTableParameters` for a consistency failure (thrown
    /// by `checkReportConsistency`), a plain `UserException` for the two argument checks.
    pub fn java_class(&self) -> &'static str {
        match self {
            AnalyzeCovariatesError::IncompatibleParameters(_) => {
                "org.broadinstitute.hellbender.exceptions.UserException$IncompatibleRecalibrationTableParameters"
            }
            _ => "org.broadinstitute.hellbender.exceptions.UserException",
        }
    }
}

/// One report and the role it was given under.
pub struct RoleReport<'a> {
    /// One of [`ROLES`]: the key `buildReportFileMap` files the report under, which is also the
    /// `Recalibration` column of every row it produces.
    pub role: &'a str,
    /// The report as `new RecalibrationReport(file)` reads it.
    pub report: &'a RecalibrationReport,
}

/// `doWork`: the checks, then the csv.
///
/// `reports` is what the arguments named, in any order; the roles are put back into the map's
/// order here, which is where the command line stops mattering.
///
/// `csv_requested` stands for `-csv` being given. The reference refuses only when neither `-csv`
/// nor `-plots` is (`pdfFile == null && csvFile == null`); the port draws no plots, and the runner
/// in `gatk-cli` documents asking for them as the port's limitation.
pub fn analyze_covariates(
    reports: &[RoleReport<'_>],
    csv_requested: bool,
) -> Result<String, AnalyzeCovariatesError> {
    if reports.is_empty() {
        return Err(AnalyzeCovariatesError::NoReport);
    }
    if !csv_requested {
        return Err(AnalyzeCovariatesError::NoOutput);
    }

    // `buildReportFileMap`'s `LinkedHashMap`, rebuilt: for each role in the map's order, the
    // report given under it, if any. `filter_map` keeps only the roles `find` answered `Some` for.
    let ordered: Vec<&RoleReport<'_>> = ROLES
        .iter()
        .filter_map(|role| reports.iter().find(|entry| entry.role == *role))
        .collect();

    check_consistency(&ordered)?;

    // `generateCsv` takes the covariates of the FIRST report and uses them for every one of them,
    // which is only safe because the check above has already run.
    let covariates = &ordered[0].report.covariates;
    // `generateCsv` writes the header once, then each report's rows in the map's order.
    let mut out = String::new();
    let _ = writeln!(out, "{}", HEADER.join(","));
    for entry in &ordered {
        write_csv(&mut out, &entry.report.tables, entry.role, covariates);
    }
    Ok(out)
}

/// `checkReportConsistency`, which compares every report against the FIRST one.
fn check_consistency(reports: &[&RoleReport<'_>]) -> Result<(), AnalyzeCovariatesError> {
    let first = reports[0];
    for other in &reports[1..] {
        let differences = compare_report_arguments(first, other);
        if !differences.is_empty() {
            return Err(AnalyzeCovariatesError::IncompatibleParameters(differences));
        }
    }
    Ok(())
}

/// `compareReportArguments`, whose fourteen comparisons are not fourteen checks.
///
/// ```java
/// compareSimpleReportArgument(result,"no_standard_covs", DO_NOT_USE_STANDARD_COVARIATES,
///         DO_NOT_USE_STANDARD_COVARIATES, thisRole, otherRole);
/// ```
///
/// (One line in the reference, broken in two here.)
///
/// Fourteen calls, of which the first four pass the **same constant** on both sides: those four
/// arguments are `static final` fields kept only for reading GATK3 reports, so they cannot differ
/// and the comparison cannot fire. And `indels_context_size` is not among the fourteen at all, so
/// two reports built with different indel contexts are combined without a word. (The `15` in
/// `new LinkedHashMap<>(15)` is the map's initial capacity and not a count of the checks.)
///
/// What is compared here is three of the ten live comparisons, in the reference's order, because
/// the message lists the differences in that order. The other seven (the three default
/// qualities, `default_platform`, `force_platform`, `quantizing_levels` and `binary_tag_name`)
/// are not compared, although `RecalibrationReport.initializeArgumentCollectionTable` reads every
/// one of them back from the report's argument table, for example:
///
/// ```java
/// else if (argument.equals("quantizing_levels"))
///     RAC.QUANTIZING_LEVELS = decodeInteger(value);
/// ```
///
/// so two reports that differ only in one of those seven are refused by the reference and combined
/// here.
fn compare_report_arguments(first: &RoleReport<'_>, other: &RoleReport<'_>) -> Vec<String> {
    let mut differences = Vec::new();
    // `compareSimpleReportArgument` as a closure: a function value that can see and push to
    // `differences`, which is why it is `let mut`. The values are compared as their decimal text,
    // which for these integers is the same test as the reference's `equals`. The message is the
    // reference's `capitalize(key) + ": " + value`, built in one step; `{{` and `}}` are a literal
    // brace in Rust's format string.
    let mut compare = |name: &str, left: String, right: String| {
        if left != right {
            differences.push(format!(
                "{}: differences between '{}' {{{}}} and '{}' {{{}}}.",
                capitalize(name),
                first.role,
                left,
                other.role,
                right
            ));
        }
    };

    let (left, right) = (&first.report.arguments, &other.report.arguments);
    // The four constant-against-constant comparisons are left out: they cannot produce a
    // difference, and reproducing them would mean comparing a literal to itself.
    compare(
        "mismatches_context_size",
        left.mismatches_context_size.to_string(),
        right.mismatches_context_size.to_string(),
    );
    compare(
        "maximum_cycle_value",
        left.maximum_cycle_value.to_string(),
        right.maximum_cycle_value.to_string(),
    );
    compare(
        "low_quality_tail",
        left.low_qual_tail.to_string(),
        right.low_qual_tail.to_string(),
    );
    differences
}

/// `AnalyzeCovariates.capitalize`: the first letter upper case, which is what the message does to
/// each key.
///
/// Rust's `to_uppercase` may answer more than one character where Java's
/// `Character.toUpperCase(char)` answers one; the keys are fixed ASCII names, where the two agree.
fn capitalize(text: &str) -> String {
    let mut characters = text.chars();
    match characters.next() {
        None => String::new(),
        Some(first) => first.to_uppercase().collect::<String>() + characters.as_str(),
    }
}

/// `RecalUtils.writeCsv` for one report: its rows, from the delta table the two loops build.
fn write_csv(
    out: &mut String,
    tables: &RecalibrationTables,
    mode: &str,
    covariates: &StandardCovariateList,
) {
    let mut delta: Vec<(DeltaKey, RecalDatum)> = Vec::new();

    // The quality score table, filed under an index one past the last covariate. `as i32` is a
    // numeric cast (here from the unsigned `usize`); `datum.borrow()` reads a datum the table holds
    // behind a shared cell (`Rc<RefCell<..>>`) that hands out access at run time.
    for (keys, datum) in tables.quality_score_table().all_leaves() {
        add_to_delta(
            &mut delta,
            DeltaKey {
                read_group: keys[0],
                covariate_index: covariates.size() as i32,
                covariate_key: keys[1],
                event: keys[2],
            },
            &datum.borrow(),
        );
    }

    // The optional covariates, whose key drops the reported quality. `.enumerate()` pairs each
    // covariate with its position, which indexes the matching additional table.
    for (index, kind) in covariates.additional_covariates().into_iter().enumerate() {
        let table = &tables.additional_tables()[index];
        for (keys, datum) in table.all_leaves() {
            add_to_delta(
                &mut delta,
                DeltaKey {
                    read_group: keys[0],
                    covariate_index: covariates.index_by_class(kind),
                    covariate_key: keys[2],
                    event: keys[3],
                },
                &datum.borrow(),
            );
        }
    }

    // `getAllLeaves` walks a nested array by index at every level, so the rows come out in KEY
    // order and not in the order the two loops filled them: the covariate index decides first, and
    // the quality score's out-of-range index puts every QualityScore row LAST, after Context and
    // Cycle, even though it was folded in first.
    delta.sort_by_key(|(key, _)| {
        (
            key.read_group,
            key.covariate_index,
            key.covariate_key,
            key.event,
        )
    });

    // `generateValuesFromKeys`, then `stringForCSV`, then the role. `iter_mut` hands out writable
    // references because reading the empirical quality caches it in the datum.
    for (key, datum) in delta.iter_mut() {
        let covariate_name = if key.covariate_index == covariates.size() as i32 {
            CovariateKind::QualityScore
        } else {
            covariates.kinds()[key.covariate_index as usize]
        };
        let value = format_key(covariates, covariate_name, key.covariate_key);
        let _ = writeln!(
            out,
            "{},{},{},{},{},{}",
            format_key(covariates, CovariateKind::ReadGroup, key.read_group),
            value,
            covariate_name.parsed_name(),
            pretty_event(key.event),
            string_for_csv(datum),
            mode
        );
    }
}

/// The four-part key of the delta table, which is a `NestedIntegerArray` in the reference and an
/// association list here. The traversal order of that array is the row order of the csv, and it is
/// **key order**, so the list is sorted before anything is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DeltaKey {
    read_group: i32,
    covariate_index: i32,
    covariate_key: i32,
    event: i32,
}

/// `addToDeltaTable`: a copy on the first key, and `combine` on every one after it.
///
/// The lookup is a linear scan, where the reference indexes its array directly; only the order
/// matters to the output, and the sort in `write_csv` restores it.
fn add_to_delta(delta: &mut Vec<(DeltaKey, RecalDatum)>, key: DeltaKey, datum: &RecalDatum) {
    match delta.iter_mut().find(|(seen, _)| *seen == key) {
        Some((_, existing)) => {
            // `combine` answers a `Result`, and `let _ =` discards it: an error the reference's
            // `combine` would throw is dropped here and the row is written regardless.
            let _ = existing.combine(datum);
        }
        None => delta.push((key, datum.clone())),
    }
}

/// `Covariate.formatKey` of the covariate `kind` names, as text.
///
/// The read group's formatter answers a `Result` and the context's a `Result` of an `Option`;
/// where they answer an error or nothing, the read group writes an empty field and the context
/// writes `null`.
fn format_key(covariates: &StandardCovariateList, kind: CovariateKind, key: i32) -> String {
    match kind {
        CovariateKind::ReadGroup => covariates
            .read_group
            .format_key(key)
            .map(|name| name.to_string())
            .unwrap_or_default(),
        CovariateKind::QualityScore => covariates.quality_score.format_key(key),
        CovariateKind::Context => covariates
            .context
            .format_key(key)
            .ok()
            .flatten()
            .unwrap_or_else(|| "null".to_string()),
        CovariateKind::Cycle => covariates.cycle.format_key(key),
    }
}

/// `EventType.prettyPrint()`, which is spelt out and not the letter the report carries.
///
/// ```java
/// BASE_SUBSTITUTION("M", "Base Substitution"),
/// BASE_INSERTION("I", "Base Insertion"),
/// BASE_DELETION("D", "Base Deletion");
/// ```
///
/// The ordinal is the key's last index, so it is always 0, 1 or 2.
fn pretty_event(ordinal: i32) -> &'static str {
    match ordinal {
        0 => "Base Substitution",
        1 => "Base Insertion",
        _ => "Base Deletion",
    }
}

/// `RecalDatum.stringForCSV`: `toString` then two more `%.2f`.
///
/// `toString` is `String.format("%d,%.2f,%.2f", ...)` over the observations, the mismatches and
/// the empirical quality; `stringForCSV` appends `getReportedQuality()` and
/// `getEmpiricalQuality() - getReportedQuality()`, each through `%.2f`.
///
/// Each value is rounded once, to two places; the accuracy column is the difference of the two
/// unrounded qualities, rounded after the subtraction, so it need not equal the difference of the
/// two rounded columns beside it.
///
/// `empirical_quality` is `&mut` because the reference caches it in the datum on first use, which
/// is why the datum is held mutably all the way here.
fn string_for_csv(datum: &mut RecalDatum) -> String {
    format!(
        "{},{},{},{},{}",
        datum.num_observations(),
        format_decimals(datum.num_mismatches(), 2),
        format_decimals(datum.empirical_quality(), 2),
        format_decimals(datum.reported_quality(), 2),
        format_decimals(datum.empirical_quality() - datum.reported_quality(), 2)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_roles_come_out_in_the_maps_order() {
        assert_eq!(ROLES, ["BQSR", "Before", "After"]);
    }

    #[test]
    fn the_first_refusal_never_closes_its_parenthesis() {
        let message = AnalyzeCovariatesError::NoReport.message();
        assert!(message.ends_with("-before or -after"), "{message}");
        assert_eq!(message.matches('(').count(), 1);
        assert_eq!(message.matches(')').count(), 0);
    }

    #[test]
    fn the_key_is_capitalised_and_the_differences_are_joined_by_a_double_slash() {
        let error = AnalyzeCovariatesError::IncompatibleParameters(vec![
            "Mismatches_context_size: differences between 'Before' {2} and 'After' {3}."
                .to_string(),
            "Low_quality_tail: differences between 'Before' {2} and 'After' {3}.".to_string(),
        ]);
        assert!(
            error.message().contains("{3}.// Low_quality_tail"),
            "{}",
            error.message()
        );
    }
}
