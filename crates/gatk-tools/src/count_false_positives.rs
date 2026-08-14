//! `CountFalsePositives`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.CountFalsePositives` (GATK 4.6.2.0).
//!
//! Two counters, a denominator and a six-column table. Each of the three hides a decision.
//!
//! # Everything that is not an indel is a SNP
//!
//! ```java
//! if (variant.isIndel()) { indelFalsePositiveCount++; } else { snpFalsePositiveCount++; }
//! ```
//!
//! `isIndel()` is `getType() == INDEL` and nothing looser, so an MNP, a symbolic allele, a mixed
//! record and a record with no alternate at all land in the `else`. The `snp` column of this table
//! is therefore "everything unfiltered that is not an indel", which is not what its name says. The
//! type itself is [`crate::remove_nearby_indels::variant_type`], measured with that tool.
//!
//! # The id is a file name
//!
//! `drivingVariantFile.getBaseName()`, which strips **one** extension, so `sample.vcf` is `sample`
//! and `sample.vcf.gz` keeps its `.vcf`. The tool's own comment says a sample name would be better
//! and that it does not know how to get one.
//!
//! # The territory is the merged intervals, in bases
//!
//! `intervalArgumentCollection.getIntervals(dictionary)` returns the sorted, merged list, so two
//! overlapping `-L` arguments contribute their union once, and `SimpleInterval.size()` is
//! `end - start + 1`. `requiresIntervals()` is true, so there is always a denominator: the missing
//! `-L` is refused by the argument parser before the tool is built.
//!
//! # The rates are per megabase, in that order
//!
//! ```java
//! return (double) snpFalsePositives / targetTerritory * 1e6;
//! ```
//!
//! The division happens first and the scaling second, which is not the same double as
//! `count * 1e6 / territory` for every input. They are written through `DataLine.set(double)`,
//! whose rounding branch is dead code ([`gatk_engine::tsv_table::java_double_to_string`]), so an
//! integral rate keeps its `.0` and the golden's `5000.0` is not `5000`.

use gatk_engine::interval::SimpleInterval;
use gatk_engine::tsv_table::{java_double_to_string, write_table};
use htsjdk_vcf::variant::VariantContext;

use crate::remove_nearby_indels::is_indel;

/// `GATKTool.getToolName()` for this tool.
pub const TOOL_NAME: &str = "GATK CountFalsePositives";

/// The table's columns, in the order `FalsePositiveTableWriter` declares them.
pub const COLUMNS: [&str; 6] = [
    "id",
    "snp",
    "indel",
    "snp_FPR",
    "indel_FPR",
    "target_territory",
];

/// What this tool refuses on its own account.
///
/// The missing `-L` is not here: `requiresIntervals()` makes the argument required, so that
/// refusal belongs to the argument parser and carries its class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CountFalsePositivesError {
    /// The output could not be opened, thrown after the whole traversal has run.
    CouldNotOpenOutput { path: String },
}

impl CountFalsePositivesError {
    /// The Java class the reference throws: a plain `UserException`, not one of its subclasses.
    pub fn class(&self) -> &'static str {
        "org.broadinstitute.hellbender.exceptions.UserException"
    }

    pub fn message(&self) -> String {
        match self {
            CountFalsePositivesError::CouldNotOpenOutput { path } => {
                format!("Encountered an IO exception while opening {path}")
            }
        }
    }
}

/// The two counters `apply` keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Counts {
    pub snp: i64,
    pub indel: i64,
}

/// `apply` over the traversal: filtered records return early, and the rest fall into two buckets.
pub fn count(variants: &[VariantContext]) -> Counts {
    let mut counts = Counts::default();
    for variant in variants {
        if variant.is_filtered() {
            continue;
        }
        if is_indel(variant) {
            counts.indel += 1;
        } else {
            counts.snp += 1;
        }
    }
    counts
}

/// The denominator: the merged intervals' bases.
pub fn target_territory(intervals: &[SimpleInterval]) -> i64 {
    intervals
        .iter()
        .map(|interval| i64::from(interval.end - interval.start + 1))
        .sum()
}

/// `FalsePositiveRecord.getSnpFalsePositiveRate`: divided first, scaled second.
pub fn false_positive_rate(count: i64, territory: i64) -> f64 {
    count as f64 / territory as f64 * 1e6
}

/// `GATKPath.getBaseName()`: the file name with **one** extension removed.
pub fn id_from_path(path: &str) -> String {
    let name = path.rsplit('/').next().unwrap_or(path);
    match name.rfind('.') {
        Some(dot) if dot > 0 => name[..dot].to_string(),
        _ => name.to_string(),
    }
}

/// The whole output file: the header line and the one record.
pub fn table(id: &str, counts: Counts, territory: i64) -> String {
    let row = vec![
        id.to_string(),
        counts.snp.to_string(),
        counts.indel.to_string(),
        java_double_to_string(false_positive_rate(counts.snp, territory)),
        java_double_to_string(false_positive_rate(counts.indel, territory)),
        territory.to_string(),
    ];
    write_table(&COLUMNS, &[row], &[])
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::allele::Allele;

    fn variant(reference: &str, alternates: &[&str], filtered: bool) -> VariantContext {
        let mut alleles = vec![Allele::create(reference.as_bytes(), true).expect("a reference")];
        for alternate in alternates {
            alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
        }
        let mut context = VariantContext::new("chr1", 100, alleles);
        context.filters = Some(if filtered {
            vec!["weak_evidence".to_string()]
        } else {
            Vec::new()
        });
        context
    }

    #[test]
    fn everything_that_is_not_an_indel_is_a_snp() {
        let variants = vec![
            variant("A", &["C"], false),        // SNP
            variant("A", &["ACC"], false),      // insertion
            variant("ACC", &["A"], false),      // deletion
            variant("AC", &["GT"], false),      // MNP
            variant("A", &["<DEL>"], false),    // symbolic
            variant("A", &["C", "ACC"], false), // mixed
            variant("A", &[], false),           // no alternate at all
            variant("A", &["C"], true),         // filtered SNP
            variant("A", &["ACC"], true),       // filtered indel
        ];
        assert_eq!(count(&variants), Counts { snp: 5, indel: 2 });
    }

    #[test]
    fn the_territory_is_the_intervals_bases() {
        let one = [SimpleInterval::new("chr1", 1, 1000).expect("valid")];
        assert_eq!(target_territory(&one), 1000);

        // What `-L chr1:1-200 -L chr1:150-400` merges to.
        let merged = [SimpleInterval::new("chr1", 1, 400).expect("valid")];
        assert_eq!(target_territory(&merged), 400);

        let disjoint = [
            SimpleInterval::new("chr1", 1, 100).expect("valid"),
            SimpleInterval::new("chr1", 900, 1000).expect("valid"),
        ];
        assert_eq!(target_territory(&disjoint), 201);
    }

    #[test]
    fn the_rate_is_per_megabase_and_keeps_its_point_zero() {
        assert_eq!(
            java_double_to_string(false_positive_rate(5, 1000)),
            "5000.0"
        );
        assert_eq!(
            java_double_to_string(false_positive_rate(1, 201)),
            "4975.124378109453"
        );
        assert_eq!(
            java_double_to_string(false_positive_rate(1, 3)),
            "333333.3333333333"
        );
        assert_eq!(java_double_to_string(false_positive_rate(0, 1000)), "0.0");
    }

    #[test]
    fn the_id_strips_one_extension() {
        assert_eq!(id_from_path("countfalsepositives-dump/mixed.vcf"), "mixed");
        assert_eq!(
            id_from_path("countfalsepositives-dump/two-extensions.vcf.vcf"),
            "two-extensions.vcf"
        );
        assert_eq!(id_from_path("plain"), "plain");
    }

    #[test]
    fn the_table_is_a_header_and_one_row() {
        assert_eq!(
            table("mixed", Counts { snp: 5, indel: 2 }, 1000),
            "id\tsnp\tindel\tsnp_FPR\tindel_FPR\ttarget_territory\nmixed\t5\t2\t5000.0\t2000.0\t1000\n"
        );
    }

    #[test]
    fn the_output_refusal_is_a_plain_user_exception() {
        let refused = CountFalsePositivesError::CouldNotOpenOutput {
            path: "countfalsepositives-dump/.".to_string(),
        };
        assert_eq!(
            refused.class(),
            "org.broadinstitute.hellbender.exceptions.UserException"
        );
        assert_eq!(
            refused.message(),
            "Encountered an IO exception while opening countfalsepositives-dump/."
        );
    }
}
