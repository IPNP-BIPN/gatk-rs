//! `VariantFiltration`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration` (GATK 4.6.2.0).
//!
//! A variant transform whose whole job is the FILTER column: JEXL expressions at the site and at
//! the genotype, a clustered-SNP test that looks at neighbours, and a mask read from a second file.
//!
//! # An empty set is not nothing
//!
//! ```java
//! //Note that making this empty set effectively converts the VC to PASS, whereas an unfiltered VC has null filters
//! if ( filters.isEmpty() ) {
//!     if (!invalidatePreviousFilters) {
//!         builder.passFilters();
//!     } else {
//!         builder.unfiltered();
//!     }
//! ```
//!
//! The comment is the reference's own. One condition, two columns: `PASS` normally, and `.` under
//! `--invalidate-previous-filters`, which also wipes whatever filter the record arrived with.
//!
//! # The order in the file is the writer's
//!
//! The tool builds a `LinkedHashSet` seeded from the record's own filters, so insertion order is
//! deliberate; the vcf writer then **sorts**, so `OldFilter` plus `LowQD` is written
//! `LowQD;OldFilter`, and two expressions given as LowQD then LowDP are written `LowDP;LowQD`. This
//! port returns the set and leaves the ordering to whoever writes, which is where it belongs.
//!
//! # What is not ported here
//!
//! The context a JEXL expression sees is [`crate::variant_filtration::Context`], built from the
//! INFO attributes and the genotype's own fields, which is what the golden exercises. GATK's full
//! `VariantJEXLContext` also exposes derived predicates (`vc.isSNP`, `isHet`, and the rest); those
//! are a brick of their own and nothing here pretends to them.

use gatk_engine::jexl::{create_expression, Expression, JexlError, Value};
use std::collections::HashMap;

/// `FILTER_DELIMITER`, which is how an existing FT is split back into names.
pub const FILTER_DELIMITER: char = ';';
/// `CLUSTERED_SNP_FILTER_NAME`.
pub const CLUSTERED_SNP_FILTER_NAME: &str = "SnpCluster";

/// One named expression, which is `JexlVCMatchExp`.
pub struct MatchExp {
    pub name: String,
    pub expression: Expression,
}

impl MatchExp {
    pub fn new(name: &str, text: &str) -> Result<MatchExp, JexlError> {
        Ok(MatchExp {
            name: name.to_string(),
            expression: create_expression(text)?,
        })
    }
}

/// The attributes an expression reads: the INFO fields, and a genotype's own fields when one is
/// given.
pub type Context = HashMap<String, String>;

/// As much of a record as the filtering reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub contig: String,
    pub start: i32,
    pub stop: i32,
    /// True where every alternate is one base against a one-base reference.
    pub is_snp: bool,
    /// The FILTER column as it arrived: `None` for `.`, empty for `PASS`.
    pub filters: Option<Vec<String>>,
    pub info: Context,
    pub genotypes: Vec<GenotypeFields>,
}

/// One sample's fields, as far as a genotype expression reads them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotypeFields {
    pub fields: Context,
    /// The FT the genotype arrived with, already split.
    pub filters: Vec<String>,
}

/// The arguments this port reads.
pub struct Arguments {
    pub cluster_size: i32,
    pub cluster_window: i32,
    pub mask_name: String,
    pub filter_records_not_in_mask: bool,
    pub invert_filter_expression: bool,
    pub invert_genotype_filter_expression: bool,
    pub missing_values_evaluate_as_failing: bool,
    pub invalidate_previous_filters: bool,
    pub set_filtered_genotypes_to_no_call: bool,
    pub mask_extension: i32,
}

impl Default for Arguments {
    fn default() -> Arguments {
        Arguments {
            cluster_size: 3,
            cluster_window: 0,
            mask_name: "Mask".to_string(),
            filter_records_not_in_mask: false,
            invert_filter_expression: false,
            invert_genotype_filter_expression: false,
            missing_values_evaluate_as_failing: false,
            invalidate_previous_filters: false,
            set_filtered_genotypes_to_no_call: false,
            mask_extension: 0,
        }
    }
}

/// What one record becomes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Filtered {
    /// `None` writes `.`, an empty list writes `PASS`, and the rest are the names, unsorted: the
    /// writer sorts them.
    pub filters: Option<Vec<String>>,
    /// One list per sample, empty where that sample is unfiltered. The FT column exists in the file
    /// only when at least one of these is non-empty.
    pub genotype_filters: Vec<Vec<String>>,
    /// Whether each call is replaced by a no-call, which `--set-filtered-genotype-to-no-call` does.
    pub no_called: Vec<bool>,
}

/// The whole traversal, which the cluster test makes non-local: a record's filters depend on its
/// neighbours.
pub fn filter_records(
    records: &[Record],
    site_expressions: &[MatchExp],
    genotype_expressions: &[MatchExp],
    mask: &[(String, i32, i32)],
    arguments: &Arguments,
) -> Vec<Filtered> {
    records
        .iter()
        .enumerate()
        .map(|(index, _record)| {
            filter_one(
                records,
                index,
                site_expressions,
                genotype_expressions,
                mask,
                arguments,
            )
        })
        .collect()
}

fn filter_one(
    records: &[Record],
    index: usize,
    site_expressions: &[MatchExp],
    genotype_expressions: &[MatchExp],
    mask: &[(String, i32, i32)],
    arguments: &Arguments,
) -> Filtered {
    let record = &records[index];

    // `invalidatePreviousFilters` empties the record BEFORE anything is added, so a filter it
    // arrived with is gone whatever happens next.
    let mut filters: Vec<String> = if arguments.invalidate_previous_filters {
        Vec::new()
    } else {
        record.filters.clone().unwrap_or_default()
    };

    // The mask is only added when it is not already there, which is what `isMaskFilterPresent`
    // guards.
    if !filters.contains(&arguments.mask_name) {
        let overlaps = mask.iter().any(|(contig, start, stop)| {
            *contig == record.contig
                && record.start <= stop + arguments.mask_extension
                && start - arguments.mask_extension <= record.stop
        });
        if overlaps != arguments.filter_records_not_in_mask {
            filters.push(arguments.mask_name.clone());
        }
    }

    if are_clustered_snps(records, index, arguments) {
        filters.push(CLUSTERED_SNP_FILTER_NAME.to_string());
    }

    for expression in site_expressions {
        if matches(
            &record.info,
            expression,
            arguments.invert_filter_expression,
            arguments.missing_values_evaluate_as_failing,
        ) {
            push_once(&mut filters, &expression.name);
        }
    }

    let mut genotype_filters = Vec::with_capacity(record.genotypes.len());
    let mut no_called = Vec::with_capacity(record.genotypes.len());
    for genotype in &record.genotypes {
        let mut names = genotype.filters.clone();
        for expression in genotype_expressions {
            let mut context = record.info.clone();
            context.extend(genotype.fields.clone());
            if matches(
                &context,
                expression,
                arguments.invert_genotype_filter_expression,
                arguments.missing_values_evaluate_as_failing,
            ) {
                push_once(&mut names, &expression.name);
            }
        }
        no_called.push(arguments.set_filtered_genotypes_to_no_call && !names.is_empty());
        genotype_filters.push(names);
    }

    Filtered {
        filters: if filters.is_empty() {
            if arguments.invalidate_previous_filters {
                None
            } else {
                Some(Vec::new())
            }
        } else {
            Some(filters)
        },
        genotype_filters,
        no_called,
    }
}

/// `areClusteredSNPs`: a window below 1 disables it, and only SNPs count, as candidate and as
/// neighbour.
fn are_clustered_snps(records: &[Record], index: usize, arguments: &Arguments) -> bool {
    if arguments.cluster_window < 1 {
        return false;
    }
    let current = &records[index];
    if !current.is_snp {
        return false;
    }
    let nearby: Vec<i32> = records
        .iter()
        .filter(|other| {
            other.is_snp
                && other.contig == current.contig
                && other.start >= current.start - arguments.cluster_window
                && other.start <= current.start + arguments.cluster_window
        })
        .map(|other| other.start)
        .collect();
    if (nearby.len() as i32) < arguments.cluster_size {
        return false;
    }

    // Any run of `clusterSize` neighbours containing this one and spanning no more than the window.
    let mut sorted = nearby;
    sorted.sort_unstable();
    let size = arguments.cluster_size as usize;
    sorted.windows(size).any(|run| {
        run.first().copied().unwrap_or(0) <= current.start
            && current.start <= run.last().copied().unwrap_or(0)
            && run[size - 1] - run[0] <= arguments.cluster_window
    })
}

/// `matchesFilter`: the match, then the inversion. A missing field is an error in the reference's
/// non-lenient engine, and `--missing-values-evaluate-as-failing` is what turns it into a match.
fn matches(
    context: &Context,
    expression: &MatchExp,
    invert: bool,
    missing_as_failing: bool,
) -> bool {
    let matched = match expression.expression.evaluate(context) {
        Ok(Value::Bool(value)) => value,
        Ok(_) => false,
        Err(_) => missing_as_failing,
    };
    matched != invert
}

/// A `LinkedHashSet` keeps the first insertion and ignores the rest.
fn push_once(names: &mut Vec<String>, name: &str) {
    if !names.iter().any(|seen| seen == name) {
        names.push(name.to_string());
    }
}

/// How the writer renders a filter set: sorted, or `PASS`, or `.`.
pub fn rendered_filters(filtered: &Filtered) -> String {
    match &filtered.filters {
        None => ".".to_string(),
        Some(names) if names.is_empty() => "PASS".to_string(),
        Some(names) => {
            let mut sorted = names.clone();
            sorted.sort();
            sorted.join(";")
        }
    }
}

/// The FT of one sample, which is `PASS` when the record has an FT column at all and this sample is
/// unfiltered, and absent when no sample of the record was filtered.
pub fn rendered_genotype_filter(filtered: &Filtered, sample: usize) -> Option<String> {
    if filtered
        .genotype_filters
        .iter()
        .all(|names| names.is_empty())
    {
        return None;
    }
    let names = &filtered.genotype_filters[sample];
    if names.is_empty() {
        return Some("PASS".to_string());
    }
    let mut sorted = names.clone();
    sorted.sort();
    Some(sorted.join(";"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(start: i32, is_snp: bool, info: &[(&str, &str)], filters: Option<&[&str]>) -> Record {
        Record {
            contig: "chr1".to_string(),
            start,
            stop: start,
            is_snp,
            filters: filters.map(|names| names.iter().map(|n| n.to_string()).collect()),
            info: info
                .iter()
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect(),
            genotypes: Vec::new(),
        }
    }

    #[test]
    fn an_empty_set_is_pass_or_a_dot() {
        let records = vec![record(100, true, &[("QD", "20.0")], None)];
        let expressions = vec![MatchExp::new("Never", "QD > 1000.0").expect("parses")];

        let kept = filter_records(&records, &expressions, &[], &[], &Arguments::default());
        assert_eq!(rendered_filters(&kept[0]), "PASS");

        let invalidated = filter_records(
            &records,
            &expressions,
            &[],
            &[],
            &Arguments {
                invalidate_previous_filters: true,
                ..Arguments::default()
            },
        );
        assert_eq!(rendered_filters(&invalidated[0]), ".");
    }

    #[test]
    fn the_writer_sorts_what_the_tool_ordered() {
        let records = vec![record(100, true, &[("QD", "1.0")], Some(&["OldFilter"]))];
        let expressions = vec![MatchExp::new("LowQD", "QD < 2.0").expect("parses")];
        let kept = filter_records(&records, &expressions, &[], &[], &Arguments::default());
        // Insertion order is OldFilter then LowQD; the file says otherwise.
        assert_eq!(
            kept[0].filters,
            Some(vec!["OldFilter".to_string(), "LowQD".to_string()])
        );
        assert_eq!(rendered_filters(&kept[0]), "LowQD;OldFilter");
    }

    #[test]
    fn the_cluster_test_ignores_indels_and_a_narrow_window() {
        let records = vec![
            record(1000, true, &[], None),
            record(1005, false, &[], None),
            record(1010, true, &[], None),
            record(1015, true, &[], None),
            record(9000, true, &[], None),
        ];
        let wide = Arguments {
            cluster_size: 3,
            cluster_window: 20,
            ..Arguments::default()
        };
        let kept = filter_records(&records, &[], &[], &[], &wide);
        let rendered: Vec<String> = kept.iter().map(rendered_filters).collect();
        assert_eq!(
            rendered,
            vec!["SnpCluster", "PASS", "SnpCluster", "SnpCluster", "PASS"]
        );

        let narrow = Arguments {
            cluster_window: 5,
            ..wide
        };
        let kept = filter_records(&records, &[], &[], &[], &narrow);
        assert!(kept.iter().all(|one| rendered_filters(one) == "PASS"));
    }
}
