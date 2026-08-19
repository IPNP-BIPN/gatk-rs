//! Conformance for `GetPileupSummaries` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/GetPileupSummariesDump.java`.
//!
//! # What this suite is for
//!
//!  * **the bounds are strict at both ends**, so a site at exactly 0.01 or exactly 0.2 is excluded
//!    under the defaults and included once they are opened;
//!  * **only the first variant at a locus is looked at**, and only if it is a biallelic SNP;
//!  * **the counts are over `ACGT` alone**, so an `N` is not counted, a deletion is skipped, and a
//!    read below the mapping quality threshold never arrives;
//!  * **the two refusals are at opposite ends of the run**, and a port that checked both up front
//!    would refuse a run the reference completes;
//!  * **and the table carries the sample as metadata**, which is what makes it readable by
//!    `CalculateContamination`.
//!
//! The sites are built here rather than traversed: the locus walker and its pileups have their own
//! suites, and what this tool decides is which of them become rows.

use gatk_corpus as corpus;
use gatk_tools::get_pileup_summaries::{
    self, SummariesError, DEFAULT_MAX_POPULATION_AF, DEFAULT_MIN_POPULATION_AF,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/get_pileup_summaries.txt.gz"),
    )
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .to_string()
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

fn allele(bases: &str, is_ref: bool) -> Allele {
    Allele::create(bases.as_bytes(), is_ref).expect("a valid allele")
}

/// One record of the dump's population VCF.
fn variant(
    position: i64,
    reference: &str,
    alternates: &[&str],
    frequency: Option<f64>,
) -> VariantContext {
    let mut alleles = vec![allele(reference, true)];
    alleles.extend(alternates.iter().map(|bases| allele(bases, false)));
    let mut record = VariantContext::new("chr1", position, alleles);
    record.filters = Some(Vec::new());
    if let Some(frequency) = frequency {
        record
            .attributes
            .push(("AF".to_string(), Value::Double(frequency)));
    }
    record
}

/// The dump's nine records, in its order.
fn variants() -> Vec<VariantContext> {
    vec![
        variant(12, "G", &["C"], Some(0.1)),
        // Exactly the default bounds, both of which are exclusive.
        variant(13, "T", &["A"], Some(0.01)),
        variant(14, "A", &["C"], Some(0.2)),
        variant(15, "C", &["G"], Some(0.005)),
        variant(16, "G", &["T"], Some(0.5)),
        // Triallelic, and therefore not summarised whatever its frequency.
        variant(17, "T", &["A", "C"], Some(0.1)),
        // An indel, for the same reason.
        variant(18, "AC", &["A"], Some(0.1)),
        // No AF at all.
        variant(19, "G", &["C"], None),
        variant(66, "C", &["A"], Some(0.15)),
    ]
}

/// The shared fixture covers every site with one read, whose base is the reference's.
fn one_read(reference_base: u8) -> [i32; 4] {
    let mut counts = [0; 4];
    let index = match reference_base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    };
    counts[index] = 1;
    counts
}

/// The sites the traversal produced, as `(variants at the locus, base counts)`.
fn sites(window: std::ops::RangeInclusive<i64>) -> Vec<(Vec<VariantContext>, [i32; 4])> {
    variants()
        .into_iter()
        .filter(|record| window.contains(&record.start))
        .map(|record| {
            let base = record.alleles[0].base_string().as_bytes()[0];
            (vec![record], one_read(base))
        })
        .collect()
}

#[test]
fn every_table_matches_the_golden() {
    let text = golden();

    // The whole of chr1 at the defaults.
    let table = get_pileup_summaries::run(
        &sites(1..=200),
        true,
        "sample1",
        DEFAULT_MIN_POPULATION_AF,
        DEFAULT_MAX_POPULATION_AF,
    )
    .expect("the run succeeds");
    assert_eq!(escape(&table), row(&text, "table", "default"));

    // The bounds opened, which brings in the two at the edges and the one above.
    let table = get_pileup_summaries::run(&sites(1..=200), true, "sample1", 0.001, 0.9)
        .expect("the run succeeds");
    assert_eq!(escape(&table), row(&text, "table", "wide-bounds"));

    // One site.
    let table = get_pileup_summaries::run(
        &sites(66..=66),
        true,
        "sample1",
        DEFAULT_MIN_POPULATION_AF,
        DEFAULT_MAX_POPULATION_AF,
    )
    .expect("the run succeeds");
    assert_eq!(escape(&table), row(&text, "table", "one-site"));

    // A window with reads and no variants at all, which still writes the header and the metadata.
    let table = get_pileup_summaries::run(
        &[],
        true,
        "sample1",
        DEFAULT_MIN_POPULATION_AF,
        DEFAULT_MAX_POPULATION_AF,
    )
    .expect("the run succeeds");
    assert_eq!(escape(&table), row(&text, "table", "no-variants"));
}

/// The stacked site is the one that exercises the counting rule.
#[test]
fn the_counts_are_over_acgt_alone() {
    let text = golden();
    // Two `G`, one `C`, one `A`; the `N`, the deletion and the low-quality read contribute nothing.
    let counts = [1, 1, 2, 0];
    let table = get_pileup_summaries::run(
        &[(vec![variant(12, "G", &["C"], Some(0.1))], counts)],
        true,
        "stacked1",
        DEFAULT_MIN_POPULATION_AF,
        DEFAULT_MAX_POPULATION_AF,
    )
    .expect("the run succeeds");
    assert_eq!(escape(&table), row(&text, "table", "stacked"));
    assert_eq!(escape(&table), row(&text, "table", "stacked-window"));
}

#[test]
fn both_refusals_match_the_golden() {
    let text = golden();

    let error = get_pileup_summaries::run(
        &sites(1..=200),
        false,
        "sample1",
        DEFAULT_MIN_POPULATION_AF,
        DEFAULT_MAX_POPULATION_AF,
    )
    .expect_err("a header without AF");
    assert_eq!(error, SummariesError::HeaderWithoutAlleleFrequency);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        row(&text, "error", "header-without-af")
    );

    // A header that declares AF and one record that never carries it.
    let error = get_pileup_summaries::run(
        &[(vec![variant(12, "G", &["C"], None)], one_read(b'G'))],
        true,
        "sample1",
        DEFAULT_MIN_POPULATION_AF,
        DEFAULT_MAX_POPULATION_AF,
    )
    .expect_err("no record with AF");
    assert_eq!(error, SummariesError::NoRecordWithAlleleFrequency);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        row(&text, "error", "records-without-af")
    );
}

/// The eleven filters are named, so a port that dropped one would fail here rather than quietly.
#[test]
fn the_default_filter_set_is_eleven() {
    assert_eq!(get_pileup_summaries::DEFAULT_READ_FILTERS.len(), 11);
    assert_eq!(
        get_pileup_summaries::DEFAULT_READ_FILTERS[0],
        "MappingQualityReadFilter"
    );
    assert_eq!(get_pileup_summaries::DEFAULT_MINIMUM_MAPPING_QUALITY, 50);
}
