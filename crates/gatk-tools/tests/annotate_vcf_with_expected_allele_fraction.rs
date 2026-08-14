//! Conformance for `AnnotateVcfWithExpectedAlleleFraction` against GATK 4.6.2.0, compared as the
//! `AF_EXP` of every written record, the header lines the tool builds, and both refusals.
//!
//! Golden from `tools/readfilter-conformance/AnnotateVcfWithExpectedAlleleFractionDump.java`.
//!
//! # What this suite is for
//!
//!  * **the weights are in column order and the fractions in sorted order**, and the two are paired
//!    by position, so a zebra-only het is annotated with alpha's fraction;
//!  * **the output carries no `##source=` line**, unlike its sibling's;
//!  * **the weight is 1.0, 0.5 or nothing**, a half-call weighing nothing;
//!  * **and the two refusals are Java's own**, one of them with no message at all.

use gatk_corpus as corpus;
use gatk_tools::annotate_vcf_with_expected_allele_fraction::{
    annotation, fractions_in_sample_order, header_lines, AF_EXP, AF_EXP_HEADER_LINE,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/annotate_vcf_with_expected_allele_fraction.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn whole(text: &str, kind: &str, label: &str) -> String {
    unescape(
        rows(text, kind)
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no {kind} {label}"))[1],
    )
}

/// The header's sample columns, in the order the file declares them.
fn header_samples(text: &str) -> Vec<String> {
    whole(text, "input", "variants")
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("a column header")
        .split('\t')
        .skip(9)
        .map(|sample| sample.to_string())
        .collect()
}

/// One mixing-fraction table, in the order its rows are written.
fn table(text: &str, label: &str) -> Vec<(String, f64)> {
    whole(text, "fractions", label)
        .lines()
        .skip(1)
        .filter(|line| !line.is_empty())
        .map(|line| {
            let (sample, value) = line.split_once('\t').expect("two columns");
            (sample.to_string(), value.parse().expect("a fraction"))
        })
        .collect()
}

/// The input's records, with their genotypes in the file's column order.
fn variants(text: &str) -> Vec<VariantContext> {
    let samples = header_samples(text);
    whole(text, "input", "variants")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let mut alleles = vec![Allele::create(field[3].as_bytes(), true).expect("a reference")];
            for alternate in field[4].split(',') {
                alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
            }
            let mut variant =
                VariantContext::new(field[0], field[1].parse().expect("a position"), alleles);
            let called = variant.alleles.clone();
            variant.genotypes = samples
                .iter()
                .enumerate()
                .map(|(index, sample)| {
                    Genotype::new(
                        sample,
                        field[9 + index]
                            .split(['/', '|'])
                            .map(|allele| match allele.parse::<usize>() {
                                Ok(at) => called[at].clone(),
                                Err(_) => Allele::no_call(),
                            })
                            .collect(),
                    )
                })
                .collect();
            variant
        })
        .collect()
}

/// The `AF_EXP` of every record of one run, taken off the golden's own output lines.
fn written(text: &str, label: &str) -> Vec<(i64, String)> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == label)
        .map(|row| unescape(row[1]))
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let value = field[7]
                .split(';')
                .filter_map(|entry| entry.split_once('='))
                .find(|(key, _)| *key == AF_EXP)
                .map(|(_, value)| value.to_string())
                .expect("every record carries the annotation");
            (field[1].parse().expect("a position"), value)
        })
        .collect()
}

#[test]
fn every_annotation_matches_the_golden() {
    let text = golden();
    let samples = header_samples(&text);
    for (run, fractions) in [("annotated", "fractions"), ("normalized", "normalized")] {
        let ordered = fractions_in_sample_order(&table(&text, fractions), &samples)
            .expect("every sample is named");
        let mine: Vec<(i64, String)> = variants(&text)
            .iter()
            .map(|variant| (variant.start, annotation(variant, &ordered)))
            .collect();
        assert_eq!(mine, written(&text, run), "the annotations of {run}");
    }
}

#[test]
fn a_zebra_only_het_is_paired_with_alphas_fraction() {
    let text = golden();
    // The columns are zebra, alpha, mike and the fractions zebra=0.3, alpha=0.2, mike=0.1.
    assert_eq!(header_samples(&text), vec!["zebra", "alpha", "mike"]);
    let first = written(&text, "annotated")
        .into_iter()
        .find(|(at, _)| *at == 20)
        .expect("the first record");
    assert_eq!(first.1, "0.100", "0.5 * 0.2, not 0.5 * 0.3");
}

#[test]
fn the_output_carries_no_source_line() {
    let text = golden();
    let written_lines: Vec<String> = rows(&text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == "annotated")
        .map(|row| unescape(row[1]))
        .collect();
    assert!(written_lines.contains(&AF_EXP_HEADER_LINE.to_string()));
    assert!(
        !written_lines
            .iter()
            .any(|line| line.starts_with("##source=")),
        "the default tool lines are added after the header is built"
    );
    assert!(
        rows(&text, "commandline").is_empty(),
        "and no command line either"
    );

    let input: Vec<String> = whole(&text, "input", "variants")
        .lines()
        .filter(|line| line.starts_with("##"))
        .map(|line| line.to_string())
        .collect();
    assert!(!header_lines(&input)
        .iter()
        .any(|line| line.starts_with("##source=")));
}

#[test]
fn a_no_call_and_a_half_call_weigh_nothing() {
    let text = golden();
    let zeroed = written(&text, "annotated")
        .into_iter()
        .find(|(at, _)| *at == 80)
        .expect("the record of no-calls");
    assert_eq!(zeroed.1, "0.00");
}

#[test]
fn both_refusals_carry_the_references_class_and_words() {
    let text = golden();
    let refusal = |label: &str| -> (String, String) {
        let row = rows(&text, "error")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no refusal {label}"));
        let (class, message) = row[1].split_once(':').expect("class and message");
        (class.to_string(), message.to_string())
    };

    let samples = header_samples(&text);
    let missing = fractions_in_sample_order(&table(&text, "missing-sample"), &samples)
        .expect_err("a sample is missing");
    let (class, message) = refusal("missing-sample");
    assert_eq!(missing.class(), class);
    assert_eq!(missing.message(), message);

    let duplicated = fractions_in_sample_order(&table(&text, "duplicate-sample"), &samples)
        .expect_err("a sample is named twice");
    let (class, message) = refusal("duplicate-sample");
    assert_eq!(duplicated.class(), class);
    assert_eq!(duplicated.message(), message);
}
