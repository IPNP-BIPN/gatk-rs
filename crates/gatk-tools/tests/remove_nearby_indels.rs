//! Conformance for `RemoveNearbyIndels` against GATK 4.6.2.0, compared as the records that reach
//! the output of every run.
//!
//! Golden from `tools/readfilter-conformance/RemoveNearbyIndelsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the distance is end to start and the test is strict**, which the same file at two
//!    neighbouring spacings shows;
//!  * **a discarded indel is still what the next one is measured against**;
//!  * **the non-indels between a discarded pair are kept**;
//!  * **a trailing indel survives on a reference comparison**;
//!  * **and a mixed site is not an indel**, so it never pairs with anything.

use gatk_corpus as corpus;
use gatk_tools::remove_nearby_indels::{is_indel, remove_nearby_indels};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/remove_nearby_indels.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

/// The reverse of the dump's `escape`, scanning once so a real backslash is never read as a tab.
fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// One `CHROM POS ID REF ALT ...` line as a variant, which is all this tool looks at.
fn parse_record(line: &str) -> VariantContext {
    let fields: Vec<&str> = line.split('\t').collect();
    let mut alleles = vec![Allele::create(fields[3].as_bytes(), true).expect("a reference")];
    for alternate in fields[4].split(',') {
        alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
    }
    VariantContext::new(fields[0], fields[1].parse().expect("a position"), alleles)
}

/// The records of one input file, in file order.
fn input(text: &str, label: &str) -> Vec<String> {
    let whole = rows(text, "input")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no input {label}"))[1]
        .to_string();
    unescape(&whole)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.to_string())
        .collect()
}

/// The records the reference kept in one run, in output order.
fn kept(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .filter(|line| !line.starts_with('#'))
        .collect()
}

/// Every run the golden holds, as its input label and its spacing.
fn runs(text: &str) -> Vec<(String, String, i32)> {
    let mut seen: Vec<(String, String, i32)> = Vec::new();
    for row in rows(text, "vcfline") {
        let run = row[0].to_string();
        if seen.iter().any(|(name, _, _)| name == &run) {
            continue;
        }
        let (label, spacing) = run.rsplit_once("-at-").expect("a run label");
        seen.push((
            run.clone(),
            label.to_string(),
            spacing.parse().expect("a spacing"),
        ));
    }
    seen
}

#[test]
fn every_run_keeps_the_records_the_reference_kept() {
    let text = golden();
    let all = runs(&text);
    assert!(all.len() >= 10, "the golden holds every run: {}", all.len());

    for (run, label, spacing) in &all {
        let records = input(&text, label);
        let variants: Vec<VariantContext> = records.iter().map(|line| parse_record(line)).collect();
        let ours: Vec<String> = remove_nearby_indels(&variants, *spacing)
            .into_iter()
            .map(|index| records[index].clone())
            .collect();
        assert_eq!(ours, kept(&text, run), "run/{run}");
    }
}

/// The same file at two neighbouring spacings, which is the strict comparison on its own.
#[test]
fn one_base_of_spacing_decides_a_pair() {
    let text = golden();
    let at_five = kept(&text, "boundary-at-5");
    let at_four = kept(&text, "boundary-at-4");
    assert_eq!(at_five.len(), 2, "one pair is dropped at 5");
    assert_eq!(at_four.len(), 4, "no pair is dropped at 4");
}

/// The mixed site reaches the output and the indel four bases away with it.
#[test]
fn a_mixed_site_pairs_with_nothing() {
    let text = golden();
    let records = input(&text, "multi-allelic");
    let variants: Vec<VariantContext> = records.iter().map(|line| parse_record(line)).collect();
    assert!(!is_indel(&variants[0]), "a snp and an insertion are MIXED");
    assert!(is_indel(&variants[1]));
    assert_eq!(
        remove_nearby_indels(&variants, 20).len(),
        records.len(),
        "nothing is dropped"
    );
}

/// A whole file lost, because the last indel is measured against one already thrown away.
#[test]
fn a_wide_spacing_reaches_past_the_indel_it_discarded() {
    let text = golden();
    let records = input(&text, "spaced");
    let variants: Vec<VariantContext> = records.iter().map(|line| parse_record(line)).collect();
    let ours: Vec<String> = remove_nearby_indels(&variants, 1000)
        .into_iter()
        .map(|index| records[index].clone())
        .collect();
    // Every snp survives and no indel does, though the last indel is 790 bases from the one before.
    assert_eq!(ours, kept(&text, "spaced-at-1000"));
    assert!(ours.iter().all(|line| !is_indel(&parse_record(line))));
}
