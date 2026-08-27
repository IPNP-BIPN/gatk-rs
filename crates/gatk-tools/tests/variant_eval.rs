//! Conformance for `VariantEval` against GATK 4.6.2.0, compared as the rows of the tables the port
//! reproduces.
//!
//! Golden from `tools/readfilter-conformance/VariantEvalDump.java`.
//!
//! The report's own formatting, the seven modules the port does not implement, and the derived
//! rates are in the golden and are not reproduced. What is compared is which strata each run has
//! and what `CountVariants` and `TiTvVariantEvaluator` count in each.
//!
//! # What this suite is for
//!
//!  * **novelty coming from dbSNP and not from `--comp`**;
//!  * **the standard stratifiers contributing three rows and one when turned off**;
//!  * **a stratifier multiplying the rows**;
//!  * **a multiallelic site counting once**;
//!  * **the two modules asking different questions of the same records**;
//!  * **and one message serving both name spaces.**

use gatk_corpus as corpus;
use gatk_tools::variant_eval::{
    check_module, count_variants, novelty, strata, ti_tv, EvalError, Novelty, Record, VariantType,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/variant_eval.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), message.to_string())
}

/// One table's data rows, whitespace-split, for the table whose first column is `name`.
fn table(text: &str, label: &str, name: &str) -> Vec<Vec<String>> {
    section(text, "out", label)
        .lines()
        .filter(|line| line.starts_with(name) && !line.starts_with("#"))
        .skip(1) // the column header, whose first cell is the table's name too
        .map(|line| {
            line.split_whitespace()
                .map(str::to_string)
                .collect::<Vec<String>>()
        })
        .collect()
}

/// The seven eval records the dump wrote.
fn records() -> Vec<Record> {
    let of = |position: i32, reference: &str, alternates: &[&str], alleles: &[i32]| Record {
        contig: "chr1".to_string(),
        position,
        reference: reference.to_string(),
        alternates: alternates.iter().map(|a| a.to_string()).collect(),
        alleles: alleles.to_vec(),
    };
    vec![
        of(1000, "A", &["G"], &[0, 1]),
        of(2000, "C", &["T"], &[0, 1]),
        of(3000, "A", &["C"], &[0, 1]),
        of(4000, "G", &["T"], &[1, 1]),
        of(5000, "A", &["ACGT"], &[0, 1]),
        of(6000, "ACGT", &["A"], &[0, 1]),
        of(7000, "A", &["C", "G"], &[1, 2]),
    ]
}

/// The three positions the comparison and dbSNP files carry.
const SHARED: &[i32] = &[1000, 3000, 5000];

#[test]
fn the_counts_match_the_golden() {
    let text = golden();
    let records = records();
    // Without dbSNP every record is novel, so `all` and `novel` carry the same counts.
    let rows = table(&text, "count-variants", "CountVariants");
    assert_eq!(rows.len(), 3, "all, known and novel");
    let counts = count_variants(&records);
    let all = &rows[0];
    assert_eq!(all[4], "all");
    assert_eq!(all[6], counts.called_loci.to_string(), "nCalledLoci");
    assert_eq!(all[11], counts.snps.to_string(), "nSNPs");
    assert_eq!(all[13], counts.insertions.to_string(), "nInsertions");
    assert_eq!(all[14], counts.deletions.to_string(), "nDeletions");
    assert_eq!(all[19], counts.hets.to_string(), "nHets");
    assert_eq!(all[21], counts.hom_var.to_string(), "nHomVar");
    // The known row is empty and the novel one is the whole set.
    assert_eq!(rows[1][6], "0", "nothing is known without dbSNP");
    assert_eq!(rows[2][6], counts.called_loci.to_string());
}

/// The same file as dbSNP splits the rows; as a comparison it does not.
#[test]
fn novelty_comes_from_dbsnp_and_not_from_comp() {
    let text = golden();
    let records = records();
    let known: Vec<&Record> = records
        .iter()
        .filter(|record| novelty(record, SHARED) == Novelty::Known)
        .collect();
    assert_eq!(known.len(), 3);

    let with_dbsnp = table(&text, "dbsnp", "CountVariants");
    assert_eq!(with_dbsnp[1][6], "3", "known");
    assert_eq!(with_dbsnp[2][6], "4", "novel");
    assert_eq!(
        with_dbsnp[1][11],
        count_variants(
            &known
                .iter()
                .map(|record| (*record).clone())
                .collect::<Vec<Record>>()
        )
        .snps
        .to_string()
    );

    // The same positions as a COMPARISON leave everything novel.
    let with_comp = table(&text, "count-variants", "CountVariants");
    assert_eq!(with_comp[1][6], "0");
    assert_eq!(with_comp[2][6], "7");
    // Which is what an empty dbSNP set gives.
    assert!(records
        .iter()
        .all(|record| novelty(record, &[]) == Novelty::Novel));
}

/// Three rows by default, one when they are off, and a stratifier multiplies them.
#[test]
fn the_stratifiers_multiply_the_rows() {
    let text = golden();
    assert_eq!(strata(true, &[]).len(), 3);
    assert_eq!(strata(false, &[]).len(), 1);
    assert_eq!(table(&text, "count-variants", "CountVariants").len(), 3);
    assert_eq!(
        table(&text, "no-standard-stratifiers", "CountVariants").len(),
        1
    );

    // VariantType has six values in this fixture's report, so three rows become eighteen.
    let stratified = table(&text, "stratify-by-type", "CountVariants");
    assert_eq!(stratified.len(), 18);
    assert_eq!(stratified.len(), 3 * 6);
    let types: Vec<String> = (0..6).map(|index| format!("t{index}")).collect();
    assert_eq!(strata(true, &[types]).len(), 18);

    // A named subset does the same, doubling them.
    assert_eq!(table(&text, "select-expression", "CountVariants").len(), 6);
    assert_eq!(
        strata(true, &[vec!["none".to_string(), "highqual".to_string()]]).len(),
        6
    );
}

/// It is a SNP and it is counted once, not once per alternate.
#[test]
fn a_multiallelic_site_counts_once() {
    let records = records();
    let multiallelic = records
        .iter()
        .find(|record| record.position == 7000)
        .expect("the multiallelic site");
    assert_eq!(multiallelic.alternates.len(), 2);
    assert_eq!(multiallelic.variant_type(), VariantType::Snp);
    // Seven records give five SNPs: four biallelic ones and this.
    assert_eq!(count_variants(&records).snps, 5);
    assert_eq!(count_variants(&records).called_loci, 7);
    // And it is a het, whose alleles are two different ALTERNATES.
    assert!(multiallelic.is_het());
    assert!(!multiallelic.is_hom_var());
}

/// The same records under a different question.
#[test]
fn the_two_modules_ask_different_questions() {
    let text = golden();
    let records = records();
    let counts = ti_tv(&records);
    assert_eq!(counts.transitions, 2, "A>G and C>T");
    assert_eq!(counts.transversions, 2, "A>C and G>T");
    assert_eq!(counts.ratio(), 1.0);
    // The multiallelic site is a SNP for CountVariants and NEITHER for TiTv, because the
    // substitution is only read for a biallelic one.
    let multiallelic = records
        .iter()
        .find(|record| record.position == 7000)
        .expect("the site");
    assert_eq!(multiallelic.variant_type(), VariantType::Snp);
    assert!(multiallelic.substitution().is_none());
    assert!(!multiallelic.is_transition() && !multiallelic.is_transversion());
    assert_eq!(counts.transitions + counts.transversions, 4);
    assert_eq!(count_variants(&records).snps, 5);

    // The golden's TiTv table carries the same two numbers.
    let rows = table(&text, "titv", "TiTvVariantEvaluator");
    assert_eq!(rows[0][5], counts.transitions.to_string());
    assert_eq!(rows[0][6], counts.transversions.to_string());

    // An indel's length is signed.
    assert_eq!(
        records
            .iter()
            .find(|r| r.position == 5000)
            .unwrap()
            .indel_length(),
        Some(3)
    );
    assert_eq!(
        records
            .iter()
            .find(|r| r.position == 6000)
            .unwrap()
            .indel_length(),
        Some(-3)
    );
    assert_eq!(
        records
            .iter()
            .find(|r| r.position == 1000)
            .unwrap()
            .indel_length(),
        None
    );
    // A ratio with no transversions is zero rather than infinite.
    assert_eq!(ti_tv(&records[..1]).ratio(), 0.0);
}

/// One message for both name spaces, and --list ends the process.
#[test]
fn one_message_serves_both_name_spaces() {
    let text = golden();
    for (label, name) in [
        ("unknown-module", "NoSuchEvaluator"),
        ("unknown-stratifier", "NoSuchStratifier"),
    ] {
        let (class, message) = refusal(&text, label);
        assert_eq!(
            class, "org.broadinstitute.barclay.argparser.CommandLineException",
            "{label}"
        );
        let produced = check_module(name, &["CountVariants".to_string()]).expect_err(label);
        assert_eq!(
            produced,
            EvalError::ModuleNotFound {
                name: name.to_string()
            }
        );
        assert_eq!(produced.message(), message, "{label}");
    }
    // A known name is accepted.
    assert!(check_module("CountVariants", &["CountVariants".to_string()]).is_ok());

    // --list ends the process: the golden carries the marker written before the call and not the
    // one written after it.
    assert!(golden().contains("none\tlist=about to run, which ends the process"));
    assert!(!golden().contains("after-list"));
}
