//! Conformance for `ApplyVQSR`'s allele-specific mode against GATK 4.6.2.0, compared as the three
//! per-allele lists and the FILTER of every written record, and as the allele-specific refusal.
//!
//! Golden from `tools/readfilter-conformance/ApplyVqsrAlleleSpecificDump.java`.
//!
//! # What this suite is for
//!
//!  * **a spanning deletion counts as a SNP** and is padded anyway, without a lookup;
//!  * **an allele of the other mode is padded rather than skipped**, so the lists keep one entry per
//!    alternate allele in either mode;
//!  * **a mixed site is left unfiltered** while its `AS_FilterStatus` already names a tranche;
//!  * **there is no site-level VQSLOD or culprit** in this mode;
//!  * **and an allele with no recal record has a refusal of its own**.

use gatk_corpus as corpus;
use gatk_engine::tranches::read_tranches;
use gatk_tools::apply_vqsr::{
    allele_specific_filtering, keep, site_filter_for_a_single_mode, AlleleSpecificSite,
    AllelicRecalRecord, RecalRecord, SiteFilteringError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/apply_vqsr_allele_specific.txt.gz"),
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

fn input(text: &str, label: &str) -> String {
    unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    )
}

/// One record of either file, decoded as far as this slice looks at it.
struct Record {
    start: i32,
    end: i32,
    reference: String,
    alternates: Vec<String>,
    info: Vec<(String, Option<String>)>,
    line: String,
}

impl RecalRecord for Record {
    fn start(&self) -> i32 {
        self.start
    }
    fn end(&self) -> i32 {
        self.end
    }
    fn lod_string(&self) -> Option<String> {
        self.attribute("VQSLOD")
    }
    fn culprit(&self) -> Option<String> {
        self.attribute("culprit")
    }
    fn has_positive_label(&self) -> bool {
        self.info
            .iter()
            .any(|(name, _)| name == "POSITIVE_TRAIN_SITE")
    }
    fn has_negative_label(&self) -> bool {
        self.info
            .iter()
            .any(|(name, _)| name == "NEGATIVE_TRAIN_SITE")
    }
}

impl AllelicRecalRecord for Record {
    fn first_alternate(&self) -> String {
        self.alternates[0].clone()
    }
}

impl Record {
    fn attribute(&self, key: &str) -> Option<String> {
        self.info
            .iter()
            .find(|(name, _)| name == key)
            .and_then(|(_, value)| value.clone())
    }

    /// `isMixed()`: more than one class of alternate allele beside the reference.
    ///
    /// The spanning deletion is not a class of its own, which is what makes `A -> *,C` a SNP site.
    fn is_mixed(&self) -> bool {
        let mut snp = false;
        let mut other = false;
        for allele in &self.alternates {
            if allele == "*" {
                continue;
            }
            if allele.len() == self.reference.len() {
                snp = true;
            } else {
                other = true;
            }
        }
        snp && other
    }

    /// `checkVariationClass(vc, mode)` at the site level: `isSNP() || isMNP()` against the rest.
    fn is_of_mode(&self, snp_mode: bool) -> bool {
        let snp_site = !self.is_mixed()
            && self
                .alternates
                .iter()
                .filter(|allele| *allele != "*")
                .all(|allele| allele.len() == self.reference.len());
        if snp_mode {
            snp_site
        } else {
            !snp_site
        }
    }
}

fn records(whole: &str) -> Vec<Record> {
    whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let start: i32 = field[1].parse().expect("a position");
            let reference = field[3].to_string();
            Record {
                start,
                end: start + reference.len() as i32 - 1,
                reference,
                alternates: field[4].split(',').map(|alt| alt.to_string()).collect(),
                info: match field[7] {
                    "." => Vec::new(),
                    list => list
                        .split(';')
                        .map(|entry| match entry.split_once('=') {
                            Some((key, value)) => (key.to_string(), Some(value.to_string())),
                            None => (entry.to_string(), None),
                        })
                        .collect(),
                },
                line: line.to_string(),
            }
        })
        .collect()
}

fn written(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

/// Every written record of one run, rebuilt from the input's columns.
fn ours(text: &str, snp_mode: bool) -> Vec<String> {
    let variants = records(&input(text, "variants"));
    let recals = records(&input(text, "recal"));
    let tranches = read_tranches("tranches", &input(text, "tranches")).expect("a good file");
    let kept = keep(&tranches, 0.0);

    variants
        .iter()
        .map(|variant| {
            let site = AlleleSpecificSite {
                start: variant.start,
                end: variant.end,
                reference: &variant.reference,
                alternates: &variant.alternates,
                record: "[VC]",
            };
            let (annotations, best) = allele_specific_filtering(&site, &recals, snp_mode, &kept)
                .expect("every allele here has its record");
            let filter = site_filter_for_a_single_mode(
                variant.is_mixed(),
                variant.is_of_mode(snp_mode),
                best,
                &kept,
            );

            // The INFO the writer produces: the three lists and whichever labels were copied, in
            // byte order of the keys.
            let mut info: Vec<String> = annotations
                .attributes()
                .into_iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect();
            if annotations.positive_label {
                info.push("POSITIVE_TRAIN_SITE".to_string());
            }
            if annotations.negative_label {
                info.push("NEGATIVE_TRAIN_SITE".to_string());
            }
            info.sort();

            let mut field: Vec<String> = variant.line.split('\t').map(|f| f.to_string()).collect();
            field[6] = filter;
            field[7] = info.join(";");
            field.join("\t")
        })
        .collect()
}

#[test]
fn every_written_record_matches_the_golden_byte_for_byte() {
    let text = golden();
    assert_eq!(ours(&text, true), written(&text, "as-snp-mode"));
    assert_eq!(ours(&text, false), written(&text, "as-indel-mode"));
}

#[test]
fn a_spanning_deletion_is_padded_without_a_lookup() {
    let text = golden();
    let line = written(&text, "as-snp-mode")
        .into_iter()
        .find(|line| line.starts_with("chr1\t300\t"))
        .expect("the spanning deletion");
    assert!(
        line.contains("AS_FilterStatus=NA,VQSRTrancheSNP90.00to99.00"),
        "{line}"
    );
    assert!(line.contains("AS_VQSLOD=NaN,2.0000"), "{line}");
    // No recal record was ever written for `*`, and the run did not refuse.
    let recals = records(&input(&text, "recal"));
    assert!(!recals
        .iter()
        .any(|recal| recal.start == 300 && recal.alternates[0] == "*"));
}

#[test]
fn an_allele_of_the_other_mode_is_padded_rather_than_skipped() {
    let text = golden();
    let snp = written(&text, "as-snp-mode")
        .into_iter()
        .find(|line| line.starts_with("chr1\t200\t"))
        .expect("the mixed site");
    let indel = written(&text, "as-indel-mode")
        .into_iter()
        .find(|line| line.starts_with("chr1\t200\t"))
        .expect("the mixed site");
    // The same two alleles, the padding on opposite sides, both lists two long.
    assert!(snp.contains("AS_VQSLOD=-3.0000,NaN"), "{snp}");
    assert!(indel.contains("AS_VQSLOD=NaN,4.0000"), "{indel}");
}

#[test]
fn a_mixed_site_is_left_unfiltered_and_carries_no_site_level_score() {
    let text = golden();
    let line = written(&text, "as-snp-mode")
        .into_iter()
        .find(|line| line.starts_with("chr1\t200\t"))
        .expect("the mixed site");
    assert_eq!(line.split('\t').nth(6), Some("."));
    assert!(
        line.contains("AS_FilterStatus=VQSRTrancheSNP99.00to100.00+"),
        "{line}"
    );
    // Nothing site-level is written in this mode, in any run.
    for run in ["as-snp-mode", "as-indel-mode"] {
        for line in written(&text, run) {
            let info = line.split('\t').nth(7).expect("an INFO");
            assert!(
                !info.contains("VQSLOD=") || info.contains("AS_VQSLOD="),
                "{line}"
            );
            assert!(!info.contains(";culprit="), "{line}");
        }
    }
}

#[test]
fn an_allele_with_no_recal_record_has_a_refusal_of_its_own() {
    let text = golden();
    let row = rows(&text, "cause")
        .into_iter()
        .find(|row| row[0] == "as-missing-allele")
        .expect("a cause");
    let (class, message) = row[1].split_once(':').expect("class and message");

    let recals = records(&input(&text, "orphan-recal"));
    let variants = records(&input(&text, "orphan"));
    let variant = &variants[0];
    let error = allele_specific_filtering(
        &AlleleSpecificSite {
            start: variant.start,
            end: variant.end,
            reference: &variant.reference,
            alternates: &variant.alternates,
            record: "",
        },
        &recals,
        true,
        &keep(
            &read_tranches("tranches", &input(&text, "tranches")).expect("a good file"),
            0.0,
        ),
    )
    .expect_err("the wrong allele");
    assert_eq!(error.class(), class);
    // The message ends with the record's own toString, which is not ported.
    assert!(unescape(message).starts_with(&error.message()), "{message}");
    assert!(matches!(error, SiteFilteringError::NoRecalAllele { .. }));
}
