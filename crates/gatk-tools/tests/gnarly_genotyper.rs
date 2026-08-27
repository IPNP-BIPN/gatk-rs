//! Conformance for `GnarlyGenotyper` against GATK 4.6.2.0, compared as the engine's verdict on
//! each site and the numbers a called one carries.
//!
//! Golden from `tools/readfilter-conformance/GnarlyGenotyperDump.java`, which asks the engine
//! directly rather than through the tool: a hand-written combined GVCF never gets past the tool's
//! reader, but `finalizeGenotype` is public and takes a variant context.
//!
//! # What this suite is for
//!
//!  * **the quality floor being the confidence argument less the log of the prior**;
//!  * **which floor applies turning on the alternates' lengths**;
//!  * **a spanning deletion not counting as a SNP**;
//!  * **QUALapprox coming from the plain key, the summed list, or nowhere**;
//!  * **a dropped site being nothing at all and a kept one being LowQual**;
//!  * **a called site losing `<NON_REF>` and gaining its annotations**;
//!  * **QD and QUAL being two different numbers**;
//!  * **MQ vanishing when its raw depth is zero**;
//!  * **and the allele-specific keys being nulled rather than removed.**

use gatk_corpus as corpus;
use gatk_tools::gnarly_genotyper::{
    called_alternates, clears_the_floor, finalize_mq, has_snp_allele, indel_quality_floor,
    is_allele_specific, outcome, parse_qual_list, phred_quality, qual_approx, quality_by_depth,
    site_prior, snp_quality_floor, Outcome, AC_ADJUSTED_KEY, ALLELE_COUNT_MISMATCH_MESSAGE,
    DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING, INDEL_HETEROZYGOSITY, LOW_QUAL_FILTER_NAME, NON_REF,
    SNP_HETEROZYGOSITY, SPAN_DEL,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/gnarly_genotyper.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// One case's finalized context, as its fields.
#[derive(Debug, Clone, PartialEq)]
struct Site {
    reference: String,
    alternates: Vec<String>,
    filter: String,
    qual: String,
    info: Vec<(String, String)>,
}

fn site(text: &str, kind: &str, name: &str) -> Option<Site> {
    let line = text
        .lines()
        .find(|line| line.starts_with(&format!("{kind}\t{name}\t")))?;
    let payload = unescape(&line[format!("{kind}\t{name}\t").len()..]);
    if payload == "null" {
        return None;
    }
    let columns: Vec<&str> = payload.split('\t').collect();
    Some(Site {
        reference: columns[2].to_string(),
        alternates: columns[3].split(',').map(str::to_string).collect(),
        filter: columns[4].to_string(),
        qual: columns[5].to_string(),
        info: columns[6]
            .split(';')
            .filter(|part| !part.is_empty())
            .map(|part| match part.split_once('=') {
                Some((k, v)) => (k.to_string(), v.to_string()),
                None => (part.to_string(), String::new()),
            })
            .collect(),
    })
}

fn attribute(site: &Site, key: &str) -> Option<String> {
    site.info
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value.clone())
}

/// Whether the golden recorded a case at all, which is how a dropped site shows.
fn is_dropped(text: &str, name: &str) -> bool {
    text.lines()
        .any(|line| line == format!("out\t{name}\tnull"))
}

/// It is 60 for a SNP and 69.03 for an indel, not the 30 the argument names.
#[test]
fn the_quality_floor_is_not_the_confidence_argument() {
    assert_eq!(DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING, 30.0);
    assert_eq!(snp_quality_floor(), 60.0);
    assert!((indel_quality_floor() - 69.0308998699).abs() < 1e-9);
    assert!(snp_quality_floor() > DEFAULT_STANDARD_CONFIDENCE_FOR_CALLING);
    // A quality of 40 is above the argument and under both floors, which is what the golden's
    // two "between the floors" runs show: both come back LowQual.
    let snp = vec!["G".to_string(), NON_REF.to_string()];
    let indel = vec!["AT".to_string(), NON_REF.to_string()];
    assert!(!clears_the_floor(40.0, "A", &snp));
    assert!(!clears_the_floor(40.0, "A", &indel));
    let text = golden();
    for name in ["between-the-floors-snp", "between-the-floors-indel"] {
        let kept = site(&text, "out", name).expect("a kept site");
        assert_eq!(kept.filter, LOW_QUAL_FILTER_NAME, "{name}");
    }
    // And 900 clears both.
    assert!(clears_the_floor(900.0, "A", &snp));
    assert!(clears_the_floor(900.0, "A", &indel));
    // The boundary is inclusive.
    assert!(clears_the_floor(60.0, "A", &snp));
    assert!(!clears_the_floor(59.9, "A", &snp));
}

/// Any alternate of the reference's own length makes it a SNP.
#[test]
fn which_floor_applies_turns_on_the_alternates() {
    let text = golden();
    assert!(has_snp_allele("A", &["G".to_string()]));
    assert!(!has_snp_allele("A", &["AT".to_string()]));
    // A mixed site has one of each and is judged as a SNP.
    let mixed = vec!["G".to_string(), "AT".to_string(), NON_REF.to_string()];
    assert!(has_snp_allele("A", &mixed));
    assert_eq!(site_prior("A", &mixed), SNP_HETEROZYGOSITY);
    // `<NON_REF>` is longer than one base, so it never makes a site a SNP on its own.
    assert!(!has_snp_allele("A", &[NON_REF.to_string()]));
    assert_eq!(
        site_prior("A", &[NON_REF.to_string()]),
        INDEL_HETEROZYGOSITY
    );
    // The golden's mixed site is under the SNP floor at 40 and comes back LowQual.
    let mixed_site = site(&text, "out", "mixed-site").expect("a kept site");
    assert_eq!(mixed_site.filter, LOW_QUAL_FILTER_NAME);
    assert!(mixed_site.alternates.contains(&"AT".to_string()));
}

/// However long the reference is.
#[test]
fn a_spanning_deletion_is_not_a_snp() {
    assert_eq!(SPAN_DEL, "*");
    assert!(!has_snp_allele("A", &[SPAN_DEL.to_string()]));
    // It is excluded by NAME before the length is looked at, so a one-base reference does not
    // make it one either.
    assert!(!has_snp_allele(
        "A",
        &[SPAN_DEL.to_string(), NON_REF.to_string()]
    ));
    assert!(has_snp_allele(
        "A",
        &[SPAN_DEL.to_string(), "G".to_string()]
    ));
    let text = golden();
    let deletion = site(&text, "out", "spanning-deletion").expect("a kept site");
    assert_eq!(deletion.filter, LOW_QUAL_FILTER_NAME);
    assert!(deletion.alternates.contains(&SPAN_DEL.to_string()));
}

/// The plain key, then the summed list, then zero.
#[test]
fn qual_approx_comes_from_one_of_three_places() {
    assert_eq!(qual_approx(Some(900), None), 900.0);
    // The plain key wins even when both are there.
    assert_eq!(qual_approx(Some(900), Some("|500|400")), 900.0);
    // The list is summed, and its leading empty field skipped.
    assert_eq!(parse_qual_list("|500|400"), vec![500, 400]);
    assert_eq!(qual_approx(None, Some("|500|400")), 900.0);
    // Neither is zero, which is under every floor.
    assert_eq!(qual_approx(None, None), 0.0);
    assert!(!clears_the_floor(0.0, "A", &["G".to_string()]));
    let text = golden();
    // The allele-specific run is called, its sum clearing the floor.
    let specific = site(&text, "out", "allele-specific-qual").expect("a called site");
    assert_eq!(specific.filter, ".");
    // The run with neither key is LowQual.
    let none = site(&text, "out", "no-qual-key").expect("a kept site");
    assert_eq!(none.filter, LOW_QUAL_FILTER_NAME);
}

/// Nothing at all when it is dropped, and a filtered stub when it is kept.
#[test]
fn a_site_under_its_floor_is_dropped_or_stubbed() {
    let text = golden();
    let snp = vec!["G".to_string(), NON_REF.to_string()];
    assert_eq!(outcome(20.0, "A", &snp, false), Outcome::Dropped);
    assert_eq!(outcome(20.0, "A", &snp, true), Outcome::LowQual);
    assert_eq!(outcome(900.0, "A", &snp, false), Outcome::Called);
    // The dropped run wrote `null`.
    assert!(is_dropped(&text, "snp-under-floor"));
    assert!(site(&text, "out", "snp-under-floor").is_none());
    // The kept one is filtered, has an adjusted count of zero, and keeps its `<NON_REF>`.
    let kept = site(&text, "out", "snp-under-floor-kept").expect("a kept site");
    assert_eq!(kept.filter, LOW_QUAL_FILTER_NAME);
    assert_eq!(attribute(&kept, AC_ADJUSTED_KEY), Some("0".to_string()));
    assert!(kept.alternates.contains(&NON_REF.to_string()));
    // It has no QUAL and none of the called annotations.
    assert_eq!(kept.qual, ".");
    for key in ["AC", "AF", "AN", "QD", "FS", "SOR", "ExcessHet"] {
        assert_eq!(attribute(&kept, key), None, "{key}");
    }
    // It does keep QUALapprox, VarDP and MQ.
    assert_eq!(attribute(&kept, "QUALapprox"), Some("20".to_string()));
    assert!(attribute(&kept, "MQ").is_some());
}

/// It loses `<NON_REF>` and gains the annotations the stub does not have.
#[test]
fn a_called_site_loses_its_non_ref() {
    let text = golden();
    let called = site(&text, "out", "snp-called").expect("a called site");
    assert_eq!(called.filter, ".");
    assert!(!called.alternates.contains(&NON_REF.to_string()));
    assert_eq!(called.alternates, vec!["G"]);
    assert_eq!(
        called_alternates(&["G".to_string(), NON_REF.to_string()]),
        vec!["G"]
    );
    for key in [
        "AC",
        "AF",
        "AN",
        "QD",
        "FS",
        "SOR",
        "ExcessHet",
        "MQ",
        "VarDP",
    ] {
        assert!(attribute(&called, key).is_some(), "{key}");
    }
    // And it has no adjusted count, which only the stub carries.
    assert_eq!(attribute(&called, AC_ADJUSTED_KEY), None);
}

/// QD is one number and the site's QUAL is another.
#[test]
fn qd_and_qual_are_two_different_numbers() {
    let text = golden();
    let called = site(&text, "out", "snp-called").expect("a called site");
    // QUALapprox 900 over a variant depth of 30 is 30.
    assert_eq!(quality_by_depth(900.0, 30), 30.0);
    assert_eq!(attribute(&called, "QD"), Some("30.0".to_string()));
    // The Phred-scaled quality is 900 plus ten times the log of the prior, which is 870.
    assert!((phred_quality(900.0, SNP_HETEROZYGOSITY) - 870.0).abs() < 1e-9);
    assert_eq!(called.qual, "870.0000");
    // At the indel prior the same QUALapprox gives a different quality.
    assert!(phred_quality(900.0, INDEL_HETEROZYGOSITY) < 870.0);
}

/// A raw depth of zero is the same absence as no raw key at all.
#[test]
fn mq_vanishes_when_its_raw_depth_is_zero() {
    assert_eq!(finalize_mq("108000,30"), Some(60.0));
    assert_eq!(finalize_mq("108000,0"), None);
    assert_eq!(finalize_mq("nonsense"), None);
    let text = golden();
    // Both runs come back with no MQ at all.
    for name in ["raw-mq-without-depth", "no-raw-mq"] {
        let without = site(&text, "out", name).expect("a called site");
        assert_eq!(attribute(&without, "MQ"), None, "{name}");
    }
    // And the ordinary run has it.
    let called = site(&text, "out", "snp-called").expect("a called site");
    assert_eq!(attribute(&called, "MQ"), Some("60.00".to_string()));
}

/// The keys stay and their values become null.
#[test]
fn the_allele_specific_keys_are_nulled_rather_than_removed() {
    assert!(is_allele_specific("AS_QD"));
    assert!(!is_allele_specific("QD"));
    let text = golden();
    let stripped = site(&text, "out", "stripped-annotations").expect("a called site");
    // The keys are present and null.
    assert_eq!(attribute(&stripped, "AS_QD"), Some("null".to_string()));
    assert_eq!(attribute(&stripped, "AS_FS"), Some("null".to_string()));
    assert_eq!(attribute(&stripped, "AS_SOR"), Some("null".to_string()));
    // The site that kept them never finished: an allele-specific list whose length does not
    // match the alleles is a refusal.
    let line = text
        .lines()
        .find(|line| line.starts_with("error\tkept-annotations\t"))
        .expect("its refusal");
    let message = &line["error\tkept-annotations\t".len()..];
    assert!(
        message.starts_with("java.lang.IllegalStateException:"),
        "{message}"
    );
    assert!(message.contains(ALLELE_COUNT_MISMATCH_MESSAGE), "{message}");
    // The allele-specific quality run, whose list has the right length, is called and carries a
    // null AS_QD too.
    let specific = site(&text, "out", "allele-specific-qual").expect("a called site");
    assert_eq!(attribute(&specific, "AS_QD"), Some("null".to_string()));
}

/// Every case the golden recorded is one of the three outcomes.
#[test]
fn every_case_is_dropped_stubbed_or_called() {
    let text = golden();
    let mut dropped = 0;
    let mut stubbed = 0;
    let mut called = 0;
    let mut refused = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("out\t") {
            let name = rest.split('\t').next().expect("a name");
            match site(&text, "out", name) {
                None => dropped += 1,
                Some(found) if found.filter == LOW_QUAL_FILTER_NAME => stubbed += 1,
                Some(_) => called += 1,
            }
        } else if line.starts_with("error\t") {
            refused += 1;
        }
    }
    assert_eq!(
        dropped, 1,
        "one site under its floor with no --keep-all-sites"
    );
    assert_eq!(stubbed, 6, "six kept as LowQual");
    assert_eq!(called, 6, "six called");
    assert_eq!(refused, 1, "one allele-specific list of the wrong length");
    assert_eq!(dropped + stubbed + called + refused, 14);
}
