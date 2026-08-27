//! Conformance for `GenotypeGVCFs` against GATK 4.6.2.0, compared as the records of every run.
//!
//! Golden from `tools/readfilter-conformance/GenotypeGVCFsDump.java`.
//!
//! The site annotations, the qualities and the finalised FORMAT fields are in the golden and are
//! not reproduced: they come from the annotation engine. What is compared is which records are
//! written, which alternates they keep and which genotype is called.
//!
//! # What this suite is for
//!
//!  * **a reference block never being written by default**;
//!  * **`<NON_REF>` and every uncarried alternate being removed**;
//!  * **the genotype being called rather than copied**, and disagreeing with the input;
//!  * **the calling threshold not deciding emission**;
//!  * **`--include-non-variant-sites` expanding a block one record per base**;
//!  * **and `--sample-ploidy` changing nothing when the arrays are diploid.**

use gatk_corpus as corpus;
use gatk_tools::genotype_gvcfs::{
    call_genotype, genotype, genotype_index, is_emitted, trim_alleles, Arguments, Record, NON_REF,
    SNP_HETEROZYGOSITY,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/genotype_gvcfs.txt.gz"),
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

/// One run's records, as position, alternates and the called genotype.
fn measured(text: &str, label: &str) -> Vec<(i32, Vec<String>, Vec<i32>)> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let alternates: Vec<String> = if columns[4] == "." {
                Vec::new()
            } else {
                columns[4].split(',').map(str::to_string).collect()
            };
            let called: Vec<i32> = columns[9]
                .split(':')
                .next()
                .expect("a genotype")
                .split(['/', '|'])
                .map(|allele| allele.parse().unwrap_or(0))
                .collect();
            (columns[1].parse().expect("a position"), alternates, called)
        })
        .collect()
}

/// The input GVCF, read as the records the genotyper sees.
fn input(text: &str) -> Vec<Record> {
    section(text, "vcf", "input")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let start: i32 = columns[1].parse().expect("a position");
            let end = columns[7]
                .split(';')
                .find_map(|part| part.strip_prefix("END="))
                .map(|value| value.parse().expect("an end"))
                .unwrap_or(start);
            let keys: Vec<&str> = columns[8].split(':').collect();
            let values: Vec<&str> = columns[9].split(':').collect();
            let field = |key: &str| {
                keys.iter()
                    .position(|name| *name == key)
                    .and_then(|at| values.get(at))
                    .copied()
            };
            Record {
                contig: columns[0].to_string(),
                start,
                end,
                reference: columns[3].to_string(),
                alternates: columns[4].split(',').map(str::to_string).collect(),
                written_alleles: values[0]
                    .split(['/', '|'])
                    .map(|allele| allele.parse().unwrap_or(0))
                    .collect(),
                likelihoods: field("PL")
                    .map(|text| {
                        text.split(',')
                            .map(|value| value.parse().expect("a likelihood"))
                            .collect()
                    })
                    .unwrap_or_default(),
                allele_depths: field("AD")
                    .map(|text| {
                        text.split(',')
                            .map(|value| value.parse().expect("a depth"))
                            .collect()
                    })
                    .unwrap_or_default(),
            }
        })
        .collect()
}

fn produced(records: &[Record], arguments: &Arguments) -> Vec<(i32, Vec<String>, Vec<i32>)> {
    genotype(records, arguments)
        .into_iter()
        .map(|record| (record.start, record.alternates, record.called))
        .collect()
}

#[test]
fn every_run_matches_the_golden() {
    let text = golden();
    let records = input(&text);
    let mut compared = 0;
    for (label, arguments) in [
        ("default", Arguments::default()),
        (
            "all-sites",
            Arguments {
                include_non_variant_sites: true,
                ..Arguments::default()
            },
        ),
        (
            "call-threshold-2",
            Arguments {
                standard_min_confidence_for_calling: 2.0,
                ..Arguments::default()
            },
        ),
        (
            "call-threshold-50",
            Arguments {
                standard_min_confidence_for_calling: 50.0,
                ..Arguments::default()
            },
        ),
        ("keep-combined", Arguments::default()),
        ("ploidy-one", Arguments::default()),
    ] {
        assert_eq!(
            produced(&records, &arguments),
            measured(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 6, "the runs the port reproduces");
}

/// It is never written by default, and every site expands it one record per base.
#[test]
fn a_reference_block_is_never_written_by_default() {
    let text = golden();
    let records = input(&text);
    let block = &records[0];
    assert!(block.is_reference_block());
    assert_eq!(block.start, 1000);
    assert_eq!(block.end, 1004, "five bases");

    let default = measured(&text, "default");
    assert!(!default
        .iter()
        .any(|(start, ..)| (1000..=1004).contains(start)));

    let all = measured(&text, "all-sites");
    let expanded: Vec<i32> = all
        .iter()
        .map(|(start, ..)| *start)
        .filter(|start| (1000..=1004).contains(start))
        .collect();
    assert_eq!(expanded, vec![1000, 1001, 1002, 1003, 1004], "one per base");
    // And each is a bare reference record.
    for (_, alternates, called) in all.iter().filter(|(start, ..)| *start == 1000) {
        assert!(alternates.is_empty());
        assert_eq!(called, &vec![0, 0]);
    }
}

/// The prior can reverse a three-point margin, so the call disagrees with the input.
#[test]
fn the_genotype_is_called_and_not_copied() {
    let text = golden();
    let records = input(&text);
    let marginal = records
        .iter()
        .find(|record| record.start == 1400)
        .expect("the marginal site");
    // The GVCF wrote it heterozygous and its best likelihood is the het.
    assert_eq!(marginal.written_alleles, vec![0, 1]);
    assert_eq!(marginal.likelihoods[genotype_index(0, 1)], 0);
    assert!(marginal.likelihoods[genotype_index(0, 0)] > 0);
    // The call is homozygous reference anyway.
    assert_eq!(call_genotype(marginal, SNP_HETEROZYGOSITY), (0, 0));
    // Which is what the golden wrote once every site was asked for.
    let all = measured(&text, "all-sites");
    let at = all
        .iter()
        .find(|(start, ..)| *start == 1400)
        .expect("the site");
    assert_eq!(at.2, vec![0, 0]);
    assert!(at.1.is_empty());

    // A confident heterozygote is called as written, so the disagreement is the margin and not
    // the rule.
    let confident = records
        .iter()
        .find(|record| record.start == 1100)
        .expect("the confident site");
    assert_eq!(call_genotype(confident, SNP_HETEROZYGOSITY), (0, 1));
    assert_eq!(confident.written_alleles, vec![0, 1]);
    // And a confident homozygote likewise.
    let homozygote = records
        .iter()
        .find(|record| record.start == 1200)
        .expect("the homozygous site");
    assert_eq!(call_genotype(homozygote, SNP_HETEROZYGOSITY), (1, 1));
}

/// Moving it changes nothing, because emission is decided by the call.
#[test]
fn the_calling_threshold_does_not_decide_emission() {
    let text = golden();
    assert_eq!(
        measured(&text, "default"),
        measured(&text, "call-threshold-2")
    );
    assert_eq!(
        measured(&text, "default"),
        measured(&text, "call-threshold-50")
    );
    // The site the threshold would have governed is absent from all three.
    for label in ["default", "call-threshold-2", "call-threshold-50"] {
        assert!(
            !measured(&text, label)
                .iter()
                .any(|(start, ..)| *start == 1400),
            "{label}"
        );
    }
    // And it is the CALL that removes it, not a quality.
    let records = input(&text);
    let marginal = records
        .iter()
        .find(|record| record.start == 1400)
        .expect("the site");
    let called = call_genotype(marginal, SNP_HETEROZYGOSITY);
    assert!(!is_emitted(marginal, called, &Arguments::default()));
    assert!(is_emitted(
        marginal,
        called,
        &Arguments {
            include_non_variant_sites: true,
            ..Arguments::default()
        }
    ));
    // --sample-ploidy changes nothing either.
    assert_eq!(measured(&text, "default"), measured(&text, "ploidy-one"));
}

/// `<NON_REF>` goes, and so does every real alternate no called genotype carries.
#[test]
fn the_uncarried_alternates_are_removed() {
    let text = golden();
    let records = input(&text);
    let two_alts = records
        .iter()
        .find(|record| record.start == 1500)
        .expect("the two-alternate site");
    assert_eq!(two_alts.alternates, vec!["C", "G", NON_REF]);
    assert_eq!(two_alts.real_alternates(), vec!["C", "G"]);
    let called = call_genotype(two_alts, SNP_HETEROZYGOSITY);
    assert_eq!(called, (0, 1), "the sample carries the first only");
    let trimmed = trim_alleles(two_alts, called);
    assert_eq!(trimmed.alternates, vec!["C"], "G and <NON_REF> both go");
    assert_eq!(
        trimmed.called,
        vec![0, 1],
        "re-indexed against what is left"
    );

    // Which is what the golden wrote.
    let at = measured(&text, "default")
        .into_iter()
        .find(|(start, ..)| *start == 1500)
        .expect("the site");
    assert_eq!(at.1, vec!["C".to_string()]);
    // Every written record has lost <NON_REF>.
    for (_, alternates, _) in measured(&text, "default") {
        assert!(!alternates.iter().any(|allele| allele == NON_REF));
    }
}

/// The VCF's own order, which the likelihood lookup depends on.
#[test]
fn the_genotype_index_is_the_vcfs_order() {
    assert_eq!(genotype_index(0, 0), 0);
    assert_eq!(genotype_index(0, 1), 1);
    assert_eq!(genotype_index(1, 1), 2);
    assert_eq!(genotype_index(0, 2), 3);
    assert_eq!(genotype_index(1, 2), 4);
    assert_eq!(genotype_index(2, 2), 5);
    // It is symmetric, so the order the pair is given in does not matter.
    assert_eq!(genotype_index(2, 1), genotype_index(1, 2));
}
