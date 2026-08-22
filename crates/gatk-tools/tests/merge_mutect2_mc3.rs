//! Conformance for `MergeMutect2CallsWithMC3` against GATK 4.6.2.0, compared as the whole set of
//! records every run writes.
//!
//! Golden from `tools/readfilter-conformance/MergeMutect2WithMC3Dump.java`.
//!
//! # What this suite is for
//!
//!  * **each of the five states writing a different record**, and one writing nothing;
//!  * **a false positive being rebuilt from scratch**, so the M2 record's ID, QUAL, FILTER column
//!    and INFO fields are dropped;
//!  * **a false negative never gaining `CENTERS`**;
//!  * **the genotype's ploidy being the number of alleles at the site**;
//!  * **the depths coming from M2 when M2 is there and from `NREF`/`NALT` when it is not**;
//!  * **and an M2 genotype without `AD` leaving the output genotype without one.**

use gatk_corpus as corpus;
use gatk_engine::concordance_walker::{concordance, ConcordanceRecord};
use gatk_tools::merge_mutect2_mc3::{concordant, merge, Merged, Variant};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/merge_mutect2_mc3.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn parse(line: &str) -> Variant {
    let fields: Vec<&str> = line.split('\t').collect();
    let mut alleles = vec![fields[3].to_string()];
    if fields[4] != "." {
        alleles.extend(fields[4].split(',').map(str::to_string));
    }
    let info = if fields[7] == "." {
        Vec::new()
    } else {
        fields[7]
            .split(';')
            .map(|entry| match entry.split_once('=') {
                Some((key, value)) => (key.to_string(), value.to_string()),
                None => (entry.to_string(), String::new()),
            })
            .collect()
    };
    let format: Vec<&str> = fields[8].split(':').collect();
    let allele_depths = format
        .iter()
        .position(|key| *key == "AD")
        .and_then(|index| fields[9].split(':').nth(index))
        .map(|depths| {
            depths
                .split(',')
                .map(|depth| depth.parse().expect("a depth"))
                .collect()
        });
    Variant {
        contig: fields[0].to_string(),
        start: fields[1].parse().expect("a start"),
        id: fields[2].to_string(),
        quality: fields[5].to_string(),
        filters: if fields[6] == "." {
            Vec::new()
        } else {
            fields[6].split(';').map(str::to_string).collect()
        },
        alleles,
        info,
        allele_depths,
    }
}

/// One written record as the VCF writer renders it: INFO sorted by key, the genotype as slashes.
fn rendered(record: &Merged) -> String {
    let mut info: Vec<String> = record
        .info
        .iter()
        .map(|(key, value)| format!("{key}={value}"))
        .collect();
    info.sort();
    let genotype: Vec<String> = record
        .genotype
        .iter()
        .map(|index| index.to_string())
        .collect();
    let (format, sample) = match &record.allele_depths {
        Some(depths) => (
            "GT:AD".to_string(),
            format!(
                "{}:{}",
                genotype.join("/"),
                depths
                    .iter()
                    .map(|depth| depth.to_string())
                    .collect::<Vec<String>>()
                    .join(",")
            ),
        ),
        None => ("GT".to_string(), genotype.join("/")),
    };
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        record.contig,
        record.start,
        record.id,
        record.alleles[0],
        record.alleles[1..].join(","),
        record.quality,
        if record.filters.is_empty() {
            ".".to_string()
        } else {
            record.filters.join(";")
        },
        if info.is_empty() {
            ".".to_string()
        } else {
            info.join(";")
        },
        format,
        sample
    )
}

fn data_lines(vcf: &str) -> Vec<String> {
    vcf.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.to_string())
        .collect()
}

/// The dictionary the walk orders by, which every fixture here shares.
const DICTIONARY: [&str; 1] = ["chr1"];

fn run(text: &str, label: &str) -> Vec<String> {
    let truth: Vec<Variant> = data_lines(&value(text, "input", &format!("{label}-truth")))
        .iter()
        .map(|line| parse(line))
        .collect();
    let eval: Vec<Variant> = data_lines(&value(text, "input", &format!("{label}-eval")))
        .iter()
        .map(|line| parse(line))
        .collect();
    // The base class drops filtered TRUTH records and keeps every eval record.
    let kept_truth: Vec<Variant> = truth
        .iter()
        .filter(|record| !record.is_filtered())
        .cloned()
        .collect();
    let dictionary: Vec<String> = DICTIONARY.iter().map(|name| name.to_string()).collect();
    let steps = concordance(&kept_truth, &eval, &dictionary, concordant);
    merge(&steps, &kept_truth, &eval)
        .iter()
        .map(rendered)
        .collect()
}

fn expected(text: &str, label: &str) -> Vec<String> {
    data_lines(&value(text, "merged", label))
}

#[test]
fn every_merged_vcf_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for label in ["every-state", "eval-without-ad", "multiallelic"] {
        assert_eq!(run(&text, label), expected(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 3, "the golden's outputs");
}

/// The M2 record's ID, QUAL, FILTER column and INFO fields are all dropped, and only CENTERS=M2
/// survives.
#[test]
fn a_false_positive_is_rebuilt_from_scratch() {
    let text = golden();
    // The eval record at 120 is `rs1 ... 99 PASS M2ONLY=dropped`.
    let source = data_lines(&value(&text, "input", "every-state-eval"))
        .into_iter()
        .find(|line| line.starts_with("chr1\t120\t"))
        .expect("the M2 record");
    assert!(source.contains("\trs1\t"));
    assert!(source.contains("M2ONLY=dropped"));

    let written = run(&text, "every-state")
        .into_iter()
        .find(|line| line.starts_with("chr1\t120\t"))
        .expect("the written record");
    assert_eq!(
        written,
        "chr1\t120\t.\tA\tC\t.\t.\tCENTERS=M2\tGT:AD\t0/1:13,14"
    );
}

/// A filtered true negative writes nothing, so the eval-only filtered call at 130 is absent.
#[test]
fn a_filtered_true_negative_writes_nothing() {
    let text = golden();
    let written = run(&text, "every-state");
    assert!(written.iter().all(|line| !line.starts_with("chr1\t130\t")));
}

/// A false negative is emitted unchanged except for its genotype, so it never gains CENTERS.
#[test]
fn a_false_negative_never_gains_centers() {
    let text = golden();
    let written = run(&text, "every-state")
        .into_iter()
        .find(|line| line.starts_with("chr1\t110\t"))
        .expect("the written record");
    assert!(!written.contains("CENTERS"));
    // And its depths come from NREF and NALT rather than from any genotype.
    assert!(written.ends_with("GT:AD\t0/1:70,30"));

    // A truth-only record with neither field defaults both to zero.
    let zeroes = run(&text, "every-state")
        .into_iter()
        .find(|line| line.starts_with("chr1\t170\t"))
        .expect("the written record");
    assert!(zeroes.ends_with("GT:AD\t0/1:0,0"));
}

/// The genotype carries one index per allele at the site, so a multiallelic false positive is
/// triploid and a biallelic record can carry three depths.
#[test]
fn the_genotype_takes_every_allele_at_the_site() {
    let text = golden();
    let written = run(&text, "multiallelic");
    let triploid = written
        .iter()
        .find(|line| line.contains("\tA\tG,T\t"))
        .expect("the multiallelic false positive");
    assert!(triploid.contains("GT:AD\t0/1/2:10,20,30"));

    // And the concordant multiallelic record is biallelic in the output while its AD is not.
    let biallelic = written
        .iter()
        .find(|line| line.starts_with("chr1\t100\t"))
        .expect("the true positive");
    assert!(biallelic.contains("\tA\tC\t"));
    assert!(biallelic.ends_with("GT:AD\t0/1:10,20,30"));
}

/// `getAD()` answering null leaves the genotype without an AD rather than throwing.
#[test]
fn an_eval_without_ad_leaves_the_genotype_without_one() {
    let text = golden();
    let written = run(&text, "eval-without-ad");
    assert_eq!(written.len(), 1);
    assert!(written[0].ends_with("\tGT\t0/1"));
    assert_eq!(written, expected(&text, "eval-without-ad"));
}
