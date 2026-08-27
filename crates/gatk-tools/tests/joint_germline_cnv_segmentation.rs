//! Conformance for `JointGermlineCNVSegmentation` against GATK 4.6.2.0, compared as which calls
//! survive the entry filter, what ploidy each site reports, and which segments the defragmenter
//! joins.
//!
//! Golden from `tools/readfilter-conformance/JointGermlineCNVSegmentationDump.java`.
//!
//! # What this suite is for
//!
//!  * **four separate reasons dropping a single-sample call**, and the quality test being strictly
//!    less than;
//!  * **ploidy coming from an argument on an autosome and from the pedigree's sex on an allosome**,
//!    with an unknown sex given the male answer on both;
//!  * **a genotype with more alleles than its ploidy being refused**;
//!  * **a single no-call allele becoming a full no-call rather than being padded**;
//!  * **the defragmenter padding by a fraction of the event, so a joined run keeps growing**;
//!  * **and more than one sample skipping defragmentation entirely.**

use gatk_corpus as corpus;
use gatk_tools::joint_germline_cnv_segmentation::{
    correct_genotype_ploidy, defragment, is_multi_sample, keeps_record, padded_interval,
    sample_ploidy, Genotype, Pedigree, Segment, SegmentationError, Sex,
};
use gatk_tools::sv_stratify::SvType;

const CONTIG_LENGTH: i32 = 199980;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/joint_germline_cnv_segmentation.txt.gz"),
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
    let line = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = line.split_once(':').expect("a class and a message");
    (class.to_string(), unescape(message))
}

fn pedigree(text: &str) -> Pedigree {
    Pedigree {
        samples: section(text, "ped", "main")
            .lines()
            .filter(|line| !line.is_empty())
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                (columns[1].to_string(), Sex::parse(columns[4]))
            })
            .collect(),
    }
}

fn parse_genotype(sample: &str, format: &str, value: &str) -> Genotype {
    let keys: Vec<&str> = format.split(':').collect();
    let values: Vec<&str> = value.split(':').collect();
    let field = |key: &str| {
        keys.iter()
            .position(|name| *name == key)
            .and_then(|at| values.get(at))
            .and_then(|text| text.parse::<i32>().ok())
    };
    Genotype {
        sample: sample.to_string(),
        alleles: values[0]
            .split(['/', '|'])
            .map(|allele| allele.parse::<i32>().ok())
            .collect(),
        copy_number: field("CN"),
        quality_some: field("QS"),
        expected_copy_number: field("ECN"),
    }
}

/// One input VCF, read as the segments the tool builds from it.
fn segments(text: &str, sample: &str) -> Vec<Segment> {
    section(text, "vcf", sample)
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Segment {
                contig: columns[0].to_string(),
                start: columns[1].parse().expect("a start"),
                end: columns[7]
                    .split(';')
                    .find_map(|part| part.strip_prefix("END="))
                    .expect("an end")
                    .parse()
                    .expect("an end"),
                sv_type: SvType::parse(columns[4].trim_start_matches('<').trim_end_matches('>'))
                    .expect("a known type"),
                genotypes: vec![parse_genotype(sample, columns[8], columns[9])],
            }
        })
        .collect()
}

/// The sites one run wrote, as contig, start, end and AN.
fn measured(text: &str, label: &str) -> Vec<(String, i32, i32, i32)> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let field = |key: &str| {
                columns[7]
                    .split(';')
                    .find_map(|part| part.strip_prefix(&format!("{key}=")))
                    .unwrap_or_else(|| panic!("{label} carries {key}"))
            };
            (
                columns[0].to_string(),
                columns[1].parse().expect("a start"),
                field("END").parse().expect("an end"),
                field("AN").parse().expect("an allele number"),
            )
        })
        .collect()
}

/// Every call the male's VCF still holds after the entry filter, in input order.
fn kept(text: &str, min_quality: i32) -> Vec<Segment> {
    segments(text, "male1")
        .into_iter()
        .filter(|segment| keeps_record(segment, min_quality))
        .collect()
}

#[test]
fn the_defragmenter_matches_the_golden() {
    let text = golden();
    let kept = kept(&text, 50);
    let mut compared = 0;
    for (label, padding) in [
        ("single", 0.25),
        ("single-no-padding", 0.0),
        ("single-wide-padding", 1.0),
        // One sample always overlaps itself, so this is the same run again.
        ("single-sample-overlap", 0.25),
    ] {
        // The defragmenter works one contig at a time.
        let mut produced: Vec<(String, i32, i32)> = Vec::new();
        for contig in ["chr1", "chrX", "chrY"] {
            let on_contig: Vec<Segment> = kept
                .iter()
                .filter(|segment| segment.contig == contig)
                .cloned()
                .collect();
            for record in defragment(&on_contig, padding, 0.0, CONTIG_LENGTH) {
                produced.push((record.contig, record.start, record.end));
            }
        }
        let expected: Vec<(String, i32, i32)> = measured(&text, label)
            .into_iter()
            .map(|(contig, start, end, _)| (contig, start, end))
            .collect();
        assert_eq!(produced, expected, "{label}");
        compared += 1;
    }
    assert_eq!(compared, 4, "the single-sample runs");
}

/// Four reasons, and the fifth call that survives because the test is strictly less than.
#[test]
fn four_reasons_drop_a_single_sample_call() {
    let text = golden();
    let all = segments(&text, "male1");
    let at = |start: i32| {
        all.iter()
            .find(|segment| segment.start == start && segment.contig == "chr1")
            .unwrap_or_else(|| panic!("a segment at {start}"))
    };
    assert!(!keeps_record(at(1000), 50), "hom-ref");
    assert!(!keeps_record(at(3000), 50), "a no-call with no CN");
    assert!(!keeps_record(at(5000), 50), "a null call");
    assert!(!keeps_record(at(7000), 50), "a quality one below");
    assert!(keeps_record(at(9000), 50), "a quality exactly at it");

    // The null call is a distinct reason: its CN of 0 is present, so it is not the no-call rule.
    assert_eq!(at(5000).genotypes[0].copy_number, Some(0));
    assert!(at(5000).genotypes[0].is_null_call());
    assert!(!at(3000).genotypes[0].is_null_call());

    // And the golden's first surviving site is the one at the threshold.
    assert_eq!(measured(&text, "single")[0].1, 9000);
    // A threshold above every quality leaves nothing at all, rather than refusing.
    assert!(measured(&text, "high-quality").is_empty());
}

/// The autosome's ploidy is an argument, the allosome's is the sex, and an unknown sex is given the
/// male answer on both.
#[test]
fn ploidy_is_an_argument_on_an_autosome_and_the_sex_on_an_allosome() {
    let text = golden();
    let pedigree = pedigree(&text);
    let of = |sample: &str, contig: &str| {
        sample_ploidy(2, Some(&pedigree), sample, contig, None).expect("a ploidy")
    };
    assert_eq!(of("male1", "chr1"), 2, "the argument's default");
    assert_eq!(
        sample_ploidy(1, Some(&pedigree), "male1", "chr1", None).expect("a ploidy"),
        1
    );
    assert_eq!(of("female1", "chrX"), 2);
    assert_eq!(of("male1", "chrX"), 1);
    assert_eq!(of("unknown1", "chrX"), 1, "not the female answer");
    assert_eq!(of("female1", "chrY"), 0);
    assert_eq!(of("male1", "chrY"), 1);
    assert_eq!(of("unknown1", "chrY"), 1, "not the female answer");

    // Which the golden shows as the site's AN, without needing a genotype read at all: chrX sums
    // to 2 + 1 + 1 and chrY to 0 + 1 + 1.
    let sites = measured(&text, "default");
    let an = |contig: &str, start: i32| {
        sites
            .iter()
            .find(|(name, at, _, _)| name == contig && *at == start)
            .unwrap_or_else(|| panic!("a site at {contig}:{start}"))
            .3
    };
    assert_eq!(an("chr1", 9000), 6);
    assert_eq!(an("chrX", 1000), 4);
    assert_eq!(an("chrY", 1000), 2);
}

/// It is read before the contig is, and it does not reach the output: the site is written from a
/// ploidy derived again for every sample.
#[test]
fn an_input_ecn_wins_the_lookup_and_still_does_not_reach_the_output() {
    let text = golden();
    let pedigree = pedigree(&text);
    let carrying = Genotype {
        sample: "male1".to_string(),
        alleles: vec![Some(0), Some(1)],
        copy_number: Some(0),
        quality_some: Some(100),
        expected_copy_number: Some(7),
    };
    assert_eq!(
        sample_ploidy(2, Some(&pedigree), "male1", "chrX", Some(&carrying)).expect("a ploidy"),
        7,
        "read before the contig"
    );
    assert_eq!(
        sample_ploidy(2, Some(&pedigree), "male1", "chrX", None).expect("a ploidy"),
        1
    );
    // And the site it produced still reports the pedigree's answer.
    let sites = measured(&text, "default");
    let site = sites
        .iter()
        .find(|(contig, start, _, _)| contig == "chrX" && *start == 5000)
        .expect("the ECN site");
    assert_eq!(site.3, 4, "2 + 1 + 1, not anything involving 7");
}

/// The alleles are padded with reference alleles, a single no-call allele is a special case, and
/// too many alleles is refused.
#[test]
fn the_genotype_is_padded_to_its_ploidy() {
    let text = golden();
    let diploid = Genotype {
        sample: "male1".to_string(),
        alleles: vec![Some(0), Some(1)],
        copy_number: Some(1),
        quality_some: Some(100),
        expected_copy_number: None,
    };
    let haploid = Genotype {
        alleles: vec![Some(1)],
        ..diploid.clone()
    };
    let single_no_call = Genotype {
        alleles: vec![None],
        ..diploid.clone()
    };
    assert_eq!(
        correct_genotype_ploidy(&haploid, 2).expect("padded"),
        vec![Some(1), Some(0)],
        "padded with a reference allele"
    );
    assert_eq!(
        correct_genotype_ploidy(&single_no_call, 2).expect("padded"),
        vec![None, None],
        "a full no-call, not one no-call and one reference"
    );
    assert_eq!(
        correct_genotype_ploidy(&diploid, 2)
            .expect("unchanged")
            .len(),
        2
    );

    // Too many alleles for the ploidy, which is what an autosomal copy number of 1 does to an
    // ordinary diploid input.
    let (class, message) = refusal(&text, "haploid-autosomes");
    assert_eq!(class, "java.lang.IllegalStateException");
    let produced = correct_genotype_ploidy(&diploid, 1).expect_err("too many alleles");
    assert_eq!(
        produced,
        SegmentationError::PloidyMismatch {
            ploidy: 1,
            alleles: 2
        }
    );
    assert_eq!(produced.message(), message);
}

/// The padding is a fraction of the event, so a joined run keeps growing and a large enough
/// fraction absorbs everything.
#[test]
fn the_padding_is_a_fraction_of_the_event() {
    let text = golden();
    assert_eq!(
        padded_interval(20000, 20500, 0.25, CONTIG_LENGTH),
        (19875, 20625)
    );
    assert_eq!(
        padded_interval(40000, 60000, 0.25, CONTIG_LENGTH),
        (35000, 65000)
    );
    assert_eq!(
        padded_interval(20000, 20500, 0.0, CONTIG_LENGTH),
        (20000, 20500)
    );
    // Clipped to the contig, and never below the first base.
    assert_eq!(padded_interval(1, 100, 1.0, CONTIG_LENGTH).0, 1);

    // The same hundred-base gap is crossed by both pairs at the default, by neither at zero, and
    // at 1.0 the run does not stop until it has taken everything above it.
    let joined = |label: &str| {
        measured(&text, label)
            .into_iter()
            .filter(|(contig, _, _, _)| contig == "chr1")
            .map(|(_, start, end, _)| (start, end))
            .collect::<Vec<(i32, i32)>>()
    };
    assert_eq!(
        joined("single"),
        vec![
            (9000, 10000),
            (20000, 21100),
            (40000, 80100),
            (90000, 91000)
        ]
    );
    assert_eq!(
        joined("single-no-padding"),
        vec![
            (9000, 10000),
            (20000, 20500),
            (20600, 21100),
            (40000, 60000),
            (60100, 80100),
            (90000, 91000)
        ]
    );
    assert_eq!(
        joined("single-wide-padding"),
        vec![(9000, 10000), (20000, 91000)],
        "everything above the first record becomes one"
    );
}

/// The three input VCFs are three samples, so no padding fraction changes anything about them.
#[test]
fn more_than_one_sample_skips_defragmentation() {
    let text = golden();
    assert!(is_multi_sample(&[
        "male1".to_string(),
        "female1".to_string(),
        "unknown1".to_string()
    ]));
    assert!(!is_multi_sample(&["male1".to_string()]));

    // The multi-sample run keeps the two pairs apart that the single-sample run at the SAME
    // padding joins, which is the whole of the difference.
    let sites = |label: &str| {
        measured(&text, label)
            .into_iter()
            .filter(|(contig, _, _, _)| contig == "chr1")
            .map(|(_, start, end, _)| (start, end))
            .collect::<Vec<(i32, i32)>>()
    };
    assert!(sites("default").contains(&(20000, 20500)));
    assert!(sites("default").contains(&(20600, 21100)));
    assert!(sites("single").contains(&(20000, 21100)));
}

/// STRICT validation refuses before a record is read.
#[test]
fn a_sample_missing_from_the_pedigree_is_refused() {
    let text = golden();
    let (class, message) = refusal(&text, "missing-from-pedigree");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException"
    );
    let short = Pedigree {
        samples: vec![("male1".to_string(), Sex::Male)],
    };
    let produced = short
        .validate(&[
            "female1".to_string(),
            "male1".to_string(),
            "unknown1".to_string(),
        ])
        .expect_err("female1 is missing");
    assert_eq!(
        produced,
        SegmentationError::MissingFromPedigree {
            sample: "female1".to_string()
        }
    );
    assert_eq!(produced.message(), message);
    // The full pedigree accepts the same three.
    assert!(pedigree(&text)
        .validate(&[
            "female1".to_string(),
            "male1".to_string(),
            "unknown1".to_string()
        ])
        .is_ok());
}
