//! Conformance for `ExcessHet`, `InbreedingCoeff`, `MQ`, `MQ0` and the tandem repeats, against the
//! oracle.
//!
//! Golden from `tools/annotation-conformance/HeterozygosityAndMqDump.java`.
//!
//! ```text
//! counts  equilibrium  rounded  4627730092099895296  4632233691727265792  4627730092099895296
//! counts  equilibrium  raw      4627730099136748619  4632233684690412495  4627730099136748607
//! eh      saturating   ExcessHet=160.0000[java.lang.String]
//! mq      all-unavailable  MQ=NaN[java.lang.String]
//! ```
//!
//! The counts are compared as **raw bits**, not as decimals: they are what the two annotations
//! divide and index by, and they come out of a chain of `Math.pow` calls this port has not proved
//! it matches. A decimal rendering would hide the last ulp, which is exactly what the suite is for.
//!
//! The `rounded` and `raw` rows on the same label are the same cohort counted twice with opposite
//! flags, so a mismatch on one and not the other localises the fault to the rounding branch.

use std::io::Read;

use gatk_annotation::heterozygosity::{
    compute_diploid_genotype_counts, excess_het, inbreeding_coeff,
};
use gatk_annotation::info_annotation::{AnnotationValue, InfoFieldAnnotation};
use gatk_annotation::mapping_quality::{MappingQualityZero, RmsMappingQuality};
use gatk_annotation::tandem_repeat::{
    find_number_of_repetitions, find_repeated_substring, TandemRepeat,
};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

const START: i64 = 105;

/// `HOM_REF`, `HET`, `HOM_VAR` and the three tied shapes, as the dump declares them.
const HOM_REF: [i32; 3] = [0, 60, 600];
const HET: [i32; 3] = [60, 0, 60];
const HOM_VAR: [i32; 3] = [600, 60, 0];
const GQ_ZERO_REF_HET: [i32; 3] = [0, 0, 60];
const GQ_ZERO_HET_VAR: [i32; 3] = [60, 0, 0];
const FLAT: [i32; 3] = [0, 0, 0];

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/heterozygosity_and_mq.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn allele(bases: &str, is_ref: bool) -> Allele {
    Allele::from_str(bases, is_ref).expect("an allele")
}

fn reference() -> Allele {
    allele("A", true)
}

fn alternate() -> Allele {
    allele("C", false)
}

fn second_alternate() -> Allele {
    allele("G", false)
}

/// The dump's `callFor`: the alleles the smallest PL implies.
fn call_for(pls: &[i32]) -> Vec<Allele> {
    let mut best = 0usize;
    for (i, pl) in pls.iter().enumerate().skip(1) {
        if *pl < pls[best] {
            best = i;
        }
    }
    match best {
        0 => vec![reference(), reference()],
        1 => vec![reference(), alternate()],
        _ => vec![alternate(), alternate()],
    }
}

fn with_genotypes(alleles: Vec<Allele>, genotypes: Vec<Genotype>) -> VariantContext {
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    vc.genotypes = genotypes;
    vc
}

fn cohort(groups: &[(&[i32], usize)]) -> VariantContext {
    let mut genotypes = Vec::new();
    let mut index = 0usize;
    for (pls, count) in groups {
        for _ in 0..*count {
            let mut genotype = Genotype::new(&format!("s{index}"), call_for(pls));
            genotype.pl = Some(pls.to_vec());
            genotypes.push(genotype);
            index += 1;
        }
    }
    with_genotypes(vec![reference(), alternate()], genotypes)
}

fn gq_only(gq: i32, count: usize) -> VariantContext {
    let genotypes = (0..count)
        .map(|i| {
            let mut genotype = Genotype::new(&format!("s{i}"), vec![reference(), reference()]);
            genotype.gq = Some(gq);
            genotype
        })
        .collect();
    with_genotypes(vec![reference(), alternate()], genotypes)
}

fn multiallelic(no_ref: bool, count: usize) -> VariantContext {
    let pls: Vec<i32> = if no_ref {
        vec![600, 300, 60, 300, 0, 60]
    } else {
        vec![60, 0, 60, 60, 60, 600]
    };
    let called = if no_ref {
        vec![alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    };
    let genotypes = (0..count)
        .map(|i| {
            let mut genotype = Genotype::new(&format!("s{i}"), called.clone());
            genotype.pl = Some(pls.clone());
            genotype
        })
        .collect();
    with_genotypes(
        vec![reference(), alternate(), second_alternate()],
        genotypes,
    )
}

fn cohort_case(label: &str) -> VariantContext {
    match label {
        "equilibrium" => cohort(&[(&HOM_REF, 25), (&HET, 50), (&HOM_VAR, 25)]),
        "all-het" => cohort(&[(&HET, 20)]),
        "all-hom-ref" => cohort(&[(&HOM_REF, 20)]),
        "all-hom-var" => cohort(&[(&HOM_VAR, 20)]),
        "excess-het-small" => cohort(&[(&HOM_REF, 2), (&HET, 8)]),
        "excess-het-large" => cohort(&[(&HOM_REF, 5), (&HET, 40), (&HOM_VAR, 5)]),
        "saturating" => cohort(&[(&HET, 200)]),
        "nine-samples" => cohort(&[(&HET, 5), (&HOM_REF, 4)]),
        "ten-samples" => cohort(&[(&HET, 5), (&HOM_REF, 5)]),
        "gq-zero-ref-het" => cohort(&[(&GQ_ZERO_REF_HET, 10)]),
        "gq-zero-het-var" => cohort(&[(&GQ_ZERO_HET_VAR, 10)]),
        "flat" => cohort(&[(&FLAT, 10)]),
        "mixed-ties" => cohort(&[(&GQ_ZERO_REF_HET, 5), (&GQ_ZERO_HET_VAR, 5)]),
        "gq-only-zero" => gq_only(0, 12),
        "gq-only-thirty" => gq_only(30, 12),
        "gq-only-ninetynine" => gq_only(99, 12),
        "multiallelic-het" => multiallelic(false, 12),
        "multiallelic-no-ref" => multiallelic(true, 12),
        "no-genotypes" => with_genotypes(vec![reference(), alternate()], Vec::new()),
        "monomorphic" => {
            let genotypes = (0..12)
                .map(|i| {
                    let mut genotype =
                        Genotype::new(&format!("s{i}"), vec![reference(), reference()]);
                    genotype.pl = Some(HOM_REF.to_vec());
                    genotype
                })
                .collect();
            with_genotypes(vec![reference()], genotypes)
        }
        other => panic!("{other} has no cohort fixture"),
    }
}

fn read(name: &str, mapping_quality: u8) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: 100,
        mapping_quality,
        cigar: htsjdk_bam::text_parse::parse_cigar("20M").expect("a cigar"),
        read_bases: vec![b'A'; 20],
        base_qualities: vec![30; 20],
        ..Default::default()
    }
}

fn matrix(mapping_qualities: &[u8]) -> AlleleLikelihoods<BamRecord> {
    let reads: Vec<BamRecord> = mapping_qualities
        .iter()
        .enumerate()
        .map(|(i, mq)| read(&format!("r{i}"), *mq))
        .collect();
    let values: Vec<Vec<f64>> = (0..2)
        .map(|a| {
            (0..reads.len())
                .map(|e| {
                    if (e % 2 == 0) == (a == 0) {
                        -1.0
                    } else {
                        -10.0
                    }
                })
                .collect()
        })
        .collect();
    AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        AlleleList::new(&[reference(), alternate()]),
        vec![reads],
        vec![values],
    )
    .expect("a matrix")
}

fn mq_case(label: &str) -> Option<AlleleLikelihoods<BamRecord>> {
    match label {
        "ordinary" => Some(matrix(&[60, 60, 60])),
        "with-zeroes" => Some(matrix(&[60, 0, 60])),
        "all-zero" => Some(matrix(&[0, 0, 0])),
        "all-unavailable" => Some(matrix(&[255, 255])),
        "mixed-unavailable" => Some(matrix(&[60, 255, 30])),
        "single-read" => Some(matrix(&[37])),
        "empty" => Some(matrix(&[])),
        "null-likelihoods" => None,
        other => panic!("{other} has no matrix fixture"),
    }
}

/// The dump's `rawRoundTrip` arguments.
fn raw_case(label: &str) -> (i64, i64) {
    match label {
        "raw-ordinary" => (60 * 60 * 3, 3),
        "raw-zero-depth" => (0, 0),
        "raw-one" => (3600, 1),
        other => panic!("{other} has no raw fixture"),
    }
}

/// The dump's `str` arguments: window start, window bases, variant start, reference, alternates.
fn str_case(label: &str) -> (i64, &'static str, i64, &'static str, &'static str) {
    match label {
        "deletion-of-one-unit" => (100, "GATCCACCACCAGTCGA", 102, "TCCA", "T"),
        "insertion-of-one-unit" => (100, "GATCCACCACCAGTCGA", 102, "T", "TCCA"),
        "not-a-repeat" => (100, "GATCCACCACCAGTCGA", 102, "TC", "T"),
        "snp-is-not-an-indel" => (100, "GATCCACCACCAGTCGA", 102, "T", "G"),
        "homopolymer" => (100, "GAAAAAAAAAAAAGTCGA", 102, "AA", "A"),
        "multiallelic-indel" => (100, "GATCCACCACCAGTCGA", 102, "TCCA", "T,TCCACCA"),
        other => panic!("{other} has no repeat fixture"),
    }
}

/// The dump's `emitMap`: `key=value[class]` joined with `;`, empty for an empty map.
fn rendered(entries: &[(String, AnnotationValue)]) -> String {
    entries
        .iter()
        .map(|(key, value)| {
            let text = value.to_java_string().expect("a renderable value");
            format!("{key}={text}[{}]", value.java_class())
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// A single `String`-valued annotation, or nothing.
fn rendered_string(key: &str, value: Option<String>) -> String {
    match value {
        Some(text) => format!("{key}={text}[java.lang.String]"),
        None => String::new(),
    }
}

#[test]
fn every_genotype_count_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("counts\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let label = fields[0];
        let rounded = match fields[1] {
            "rounded" => true,
            "raw" => false,
            other => panic!("unknown rounding {other}"),
        };
        let vc = cohort_case(label);
        let genotypes: Vec<&Genotype> = vc.genotypes.iter().collect();
        let ours = compute_diploid_genotype_counts(&vc, &genotypes, rounded);
        for (index, name) in ["refs", "hets", "homs"].iter().enumerate() {
            let want: u64 = fields[2 + index].parse::<i64>().expect("bits") as u64;
            let got = match index {
                0 => ours.refs,
                1 => ours.hets,
                _ => ours.homs,
            };
            assert_eq!(
                got.to_bits(),
                want,
                "{name} of {label} ({}) = {got}, reference {}",
                fields[1],
                f64::from_bits(want)
            );
        }
        count += 1;
    }
    assert!(count > 0, "the golden carries no count rows");
    println!("{count} genotype-count triples bit-identical");
}

#[test]
fn every_cohort_annotation_answers_as_the_reference_answers() {
    let text = golden();
    let mut excess = 0;
    let mut inbreeding = 0;
    for line in text.lines() {
        let (key, kind, rest) = if let Some(rest) = line.strip_prefix("eh\t") {
            ("ExcessHet", "eh", rest)
        } else if let Some(rest) = line.strip_prefix("ic\t") {
            ("InbreedingCoeff", "ic", rest)
        } else {
            continue;
        };
        let (label, expected) = match rest.split_once('\t') {
            Some((label, expected)) => (label, expected),
            None => (rest, ""),
        };
        let vc = cohort_case(label);
        let ours = if kind == "eh" {
            excess += 1;
            rendered_string(key, excess_het(&vc))
        } else {
            inbreeding += 1;
            rendered_string(key, inbreeding_coeff(&vc))
        };
        assert_eq!(ours, expected, "{key} on {label}");
    }
    assert!(
        excess > 0 && inbreeding > 0,
        "the golden carries no cohort rows"
    );
    println!("{excess} ExcessHet and {inbreeding} InbreedingCoeff answers identical");
}

#[test]
fn every_mapping_quality_answers_as_the_reference_answers() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let kind = ["mq", "mq0", "rawmq", "finalized"]
            .iter()
            .find(|kind| line.starts_with(&format!("{kind}\t")));
        let Some(kind) = kind else { continue };
        let rest = &line[kind.len() + 1..];
        let (label, expected) = match rest.split_once('\t') {
            Some((label, expected)) => (label, expected),
            None => (rest, ""),
        };
        let mut vc = VariantContext::new("chr1", START, vec![reference(), alternate()]);
        vc.stop = START;
        let ours = match *kind {
            "mq" => rendered(&RmsMappingQuality.annotate(None, &vc, mq_case(label).as_ref())),
            "mq0" => rendered(&MappingQualityZero.annotate(None, &vc, mq_case(label).as_ref())),
            "rawmq" => rendered(&RmsMappingQuality::annotate_raw_data(
                mq_case(label).as_ref(),
            )),
            _ => {
                let (square_sum, depth) = raw_case(label);
                let raw =
                    gatk_annotation::mapping_quality::raw_annotation_string(square_sum, depth);
                rendered_string(
                    "MQ",
                    Some(RmsMappingQuality::finalize_raw_data(&raw).expect("a well-formed tuple")),
                )
            }
        };
        assert_eq!(ours, expected, "{kind} on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no mapping-quality rows");
    println!("{count} mapping-quality answers identical");
}

#[test]
fn every_repeat_measurement_matches_the_reference() {
    let text = golden();
    let mut units = 0;
    let mut repetitions = 0;
    let mut annotations = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("repunit\t") {
            let fields: Vec<&str> = rest.split('\t').collect();
            let bases = fields[0];
            let want: usize = fields[1].parse().expect("a length");
            assert_eq!(
                find_repeated_substring(bases.as_bytes()),
                want,
                "repeat unit of {bases:?}"
            );
            units += 1;
        } else if let Some(rest) = line.strip_prefix("reps\t") {
            let fields: Vec<&str> = rest.split('\t').collect();
            let (unit, test) = (fields[0], fields[1]);
            let leading = fields[2] == "true";
            let want: i32 = fields[3].parse().expect("a count");
            assert_eq!(
                find_number_of_repetitions(unit.as_bytes(), test.as_bytes(), leading),
                Some(want),
                "{unit:?} in {test:?} leading={leading}"
            );
            repetitions += 1;
        } else if let Some(rest) = line.strip_prefix("str\t") {
            let (label, expected) = match rest.split_once('\t') {
                Some((label, expected)) => (label, expected),
                None => (rest, ""),
            };
            let (window_start, bases, variant_start, ref_allele, alts) = str_case(label);
            let mut alleles = vec![allele(ref_allele, true)];
            alleles.extend(alts.split(',').map(|alt| allele(alt, false)));
            let mut vc = VariantContext::new("chr1", variant_start, alleles);
            vc.stop = variant_start + ref_allele.len() as i64 - 1;
            let ours = rendered(&TandemRepeat::local_annotate(
                window_start,
                bases.as_bytes(),
                &vc,
            ));
            assert_eq!(ours, expected, "TandemRepeat on {label}");
            annotations += 1;
        }
    }
    assert!(
        units > 0 && repetitions > 0 && annotations > 0,
        "the golden carries no repeat rows"
    );
    println!(
        "{units} repeat units, {repetitions} repetition counts and {annotations} STR annotations identical"
    );
}
