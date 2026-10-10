//! `PairHMMLikelihoodCalculationEngine.computeReadLikelihoods`, against the golden
//! `tools/readfilter-conformance/PairHmmLikelihoodEngineDump.java` wrote in the pinned container:
//! the PCR error model's qualities, the tandem repeats it reads, and for each case the prepared
//! reads, the evidence left after filtering, every normalized likelihood as raw bits, and the
//! filtered reads.

use gatk_corpus as corpus;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::pair_hmm_likelihood_engine::{
    error_model_adjusted_qual, find_tandem_repeat_units, LikelihoodEngineArguments,
    LikelihoodEngineError, PairHmmLikelihoodEngine, PcrErrorModel, MAX_REPEAT_LENGTH,
};
use gatk_engine::tsv_table::java_double_to_string;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const REF: &str = "ACGTTGCATGTCAAAAAAAGATGCACACACAGAGCTCAGTCTAGGCTTAC";

/// A read: sample, name, start, bases, phred+33 qualities, cigar, mapping quality, and the BI/BD
/// strings when the read carries them.
type Read = (
    &'static str,
    &'static str,
    i32,
    String,
    String,
    &'static str,
    u8,
    Option<(&'static str, &'static str)>,
);

fn q(length: usize, c: char) -> String {
    c.to_string().repeat(length)
}

fn reads() -> Vec<Read> {
    let snp = format!("{}G{}", &REF[..25], &REF[26..]);
    let hom_del = format!("{}{}", &REF[..12], &REF[13..]);
    let ca_ins = format!("{}CA{}", &REF[..23], &REF[23..]);
    vec![
        (
            "s1",
            "exact",
            6,
            REF[5..45].to_string(),
            q(40, '?'),
            "40M",
            60,
            None,
        ),
        (
            "s1",
            "snp",
            6,
            snp[5..45].to_string(),
            q(40, '?'),
            "40M",
            60,
            None,
        ),
        (
            "s1",
            "lowqual",
            6,
            REF[5..45].to_string(),
            "?????#####((((((5555)))))+++++????????&&".to_string(),
            "40M",
            60,
            None,
        ),
        (
            "s1",
            "mapq20",
            6,
            snp[5..45].to_string(),
            q(40, 'I'),
            "40M",
            20,
            None,
        ),
        (
            "s1",
            "softclip",
            6,
            format!("GGGGG{}TTTTT", &REF[10..40]),
            q(40, '?'),
            "5S30M5S",
            60,
            None,
        ),
        (
            "s1",
            "indelquals",
            6,
            hom_del[5..40].to_string(),
            q(35, '?'),
            "35M",
            60,
            Some((
                "####%%%%''''))))++++----////1111333",
                "DDDD::::0000&&&&!!!!######$$$$%%%&&",
            )),
        ),
        (
            "s2",
            "garbage",
            6,
            "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGG".to_string(),
            q(40, 'I'),
            "40M",
            60,
            None,
        ),
        (
            "s2",
            "homdel",
            6,
            hom_del[5..45].to_string(),
            q(40, '?'),
            "40M",
            60,
            None,
        ),
        (
            "s2",
            "ins",
            18,
            ca_ins[17..47].to_string(),
            q(30, '?'),
            "30M",
            60,
            None,
        ),
        (
            "s2",
            "short",
            11,
            REF[10..20].to_string(),
            q(10, '5'),
            "10M",
            60,
            None,
        ),
        (
            "s2",
            "mapq0",
            6,
            REF[5..45].to_string(),
            q(40, '?'),
            "40M",
            0,
            None,
        ),
        (
            "s2",
            "leftclip",
            1,
            format!("CCCCCCC{}", &REF[..25]),
            q(32, '?'),
            "7S25M",
            60,
            None,
        ),
    ]
}

fn phred(s: &str) -> Vec<u8> {
    s.bytes().map(|b| b - 33).collect()
}

fn record(r: &Read) -> BamRecord {
    let mut record = BamRecord {
        read_name: r.1.to_string(),
        reference_index: 0,
        alignment_start: r.2,
        mapping_quality: r.6,
        read_bases: r.3.as_bytes().to_vec(),
        base_qualities: phred(&r.4),
        cigar: htsjdk_bam::text_parse::parse_cigar(r.5).expect("a cigar"),
        ..Default::default()
    };
    if let Some((ins, del)) = r.7 {
        // `ReadUtils.setInsertionBaseQualities`: phred+33 text in BI and BD.
        record
            .tags
            .insert(Tag::new(b"BI"), TagValue::Str(ins.to_string()));
        record
            .tags
            .insert(Tag::new(b"BD"), TagValue::Str(del.to_string()));
    }
    record
}

struct Case {
    label: &'static str,
    arguments: LikelihoodEngineArguments,
    haplotypes: Vec<String>,
    samples: Vec<&'static str>,
}

fn pcr_name(model: PcrErrorModel) -> &'static str {
    match model {
        PcrErrorModel::None => "NONE",
        PcrErrorModel::Hostile => "HOSTILE",
        PcrErrorModel::Aggressive => "AGGRESSIVE",
        PcrErrorModel::Conservative => "CONSERVATIVE",
    }
}

fn cases() -> Vec<Case> {
    let snp = format!("{}G{}", &REF[..25], &REF[26..]);
    let hom_del = format!("{}{}", &REF[..12], &REF[13..]);
    let ca_ins = format!("{}CA{}", &REF[..23], &REF[23..]);
    let haps = vec![REF.to_string(), snp.clone(), hom_del, ca_ins];
    let both = vec!["s1", "s2"];
    let d = LikelihoodEngineArguments::default;
    let case = |label, arguments, haplotypes: &Vec<String>, samples: &Vec<&'static str>| Case {
        label,
        arguments,
        haplotypes: haplotypes.clone(),
        samples: samples.clone(),
    };
    vec![
        case("default", d(), &haps, &both),
        case(
            "pcr-none",
            LikelihoodEngineArguments {
                pcr_error_model: PcrErrorModel::None,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "pcr-hostile",
            LikelihoodEngineArguments {
                pcr_error_model: PcrErrorModel::Hostile,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "pcr-aggressive",
            LikelihoodEngineArguments {
                pcr_error_model: PcrErrorModel::Aggressive,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "dynamic-1",
            LikelihoodEngineArguments {
                dynamic_disqualification: true,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "dynamic-0.5",
            LikelihoodEngineArguments {
                dynamic_disqualification: true,
                read_disqualification_scale: 0.5,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "error-0.05",
            LikelihoodEngineArguments {
                expected_error_rate_per_base: 0.05,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "asymmetric",
            LikelihoodEngineArguments {
                symmetrically_normalize_alleles_to_reference: false,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "disable-cap",
            LikelihoodEngineArguments {
                disable_cap_read_qualities_to_map_q: true,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "keep-softclips",
            LikelihoodEngineArguments {
                modify_soft_clipped_bases: true,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "no-mismapping-cap",
            LikelihoodEngineArguments {
                log10_global_read_mismapping_rate: f64::NEG_INFINITY,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "mismapping-1",
            LikelihoodEngineArguments {
                log10_global_read_mismapping_rate: -1.0,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "gcp-5",
            LikelihoodEngineArguments {
                gap_continuation_penalty: 5,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "bq-6",
            LikelihoodEngineArguments {
                base_quality_score_threshold: 6,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "bq-25",
            LikelihoodEngineArguments {
                base_quality_score_threshold: 25,
                ..d()
            },
            &haps,
            &both,
        ),
        case("one-haplotype", d(), &vec![REF.to_string()], &both),
        case("ref-second", d(), &vec![snp, REF.to_string()], &both),
        case("empty-sample", d(), &haps, &vec!["s1", "s3"]),
        case(
            "negative-gcp",
            LikelihoodEngineArguments {
                gap_continuation_penalty: -1,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "positive-mismapping",
            LikelihoodEngineArguments {
                log10_global_read_mismapping_rate: 0.5,
                ..d()
            },
            &haps,
            &both,
        ),
        case(
            "bq-5",
            LikelihoodEngineArguments {
                base_quality_score_threshold: 5,
                ..d()
            },
            &haps,
            &both,
        ),
    ]
}

fn join(values: &[u8]) -> String {
    values
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn describe(error: &LikelihoodEngineError) -> String {
    match error {
        LikelihoodEngineError::NegativeGapContinuationPenalty => {
            "IllegalArgumentException: gap continuation penalty must be non-negative".to_string()
        }
        LikelihoodEngineError::PositiveMismappingRate => {
            "IllegalArgumentException: log10globalReadMismappingRate must be negative".to_string()
        }
        LikelihoodEngineError::BaseQualityThresholdTooLow => "IllegalArgumentException: baseQualityScoreThreshold must be greater than or equal to 6 (QualityUtils.MIN_USABLE_Q_SCORE)".to_string(),
        other => format!("{other:?}"),
    }
}

fn render(c: &Case) -> Vec<String> {
    let a = &c.arguments;
    let mut out = vec![format!(
        "case\t{}\tgcp={}\tmismapping={}\tpcr={}\tbq={}\tdynamic={}\tscale={}\terror={}\tsymmetric={}\tdisablecap={}\tkeepsoftclips={}",
        c.label,
        a.gap_continuation_penalty,
        java_double_to_string(a.log10_global_read_mismapping_rate),
        pcr_name(a.pcr_error_model),
        a.base_quality_score_threshold,
        a.dynamic_disqualification,
        java_double_to_string(a.read_disqualification_scale),
        java_double_to_string(a.expected_error_rate_per_base),
        a.symmetrically_normalize_alleles_to_reference,
        a.disable_cap_read_qualities_to_map_q,
        a.modify_soft_clipped_bases
    )];
    let engine = match PairHmmLikelihoodEngine::new(a.clone()) {
        Ok(engine) => engine,
        Err(error) => {
            out.push(format!("error\t{}\t{}", c.label, describe(&error)));
            return out;
        }
    };
    let all = reads();
    let samples: Vec<String> = c.samples.iter().map(|s| s.to_string()).collect();
    let per_sample: Vec<Vec<BamRecord>> = c
        .samples
        .iter()
        .map(|s| all.iter().filter(|r| r.0 == *s).map(record).collect())
        .collect();
    let mut hmm_length = std::collections::HashMap::new();
    for (sample, sample_reads) in c.samples.iter().zip(&per_sample) {
        let processed = engine.modify_read_qualities(sample_reads).unwrap();
        for (read, p) in sample_reads.iter().zip(&processed) {
            out.push(format!(
                "processed\t{}\t{sample}\t{}\t{}\t{}\t{}\t{}",
                c.label,
                read.read_name,
                String::from_utf8_lossy(&p.bases),
                join(&p.qualities),
                join(&p.insertion_qualities),
                join(&p.deletion_qualities)
            ));
            hmm_length.insert(
                (sample.to_string(), read.read_name.clone()),
                p.qualities.len(),
            );
        }
    }
    let haplotypes: Vec<Haplotype> = c
        .haplotypes
        .iter()
        .map(|h| Haplotype::new(h.as_bytes(), h == REF).unwrap())
        .collect();
    let result = engine
        .compute_read_likelihoods(&haplotypes, &samples, &per_sample)
        .unwrap();
    for (s, sample) in samples.iter().enumerate() {
        let evidence = result.sample_evidence(s).unwrap();
        for (r, read) in evidence.iter().enumerate() {
            out.push(format!(
                "evidence\t{}\t{sample}\t{r}\t{}\t{}",
                c.label,
                read.read_name,
                hmm_length[&(sample.clone(), read.read_name.clone())]
            ));
        }
        for allele in 0..result.number_of_alleles() {
            for r in 0..evidence.len() {
                let v = result.value(s, allele, r);
                out.push(format!(
                    "lk\t{}\t{sample}\t{allele}\t{r}\t{:016x}\t{}",
                    c.label,
                    v.to_bits(),
                    java_double_to_string(v)
                ));
            }
        }
    }
    for (s, sample) in samples.iter().enumerate() {
        for read in result.filtered_evidence(s) {
            out.push(format!(
                "filtered\t{}\t{sample}\t{}",
                c.label, read.read_name
            ));
        }
    }
    out
}

fn preamble() -> Vec<String> {
    let mut out = Vec::new();
    for model in [
        PcrErrorModel::Hostile,
        PcrErrorModel::Aggressive,
        PcrErrorModel::Conservative,
    ] {
        for i in 0..=MAX_REPEAT_LENGTH {
            out.push(format!(
                "pcr\t{}\t{i}\t{}",
                java_double_to_string(model.rate_factor()),
                error_model_adjusted_qual(i, model.rate_factor())
            ));
        }
    }
    for s in [
        REF,
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "ACACACACACACAC",
        "TTCTTCCCCAGT",
        "ACGT",
        "A",
        "AACAACAACAACGTGTGT",
        "GATATATCGCGCGCGCGCGA",
    ] {
        for offset in 0..s.len() {
            let (unit, length) = find_tandem_repeat_units(s.as_bytes(), offset);
            out.push(format!(
                "repeat\t{s}\t{offset}\t{}\t{length}",
                String::from_utf8_lossy(&unit)
            ));
        }
    }
    out
}

#[test]
fn every_case_scores_reads_as_the_reference_does() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/pair_hmm_likelihood_engine.txt.gz"),
    );
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();
    let mut actual = preamble();
    actual.extend(cases().iter().flat_map(render));
    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(got, want, "line {i}");
    }
    assert_eq!(actual.len(), expected.len(), "line count");
}
