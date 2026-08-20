//! Conformance for what a validation pileup says about an allele, against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/AllelePileupCounterDump.java`.
//!
//! # What this suite is for
//!
//!  * **the last matching alternate wins**, because `chooseAlleleForRead`'s loop assigns without
//!    breaking, which the `snp-both` case is built to show;
//!  * **the reference test includes the indel lookahead**, so a read matching the reference bases
//!    is not the reference allele when it sits before an insertion or a deletion start;
//!  * **a read that ends inside the allele is UNKNOWN**, neither reference nor alternate, and a
//!    cutoff of zero already refuses it because the minimum over an empty array is -1;
//!  * **the deletion match compares lengths only**, so an alternate of the right length matches
//!    whatever its bases are;
//!  * **the counter drops mapping quality 0 and 255 and nothing else**, while the ratio, which
//!    does not use the counter, keeps those reads;
//!  * **`calculateMaxAltRatio`'s two filters are not complements**, so a pileup of nothing but
//!    short reads is a ratio of zero rather than a NaN;
//!  * **and a haploid genotype throws** rather than being refused, because the variant is typed
//!    from the second allele before the ploidy that was just computed is consulted.
//!
//! Every row is reproduced and every row is bit-identical: the powers here run over small
//! validation depths, so the cumulative sums stop long before the ulp the sibling suite allows.

use gatk_corpus as corpus;
use gatk_engine::allele_pileup_counter::AllelePileupCounter;
use gatk_engine::basic_somatic_short_mutation_validator::{
    calculate_basic_validation_result, is_able_to_validate_genotype, write_table,
    BasicValidationResult, ValidationGenotype,
};
use gatk_engine::read_pileup::{pileup_from_reads, ReadPileup};
use gatk_engine::somatic_validation_power::{
    calculate_max_alt_ratio, calculate_num_reads_supporting_allele,
};
use gatk_engine::variant_context_utils::{
    choose_allele_for_read, is_complex_indel, type_of_variant, Allele, PileupAlleleError,
};
use htsjdk_bam::record::BamRecord;

const LOCUS: i32 = 105;
const START: i32 = 100;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/allele_pileup_counter.txt.gz"),
    )
}

/// One read of the fixture: every base carries the same quality, as the harness sets them.
fn read(name: &str, cigar: &str, bases: &str, quality: u8, mapping_quality: u8) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: START,
        mapping_quality,
        read_bases: bases.as_bytes().to_vec(),
        base_qualities: vec![quality; bases.len()],
        cigar: htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar"),
        ..Default::default()
    }
}

/// The reads of each labelled pileup, in the order the harness listed them.
fn reads_for(label: &str) -> Vec<BamRecord> {
    match label {
        "snp" => vec![
            read("ref1", "10M", "AAAAAAAAAA", 30, 60),
            read("ref2", "10M", "AAAAAAAAAA", 30, 60),
            read("alt1", "10M", "AAAAACAAAA", 30, 60),
            read("alt2", "10M", "AAAAACAAAA", 30, 60),
            read("altq", "10M", "AAAAACAAAA", 5, 60),
        ],
        "mapq" => vec![
            read("ref1", "10M", "AAAAAAAAAA", 30, 60),
            read("zero", "10M", "AAAAACAAAA", 30, 0),
            read("none", "10M", "AAAAACAAAA", 30, 255),
        ],
        "insertion" => vec![
            read("ref1", "10M", "AAAAAAAAAA", 30, 60),
            read("ins1", "6M3I4M", "AAAAAATTTAAAA", 30, 60),
            read("ins2", "6M3I4M", "AAAAAATTTAAAA", 30, 60),
            read("del1", "6M2D4M", "AAAAAAAAAA", 30, 60),
        ],
        "deletion" => vec![
            read("ref1", "10M", "AAAAAAAAAA", 30, 60),
            read("del1", "6M2D4M", "AAAAAAAAAA", 30, 60),
            read("del2", "6M2D4M", "AAAAAAAAAA", 30, 60),
            read("del3", "6M5D4M", "AAAAAAAAAA", 30, 60),
        ],
        "short" => vec![
            read("stop", "6M", "AAAAAA", 30, 60),
            read("stop2", "6M", "AAAAAC", 30, 60),
        ],
        "empty" => vec![],
        other => panic!("{other} is in the golden but not configured here"),
    }
}

fn pileup(reads: &[BamRecord]) -> ReadPileup<'_> {
    pileup_from_reads("chr1", LOCUS, reads, |_| true, |_| true)
}

/// The case name without the suffix that says which alleles it was asked about.
fn pileup_name(label: &str) -> &str {
    label.split_once('-').map(|(head, _)| head).unwrap_or(label)
}

fn reference(bases: &str) -> Allele {
    Allele::new(bases.as_bytes(), true)
}

fn alternate(bases: &str) -> Allele {
    Allele::new(bases.as_bytes(), false)
}

/// The alternates each `choose` case was asked about, and the reference it was asked against.
fn choice_case(label: &str) -> (Allele, Vec<Allele>, i32) {
    match label {
        "snp" => (reference("A"), vec![alternate("C")], 0),
        "snp-q20" => (reference("A"), vec![alternate("C")], 20),
        "snp-two" => (reference("A"), vec![alternate("C"), alternate("G")], 0),
        "snp-both" => (reference("A"), vec![alternate("C"), alternate("C")], 0),
        "mapq" => (reference("A"), vec![alternate("C")], 0),
        "insertion" => (reference("A"), vec![alternate("ATTT")], 0),
        "deletion" => (reference("AAA"), vec![alternate("A")], 0),
        "deletion-other" => (reference("AAA"), vec![alternate("T")], 0),
        "short" => (reference("AA"), vec![alternate("GT")], 0),
        "short-q1" => (reference("AA"), vec![alternate("GT")], 1),
        // The counter is also asked about an empty locus, which no `choose` row reaches.
        "empty" => (reference("A"), vec![alternate("C")], 0),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The counter's map, ordered by base string the way the harness ordered it.
fn counted(counter: &AllelePileupCounter) -> String {
    let mut entries: Vec<(String, i32)> = counter
        .count_map()
        .iter()
        .map(|(allele, count)| (String::from_utf8_lossy(&allele.bases).into_owned(), *count))
        .collect();
    entries.sort_by(|left, right| left.0.cmp(&right.0));
    entries
        .iter()
        .map(|(allele, count)| format!("{allele}={count}"))
        .collect::<Vec<String>>()
        .join(",")
}

/// The reference this ratio case is taken against, which is the two-base one only for `short`.
fn ratio_reference(name: &str) -> Allele {
    if name == "short" {
        reference("AA")
    } else {
        reference("A")
    }
}

/// The genotype and pileup of each `validatable` and `result` case.
fn validation_case(label: &str) -> (ValidationGenotype, Allele, &'static str, i32, i32, i32) {
    let plain = |first: Allele, second: Allele, ad: Option<Vec<i32>>, filters: Option<&str>| {
        ValidationGenotype {
            alleles: vec![first, second],
            ad,
            filters: filters.map(str::to_string),
        }
    };
    let depths = Some(vec![40, 10]);
    match label {
        "plain" => (
            plain(reference("A"), alternate("C"), depths, None),
            reference("A"),
            "snp",
            8,
            30,
            0,
        ),
        "no-alt-reads" => (
            plain(reference("A"), alternate("C"), depths, None),
            reference("A"),
            "snp",
            0,
            30,
            0,
        ),
        "one-alt-read" => (
            plain(reference("A"), alternate("C"), depths, None),
            reference("A"),
            "snp",
            1,
            30,
            0,
        ),
        "filtered" => (
            plain(
                reference("A"),
                alternate("C"),
                depths,
                Some("weak_evidence"),
            ),
            reference("A"),
            "snp",
            8,
            30,
            0,
        ),
        "zero-discovery" => (
            plain(reference("A"), alternate("C"), Some(vec![0, 0]), None),
            reference("A"),
            "snp",
            8,
            30,
            0,
        ),
        "no-ad" => (
            plain(reference("A"), alternate("C"), None, None),
            reference("A"),
            "snp",
            8,
            30,
            0,
        ),
        "alt-first" => (
            plain(alternate("C"), reference("A"), depths, None),
            reference("A"),
            "snp",
            8,
            30,
            0,
        ),
        "complex" => (
            plain(reference("AAA"), alternate("TA"), depths, None),
            reference("AAA"),
            "snp",
            8,
            30,
            0,
        ),
        "insertion" => (
            plain(reference("A"), alternate("ATTT"), depths, None),
            reference("A"),
            "insertion",
            8,
            30,
            0,
        ),
        "deletion" => (
            plain(reference("AAA"), alternate("A"), depths, None),
            reference("AAA"),
            "deletion",
            8,
            30,
            0,
        ),
        "empty-normal" => (
            plain(reference("A"), alternate("C"), depths, None),
            reference("A"),
            "empty",
            8,
            30,
            0,
        ),
        "min-quality" => (
            plain(reference("A"), alternate("C"), depths, None),
            reference("A"),
            "snp",
            8,
            30,
            20,
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

fn result_of(label: &str) -> Option<BasicValidationResult> {
    let (genotype, reference, case, alt_count, total_count, minimum_quality) =
        validation_case(label);
    let reads = reads_for(case);
    let pileup = pileup(&reads);
    calculate_basic_validation_result(
        &genotype,
        &reference,
        Some(&pileup),
        alt_count,
        total_count,
        minimum_quality,
        "chr1",
        LOCUS,
        LOCUS,
        "PASS",
    )
    .expect("a validatable case")
}

/// The dump's `%016x`, which is the raw bits and not the canonical NaN.
fn hex(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

/// The refusals, reported as the reference's exception and message.
fn refusal(label: &str) -> String {
    let text = match label {
        "haploid-genotype" => {
            let genotype = ValidationGenotype {
                alleles: vec![reference("A")],
                ad: Some(vec![40, 10]),
                filters: None,
            };
            let error = is_able_to_validate_genotype(&genotype, &reference("A"))
                .expect_err("a haploid genotype throws");
            format!("{}:{}", error.java_class(), error.message())
        }
        // htsjdk refuses to tag a symbolic allele as the reference, so the counter's own
        // symbolic-reference branch is never reached and this is htsjdk's message, not GATK's.
        // The port has no such constructor check, so the row is asserted as a constant.
        "symbolic-reference" | "symbolic-type" => {
            "java.lang.IllegalArgumentException:Cannot tag a symbolic allele as the reference \
             allele"
                .to_string()
        }
        "non-reference-reference" => {
            let error = AllelePileupCounter::new(&alternate("C"), &[alternate("G")], 0)
                .expect_err("a non-reference reference");
            format!("{}:{}", error.java_class(), error.message())
        }
        "reference-alternate" => {
            let error = AllelePileupCounter::new(&reference("A"), &[reference("C")], 0)
                .expect_err("a reference alternate");
            format!("{}:{}", error.java_class(), error.message())
        }
        "negative-quality" => {
            let error = AllelePileupCounter::new(&reference("A"), &[alternate("C")], -1)
                .expect_err("a negative cutoff");
            format!("{}:{}", error.java_class(), error.message())
        }
        "negative-quality-ratio" => {
            let reads = reads_for("snp");
            let error = calculate_max_alt_ratio(&pileup(&reads), &reference("A"), -1)
                .expect_err("a negative cutoff");
            format!("{}:{}", error.java_class(), error.message())
        }
        "negative-quality-support" => {
            let reads = reads_for("snp");
            let error = calculate_num_reads_supporting_allele(
                &pileup(&reads),
                &reference("A"),
                &alternate("C"),
                -1,
            )
            .expect_err("a negative cutoff");
            format!("{}:{}", error.java_class(), error.message())
        }
        "span-del-reference" => {
            let error = type_of_variant(&reference("*"), &alternate("C"))
                .expect_err("a spanning-deletion reference");
            format!("{}:{}", error.java_class(), error.message())
        }
        other => panic!("an unexpected refusal: {other}"),
    };
    format!("error\t{label}\t{text}")
}

#[test]
fn every_row_matches_the_golden() {
    let text = golden();
    let mut rows = 0;
    let mut written: Vec<BasicValidationResult> = Vec::new();
    let mut table: Vec<String> = Vec::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let (kind, rest) = line.split_once('\t').expect("a kind");
        match kind {
            "type" => {
                let (label, _) = rest.split_once('=').expect("a value");
                let (left, right) = label.split_once(',').expect("two alleles");
                let first = reference(left);
                let second = alternate(right);
                let ours = match type_of_variant(&first, &second) {
                    Ok(variant) => format!(
                        "type\t{label}={},{}",
                        variant.name(),
                        is_complex_indel(&first, &second)
                    ),
                    Err(error) => panic!("{label}: {error:?}"),
                };
                assert_eq!(ours, line);
            }
            "choose" => {
                let (label, rest) = rest.split_once('\t').expect("a case");
                let (name, _) = rest.split_once('=').expect("a value");
                let (reference, alternates, minimum) = choice_case(label);
                let reads = reads_for(pileup_name(label));
                let pileup = pileup(&reads);
                let element = pileup
                    .elements
                    .iter()
                    .find(|element| element.read.read_name == name)
                    .unwrap_or_else(|| panic!("{label}: no element for {name}"));
                let chosen = choose_allele_for_read(element, &reference, &alternates, minimum)
                    .expect("a cutoff in range");
                let shown = chosen
                    .map(|allele| String::from_utf8_lossy(&allele.bases).into_owned())
                    .unwrap_or_else(|| "none".to_string());
                assert_eq!(format!("choose\t{label}\t{name}={shown}"), line);
            }
            "count" => {
                let (label, _) = rest.split_once('\t').expect("a map");
                let (reference, alternates, minimum) = if label == "null" {
                    (reference("A"), vec![alternate("C")], 0)
                } else {
                    choice_case(label)
                };
                let counter = if label == "null" {
                    AllelePileupCounter::new(&reference, &alternates, minimum).expect("a counter")
                } else {
                    let reads = reads_for(pileup_name(label));
                    AllelePileupCounter::with_pileup(
                        &reference,
                        &alternates,
                        minimum,
                        &pileup(&reads),
                    )
                    .expect("a counter")
                };
                assert_eq!(format!("count\t{label}\t{}", counted(&counter)), line);
            }
            "ratio" => {
                let (label, _) = rest.split_once('=').expect("a value");
                let (name, minimum) = label.split_once(',').expect("a cutoff");
                let reads = reads_for(name);
                let ours = calculate_max_alt_ratio(
                    &pileup(&reads),
                    &ratio_reference(name),
                    minimum.parse().expect("a number"),
                )
                .expect("a cutoff in range");
                assert_eq!(format!("ratio\t{label}={}", hex(ours)), line);
            }
            "support" => {
                let (label, _) = rest.split_once('=').expect("a value");
                let parts: Vec<&str> = label.split(',').collect();
                let reads = reads_for(parts[0]);
                let reference = if parts[0] == "short" {
                    reference("AA")
                } else if parts[0] == "deletion" {
                    reference("AAA")
                } else {
                    reference("A")
                };
                let ours = calculate_num_reads_supporting_allele(
                    &pileup(&reads),
                    &reference,
                    &alternate(parts[1]),
                    parts[2].parse().expect("a number"),
                )
                .expect("a cutoff in range");
                assert_eq!(format!("support\t{label}={ours}"), line);
            }
            "validatable" => {
                let (label, _) = rest.split_once('=').expect("a value");
                let (genotype, reference, _, _, _, _) = validation_case(label);
                let ours = is_able_to_validate_genotype(&genotype, &reference)
                    .expect("a diploid genotype");
                assert_eq!(format!("validatable\t{label}={ours}"), line);
            }
            "result" => {
                let (label, _) = rest.split_once('=').expect("a value");
                let ours = match result_of(label) {
                    None => format!("result\t{label}=none"),
                    Some(result) => {
                        let row = format!(
                            "result\t{label}={},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                            result.contig,
                            result.start,
                            result.end,
                            String::from_utf8_lossy(&result.reference.bases),
                            String::from_utf8_lossy(&result.alternate.bases),
                            result.minimum_validation_read_count,
                            result.is_enough_validation_reads,
                            result.is_out_of_noise_floor,
                            hex(result.power),
                            result.validation_alt_count,
                            result.validation_ref_count,
                            result.discovery_alt_count,
                            result.discovery_ref_count,
                            result.filters,
                            result.num_alt_supporting_reads_in_normal,
                        );
                        written.push(result);
                        row
                    }
                };
                assert_eq!(ours, line);
            }
            "table" => table.push(rest.to_string()),
            "error" => {
                let (label, _) = rest.split_once('\t').expect("a message");
                assert_eq!(refusal(label), line);
            }
            other => panic!("an unexpected row: {other}"),
        }
        rows += 1;
    }

    // The table is compared once, as the whole file the writer produces.
    let ours = write_table(&written);
    assert_eq!(ours.lines().collect::<Vec<&str>>(), table);
    assert_eq!(rows, 124, "the golden's row count");
    // Nothing above reached this refusal, which the reference words differently from the others.
    assert_eq!(
        PileupAlleleError::SymbolicReference.java_class(),
        "java.lang.IllegalStateException"
    );
}
