//! Conformance for the three allele-specific rank sums and the histogram they travel as, against
//! the oracle.
//!
//! Golden from `tools/annotation-conformance/AlleleSpecificRankSumDump.java`.
//!
//! ```text
//! bin     -1.23   -1.3,1  null
//! asraw   AS_BaseQualityRankSumTest  ref-only  AS_RAW_BaseQRankSum=|NaN[java.lang.String]
//! asfinal AS_BaseQualityRankSumTest  empty-slots  AS_BaseQRankSum=.,.[...];AS_RAW_...=||[...]
//! ```
//!
//! The `bin` rows are the binning rule one value at a time, which is what makes the claim that a
//! negative score bins **away** from zero a measurement rather than a reading of the source.

use std::io::Read;

use gatk_annotation::allele_specific_rank_sum::{
    annotate_raw_data, combine_raw_data, finalize_raw_data, AsRankSum, AsRankSumError,
};
use gatk_annotation::rank_sum;
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::histogram::{CompressedDataList, Histogram};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

const START: i64 = 105;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/allele_specific_rank_sum.txt.gz");
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

fn histogram_values(label: &str) -> Vec<f64> {
    match label {
        "empty" => vec![],
        "one-value" => vec![0.5],
        "two-values" => vec![0.1, 0.2],
        "three-values" => vec![0.1, 0.2, 0.3],
        "four-values" => vec![0.1, 0.2, 0.3, 0.4],
        "negative" => vec![-1.5, -2.5, -3.5],
        "straddling-zero" => vec![-0.2, -0.1, 0.1, 0.2],
        "on-a-boundary" => vec![0.1, 0.2, 0.300_000_000_000_000_04],
        "just-below-a-boundary" => vec![0.0999, 0.1999],
        "repeated" => vec![0.5, 0.5, 0.5, 0.5, 0.5],
        "wide-bins" => vec![1.0, 2.0, 3.0],
        other => panic!("{other} has no histogram fixture"),
    }
}

fn compressed_values(label: &str) -> Vec<i32> {
    match label {
        "empty" => vec![],
        "ascending" => vec![1, 2, 3],
        "descending" => vec![3, 2, 1],
        "repeated" => vec![2, 2, 2, 5, 5],
        "negative" => vec![-3, 0, 3],
        other => panic!("{other} has no list fixture"),
    }
}

fn read(name: &str, mapping_quality: u8, base_quality: u8) -> BamRecord {
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: 100,
        mapping_quality: mapping_quality.min(60),
        cigar: htsjdk_bam::text_parse::parse_cigar("20M").expect("a cigar"),
        read_bases: vec![b'A'; 20],
        base_qualities: vec![base_quality.min(60); 20],
        ..Default::default()
    }
}

fn variant_context(label: &str) -> VariantContext {
    let alleles = if label == "multiallelic" {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    };
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    match label {
        "no-genotypes" => {}
        "two-samples" => {
            vc.genotypes
                .push(Genotype::new("s1", vec![reference(), alternate()]));
            vc.genotypes
                .push(Genotype::new("s2", vec![reference(), alternate()]));
        }
        _ => vc
            .genotypes
            .push(Genotype::new("s1", vec![reference(), alternate()])),
    }
    vc
}

fn likelihoods(label: &str) -> Option<AlleleLikelihoods<BamRecord>> {
    if label == "null-likelihoods" {
        return None;
    }
    let count = if label == "single-read" { 1 } else { 12 };
    let reads: Vec<BamRecord> = (0..count)
        .map(|i| read(&format!("r{i}"), (20 + i * 3) as u8, (30 + i) as u8))
        .collect();
    let alleles = if label == "multiallelic" {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    };
    let values: Vec<Vec<f64>> = (0..alleles.len())
        .map(|a| {
            (0..reads.len())
                .map(|e| {
                    let best = match label {
                        "ref-only" => 0,
                        "alt-only" => 1,
                        "multiallelic" => e % 3,
                        _ => {
                            if e < reads.len() / 2 {
                                0
                            } else {
                                1
                            }
                        }
                    };
                    let strong = if label == "overlapping" {
                        -1.0 - (e as f64 * 0.1)
                    } else {
                        -1.0
                    };
                    let weak = if label == "overlapping" {
                        -5.0 - (e as f64 * 0.1)
                    } else {
                        -10.0
                    };
                    if a == best {
                        strong
                    } else {
                        weak
                    }
                })
                .collect()
        })
        .collect();
    Some(
        AlleleLikelihoods::new(
            SampleList::new(&["s1".to_string()]),
            AlleleList::new(&alleles),
            vec![reads],
            vec![values],
        )
        .expect("a matrix"),
    )
}

fn annotation(name: &str) -> AsRankSum {
    match name {
        "AS_BaseQualityRankSumTest" => AsRankSum::BaseQuality,
        "AS_MappingQualityRankSumTest" => AsRankSum::MappingQuality,
        "AS_ReadPosRankSumTest" => AsRankSum::ReadPosition,
        other => panic!("unknown annotation {other}"),
    }
}

/// The dump's `combineAndFinalize` inputs.
fn combine_case(label: &str) -> Vec<String> {
    let raws: Vec<&str> = match label {
        "one-source" => vec!["|-1.3,1|0.4,1"],
        "two-sources" => vec!["|-1.3,1|0.4,1", "|-1.3,1|0.5,2"],
        "even-count" => vec!["|1.0,1", "|2.0,1"],
        "odd-count" => vec!["|1.0,1", "|2.0,2"],
        "empty-slots" => vec!["||"],
        "bracketed" => vec!["[|-1.3,1|0.4,1]"],
        "missing-alt" => vec!["|-1.3,1"],
        other => panic!("{other} has no combine fixture"),
    };
    raws.into_iter().map(|s| s.to_string()).collect()
}

/// The dump's allele count, taken from the first raw string exactly as the dump takes it.
fn combine_alleles(label: &str) -> Vec<Allele> {
    let first = combine_case(label).remove(0);
    let stripped = first.replace(['[', ']'], "");
    // `split("\\|", -1)`, which keeps trailing empties.
    let slots = stripped.split('|').count();
    if slots > 2 {
        vec![reference(), alternate(), second_alternate()]
    } else {
        vec![reference(), alternate()]
    }
}

/// `key=value[class]` joined with `;`, keys sorted, as the dump emits them.
fn rendered(entries: &[(String, String)]) -> String {
    let mut sorted: Vec<&(String, String)> = entries.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));
    sorted
        .iter()
        .map(|(key, value)| format!("{key}={value}[java.lang.String]"))
        .collect::<Vec<_>>()
        .join(";")
}

#[test]
fn every_histogram_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("hist\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let label = fields[0];
        let mut histogram = if label == "wide-bins" {
            Histogram::with_bin_size(0.01)
        } else {
            Histogram::new()
        };
        for value in histogram_values(label) {
            histogram.add(value).expect("a bin");
        }
        assert_eq!(histogram.to_string(), fields[1], "rendering of {label}");
        let median = match histogram.median() {
            None => "null".to_string(),
            // `Double.toString`, which for every value here is the shortest round-tripping
            // decimal and matches Rust's `{}` for a finite non-integral double.
            Some(value) => {
                if value == value.trunc() && value.abs() < 1e7 {
                    format!("{value:.1}")
                } else {
                    format!("{value}")
                }
            }
        };
        assert_eq!(median, fields[2], "median of {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no histogram rows");
    println!("{count} histograms identical");
}

#[test]
fn every_compressed_list_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("cdl\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.splitn(3, '\t').collect();
        let label = fields[0];
        let mut list = CompressedDataList::new();
        for value in compressed_values(label) {
            list.add(value);
        }
        assert_eq!(list.to_string(), fields[1], "rendering of {label}");
        let iterated = list
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        assert_eq!(
            iterated,
            fields.get(2).copied().unwrap_or(""),
            "iteration of {label}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no list rows");
    println!("{count} compressed lists identical");
}

#[test]
fn every_bin_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("bin\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let value: f64 = fields[0].parse().expect("a double");
        let mut histogram = Histogram::new();
        histogram.add(value).expect("a bin");
        assert_eq!(histogram.to_string(), fields[1], "bin of {}", fields[0]);
        let at_zero = match histogram.get(0.0).expect("a valid key") {
            None => "null".to_string(),
            Some(count) => count.to_string(),
        };
        assert_eq!(at_zero, fields[2], "count at bin zero for {}", fields[0]);
        count += 1;
    }
    assert!(count > 0, "the golden carries no bin rows");
    println!("{count} binnings identical");
}

#[test]
fn every_allele_specific_answer_matches_the_reference() {
    let text = golden();
    let mut direct = 0;
    let mut raw = 0;
    for line in text.lines() {
        let kind = if line.starts_with("as\t") {
            "as"
        } else if line.starts_with("asraw\t") {
            "asraw"
        } else {
            continue;
        };
        let rest = &line[kind.len() + 1..];
        let fields: Vec<&str> = rest.splitn(3, '\t').collect();
        let name = fields[0];
        let label = fields[1];
        let expected = fields.get(2).copied().unwrap_or("");
        let annotation = annotation(name);
        let vc = variant_context(label);
        let matrix = likelihoods(label);

        let ours = if kind == "as" {
            direct += 1;
            rank_sum::annotate(&annotation, None, &vc, matrix.as_ref())
                .iter()
                .map(|(key, value)| {
                    let text = value.to_java_string().expect("a renderable value");
                    format!("{key}={text}[java.lang.String]")
                })
                .collect::<Vec<_>>()
                .join(";")
        } else {
            raw += 1;
            match annotate_raw_data(annotation, &vc, matrix.as_ref()) {
                Ok(None) => String::new(),
                Ok(Some((key, value))) => rendered(&[(key, value)]),
                // The reference's `IllegalStateException` when the sample count is not one.
                Err(AsRankSumError::NotExactlyOneSample { .. }) => {
                    "E:java.lang.IllegalStateException".to_string()
                }
                Err(other) => panic!("{name} on {label}: {other:?}"),
            }
        };
        assert_eq!(ours, expected, "{kind} {name} on {label}");
    }
    assert!(
        direct > 0 && raw > 0,
        "the golden carries no annotation rows"
    );
    println!("{direct} direct and {raw} raw allele-specific answers identical");
}

#[test]
fn every_combination_and_finalisation_matches_the_reference() {
    let text = golden();
    let mut combined = 0;
    let mut finalised = 0;
    for line in text.lines() {
        let kind = if line.starts_with("ascombine\t") {
            "ascombine"
        } else if line.starts_with("asfinal\t") {
            "asfinal"
        } else {
            continue;
        };
        let rest = &line[kind.len() + 1..];
        let fields: Vec<&str> = rest.splitn(3, '\t').collect();
        let label = fields[1];
        let expected = fields.get(2).copied().unwrap_or("");
        let annotation = AsRankSum::BaseQuality;
        let alleles = combine_alleles(label);
        let text_combined =
            combine_raw_data(&alleles, &combine_case(label)).expect("a combination");

        let ours = if kind == "ascombine" {
            combined += 1;
            rendered(&[(annotation.raw_key().to_string(), text_combined)])
        } else {
            finalised += 1;
            let alternates: Vec<Allele> = alleles
                .iter()
                .filter(|a| !a.is_reference())
                .cloned()
                .collect();
            match finalize_raw_data(annotation, &alternates, &alleles, Some(&text_combined))
                .expect("a finalisation")
            {
                None => String::new(),
                Some((key, reduced, raw_key, raw_value)) => {
                    rendered(&[(key, reduced), (raw_key, raw_value)])
                }
            }
        };
        assert_eq!(ours, expected, "{kind} on {label}");
    }
    assert!(
        combined > 0 && finalised > 0,
        "the golden carries no combination rows"
    );
    println!("{combined} combinations and {finalised} finalisations identical");
}
