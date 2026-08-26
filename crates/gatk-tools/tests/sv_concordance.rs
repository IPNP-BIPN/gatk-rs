//! Conformance for `SVConcordance` against GATK 4.6.2.0, compared as every concordance field of
//! every eval record of every run.
//!
//! Golden from `tools/readfilter-conformance/SVConcordanceDump.java`.
//!
//! # What this suite is for
//!
//!  * **each eval record taking one truth record**, the closest by total breakend distance;
//!  * **the tiebreaker being the closest single breakend**;
//!  * **the CNV override being asymmetric**, so an eval DEL matches a truth CNV and not the other
//!    way round;
//!  * **an unmatched eval record still being written**, with no truth allele counts at all;
//!  * **a multiallelic CNV being scored on copy state instead**;
//!  * **concordance being computed on the common samples only**;
//!  * **a truth record's allele counts being copied when it has them and recounted when it does
//!    not**, which the formatting gives away;
//!  * **and --do-not-sort changing nothing for an in-order input.**

use gatk_corpus as corpus;
use gatk_tools::sv_cluster::{CallRecord, Linkage};
use gatk_tools::sv_concordance::{
    annotate, are_clusterable, closest, count_alleles, min_distance, run, total_distance,
    AlleleCounts, Annotation, Genotype, Record,
};
use gatk_tools::sv_stratify::SvType;
use htsjdk_vcf::variant::format_vcf_double;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sv_concordance.txt.gz"),
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

/// One input VCF, read as the records the tool builds from it.
fn records(text: &str, which: &str) -> Vec<Record> {
    let vcf = section(text, "vcf", which);
    let mut samples: Vec<String> = Vec::new();
    let mut out = Vec::new();
    for line in vcf.lines() {
        if line.starts_with("#CHROM") {
            samples = line.split('\t').skip(9).map(str::to_string).collect();
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        let info: Vec<(&str, &str)> = columns[7]
            .split(';')
            .filter_map(|part| part.split_once('='))
            .collect();
        let field = |key: &str| {
            info.iter()
                .find(|(name, _)| *name == key)
                .map(|(_, value)| *value)
        };
        let sv_type = SvType::parse(field("SVTYPE").expect("a type")).expect("a known type");
        let end: i32 = field("END").expect("an end").parse().expect("an end");
        let keys: Vec<&str> = columns[8].split(':').collect();
        let index_of = |key: &str| keys.iter().position(|name| *name == key);
        let genotype_index = index_of("GT").expect("a GT");
        let copy_number_index = index_of("CN");
        let genotypes: Vec<Genotype> = samples
            .iter()
            .enumerate()
            .map(|(index, sample)| {
                let parts: Vec<&str> = columns[9 + index].split(':').collect();
                Genotype {
                    sample: sample.clone(),
                    alleles: parts[genotype_index]
                        .split(['/', '|'])
                        .map(|allele| allele.parse::<i32>().ok())
                        .collect(),
                    copy_number: copy_number_index
                        .and_then(|at| parts.get(at))
                        .and_then(|value| value.parse::<i32>().ok()),
                }
            })
            .collect();
        // All three must be present for the record to count as carrying them, and the TEXT is kept
        // because a copied value is written back exactly as it was read.
        let allele_counts = match (field("AC"), field("AF"), field("AN")) {
            (Some(count), Some(frequency), Some(number)) => Some(AlleleCounts {
                count: count
                    .split(',')
                    .map(|v| v.parse().expect("a count"))
                    .collect(),
                frequency: frequency
                    .split(',')
                    .map(|v| v.parse().expect("a frequency"))
                    .collect(),
                number: number.parse().expect("a number"),
                verbatim: Some((count.to_string(), frequency.to_string(), number.to_string())),
            }),
            _ => None,
        };
        out.push(Record {
            call: CallRecord {
                id: columns[2].to_string(),
                sv_type,
                contig_a: columns[0].to_string(),
                position_a: columns[1].parse().expect("a position"),
                contig_b: columns[0].to_string(),
                position_b: end,
                strand_a: None,
                strand_b: None,
                length: match sv_type {
                    SvType::Ins | SvType::Bnd | SvType::Ctx => None,
                    _ => Some(field("SVLEN").expect("a length").parse().expect("a length")),
                },
                algorithms: field("ALGORITHMS")
                    .expect("an algorithms field")
                    .split(',')
                    .map(str::to_string)
                    .collect(),
                carriers: genotypes
                    .iter()
                    .filter(|genotype| {
                        genotype
                            .alleles
                            .iter()
                            .any(|allele| matches!(allele, Some(index) if *index > 0))
                    })
                    .map(|genotype| genotype.sample.clone())
                    .collect(),
            },
            genotypes,
            allele_counts,
        });
    }
    out
}

/// The keys the annotator writes, which are the ones this suite compares. Everything else in the
/// INFO column came straight off the input record.
const KEYS: &[&str] = &[
    "AC",
    "AF",
    "AN",
    "STATUS",
    "TRUTH_VID",
    "TRUTH_AC",
    "TRUTH_AF",
    "TRUTH_AN",
    "GENOTYPE_CONCORDANCE",
    "NON_REF_GENOTYPE_CONCORDANCE",
    "HET_PPV",
    "HET_SENSITIVITY",
    "HOMVAR_PPV",
    "HOMVAR_SENSITIVITY",
    "VAR_PPV",
    "VAR_SENSITIVITY",
    "VAR_SPECIFICITY",
    "CNV_CONCORDANCE",
];

/// One eval record's measured output: the concordance INFO fields, and the per-sample values.
#[derive(Debug, PartialEq)]
struct Measured {
    id: String,
    info: Vec<(String, String)>,
    contingency: Vec<(String, Option<String>)>,
    copy_number_equal: Vec<(String, Option<String>)>,
}

fn measured(text: &str, label: &str) -> Vec<Measured> {
    let body = section(text, "out", label);
    let mut lines = body.lines();
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let samples: Vec<String> = header[9..].iter().map(|s| s.to_string()).collect();
    lines
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let mut info: Vec<(String, String)> = columns[7]
                .split(';')
                .filter_map(|part| part.split_once('='))
                .filter(|(key, _)| KEYS.contains(key))
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .collect();
            info.sort();
            let keys: Vec<&str> = columns[8].split(':').collect();
            let at = |key: &str| keys.iter().position(|name| *name == key);
            let read = |index: Option<usize>, sample: usize| -> Option<String> {
                index
                    .and_then(|position| columns[9 + sample].split(':').nth(position))
                    .map(str::to_string)
            };
            let contingency_at = at("CONC_ST");
            let copy_number_at = at("TRUTH_CN_EQUAL");
            Measured {
                id: columns[2].to_string(),
                info,
                contingency: samples
                    .iter()
                    .enumerate()
                    .map(|(index, sample)| {
                        (
                            sample.clone(),
                            read(contingency_at, index).filter(|value| value != "."),
                        )
                    })
                    .collect(),
                copy_number_equal: samples
                    .iter()
                    .enumerate()
                    .map(|(index, sample)| {
                        (
                            sample.clone(),
                            read(copy_number_at, index).filter(|value| value != "."),
                        )
                    })
                    .collect(),
            }
        })
        .collect()
}

/// A metric is dropped from the INFO column when it is NaN.
fn push_metric(info: &mut Vec<(String, String)>, key: &str, value: f64) {
    if !value.is_nan() {
        info.push((key.to_string(), format_vcf_double(value)));
    }
}

fn push_counts(info: &mut Vec<(String, String)>, prefix: &str, counts: &AlleleCounts) {
    let (count, frequency, number) = match &counts.verbatim {
        Some(text) => text.clone(),
        None => (
            counts
                .count
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<String>>()
                .join(","),
            counts
                .frequency
                .iter()
                .map(|v| format_vcf_double(*v))
                .collect::<Vec<String>>()
                .join(","),
            counts.number.to_string(),
        ),
    };
    info.push((format!("{prefix}AC"), count));
    info.push((format!("{prefix}AF"), frequency));
    info.push((format!("{prefix}AN"), number));
}

/// One annotation rendered the way the writer renders it.
fn rendered(annotation: &Annotation) -> Measured {
    let mut info: Vec<(String, String)> =
        vec![("STATUS".to_string(), annotation.status.to_string())];
    if let Some(id) = &annotation.truth_variant_id {
        info.push(("TRUTH_VID".to_string(), id.clone()));
    }
    if let Some(counts) = &annotation.allele_counts {
        push_counts(&mut info, "", counts);
    }
    if let Some(counts) = &annotation.truth_allele_counts {
        push_counts(&mut info, "TRUTH_", counts);
    }
    if let Some(concordance) = annotation.copy_number_concordance {
        push_metric(&mut info, "CNV_CONCORDANCE", concordance);
    }
    if let Some(metrics) = &annotation.metrics {
        push_metric(
            &mut info,
            "GENOTYPE_CONCORDANCE",
            metrics.genotype_concordance,
        );
        push_metric(
            &mut info,
            "NON_REF_GENOTYPE_CONCORDANCE",
            metrics.non_ref_genotype_concordance,
        );
        push_metric(&mut info, "HET_PPV", metrics.het_ppv);
        push_metric(&mut info, "HET_SENSITIVITY", metrics.het_sensitivity);
        push_metric(&mut info, "HOMVAR_PPV", metrics.homvar_ppv);
        push_metric(&mut info, "HOMVAR_SENSITIVITY", metrics.homvar_sensitivity);
        push_metric(&mut info, "VAR_PPV", metrics.var_ppv);
        push_metric(&mut info, "VAR_SENSITIVITY", metrics.var_sensitivity);
        push_metric(&mut info, "VAR_SPECIFICITY", metrics.var_specificity);
    }
    info.sort();
    Measured {
        id: annotation.id.clone(),
        info,
        contingency: annotation.contingency.clone(),
        copy_number_equal: annotation
            .truth_copy_number_equal
            .iter()
            .map(|(sample, value)| {
                (
                    sample.clone(),
                    value.map(|matched| if matched { "1" } else { "0" }.to_string()),
                )
            })
            .collect(),
    }
}

fn common() -> Vec<String> {
    vec!["s2".to_string(), "s3".to_string()]
}

#[test]
fn every_annotation_matches_the_golden() {
    let text = golden();
    let eval = records(&text, "eval");
    let truth = records(&text, "truth");
    let mut compared = 0;
    for (label, overlap) in [
        ("default", 0.8),
        ("do-not-sort", 0.8),
        ("overlap-high", 0.999),
    ] {
        let mut linkage = Linkage::default();
        linkage.depth.reciprocal_overlap = overlap;
        linkage.mixed.reciprocal_overlap = overlap;
        linkage.pesr.reciprocal_overlap = overlap;
        let produced: Vec<Measured> = run(&linkage, &eval, &truth, &common())
            .iter()
            .map(rendered)
            .collect();
        assert_eq!(produced, measured(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 3, "the runs that wrote records");
}

/// The comment says CNV/DEL matching is not allowed. The condition only refuses it when the CNV is
/// the EVAL record, so the same pair of types matches one way round and not the other.
#[test]
fn the_cnv_override_is_asymmetric() {
    let text = golden();
    let eval = records(&text, "eval");
    let truth = records(&text, "truth");
    let of = |records: &[Record], id: &str| {
        records
            .iter()
            .find(|record| record.call.id == id)
            .unwrap_or_else(|| panic!("{id}"))
            .clone()
    };
    let eval_cnv = of(&eval, "eval-cnv");
    let truth_del = of(&truth, "t-eval-cnv");
    let eval_del = of(&eval, "eval-del");
    let truth_cnv = of(&truth, "t-eval-del");
    let linkage = Linkage::default();

    // The eval CNV against a truth DEL is refused.
    assert_eq!(eval_cnv.call.sv_type, SvType::Cnv);
    assert_eq!(truth_del.call.sv_type, SvType::Del);
    assert!(!are_clusterable(&linkage, &eval_cnv.call, &truth_del.call));

    // The eval DEL against a truth CNV is not.
    assert_eq!(eval_del.call.sv_type, SvType::Del);
    assert_eq!(truth_cnv.call.sv_type, SvType::Cnv);
    assert!(are_clusterable(&linkage, &eval_del.call, &truth_cnv.call));

    // And the base linkage would have allowed BOTH, so the override is what makes the difference
    // in one direction and nothing in the other.
    assert!(linkage.are_clusterable(&eval_cnv.call, &truth_del.call));
    assert!(linkage.are_clusterable(&eval_del.call, &truth_cnv.call));

    // Which the golden shows as a false positive on one and a true positive on the other.
    let out = measured(&text, "default");
    let status = |id: &str| {
        out.iter()
            .find(|record| record.id == id)
            .expect(id)
            .info
            .iter()
            .find(|(key, _)| key == "STATUS")
            .expect("a status")
            .1
            .clone()
    };
    assert_eq!(status("eval-cnv"), "FP");
    assert_eq!(status("eval-del"), "TP");
}

/// Two truth records at the same total distance are separated by the smaller of their two ends.
#[test]
fn the_tiebreaker_is_the_closest_single_breakend() {
    let text = golden();
    let eval = records(&text, "eval");
    let truth = records(&text, "truth");
    let of = |records: &[Record], id: &str| {
        records
            .iter()
            .find(|record| record.call.id == id)
            .unwrap_or_else(|| panic!("{id}"))
            .clone()
    };
    let tie = of(&eval, "tie");
    let both = of(&truth, "t-tie-both");
    let one = of(&truth, "t-tie-one");
    assert_eq!(
        total_distance(&tie.call, &both.call),
        total_distance(&tie.call, &one.call),
        "the totals are equal, which is what makes it a tie"
    );
    assert_eq!(min_distance(&tie.call, &one.call), 0);
    assert_eq!(min_distance(&tie.call, &both.call), 100);
    assert_eq!(
        closest(&Linkage::default(), &tie, &truth)
            .expect("a match")
            .call
            .id,
        "t-tie-one"
    );
}

/// A sample present in one VCF and not the other is annotated on neither side.
#[test]
fn concordance_is_computed_on_the_common_samples_only() {
    let text = golden();
    let out = measured(&text, "default");
    for record in &out {
        let (sample, value) = &record.contingency[0];
        assert_eq!(sample, "s1", "the eval-only sample comes first");
        assert_eq!(*value, None, "{}", record.id);
    }
    // And it is not counted either: annotating with s1 in the common set moves the numbers.
    let eval = records(&text, "eval");
    let truth = records(&text, "truth");
    fn of<'a>(records: &'a [Record], id: &str) -> &'a Record {
        records
            .iter()
            .find(|record| record.call.id == id)
            .unwrap_or_else(|| panic!("{id}"))
    }
    let exact = of(&eval, "exact");
    let matched = Some(of(&truth, "t-exact"));
    let narrow = annotate(exact, matched, &common());
    let wide = annotate(
        exact,
        matched,
        &["s1".to_string(), "s2".to_string(), "s3".to_string()],
    );
    assert_ne!(
        narrow.metrics.as_ref().map(|m| m.genotype_concordance),
        wide.metrics.as_ref().map(|m| m.genotype_concordance)
    );
}

/// A truth record carrying AC, AF and AN has them copied verbatim, and one without has them
/// recounted. The formatting gives it away: the copied value keeps the input's own text.
#[test]
fn a_truth_records_allele_counts_are_copied_or_recounted() {
    let text = golden();
    let out = measured(&text, "default");
    let field = |id: &str, key: &str| {
        out.iter()
            .find(|record| record.id == id)
            .expect(id)
            .info
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.clone())
    };
    assert_eq!(field("af-given", "TRUTH_AF").as_deref(), Some("0.5"));
    assert_eq!(field("af-missing", "TRUTH_AF").as_deref(), Some("0.500"));
    assert_eq!(
        format_vcf_double(0.5),
        "0.500",
        "the recounted value goes through the writer's formatter"
    );

    // And an unmatched record has no truth counts at all rather than empty ones.
    assert_eq!(field("nomatch", "TRUTH_AC"), None);
    assert_eq!(field("nomatch", "TRUTH_AF"), None);
    assert_eq!(field("nomatch", "TRUTH_AN"), None);

    // The recount is over the TRUTH genotypes but against the EVAL record's alternate alleles, so
    // a truth record with no called allele at all yields a count of zero over a number of zero.
    assert_eq!(field("eval-del", "TRUTH_AN").as_deref(), Some("0"));
    assert_eq!(field("eval-del", "TRUTH_AF").as_deref(), Some("NaN"));
    let none_called = count_alleles(
        1,
        &[Genotype {
            sample: "s2".to_string(),
            alleles: vec![None, None],
            copy_number: Some(2),
        }],
    );
    assert_eq!(none_called.number, 0);
    assert!(none_called.frequency[0].is_nan());
}

/// It gets a copy-state answer per sample and one number over the record, and no genotype
/// concordance at all.
#[test]
fn a_multiallelic_cnv_is_scored_on_copy_state_instead() {
    let text = golden();
    let out = measured(&text, "default");
    let record = out.iter().find(|r| r.id == "cnv-pair").expect("cnv-pair");
    assert_eq!(
        record.copy_number_equal,
        vec![
            ("s1".to_string(), None),
            ("s2".to_string(), Some("1".to_string())),
            ("s3".to_string(), Some("0".to_string())),
        ]
    );
    assert!(record
        .info
        .contains(&("CNV_CONCORDANCE".to_string(), "0.500".to_string())));
    for key in ["GENOTYPE_CONCORDANCE", "AC", "AF", "AN"] {
        assert!(
            !record.info.iter().any(|(name, _)| name == key),
            "a CNV carries no {key}"
        );
    }
    // The unmatched CNV gets no copy-state answer on any sample, because there is nothing to
    // compare against.
    let unmatched = out.iter().find(|r| r.id == "eval-cnv").expect("eval-cnv");
    assert!(unmatched
        .copy_number_equal
        .iter()
        .all(|(_, value)| value.is_none()));
    assert!(!unmatched
        .info
        .iter()
        .any(|(name, _)| name == "CNV_CONCORDANCE"));
}

/// The flag removes a buffer that sorts by position, and the records complete in the order they
/// were added, so for an in-order input it changes nothing.
#[test]
fn do_not_sort_changes_nothing_for_an_in_order_input() {
    let text = golden();
    assert_eq!(measured(&text, "default"), measured(&text, "do-not-sort"));
    // The control that says the comparison can see a difference at all.
    assert_ne!(measured(&text, "default"), measured(&text, "overlap-high"));
}

/// Raising it to 0.999 breaks exactly one match: the tie pair is the only one whose truth record
/// is not an exact interval.
#[test]
fn a_higher_reciprocal_overlap_breaks_one_match() {
    let text = golden();
    let default = measured(&text, "default");
    let high = measured(&text, "overlap-high");
    let changed: Vec<&String> = default
        .iter()
        .zip(&high)
        .filter(|(a, b)| a != b)
        .map(|(a, _)| &a.id)
        .collect();
    assert_eq!(changed, vec!["tie"]);
}

/// The dictionary comes from the argument, not from either VCF, and the finder refuses an input
/// that is not in position order.
#[test]
fn the_two_refusals() {
    let text = golden();
    assert_eq!(
        refusal(&text, "no-dictionary"),
        (
            "org.broadinstitute.hellbender.exceptions.UserException".to_string(),
            "Reference sequence dictionary required".to_string()
        )
    );
    assert_eq!(
        refusal(&text, "unsorted"),
        (
            "java.lang.IllegalArgumentException".to_string(),
            "Items must be added in dictionary-sorted order".to_string()
        )
    );
}
