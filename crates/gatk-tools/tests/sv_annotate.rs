//! Conformance for `SVAnnotate` against GATK 4.6.2.0, compared as every `PREDICTED_` field of
//! every record of every run.
//!
//! Golden from `tools/readfilter-conformance/SVAnnotateDump.java`.
//!
//! # What this suite is for
//!
//!  * **spanning and overlapping being different questions**, which is the whole duplication rule;
//!  * **one breakpoint in coding sequence with one in a UTR being LOF**;
//!  * **the start site being the transcript's END on the minus strand**;
//!  * **a breakend using the simple rule and then downgrading its LOF**, and producing two
//!    segments;
//!  * **a CNV being annotated as a duplication and then reclassified**;
//!  * **the promoter being inferred from a window**, and a promoter also being intergenic;
//!  * **two arguments being accepted and then crashing**;
//!  * **and the non-coding BED refusing the header its own documentation asks for.**

use gatk_corpus as corpus;
use gatk_tools::sv_annotate::{
    annotate_breakend, annotate_copy_number_variant, annotate_duplication,
    annotate_structural_variant, annotation_type_for_breakend, count_breakends_inside_feature,
    promoter_interval, sv_segments, AnnotateError, Annotation, ComplexSubtype, Feature,
    FeatureKind, Interval, NonCodingElement, SvType, Transcript, Variant,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sv_annotate.txt.gz"),
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

/// The GTF the golden carries, read as one transcript per `transcript` line with the features that
/// follow it.
fn transcripts(text: &str) -> Vec<Transcript> {
    let mut out: Vec<Transcript> = Vec::new();
    for line in section(text, "gtf", "main").lines() {
        if line.is_empty() {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        let kind = match columns[2] {
            "transcript" => FeatureKind::Transcript,
            "exon" => FeatureKind::Exon,
            "CDS" => FeatureKind::Cds,
            "start_codon" => FeatureKind::StartCodon,
            "stop_codon" => FeatureKind::StopCodon,
            "UTR" => FeatureKind::Utr,
            // The gene record, which owns no features of its own.
            _ => continue,
        };
        let start: i32 = columns[3].parse().expect("a start");
        let end: i32 = columns[4].parse().expect("an end");
        if kind == FeatureKind::Transcript {
            let attributes = columns[8];
            let gene_name = attributes
                .split("gene_name \"")
                .nth(1)
                .and_then(|rest| rest.split('"').next())
                .expect("a gene name")
                .to_string();
            out.push(Transcript {
                gene_name,
                contig: columns[0].to_string(),
                start,
                end,
                negative_strand: columns[6] == "-",
                // `getAllFeatures` includes the transcript record itself.
                features: vec![Feature { kind, start, end }],
            });
        } else {
            out.last_mut()
                .expect("a transcript first")
                .features
                .push(Feature { kind, start, end });
        }
    }
    out
}

/// The BED the golden carries. BED is half-open, so a start of 10399 is base 10400.
fn non_coding(text: &str) -> Vec<NonCodingElement> {
    section(text, "bed", "main")
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            NonCodingElement {
                name: columns[3].to_string(),
                interval: Interval {
                    contig: columns[0].to_string(),
                    start: columns[1].parse::<i32>().expect("a start") + 1,
                    end: columns[2].parse().expect("an end"),
                },
            }
        })
        .collect()
}

fn variants(text: &str) -> Vec<Variant> {
    section(text, "vcf", "input")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
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
            Variant {
                id: columns[2].to_string(),
                contig: columns[0].to_string(),
                position: columns[1].parse().expect("a position"),
                end: field("END").expect("an end").parse().expect("an end"),
                sv_type: match field("SVTYPE").expect("a type") {
                    "DEL" => SvType::Del,
                    "DUP" => SvType::Dup,
                    "CNV" => SvType::Cnv,
                    "INS" => SvType::Ins,
                    "INV" => SvType::Inv,
                    "BND" => SvType::Bnd,
                    "CTX" => SvType::Ctx,
                    "CPX" => SvType::Cpx,
                    other => panic!("an unexpected type {other}"),
                },
                sv_length: field("SVLEN").expect("a length").parse().expect("a length"),
                contig2: field("CHR2").map(str::to_string),
                end2: field("END2").map(|value| value.parse().expect("a second position")),
                strands: field("STRANDS").map(str::to_string),
                complex_type: None,
                complex_intervals: Vec::new(),
            }
        })
        .collect()
}

/// One record's measured annotation: the id, its `PREDICTED_` fields, and whether the intergenic
/// flag was written.
#[derive(Debug, PartialEq, Eq)]
struct Measured {
    id: String,
    consequences: Vec<(String, Vec<String>)>,
    intergenic: bool,
}

fn measured(text: &str, label: &str) -> Vec<Measured> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let mut consequences = Vec::new();
            let mut intergenic = false;
            for part in columns[7].split(';') {
                if part == "PREDICTED_INTERGENIC" {
                    intergenic = true;
                    continue;
                }
                if let Some((key, value)) = part.split_once('=') {
                    if key.starts_with("PREDICTED_") {
                        consequences.push((
                            key.to_string(),
                            value.split(',').map(str::to_string).collect(),
                        ));
                    }
                }
            }
            Measured {
                id: columns[2].to_string(),
                consequences,
                intergenic,
            }
        })
        .collect()
}

fn rendered(annotation: &Annotation) -> Measured {
    Measured {
        id: annotation.id.clone(),
        consequences: annotation.consequences.clone(),
        intergenic: annotation.intergenic == Some(true),
    }
}

/// label, promoter window, maximum breakend length, GTF given, BED given.
fn runs() -> Vec<(&'static str, i32, i32, bool, bool)> {
    vec![
        ("default", 1000, -1, true, true),
        ("wide-promoter", 5000, -1, true, true),
        ("breakend-as-cnv", 1000, 1_000_000, true, true),
        ("gtf-only", 1000, -1, true, false),
        ("bed-only", 1000, -1, false, true),
        ("neither", 1000, -1, false, false),
    ]
}

#[test]
fn every_annotation_matches_the_golden() {
    let text = golden();
    let transcripts = transcripts(&text);
    let non_coding = non_coding(&text);
    let variants = variants(&text);
    let mut compared = 0;
    for (label, window, max_breakend, have_gtf, have_bed) in runs() {
        let produced: Vec<Measured> = variants
            .iter()
            .map(|variant| {
                rendered(
                    &annotate_structural_variant(
                        variant,
                        &transcripts,
                        &non_coding,
                        have_gtf,
                        have_bed,
                        window,
                        max_breakend,
                    )
                    .unwrap_or_else(|error| panic!("{label}/{}: {}", variant.id, error.message())),
                )
            })
            .collect();
        assert_eq!(produced, measured(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 6, "the runs that wrote records");
}

/// The four duplication answers over the same gene, differing only in where the two breakpoints
/// land.
#[test]
fn spanning_and_overlapping_are_different_questions() {
    let text = golden();
    let out = measured(&text, "default");
    let of = |id: &str| {
        out.iter()
            .find(|record| record.id == id)
            .expect(id)
            .consequences
            .iter()
            .find(|(key, _)| {
                key.starts_with("PREDICTED_")
                    && key != "PREDICTED_NONCODING_SPAN"
                    && key != "PREDICTED_NONCODING_BREAKPOINT"
            })
            .expect("a coding consequence")
            .0
            .clone()
    };
    assert_eq!(of("dup-span"), "PREDICTED_COPY_GAIN");
    assert_eq!(of("dup-int-exon"), "PREDICTED_INTRAGENIC_EXON_DUP");
    assert_eq!(of("dup-partial-exon"), "PREDICTED_PARTIAL_EXON_DUP");
    assert_eq!(of("dup-lof"), "PREDICTED_LOF");
    assert_eq!(of("dup-partial"), "PREDICTED_DUP_PARTIAL");

    // And the count that separates them: a variant CONTAINING a feature has no breakpoint inside
    // it, however much they overlap.
    let feature = Interval::new("chr1", 10300, 10400).expect("a feature");
    assert_eq!(
        count_breakends_inside_feature(&Interval::new("chr1", 10250, 10650).unwrap(), &feature),
        0,
        "spanned, so none inside"
    );
    assert_eq!(
        count_breakends_inside_feature(&Interval::new("chr1", 10350, 10650).unwrap(), &feature),
        1
    );
    assert_eq!(
        count_breakends_inside_feature(&Interval::new("chr1", 10320, 10380).unwrap(), &feature),
        2
    );
}

/// The same interval gives a duplication one answer and a CNV another, and an intronic one is left
/// alone because it is not in the reclassified set.
#[test]
fn a_cnv_is_annotated_as_a_duplication_and_then_reclassified() {
    let text = golden();
    let transcripts = transcripts(&text);
    let alpha = transcripts
        .iter()
        .find(|t| t.gene_name == "ALPHA")
        .expect("ALPHA");
    let exonic = Interval::new("chr1", 10350, 10450).expect("an interval");
    let intronic = Interval::new("chr1", 10420, 10480).expect("an interval");
    assert_eq!(
        annotate_duplication(&exonic, alpha, false),
        "PREDICTED_PARTIAL_EXON_DUP"
    );
    assert_eq!(
        annotate_copy_number_variant(&exonic, alpha),
        "PREDICTED_MSV_EXON_OVERLAP"
    );
    assert_eq!(
        annotate_duplication(&intronic, alpha, false),
        "PREDICTED_INTRONIC"
    );
    assert_eq!(
        annotate_copy_number_variant(&intronic, alpha),
        "PREDICTED_INTRONIC",
        "not in the reclassified set"
    );
}

/// The GTF codec puts the smaller coordinate in `start` whatever the strand, so the start site of a
/// minus-strand transcript is its END.
#[test]
fn the_start_site_is_the_end_on_the_minus_strand() {
    let text = golden();
    let transcripts = transcripts(&text);
    let beta = transcripts
        .iter()
        .find(|t| t.gene_name == "BETA")
        .expect("BETA");
    assert!(beta.negative_strand);
    assert_eq!(beta.transcription_start_site(), beta.end);
    assert_eq!(beta.transcription_start_site(), 21000);
    let alpha = transcripts
        .iter()
        .find(|t| t.gene_name == "ALPHA")
        .expect("ALPHA");
    assert_eq!(alpha.transcription_start_site(), alpha.start);

    // Which the golden shows as LOF at BETA's larger coordinate and only UTR at its smaller one.
    let out = measured(&text, "default");
    let keys = |id: &str| {
        out.iter()
            .find(|record| record.id == id)
            .expect(id)
            .consequences
            .iter()
            .map(|(key, _)| key.clone())
            .collect::<Vec<String>>()
    };
    assert_eq!(keys("del-tss-minus"), vec!["PREDICTED_LOF"]);
    assert_eq!(keys("del-end-minus"), vec!["PREDICTED_UTR"]);
}

/// It uses the SIMPLE rule rather than the deletion one, so it never picks up the start-site LOF,
/// and the LOF it does find is downgraded. It also produces two segments.
#[test]
fn a_breakend_downgrades_its_lof_and_produces_two_segments() {
    let text = golden();
    let transcripts = transcripts(&text);
    let alpha = transcripts
        .iter()
        .find(|t| t.gene_name == "ALPHA")
        .expect("ALPHA");
    let coding = Interval::new("chr1", 10350, 10350).expect("an interval");
    assert_eq!(
        annotate_breakend(&coding, alpha),
        "PREDICTED_BREAKEND_EXONIC"
    );

    // Two segments, which is how one record carries two consequences for one gene.
    let bnd = variants(&text)
        .into_iter()
        .find(|variant| variant.id == "bnd-dup")
        .expect("bnd-dup");
    assert_eq!(sv_segments(&bnd, -1).expect("segments").len(), 2);
    let out = measured(&text, "default");
    let record = out.iter().find(|r| r.id == "bnd-dup").expect("bnd-dup");
    assert_eq!(
        record
            .consequences
            .iter()
            .map(|(key, _)| key.clone())
            .collect::<Vec<String>>(),
        vec!["PREDICTED_BREAKEND_EXONIC", "PREDICTED_INTRONIC"]
    );

    // And with a maximum length given, the same record is read as a duplication instead.
    assert_eq!(annotation_type_for_breakend(&bnd, 1_000_000), SvType::Dup);
    assert_eq!(annotation_type_for_breakend(&bnd, -1), SvType::Bnd);
}

/// It is the window upstream of the start site, and it is only added for a gene with no coding
/// consequence. A promoter is also intergenic.
#[test]
fn the_promoter_is_inferred_from_a_window() {
    let text = golden();
    let transcripts = transcripts(&text);
    let alpha = transcripts
        .iter()
        .find(|t| t.gene_name == "ALPHA")
        .expect("ALPHA");
    assert_eq!(
        promoter_interval(alpha, 1000).expect("a promoter"),
        Interval::new("chr1", 9000, 9999).expect("an interval")
    );
    let beta = transcripts
        .iter()
        .find(|t| t.gene_name == "BETA")
        .expect("BETA");
    assert_eq!(
        promoter_interval(beta, 1000).expect("a promoter"),
        Interval::new("chr1", 21001, 22000).expect("an interval"),
        "upstream is the other way on the minus strand"
    );

    // A promoter overlap is also intergenic, and the window's width is what moves `far-upstream`.
    let default = measured(&text, "default");
    let wide = measured(&text, "wide-promoter");
    let record = |out: &[Measured], id: &str| {
        out.iter()
            .find(|record| record.id == id)
            .expect(id)
            .consequences
            .iter()
            .map(|(key, _)| key.clone())
            .collect::<Vec<String>>()
    };
    assert!(
        default
            .iter()
            .find(|r| r.id == "promoter")
            .expect("promoter")
            .intergenic
    );
    assert_eq!(record(&default, "promoter"), vec!["PREDICTED_PROMOTER"]);
    assert_eq!(
        record(&default, "far-upstream"),
        vec!["PREDICTED_NEAREST_TSS"]
    );
    assert_eq!(record(&wide, "far-upstream"), vec!["PREDICTED_PROMOTER"]);
}

/// Both declare themselves valid and then build an interval that ends one base before it starts.
#[test]
fn two_arguments_are_accepted_and_then_crash() {
    let text = golden();

    // --promoter-window-length 0, which passes its own minValue of 0.
    let transcripts = transcripts(&text);
    let alpha = transcripts
        .iter()
        .find(|t| t.gene_name == "ALPHA")
        .expect("ALPHA");
    let (class, message) = refusal(&text, "zero-promoter");
    assert_eq!(class, "java.lang.IllegalArgumentException");
    let produced = promoter_interval(alpha, 0).expect_err("an empty window");
    assert_eq!(produced.message(), message);
    assert_eq!(
        produced,
        AnnotateError::InvalidInterval {
            contig: "chr1".to_string(),
            start: 10000,
            end: 9999
        }
    );
    // A window of one is the smallest that works, which is what makes zero the boundary.
    assert!(promoter_interval(alpha, 1).is_ok());

    // A breakend carrying the conventional SVLEN of -1, which is <= any maximum.
    let (class, message) = refusal(&text, "bnd-no-length");
    assert_eq!(class, "java.lang.IllegalArgumentException");
    let bnd = Variant {
        id: "bnd-no-length".to_string(),
        contig: "chr1".to_string(),
        position: 10350,
        end: 10350,
        sv_type: SvType::Bnd,
        sv_length: -1,
        contig2: Some("chr1".to_string()),
        end2: Some(15000),
        strands: Some("+-".to_string()),
        complex_type: None,
        complex_intervals: Vec::new(),
    };
    assert_eq!(
        annotation_type_for_breakend(&bnd, 1_000_000),
        SvType::Del,
        "a length of -1 is under the maximum"
    );
    assert_eq!(
        sv_segments(&bnd, 1_000_000)
            .expect_err("an empty interval")
            .message(),
        message
    );
    // Without a maximum the same record is fine, which is what makes the argument the cause.
    assert!(sv_segments(&bnd, -1).is_ok());
}

/// A complex variant needs both of its fields and a translocation needs its second contig.
#[test]
fn the_three_missing_field_refusals() {
    let text = golden();
    let base = Variant {
        id: "cpx".to_string(),
        contig: "chr1".to_string(),
        position: 10300,
        end: 10400,
        sv_type: SvType::Cpx,
        sv_length: 101,
        contig2: None,
        end2: None,
        strands: None,
        complex_type: None,
        complex_intervals: Vec::new(),
    };
    for (label, variant, expected) in [
        (
            "cpx-no-intervals",
            Variant {
                complex_type: Some(ComplexSubtype::DelInv),
                ..base.clone()
            },
            AnnotateError::CpxWithoutIntervals,
        ),
        (
            "cpx-no-type",
            Variant {
                complex_intervals: vec!["DEL_chr1:10300-10400".to_string()],
                ..base.clone()
            },
            AnnotateError::CpxWithoutType,
        ),
        (
            "ctx-no-chr2",
            Variant {
                sv_type: SvType::Ctx,
                position: 10450,
                end: 10450,
                ..base.clone()
            },
            AnnotateError::CtxWithoutContig2,
        ),
    ] {
        let (class, message) = refusal(&text, label);
        assert_eq!(
            class, "org.broadinstitute.hellbender.exceptions.UserException",
            "{label}"
        );
        let produced = sv_segments(&variant, -1).expect_err(label);
        assert_eq!(produced, expected, "{label}");
        assert_eq!(produced.message(), message, "{label}");
    }
}

/// The argument documents itself as taking a BED "with header" and a header row is what it cannot
/// read.
#[test]
fn the_non_coding_bed_refuses_its_documented_header() {
    let text = golden();
    let (class, message) = refusal(&text, "bed-header");
    assert_eq!(class, "java.lang.NumberFormatException");
    assert_eq!(
        AnnotateError::NumberFormat {
            text: "start".to_string()
        }
        .message(),
        message
    );
    // The header row's second column is the one that fails, and the file underneath it parses.
    assert_eq!(non_coding(&text).len(), 3);
}

/// Nothing is written at all when neither reference file is given, not even the intergenic flag.
#[test]
fn without_a_gtf_there_is_no_intergenic_flag() {
    let text = golden();
    for record in measured(&text, "neither") {
        assert!(record.consequences.is_empty(), "{}", record.id);
        assert!(!record.intergenic, "{}", record.id);
    }
    for record in measured(&text, "bed-only") {
        assert!(!record.intergenic, "{}", record.id);
        assert!(
            record
                .consequences
                .iter()
                .all(|(key, _)| key.starts_with("PREDICTED_NONCODING_")),
            "{}",
            record.id
        );
    }
}
