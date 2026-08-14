//! Conformance for `ApplyVQSR`'s per-record filtering against GATK 4.6.2.0, compared as the FILTER
//! and INFO of every written record and as all three refusals.
//!
//! Golden from `tools/readfilter-conformance/ApplyVqsrSiteFilteringDump.java`.
//!
//! # What this suite is for
//!
//!  * **the bands are the intervals the header lines describe**, the last tranche meaning `PASS`
//!    and a LOD on a boundary belonging to the tranche below it;
//!  * **the recal record must agree on both ends**, one of the two candidates at `chr1:800` being
//!    skipped for its end alone;
//!  * **every negative VQSLOD is written in scientific notation**;
//!  * **a record emitted untouched survives `--exclude-filtered`**;
//!  * **and the three refusals are told apart only by their cause**.
//!
//! The tranches file is read with [`gatk_engine::tranches`], measured by the `apply-vqsr-tranches`
//! suite.

use gatk_corpus as corpus;
use gatk_engine::tranches::read_tranches;
use gatk_tools::apply_vqsr::{
    filter_string, filter_string_by_cutoff, keep, site_specific_filtering, writes_out, Annotation,
    RecalRecord, SiteFilteringError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/apply_vqsr_site_filtering.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn input(text: &str, label: &str) -> String {
    unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    )
}

/// One record of either file, decoded as far as this slice looks at it.
struct Record {
    start: i32,
    end: i32,
    reference: String,
    alternates: Vec<String>,
    filters: Vec<String>,
    info: Vec<(String, Option<String>)>,
    /// The whole line, which is what the output is compared against.
    line: String,
}

impl RecalRecord for Record {
    fn start(&self) -> i32 {
        self.start
    }
    fn end(&self) -> i32 {
        self.end
    }
    fn lod_string(&self) -> Option<String> {
        self.attribute("VQSLOD")
    }
    fn culprit(&self) -> Option<String> {
        self.attribute("culprit")
    }
    fn has_positive_label(&self) -> bool {
        self.has("POSITIVE_TRAIN_SITE")
    }
    fn has_negative_label(&self) -> bool {
        self.has("NEGATIVE_TRAIN_SITE")
    }
}

impl Record {
    fn attribute(&self, key: &str) -> Option<String> {
        self.info
            .iter()
            .find(|(name, _)| name == key)
            .and_then(|(_, value)| value.clone())
    }

    fn has(&self, key: &str) -> bool {
        self.info.iter().any(|(name, _)| name == key)
    }

    /// `isSNP()` for the biallelic records this fixture carries.
    fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternates.iter().all(|alt| alt.len() == 1)
    }
}

fn records(whole: &str) -> Vec<Record> {
    whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let start: i32 = field[1].parse().expect("a position");
            let reference = field[3].to_string();
            Record {
                start,
                // htsjdk's stop for a record with no END attribute.
                end: start + reference.len() as i32 - 1,
                reference,
                alternates: field[4].split(',').map(|alt| alt.to_string()).collect(),
                filters: match field[6] {
                    "." | "PASS" => Vec::new(),
                    list => list.split(';').map(|name| name.to_string()).collect(),
                },
                info: match field[7] {
                    "." => Vec::new(),
                    list => list
                        .split(';')
                        .map(|entry| match entry.split_once('=') {
                            Some((key, value)) => (key.to_string(), Some(value.to_string())),
                            None => (entry.to_string(), None),
                        })
                        .collect(),
                },
                line: line.to_string(),
            }
        })
        .collect()
}

/// The record lines of one run's output, in file order.
fn written(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

/// The cause of one run's refusal, which is where the tool's own wording is.
fn cause(text: &str, run: &str) -> (String, String) {
    let row = rows(text, "cause")
        .into_iter()
        .find(|row| row[0] == run)
        .unwrap_or_else(|| panic!("no cause {run}"));
    let (class, message) = row[1].split_once(':').expect("class and message");
    (class.to_string(), unescape(message))
}

/// How a run was configured: which files, which mode, and the four flags.
struct Run {
    label: &'static str,
    variants: &'static str,
    recal: &'static str,
    snp_mode: bool,
    level: Option<f64>,
    cutoff: Option<f64>,
    ignore_all_filters: bool,
    ignored_filters: &'static [&'static str],
    exclude_filtered: bool,
}

const RUNS: [Run; 7] = [
    Run {
        label: "snp-mode",
        variants: "variants",
        recal: "recal",
        snp_mode: true,
        level: Some(0.0),
        cutoff: None,
        ignore_all_filters: false,
        ignored_filters: &[],
        exclude_filtered: false,
    },
    Run {
        label: "indel-mode",
        variants: "variants",
        recal: "recal",
        snp_mode: false,
        level: Some(0.0),
        cutoff: None,
        ignore_all_filters: false,
        ignored_filters: &[],
        exclude_filtered: false,
    },
    Run {
        label: "ignore-all-filters",
        variants: "variants",
        recal: "recal",
        snp_mode: true,
        level: Some(0.0),
        cutoff: None,
        ignore_all_filters: true,
        ignored_filters: &[],
        exclude_filtered: false,
    },
    Run {
        label: "ignore-named-filter",
        variants: "variants",
        recal: "recal",
        snp_mode: true,
        level: Some(0.0),
        cutoff: None,
        ignore_all_filters: false,
        ignored_filters: &["weak"],
        exclude_filtered: false,
    },
    Run {
        label: "exclude-filtered",
        variants: "variants",
        recal: "recal",
        snp_mode: true,
        level: Some(0.0),
        cutoff: None,
        ignore_all_filters: false,
        ignored_filters: &[],
        exclude_filtered: true,
    },
    Run {
        label: "lod-cutoff",
        variants: "variants",
        recal: "recal",
        snp_mode: true,
        level: None,
        cutoff: Some(1.0),
        ignore_all_filters: false,
        ignored_filters: &[],
        exclude_filtered: false,
    },
    Run {
        label: "ends",
        variants: "ends",
        recal: "ends-recal",
        snp_mode: true,
        level: Some(0.0),
        cutoff: None,
        ignore_all_filters: false,
        ignored_filters: &[],
        exclude_filtered: false,
    },
];

/// The FILTER and INFO the tool gives each record, as the output line the writer would produce.
///
/// The record's other columns are copied from the input, since `VariantContextBuilder(vc)` keeps
/// them: only the seventh and eighth are this slice's.
fn ours(text: &str, run: &Run) -> Vec<String> {
    let variants = records(&input(text, run.variants));
    let recals = records(&input(text, run.recal));
    let tranches = read_tranches("tranches", &input(text, "tranches")).expect("a good file");
    let ignored: Vec<String> = run
        .ignored_filters
        .iter()
        .map(|filter| filter.to_string())
        .collect();

    let mut lines = Vec::new();
    for variant in &variants {
        let of_this_mode = variant.is_snp() == run.snp_mode;
        let recalibrated = gatk_tools::apply_vqsr::recalibrates(
            of_this_mode,
            &variant.filters,
            run.ignore_all_filters,
            &ignored,
        );
        if !recalibrated {
            // Emitted exactly as it came in, and `--exclude-filtered` never reaches it.
            assert!(writes_out(false, "", run.exclude_filtered));
            lines.push(variant.line.clone());
            continue;
        }
        let (annotation, lod) =
            site_specific_filtering(variant.start, variant.end, &recals, "[VC]")
                .expect("every recalibrated record here has a partner");
        let filter = match (run.level, run.cutoff) {
            (Some(level), _) => filter_string(&keep(&tranches, level), lod),
            (None, Some(cutoff)) => filter_string_by_cutoff(lod, cutoff),
            (None, None) => unreachable!("one of the two"),
        };
        if !writes_out(true, &filter, run.exclude_filtered) {
            continue;
        }
        lines.push(rendered(variant, &filter, &annotation));
    }
    lines
}

/// The output line: the input's columns with the tool's FILTER and its INFO, which htsjdk writes in
/// byte order of the keys.
fn rendered(variant: &Record, filter: &str, annotation: &Annotation) -> String {
    let mut info: Vec<String> = Vec::new();
    if annotation.negative_label {
        info.push("NEGATIVE_TRAIN_SITE".to_string());
    }
    if annotation.positive_label {
        info.push("POSITIVE_TRAIN_SITE".to_string());
    }
    info.push(format!("VQSLOD={}", annotation.vqslod));
    info.push(format!("culprit={}", annotation.culprit));
    info.sort();
    let mut field: Vec<String> = variant.line.split('\t').map(|f| f.to_string()).collect();
    field[6] = filter.to_string();
    field[7] = info.join(";");
    field.join("\t")
}

#[test]
fn every_written_record_matches_the_golden_byte_for_byte() {
    let text = golden();
    for run in &RUNS {
        assert_eq!(ours(&text, run), written(&text, run.label), "{}", run.label);
    }
}

#[test]
fn the_bands_are_the_intervals_the_header_lines_describe() {
    let text = golden();
    let filters: Vec<String> = written(&text, "snp-mode")
        .into_iter()
        .map(|line| line.split('\t').nth(6).expect("a FILTER").to_string())
        .collect();
    // 5.0, 2.0, 0.0, -3.0, an indel, an already-filtered record, and 1.5 on the boundary.
    assert_eq!(
        filters,
        vec![
            "PASS",
            "VQSRTrancheSNP90.00to99.00",
            "VQSRTrancheSNP99.00to100.00",
            "VQSRTrancheSNP99.00to100.00+",
            ".",
            "weak",
            "VQSRTrancheSNP90.00to99.00"
        ]
    );
}

#[test]
fn every_negative_lod_is_written_in_scientific_notation() {
    let text = golden();
    let line = written(&text, "snp-mode")
        .into_iter()
        .find(|line| line.starts_with("chr1\t400\t"))
        .expect("the negative one");
    assert!(line.contains("VQSLOD=-3.000e+00"), "{line}");
    // And the culprit the recal record did not carry.
    assert!(line.contains("culprit=."), "{line}");
}

#[test]
fn the_recal_record_must_agree_on_both_ends() {
    let text = golden();
    let line = written(&text, "ends").first().expect("one record").clone();
    // Two candidates start at 800 and the first does not end there.
    assert!(line.contains("culprit=TAKEN"), "{line}");
    let recals = records(&input(&text, "ends-recal"));
    assert_eq!(recals.len(), 3);
    assert_eq!(
        site_specific_filtering(850, 850, &recals, "[VC]").expect_err("nothing ends at 850"),
        SiteFilteringError::NoRecalRecord {
            record: "[VC]".to_string()
        }
    );
}

#[test]
fn a_record_emitted_untouched_survives_exclude_filtered() {
    let text = golden();
    let lines = written(&text, "exclude-filtered");
    // Everything the tool filtered is gone, and the two it never touched are not.
    assert_eq!(lines.len(), 3);
    assert!(lines.iter().any(|line| line.contains("\tweak\t")));
    assert!(lines.iter().any(|line| line.starts_with("chr1\t500\t")));
    assert!(lines.iter().any(|line| line.contains("\tPASS\t")));
}

#[test]
fn the_three_refusals_are_told_apart_only_by_their_cause() {
    let text = golden();
    for (run, error) in [
        (
            "ends-mismatch",
            SiteFilteringError::NoRecalRecord {
                record: String::new(),
            },
        ),
        (
            "no-lod",
            SiteFilteringError::NoLod {
                record: String::new(),
            },
        ),
        (
            "bad-lod",
            SiteFilteringError::UnreadableLod {
                record: String::new(),
            },
        ),
    ] {
        let (class, message) = cause(&text, run);
        assert_eq!(error.class(), class, "{run}");
        // The message ends with the record's own toString, which is not ported: what is compared is
        // everything the tool wrote in front of it.
        assert!(message.starts_with(&error.message()), "{run}: {message}");
    }
    // And the wrapper is the same class and shape for all three.
    let rows = rows(&text, "error");
    for run in ["ends-mismatch", "no-lod", "bad-lod"] {
        let row = rows.iter().find(|row| row[0] == run).expect("an error row");
        assert!(row[1].starts_with(
            "org.broadinstitute.hellbender.exceptions.GATKException:Exception thrown at chr1:"
        ));
    }
}
