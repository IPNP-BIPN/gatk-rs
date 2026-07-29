//! Conformance for the `IntervalWalker` traversal against GATK 4.6.2.0.
//!
//! The golden was produced by a real `IntervalWalker` subclass run through the real command line
//! over twenty-four argument combinations, so what is replayed here is the interval list a tool
//! gets out of `-L`, `-XL`, `--interval-padding`, `--interval-exclusion-padding`,
//! `--interval-set-rule` and `--interval-merging-rule`, along with the reads and reference each
//! interval arrived with.
//!
//! The arguments live in a table here rather than in the golden: a label is a configuration, and
//! the row carries the *result*, which is the whole point of comparing them.

use gatk_corpus as corpus;
use gatk_engine::interval::MergingRule;
use gatk_engine::interval_args::{IntervalArgumentError, IntervalArguments, SetRule};
use gatk_engine::reads::ReadsDataSource;
use gatk_engine::reference::ReferenceFileSource;
use gatk_readfilter::with_header;
use gatk_tools::interval_walker::{self, TraversalError};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/interval_walker.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => out.push(other),
            None => {}
        }
    }
    out
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .unwrap_or_else(|| panic!("the golden carries no {kind} row"))
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|v| v.to_string()).collect()
}

/// The arguments each labelled run was given, and whether it had a reference.
fn configuration(label: &str) -> (IntervalArguments, bool) {
    let base = |include: &[&str]| IntervalArguments {
        include: strings(include),
        ..IntervalArguments::default()
    };
    let arguments = match label {
        "chr1" => base(&["chr1"]),
        "chr1:100-160" => base(&["chr1:100-160"]),
        "chr2-then-chr1" => base(&["chr2", "chr1"]),
        "abutting-all" => base(&["chr1:1-100", "chr1:101-200"]),
        "abutting-overlapping-only" => IntervalArguments {
            merging_rule: MergingRule::OverlappingOnly,
            ..base(&["chr1:1-100", "chr1:101-200"])
        },
        "overlapping-overlapping-only" => IntervalArguments {
            merging_rule: MergingRule::OverlappingOnly,
            ..base(&["chr1:1-100", "chr1:50-150"])
        },
        "one-base-gap" => base(&["chr1:1-10", "chr1:12-20"]),
        "padded" => IntervalArguments {
            padding: 20,
            ..base(&["chr1:50-60"])
        },
        "padded-clamped-start" => IntervalArguments {
            padding: 20,
            ..base(&["chr1:1-5"])
        },
        "padded-clamped-end" => IntervalArguments {
            padding: 20,
            ..base(&["chr1:195-200"])
        },
        "padded-overlapping-only" => IntervalArguments {
            padding: 5,
            merging_rule: MergingRule::OverlappingOnly,
            ..base(&["chr1:1-50", "chr1:60-100"])
        },
        "intersection" => IntervalArguments {
            set_rule: SetRule::Intersection,
            ..base(&["chr1:1-100", "chr1:50-150"])
        },
        "intersection-single" => IntervalArguments {
            set_rule: SetRule::Intersection,
            ..base(&["chr1:1-100"])
        },
        "intersection-empty" => IntervalArguments {
            set_rule: SetRule::Intersection,
            ..base(&["chr1", "chr2"])
        },
        "exclude-only" => IntervalArguments {
            exclude: strings(&["chr1"]),
            ..IntervalArguments::default()
        },
        "exclude-middle" => IntervalArguments {
            exclude: strings(&["chr1:50-100"]),
            ..base(&["chr1"])
        },
        "exclude-head" => IntervalArguments {
            exclude: strings(&["chr1:1-50"]),
            ..base(&["chr1"])
        },
        "exclude-tail" => IntervalArguments {
            exclude: strings(&["chr1:150-200"]),
            ..base(&["chr1"])
        },
        "exclude-everything" => IntervalArguments {
            exclude: strings(&["chr1"]),
            ..base(&["chr1"])
        },
        "exclude-other-contig" => IntervalArguments {
            exclude: strings(&["chr2"]),
            ..base(&["chr1", "chr2"])
        },
        "exclude-padded" => IntervalArguments {
            exclude: strings(&["chr1:40-60"]),
            exclusion_padding: 10,
            ..base(&["chr1:1-100"])
        },
        "unmapped" => base(&["unmapped"]),
        "unmapped-and-chr1" => base(&["unmapped", "chr1:1-20"]),
        "chr1-noref" => base(&["chr1:100-160"]),
        other => panic!("{other} is in the golden but not configured here"),
    };
    (arguments, label != "chr1-noref")
}

#[test]
fn every_traversal_hands_apply_what_the_reference_hands_it() {
    let text = golden();

    let dir = std::env::temp_dir().join(format!("gatk-rs-intervalwalker-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bam = dir.join("reads.bam");
    let bai = dir.join("reads.bai");
    let fasta = dir.join("ref.fasta");
    std::fs::write(&bam, corpus::decode_base64(field(&text, "bam"))).unwrap();
    std::fs::write(&bai, corpus::decode_base64(field(&text, "bai"))).unwrap();
    std::fs::write(&fasta, unescape(field(&text, "fasta"))).unwrap();
    std::fs::write(dir.join("ref.fasta.fai"), unescape(field(&text, "fai"))).unwrap();

    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();

    // The labels in the order the harness ran them, with their apply rows and their outcome.
    let mut labels: Vec<String> = Vec::new();
    let mut applied_rows: std::collections::HashMap<String, Vec<String>> = Default::default();
    let mut outcomes: std::collections::HashMap<String, String> = Default::default();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("apply\t") {
            let mut parts = rest.splitn(3, '\t');
            let label = parts.next().expect("label").to_string();
            let _index = parts.next();
            let payload = parts.next().unwrap_or("").to_string();
            applied_rows.entry(label).or_default().push(payload);
        } else if let Some(rest) = line.strip_prefix("summary\t") {
            let (label, outcome) = rest.split_once('\t').expect("a summary row has an outcome");
            labels.push(label.to_string());
            outcomes.insert(label.to_string(), outcome.to_string());
        }
    }
    assert!(!labels.is_empty(), "the golden carries no summary rows");

    let mut compared = 0;
    let mut failures = 0;
    for label in &labels {
        let (arguments, has_reference) = configuration(label);
        let expected_rows = applied_rows.get(label).cloned().unwrap_or_default();
        let outcome = &outcomes[label];

        let header_for_filter = header.clone();
        let filter: Box<dyn Fn(&BamRecord) -> bool> =
            Box::new(move |read: &BamRecord| with_header::wellformed(read, &header_for_filter));

        let mut reference = if has_reference {
            Some(ReferenceFileSource::open(&fasta).expect("the fixture reference opens"))
        } else {
            None
        };
        let result =
            interval_walker::traverse(&source, reference.as_mut(), &arguments, filter.as_ref());

        // The reference reports the exception's class, and Barclay collapses two different causes
        // into `BadArgumentValue`: an empty intersection and an exclusion that removed everything.
        // Both are listed, because a port that refused for the other reason would still be wrong.
        let allowed: Vec<TraversalError> = match outcome.strip_prefix("E:") {
            None => Vec::new(),
            Some("org.broadinstitute.barclay.argparser.CommandLineException$MissingArgument") => {
                vec![TraversalError::MissingIntervalArgument]
            }
            Some("org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue") => {
                vec![
                    TraversalError::Intervals(IntervalArgumentError::EmptyIntersection),
                    TraversalError::Intervals(IntervalArgumentError::ExclusionRemovedEverything),
                ]
            }
            Some(other) => panic!("{label}: the golden carries an unhandled exception {other}"),
        };

        match (&result, outcome.as_str()) {
            (Ok(_), "ok") => {}
            (Err(error), outcome) if outcome.starts_with("E:") => {
                assert!(
                    allowed.contains(error),
                    "{label}: the reference raised {outcome}, the port raised {error:?}"
                );
                failures += 1;
                continue;
            }
            (Ok(_), outcome) => {
                panic!("{label}: the reference failed with {outcome}, the port did not")
            }
            (Err(error), _) => {
                panic!("{label}: the port failed with {error:?}, the reference did not")
            }
        }

        let ours: Vec<String> = result
            .expect("checked above")
            .into_iter()
            .map(|mut applied| {
                let window = match applied.context.window() {
                    None => "null".to_string(),
                    Some(window) => format!("{}:{}-{}", window.contig, window.start, window.end),
                };
                let bases = match reference.as_mut() {
                    Some(source) => applied.context.bases(source).expect("bases"),
                    None => Vec::new(),
                };
                format!(
                    "{}:{}-{}|{}|{}|{}",
                    applied.interval.contig,
                    applied.interval.start,
                    applied.interval.end,
                    applied.reads.len(),
                    window,
                    String::from_utf8(bases).expect("ASCII bases"),
                )
            })
            .collect();

        assert_eq!(&ours, &expected_rows, "{label}");
        compared += ours.len();
    }

    std::fs::remove_dir_all(&dir).ok();
    println!(
        "{compared} apply calls over {} traversals, {failures} refused exactly where the \
         reference refused",
        labels.len()
    );
}
