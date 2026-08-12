//! Conformance for `ApplyBQSR` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/ApplyBqsrDump.java`. The output BAMs and their indexes
//! travel in full, base64, as every record-transform tool's do, and so do the two input fixtures,
//! their indexes and the two recalibration tables: a test that rebuilt any of them would be
//! inventing part of the input.
//!
//! # What this suite is for
//!
//! The eighth whole tool of the largest archetype, and the one every BQSR brick was built for. Its
//! own body is four lines, so what this pins is where the transformer is hooked and what the
//! traversal hands over:
//!
//!  * **the transformer runs after the read filters**, so a read the filters drop is never
//!    recalibrated and never written. The `filtered-read` and `filters-disabled` runs are the same
//!    file with and without `WellformedReadFilter`;
//!  * **the recalibration file decides the read group keys**, so a BAM whose group the report does
//!    not name is refused, and `--allow-missing-read-group` quantizes those reads without
//!    recalibrating them;
//!  * **the quantization arguments are two separate mechanisms** that both reach the output.

use gatk_corpus as corpus;
use gatk_engine::bqsr_transformer::ApplyBqsrArguments;
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::with_header;
use gatk_tools::apply_bqsr::{self, ApplyBqsrError};
use gatk_tools::sam_output::Options;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/apply_bqsr.txt.gz"),
    )
}

/// The escaping the harness applies to text that carries newlines.
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
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// Rows of one kind that carry `<label>\t<value>`.
fn pairs<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str)> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .filter_map(|rest| rest.split_once('\t'))
        .collect()
}

fn value<'a>(text: &'a str, kind: &str, label: &str) -> Option<&'a str> {
    pairs(text, kind)
        .into_iter()
        .find(|(name, _)| *name == label)
        .map(|(_, value)| value)
}

/// The arguments and the fixture each labelled run was given.
struct Configuration {
    fixture: &'static str,
    recal: &'static str,
    arguments: ApplyBqsrArguments,
    /// Whether the default `WellformedReadFilter` was left in place.
    wellformed: bool,
}

fn configuration(label: &str) -> Configuration {
    let base = Configuration {
        fixture: "input",
        recal: "one-group",
        arguments: ApplyBqsrArguments::default(),
        wellformed: true,
    };
    match label {
        "plain" => base,
        "emit-original-quals" => Configuration {
            arguments: ApplyBqsrArguments {
                emit_original_quals: true,
                ..ApplyBqsrArguments::default()
            },
            ..base
        },
        "quantize-4" => Configuration {
            arguments: ApplyBqsrArguments {
                quantization_levels: 4,
                ..ApplyBqsrArguments::default()
            },
            ..base
        },
        "no-quantization" => base,
        "static-quals" => Configuration {
            arguments: ApplyBqsrArguments {
                static_quantization_quals: vec![10, 30],
                ..ApplyBqsrArguments::default()
            },
            ..base
        },
        "static-quals-round-down" => Configuration {
            arguments: ApplyBqsrArguments {
                static_quantization_quals: vec![10, 30],
                round_down: true,
                ..ApplyBqsrArguments::default()
            },
            ..base
        },
        "preserve-25" => Configuration {
            arguments: ApplyBqsrArguments {
                preserve_qscores_less_than: 25,
                ..ApplyBqsrArguments::default()
            },
            ..base
        },
        "global-prior" => Configuration {
            arguments: ApplyBqsrArguments {
                global_qscore_prior: 20.0,
                ..ApplyBqsrArguments::default()
            },
            ..base
        },
        "missing-read-group" => Configuration {
            recal: "elsewhere",
            ..base
        },
        "missing-read-group-allowed" => Configuration {
            recal: "elsewhere",
            arguments: ApplyBqsrArguments {
                allow_missing_read_groups: true,
                ..ApplyBqsrArguments::default()
            },
            ..base
        },
        "filtered-read" => Configuration {
            fixture: "filtered",
            ..base
        },
        "filters-disabled" => Configuration {
            fixture: "filtered",
            wellformed: false,
            ..base
        },
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The two fixtures and their indexes, written out so the port's reader opens the same bytes.
fn install_fixtures(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    for (label, encoded) in pairs(text, "fixture") {
        std::fs::write(
            dir.join(format!("{label}.bam")),
            corpus::decode_base64(encoded),
        )
        .expect("the fixture bam");
    }
    for (label, encoded) in pairs(text, "fixtureindex") {
        std::fs::write(
            dir.join(format!("{label}.bai")),
            corpus::decode_base64(encoded),
        )
        .expect("the fixture index");
    }
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-applybqsr-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let outputs = pairs(&text, "output");
    let errors = pairs(&text, "error");
    assert_eq!(outputs.len(), 11, "eleven runs finished");
    assert_eq!(errors.len(), 1, "and one was refused");

    let mut compared = 0usize;
    for (label, expected_base64) in &outputs {
        let config = configuration(label);
        let bam = dir.join(format!("{}.bam", config.fixture));
        let bai = dir.join(format!("{}.bai", config.fixture));
        let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
        let header = source.header().clone();

        // `WellformedReadFilter` is this tool's default, taken from `GATKTool` and not overridden.
        let header_for_filter = header.clone();
        let filter: Box<dyn Fn(&BamRecord) -> bool> = if config.wellformed {
            Box::new(move |read: &BamRecord| with_header::wellformed(read, &header_for_filter))
        } else {
            Box::new(|_: &BamRecord| true)
        };

        let options = Options {
            command_line: value(&text, "commandline", label)
                .expect("the golden carries the command line"),
            ..Options::default()
        };
        let recal = unescape(
            value(&text, "recal", config.recal).expect("the golden carries the recal table"),
        );

        let (out, index) =
            apply_bqsr::apply_bqsr(&source, &recal, &config.arguments, &options, &filter)
                .unwrap_or_else(|error| panic!("{label}: {}", error.message()));

        let expected = corpus::decode_base64(expected_base64);
        assert_eq!(out.len(), expected.len(), "{label}: output length differs");
        if out != expected {
            let at = out
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            panic!("{label}: first byte difference at offset {at}");
        }
        let expected_index = value(&text, "index", label).expect("an index row");
        match index {
            Some(index) => assert_eq!(
                index,
                corpus::decode_base64(expected_index),
                "{label}: the .bai"
            ),
            None => assert_eq!(expected_index, "absent", "{label}: the index"),
        }
        compared += 1;
    }
    println!("apply-bqsr: {compared} outputs compared byte for byte");
}

/// The one run the reference refuses, and the words it refuses it in.
#[test]
fn a_read_group_the_table_does_not_name_aborts_the_run() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-applybqsr-err-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let (label, message) = pairs(&text, "error")
        .into_iter()
        .map(|(label, rest)| {
            let (_exception, message) = rest.split_once('\t').expect("an error row has a message");
            (label, message)
        })
        .next()
        .expect("the golden lost the run that aborts");
    assert_eq!(label, "missing-read-group");

    let config = configuration(label);
    let source = ReadsDataSource::open(
        &dir.join(format!("{}.bam", config.fixture)),
        &dir.join(format!("{}.bai", config.fixture)),
    )
    .expect("the fixture opens");
    let header = source.header().clone();
    let filter = move |read: &BamRecord| with_header::wellformed(read, &header);
    let recal = unescape(value(&text, "recal", config.recal).unwrap());

    let error = apply_bqsr::apply_bqsr(
        &source,
        &recal,
        &config.arguments,
        &Options::default(),
        &filter,
    )
    .unwrap_err();
    assert!(matches!(error, ApplyBqsrError::Transform(_)), "{error:?}");
    assert_eq!(error.message(), message);
}

/// The filters run before the transformer, so a dropped read is not in the output at all.
#[test]
fn a_filtered_read_is_never_recalibrated_and_never_written() {
    let text = golden();
    // The two runs are the same file; only the filter differs, and the outputs differ in length.
    let filtered = value(&text, "output", "filtered-read").expect("the filtered run");
    let disabled = value(&text, "output", "filters-disabled").expect("the unfiltered run");
    assert_ne!(filtered, disabled);
    assert!(
        corpus::decode_base64(disabled).len() > corpus::decode_base64(filtered).len(),
        "the unfiltered run writes one read more"
    );
}
