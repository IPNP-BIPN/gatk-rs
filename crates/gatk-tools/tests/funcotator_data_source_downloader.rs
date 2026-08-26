//! Conformance for `FuncotatorDataSourceDownloader` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/FuncotatorDataSourceDownloaderDump.java`.
//!
//! # What this suite is for
//!
//!  * **the four bucket paths**, and the version string that carries the reference inside it;
//!  * **the checksum file being cleaned before comparison**, in four spellings that all match;
//!  * **the mismatch message printing an upper-case sum against a lower-case one**, and the
//!    corrupt file staying on disk;
//!  * **the destination with no `-O` being a relative path**;
//!  * **and the startup checks running in an order that misreports which argument was missing**.

use gatk_corpus as corpus;
use gatk_tools::funcotator_data_source_downloader::{
    checksum_path, clean_expected_sha256, data_sources_path, description, expected_sha256,
    is_dest_file_valid, max_version_string, min_version_string, output_location, print_hex_binary,
    startup, validate_integrity, Arguments, DataSourceKind, DownloadError,
};
use std::path::Path;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/funcotator_data_source_downloader.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn constant(text: &str, name: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("constant\t{name}=")))
        .unwrap_or_else(|| panic!("the golden carries constant/{name}"))
        .to_string()
}

fn fixture(text: &str, name: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("fixture\t{name}=")))
        .unwrap_or_else(|| panic!("the golden carries fixture/{name}"))
        .to_string()
}

/// The two paths one bundle is fetched from.
fn paths(text: &str, kind: &str, reference: &str) -> (String, String) {
    let prefix = format!("path\t{kind}\t{reference}\t");
    let line = text
        .lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries path/{kind}/{reference}"));
    let (data, checksum) = line.split_once('\t').expect("two paths");
    (data.to_string(), checksum.to_string())
}

/// The raw contents of one checksum file, as the run was given it.
fn checksum_file(text: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("clean\t{name}\tin=")))
            .unwrap_or_else(|| panic!("the golden carries clean/{name}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    text.lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"))
        .to_string()
}

fn ran(text: &str, label: &str) -> bool {
    text.lines()
        .any(|line| line.starts_with(&format!("run\t{label}\tcode=true")))
}

/// The sha256 of whatever landed at a destination.
fn landed(text: &str, label: &str) -> String {
    text.lines()
        .find_map(|line| {
            line.strip_prefix(&format!("file\t{label}\t"))
                .and_then(|rest| rest.split_once('=').map(|(_, sum)| sum.to_string()))
        })
        .unwrap_or_else(|| panic!("the golden carries file/{label}"))
}

#[test]
fn the_four_bucket_paths_match_the_golden() {
    let text = golden();
    assert_eq!(max_version_string(38), constant(&text, "max-version-38"));
    assert_eq!(max_version_string(19), constant(&text, "max-version-19"));
    assert_eq!(min_version_string(), constant(&text, "min-version"));
    // The current maximum IS the hg38 one, and the minimum carries no reference at all.
    assert_eq!(
        constant(&text, "current-maximum"),
        constant(&text, "max-version-38")
    );
    assert_eq!(
        constant(&text, "current-minimum"),
        constant(&text, "min-version")
    );

    let mut compared = 0;
    for (kind, label) in [
        (DataSourceKind::Somatic, "somatic"),
        (DataSourceKind::Germline, "germline"),
    ] {
        for reference in [38, 19] {
            let (data, checksum) = paths(&text, label, &format!("hg{reference}"));
            assert_eq!(data_sources_path(kind, reference), data);
            assert_eq!(checksum_path(kind, reference), checksum);
            compared += 1;
        }
    }
    assert_eq!(compared, 4);

    assert_eq!(description(DataSourceKind::Somatic, 38), "HG38_Somatic");
    assert_eq!(description(DataSourceKind::Germline, 19), "HG19_Germline");
}

/// The reference number sits INSIDE the version string, between the minor version and the date,
/// and the bundle modifier is a single letter glued to the end of it.
#[test]
fn the_reference_is_part_of_the_version_and_not_of_the_name() {
    let text = golden();
    let (somatic, _) = paths(&text, "somatic", "hg38");
    let (germline, _) = paths(&text, "germline", "hg38");
    assert!(somatic.contains("funcotator_dataSources.v1.8.hg38.20230908s.tar.gz"));
    assert!(germline.contains("funcotator_dataSources.v1.8.hg38.20230908g.tar.gz"));
    // The two differ by exactly one character.
    assert_eq!(somatic.len(), germline.len());
    assert_eq!(
        somatic
            .chars()
            .zip(germline.chars())
            .filter(|(left, right)| left != right)
            .count(),
        1
    );
}

/// Four spellings of the same sum all reduce to it, and the tab test runs on what the space test
/// already shortened.
#[test]
fn the_checksum_file_is_cleaned_before_it_is_compared() {
    let text = golden();
    let expected = fixture(&text, "archive-sha256");
    for name in [
        "plain.sha256",
        "with-name.sha256",
        "tabbed.sha256",
        "upper.sha256",
    ] {
        let contents = checksum_file(&text, name);
        assert_eq!(
            expected_sha256(&contents, "file:///unused").expect("a checksum"),
            expected,
            "{name}"
        );
        assert!(ran(
            &text,
            name.strip_suffix(".sha256").map(run_label).unwrap()
        ));
    }
    // The space is cut first, so a tab after it is never reached.
    assert_eq!(clean_expected_sha256("abc def\tghi"), "abc");
    assert_eq!(clean_expected_sha256("abc\tdef ghi"), "abc");
    assert_eq!(clean_expected_sha256("  ABC  \n"), "abc");
}

/// The run label a checksum file was used under.
fn run_label(name: &str) -> &'static str {
    match name {
        "plain" => "validate-plain",
        "with-name" => "validate-with-name",
        "tabbed" => "validate-tabbed",
        "upper" => "validate-upper",
        other => panic!("no run for {other}"),
    }
}

/// A file with no first line at all is a refusal naming the path as a URI.
#[test]
fn an_empty_checksum_file_is_a_refusal() {
    let text = golden();
    assert_eq!(checksum_file(&text, "empty.sha256"), "");
    let error = expected_sha256("", "file://<dir>/empty.sha256").expect_err("no checksum");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "validate-empty")
    );
}

/// The computed sum is upper-case hex, the expected one was lower-cased, the comparison ignores
/// the difference and the MESSAGE does not.
#[test]
fn a_mismatch_names_an_upper_case_sum_against_a_lower_case_one() {
    let text = golden();
    let computed = fixture(&text, "archive-sha256").to_uppercase();
    let wrong = clean_expected_sha256(&checksum_file(&text, "wrong.sha256"));
    assert!(is_dest_file_valid(
        &computed,
        &fixture(&text, "archive-sha256")
    ));

    let error = validate_integrity(&computed, &wrong).expect_err("a mismatch");
    assert_eq!(
        error,
        DownloadError::Corrupt {
            checksum: computed,
            expected: wrong
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "validate-wrong")
    );

    // And the file that failed its checksum is still there, byte for byte the archive.
    assert_eq!(
        landed(&text, "validate-wrong"),
        fixture(&text, "archive-sha256")
    );
}

/// Upper-case hex with no separator, which is what `printHexBinary` produces.
#[test]
fn the_computed_checksum_is_upper_case_hex() {
    assert_eq!(print_hex_binary(&[0xd8, 0x1c, 0x87, 0x84]), "D81C8784");
    assert_eq!(print_hex_binary(&[0x00, 0x0f, 0xff]), "000FFF");
}

/// With no `-O` the destination is the source's file name as a RELATIVE path, so the copy lands in
/// the working directory rather than beside the source.
#[test]
fn the_default_destination_is_relative() {
    let text = golden();
    let source = Path::new("/somewhere/else/funcotator_dataSources.tar.gz");
    let destination = output_location(None, source);
    assert_eq!(destination, Path::new("funcotator_dataSources.tar.gz"));
    assert!(destination.is_relative());
    assert_eq!(
        output_location(Some(Path::new("/out/copied.tar.gz")), source),
        Path::new("/out/copied.tar.gz")
    );
    // The golden's default run landed beside the harness, not beside the archive.
    assert!(text.contains("file\tdefault-output\t<cwd>/funcotator_dataSources.tar.gz="));
}

/// The startup checks run in their own order, so the sha256 override alone is reported as a
/// missing data source rather than as an incomplete pair.
#[test]
fn the_startup_checks_misreport_which_argument_was_missing() {
    let text = golden();

    let cases = [
        (
            Arguments {
                hg38: true,
                ..Arguments::default()
            },
            DownloadError::NoDataSource,
            "no-source",
        ),
        (
            Arguments {
                somatic: true,
                ..Arguments::default()
            },
            DownloadError::NoReference,
            "no-reference",
        ),
        (
            Arguments {
                testing_data_sources_path: Some("<dir>/archive.tar.gz".to_string()),
                ..Arguments::default()
            },
            DownloadError::IncompleteTestingArguments,
            "testing-path-alone",
        ),
        (
            // The argument that is actually incomplete is the sha256 one, and the message names a
            // data source.
            Arguments {
                testing_sha256_path: Some("<dir>/plain.sha256".to_string()),
                ..Arguments::default()
            },
            DownloadError::NoDataSource,
            "testing-sha-alone",
        ),
    ];
    for (arguments, expected, label) in cases {
        let error = startup(&arguments).expect_err(label);
        assert_eq!(error, expected, "{label}");
        assert_eq!(
            format!("{}:{}", error.java_class(), error.message()),
            refusal(&text, label),
            "{label}"
        );
    }

    // Both testing arguments together are accepted, and are what every measured run used.
    assert!(startup(&Arguments {
        testing_data_sources_path: Some("a".to_string()),
        testing_sha256_path: Some("b".to_string()),
        ..Arguments::default()
    })
    .is_ok());
}

/// Overwriting is the only thing that lets a copy land on an existing file, and until it does the
/// old contents are still there.
#[test]
fn an_existing_destination_is_refused_until_overwrite() {
    let text = golden();
    let archive = fixture(&text, "archive-sha256");
    // Refused: what is on disk afterwards is still the file that was there before.
    assert_ne!(landed(&text, "existing-refused"), archive);
    assert!(refusal(&text, "existing-refused").contains("Output data sources file already exists!"));
    // Overwritten: the archive.
    assert!(ran(&text, "existing-overwritten"));
    assert_eq!(landed(&text, "existing-overwritten"), archive);
}
