//! Conformance for `GetSampleName` against GATK 4.6.2.0, compared as text **and length**.
//!
//! Golden from `tools/readfilter-conformance/GetSampleNameDump.java`. The seven input BAMs and their
//! indexes travel in full, base64, so the port reads the same headers.
//!
//! # What this suite is for
//!
//!  * **there is no trailing newline**, which only the byte count shows;
//!  * **the order is the header's and repeats collapse**, because `distinct()` keeps the first;
//!  * **a header with no read groups takes the second refusal**, not the first;
//!  * **a read group with no `SM` writes the four letters `null`** rather than refusing;
//!  * **URL encoding turns a space into `+`**, not `%20`;
//!  * **and the tool reads only the header**, so a BAM with no records still writes its sample.

use gatk_corpus as corpus;
use gatk_engine::reads::ReadsDataSource;
use gatk_tools::get_sample_name::{get_sample_name, GetSampleNameError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/get_sample_name.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn install(text: &str, dir: &std::path::Path) {
    std::fs::create_dir_all(dir).expect("a scratch directory");
    for row in rows(text, "fixture") {
        std::fs::write(
            dir.join(format!("{}.bam", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture bam");
    }
    for row in rows(text, "fixtureindex") {
        if row[1] == "absent" {
            continue;
        }
        std::fs::write(
            dir.join(format!("{}.bai", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture index");
    }
}

/// The fixture each labelled run used, and whether it asked for URL encoding.
fn configuration(label: &str) -> (&str, bool) {
    match label {
        "single-encoded" => ("single", true),
        "special-encoded" => ("special", true),
        "two-encoded" => ("two", true),
        // The run with the argument left out takes the field's own default.
        "no-encoding-argument" => ("single", false),
        other => (other, false),
    }
}

fn run(dir: &std::path::Path, label: &str) -> Result<String, GetSampleNameError> {
    let (fixture, url_encode) = configuration(label);
    let bam = dir.join(format!("{fixture}.bam"));
    let bai = dir.join(format!("{fixture}.bai"));
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    get_sample_name(&source, url_encode)
}

#[test]
fn every_written_file_is_the_reference() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-getsample-{}", std::process::id()));
    install(&text, &dir);

    let samples = rows(&text, "sample");
    let lengths = rows(&text, "bytes");
    assert_eq!(samples.len(), 10, "ten runs finish");

    for row in &samples {
        let label = row[0];
        let ours = run(&dir, label).unwrap_or_else(|error| {
            panic!("{label} was refused: {}", error.message());
        });
        assert_eq!(ours, unescape(row[1]), "sample/{label}");

        // The byte count is the only place the missing trailing newline shows.
        let expected: usize = lengths
            .iter()
            .find(|other| other[0] == label)
            .expect("every finished run dumps its length")[1]
            .parse()
            .expect("a number");
        assert_eq!(ours.len(), expected, "bytes/{label}");
    }
}

#[test]
fn the_refusal_is_the_reference() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-getsample-error-{}", std::process::id()));
    install(&text, &dir);

    let errors = rows(&text, "error");
    assert_eq!(errors.len(), 1, "one refusal");
    for row in errors {
        let error = run(&dir, row[0]).expect_err("this run is refused");
        assert_eq!(
            format!(
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput:Bad input: {}",
                error.message()
            ),
            row[1],
            "error/{}",
            row[0]
        );
    }
}
