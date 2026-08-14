//! Conformance for `CountVariants` against GATK 4.6.2.0, compared as the count of every run, the
//! bytes of every output file, and the class and message of every refusal.
//!
//! Golden from `tools/readfilter-conformance/CountVariantsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the count reaches no stream without `-O`**, whatever the class documentation says: the
//!    golden's `no-output-argument` row is a file that was never written, not an empty one;
//!  * **the file has no trailing newline and truncates what was there**, so a count of 5 over ten
//!    bytes leaves one byte;
//!  * **a record is selected by its whole span**, `END` or the length of `REF`, and not by `POS`;
//!  * **a record spanning two intervals is counted once**;
//!  * **and the refusal for an unwritable `-O` is the path alone**, the overload the call site
//!    reaches carrying no reason at all.
//!
//! # What is compared, and what is not
//!
//! The `interval-off-the-dictionary` row is asserted as a message rather than replayed: that
//! refusal belongs to the interval argument layer, which has its own suite, and nothing in this
//! tool decides it.

use gatk_corpus as corpus;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::variant_source::Located;
use gatk_tools::count_variants::{count, output_bytes, CountVariantsError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/count_variants.txt.gz"),
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
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
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

/// One record, decoded as far as the traversal looks at it.
struct Variant {
    contig: String,
    start: i32,
    stop: i32,
}

impl Located for Variant {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn start(&self) -> i32 {
        self.start
    }
    fn stop(&self) -> i32 {
        self.stop
    }
}

/// The whole input file, as the dump escaped it.
fn input(text: &str, label: &str) -> String {
    unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no input {label}"))[1],
    )
}

/// The records of one input, whose stop is `END` when there is one and `start + len(REF) - 1`
/// otherwise. That is the span a Tribble query matches, and it is why `-L` reaches records whose
/// position it does not.
fn records(text: &str, label: &str) -> Vec<Variant> {
    input(text, label)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let start: i32 = field[1].parse().expect("a position");
            let end = field[7]
                .split(';')
                .filter_map(|entry| entry.split_once('='))
                .find(|(key, _)| *key == "END")
                .and_then(|(_, value)| value.parse::<i32>().ok());
            Variant {
                contig: field[0].to_string(),
                start,
                stop: end.unwrap_or(start + field[3].len() as i32 - 1),
            }
        })
        .collect()
}

/// The contig lengths the header declares, which is what `-L chr2` expands to.
fn contig_length(text: &str, label: &str, contig: &str) -> i32 {
    input(text, label)
        .lines()
        .filter(|line| line.starts_with("##contig=<ID="))
        .find_map(|line| {
            let body = line.trim_start_matches("##contig=<").trim_end_matches('>');
            let mut id = None;
            let mut length = None;
            for entry in body.split(',') {
                match entry.split_once('=') {
                    Some(("ID", value)) => id = Some(value),
                    Some(("length", value)) => length = value.parse::<i32>().ok(),
                    _ => {}
                }
            }
            match (id, length) {
                (Some(name), Some(length)) if name == contig => Some(length),
                _ => None,
            }
        })
        .unwrap_or_else(|| panic!("no contig {contig}"))
}

fn interval(contig: &str, start: i32, end: i32) -> SimpleInterval {
    SimpleInterval::new(contig, start, end).expect("a valid interval")
}

/// Every run: which file it reads, the merged `-L` list, whether the input is indexed, and whether
/// `-O` was given.
struct Run {
    file: &'static str,
    intervals: Vec<SimpleInterval>,
    indexed: bool,
    has_output: bool,
}

fn setup(text: &str, label: &str) -> Run {
    let plain = |file: &'static str, has_output: bool| Run {
        file,
        intervals: Vec::new(),
        indexed: true,
        has_output,
    };
    match label {
        "plain-no-output" => plain("plain", false),
        "plain" => plain("plain", true),
        "filtered-only" => plain("filtered-only", true),
        "empty" => plain("empty", true),
        "overwrite" => plain("plain", true),
        "span-by-end" => Run {
            intervals: vec![interval("chr1", 300, 310)],
            ..plain("spanning", true)
        },
        "span-by-ref-length" => Run {
            intervals: vec![interval("chr1", 605, 606)],
            ..plain("spanning", true)
        },
        "span-missed" => Run {
            intervals: vec![interval("chr1", 500, 510)],
            ..plain("spanning", true)
        },
        "two-intervals-one-record" => Run {
            intervals: vec![interval("chr1", 150, 160), interval("chr1", 350, 360)],
            ..plain("spanning", true)
        },
        "interval-matches-nothing" => Run {
            intervals: vec![interval("chr1", 900, 950)],
            ..plain("plain", true)
        },
        "interval-selects-contig" => Run {
            intervals: vec![interval(
                "chr2",
                1,
                contig_length(text, "two-contigs", "chr2"),
            )],
            ..plain("two-contigs", true)
        },
        "interval-without-index" => Run {
            intervals: vec![interval("chr1", 100, 200)],
            indexed: false,
            ..plain("unindexed", true)
        },
        other => panic!("no setup for {other}"),
    }
}

/// The path the dump gave `-O`, which is what the second refusal's message is.
const OUTPUT_IS_A_DIRECTORY: &str = "countvariants-dump/.";

#[test]
fn every_counted_run_matches_the_golden() {
    let text = golden();
    let counted = rows(&text, "count");
    assert_eq!(counted.len(), 11, "eleven of the fourteen runs finish");

    for row in counted {
        let (label, expected, class) = (row[0], row[1], row[2]);
        let run = setup(&text, label);
        let intervals = if run.intervals.is_empty() {
            None
        } else {
            Some(run.intervals.as_slice())
        };
        let mine = count(
            &records(&text, run.file),
            intervals,
            run.indexed,
            &format!("countvariants-dump/{}.vcf", run.file),
        )
        .unwrap_or_else(|_| panic!("{label} finishes"));
        assert_eq!(mine.to_string(), expected, "the count of {label}");
        assert_eq!(class, "java.lang.Long", "the count of {label} is a Long");
    }
}

#[test]
fn every_output_file_matches_the_golden_byte_for_byte() {
    let text = golden();
    for row in rows(&text, "file") {
        let label = row[0];
        // The file written before the overwriting run, which is not a run of the tool.
        if label == "before-overwrite" {
            assert_eq!(row[1], "present");
            assert_eq!(row[2], "10");
            continue;
        }
        let run = setup(&text, label);
        let intervals = if run.intervals.is_empty() {
            None
        } else {
            Some(run.intervals.as_slice())
        };
        let Ok(counted) = count(
            &records(&text, run.file),
            intervals,
            run.indexed,
            &format!("countvariants-dump/{}.vcf", run.file),
        ) else {
            panic!("{label} finishes");
        };

        match row[1] {
            // A run with no -O writes nothing at all, which is not the same as an empty file.
            "no-output-argument" => {
                assert!(!run.has_output);
                assert_eq!(output_bytes(counted, run.has_output), None, "{label}");
            }
            "present" => {
                let bytes = output_bytes(counted, run.has_output).expect("a file");
                assert_eq!(bytes.len().to_string(), row[2], "the length of {label}");
                assert_eq!(String::from_utf8(bytes).expect("utf-8"), row[3], "{label}");
            }
            other => panic!("no such file state {other}"),
        }
    }
}

#[test]
fn the_output_has_no_trailing_newline() {
    let text = golden();
    let plain = rows(&text, "file")
        .into_iter()
        .find(|row| row[0] == "plain")
        .expect("the plain run");
    assert_eq!(plain[2], "1", "one byte for a count of five");
    assert_eq!(plain[3], "5");
}

#[test]
fn the_overwriting_run_leaves_one_byte_of_ten() {
    let text = golden();
    let before = rows(&text, "file")
        .into_iter()
        .find(|row| row[0] == "before-overwrite")
        .expect("the file written first");
    let after = rows(&text, "file")
        .into_iter()
        .find(|row| row[0] == "overwrite")
        .expect("the overwriting run");
    assert_eq!(before[2], "10");
    assert_eq!(after[2], "1");
    assert_eq!(after[3], "5");
}

#[test]
fn both_refusals_carry_the_references_class_and_words() {
    let text = golden();
    let errors = rows(&text, "error");
    assert_eq!(errors.len(), 3, "three refusals");

    let refusal = |label: &str| -> (String, String) {
        let row = errors
            .iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no refusal {label}"));
        let whole = row[1..].join("\t");
        let (class, message) = whole.split_once(':').expect("class and message");
        (class.to_string(), message.to_string())
    };

    let run = setup(&text, "interval-without-index");
    let refused = count(
        &records(&text, run.file),
        Some(run.intervals.as_slice()),
        run.indexed,
        "countvariants-dump/unindexed.vcf",
    )
    .expect_err("no index");
    let (class, message) = refusal("interval-without-index");
    assert_eq!(refused.class(), class);
    assert_eq!(refused.message(), message);

    let unwritable = CountVariantsError::CouldNotCreateOutputFile {
        path: OUTPUT_IS_A_DIRECTORY.to_string(),
    };
    let (class, message) = refusal("output-is-a-directory");
    assert_eq!(unwritable.class(), class);
    assert_eq!(unwritable.message(), message, "the path and nothing else");
    assert_eq!(message, OUTPUT_IS_A_DIRECTORY);
}

/// The third refusal is the interval argument layer's, and it is asserted rather than replayed.
#[test]
fn an_interval_off_the_dictionary_is_refused_before_this_tool() {
    let text = golden();
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "interval-off-the-dictionary")
        .expect("the third refusal");
    let whole = row[1..].join("\t");
    assert!(
        whole.starts_with(
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedGenomeLoc:"
        ),
        "{whole}"
    );
    assert!(whole.contains(r#"Query interval "chr3:1-10" is not valid for this input."#));
}
