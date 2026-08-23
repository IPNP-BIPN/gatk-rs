//! Conformance for `PrintSVEvidence` against GATK 4.6.2.0, compared as the whole merged file of
//! every run.
//!
//! Golden from `tools/readfilter-conformance/PrintSVEvidenceDump.java`.
//!
//! # What this suite is for
//!
//!  * **merging being only a widening**, and the refusal when it is not;
//!  * **the sample list being alphabetical rather than in file order**;
//!  * **`--sample-names` subsetting and reordering**, and an unknown name becoming a column of -1;
//!  * **the output's header being the header of no input**;
//!  * **and the two type refusals**, one of which fires before the tool's own check.

use gatk_corpus as corpus;
use gatk_tools::print_sv_evidence::{
    check_types, run, sample_names, write, DepthEvidence, EvidenceFile, PrintError,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/print_sv_evidence.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn input(text: &str, label: &str, name: &str) -> String {
    let prefix = format!("input\t{label}\t{name}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries input/{label}/{name}")),
    )
}

fn merged(text: &str, label: &str) -> String {
    let prefix = format!("merged\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries merged/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

/// One `.rd.txt` file read back: the header's samples and the records, one-based inside.
fn parse(text: &str) -> EvidenceFile {
    let mut samples = Vec::new();
    let mut records = Vec::new();
    for line in text.lines().filter(|line| !line.is_empty()) {
        let columns: Vec<&str> = line.split('\t').collect();
        if line.starts_with("#Chr") {
            samples = columns[3..].iter().map(|name| name.to_string()).collect();
            continue;
        }
        records.push(DepthEvidence {
            contig: columns[0].to_string(),
            start: columns[1].parse::<i32>().expect("a start") + 1,
            end: columns[2].parse().expect("an end"),
            counts: columns[3..]
                .iter()
                .map(|count| count.parse().expect("a count"))
                .collect(),
        });
    }
    EvidenceFile { samples, records }
}

/// The order the walker hands the records over: every file's records, merged by locus, ties in the
/// order the heap gives them, which for these fixtures is the order the files were named.
fn walked(files: &[EvidenceFile]) -> Vec<(usize, DepthEvidence)> {
    let mut all: Vec<(usize, DepthEvidence)> = files
        .iter()
        .enumerate()
        .flat_map(|(index, file)| {
            file.records
                .iter()
                .map(move |record| (index, record.clone()))
        })
        .collect();
    all.sort_by_key(|(index, record)| (record.start, record.end, *index));
    all
}

fn produced(text: &str, label: &str, names: &[&str], sources: &[&str]) -> String {
    let files: Vec<EvidenceFile> = sources
        .iter()
        .map(|name| parse(&input(text, label, name)))
        .collect();
    let requested: Vec<String> = names.iter().map(|name| name.to_string()).collect();
    let (samples, records) = run(&files, &walked(&files), &requested).expect("a run");
    write(&samples, &records)
}

#[test]
fn every_merged_file_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, names, sources) in [
        ("widening", vec![], vec!["a", "b"]),
        ("subset-one", vec!["beta"], vec!["a", "b"]),
        ("unknown-sample", vec!["gamma"], vec!["a", "b"]),
        ("reordered", vec!["beta", "alpha"], vec!["a", "b"]),
        ("disjoint-bins", vec![], vec!["a", "b"]),
        ("three-samples", vec![], vec!["a", "b", "c"]),
    ] {
        assert_eq!(
            produced(&text, label, &names, &sources),
            merged(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 6, "the golden's outputs");
}

/// Two files that both report a sample at one bin are refused, and the message names the sample
/// one-based and the bin one-based while the file writes it zero-based.
#[test]
fn merging_the_same_sample_twice_is_refused() {
    let text = golden();
    let files: Vec<EvidenceFile> = ["a", "b"]
        .iter()
        .map(|name| parse(&input(&text, "same-sample-twice", name)))
        .collect();
    let error = run(&files, &walked(&files), &[]).expect_err("a refused run");
    assert_eq!(
        error,
        PrintError::MultipleSources {
            sample: 1,
            contig: "chr1".to_string(),
            start: 1,
            end: 100
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "same-sample-twice")
    );
    // The file writes that same bin as `0 100`.
    assert!(input(&text, "same-sample-twice", "a").contains("chr1\t0\t100\t"));
}

/// The union is a TreeSet, so the columns are alphabetical whatever order the files came in.
#[test]
fn the_sample_list_is_alphabetical_not_file_order() {
    let text = golden();
    let files: Vec<EvidenceFile> = ["a", "b", "c"]
        .iter()
        .map(|name| parse(&input(&text, "three-samples", name)))
        .collect();
    assert_eq!(
        files
            .iter()
            .flat_map(|file| file.samples.clone())
            .collect::<Vec<String>>(),
        vec!["zulu", "alpha", "mike"]
    );
    assert_eq!(sample_names(&[], &files), vec!["alpha", "mike", "zulu"]);
    assert!(merged(&text, "three-samples").starts_with("#Chr\tStart\tEnd\talpha\tmike\tzulu\n"));
}

/// A name no file knows is a column of -1, not a refusal.
#[test]
fn an_unknown_sample_is_a_column_of_missing_data() {
    let text = golden();
    let produced = produced(&text, "unknown-sample", &["gamma"], &["a", "b"]);
    assert!(produced.contains("chr1\t0\t100\t-1\n"));
    assert_eq!(produced, merged(&text, "unknown-sample"));
}

#[test]
fn the_two_type_refusals_match_the_golden() {
    let text = golden();

    // An extension whose codec produces another SV feature type.
    let error = check_types(
        "<dir>/merged.baf.txt",
        &["<dir>/wrong-sv-type-a.rd.txt".to_string()],
    )
    .expect_err("an incompatible input");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "wrong-sv-type")
    );

    // And one that names no feature type at all, which fails first.
    let error = check_types(
        "<dir>/merged.vcf",
        &["<dir>/not-an-sv-type-a.rd.txt".to_string()],
    )
    .expect_err("no output codec");
    assert_eq!(
        error,
        PrintError::NoOutputCodec {
            path: "<dir>/merged.vcf".to_string()
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "not-an-sv-type")
    );
}
