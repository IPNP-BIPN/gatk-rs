//! Conformance for `PrintReadCounts` against GATK 4.6.2.0, compared as the whole of every file the
//! run left behind.
//!
//! Golden from `tools/readfilter-conformance/PrintReadCountsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the built header**, which is the dictionary plus one read group always called
//!    `GATKCopyNumber`, so a counts file fed back in loses everything else its header carried;
//!  * **the change of base**, a bin written `0 100` coming back `1 100`;
//!  * **the raw concatenation of the output names**, and the same names in the output file list;
//!  * **the two writers one duplicated column name opens over one path**;
//!  * **the half header a crash leaves behind**, the SAM header flushed and the column line lost;
//!  * **and the refusals**, including the two the feature reader wraps twice.
//!
//! The `counts-header-late` run is not replayed here: it fails inside `SimpleCountCodec`'s own
//! header reader, which is the codec's rule rather than the tool's, and nothing of the tool runs.

use gatk_corpus as corpus;
use gatk_tools::print_read_counts::{
    decode_depth, depth_header_samples, output_name, read_sample_name, reader_error_chain, run,
    CountsFile, DepthFile, Disk, Input, Interval, PrintError, SimpleCount,
};
use htsjdk_bam::header::{ProgramRecord, ReadGroup, SamHeader, SequenceRecord};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/print_read_counts.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, prefix: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix))
            .unwrap_or_else(|| panic!("the golden carries {prefix}")),
    )
}

fn input(text: &str, label: &str, name: &str) -> String {
    field(text, &format!("input\t{label}\t{name}="))
}

fn out(text: &str, label: &str, name: &str) -> String {
    field(text, &format!("out\t{label}\t{name}="))
}

fn list(text: &str, label: &str) -> String {
    field(text, &format!("list\t{label}="))
}

/// Every `out` line of one run, in the order the dump printed them, which is by name.
fn outputs(text: &str, label: &str) -> Vec<(String, String)> {
    let prefix = format!("out\t{label}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|rest| {
            let (name, content) = rest.split_once('=').expect("an out line names its file");
            (name.to_string(), unescape(content))
        })
        .collect()
}

/// The `error` line, as `<class>:<message>`.
fn refusal(text: &str, label: &str) -> String {
    field(text, &format!("error\t{label}\t"))
}

/// The `error` line and every `cause` line under it, outermost first.
fn chain(text: &str, label: &str) -> Vec<String> {
    let mut chain = vec![refusal(text, label)];
    let prefix = format!("cause\t{label}\t");
    chain.extend(
        text.lines()
            .filter_map(|line| line.strip_prefix(prefix.as_str()))
            .map(unescape),
    );
    chain
}

/// A SAM text header read back into the record classes, which is what the counts codec hands the
/// tool and what the dictionary file gives the `.rd.txt` runs.
fn parse_header(text: &str) -> SamHeader {
    let mut header = SamHeader::new();
    header.attributes = Default::default();
    for line in text.lines().filter(|line| line.starts_with('@')) {
        let mut fields = line.split('\t');
        let tag = fields.next().expect("a line tag");
        let rest: Vec<&str> = fields.collect();
        let attribute = |field: &str| {
            let (key, value) = field.split_once(':').expect("a tagged field");
            (key.to_string(), value.to_string())
        };
        match tag {
            "@HD" => {
                for pair in &rest {
                    let (key, value) = attribute(pair);
                    header.attributes.set(&key, &value);
                }
            }
            "@SQ" => {
                let mut record = SequenceRecord::new("", 0);
                for pair in &rest {
                    let (key, value) = attribute(pair);
                    match key.as_str() {
                        "SN" => record.name = value,
                        "LN" => record.length = value.parse().expect("a length"),
                        _ => record.attributes.set(&key, &value),
                    }
                }
                header.sequences.push(record);
            }
            "@RG" | "@PG" => {
                let mut id = String::new();
                let mut attributes = Vec::new();
                for pair in &rest {
                    let (key, value) = attribute(pair);
                    if key == "ID" {
                        id = value;
                    } else {
                        attributes.push((key, value));
                    }
                }
                if tag == "@RG" {
                    let mut group = ReadGroup::new(&id);
                    for (key, value) in attributes {
                        group.attributes.set(&key, &value);
                    }
                    header.read_groups.push(group);
                } else {
                    let mut program = ProgramRecord::new(&id);
                    for (key, value) in attributes {
                        program.attributes.set(&key, &value);
                    }
                    header.programs.push(program);
                }
            }
            "@CO" => header.add_comment(line),
            _ => panic!("the golden's headers hold no {tag}"),
        }
    }
    header
}

/// The dictionary every `.rd.txt` run is handed, read from the `.dict` the dump printed.
fn dictionary(text: &str) -> Vec<SequenceRecord> {
    parse_header(&field(text, "dict\tcounts=")).sequences
}

fn depth_file(text: &str) -> DepthFile {
    let mut samples = Vec::new();
    let mut records = Vec::new();
    for line in text.lines().filter(|line| !line.is_empty()) {
        if line.starts_with("#Chr") {
            samples = depth_header_samples(line);
            continue;
        }
        records.push(decode_depth(line).expect("a record"));
    }
    DepthFile { samples, records }
}

fn counts_file(text: &str) -> CountsFile {
    let header = parse_header(text);
    let records = text
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('@') && !line.starts_with("CONTIG"))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            SimpleCount {
                contig: columns[0].to_string(),
                start: columns[1].parse().expect("a start"),
                end: columns[2].parse().expect("an end"),
                count: columns[3].parse().expect("a count"),
            }
        })
        .collect();
    CountsFile { header, records }
}

/// The prefix the dump handed each run: its own directory, so the names in the output file list
/// are absolute and the names the dump printed are the last component.
fn prefix(label: &str, glued: &str) -> String {
    format!("<dir>/{label}/{glued}")
}

/// Every `.counts.tsv` the run left, keyed by its last component, which is what the dump printed.
fn written(disk: &Disk, label: &str) -> Vec<(String, String)> {
    let directory = format!("<dir>/{label}/");
    disk.files()
        .into_iter()
        .filter(|(name, _)| name.ends_with(".counts.tsv"))
        .map(|(name, content)| (name.trim_start_matches(&directory).to_string(), content))
        .collect()
}

fn depth_run(
    text: &str,
    label: &str,
    glued: &str,
    list_path: bool,
    intervals: &[Interval],
) -> Disk {
    let file = depth_file(&input(text, label, "input.rd.txt"));
    let list = format!("<dir>/{label}/outputs.list");
    run(
        &Input::Depth(file),
        Some(&dictionary(text)),
        &prefix(label, glued),
        list_path.then_some(list.as_str()),
        intervals,
    )
    .disk
}

#[test]
fn two_samples_become_two_files_one_based() {
    let text = golden();
    let disk = depth_run(&text, "rd-two-samples", "", true, &[]);
    assert_eq!(
        written(&disk, "rd-two-samples"),
        outputs(&text, "rd-two-samples")
    );
    // The bins were written 0 100 and 100 200; both files carry 1 100 and 101 200.
    assert!(out(&text, "rd-two-samples", "alpha.counts.tsv").contains("chr1\t1\t100\t11\n"));
    assert_eq!(
        disk.read("<dir>/rd-two-samples/outputs.list"),
        Some(list(&text, "rd-two-samples"))
    );
}

#[test]
fn an_rd_header_carries_no_dictionary() {
    let text = golden();
    let file = depth_file(&input(&text, "rd-no-dictionary", "input.rd.txt"));
    let outcome = run(
        &Input::Depth(file),
        None,
        "<dir>/rd-no-dictionary/",
        None,
        &[],
    );
    let error = outcome.error.expect("a refusal");
    assert_eq!(error, PrintError::NoDictionary);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "rd-no-dictionary")
    );
    assert!(outcome.disk.files().is_empty());
}

#[test]
fn the_prefix_is_concatenated_raw() {
    let text = golden();
    let disk = depth_run(&text, "rd-glued-prefix", "sample-", false, &[]);
    assert_eq!(
        written(&disk, "rd-glued-prefix"),
        outputs(&text, "rd-glued-prefix")
    );
    assert_eq!(
        output_name("<dir>/rd-glued-prefix/sample-", "alpha"),
        "<dir>/rd-glued-prefix/sample-alpha.counts.tsv"
    );
}

#[test]
fn one_sample_is_one_whole_file() {
    let text = golden();
    let disk = depth_run(&text, "rd-one-sample", "", false, &[]);
    assert_eq!(
        written(&disk, "rd-one-sample"),
        outputs(&text, "rd-one-sample")
    );
}

#[test]
fn a_duplicated_column_name_opens_two_writers_over_one_path() {
    let text = golden();
    let disk = depth_run(&text, "rd-duplicate-sample", "", true, &[]);
    assert_eq!(
        written(&disk, "rd-duplicate-sample"),
        outputs(&text, "rd-duplicate-sample")
    );
    // One file, carrying the second writer's count, and a list that names it twice.
    assert_eq!(written(&disk, "rd-duplicate-sample").len(), 1);
    assert_eq!(
        disk.read("<dir>/rd-duplicate-sample/outputs.list"),
        Some(list(&text, "rd-duplicate-sample"))
    );
}

#[test]
fn a_short_record_crashes_and_leaves_half_a_header() {
    let text = golden();
    let file = depth_file(&input(&text, "rd-short-record", "input.rd.txt"));
    let outcome = run(
        &Input::Depth(file),
        Some(&dictionary(&text)),
        &prefix("rd-short-record", ""),
        None,
        &[],
    );
    let error = outcome.error.expect("a refusal");
    assert_eq!(
        error,
        PrintError::IndexOutOfBounds {
            index: 1,
            length: 1
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "rd-short-record")
    );
    let written = written(&outcome.disk, "rd-short-record");
    assert_eq!(written, outputs(&text, "rd-short-record"));
    // The column line never reached the disk, though the SAM header did.
    assert!(written.iter().all(|(_, text)| !text.contains("CONTIG")));
}

#[test]
fn a_long_record_drops_the_extra_count() {
    let text = golden();
    let disk = depth_run(&text, "rd-long-record", "", false, &[]);
    assert_eq!(
        written(&disk, "rd-long-record"),
        outputs(&text, "rd-long-record")
    );
}

#[test]
fn intervals_subset_the_records_not_the_files() {
    let text = golden();
    let disk = depth_run(
        &text,
        "rd-intervals",
        "",
        false,
        &[Interval {
            contig: "chr2".to_string(),
            start: 1,
            end: 1040,
        }],
    );
    assert_eq!(
        written(&disk, "rd-intervals"),
        outputs(&text, "rd-intervals")
    );
    assert_eq!(written(&disk, "rd-intervals").len(), 2);
}

#[test]
fn a_counts_file_keeps_its_records_and_loses_its_header() {
    let text = golden();
    let file = counts_file(&input(&text, "counts-round-trip", "input.counts.tsv"));
    // What the input carried and the output will not.
    assert_eq!(file.header.read_groups[0].id, "not-the-cnv-id");
    assert_eq!(file.header.programs.len(), 1);
    assert_eq!(file.header.comments.len(), 1);
    let list_path = "<dir>/counts-round-trip/outputs.list";
    let outcome = run(
        &Input::Counts(file),
        None,
        &prefix("counts-round-trip", ""),
        Some(list_path),
        &[],
    );
    assert!(outcome.error.is_none());
    assert_eq!(
        written(&outcome.disk, "counts-round-trip"),
        outputs(&text, "counts-round-trip")
    );
    assert_eq!(
        outcome.disk.read(list_path),
        Some(list(&text, "counts-round-trip"))
    );
}

#[test]
fn two_read_groups_naming_two_samples_are_refused_through_the_reader() {
    let text = golden();
    let file = counts_file(&input(&text, "counts-two-samples", "input.counts.tsv"));
    let inner = read_sample_name(&file.header).expect_err("a refusal");
    assert_eq!(
        inner,
        PrintError::ManySampleNames(vec!["gamma".to_string(), "delta".to_string()])
    );
    let path = "<dir>/counts-two-samples/input.counts.tsv";
    assert_eq!(
        reader_error_chain(path, &inner)
            .into_iter()
            .map(|(class, message)| format!("{class}:{message}"))
            .collect::<Vec<_>>(),
        chain(&text, "counts-two-samples")
    );
}

#[test]
fn a_read_group_with_no_sample_is_refused_by_the_emptiness_test() {
    let text = golden();
    let file = counts_file(&input(&text, "counts-no-sample", "input.counts.tsv"));
    // readSampleName itself is content with the lone null it found.
    assert_eq!(read_sample_name(&file.header), Ok(None));
    let outcome = run(
        &Input::Counts(file),
        None,
        &prefix("counts-no-sample", ""),
        None,
        &[],
    );
    let inner = outcome.error.expect("a refusal");
    assert_eq!(inner, PrintError::NullSampleName);
    let path = "<dir>/counts-no-sample/input.counts.tsv";
    assert_eq!(
        reader_error_chain(path, &inner)
            .into_iter()
            .map(|(class, message)| format!("{class}:{message}"))
            .collect::<Vec<_>>(),
        chain(&text, "counts-no-sample")
    );
}

#[test]
fn any_other_feature_type_never_reaches_the_tool() {
    let text = golden();
    let path = "<dir>/wrong-feature-type/input.baf.txt";
    let error = PrintError::WrongFeatureType {
        path: path.to_string(),
    };
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "wrong-feature-type")
    );
}
