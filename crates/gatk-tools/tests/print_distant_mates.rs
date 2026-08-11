//! Conformance for `PrintDistantMates` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/PrintDistantMatesDump.java`. The five output BAMs and
//! their indexes travel in full, base64, as the rest of this archetype's do, and so does the input
//! fixture with its index: a test that built its own index would be inventing part of the input.
//!
//! # What this suite is for
//!
//! The sixth whole tool of the record-transform archetype, and three of the four things it pins
//! are things the archetype hides.
//!
//!  * **the writer is opened with `preSorted = false`.** Every other tool here passes `true` and
//!    writes the traversal order. This one moves each read onto its mate, so the traversal order
//!    is wrong by construction, and the fixture is built so the two orders cannot agree by
//!    accident: three reads leave the traversal for `chr2:600`, `chr1:2500`, `chr2:150` and are
//!    written `chr1:2500`, `chr2:150`, `chr2:600`;
//!  * **the default read filters are extended**, not taken and not replaced. Both lists are in the
//!    golden, so the difference between `PrintReads`'s default and this tool's is read rather than
//!    argued;
//!  * **the same `OA` tag is spelled two ways.** One read with no `NM` goes through both tools and
//!    the golden carries `chr1,200,+,10M,60,;` from this one against `chr1,200,+,10M,60,null;`
//!    from `AddOriginalAlignmentTags`;
//!  * **`undoDistantMateAlterations` is the inverse**, tag block included, and its guard is
//!    narrower than the `catch` suggests: an `OA` naming a contig the dictionary does not hold
//!    escapes it.
//!
//! The command line lands in the `@PG` record's `CL`, so it is read out of the golden and handed
//! to the port rather than reconstructed: it carries the paths of the run that produced it.

use gatk_corpus as corpus;
use gatk_engine::interval::{self, SimpleInterval};
use gatk_engine::reads::ReadsDataSource;
use gatk_tools::add_original_alignment_tags as add_oa;
use gatk_tools::print_distant_mates as tool;
use gatk_tools::sam_output::Options;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/print_distant_mates.txt.gz"),
    )
}

/// Rows of one kind, split on tabs, with the kind dropped.
fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// Rows of one kind that carry `<label>\t<value>`.
fn pairs<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str)> {
    rows(text, kind)
        .into_iter()
        .map(|row| (row[0], row[1]))
        .collect()
}

/// Rows of one kind that carry `<tool>\t<label>\t<value>`.
fn triples<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str, &'a str)> {
    rows(text, kind)
        .into_iter()
        .map(|row| (row[0], row[1], row[2]))
        .collect()
}

fn lookup<'a>(rows: &[(&'a str, &'a str, &'a str)], tool: &str, label: &str) -> Option<&'a str> {
    rows.iter()
        .find(|(t, l, _)| *t == tool && *l == label)
        .map(|(_, _, value)| *value)
}

/// What each labelled run was given. A label is a configuration and the row carries nothing to
/// derive it from, so it is written here beside the dump that produced it.
struct Configuration {
    intervals: &'static [&'static str],
    /// `--mate-too-distant-length`, `MateDistantReadFilter`'s one argument.
    mate_too_distant_length: i32,
    create_index: bool,
    program_record: bool,
}

fn configuration(label: &str) -> Configuration {
    let base = Configuration {
        intervals: &[],
        mate_too_distant_length: tool::DEFAULT_MATE_TOO_DISTANT_THRESHOLD,
        create_index: true,
        program_record: true,
    };
    match label {
        "all" => base,
        // Two of the three survivors, and the one it drops is the one that would have been written
        // second: the interval changes the output order as well as its contents.
        "chr1head" => Configuration {
            intervals: &["chr1:1-250"],
            ..base
        },
        "noindex" => Configuration {
            create_index: false,
            ..base
        },
        "nopg" => Configuration {
            program_record: false,
            ..base
        },
        // At 5, the read whose mate is ten bases away survives, and it is the read with no NM.
        "distant5" => Configuration {
            mate_too_distant_length: 5,
            ..base
        },
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The fixture, written out so the port's reader can open it.
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

fn source(dir: &std::path::Path) -> ReadsDataSource {
    ReadsDataSource::open(&dir.join("plain.bam"), &dir.join("plain.bai"))
        .expect("the fixture opens")
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-distantmates-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let outputs = triples(&text, "output");
    let indexes = triples(&text, "index");
    let command_lines = triples(&text, "commandline");
    assert_eq!(outputs.len(), command_lines.len());

    let mut compared = 0usize;
    for (tool_name, label, expected_base64) in &outputs {
        // The other tool's run is here for its OA spelling, not for its bytes: it has a suite of
        // its own, and comparing it twice would say nothing new.
        if *tool_name != "PrintDistantMates" {
            continue;
        }
        let config = configuration(label);
        let source = source(&dir);
        let header = source.header().clone();

        let intervals: Vec<SimpleInterval> = config
            .intervals
            .iter()
            .map(|text| interval::parse_interval(text, &header).expect("a parsable interval"))
            .collect();

        let threshold = config.mate_too_distant_length;
        let filter = move |read: &BamRecord| tool::default_read_filter(read, &header, threshold);

        let command_line =
            lookup(&command_lines, tool_name, label).expect("a command line for every output");
        let options = Options {
            intervals,
            create_output_bam_index: config.create_index,
            add_output_sam_program_record: config.program_record,
            command_line,
            ..Options::default()
        };

        let (ours, our_index) =
            tool::print_distant_mates(&source, &options, &filter).expect("the source reads");

        let expected = corpus::decode_base64(expected_base64);
        assert_eq!(
            ours.len(),
            expected.len(),
            "{tool_name}/{label}: output length differs"
        );
        if ours != expected {
            let at = ours
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            panic!("{tool_name}/{label}: first byte difference at offset {at}");
        }

        let expected_index =
            lookup(&indexes, tool_name, label).expect("an index row for every output");
        match (our_index, expected_index) {
            (None, "absent") => {}
            (Some(_), "absent") => {
                panic!("{tool_name}/{label}: the reference wrote no index and the port did")
            }
            (None, _) => {
                panic!("{tool_name}/{label}: the reference wrote an index and the port did not")
            }
            (Some(ours), expected) => {
                assert_eq!(
                    ours,
                    corpus::decode_base64(expected),
                    "{tool_name}/{label}: the .bai"
                );
            }
        }
        compared += 1;
    }

    assert_eq!(compared, 5, "five runs of this tool");
    println!("print-distant-mates: {compared} output files byte-identical");
}

/// `preSorted = false`, which is the one argument to the writer this archetype does not share.
///
/// A port that wrote the traversal order would produce a file with the same reads in it and a
/// different order, and an index built over that order. The golden's read rows are the output
/// order, so this asserts the order rather than the set.
#[test]
fn the_output_is_written_in_coordinate_order_not_traversal_order() {
    let text = golden();
    let dir =
        std::env::temp_dir().join(format!("gatk-rs-distantmates-order-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let expected: Vec<&str> = rows(&text, "reads")
        .into_iter()
        .filter(|row| row[0] == "PrintDistantMates" && row[1] == "all")
        .map(|row| row[2])
        .collect();
    assert_eq!(
        expected,
        ["r4", "r5", "r0"],
        "the golden lost the order this suite is for"
    );

    let source = source(&dir);
    let header = source.header().clone();
    let filter = {
        let header = header.clone();
        move |read: &BamRecord| {
            tool::default_read_filter(read, &header, tool::DEFAULT_MATE_TOO_DISTANT_THRESHOLD)
        }
    };
    let traversal = gatk_tools::read_walker::traverse(&source, &[], &filter).expect("the reads");
    let names: Vec<&str> = traversal.iter().map(|r| r.read_name.as_str()).collect();
    assert_eq!(
        names,
        ["r0", "r4", "r5"],
        "the traversal order is the input order"
    );

    let mut altered: Vec<BamRecord> = traversal
        .iter()
        .map(|read| tool::do_distant_mate_alterations(read, &header))
        .collect();
    tool::output_order(&mut altered);
    let sorted: Vec<&str> = altered.iter().map(|r| r.read_name.as_str()).collect();
    assert_eq!(sorted, expected);
}

/// The four filters this tool adds to the list it inherits, in the order it adds them.
#[test]
fn the_default_filters_extend_the_engine_s_rather_than_replacing_them() {
    let text = golden();
    let listed = |tool_name: &str| -> Vec<String> {
        rows(&text, "filters")
            .into_iter()
            .filter(|row| row[0] == tool_name)
            .map(|row| row[2].to_string())
            .collect()
    };
    let inherited = listed("PrintReads");
    let ours = listed("PrintDistantMates");

    assert_eq!(inherited, ["WellformedReadFilter"]);
    assert_eq!(ours, tool::DEFAULT_READ_FILTERS);
    assert!(
        ours.starts_with(&inherited),
        "the list is built from super's, so super's is its prefix"
    );
    assert_eq!(ours.len(), inherited.len() + 4);
}

/// Every read the reference wrote, rebuilt from the same input: the position it was moved to, the
/// flags, the cleared cigar and mapping quality, and both tags.
#[test]
fn every_read_is_moved_onto_its_mate() {
    let text = golden();
    let dir =
        std::env::temp_dir().join(format!("gatk-rs-distantmates-reads-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let source = source(&dir);
    let header = source.header().clone();
    let filter = {
        let header = header.clone();
        move |read: &BamRecord| {
            tool::default_read_filter(read, &header, tool::DEFAULT_MATE_TOO_DISTANT_THRESHOLD)
        }
    };
    let traversal = gatk_tools::read_walker::traverse(&source, &[], &filter).expect("the reads");
    let mut altered: Vec<BamRecord> = traversal
        .iter()
        .map(|read| tool::do_distant_mate_alterations(read, &header))
        .collect();
    tool::output_order(&mut altered);

    let expected: Vec<Vec<&str>> = rows(&text, "reads")
        .into_iter()
        .filter(|row| row[0] == "PrintDistantMates" && row[1] == "all")
        .collect();
    assert_eq!(altered.len(), expected.len());

    for (ours, row) in altered.iter().zip(&expected) {
        let name = row[2];
        assert_eq!(ours.read_name, name);
        assert_eq!(ours.flags.to_string(), row[3], "flags of {name}");
        assert_eq!(
            header.sequences[ours.reference_index as usize].name, row[4],
            "contig of {name}"
        );
        assert_eq!(ours.alignment_start.to_string(), row[5], "start of {name}");
        // `SAMRecord.NO_ALIGNMENT_CIGAR`, which prints as a star.
        assert_eq!(row[6], "*", "the reference cleared the cigar of {name}");
        assert!(ours.cigar.is_empty(), "cigar of {name}");
        assert_eq!(ours.mapping_quality.to_string(), row[7], "mapq of {name}");

        let oa = match ours.tags.get(htsjdk_bam::tag::Tag::new(b"OA")) {
            Some(htsjdk_bam::tag::TagValue::Str(value)) => value.clone(),
            other => panic!("{name} has no OA: {other:?}"),
        };
        assert_eq!(oa, row[8], "OA of {name}");
        // The DM tag's value is the empty string: its presence is the message.
        assert_eq!(row[9], "", "DM of {name}");
        assert!(tool::is_distant_mate(ours), "DM of {name}");
        // The NM was moved into the OA, so it is gone from the record. The dump prints an absent
        // attribute as the word null.
        assert_eq!(row[10], "null", "NM of {name}");
        assert!(ours.tags.get(htsjdk_bam::tag::Tag::new(b"NM")).is_none());
    }
}

/// One tag, one missing field, two spellings, on the same read of the same file.
#[test]
fn a_missing_nm_is_spelled_differently_by_the_other_tool() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-distantmates-oa-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let source = source(&dir);
    let header = source.header().clone();
    let all =
        gatk_tools::read_walker::traverse(&source, &[], &|_: &BamRecord| true).expect("the reads");

    // r4 is the read with no NM that both tools reach: distant on its own contig, and on chr1, so
    // AddOriginalAlignmentTags' interval keeps it.
    let read = all
        .iter()
        .find(|read| read.read_name == "r4")
        .expect("the fixture lost the read with no NM");

    let ours = tool::original_alignment_value(read, &header);
    let theirs = add_oa::original_alignment_value(read, &header);
    assert_eq!(ours, "chr1,200,+,10M,60,;");
    assert_eq!(theirs, "chr1,200,+,10M,60,null;");

    // And both are in the golden, off the two outputs, rather than only here.
    let oa_of = |tool_name: &str, label: &str| -> String {
        rows(&text, "reads")
            .into_iter()
            .find(|row| row[0] == tool_name && row[1] == label && row[2] == "r4")
            .map(|row| row[8].to_string())
            .expect("a row for r4")
    };
    assert_eq!(oa_of("PrintDistantMates", "all"), ours);
    assert_eq!(oa_of("AddOriginalAlignmentTags", "chr1"), theirs);
}

/// `undoDistantMateAlterations` over every read of the fixture, against the reference's verdict.
///
/// The reference compares the whole SAM text, tag block included, and says `same` for every read
/// it altered. The port compares the records, which is the same claim in the terms it has.
#[test]
fn the_round_trip_returns_the_record_it_started_from() {
    let text = golden();
    let dir =
        std::env::temp_dir().join(format!("gatk-rs-distantmates-undo-{}", std::process::id()));
    install_fixtures(&text, &dir);

    let source = source(&dir);
    let header = source.header().clone();
    let all =
        gatk_tools::read_walker::traverse(&source, &[], &|_: &BamRecord| true).expect("the reads");

    let verdicts = rows(&text, "roundtrip");
    assert_eq!(verdicts.len(), all.len(), "one verdict per fixture read");

    let mut compared = 0;
    for row in &verdicts {
        let (name, verdict) = (row[0], row[1]);
        let read = all
            .iter()
            .find(|read| read.read_name == name)
            .expect("a read for every verdict");
        if verdict == "unaltered" {
            // The unpaired read: the transform reads a mate contig it does not have, so the
            // reference never offered it one.
            assert!(read.mate_reference_index < 0 || read.flags & 0x1 == 0);
            continue;
        }
        assert_eq!(verdict, "same", "the reference's verdict for {name}");
        let moved = tool::do_distant_mate_alterations(read, &header);
        let back = tool::undo_distant_mate_alterations(&moved, &header).expect("it recovers");
        assert_eq!(&back, read, "the round trip of {name}");
        compared += 1;
    }
    assert_eq!(compared, 7, "seven of the eight fixture reads are altered");

    // And the tag it moves comes back in its old place, which is what makes the SAM texts equal
    // rather than merely equivalent.
    let isdistant = rows(&text, "isdistant");
    assert_eq!(isdistant.len(), compared);
    for row in &isdistant {
        assert_eq!((row[1], row[2], row[3]), ("false", "true", "false"));
    }
}

/// What `undo` refuses, and the one refusal that is not its own.
#[test]
fn an_oa_it_cannot_use_is_refused_in_two_different_ways() {
    let text = golden();
    let header = {
        let dir = std::env::temp_dir().join(format!(
            "gatk-rs-distantmates-refuse-{}",
            std::process::id()
        ));
        install_fixtures(&text, &dir);
        source(&dir).header().clone()
    };

    let refusals = pairs(&text, "undo");
    // The first row is not an OA at all: a read with no OA comes back as the same object.
    assert_eq!(refusals[0], ("nooa", "same object"));

    let mut user_exceptions = 0;
    let mut deferred = 0;
    for (oa, expected) in refusals.iter().skip(1) {
        let mut read = BamRecord {
            read_name: "oa".to_string(),
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            ..BamRecord::default()
        };
        read.tags.insert(
            htsjdk_bam::tag::Tag::new(b"OA"),
            htsjdk_bam::tag::TagValue::Str((*oa).to_string()),
        );
        let error = tool::undo_distant_mate_alterations(&read, &header)
            .expect_err("every one of these is refused");
        let (class, message) = expected.split_once(':').expect("a class and a message");
        assert_eq!(class, error.class(), "the class for {oa}");
        assert_eq!(message.trim(), error.message(), "the message for {oa}");
        match error {
            tool::UndoError::Unrecoverable(_) => user_exceptions += 1,
            // The reference does not raise this one here at all: it keeps the name and resolves it
            // lazily, so the failure surfaces out of whatever first asks for the index.
            tool::UndoError::UnknownContig(_) => deferred += 1,
        }
    }
    assert_eq!(user_exceptions, 3, "three reach the tool's own catch");
    assert_eq!(deferred, 1, "and one escapes it");
}
