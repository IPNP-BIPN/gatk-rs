//! Conformance for `ClipReads` against GATK 4.6.2.0, compared as **bytes**.
//!
//! Golden from `tools/readfilter-conformance/ClipReadsDump.java`. The output BAMs and their indexes
//! travel in full, base64, as the rest of this archetype's do, and so does the fixture with its
//! index. The statistics file travels beside them as escaped text, which makes this the first tool
//! here with a second output to compare, and the first whose second output is Java text formatting
//! rather than htsjdk bytes.
//!
//! # What this suite is for
//!
//! The ninth whole tool of the archetype. `ReadClipper` has its own suite; what is here is what the
//! tool builds and hands to it:
//!
//!  * **the quality clipper reads the read in machine-cycle order.** The loop counts down but the
//!    index is flipped for a reverse-strand read, so the walk goes up the array on one strand and
//!    down it on the other, and the op is flipped to match. The fixture's two quality reads carry
//!    the bad qualities at opposite ends of the array, so a port that walked the array backwards on
//!    both would find no clip point at all on the reverse one;
//!  * **the representation decides whether the writer sorts.** Only `WRITE_NS`, `WRITE_NS_Q0S` and
//!    `WRITE_Q0S` are written presorted; the other three get a sorting writer, and the golden's
//!    read order under `SOFTCLIP_BASES` and `HARDCLIP_BASES` is not its read order under `WRITE_NS`;
//!  * **`HARDCLIP_BASES` and `REVERT_SOFTCLIPPED_BASES` revert first**, so `3S7M` at 6 is `10M` at
//!    3 before any op is built;
//!  * **a sequence clips a read as often as it matches**, case-insensitively, against the reverse
//!    complement on a reverse-strand read, and the per-sequence counts print in ASCII order of the
//!    argument as typed;
//!  * **the adapter clipper does not flip for strand**, counts `xf` rather than `xf - 1`, and is
//!    the one place the tool writes a tag;
//!  * **`--read` drops what it does not name**, so naming nothing gives an empty BAM and a
//!    statistics file whose percentages are `NaN`.
//!
//! The command line lands in the `@PG` record's `CL`, so it is read out of the golden and handed to
//! the port rather than reconstructed: it carries the paths of the run that produced it.

use gatk_corpus as corpus;
use gatk_engine::clipping::ClippingRepresentation;
use gatk_engine::reads::ReadsDataSource;
use gatk_readfilter::with_header;
use gatk_tools::clip_reads::{self as tool, ClipArguments};
use gatk_tools::sam_output::Options;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/clip_reads.txt.gz"),
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
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// A row of one kind that carries a single field.
fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| {
            line.strip_prefix(kind)
                .and_then(|rest| rest.strip_prefix('\t'))
        })
        .unwrap_or_else(|| panic!("the golden lost its {kind} row"))
}

fn of_run<'a>(text: &'a str, kind: &str, label: &str) -> Vec<Vec<&'a str>> {
    rows(text, kind)
        .into_iter()
        .filter(|row| row[0] == label)
        .collect()
}

/// What each labelled run was given. A label is a configuration and the row carries nothing to
/// derive it from, so it is written here beside the dump that produced it.
fn arguments(label: &str, clip_fasta: &str) -> ClipArguments {
    let quality = ClipArguments {
        q_trimming_threshold: 10,
        ..ClipArguments::default()
    };
    match label {
        "qt" => quality,
        "cycles" => ClipArguments {
            cycles_to_clip: Some("1-3,8-12".to_string()),
            ..ClipArguments::default()
        },
        "seq" => ClipArguments {
            clip_sequences: vec!["GGGGG".to_string(), "acgt".to_string()],
            ..ClipArguments::default()
        },
        "seqfile" => ClipArguments {
            clip_sequence_file: tool::parse_clip_sequence_file(clip_fasta),
            ..ClipArguments::default()
        },
        "combo" => ClipArguments {
            cycles_to_clip: Some("1-2".to_string()),
            clip_sequences: vec!["GGGGG".to_string()],
            ..quality
        },
        "q0s" => ClipArguments {
            clipping_representation: ClippingRepresentation::WriteQ0s,
            ..quality
        },
        "nsq0s" => ClipArguments {
            clipping_representation: ClippingRepresentation::WriteNsQ0s,
            ..quality
        },
        "soft" => ClipArguments {
            clipping_representation: ClippingRepresentation::SoftclipBases,
            ..quality
        },
        "hard" => ClipArguments {
            clipping_representation: ClippingRepresentation::HardclipBases,
            ..quality
        },
        "revert" => ClipArguments {
            clipping_representation: ClippingRepresentation::RevertSoftclippedBases,
            ..quality
        },
        "adapter" => ClipArguments {
            clip_adapter: true,
            ..ClipArguments::default()
        },
        "onlyread" => ClipArguments {
            only_do_read: Some("r0".to_string()),
            ..quality
        },
        "noread" => ClipArguments {
            only_do_read: Some("nosuchread".to_string()),
            ..quality
        },
        "minlen" => ClipArguments {
            min_read_length: 8,
            ..quality
        },
        // The quality clipper still runs, at its default threshold of -1, and clips nothing.
        "noclip" => ClipArguments::default(),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The fixture, written out so the port can open it.
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
        std::fs::write(
            dir.join(format!("{}.bai", row[0])),
            corpus::decode_base64(row[1]),
        )
        .expect("the fixture index");
    }
}

#[test]
fn every_output_file_is_byte_identical() {
    let text = golden();
    let dir = std::env::temp_dir().join(format!("gatk-rs-clipreads-{}", std::process::id()));
    install(&text, &dir);
    let clip_fasta = unescape(field(&text, "clipfasta"));

    let outputs = rows(&text, "output");
    let indexes = rows(&text, "index");
    let statistics = rows(&text, "stats");
    assert_eq!(outputs.len(), 15, "fifteen runs, none refused");

    let mut compared = 0usize;
    for row in &outputs {
        let (label, expected_base64) = (row[0], row[1]);
        let source = ReadsDataSource::open(&dir.join("plain.bam"), &dir.join("plain.bai"))
            .expect("the fixture opens");
        let header = source.header().clone();

        // This tool does not override its read filters, so its default is GATKTool's.
        let filter = move |read: &BamRecord| with_header::wellformed(read, &header);

        let command_line = of_run(&text, "commandline", label)
            .first()
            .map(|row| row.get(1).copied().unwrap_or(""))
            .unwrap_or("");
        let options = Options {
            command_line,
            ..Options::default()
        };

        let (ours, our_index, our_stats) =
            tool::clip_reads(&source, &options, &arguments(label, &clip_fasta), &filter)
                .expect("the source reads")
                .expect("no label is refused");

        let expected = corpus::decode_base64(expected_base64);
        assert_eq!(ours.len(), expected.len(), "{label}: output length differs");
        if ours != expected {
            let at = ours
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            panic!("{label}: first byte difference at offset {at}");
        }

        let expected_index = indexes
            .iter()
            .find(|index| index[0] == label)
            .map(|index| index[1])
            .expect("an index row for every output");
        match (our_index, expected_index) {
            (None, "absent") => {}
            (Some(_), "absent") => panic!("{label}: the reference wrote no index and the port did"),
            (None, _) => panic!("{label}: the reference wrote an index and the port did not"),
            (Some(ours), expected) => {
                assert_eq!(ours, corpus::decode_base64(expected), "{label}: the .bai");
            }
        }

        let expected_stats = statistics
            .iter()
            .find(|stats| stats[0] == label)
            .map(|stats| unescape(stats[1]))
            .expect("a statistics row for every output");
        assert_eq!(our_stats, expected_stats, "{label}: the statistics file");

        compared += 1;
    }

    assert_eq!(compared, 15);
    println!("clip-reads: {compared} runs, output, index and statistics byte-identical");
}

/// The scan is in machine-cycle order, which is the opposite of array order on one strand.
///
/// A port that walked the array backwards on both strands would clip the forward read correctly and
/// find no clip point at all on the reverse one, which is a healthier-looking file and a wrong one.
#[test]
fn the_quality_scan_reaches_opposite_ends_of_the_two_strands() {
    let text = golden();
    let clipped = of_run(&text, "reads", "qt");
    assert_eq!(clipped.len(), 7);

    let by_name = |name: &str| -> Vec<String> {
        clipped
            .iter()
            .find(|row| row[1] == name)
            .map(|row| row.iter().map(|field| field.to_string()).collect())
            .unwrap_or_else(|| panic!("the golden lost {name}"))
    };

    // Both went in with `IIIII#####`-shaped qualities in array order; r1 is the reverse-strand one.
    assert_eq!(by_name("r0")[6], "ACGTANNNNN");
    assert_eq!(by_name("r1")[6], "NNNNNACGTA");
    // And nothing else in the fixture has a quality low enough to reach.
    assert_eq!(by_name("r5")[6], "GGGGGGGGGG");
}

/// The three representations that can move a read get a sorting writer, and the order changes.
#[test]
fn the_representation_decides_whether_the_output_is_in_traversal_order() {
    let text = golden();
    let order = |label: &str| -> Vec<String> {
        of_run(&text, "reads", label)
            .iter()
            .map(|row| row[1].to_string())
            .collect()
    };

    let traversal = order("qt");
    assert_eq!(traversal, ["r0", "r1", "r2", "r3", "r4", "r5", "r6"]);
    // A front soft clip moves r1 from 5 to 10, behind r2.
    assert_eq!(order("soft"), ["r0", "r2", "r1", "r3", "r4", "r5", "r6"]);
    // And a reverted soft clip moves r2 from 6 back to 3, in front of r1, which the hard clip has
    // moved to 10.
    assert_eq!(order("hard"), ["r0", "r2", "r1", "r3", "r4", "r5", "r6"]);
    assert_eq!(order("revert"), ["r0", "r2", "r1", "r3", "r4", "r5", "r6"]);
}

/// What each clipper did to a read, which is the finding rather than the byte count.
#[test]
fn each_clipper_writes_over_what_the_reference_says_it_does() {
    let text = golden();
    let bases = |label: &str, name: &str| -> String {
        of_run(&text, "reads", label)
            .iter()
            .find(|row| row[1] == name)
            .map(|row| row[6].to_string())
            .unwrap_or_else(|| panic!("the golden lost {label}/{name}"))
    };

    // 1-3 and 8-12 on a ten-base read: the second range is clamped, not dropped.
    assert_eq!(bases("cycles", "r0"), "NNNTAGGNNN");
    // The same on a five-base read: the second range never starts.
    assert_eq!(bases("cycles", "r6"), "NNNTA");
    // `GGGGG` matches twice on a ten-G read, so the find() loop goes round twice.
    assert_eq!(bases("seq", "r5"), "NNNNNNNNNN");
    // On a reverse read the pattern is reverse-complemented: `GGGGG` becomes `CCCCC` and misses,
    // while `acgt` is its own reverse complement and matches twice, case-insensitively.
    assert_eq!(bases("seq", "r3"), "NNNNNNNNAC");
    // `GGTACC` out of the -XF file, which is a palindrome, on the soft-clipped read.
    assert_eq!(bases("seqfile", "r2"), "TTTNNNNNNA");
    // XF=3, XT=8 on a reverse-strand read, used as written rather than flipped.
    assert_eq!(bases("adapter", "r3"), "NNGTACGNNN");
    // Both adapter tags zero clips the whole read.
    assert_eq!(bases("adapter", "r4"), "NNNNNNNNNN");
    // And with no clipping argument the tool runs and changes nothing.
    assert_eq!(bases("noclip", "r0"), "ACGTAGGTAC");
}

/// The statistics file is the tool's second output, and its text is Java's rather than htsjdk's.
#[test]
fn the_statistics_file_is_java_text_formatting() {
    let text = golden();
    let stats = |label: &str| -> String {
        of_run(&text, "stats", label)
            .first()
            .map(|row| unescape(row[1]))
            .unwrap_or_else(|| panic!("the golden lost the {label} statistics"))
    };

    // Nothing examined, so both percentages are a division by zero, which Java prints as NaN.
    let empty = stats("noread");
    assert!(
        empty.contains("Number of examined reads              0\n"),
        "{empty}"
    );
    assert!(
        empty.contains("Percent of clipped reads              NaN\n"),
        "{empty}"
    );
    assert!(
        empty.contains("Percent of clipped bases              NaN\n"),
        "{empty}"
    );

    // The per-sequence rows are a TreeMap, so upper case sorts before lower case, and the `%8d`
    // width is the reference's.
    let sequences = stats("seq");
    let rows: Vec<&str> = sequences
        .lines()
        .filter(|line| line.contains("clip sites matching"))
        .collect();
    assert_eq!(
        rows,
        vec![
            "        10 clip sites matching GGGGG",
            "        20 clip sites matching acgt",
        ]
    );

    // The adapter row exists only under -CA, and counts `xf` rather than `xf - 1`: three bases from
    // r3's 3' op, three counted for its 5' op, and ten for r4's whole-read clip.
    assert!(stats("adapter").contains("Number of adapter clipped bases       16\n"));
    assert!(!stats("qt").contains("adapter"));
}

/// `--read` drops what it does not name, and the length filter drops what is too short.
#[test]
fn two_arguments_decide_which_reads_are_written_at_all() {
    let text = golden();
    let names = |label: &str| -> Vec<String> {
        of_run(&text, "reads", label)
            .iter()
            .map(|row| row[1].to_string())
            .collect()
    };

    assert_eq!(
        names("onlyread"),
        ["r0"],
        "the other six are gone, not passed through"
    );
    assert!(names("noread").is_empty());
    // r6 is five bases long and every other read is ten.
    assert_eq!(names("minlen"), ["r0", "r1", "r2", "r3", "r4", "r5"]);
}
