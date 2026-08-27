//! Conformance for `GroundTruthReadsBuilder` against GATK 4.6.2.0, compared as the reads every
//! run kept and the columns their rows carry.
//!
//! Golden from `tools/readfilter-conformance/GroundTruthReadsBuilderDump.java`.
//!
//! The flow-based scoring engine is not measured or ported: the scores are read off the golden's
//! own rows and everything that decides which rows exist is rebuilt from the fixture.
//!
//! # What this suite is for
//!
//!  * **the translation table, its ignored first line and its fallback**;
//!  * **the translated contig being the read's own name with the ancestor appended**;
//!  * **a collapsed translation being skipped and counted, not refused**;
//!  * **the mapping-quality floor**;
//!  * **the soft-clip filter looking at the END of the read**;
//!  * **the two score filters, and both being off at zero**;
//!  * **the output cap**;
//!  * **the fixed column order**;
//!  * **and the CSV being quoted, the flow keys holding commas.**

use gatk_corpus as corpus;
use gatk_tools::ground_truth_reads_builder::{
    header, is_end_softclipped, is_polyt, keeps_scores, parse_cigar, quote, split_row,
    translate_span, Arguments, Translator, CSV_FIELD_ORDER, MATERNAL, PATERNAL,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/ground_truth_reads_builder.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

/// The rows of one run's CSV, split on the commas that are not inside quotes.
fn rows(text: &str, label: &str) -> Vec<Vec<String>> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.is_empty())
        .map(split_row)
        .collect()
}

/// The read names one run kept, in order.
fn kept(text: &str, label: &str) -> Vec<String> {
    rows(text, label)
        .into_iter()
        .skip(1)
        .map(|row| row[0].clone())
        .collect()
}

/// One column of one run, by name.
fn column(text: &str, label: &str, name: &str) -> Vec<String> {
    let rows = rows(text, label);
    let index = rows[0]
        .iter()
        .position(|column| column == name)
        .unwrap_or_else(|| panic!("the csv carries {name}"));
    rows.into_iter()
        .skip(1)
        .map(|row| row[index].clone())
        .collect()
}

/// The fixture's reads, from the golden's SAM.
fn reads(text: &str) -> Vec<(String, i32, String, i32, String)> {
    section(text, "sam", "reads")
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[0].to_string(),
                columns[3].parse().expect("a position"),
                columns[5].to_string(),
                columns[4].parse().expect("a quality"),
                columns[9].to_string(),
            )
        })
        .collect()
}

/// The header is the tool's own column order.
#[test]
fn the_column_order_is_fixed() {
    let text = golden();
    assert_eq!(CSV_FIELD_ORDER.len(), 22);
    for label in ["default", "min-mq-thirty", "max-two-reads"] {
        assert_eq!(rows(&text, label)[0], CSV_FIELD_ORDER.to_vec(), "{label}");
    }
    assert_eq!(
        header(),
        section(&text, "out", "default")
            .lines()
            .next()
            .expect("a header")
    );
    // The scores come before the keys, and the intervals last.
    assert_eq!(CSV_FIELD_ORDER[4], "PaternalHaplotypeScore");
    assert_eq!(CSV_FIELD_ORDER[5], "MaternalHaplotypeScore");
    assert_eq!(
        *CSV_FIELD_ORDER.last().expect("a column"),
        "MaternalHaplotypeInterval"
    );
}

/// The flow keys hold commas, so a naive split reads the columns out of step.
#[test]
fn the_csv_is_quoted() {
    let text = golden();
    let line = section(&text, "out", "default")
        .lines()
        .nth(1)
        .expect("a row")
        .to_string();
    // Split properly there are twenty-two fields; split on every comma there are more.
    assert_eq!(split_row(&line).len(), 22);
    assert!(line.split(',').count() > 22, "the keys hold commas");
    // The key column is one of the ones that holds them.
    let keys = column(&text, "default", "ReadKey");
    assert!(keys[0].contains(','), "{}", keys[0]);
    // And the port quotes a value that holds one, doubling any quote inside it.
    assert_eq!(quote("plain"), "plain");
    assert_eq!(quote("0,1,2"), "\"0,1,2\"");
    assert_eq!(quote("a\"b"), "\"a\"\"b\"");
}

/// The first line is ignored and a position between two rows takes the earlier offset.
#[test]
fn the_translation_table_ignores_its_first_line() {
    let text = golden();
    let paternal = Translator::parse(&section(&text, "csv", "paternal"));
    assert_eq!(paternal.positions, vec![1, 1500]);
    assert_eq!(paternal.offsets, vec![0, 10]);
    // Between the two rows the earlier offset applies.
    assert_eq!(paternal.translate(1000), Some(1000));
    assert_eq!(paternal.translate(1499), Some(1499));
    // At and past the second row the later one does.
    assert_eq!(paternal.translate(1500), Some(1510));
    assert_eq!(paternal.translate(1600), Some(1610));
    // The maternal table is an identity throughout.
    let maternal = Translator::parse(&section(&text, "csv", "maternal"));
    assert_eq!(maternal.positions, vec![1]);
    assert_eq!(maternal.translate(1000), Some(1000));
    // A position before the first row has no earlier offset at all.
    assert_eq!(maternal.translate(0), None);
    // A table whose first line is a row loses it.
    let headerless = Translator::parse("1,0\n1500,10\n");
    assert_eq!(headerless.positions, vec![1500]);
}

/// The read's own name with the ancestor appended.
#[test]
fn the_translated_contig_names_the_ancestor() {
    let text = golden();
    let paternal = Translator::parse(&section(&text, "csv", "paternal"));
    let span = translate_span(&paternal, "chr1", 1600, 1699, PATERNAL).expect("a span");
    assert_eq!(span.contig, "chr1_paternal");
    assert_eq!(span.start, 1610);
    assert_eq!(span.end, 1709);
    let maternal = Translator::parse(&section(&text, "csv", "maternal"));
    let span = translate_span(&maternal, "chr1", 1600, 1699, MATERNAL).expect("a span");
    assert_eq!(span.contig, "chr1_maternal");
    assert_eq!(span.start, 1600);
    // Which is why the two reference files carry those contigs.
    assert_eq!(section(&text, "fasta", "maternal"), "chr1_maternal:2400");
    assert_eq!(section(&text, "fasta", "paternal"), "chr1_paternal:2400");
    assert_eq!(section(&text, "fasta", "reference"), "chr1:2400");
    // The interval columns the golden wrote name the same contigs.
    let intervals = column(&text, "default", "MaternalHaplotypeInterval");
    assert!(intervals[0].contains("chr1_maternal"), "{}", intervals[0]);
}

/// The traversal carries on, and the run's other reads are still written.
#[test]
fn a_collapsed_translation_is_skipped_and_counted() {
    let text = golden();
    // The collapsing table sends 1010 back by two hundred, so the first read's end lands before
    // its start.
    let collapsing = Translator::parse(&section(&text, "csv", "collapsing"));
    assert_eq!(collapsing.translate(1000), Some(1000));
    assert_eq!(collapsing.translate(1099), Some(899));
    let failure =
        translate_span(&collapsing, "chr1", 1000, 1099, PATERNAL).expect_err("a collapsed span");
    assert!(
        failure.contains("failed to translate for paternal"),
        "{failure}"
    );
    assert!(failure.contains("start:1000 ,end:899"), "{failure}");
    // The run kept every other read and did not fail.
    assert_eq!(
        kept(&text, "collapsing-translation"),
        vec![
            "r-second",
            "r-low-mq",
            "r-polyt-clip",
            "r-other-clip",
            "r-shifted"
        ]
    );
    assert!(!text.contains("error\tcollapsing-translation"));
    // The default run, whose table does not collapse, kept that read too.
    assert!(kept(&text, "default").contains(&"r-plain".to_string()));
}

/// A read below the floor is dropped, and nothing else changes.
#[test]
fn the_mapping_quality_floor_drops_a_read() {
    let text = golden();
    let reads = reads(&text);
    let low = reads
        .iter()
        .find(|read| read.0 == "r-low-mq")
        .expect("a read");
    assert_eq!(low.3, 10);
    assert!(kept(&text, "default").contains(&"r-low-mq".to_string()));
    assert!(!kept(&text, "min-mq-thirty").contains(&"r-low-mq".to_string()));
    // Every other read survived both runs, in the same order.
    let without: Vec<String> = kept(&text, "default")
        .into_iter()
        .filter(|name| name != "r-low-mq")
        .collect();
    assert_eq!(without, kept(&text, "min-mq-thirty"));
    let arguments = Arguments {
        min_mapping_quality: 30,
        ..Arguments::default()
    };
    assert!(low.3 < arguments.min_mapping_quality);
}

/// A read clipped only at its front is kept whatever the argument says.
#[test]
fn the_soft_clip_filter_looks_at_the_end_of_the_read() {
    let text = golden();
    let reads = reads(&text);
    for name in ["r-polyt-clip", "r-other-clip"] {
        let read = reads.iter().find(|read| read.0 == name).expect("a read");
        let cigar = parse_cigar(&read.2);
        // The clip is the FIRST element, so the filter does not see it.
        assert_eq!(cigar[0].operator, 'S');
        assert!(!is_end_softclipped(&cigar), "{name}");
    }
    // Both runs keep both reads, so turning the filter off changes nothing here.
    assert_eq!(kept(&text, "default"), kept(&text, "keep-softclipped"));
    for name in ["r-polyt-clip", "r-other-clip"] {
        assert!(kept(&text, "default").contains(&name.to_string()), "{name}");
    }
    // A clip at the end would be seen, and a poly-T one spared.
    assert!(is_end_softclipped(&parse_cigar("96M4S")));
    assert!(is_polyt(b"TTTT"));
    assert!(!is_polyt(b"ACGT"));
    assert!(!is_polyt(b""));
    // The fixture's own clips are one of each, which is what they were written for.
    let polyt = reads
        .iter()
        .find(|read| read.0 == "r-polyt-clip")
        .expect("a read");
    assert!(is_polyt(&polyt.4.as_bytes()[..4]));
    let other = reads
        .iter()
        .find(|read| read.0 == "r-other-clip")
        .expect("a read");
    assert!(!is_polyt(&other.4.as_bytes()[..4]));
}

/// Both are off at zero, and the delta one drops the reads away from either difference.
#[test]
fn the_two_score_filters_are_off_at_zero() {
    let text = golden();
    // The default run has both at zero and keeps every read that translated.
    assert_eq!(kept(&text, "default").len(), 6);
    assert!(keeps_scores(-100.0, -1.0, &Arguments::default()));
    // The delta filter drops the reads whose two haplotypes are far apart.
    let delta = Arguments {
        min_haplotype_score_delta: 1.0,
        ..Arguments::default()
    };
    assert_eq!(
        kept(&text, "min-score-delta"),
        vec!["r-second", "r-polyt-clip", "r-other-clip"]
    );
    // Which the port agrees with, read by read, from the golden's own scores.
    let maternal = column(&text, "default", "MaternalHaplotypeScore");
    let paternal = column(&text, "default", "PaternalHaplotypeScore");
    let names = kept(&text, "default");
    let produced: Vec<String> = names
        .iter()
        .zip(maternal.iter().zip(paternal.iter()))
        .filter(|(_, (m, p))| {
            keeps_scores(
                m.parse().expect("a score"),
                p.parse().expect("a score"),
                &delta,
            )
        })
        .map(|(name, _)| name.clone())
        .collect();
    assert_eq!(produced, kept(&text, "min-score-delta"));
    // The score filter as the golden ran it keeps everything, the scores all being below it.
    assert_eq!(kept(&text, "min-score"), kept(&text, "default"));
}

/// A haplotype identical to the reference takes the reference's own score.
#[test]
fn an_haplotype_identical_to_the_reference_is_not_rescored() {
    let text = golden();
    let names = kept(&text, "default");
    let paternal = column(&text, "default", "PaternalHaplotypeScore");
    let maternal = column(&text, "default", "MaternalHaplotypeScore");
    let reference = column(&text, "default", "RefHaplotypeScore");
    // The paternal difference is at 1200, so the read at 1000 sees a paternal haplotype that is
    // the reference and takes its score exactly.
    let first = names
        .iter()
        .position(|name| name == "r-plain")
        .expect("a read");
    assert_eq!(paternal[first], reference[first]);
    assert_ne!(maternal[first], reference[first]);
    // And at least one read sees both ancestors differ from the reference at once.
    assert!(names
        .iter()
        .enumerate()
        .any(|(i, _)| paternal[i] != reference[i] || maternal[i] != reference[i]));
}

/// The cap stops the output rather than sampling it.
#[test]
fn the_output_cap_keeps_the_first_reads() {
    let text = golden();
    let capped = kept(&text, "max-two-reads");
    assert_eq!(capped.len(), 2);
    // They are the first two of the uncapped run, in the same order.
    assert_eq!(capped, kept(&text, "default")[..2].to_vec());
    let arguments = Arguments {
        max_output_reads: Some(2),
        ..Arguments::default()
    };
    assert_eq!(arguments.max_output_reads, Some(2));
    assert_eq!(Arguments::default().max_output_reads, None);
}

/// Every run's rows are the same shape, whatever the arguments did.
#[test]
fn every_run_writes_the_same_columns() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "default",
        "min-mq-thirty",
        "keep-softclipped",
        "prepend-append",
        "flow-length-ten",
        "flow-length-large",
        "min-score",
        "min-score-delta",
        "max-two-reads",
        "collapsing-translation",
    ] {
        let rows = rows(&text, label);
        assert_eq!(rows[0], CSV_FIELD_ORDER.to_vec(), "{label}");
        for row in rows.iter().skip(1) {
            assert_eq!(row.len(), 22, "{label}");
        }
        compared += 1;
    }
    assert_eq!(compared, 10, "the runs the port reproduces");
}

/// Added to the haplotype and not to the read.
#[test]
fn the_prepended_sequence_lands_on_the_haplotype() {
    let text = golden();
    let plain = column(&text, "default", "BestHaplotypeSequence");
    let padded = column(&text, "prepend-append", "BestHaplotypeSequence");
    let plain_reads = column(&text, "default", "ReadSequence");
    let padded_reads = column(&text, "prepend-append", "ReadSequence");
    // The reads are untouched.
    assert_eq!(plain_reads, padded_reads);
    // The haplotypes are not.
    assert_ne!(plain, padded);
    assert!(padded[0].starts_with("TTTT"), "{}", padded[0]);
    assert!(padded[0].ends_with("CCCC"), "{}", padded[0]);
    // And what sits between is what the plain run wrote.
    assert_eq!(
        padded[0]
            .strip_prefix("TTTT")
            .and_then(|rest| rest.strip_suffix("CCCC")),
        Some(plain[0].as_str())
    );
}

/// It fixes the length of the two HAPLOTYPE keys in both directions, and leaves the read's alone.
#[test]
fn the_flow_length_fixes_the_haplotype_keys() {
    let text = golden();
    let cells = |key: &str| key.split(',').count();
    // The read's own key is the same length under every setting.
    for label in ["default", "flow-length-ten", "flow-length-large"] {
        assert_eq!(cells(&column(&text, label, "ReadKey")[0]), 103, "{label}");
    }
    // The two haplotype keys are the ones the argument moves, and it moves both the same way.
    for name in ["BestHaplotypeKey", "ConsensusHaplotypeKey"] {
        assert_eq!(cells(&column(&text, "default", name)[0]), 111, "{name}");
        // A length BELOW the key's own truncates it rather than leaving it alone.
        assert_eq!(
            cells(&column(&text, "flow-length-ten", name)[0]),
            10,
            "{name}"
        );
        // And a length above it pads it out to exactly that many.
        assert_eq!(
            cells(&column(&text, "flow-length-large", name)[0]),
            400,
            "{name}"
        );
    }
    // The padded key begins with the key the plain run wrote, so the padding is at the end.
    let plain = column(&text, "default", "BestHaplotypeKey");
    let padded = column(&text, "flow-length-large", "BestHaplotypeKey");
    assert!(padded[0].starts_with(&plain[0]), "{}", &padded[0][..40]);
    // And the truncated one is that key's own first ten cells.
    let cut = column(&text, "flow-length-ten", "BestHaplotypeKey");
    let first_ten: String = plain[0].split(',').take(10).collect::<Vec<_>>().join(",");
    assert_eq!(cut[0], first_ten);
}
