//! Conformance for `FlowFeatureMapper` against GATK 4.6.2.0, compared as the records of every run.
//!
//! Golden from `tools/readfilter-conformance/FlowFeatureMapperDump.java`.
//!
//! The flow matrix's score is not measured or ported: the scores are read off the golden's own
//! records and everything else is rebuilt from the fixture's BAM and reference.
//!
//! # What this suite is for
//!
//!  * **a feature being a base surrounded by bases that match**;
//!  * **the surround being two arguments, the second defaulting to the first**;
//!  * **an element shorter than the surround being skipped whole**;
//!  * **an `N` in the reference not being a mismatch**;
//!  * **X_FC1 being the mismatch count and X_FC2 the feature count**;
//!  * **X_INDEX being the offset in the whole read and X_LENGTH the unclipped length**;
//!  * **an interval selecting reads and not features**;
//!  * **the two score bounds**;
//!  * **a duplicate read being dropped**;
//!  * **and --copy-attr carrying a type and a description of its own.**

use gatk_corpus as corpus;
use gatk_tools::flow_feature_mapper::{
    edit_distance, features, info_column, is_surrounded, keeps_read, keeps_score, mismatch_count,
    parse_cigar, Arguments, CopyAttribute, Read, Surround,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/flow_feature_mapper.txt.gz"),
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

/// The fixture's reference: `TGCA` repeated, with `N` from 1500 to 1509.
fn reference_base(position: i32) -> u8 {
    if (1500..1510).contains(&position) {
        return b'N';
    }
    b"TGCA"[((position - 1) % 4) as usize]
}

fn reference_bases(start: i32, length: i32) -> Vec<u8> {
    (0..length).map(|i| reference_base(start + i)).collect()
}

/// The fixture's reads, read back from the golden's SAM text.
fn reads(text: &str) -> Vec<Read> {
    section(text, "sam", "reads")
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Read {
                name: columns[0].to_string(),
                contig: columns[2].to_string(),
                start: columns[3].parse().expect("a position"),
                cigar: parse_cigar(columns[5]),
                bases: columns[9].as_bytes().to_vec(),
                flags: columns[1].parse().expect("a flag"),
                mapping_quality: columns[4].parse().expect("a quality"),
            }
        })
        .collect()
}

/// The reference the walker hands a read: its own span, soft clips excluded.
fn reference_for(read: &Read) -> Vec<u8> {
    let length: i32 = read
        .cigar
        .iter()
        .filter(|element| element.consumes_reference_bases())
        .map(|element| element.length)
        .sum();
    reference_bases(read.start, length)
}

/// One record the golden wrote: its position, alleles and INFO fields.
#[derive(Debug, Clone, PartialEq)]
struct Record {
    position: i32,
    reference: String,
    alternate: String,
    info: Vec<(String, String)>,
}

impl Record {
    fn get(&self, key: &str) -> Option<&str> {
        self.info
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }
}

fn records(text: &str, label: &str) -> Vec<Record> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Record {
                position: columns[1].parse().expect("a position"),
                reference: columns[3].to_string(),
                alternate: columns[4].to_string(),
                info: columns[7]
                    .split(';')
                    .map(|part| match part.split_once('=') {
                        Some((key, value)) => (key.to_string(), value.to_string()),
                        None => (part.to_string(), String::new()),
                    })
                    .collect(),
            }
        })
        .collect()
}

/// The features one run should have produced, as `(position, ref, alt, read name, index)`.
fn produced(text: &str, arguments: &Arguments) -> Vec<(i32, String, String, String, i32)> {
    let mut rows = Vec::new();
    for read in reads(text) {
        if !keeps_read(&read, arguments) {
            continue;
        }
        let reference = reference_for(&read);
        for feature in features(&read, &reference, arguments.surround) {
            rows.push((
                feature.start,
                (feature.reference_base as char).to_string(),
                (feature.read_base as char).to_string(),
                read.name.clone(),
                feature.index,
            ));
        }
    }
    // The records come out in reference order, the reads being held in a queue behind the walk.
    rows.sort_by_key(|row| row.0);
    rows
}

fn measured(text: &str, label: &str) -> Vec<(i32, String, String, String, i32)> {
    records(text, label)
        .into_iter()
        .map(|record| {
            (
                record.position,
                record.reference.clone(),
                record.alternate.clone(),
                record.get("X_RN").expect("a read name").to_string(),
                record
                    .get("X_INDEX")
                    .expect("an index")
                    .parse()
                    .expect("a number"),
            )
        })
        .collect()
}

/// label, the arguments it was run with.
fn runs() -> Vec<(&'static str, Arguments)> {
    vec![
        ("default", Arguments::default()),
        (
            "identical-two",
            Arguments {
                surround: Surround::new(2, 0),
                ..Arguments::default()
            },
        ),
        (
            "identical-three-after-one",
            Arguments {
                surround: Surround::new(3, 1),
                ..Arguments::default()
            },
        ),
        (
            "include-duplicates",
            Arguments {
                include_duplicate_reads: true,
                ..Arguments::default()
            },
        ),
        ("copy-attr", Arguments::default()),
    ]
}

/// Every run's features, in the order and with the fields the golden wrote them.
#[test]
fn every_run_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, arguments) in runs() {
        assert_eq!(
            produced(&text, &arguments),
            measured(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "the runs the port reproduces");
}

/// A mismatch on the first base of its element has nothing before it.
#[test]
fn a_feature_must_be_surrounded_by_matching_bases() {
    let text = golden();
    let reads = reads(&text);
    let edge = reads.iter().find(|r| r.name == "r-edge").expect("a read");
    let reference = reference_for(edge);
    // Its two mismatches are at offsets 0 and 20, and only the second becomes a feature.
    assert_ne!(edge.bases[0], reference[0]);
    assert_ne!(edge.bases[20], reference[20]);
    let found = features(edge, &reference, Surround::default());
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].index, 20);
    assert!(!is_surrounded(
        &edge.bases,
        &reference,
        0,
        0,
        Surround::default()
    ));
    assert!(is_surrounded(
        &edge.bases,
        &reference,
        20,
        20,
        Surround::default()
    ));
    // Which is what the golden wrote: one record for that read, at 1220.
    let from_edge: Vec<i32> = records(&text, "default")
        .into_iter()
        .filter(|record| record.get("X_RN") == Some("r-edge"))
        .map(|record| record.position)
        .collect();
    assert_eq!(from_edge, vec![1220]);
}

/// Two mismatches two apart survive a surround of one and go at two.
#[test]
fn the_surround_is_two_arguments() {
    let text = golden();
    let at = |label: &str, name: &str| -> Vec<i32> {
        records(&text, label)
            .into_iter()
            .filter(|record| record.get("X_RN") == Some(name))
            .map(|record| record.position)
            .collect()
    };
    // The base between them matches, so both are surrounded at one.
    assert_eq!(at("default", "r-adjacent"), vec![1110, 1112]);
    assert_eq!(at("identical-two", "r-adjacent"), Vec::<i32>::new());
    // Three before and one after keeps the first and drops the second, the second needing three
    // matching bases before it and having a mismatch two back.
    assert_eq!(at("identical-three-after-one", "r-adjacent"), vec![1110]);
    // The second argument defaults to the first rather than to zero.
    assert_eq!(
        Surround::new(2, 0),
        Surround {
            before: 2,
            after: 2
        }
    );
    assert_eq!(
        Surround::new(3, 1),
        Surround {
            before: 3,
            after: 1
        }
    );
    assert_eq!(
        Surround::default(),
        Surround {
            before: 1,
            after: 1
        }
    );
}

/// An element shorter than the surround plus one is skipped whole.
#[test]
fn a_short_cigar_element_is_skipped_whole() {
    let text = golden();
    let reads = reads(&text);
    let short = reads
        .iter()
        .find(|r| r.name == "r-short-elements")
        .expect("a read");
    assert_eq!(short.cigar, parse_cigar("3M1I3M1D33M"));
    assert_eq!(Surround::default().minimum_element_length(), 3);
    assert_eq!(Surround::new(2, 0).minimum_element_length(), 5);
    let at = |label: &str| -> Vec<i32> {
        records(&text, label)
            .into_iter()
            .filter(|record| record.get("X_RN") == Some("r-short-elements"))
            .map(|record| record.position)
            .collect()
    };
    // The mismatch at offset 1 sits in the first 3M, which a surround of one looks inside and a
    // surround of two does not.
    assert_eq!(at("default"), vec![1301, 1320]);
    assert_eq!(at("identical-two"), vec![1320]);
}

/// It does count towards the edit distance, though.
#[test]
fn an_n_in_the_reference_is_not_a_mismatch() {
    let text = golden();
    let reads = reads(&text);
    let over_n = reads.iter().find(|r| r.name == "r-over-n").expect("a read");
    let reference = reference_for(over_n);
    // Its two changed bases are at offsets 10 and 25, and the second sits under the `N` run.
    assert_eq!(reference[10], reference_base(1490));
    assert_eq!(reference[25], b'N');
    let found = features(over_n, &reference, Surround::default());
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].start, 1490);
    // The mismatch count agrees, and the edit distance does not.
    assert_eq!(mismatch_count(over_n, &reference), 1);
    let distance = edit_distance(over_n.aligned_bases(), &reference);
    assert_eq!(distance, 11, "ten Ns and one mismatch");
    let record = records(&text, "default")
        .into_iter()
        .find(|record| record.get("X_RN") == Some("r-over-n"))
        .expect("its record");
    assert_eq!(record.get("X_EDIST"), Some("11"));
    assert_eq!(record.get("X_FC1"), Some("1"));
}

/// The two differ whenever a mismatch failed the surround test.
#[test]
fn fc1_is_the_mismatch_count_and_fc2_the_feature_count() {
    let text = golden();
    let reads = reads(&text);
    for (name, mismatches, count) in [
        ("r-three", 3, 3),
        ("r-adjacent", 2, 2),
        // Two mismatches, one of them on the element's first base: one feature.
        ("r-edge", 2, 1),
        ("r-over-n", 1, 1),
    ] {
        let read = reads.iter().find(|r| r.name == name).expect("a read");
        let reference = reference_for(read);
        assert_eq!(mismatch_count(read, &reference), mismatches, "{name}");
        assert_eq!(
            features(read, &reference, Surround::default()).len(),
            count,
            "{name}"
        );
        let record = records(&text, "default")
            .into_iter()
            .find(|record| record.get("X_RN") == Some(name))
            .expect("a record");
        assert_eq!(
            record.get("X_FC1"),
            Some(mismatches.to_string().as_str()),
            "{name}"
        );
        assert_eq!(
            record.get("X_FC2"),
            Some(count.to_string().as_str()),
            "{name}"
        );
    }
    // X_FC2 moves with the surround arguments and X_FC1 does not.
    let adjacent = records(&text, "identical-three-after-one")
        .into_iter()
        .find(|record| record.get("X_RN") == Some("r-adjacent"))
        .expect("a record");
    assert_eq!(adjacent.get("X_FC1"), Some("2"));
    assert_eq!(adjacent.get("X_FC2"), Some("1"));
}

/// The soft clip counts in both.
#[test]
fn the_index_and_the_length_include_the_soft_clip() {
    let text = golden();
    let reads = reads(&text);
    let clipped = reads
        .iter()
        .find(|r| r.name == "r-clipped")
        .expect("a read");
    assert_eq!(clipped.cigar, parse_cigar("5S35M"));
    assert_eq!(clipped.bases.len(), 40);
    assert_eq!(clipped.aligned_bases().len(), 35);
    assert_eq!(clipped.unclipped_length(), 40);
    let reference = reference_for(clipped);
    assert_eq!(reference.len(), 35);
    let found = features(clipped, &reference, Surround::default());
    // The indices are into the whole read, so they are five more than the aligned offsets.
    assert_eq!(
        found.iter().map(|f| f.index).collect::<Vec<_>>(),
        vec![10, 20]
    );
    assert_eq!(
        found.iter().map(|f| f.start).collect::<Vec<_>>(),
        vec![1705, 1715]
    );
    for record in records(&text, "default")
        .into_iter()
        .filter(|record| record.get("X_RN") == Some("r-clipped"))
    {
        assert_eq!(record.get("X_LENGTH"), Some("40"));
        assert_eq!(record.get("X_CIGAR"), Some("5S35M"));
    }
}

/// A read that starts inside it contributes every feature it carries.
#[test]
fn an_interval_selects_reads_and_not_features() {
    let text = golden();
    let inside = records(&text, "one-interval");
    // chr1:1000-1100 keeps the two reads starting at 1000 and 1100.
    let names: Vec<String> = inside
        .iter()
        .map(|record| record.get("X_RN").expect("a name").to_string())
        .collect();
    assert!(names
        .iter()
        .all(|name| name == "r-three" || name == "r-adjacent"));
    // And it keeps the features past the interval's end, at 1110 and 1112.
    let positions: Vec<i32> = inside.iter().map(|record| record.position).collect();
    assert_eq!(positions, vec![1010, 1020, 1030, 1110, 1112]);
    assert!(positions.iter().any(|position| *position > 1100));
    // Which is exactly what those two reads contribute to the unrestricted run.
    let whole: Vec<i32> = records(&text, "default")
        .into_iter()
        .filter(|record| matches!(record.get("X_RN"), Some("r-three") | Some("r-adjacent")))
        .map(|record| record.position)
        .collect();
    assert_eq!(positions, whole);
}

/// Each bound cuts the set somewhere.
#[test]
fn the_two_score_bounds_drop_the_records_outside() {
    let text = golden();
    let scores = |label: &str| -> Vec<f64> {
        records(&text, label)
            .into_iter()
            .map(|record| {
                record
                    .get("X_SCORE")
                    .expect("a score")
                    .parse()
                    .expect("a number")
            })
            .collect()
    };
    let whole = scores("default");
    assert_eq!(whole.len(), 11);
    // The fixture produces three distinct scores.
    let mut distinct = whole.clone();
    distinct.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    distinct.dedup();
    assert_eq!(distinct.len(), 3);
    // The upper bound drops the highest and the lower bound drops the lowest.
    let upper = Arguments {
        maximum_score: 6.0,
        ..Arguments::default()
    };
    let lower = Arguments {
        minimum_score: 5.5,
        ..Arguments::default()
    };
    assert_eq!(
        whole.iter().filter(|s| keeps_score(**s, &upper)).count(),
        scores("max-score-6").len()
    );
    assert_eq!(
        whole.iter().filter(|s| keeps_score(**s, &lower)).count(),
        scores("min-score-5-5").len()
    );
    assert!(scores("max-score-6").iter().all(|s| *s <= 6.0));
    assert!(scores("min-score-5-5").iter().all(|s| *s >= 5.5));
    // A bound below every score empties the file.
    assert!(scores("max-score-0").is_empty());
    // A NaN is kept unless it is asked about.
    assert!(keeps_score(f64::NAN, &Arguments::default()));
    assert!(!keeps_score(
        f64::NAN,
        &Arguments {
            exclude_nan_scores: true,
            ..Arguments::default()
        }
    ));
}

/// Unless --include-dup-reads asks for it, and then its records carry the flag.
#[test]
fn a_duplicate_read_is_dropped() {
    let text = golden();
    let reads = reads(&text);
    let duplicate = reads
        .iter()
        .find(|r| r.name == "r-duplicate")
        .expect("a read");
    assert!(duplicate.is_duplicate());
    assert!(!keeps_read(duplicate, &Arguments::default()));
    assert!(keeps_read(
        duplicate,
        &Arguments {
            include_duplicate_reads: true,
            ..Arguments::default()
        }
    ));
    assert!(!section(&text, "out", "default").contains("r-duplicate"));
    let included: Vec<i32> = records(&text, "include-duplicates")
        .into_iter()
        .filter(|record| record.get("X_RN") == Some("r-duplicate"))
        .map(|record| record.position)
        .collect();
    assert_eq!(included, vec![1610, 1620, 1630]);
    // And the flag reaches the record.
    let record = records(&text, "include-duplicates")
        .into_iter()
        .find(|record| record.get("X_RN") == Some("r-duplicate"))
        .expect("a record");
    assert_eq!(record.get("X_FLAGS"), Some("1024"));
}

/// With a type and a description read out of the argument itself.
#[test]
fn copy_attr_carries_its_own_type_and_description() {
    let text = golden();
    let parsed = CopyAttribute::parse("za,Integer,a number");
    assert_eq!(parsed.name, "za");
    assert_eq!(parsed.kind, "Integer");
    assert_eq!(parsed.description, "a number");
    assert_eq!(parsed.key("P_"), "P_za");
    // The header line the golden wrote says the same.
    let line = text
        .lines()
        .find(|line| line.starts_with("header\tcopy-attr\t##INFO=<ID=P_za"))
        .expect("its header line");
    assert!(line.contains("Type=Integer"), "{line}");
    assert!(line.contains("Description=\"a number\""), "{line}");
    // Every record carries the tag, whose value is the read's own.
    for record in records(&text, "copy-attr") {
        let name = record.get("X_RN").expect("a name");
        let read = reads(&text)
            .into_iter()
            .find(|read| read.name == name)
            .expect("a read");
        // The fixture sets `za` to the read's alignment start.
        assert_eq!(
            record.get("P_za"),
            Some(read.start.to_string().as_str()),
            "{name}"
        );
    }
    // A spec with no type is a String described after itself, and a description may hold commas.
    let bare = CopyAttribute::parse("zb");
    assert_eq!(bare.kind, "String");
    assert_eq!(bare.description, "copy-attr: zb");
    let comma = CopyAttribute::parse("zc,Float,one, two, three");
    assert_eq!(comma.description, "one, two, three");
    // The prefix may be empty, which is the default.
    assert_eq!(parsed.key(""), "za");
}

/// The INFO keys are written in sorted order, so a copied tag lands where its name puts it.
#[test]
fn the_info_keys_are_sorted() {
    let text = golden();
    for label in ["default", "copy-attr"] {
        for record in records(&text, label) {
            let keys: Vec<String> = record.info.iter().map(|(key, _)| key.clone()).collect();
            let mut sorted = keys.clone();
            sorted.sort();
            assert_eq!(keys, sorted, "{label}");
        }
    }
    // The prefixed tag sorts before every X_ key, which is where the golden has it.
    let first = records(&text, "copy-attr")[0].info[0].clone();
    assert_eq!(first.0, "P_za");
    assert_eq!(
        info_column(&[
            ("X_RN".to_string(), "r".to_string()),
            ("P_za".to_string(), "1".to_string())
        ]),
        "P_za=1;X_RN=r"
    );
}
