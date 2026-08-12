//! Conformance for `BaseRecalibrationEngine` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/BaseRecalibrationEngineDump.java`. BQSR's counting
//! pass: the other half of the cycle `ApplyBQSR` closes, and the last brick before
//! `BaseRecalibrator`.
//!
//! # What this suite is for
//!
//!  * **BAQ is off by default**, so the ordinary run uses a flat array and skips the model;
//!  * **an indel is marked at a different base on each strand**, and clamped away at the ends;
//!  * **the fractional error array spreads each error over a block** that reaches one base before
//!    the uncertain run and includes the base that closes it;
//!  * **a known site is skipped by read offset**, through a conversion whose deletion case steps
//!    back one base;
//!  * **the read group table is collapsed, not counted**;
//!  * **and the table is rounded with the ulp added before the rounding**.

use gatk_corpus as corpus;
use gatk_engine::base_recalibration_engine::{
    calculate_fractional_error_array, calculate_is_snp_or_indel, round_to_n_decimal_places,
    BaseRecalibrationEngine, EngineArguments,
};
use gatk_engine::interval::SimpleInterval;
use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/base_recalibration_engine.txt.gz"),
    )
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

fn constant(text: &str, name: &str) -> String {
    rows(text, "const")
        .into_iter()
        .find(|row| row[0] == name)
        .unwrap_or_else(|| panic!("the golden has no constant {name}"))[1]
        .to_string()
}

fn ints(values: &[i32]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn parse_cigar(text: &str) -> Vec<CigarElement> {
    let mut elements = Vec::new();
    let mut length = 0u32;
    for c in text.chars() {
        if c.is_ascii_digit() {
            length = length * 10 + c.to_digit(10).unwrap();
            continue;
        }
        let op = match c {
            'M' => Op::M,
            'I' => Op::I,
            'D' => Op::D,
            'N' => Op::N,
            'S' => Op::S,
            'H' => Op::H,
            'P' => Op::P,
            '=' => Op::Eq,
            'X' => Op::X,
            other => panic!("no cigar operator {other}"),
        };
        elements.push(CigarElement { length, op });
        length = 0;
    }
    elements
}

/// The header the dump built: one contig, one read group with a platform unit.
fn header(reference: &str) -> SamHeader {
    let mut header = SamHeader::new();
    header
        .sequences
        .push(SequenceRecord::new("chr1", reference.len() as i32));
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "s1");
    group.attributes.set("PL", "ILLUMINA");
    group.attributes.set("PU", "unit-rg1");
    header.read_groups.push(group);
    header
}

/// The reads the dump built, with the quality gradient it gave them.
fn read(name: &str, cigar: &str, start: i32, bases: &str, flags: u16) -> BamRecord {
    let qualities: Vec<u8> = (0..bases.len()).map(|i| (2 + i * 4) as u8).collect();
    let mut record = BamRecord {
        read_name: name.to_string(),
        flags,
        reference_index: 0,
        alignment_start: start,
        cigar: Cigar::new(parse_cigar(cigar)),
        read_bases: bases.as_bytes().to_vec(),
        base_qualities: qualities,
        mapping_quality: 60,
        ..BamRecord::default()
    };
    record
        .tags
        .insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
    record
}

/// The reference bases the read covers, which is what `queryAndPrefetch` hands the counting loop.
fn reference_for(reference: &str, record: &BamRecord) -> Vec<u8> {
    let start = gatk_engine::read_utils::start(record) as usize - 1;
    let end = gatk_engine::read_utils::end(record) as usize;
    reference.as_bytes()[start..end.min(reference.len())].to_vec()
}

#[test]
fn the_constants_are_the_references() {
    let text = golden();
    let arguments = EngineArguments::default();
    assert_eq!(
        constant(&text, "enableBAQ"),
        arguments.enable_baq.to_string()
    );
    assert_eq!(
        constant(&text, "PRESERVE_QSCORES_LESS_THAN"),
        arguments.preserve_qscores_less_than.to_string()
    );
    assert_eq!(
        constant(&text, "defaultBaseQualities"),
        arguments.default_base_qualities.to_string()
    );
    assert_eq!(
        constant(&text, "useOriginalBaseQualities"),
        arguments.use_original_base_qualities.to_string()
    );
}

/// Every event the counting loop marks, over the cigars that move the mark.
#[test]
fn every_marked_event_is_the_reference() {
    let text = golden();
    let reference = constant(&text, "reference");
    let cases: Vec<(&str, &str, i32, &str, u16)> = vec![
        ("exact", "10M", 1, "ACGTACGTAC", 0),
        ("one-mismatch", "10M", 1, "ACGTTCGTAC", 0),
        ("one-mismatch-reverse", "10M", 1, "ACGTTCGTAC", 16),
        ("deletion", "4M2D6M", 1, "ACGTACGTAC", 0),
        ("deletion-reverse", "4M2D6M", 1, "ACGTACGTAC", 16),
        ("insertion", "4M2I4M", 1, "ACGTACGTAC", 0),
        ("insertion-reverse", "4M2I4M", 1, "ACGTACGTAC", 16),
        ("leading-deletion", "1D9M", 1, "ACGTACGTAC", 0),
        ("trailing-deletion", "9M1D", 1, "ACGTACGTAC", 0),
        ("leading-insertion", "1I9M", 1, "ACGTACGTAC", 0),
        ("trailing-insertion", "9M1I", 1, "ACGTACGTAC", 0),
        ("soft-clipped", "3S7M", 1, "ACGTACGTAC", 0),
        ("skipped", "4M2N6M", 1, "ACGTACGTAC", 0),
        ("n-in-read", "10M", 1, "ACGTNCGTAC", 0),
        ("many-mismatches", "10M", 13, "ACGTACGTAC", 0),
    ];

    for (label, cigar, start, bases, flags) in cases {
        let record = read(label, cigar, start, bases, flags);
        let bases_for_read = reference_for(&reference, &record);
        let length = record.read_bases.len();
        let mut snp = vec![0i32; length];
        let mut ins = vec![0i32; length];
        let mut del = vec![0i32; length];
        let events =
            calculate_is_snp_or_indel(&record, &bases_for_read, &mut snp, &mut ins, &mut del)
                .unwrap();
        let theirs = rows(&text, "events")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no events row {label}"));
        assert_eq!(events.to_string(), theirs[1], "{label}: nErrors");
        assert_eq!(ints(&snp), theirs[2], "{label}: isSNP");
        assert_eq!(ints(&ins), theirs[3], "{label}: isInsertion");
        assert_eq!(ints(&del), theirs[4], "{label}: isDeletion");
    }
}

/// The fractional error array, bit for bit, over every error-and-BAQ pair the dump ran.
#[test]
fn every_fractional_error_array_is_the_reference() {
    let text = golden();
    let errors: Vec<Vec<i32>> = vec![
        vec![0, 0, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0],
        vec![1, 0, 0, 0, 0, 1],
        vec![1, 1, 1, 1, 1, 1],
    ];
    let baqs: Vec<Vec<u8>> = vec![
        vec![64, 64, 64, 64, 64, 64],
        vec![64, 60, 60, 64, 64, 64],
        vec![60, 60, 60, 60, 60, 60],
        vec![64, 64, 64, 64, 64, 60],
        vec![60, 64, 64, 64, 64, 64],
    ];
    let mut compared = 0;
    for (e, error) in errors.iter().enumerate() {
        for (b, baq) in baqs.iter().enumerate() {
            let ours = calculate_fractional_error_array(error, baq).unwrap();
            let text_of: Vec<String> = ours
                .iter()
                .map(|value| format!("{:x}", value.to_bits()))
                .collect();
            let label = format!("e{e}-b{b}");
            let theirs = rows(&text, "fractional")
                .into_iter()
                .find(|row| row[0] == label)
                .unwrap_or_else(|| panic!("no fractional row {label}"))[1]
                .to_string();
            assert_eq!(text_of.join(","), theirs, "{label}");
            compared += 1;
        }
    }
    assert_eq!(compared, 20);

    // The length check, which the reference words exactly.
    let message = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "fractional-length-mismatch")
        .unwrap()[2]
        .to_string();
    assert_eq!(
        calculate_fractional_error_array(&[0; 3], &[64; 4])
            .unwrap_err()
            .message(),
        message
    );
}

/// The rounding, whose ulp is added before the round.
#[test]
fn the_rounding_is_the_reference_bit_for_bit() {
    let text = golden();
    for row in rows(&text, "round") {
        let value: f64 = row[0].parse().unwrap();
        let places: i32 = row[1].parse().unwrap();
        let ours = round_to_n_decimal_places(value, places).unwrap();
        assert_eq!(
            format!("{:x}", ours.to_bits()),
            row[2],
            "roundToNDecimalPlaces({value}, {places})"
        );
    }
    let message = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "round-zero-places")
        .unwrap()[2]
        .to_string();
    assert_eq!(
        round_to_n_decimal_places(1.0, 0).unwrap_err().message(),
        message
    );
}

/// The whole engine, over the same corpus and the same known sites.
#[test]
fn every_counted_datum_is_the_reference() {
    let text = golden();
    let reference = constant(&text, "reference");
    let header = header(&reference);

    let corpus: Vec<BamRecord> = vec![
        read("r0", "10M", 1, "ACGTACGTAC", 0),
        read("r1", "10M", 1, "ACGTTCGTAC", 0),
        read("r2", "4M2D6M", 5, "ACGTACGTAC", 0),
        read("r3", "4M2I4M", 9, "ACGTACGTAC", 16),
        read("r4", "10M", 13, "ACGTACGTAC", 0),
    ];

    let runs: Vec<(&str, Vec<SimpleInterval>, bool)> = vec![
        ("plain", vec![], false),
        (
            "known-site",
            vec![SimpleInterval {
                contig: "chr1".to_string(),
                start: 3,
                end: 5,
            }],
            false,
        ),
        (
            "known-everything",
            vec![SimpleInterval {
                contig: "chr1".to_string(),
                start: 1,
                end: 50,
            }],
            false,
        ),
        ("baq-enabled", vec![], true),
    ];

    let mut compared = 0;
    for (label, known_sites, enable_baq) in runs {
        let arguments = EngineArguments {
            enable_baq,
            ..EngineArguments::default()
        };
        let mut engine = BaseRecalibrationEngine::new(arguments, &header).unwrap();
        for record in &corpus {
            engine
                .process_read(record, &header, reference.as_bytes(), &known_sites)
                .unwrap_or_else(|error| {
                    panic!("{label}/{}: {}", record.read_name, error.message())
                });
        }
        let expected_reads = rows(&text, "reads")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("no reads row {label}"))[1]
            .to_string();
        assert_eq!(
            engine.num_reads_processed().to_string(),
            expected_reads,
            "{label}: reads processed"
        );
        engine.finalize_data().unwrap();

        let expected: Vec<Vec<&str>> = rows(&text, "datum")
            .into_iter()
            .filter(|row| row[0] == label)
            .collect();
        let mut ours: Vec<Vec<String>> = Vec::new();
        for (index, table) in engine.tables.all_tables.iter().enumerate() {
            for (keys, datum) in table.all_leaves() {
                ours.push(vec![
                    index.to_string(),
                    ints(&keys),
                    datum.borrow().num_observations().to_string(),
                    format!("{:.2}", datum.borrow().num_mismatches()),
                    format!("{:x}", datum.borrow().reported_quality().to_bits()),
                ]);
            }
        }
        assert_eq!(ours.len(), expected.len(), "{label}: datum count");
        for (ours, theirs) in ours.iter().zip(&expected) {
            assert_eq!(ours[0], theirs[1], "{label}: table index");
            assert_eq!(ours[1], theirs[2], "{label}: keys");
            assert_eq!(ours[2], theirs[3], "{label}: observations");
            assert_eq!(ours[3], theirs[4], "{label}: errors");
            assert_eq!(ours[4], theirs[5], "{label}: reported quality");
            compared += 1;
        }
    }
    println!("base-recalibration-engine: {compared} datums compared");
}
