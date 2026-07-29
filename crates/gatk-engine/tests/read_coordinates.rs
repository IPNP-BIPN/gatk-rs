//! Conformance for `ReadUtils`' coordinate mapping against GATK 4.6.2.0.
//!
//! The golden is produced by `tools/readfilter-conformance/ReadCoordinateDump.java` over the same
//! corpus as the read filters, and probes every reference position from three before each read's
//! start to three past its end, so the not-found answers at both edges are covered rather than
//! assumed.
//!
//! `E` in the golden is the reference throwing, and it is kept apart from `.`, which is the
//! reference answering "absent". Both appear here: a read with no cigar throws when asked for its
//! last insertion offset, a read of `*` bases throws when reverse-complemented, and a read with no
//! qualities throws for a quality it answers a base for.

use gatk_corpus as corpus;
use gatk_engine::read_utils::{self, BaseAt};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/read_coordinates.txt.gz"),
    )
}

fn show(value: Option<impl ToString>) -> String {
    value.map_or_else(|| "E".to_string(), |v| v.to_string())
}

fn show_base(value: BaseAt) -> String {
    match value {
        BaseAt::Absent => ".".to_string(),
        // The dump prints a byte as a signed decimal, which is what String.valueOf((int)(byte) b)
        // gives, so a quality of 30 is "30" and a base 'A' is "65".
        BaseAt::Present(byte) => (byte as i8).to_string(),
        BaseAt::Threw => "E".to_string(),
    }
}

#[test]
fn every_coordinate_maps_the_way_the_reference_maps_it() {
    let text = golden();
    let records = corpus::records(&text);

    let mut soft_rows = 0;
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts[0] != "soft" {
            continue;
        }
        let index: usize = parts[1].parse().unwrap();
        let record = &records[index];
        let name = &record.read_name;
        assert_eq!(
            read_utils::soft_start(record).to_string(),
            parts[2],
            "{name}: soft start"
        );
        assert_eq!(
            read_utils::soft_end(record).to_string(),
            parts[3],
            "{name}: soft end"
        );
        assert_eq!(
            show(read_utils::last_insertion_offset(record)),
            parts[4],
            "{name}: last insertion offset"
        );
        assert_eq!(
            show(read_utils::bases_reverse_complement(record)),
            parts[5],
            "{name}: reverse complement"
        );
        soft_rows += 1;
    }
    assert_eq!(
        soft_rows,
        records.len(),
        "the golden has {soft_rows} soft rows for {} records",
        records.len()
    );

    let mut coord_rows = 0;
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts[0] != "coord" {
            continue;
        }
        let index: usize = parts[1].parse().unwrap();
        let ref_coord: i32 = parts[2].parse().unwrap();
        let record = &records[index];
        let at = format!("{} at {ref_coord}", record.read_name);

        let (read_index, operator) = read_utils::read_index_for_read(record, ref_coord);
        assert_eq!(read_index.to_string(), parts[3], "{at}: read index");
        // CigarOperator's toString is its character, and `=` is what Eq prints.
        let operator = operator.map_or(".".to_string(), |op| (op.to_char() as char).to_string());
        assert_eq!(operator, parts[4], "{at}: operator");

        assert_eq!(
            show_base(read_utils::read_base_at_reference_coordinate(
                record, ref_coord
            )),
            parts[5],
            "{at}: base"
        );
        assert_eq!(
            show_base(read_utils::read_base_quality_at_reference_coordinate(
                record, ref_coord
            )),
            parts[6],
            "{at}: base quality"
        );
        assert_eq!(
            if read_utils::is_inside_read(record, ref_coord) {
                "1"
            } else {
                "0"
            },
            parts[7],
            "{at}: inside read"
        );
        coord_rows += 1;
    }
    assert!(coord_rows > 0, "the golden carries no coordinate rows");
    println!(
        "{coord_rows} coordinates over {} records, all identical",
        records.len()
    );
}
