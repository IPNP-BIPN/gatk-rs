//! Conformance for `ReadClipper` against GATK 4.6.2.0.
//!
//! The golden carries the **whole clipped read** for every entry point over the read-filter
//! corpus: start, flags, cigar, bases, qualities, mapping quality and tags. Comparing one field
//! would hide the others, and clipping changes all of them at once.
//!
//! `E` is the reference throwing. It is common here: clipping the middle of a read, soft-clipping
//! an unmapped one, and asking for a reference coordinate the read does not cover are exceptions
//! rather than reads.

use gatk_corpus as corpus;
use gatk_engine::clipping::{self, ClippingRepresentation};
use gatk_engine::read_utils;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_clipping.txt.gz"),
    )
}

/// The read as the dump renders it.
fn render(read: &BamRecord) -> String {
    let quals: Vec<String> = read.base_qualities.iter().map(|q| q.to_string()).collect();
    let tags: Vec<String> = read
        .tags
        .iter()
        .map(|(tag, value)| {
            let text = match value {
                htsjdk_bam::tag::TagValue::Str(text) => text.clone(),
                htsjdk_bam::tag::TagValue::ByteArray { values, .. } => values
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                htsjdk_bam::tag::TagValue::Int(value) => value.to_string(),
                htsjdk_bam::tag::TagValue::Float(value) => value.to_string(),
                other => format!("{other:?}"),
            };
            format!("{tag}={text}")
        })
        .collect();
    format!(
        "{}|{}|{}|{}|{}|{}|{}",
        read_utils::start(read),
        read.flags,
        read.cigar.to_text(),
        String::from_utf8_lossy(&read.read_bases),
        quals.join(","),
        read.mapping_quality,
        tags.join(";")
    )
}

fn run(
    operation: &str,
    read: &BamRecord,
    header: &SamHeader,
) -> Result<BamRecord, gatk_engine::clipping::ClipError> {
    let header = Some(header);
    let (name, argument) = operation.split_once('@').unwrap_or((operation, ""));
    match name {
        "leftTail" => clipping::hard_clip_by_reference_coordinates_left_tail(
            read,
            header,
            argument.parse().unwrap(),
        ),
        "rightTail" => clipping::hard_clip_by_reference_coordinates_right_tail(
            read,
            header,
            argument.parse().unwrap(),
        ),
        "toRegion" => {
            let (start, stop) = argument.split_once(',').unwrap();
            clipping::hard_clip_to_region(
                read,
                header,
                start.parse().unwrap(),
                stop.parse().unwrap(),
            )
        }
        "revertSoftClipped" => clipping::revert_soft_clipped_bases(read, header),
        _ => {
            let low_qual = argument.parse().unwrap();
            let representation = match name {
                "hardClipLowQual" => ClippingRepresentation::HardclipBases,
                "softClipLowQual" => ClippingRepresentation::SoftclipBases,
                "writeNsLowQual" => ClippingRepresentation::WriteNs,
                "writeQ0sLowQual" => ClippingRepresentation::WriteQ0s,
                "writeNsQ0sLowQual" => ClippingRepresentation::WriteNsQ0s,
                other => panic!("{other} is in the golden but not ported"),
            };
            clipping::clip_low_qual_ends(read, header, low_qual, representation)
        }
    }
}

#[test]
fn every_clipped_read_is_the_read_the_reference_produced() {
    let text = golden();
    let records = corpus::records(&text);
    let header = corpus::header(&text);

    let mut compared = 0;
    for line in text.lines() {
        let parts: Vec<&str> = line.splitn(4, '\t').collect();
        if parts[0] != "clipped" {
            continue;
        }
        let index: usize = parts[1].parse().unwrap();
        let operation = parts[2];
        let expected = parts[3];
        let record = &records[index];

        let ours = run(operation, record, &header).map_or_else(|_| "E".to_string(), |r| render(&r));
        assert_eq!(ours, expected, "{}: {operation}", records[index].read_name);
        compared += 1;
    }
    assert!(compared > 0, "the golden carries no clipped reads");
    println!("{compared} clipped reads, all identical");
}
