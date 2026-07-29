//! Conformance for `ReadPileup` against GATK 4.6.2.0.
//!
//! Three pileups and twenty-four overlap-fix pairs.
//!
//! The `mixed` pileup is the one that carries the suite's point: its bases are `ACDN*:` and its
//! base counts are `A=2`, because the wildcard `*` maps to `A` in `BaseUtils.baseIndexMap` while
//! `N` and `:` both fall out through the `-1` check and the deletion falls out through an explicit
//! test. Three exclusions, three different mechanisms, one row.

use gatk_corpus as corpus;
use gatk_engine::read_pileup;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const CONTIG_LENGTH: i32 = 200;
const LOCUS: i32 = 105;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_pileup.txt.gz"),
    )
}

/// The header the harness built: two read groups with samples, one with none.
fn header() -> SamHeader {
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", CONTIG_LENGTH));
    for (id, sample) in [
        ("rg1", Some("sampleA")),
        ("rg2", Some("sampleB")),
        ("rg3", None),
    ] {
        let mut group = ReadGroup::new(id);
        if let Some(sample) = sample {
            group.attributes.set("SM", sample);
        }
        group.attributes.set("PL", "ILLUMINA");
        header.read_groups.push(group);
    }
    header
}

/// One read of the fixture: qualities are `20 + offset`, as the harness sets them.
fn read(name: &str, group: &str, cigar: &str, start: i32, bases: &str) -> BamRecord {
    let mut tags = htsjdk_bam::tag::Tags::new();
    tags.insert(Tag::new(b"RG"), TagValue::Str(group.to_string()));
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        read_bases: bases.as_bytes().to_vec(),
        base_qualities: (0..bases.len()).map(|i| 20 + i as u8).collect(),
        cigar: htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar"),
        tags,
        ..Default::default()
    }
}

/// The reads of each labelled pileup, in the order the harness listed them.
fn reads_for(label: &str) -> Vec<BamRecord> {
    match label {
        "mixed" => vec![
            read("r1", "rg1", "10M", 101, "ACGTACGTAC"),
            read("r2", "rg1", "10M", 101, "CCCCCCCCCC"),
            read("r3", "rg2", "2M5D3M", 101, "ACGTA"),
            read("r4", "rg2", "10M", 101, "NNNNNNNNNN"),
            read("r5", "rg1", "10M", 101, "AAAA*AAAAA"),
            read("r6", "rg2", "10M", 101, "::::::::::"),
        ],
        "staggered" => vec![
            read("s3", "rg1", "10M", 100, "ACGTACGTAC"),
            read("s1", "rg1", "10M", 96, "ACGTACGTAC"),
            read("s2", "rg2", "10M", 100, "TTTTTTTTTT"),
            read("s4", "rg1", "10M", 98, "GGGGGGGGGG"),
        ],
        "nosample" => vec![
            read("n1", "rg1", "10M", 101, "ACGTACGTAC"),
            read("n2", "rg3", "10M", 101, "TTTTTTTTTT"),
        ],
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_pileup_answers_what_the_reference_answers() {
    let text = golden();
    let header = header();

    let mut pileups = 0;
    let mut overlaps = 0;
    for line in text.lines() {
        let Some((kind, rest)) = line.split_once('\t') else {
            continue;
        };
        match kind {
            "pileup" => {
                let mut parts = rest.split('\t');
                let label = parts.next().expect("a label");
                let size: usize = parts.next().expect("a size").parse().expect("a number");
                let bases = parts.next().expect("bases");
                let quals = parts.next().expect("quals");
                let counts = parts.next().expect("counts");
                let offsets = parts.next().expect("offsets");
                let string = parts.next().expect("a pileup string");

                let reads = reads_for(label);
                // No filter is disabled in the harness, so both default predicates pass every
                // fixture read; they are still applied, because the constructor applies them.
                let pileup =
                    read_pileup::pileup_from_reads("chr1", LOCUS, &reads, |_| true, |_| true);
                assert_eq!(pileup.size(), size, "{label}: size");
                assert_eq!(
                    String::from_utf8_lossy(&pileup.bases()),
                    bases,
                    "{label}: bases"
                );
                assert_eq!(
                    pileup
                        .quals()
                        .iter()
                        .map(|q| q.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                    quals,
                    "{label}: quals"
                );
                let ours = pileup.base_counts();
                assert_eq!(
                    format!("{},{},{},{}", ours[0], ours[1], ours[2], ours[3]),
                    counts,
                    "{label}: base counts"
                );
                assert_eq!(
                    format!(
                        "[{}]",
                        pileup
                            .offsets()
                            .iter()
                            .map(|o| o.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    ),
                    offsets,
                    "{label}: offsets"
                );
                assert_eq!(pileup.pileup_string('A'), string, "{label}: pileup string");
                pileups += 1;
            }
            "sorted" => {
                let (label, expected) = rest.split_once('\t').expect("a label and an order");
                let reads = reads_for(label);
                let pileup =
                    read_pileup::pileup_from_reads("chr1", LOCUS, &reads, |_| true, |_| true);
                let ours = pileup
                    .sorted()
                    .iter()
                    .map(|e| format!("{}@{}", e.read.read_name, e.read.alignment_start))
                    .collect::<Vec<_>>()
                    .join("|");
                assert_eq!(ours, expected, "{label}: sortedIterator order");
            }
            "sample" => {
                let (label, expected) = rest.split_once('\t').expect("a label and samples");
                let reads = reads_for(label);
                let pileup =
                    read_pileup::pileup_from_reads("chr1", LOCUS, &reads, |_| true, |_| true);
                // The harness sorted the names, with a missing sample rendered as "null".
                let mut names: Vec<String> = pileup
                    .samples(&header)
                    .into_iter()
                    .map(|s| s.unwrap_or_else(|| "null".to_string()))
                    .collect();
                names.sort();
                let ours = names
                    .iter()
                    .map(|name| {
                        let wanted = if name == "null" {
                            None
                        } else {
                            Some(name.as_str())
                        };
                        format!(
                            "{name}={}",
                            pileup.pileup_for_sample(wanted, &header).size()
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("|");
                assert_eq!(ours, expected, "{label}: samples");
            }
            "split" => {
                let (label, expected) = rest.split_once('\t').expect("a label and an outcome");
                let reads = reads_for(label);
                let pileup =
                    read_pileup::pileup_from_reads("chr1", LOCUS, &reads, |_| true, |_| true);
                let ours = match pileup.split_by_sample(&header, None) {
                    Ok(split) => format!("ok:{}", split.len()),
                    Err(_) => {
                        "E:org.broadinstitute.hellbender.exceptions.UserException$ReadMissingReadGroup"
                            .to_string()
                    }
                };
                assert_eq!(ours, expected, "{label}: splitBySample");
            }
            "overlap" => {
                let mut parts = rest.split('\t');
                let same = parts.next().expect("same or differ") == "same";
                let first: u8 = parts.next().expect("q1").parse().expect("a number");
                let second: u8 = parts.next().expect("q2").parse().expect("a number");
                let expected = parts.next().expect("the new pair");
                let (a, b) = read_pileup::fix_pair_overlapping_qualities(
                    b'A',
                    first,
                    if same { b'A' } else { b'C' },
                    second,
                );
                assert_eq!(format!("{a},{b}"), expected, "overlap {first},{second}");
                overlaps += 1;
            }
            _ => {}
        }
    }

    assert!(pileups > 0, "the golden carries no pileup rows");
    println!("{pileups} pileups and {overlaps} overlap fixes, all identical");
}
