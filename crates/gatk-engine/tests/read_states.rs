//! Conformance for `ReadStateManager` and its per-sample managers against GATK 4.6.2.0.
//!
//! One row per step of a traversal: the total state count, and the read names each sample holds
//! with the genome position each is at.
//!
//! Four rows carry the suite:
//!
//!  * `two-samples` declares `sampleB` before `sampleA` while the header has them the other way,
//!    and every row lists B first. The order is the declaration's, because the reference's map is
//!    a `LinkedHashMap` and says so in capitals;
//!  * `all-clipped` holds three reads and reports two: `5S5I` never enters, because its first
//!    `stepForwardOnGenome` returns null and `addReadsToSample` keeps only non-null ones;
//!  * `deletion-boundary` admits `d2` at position 105 while `d1` is part-way through a `10D`,
//!    which is where `d1` *is* rather than where it started;
//!  * `undeclared-sample` is an `IllegalStateException` rather than a new bucket.

use gatk_corpus as corpus;
use gatk_engine::read_states::{DownsamplingInfo, ReadStateError, ReadStateManager};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use std::collections::VecDeque;

const CONTIG_LENGTH: i32 = 300;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_states.txt.gz"),
    )
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for name in ["chr1", "chr2"] {
        header
            .sequences
            .push(SequenceRecord::new(name, CONTIG_LENGTH));
    }
    for (id, sample) in [("rg1", "sampleA"), ("rg2", "sampleB")] {
        let mut group = ReadGroup::new(id);
        group.attributes.set("SM", sample);
        group.attributes.set("PL", "ILLUMINA");
        header.read_groups.push(group);
    }
    header
}

fn read(name: &str, group: Option<&str>, cigar: &str, start: i32) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar");
    let length = cigar.read_length() as usize;
    let mut tags = htsjdk_bam::tag::Tags::new();
    if let Some(group) = group {
        tags.insert(Tag::new(b"RG"), TagValue::Str(group.to_string()));
    }
    BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        read_bases: (0..length).map(|i| b"ACGT"[i % 4]).collect(),
        base_qualities: vec![30; length],
        cigar,
        tags,
        ..Default::default()
    }
}

/// The samples each labelled run declared, in the order it declared them, and its reads.
fn configuration(label: &str) -> (Vec<Option<String>>, Vec<BamRecord>) {
    let samples = |names: &[Option<&str>]| -> Vec<Option<String>> {
        names.iter().map(|n| n.map(|s| s.to_string())).collect()
    };
    match label {
        "two-samples" => (
            samples(&[Some("sampleB"), Some("sampleA")]),
            vec![
                read("a1", Some("rg1"), "10M", 101),
                read("b1", Some("rg2"), "10M", 101),
                read("a2", Some("rg1"), "10M", 103),
                read("b2", Some("rg2"), "10M", 105),
            ],
        ),
        "all-clipped" => (
            samples(&[Some("sampleA")]),
            vec![
                read("c1", Some("rg1"), "10M", 101),
                read("c2", Some("rg1"), "5S5I", 101),
                read("c3", Some("rg1"), "10M", 101),
            ],
        ),
        "deletion-boundary" => (
            samples(&[Some("sampleA")]),
            vec![
                read("d1", Some("rg1"), "2M10D8M", 101),
                read("d2", Some("rg1"), "10M", 105),
                read("d3", Some("rg1"), "10M", 113),
            ],
        ),
        "undeclared-sample" => (
            samples(&[Some("sampleA")]),
            vec![
                read("f1", Some("rg1"), "10M", 101),
                read("f2", Some("rg2"), "10M", 101),
            ],
        ),
        "null-sample" => (
            samples(&[Some("sampleA"), None]),
            vec![
                read("g1", Some("rg1"), "10M", 101),
                read("g2", None, "10M", 101),
            ],
        ),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

#[test]
fn every_step_holds_what_the_reference_holds() {
    let text = golden();
    let header = header();

    // The rows of each run, in order.
    let mut labels: Vec<String> = Vec::new();
    let mut steps: std::collections::HashMap<String, Vec<String>> = Default::default();
    let mut errors: std::collections::HashMap<String, String> = Default::default();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("step\t") {
            let mut parts = rest.splitn(3, '\t');
            let label = parts.next().expect("a label").to_string();
            let _index = parts.next();
            if !labels.contains(&label) {
                labels.push(label.clone());
            }
            steps
                .entry(label)
                .or_default()
                .push(parts.next().unwrap_or("").to_string());
        } else if let Some(rest) = line.strip_prefix("error\t") {
            let (label, class) = rest.split_once('\t').expect("a label and a class");
            if !labels.contains(&label.to_string()) {
                labels.push(label.to_string());
            }
            errors.insert(label.to_string(), class.to_string());
        }
    }
    assert!(!labels.is_empty(), "the golden carries no rows");

    let mut compared = 0;
    for label in &labels {
        let (samples, reads) = configuration(label);
        let mut manager = ReadStateManager::new(samples.clone(), DownsamplingInfo::NONE)
            .expect("no downsampling");
        let mut pending: VecDeque<&BamRecord> = reads.iter().collect();

        if let Some(class) = errors.get(label) {
            // The reference refuses at the first collect, when the undeclared sample is submitted.
            let result = manager.collect_pending_reads(&mut pending, &header);
            assert!(
                matches!(result, Err(ReadStateError::UndeclaredSample(_))),
                "{label}: the reference raised {class}, the port gave {result:?}"
            );
            continue;
        }

        let expected = &steps[label];
        for (index, row) in expected.iter().enumerate() {
            manager
                .collect_pending_reads(&mut pending, &header)
                .unwrap_or_else(|e| panic!("{label} step {index}: {e:?}"));

            let (total, contents) = row.split_once('\t').expect("a total and contents");
            assert_eq!(
                manager.size().to_string(),
                total,
                "{label} step {index}: total states"
            );

            let ours = manager
                .samples
                .iter()
                .zip(&manager.by_sample)
                .map(|(sample, per_sample)| {
                    let names = if per_sample.states.is_empty() {
                        "-".to_string()
                    } else {
                        per_sample
                            .states
                            .iter()
                            .map(|s| {
                                format!("{}@{}", s.read.read_name, s.machine.genome_position())
                            })
                            .collect::<Vec<_>>()
                            .join(",")
                    };
                    format!(
                        "{}={names}",
                        sample.clone().unwrap_or_else(|| "null".to_string())
                    )
                })
                .collect::<Vec<_>>()
                .join("|");
            assert_eq!(ours, contents, "{label} step {index}");
            compared += 1;

            manager
                .update_read_states()
                .unwrap_or_else(|e| panic!("{label} step {index}: {e:?}"));
        }
    }

    println!(
        "{compared} traversal steps over {} runs, all identical",
        labels.len()
    );
}
