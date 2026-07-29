//! Conformance for `LocusIteratorByState` against GATK 4.6.2.0.
//!
//! Twelve runs, three fixtures crossed with the four `(includeDeletions, includeNs)` settings, and
//! one row per yielded pileup.
//!
//! Two families of rows carry the suite. `indels-*` shows the two exclusions being independent: at
//! `includeDeletions=false` the deleted positions drop to one element, while the `N` read is
//! absent from every setting where `includeNs` is false. `adaptor-*` shows the per-base test: `p1`
//! is missing from chr1:101 to chr1:105 and present from 106, because its adaptor boundary is 105
//! and the reverse-strand comparison is `position <= boundary`.

use gatk_corpus as corpus;
use gatk_engine::locus_iterator::{self, LocusIteratorOptions};
use gatk_engine::read_states::{DownsamplingInfo, ReadStateManager};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const CONTIG_LENGTH: i32 = 400;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/locus_iterator.txt.gz"),
    )
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", CONTIG_LENGTH));
    for (id, sample) in [("rg1", "sampleA"), ("rg2", "sampleB")] {
        let mut group = ReadGroup::new(id);
        group.attributes.set("SM", sample);
        group.attributes.set("PL", "ILLUMINA");
        header.read_groups.push(group);
    }
    header
}

/// One fixture read. A non-zero fragment length makes it a reverse-strand read whose mate starts
/// five bases in, which is what puts its adaptor boundary inside its own span.
fn read(name: &str, group: &str, cigar: &str, start: i32, fragment: i32) -> BamRecord {
    let cigar = htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar");
    let length = cigar.read_length() as usize;
    let mut tags = htsjdk_bam::tag::Tags::new();
    tags.insert(Tag::new(b"RG"), TagValue::Str(group.to_string()));
    let mut record = BamRecord {
        read_name: name.to_string(),
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        read_bases: (0..length).map(|i| b"ACGT"[i % 4]).collect(),
        base_qualities: vec![30; length],
        cigar,
        tags,
        ..Default::default()
    };
    if fragment != 0 {
        // paired | proper pair | reverse strand, with a mapped forward mate.
        record.flags = 0x1 | 0x2 | 0x10;
        record.mate_reference_index = 0;
        record.mate_alignment_start = start + 5;
        record.inferred_insert_size = -fragment;
    }
    record
}

fn reads_for(label: &str) -> Vec<BamRecord> {
    match label.split('-').next().expect("a fixture name") {
        "plain" => vec![
            read("a1", "rg1", "10M", 101, 0),
            read("b1", "rg2", "10M", 104, 0),
            read("a2", "rg1", "10M", 108, 0),
        ],
        "indels" => vec![
            read("d1", "rg1", "3M4D3M", 101, 0),
            read("n1", "rg1", "3M4N3M", 101, 0),
            read("m1", "rg2", "10M", 101, 0),
        ],
        "adaptor" => vec![
            read("p1", "rg1", "10M", 101, 60),
            read("q1", "rg1", "10M", 101, 0),
        ],
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The label's suffix encodes the two flags: upper case is on.
fn options_for(label: &str) -> LocusIteratorOptions {
    let suffix = label.rsplit('-').next().expect("a suffix");
    let mut chars = suffix.chars();
    LocusIteratorOptions {
        include_deletions: chars.next() == Some('D'),
        include_ns: chars.next() == Some('N'),
    }
}

#[test]
fn every_pileup_matches_the_reference() {
    let text = golden();
    let header = header();

    let mut labels: Vec<String> = Vec::new();
    let mut rows: std::collections::HashMap<String, Vec<String>> = Default::default();
    let mut counts: std::collections::HashMap<String, usize> = Default::default();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("ctx\t") {
            let mut parts = rest.splitn(3, '\t');
            let label = parts.next().expect("a label").to_string();
            let _index = parts.next();
            rows.entry(label)
                .or_default()
                .push(parts.next().unwrap_or("").to_string());
        } else if let Some(rest) = line.strip_prefix("count\t") {
            let (label, count) = rest.split_once('\t').expect("a label and a count");
            labels.push(label.to_string());
            counts.insert(label.to_string(), count.parse().expect("a number"));
        }
    }
    assert!(!labels.is_empty(), "the golden carries no count rows");

    let mut compared = 0;
    for label in &labels {
        let reads = reads_for(label);
        let samples = vec![Some("sampleA".to_string()), Some("sampleB".to_string())];
        let states = ReadStateManager::new(samples.clone(), DownsamplingInfo::NONE)
            .expect("no downsampling");
        let contexts =
            locus_iterator::contexts(&reads, samples, &header, options_for(label), states)
                .unwrap_or_else(|e| panic!("{label}: {e:?}"));

        assert_eq!(contexts.len(), counts[label], "{label}: context count");
        let expected = rows.get(label).cloned().unwrap_or_default();
        for (index, context) in contexts.iter().enumerate() {
            let bases: String = context.pileup.bases().iter().map(|b| *b as char).collect();
            let names: Vec<String> = context
                .pileup
                .elements
                .iter()
                .map(|e| e.read.read_name.clone())
                .collect();
            let ours = format!(
                "{}:{}\t{}\t{}\t{}",
                context.contig,
                context.position,
                context.pileup.size(),
                bases,
                names.join(",")
            );
            assert_eq!(ours, expected[index], "{label} context {index}");
            compared += 1;
        }
    }

    println!(
        "{compared} pileups over {} runs, all identical",
        labels.len()
    );
}
