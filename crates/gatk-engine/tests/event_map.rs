//! Conformance for `EventMap` and `Event` against the oracle.
//!
//! Golden from `tools/readfilter-conformance/EventMapDump.java`. Each `case` row carries the
//! reference, the haplotype, its CIGAR, its alignment start and the MNP distance; this test builds
//! the map from them and renders the `event`, `overlap` and `error` rows that follow. The `direct`
//! rows are `new Event` and `makeCompoundEvents` on the inputs the dump lists, repeated here.

use gatk_corpus as corpus;
use gatk_engine::event_map::{
    event_start_positions, make_compound_events, Event, EventError, EventMap,
};
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use htsjdk_vcf::allele::Allele;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/event_map.txt.gz"),
    )
}

fn error(label: &str, e: &EventError) -> String {
    format!("error\t{label}\t{}: {}", e.class(), e.message())
}

/// The rows one `case` produces, and the haplotype if it was built.
fn run(fields: &[&str]) -> (Vec<String>, Option<Haplotype>) {
    let label = fields[1];
    let reference = fields[2].as_bytes();
    let ref_start: i32 = fields[3].parse().unwrap();
    let cigar = htsjdk_bam::text_parse::parse_cigar(fields[5]).unwrap();
    let alignment_start: i32 = fields[6].parse().unwrap();
    let max_mnp: i32 = fields[7].parse().unwrap();
    let mut out = Vec::new();
    let mut h = match Haplotype::new(fields[4].as_bytes(), false) {
        Ok(h) => h,
        Err(e) => {
            out.push(error(label, &EventError::Allele(e)));
            return (out, None);
        }
    };
    h.set_cigar(&cigar).unwrap();
    h.set_alignment_start_hap_wrt_ref(alignment_start);
    let ref_loc =
        SimpleInterval::new("chr1", ref_start, ref_start + reference.len() as i32 - 1).unwrap();
    let map = match EventMap::from_haplotype(&h, reference, &ref_loc, max_mnp) {
        Ok(map) => map,
        Err(e) => {
            out.push(error(label, &e));
            return (out, None);
        }
    };
    let events: Vec<&Event> = map.events().collect();
    for (i, e) in events.iter().enumerate() {
        out.push(format!(
            "event\t{label}\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.start(),
            e.end(),
            e.ref_allele().base_string(),
            e.alt_allele().base_string(),
            e.is_snp(),
            e.is_indel(),
            e.is_simple_insertion(),
            e.is_simple_deletion(),
            e.is_mnp()
        ));
    }
    for locus in ref_start - 1..=ref_start + reference.len() as i32 {
        let overlapping = map.overlapping_events(locus);
        if !overlapping.is_empty() {
            let indices: Vec<String> = overlapping
                .iter()
                .map(|o| events.iter().position(|e| e == o).unwrap().to_string())
                .collect();
            out.push(format!("overlap\t{label}\t{locus}\t{}", indices.join(",")));
        }
    }
    h.set_event_map(map);
    (out, Some(h))
}

fn event(start: i32, reference: &str, alt: &str) -> Event {
    Event::new(
        "chr1",
        start,
        Allele::create(reference.as_bytes(), true).unwrap(),
        Allele::create(alt.as_bytes(), false).unwrap(),
    )
    .unwrap()
}

fn render(label: &str, result: Result<Event, EventError>) -> String {
    match result {
        Ok(e) => format!(
            "direct\t{label}\t{}\t{}\t{}\t{}",
            e.start(),
            e.end(),
            e.ref_allele().base_string(),
            e.alt_allele().base_string()
        ),
        Err(e) => error(label, &e),
    }
}

/// The dump's `direct` and `compound` calls, in its order.
fn direct_rows() -> Vec<String> {
    let mut out = Vec::new();
    for (label, reference, alt) in [
        ("trim-two", "ACGT", "AGT"),
        ("trim-to-mnp", "ACGTT", "AGCTT"),
        ("no-trim-length-one", "A", "AGT"),
        ("no-trim-last-base-differs", "ACG", "AT"),
        ("identical", "ACG", "ACG"),
        ("lowercase-shared", "ACgt", "AGgt"),
    ] {
        let result = Event::new(
            "chr1",
            1019,
            Allele::create(reference.as_bytes(), true).unwrap(),
            Allele::create(alt.as_bytes(), false).unwrap(),
        );
        out.push(render(label, result));
    }
    let pairs = [
        ("snp+snp", event(1019, "G", "T"), event(1019, "G", "C")),
        (
            "snp+insertion",
            event(1019, "G", "T"),
            event(1019, "G", "GCC"),
        ),
        (
            "insertion+snp",
            event(1019, "G", "GCC"),
            event(1019, "G", "T"),
        ),
        (
            "snp+deletion",
            event(1019, "G", "T"),
            event(1019, "GAT", "G"),
        ),
        (
            "deletion+insertion",
            event(1019, "GAT", "G"),
            event(1019, "G", "GCT"),
        ),
        (
            "insertion+deletion",
            event(1019, "G", "GCC"),
            event(1019, "GAT", "G"),
        ),
        (
            "insertion+insertion",
            event(1019, "G", "GCC"),
            event(1019, "G", "GA"),
        ),
        (
            "different-starts",
            event(1019, "G", "T"),
            event(1020, "A", "AC"),
        ),
    ];
    for (label, e1, e2) in pairs {
        out.push(render(label, make_compound_events(&e1, &e2)));
    }
    let result = Event::new(
        "chr1",
        1019,
        Allele::create(b"G", false).unwrap(),
        Allele::create(b"T", false).unwrap(),
    );
    out.push(render("non-reference-ref", result));
    out
}

#[test]
fn every_row_matches_the_reference() {
    let golden = golden();
    let mut expected: Vec<&str> = Vec::new();
    let mut got: Vec<String> = Vec::new();
    let mut built: Vec<Haplotype> = Vec::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let fields: Vec<&str> = line.split('\t').collect();
        match fields[0] {
            "case" => {
                let (rows, h) = run(&fields);
                got.extend(rows);
                built.extend(h);
            }
            _ => expected.push(line),
        }
    }
    got.extend(direct_rows());
    let starts = event_start_positions(&built).unwrap();
    got.push(format!(
        "starts\tall\t{}",
        starts
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    ));
    assert_eq!(got.len(), expected.len(), "{}", got.join("\n"));
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g, e);
    }
}
