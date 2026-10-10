//! Conformance for `AssemblyRegionTrimmer` against the oracle.
//!
//! Golden from `tools/readfilter-conformance/AssemblyRegionTrimmerDump.java`. The golden carries
//! the contig, the reads and each case's region, arguments and events; this test rebuilds them,
//! trims, and renders the `result`, `region` and `error` rows, which must be the golden's in order.

use gatk_corpus as corpus;
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::assembly_region_trimmer::{
    AssemblyRegionTrimmer, TrimResult, TrimmerArguments, TrimmerError,
};
use gatk_engine::event_map::Event;
use gatk_engine::interval::SimpleInterval;
use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;

fn render(i: Option<&SimpleInterval>) -> String {
    i.map_or("null".to_string(), |i| {
        format!("{}:{}-{}", i.contig, i.start, i.end)
    })
}

fn render_error(e: &TrimmerError) -> String {
    match e {
        TrimmerError::NoVariation => {
            "IllegalStateException: There is no variation present.".to_string()
        }
        TrimmerError::NullPointer => "NullPointerException: Cannot invoke \"org.broadinstitute.hellbender.utils.SimpleInterval.getEnd()\" because \"this.variantSpan\" is null".to_string(),
        TrimmerError::BadArgument(m) => format!("BadArgumentValue: {m}"),
        other => format!("{other:?}"),
    }
}

fn region_row(
    label: &str,
    which: &str,
    result: Result<Option<AssemblyRegion>, TrimmerError>,
) -> String {
    match result {
        Ok(None) => format!("region\t{label}\t{which}\tempty"),
        Ok(Some(r)) => format!(
            "region\t{label}\t{which}\t{}\t{}\t{}",
            render(Some(r.span())),
            render(Some(r.padded_span())),
            r.reads()
                .iter()
                .map(|read| read.read_name.as_str())
                .collect::<Vec<_>>()
                .join(",")
        ),
        Err(e) => format!("error\t{label}\t{which}\t{}", render_error(&e)),
    }
}

#[test]
fn every_trim_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/assembly_region_trimmer.txt.gz"),
    );
    let mut contig = String::new();
    let mut reads: Vec<BamRecord> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();
    let mut cases: Vec<(Vec<&str>, Vec<Vec<&str>>)> = Vec::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "contig" => contig = f[1].to_string(),
            "read" => reads.push(BamRecord {
                read_name: f[1].to_string(),
                reference_index: 0,
                alignment_start: f[2].parse().unwrap(),
                mapping_quality: 60,
                cigar: htsjdk_bam::text_parse::parse_cigar(f[3]).unwrap(),
                read_bases: vec![b'A'; 50],
                base_qualities: vec![30; 50],
                ..Default::default()
            }),
            "case" => cases.push((f, Vec::new())),
            "event" => cases.last_mut().unwrap().1.push(f),
            _ => expected.push(line),
        }
    }
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", contig.len() as _));

    let mut got: Vec<String> = Vec::new();
    for (case, event_rows) in &cases {
        let label = case[1];
        let n = |i: usize| -> i32 { case[i].parse().unwrap() };
        let active = SimpleInterval::new("chr1", n(2), n(3)).unwrap();
        let mut region = AssemblyRegion::with_padding(active, true, n(4), &header).unwrap();
        let padded = region.padded_span().clone();
        for read in &reads {
            let end = read.alignment_start + 49;
            if padded.overlaps("chr1", read.alignment_start, end) {
                region.add(read.clone(), &header).unwrap();
            }
        }
        let arguments = TrimmerArguments {
            snp_padding_for_genotyping: n(6),
            indel_padding_for_genotyping: n(7),
            str_padding_for_genotyping: n(8),
            max_extension_into_region_padding: n(9),
            assembly_region_padding: n(10),
            enable_legacy_assembly_region_trimming: case[5] == "true",
        };
        let events: Vec<Event> = event_rows
            .iter()
            .map(|e| {
                Event::new(
                    "chr1",
                    e[2].parse().unwrap(),
                    Allele::create(e[3].as_bytes(), true).unwrap(),
                    Allele::create(e[4].as_bytes(), false).unwrap(),
                )
                .unwrap()
            })
            .collect();
        let reference = &contig.as_bytes()[padded.start as usize - 1..padded.end as usize];
        let trimmed = AssemblyRegionTrimmer::new(arguments, contig.len() as i32).and_then(|t| {
            t.trim(&region, &events, reference, padded.start)
                .map(|r| (t, r))
        });
        let (trimmer, result): (AssemblyRegionTrimmer, TrimResult) = match trimmed {
            Ok(x) => x,
            Err(e) => {
                got.push(format!("error\t{label}\ttrim\t{}", render_error(&e)));
                continue;
            }
        };
        got.push(format!(
            "result\t{label}\t{}\t{}",
            render(result.variant_span.as_ref()),
            render(result.padded_span.as_ref())
        ));
        got.push(region_row(
            label,
            "variant",
            trimmer.variant_region(&result, &region, &header).map(Some),
        ));
        got.push(region_row(
            label,
            "left",
            trimmer.non_variant_left_flank_region(&result, &region, &header),
        ));
        got.push(region_row(
            label,
            "right",
            trimmer.non_variant_right_flank_region(&result, &region, &header),
        ));
    }
    assert_eq!(got.len(), expected.len(), "{}", got.join("\n"));
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g, e);
    }
}
