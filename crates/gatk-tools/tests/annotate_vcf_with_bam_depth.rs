//! Conformance for `AnnotateVcfWithBamDepth` against GATK 4.6.2.0, compared as the `BAM_DEPTH` of
//! every written record and as the header lines the tool builds.
//!
//! Golden from `tools/readfilter-conformance/AnnotateVcfWithBamDepthDump.java`.
//!
//! # What this suite is for
//!
//!  * **the output header declares `BAM_DEPTH` twice** when the input already declared it;
//!  * **a read one base long is never counted**;
//!  * **the read must contain the variant's whole span**, `END` included;
//!  * **duplicates, vendor-failed and unmapped reads are excluded by the tool itself**;
//!  * **and the annotation is written at zero and overwrites what the record carried**.
//!
//! # What is compared, and what is not
//!
//! The `##GATKCommandLine` line is the dump's own row and is not replayed: it carries the whole
//! command line and a masked date, neither of which this port produces. Everything else in the
//! output, header lines and records alike, is compared.

use gatk_corpus as corpus;
use gatk_tools::annotate_vcf_with_bam_depth::{
    annotate, bam_depth, header_lines, Read, BAM_DEPTH, BAM_DEPTH_HEADER_LINE,
};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/annotate_vcf_with_bam_depth.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The reads of the fixture BAM, as the five conditions see them.
fn reads(text: &str) -> Vec<BamRecord> {
    let encoded = rows(text, "fixture")
        .into_iter()
        .find(|row| row[0] == "reads")
        .expect("the bam")[1]
        .to_string();
    let bytes = corpus::decode_base64(&encoded);
    let decompressed = htsjdk_bgzf::read::decompress_all(&bytes).expect("the fixture is BGZF");
    let reader = htsjdk_bam::reader::BamReader::new(&decompressed).expect("the fixture opens");
    reader.map(|record| record.expect("a record")).collect()
}

fn as_read(record: &BamRecord) -> Read<'_> {
    Read {
        contig: "chr1",
        start: record.alignment_start,
        end: record.alignment_start + record.cigar.reference_length() as i32 - 1,
        flags: record.flags,
    }
}

/// The input's records, decoded as far as the conditions and the annotation look at them.
fn variants(text: &str) -> Vec<VariantContext> {
    let whole = unescape(
        rows(text, "input")
            .into_iter()
            .find(|row| row[0] == "variants")
            .expect("the input")[1],
    );
    whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let mut alleles = vec![Allele::create(field[3].as_bytes(), true).expect("a reference")];
            for alternate in field[4].split(',') {
                alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
            }
            let start: i64 = field[1].parse().expect("a position");
            let mut variant = VariantContext::new(field[0], start, alleles);
            variant.stop = start + field[3].len() as i64 - 1;
            if field[7] != "." {
                for entry in field[7].split(';') {
                    if let Some((key, value)) = entry.split_once('=') {
                        if key == "END" {
                            variant.stop = value.parse().expect("an end");
                        }
                        variant
                            .attributes
                            .push((key.to_string(), Value::Str(value.to_string())));
                    }
                }
            }
            variant
        })
        .collect()
}

/// The `BAM_DEPTH` of every record of one run, taken off the golden's own output lines.
fn written_depths(text: &str, label: &str) -> Vec<(i64, i32)> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == label)
        .map(|row| unescape(row[1]))
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let depth = field[7]
                .split(';')
                .filter_map(|entry| entry.split_once('='))
                .find(|(key, _)| *key == BAM_DEPTH)
                .map(|(_, value)| value.parse::<i32>().expect("an integer"))
                .expect("every record carries the annotation");
            (field[1].parse().expect("a position"), depth)
        })
        .collect()
}

#[test]
fn every_depth_matches_the_golden() {
    let text = golden();
    let records = reads(&text);
    let pooled: Vec<Read> = records.iter().map(as_read).collect();
    let mine: Vec<(i64, i32)> = variants(&text)
        .iter()
        .map(|variant| (variant.start, bam_depth(&pooled, variant)))
        .collect();
    assert_eq!(mine, written_depths(&text, "annotated"));
}

#[test]
fn without_reads_every_record_is_zero() {
    let text = golden();
    let mine: Vec<(i64, i32)> = variants(&text)
        .iter()
        .map(|variant| (variant.start, bam_depth(&[], variant)))
        .collect();
    assert_eq!(mine, written_depths(&text, "no-reads"));
    assert!(mine.iter().all(|(_, depth)| *depth == 0));
}

#[test]
fn a_one_base_read_and_a_partial_cover_are_both_zero() {
    let text = golden();
    let depths = written_depths(&text, "annotated");
    // The site only a 1M read sits on.
    assert_eq!(
        depths.iter().find(|(at, _)| *at == 40).expect("the site").1,
        0
    );
    // The deletion, covered wholly by one of its two reads.
    assert_eq!(
        depths
            .iter()
            .find(|(at, _)| *at == 60)
            .expect("the deletion")
            .1,
        1
    );
    // The block carrying END, which no read contains.
    assert_eq!(
        depths
            .iter()
            .find(|(at, _)| *at == 80)
            .expect("the block")
            .1,
        0
    );
}

#[test]
fn the_annotation_overwrites_what_the_record_carried() {
    let text = golden();
    let carried = variants(&text)
        .into_iter()
        .find(|variant| variant.start == 160)
        .expect("the record already carrying BAM_DEPTH");
    assert!(carried
        .attributes
        .iter()
        .any(|(key, value)| key == BAM_DEPTH && *value == Value::Str("99".to_string())));

    let annotated = annotate(&carried, 0);
    assert_eq!(
        annotated
            .attributes
            .iter()
            .filter(|(key, _)| key == BAM_DEPTH)
            .collect::<Vec<_>>(),
        vec![&(BAM_DEPTH.to_string(), Value::Int(0))]
    );
    assert_eq!(
        written_depths(&text, "annotated")
            .iter()
            .find(|(at, _)| *at == 160)
            .expect("the record")
            .1,
        0
    );
}

#[test]
fn the_output_header_declares_bam_depth_twice() {
    let text = golden();
    let written: Vec<String> = rows(&text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == "annotated")
        .map(|row| unescape(row[1]))
        .filter(|line| line.starts_with("##INFO=<ID=BAM_DEPTH"))
        .collect();
    assert_eq!(written.len(), 2, "the input's line and the tool's");
    assert!(written.contains(&BAM_DEPTH_HEADER_LINE.to_string()));

    let input: Vec<String> = unescape(
        rows(&text, "input")
            .into_iter()
            .find(|row| row[0] == "variants")
            .expect("the input")[1],
    )
    .lines()
    .filter(|line| line.starts_with("##"))
    .map(|line| line.to_string())
    .collect();
    let mine = header_lines(&input);
    assert_eq!(
        mine.iter()
            .filter(|line| line.starts_with("##INFO=<ID=BAM_DEPTH"))
            .count(),
        2
    );
}

#[test]
fn the_refusal_is_htsjdks_own() {
    let text = golden();
    let row = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "output-is-a-directory")
        .expect("the refusal");
    assert_eq!(
        row[1],
        "htsjdk.samtools.util.RuntimeIOException:File not found: annotatevcfwithbamdepth-dump/."
    );
}
