//! Conformance for `ValidateBasicSomaticShortMutations` against GATK 4.6.2.0, compared as every
//! validation table, every concordance summary and every record of every annotated VCF.
//!
//! Golden from `tools/readfilter-conformance/ValidateBasicSomaticShortMutationsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the validation depth is the counter's reference plus first alternate**, not the pileup's
//!    size, and the fixture's reads overlap their neighbours' loci so the two differ;
//!  * **a validatable genotype whose result is null kills the run**, which the `zero-ad` row is;
//!  * **the artifact test is strictly greater**, so the record with exactly one alternate read in
//!    the validation normal still validates while the one with three does not;
//!  * **an artifact is powered whatever its power is**, so it is a false positive either way;
//!  * **the judgment and the table's `validated` column disagree for an artifact**, the column
//!    being `isOutOfNoiseFloor`, which knows nothing about the normal;
//!  * **a genotype that cannot be validated is written to the annotated VCF and nowhere else**;
//!  * **and a missing validation control sample skips the record in silence**, while a missing
//!    case sample is a null pileup and a table row of zeros.
//!
//! # What is compared, and what is not
//!
//! The `##GATKCommandLine` line is not in the golden: it carries the whole command line and a date.
//! The header lines this port builds are compared against the rest of the golden's header, with one
//! caveat carried in [`gatk_tools::validate_basic_somatic_short_mutations::header_lines`]: the `FT`
//! line the output declares is htsjdk's reserved definition rather than the one the input file
//! wrote, so the input lines fed here are the ones htsjdk holds after reading.

use gatk_corpus as corpus;
use gatk_engine::basic_somatic_short_mutation_validator::{
    write_table, BasicValidationResult, ValidationGenotype,
};
use gatk_engine::read_pileup::{pileup_from_reads, ReadPileup};
use gatk_engine::variant_context_utils::Allele;
use gatk_tools::concordance::Summary;
use gatk_tools::validate_basic_somatic_short_mutations::{
    apply, count_towards_summary, header_lines, Applied, Arguments, ToolError,
};
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const CONTIG_LENGTH: i32 = 200;
const CASE: &str = "valcase";
const CONTROL: &str = "valcontrol";

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/validate_basic_somatic_short_mutations.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The rows of one kind, as `(label, rest)`.
fn rows<'a>(text: &'a str, kind: &str) -> Vec<(&'a str, &'a str)> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split_once('\t').expect("a label"))
        .collect()
}

fn one(text: &str, kind: &str, label: &str) -> Option<String> {
    rows(text, kind)
        .into_iter()
        .find(|(row, _)| *row == label)
        .map(|(_, rest)| unescape(rest))
}

/// The reference is `ACGT` repeating, so a position's base is fixed by its offset.
fn reference_base(position: i32) -> u8 {
    b"ACGT"[((position - 1) % 4) as usize]
}

/// The ten reference bases a `10M` read starting five before the locus carries.
fn window(locus: i32) -> Vec<u8> {
    (locus - 5..=locus + 4).map(reference_base).collect()
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for contig in ["chr1", "chr2"] {
        header
            .sequences
            .push(SequenceRecord::new(contig, CONTIG_LENGTH));
    }
    for (id, sample) in [("rgcase", CASE), ("rgcontrol", CONTROL)] {
        let mut group = ReadGroup::new(id);
        group.attributes.set("SM", sample);
        group.attributes.set("PL", "ILLUMINA");
        header.read_groups.push(group);
    }
    header
}

const DUPLICATE: u16 = 0x400;
const VENDOR_FAILED: u16 = 0x200;

fn record(
    sample: &str,
    name: &str,
    locus: i32,
    cigar: &str,
    bases: Vec<u8>,
    quality: u8,
) -> BamRecord {
    let mut tags = htsjdk_bam::tag::Tags::new();
    let group = if sample == CASE {
        "rgcase"
    } else {
        "rgcontrol"
    };
    tags.insert(Tag::new(b"RG"), TagValue::Str(group.to_string()));
    BamRecord {
        read_name: format!("{sample}-{locus}-{name}"),
        reference_index: 0,
        alignment_start: locus - 5,
        mapping_quality: 60,
        base_qualities: vec![quality; bases.len()],
        read_bases: bases,
        cigar: htsjdk_bam::text_parse::parse_cigar(cigar).expect("a cigar"),
        tags,
        ..Default::default()
    }
}

fn substitution(sample: &str, name: &str, locus: i32, base: u8, quality: u8) -> BamRecord {
    let mut bases = window(locus);
    bases[5] = base;
    record(sample, name, locus, "10M", bases, quality)
}

#[allow(
    clippy::too_many_arguments,
    reason = "the harness's own read specification"
)]
fn substitutions(
    into: &mut Vec<BamRecord>,
    sample: &str,
    locus: i32,
    reference: u8,
    reference_count: usize,
    alternate: u8,
    alternate_count: usize,
    quality: u8,
) {
    for i in 0..reference_count {
        into.push(substitution(
            sample,
            &format!("r{i}"),
            locus,
            reference,
            quality,
        ));
    }
    for i in 0..alternate_count {
        into.push(substitution(
            sample,
            &format!("a{quality}-{i}"),
            locus,
            alternate,
            quality,
        ));
    }
}

fn insertions(into: &mut Vec<BamRecord>, sample: &str, locus: i32, inserted: &str, count: usize) {
    let reference = window(locus);
    let mut bases = reference[..6].to_vec();
    bases.extend_from_slice(inserted.as_bytes());
    bases.extend_from_slice(&reference[6..]);
    for i in 0..count {
        into.push(record(
            sample,
            &format!("i{i}"),
            locus,
            &format!("6M{}I4M", inserted.len()),
            bases.clone(),
            30,
        ));
    }
}

fn deletions(into: &mut Vec<BamRecord>, sample: &str, locus: i32, length: i32, count: usize) {
    let mut bases: Vec<u8> = (locus - 5..=locus).map(reference_base).collect();
    bases.extend((locus + 1 + length..=locus + 4 + length).map(reference_base));
    for i in 0..count {
        into.push(record(
            sample,
            &format!("d{i}"),
            locus,
            &format!("6M{length}D4M"),
            bases.clone(),
            30,
        ));
    }
}

/// The fixture the harness wrote, read for read.
fn fixture() -> Vec<BamRecord> {
    let mut reads = Vec::new();
    substitutions(&mut reads, CASE, 20, b'T', 22, b'A', 8, 30);
    substitutions(&mut reads, CASE, 20, b'T', 0, b'A', 4, 5);
    substitutions(&mut reads, CONTROL, 20, b'T', 20, b'A', 0, 30);
    substitutions(&mut reads, CASE, 24, b'T', 30, b'G', 0, 30);
    substitutions(&mut reads, CONTROL, 24, b'T', 20, b'G', 0, 30);
    substitutions(&mut reads, CASE, 30, b'C', 20, b'A', 10, 30);
    substitutions(&mut reads, CONTROL, 30, b'C', 20, b'A', 3, 30);
    substitutions(&mut reads, CASE, 34, b'C', 20, b'C', 0, 30);
    insertions(&mut reads, CASE, 34, "GGG", 6);
    substitutions(&mut reads, CONTROL, 34, b'C', 20, b'C', 0, 30);
    substitutions(&mut reads, CASE, 40, b'T', 20, b'T', 0, 30);
    deletions(&mut reads, CASE, 40, 1, 6);
    substitutions(&mut reads, CONTROL, 40, b'T', 20, b'T', 0, 30);
    substitutions(&mut reads, CASE, 44, b'T', 20, b'A', 5, 30);
    substitutions(&mut reads, CONTROL, 44, b'T', 20, b'A', 0, 30);
    substitutions(&mut reads, CASE, 48, b'T', 10, b'A', 0, 30);
    substitutions(&mut reads, CONTROL, 48, b'T', 10, b'A', 0, 30);
    substitutions(&mut reads, CASE, 56, b'T', 20, b'A', 8, 30);
    substitutions(&mut reads, CONTROL, 56, b'T', 20, b'A', 1, 30);

    // The three reads the tool's own filters drop, all alternate at position 20.
    let mut duplicate = substitution(CASE, "dup", 20, b'A', 30);
    duplicate.flags |= DUPLICATE;
    reads.push(duplicate);
    let mut vendor = substitution(CASE, "vendor", 20, b'A', 30);
    vendor.flags |= VENDOR_FAILED;
    reads.push(vendor);
    let mut unmapped = substitution(CASE, "mapq0", 20, b'A', 30);
    unmapped.mapping_quality = 0;
    reads.push(unmapped);
    reads
}

/// One record of the discovery VCF, as much of it as the walker reads.
struct Record {
    start: i32,
    stop: i32,
    reference: Allele,
    alternates: Vec<Allele>,
    filters: Vec<String>,
    genotype: ValidationGenotype,
    /// `VariantContext.isSNP()`: biallelic, and both alleles one base.
    is_snp: bool,
}

fn allele(bases: &str, is_reference: bool) -> Allele {
    Allele::new(bases.as_bytes(), is_reference)
}

fn record_of(
    start: i32,
    reference: &str,
    alternates: &[&str],
    filters: &[&str],
    ad: Option<Vec<i32>>,
    genotype_filters: Option<&str>,
) -> Record {
    let reference_allele = allele(reference, true);
    let alternate_alleles: Vec<Allele> = alternates
        .iter()
        .map(|bases| allele(bases, false))
        .collect();
    let is_snp = alternate_alleles.len() == 1
        && reference_allele.len() == 1
        && alternate_alleles[0].len() == 1;
    Record {
        start,
        stop: start + reference.len() as i32 - 1,
        genotype: ValidationGenotype {
            alleles: vec![reference_allele.clone(), alternate_alleles[0].clone()],
            ad,
            filters: genotype_filters.map(str::to_string),
        },
        reference: reference_allele,
        alternates: alternate_alleles,
        filters: filters.iter().map(|f| f.to_string()).collect(),
        is_snp,
    }
}

/// The dump's eight records, in its order.
fn variants() -> Vec<Record> {
    vec![
        record_of(20, "T", &["A"], &[], Some(vec![40, 10]), None),
        record_of(24, "T", &["G"], &[], Some(vec![40, 10]), None),
        record_of(30, "C", &["A"], &[], Some(vec![40, 10]), None),
        record_of(34, "C", &["CGGG"], &[], Some(vec![40, 10]), None),
        record_of(40, "TA", &["T"], &[], Some(vec![40, 10]), None),
        record_of(44, "T", &["A", "C"], &[], Some(vec![30, 10, 5]), None),
        record_of(48, "T", &["A"], &[], None, None),
        record_of(
            56,
            "T",
            &["A"],
            &["weak_evidence", "strand_bias"],
            Some(vec![40, 10]),
            Some("base_qual"),
        ),
    ]
}

/// The one record whose result is null while its genotype is validatable.
fn zero_ad() -> Vec<Record> {
    vec![record_of(20, "T", &["A"], &[], Some(vec![0, 0]), None)]
}

/// The reads the walker's own filters keep: a mapping quality of zero is dropped here, and the
/// duplicate and vendor-failed reads are dropped again by the pileup constructor.
fn kept(reads: &[BamRecord]) -> Vec<BamRecord> {
    reads
        .iter()
        .filter(|read| read.mapping_quality != 0)
        .cloned()
        .collect()
}

/// The pileup of one sample at one record's start.
fn pileup_for<'a>(
    reads: &'a [BamRecord],
    header: &SamHeader,
    start: i32,
    sample: &str,
) -> Option<ReadPileup<'a>> {
    let whole = pileup_from_reads(
        "chr1",
        start,
        reads,
        |read| read.flags & VENDOR_FAILED == 0,
        |read| read.flags & DUPLICATE == 0,
    );
    if whole.elements.is_empty() {
        // `pileupsBySample.isEmpty()`, which is the record being skipped before either sample is
        // looked up.
        return None;
    }
    whole
        .split_by_sample(header, Some("__UNKNOWN__"))
        .expect("every read has a read group")
        .into_iter()
        .find(|(name, _)| name == sample)
        .map(|(_, pileup)| pileup)
}

/// One whole run: the table, the summary and the INFO field of every record written.
#[derive(Debug)]
struct Run {
    table: String,
    summary: String,
    info: Vec<(i32, Option<String>)>,
}

fn run(
    records: &[Record],
    arguments: &Arguments,
    window: Option<(i32, i32)>,
) -> Result<Run, ToolError> {
    let reads = kept(&fixture());
    let header = header();
    let mut results: Vec<BasicValidationResult> = Vec::new();
    let mut summary = Summary::default();
    let mut info = Vec::new();
    for record in records {
        if let Some((from, to)) = window {
            if record.start < from || record.start > to {
                continue;
            }
        }
        let case = pileup_for(
            &reads,
            &header,
            record.start,
            &arguments.validation_case_name,
        );
        let control = pileup_for(
            &reads,
            &header,
            record.start,
            &arguments.validation_control_name,
        );
        let applied: Option<Applied> = apply(
            "chr1",
            record.start,
            record.stop,
            &record.reference,
            &record.alternates,
            &record.filters,
            &record.genotype,
            case.as_ref(),
            control.as_ref(),
            arguments,
        )?;
        let Some(applied) = applied else {
            continue;
        };
        if let Some(result) = &applied.result {
            results.push(result.clone());
        }
        count_towards_summary(&mut summary, &applied, record.is_snp);
        info.push((record.start, Some(applied.info())));
    }
    Ok(Run {
        table: write_table(&results),
        summary: summary.table(),
        info,
    })
}

/// The arguments of each labelled run.
fn arguments(label: &str) -> Arguments {
    let base = Arguments {
        discovery_sample: "discovery".to_string(),
        validation_case_name: CASE.to_string(),
        validation_control_name: CONTROL.to_string(),
        ..Arguments::default()
    };
    match label {
        "default" | "zero-ad" | "no-records" => base,
        "min-power-zero" => Arguments {
            min_power: 0.0,
            ..base
        },
        "normal-count-three" => Arguments {
            max_validation_normal_count: 3,
            ..base
        },
        "cutoff-zero" => Arguments {
            min_bq_cutoff: 0,
            ..base
        },
        "cutoff-fifty" => Arguments {
            min_bq_cutoff: 50,
            ..base
        },
        "missing-control" => Arguments {
            validation_control_name: "absent".to_string(),
            ..base
        },
        "missing-case" => Arguments {
            validation_case_name: "absent".to_string(),
            ..base
        },
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The INFO fields the golden's annotated VCF carries, by position.
fn golden_info(text: &str, label: &str) -> Vec<(i32, Option<String>)> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|(row, _)| *row == label)
        .map(|(_, rest)| unescape(rest))
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            (
                field[1].parse().expect("a position"),
                Some(field[7].to_string()),
            )
        })
        .collect()
}

#[test]
fn every_run_writes_what_the_reference_writes() {
    let text = golden();
    let mut compared = 0;
    for label in [
        "default",
        "min-power-zero",
        "normal-count-three",
        "cutoff-zero",
        "cutoff-fifty",
        "missing-control",
        "missing-case",
        "no-records",
    ] {
        let window = if label == "no-records" {
            Some((150, 160))
        } else {
            None
        };
        let ours = run(&variants(), &arguments(label), window).expect("a run that finishes");
        assert_eq!(
            ours.table,
            one(&text, "table", label).expect("the golden carries the table"),
            "{label}: the validation table"
        );
        assert_eq!(
            ours.summary,
            one(&text, "summary", label).expect("the golden carries the summary"),
            "{label}: the concordance summary"
        );
        assert_eq!(
            ours.info,
            golden_info(&text, label),
            "{label}: the INFO fields"
        );
        compared += 1;
    }
    assert_eq!(compared, 8, "the golden's runs");

    // The one run that does not finish.
    let error = run(&zero_ad(), &arguments("zero-ad"), None).expect_err("a null result");
    assert_eq!(
        format!("error\tzero-ad\t{}:{}", error.java_class(), error.message()),
        text.lines()
            .find(|line| line.starts_with("error\tzero-ad"))
            .expect("the golden carries the refusal")
    );
}

#[test]
fn the_annotated_header_is_the_input_plus_three_info_lines_and_a_source() {
    let text = golden();
    let input = [
        "##fileformat=VCFv4.2",
        "##FILTER=<ID=strand_bias,Description=\"strand\">",
        "##FILTER=<ID=weak_evidence,Description=\"weak\">",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
        // htsjdk's own reserved definition, which replaces what the input file declared.
        "##FORMAT=<ID=FT,Number=.,Type=String,Description=\"Genotype-level filter\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##contig=<ID=chr1,length=200>",
        "##contig=<ID=chr2,length=200>",
    ]
    .map(str::to_string);
    let theirs: Vec<String> = rows(&text, "vcfline")
        .into_iter()
        .filter(|(label, _)| *label == "default")
        .map(|(_, rest)| unescape(rest))
        .filter(|line| line.starts_with("##"))
        .collect();
    assert_eq!(header_lines(&input), theirs);
}
