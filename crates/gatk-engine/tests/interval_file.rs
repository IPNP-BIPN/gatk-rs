//! Conformance for `IntervalUtils.parseIntervalArguments` against GATK 4.6.2.0.
//!
//! One row per `-L` argument: what it resolved to, or which exception refused it.
//!
//! The golden corrected an assumption the port was written on. `FeatureManager.isFeatureFile`
//! asks every codec whether it `canDecode` the path, and the codecs answer by **extension** rather
//! than by content, so a `.list` file holding a BED body is not a Feature file: it falls through
//! to the interval-file reader and dies parsing `chr1\t0\t10` as a genome location. The row named
//! `bed-contents-list-extension` is that measurement.
//!
//! Every case is compared. The two that waited for codecs, `.bed` and `.interval_list`, are
//! resolved by [`gatk_engine::feature_intervals::RegisteredCodecs`], and the interval list brought
//! four cases of its own: it is the one format here that validates against **two** dictionaries,
//! its own `@SQ` lines and the reference's, and the two disagreements have different outcomes.
//!
//! ```text
//! case  interval-list-bad-strand                    E:...UserException$MalformedFile
//! case  interval-list-short-record                  E:htsjdk.tribble.TribbleException
//! case  interval-list-unknown-contig                ok  1  chr1:1-10
//! case  interval-list-contig-absent-from-reference  E:...UserException$MalformedGenomeLoc
//! ```
//!
//! One malformed file, two exception classes: `featureFileToIntervals` catches
//! `IllegalArgumentException` and nothing else, so the bad strand is wrapped as a `MalformedFile`
//! and the short record leaves the engine as the `TribbleException` it was. And the two
//! dictionaries fail in opposite directions: a contig the file does not declare costs one line,
//! a contig the reference does not hold costs the argument.

use gatk_corpus as corpus;
use gatk_engine::interval_args::{self, IntervalArgumentError};
use htsjdk_bam::header::{SamHeader, SequenceRecord};

const CONTIG_LENGTH: i32 = 200;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/interval_file.txt.gz"),
    )
}

fn header() -> SamHeader {
    let mut header = SamHeader::default();
    for name in ["chr1", "chr2"] {
        header
            .sequences
            .push(SequenceRecord::new(name, CONTIG_LENGTH));
    }
    header
}

/// The files the harness wrote, by the label of the case that used them. The contents are the
/// input, so they are here rather than in the golden.
fn fixture(label: &str) -> Option<(&'static str, &'static str)> {
    // (file name, contents)
    match label {
        "list" => Some(("a.list", "chr1:1-10\nchr1:50-60\nchr2\n")),
        "intervals" => Some(("b.intervals", "chr1:1-10\nchr1:50-60\nchr2\n")),
        "whitespace" => Some(("ws.list", "\n  chr1:1-10  \n\n\tchr2:5-6\n\n")),
        "blank-only" => Some(("blank.list", "\n\n   \n")),
        "empty" => Some(("empty.list", "")),
        "uppercase-extension" => Some(("c.LIST", "chr1:1-5\n")),
        "unknown-extension" => Some(("d.txt", "chr1:1-10\n")),
        "bed" => Some(("f.bed", "chr1\t0\t10\nchr2\t4\t6\n")),
        "bed-contents-list-extension" => Some(("g.list", "chr1\t0\t10\nchr2\t4\t6\n")),
        // Written by htsjdk's IntervalList writer, header included.
        "picard-interval-list" => Some((
            "e.interval_list",
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\n@SQ\tSN:chr2\tLN:200\n\
             chr1\t1\t10\t+\t.\nchr2\t5\t6\t+\t.\n",
        )),
        // A strand of `.`, which the BED codec in the same package accepts and this one refuses.
        "interval-list-bad-strand" => Some((
            "strand.interval_list",
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\nchr1\t1\t10\t.\t.\n",
        )),
        // Four fields, because the codec counts exactly five.
        "interval-list-short-record" => Some((
            "short.interval_list",
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\nchr1\t1\t10\t+\n",
        )),
        // A contig the file's own header does not declare: dropped, and the file still loads.
        "interval-list-unknown-contig" => Some((
            "unknown.interval_list",
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\nchr1\t1\t10\t+\t.\n\
             chr3\t1\t10\t+\t.\n",
        )),
        // A contig the file declares and the reference dictionary does not.
        "interval-list-contig-absent-from-reference" => Some((
            "foreign.interval_list",
            "@HD\tVN:1.6\n@SQ\tSN:chr9\tLN:200\nchr9\t1\t10\t+\t.\n",
        )),
        // VCF: the one codec that decides by content, so the same body under a `.list` name is
        // still a Feature file.
        "vcf" => Some(("h.vcf", VCF_BODY)),
        "vcf-list-extension" => Some(("i.list", VCF_BODY)),
        "vcf-symbolic-end" => Some((
            "j.vcf",
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=200>\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
             ##ALT=<ID=DEL,Description=\"Deletion\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t20\t.\tA\t<DEL>\t.\t.\tEND=80\n",
        )),
        "vcf-malformed-record" => Some((
            "k.vcf",
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=200>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\tNOTANUMBER\t.\tA\tC\t.\t.\t.\n",
        )),
        "vcf-unknown-contig" => Some((
            "l.vcf",
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr9,length=200>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr9\t10\t.\tA\tC\t.\t.\t.\n",
        )),
        "vcf-magic-only" => Some(("m.vcf", "##fileformat=VCFv4.2\n")),
        _ => None,
    }
}

/// The VCF body two cases share, under two different extensions.
const VCF_BODY: &str = "##fileformat=VCFv4.2\n\
             ##contig=<ID=chr1,length=200>\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t10\t.\tA\tC\t.\t.\t.\n\
             chr2\t50\t.\tACGT\tA\t.\t.\t.\n";

/// The argument each case passed, given the directory the fixtures live in.
fn argument(label: &str, dir: &std::path::Path) -> String {
    match label {
        "missing-list" => dir.join("absent.list").display().to_string(),
        "missing-unknown-extension" => dir.join("absent.txt").display().to_string(),
        "semicolon" => "chr1:1-10;chr2:1-10".to_string(),
        "literal" => "chr1:1-10".to_string(),
        "literal-whole-contig" => "chr2".to_string(),
        other => {
            let (name, _) = fixture(other).unwrap_or_else(|| panic!("{other} has no fixture"));
            dir.join(name).display().to_string()
        }
    }
}

/// `GenomeLoc.toString`, which is not `contig:start-end` for every locus.
///
/// A one-base locus prints as `contig:start`, with no range at all, and the golden's first VCF row
/// is one: a SNV at chr1:10 prints `chr1:10` where a four-base reference allele prints
/// `chr2:50-53`. The whole-contig form (`contig` alone) needs a stop of `Integer.MAX_VALUE`, which
/// nothing here produces, so it is not reached.
fn genome_loc_to_string(interval: &gatk_engine::interval::SimpleInterval) -> String {
    if interval.start == interval.end {
        format!("{}:{}", interval.contig, interval.start)
    } else {
        format!("{}:{}-{}", interval.contig, interval.start, interval.end)
    }
}

/// The exception class the reference raised, for the refusal the port produced.
fn class_of(error: &IntervalArgumentError) -> &'static str {
    match error {
        IntervalArgumentError::IntervalFileEmpty => {
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile"
        }
        IntervalArgumentError::IntervalFileMissing(_)
        | IntervalArgumentError::FileIsNeitherFeaturesNorIntervals(_) => {
            "org.broadinstitute.hellbender.exceptions.UserException$CouldNotReadInputFile"
        }
        IntervalArgumentError::LegacySemicolonSyntax(_) => {
            "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue"
        }
        // `featureFileToIntervals` catches IllegalArgumentException and nothing else, so the same
        // malformed interval list raises two different classes depending on which line is wrong.
        IntervalArgumentError::FeatureFileMalformed(_) => {
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedFile"
        }
        IntervalArgumentError::FeatureCodecRefused(_) => "htsjdk.tribble.TribbleException",
        IntervalArgumentError::FeatureSourceFailed(_) => {
            "org.broadinstitute.hellbender.exceptions.GATKException"
        }
        // Every parse failure surfaces as one class: an unknown contig and malformed positions
        // are different messages of the same exception.
        IntervalArgumentError::Parse(_) => {
            "org.broadinstitute.hellbender.exceptions.UserException$MalformedGenomeLoc"
        }
        other => panic!("{other:?} has no reference class"),
    }
}

#[test]
fn every_argument_resolves_the_way_the_reference_resolves_it() {
    let text = golden();
    let header = header();

    let dir = std::env::temp_dir().join(format!("gatk-rs-intervalfile-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for label in [
        "list",
        "intervals",
        "whitespace",
        "blank-only",
        "empty",
        "uppercase-extension",
        "unknown-extension",
        "bed",
        "bed-contents-list-extension",
        "picard-interval-list",
        "interval-list-bad-strand",
        "interval-list-short-record",
        "interval-list-unknown-contig",
        "interval-list-contig-absent-from-reference",
        "vcf",
        "vcf-list-extension",
        "vcf-symbolic-end",
        "vcf-malformed-record",
        "vcf-unknown-contig",
        "vcf-magic-only",
    ] {
        let (name, contents) = fixture(label).expect("a fixture");
        std::fs::write(dir.join(name), contents).unwrap();
    }

    let mut compared = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("case\t") else {
            continue;
        };
        let mut parts = rest.split('\t');
        let label = parts.next().expect("a label");
        let outcome = parts.next().expect("an outcome");
        let count: usize = parts.next().expect("a count").parse().expect("a number");
        let expected = parts.next().unwrap_or("");

        let result = interval_args::parse_interval_arguments(
            &argument(label, &dir),
            &header,
            &gatk_engine::feature_intervals::RegisteredCodecs,
        );

        match (result, outcome) {
            (Ok(intervals), "ok") => {
                let ours: Vec<String> = intervals.iter().map(genome_loc_to_string).collect();
                assert_eq!(ours.len(), count, "{label}: interval count");
                assert_eq!(ours.join("|"), expected, "{label}");
            }
            (Err(error), outcome) if outcome.starts_with("E:") => {
                assert_eq!(
                    format!("E:{}", class_of(&error)),
                    outcome,
                    "{label}: the wrong refusal"
                );
            }
            (Ok(_), outcome) => panic!("{label}: the reference raised {outcome}, the port did not"),
            (Err(error), _) => {
                panic!("{label}: the port raised {error:?}, the reference did not")
            }
        }
        compared += 1;
    }

    std::fs::remove_dir_all(&dir).ok();
    assert!(compared > 0, "the golden carries no case rows");
    println!("{compared} -L arguments resolved identically, none pending");
}
