//! Conformance for `HaplotypeCaller` against GATK 4.6.2.0, every run of the dump replayed through
//! the command line and compared as the records the reference wrote.
//!
//! Golden from `tools/readfilter-conformance/HcEndToEndDump.java`, which writes two random contigs
//! and two samples' reads (SNPs, a deletion, an insertion, SNPs in cis and in trans, reads at low
//! mapping quality, duplicates, soft clips, IUPAC bases and one over-deep start) as an indexed BAM,
//! and runs the tool seven times in ONE JVM. The golden carries the contigs, the SAM text the BAM
//! was written from, and each output from its `#CHROM` line on. This test writes the same BAM and
//! index from that SAM text.
//!
//! The runs go in the dump's order because the generator a `QD` above 35 draws from is one per
//! process on both sides; this file holds the only test that drives the tool.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_end_to_end.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

fn scratch() -> std::path::PathBuf {
    let dir = std::env::temp_dir().join("gatk-cli-haplotype-caller");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    dir
}

/// `writeReference`: each contig in lines of sixty, with the index and dictionary.
fn write_reference(dir: &std::path::Path, contigs: &[(&str, String)]) -> String {
    let mut fasta = String::new();
    let mut fai = String::new();
    let mut dict = String::from("@HD\tVN:1.6\n");
    for (name, bases) in contigs {
        fasta.push_str(&format!(">{name}\n"));
        let offset = fasta.len();
        for line in bases.as_bytes().chunks(60) {
            fasta.push_str(std::str::from_utf8(line).expect("ASCII"));
            fasta.push('\n');
        }
        fai.push_str(&format!("{name}\t{}\t{offset}\t60\t61\n", bases.len()));
        dict.push_str(&format!("@SQ\tSN:{name}\tLN:{}\n", bases.len()));
    }
    std::fs::write(dir.join("reference.fasta"), fasta).expect("the reference");
    std::fs::write(dir.join("reference.fasta.fai"), fai).expect("the index");
    std::fs::write(dir.join("reference.dict"), dict).expect("the dictionary");
    dir.join("reference.fasta").to_string_lossy().to_string()
}

/// The BAM and its index, written from the golden's SAM text.
fn write_bam(dir: &std::path::Path, sam: &str) -> String {
    let (header, records) = htsjdk_bam::sam_file::read_sam(sam).expect("the SAM text");
    let mut writer = htsjdk_bam::writer::BamWriter::new(Vec::new(), &header)
        .expect("a writer")
        .with_index();
    for record in &records {
        writer.write(record).expect("a record");
    }
    let (bam, bai) = writer.finish_with_index().expect("the BAM");
    std::fs::write(dir.join("input.bam"), bam).expect("the BAM file");
    std::fs::write(dir.join("input.bai"), bai).expect("the index file");
    dir.join("input.bam").to_string_lossy().to_string()
}

#[test]
fn every_run_matches_the_golden_in_the_dumps_order() {
    let text = golden();
    let dir = scratch();
    let contigs: Vec<(&str, String)> = ["chr1", "chr2"]
        .iter()
        .map(|name| (*name, section(&text, "fasta", name)))
        .collect();
    let reference = write_reference(&dir, &contigs);
    let bam = write_bam(&dir, &section(&text, "sam", "input"));
    let runs: [(&str, &[&str]); 7] = [
        ("default", &[]),
        ("interval", &["-L", "chr1:200-600", "-L", "chr2"]),
        ("no-phasing", &["--do-not-run-physical-phasing", "true"]),
        (
            "call-threshold-500",
            &["--standard-min-confidence-threshold-for-calling", "500"],
        ),
        ("no-downsampling", &["--max-reads-per-alignment-start", "0"]),
        (
            "low-mapq-kept",
            &["--disable-read-filter", "MappingQualityReadFilter"],
        ),
        ("min-base-quality-30", &["--min-base-quality-score", "30"]),
    ];
    for (case, extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            "HaplotypeCaller",
            "-I",
            &bam,
            "-O",
            &output.to_string_lossy(),
            "-R",
            &reference,
        ]
        .iter()
        .map(|arg| arg.to_string())
        .collect();
        args.extend(extra.iter().map(|arg| arg.to_string()));
        let outcome = gatk_cli::run(&args);
        assert_eq!(outcome.status, 0, "{case}: {}", outcome.stderr);
        let written = std::fs::read_to_string(&output).expect("the output");
        let body: String = written
            .lines()
            .skip_while(|line| !line.starts_with("#CHROM"))
            .filter(|line| !line.is_empty())
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(body, section(&text, "out", case), "{case}");
    }
}
