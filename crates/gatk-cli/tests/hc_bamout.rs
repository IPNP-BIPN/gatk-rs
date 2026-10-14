//! Conformance for the callers' --bam-output against GATK 4.6.2.0, every run of the dump replayed
//! through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcBamoutDump.java`, over `HcEndToEndDump`'s reads:
//! HaplotypeCaller with each writer type, in GVCF mode, with --disable-optimizations and with the
//! PDHMM, RampedHaplotypeCaller, and Mutect2 tumour-only and tumour-normal. Each VCF is compared
//! from its ##INFO and ##FORMAT lines to its last record; each bamout as SAM text, header included,
//! and whether an index was written beside it.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_bamout.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-bamout");
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
    let set = |name: &str, contigs: &[&str]| -> (String, String) {
        let sub = dir.join(name);
        std::fs::create_dir_all(&sub).expect("a set directory");
        let contigs: Vec<(&str, String)> = contigs
            .iter()
            .map(|&contig| (contig, section(&text, "fasta", &format!("{name}:{contig}"))))
            .collect();
        let reference = write_reference(&sub, &contigs);
        let bam = write_bam(&sub, &section(&text, "sam", name));
        (reference, bam)
    };
    let (e2e_reference, e2e_bam) = set("e2e", &["chr1", "chr2"]);
    let e2e = (e2e_bam.as_str(), e2e_reference.as_str());
    let runs: Vec<Run<'_>> = vec![
        ("hc-bamout", "HaplotypeCaller", e2e, vec![]),
        (
            "hc-bamout-all",
            "HaplotypeCaller",
            e2e,
            vec!["--bam-writer-type", "ALL_POSSIBLE_HAPLOTYPES"],
        ),
        (
            "hc-bamout-no-haplotypes",
            "HaplotypeCaller",
            e2e,
            vec![
                "--bam-writer-type",
                "NO_HAPLOTYPES",
                "--create-output-bam-index",
                "false",
            ],
        ),
        (
            "hc-bamout-no-reads",
            "HaplotypeCaller",
            e2e,
            vec!["--bam-writer-type", "CALLED_HAPLOTYPES_NO_READS"],
        ),
        (
            "hc-bamout-gvcf",
            "HaplotypeCaller",
            e2e,
            vec!["-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"],
        ),
        (
            "hc-bamout-disable-optimizations",
            "HaplotypeCaller",
            e2e,
            vec!["--disable-optimizations", "-L", "chr2"],
        ),
        (
            "hc-bamout-pdhmm",
            "HaplotypeCaller",
            e2e,
            vec![
                "--use-pdhmm",
                "--bam-writer-type",
                "ALL_POSSIBLE_HAPLOTYPES",
                "-L",
                "chr1:1100-1350",
            ],
        ),
        ("ramped-bamout", "RampedHaplotypeCaller", e2e, vec![]),
        ("mutect2-bamout", "Mutect2", e2e, vec![]),
        (
            "mutect2-bamout-all-tn",
            "Mutect2",
            e2e,
            vec![
                "--bam-writer-type",
                "ALL_POSSIBLE_HAPLOTYPES",
                "--normal-sample",
                "s2",
            ],
        ),
    ];
    for (case, tool, (bam, reference), extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let bamout = dir.join(format!("bamout-{case}.bam"));
        let mut args: Vec<String> = [
            tool,
            "-I",
            bam,
            "-O",
            &output.to_string_lossy(),
            "-R",
            reference,
        ]
        .iter()
        .map(|arg| arg.to_string())
        .collect();
        args.extend(extra.iter().map(|arg| arg.to_string()));
        // The JDK deflater: the bytes are not compared, only the SAM text, and GKL's level two is
        // ISA-L's, which a build without the `gkl-igzip` feature cannot reach.
        args.extend([
            "--bam-output".to_string(),
            bamout.to_string_lossy().to_string(),
            "--use-jdk-deflater".to_string(),
        ]);
        let outcome = gatk_cli::run(&args);
        let refusal = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{case}\t")));
        if let Some(refusal) = refusal {
            let message = refusal.split_once(':').expect("a class and a message").1;
            assert_ne!(outcome.status, 0, "{case} should be refused");
            assert!(
                outcome.stderr.contains(message),
                "{case}: {}",
                outcome.stderr
            );
            continue;
        }
        assert_eq!(outcome.status, 0, "{case}: {}", outcome.stderr);
        let written = std::fs::read_to_string(&output).expect("the output");
        let mut in_body = false;
        let body: String = written
            .lines()
            .filter(|line| {
                in_body |= line.starts_with("#CHROM");
                (in_body || line.starts_with("##INFO") || line.starts_with("##FORMAT"))
                    && !line.is_empty()
            })
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(body, section(&text, "out", case), "{case}");
        // The bamout, as SAM text, and its index.
        let bytes = std::fs::read(&bamout).expect("the bamout");
        let raw = htsjdk_bgzf::read::decompress_all(&bytes).expect("a BGZF file");
        let mut reader = htsjdk_bam::reader::BamReader::new(&raw).expect("a BAM");
        let header = reader.header.text.clone();
        let records: Vec<_> = reader.by_ref().map(|r| r.expect("a record")).collect();
        let sam = htsjdk_bam::sam_file::write_sam(&header, &records).expect("SAM text");
        assert_eq!(sam, section(&text, "bam", case), "{case} bamout");
        let indexed = dir.join(format!("bamout-{case}.bai")).exists();
        assert_eq!(
            indexed.to_string(),
            section(&text, "bai", case),
            "{case} index"
        );
    }
}
