//! Conformance for the callers force-calling against GATK 4.6.2.0, every run of the dump replayed
//! through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcAllelesDump.java`, over `HcEndToEndDump`'s reads and
//! `HcPileupDump`'s paired set, each with a VCF of given alleles: HaplotypeCaller with --alleles,
//! filtered records forced, GVCF mode, no assembly, beside pileup detection and the PDHMM, under
//! intervals, RampedHaplotypeCaller, and Mutect2 tumour-only, tumour-normal and in GVCF mode. Each
//! output is compared from its ##INFO and ##FORMAT lines to its last record.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_alleles.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-alleles");
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
    let (paired_reference, paired_bam) = set("paired", &["chr1"]);
    let given = |name: &str| {
        let path = dir.join(name).join("given.vcf");
        std::fs::write(&path, section(&text, "vcf", name)).expect("the given alleles");
        path.to_string_lossy().to_string()
    };
    let e2e_given = given("e2e");
    let paired_given = given("paired");
    let e2e = (e2e_bam.as_str(), e2e_reference.as_str());
    let paired = (paired_bam.as_str(), paired_reference.as_str());
    let alleles = "--alleles";
    let eg = e2e_given.as_str();
    let pg = paired_given.as_str();
    let runs: Vec<Run<'_>> = vec![
        ("hc-alleles", "HaplotypeCaller", e2e, vec![alleles, eg]),
        (
            "hc-alleles-filtered",
            "HaplotypeCaller",
            e2e,
            vec![alleles, eg, "--force-call-filtered-alleles"],
        ),
        (
            "hc-alleles-gvcf",
            "HaplotypeCaller",
            e2e,
            vec![
                alleles,
                eg,
                "-ERC",
                "GVCF",
                "--sample-name",
                "s1",
                "-L",
                "chr1:250-750",
            ],
        ),
        (
            "hc-alleles-no-assembly",
            "HaplotypeCaller",
            e2e,
            vec![
                alleles,
                eg,
                "--kmer-size",
                "101",
                "--dont-increase-kmer-sizes-for-cycles",
            ],
        ),
        (
            "hc-alleles-pileup-pdhmm",
            "HaplotypeCaller",
            e2e,
            vec![alleles, eg, "--pileup-detection", "--use-pdhmm"],
        ),
        (
            "hc-alleles-intervals",
            "HaplotypeCaller",
            e2e,
            vec![alleles, eg, "-L", "chr1:400-470", "-L", "chr2:100-140"],
        ),
        (
            "ramped-alleles",
            "RampedHaplotypeCaller",
            e2e,
            vec![alleles, eg],
        ),
        ("mutect2-alleles", "Mutect2", e2e, vec![alleles, eg]),
        (
            "mutect2-alleles-filtered-tn",
            "Mutect2",
            e2e,
            vec![
                alleles,
                eg,
                "--force-call-filtered-alleles",
                "--normal-sample",
                "s2",
            ],
        ),
        (
            "mutect2-alleles-gvcf",
            "Mutect2",
            e2e,
            vec![alleles, eg, "-ERC", "GVCF", "-L", "chr2"],
        ),
        (
            "paired-hc-alleles",
            "HaplotypeCaller",
            paired,
            vec![alleles, pg],
        ),
        (
            "paired-mutect2-alleles",
            "Mutect2",
            paired,
            vec![alleles, pg],
        ),
    ];
    for (case, tool, (bam, reference), extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
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
    }
}
