//! Conformance for the callers with the flow-based likelihood engine against GATK 4.6.2.0, every
//! run of the dump replayed through the command line.
//!
//! Golden from `tools/readfilter-conformance/HcFlowDump.java`, over one contig and one sample of
//! Ultima reads (a het SNP, a het deletion and a het insertion in homopolymers): HaplotypeCaller
//! with `--likelihood-calculation-engine FlowBased`, its optimized computation, exact matching,
//! untrimmed reads, the flow matrix arguments, dynamic read disqualification and GVCF mode, and
//! Mutect2 with the same engine. Each output is compared from its header lines to its last record.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_flow.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-flow");
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
    let contigs = vec![("chr1", section(&text, "fasta", "chr1"))];
    let reference = write_reference(&dir, &contigs);
    let bam = write_bam(&dir, &section(&text, "sam", "input"));
    let flow = ["--likelihood-calculation-engine", "FlowBased"];
    let collapse = "--flow-assembly-collapse-hmer-size";
    let hmm = ["--likelihood-calculation-engine", "FlowBasedHMM"];
    let stepwise = "--use-flow-aligner-for-stepwise-hc-filtering";
    let annotations = [
        "-A",
        "CycleSkipStatus",
        "-A",
        "HmerIndelLength",
        "-A",
        "HmerMotifs",
        "-A",
        "VariantType",
        "-A",
        "GcContent",
        "-A",
        "IndelClassify",
    ];
    let with = |extra: &[&'static str]| -> Vec<&'static str> {
        flow.iter().chain(extra).copied().collect()
    };
    let runs: Vec<(&str, &str, Vec<&str>)> = vec![
        ("hc-flow", "HaplotypeCaller", with(&[])),
        (
            "hc-flow-optimized",
            "HaplotypeCaller",
            with(&["--flow-likelihood-optimized-comp"]),
        ),
        (
            "hc-flow-exact",
            "HaplotypeCaller",
            with(&["--exact-matching", "--trim-to-haplotype", "false"]),
        ),
        (
            "hc-flow-untrimmed",
            "HaplotypeCaller",
            with(&["--trim-to-haplotype", "false"]),
        ),
        (
            "hc-flow-matrix",
            "HaplotypeCaller",
            with(&[
                "--flow-use-t0-tag",
                "--flow-fill-empty-bins-value",
                "0.0001",
                "--keep-boundary-flows",
                "--flow-symmetric-indel-probs",
                "--flow-lump-probs",
            ]),
        ),
        (
            "hc-flow-dynamic",
            "HaplotypeCaller",
            with(&[
                "--enable-dynamic-read-disqualification-for-genotyping",
                "--dynamic-read-disqualification-threshold",
                "10",
            ]),
        ),
        (
            "hc-flow-gvcf",
            "HaplotypeCaller",
            with(&["-ERC", "GVCF", "-L", "chr1:550-700"]),
        ),
        ("mutect2-flow", "Mutect2", with(&[])),
        ("hc-collapse-4", "HaplotypeCaller", vec![collapse, "4"]),
        (
            "hc-collapse-3-partial",
            "HaplotypeCaller",
            vec![collapse, "3", "--flow-assembly-collapse-partial-mode"],
        ),
        (
            "hc-collapse-determined",
            "HaplotypeCaller",
            vec![collapse, "-1"],
        ),
        (
            "hc-flow-collapse-4",
            "HaplotypeCaller",
            with(&[collapse, "4"]),
        ),
        (
            "ramped-collapse-4",
            "RampedHaplotypeCaller",
            vec![collapse, "4"],
        ),
        ("mutect2-collapse-4", "Mutect2", vec![collapse, "4"]),
        (
            "hc-flow-mode-standard",
            "HaplotypeCaller",
            vec!["--flow-mode", "STANDARD"],
        ),
        (
            "hc-flow-mode-advanced",
            "HaplotypeCaller",
            vec!["--flow-mode", "ADVANCED"],
        ),
        (
            "hc-flow-mode-standard-sor-1",
            "HaplotypeCaller",
            vec![
                "--flow-mode",
                "STANDARD",
                "--flow-filter-alleles-sor-threshold",
                "1",
            ],
        ),
        (
            "ramped-flow-mode-advanced",
            "RampedHaplotypeCaller",
            vec!["--flow-mode", "ADVANCED"],
        ),
        (
            "mutect2-flow-mode-standard",
            "Mutect2",
            vec!["--flow-mode", "STANDARD"],
        ),
        (
            "mutect2-flow-mode-advanced",
            "Mutect2",
            vec!["--flow-mode", "ADVANCED"],
        ),
        ("hc-flow-hmm", "HaplotypeCaller", hmm.to_vec()),
        (
            "hc-flow-hmm-penalties",
            "HaplotypeCaller",
            hmm.iter()
                .copied()
                .chain([
                    "--flow-hmm-engine-min-indel-adjust",
                    "20",
                    "--flow-hmm-engine-flat-insertion-penatly",
                    "30",
                    "--flow-hmm-engine-flat-deletion-penatly",
                    "35",
                    "--pcr-indel-model",
                    "AGGRESSIVE",
                ])
                .collect(),
        ),
        (
            "hc-flow-hmm-no-pcr",
            "HaplotypeCaller",
            hmm.iter()
                .copied()
                .chain([
                    "--pcr-indel-model",
                    "NONE",
                    "--enable-dynamic-read-disqualification-for-genotyping",
                ])
                .collect(),
        ),
        (
            "hc-flow-hmm-cycle-8",
            "HaplotypeCaller",
            hmm.iter()
                .copied()
                .chain(["--flow-order-cycle-length", "8"])
                .collect(),
        ),
        ("mutect2-flow-hmm", "Mutect2", hmm.to_vec()),
        (
            "hc-stepwise",
            "HaplotypeCaller",
            vec![
                stepwise,
                "--flow-filter-alleles",
                "--flow-filter-alleles-sor-threshold",
                "1",
            ],
        ),
        (
            "hc-stepwise-hmm",
            "HaplotypeCaller",
            hmm.iter()
                .copied()
                .chain([stepwise, "--flow-filter-alleles"])
                .collect(),
        ),
        (
            "ramped-stepwise",
            "RampedHaplotypeCaller",
            vec![stepwise, "--flow-filter-alleles"],
        ),
        ("hc-stepwise-no-filter", "HaplotypeCaller", vec![stepwise]),
        (
            "hc-stepwise-flow-based",
            "HaplotypeCaller",
            with(&[stepwise, "--flow-filter-alleles"]),
        ),
        (
            "hc-flow-annotations",
            "HaplotypeCaller",
            annotations.to_vec(),
        ),
        (
            "hc-flow-annotations-order",
            "HaplotypeCaller",
            annotations
                .iter()
                .copied()
                .chain(["--flow-order-for-annotations", "ACGT"])
                .collect(),
        ),
        (
            "hc-flow-annotations-group-order",
            "HaplotypeCaller",
            annotations
                .iter()
                .copied()
                .chain([
                    "--flow-order-for-annotations",
                    "rgF:CATG",
                    "--flow-order-for-annotations",
                    "GTAC",
                ])
                .collect(),
        ),
        (
            "hc-flow-annotations-only-group",
            "HaplotypeCaller",
            annotations
                .iter()
                .copied()
                .chain(["--flow-order-for-annotations", "rgF:CATG"])
                .collect(),
        ),
        (
            "hc-flow-order-without-annotation",
            "HaplotypeCaller",
            vec!["--flow-order-for-annotations", "ACGT"],
        ),
        (
            "mutect2-flow-annotations-order",
            "Mutect2",
            annotations
                .iter()
                .copied()
                .chain(["--flow-order-for-annotations", "TACG"])
                .collect(),
        ),
    ];
    for (case, tool, extra) in runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            tool,
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
