//! Conformance for `--flow-filter-alleles` in the callers against GATK 4.6.2.0, every run of the
//! dump replayed through the command line, in its order.
//!
//! Golden from `tools/readfilter-conformance/HcAlleleFilteringDump.java`, over `HcEndToEndDump`'s
//! reads: the default thresholds, a quality threshold of 300 and of 2000, a strand odds ratio
//! threshold of 0.5, lone alleles removed, RampedHaplotypeCaller, and GVCF mode on one sample.
//! Each output is compared from its ##INFO and ##FORMAT lines to its last record.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_allele_filtering.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-allele-filtering");
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
    let contigs = vec![
        ("chr1", section(&text, "fasta", "chr1")),
        ("chr2", section(&text, "fasta", "chr2")),
    ];
    let reference = write_reference(&dir, &contigs);
    let bam = write_bam(&dir, &section(&text, "sam", "input"));
    let filter = "--flow-filter-alleles";
    let qual = "--flow-filter-alleles-qual-threshold";
    let runs: Vec<(&str, &str, Vec<&str>)> = vec![
        ("hc-filter", "HaplotypeCaller", vec![filter]),
        (
            "hc-filter-qual-300",
            "HaplotypeCaller",
            vec![filter, qual, "300"],
        ),
        (
            "hc-filter-qual-2000",
            "HaplotypeCaller",
            vec![filter, qual, "2000"],
        ),
        (
            "hc-filter-sor",
            "HaplotypeCaller",
            vec![filter, "--flow-filter-alleles-sor-threshold", "0.5"],
        ),
        (
            "hc-filter-lone",
            "HaplotypeCaller",
            vec![filter, "--flow-filter-lone-alleles", qual, "300"],
        ),
        (
            "ramped-filter",
            "RampedHaplotypeCaller",
            vec![filter, qual, "300"],
        ),
        (
            "hc-filter-gvcf",
            "HaplotypeCaller",
            vec![filter, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"],
        ),
        ("mutect2-filter", "Mutect2", vec![filter]),
        (
            "mutect2-filter-qual-300",
            "Mutect2",
            vec![filter, qual, "300", "--normal-sample", "s2"],
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
