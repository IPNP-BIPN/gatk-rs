//! Conformance for `--flow-filter-alleles-debug-graphs` in HaplotypeCaller, RampedHaplotypeCaller
//! and Mutect2 against GATK 4.6.2.0, every run of the dump replayed through the command line.
//!
//! Golden from `tools/readfilter-conformance/HcAlleleGraphsDump.java`, over reads tiling four
//! haplotypes: a linked group of two SNVs and one of three. Each output is compared from its
//! ##INFO and ##FORMAT lines to its last record, and each interaction graph the run leaves in the
//! working directory as written.

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_allele_graphs.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-allele-graphs");
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
    let reference = write_reference(&dir, &[("chr1", section(&text, "fasta", "chr1"))]);
    let bam = write_bam(&dir, &section(&text, "sam", "input"));
    // The graphs land in the working directory, as the reference writes them.
    std::env::set_current_dir(&dir).expect("the scratch directory");
    let filter = "--flow-filter-alleles";
    let graphs = "--flow-filter-alleles-debug-graphs";
    let runs: Vec<(&str, &str, Vec<&str>)> = vec![
        ("hc-graphs", "HaplotypeCaller", vec![filter, graphs]),
        (
            "hc-graphs-qual-300",
            "HaplotypeCaller",
            vec![
                filter,
                graphs,
                "--flow-filter-alleles-qual-threshold",
                "300",
            ],
        ),
        ("hc-graphs-unfiltered", "HaplotypeCaller", vec![graphs]),
        (
            "ramped-graphs",
            "RampedHaplotypeCaller",
            vec![filter, graphs],
        ),
        ("mutect2-graphs", "Mutect2", vec![filter, graphs]),
    ];
    for (case, tool, extra) in &runs {
        let output = dir.join(format!("out-{case}.vcf"));
        let mut args: Vec<String> = [
            *tool,
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
        let mut found = std::collections::BTreeMap::new();
        for entry in std::fs::read_dir(&dir).expect("the directory") {
            let path = entry.expect("an entry").path();
            let name = path
                .file_name()
                .expect("a name")
                .to_string_lossy()
                .to_string();
            if name.starts_with("allele.interaction.") && name.ends_with(".dot") {
                found.insert(name, std::fs::read_to_string(&path).expect("a graph"));
                std::fs::remove_file(&path).expect("removed");
            }
        }
        let expected: std::collections::BTreeMap<String, String> = text
            .lines()
            .filter_map(|line| line.strip_prefix(&format!("graph\t{case}\t")))
            .map(|rest| {
                let (name, content) = rest.split_once('=').expect("a graph");
                (name.to_string(), unescape(content))
            })
            .collect();
        assert_eq!(found, expected, "{case}");
    }
}
