//! Conformance for the callers with the linked de Bruijn graph against GATK 4.6.2.0, every run of
//! the dump replayed through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcLinkedGraphDump.java`, over `HcEndToEndDump`'s reads
//! and `HcDragenEndToEndDump`'s: HaplotypeCaller, RampedHaplotypeCaller and Mutect2 with
//! --linked-de-bruijn-graph, with a fixed k-mer size, every dangling branch recovered and in GVCF
//! mode, and --disable-artificial-haplotype-recovery, which the reference refuses in that mode.
//! Each output is compared from its ##INFO and ##FORMAT lines to its last record.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_linked_graph.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-linked-graph");
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
    let (str_reference, str_bam) = set("str", &["chr1"]);
    let linked = "--linked-de-bruijn-graph";
    let e2e = (e2e_bam.as_str(), e2e_reference.as_str());
    let strs = (str_bam.as_str(), str_reference.as_str());
    let runs: Vec<Run<'_>> = vec![
        ("hc-linked", "HaplotypeCaller", e2e, vec![linked]),
        (
            "hc-linked-no-recovery",
            "HaplotypeCaller",
            e2e,
            vec![linked, "--disable-artificial-haplotype-recovery"],
        ),
        (
            "hc-recovery-without-linked",
            "HaplotypeCaller",
            e2e,
            vec!["--disable-artificial-haplotype-recovery"],
        ),
        (
            "hc-linked-kmer-15",
            "HaplotypeCaller",
            e2e,
            vec![
                linked,
                "--kmer-size",
                "15",
                "--dont-increase-kmer-sizes-for-cycles",
            ],
        ),
        (
            "hc-linked-recover-all",
            "HaplotypeCaller",
            e2e,
            vec![linked, "--recover-all-dangling-branches"],
        ),
        (
            "hc-linked-gvcf",
            "HaplotypeCaller",
            e2e,
            vec![linked, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"],
        ),
        ("ramped-linked", "RampedHaplotypeCaller", e2e, vec![linked]),
        ("mutect2-linked", "Mutect2", e2e, vec![linked]),
        (
            "mutect2-linked-tumor-normal",
            "Mutect2",
            e2e,
            vec![linked, "--normal-sample", "s2"],
        ),
        ("str-hc-linked", "HaplotypeCaller", strs, vec![linked]),
        (
            "str-hc-linked-no-recovery",
            "HaplotypeCaller",
            strs,
            vec![linked, "--disable-artificial-haplotype-recovery"],
        ),
        ("str-mutect2-linked", "Mutect2", strs, vec![linked]),
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
            // `PedReader`'s source is an identity hash, masked in the golden.
            for part in message.split("<hash>") {
                assert!(outcome.stderr.contains(part), "{case}: {}", outcome.stderr);
            }
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
