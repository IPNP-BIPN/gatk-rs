//! Conformance for HaplotypeCaller's variant resources against GATK 4.6.2.0, every run of the dump
//! replayed through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcResourcesDump.java`, over `HcEndToEndDump`'s reads
//! and resource VCFs written from a first plain run's calls: `--dbsnp` alone, beside comparisons
//! (untagged, tagged, the dbSNP file itself, and tagged `DB`, which is refused), in GVCF mode, and
//! `--population-callset` in VCF, GVCF and base-pair-resolution mode, with and without reference
//! samples assumed, and RampedHaplotypeCaller with all three. Each output is compared from its
//! ##INFO and ##FORMAT lines to its last record; an untagged comparison is named by its absolute
//! path, so the oracle's directory is read as this test's.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_resources.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-resources");
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
    let e2e_contigs: Vec<(&str, String)> = ["chr1", "chr2"]
        .iter()
        .map(|&contig| (contig, section(&text, "fasta", &format!("e2e:{contig}"))))
        .collect();
    let sub = dir.join("e2e");
    std::fs::create_dir_all(&sub).expect("a set directory");
    let e2e_reference = write_reference(&sub, &e2e_contigs);
    let e2e_bam = write_bam(&sub, &section(&text, "sam", "e2e"));
    let e2e = (e2e_bam.as_str(), e2e_reference.as_str());
    let resource = |name: &str| {
        let path = dir.join(format!("{name}.vcf"));
        std::fs::write(&path, section(&text, "vcf", name)).expect("a resource");
        path.to_string_lossy().to_string()
    };
    let (db, co, po) = (resource("dbsnp"), resource("comp"), resource("population"));
    let (db, co, po) = (db.as_str(), co.as_str(), po.as_str());
    let gvcf = ["-ERC", "GVCF", "--sample-name", "s1", "-L", "chr1"];
    let runs: Vec<Run<'_>> = vec![
        ("hc-plain", "HaplotypeCaller", e2e, vec![]),
        ("hc-dbsnp", "HaplotypeCaller", e2e, vec!["--dbsnp", db]),
        (
            "hc-dbsnp-comp",
            "HaplotypeCaller",
            e2e,
            vec!["--dbsnp", db, "--comparison", co, "--comparison:truth", db],
        ),
        (
            "hc-comp-is-dbsnp",
            "HaplotypeCaller",
            e2e,
            vec!["--dbsnp", db, "--comparison", db],
        ),
        (
            "hc-comp-tagged-db",
            "HaplotypeCaller",
            e2e,
            vec!["--comparison:DB", co],
        ),
        (
            "hc-dbsnp-gvcf",
            "HaplotypeCaller",
            e2e,
            [&["--dbsnp", db][..], &gvcf[..]].concat(),
        ),
        (
            "hc-population-vcf",
            "HaplotypeCaller",
            e2e,
            vec!["--population-callset", po],
        ),
        (
            "hc-population-gvcf",
            "HaplotypeCaller",
            e2e,
            [&["--population-callset", po][..], &gvcf[..]].concat(),
        ),
        (
            "hc-population-gvcf-numref",
            "HaplotypeCaller",
            e2e,
            vec![
                "--population-callset",
                po,
                "--num-reference-samples-if-no-call",
                "5",
                "-ERC",
                "GVCF",
                "--sample-name",
                "s2",
                "-L",
                "chr1",
            ],
        ),
        (
            "hc-population-bp",
            "HaplotypeCaller",
            e2e,
            vec![
                "--population-callset",
                po,
                "-ERC",
                "BP_RESOLUTION",
                "--sample-name",
                "s1",
                "-L",
                "chr1:250-450",
            ],
        ),
        (
            "ramped-resources",
            "RampedHaplotypeCaller",
            e2e,
            [
                &[
                    "--dbsnp",
                    db,
                    "--comparison:truth",
                    co,
                    "--population-callset",
                    po,
                ][..],
                &gvcf[..],
            ]
            .concat(),
        ),
    ];
    let oracle_dir = "/work/hc-resources-dump/";
    let our_dir = format!("{}/", dir.to_string_lossy());
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
        assert_eq!(
            body,
            section(&text, "out", case).replace(oracle_dir, &our_dir),
            "{case}"
        );
    }
}
