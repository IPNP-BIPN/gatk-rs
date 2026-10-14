//! Conformance for the callers with the partially determined PairHMM against GATK 4.6.2.0, every
//! run of the dump replayed through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcPdhmmDump.java`, over `HcPileupDump`'s read sets:
//! --use-pdhmm with and without pileup detection, the overlap optimisation, determined haplotypes,
//! the GGA fallback, GVCF mode (which the reference refuses for want of a reference haplotype),
//! --dragen-378-concordance-mode, RampedHaplotypeCaller and Mutect2 (which make no PD haplotypes)
//! the results file and the debug status printed. Each output is compared from its ##INFO and ##FORMAT lines to its last
//! record, and each results file whole.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_pdhmm.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-pdhmm");
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
    let (paired_reference, paired_bam) = set("paired", &["chr1"]);
    let e2e = (e2e_bam.as_str(), e2e_reference.as_str());
    let strs = (str_bam.as_str(), str_reference.as_str());
    let paired = (paired_bam.as_str(), paired_reference.as_str());
    let on = "--pileup-detection";
    let indels = "--pileup-detection-enable-indel-pileup-calling";
    let pdhmm = "--use-pdhmm";
    let dragen = "--dragen-378-concordance-mode";
    let status = "--print-pileupcalling-status";
    let no_assembly = [
        "--kmer-size",
        "101",
        "--dont-increase-kmer-sizes-for-cycles",
    ];
    let with = |parts: &[&[&'static str]]| -> Vec<&'static str> { parts.concat() };
    let results = |case: &str| {
        dir.join(format!("results-{case}.txt"))
            .to_string_lossy()
            .to_string()
    };
    let first_results = results("hc-pdhmm-results");
    let second_results = results("hc-pdhmm-overlap-results");
    let runs: Vec<Run<'_>> = vec![
        ("hc-pdhmm", "HaplotypeCaller", e2e, vec![pdhmm]),
        (
            "hc-pdhmm-pileup",
            "HaplotypeCaller",
            e2e,
            vec![pdhmm, on, indels],
        ),
        (
            "hc-pdhmm-overlap",
            "HaplotypeCaller",
            e2e,
            vec![pdhmm, on, indels, "--use-pdhmm-overlap-optimization"],
        ),
        (
            "hc-pdhmm-determined",
            "HaplotypeCaller",
            e2e,
            vec![pdhmm, on, indels, "--make-determined-haps-from-pd-code"],
        ),
        (
            "hc-pdhmm-no-assembly",
            "HaplotypeCaller",
            e2e,
            with(&[&[pdhmm, on, indels], &no_assembly]),
        ),
        (
            "hc-pdhmm-no-fallback",
            "HaplotypeCaller",
            e2e,
            with(&[
                &[pdhmm, on, indels, "--fallback-gga-if-pdhmm-fails", "false"],
                &no_assembly,
            ]),
        ),
        (
            "hc-pdhmm-gvcf",
            "HaplotypeCaller",
            e2e,
            vec![pdhmm, "-ERC", "GVCF", "--sample-name", "s1", "-L", "chr2"],
        ),
        ("hc-dragen-378", "HaplotypeCaller", e2e, vec![dragen]),
        (
            "hc-dragen-378-gvcf",
            "HaplotypeCaller",
            e2e,
            vec![
                dragen,
                "-ERC",
                "GVCF",
                "--sample-name",
                "s2",
                "-L",
                "chr1:1150-1300",
            ],
        ),
        (
            "hc-dragen-378-with-3412",
            "HaplotypeCaller",
            e2e,
            vec![dragen, "--dragen-mode"],
        ),
        (
            "ramped-pdhmm",
            "RampedHaplotypeCaller",
            e2e,
            vec![pdhmm, on, indels],
        ),
        ("mutect2-pdhmm", "Mutect2", e2e, vec![pdhmm, on, indels]),
        (
            "hc-pdhmm-results",
            "HaplotypeCaller",
            e2e,
            vec![
                pdhmm,
                on,
                indels,
                "-L",
                "chr2:360-440",
                "--pdhmm-results-file",
                &first_results,
            ],
        ),
        (
            "hc-pdhmm-overlap-results",
            "HaplotypeCaller",
            e2e,
            vec![
                pdhmm,
                on,
                indels,
                "--use-pdhmm-overlap-optimization",
                "-L",
                "chr1:1190-1240",
                "--pdhmm-results-file",
                &second_results,
            ],
        ),
        (
            "str-hc-pdhmm",
            "HaplotypeCaller",
            strs,
            vec![pdhmm, on, indels],
        ),
        ("str-hc-dragen-378", "HaplotypeCaller", strs, vec![dragen]),
        (
            "str-hc-pdhmm-determined",
            "HaplotypeCaller",
            strs,
            vec![pdhmm, on, indels, "--make-determined-haps-from-pd-code"],
        ),
        (
            "paired-hc-pdhmm",
            "HaplotypeCaller",
            paired,
            vec![pdhmm, on, indels],
        ),
        (
            "paired-hc-dragen-378",
            "HaplotypeCaller",
            paired,
            vec![dragen],
        ),
        (
            "paired-hc-dragen-378-no-assembly",
            "HaplotypeCaller",
            paired,
            with(&[&[dragen], &no_assembly]),
        ),
        (
            "hc-pdhmm-status",
            "HaplotypeCaller",
            e2e,
            vec![pdhmm, on, indels, status],
        ),
        (
            "hc-pdhmm-determined-status",
            "HaplotypeCaller",
            e2e,
            vec![
                pdhmm,
                on,
                indels,
                "--make-determined-haps-from-pd-code",
                status,
            ],
        ),
        (
            "hc-pileup-status",
            "HaplotypeCaller",
            e2e,
            with(&[&[on, indels, status], &no_assembly]),
        ),
        (
            "str-hc-pdhmm-status",
            "HaplotypeCaller",
            strs,
            vec![pdhmm, on, indels, status],
        ),
        (
            "paired-hc-dragen-378-status",
            "HaplotypeCaller",
            paired,
            vec![dragen, status],
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
        // What `--print-pileupcalling-status` printed, where the dump captured it.
        if let Some(expected) = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("stdout\t{case}=")))
        {
            assert_eq!(outcome.stdout, unescape(expected), "{case} stdout");
        }
        // The results file, whole, where the run wrote one.
        if let Some(expected) = text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("file\t{case}=")))
        {
            let written = std::fs::read_to_string(results(case)).expect("the results file");
            assert_eq!(written, unescape(expected), "{case} results");
        }
    }
}
