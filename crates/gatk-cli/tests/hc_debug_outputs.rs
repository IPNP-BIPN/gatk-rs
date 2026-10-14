//! Conformance for the callers' debug outputs and switches against GATK 4.6.2.0, every run of the
//! dump replayed through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcDebugOutputsDump.java`, over `HcEndToEndDump`'s
//! reads: --haplotype-debug-histogram-output, --debug-assembly-variants-out, --alm-path,
//! --pair-hmm-results-file,
//! --debug-assembly-region-state and --assembly-region-out in each caller, --keep-rg, --just-determine-active-regions and
//! --dont-genotype. Each VCF is compared from its ##INFO and ##FORMAT lines to its last record, and
//! each debug file whole.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_debug_outputs.txt.gz"),
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

/// The dot files a run wrote in the working directory against the golden's, then removed. Vertex
/// and tree-node numbers are identities in the reference, so they are masked; a linked graph's
/// junction trees come out of identity-keyed hash maps, so a file holding them is compared as a
/// set of lines.
fn compare_dot_files(text: &str, case: &str, dir: &std::path::Path) {
    let prefix = format!("dot\t{case}\t");
    let expected: Vec<(&str, String)> = text
        .lines()
        .filter_map(|line| line.strip_prefix(&prefix))
        .map(|rest| {
            let (name, body) = rest.split_once('=').expect("a name and a body");
            (name, unescape(body))
        })
        .collect();
    let mut written: Vec<String> = std::fs::read_dir(dir)
        .expect("the scratch directory")
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.file_name().to_string_lossy().to_string())
        .filter(|name| name.ends_with(".dot") && !name.starts_with("debug-"))
        .collect();
    written.sort();
    let names: Vec<&str> = expected.iter().map(|(name, _)| *name).collect();
    assert_eq!(written, names, "{case} dot files");
    for (name, body) in expected {
        let path = dir.join(name);
        let ours = std::fs::read_to_string(&path).expect("a dot file");
        assert_dot_eq(&ours, &body, case, name);
        std::fs::remove_file(path).expect("removable");
    }
}

/// Two dot files, identities masked, compared as sets of lines where junction trees are printed.
fn assert_dot_eq(ours: &str, theirs: &str, case: &str, name: &str) {
    let (ours, theirs) = (mask_dot(ours), mask_dot(theirs));
    if theirs.contains("TreadingNode_") {
        let mut ours: Vec<&str> = ours.lines().collect();
        let mut theirs: Vec<&str> = theirs.lines().collect();
        ours.sort_unstable();
        theirs.sort_unstable();
        assert_eq!(ours, theirs, "{case} {name}");
    } else {
        assert_eq!(ours, theirs, "{case} {name}");
    }
}

/// Masks `_id_<n>_`, `(<n>)` and `TreadingNode_<hex>`.
fn mask_dot(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    loop {
        let next = [("_id_", 0), ("(", 1), ("TreadingNode_", 2)]
            .iter()
            .filter_map(|&(marker, kind)| rest.find(marker).map(|at| (at, marker, kind)))
            .min();
        let Some((at, marker, kind)) = next else {
            out.push_str(rest);
            return out;
        };
        out.push_str(&rest[..at + marker.len()]);
        rest = &rest[at + marker.len()..];
        let digits = match kind {
            2 => {
                rest.len()
                    - rest
                        .trim_start_matches(|c: char| c.is_ascii_hexdigit())
                        .len()
            }
            _ => {
                let unsigned = rest.strip_prefix('-').unwrap_or(rest);
                let count = unsigned.len()
                    - unsigned
                        .trim_start_matches(|c: char| c.is_ascii_digit())
                        .len();
                let sign = rest.len() - unsigned.len();
                // Only a number closed by the marker's end: `_id_<n>_` and `(<n>)`.
                let close = match kind {
                    0 => '_',
                    _ => ')',
                };
                if count > 0 && unsigned[count..].starts_with(close) {
                    sign + count
                } else {
                    0
                }
            }
        };
        if digits > 0 {
            out.push('N');
        }
        rest = &rest[digits..];
    }
}

/// FRD's debug lines print three `double[]` by identity (`[D@1b2c3d4`), which no two JVM runs
/// share; only the prefix is compared.
fn mask_identities(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.find("[D@") {
        out.push_str(&rest[..at + 3]);
        rest = rest[at + 3..].trim_start_matches(|c: char| c.is_ascii_hexdigit());
    }
    out.push_str(rest);
    out
}

fn scratch() -> std::path::PathBuf {
    let dir = std::env::temp_dir().join("gatk-cli-hc-debug-outputs");
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
    // `{file}` stands for the run's debug file.
    let region_out = "--assembly-region-out";
    let histogram = "--haplotype-debug-histogram-output";
    let variants_out = "--debug-assembly-variants-out";
    let alm = "--alm-path";
    let hmm_results = "--pair-hmm-results-file";
    let state = "--debug-assembly-region-state";
    let genotyper = "--debug-genotyper-output";
    let graph = "--graph-output";
    let transformations = "--debug-graph-transformations";
    let disqualify = ["--expected-mismatch-rate-for-read-disqualification", "0.0"];
    let gvcf = ["-ERC", "GVCF", "--sample-name", "s1", "-L", "chr2"];
    let runs: Vec<Run<'_>> = vec![
        (
            "hc-region-out",
            "HaplotypeCaller",
            e2e,
            vec![region_out, "{file}"],
        ),
        (
            "hc-region-out-gvcf",
            "HaplotypeCaller",
            e2e,
            [&gvcf[..], &[region_out, "{file}"]].concat(),
        ),
        (
            "ramped-region-out",
            "RampedHaplotypeCaller",
            e2e,
            vec!["-L", "chr2", region_out, "{file}"],
        ),
        (
            "mutect2-region-out",
            "Mutect2",
            e2e,
            vec![region_out, "{file}"],
        ),
        (
            "hc-histogram",
            "HaplotypeCaller",
            e2e,
            vec![histogram, "{file}"],
        ),
        (
            "hc-histogram-linked",
            "HaplotypeCaller",
            e2e,
            vec!["--linked-de-bruijn-graph", histogram, "{file}"],
        ),
        (
            "hc-histogram-empty",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr2:10-20", histogram, "{file}"],
        ),
        (
            "ramped-histogram",
            "RampedHaplotypeCaller",
            e2e,
            vec![histogram, "{file}"],
        ),
        (
            "mutect2-histogram",
            "Mutect2",
            e2e,
            vec![histogram, "{file}"],
        ),
        (
            "hc-assembly-variants",
            "HaplotypeCaller",
            e2e,
            vec![variants_out, "{file}"],
        ),
        (
            "hc-assembly-variants-chr2",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr2", variants_out, "{file}"],
        ),
        (
            "ramped-assembly-variants",
            "RampedHaplotypeCaller",
            e2e,
            vec!["-L", "chr2", variants_out, "{file}"],
        ),
        (
            "mutect2-assembly-variants",
            "Mutect2",
            e2e,
            vec![variants_out, "{file}"],
        ),
        (
            "mutect2-assembly-variants-chr2",
            "Mutect2",
            e2e,
            vec!["-L", "chr2", variants_out, "{file}"],
        ),
        (
            "hc-alm",
            "HaplotypeCaller",
            e2e,
            vec!["--alm-interval", "chr2", "-L", "chr2", alm, "{file}"],
        ),
        (
            "hc-alm-window",
            "HaplotypeCaller",
            e2e,
            vec!["--alm-interval", "chr1:1000-1400", alm, "{file}"],
        ),
        (
            "hc-alm-no-interval",
            "HaplotypeCaller",
            e2e,
            vec![alm, "{file}"],
        ),
        (
            "hc-alm-pdhmm",
            "HaplotypeCaller",
            e2e,
            vec![
                "--alm-interval",
                "chr1:1100-1350",
                "--use-pdhmm",
                "-L",
                "chr1:1150-1300",
                alm,
                "{file}",
            ],
        ),
        (
            "ramped-alm",
            "RampedHaplotypeCaller",
            e2e,
            vec!["--alm-interval", "chr2", "-L", "chr2", alm, "{file}"],
        ),
        (
            "mutect2-alm",
            "Mutect2",
            e2e,
            vec!["--alm-interval", "chr2", "-L", "chr2", alm, "{file}"],
        ),
        (
            "hc-pairhmm-results",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr2:380-420", hmm_results, "{file}"],
        ),
        (
            "ramped-pairhmm-results",
            "RampedHaplotypeCaller",
            e2e,
            vec!["-L", "chr2:380-420", hmm_results, "{file}"],
        ),
        (
            "mutect2-pairhmm-results",
            "Mutect2",
            e2e,
            vec!["-L", "chr2:380-420", hmm_results, "{file}"],
        ),
        (
            "hc-region-state",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr2", state, "{file}"],
        ),
        (
            "hc-region-state-filter",
            "HaplotypeCaller",
            e2e,
            vec![
                "-L",
                "chr1:1100-1300",
                "--flow-filter-alleles",
                state,
                "{file}",
            ],
        ),
        (
            "ramped-region-state",
            "RampedHaplotypeCaller",
            e2e,
            vec!["-L", "chr2", state, "{file}"],
        ),
        (
            "hc-genotyper",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr2", genotyper, "{file}"],
        ),
        (
            "hc-genotyper-cis",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr1:1150-1300", genotyper, "{file}"],
        ),
        (
            "hc-genotyper-dynamic",
            "HaplotypeCaller",
            e2e,
            vec![
                "-L",
                "chr1:900-1050",
                "--enable-dynamic-read-disqualification-for-genotyping",
                genotyper,
                "{file}",
            ],
        ),
        (
            "hc-genotyper-legacy",
            "HaplotypeCaller",
            e2e,
            vec![
                "-L",
                "chr1:700-800",
                "--enable-legacy-assembly-region-trimming",
                genotyper,
                "{file}",
            ],
        ),
        (
            "hc-genotyper-dragen",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr1:900-1050", "--dragen-mode", genotyper, "{file}"],
        ),
        (
            "hc-genotyper-pdhmm",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr1:1150-1300", "--use-pdhmm", genotyper, "{file}"],
        ),
        (
            "hc-genotyper-dragen-snp",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr1:1150-1300", "--dragen-mode", genotyper, "{file}"],
        ),
        (
            "hc-genotyper-disqualify",
            "HaplotypeCaller",
            e2e,
            [
                &["-L", "chr1:900-1050"][..],
                &disqualify[..],
                &[genotyper, "{file}"],
            ]
            .concat(),
        ),
        (
            "hc-genotyper-disqualify-dynamic",
            "HaplotypeCaller",
            e2e,
            [
                &["-L", "chr1:900-1050"][..],
                &disqualify[..],
                &[
                    "--enable-dynamic-read-disqualification-for-genotyping",
                    "--dynamic-read-disqualification-threshold",
                    "0.02",
                    genotyper,
                    "{file}",
                ],
            ]
            .concat(),
        ),
        (
            "ramped-genotyper",
            "RampedHaplotypeCaller",
            e2e,
            vec!["-L", "chr2", genotyper, "{file}"],
        ),
        ("hc-graph", "HaplotypeCaller", e2e, vec![graph, "{file}"]),
        (
            "hc-graph-debug",
            "HaplotypeCaller",
            e2e,
            vec!["-L", "chr1:900-1300", transformations, graph, "{file}"],
        ),
        (
            "hc-graph-adaptive",
            "HaplotypeCaller",
            e2e,
            vec![
                "-L",
                "chr1:900-1300",
                "--adaptive-pruning",
                transformations,
                graph,
                "{file}",
            ],
        ),
        (
            "hc-graph-linked",
            "HaplotypeCaller",
            e2e,
            vec![
                "-L",
                "chr1:900-1300",
                "--linked-de-bruijn-graph",
                "--assembly-region-padding",
                "30",
                transformations,
                graph,
                "{file}",
            ],
        ),
        (
            "hc-graph-linked-long-name",
            "HaplotypeCaller",
            e2e,
            vec![
                "-L",
                "chr2",
                "--linked-de-bruijn-graph",
                transformations,
                graph,
                "{file}",
            ],
        ),
        (
            "mutect2-graph",
            "Mutect2",
            e2e,
            vec!["-L", "chr1:1150-1300", transformations, graph, "{file}"],
        ),
        (
            "ramped-graph",
            "RampedHaplotypeCaller",
            e2e,
            vec!["-L", "chr2", transformations, graph, "{file}"],
        ),
        (
            "hc-keep-rg",
            "HaplotypeCaller",
            e2e,
            vec!["--keep-rg", "rg1"],
        ),
        (
            "ramped-keep-rg",
            "RampedHaplotypeCaller",
            e2e,
            vec!["--keep-rg", "rg2"],
        ),
        (
            "hc-just-active",
            "HaplotypeCaller",
            e2e,
            vec!["--just-determine-active-regions"],
        ),
        (
            "hc-just-active-gvcf",
            "HaplotypeCaller",
            e2e,
            [&["--just-determine-active-regions"][..], &gvcf[..]].concat(),
        ),
        (
            "ramped-just-active-gvcf",
            "RampedHaplotypeCaller",
            e2e,
            [&["--just-determine-active-regions"][..], &gvcf[..]].concat(),
        ),
        (
            "hc-dont-genotype",
            "HaplotypeCaller",
            e2e,
            vec!["--dont-genotype"],
        ),
        (
            "ramped-dont-genotype",
            "RampedHaplotypeCaller",
            e2e,
            vec!["--dont-genotype"],
        ),
        (
            "hc-dont-genotype-gvcf",
            "HaplotypeCaller",
            e2e,
            vec!["--dont-genotype", "-ERC", "GVCF"],
        ),
    ];
    // `--debug-graph-transformations` writes in the working directory, as the reference does.
    std::env::set_current_dir(&dir).expect("the scratch directory");
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
        let file = dir.join(format!("debug-{case}.txt"));
        args.extend(extra.iter().map(|arg| {
            if *arg == "{file}" {
                file.to_string_lossy().to_string()
            } else {
                arg.to_string()
            }
        }));
        let outcome = gatk_cli::run(&args);
        compare_dot_files(&text, case, &dir);
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
        match text
            .lines()
            .find_map(|line| line.strip_prefix(&format!("file\t{case}=")))
        {
            Some(expected) => {
                let written = std::fs::read_to_string(&file).expect("the debug file");
                if case.contains("graph") {
                    assert_dot_eq(&written, &unescape(expected), case, "graph output");
                } else {
                    assert_eq!(
                        mask_identities(&written),
                        mask_identities(&unescape(expected)),
                        "{case} file"
                    );
                }
            }
            None => assert!(!file.exists(), "{case} wrote a debug file"),
        }
    }
}
