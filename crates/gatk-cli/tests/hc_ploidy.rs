//! Conformance for `--ploidy-regions` in HaplotypeCaller and RampedHaplotypeCaller against
//! GATK 4.6.2.0, every run of the dump replayed through the command line, in the dump's order.
//!
//! Golden from `tools/readfilter-conformance/HcPloidyDump.java`, over `HcEndToEndDump`'s reads and
//! ploidy BEDs: haploid, triploid and tetraploid intervals, overlapping ones (the largest wins), an
//! interval of ploidy zero, one overlapping a region's span only in part, GVCF and base-pair
//! resolution, a `--sample-ploidy` beside them, and the refusals of an interval off the dictionary,
//! past a contig's end, of a name that is no integer and of a negative ploidy. Each output is
//! compared from its ##INFO and ##FORMAT lines to its last record.

use gatk_corpus as corpus;

/// A run: its label, the tool, the BAM and the reference, and the extra arguments.
type Run<'a> = (&'a str, &'a str, (&'a str, &'a str), Vec<&'a str>);

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gatk-tools/tests/data/hc_ploidy.txt.gz"),
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
    let dir = std::env::temp_dir().join("gatk-cli-hc-ploidy");
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
    let bed = |name: &str| {
        let path = dir.join(format!("{name}.bed"));
        std::fs::write(&path, section(&text, "bed", name)).expect("a BED");
        path.to_string_lossy().to_string()
    };
    let names = [
        "haploid",
        "triploid",
        "overlapping",
        "zero",
        "partial",
        "off-contig",
        "past-end",
        "not-a-number",
        "negative",
    ];
    let beds: Vec<String> = names.iter().map(|name| bed(name)).collect();
    let path = |name: &str| beds[names.iter().position(|n| *n == name).expect("a BED")].as_str();
    let regions = "--ploidy-regions";
    let gvcf = |bed: &'static str, sample: &'static str, span: &'static str, mode: &'static str| {
        vec![
            regions,
            path(bed),
            "-ERC",
            mode,
            "--sample-name",
            sample,
            "-L",
            span,
        ]
    };
    let runs: Vec<Run<'_>> = vec![
        (
            "hc-haploid",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("haploid")],
        ),
        (
            "hc-triploid",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("triploid")],
        ),
        (
            "hc-overlapping",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("overlapping")],
        ),
        (
            "hc-zero",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("zero"), "-L", "chr2"],
        ),
        (
            "hc-partial",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("partial"), "-L", "chr1:900-1300"],
        ),
        (
            "hc-haploid-gvcf",
            "HaplotypeCaller",
            e2e,
            gvcf("haploid", "s1", "chr1:200-900", "GVCF"),
        ),
        (
            "hc-haploid-bp",
            "HaplotypeCaller",
            e2e,
            gvcf("haploid", "s1", "chr1:200-900", "BP_RESOLUTION"),
        ),
        (
            "hc-triploid-bp",
            "HaplotypeCaller",
            e2e,
            gvcf("triploid", "s2", "chr1:950-1050", "BP_RESOLUTION"),
        ),
        (
            "hc-triploid-sample-ploidy",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("triploid"), "--sample-ploidy", "1"],
        ),
        (
            "hc-off-contig",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("off-contig")],
        ),
        (
            "hc-past-end",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("past-end")],
        ),
        (
            "hc-not-a-number",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("not-a-number")],
        ),
        (
            "hc-negative",
            "HaplotypeCaller",
            e2e,
            vec![regions, path("negative")],
        ),
        (
            "ramped-triploid",
            "RampedHaplotypeCaller",
            e2e,
            vec![regions, path("triploid")],
        ),
        (
            "ramped-haploid-gvcf",
            "RampedHaplotypeCaller",
            e2e,
            gvcf("haploid", "s1", "chr1:200-900", "GVCF"),
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
