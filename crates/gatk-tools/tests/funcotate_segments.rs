//! Conformance for `FuncotateSegments` against GATK 4.6.2.0, compared as the refusals and the
//! outputs of every run over one folder of data sources.
//!
//! Golden from `tools/readfilter-conformance/FuncotateSegmentsDump.java`.
//!
//! Reading a GTF, reading a reference and the annotation itself are not measured or ported. What
//! is measured is the folder: the version gate, the three-level walk, the config keys, and what a
//! segment has to look like to be one.
//!
//! # What this suite is for
//!
//!  * **the manifest's version being MISREAD**, so the refusal quotes `1.2.382015101`;
//!  * **the date range being unreachable**, because the throw it causes is swallowed;
//!  * **a missing manifest being accepted**, so the file that guards the range is optional;
//!  * **a reference version with no directory being refused**, not passed over;
//!  * **the config keys depending on the `type`**, and a missing one being named back;
//!  * **a GENCODE source being required, under that name in that case**;
//!  * **a segment having to EXCEED 150 bases**;
//!  * **and the gene list holding only the genes a segment covers.**

use gatk_corpus as corpus;
use gatk_tools::funcotate_segments::{
    acceptable, allele_of, call_of, check_config, column_prefix, gene_rows, is_segment,
    parse_manifest_version, resolve, validate_version, version_refusal, Call, ConfigError, GeneRow,
    ManifestVersion, ResolveError, SegmentGenes, SourceConfig, SourceFolder, SourceType,
    VersionVerdict, MIN_BASES_FOR_VALID_SEGMENT, UNKNOWN,
};
use std::collections::BTreeMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/funcotate_segments.txt.gz"),
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

/// A run the golden refused, as its exception class and its message.
fn refusal(text: &str, label: &str) -> (String, String) {
    let row = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = row.split_once(':').expect("a class and a message");
    (class.to_string(), unescape(message))
}

fn properties(pairs: &[(&str, &str)]) -> BTreeMap<String, String> {
    pairs
        .iter()
        .map(|(key, value)| ((*key).to_string(), (*value).to_string()))
        .collect()
}

/// The locatable XSV source the fixture writes beside GENCODE.
fn regions(reference: &str, folder: &str) -> SourceConfig {
    SourceConfig {
        path: format!("file://<dir>/{folder}/regions/{reference}/regions.config"),
        reference: reference.to_string(),
        properties: properties(&[
            ("name", "regions"),
            ("version", "1"),
            ("src_file", "regions.tsv"),
            ("origin_location", "test"),
            ("preprocessing_script", ""),
            ("type", "locatableXSV"),
            ("contig_column", "0"),
            ("start_column", "1"),
            ("end_column", "2"),
            ("xsv_delimiter", "\t"),
        ]),
    }
}

/// The GENCODE source the fixture writes, which is the one the tool insists on.
fn gencode(reference: &str, folder: &str) -> SourceConfig {
    SourceConfig {
        path: format!("file://<dir>/{folder}/gencode/{reference}/gencode.config"),
        reference: reference.to_string(),
        properties: properties(&[
            ("name", "Gencode"),
            ("version", "1"),
            ("src_file", "gencode.gtf"),
            ("origin_location", "test"),
            ("preprocessing_script", ""),
            ("type", "gencode"),
            ("gencode_fasta_path", "gencode.transcripts.fasta"),
            ("ncbi_build_version", reference),
        ]),
    }
}

fn folder(manifest: Option<&str>, reference: &str, name: &str) -> SourceFolder {
    SourceFolder {
        manifest: manifest.and_then(parse_manifest_version),
        sources: vec![regions(reference, name), gencode(reference, name)],
    }
}

/// The four folders the golden ran against, by the label of the run that used them.
fn folders() -> Vec<(&'static str, SourceFolder, &'static str)> {
    vec![
        (
            "annotated",
            folder(Some("Version: 1.7.hg38.20220101"), "hg38", "ds-good"),
            "hg38",
        ),
        (
            "old-version",
            folder(Some("Version: 1.2.hg38.20150101"), "hg38", "ds-old"),
            "hg38",
        ),
        ("no-manifest", folder(None, "hg38", "ds-none"), "hg38"),
        (
            "wrong-ref-version",
            folder(Some("Version: 1.7.hg38.20220101"), "hg19", "ds-hg19"),
            "hg38",
        ),
        (
            "hg19",
            folder(Some("Version: 1.7.hg38.20220101"), "hg19", "ds-hg19"),
            "hg19",
        ),
    ]
}

/// Every run of the golden either produced an output or a refusal, and the port agrees on which.
#[test]
fn every_run_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, source_folder, reference) in folders() {
        let produced = resolve(&source_folder, reference);
        let refused = text
            .lines()
            .any(|line| line.starts_with(&format!("error\t{label}\t")));
        assert_eq!(produced.is_err(), refused, "{label}");
        if let Err(error) = &produced {
            let (_, message) = refusal(&text, label);
            assert_eq!(error.message(), message, "{label}");
        } else {
            // A run that resolved wrote both files, and the second one holds the genes.
            assert!(!section(&text, "out", label).is_empty());
            assert!(!section(&text, "genes", label).is_empty());
        }
        compared += 1;
    }
    assert_eq!(compared, 5, "the folders the port reproduces");
}

/// Seven groups read as six: the `hg` number becomes the year and the day becomes the decorator.
#[test]
fn the_manifest_version_is_misread() {
    let text = golden();
    let version = parse_manifest_version("Version: 1.2.hg38.20150101").expect("a version");
    assert_eq!(
        version,
        ManifestVersion {
            major: 1,
            minor: 2,
            year: 38,
            month: 2015,
            day: 1,
            decorator: "01".to_string(),
        }
    );
    // The refusal quotes the scrambled string, which is how the misreading is visible at all.
    assert_eq!(version.display(), "1.2.382015101");
    let (class, message) = refusal(&text, "old-version");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException"
    );
    assert_eq!(version_refusal(Some(&version)), message);
    assert!(message.contains("1.2.382015101"), "{message}");
    // A line that does not match the pattern yields nothing rather than a refusal.
    assert_eq!(parse_manifest_version("Version: 1.7.20220101"), None);
    assert_eq!(parse_manifest_version("Source: test"), None);
}

/// The written year lands in the month, so the date can never be built and its range never runs.
#[test]
fn the_date_range_is_unreachable() {
    let good = parse_manifest_version("Version: 1.7.hg38.20220101").expect("a version");
    assert_eq!(
        validate_version(&good),
        VersionVerdict::DateUnrepresentable,
        "month 2022 is no month"
    );
    // The throw is swallowed and the flag keeps its initial value, so the folder is accepted.
    assert!(acceptable(Some(&good)));
    // Every plausible manifest lands the same way, whatever its date.
    for line in [
        "Version: 1.6.hg19.20190124",
        "Version: 1.7.hg38.20220101",
        "Version: 1.8.hg38.20230908",
        "Version: 1.8.hg38.19000101",
    ] {
        let version = parse_manifest_version(line).expect("a version");
        assert_eq!(
            validate_version(&version),
            VersionVerdict::DateUnrepresentable,
            "{line}"
        );
        assert!(acceptable(Some(&version)), "{line}");
    }
    // Only the major and minor numbers can turn a folder away.
    let old = parse_manifest_version("Version: 1.2.hg38.20150101").expect("a version");
    assert_eq!(validate_version(&old), VersionVerdict::Refused);
    let new = parse_manifest_version("Version: 1.9.hg38.20220101").expect("a version");
    assert_eq!(validate_version(&new), VersionVerdict::Refused);
    let other_major = parse_manifest_version("Version: 2.7.hg38.20220101").expect("a version");
    assert_eq!(validate_version(&other_major), VersionVerdict::Refused);
}

/// The one file that guards the range is the one file the folder can do without.
#[test]
fn a_missing_manifest_is_accepted() {
    let text = golden();
    assert!(acceptable(None));
    // The run with no manifest produced exactly what the run with a good one produced.
    assert_eq!(
        section(&text, "out", "no-manifest"),
        section(&text, "out", "annotated")
    );
    assert_eq!(
        section(&text, "genes", "no-manifest"),
        section(&text, "genes", "annotated")
    );
    let no_manifest = folder(None, "hg38", "ds-none");
    assert!(resolve(&no_manifest, "hg38").is_ok());
}

/// The whole folder coming up empty is refused; the individual source is passed over in silence.
#[test]
fn a_reference_version_with_no_directory_is_refused() {
    let text = golden();
    let hg19_only = folder(Some("Version: 1.7.hg38.20220101"), "hg19", "ds-hg19");
    let produced = resolve(&hg19_only, "hg38").expect_err("no such reference");
    assert_eq!(
        produced,
        ResolveError::NoSources {
            reference: "hg38".to_string()
        }
    );
    let (class, message) = refusal(&text, "wrong-ref-version");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException"
    );
    assert_eq!(produced.message(), message);
    // The same folder answers normally for the version it does have.
    assert_eq!(resolve(&hg19_only, "hg19").expect("sources").len(), 2);
    // One source missing that version is not enough on its own: the GENCODE one is still there.
    let mut mixed = hg19_only.clone();
    mixed.sources[1] = gencode("hg38", "ds-hg19");
    assert_eq!(resolve(&mixed, "hg38").expect("sources").len(), 1);
}

/// The universal keys first, then the type's own, and the missing one is named back.
#[test]
fn the_config_keys_depend_on_the_type() {
    let text = golden();
    let path = "file://<dir>/ds-broken/regions/hg38/regions.config";
    let mut broken = regions("hg38", "ds-broken").properties;
    broken.remove("end_column");
    let produced = check_config(path, &broken).expect_err("no end column");
    assert_eq!(
        produced,
        ConfigError::MissingKey {
            path: path.to_string(),
            key: "end_column".to_string()
        }
    );
    let (class, message) = refusal(&text, "missing-config-key");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    );
    assert_eq!(produced.message(), message);
    // A VCF source needs none of the type-specific keys, so the same removal leaves it valid.
    let mut as_vcf = broken.clone();
    as_vcf.insert("type".to_string(), "vcf".to_string());
    assert_eq!(check_config(path, &as_vcf), Ok(SourceType::Vcf));
    // The universal keys are checked before the type, so a config missing both is refused for the
    // universal one.
    let mut no_type = broken.clone();
    no_type.remove("type");
    assert_eq!(
        check_config(path, &no_type),
        Err(ConfigError::MissingKey {
            path: path.to_string(),
            key: "type".to_string()
        })
    );
    // And a type that is not one of the five is refused by value.
    let mut odd = broken.clone();
    odd.insert("type".to_string(), "spreadsheet".to_string());
    assert!(matches!(
        check_config(path, &odd),
        Err(ConfigError::UnknownType { .. })
    ));
    // The name is matched without regard to case, which the directory names are not.
    let mut upper = regions("hg38", "ds-good").properties;
    upper.insert("type".to_string(), "LOCATABLEXSV".to_string());
    assert_eq!(check_config(path, &upper), Ok(SourceType::LocatableXsv));
}

/// Whatever else the folder holds, and its columns are prefixed with its `name` key.
#[test]
fn a_gencode_source_is_required() {
    let text = golden();
    let without = SourceFolder {
        manifest: parse_manifest_version("Version: 1.7.hg38.20220101"),
        sources: vec![regions("hg38", "ds-good")],
    };
    let produced = resolve(&without, "hg38").expect_err("no gencode");
    assert_eq!(produced, ResolveError::NoGencode);
    assert_eq!(
        produced.message(),
        "ERROR: a Gencode datasource is required!"
    );
    // The prefix is the source's own name, capital G included, which is what the output carries.
    let source = gencode("hg38", "ds-good");
    assert_eq!(column_prefix(&source.properties), "Gencode_1");
    assert!(section(&text, "out", "annotated").contains("Gencode_1_genes"));
    // A source named in lower case produces a prefix the renderers do not look for.
    let mut lower = source.properties.clone();
    lower.insert("name".to_string(), "gencode".to_string());
    assert_eq!(column_prefix(&lower), "gencode_1");
    // The empty check runs before the gencode one, so no sources at all is the other message.
    assert_eq!(
        resolve(&without, "hg38").expect_err("no gencode"),
        ResolveError::NoGencode
    );
    assert!(matches!(
        resolve(&without, "hg19").expect_err("no sources"),
        ResolveError::NoSources { .. }
    ));
}

/// Exceed, not reach: a segment of exactly the minimum is not one.
#[test]
fn a_segment_must_exceed_the_minimum_length() {
    let text = golden();
    assert_eq!(MIN_BASES_FOR_VALID_SEGMENT, 150);
    let insertion = vec!["<INS>".to_string()];
    // The refused run's segment: 1200 to 1300 is 101 bases.
    assert!(!is_segment(
        &insertion,
        1200,
        1300,
        MIN_BASES_FOR_VALID_SEGMENT
    ));
    let (class, message) = refusal(&text, "short-segment");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    );
    assert!(message.contains("chr1:1200-1300"), "{message}");
    // The three segments the good run carried are all long enough.
    for (start, end) in [(1200, 1800), (5000, 5500), (9400, 9900)] {
        assert!(is_segment(
            &insertion,
            start,
            end,
            MIN_BASES_FOR_VALID_SEGMENT
        ));
    }
    // The boundary is strict: 150 bases is not enough, 151 is.
    assert!(!is_segment(&insertion, 1, 150, MIN_BASES_FOR_VALID_SEGMENT));
    assert!(is_segment(&insertion, 1, 151, MIN_BASES_FOR_VALID_SEGMENT));
    assert!(!is_segment(
        &["<SNP>".to_string()],
        1200,
        1800,
        MIN_BASES_FOR_VALID_SEGMENT
    ));
}

/// The first name present decides, and an amplification becomes an insertion.
#[test]
fn a_call_becomes_a_symbolic_allele() {
    let text = golden();
    // The fixture's segment file writes its calls under `CALL`.
    let segments = section(&text, "tsv", "segments");
    let mut lines = segments.lines();
    assert_eq!(lines.next().expect("a header"), "CONTIG\tSTART\tEND\tCALL");
    let calls: Vec<Option<Call>> = lines
        .filter(|line| !line.is_empty())
        .map(|line| {
            let column = line.split('\t').nth(3).expect("a call");
            call_of(&properties(&[("CALL", column)]))
        })
        .collect();
    assert_eq!(
        calls,
        vec![
            Some(Call::Amplification),
            Some(Call::Neutral),
            Some(Call::Deletion)
        ]
    );
    let alleles: Vec<&str> = calls.iter().map(|call| allele_of(*call)).collect();
    assert_eq!(alleles, vec!["<INS>", "<COPY_NEUTRAL>", "<DEL>"]);
    // The refused short segment is an amplification, and the message shows its `<INS>`.
    let (_, message) = refusal(&text, "short-segment");
    assert!(message.contains("<INS>"), "{message}");
    // A file with no call column at all yields the unspecified allele.
    assert_eq!(call_of(&properties(&[])), None);
    assert_eq!(allele_of(None), "<*>");
    // The first NAME present decides even when its value is not a call, so a later column that
    // does carry one does not rescue the segment.
    assert_eq!(
        call_of(&properties(&[("CALL", "NA"), ("Segment_Call", "+")])),
        None
    );
    assert_eq!(
        call_of(&properties(&[("Segment_Call", "+")])),
        Some(Call::Amplification)
    );
    assert_eq!(call_of(&properties(&[("Call", "-")])), Some(Call::Deletion));
}

/// Only the genes a segment covers, and a gene can hold two rows.
#[test]
fn the_gene_list_holds_only_the_covered_genes() {
    let text = golden();
    let list = section(&text, "genes", "annotated");
    let mut lines = list.lines();
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    assert_eq!(header[0], "gene");
    assert_eq!(header[1], "exon");
    let measured: Vec<GeneRow> = lines
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            GeneRow {
                gene: columns[0].to_string(),
                exon: columns[1].to_string(),
            }
        })
        .collect();
    // The three segments: the first covers ALPHA and starts inside its one exon, the other two
    // cover nothing at all.
    let produced = gene_rows(&[
        SegmentGenes {
            genes: vec!["ALPHA".to_string()],
            start: Some(("ALPHA".to_string(), "0+".to_string())),
            end: None,
        },
        SegmentGenes {
            genes: Vec::new(),
            start: None,
            end: None,
        },
        SegmentGenes {
            genes: Vec::new(),
            start: None,
            end: None,
        },
    ]);
    assert_eq!(produced, measured);
    assert_eq!(produced.len(), 2, "one gene, two rows");
    // The output beside it carried a row per segment, so the two files disagree on length.
    let rows = section(&text, "out", "annotated").lines().count() - 1;
    assert_eq!(rows, 3);
    // A run over segments that cover nothing yields an empty list rather than no file.
    assert!(gene_rows(&[SegmentGenes {
        genes: Vec::new(),
        start: None,
        end: None
    }])
    .is_empty());
}

/// A column the segment file did not carry, which is not the same absence as an empty one.
#[test]
fn an_absent_segment_column_is_unknown() {
    let text = golden();
    let out = section(&text, "out", "annotated");
    let mut lines = out.lines();
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let first: Vec<&str> = lines.next().expect("a row").split('\t').collect();
    let at = |name: &str| first[header.iter().position(|c| *c == name).expect(name)];
    // The fixture's segment file has four columns, so the rest of the SEG columns are absent.
    assert_eq!(at("ref_allele"), UNKNOWN);
    assert_eq!(at("alt_allele"), UNKNOWN);
    assert_eq!(at("start_gene"), UNKNOWN);
    // The columns it did carry are its own values.
    assert_eq!(at("chr"), "chr1");
    assert_eq!(at("start"), "1200");
    assert_eq!(at("end"), "1800");
    assert_eq!(at("Segment_Call"), "+");
    // A column a SOURCE did not fill is empty rather than unknown, which is the other absence.
    assert_eq!(at("Gencode_1_ref_allele"), "");
    assert_eq!(at("Gencode_1_genes"), "ALPHA");
    // And the locatable XSV source beside GENCODE contributes no column at all.
    assert!(!header.iter().any(|column| column.starts_with("regions_")));
}
