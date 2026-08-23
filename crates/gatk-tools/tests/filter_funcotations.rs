//! Conformance for `FilterFuncotations` against GATK 4.6.2.0, compared as the CLINSIG value and
//! the FILTER value of every variant.
//!
//! Golden from `tools/readfilter-conformance/FilterFuncotationsDump.java`.
//!
//! # What this suite is for
//!
//!  * **every rule of the four filters**, and the fact that a filter needs all of its rules;
//!  * **the three kinds of missing ExAC data that all read as a frequency of zero**;
//!  * **the gnomAD path's present-dataset rule**, and its refusal to catch anything;
//!  * **the joined CLINSIG value**, whose order is a `HashSet`'s and whose two autosomal recessive
//!    filters share one name;
//!  * **and the second pass**, where a gene needs more than one het variant.
//!
//! The `no-funcotation-header` run is not replayed: it fails in the VCF header, which is htsjdk's
//! parsing rather than the tool's rules.

use gatk_corpus as corpus;
use gatk_tools::filter_funcotations::{
    clinsig, compound_het_variants, gnomad_max_maf, matching_filters, AlleleFrequencySource,
    Funcotations, NumberFormatError, Reference, Variant,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/filter_funcotations.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .expect("the golden carries the refusal")
        .to_string()
}

/// The funcotation keys the header declares, in the order the values are written.
fn keys(vcf: &str) -> Vec<String> {
    let line = vcf
        .lines()
        .find(|line| line.starts_with("##INFO=<ID=FUNCOTATION,"))
        .expect("a FUNCOTATION header line");
    let description = line
        .split_once("Funcotation fields are: ")
        .expect("the preamble")
        .1
        .trim_end_matches("\">");
    description.split('|').map(str::to_string).collect()
}

/// Every variant of a VCF, with its transcripts' funcotations pruned of empty values.
fn variants(vcf: &str) -> Vec<(Variant, Vec<Funcotations>)> {
    let names = keys(vcf);
    vcf.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let genotype = fields[9];
            let variant = Variant {
                contig: fields[0].to_string(),
                start: fields[1].parse().expect("a position"),
                end: fields[1].parse().expect("a position"),
                reference_allele: fields[3].to_string(),
                alternate_alleles: fields[4].split(',').map(str::to_string).collect(),
                het_count: i32::from(genotype == "0/1"),
                hom_var_count: i32::from(genotype == "1/1"),
            };
            let attribute = fields[7]
                .split(';')
                .find_map(|entry| entry.strip_prefix("FUNCOTATION="))
                .expect("a FUNCOTATION attribute");
            let transcripts = attribute
                .split('#')
                .map(|transcript| {
                    let values: Vec<&str> = transcript
                        .trim_start_matches('[')
                        .trim_end_matches(']')
                        .split('|')
                        .collect();
                    let entries: Vec<(&str, &str)> = names
                        .iter()
                        .enumerate()
                        .map(|(index, name)| {
                            (name.as_str(), values.get(index).copied().unwrap_or(""))
                        })
                        .collect();
                    Funcotations::new(&entries)
                })
                .collect();
            (variant, transcripts)
        })
        .collect()
}

/// What the tool wrote for each variant: its position, its FILTER and its CLINSIG.
fn written(vcf: &str) -> Vec<(i32, String, String)> {
    vcf.lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let clinsig = fields[7]
                .split(';')
                .find_map(|entry| entry.strip_prefix("CLINSIG="))
                .expect("a CLINSIG attribute")
                .to_string();
            (
                fields[1].parse().expect("a position"),
                fields[6].to_string(),
                clinsig,
            )
        })
        .collect()
}

fn check(text: &str, label: &str, reference: Reference, source: AlleleFrequencySource) {
    let input = variants(&section(text, "input", label));
    let compound = compound_het_variants(&input, reference);
    let ours: Vec<(i32, String, String)> = input
        .iter()
        .map(|(variant, transcripts)| {
            let matching = matching_filters(transcripts, variant, reference, source, &compound)
                .expect("no unparseable frequency in this run");
            let (value, passes) = clinsig(&matching);
            (
                variant.start,
                if passes { "PASS" } else { "NOT_CLINSIG" }.to_string(),
                value,
            )
        })
        .collect();
    assert_eq!(ours, written(&section(text, "output", label)), "{label}");
}

#[test]
fn every_exac_rule() {
    let text = golden();
    check(&text, "exac", Reference::Hg19, AlleleFrequencySource::Exac);
}

#[test]
fn the_gnomad_path_has_a_rule_of_its_own() {
    let text = golden();
    check(
        &text,
        "gnomad",
        Reference::Hg19,
        AlleleFrequencySource::Gnomad,
    );
}

#[test]
fn an_unparseable_gnomad_frequency_is_not_caught() {
    let text = golden();
    let input = variants(&section(&text, "input", "gnomad-unparseable"));
    let error = gnomad_max_maf(&input[0].1[0]).expect_err("a refusal");
    assert_eq!(
        error,
        NumberFormatError {
            value: "many".to_string()
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "gnomad-unparseable")
    );
}

#[test]
fn the_keys_carry_the_gencode_version_of_the_reference() {
    let text = golden();
    check(
        &text,
        "hg38-keys",
        Reference::Hg38,
        AlleleFrequencySource::Exac,
    );
    assert_eq!(Reference::Hg38.gencode_version(), 27);
    assert_eq!(Reference::B37.gencode_version(), 19);
    // The same funcotations under hg19 would have matched, which is what the run shows.
    let input = variants(&section(&text, "input", "hg38-keys"));
    let compound = compound_het_variants(&input, Reference::Hg19);
    let matching = matching_filters(
        &input[0].1,
        &input[0].0,
        Reference::Hg19,
        AlleleFrequencySource::Exac,
        &compound,
    )
    .expect("no unparseable frequency");
    assert_eq!(matching, vec!["LOF".to_string()]);
}

#[test]
fn the_two_autosomal_recessive_filters_share_one_name() {
    let text = golden();
    let input = variants(&section(&text, "input", "exac"));
    let compound = compound_het_variants(&input, Reference::Hg19);
    // A hom-var and a compound het in the same gene both answer to AR, and the variant that is
    // both contributes one entry.
    let everything = input
        .iter()
        .find(|(variant, _)| variant.start == 2100)
        .expect("the variant that matches everything");
    let matching = matching_filters(
        &everything.1,
        &everything.0,
        Reference::Hg19,
        AlleleFrequencySource::Exac,
        &compound,
    )
    .expect("no unparseable frequency");
    assert_eq!(matching.len(), 4);
    assert_eq!(clinsig(&matching).0, "AR,CLINVAR,LOF,LMM");
}

#[test]
fn a_lone_het_in_an_autosomal_recessive_gene_is_not_compound() {
    let text = golden();
    let input = variants(&section(&text, "input", "exac"));
    let compound = compound_het_variants(&input, Reference::Hg19);
    // Two hets in MUTYH, and one in ATP7B whose only company is a hom-var.
    assert_eq!(compound.len(), 2);
    assert!(compound
        .iter()
        .all(|key| key.start == 1900 || key.start == 2000));
}
