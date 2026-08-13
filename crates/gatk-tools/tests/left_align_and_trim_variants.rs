//! Conformance for `LeftAlignAndTrimVariants` against GATK 4.6.2.0, compared as the records every
//! run writes.
//!
//! Golden from `tools/readfilter-conformance/LeftAlignAndTrimVariantsDump.java`.
//!
//! # What this suite is for
//!
//!  * **the window is bounded by the previous record as written**, so aligning one frees the next;
//!  * **a skipped indel still bounds the record after it**;
//!  * **the pieces of a split record bound each other**;
//!  * **and a contig boundary lifts the bound entirely**.

use gatk_corpus as corpus;
use gatk_engine::variant_context_utils::{Allele, Variant};
use gatk_tools::left_align_and_trim_variants::{
    left_align_and_trim_variants, Arguments, DEFAULT_MAX_INDEL_SIZE, DEFAULT_MAX_LEADING_BASES,
};
use std::collections::HashMap;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/left_align_and_trim_variants.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

/// The reverse of the dump's `escape`, scanning once so a real backslash is never read as a tab.
fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut characters = text.chars();
    while let Some(character) = characters.next() {
        if character != '\\' {
            out.push(character);
            continue;
        }
        match characters.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// The reference the tool aligned against, contig by contig.
fn reference(text: &str) -> HashMap<String, Vec<u8>> {
    rows(text, "reference")
        .into_iter()
        .map(|row| (row[0].to_string(), row[1].as_bytes().to_vec()))
        .collect()
}

/// The arguments each run was made with, which the golden does not carry.
fn arguments(label: &str) -> Arguments {
    let base = Arguments {
        dont_trim_alleles: false,
        split_multiallelics: false,
        max_indel_size: DEFAULT_MAX_INDEL_SIZE,
        max_leading_bases: DEFAULT_MAX_LEADING_BASES,
    };
    match label {
        "long-indel-allowed" => Arguments {
            max_indel_size: 500,
            ..base
        },
        "long-indel-skipped" => Arguments {
            max_indel_size: 5,
            ..base
        },
        "multiallelic-split" => Arguments {
            split_multiallelics: true,
            ..base
        },
        "untrimmed-no-trim" => Arguments {
            dont_trim_alleles: true,
            ..base
        },
        "narrow-window" => Arguments {
            max_leading_bases: 2,
            ..base
        },
        "zero-window" => Arguments {
            max_leading_bases: 0,
            ..base
        },
        _ => base,
    }
}

/// Which input file each run was given.
fn input_of(label: &str) -> &str {
    match label {
        "long-indel-allowed" | "long-indel-skipped" => "long-indel",
        "multiallelic-split" => "multiallelic",
        "untrimmed-no-trim" => "untrimmed",
        "narrow-window" | "zero-window" => "apart",
        other => other,
    }
}

/// One `CHROM POS ID REF ALT ...` line as a record.
fn parse_record(line: &str) -> Variant {
    let fields: Vec<&str> = line.split('\t').collect();
    let mut alleles = vec![Allele::new(fields[3].as_bytes(), true)];
    for alternate in fields[4].split(',') {
        alleles.push(Allele::new(alternate.as_bytes(), false));
    }
    let start: i32 = fields[1].parse().expect("a position");
    Variant {
        contig: fields[0].to_string(),
        start,
        stop: start + fields[3].len() as i32 - 1,
        alleles,
        genotypes: Vec::new(),
        attributes: Vec::new(),
    }
}

/// The records of one input vcf, in file order.
fn input(text: &str, label: &str) -> Vec<Variant> {
    let whole = rows(text, "input")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no input {label}"))[1]
        .to_string();
    unescape(&whole)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(parse_record)
        .collect()
}

/// The records one run wrote, as the dump printed them.
fn written(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .filter(|line| !line.starts_with('#'))
        .collect()
}

/// A record as a vcf line, so the two can be compared as text.
fn rendered(variant: &Variant) -> String {
    let alternates: Vec<String> = variant.alleles[1..]
        .iter()
        .map(|allele| String::from_utf8_lossy(&allele.bases).to_string())
        .collect();
    format!(
        "{}\t{}\t.\t{}\t{}\t.\t.\t.",
        variant.contig,
        variant.start,
        String::from_utf8_lossy(&variant.alleles[0].bases),
        alternates.join(",")
    )
}

fn runs(text: &str) -> Vec<String> {
    let mut seen: Vec<String> = Vec::new();
    for row in rows(text, "vcfline") {
        if !seen.iter().any(|name| name == row[0]) {
            seen.push(row[0].to_string());
        }
    }
    seen
}

#[test]
fn every_run_writes_the_records_the_reference_wrote() {
    let text = golden();
    let bases = reference(&text);
    let lookup = |contig: &str| bases.get(contig).cloned();
    let all = runs(&text);
    assert!(all.len() >= 13, "every run is in the golden: {}", all.len());

    for run in &all {
        let variants = input(&text, input_of(run));
        let ours = left_align_and_trim_variants(&variants, &lookup, arguments(run))
            .unwrap_or_else(|error| panic!("{run}: {}", error.message()));
        let ours: Vec<String> = ours.iter().map(rendered).collect();
        assert_eq!(ours, written(&text, run), "run/{run}");
    }
}

/// Two indels a base apart in the input, both of which move.
#[test]
fn aligning_one_record_relaxes_the_bound_on_the_next() {
    let text = golden();
    let adjacent = written(&text, "adjacent");
    assert_eq!(adjacent[0], "chr1\t10\t.\tGA\tG\t.\t.\t.");
    assert_eq!(adjacent[1], "chr1\t18\t.\tAG\tA\t.\t.\t.");
    // One base further apart in the input, and the same two records out.
    assert_eq!(written(&text, "nearly-adjacent"), adjacent);
}

/// The same file, two lengths, and the record that was skipped stays where it was.
#[test]
fn a_skipped_indel_is_written_untouched() {
    let text = golden();
    let aligned = written(&text, "long-indel");
    let skipped = written(&text, "long-indel-skipped");
    assert!(aligned[0].starts_with("chr1\t30\t"));
    assert!(skipped[0].starts_with("chr1\t34\t"));
    // And the record after it lands in the same place either way.
    assert_eq!(aligned[1], skipped[1]);
}

/// One record in, two out, at two positions, because the first piece bounds the second.
#[test]
fn the_pieces_of_a_split_record_bound_each_other() {
    let text = golden();
    let split = written(&text, "multiallelic-split");
    assert_eq!(split.len(), 3);
    assert_eq!(split[0], "chr1\t10\t.\tGA\tG\t.\t.\t.");
    assert_eq!(split[1], "chr1\t12\t.\tA\tAA\t.\t.\t.");
    // Unsplit, the same input gives one record carrying both alternates.
    assert_eq!(
        written(&text, "multiallelic")[0],
        "chr1\t10\t.\tGA\tG,GAA\t.\t.\t."
    );
}
