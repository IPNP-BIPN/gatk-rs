//! Conformance for `CheckReferenceCompatibility` against GATK 4.6.2.0, compared as the whole
//! output table of every run.
//!
//! Golden from `tools/readfilter-conformance/CheckReferenceCompatibilityDump.java`.
//!
//! # What this suite is for
//!
//!  * **the verdict and its summary**, which differ between the two paths for the same relation;
//!  * **a VCF never reaching the MD5 path**, however many `M5` fields its header carries;
//!  * **`COMPATIBLE_SUBSET` needing the flag set to be exactly `{SUBSET}`**;
//!  * **the missing sequences**, listed from the reference's side as a Java list;
//!  * **and the two refusals**, a BAM with a VCF and neither.

use gatk_corpus as corpus;
use gatk_tools::check_reference_compatibility::{
    check_input, evaluate_with_md5, evaluate_without_md5, md5s_present, write_table, Compatibility,
    DictionaryCompatibility, InputError,
};
use gatk_tools::compare_references::{
    build, compare_all, Md5Mode, Pair, Reference, Sequence, Status,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/check_reference_compatibility.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn value(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}=");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(prefix.as_str()))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

/// A reference read out of the golden's dictionary, whose `M5` is the truth.
fn reference(text: &str, label: &str) -> Reference {
    let sequences = value(text, "dict", label)
        .lines()
        .filter(|line| line.starts_with("@SQ\t"))
        .map(|line| {
            let field = |tag: &str| {
                line.split('\t')
                    .find_map(|part| part.strip_prefix(&format!("{tag}:")))
                    .map(str::to_string)
            };
            let md5 = field("M5");
            Sequence {
                name: field("SN").expect("a name"),
                length: field("LN").and_then(|v| v.parse().ok()).expect("a length"),
                calculated_md5: md5.clone().unwrap_or_default(),
                md5,
            }
        })
        .collect();
    Reference {
        column: format!("{label}.fasta"),
        sequences,
    }
}

/// The BAM's dictionary, which is the base reference's with or without the digests.
fn bam_input(text: &str, with_md5: bool) -> Reference {
    let mut input = reference(text, "base");
    input.column = if with_md5 {
        "with-md5.bam".to_string()
    } else {
        "without-md5.bam".to_string()
    };
    if !with_md5 {
        for sequence in &mut input.sequences {
            sequence.md5 = None;
        }
    }
    input
}

/// The pair the MD5 path analyses: the input against one reference, in that order.
fn pair_against(input: &Reference, reference: &Reference) -> Pair {
    let references = vec![input.clone(), reference.clone()];
    let table = build(&references, Md5Mode::UseDict).expect("a table");
    compare_all(&table, &references)
        .expect("an analysis")
        .into_iter()
        .next()
        .expect("one pair")
}

/// `getMissingSequencesIfSubset`: the reference's names the input does not have.
fn missing(input: &Reference, reference: &Reference) -> Vec<String> {
    reference
        .sequences
        .iter()
        .map(|sequence| sequence.name.clone())
        .filter(|name| !input.sequences.iter().any(|other| other.name == *name))
        .collect()
}

#[test]
fn the_md5_path_matches_the_golden() {
    let text = golden();
    let input = bam_input(&text, true);
    assert!(md5s_present(
        &input
            .sequences
            .iter()
            .map(|sequence| sequence.md5.clone())
            .collect::<Vec<Option<String>>>()
    ));
    let mut compared = 0;
    for (label, name) in [
        ("bam-md5-exact", "base"),
        ("bam-md5-altered", "altered"),
        ("bam-md5-subset", "extra"),
    ] {
        let against = reference(&text, name);
        let record = evaluate_with_md5(&pair_against(&input, &against), &missing(&input, &against));
        assert_eq!(
            write_table(&input.column, &[record]),
            value(&text, "table", label),
            "{label}"
        );
        compared += 1;
    }

    // Every reference at once is one row each, in the order they were given.
    let records: Vec<_> = ["base", "altered", "extra"]
        .iter()
        .map(|name| {
            let against = reference(&text, name);
            evaluate_with_md5(&pair_against(&input, &against), &missing(&input, &against))
        })
        .collect();
    assert_eq!(
        write_table(&input.column, &records),
        value(&text, "table", "bam-md5-all")
    );
    compared += 1;
    assert_eq!(compared, 4, "the golden's MD5 runs");
}

/// The other path, which every VCF run takes whatever its header says.
#[test]
fn the_name_and_length_path_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, input_name, status, missing_names) in [
        (
            "no-md5-exact",
            "without-md5.vcf",
            DictionaryCompatibility::Identical,
            vec![],
        ),
        (
            // The altered reference has the same names and lengths, so this path calls it
            // compatible where the MD5 path calls it not.
            "no-md5-altered",
            "without-md5.vcf",
            DictionaryCompatibility::Identical,
            vec![],
        ),
        (
            "no-md5-subset",
            "without-md5.vcf",
            DictionaryCompatibility::Superset,
            vec!["chr3".to_string()],
        ),
        // A VCF whose header carries M5 for every contig still lands here.
        (
            "md5-exact",
            "with-md5.vcf",
            DictionaryCompatibility::Identical,
            vec![],
        ),
        // Including the one whose M5 is a lie.
        (
            "md5-lying",
            "lying.vcf",
            DictionaryCompatibility::Identical,
            vec![],
        ),
    ] {
        let reference_name = match label {
            "no-md5-altered" => "altered.fasta",
            "no-md5-subset" => "extra.fasta",
            _ => "base.fasta",
        };
        let record = evaluate_without_md5(reference_name, input_name, status, &missing_names);
        assert_eq!(
            write_table(input_name, &[record]),
            value(&text, "table", label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "the golden's name-and-length runs");
}

/// The two paths disagree about the altered reference, which is the point of having both.
#[test]
fn the_two_paths_disagree_about_the_same_pair() {
    let text = golden();
    assert!(value(&text, "table", "bam-md5-altered").contains("NOT_COMPATIBLE"));
    assert!(value(&text, "table", "no-md5-altered").contains("COMPATIBLE\t"));
}

/// `COMPATIBLE_SUBSET` needs the flag set to be exactly `{SUBSET}`.
#[test]
fn a_subset_with_another_flag_is_not_compatible() {
    let mut pair = Pair {
        first: "input.bam".to_string(),
        second: "reference.fasta".to_string(),
        analysis: Default::default(),
    };
    pair.analysis.insert(Status::Subset);
    assert_eq!(
        evaluate_with_md5(&pair, &["chr3".to_string()]).compatibility,
        Compatibility::CompatibleSubset
    );
    pair.analysis.insert(Status::DifferInSequenceNames);
    let record = evaluate_with_md5(&pair, &["chr3".to_string()]);
    assert_eq!(record.compatibility, Compatibility::NotCompatible);
    assert!(record
        .summary
        .starts_with("Status: [DIFFER_IN_SEQUENCE_NAMES, SUBSET]."));
}

#[test]
fn the_two_refusals_match_the_golden() {
    let text = golden();
    let error = check_input(true, 1, true).expect_err("both inputs");
    assert_eq!(error, InputError::BothBamAndVcf);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "both-inputs")
    );

    let error = check_input(false, 0, false).expect_err("no input");
    assert_eq!(error, InputError::NoInput);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-input")
    );

    assert_eq!(check_input(true, 1, false), Ok(()));
}
