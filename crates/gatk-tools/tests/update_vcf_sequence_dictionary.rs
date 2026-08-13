//! Conformance for `UpdateVCFSequenceDictionary` against GATK 4.6.2.0, compared as the records
//! each run put on disk and as the refusals, message and Java class both.
//!
//! Golden from `tools/readfilter-conformance/UpdateVCFSequenceDictionaryDump.java`.
//!
//! # What this suite is for
//!
//!  * **five refusals in five Java classes**, for one tool's arguments and records;
//!  * **the validation is per record and the writer is already open**, so the golden's partial
//!    files are what a refusal leaves behind;
//!  * **the end is `vc.getEnd()`**, which an INFO END overrides;
//!  * **and the contig lines are replaced, not merged**, so a contig the input had and the
//!    dictionary lacks is gone from the header.

use gatk_corpus as corpus;
use gatk_tools::update_vcf_sequence_dictionary::{
    best_available_dictionary, check_replace, check_variant, update_dictionary,
    UpdateDictionaryError,
};
use htsjdk_bam::header::SequenceRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/update_vcf_sequence_dictionary.txt.gz"),
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

fn labelled(text: &str, kind: &str, label: &str) -> String {
    rows(text, kind)
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no {kind} row for {label}"))
        .get(1)
        .copied()
        .unwrap_or("")
        .to_string()
}

/// The dictionary of one `.dict` the golden holds, which is a SAM header.
fn dictionary(text: &str, label: &str) -> Vec<SequenceRecord> {
    unescape(&labelled(text, "dict", label))
        .lines()
        .filter(|line| line.starts_with("@SQ"))
        .map(|line| {
            let mut name = String::new();
            let mut length = 0;
            for field in line.split('\t') {
                if let Some(value) = field.strip_prefix("SN:") {
                    name = value.to_string();
                }
                if let Some(value) = field.strip_prefix("LN:") {
                    length = value.parse().expect("a length");
                }
            }
            SequenceRecord::new(&name, length)
        })
        .collect()
}

/// One `CHROM POS ID REF ALT QUAL FILTER INFO` line, with END applied where it is present.
fn parse_record(line: &str) -> VariantContext {
    let fields: Vec<&str> = line.split('\t').collect();
    let mut alleles = vec![Allele::create(fields[3].as_bytes(), true).expect("a reference")];
    for alternate in fields[4].split(',') {
        alleles.push(Allele::create(alternate.as_bytes(), false).expect("an alternate"));
    }
    let mut variant =
        VariantContext::new(fields[0], fields[1].parse().expect("a position"), alleles);
    variant.id = fields[2].to_string();
    for entry in fields[7].split(';') {
        if let Some(value) = entry.strip_prefix("END=") {
            // `getEnd()` is the END attribute where there is one, with no check at all.
            variant.stop = value.parse().expect("an end");
        }
    }
    variant
}

/// The records of one input vcf, in file order.
fn input(text: &str, label: &str) -> Vec<String> {
    unescape(&labelled(text, "input", label))
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.to_string())
        .collect()
}

/// The records one run left on disk, from its output or its partial output.
fn written(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .filter(|line| !line.starts_with('#'))
        .collect()
}

/// The refusal the golden holds. It is unescaped, because the missing-lengths message carries a
/// real newline that the dump escaped along with everything else.
fn refusal(text: &str, label: &str) -> String {
    unescape(&labelled(text, "error", label))
}

/// The message as the golden holds it, which is the class, a colon, and whatever prefix the class
/// adds to `getMessage`.
fn rendered(error: &UpdateDictionaryError) -> String {
    let prefix = match error {
        UpdateDictionaryError::TwoDictionaries
        | UpdateDictionaryError::MissingContigLengths { .. } => String::new(),
        UpdateDictionaryError::NoDictionary => {
            "Argument source-dictionary was missing: ".to_string()
        }
        _ => "Illegal argument value: ".to_string(),
    };
    format!("{}:{}{}", error.java_class(), prefix, error.message())
}

#[test]
fn every_record_refusal_is_the_reference() {
    let text = golden();
    let good = dictionary(&text, "dictionary");
    let partial = dictionary(&text, "partial-dictionary");

    for (run, label, dictionary) in [
        ("unknown-contig", "unknown-contig", &good),
        ("past-end", "past-end", &good),
        ("end-attribute", "end-attribute", &good),
        ("partial", "bare", &partial),
    ] {
        let records = input(&text, label);
        let variants: Vec<VariantContext> = records.iter().map(|line| parse_record(line)).collect();
        let (kept, error) = update_dictionary(dictionary, &variants);
        let error = error.unwrap_or_else(|| panic!("run {run} is refused"));
        assert_eq!(rendered(&error), refusal(&text, run), "error/{run}");

        // What the reference left on disk before it threw.
        let ours: Vec<String> = kept
            .into_iter()
            .map(|index| records[index].clone())
            .collect();
        assert_eq!(
            ours,
            written(&text, &format!("{run}-partial")),
            "partial/{run}"
        );
    }
}

#[test]
fn every_argument_refusal_is_the_reference() {
    let text = golden();
    let good = dictionary(&text, "dictionary");
    let empty = dictionary(&text, "empty-dictionary");
    let no_length = dictionary(&text, "no-length");

    let two = best_available_dictionary(Some(("d", &good)), Some(&good), None, true)
        .expect_err("two dictionaries");
    assert_eq!(rendered(&two), refusal(&text, "both-dictionaries"));

    let none = best_available_dictionary(None, None, None, true).expect_err("no dictionary");
    assert_eq!(rendered(&none), refusal(&text, "no-dictionary"));

    let source = "updatevcfdictionary-dump/empty-dictionary.dict";
    let empty_source =
        best_available_dictionary(Some((source, &empty)), None, None, true).expect_err("empty");
    assert_eq!(rendered(&empty_source), refusal(&text, "empty-dictionary"));

    let source = "updatevcfdictionary-dump/no-length.dict";
    let missing = best_available_dictionary(Some((source, &no_length)), None, None, true)
        .expect_err("a sequence with no length");
    assert_eq!(rendered(&missing), refusal(&text, "no-length"));
}

/// The refusal that needs no records at all, and the argument that turns it off.
#[test]
fn an_input_with_a_dictionary_needs_replace() {
    let text = golden();
    let input_dictionary = dictionary(&text, "dictionary");
    let error = check_replace(&input_dictionary, false).expect_err("already has one");
    assert_eq!(rendered(&error), refusal(&text, "with-contigs-refused"));

    assert!(check_replace(&input_dictionary, true).is_ok());
    // The bare input has no dictionary of its own, so it never needs the argument.
    assert!(check_replace(&[], false).is_ok());
}

/// Every run that succeeded wrote every record it was given.
#[test]
fn every_accepted_run_writes_every_record() {
    let text = golden();
    let good = dictionary(&text, "dictionary");
    for (run, label) in [("bare", "bare"), ("with-contigs-replaced", "with-contigs")] {
        let records = input(&text, label);
        let variants: Vec<VariantContext> = records.iter().map(|line| parse_record(line)).collect();
        let (kept, error) = update_dictionary(&good, &variants);
        assert!(error.is_none(), "run {run} is accepted");
        let ours: Vec<String> = kept
            .into_iter()
            .map(|index| records[index].clone())
            .collect();
        assert_eq!(ours, written(&text, run), "written/{run}");
    }
}

/// The end is the END attribute where there is one, and the record is one base long.
#[test]
fn an_end_attribute_is_what_is_checked() {
    let text = golden();
    let good = dictionary(&text, "dictionary");
    let record = input(&text, "end-attribute")[0].clone();
    let variant = parse_record(&record);
    assert_eq!(variant.start, 100);
    assert_eq!(variant.stop, 250_000_001);
    let error = check_variant(&good, &variant).expect_err("past the end of chr1");
    assert!(
        error.message().contains("ends at a position (250000001)"),
        "{}",
        error.message()
    );
}
