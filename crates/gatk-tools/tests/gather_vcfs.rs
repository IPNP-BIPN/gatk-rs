//! Conformance for `GatherVcfsCloud` against GATK 4.6.2.0, compared as the records every gather
//! writes and as the class and message of every refusal.
//!
//! Golden from `tools/readfilter-conformance/GatherVcfsDump.java`.
//!
//! # What this suite is for
//!
//!  * **two order checks, two exception classes**, one on first records and one on the last written;
//!  * **`--disable-contig-ordering-check` writes an unordered file** rather than refusing;
//!  * **`--ignore-safety-checks` relabels a genotype** rather than dropping it;
//!  * **and a missing dictionary is the indexer's refusal**, not the validation's.

use gatk_corpus as corpus;
use gatk_tools::gather_vcfs::{gather, Arguments, GatherError, GatherType, Shard};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gather_vcfs.txt.gz"),
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

/// One input file of the golden as a shard, named as the reference names it.
fn shard(text: &str, label: &str) -> Shard {
    let whole = rows(text, "input")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no input {label}"))[1]
        .to_string();
    let content = unescape(&whole);
    let mut dictionary = Vec::new();
    let mut samples = Vec::new();
    let mut records = Vec::new();
    for line in content.lines() {
        if let Some(rest) = line.strip_prefix("##contig=<ID=") {
            dictionary.push(rest.split(',').next().expect("a contig").to_string());
        } else if let Some(rest) = line.strip_prefix("#CHROM") {
            samples = rest
                .split('\t')
                .skip(9)
                .map(|name| name.to_string())
                .collect();
        } else if !line.starts_with('#') {
            let field: Vec<&str> = line.split('\t').collect();
            records.push((field[0].to_string(), field[1].parse().expect("a position")));
        }
    }
    Shard {
        name: format!("file:///work/gathervcfs-dump/{label}.vcf"),
        dictionary,
        samples,
        records,
    }
}

/// The records one gather wrote, as `contig:position`.
fn written(text: &str, run: &str) -> Vec<String> {
    rows(text, "vcfline")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            format!("{}:{}", field[0], field[1])
        })
        .collect()
}

/// How the golden holds a refusal: the class, a colon, and whatever prefix the class adds to
/// `getMessage`. `UserException$BadInput` prepends "Bad input: ", which is the class's doing and
/// not the tool's, so it belongs here rather than in the port's message.
fn rendered(error: &GatherError) -> String {
    let prefix = match error {
        GatherError::BlockCopyImpossible => "Bad input: ",
        _ => "",
    };
    format!("{}:{}{}", error.java_class(), prefix, error.message())
}

fn refusal(text: &str, run: &str) -> Option<String> {
    rows(text, "error")
        .into_iter()
        .find(|row| row[0] == run)
        .map(|row| unescape(row[1]))
}

/// The inputs and arguments of each run, which the golden does not carry.
fn setup(run: &str) -> (Vec<&'static str>, Arguments) {
    let base = Arguments::default;
    match run {
        "three-shards" | "conventional" => (vec!["first", "second", "third"], base()),
        "one-shard" => (vec!["first"], base()),
        "block-refused" => (
            vec!["first", "second", "third"],
            Arguments {
                gather_type: GatherType::Block,
                ..base()
            },
        ),
        "out-of-order" => (vec!["second", "first"], base()),
        "out-of-order-check-disabled" => (
            vec!["second", "first"],
            Arguments {
                disable_contig_ordering_check: true,
                ..base()
            },
        ),
        "contig-out-of-order" => (vec!["third", "first"], base()),
        "contig-out-of-order-check-disabled" => (
            vec!["third", "first"],
            Arguments {
                disable_contig_ordering_check: true,
                ..base()
            },
        ),
        "different-samples" => (vec!["first", "other-sample"], base()),
        "different-samples-ignored" => (
            vec!["first", "other-sample"],
            Arguments {
                ignore_safety_checks: true,
                ..base()
            },
        ),
        "no-dictionary" => (vec!["no-dictionary", "second"], base()),
        "overlapping-records" => (vec!["first", "overlapping"], base()),
        other => panic!("no setup for {other}"),
    }
}

const RUNS: [&str; 12] = [
    "three-shards",
    "one-shard",
    "conventional",
    "block-refused",
    "out-of-order",
    "out-of-order-check-disabled",
    "contig-out-of-order",
    "contig-out-of-order-check-disabled",
    "different-samples",
    "different-samples-ignored",
    "no-dictionary",
    "overlapping-records",
];

#[test]
fn every_gather_is_the_reference_s() {
    let text = golden();
    for run in RUNS {
        let (labels, arguments) = setup(run);
        let shards: Vec<Shard> = labels.iter().map(|label| shard(&text, label)).collect();
        let result = gather(&shards, &arguments);

        match refusal(&text, run) {
            Some(expected) => {
                let error = result.expect_err(run);
                assert_eq!(rendered(&error), expected, "error/{run}");
            }
            None => {
                let ours: Vec<String> = result
                    .unwrap_or_else(|error| panic!("{run}: {}", error.message()))
                    .into_iter()
                    .map(|(shard, record)| {
                        let (contig, at) = &shards[shard].records[record];
                        format!("{contig}:{at}")
                    })
                    .collect();
                assert_eq!(ours, written(&text, run), "written/{run}");
            }
        }
    }
}

/// The same pair of files, two checks, two classes.
#[test]
fn the_two_order_checks_are_two_different_refusals() {
    let text = golden();
    let first = refusal(&text, "out-of-order").expect("a refusal");
    let second = refusal(&text, "overlapping-records").expect("a refusal");
    assert!(first.starts_with("java.lang.IllegalArgumentException:"));
    assert!(second.starts_with("java.lang.IllegalStateException:"));
    // The second names both positions, which the first does not.
    assert!(second.contains("is at chr1:150 but last variant"));
}

/// The flag writes a file whose records are not in dictionary order.
#[test]
fn the_weakened_check_writes_an_unordered_file() {
    let text = golden();
    let records = written(&text, "contig-out-of-order-check-disabled");
    assert_eq!(records, vec!["chr2:100", "chr1:100", "chr1:200"]);

    // And our port produces the same order rather than refusing.
    let (labels, arguments) = setup("contig-out-of-order-check-disabled");
    let shards: Vec<Shard> = labels.iter().map(|label| shard(&text, label)).collect();
    let ours = gather(&shards, &arguments).expect("accepted");
    assert_eq!(ours.len(), 3);
}

/// Ignoring the checks writes the other sample's record under the first file's sample name.
#[test]
fn ignoring_the_safety_checks_relabels_a_genotype() {
    let text = golden();
    let refused = refusal(&text, "different-samples").expect("a refusal");
    assert!(refused.contains("Samples unique to first file: [s0]"));
    assert!(refused.contains("[s1]"));

    // The same pair with the flag writes three records under a header declaring only s0.
    assert_eq!(
        written(&text, "different-samples-ignored"),
        vec!["chr1:100", "chr1:200", "chr1:400"]
    );
    let (labels, arguments) = setup("different-samples-ignored");
    let shards: Vec<Shard> = labels.iter().map(|label| shard(&text, label)).collect();
    assert_eq!(shards[1].samples, vec!["s1".to_string()]);
    assert!(gather(&shards, &arguments).is_ok());
}

/// The dictionary refusal comes from the indexer, after the validation has passed.
#[test]
fn a_missing_dictionary_is_the_indexers_refusal() {
    let text = golden();
    let expected = refusal(&text, "no-dictionary").expect("a refusal");
    assert!(expected.starts_with("org.broadinstitute.hellbender.exceptions.UserException:"));
    assert!(expected.contains("##contig lines"));

    let (labels, arguments) = setup("no-dictionary");
    let shards: Vec<Shard> = labels.iter().map(|label| shard(&text, label)).collect();
    let error = gather(&shards, &arguments).unwrap_err();
    assert_eq!(error, GatherError::NoDictionary);
}
