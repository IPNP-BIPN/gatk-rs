//! Conformance for `LocalAssembler` against GATK 4.6.2.0, compared as the graph it writes and the
//! sequences it reads off it.
//!
//! Golden from `tools/readfilter-conformance/LocalAssemblerDump.java`, whose reference is drawn by
//! a linear congruential generator: a periodic one collapses the assembly into a single forty-base
//! contig looping on itself.
//!
//! The assembly is not ported. What the port is asked for is everything the two files are made of:
//! given the GFA's own segments and paths, it must spell out the FASTA the tool wrote.
//!
//! # What this suite is for
//!
//!  * **a traversal's sequence being its contigs joined on a thirty-base overlap**;
//!  * **the FASTA naming a traversal `<assembly>_t<n>` with a one-based counter**;
//!  * **the path notation using `+` and `RC` while the GFA uses `+` and `-`**;
//!  * **a traversal being written in either orientation**;
//!  * **the edge coordinates being the overlap's place in each contig**;
//!  * **a case that assembles nothing writing an empty FASTA rather than none**;
//!  * **and the thin-observation floor and the scaffolding changing nothing here.**

use gatk_corpus as corpus;
use gatk_tools::local_assembler::{
    fasta_header, gfa_edge, gfa_path, parse_gfa_path, parse_traversal_name, reverse_complement,
    reverse_traversal, traversal_name, traversal_sequence, traversal_sequence_length, Step,
    GFA_HEADER, KMER_SIZE, MIN_THIN_OBS_DEFAULT,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/local_assembler.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn payload(text: &str, kind: &str, label: &str) -> Option<String> {
    let prefix = format!("{kind}\t{label}=");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// The GFA's segments of one case: the contig id, its sequence and its three counts.
fn segments(text: &str, label: &str) -> Vec<(String, String, i32)> {
    payload(text, "gfa", label)
        .unwrap_or_else(|| panic!("{label} wrote a gfa"))
        .lines()
        .filter(|line| line.starts_with("S\t"))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[1].to_string(),
                columns[3].to_string(),
                columns[4]
                    .strip_prefix("MO:i:")
                    .expect("a max observation")
                    .parse()
                    .expect("a number"),
            )
        })
        .collect()
}

/// The GFA's paths of one case, as their steps.
fn paths(text: &str, label: &str) -> Vec<Vec<Step>> {
    payload(text, "gfa", label)
        .unwrap_or_else(|| panic!("{label} wrote a gfa"))
        .lines()
        .filter(|line| line.starts_with("O\t"))
        .map(|line| parse_gfa_path(line.split('\t').nth(2).expect("the references")))
        .collect()
}

/// The FASTA's records of one case: the header's name and the sequence.
fn records(text: &str, label: &str) -> Vec<(String, String, String)> {
    let fasta = payload(text, "fasta", label).unwrap_or_else(|| panic!("{label} wrote a fasta"));
    let lines: Vec<&str> = fasta.lines().filter(|line| !line.is_empty()).collect();
    lines
        .chunks(2)
        .map(|pair| {
            let header = pair[0].strip_prefix('>').expect("a header");
            let (identifier, name) = header.split_once(' ').expect("a name");
            (
                identifier.to_string(),
                name.to_string(),
                pair[1].to_string(),
            )
        })
        .collect()
}

const ASSEMBLED: &[&str] = &[
    "one-read",
    "overlapping-reads",
    "disjoint-reads",
    "bubble",
    "n-in-the-middle",
    "thin-observations-ten",
    "no-scaffolding",
];

/// Every assembled case's FASTA sequences are spelled out from its own GFA segments, and the
/// traversal named in each header walks those segments.
#[test]
fn every_fasta_is_spelled_out_of_its_gfa() {
    let text = golden();
    for label in ASSEMBLED {
        let segments = segments(&text, label);
        let sequence_of = |contig: &str| {
            segments
                .iter()
                .find(|(id, _, _)| id == contig)
                .map(|(_, sequence, _)| sequence.clone())
        };
        for (index, (identifier, name, sequence)) in records(&text, label).iter().enumerate() {
            assert_eq!(
                *identifier,
                format!("{label}_t{}", index + 1),
                "{label} header"
            );
            let steps = parse_traversal_name(name);
            assert_eq!(traversal_name(&steps), *name, "{label} name round trip");
            assert_eq!(
                traversal_sequence(&steps, sequence_of),
                *sequence,
                "{label} sequence"
            );
        }
    }
}

/// The FASTA header is the assembly name, an underscore, `t` and a ONE-based counter.
#[test]
fn the_fasta_header_counts_from_one() {
    let text = golden();
    assert_eq!(fasta_header(Some("one-read"), 1, "c1"), ">one-read_t1 c1");
    assert_eq!(fasta_header(None, 2, "c1+c2"), ">t2 c1+c2");
    let first = &records(&text, "disjoint-reads")[0];
    assert_eq!(first.0, "disjoint-reads_t1");
    let second = &records(&text, "disjoint-reads")[1];
    assert_eq!(second.0, "disjoint-reads_t2");
}

/// The two notations differ: the FASTA joins with `+` and marks a reverse-complemented contig
/// `RC`, while the GFA separates with spaces and marks it `-`.
#[test]
fn the_two_notations_differ() {
    let steps = parse_traversal_name("c1+c3RC+c4RC");
    assert_eq!(traversal_name(&steps), "c1+c3RC+c4RC");
    assert_eq!(gfa_path(&steps), "O\t*\tc1+ c3- c4-");
    assert_eq!(parse_gfa_path("c1+ c3- c4-"), steps);
    let text = golden();
    let gfa = payload(&text, "gfa", "bubble").expect("the bubble gfa");
    assert!(gfa.contains("O\t*\tc1+ c3- c4-"), "{gfa}");
    let names: Vec<String> = records(&text, "bubble")
        .into_iter()
        .map(|(_, name, _)| name)
        .collect();
    assert!(names.contains(&"c1+c3RC+c4RC".to_string()), "{names:?}");
}

/// A traversal may be written in either orientation: the bubble's second FASTA record walks the
/// reverse of the graph's own second path.
#[test]
fn a_traversal_may_be_written_reversed() {
    let text = golden();
    let paths = paths(&text, "bubble");
    let names: Vec<String> = records(&text, "bubble")
        .into_iter()
        .map(|(_, name, _)| name)
        .collect();
    assert_eq!(names.len(), 2);
    assert_eq!(paths.len(), 2);
    // One record walks a path as the graph wrote it, the other walks its reverse.
    let forwards: Vec<String> = paths.iter().map(|steps| traversal_name(steps)).collect();
    let backwards: Vec<String> = paths
        .iter()
        .map(|steps| traversal_name(&reverse_traversal(steps)))
        .collect();
    for name in &names {
        assert!(
            forwards.contains(name) || backwards.contains(name),
            "{name} is neither {forwards:?} nor {backwards:?}"
        );
    }
    assert!(names.iter().any(|name| backwards.contains(name)));
    // Reversing twice is the identity.
    let steps = parse_gfa_path("c1+ c2- c4-");
    assert_eq!(reverse_traversal(&reverse_traversal(&steps)), steps);
}

/// The contigs overlap by thirty bases, which is the kmer size less one, so three of 59, 61 and 60
/// spell out 120 rather than 180.
#[test]
fn the_contigs_overlap_by_thirty() {
    let text = golden();
    assert_eq!(KMER_SIZE, 31);
    let segments = segments(&text, "bubble");
    let lengths: Vec<usize> = ["c1", "c3", "c4"]
        .iter()
        .map(|id| {
            segments
                .iter()
                .find(|(name, _, _)| name == id)
                .expect("a segment")
                .1
                .len()
        })
        .collect();
    assert_eq!(lengths, vec![59, 61, 60]);
    let spelled: usize = lengths[0] + (lengths[1] - 30) + (lengths[2] - 30);
    assert_eq!(spelled, 120);
    let record = records(&text, "bubble")
        .into_iter()
        .find(|(_, name, _)| name == "c1+c3RC+c4RC")
        .expect("the record");
    assert_eq!(record.2.len(), 120);
    // Which is what counting kmers and adding the overlap back once gives.
    let kmer_counts: Vec<usize> = lengths
        .iter()
        .map(|length| length - KMER_SIZE + 1)
        .collect();
    assert_eq!(traversal_sequence_length(&kmer_counts), 120);
}

/// The edge coordinates are the overlap's place in each contig, the first written with a `$`.
#[test]
fn the_edge_names_the_overlap_in_both_contigs() {
    let text = golden();
    let from = Step {
        contig: "c1".to_string(),
        reverse_complemented: false,
    };
    let to = Step {
        contig: "c2".to_string(),
        reverse_complemented: true,
    };
    assert_eq!(
        gfa_edge(&from, &to, 59),
        "E\t*\tc1+\tc2-\t29\t59$\t0\t30\t30M"
    );
    let gfa = payload(&text, "gfa", "bubble").expect("the bubble gfa");
    assert!(gfa.contains("E\t*\tc1+\tc2-\t29\t59$\t0\t30\t30M"), "{gfa}");
    assert!(gfa.starts_with(GFA_HEADER), "{gfa}");
}

/// A case that assembles nothing writes an empty FASTA rather than none at all.
#[test]
fn a_case_that_assembles_nothing_writes_an_empty_fasta() {
    let text = golden();
    for label in ["read-shorter-than-k", "low-quality-base", "empty-interval"] {
        let fasta = payload(&text, "fasta", label).unwrap_or_else(|| panic!("{label}"));
        assert!(fasta.trim().is_empty(), "{label}: {fasta}");
    }
}

/// Raising the thin-observation floor to ten changes nothing when the reads already meet it, and
/// neither does turning the scaffolding off: both cases match the plain overlapping run.
#[test]
fn the_floor_and_the_scaffolding_change_nothing_here() {
    let text = golden();
    assert_eq!(MIN_THIN_OBS_DEFAULT, 4);
    let plain = records(&text, "overlapping-reads");
    for label in ["thin-observations-ten", "no-scaffolding"] {
        let other = records(&text, label);
        assert_eq!(
            other.iter().map(|record| &record.2).collect::<Vec<_>>(),
            plain.iter().map(|record| &record.2).collect::<Vec<_>>(),
            "{label}"
        );
    }
}

/// The reverse complement, which a reverse-complemented step reads its contig through.
#[test]
fn the_reverse_complement_turns_the_sequence_round() {
    assert_eq!(reverse_complement("ACGT"), "ACGT");
    assert_eq!(reverse_complement("AAAA"), "TTTT");
    assert_eq!(reverse_complement("ACGTTT"), "AAACGT");
    assert_eq!(reverse_complement(""), "");
    // A base the alphabet does not name is carried through untouched.
    assert_eq!(reverse_complement("ANT"), "ANT");
}
