//! Conformance for `CollectSVEvidence` against GATK 4.6.2.0, compared as all four evidence files
//! of every run and every refusal.
//!
//! Golden from `tools/readfilter-conformance/CollectSVEvidenceDump.java`.
//!
//! # What this suite is for
//!
//!  * **a discordant pair being written once, by one of its two reads**, and at an equal start by
//!    the first of the two SEEN, tracked by name;
//!  * **a read being split only if exactly one end is soft-clipped**, and which end deciding both
//!    the direction and the position;
//!  * **the match length summing every reference-consuming operator**;
//!  * **split positions being counted, and the two directions staying apart**;
//!  * **the site depth reading only biallelic SNPs at new loci**;
//!  * **three quality floors with three different defaults**;
//!  * **an interval no read reaches still being written**;
//!  * **an unpaired read crashing the discordant writer**;
//!  * **and each writer refusing a file name it could not read back.**

use gatk_corpus as corpus;
use gatk_tools::collect_sv_evidence::{
    bad_name_message, baf_sites, crashes_on_an_unpaired_read, depth_evidence, discordant_pairs,
    empty_intervals_message, is_soft_clipped, site_depths, split_position, split_reads,
    DepthEvidence, DiscordantPair, Read, Side, Site, SiteDepth, DEFAULT_DEPTH_EVIDENCE_MIN_MAPQ,
    DEFAULT_SITE_DEPTH_MIN_BASEQ, DEFAULT_SITE_DEPTH_MIN_MAPQ, NO_OUTPUT_MESSAGE,
    UNPAIRED_MATE_MESSAGE,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/collect_sv_evidence.txt.gz"),
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

fn refusal(text: &str, label: &str) -> (String, String) {
    let line = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = line.split_once(':').expect("a class and a message");
    (class.to_string(), unescape(message))
}

fn parse_cigar(text: &str) -> Vec<(char, i32)> {
    let mut out = Vec::new();
    let mut digits = String::new();
    for character in text.chars() {
        if character.is_ascii_digit() {
            digits.push(character);
        } else {
            out.push((character, digits.parse().expect("a length")));
            digits.clear();
        }
    }
    out
}

const CONTIGS: &[&str] = &["chr1", "chr2"];

/// The reads the golden reports, as the four writers see them. The SAM flags carry every boolean
/// the rules ask about.
fn reads(text: &str) -> Vec<Read> {
    section(text, "bam", "reads")
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let flags: u32 = columns[1].parse().expect("flags");
            let mate_contig = columns[6];
            Read {
                name: columns[0].to_string(),
                contig_index: CONTIGS
                    .iter()
                    .position(|name| *name == columns[2])
                    .expect("a known contig"),
                contig: columns[2].to_string(),
                start: columns[3].parse().expect("a start"),
                mapping_quality: columns[4].parse().expect("a mapping quality"),
                cigar: parse_cigar(columns[5]),
                paired: flags & 0x1 != 0,
                properly_paired: flags & 0x2 != 0,
                mate_unmapped: flags & 0x8 != 0,
                mate_contig_index: (mate_contig != ".")
                    .then(|| CONTIGS.iter().position(|name| *name == mate_contig))
                    .flatten(),
                mate_contig: (mate_contig != ".").then(|| mate_contig.to_string()),
                mate_start: (mate_contig != ".").then(|| columns[7].parse().expect("a mate start")),
                reverse_strand: flags & 0x10 != 0,
                mate_reverse_strand: flags & 0x20 != 0,
                secondary: flags & 0x100 != 0,
                duplicate: flags & 0x400 != 0,
                supplementary: flags & 0x800 != 0,
                unmapped: flags & 0x4 != 0,
                bases: columns[8].as_bytes().to_vec(),
                base_qualities: columns[9]
                    .split(',')
                    .map(|value| value.parse().expect("a base quality"))
                    .collect(),
            }
        })
        .collect()
}

fn sites(text: &str) -> Vec<Site> {
    section(text, "vcf", "sites")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Site {
                contig: columns[0].to_string(),
                position: columns[1].parse().expect("a position"),
                reference: columns[3].to_string(),
                alternates: columns[4].split(',').map(str::to_string).collect(),
            }
        })
        .collect()
}

fn intervals(text: &str) -> Vec<(String, i32, i32)> {
    section(text, "bed", "intervals")
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[0].to_string(),
                columns[1].parse().expect("a start"),
                columns[2].parse().expect("an end"),
            )
        })
        .collect()
}

fn measured_pe(text: &str, label: &str) -> Vec<DiscordantPair> {
    section(text, "out", &format!("{label}.pe"))
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            DiscordantPair {
                contig: columns[0].to_string(),
                position: columns[1].parse().expect("a position"),
                strand: columns[2] == "+",
                mate_contig: columns[3].to_string(),
                mate_position: columns[4].parse().expect("a mate position"),
                mate_strand: columns[5] == "+",
            }
        })
        .collect()
}

fn measured_sr(text: &str, label: &str) -> Vec<(String, i32, Side, i32)> {
    section(text, "out", &format!("{label}.sr"))
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[0].to_string(),
                columns[1].parse().expect("a position"),
                if columns[2] == "right" {
                    Side::Right
                } else {
                    Side::Left
                },
                columns[3].parse().expect("a count"),
            )
        })
        .collect()
}

fn measured_sd(text: &str, label: &str) -> Vec<SiteDepth> {
    section(text, "out", &format!("{label}.sd"))
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            SiteDepth {
                contig: columns[0].to_string(),
                position: columns[1].parse().expect("a position"),
                counts: [
                    columns[3].parse().expect("A"),
                    columns[4].parse().expect("C"),
                    columns[5].parse().expect("G"),
                    columns[6].parse().expect("T"),
                ],
            }
        })
        .collect()
}

fn measured_rd(text: &str, label: &str) -> Vec<DepthEvidence> {
    section(text, "out", &format!("{label}.rd"))
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            DepthEvidence {
                contig: columns[0].to_string(),
                start: columns[1].parse().expect("a start"),
                end: columns[2].parse().expect("an end"),
                count: columns[3].parse().expect("a count"),
            }
        })
        .collect()
}

#[test]
fn every_evidence_file_matches_the_golden() {
    let text = golden();
    let reads = reads(&text);
    let sites = sites(&text);
    let intervals = intervals(&text);
    let mut compared = 0;

    // The discordant pairs, from the two runs that wrote them.
    for label in ["default", "pe-only"] {
        assert_eq!(
            discordant_pairs(&reads),
            measured_pe(&text, label),
            "{label}"
        );
        compared += 1;
    }
    // The split reads, from the two runs that wrote them.
    for label in ["default", "sr-only"] {
        let produced: Vec<(String, i32, Side, i32)> = split_reads(&reads)
            .into_iter()
            .map(|record| (record.contig, record.position, record.side, record.count))
            .collect();
        assert_eq!(produced, measured_sr(&text, label), "{label}");
        compared += 1;
    }
    // The site depths, at the three quality settings that were run.
    for (label, min_mapq, min_baseq) in [
        (
            "default",
            DEFAULT_SITE_DEPTH_MIN_MAPQ,
            DEFAULT_SITE_DEPTH_MIN_BASEQ,
        ),
        ("high-mapq", 30, DEFAULT_SITE_DEPTH_MIN_BASEQ),
        ("low-baseq", DEFAULT_SITE_DEPTH_MIN_MAPQ, 0),
    ] {
        assert_eq!(
            site_depths(&reads, &sites, min_mapq, min_baseq),
            measured_sd(&text, label),
            "{label}"
        );
        compared += 1;
    }
    // The interval depths, at the two mapping-quality floors that were run.
    for (label, min_mapq) in [
        ("default", DEFAULT_DEPTH_EVIDENCE_MIN_MAPQ),
        ("high-mapq", 30),
    ] {
        assert_eq!(
            depth_evidence(&reads, &intervals, min_mapq),
            measured_rd(&text, label),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 9, "the evidence files that were written");
}

/// One of the two reads writes it, and at an equal start the FIRST one seen. Two different pairs
/// at one position are still two records, which is what makes the rule about names.
#[test]
fn a_discordant_pair_is_written_once() {
    let text = golden();
    let produced = discordant_pairs(&reads(&text));
    // The pair whose two reads are 500 apart: one record, at the smaller start.
    assert_eq!(
        produced
            .iter()
            .filter(|record| record.position == 1999)
            .count(),
        1
    );
    assert!(!produced.iter().any(|record| record.position == 2499));
    // The pair at one position: one record, and its strands are its own read's.
    let same: Vec<&DiscordantPair> = produced
        .iter()
        .filter(|record| record.position == 2999)
        .collect();
    assert_eq!(same.len(), 1);
    assert!(same[0].strand, "the forward read wrote it");
    assert!(!same[0].mate_strand, "and its mate is reverse");
    // Two DIFFERENT pairs at one position: two records.
    assert_eq!(
        produced
            .iter()
            .filter(|record| record.position == 3999)
            .count(),
        2
    );
    // Across contigs, written by the read on the smaller index only.
    assert_eq!(
        produced
            .iter()
            .filter(|record| record.mate_contig == "chr2")
            .count(),
        1
    );
}

/// Exactly one end, and which end decides both the direction and the position. The match length
/// counts deletions.
#[test]
fn a_read_is_split_only_if_exactly_one_end_is_clipped() {
    let text = golden();
    let reads = reads(&text);
    let of = |name: &str| {
        reads
            .iter()
            .find(|read| read.name == name)
            .unwrap_or_else(|| panic!("{name}"))
    };
    assert!(is_soft_clipped(of("clipRight")));
    assert!(is_soft_clipped(of("clipLeft")));
    assert!(!is_soft_clipped(of("clipBoth")), "both ends is not a split");
    assert_eq!(split_position(of("clipRight")), (10020, Side::Right));
    assert_eq!(split_position(of("clipLeft")), (11000, Side::Left));
    // 10M5D10M5S: the deletion is counted, so the position is 25 out rather than 20.
    assert_eq!(split_position(of("withDeletion")), (15025, Side::Right));

    // And the golden shows the read clipped at both ends contributing nothing.
    let produced = measured_sr(&text, "default");
    assert!(!produced
        .iter()
        .any(|(_, position, _, _)| (13999..14030).contains(position)));
    // Two reads clipped at one place are a count of two.
    assert!(produced.contains(&("chr1".to_string(), 12019, Side::Right, 1 + 1)));
    // The same position from the two directions stays two records.
    assert!(produced.contains(&("chr1".to_string(), 13019, Side::Left, 1)));
    assert!(produced.contains(&("chr1".to_string(), 13019, Side::Right, 1)));
}

/// A repeated position, a triallelic site and an indel are all walked past.
#[test]
fn the_site_depth_reads_only_biallelic_snps_at_new_loci() {
    let text = golden();
    let all = sites(&text);
    assert_eq!(all.len(), 6);
    let taken = baf_sites(&all);
    assert_eq!(
        taken.iter().map(|site| site.position).collect::<Vec<i32>>(),
        vec![20005, 20012, 50000]
    );
    // Which is exactly what the golden wrote, one record each, zero-based.
    assert_eq!(
        measured_sd(&text, "default")
            .iter()
            .map(|record| record.position)
            .collect::<Vec<i32>>(),
        vec![20004, 20011, 49999]
    );
    // The locus no read covers is still written, with four zeros.
    assert_eq!(measured_sd(&text, "default")[2].counts, [0, 0, 0, 0]);
}

/// Three floors, three defaults, and one filters a read where another filters a base.
#[test]
fn the_three_quality_floors_have_three_defaults() {
    let text = golden();
    assert_eq!(DEFAULT_SITE_DEPTH_MIN_MAPQ, 30);
    assert_eq!(DEFAULT_SITE_DEPTH_MIN_BASEQ, 20);
    assert_eq!(DEFAULT_DEPTH_EVIDENCE_MIN_MAPQ, 0);

    // The low-mapping-quality read is outside the site depth by default and inside the read depth,
    // so raising the depth floor to 30 is what removes it.
    assert_eq!(measured_rd(&text, "default")[0].count, 4);
    assert_eq!(measured_rd(&text, "high-mapq")[0].count, 3);
    // Raising the SITE depth floor to the same 30 changes nothing, because 30 is already the
    // default there.
    assert_eq!(
        measured_sd(&text, "default"),
        measured_sd(&text, "high-mapq")
    );
    // The low-base-quality read needs the floor LOWERED before its base appears at all.
    assert_eq!(measured_sd(&text, "default")[0].counts, [0, 2, 0, 0]);
    assert_eq!(measured_sd(&text, "low-baseq")[0].counts, [0, 3, 0, 0]);
}

/// An interval nothing reaches is written with a count of zero.
#[test]
fn an_interval_no_read_reaches_is_still_written() {
    let text = golden();
    let produced = measured_rd(&text, "default");
    assert_eq!(produced.len(), 3);
    assert_eq!(produced[2].count, 0);
    assert_eq!(produced[2].start, 99999);
}

/// It reports isProperlyPaired as false, so the writer asks it for a mate it does not have.
#[test]
fn an_unpaired_read_crashes_the_discordant_writer() {
    let text = golden();
    let (class, message) = refusal(&text, "unpaired");
    assert_eq!(class, "java.lang.IllegalStateException");
    assert_eq!(message, UNPAIRED_MATE_MESSAGE);
    let lone = Read {
        name: "lone".to_string(),
        contig_index: 0,
        contig: "chr1".to_string(),
        start: 1000,
        mapping_quality: 60,
        cigar: vec![('M', 20)],
        paired: false,
        properly_paired: false,
        mate_unmapped: false,
        mate_contig_index: None,
        mate_contig: None,
        mate_start: None,
        reverse_strand: false,
        mate_reverse_strand: false,
        supplementary: false,
        secondary: false,
        duplicate: false,
        unmapped: false,
        bases: b"ACGTACGTACGTACGTACGT".to_vec(),
        base_qualities: vec![30; 20],
    };
    assert!(crashes_on_an_unpaired_read(&lone));
    // Every read in the main fixture is properly paired, which is why the other runs survive.
    assert!(!reads(&text).iter().any(crashes_on_an_unpaired_read));
}

/// No output at all, four wrong file names each with its own wording, and an empty interval file.
#[test]
fn the_seven_refusals() {
    let text = golden();
    assert_eq!(refusal(&text, "no-output").1, NO_OUTPUT_MESSAGE);
    for (label, kind) in [
        ("bad-pe-name", "pe"),
        ("bad-sr-name", "sr"),
        ("bad-sd-name", "sd"),
        ("bad-rd-name", "rd"),
    ] {
        let (class, message) = refusal(&text, label);
        assert_eq!(
            class, "org.broadinstitute.hellbender.exceptions.UserException",
            "{label}"
        );
        assert_eq!(
            message,
            bad_name_message(kind, "<dir>/wrong.txt"),
            "{label}"
        );
    }
    // The four messages are four different messages, not one wording repeated.
    let wordings: Vec<String> = ["pe", "sr", "sd", "rd"]
        .iter()
        .map(|kind| bad_name_message(kind, "x"))
        .collect();
    let mut unique = wordings.clone();
    unique.sort();
    unique.dedup();
    assert_eq!(unique.len(), 4);

    assert_eq!(
        refusal(&text, "empty-intervals").1,
        empty_intervals_message("<dir>/empty.bed")
    );
}
