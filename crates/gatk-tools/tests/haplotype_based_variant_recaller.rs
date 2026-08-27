//! Conformance for `HaplotypeBasedVariantRecaller` against GATK 4.6.2.0, compared as the whole
//! matrix file of every run.
//!
//! Golden from `tools/readfilter-conformance/HaplotypeBasedVariantRecallerDump.java`.
//!
//! The PairHMM is not measured or ported. The likelihoods are read off the golden's own matrix
//! and everything around them is rebuilt: the header line, the bases column, the unclipped
//! offset, the sort and the lines that never appear.
//!
//! # What this suite is for
//!
//!  * **a haplotype being a record named `HC_`**;
//!  * **haplotypes grouping by identical span**;
//!  * **the group chosen being the one that centres the variant, ties going to the first**;
//!  * **the header line omitting the end for a one-base variant and for a MIXED one**;
//!  * **the lines being sorted by the LAST likelihood and not by the best**;
//!  * **a read that does not span the whole variant being dropped**;
//!  * **a variant inside a deletion NOT being dropped but coming out from the wrong base**;
//!  * **the offset being counted from the other end on the reverse strand**;
//!  * **and a duplicate read never reaching the matrix at all.**

use gatk_corpus as corpus;
use gatk_tools::haplotype_based_variant_recaller::{
    best_group, fitness_score, groups, header_line, is_haplotype_record, matrix_line,
    offset_on_read, parse_cigar, sorted_lines, variant_bases, variant_block, HaplotypeRecord,
    MatrixLine, Read, Span, VariantKind, HAPLOTYPE_NAME_PREFIX,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/haplotype_based_variant_recaller.txt.gz"),
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

/// The reads of one of the golden's BAMs, read back from its SAM text.
fn reads(text: &str, name: &str) -> Vec<Read> {
    section(text, "sam", name)
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let start: i32 = columns[3].parse().expect("a position");
            let cigar = parse_cigar(columns[5]);
            let reference_length: i32 = cigar
                .iter()
                .filter(|element| element.consumes_reference_bases())
                .map(|element| element.length)
                .sum();
            let flags: i32 = columns[1].parse().expect("a flag");
            let bases = columns[9].as_bytes().to_vec();
            Read {
                name: columns[0].to_string(),
                span: Span::new(columns[2], start, start + reference_length - 1),
                cigar,
                bases,
                is_duplicate: flags & 1024 != 0,
                is_reverse: flags & 16 != 0,
                mapping_quality: columns[4].parse().expect("a quality"),
                // Nothing in this fixture is flow-based.
                key_length: 0,
                sample: "sm1".to_string(),
                // The reads are trimmed to the haplotype before the line is built, and a trimmed
                // read's unclipped ends are its own: see the soft-clipped read's own test.
                unclipped_start: start,
                unclipped_end: start + reference_length - 1,
            }
        })
        .collect()
}

/// The haplotype records of the golden's haplotype BAM.
fn haplotypes(text: &str) -> Vec<HaplotypeRecord> {
    reads(text, "haplotypes")
        .into_iter()
        .map(|read| HaplotypeRecord {
            name: read.name,
            span: read.span,
        })
        .collect()
}

/// One variant block the golden carries: its header line and its lines.
#[derive(Debug, Clone)]
struct Block {
    header: String,
    lines: Vec<String>,
}

fn blocks(text: &str, label: &str) -> Vec<Block> {
    let mut blocks: Vec<Block> = Vec::new();
    for line in section(text, "out", label).lines() {
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            blocks.push(Block {
                header: line.to_string(),
                lines: Vec::new(),
            });
        } else {
            blocks
                .last_mut()
                .expect("a header before a line")
                .lines
                .push(line.to_string());
        }
    }
    blocks
}

/// The likelihood columns of one written line, and the read it names.
fn likelihoods(line: &str) -> (String, Vec<f64>) {
    let fields: Vec<&str> = line.split(' ').collect();
    // name, keyLength, duplicate, reverse, mappingQuality, ...likelihoods..., bases, offset,
    // sample.
    let values = fields[5..fields.len() - 3]
        .iter()
        .map(|value| value.parse().expect("a likelihood"))
        .collect();
    (fields[0].to_string(), values)
}

/// The start one header line names, which is the one field it never leaves out.
fn header_start(header: &str) -> i32 {
    let place = header
        .trim_start_matches('#')
        .split(' ')
        .next()
        .expect("a place");
    let (_, rest) = place.split_once(':').expect("a contig");
    rest.split('-')
        .next()
        .expect("a start")
        .parse()
        .expect("a number")
}

/// The spans the ALLELES VCF declares, by start.
///
/// The header line is not enough to recover them: a MIXED site's end is left out of it, so a
/// two-base mixed variant reads exactly like a one-base one. The bases column is cut from the
/// VCF's span rather than from the header's, which is why the two can disagree.
fn variant_spans(text: &str) -> Vec<(i32, Span)> {
    section(text, "vcf", "alleles")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let start: i32 = columns[1].parse().expect("a position");
            let end = start + columns[3].len() as i32 - 1;
            (start, Span::new(columns[0], start, end))
        })
        .collect()
}

/// Every run's whole matrix file, rebuilt from the golden's own likelihoods.
///
/// The read that was soft-clipped is left out of the rebuild: it is the one read whose unclipped
/// offset the trimming changed, and its own test covers that.
#[test]
fn every_matrix_file_matches_the_golden() {
    let text = golden();
    let reads = reads(&text, "reads");
    let mut compared = 0;
    for label in [
        "whole-reference",
        "one-interval",
        "empty-interval",
        "mapping-quality-30",
    ] {
        let mut rebuilt = String::new();
        let spans = variant_spans(&text);
        for block in blocks(&text, label) {
            let start = header_start(&block.header);
            let variant = spans
                .iter()
                .find(|(at, _)| *at == start)
                .map(|(_, span)| span.clone())
                .unwrap_or_else(|| panic!("the alleles vcf carries {start}"));
            let lines: Vec<MatrixLine> = block
                .lines
                .iter()
                .filter(|line| !line.starts_with("r-clipped "))
                .map(|line| {
                    let (name, values) = likelihoods(line);
                    let read = reads
                        .iter()
                        .find(|read| read.name == name)
                        .unwrap_or_else(|| panic!("the fixture carries {name}"));
                    matrix_line(read, &variant, &values)
                        .unwrap_or_else(|| panic!("{label}: {name} at {variant}"))
                })
                .collect();
            let expected: Vec<String> = block
                .lines
                .iter()
                .filter(|line| !line.starts_with("r-clipped "))
                .cloned()
                .collect();
            rebuilt.push_str(&variant_block(&block.header, &lines));
            assert_eq!(sorted_lines(&lines), expected, "{label} at {variant}");
        }
        // The file is those blocks and nothing else.
        let written: String = section(&text, "out", label)
            .lines()
            .filter(|line| !line.starts_with("r-clipped "))
            .map(|line| format!("{line}\n"))
            .collect();
        assert_eq!(rebuilt, written, "{label}");
        compared += 1;
    }
    assert_eq!(compared, 4, "the runs the port reproduces");
}

/// Any other record in the haplotype BAM is passed over however well it fits.
#[test]
fn a_haplotype_is_a_record_named_hc() {
    let text = golden();
    let records = haplotypes(&text);
    assert_eq!(HAPLOTYPE_NAME_PREFIX, "HC_");
    // The fixture holds one record that is not a haplotype, over the same span as a group.
    let intruder = records
        .iter()
        .find(|record| record.name == "not_a_haplotype")
        .expect("the intruder");
    assert!(!is_haplotype_record(&intruder.name));
    assert!(
        records
            .iter()
            .filter(|r| is_haplotype_record(&r.name))
            .count()
            >= 6
    );
    // It sits over the same span as the centred group and is still left out of it.
    let grouped = groups(&records);
    let centred = grouped
        .iter()
        .find(|group| group[0].span == intruder.span)
        .expect("the centred group");
    assert!(!centred
        .iter()
        .any(|record| record.name == "not_a_haplotype"));
    assert_eq!(centred.len(), 2);
}

/// A record whose span differs closes the group and opens a new one.
#[test]
fn haplotypes_group_by_identical_span() {
    let text = golden();
    let grouped = groups(&haplotypes(&text));
    let spans: Vec<String> = grouped
        .iter()
        .map(|group| group[0].span.to_string())
        .collect();
    assert_eq!(
        spans,
        vec!["chr1:1000-1099", "chr1:1040-1239", "chr1:1250-1499"]
    );
    for group in &grouped {
        assert!(group.iter().all(|record| record.span == group[0].span));
    }
    // Two runs of the same span separated by a third are two groups, not one.
    let split = vec![
        HaplotypeRecord {
            name: "HC_a".to_string(),
            span: Span::new("chr1", 1, 10),
        },
        HaplotypeRecord {
            name: "HC_b".to_string(),
            span: Span::new("chr1", 5, 15),
        },
        HaplotypeRecord {
            name: "HC_c".to_string(),
            span: Span::new("chr1", 1, 10),
        },
    ];
    assert_eq!(groups(&split).len(), 3);
}

/// The one that centres the variant best, ties going to the first.
#[test]
fn the_group_chosen_centres_the_variant() {
    let text = golden();
    let grouped = groups(&haplotypes(&text));
    // The variant at 1050 falls in two groups, and the centred one wins.
    let variant = Span::new("chr1", 1050, 1050);
    let best = best_group(&variant, &grouped).expect("a group");
    assert_eq!(best[0].span, Span::new("chr1", 1000, 1099));
    assert!(fitness_score(&variant, best) > fitness_score(&variant, &grouped[1]));
    // Which is the span the golden's own header line names.
    let header = &blocks(&text, "whole-reference")[0].header;
    assert!(header.contains("chr1:1000-1099"), "{header}");
    // A variant flush against an end scores above zero rather than dividing by zero.
    let flush = Span::new("chr1", 1000, 1000);
    let score = fitness_score(&flush, &grouped[0]);
    assert!(score > 0.0 && score < 0.1, "{score}");
    // A perfectly centred variant scores one.
    let centred = Span::new("chr1", 1050, 1049);
    assert_eq!(fitness_score(&centred, &grouped[0]), 1.0);
    // The comparison is strict, so a later group of equal fitness does not displace the first.
    let twins = vec![grouped[0].clone(), grouped[0].clone()];
    assert!(std::ptr::eq(
        best_group(&variant, &twins).expect("a group"),
        &twins[0]
    ));
    // An empty list yields nothing.
    assert!(best_group(&variant, &[]).is_none());
}

/// For a one-base variant, and for a MIXED one however long it is.
#[test]
fn the_header_omits_the_end_for_a_mixed_variant() {
    let text = golden();
    let headers: Vec<String> = blocks(&text, "whole-reference")
        .into_iter()
        .map(|block| block.header)
        .collect();
    assert_eq!(headers.len(), 4);
    // The SNP and the insertion are one base, so neither carries an end.
    assert!(headers[0].starts_with("#chr1:1050 "), "{}", headers[0]);
    assert!(headers[1].starts_with("#chr1:1199 "), "{}", headers[1]);
    // The deletion is four bases and carries one.
    assert!(headers[2].starts_with("#chr1:1299-1302 "), "{}", headers[2]);
    // The mixed site is two bases and carries none.
    assert!(headers[3].starts_with("#chr1:1400 "), "{}", headers[3]);
    let span = Span::new("chr1", 1250, 1499);
    assert_eq!(
        header_line(
            "chr1",
            1400,
            1401,
            VariantKind::Mixed,
            &span,
            &["T*".into(), "A".into()]
        ),
        headers[3]
    );
    // The same site as anything but MIXED does carry its end.
    assert!(header_line(
        "chr1",
        1400,
        1401,
        VariantKind::Other,
        &span,
        &["T*".into(), "A".into()]
    )
    .starts_with("#chr1:1400-1401 "));
}

/// The key is overwritten by each allele in turn, so the best column does not decide.
#[test]
fn the_lines_are_sorted_by_the_last_likelihood() {
    let text = golden();
    let block = blocks(&text, "whole-reference")
        .into_iter()
        .find(|block| block.header.starts_with("#chr1:1400 "))
        .expect("the mixed site");
    let written: Vec<(String, Vec<f64>)> =
        block.lines.iter().map(|line| likelihoods(line)).collect();
    // Its first line has the WORSE first column and the BETTER last one.
    assert!(written[0].1[0] < written[1].1[0], "{written:?}");
    assert!(written[0].1[1] > written[1].1[1], "{written:?}");
    // The last column is what the order follows, descending.
    let last: Vec<f64> = written.iter().map(|(_, values)| values[1]).collect();
    for i in 1..last.len() {
        assert!(last[i - 1] >= last[i], "{last:?}");
    }
    // And the sort is stable: two reads with the same last column keep their order.
    let block = blocks(&text, "whole-reference")
        .into_iter()
        .next()
        .expect("the first site");
    let names: Vec<String> = block.lines.iter().map(|line| likelihoods(line).0).collect();
    let spans_all = names
        .iter()
        .position(|n| n == "r-spans-all")
        .expect("a line");
    let reverse = names.iter().position(|n| n == "r-reverse").expect("a line");
    assert!(spans_all < reverse);
    assert_eq!(
        likelihoods(&block.lines[spans_all]).1,
        likelihoods(&block.lines[reverse]).1
    );
}

/// The bases column comes out empty and the line is never added.
#[test]
fn a_read_that_does_not_span_the_variant_is_dropped() {
    let text = golden();
    let reads = reads(&text, "reads");
    let short = reads.iter().find(|r| r.name == "r-short").expect("a read");
    assert_eq!(short.span, Span::new("chr1", 1000, 1059));
    // It spans the SNP and is written.
    assert!(variant_bases(short, &Span::new("chr1", 1050, 1050)).is_some());
    assert!(blocks(&text, "whole-reference")[0]
        .lines
        .iter()
        .any(|line| line.starts_with("r-short ")));
    // It does not span the insertion and is absent from every later block.
    assert!(variant_bases(short, &Span::new("chr1", 1199, 1199)).is_none());
    for block in blocks(&text, "whole-reference").iter().skip(1) {
        assert!(
            !block.lines.iter().any(|line| line.starts_with("r-short ")),
            "{}",
            block.header
        );
    }
    // A line whose likelihoods are all negative infinity is dropped for the other reason.
    assert!(matrix_line(
        short,
        &Span::new("chr1", 1050, 1050),
        &[f64::NEG_INFINITY, f64::NEG_INFINITY]
    )
    .is_none());
}

/// The offset goes negative inside the element that follows, and the base is read early.
#[test]
fn a_variant_inside_a_deletion_is_not_dropped_but_comes_out_wrong() {
    let text = golden();
    let reads = reads(&text, "reads");
    let deleted = reads
        .iter()
        .find(|r| r.name == "r-deletion")
        .expect("a read");
    // Its deletion covers 1045 to 1054, so the variant at 1050 has no base of its own.
    assert_eq!(deleted.cigar, parse_cigar("45M10D65M"));
    // The walk returns an offset all the same, and it is ten less than the honest one.
    let offset = offset_on_read(&deleted.cigar, 50).expect("an offset");
    assert_eq!(offset, 40);
    assert_eq!(offset_on_read(&deleted.cigar, 55).expect("an offset"), 45);
    // Which is the base the golden wrote, taken ten positions early.
    let (bases, first) = variant_bases(deleted, &Span::new("chr1", 1050, 1050)).expect("bases");
    assert_eq!(first, 40);
    assert_eq!(bases, "T");
    assert!(blocks(&text, "whole-reference")[0]
        .lines
        .iter()
        .any(|line| line.starts_with("r-deletion 0 0 0 60 ") && line.contains(" T 40 sm1")));
    // Only running off the end of the read returns nothing at all.
    assert_eq!(offset_on_read(&parse_cigar("10M"), 20), None);
    assert_eq!(offset_on_read(&parse_cigar("10M"), 9), Some(9));
}

/// Counted from the other end of the read.
#[test]
fn the_offset_is_reversed_on_the_reverse_strand() {
    let text = golden();
    let reads = reads(&text, "reads");
    let forward = reads
        .iter()
        .find(|r| r.name == "r-spans-all")
        .expect("a read");
    let reverse = reads
        .iter()
        .find(|r| r.name == "r-reverse")
        .expect("a read");
    assert!(!forward.is_reverse);
    assert!(reverse.is_reverse);
    let variant = Span::new("chr1", 1050, 1050);
    assert_eq!(variant_bases(forward, &variant).expect("bases").1, 50);
    // 120 bases, so the offset from the other end is 120 - 50 - 1.
    assert_eq!(variant_bases(reverse, &variant).expect("bases").1, 69);
    assert_eq!(forward.bases.len(), 120);
    // The bases themselves are the same: only the offset is mirrored.
    assert_eq!(
        variant_bases(forward, &variant).expect("bases").0,
        variant_bases(reverse, &variant).expect("bases").0
    );
    // Over a multi-base variant the offset is the LAST base's, not the first's.
    let deletion = Span::new("chr1", 1299, 1302);
    let second = reads
        .iter()
        .find(|r| r.name == "r-second-reverse")
        .expect("a read");
    assert_eq!(variant_bases(second, &deletion).expect("bases").1, 167);
    assert_eq!(variant_bases(second, &deletion).expect("bases").0, "GTAC");
}

/// The default read filters take it out before the matrix is built.
#[test]
fn a_duplicate_read_never_reaches_the_matrix() {
    let text = golden();
    let reads = reads(&text, "reads");
    let duplicate = reads
        .iter()
        .find(|r| r.name == "r-duplicate")
        .expect("a read");
    assert!(duplicate.is_duplicate);
    // It spans the SNP as well as its twin does, and is in none of the four runs.
    assert!(variant_bases(duplicate, &Span::new("chr1", 1050, 1050)).is_some());
    for label in ["whole-reference", "one-interval", "mapping-quality-30"] {
        assert!(
            !section(&text, "out", label).contains("r-duplicate"),
            "{label}"
        );
    }
    // So the duplicate column is 0 on every line the golden wrote.
    for block in blocks(&text, "whole-reference") {
        for line in block.lines {
            assert_eq!(line.split(' ').nth(2), Some("0"), "{line}");
        }
    }
    // The column is still written from the flag: a duplicate that survived would read 1.
    let line =
        matrix_line(duplicate, &Span::new("chr1", 1050, 1050), &[-1.0, -2.0]).expect("a line");
    assert!(
        line.text.starts_with("r-duplicate 0 1 0 60 "),
        "{}",
        line.text
    );
}

/// A mapping-quality filter takes reads out, and an interval with no allele writes an empty file.
#[test]
fn the_runs_differ_only_in_which_reads_and_alleles_they_keep() {
    let text = golden();
    // The filter at 30 drops the read of quality 25 and the one of quality 20.
    let plain = section(&text, "out", "whole-reference");
    let filtered = section(&text, "out", "mapping-quality-30");
    assert!(plain.contains("r-second-low-mq"));
    assert!(!filtered.contains("r-second-low-mq"));
    assert!(plain.contains("r-clipped"));
    assert!(!filtered.contains("r-clipped"));
    // The interval limits the walk rather than the VCF: chr1:1000-1250 keeps the two alleles
    // that fall in it and drops the two that do not.
    assert_eq!(blocks(&text, "one-interval").len(), 2);
    assert_eq!(blocks(&text, "whole-reference").len(), 4);
    for i in 0..2 {
        assert_eq!(
            blocks(&text, "one-interval")[i].header,
            blocks(&text, "whole-reference")[i].header
        );
    }
    // An interval with no allele in it writes an empty file rather than none.
    assert_eq!(section(&text, "out", "empty-interval"), "");
    assert!(!text.contains("none\tempty-interval"));
}

/// The read that was soft-clipped is the one the trimming changed.
#[test]
fn the_reads_are_trimmed_before_the_line_is_built() {
    let text = golden();
    let reads = reads(&text, "reads");
    let clipped = reads
        .iter()
        .find(|r| r.name == "r-clipped")
        .expect("a read");
    assert_eq!(clipped.cigar, parse_cigar("10S110M"));
    // The cigar walk is unaffected: the offset into the read is still 60.
    assert_eq!(offset_on_read(&clipped.cigar, 50), Some(60));
    // And the golden wrote 60 as its UNCLIPPED offset, which the ten clipped bases would have
    // made 70 had the read reached the writer as it was read from the file.
    let line = blocks(&text, "whole-reference")[0]
        .lines
        .iter()
        .find(|line| line.starts_with("r-clipped "))
        .expect("its line")
        .clone();
    assert!(line.ends_with(" T 60 sm1"), "{line}");
    let untrimmed = Read {
        unclipped_start: clipped.span.start - 10,
        ..clipped.clone()
    };
    assert_eq!(
        variant_bases(&untrimmed, &Span::new("chr1", 1050, 1050))
            .expect("bases")
            .1,
        70,
        "what the untrimmed read would have written"
    );
    // The trimmed read, whose unclipped start is its own start, writes what the golden has.
    assert_eq!(
        variant_bases(clipped, &Span::new("chr1", 1050, 1050))
            .expect("bases")
            .1,
        60
    );
}
