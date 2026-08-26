//! Conformance for `GeneExpressionEvaluation` against GATK 4.6.2.0, compared as the whole count
//! table of every run.
//!
//! Golden from `tools/readfilter-conformance/GeneExpressionEvaluationDump.java`.
//!
//! # What this suite is for
//!
//!  * **a single-end read being a refusal**, not a skipped read;
//!  * **PROPORTIONAL counting the bases a read does not cover**, where EQUAL does not;
//!  * **EQUAL multi-mapping moving the mapping quality filter** as a side effect;
//!  * **an unstranded feature emitting one row**, and everything over it counting as sense;
//!  * **unspliced mode swallowing the intron**;
//!  * **and a count of exactly one being written `1.0`**, because the branch that would write `1`
//!    is overwritten by the line after it.

use gatk_corpus as corpus;
use gatk_tools::gene_expression_evaluation::{
    count, in_good_pair, write_counts, BaseData, CountError, FeatureLabel, GroupingFeature,
    Interval, MultiMapMethod, MultiOverlapMethod, Read, ReadStrands, Settings, Strand,
};
use std::collections::BTreeMap;

const SAMPLE: &str = "sm1";
const BAM: &str = "<dir>/reads.bam";
const GFF: &str = "<dir>/features.gff3";

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/gene_expression_evaluation.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field<'a>(line: &'a str, name: &str) -> &'a str {
    line.split('\t')
        .find_map(|part| part.strip_prefix(&format!("{name}=")))
        .unwrap_or_else(|| panic!("the row carries {name}"))
}

/// The gff3 the run was given, verbatim.
fn gff(text: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix("gff\t"))
            .expect("the golden carries the gff3"),
    )
}

/// The gff3 attributes column, over the `key=value;key=value` shape this fixture uses. htsjdk's
/// codec also URL-decodes and splits on commas, neither of which is exercised here.
fn attributes(column: &str) -> BTreeMap<String, Vec<String>> {
    column
        .split(';')
        .filter(|part| !part.is_empty())
        .filter_map(|part| part.split_once('='))
        .map(|(key, value)| (key.to_string(), vec![value.to_string()]))
        .collect()
}

/// The grouping features and the descendants that give them their overlap intervals.
///
/// `shrinkBaseData` drops every attribute but the label's before the feature becomes a key, which
/// is what this reproduces: only the label survives.
fn features(
    text: &str,
    grouping: &str,
    overlap: &str,
    label: FeatureLabel,
) -> Vec<GroupingFeature> {
    let file = gff(text);
    let rows: Vec<Vec<String>> = file
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();

    let mut out = Vec::new();
    for row in rows.iter().filter(|row| row[2] == grouping) {
        let all = attributes(&row[8]);
        let id = all.get("ID").and_then(|values| values.first()).cloned();
        let mut kept = BTreeMap::new();
        if let Some(values) = all.get(label.key()) {
            kept.insert(label.key().to_string(), values.clone());
        }
        let overlaps = rows
            .iter()
            .filter(|child| child[2] == overlap)
            .filter(|child| {
                attributes(&child[8])
                    .get("Parent")
                    .and_then(|values| values.first())
                    .cloned()
                    == id
            })
            .map(|child| Interval {
                contig: child[0].clone(),
                start: child[3].parse().expect("a start"),
                end: child[4].parse().expect("an end"),
            })
            .collect();
        out.push(GroupingFeature {
            base: BaseData {
                contig: row[0].clone(),
                source: row[1].clone(),
                kind: row[2].clone(),
                start: row[3].parse().expect("a start"),
                end: row[4].parse().expect("an end"),
                strand: Strand::decode(&row[6]),
                attributes: kept,
            },
            overlaps,
        });
    }
    out
}

/// `SAMUtils.getAlignmentBlocks`, over the operators this fixture uses: M covers reference and
/// read, N and D cover reference alone.
fn blocks(contig: &str, start: i32, cigar: &str) -> Vec<Interval> {
    let mut out = Vec::new();
    let mut position = start;
    let mut length = String::new();
    for character in cigar.chars() {
        if character.is_ascii_digit() {
            length.push(character);
            continue;
        }
        let count: i32 = length.parse().expect("a cigar length");
        length.clear();
        match character {
            'M' | '=' | 'X' => {
                out.push(Interval {
                    contig: contig.to_string(),
                    start: position,
                    end: position + count - 1,
                });
                position += count;
            }
            'N' | 'D' => position += count,
            _ => {}
        }
    }
    out
}

fn reference_end(start: i32, cigar: &str) -> i32 {
    blocks("chr1", start, cigar)
        .last()
        .expect("an aligned block")
        .end
}

/// Every read of the measured BAM, in the order the walker saw them.
fn reads(text: &str) -> Vec<Read> {
    text.lines()
        .filter(|line| line.starts_with("read\t"))
        .map(|line| {
            let name = line.split('\t').nth(1).expect("a name").to_string();
            let start: i32 = field(line, "start").parse().expect("a start");
            let cigar = field(line, "cigar");
            let optional = |part: &str| match field(line, part) {
                "none" => None,
                value => Some(value.to_string()),
            };
            let paired = field(line, "paired") == "true";
            let mate_start = optional("mate-start").map(|value| value.parse().expect("a start"));
            let mate_cigar = optional("mate-cigar");
            Read {
                name,
                contig: field(line, "contig").to_string(),
                start,
                blocks: blocks("chr1", start, cigar),
                end: reference_end(start, cigar),
                reverse: field(line, "reverse") == "true",
                paired,
                proper_pair: field(line, "proper") == "true",
                first_of_pair: field(line, "first") == "true",
                mate_unmapped: false,
                mate_contig: paired.then(|| "chr1".to_string()),
                mate_start,
                mate_blocks: mate_cigar
                    .as_ref()
                    .map(|cigar| blocks("chr1", mate_start.expect("a mate start"), cigar)),
                mate_reverse: field(line, "mate-reverse") == "true",
                mate_quality: optional("mate-mq").map(|value| value.parse().expect("a quality")),
                mapping_quality: 60,
                hits: optional("nh").map(|value| value.parse().expect("a hit count")),
                fragment_length: field(line, "fragment").parse().expect("a fragment length"),
            }
        })
        .collect()
}

fn table(text: &str, label: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("counts\t{label}=")))
            .unwrap_or_else(|| panic!("the golden carries counts/{label}")),
    )
}

fn refusal(text: &str, label: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
            .unwrap_or_else(|| panic!("the golden carries error/{label}")),
    )
}

fn defaults() -> Settings {
    Settings {
        multi_overlap_method: MultiOverlapMethod::Proportional,
        multi_map_method: MultiMapMethod::Ignore,
        read_strands: ReadStrands::ForwardReverse,
        unspliced: false,
        feature_label: FeatureLabel::Name,
        minimum_mapping_quality: 10,
    }
}

/// label, settings.
fn runs() -> Vec<(&'static str, Settings)> {
    vec![
        ("default", defaults()),
        (
            "equal-overlap",
            Settings {
                multi_overlap_method: MultiOverlapMethod::Equal,
                ..defaults()
            },
        ),
        (
            "equal-multimap",
            Settings {
                multi_map_method: MultiMapMethod::Equal,
                ..defaults()
            },
        ),
        (
            "unspliced",
            Settings {
                unspliced: true,
                ..defaults()
            },
        ),
        (
            "forward-forward",
            Settings {
                read_strands: ReadStrands::ForwardForward,
                ..defaults()
            },
        ),
        (
            "reverse-forward",
            Settings {
                read_strands: ReadStrands::ReverseForward,
                ..defaults()
            },
        ),
        (
            "label-id",
            Settings {
                feature_label: FeatureLabel::Id,
                ..defaults()
            },
        ),
    ]
}

fn produced(text: &str, settings: &Settings) -> String {
    let features = features(text, "gene", "exon", settings.feature_label);
    let coverages = count(&features, &reads(text), settings).expect("a counted run");
    write_counts(
        &features,
        &coverages,
        SAMPLE,
        settings.feature_label,
        &[BAM.to_string()],
        GFF,
    )
}

#[test]
fn every_count_table_matches_the_golden() {
    let text = golden();
    let mut compared = 0;
    for (label, settings) in runs() {
        assert_eq!(produced(&text, &settings), table(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 7, "the golden's counted runs");
}

/// The uncovered part of a read is added to the denominator, so a gene a read half misses gets
/// less than a whole count. EQUAL never looks at a length and gives it one.
#[test]
fn proportional_charges_a_read_for_the_bases_it_does_not_cover() {
    let text = golden();
    let proportional = table(&text, "default");
    let equal = table(&text, "equal-overlap");
    assert!(proportional.contains("geneA\tchr1\t100\t299\t+\tsense\t3.3\n"));
    assert!(equal.contains("geneA\tchr1\t100\t299\t+\tsense\t3.5\n"));
    assert!(proportional.contains("geneB\tchr1\t260\t400\t+\tsense\t1.2\n"));
    assert!(equal.contains("geneB\tchr1\t260\t400\t+\tsense\t1.5\n"));
    // And the port produces both from the same reads.
    assert_eq!(produced(&text, &defaults()), proportional);
}

/// Choosing EQUAL for multi-mapped reads drops the mapping quality filter to zero as a side
/// effect, and admits the NH=3 read that IGNORE discards.
#[test]
fn equal_multimapping_admits_the_read_ignore_discards() {
    let text = golden();
    let ignore = defaults();
    let equal = Settings {
        multi_map_method: MultiMapMethod::Equal,
        ..defaults()
    };
    assert_eq!(ignore.effective_minimum_mapping_quality(), 10);
    assert_eq!(equal.effective_minimum_mapping_quality(), 0);
    assert!(table(&text, "equal-multimap").contains("geneA\tchr1\t100\t299\t+\tsense\t3.63\n"));
    // A third of a count more than IGNORE, which is the multi-mapped read's share.
    assert!(table(&text, "default").contains("geneA\tchr1\t100\t299\t+\tsense\t3.3\n"));
    assert!(reads(&text).iter().any(|read| read.hits == Some(3)));
}

/// An unstranded feature has no antisense row at all, and a read over it is sense whichever way it
/// points.
#[test]
fn an_unstranded_feature_emits_one_row() {
    let text = golden();
    let table = table(&text, "default");
    assert!(table.contains("geneD\tchr1\t2000\t2199\t.\tsense\t1.0\n"));
    assert!(!table.contains("geneD\tchr1\t2000\t2199\t.\tantisense"));
    // Every other feature has both rows.
    for gene in ["geneA", "geneB", "geneC"] {
        assert!(table.contains(&format!("{gene}\tchr1")));
        assert_eq!(
            table.lines().filter(|line| line.starts_with(gene)).count(),
            2,
            "{gene}"
        );
    }
}

/// Unspliced replaces the alignment blocks with one interval, so the intron the spliced read
/// straddles counts as covered and the proportions move.
#[test]
fn unspliced_swallows_the_intron() {
    let text = golden();
    assert!(table(&text, "default").contains("geneA\tchr1\t100\t299\t+\tsense\t3.3\n"));
    assert!(table(&text, "unspliced").contains("geneA\tchr1\t100\t299\t+\tsense\t2.62\n"));
    let spliced = reads(&text)
        .into_iter()
        .find(|read| read.name == "spliced")
        .expect("the spliced read");
    assert_eq!(
        spliced.blocks.len(),
        2,
        "two blocks, an intron between them"
    );
    assert_eq!(spliced.blocks[0].end + 1 + 100, spliced.blocks[1].start);
}

/// A count of exactly one is written `1.0`: the branch that writes `1` is overwritten by the line
/// that follows it.
#[test]
fn an_integral_count_is_written_with_its_decimal() {
    let text = golden();
    let table = table(&text, "default");
    assert!(table.contains("\t1.0\n"), "a whole count");
    assert!(table.contains("\t0.0\n"), "and a zero one");
    assert!(!table.contains("\t1\n"));
    assert!(!table.contains("\t0\n"));
}

/// The three refusals, each of which stops the run rather than skipping a read.
#[test]
fn the_three_refusals_match_the_golden() {
    let text = golden();

    // A read that is not paired at all, asked about its mate.
    let mut unpaired = reads(&text)[0].clone();
    unpaired.paired = false;
    let error = in_good_pair(&unpaired, 10, ReadStrands::ForwardReverse).expect_err("unpaired");
    assert_eq!(error, CountError::UnpairedRead);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "unpaired")
    );

    // A good pair whose MQ tag was never written.
    let mut no_quality = reads(&text)
        .into_iter()
        .find(|read| read.proper_pair && read.first_of_pair)
        .expect("a proper pair");
    no_quality.mate_quality = None;
    let error =
        in_good_pair(&no_quality, 10, ReadStrands::ForwardReverse).expect_err("no mate quality");
    assert_eq!(error, CountError::MissingMateQuality);
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-mate-quality")
    );

    // Grouping by a type that carries no label, which names the first feature it reaches.
    let exons = features(&text, "exon", "exon", FeatureLabel::Name);
    let error = count(&exons, &reads(&text), &defaults()).expect_err("no label");
    assert_eq!(
        error,
        CountError::NoLabel {
            key: "NAME".to_string(),
            contig: "chr1".to_string(),
            start: 100,
            end: 149
        }
    );
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "group-exon")
    );
}
