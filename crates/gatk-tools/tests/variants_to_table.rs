//! Conformance for `VariantsToTable` against GATK 4.6.2.0, compared as every line of every table.
//!
//! Golden from `tools/readfilter-conformance/VariantsToTableDump.java`.
//!
//! # What this suite is for
//!
//!  * **the wildcard sorts the values**, not the keys, and renders a list as Java prints one;
//!  * **`-ASF` leaves the space behind**, so a split column reads ` 5`;
//!  * **splitting spreads a value only when its length matches**;
//!  * **and a missing field is the string `NA`**.

use gatk_corpus as corpus;
use gatk_tools::variants_to_table::{extract_fields, header, Arguments, CountTypes, Record, Value};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/variants_to_table.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.splitn(2, '\t').collect())
        .collect()
}

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

/// The header's declarations: which INFO and FORMAT fields exist, in order, and which are `R`.
struct Declarations {
    info: Vec<String>,
    format: Vec<String>,
    per_allele: CountTypes,
}

fn declarations(whole: &str) -> Declarations {
    let mut info = Vec::new();
    let mut format = Vec::new();
    let mut per_allele = CountTypes::new();
    for line in whole.lines() {
        let kind = if line.starts_with("##INFO=<") {
            "INFO"
        } else if line.starts_with("##FORMAT=<") {
            "FORMAT"
        } else {
            continue;
        };
        let id = line
            .split_once("ID=")
            .and_then(|(_, rest)| rest.split(',').next())
            .expect("an ID")
            .to_string();
        per_allele.insert(id.clone(), line.contains("Number=R,"));
        if kind == "INFO" {
            info.push(id);
        } else {
            format.push(id);
        }
    }
    Declarations {
        info,
        format,
        per_allele,
    }
}

/// A value as the table sees it: a list where the file wrote commas AND the header says so.
fn value(text: &str, is_list: bool) -> Value {
    if is_list && text.contains(',') {
        Value::Many(text.split(',').map(|part| part.to_string()).collect())
    } else {
        Value::One(text.to_string())
    }
}

fn input(text: &str) -> (Vec<String>, Vec<Record>, Declarations) {
    let whole = unescape(rows(text, "input").first().expect("an input")[1]);
    let declared = declarations(&whole);
    let samples: Vec<String> = whole
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("a header")
        .split('\t')
        .skip(9)
        .map(|name| name.to_string())
        .collect();
    // Every INFO and FORMAT field of this file whose number is not 1 is a list.
    let is_list = |key: &str| -> bool {
        whole
            .lines()
            .any(|line| line.contains(&format!("ID={key},")) && !line.contains("Number=1,"))
    };
    let records = whole
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let field: Vec<&str> = line.split('\t').collect();
            let info = field[7]
                .split(';')
                .filter_map(|entry| entry.split_once('='))
                .map(|(key, text)| (key.to_string(), value(text, is_list(key))))
                .collect();
            let keys: Vec<&str> = field[8].split(':').collect();
            let genotypes = (0..samples.len())
                .map(|index| {
                    field[9 + index]
                        .split(':')
                        .enumerate()
                        .map(|(at, text)| {
                            let key = keys[at];
                            (key.to_string(), value(text, is_list(key)))
                        })
                        .collect()
                })
                .collect();
            Record {
                contig: field[0].to_string(),
                start: field[1].parse().expect("a position"),
                id: field[2].to_string(),
                reference: field[3].to_string(),
                alternates: field[4].split(',').map(|alt| alt.to_string()).collect(),
                qual: field[5].parse().ok(),
                filters: match field[6] {
                    "." | "PASS" => Vec::new(),
                    names => names.split(';').map(|name| name.to_string()).collect(),
                },
                info,
                genotypes,
            }
        })
        .collect();
    (samples, records, declared)
}

fn names(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

/// The arguments of each run, plus whether it asks for genotype columns at all.
fn setup(run: &str, declared: &Declarations) -> Arguments {
    let base = Arguments::default;
    match run {
        // No field at all takes the mandatory columns, then the header's INFO, then its FORMAT
        // with GT first.
        "no-fields" => {
            let mut fields = names(&["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]);
            fields.extend(declared.info.clone());
            let mut genotype_fields: Vec<String> = Vec::new();
            for field in &declared.format {
                if field == "GT" {
                    genotype_fields.insert(0, field.clone());
                } else {
                    genotype_fields.push(field.clone());
                }
            }
            Arguments {
                fields,
                genotype_fields,
                ..base()
            }
        }
        "mandatory" => Arguments {
            fields: names(&["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]),
            ..base()
        },
        "derived" => Arguments {
            fields: names(&[
                "TYPE",
                "EVENTLENGTH",
                "TRANSITION",
                "HET",
                "HOM-REF",
                "HOM-VAR",
                "NO-CALL",
                "NSAMPLES",
                "NCALLED",
                "MULTI-ALLELIC",
            ]),
            ..base()
        },
        "info" => Arguments {
            fields: names(&["POS", "DP", "AC"]),
            ..base()
        },
        "missing" => Arguments {
            fields: names(&["POS", "ZZ"]),
            ..base()
        },
        "wildcard" => Arguments {
            fields: names(&["POS", "AS_*"]),
            ..base()
        },
        "genotype" => Arguments {
            fields: names(&["POS"]),
            genotype_fields: names(&["GT", "GQ"]),
            ..base()
        },
        "genotype-filter" => Arguments {
            fields: names(&["POS"]),
            genotype_fields: names(&["FT"]),
            ..base()
        },
        "genotype-depths" => Arguments {
            fields: names(&["POS"]),
            genotype_fields: names(&["AD"]),
            ..base()
        },
        "show-filtered" => Arguments {
            fields: names(&["POS", "FILTER"]),
            show_filtered: true,
            ..base()
        },
        "split" => Arguments {
            fields: names(&["POS", "ALT", "AC", "AS_SB"]),
            split_multi_allelic: true,
            ..base()
        },
        "split-allele-specific" => Arguments {
            fields: names(&["POS", "ALT"]),
            allele_specific_fields: names(&["AS_SB", "AS_QD"]),
            split_multi_allelic: true,
            ..base()
        },
        "allele-specific-unsplit" => Arguments {
            fields: names(&["POS"]),
            allele_specific_fields: names(&["AS_SB"]),
            ..base()
        },
        "split-genotype-depths" => Arguments {
            fields: names(&["POS"]),
            allele_specific_genotype_fields: names(&["AD"]),
            split_multi_allelic: true,
            ..base()
        },
        other => panic!("no setup for {other}"),
    }
}

const RUNS: [&str; 14] = [
    "no-fields",
    "mandatory",
    "derived",
    "info",
    "missing",
    "wildcard",
    "genotype",
    "genotype-filter",
    "genotype-depths",
    "show-filtered",
    "split",
    "split-allele-specific",
    "allele-specific-unsplit",
    "split-genotype-depths",
];

fn table(text: &str, run: &str) -> Vec<String> {
    rows(text, "line")
        .into_iter()
        .filter(|row| row[0] == run)
        .map(|row| unescape(row[1]))
        .collect()
}

/// One whole run, header line included.
fn produced(
    records: &[Record],
    samples: &[String],
    declared: &Declarations,
    run: &str,
) -> Vec<String> {
    let arguments = setup(run, declared);
    // Asking for no genotype field empties the sample list, so the columns are site columns alone.
    let samples: Vec<String> = if arguments.genotype_fields.is_empty()
        && arguments.allele_specific_genotype_fields.is_empty()
    {
        Vec::new()
    } else {
        samples.to_vec()
    };

    let mut lines = vec![header(&samples, &arguments).join("\t")];
    for record in records {
        if !arguments.show_filtered && !record.filters.is_empty() {
            continue;
        }
        let rows = extract_fields(record, &samples, &arguments, &declared.per_allele)
            .expect("no missing field");
        for row in rows {
            lines.push(row.join("\t"));
        }
    }
    lines
}

#[test]
fn every_table_is_the_reference_s() {
    let text = golden();
    let (samples, records, declared) = input(&text);
    for run in RUNS {
        assert_eq!(
            produced(&records, &samples, &declared, run),
            table(&text, run),
            "table/{run}"
        );
    }
}

/// The wildcard's order is its values', and its rendering is Java's.
#[test]
fn the_wildcard_sorts_the_values_and_prints_a_list_as_java_does() {
    let text = golden();
    let (samples, records, declared) = input(&text);
    let ours = produced(&records, &samples, &declared, "wildcard");
    // AS_QD comes before AS_SB because `2` sorts before `[`, not because Q sorts before S.
    assert_eq!(ours[1], "100\t20.0,[10, 5]");
    assert_eq!(ours, table(&text, "wildcard"));
}

/// The space Java leaves in a list survives into the cell.
#[test]
fn an_allele_specific_field_keeps_the_space() {
    let text = golden();
    let (samples, records, declared) = input(&text);
    let split = produced(&records, &samples, &declared, "split-allele-specific");
    assert!(split[1].ends_with("\t 5\t20.0"), "{}", split[1]);
    let unsplit = produced(&records, &samples, &declared, "allele-specific-unsplit");
    assert_eq!(unsplit[1], "100\t10, 5");
}

/// A list of the wrong length is repeated whole rather than spread.
#[test]
fn splitting_spreads_only_a_list_of_the_right_length() {
    let text = golden();
    let (samples, records, declared) = input(&text);
    let ours = produced(&records, &samples, &declared, "split");
    // The multi-allelic record's two rows: AC is Number=A and splits, AS_SB is Number=R and does
    // not, so it lands whole in both.
    assert_eq!(ours[2], "200\tC\t1\t10,5,3");
    assert_eq!(ours[3], "200\tG\t2\t10,5,3");
}

/// The missing field, and the record whose FILTER column was a dot.
#[test]
fn a_missing_field_is_na_and_an_unfiltered_record_is_pass() {
    let text = golden();
    let (samples, records, declared) = input(&text);
    let missing = produced(&records, &samples, &declared, "missing");
    assert!(missing[1..].iter().all(|line| line.ends_with("\tNA")));

    let filtered = produced(&records, &samples, &declared, "show-filtered");
    assert_eq!(filtered[1], "100\tPASS");
    // The filtered record is only there because the flag is.
    assert!(filtered.iter().any(|line| line == "300\tLowQD"));
    let without = produced(&records, &samples, &declared, "mandatory");
    assert!(!without.iter().any(|line| line.starts_with("chr1\t300")));
}
