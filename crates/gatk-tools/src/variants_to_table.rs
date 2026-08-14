//! `VariantsToTable`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.variantutils.VariantsToTable` (GATK 4.6.2.0).
//!
//! A VCF in, a tab-separated table out, where almost no column is a copy of the file.
//!
//! # Two cells carry Java's rendering rather than the file's
//!
//! `-F AS_*` is a wildcard over the INFO keys, and what it collects into is a `TreeSet<String>` of
//! the **values**. So the order is the values' order and not the keys': a record whose `AS_QD` is
//! `20.0` and whose `AS_SB` is a list comes out `20.0,[10, 5]`, the AS_QD first because `2` sorts
//! before `[`, and the list rendered as Java prints a `List`, brackets and space included.
//!
//! `-ASF` then strips the brackets and splits on the comma without trimming, so the second entry of
//! `[10, 5]` is **` 5`**, a leading space and all, and that is what lands in the cell. The same
//! field read with `-F` gives `10, 5`. Both are reproduced: a port that trimmed would produce a
//! tidier table than the reference does.
//!
//! # What the defaults do
//!
//! No field at all means every field the header declares: the mandatory columns except INFO, then
//! every INFO line, then every FORMAT line with `GT` forced to the front. And asking for no
//! genotype field empties the sample list, so `-F CHROM` produces one column rather than one per
//! sample.
//!
//! # Splitting
//!
//! `--split-multi-allelic` produces one row per alternate, and `addFieldValue` spreads a value
//! across those rows **only when it is a list of exactly the right length**. A `Number=R` field has
//! one entry too many, so it lands whole in every row; a `Number=A` field is split. `-ASF` is what
//! makes an R-type field split properly, by dropping its first entry, which is the reference's.
//!
//! `-ASGF AD` is special-cased to `<reference depth>,<alternate depth>` per row, which no other
//! field is.

use std::collections::BTreeSet;

/// A value as the table sees it: either one string or a list, which is what decides splitting.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Value {
    One(String),
    /// A list, whose Java `toString` is `[a, b]` and which splits across rows when its length
    /// matches the number of alternates.
    Many(Vec<String>),
}

impl Value {
    /// `prettyPrintObject`, which joins a list with commas and no brackets.
    fn pretty(&self) -> String {
        match self {
            Value::One(value) => value.clone(),
            Value::Many(values) => values.join(","),
        }
    }

    /// `Object.toString`, which is what the wildcard and the allele-specific paths see.
    fn java_string(&self) -> String {
        match self {
            Value::One(value) => value.clone(),
            Value::Many(values) => format!("[{}]", values.join(", ")),
        }
    }
}

/// As much of a record as the table reads.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub contig: String,
    pub start: i32,
    pub id: String,
    pub reference: String,
    pub alternates: Vec<String>,
    /// The QUAL column, absent where the file said `.`.
    pub qual: Option<f64>,
    /// The FILTER column: empty for `.` or `PASS`.
    pub filters: Vec<String>,
    /// The INFO fields in the order the record wrote them.
    pub info: Vec<(String, Value)>,
    /// One per sample, keyed by FORMAT field.
    pub genotypes: Vec<Vec<(String, Value)>>,
}

/// The arguments that decide the shape of the table.
#[derive(Debug, Clone, Default)]
pub struct Arguments {
    pub fields: Vec<String>,
    pub genotype_fields: Vec<String>,
    pub allele_specific_fields: Vec<String>,
    pub allele_specific_genotype_fields: Vec<String>,
    pub split_multi_allelic: bool,
    pub show_filtered: bool,
    pub moltenize: bool,
    pub error_if_missing_data: bool,
}

/// What `--error-if-missing-data` raises.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MissingField {
    pub field: String,
}

/// Whether a field's declared count is `R`, which is the only thing `-ASF` treats differently.
pub type CountTypes = std::collections::HashMap<String, bool>;

const MISSING_DATA: &str = "NA";

/// `extractFields`: the rows one record produces.
pub fn extract_fields(
    record: &Record,
    samples: &[String],
    arguments: &Arguments,
    counts_are_per_allele: &CountTypes,
) -> Result<Vec<Vec<String>>, MissingField> {
    let rows = if arguments.split_multi_allelic {
        record.alternates.len()
    } else {
        1
    };
    let mut records: Vec<Vec<String>> = vec![Vec::new(); rows];

    for field in &arguments.fields {
        if arguments.split_multi_allelic && field == "ALT" {
            add(&mut records, &Value::Many(record.alternates.clone()));
        } else if let Some(value) = getter(record, field) {
            add(&mut records, &Value::One(value));
        } else if let Some(value) = attribute(record, field) {
            add(&mut records, value);
        } else if let Some(prefix) = field.strip_suffix('*') {
            // A `TreeSet<String>` of the VALUES, so the order is theirs and not the keys'.
            let matched: BTreeSet<String> = record
                .info
                .iter()
                .filter(|(key, _)| key.starts_with(prefix))
                .map(|(_, value)| value.java_string())
                .collect();
            let joined = if matched.is_empty() {
                MISSING_DATA.to_string()
            } else {
                matched.into_iter().collect::<Vec<_>>().join(",")
            };
            add(&mut records, &Value::One(joined));
        } else {
            missing(&mut records, field, arguments)?;
        }
    }

    for field in &arguments.allele_specific_fields {
        match attribute(record, field) {
            Some(value) => {
                // `getAttributeAsString` then `replace("[","").replace("]","")`, with no trim: the
                // space after the comma survives into the cell.
                let text = value.java_string();
                let stripped = text.replace(['[', ']'], "");
                if arguments.split_multi_allelic {
                    let parts: Vec<String> =
                        stripped.split(',').map(|part| part.to_string()).collect();
                    let per_allele = if counts_are_per_allele.get(field) == Some(&true) {
                        // An R-type field drops the reference's entry, which is what makes the
                        // rest line up with the alternates.
                        parts[1..].to_vec()
                    } else {
                        parts
                    };
                    add(&mut records, &Value::Many(per_allele));
                } else {
                    add(&mut records, &Value::One(stripped));
                }
            }
            None => missing(&mut records, field, arguments)?,
        }
    }

    if !arguments.genotype_fields.is_empty()
        || !arguments.allele_specific_genotype_fields.is_empty()
    {
        for (index, _) in samples.iter().enumerate() {
            let genotype = record.genotypes.get(index);
            for field in &arguments.genotype_fields {
                match genotype.and_then(|fields| present(fields, field)) {
                    // `getGenotypeString(true)` writes the BASES rather than the indices, so a
                    // `0/1` on an `A -> C` record is `A/C` and on an indel is `ACGT/A`.
                    Some(_) if field == "GT" => {
                        let call = genotype_string(record, genotype.expect("a genotype"));
                        add(&mut records, &Value::One(call));
                    }
                    Some(value) => add(&mut records, value),
                    None => missing(&mut records, field, arguments)?,
                }
            }
            for field in &arguments.allele_specific_genotype_fields {
                match genotype.and_then(|fields| find(fields, field)) {
                    Some(value) => {
                        if arguments.split_multi_allelic && field == "AD" {
                            // The one field with its own rule: the reference depth is repeated
                            // beside each alternate's.
                            let depths = match value {
                                Value::Many(values) => values.clone(),
                                Value::One(value) => vec![value.clone()],
                            };
                            let pairs: Vec<String> = depths[1..]
                                .iter()
                                .map(|depth| format!("{},{depth}", depths[0]))
                                .collect();
                            add(&mut records, &Value::Many(pairs));
                        } else if arguments.split_multi_allelic {
                            let text = value.java_string();
                            let parts: Vec<String> = text
                                .replace(['[', ']'], "")
                                .split(',')
                                .map(|part| part.to_string())
                                .collect();
                            let per_allele = if counts_are_per_allele.get(field) == Some(&true) {
                                parts[1..].to_vec()
                            } else {
                                parts
                            };
                            add(&mut records, &Value::Many(per_allele));
                        } else {
                            let text = value.java_string();
                            let rendered = if field == "AD" {
                                text.replace(['[', ']'], "").replace(' ', "")
                            } else {
                                text
                            };
                            add(&mut records, &Value::One(rendered));
                        }
                    }
                    None => missing(&mut records, field, arguments)?,
                }
            }
        }
    }

    Ok(records)
}

/// `addFieldValue`: a list of exactly the right length is spread, anything else is repeated.
fn add(records: &mut [Vec<String>], value: &Value) {
    if records.len() == 1 {
        records[0].push(value.pretty());
        return;
    }
    if let Value::Many(values) = value {
        if values.len() == records.len() {
            for (row, entry) in records.iter_mut().zip(values.iter()) {
                row.push(entry.clone());
            }
            return;
        }
    }
    let rendered = value.pretty();
    for row in records.iter_mut() {
        row.push(rendered.clone());
    }
}

fn missing(
    records: &mut [Vec<String>],
    field: &str,
    arguments: &Arguments,
) -> Result<(), MissingField> {
    if arguments.error_if_missing_data {
        return Err(MissingField {
            field: field.to_string(),
        });
    }
    add(records, &Value::One(MISSING_DATA.to_string()));
    Ok(())
}

fn attribute<'a>(record: &'a Record, field: &str) -> Option<&'a Value> {
    find(&record.info, field)
}

fn find<'a>(fields: &'a [(String, Value)], key: &str) -> Option<&'a Value> {
    fields
        .iter()
        .find(|(name, _)| name == key)
        .map(|(_, value)| value)
}

/// The same lookup, minus the values the reference reads as absent.
///
/// `getAnyAttribute("FT")` returns null for a genotype that is not filtered, and `PASS` is not
/// filtered, so a `PASS` in the file becomes an `NA` in the table. The reference's own comment
/// calls this out as a bug in `hasAnyAttribute`; the branch that follows it is what produces the
/// `NA`, and this is that branch.
fn present<'a>(fields: &'a [(String, Value)], key: &str) -> Option<&'a Value> {
    let value = find(fields, key)?;
    if key == "FT" && matches!(value, Value::One(text) if text == "PASS" || text == ".") {
        return None;
    }
    Some(value)
}

/// `Genotype.getGenotypeString(true)`: the alleles, spelled out, separated by `/`.
fn genotype_string(record: &Record, fields: &[(String, Value)]) -> String {
    calls_of(fields)
        .into_iter()
        .map(|allele| match allele {
            Some(0) => record.reference.clone(),
            Some(index) => record
                .alternates
                .get(index - 1)
                .cloned()
                .unwrap_or_else(|| ".".to_string()),
            None => ".".to_string(),
        })
        .collect::<Vec<_>>()
        .join("/")
}

/// The columns computed from the record rather than read off it.
fn getter(record: &Record, field: &str) -> Option<String> {
    Some(match field {
        "CHROM" => record.contig.clone(),
        "POS" => record.start.to_string(),
        "REF" => record.reference.clone(),
        "ALT" => {
            if record.alternates.is_empty() {
                ".".to_string()
            } else {
                record.alternates.join(",")
            }
        }
        "ID" => record.id.clone(),
        // `Double.toString` of the phred score, so an integer quality is written with a `.0`.
        "QUAL" => gatk_engine::tsv_table::java_double_to_string(record.qual.unwrap_or(0.0)),
        // `PASS` for a record whose column said `.`, which is not what the file holds.
        "FILTER" => {
            if record.filters.is_empty() {
                "PASS".to_string()
            } else {
                record.filters.join(",")
            }
        }
        "TYPE" => variant_type(record),
        "EVENTLENGTH" => {
            let mut longest = 0i32;
            for alternate in &record.alternates {
                let length = alternate.len() as i32 - record.reference.len() as i32;
                if length.abs() > longest.abs() {
                    longest = length;
                }
            }
            longest.to_string()
        }
        "TRANSITION" => {
            if variant_type(record) == "SNP" && record.alternates.len() == 1 {
                if is_transition(&record.reference, &record.alternates[0]) {
                    "1".to_string()
                } else {
                    "0".to_string()
                }
            } else {
                "-1".to_string()
            }
        }
        "HET" => count(record, |calls| {
            calls.len() > 1 && calls.iter().any(|a| a != &calls[0])
        }),
        "HOM-REF" => count(record, |calls| calls.iter().all(|a| a == &Some(0))),
        "HOM-VAR" => count(record, |calls| {
            calls.iter().all(|a| a.is_some_and(|index| index > 0))
                && calls.iter().all(|a| a == &calls[0])
        }),
        "NO-CALL" => count(record, |calls| calls.iter().all(Option::is_none)),
        "NSAMPLES" => record.genotypes.len().to_string(),
        "NCALLED" => count(record, |calls| calls.iter().any(Option::is_some)),
        "MULTI-ALLELIC" => (record.alternates.len() > 1).to_string(),
        _ => return None,
    })
}

/// `determineType` as far as the TYPE column needs it.
fn variant_type(record: &Record) -> String {
    if record.alternates.is_empty() {
        return "NO_VARIATION".to_string();
    }
    let mut kind: Option<&str> = None;
    for alternate in &record.alternates {
        let this = if alternate.len() == record.reference.len() {
            if record.reference.len() == 1 {
                "SNP"
            } else {
                "MNP"
            }
        } else {
            "INDEL"
        };
        match kind {
            None => kind = Some(this),
            Some(seen) if seen != this => return "MIXED".to_string(),
            Some(_) => {}
        }
    }
    kind.unwrap_or("NO_VARIATION").to_string()
}

fn is_transition(reference: &str, alternate: &str) -> bool {
    matches!(
        (reference, alternate),
        ("A", "G") | ("G", "A") | ("C", "T") | ("T", "C")
    )
}

/// The genotype counts, over the calls this port models as allele indices.
fn count(record: &Record, predicate: impl Fn(&[Option<usize>]) -> bool) -> String {
    record
        .genotypes
        .iter()
        .filter(|fields| {
            let calls = calls_of(fields);
            !calls.is_empty() && predicate(&calls)
        })
        .count()
        .to_string()
}

fn calls_of(fields: &[(String, Value)]) -> Vec<Option<usize>> {
    match find(fields, "GT") {
        Some(Value::One(call)) => call
            .split(['/', '|'])
            .map(|allele| allele.parse::<usize>().ok())
            .collect(),
        _ => Vec::new(),
    }
}

/// The header line, and the sample-qualified genotype columns under it.
pub fn header(samples: &[String], arguments: &Arguments) -> Vec<String> {
    if arguments.moltenize {
        return vec![
            "RecordID".to_string(),
            "Sample".to_string(),
            "Variable".to_string(),
            "Value".to_string(),
        ];
    }
    let mut fields = arguments.fields.clone();
    fields.extend(arguments.allele_specific_fields.clone());
    for sample in samples {
        // A space in a sample name is legal and wrecks an R data frame, so it becomes an
        // underscore here and nowhere else.
        let name = sample.replace(' ', "_");
        for field in arguments
            .genotype_fields
            .iter()
            .chain(arguments.allele_specific_genotype_fields.iter())
        {
            fields.push(format!("{name}.{field}"));
        }
    }
    fields
}
