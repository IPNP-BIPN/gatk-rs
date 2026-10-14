//! Ported from `org.broadinstitute.hellbender.utils.samples` (`PedReader`, `Sample`, `SampleDB`,
//! `SampleDBBuilder`, `Trio`) (GATK 4.6.2.0): a pedigree file read into samples, and the trios and
//! founders the pedigree annotations ask it for.
//!
//! A PED line is six whitespace-separated fields: family, individual, father, mother, sex and
//! phenotype, a parent or family of `0` or `-9` meaning none. A parent named but given no line of
//! its own is added as a sample of the same family with the sex its role implies.
//!
//! What a reader would not guess:
//!
//!  * a line is split with Java's `split("\\s+")`, which keeps the empty field a leading space
//!    makes, so an indented line has seven fields and is refused;
//!  * the phenotype column decides the type of the whole file: if any line holds a value other
//!    than `-9`, `0`, `1` or `2`, every phenotype is quantitative, and only `-9` is then missing;
//!  * a sample seen twice is merged field by field, a missing value giving way to a real one and
//!    two different real values refused;
//!  * the trios come family by family in family-ID order, and within a family in sample order
//!    (ID first), each child with both parents in the database giving one trio; two children of
//!    the same parents give two trios, since `Trio` has no equality of its own;
//!  * a founder is a sample with neither parent in the database, in the order the samples were
//!    first added.

use crate::java_hash::compare_strings;

/// `Sex`, in declaration order, which `Sample.compareTo` sorts by.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Sex {
    Male,
    Female,
    Unknown,
}

/// `Affection`, in declaration order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Affection {
    Unknown,
    Affected,
    Unaffected,
    Other,
}

/// `Sample`.
#[derive(Debug, Clone, PartialEq)]
pub struct Sample {
    pub id: String,
    pub family_id: Option<String>,
    pub paternal_id: Option<String>,
    pub maternal_id: Option<String>,
    pub sex: Sex,
    pub affection: Affection,
}

/// What reading or merging a pedigree refuses, as the reference's exception and message.
#[derive(Debug, Clone, PartialEq)]
pub struct PedigreeError {
    pub class: &'static str,
    pub message: String,
}

impl PedigreeError {
    fn user(message: String) -> PedigreeError {
        PedigreeError {
            class: "org.broadinstitute.hellbender.exceptions.UserException",
            message,
        }
    }
}

/// `String.compareTo` on two optional strings, a null sorting first.
fn compare_optional(a: &Option<String>, b: &Option<String>) -> std::cmp::Ordering {
    match (a, b) {
        (None, None) => std::cmp::Ordering::Equal,
        (None, Some(_)) => std::cmp::Ordering::Less,
        (Some(_), None) => std::cmp::Ordering::Greater,
        (Some(a), Some(b)) => compare_strings(a, b),
    }
}

impl Sample {
    /// `Sample.compareTo`: the ID, then the family, the parents, the sex and the affection.
    fn compare(&self, other: &Sample) -> std::cmp::Ordering {
        compare_strings(&self.id, &other.id)
            .then_with(|| compare_optional(&self.family_id, &other.family_id))
            .then_with(|| compare_optional(&self.paternal_id, &other.paternal_id))
            .then_with(|| compare_optional(&self.maternal_id, &other.maternal_id))
            .then_with(|| self.sex.cmp(&other.sex))
            .then_with(|| self.affection.cmp(&other.affection))
    }

    /// `mergeSamples`: a missing (or the empty) value gives way to a real one, and two different
    /// real values are refused.
    fn merge(&self, newer: &Sample) -> Result<Sample, PedigreeError> {
        if self == newer {
            return Ok(newer.clone());
        }
        fn merge_text(
            id: &str,
            field: &str,
            a: &Option<String>,
            b: &Option<String>,
        ) -> Result<Option<String>, PedigreeError> {
            match (a, b) {
                (None, _) => Ok(b.clone()),
                (Some(_), None) => Ok(a.clone()),
                (Some(x), Some(y)) if x == y => Ok(a.clone()),
                (Some(x), Some(y)) => Err(PedigreeError::user(format!(
                    "Inconsistent values detected for {id} for field {field} value1 {x} value2 {y}"
                ))),
            }
        }
        fn merge_enum<T: Copy + PartialEq + std::fmt::Debug>(
            id: &str,
            field: &str,
            a: T,
            b: T,
            empty: T,
            show: fn(T) -> &'static str,
        ) -> Result<T, PedigreeError> {
            if a == empty {
                Ok(b)
            } else if b == empty || a == b {
                Ok(a)
            } else {
                Err(PedigreeError::user(format!(
                    "Inconsistent values detected for {id} for field {field} value1 {} value2 {}",
                    show(a),
                    show(b)
                )))
            }
        }
        Ok(Sample {
            id: self.id.clone(),
            family_id: merge_text(&self.id, "FamilyID", &self.family_id, &newer.family_id)?,
            paternal_id: merge_text(
                &self.id,
                "PaternalID",
                &self.paternal_id,
                &newer.paternal_id,
            )?,
            maternal_id: merge_text(
                &self.id,
                "MaterialID",
                &self.maternal_id,
                &newer.maternal_id,
            )?,
            sex: merge_enum(
                &self.id,
                "Gender",
                self.sex,
                newer.sex,
                Sex::Unknown,
                |s| match s {
                    Sex::Male => "MALE",
                    Sex::Female => "FEMALE",
                    Sex::Unknown => "UNKNOWN",
                },
            )?,
            affection: merge_enum(
                &self.id,
                "Affection",
                self.affection,
                newer.affection,
                Affection::Unknown,
                |a| match a {
                    Affection::Unknown => "UNKNOWN",
                    Affection::Affected => "AFFECTED",
                    Affection::Unaffected => "UNAFFECTED",
                    Affection::Other => "OTHER",
                },
            )?,
        })
    }
}

/// `Trio`: a child and both its parents, by ID.
#[derive(Debug, Clone, PartialEq)]
pub struct Trio {
    pub family_id: Option<String>,
    pub maternal_id: String,
    pub paternal_id: String,
    pub child_id: String,
}

/// `SampleDB`: the samples in the order their ID was first added, a later sighting merged into
/// the first.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SampleDb {
    samples: Vec<Sample>,
}

impl SampleDb {
    /// `addSample`.
    fn add_sample(&mut self, sample: Sample) -> Result<(), PedigreeError> {
        match self.samples.iter_mut().find(|s| s.id == sample.id) {
            Some(existing) => {
                *existing = existing.merge(&sample)?;
            }
            None => self.samples.push(sample),
        }
        Ok(())
    }

    pub fn sample(&self, id: &str) -> Option<&Sample> {
        self.samples.iter().find(|s| s.id == id)
    }

    /// `getParents`: the mother then the father, each only when it is in the database.
    fn parent_count(&self, sample: &Sample) -> usize {
        [&sample.maternal_id, &sample.paternal_id]
            .iter()
            .filter(|id| id.as_ref().is_some_and(|id| self.sample(id).is_some()))
            .count()
    }

    /// `getTrios()`: per family in family-ID order, per member in sample order, a trio for each
    /// child with both parents in the database.
    pub fn trios(&self) -> Vec<Trio> {
        let mut families: Vec<&str> = self
            .samples
            .iter()
            .filter_map(|s| s.family_id.as_deref())
            .collect();
        families.sort_by(|a, b| compare_strings(a, b));
        families.dedup();
        let mut trios = Vec::new();
        for family in families {
            let mut members: Vec<&Sample> = self
                .samples
                .iter()
                .filter(|s| s.family_id.as_deref() == Some(family))
                .collect();
            members.sort_by(|a, b| a.compare(b));
            members.dedup_by(|a, b| a.compare(b) == std::cmp::Ordering::Equal);
            for sample in members {
                if self.parent_count(sample) == 2 {
                    trios.push(Trio {
                        family_id: sample.family_id.clone(),
                        maternal_id: sample.maternal_id.clone().expect("a mother"),
                        paternal_id: sample.paternal_id.clone().expect("a father"),
                        child_id: sample.id.clone(),
                    });
                }
            }
        }
        trios
    }

    /// `getFounderIds()`: the samples with no parent in the database, in insertion order.
    pub fn founder_ids(&self) -> Vec<String> {
        self.samples
            .iter()
            .filter(|s| self.parent_count(s) == 0)
            .map(|s| s.id.clone())
            .collect()
    }
}

/// `String.split("\\s+")`: Java keeps a leading empty field and drops trailing ones.
fn java_split_whitespace(line: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = Vec::new();
    let mut start = 0;
    let mut in_whitespace = false;
    for (i, c) in line.char_indices() {
        if c.is_whitespace() || c == '\u{000B}' {
            if !in_whitespace {
                fields.push(&line[start..i]);
                in_whitespace = true;
            }
        } else if in_whitespace {
            start = i;
            in_whitespace = false;
        }
    }
    if !in_whitespace {
        fields.push(&line[start..]);
    }
    while fields.last() == Some(&"") {
        fields.pop();
    }
    fields
}

const MISSING_VALUE1: &str = "-9";
const MISSING_VALUE2: &str = "0";

fn maybe_missing(text: &str) -> Option<String> {
    (text != MISSING_VALUE1 && text != MISSING_VALUE2).then(|| text.to_string())
}

/// `SampleDBBuilder(STRICT).addSamplesFromPedigreeFiles([file]).getFinalSampleDB()`, over the
/// file's text, with all six fields expected. `source` is what the refusal of a short line names
/// where the reference names its reader.
pub fn read_pedigree(text: &str, source: &str) -> Result<SampleDb, PedigreeError> {
    let mut lines: Vec<Vec<&str>> = Vec::new();
    let mut line_no = 1;
    for line in text.lines() {
        // `XReadLines` skips the comment lines, and `splitLines` the blank ones.
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts = java_split_whitespace(line);
        if parts.len() != 6 {
            return Err(PedigreeError::user(format!(
                "{source}Bad PED line {line_no}: wrong number of fields"
            )));
        }
        lines.push(parts);
        line_no += 1;
    }
    let categorical = ["-9", "0", "1", "2"];
    let quantitative = lines.iter().any(|parts| !categorical.contains(&parts[5]));

    let mut db = SampleDb::default();
    let mut read: Vec<Sample> = Vec::new();
    for (index, parts) in lines.iter().enumerate() {
        let sex = match parts[4] {
            "1" => Sex::Male,
            "2" => Sex::Female,
            _ => Sex::Unknown,
        };
        let affection = if quantitative {
            if parts[5] == MISSING_VALUE1 {
                Affection::Unknown
            } else {
                Affection::Other
            }
        } else {
            match parts[5] {
                "-9" | "0" => Affection::Unknown,
                "1" => Affection::Unaffected,
                "2" => Affection::Affected,
                other => {
                    return Err(PedigreeError {
                        class: "org.broadinstitute.hellbender.exceptions.GATKException",
                        message: format!("Unexpected phenotype type {other} at line {}", index + 1),
                    })
                }
            }
        };
        let sample = Sample {
            id: parts[1].to_string(),
            family_id: maybe_missing(parts[0]),
            paternal_id: maybe_missing(parts[2]),
            maternal_id: maybe_missing(parts[3]),
            sex,
            affection,
        };
        read.push(sample.clone());
        db.add_sample(sample)?;
    }
    // `maybeAddImplicitSample`: a named parent with no line of its own.
    for sample in &read {
        for (parent, sex) in [
            (&sample.paternal_id, Sex::Male),
            (&sample.maternal_id, Sex::Female),
        ] {
            if let Some(id) = parent {
                if db.sample(id).is_none() {
                    db.add_sample(Sample {
                        id: id.clone(),
                        family_id: sample.family_id.clone(),
                        paternal_id: None,
                        maternal_id: None,
                        sex,
                        affection: Affection::Unknown,
                    })?;
                }
            }
        }
    }
    Ok(db)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_leading_space_is_a_field() {
        assert_eq!(java_split_whitespace(" a b"), vec!["", "a", "b"]);
        assert_eq!(java_split_whitespace("a  b\t"), vec!["a", "b"]);
    }

    #[test]
    fn a_trio_needs_both_parents_in_the_database() {
        let db = read_pedigree("f1 kid dad mom 1 1\nf1 dad 0 0 1 1\nf1 mom 0 0 2 1\n", "")
            .expect("a pedigree");
        assert_eq!(db.trios().len(), 1);
        assert_eq!(db.founder_ids(), vec!["dad", "mom"]);
        let implicit = read_pedigree("f1 kid dad mom 1 1\n", "").expect("a pedigree");
        assert_eq!(implicit.founder_ids(), vec!["dad", "mom"]);
        assert_eq!(implicit.trios().len(), 1);
    }
}
