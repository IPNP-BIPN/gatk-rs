//! `CheckReferenceCompatibility`, ported from
//! `org.broadinstitute.hellbender.tools.reference.CheckReferenceCompatibility` (GATK 4.6.2.0).
//!
//! A BAM or a VCF checked against several references, on the table
//! [`crate::compare_references`] already builds.
//!
//! # One property of the input decides the whole algorithm
//!
//! ```java
//! if(md5sPresent){ ... compareAgainstKeyReference ... } else { ... compareDictionaries ... }
//! ```
//!
//! With MD5s the sequences are compared by their bases; without, by name and length alone, and the
//! summaries say so outright. `dictionaryHasMD5s` needs EVERY sequence to have one, so a single
//! missing `M5` sends the whole run down the other path.
//!
//! # A VCF can never take the MD5 path
//!
//! `VCFContigHeaderLine.getSAMSequenceRecord` copies the ID, the length and `assembly`, and drops
//! `M5`. So a VCF header carrying an MD5 for every contig still produces a dictionary with none:
//! only a BAM or a CRAM reaches the MD5 branch. The golden holds five VCF runs to prove it,
//! including one whose `M5` is a lie and which is nevertheless called compatible.
//!
//! # The two paths disagree on what a subset is
//!
//! The MD5 path reads `ReferencePair`'s `SUBSET`, which means the INPUT is contained in the
//! reference. The other path reads `SequenceDictionaryUtils`' `SUPERSET`, which is the same
//! relation named from the reference's side. Both produce `COMPATIBLE_SUBSET`.

use crate::compare_references::{Pair, Status};

/// The verdict one reference earns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compatibility {
    Compatible,
    CompatibleSubset,
    NotCompatible,
}

impl Compatibility {
    pub fn name(&self) -> &'static str {
        match self {
            Compatibility::Compatible => "COMPATIBLE",
            Compatibility::CompatibleSubset => "COMPATIBLE_SUBSET",
            Compatibility::NotCompatible => "NOT_COMPATIBLE",
        }
    }
}

/// One row of the output table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub reference: String,
    pub compatibility: Compatibility,
    pub summary: String,
}

/// What the tool refuses before it reads anything.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputError {
    BothBamAndVcf,
    NoInput,
    ManyReadInputs,
}

impl InputError {
    pub fn java_class(&self) -> &str {
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    }

    pub fn message(&self) -> String {
        match self {
            InputError::BothBamAndVcf => {
                "Bad input: Both BAM and VCF specified. Tool analyzes one input at a time."
                    .to_string()
            }
            InputError::NoInput => "Bad input: No input provided.".to_string(),
            InputError::ManyReadInputs => {
                "Bad input: Tool analyzes one reads input at a time.".to_string()
            }
        }
    }
}

/// `initializeSequenceDictionaryForInput`'s two refusals, in the order it makes them.
pub fn check_input(has_reads: bool, read_inputs: usize, has_vcf: bool) -> Result<(), InputError> {
    if has_reads && has_vcf {
        return Err(InputError::BothBamAndVcf);
    }
    if has_reads && read_inputs > 1 {
        return Err(InputError::ManyReadInputs);
    }
    if !has_reads && !has_vcf {
        return Err(InputError::NoInput);
    }
    Ok(())
}

/// `dictionaryHasMD5s`: every sequence, or the run takes the other path.
pub fn md5s_present(md5s: &[Option<String>]) -> bool {
    !md5s.is_empty() && md5s.iter().all(|md5| md5.is_some())
}

/// `evaluateCompatibilityWithMD5Table`, whose verdict is read off the pair's flags alone.
///
/// `missing` is what `getMissingSequencesIfSubset` answers: the reference's sequence names the
/// input does not have, in the reference's own order.
pub fn evaluate_with_md5(pair: &Pair, missing: &[String]) -> Record {
    let compatibility = if pair.analysis.contains(&Status::ExactMatch) {
        Compatibility::Compatible
    } else if pair.analysis.contains(&Status::Subset) && pair.analysis.len() == 1 {
        Compatibility::CompatibleSubset
    } else {
        Compatibility::NotCompatible
    };
    let summary = match compatibility {
        Compatibility::Compatible => "The sequence dictionaries exactly match".to_string(),
        Compatibility::CompatibleSubset => format!(
            "The sequence dictionary in {} is a subset of the {} reference sequence dictionary. \
             Missing sequence(s): {}",
            pair.first,
            pair.second,
            rendered_list(missing)
        ),
        Compatibility::NotCompatible => format!(
            "Status: {}. Run CompareReferences tool for more information on reference differences.",
            rendered_set(pair)
        ),
    };
    Record {
        reference: pair.second.clone(),
        compatibility,
        summary,
    }
}

/// `SequenceDictionaryUtils.SequenceDictionaryCompatibility`, for the three answers this tool
/// reads. The name is the reference's view: `SUPERSET` means the reference has sequences the input
/// does not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DictionaryCompatibility {
    Identical,
    Superset,
    Other(&'static str),
}

/// `evaluateCompatibilityWithoutMD5`.
pub fn evaluate_without_md5(
    reference_name: &str,
    input_name: &str,
    status: DictionaryCompatibility,
    missing: &[String],
) -> Record {
    match status {
        DictionaryCompatibility::Identical => Record {
            reference: reference_name.to_string(),
            compatibility: Compatibility::Compatible,
            summary: "All sequence names and lengths match in the sequence dictionaries. Since \
                      the MD5s are lacking, we can't confirm there aren't mismatching bases in \
                      the references."
                .to_string(),
        },
        DictionaryCompatibility::Superset => Record {
            reference: reference_name.to_string(),
            compatibility: Compatibility::CompatibleSubset,
            summary: format!(
                "All sequence names and lengths present in the sequence dictionaries match, but \
                 {input_name} is a subset of {reference_name}. Missing sequence(s): {}. Since the \
                 MD5s are lacking, we can't confirm there aren't mismatching bases in the \
                 references.",
                rendered_list(missing)
            ),
        },
        DictionaryCompatibility::Other(name) => Record {
            reference: reference_name.to_string(),
            compatibility: Compatibility::NotCompatible,
            summary: format!(
                "Status: {name}. Run CompareReferences tool for more information on reference \
                 differences."
            ),
        },
    }
}

/// A Java `List.toString`.
fn rendered_list(values: &[String]) -> String {
    format!("[{}]", values.join(", "))
}

/// A Java `EnumSet.toString`, which is the enum's declaration order.
fn rendered_set(pair: &Pair) -> String {
    format!(
        "[{}]",
        pair.analysis
            .iter()
            .map(|status| status.name())
            .collect::<Vec<&str>>()
            .join(", ")
    )
}

/// `writeOutput`: the comment line, the header and one row per reference.
pub fn write_table(input_name: &str, records: &[Record]) -> String {
    let mut out = format!("#Current Reference: {input_name}\n");
    out.push_str("Reference\tCompatibility\tSummary\n");
    for record in records {
        out.push_str(&format!(
            "{}\t{}\t{}\n",
            record.reference,
            record.compatibility.name(),
            record.summary
        ));
    }
    out
}
