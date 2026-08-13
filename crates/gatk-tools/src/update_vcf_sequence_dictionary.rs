//! `UpdateVCFSequenceDictionary`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.variantutils.UpdateVCFSequenceDictionary`
//! (GATK 4.6.2.0).
//!
//! The second member of the variant-transform archetype: it replaces a vcf's contig lines with a
//! dictionary taken from elsewhere, and refuses in five different Java classes.
//!
//! # The header order is the writer's, not the tool's
//!
//! ```java
//! new VCFHeader(inputHeader.getMetaDataInInputOrder(), inputHeader.getGenotypeSamples())
//! ```
//!
//! [`crate::remove_nearby_indels`] uses `getMetaDataInSortedOrder()` on the same line, and the
//! difference never reaches the file: the writer emits the lines in its own order, so an input
//! carrying `INFO` before `ALT` comes out `ALT` before `INFO` from either tool. The choice is
//! therefore not part of this port; what is ported is the dictionary and the checks.
//!
//! # The check that needs an argument reads the file's own header
//!
//! ```java
//! SAMSequenceDictionary oldDictionary = inputHeader == null ? null : inputHeader.getSequenceDictionary();
//! ```
//!
//! and not the engine's best dictionary, "since it might dig one up from an index", says the
//! comment. The refusal names the **feature input's name**, which is `drivingVariantFile`, and not
//! the path the user typed.
//!
//! # Validation is per record, after the file is already open
//!
//! `apply` checks each variant against the dictionary and throws, so a header and every earlier
//! record are on disk when the refusal comes. A contig the input's header had and the dictionary
//! lacks is dropped from the output header by `setSequenceDictionary` and only refused when the
//! traversal reaches a record on it.

use htsjdk_bam::header::SequenceRecord;
use htsjdk_vcf::variant::VariantContext;

/// `SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH`, which is what makes `LN:0` a refusal.
pub const UNKNOWN_SEQUENCE_LENGTH: i32 = 0;

/// The name a `FeatureInput` reports when it was given no logical name, which is what the
/// already-has-a-dictionary refusal quotes.
pub const DRIVING_VARIANTS_NAME: &str = "drivingVariantFile";

/// What the tool refuses, each with the Java class the reference throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UpdateDictionaryError {
    /// The input already had a dictionary and `--replace` was not given.
    AlreadyHasDictionary { input_name: String },
    /// A record on a contig the dictionary does not have.
    UnknownSequence { id: String, contig: String },
    /// A record whose end runs past the sequence's length.
    PastSequenceEnd {
        id: String,
        contig: String,
        end: i64,
        length: i32,
    },
    /// The dictionary source had no sequences at all.
    EmptyDictionarySource { source: String },
    /// A sequence with no length, which is `LN:0`. The message names the source and every
    /// sequence that lacked one.
    MissingContigLengths {
        source: String,
        missing: Vec<String>,
    },
    /// Both `--source-dictionary` and `--sequence-dictionary`.
    TwoDictionaries,
    /// Neither, and no reference either.
    NoDictionary,
}

impl UpdateDictionaryError {
    /// The message the reference carries, without the prefix its exception class adds.
    pub fn message(&self) -> String {
        match self {
            UpdateDictionaryError::AlreadyHasDictionary { input_name } => format!(
                "The input variant file {input_name} already contains a sequence dictionary. \
                 Use replace to force the dictionary to be replaced."
            ),
            UpdateDictionaryError::UnknownSequence { id, contig } => format!(
                "The input variant file contains a variant (ID: \"{id}\") with a reference to a \
                 sequence (\"{contig}\") that is not present in the provided dictionary"
            ),
            UpdateDictionaryError::PastSequenceEnd {
                id,
                contig,
                end,
                length,
            } => format!(
                "The input variant file contains a variant (ID: \"{id}\") with a reference to a \
                 sequence (\"{contig}\") that ends at a position ({end}) that exceeds the length \
                 of that sequence ({length}) in the provided dictionary"
            ),
            UpdateDictionaryError::EmptyDictionarySource { source } => format!(
                "The specified dictionary source has an empty or invalid sequence dictionary: \
                 {source}"
            ),
            // Two spaces before the newline, and the missing names listed at the end, which is
            // how the reference builds it.
            UpdateDictionaryError::MissingContigLengths { source, missing } => format!(
                "GATK SequenceDictionaryValidation requires all contigs in the dictionary to have \
                 lengths associated with them.  \nOne or more contigs in the dictionary from \
                 {source} are missing contig lengths.\nThe following contigs are missing lengths: \
                 {}",
                missing.join(", ")
            ),
            UpdateDictionaryError::TwoDictionaries => {
                "Only one of sequence-dictionary or source-dictionary may be specified on the \
                 command line"
                    .to_string()
            }
            UpdateDictionaryError::NoDictionary => {
                "A dictionary source file or reference file must be provided".to_string()
            }
        }
    }

    /// The Java class, which is a different one for four of the seven.
    pub fn java_class(&self) -> &'static str {
        match self {
            UpdateDictionaryError::TwoDictionaries => {
                "org.broadinstitute.barclay.argparser.CommandLineException"
            }
            UpdateDictionaryError::NoDictionary => {
                "org.broadinstitute.barclay.argparser.CommandLineException$MissingArgument"
            }
            UpdateDictionaryError::MissingContigLengths { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$SequenceDictionaryIsMissingContigLengths"
            }
            _ => "org.broadinstitute.barclay.argparser.CommandLineException$BadArgumentValue",
        }
    }
}

/// Which dictionary the tool works from, which is `getBestAvailableSequenceDictionary` overridden.
///
/// The override exists so that the NEW dictionary, and not the vcf's own, is what every caller
/// sees, "otherwise the wrong dictionary would be used when writing the index for the output vcf".
pub fn best_available_dictionary(
    source: Option<(&str, &[SequenceRecord])>,
    master: Option<&[SequenceRecord]>,
    reference: Option<&[SequenceRecord]>,
    validate_lengths: bool,
) -> Result<Vec<SequenceRecord>, UpdateDictionaryError> {
    let dictionary = match (source, master) {
        (None, Some(master)) => master.to_vec(),
        (None, None) => match reference {
            Some(reference) => reference.to_vec(),
            None => return Err(UpdateDictionaryError::NoDictionary),
        },
        (Some(_), Some(_)) => return Err(UpdateDictionaryError::TwoDictionaries),
        (Some((name, sequences)), None) => {
            if sequences.is_empty() {
                return Err(UpdateDictionaryError::EmptyDictionarySource {
                    source: name.to_string(),
                });
            }
            sequences.to_vec()
        }
    };

    // The length check runs on whichever dictionary won, and names the SOURCE argument even when
    // the dictionary came from somewhere else, which is why the source is threaded through.
    let missing: Vec<String> = dictionary
        .iter()
        .filter(|record| record.length == UNKNOWN_SEQUENCE_LENGTH)
        .map(|record| record.name.clone())
        .collect();
    if validate_lengths && !missing.is_empty() {
        return Err(UpdateDictionaryError::MissingContigLengths {
            source: source.map(|(name, _)| name.to_string()).unwrap_or_default(),
            missing,
        });
    }
    Ok(dictionary)
}

/// `onTraversalStart`'s check: an input that already has a dictionary needs `--replace`.
///
/// `input_dictionary` is what the FILE's header carried, not what the engine would hand back.
pub fn check_replace(
    input_dictionary: &[SequenceRecord],
    replace: bool,
) -> Result<(), UpdateDictionaryError> {
    if !replace && !input_dictionary.is_empty() {
        return Err(UpdateDictionaryError::AlreadyHasDictionary {
            input_name: DRIVING_VARIANTS_NAME.to_string(),
        });
    }
    Ok(())
}

/// `apply`: one record against the dictionary, refusing on the contig or on the end.
pub fn check_variant(
    dictionary: &[SequenceRecord],
    variant: &VariantContext,
) -> Result<(), UpdateDictionaryError> {
    let record = dictionary
        .iter()
        .find(|record| record.name == variant.contig);
    match record {
        None => Err(UpdateDictionaryError::UnknownSequence {
            id: variant.id.clone(),
            contig: variant.contig.clone(),
        }),
        // `vc.getEnd()`, which the INFO field END overrides, so a one-base record can end anywhere.
        Some(record) if variant.stop > i64::from(record.length) => {
            Err(UpdateDictionaryError::PastSequenceEnd {
                id: variant.id.clone(),
                contig: variant.contig.clone(),
                end: variant.stop,
                length: record.length,
            })
        }
        Some(_) => Ok(()),
    }
}

/// The whole traversal: the records that reach the output, and the refusal if one comes.
///
/// The two are returned together because the reference writes as it goes: when a record is
/// refused, everything before it is already on disk.
pub fn update_dictionary(
    dictionary: &[SequenceRecord],
    variants: &[VariantContext],
) -> (Vec<usize>, Option<UpdateDictionaryError>) {
    let mut written = Vec::new();
    for (index, variant) in variants.iter().enumerate() {
        if let Err(error) = check_variant(dictionary, variant) {
            return (written, Some(error));
        }
        written.push(index);
    }
    (written, None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::allele::Allele;

    fn dictionary() -> Vec<SequenceRecord> {
        vec![
            SequenceRecord::new("chr1", 250_000_000),
            SequenceRecord::new("chr2", 240_000_000),
        ]
    }

    fn variant(contig: &str, start: i64, id: &str, reference: &str) -> VariantContext {
        let alleles = vec![
            Allele::create(reference.as_bytes(), true).expect("a reference"),
            Allele::create(b"A", false).expect("an alternate"),
        ];
        let mut variant = VariantContext::new(contig, start, alleles);
        variant.id = id.to_string();
        variant
    }

    #[test]
    fn a_record_with_no_id_is_quoted_as_a_dot() {
        let error = check_variant(&dictionary(), &variant("chrUn", 1, ".", "A")).unwrap_err();
        assert_eq!(
            error.message(),
            "The input variant file contains a variant (ID: \".\") with a reference to a sequence \
             (\"chrUn\") that is not present in the provided dictionary"
        );
    }

    #[test]
    fn the_end_is_what_is_checked_and_an_end_attribute_moves_it() {
        let mut symbolic = variant("chr1", 100, "symbolic", "A");
        assert!(check_variant(&dictionary(), &symbolic).is_ok());
        // What INFO END does to the record, which is to move its end and nothing else.
        symbolic.stop = 250_000_001;
        let error = check_variant(&dictionary(), &symbolic).unwrap_err();
        assert!(
            error.message().contains("ends at a position (250000001)"),
            "{}",
            error.message()
        );
    }

    #[test]
    fn every_record_before_the_refusal_is_already_written() {
        let variants = vec![
            variant("chr1", 100, ".", "A"),
            variant("chrUn", 1, "bad", "A"),
            variant("chr1", 200, ".", "A"),
        ];
        let (written, error) = update_dictionary(&dictionary(), &variants);
        assert_eq!(written, vec![0]);
        assert!(matches!(
            error,
            Some(UpdateDictionaryError::UnknownSequence { .. })
        ));
    }

    #[test]
    fn the_four_argument_refusals_are_four_classes() {
        let empty: Vec<SequenceRecord> = Vec::new();
        let good = dictionary();

        let two =
            best_available_dictionary(Some(("d", &good)), Some(&good), None, true).unwrap_err();
        assert_eq!(
            two.java_class(),
            "org.broadinstitute.barclay.argparser.CommandLineException"
        );
        let none = best_available_dictionary(None, None, None, true).unwrap_err();
        assert_eq!(
            none.java_class(),
            "org.broadinstitute.barclay.argparser.CommandLineException$MissingArgument"
        );
        let empty_source =
            best_available_dictionary(Some(("d", &empty)), None, None, true).unwrap_err();
        assert!(matches!(
            empty_source,
            UpdateDictionaryError::EmptyDictionarySource { .. }
        ));
        let no_length = vec![
            SequenceRecord::new("chr1", 250_000_000),
            SequenceRecord::new("chr2", UNKNOWN_SEQUENCE_LENGTH),
        ];
        let missing =
            best_available_dictionary(Some(("d", &no_length)), None, None, true).unwrap_err();
        assert_eq!(
            missing.java_class(),
            "org.broadinstitute.hellbender.exceptions.UserException$SequenceDictionaryIsMissingContigLengths"
        );
        // And with validation off, the same dictionary is accepted.
        assert!(best_available_dictionary(Some(("d", &no_length)), None, None, false).is_ok());
    }

    #[test]
    fn the_replace_check_names_the_feature_input() {
        let error = check_replace(&dictionary(), false).unwrap_err();
        assert_eq!(
            error.message(),
            "The input variant file drivingVariantFile already contains a sequence dictionary. \
             Use replace to force the dictionary to be replaced."
        );
        assert!(check_replace(&dictionary(), true).is_ok());
        assert!(check_replace(&[], false).is_ok());
    }
}
