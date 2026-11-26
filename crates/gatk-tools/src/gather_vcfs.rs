//! `GatherVcfsCloud`, ported from `org.broadinstitute.hellbender.tools.GatherVcfsCloud`
//! (GATK 4.7.0.0; `GatherVcfsCloud.java` is unchanged since 4.6.2.0).
//!
//! Concatenating vcfs that are already in order, which sounds like `cat` and is not: two gathering
//! modes with different code paths, and a validation that runs before anything is written.
//!
//! # Two order checks, not one
//!
//! `assertSameSamplesAndValidOrdering` compares the **first** record of each file and throws an
//! `IllegalArgumentException`. Then `gatherConventionally` compares the **last** record it wrote
//! against the next file's first, and throws an `IllegalStateException` naming both positions. A
//! pair of files that overlap in the middle passes the first and is caught by the second, which is
//! why this port returns both refusals as distinct variants rather than one.
//!
//! # Two flags, two broken files
//!
//! `--disable-contig-ordering-check` does not disable the check, it **weakens** it to comparing
//! positions within a contig, so a chr2 shard gathered before a chr1 shard is accepted and the
//! output holds `chr2:100` followed by `chr1:100`: a vcf whose records are not in dictionary order.
//!
//! `--ignore-safety-checks` skips the sample-list comparison and writes anyway, so a record
//! belonging to `s1` comes out under a header declaring only `s0`. The genotype is not dropped, it
//! is relabelled. Both are reproduced here, because a port that refused where the reference writes
//! would not be the reference.
//!
//! # What is modelled, and what is not
//!
//! A shard here is its header's contig names and samples and its records' `(contig, position)`;
//! the bytes are the caller's business. `--create-output-variant-index` is not an argument: the
//! reference's default, `createIndex = true`, is assumed throughout. The validation's
//! `dict.assertSameDictionary(header.getSequenceDictionary())` is not ported.
//!
//! # Reading this file without reading Rust
//!
//! - `&[Shard]` is a borrowed, read-only view of a list; the caller keeps it;
//! - `Option<T>` is a value or nothing; `Result<T, GatherError>` is a value or one of the
//!   refusals, and `?` hands a refusal straight back to the caller;
//! - `if let Some(x) = y` runs its block only when `y` holds a value, and
//!   `if let (Some(a), Some(b)) = (..)` only when both do;
//! - `match` must cover every variant, and `_` is "every other one";
//! - `|x| ..` is a closure (Java's lambda); `.iter().filter(..).cloned().collect()` is a stream
//!   that copies each kept element into a new list;
//! - `..Arguments::default()` (in the tests) fills the fields not named from the defaults.

/// One input file, as far as gathering reads it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Shard {
    /// The name the messages use, which is a URI in the reference.
    pub name: String,
    /// The sequence dictionary from the header: empty where the file has no `##contig` line.
    pub dictionary: Vec<String>,
    /// `header.getGenotypeSamples()`, in header order.
    pub samples: Vec<String>,
    /// `(contig, position)` for each record, in file order.
    pub records: Vec<(String, i32)>,
}

/// How to gather, which `AUTOMATIC` resolves from the file names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GatherType {
    /// Copy the inputs' compressed blocks; in the reference the writer's overlap check does not
    /// run.
    Block,
    /// Decode and re-encode each record, checking each file's first against the last written.
    Conventional,
    /// `BLOCK` when every input and the output are block compressed, else `CONVENTIONAL`.
    Automatic,
}

/// The tool's arguments that change what is gathered or refused.
pub struct Arguments {
    pub gather_type: GatherType,
    pub ignore_safety_checks: bool,
    pub disable_contig_ordering_check: bool,
    /// Whether the output path is block compressed, which is half of what `AUTOMATIC` asks.
    pub output_is_block_compressed: bool,
    /// Whether every input is block compressed, which is the other half.
    pub inputs_are_block_compressed: bool,
}

impl Default for Arguments {
    fn default() -> Arguments {
        Arguments {
            gather_type: GatherType::Automatic,
            ignore_safety_checks: false,
            disable_contig_ordering_check: false,
            output_is_block_compressed: false,
            inputs_are_block_compressed: false,
        }
    }
}

/// What the gather refuses, each with the class the reference throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GatherError {
    /// The first file has no `##contig` line, so the output cannot be indexed: `UserException`.
    ///
    /// The reference throws this BEFORE the validation, and the port checks it after (see
    /// [`gather`]).
    NoDictionary,
    /// Two files whose sample lists differ, named in both directions.
    ///
    /// Each list is the reference's `TreeSet` after `removeAll`: sorted, and printed by the
    /// set's `toString` as `[a, b]`.
    DifferentSamples {
        first_only: Vec<String>,
        other: String,
        other_only: Vec<String>,
    },
    /// The validation's check, on first records.
    FirstRecordNotAfter { file: String, previous: String },
    /// The writer's check, on the last record written against the next file's first.
    RecordsOverlap {
        file: String,
        at: (String, i32),
        previous: String,
        last: (String, i32),
    },
    /// The same check under --disable-contig-ordering-check, which compares positions only and
    /// says so in different words.
    PositionsOverlap {
        file: String,
        at: i32,
        previous: String,
        last: i32,
    },
    /// `BLOCK` asked for where something is not bgzipped.
    BlockCopyImpossible,
}

impl GatherError {
    /// The exception's message, as the reference builds it from the file URIs and positions.
    pub fn message(&self) -> String {
        match self {
            GatherError::NoDictionary => {
                "In order to index the resulting VCF, the input VCFs must contain ##contig lines."
                    .to_string()
            }
            GatherError::DifferentSamples {
                first_only,
                other,
                other_only,
            } => format!(
                "VCFs do not have identical sample lists. Samples unique to first file: [{}]. \
                 Samples unique to {other}: [{}].",
                first_only.join(", "),
                other_only.join(", ")
            ),
            GatherError::FirstRecordNotAfter { file, previous } => format!(
                "First record in file {file} is not after first record in previous file {previous}"
            ),
            GatherError::RecordsOverlap {
                file,
                at,
                previous,
                last,
            } => format!(
                "First variant in file {file} is at {}:{} but last variant in earlier file \
                 {previous} is at {}:{}",
                at.0, at.1, last.0, last.1
            ),
            GatherError::PositionsOverlap {
                file,
                at,
                previous,
                last,
            } => format!(
                "First variant in file {file} is at start position {at} but last variant in \
                 earlier file {previous} is at start position {last}"
            ),
            GatherError::BlockCopyImpossible => {
                "Requested block copy but some files are not bgzipped, all inputs and the output \
                 must be bgzipped to block copy"
                    .to_string()
            }
        }
    }

    /// The fully qualified class of the exception the reference throws for this refusal.
    pub fn java_class(&self) -> &'static str {
        match self {
            GatherError::NoDictionary => "org.broadinstitute.hellbender.exceptions.UserException",
            GatherError::BlockCopyImpossible => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
            GatherError::RecordsOverlap { .. } | GatherError::PositionsOverlap { .. } => {
                "java.lang.IllegalStateException"
            }
            _ => "java.lang.IllegalArgumentException",
        }
    }
}

/// The gather: the records written, in order, or the refusal.
///
/// The records are returned as `(shard, index)` rather than copied, since block gathering copies
/// bytes and conventional gathering re-encodes: what both agree on is which record goes where.
///
/// `GatherVcfsCloud.doWork` checks in this order: the missing dictionary
/// (`if (createIndex && sequenceDictionary == null)`), then the validation unless
/// `--ignore-safety-checks`, then the block-copy refusal. This function checks the block copy
/// first and the dictionary last, so an input that fails two of them gets a different refusal
/// from the reference's. The writer's check runs whatever the
/// gather type, where the reference runs it only in `gatherConventionally`; the caller refuses
/// a block gather on its own afterwards.
pub fn gather(shards: &[Shard], arguments: &Arguments) -> Result<Vec<(usize, usize)>, GatherError> {
    if arguments.gather_type == GatherType::Block
        && !(arguments.inputs_are_block_compressed && arguments.output_is_block_compressed)
    {
        return Err(GatherError::BlockCopyImpossible);
    }

    if !arguments.ignore_safety_checks {
        validate(shards, arguments)?;
    }

    // The missing-dictionary refusal, a `UserException` that the reference raises before the
    // validation.
    // `is_some_and` is false for an empty list of shards.
    if shards
        .first()
        .is_some_and(|first| first.dictionary.is_empty())
    {
        return Err(GatherError::NoDictionary);
    }

    let mut written: Vec<(usize, usize)> = Vec::new();
    // `lastContext`: the last record written so far, `None` until a shard with records is seen.
    let mut last: Option<(String, i32)> = None;
    // `.enumerate()` pairs each shard with its index, used to name the previous file.
    for (index, shard) in shards.iter().enumerate() {
        if let (Some(last), Some(first)) = (&last, shard.records.first()) {
            // The writer's own check, on the LAST record written, which honours the same flag as
            // the validation and says something different when it is set.
            // Under the flag the reference tests only `vc.getStart() <= lastContext.getStart()` on
            // one contig; without it, `comparator.compare(vc, lastContext) <= 0`. The previous
            // file named is `lastFile`, which the reference sets after every file, records or not.
            if arguments.disable_contig_ordering_check {
                if last.0 == first.0 && first.1 <= last.1 {
                    return Err(GatherError::PositionsOverlap {
                        file: shard.name.clone(),
                        at: first.1,
                        previous: shards[index - 1].name.clone(),
                        last: last.1,
                    });
                }
            } else if !after(first, last, &dictionary_of(shards)) {
                return Err(GatherError::RecordsOverlap {
                    file: shard.name.clone(),
                    at: first.clone(),
                    previous: shards[index - 1].name.clone(),
                    last: last.clone(),
                });
            }
        }
        for record in 0..shard.records.len() {
            written.push((index, record));
        }
        if let Some(record) = shard.records.last() {
            last = Some(record.clone());
        }
    }
    Ok(written)
}

/// `assertSameSamplesAndValidOrdering`: the samples, then the first records.
///
/// Each file's samples must equal the first file's, list for list (order included, as
/// `samples.equals(theseSamples)` compares). Each file's first record must be strictly after the
/// previous non-empty file's first record. The reference's refusal of a first file without a
/// dictionary (`UserException.BadInput`) is not here, since with `createIndex` true the
/// missing-dictionary refusal has already been thrown.
fn validate(shards: &[Shard], arguments: &Arguments) -> Result<(), GatherError> {
    // No shards, nothing to validate (the tool requires at least one input).
    let first = match shards.first() {
        None => return Ok(()),
        Some(first) => first,
    };
    let dictionary = first.dictionary.clone();

    let mut previous: Option<(&Shard, (String, i32))> = None;
    for shard in shards {
        if shard.samples != first.samples {
            // The names in `left` and not in `right`, sorted, as `TreeSet` then `removeAll`.
            let unique = |left: &[String], right: &[String]| {
                let mut names: Vec<String> = left
                    .iter()
                    .filter(|name| !right.contains(name))
                    .cloned()
                    .collect();
                names.sort();
                names
            };
            return Err(GatherError::DifferentSamples {
                first_only: unique(&first.samples, &shard.samples),
                other: shard.name.clone(),
                other_only: unique(&shard.samples, &first.samples),
            });
        }

        if let Some(record) = shard.records.first() {
            if let Some((previous_shard, previous_record)) = &previous {
                let out_of_order = if arguments.disable_contig_ordering_check {
                    // The weakened check: positions only, and only within one contig. A different
                    // contig is not compared at all, which is what lets an unordered file out.
                    previous_record.0 == record.0 && previous_record.1 >= record.1
                } else {
                    !after(record, previous_record, &dictionary)
                };
                if out_of_order {
                    return Err(GatherError::FirstRecordNotAfter {
                        file: shard.name.clone(),
                        previous: previous_shard.name.clone(),
                    });
                }
            }
            previous = Some((shard, record.clone()));
        }
    }
    Ok(())
}

/// The first shard's contig names, which `gatherConventionally` builds its comparator from
/// (`new VariantContextComparator(firstHeader.getContigLines())`); empty for no shards.
fn dictionary_of(shards: &[Shard]) -> Vec<String> {
    shards
        .first()
        .map(|shard| shard.dictionary.clone())
        .unwrap_or_default()
}

/// `VariantContextComparator`: is `left` strictly after `right`? Contig order is the DICTIONARY'S,
/// then position.
///
/// ```java
/// final int contigCompare = this.contigIndexLookup.get(firstVariantContext.getContig())
///         - this.contigIndexLookup.get(secondVariantContext.getContig());
/// ```
///
/// and on a tie, the start positions are subtracted.
///
/// A contig the dictionary lacks gets index -1 here, so two such contigs compare by position.
/// The reference has no such case: `contigIndexLookup` is a map, `get` returns null for a
/// missing contig, and unboxing it throws a `NullPointerException` (htsjdk 5.0.0's own comment:
/// "Will throw NullPointerException -- happily").
fn after(left: &(String, i32), right: &(String, i32), dictionary: &[String]) -> bool {
    // `position` is `indexOf`: `Some(i)` for the first match, `None` for none; `map_or` turns
    // `None` into -1 and widens the index to `i64`.
    let index = |contig: &String| {
        dictionary
            .iter()
            .position(|name| name == contig)
            .map_or(-1i64, |at| at as i64)
    };
    let (first, second) = (index(&left.0), index(&right.0));
    if first != second {
        first > second
    } else {
        left.1 > right.1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shard(name: &str, samples: &[&str], records: &[(&str, i32)]) -> Shard {
        Shard {
            name: name.to_string(),
            dictionary: vec!["chr1".to_string(), "chr2".to_string()],
            samples: samples.iter().map(|name| name.to_string()).collect(),
            records: records
                .iter()
                .map(|(contig, at)| (contig.to_string(), *at))
                .collect(),
        }
    }

    #[test]
    fn the_two_order_checks_are_two_different_refusals() {
        // Out of order at the first record: the validation catches it.
        let shards = vec![
            shard("second", &["s0"], &[("chr1", 300)]),
            shard("first", &["s0"], &[("chr1", 100), ("chr1", 200)]),
        ];
        let error = gather(&shards, &Arguments::default()).unwrap_err();
        assert_eq!(error.java_class(), "java.lang.IllegalArgumentException");

        // Overlapping in the middle: the first records are in order, the writer catches it.
        let shards = vec![
            shard("first", &["s0"], &[("chr1", 100), ("chr1", 200)]),
            shard("overlapping", &["s0"], &[("chr1", 150), ("chr1", 250)]),
        ];
        let error = gather(&shards, &Arguments::default()).unwrap_err();
        assert_eq!(error.java_class(), "java.lang.IllegalStateException");
        assert!(error.message().contains("is at chr1:150 but last variant"));
    }

    #[test]
    fn the_weakened_check_writes_an_unordered_file() {
        let shards = vec![
            shard("third", &["s0"], &[("chr2", 100)]),
            shard("first", &["s0"], &[("chr1", 100), ("chr1", 200)]),
        ];
        assert!(gather(&shards, &Arguments::default()).is_err());

        let written = gather(
            &shards,
            &Arguments {
                disable_contig_ordering_check: true,
                ..Arguments::default()
            },
        )
        .expect("accepted");
        // chr2 before chr1, which is not the dictionary's order.
        assert_eq!(written, vec![(0, 0), (1, 0), (1, 1)]);
    }

    #[test]
    fn ignoring_the_safety_checks_relabels_a_genotype() {
        let shards = vec![
            shard("first", &["s0"], &[("chr1", 100)]),
            shard("other", &["s1"], &[("chr1", 400)]),
        ];
        let error = gather(&shards, &Arguments::default()).unwrap_err();
        assert!(error
            .message()
            .contains("Samples unique to first file: [s0]"));

        // With the flag, the record is written under the first file's sample.
        let written = gather(
            &shards,
            &Arguments {
                ignore_safety_checks: true,
                ..Arguments::default()
            },
        )
        .expect("accepted");
        assert_eq!(written, vec![(0, 0), (1, 0)]);
    }
}
