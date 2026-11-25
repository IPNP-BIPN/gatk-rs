//! `SequenceDictionaryUtils`: whether two sequence dictionaries may be used together (GATK
//! 4.7.0.0; `SequenceDictionaryUtils.java` is unchanged since 4.6.2.0, and the two exception
//! classes and `prettyPrintSequenceRecords` lie outside every hunk of the 4.6.2.0 to 4.7.0.0
//! diffs of `UserException.java` and `ReadUtils.java`).
//!
//! A walker given both a reads input and a feature input compares their dictionaries **before the
//! traversal**, so a pair that shares no contig is refused rather than counted. That is what four
//! rows of `CountVariants`' covering array disagreed on, and nothing else in them (#1038).
//!
//! # Eight outcomes, three of which need to be asked for
//!
//! `NON_CANONICAL_HUMAN_ORDER`, `OUT_OF_ORDER` and `DIFFERENT_INDICES` are only ever returned when
//! the caller passed `check_contig_ordering`. Without it the same pairs fall through to
//! `IDENTICAL`, `SUPERSET` or `COMMON_SUBSET`: every such pair in the golden (`reversed`,
//! `different-indices` and the three `lexicographic-*` rows of
//! `sequence_dictionary_validation.txt.gz`) comes back `SUPERSET`, and a pair whose second
//! dictionary also has a contig the first lacks would come back `COMMON_SUBSET` instead. So an
//! argument decides whether a dictionary pair is accepted, not the dictionaries.
//!
//! # A length of zero is equivalent to any length
//!
//! `sequenceRecordsAreEquivalent` compares lengths only when both are non-zero, so a contig
//! declared without one agrees with the same contig declared with one, and two dictionaries that
//! differ in exactly that way are `IDENTICAL`.
//!
//! # Two empty dictionaries share no contigs
//!
//! `NO_COMMON_CONTIGS` is the empty intersection, and an empty set intersected with an empty set is
//! empty. A tool handed two dictionaries with nothing in them refuses them as incompatible.
//!
//! # Two of the messages end in two full stops
//!
//! `IncompatibleSequenceDictionaries` appends `.` to whatever reason it is handed, and the reasons
//! for `OUT_OF_ORDER` and `DIFFERENT_INDICES` already end in one. The reference prints `files..`
//! and `dictionaries..`, and a port that tidied that up would differ from it on every such run.
//!
//! # The human check is a length as much as a name
//!
//! `nonCanonicalHumanContigOrder` looks for chr1, chr2 and chr10 by name **and** by length, under
//! four assemblies' spellings. The same three names in the same wrong order with lengths that
//! none of hg18, hg19, b36 and b37 has are not human at all, and come back `OUT_OF_ORDER` instead
//! (golden row `lexicographic-not-human`, lengths 100, 300 and 200).
//!
//! Ported from `org.broadinstitute.hellbender.utils.SequenceDictionaryUtils`,
//! `org.broadinstitute.hellbender.exceptions.UserException$IncompatibleSequenceDictionaries` and
//! `org.broadinstitute.hellbender.utils.read.ReadUtils.prettyPrintSequenceRecords`.
//!
//! # Lookups are by primary name
//!
//! `SAMSequenceDictionary.getSequence(name)` also answers a contig's `AN` aliases, and the port's
//! `find` does not. That cannot change an answer here: the common contigs come from
//! `getContigNames`, which is primary names only, and the one lookup that could land on an
//! alias (`supersets`) is followed by `sequenceRecordsAreEquivalent`, which requires the two
//! names to be equal.
//!
//! # Reading this file without reading Rust
//!
//! - `&[SequenceRecord]` is a borrowed, read-only view of a dictionary's records in order; the
//!   caller keeps ownership. A record's position in the slice is its `getSequenceIndex()`.
//! - `Option<T>` is a value or nothing, and `Result<(), DictionaryRefusal>` is "returned
//!   normally, or threw this": `Ok(())` is the Java's plain `return`, `Err(..)` the `throw`.
//! - `.iter().filter(..).map(..).collect()`, `.find`, `.any`, `.all`, `.position` are Java stream
//!   operations; `|x| ...` is a closure, Java's lambda.
//! - `match` is a `switch` the compiler checks covers every case, which is why `validate` needs
//!   no `default:` branch for the Java's `GATKException`.
//! - A function's last expression is its return value; `return` is only written for early exits.

use htsjdk_bam::header::SequenceRecord;

/// `SequenceDictionaryCompatibility`, the variants in the Java's declaration order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compatibility {
    Identical,
    CommonSubset,
    Superset,
    NoCommonContigs,
    UnequalCommonContigs,
    NonCanonicalHumanOrder,
    OutOfOrder,
    DifferentIndices,
}

impl Compatibility {
    /// The enum constant's own name, which is what the golden carries.
    pub fn name(self) -> &'static str {
        match self {
            Compatibility::Identical => "IDENTICAL",
            Compatibility::CommonSubset => "COMMON_SUBSET",
            Compatibility::Superset => "SUPERSET",
            Compatibility::NoCommonContigs => "NO_COMMON_CONTIGS",
            Compatibility::UnequalCommonContigs => "UNEQUAL_COMMON_CONTIGS",
            Compatibility::NonCanonicalHumanOrder => "NON_CANONICAL_HUMAN_ORDER",
            Compatibility::OutOfOrder => "OUT_OF_ORDER",
            Compatibility::DifferentIndices => "DIFFERENT_INDICES",
        }
    }
}

/// What `validateDictionaries` throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DictionaryRefusal {
    /// `UserException$IncompatibleSequenceDictionaries`, which wraps every message but one.
    Incompatible { message: String },
    /// `UserException$LexicographicallySortedSequenceDictionary`, which names one dictionary.
    LexicographicallySorted { name: String, contigs: String },
}

impl DictionaryRefusal {
    /// The class the reference throws, as the golden prints it before the message.
    pub fn java_class(&self) -> &'static str {
        match self {
            DictionaryRefusal::Incompatible { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$IncompatibleSequenceDictionaries"
            }
            DictionaryRefusal::LexicographicallySorted { .. } => {
                "org.broadinstitute.hellbender.exceptions.UserException$LexicographicallySortedSequenceDictionary"
            }
        }
    }

    /// The exception's message.
    ///
    /// For the lexicographic case this is the `LexicographicallySortedSequenceDictionary`
    /// constructor's format string, with
    /// `HelpConstants.forumPost("discussion/58/companion-utilities-reordersam")` already
    /// expanded: `LEGACY_GATK_FORUM_URL` is `"http://gatkforums.broadinstitute.org/"`. A `\` at
    /// the end of a line inside a Rust string joins the next line and drops its leading spaces,
    /// so the spaces that matter are the ones written before each `\`.
    pub fn message(&self) -> String {
        match self {
            DictionaryRefusal::Incompatible { message } => message.clone(),
            DictionaryRefusal::LexicographicallySorted { name, contigs } => format!(
                "Lexicographically sorted human genome sequence detected in {name}.\n\
                 For safety's sake the GATK requires human contigs in karyotypic order: 1, 2, ..., \
                 10, 11, ..., 20, 21, 22, X, Y with M either leading or trailing these contigs.\n\
                 This is because all distributed GATK resources are sorted in karyotypic order, \
                 and your processing will fail when you need to use these files.\n\
                 You can use the ReorderSam utility to fix this problem: \
                 http://gatkforums.broadinstitute.org/discussion/58/companion-utilities-reordersam\n  \
                 {name} contigs = {contigs}"
            ),
        }
    }
}

/// `ReadUtils.prettyPrintSequenceRecords`, which is `Arrays.deepToString` on the NAMES.
///
/// So `[chr1, chr2]`: square brackets, comma and space, no lengths, and `[]` for an empty
/// dictionary.
pub fn pretty_print(dictionary: &[SequenceRecord]) -> String {
    format!(
        "[{}]",
        dictionary
            .iter()
            .map(|record| record.name.as_str())
            .collect::<Vec<_>>()
            .join(", ")
    )
}

/// `sequenceRecordsAreEquivalent`: the same name, and the same length unless one of them is zero.
///
/// The Java refuses only when `length1 != length2` and neither is
/// `SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH`, which is 0, then compares names. Attributes such
/// as `M5` are not compared, unlike `SAMSequenceRecord.isSameSequence`; the Java's `supersets`
/// says why: `Cannot rely on SAMSequenceRecord.equals() as it's too strict`.
fn equivalent(a: &SequenceRecord, b: &SequenceRecord) -> bool {
    a.name == b.name && (a.length == 0 || b.length == 0 || a.length == b.length)
}

/// `SAMSequenceDictionary.getSequence(name)`, by primary name (see the module note on aliases).
/// The lifetime `'a` says the record returned is borrowed from `dictionary`.
fn find<'a>(dictionary: &'a [SequenceRecord], name: &str) -> Option<&'a SequenceRecord> {
    dictionary.iter().find(|record| record.name == name)
}

/// `getCommonContigsByName`, in the FIRST dictionary's order.
///
/// The order matters for the message a disequal pair produces: the reference reports the first
/// mismatch it walks into, and it walks the intersection the same way every time. The Java's set
/// is a `LinkedHashSet` of the first dictionary's names, then
/// `intersectingSequenceNames.retainAll(getContigNames(dict2));`, which keeps that order.
fn common_contigs(a: &[SequenceRecord], b: &[SequenceRecord]) -> Vec<String> {
    a.iter()
        .filter(|record| find(b, &record.name).is_some())
        .map(|record| record.name.clone())
        .collect()
}

/// `hg18`, `hg19`, `b36` and `b37`'s chr1, chr2 and chr10, by the names and lengths the check uses.
///
/// The Java's `CHR1_HG18`, `CHR1_HG19`, `CHR1_B36` and `CHR1_B37` (and the same for 2 and 10),
/// in that order. The b assemblies reuse the hg lengths under bare names, so b36 is hg18's
/// lengths and b37 is hg19's.
const HUMAN_CHR1: [(&str, i32); 4] = [
    ("chr1", 247249719),
    ("chr1", 249250621),
    ("1", 247249719),
    ("1", 249250621),
];
const HUMAN_CHR2: [(&str, i32); 4] = [
    ("chr2", 242951149),
    ("chr2", 243199373),
    ("2", 242951149),
    ("2", 243199373),
];
const HUMAN_CHR10: [(&str, i32); 4] = [
    ("chr10", 135374737),
    ("chr10", 135534747),
    ("10", 135374737),
    ("10", 135534747),
];

/// `isHumanSeqRecord(elt, recs...)`: an exact name and an exact length, both.
fn is_human(record: &SequenceRecord, candidates: &[(&str, i32); 4]) -> bool {
    candidates
        .iter()
        .any(|(name, length)| record.name == *name && record.length == *length)
}

/// `nonCanonicalHumanContigOrder`: chr1 before chr2 before chr10, when all three are recognised.
///
/// The Java walks the whole dictionary and assigns `chr1 = elt` on every match, so a dictionary
/// with two records that both pass (`chr1` and `1` with one length) is judged by the LAST.
/// `.next_back()` takes the last of the matching indices for the same reason.
pub fn non_canonical_human_order(dictionary: &[SequenceRecord]) -> bool {
    let index = |candidates: &[(&str, i32); 4]| {
        dictionary
            .iter()
            .enumerate()
            .filter(|(_, record)| is_human(record, candidates))
            .map(|(index, _)| index)
            .next_back()
    };
    match (index(&HUMAN_CHR1), index(&HUMAN_CHR2), index(&HUMAN_CHR10)) {
        (Some(one), Some(two), Some(ten)) => !(one < two && two < ten),
        // A dictionary missing any of the three is not judged: the check has nothing to go on.
        _ => false,
    }
}

/// `supersets(dict1, dict2)`: whether the first dictionary holds an equivalent record for every
/// one of the second's.
fn supersets(a: &[SequenceRecord], b: &[SequenceRecord]) -> bool {
    b.iter()
        .all(|record| find(a, &record.name).is_some_and(|theirs| equivalent(record, theirs)))
}

/// `findDisequalCommonContigs`: the first common contig whose lengths disagree, as the pair the
/// message names. The Java returns `null` for none, here `None`; `find_map` stops at the first
/// closure that answers `Some`, and the `?` inside it skips a name either side lacks (which
/// cannot happen for a common contig).
fn disequal_common_contig<'a>(
    common: &[String],
    a: &'a [SequenceRecord],
    b: &'a [SequenceRecord],
) -> Option<(&'a SequenceRecord, &'a SequenceRecord)> {
    common.iter().find_map(|name| {
        let (one, two) = (find(a, name)?, find(b, name)?);
        (!equivalent(one, two)).then_some((one, two))
    })
}

/// `commonContigsAreInSameRelativeOrder`: whether the common contigs appear in the same relative
/// order in both dictionaries.
///
/// The Java looks each common name up, sorts the two lists by `getSequenceIndex()` and compares
/// names pairwise. Filtering each dictionary in its own order down to the common names is the
/// same list without the sort.
fn same_relative_order(common: &[String], a: &[SequenceRecord], b: &[SequenceRecord]) -> bool {
    let order = |dictionary: &[SequenceRecord]| -> Vec<String> {
        dictionary
            .iter()
            .filter(|record| common.contains(&record.name))
            .map(|record| record.name.clone())
            .collect()
    };
    order(a) == order(b)
}

/// `commonContigsAreAtSameIndices`: whether every common contig sits at the same absolute index in
/// both dictionaries, the index being the record's position.
fn same_indices(common: &[String], a: &[SequenceRecord], b: &[SequenceRecord]) -> bool {
    let index = |dictionary: &[SequenceRecord], name: &str| {
        dictionary.iter().position(|record| record.name == name)
    };
    common.iter().all(|name| index(a, name) == index(b, name))
}

/// `compareDictionaries`, whose branch order is what decides which of two true things is reported.
///
/// The order, exactly the Java's: the human check (only when asked), no common contigs, unequal
/// common contigs, then out of order (only when asked), identical, different indices (only when
/// asked), superset, and common subset as what is left. So `IDENTICAL` needs the same relative
/// order and nothing outside the intersection, and is tested before `DIFFERENT_INDICES`, which
/// therefore needs a contig outside the intersection to be reachable (golden row
/// `different-indices`: a leading `chrX`).
pub fn compare(
    a: &[SequenceRecord],
    b: &[SequenceRecord],
    check_contig_ordering: bool,
) -> Compatibility {
    if check_contig_ordering && (non_canonical_human_order(a) || non_canonical_human_order(b)) {
        return Compatibility::NonCanonicalHumanOrder;
    }

    let common = common_contigs(a, b);
    if common.is_empty() {
        // Two EMPTY dictionaries land here as well: an empty intersection is an empty
        // intersection, whether or not there was anything to intersect.
        return Compatibility::NoCommonContigs;
    }
    if disequal_common_contig(&common, a, b).is_some() {
        return Compatibility::UnequalCommonContigs;
    }

    let same_order = same_relative_order(&common, a, b);
    if check_contig_ordering && !same_order {
        Compatibility::OutOfOrder
    } else if same_order && common.len() == a.len() && common.len() == b.len() {
        Compatibility::Identical
    } else if check_contig_ordering && !same_indices(&common, a, b) {
        Compatibility::DifferentIndices
    } else if supersets(a, b) {
        Compatibility::Superset
    } else {
        Compatibility::CommonSubset
    }
}

/// `new UserException.IncompatibleSequenceDictionaries(...)`, which wraps a reason in the names
/// and both contig lists.
fn incompatible(
    reason: &str,
    name1: &str,
    a: &[SequenceRecord],
    name2: &str,
    b: &[SequenceRecord],
) -> DictionaryRefusal {
    DictionaryRefusal::Incompatible {
        message: incompatible_message(reason, name1, a, name2, b),
    }
}

/// The message of `UserException.IncompatibleSequenceDictionaries`.
///
/// ```java
/// "Input files %s and %s have incompatible contigs: %s.\n  %s contigs = %s\n  %s contigs = %s",
/// ```
///
/// The `.` after the reason is the format's, which is where the doubled full stops come from.
pub fn incompatible_message(
    reason: &str,
    name1: &str,
    a: &[SequenceRecord],
    name2: &str,
    b: &[SequenceRecord],
) -> String {
    format!(
        "Input files {name1} and {name2} have incompatible contigs: {reason}.\n  \
         {name1} contigs = {}\n  {name2} contigs = {}",
        pretty_print(a),
        pretty_print(b)
    )
}

/// `validateDictionaries`: `Ok` where the reference returns, and the refusal where it throws.
///
/// The Java's `Utils.nonNull` checks on the two dictionaries have no counterpart: a slice cannot
/// be null. `require_superset` only changes `COMMON_SUBSET`, and `check_contig_ordering` only
/// what `compare` can return.
pub fn validate(
    name1: &str,
    a: &[SequenceRecord],
    name2: &str,
    b: &[SequenceRecord],
    require_superset: bool,
    check_contig_ordering: bool,
) -> Result<(), DictionaryRefusal> {
    match compare(a, b, check_contig_ordering) {
        Compatibility::Identical | Compatibility::Superset => Ok(()),
        Compatibility::CommonSubset => {
            if !require_superset {
                return Ok(());
            }
            // Names of the second dictionary absent from the first, in the second's order, as
            // the Java's stream over `dict2.getSequences()` filtered by a set of `dict1`'s names.
            let missing: Vec<&str> = b
                .iter()
                .filter(|record| find(a, &record.name).is_none())
                .map(|record| record.name.as_str())
                .collect();
            // The two spaces after the full stop and the spaces around the newlines are the
            // format string's, not a rendering choice.
            Err(incompatible(
                &format!(
                    "Dictionary {name1} is missing contigs found in dictionary {name2}.  \
                     Missing contigs: \n {} \n",
                    missing.join(", ")
                ),
                name1,
                a,
                name2,
                b,
            ))
        }
        Compatibility::NoCommonContigs => Err(incompatible(
            "No overlapping contigs found",
            name1,
            a,
            name2,
            b,
        )),
        Compatibility::UnequalCommonContigs => {
            // The Java recomputes the pair here rather than carrying it out of `compare`, and so
            // does this. `.expect` would panic if there were none, which `compare` rules out.
            let common = common_contigs(a, b);
            let (one, two) = disequal_common_contig(&common, a, b)
                .expect("an unequal common contig, which is what this case is");
            Err(incompatible(
                &format!(
                    "Found contigs with the same name but different lengths:\n  contig {name1} = \
                     {} / {}\n  contig {name2} = {} / {}",
                    one.name, one.length, two.name, two.length
                ),
                name1,
                a,
                name2,
                b,
            ))
        }
        Compatibility::NonCanonicalHumanOrder => {
            // Whichever dictionary is the lexicographic one is the one named, and the first is
            // tested first.
            let (name, dictionary) = if non_canonical_human_order(a) {
                (name1, a)
            } else {
                (name2, b)
            };
            Err(DictionaryRefusal::LexicographicallySorted {
                name: name.to_string(),
                contigs: pretty_print(dictionary),
            })
        }
        Compatibility::OutOfOrder => Err(incompatible(
            &format!(
                "The relative ordering of the common contigs in {name1} and {name2} is not the \
                 same; to fix this please see: \
                 (https://www.broadinstitute.org/gatk/guide/article?id=1328),  which describes \
                 reordering contigs in BAM and VCF files."
            ),
            name1,
            a,
            name2,
            b,
        )),
        Compatibility::DifferentIndices => Err(incompatible(
            "One or more contigs common to both dictionaries have different indices (ie., \
             absolute positions) in each dictionary. Code that is sensitive to contig ordering \
             can fail when this is the case. You should fix the sequence dictionaries so that all \
             shared contigs occur at the same absolute positions in both dictionaries.",
            name1,
            a,
            name2,
            b,
        )),
    }
}
