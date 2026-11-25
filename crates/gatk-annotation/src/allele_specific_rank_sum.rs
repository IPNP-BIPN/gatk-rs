//! `AS_BaseQualityRankSumTest`, `AS_MappingQualityRankSumTest` and `AS_ReadPosRankSumTest`, and the
//! `AS_RankSumTest` machinery under them, from GATK 4.7.0.0; `AS_RankSumTest.java` and its three
//! subclasses are unchanged since 4.6.2.0.
//!
//! `AS_BaseQRankSum`, `AS_MQRankSum` and `AS_ReadPosRankSum`, plus the three `AS_RAW_` keys the
//! gVCF path writes so that combining two samples is adding two histograms.
//!
//! "The golden" below is `crates/gatk-annotation/tests/data/allele_specific_rank_sum.txt.gz`;
//! a row is cited by its kind and label, the first two or three tab-separated columns.
//!
//! # The direct path of an allele-specific annotation is **not** allele-specific
//!
//! ```java
//! final MannWhitneyU.Result result = mannWhitneyU.test(
//!         Doubles.toArray(altQuals), Doubles.toArray(refQuals),
//!         MannWhitneyU.TestType.FIRST_DOMINATES);
//! ```
//!
//! (One line in `AS_RankSumTest.java`, wrapped here.) `AS_RankSumTest.annotate` pools every
//! alternate allele's reads into one series and reports one Z score, exactly as its
//! non-allele-specific parent does. The two differ only in the key they write it under, so
//! `AS_MQRankSum` and `MQRankSum` on the same site carry the same number. The golden's
//! `multiallelic` rows show the pooling: one direct value, `AS_BaseQRankSum=0.651`, for a site
//! whose raw string carries two alternates, `|0.5,1|0.5,1`.
//!
//! The allele-specific part lives entirely in the raw path. [`annotate_raw_data`] runs only when
//! the annotator engine was built to write raw annotations
//! (`if (useRawAnnotations && annotationType instanceof ReducibleAnnotation)` in
//! `VariantAnnotatorEngine.addInfoAnnotations`), as the HaplotypeCaller's gVCF mode is;
//! [`combine_raw_data`] is the gVCF merge; [`finalize_raw_data`] is reached from GenotypeGVCFs
//! (through `VariantAnnotatorEngine.finalizeAnnotations`) and from GnarlyGenotyper.
//!
//! # The raw string starts with its delimiter, and the parser depends on that
//!
//! ```java
//! if (!vcAlleles.get(i).isReference()) {
//!     if (i != 0) { //strings will always start with a printDelim ...
//!         annotationString += AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM;
//!     }
//! ```
//!
//! The reference allele is skipped but its **slot is not**: the writer emits a leading `|`, so
//! splitting the string yields an empty first token that the parser assigns back to the reference.
//! Position, not content, is what carries the allele identity through a gVCF, and a raw string
//! whose leading delimiter were trimmed would silently shift every allele by one.
//!
//! # A site with no unambiguous reference read produces **no** rank sum
//!
//! ```java
//! if (perAlleleValues.get(ref).isEmpty())
//!     return perAltRankSumResults;
//! ```
//!
//! Not a missing value per allele: an empty map, which makes the raw string a bare run of
//! delimiters. The raw key is still written: the golden's `alt-only` and `single-read` rows carry
//! `AS_RAW_BaseQRankSum=|`. The finalising step then has nothing to take a median of, and writes
//! the missing value for each alternate (the `empty-slots` row: `AS_BaseQRankSum=.,.`).
//!
//! # The opposite case writes the four characters `NaN` into the record
//!
//! A site with reference reads and **no** alternate read takes a rank sum of an empty series
//! against a non-empty one, which is `NaN`. `Histogram.add` drops a `NaN` silently, so the
//! histogram is empty, and an empty histogram renders as `Double.toString(Double.NaN)`. The raw
//! field is then `|NaN`, and the golden carries that row (`ref-only`). The guard in
//! [`make_combined_annotation_string`] keeps those four characters out of the *combined* string,
//! but there is no such guard on the way in.
//!
//! # Each Z score is stored as a one-entry histogram, so the raw form is lossy
//!
//! ```java
//! public String outputSingletonValueAsHistogram(final Double rankSumValue) {
//!     Histogram h = new Histogram();
//!     h.add(rankSumValue);
//!     return h.toString();
//! }
//! ```
//!
//! The Z score is binned to a tenth before it is written, so `-1.23` and `-1.27` both come back as
//! `-1.3`: the bin index is a **floor**, so a negative score bins away from zero and the binning is
//! not symmetric about it. The golden's `bin` rows show it: `-1.23` and `-1.27` give `-1.3,1`,
//! `1.23` and `1.27` give `1.2,1`. Combining gVCFs then takes the **median of the bins**, which is
//! why `AS_MQRankSum` from a joint call is not the `MQRankSum` of the pooled reads and cannot be.
//!
//! # Reading this file without reading Rust
//!
//! - Java's `Map<Allele, X>` is modelled as `Vec<(Allele, X)>`, a list of pairs looked up with
//!   `.iter().find(|(a, _)| a == allele)`. `find` returns the first pair whose allele matches, or
//!   nothing. The order of the list never reaches the output: every string is written by walking an
//!   allele list, as the Java walks `vcAlleles`.
//! - `Option<T>` is "a `T` or nothing", Java's nullable reference with the null check enforced.
//!   `Some(x)` is present, `None` is null. `if let Some(x) = y { .. }` runs the block only when `y`
//!   is present, and `let Some(x) = y else { continue };` takes the early exit when it is not.
//! - `Result<T, AsRankSumError>` is "a `T`, or one of the exceptions listed on [`AsRankSumError`]".
//!   `Err(..)` is the throw. A `?` after a call rethrows that call's error to this function's
//!   caller, and `.map_err(AsRankSumError::Histogram)?` wraps the histogram's own error first.
//! - `&[Allele]` is a borrowed, read-only view of an array the caller keeps.
//! - A function's last expression, written without `return` or `;`, is its value.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::histogram::{CompressedDataList, Histogram, HistogramError};
use gatk_engine::mann_whitney::{self, TestType};
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

use crate::rank_sum::{
    format_decimals, BaseQualityRankSumTest, MappingQualityRankSumTest, RankSumTest,
    ReadPosRankSumTest, INVALID_ELEMENT_FROM_READ,
};

/// `GATKVCFConstants.AS_BASE_QUAL_RANK_SUM_KEY` and its raw form.
pub const AS_BASE_QUAL_RANK_SUM_KEY: &str = "AS_BaseQRankSum";
pub const AS_RAW_BASE_QUAL_RANK_SUM_KEY: &str = "AS_RAW_BaseQRankSum";
/// `GATKVCFConstants.AS_MAP_QUAL_RANK_SUM_KEY` and its raw form.
pub const AS_MAP_QUAL_RANK_SUM_KEY: &str = "AS_MQRankSum";
pub const AS_RAW_MAP_QUAL_RANK_SUM_KEY: &str = "AS_RAW_MQRankSum";
/// `GATKVCFConstants.AS_READ_POS_RANK_SUM_KEY` and its raw form.
pub const AS_READ_POS_RANK_SUM_KEY: &str = "AS_ReadPosRankSum";
pub const AS_RAW_READ_POS_RANK_SUM_KEY: &str = "AS_RAW_ReadPosRankSum";

/// `AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM`.
pub const ALLELE_SPECIFIC_RAW_DELIM: char = '|';
/// `AnnotationUtils.ALLELE_SPECIFIC_REDUCED_DELIM`.
pub const ALLELE_SPECIFIC_REDUCED_DELIM: char = ',';
/// `VCFConstants.MISSING_VALUE_v4`.
pub const MISSING_VALUE: &str = ".";
/// `AS_RankSumTest.RAW_DELIM`: the comma **inside** one allele's histogram, between a bin and its
/// count and between one pair and the next. The pipe separates alleles.
const RAW_DELIM: char = ',';

/// What this family refuses, one variant per reference exception. The port returns these where the
/// Java throws; which of them a given input reaches first follows the Java's order of checks.
#[derive(Debug, Clone, PartialEq)]
pub enum AsRankSumError {
    /// `IllegalStateException`: the raw path needs a variant context with exactly one sample, "as
    /// in a gVCF".
    NotExactlyOneSample { found: usize },
    /// `IllegalArgumentException` from `AlleleSpecificAnnotationData`: no reference allele, or more
    /// than one.
    ReferenceAlleleCount { found: usize },
    /// `ArrayIndexOutOfBoundsException`: the raw string carries more allele slots than the variant
    /// has alleles.
    TooManyAlleleSlots { slots: usize, alleles: usize },
    /// A bin index the histogram refuses: `Histogram`'s `GATKException`, "Histogram values are
    /// suspiciously extreme", or a non-positive count.
    Histogram(HistogramError),
}

/// `AS_BaseQualityRankSumTest`, `AS_MappingQualityRankSumTest` and `AS_ReadPosRankSumTest`.
///
/// The three members, which differ only in their keys and in the element they take from a read.
/// Java has three subclasses of an abstract class; here one `enum` (a closed set of named cases,
/// like a Java `enum`) and a `match` on it at each point where the subclasses differ. A `match`
/// must name every case, so a fourth member cannot be added and forgotten at one of those points.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsRankSum {
    BaseQuality,
    MappingQuality,
    ReadPosition,
}

impl AsRankSum {
    /// `getPrimaryRawKey()`.
    pub fn raw_key(&self) -> &'static str {
        match self {
            AsRankSum::BaseQuality => AS_RAW_BASE_QUAL_RANK_SUM_KEY,
            AsRankSum::MappingQuality => AS_RAW_MAP_QUAL_RANK_SUM_KEY,
            AsRankSum::ReadPosition => AS_RAW_READ_POS_RANK_SUM_KEY,
        }
    }

    /// `getEmptyRawValue()`, which is the empty string and not the missing value.
    pub fn empty_raw_value(&self) -> &'static str {
        ""
    }
}

/// The members' overrides of `RankSumTest`. Each `AS_` subclass's `getElementForRead` calls the
/// same code its non-allele-specific sibling uses (the static helpers
/// `BaseQualityRankSumTest.getReadBaseQuality` and `ReadPosRankSumTest.getReadPosition`, and the
/// same `OptionalDouble.of(read.getMappingQuality())`), so delegating to the sibling's port is the
/// reference's own reuse, not a shortcut.
impl RankSumTest for AsRankSum {
    /// `getKeyNames().get(0)`: the reduced key, `AS_BaseQRankSum` and so on.
    fn vcf_key(&self) -> &'static str {
        match self {
            AsRankSum::BaseQuality => AS_BASE_QUAL_RANK_SUM_KEY,
            AsRankSum::MappingQuality => AS_MAP_QUAL_RANK_SUM_KEY,
            AsRankSum::ReadPosition => AS_READ_POS_RANK_SUM_KEY,
        }
    }

    /// `getElementForRead(read, vc)`. `None` is `OptionalDouble.empty()`. The three-argument
    /// overload the raw path calls, `getElementForRead(read, vc, bestAllele)`, is `RankSumTest`'s
    /// default, which drops `bestAllele` and calls this one; no member overrides it.
    fn element_for_read(&self, read: &BamRecord, vc: &VariantContext) -> Option<f64> {
        match self {
            AsRankSum::BaseQuality => BaseQualityRankSumTest.element_for_read(read, vc),
            AsRankSum::MappingQuality => MappingQualityRankSumTest.element_for_read(read, vc),
            AsRankSum::ReadPosition => ReadPosRankSumTest.element_for_read(read, vc),
        }
    }

    /// `isUsableRead(read, vc)`.
    fn is_usable_read(&self, read: &BamRecord, vc: &VariantContext) -> bool {
        match self {
            // Only `AS_ReadPosRankSumTest` overrides the usability test, and it overrides it with
            // the same one its non-allele-specific sibling uses.
            AsRankSum::ReadPosition => ReadPosRankSumTest.is_usable_read(read, vc),
            AsRankSum::BaseQuality => BaseQualityRankSumTest.is_usable_read(read, vc),
            AsRankSum::MappingQuality => MappingQualityRankSumTest.is_usable_read(read, vc),
        }
    }
}

/// `AlleleSpecificAnnotationData`'s reference-allele check, which runs before any data is stored.
///
/// `checkRefAlleles` in the constructor: "multiple reference alleles found" for two or more, "no
/// reference alleles found" for none, both `IllegalArgumentException`. The returned `&Allele`
/// borrows from `alleles` rather than copying it, which is why callers that keep it write
/// `.clone()`.
fn reference_allele(alleles: &[Allele]) -> Result<&Allele, AsRankSumError> {
    // `.filter(..)` keeps the alleles the closure `|a| a.is_reference()` accepts (Java's lambda
    // `a -> a.isReference()`), and `.collect()` gathers them into a new list.
    let references: Vec<&Allele> = alleles.iter().filter(|a| a.is_reference()).collect();
    match references.len() {
        1 => Ok(references[0]),
        found => Err(AsRankSumError::ReferenceAlleleCount { found }),
    }
}

/// `AS_RankSumTest.calculateRawData`: one run-length list of **truncated** values per allele.
///
/// ```java
/// perAlleleValues.get(bestAllele.allele).add((int) value.getAsDouble());
/// ```
///
/// The cast to `int` is a truncation towards zero. For these three members it changes nothing,
/// because every value they produce is already integral: a base quality byte, a mapping quality,
/// and `Math.min(leftDistance, rightDistance)` of two `int`s for the read position. (The clipping
/// rank sum, which has no allele-specific member, is integral too:
/// `OptionalDouble.of(AlignmentUtils.getNumHardClippedBases(read))`, an `int`.)
///
/// The lists start out holding one empty entry per variant allele, which is
/// `initializeNewRawAnnotationData`; `None` likelihoods return them empty, which is the Java's
/// `if(likelihoods == null) {` exit **after** the one-sample check, so the check still fires.
pub fn calculate_raw_data(
    annotation: AsRankSum,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Result<Vec<(Allele, CompressedDataList)>, AsRankSumError> {
    // `vc.getGenotypes().getSampleNames().size() != 1`: a set of names, so two genotypes under one
    // name count once. The golden's `no-genotypes` and `two-samples` rows both throw here
    // (`E:java.lang.IllegalStateException`) while the direct `as` row for `two-samples` succeeds.
    let samples: std::collections::BTreeSet<&str> = vc
        .genotypes
        .iter()
        .map(|g| g.sample_name.as_str())
        .collect();
    if samples.len() != 1 {
        return Err(AsRankSumError::NotExactlyOneSample {
            found: samples.len(),
        });
    }
    let mut per_allele: Vec<(Allele, CompressedDataList)> = vc
        .alleles
        .iter()
        .map(|allele| (allele.clone(), CompressedDataList::new()))
        .collect();
    let Some(likelihoods) = likelihoods else {
        return Ok(per_allele);
    };

    // `likelihoods.bestAllelesBreakingTies()`, every sample's reads in the matrix's own order. The
    // Java's single `if (a && b && c ...)` is split into one early `continue` per condition, in the
    // same order, so the same reads are skipped.
    for best in likelihoods.best_alleles_breaking_ties(None) {
        if !best.is_informative() {
            continue;
        }
        // `bestAllele.evidence`: the Java object holds the read itself, the port holds its sample
        // and index, so the read is fetched back from the matrix. `.and_then(..)` chains a second
        // lookup that may also come back empty; either miss skips the read. The sample name came
        // from this matrix, so the `.unwrap_or(0)` fallback for an unknown name is never taken.
        let Some(read) = likelihoods
            .sample_evidence(likelihoods.index_of_sample(&best.sample).unwrap_or(0))
            .and_then(|reads| reads.get(best.evidence_index))
        else {
            continue;
        };
        if !annotation.is_usable_read(read, vc) {
            continue;
        }
        // `value.isPresent()`.
        let Some(value) = annotation.element_for_read(read, vc) else {
            continue;
        };
        if value == INVALID_ELEMENT_FROM_READ {
            continue;
        }
        let Some(allele) = best.allele else { continue };
        // `perAlleleValues.containsKey(bestAllele.allele)`: an allele the matrix holds and the
        // variant does not is dropped rather than counted into a new bucket.
        // `.iter_mut()` walks the list with write access; `slot.1` is the second half of the pair.
        // `value as i32` is the Java `(int)` cast: a float-to-integer `as` truncates towards zero
        // and saturates out of range, which is Java's narrowing rule for `double` to `int` too.
        if let Some(slot) = per_allele.iter_mut().find(|(a, _)| *a == allele) {
            slot.1.add(value as i32);
        }
    }
    Ok(per_allele)
}

/// `AS_RankSumTest.calculateRankSum`: one Z score per alternate, or nothing at all.
///
/// Each alternate's reads are tested against the reference's with `TestType.FIRST_DOMINATES`,
/// alternates as the first series, exactly as the Java loads "alts (series 1)" and "refs (series
/// 2)". An alternate with no reads still gets a test, whose Z is `NaN` (see the module note). The
/// Java iterates a `HashMap`'s keys; the order cannot matter because each result goes into a map
/// and every test is independent.
pub fn calculate_rank_sum(
    per_allele: &[(Allele, CompressedDataList)],
    reference: &Allele,
) -> Vec<(Allele, f64)> {
    // `(double) qual.intValue()` for each stored value. `.map(..).unwrap_or_default()` turns a
    // missing reference entry into an empty list; the Java would throw there instead, but the
    // entries are built from the same allele list, so the reference is always present.
    let ref_values: Vec<f64> = per_allele
        .iter()
        .find(|(a, _)| a == reference)
        .map(|(_, list)| list.iter().map(|v| v as f64).collect())
        .unwrap_or_default();
    // "shortcut to not try to calculate rank sum if there are no reads that unambiguously support
    // the ref": every alternate is dropped, not just the ones with no reads.
    if ref_values.is_empty() {
        return Vec::new();
    }
    // `!alt.equals(ref, false)`: the `false` means the reference flag is **not** ignored, which is
    // exactly the port's `!=` on `Allele` (bases and reference flag both compared).
    per_allele
        .iter()
        .filter(|(allele, _)| allele != reference)
        .map(|(allele, list)| {
            let alts: Vec<f64> = list.iter().map(|v| v as f64).collect();
            let result = mann_whitney::test(&alts, &ref_values, TestType::FirstDominates);
            (allele.clone(), result.z)
        })
        .collect()
}

/// `AS_RankSumTest.outputSingletonValueAsHistogram`: a Z score binned to a tenth.
///
/// A `NaN` Z adds nothing and the empty histogram renders as `NaN`; a finite one renders as
/// `bin,1`. The only refusal is a bin index outside `int`, "Histogram values are suspiciously
/// extreme".
pub fn singleton_histogram(z: f64) -> Result<String, AsRankSumError> {
    let mut histogram = Histogram::new();
    histogram.add(z).map_err(AsRankSumError::Histogram)?;
    Ok(histogram.to_string())
}

/// `AS_RankSumTest.makeRawAnnotationString`: the leading delimiter is deliberate. See the module
/// note.
///
/// The delimiter test is on the **position** `i != 0`, not on whether anything has been written,
/// so with the reference first every alternate, including the first, is preceded by a pipe. A
/// reference that was not first would move the pipes, in the port as in the Java.
pub fn make_raw_annotation_string(
    alleles: &[Allele],
    per_allele: &[(Allele, f64)],
) -> Result<String, AsRankSumError> {
    let mut out = String::new();
    // `.enumerate()` pairs each allele with its index, Java's `for (int i = 0; ...)` with
    // `vcAlleles.get(i)`.
    for (i, allele) in alleles.iter().enumerate() {
        if allele.is_reference() {
            continue;
        }
        if i != 0 {
            out.push(ALLELE_SPECIFIC_RAW_DELIM);
        }
        // Null when `calculateRankSum` returned nothing, in which case only the delimiters are
        // written and the field is a run of pipes.
        if let Some((_, z)) = per_allele.iter().find(|(a, _)| a == allele) {
            out.push_str(&singleton_histogram(*z)?);
        }
    }
    Ok(out)
}

/// `AS_RankSumTest.annotateRawData`: the raw key and its string, or `None` for the Java's empty
/// map.
///
/// The order of the steps is the Java's and decides which refusal wins: the null-matrix exit, then
/// the reference-allele check (`initializeNewRawAnnotationData` constructs an
/// `AlleleSpecificAnnotationData`), then the one-sample check inside `calculateRawData`.
pub fn annotate_raw_data(
    annotation: AsRankSum,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Result<Option<(String, String)>, AsRankSumError> {
    // `if ( likelihoods == null) {` then `return Collections.emptyMap();` happens **before** the
    // one-sample check, so a null matrix is an empty answer and not an exception (the golden's
    // `null-likelihoods` rows).
    if likelihoods.is_none() {
        return Ok(None);
    }
    let reference = reference_allele(&vc.alleles)?.clone();
    let per_allele = calculate_raw_data(annotation, vc, likelihoods)?;
    let rank_sums = calculate_rank_sum(&per_allele, &reference);
    let text = make_raw_annotation_string(&vc.alleles, &rank_sums)?;
    Ok(Some((annotation.raw_key().to_string(), text)))
}

/// `AS_RankSumTest.parseRawDataString`: positional, over the alleles the caller supplies.
///
/// The raw string's tokens line up with the allele list index for index, the reference's token
/// being the empty one the leading delimiter creates. More tokens than alleles is the reference's
/// `ArrayIndexOutOfBoundsException`.
///
/// Java's `parseRawDataString` takes an `AlleleSpecificAnnotationData` whose constructor has
/// already checked the reference allele; the port takes the allele list and runs that check first.
///
/// Malformed input the Java throws on is not all refused here. The Java reads
/// `rawListEntriesAsStringVector[j + 1]` without a bounds test and calls `Integer.parseInt` on it;
/// the port skips a bin with no count and an unparsable bin, and reads an unparsable count as zero,
/// which `add_count` then refuses as non-positive (a refusal, but not the Java's
/// `NumberFormatException`). `rawDataString.charAt(0)` throws on an empty string; `starts_with`
/// does not. No writer of this field produces any of these.
pub fn parse_raw_data_string(
    alleles: &[Allele],
    raw: &str,
) -> Result<Vec<(Allele, Histogram)>, AsRankSumError> {
    // The value is discarded: only the refusal matters here.
    reference_allele(alleles)?;
    let mut per_allele: Vec<(Allele, Histogram)> = alleles
        .iter()
        .map(|allele| (allele.clone(), Histogram::new()))
        .collect();

    // "Map gives back list with []": only a **leading** bracket is tested, and then both ends are
    // trimmed, so a string with a trailing bracket and no leading one keeps it. `&raw[1..len - 1]`
    // is `substring(1, length() - 1)`, a borrowed slice of the same text rather than a copy. The
    // golden's `bracketed` rows parse the same as the unbracketed ones.
    let without_brackets = if raw.starts_with('[') {
        &raw[1..raw.len() - 1]
    } else {
        raw
    };

    // `rawDataNoBrackets.split(AnnotationUtils.ALLELE_SPECIFIC_SPLIT_REGEX)`, no limit argument.
    let tokens = split_dropping_trailing_empties(without_brackets, ALLELE_SPECIFIC_RAW_DELIM);
    if tokens.len() > alleles.len() {
        return Err(AsRankSumError::TooManyAlleleSlots {
            slots: tokens.len(),
            alleles: alleles.len(),
        });
    }
    for (index, token) in tokens.iter().enumerate() {
        // `alleleData.split(RAW_DELIM)`. Rust's `split` keeps trailing empty pieces where Java's
        // drops them; the extra pieces are empty and the `is_empty` test below skips them.
        let entries: Vec<&str> = token.split(RAW_DELIM).collect();
        // Pairs of `bin,count`: `j` steps by two, Java's `j+=2`.
        let mut j = 0usize;
        while j < entries.len() {
            if !entries[j].is_empty() {
                let Ok(value) = entries[j].trim().parse::<f64>() else {
                    // `Double.parseDouble` throws a NumberFormatException, which the reference does
                    // not catch. No writer of this field produces one.
                    j += 2;
                    continue;
                };
                // A NaN bin is skipped here, which is the only reason `Histogram.add(d, count)`
                // never sees one despite having no guard of its own.
                if !value.is_nan() && j + 1 < entries.len() && !entries[j + 1].is_empty() {
                    let count: i32 = entries[j + 1].trim().parse().unwrap_or(0);
                    per_allele[index]
                        .1
                        .add_count(value, count)
                        .map_err(AsRankSumError::Histogram)?;
                }
            }
            j += 2;
        }
    }
    Ok(per_allele)
}

/// `String.split(regex)`, which drops trailing empty tokens but keeps leading and interior ones.
///
/// So `"|"` splits to no tokens at all and `"||x"` to three. Rust's own `split` keeps every empty
/// piece, which is why the trailing ones are popped off by hand. The leading empty token must
/// survive: it is the reference's slot.
fn split_dropping_trailing_empties(text: &str, delimiter: char) -> Vec<String> {
    if text.is_empty() {
        // Java's split on an empty input yields one empty token, not none.
        return vec![String::new()];
    }
    let mut parts: Vec<String> = text.split(delimiter).map(|s| s.to_string()).collect();
    // `.last()` is the final element or nothing; `.is_some_and(..)` is true only when there is one
    // and the closure accepts it.
    while parts.last().is_some_and(|part| part.is_empty()) {
        parts.pop();
    }
    parts
}

/// `AS_RankSumTest.combineRawData`, with `combineAttributeMap` folded in: the histograms of every
/// source summed bin by bin.
///
/// One allele list serves every source here. In the Java each source is an
/// `AlleleSpecificAnnotationData` built by `ReferenceConfidenceVariantContextMerger` over that
/// source's own `vcPair.getNewAlleles()`, and only the merged list is `vcAlleles`. A caller whose
/// sources carry different allele lists than the merged one would parse them positionally against
/// the wrong alleles; the golden's `ascombine` rows all use one list.
pub fn combine_raw_data(
    alleles: &[Allele],
    raw_strings: &[String],
) -> Result<String, AsRankSumError> {
    let mut combined: Vec<(Allele, Histogram)> = alleles
        .iter()
        .map(|allele| (allele.clone(), Histogram::new()))
        .collect();
    for raw in raw_strings {
        let parsed = parse_raw_data_string(alleles, raw)?;
        // `if (toAdd.hasAttribute(a))`: an allele the source lacks adds nothing.
        for (allele, histogram) in &parsed {
            if let Some(slot) = combined.iter_mut().find(|(a, _)| a == allele) {
                slot.1
                    .add_histogram(histogram)
                    .map_err(AsRankSumError::Histogram)?;
            }
        }
    }
    make_combined_annotation_string(alleles, &combined)
}

/// `AS_RankSumTest.makeCombinedAnnotationString`.
///
/// The same position-keyed delimiter as [`make_raw_annotation_string`], with one difference: an
/// empty histogram writes nothing rather than `NaN`, so `||` (the golden's `empty-slots` row)
/// rather than `|NaN|NaN`. Where the Java's `perAlleleValues.get(..)` would be null and
/// `alleleValue.isEmpty()` would throw a `NullPointerException`, the port writes nothing; see
/// [`finalize_raw_data`] for the one caller that can reach it.
pub fn make_combined_annotation_string(
    alleles: &[Allele],
    per_allele: &[(Allele, Histogram)],
) -> Result<String, AsRankSumError> {
    let mut out = String::new();
    for (i, allele) in alleles.iter().enumerate() {
        if allele.is_reference() {
            continue;
        }
        if i != 0 {
            out.push(ALLELE_SPECIFIC_RAW_DELIM);
        }
        if let Some((_, histogram)) = per_allele.iter().find(|(a, _)| a == allele) {
            // An empty histogram would render as `NaN`, so the guard is what keeps those four
            // characters out of the field.
            if !histogram.is_empty() {
                out.push_str(&histogram.to_string());
            }
        }
    }
    Ok(out)
}

/// `AS_RankSumTest.finalizeRawData`: the per-allele medians, and the combined raw string.
///
/// `None` when the raw key is absent. The reduced string is written over the **current** variant's
/// alternates while the histograms were parsed over the **original** variant's alleles, so an
/// alternate that survived trimming under a different representation is written as the missing
/// value rather than dropped, keeping the field's arity equal to the alternate count. The combined
/// raw string is written over the current alleles as well, `makeCombinedAnnotationString(
/// vc.getAlleles(), ...)`, so an alternate the genotyper dropped leaves it.
///
/// That second walk is where the two languages part on the trimming case. The Java's
/// `makeCombinedAnnotationString` calls `isEmpty()` on the map's value for each current alternate,
/// and an alternate absent from the original alleles has no entry, so the Java throws a
/// `NullPointerException` after having written the missing value for it in the reduced string.
/// The port writes nothing for that alternate and returns both keys.
///
/// The four strings returned are the reduced key and value, then the raw key and value, the two
/// entries the Java puts in its map, grouped as one tuple (a fixed-size group of values written
/// `(a, b, c, d)`).
pub fn finalize_raw_data(
    annotation: AsRankSum,
    vc_alternates: &[Allele],
    original_alleles: &[Allele],
    raw: Option<&str>,
) -> Result<Option<(String, String, String, String)>, AsRankSumError> {
    // `!vc.hasAttribute(getPrimaryRawKey())` and the null string both come in as `None`.
    let Some(raw) = raw else { return Ok(None) };
    let reference = reference_allele(original_alleles)?.clone();
    let per_allele = parse_raw_data_string(original_alleles, raw)?;

    // `calculateReducedData`: `Histogram.median()` per alternate, `None` for an empty histogram,
    // which is the Java's `return null` after a loop over no bins. The golden's `hist` rows cover
    // the median itself: `empty` gives `null`, `two-values` gives `0.15000000000000002`.
    let medians: Vec<(Allele, Option<f64>)> = per_allele
        .iter()
        .filter(|(allele, _)| *allele != reference)
        .map(|(allele, histogram)| (allele.clone(), histogram.median()))
        .collect();
    // "shortcut for no ref values": an empty map here is an empty answer, both keys absent.
    if medians.is_empty() {
        return Ok(None);
    }

    // `makeReducedAnnotationString`. The delimiter is written before the lookup, so every
    // alternate takes a slot whatever the lookup finds: the arity is always the alternate count.
    let mut reduced = String::new();
    for allele in vc_alternates {
        if !reduced.is_empty() {
            reduced.push(ALLELE_SPECIFIC_REDUCED_DELIM);
        }
        // `match` on the lookup names all three outcomes: no entry, an entry holding null, an
        // entry holding a median. `String.format("%.3f", ..)` is `format_decimals(.., 3)`.
        match medians.iter().find(|(a, _)| a == allele) {
            // "VC allele not found in annotation alleles -- maybe there was trimming?"
            None => reduced.push_str(MISSING_VALUE),
            Some((_, None)) => reduced.push_str(MISSING_VALUE),
            Some((_, Some(median))) => reduced.push_str(&format_decimals(*median, 3)),
        }
    }
    // Over the CURRENT alleles too, so an alternate `--max-alternate-alleles` dropped leaves the
    // raw key `--keep-combined-raw-annotations` writes. The Java passes `vc.getAlleles()`; the
    // port rebuilds that list from the original reference and the current alternates, which
    // writes the same string because the reference's slot is skipped either way. `once(x)` is a
    // one-element sequence and `.chain(..)` appends the alternates after it.
    let current: Vec<Allele> = std::iter::once(reference.clone())
        .chain(vc_alternates.iter().cloned())
        .collect();
    let combined = make_combined_annotation_string(&current, &per_allele)?;
    Ok(Some((
        annotation.vcf_key().to_string(),
        reduced,
        annotation.raw_key().to_string(),
        combined,
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles() -> Vec<Allele> {
        vec![
            Allele::from_str("A", true).expect("an allele"),
            Allele::from_str("C", false).expect("an allele"),
            Allele::from_str("G", false).expect("an allele"),
        ]
    }

    #[test]
    fn a_dropped_alternate_leaves_the_combined_raw_string() {
        let alleles = alleles();
        let finalized = finalize_raw_data(
            AsRankSum::BaseQuality,
            &alleles[1..2],
            &alleles,
            Some("|-0.5,1|-0.5,1"),
        )
        .expect("a finalized pair")
        .expect("both keys");
        assert_eq!(finalized.3, "|-0.5,1");
    }

    #[test]
    fn the_raw_string_starts_with_its_delimiter_and_round_trips() {
        let alleles = alleles();
        let raw = make_raw_annotation_string(
            &alleles,
            &[(alleles[1].clone(), -1.23), (alleles[2].clone(), 0.4)],
        )
        .expect("a raw string");
        // Binned to a tenth by a **floor**, so a negative Z score bins away from zero: -1.23
        // becomes -1.3 and not -1.2. The binning is not symmetric about zero.
        assert_eq!(raw, "|-1.3,1|0.4,1");
        let parsed = parse_raw_data_string(&alleles, &raw).expect("a parse");
        assert!(parsed[0].1.is_empty(), "the reference keeps an empty slot");
        let median = parsed[1].1.median().expect("a median");
        assert!((median + 1.3).abs() < 1e-12, "{median}");
    }

    #[test]
    fn a_site_with_no_reference_read_produces_only_delimiters() {
        let alleles = alleles();
        let raw = make_raw_annotation_string(&alleles, &[]).expect("a raw string");
        assert_eq!(raw, "||");
    }

    #[test]
    fn the_reference_allele_count_is_checked_before_anything_is_stored() {
        let two_references = vec![
            Allele::from_str("A", true).expect("an allele"),
            Allele::from_str("C", true).expect("an allele"),
        ];
        assert_eq!(
            parse_raw_data_string(&two_references, "|"),
            Err(AsRankSumError::ReferenceAlleleCount { found: 2 })
        );
    }
}
