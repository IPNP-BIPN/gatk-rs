//! `AS_FisherStrand`, `AS_StrandOddsRatio` and the `AS_StrandBiasTest` / `StrandBiasUtils`
//! machinery under them, from GATK 4.7.0.0; the four `.java` files are unchanged since 4.6.2.0.
//!
//! `AS_FS`, `AS_SOR` and the `AS_SB_TABLE` they are both derived from: one forward/reverse pair per
//! allele, summed across samples, carried through a gVCF as a pipe-separated list.
//!
//! "The golden" below is `crates/gatk-annotation/tests/data/allele_specific_strand_bias.txt.gz`;
//! a row is cited by its kind and label, the first two or three tab-separated columns.
//!
//! # The raw string has an entry for the **reference**, unlike the rank sums'
//!
//! ```java
//! for (final Allele a : vcAlleles) {
//!     if (!annotationString.isEmpty()) {
//!         annotationString += AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM;
//!     }
//! ```
//!
//! `AS_SB_TABLE` writes every allele including the reference and puts the delimiter **between**
//! entries, so it does not start with one. The rank sums skip the reference's value but keep its
//! slot, so theirs does. Two allele-specific families, two incompatible conventions for the same
//! delimiter, and the parser here refuses a count mismatch outright:
//!
//! ```java
//! throw new IllegalStateException(
//! "Number of alleles and number of allele-specific entries do not match.  " +
//! "Allele-specific annotations should have an entry for each allele including the reference.");
//! ```
//!
//! (Two lines in `AS_StrandBiasTest.parseRawDataString`, rewrapped here to fit.) The golden's
//! `wrong-count` row is that exception.
//!
//! # A sample contributes nothing at all unless it has **more than two** informative reads
//!
//! ```java
//! return readCount > minCount;
//! ```
//!
//! with `minCount` passed as `MIN_COUNT`, which is 2. Strictly greater, and the threshold is on the
//! sample's whole table rather than on any one allele. A sample with two informative reads is
//! dropped entirely, so its strand counts do not appear in `AS_SB_TABLE` even for the allele that
//! had both of them. The golden shows the edge: `two-reads` gives `AS_SB_TABLE=0,0|0,0` and
//! `three-reads` gives `2,1|0,0`.
//!
//! # `AS_SOR` computes a value for the reference allele and then never prints it
//!
//! `AS_StrandOddsRatio.calculateReducedData`, no filter:
//!
//! ```java
//! for (final Allele a : perAlleleData.keySet()) {
//!     List<Integer> altStrandCounts = perAlleleData.get(a);
//! ```
//!
//! `AS_FisherStrand.calculateReducedData`, filtered:
//!
//! ```java
//! for (final Allele a : perAlleleData.keySet()) {
//!     if(!a.equals(combinedData.getRefAllele(),true)) {
//! ```
//!
//! `AS_StrandOddsRatio.calculateReducedData` iterates every allele, so it takes a symmetric odds
//! ratio of the reference against itself. `makeReducedAnnotationString` then walks only the
//! alternates, so that value is computed and discarded. `AS_FisherStrand` filters the reference out
//! first. The two siblings differ in exactly that line.
//!
//! # `AS_FS` floors its p-value at `1e-320`, which is a subnormal
//!
//! ```java
//! QualityUtils.phredScaleErrorRate(
//!         Math.max(FisherStrand.pValueForContingencyTable(refAltTable), MIN_PVALUE))
//! ```
//!
//! (A fragment of one line of `AS_FisherStrand.calculateReducedData`, wrapped here.) `1.0E-320` is
//! below `Double.MIN_NORMAL`, so the floor is a denormalised double and the phred
//! scale of it is just above 3200: the golden's `phred 1.0E-320` row holds the bits of
//! `3200.00004834948`, the excess being the subnormal's rounding away from exactly `1e-320`. The
//! clamp inside `phredScaleLog10ErrorRate` is at `log10(Double.MIN_VALUE)`, about -323.3, which is
//! lower still, so the `MIN_PVALUE` floor is the one that binds and the value never reaches 3233
//! (the `phred 4.9E-324` and `phred 0.0` rows, `3233.062153431158`, are that clamp).
//!
//! # Reading this file without reading Rust
//!
//! - Java's `Map<Allele, List<Integer>>` is modelled as `Vec<(Allele, Option<Vec<i32>>)>`: a list
//!   of pairs, one per allele, found with `.iter().find(|(a, _)| a == allele)`. `Option` is the
//!   map's null: `None` is an allele whose value is null, `Some(vec![])` is one whose list is
//!   present but empty, and the file turns on the difference (`0,0` against an empty entry).
//! - `Result<T, AsStrandBiasError>` is "a `T`, or one of the exceptions on [`AsStrandBiasError`]";
//!   `?` after a call rethrows that call's error to this function's caller.
//! - `let Some(x) = y else { continue };` unwraps `y` or skips to the next loop iteration.
//! - `&[Allele]` is a borrowed, read-only view of an array the caller keeps; `&mut` is a writable
//!   one. `.clone()` makes an independent copy where the Java would share one object.
//! - A function's last expression, written without `return` or `;`, is its value.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

use crate::info_annotation::AnnotationValue;
use crate::rank_sum::format_decimals;
use crate::strand_bias::{
    calculate_sor, contingency_table, p_value_for_contingency_table, STRAND_BIAS_BY_SAMPLE_KEY,
};

/// `GATKVCFConstants.AS_SB_TABLE_KEY`.
pub const AS_SB_TABLE_KEY: &str = "AS_SB_TABLE";
/// `GATKVCFConstants.AS_FISHER_STRAND_KEY`.
pub const AS_FISHER_STRAND_KEY: &str = "AS_FS";
/// `GATKVCFConstants.AS_STRAND_ODDS_RATIO_KEY`.
pub const AS_STRAND_ODDS_RATIO_KEY: &str = "AS_SOR";

/// `StrandBiasUtils.FORWARD` and `REVERSE`, which are indices into a two-element list.
pub const FORWARD: usize = 0;
pub const REVERSE: usize = 1;
/// `StrandBiasUtils.MIN_COUNT`, equal to `AS_StrandBiasTest.MIN_COUNT`: both are 2.
pub const MIN_COUNT: i32 = 2;
/// `AS_StrandBiasTest.MIN_PVALUE`, the same `1E-320` as `FisherStrand.MIN_PVALUE`, which is the
/// floor the direct `AS_FS` path reaches through `FisherStrand.makeValueObjectForAnnotation`.
pub const MIN_PVALUE: f64 = 1.0E-320;
/// `AnnotationUtils.ALLELE_SPECIFIC_RAW_DELIM`, which `StrandBiasUtils.makeRawAnnotationString`
/// writes; `AS_StrandBiasTest.PRINT_DELIM` is the same `"|"`.
const RAW_DELIM: char = '|';
/// `AS_StrandBiasTest.REDUCED_DELIM`.
const REDUCED_DELIM: char = ',';
/// `VCFConstants.MISSING_VALUE_v4`.
const MISSING_VALUE: &str = ".";

/// What this family refuses, one variant per reference exception. The port returns these where the
/// Java throws.
#[derive(Debug, Clone, PartialEq)]
pub enum AsStrandBiasError {
    /// `IllegalStateException`: "Number of alleles and number of allele-specific entries do not
    /// match. Allele-specific annotations should have an entry for each allele including the
    /// reference."
    AlleleCountMismatch { entries: usize, alleles: usize },
    /// `NumberFormatException` from `Integer.parseInt` on a non-integer count.
    MalformedCount { raw: String },
    /// `IllegalArgumentException` from the `AlleleSpecificAnnotationData` constructor: no reference
    /// allele, or more than one.
    ReferenceAlleleCount { found: usize },
}

/// `AS_FisherStrand` and `AS_StrandOddsRatio`, the two subclasses of `AS_StrandBiasTest`.
///
/// The two members, which differ in their key and in whether they filter the reference out. One
/// `enum` (a closed set of named cases) stands for the two subclasses, and a `match` on it marks
/// each place they differ.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AsStrandBias {
    Fisher,
    OddsRatio,
}

impl AsStrandBias {
    /// `getKeyNames().get(0)`: `AS_FS` or `AS_SOR`.
    pub fn vcf_key(&self) -> &'static str {
        match self {
            AsStrandBias::Fisher => AS_FISHER_STRAND_KEY,
            AsStrandBias::OddsRatio => AS_STRAND_ODDS_RATIO_KEY,
        }
    }

    /// `getPrimaryRawKey()`, which is the **same** key for both: they share one raw table.
    pub fn raw_key(&self) -> &'static str {
        AS_SB_TABLE_KEY
    }

    /// `getEmptyRawValue()`, which is a pair of zeros rather than the empty string the rank sums
    /// use.
    pub fn empty_raw_value(&self) -> &'static str {
        "0,0"
    }
}

/// `AlleleSpecificAnnotationData`'s reference-allele check (`checkRefAlleles`), which runs in its
/// constructor before any data is parsed.
fn reference_allele(alleles: &[Allele]) -> Result<&Allele, AsStrandBiasError> {
    let references: Vec<&Allele> = alleles.iter().filter(|a| a.is_reference()).collect();
    match references.len() {
        1 => Ok(references[0]),
        found => Err(AsStrandBiasError::ReferenceAlleleCount { found }),
    }
}

/// `StrandBiasUtils.getStrandCountsFromLikelihoodMap`: one forward/reverse pair per allele.
///
/// Absent rather than zero: an allele no read was assigned to has **no** entry, which is what makes
/// the difference between `0,0` and an empty entry in the raw string.
///
/// **How**: per sample, a fresh table over the variant's alleles; each informative read adds one to
/// its best allele's forward or reverse count (`updateTable`); the sample's table is added into
/// the running one only if its total passes `readCount > minCount` (`passesMinimumThreshold`, then
/// `combineAttributeMap`). The golden's `two-samples-each-small` row gives `0,0|0,0` and
/// `two-samples-one-large` gives `5,5|5,5`: the threshold is per sample, not on the total.
///
/// One difference from the Java that no golden row reaches: `updateTable` matches the reference
/// with `bestAllele.equals(ref, true)`, which ignores the reference flag, while the port's `==`
/// compares it. A best allele with the reference's bases but no reference flag would count towards
/// the Java's sample threshold (under its own key, which `combineAttributeMap` then never copies)
/// and is dropped here. It takes a matrix allele with the reference's bases and no reference flag.
pub fn strand_counts(
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
    min_count: i32,
) -> Vec<(Allele, Option<Vec<i32>>)> {
    let mut combined: Vec<(Allele, Option<Vec<i32>>)> = vc
        .alleles
        .iter()
        .map(|allele| (allele.clone(), None))
        .collect();
    // `if( likelihoods == null || vc == null ) {`, then `return;`: the empty table.
    let Some(likelihoods) = likelihoods else {
        return combined;
    };
    // `vc.getReference()`, as an `Option` because the port's variant may lack one.
    let reference = vc.alleles.iter().find(|a| a.is_reference());

    // `for (final String sample : likelihoods.samples())`, by index rather than by name.
    for sample_index in 0..likelihoods.number_of_samples() {
        let mut sample_table: Vec<(Allele, Option<Vec<i32>>)> = vc
            .alleles
            .iter()
            .map(|allele| (allele.clone(), None))
            .collect();
        for best in likelihoods.best_alleles_breaking_ties_for_sample(sample_index, None) {
            if !best.is_informative() {
                continue;
            }
            let Some(allele) = best.allele.as_ref() else {
                continue;
            };
            // "can happen if a read's most likely allele has been removed when
            // --max_alternate_alleles is exceeded": neither the reference nor a declared
            // alternate, so the read is dropped rather than counted somewhere.
            let matches_reference = reference == Some(allele);
            let matches_alternate = vc
                .alternate_alleles()
                .iter()
                .any(|alternate| alternate == allele);
            if !(matches_reference || matches_alternate) {
                continue;
            }
            let Some(read) = likelihoods
                .sample_evidence(sample_index)
                .and_then(|reads| reads.get(best.evidence_index))
            else {
                continue;
            };
            let strand = if is_reverse_strand(read) {
                REVERSE
            } else {
                FORWARD
            };
            // `get_or_insert_with` returns the slot's list, first storing `[0, 0]` if the slot was
            // `None`: the Java's "existing list, or a new `0,0` one" followed by `putAttribute`.
            if let Some(slot) = sample_table.iter_mut().find(|(a, _)| a == allele) {
                let counts = slot.1.get_or_insert_with(|| vec![0, 0]);
                counts[strand] += 1;
            }
        }
        // The threshold is on the sample's whole table, and it is strictly greater than.
        // `.filter_map(..)` keeps only the present lists (Java's `filter(alleleValues ->
        // alleleValues != null)`), and `.sum()` adds the per-allele totals.
        let read_count: i32 = sample_table
            .iter()
            .filter_map(|(_, counts)| counts.as_ref())
            .map(|counts| counts[FORWARD] + counts[REVERSE])
            .sum();
        if read_count > min_count {
            // `combineAttributeMap(sampleTable, perAlleleValues)`: add into an existing pair, or
            // store a copy when the running table has none. A sample's list is never empty here,
            // so the Java's empty-list branch cannot arise on this path.
            for (allele, counts) in &sample_table {
                let Some(counts) = counts else { continue };
                if let Some(slot) = combined.iter_mut().find(|(a, _)| a == allele) {
                    // `&mut slot.1` borrows the slot's `Option` writably, so the `Some` arm can
                    // add into the list in place.
                    match &mut slot.1 {
                        Some(existing) => {
                            existing[FORWARD] += counts[FORWARD];
                            existing[REVERSE] += counts[REVERSE];
                        }
                        None => slot.1 = Some(counts.clone()),
                    }
                }
            }
        }
    }
    combined
}

/// `GATKRead.isReverseStrand()`: the `0x10` flag. `&` is the bitwise and.
fn is_reverse_strand(read: &BamRecord) -> bool {
    read.flags & 0x10 != 0
}

/// `StrandBiasUtils.makeRawAnnotationString`: the delimiter goes **between** entries, and a missing
/// allele is written as the shared `ZERO_LIST`.
///
/// An entry that is **present but empty** is not the same as an absent one: `encode` renders it as
/// the empty string, so the field can read `10,8|` with nothing after the delimiter. That state is
/// reachable, because `combineAttributeMap` stores an empty list whenever the source's own entry
/// was empty. The golden's `empty-entry` row writes `AS_SB_TABLE=10,8|`; `zero-entry` writes
/// `10,8|0,0`.
pub fn make_raw_annotation_string(per_allele: &[(Allele, Option<Vec<i32>>)]) -> String {
    let mut out = String::new();
    for (_, counts) in per_allele {
        if !out.is_empty() {
            out.push(RAW_DELIM);
        }
        match counts {
            // `ZERO_LIST` is a `static` mutable `ArrayList` shared by every call site; nothing
            // writes through it, so the sharing is safe and remains one edit away from not being.
            None => out.push_str("0,0"),
            Some(counts) => out.push_str(&encode(counts)),
        }
    }
    out
}

/// `StrandBiasUtils.encode`: the values joined with commas, and nothing at all when there are none.
///
/// The Java appends a comma after every value but the last; `.join(",")` writes the same string.
/// `.collect::<Vec<_>>()` gathers the rendered values into a list, `_` letting the compiler infer
/// the element type.
pub fn encode(values: &[i32]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// `AS_StrandBiasTest.annotateRawData`, which is `StrandBiasUtils.computeSBAnnotation`: the
/// `AS_SB_TABLE` key and string, or `None` for the Java's empty map. The same for both members.
pub fn annotate_raw_data(
    annotation: AsStrandBias,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Option<(String, String)> {
    // "for allele-specific annotations we only call from HC and we only use likelihoods".
    // `likelihoods?;` returns `None` from this function when there is no matrix, and otherwise
    // does nothing: `?` on an `Option` is the early exit for null.
    likelihoods?;
    let counts = strand_counts(vc, likelihoods, MIN_COUNT);
    Some((
        annotation.raw_key().to_string(),
        make_raw_annotation_string(&counts),
    ))
}

/// `StrandBiasTest.annotate` as the allele-specific members inherit it, which is **not** what
/// their non-allele-specific siblings do.
///
/// Two differences, both invisible in the subclasses:
///
/// - the minimum count is `AS_StrandBiasTest.MIN_COUNT`, which is 2 for **both** members, while the
///   plain `StrandOddsRatio` uses 0. So a sample with one or two informative reads contributes to
///   `SOR` and not to `AS_SOR`, and the two annotations on such a site disagree;
/// - `calculateAnnotationFromGTfield` is overridden to return nothing, because "Allele-specific
///   annotations cannot be called from walkers other than HaplotypeCaller". A site whose genotypes
///   carry `SB` therefore gets `SOR` from that field and no `AS_SOR` at all.
///
/// And the value is still the **pooled** table over every alternate, so this path is no more
/// allele-specific than the rank sums' is.
///
/// The steps are `StrandBiasTest.annotate`'s, in its order: not a variant, empty; any genotype
/// with `SB`, the overridden empty map; no matrix, empty; otherwise the members'
/// `calculateAnnotationFromLikelihoods`. `AS_FisherStrand`'s is
/// `FisherStrand.makeValueObjectForAnnotation(pValue)`, which is
/// `String.format("%.3f", QualityUtils.phredScaleErrorRate(Math.max(pValue, MIN_PVALUE)))`;
/// `AS_StrandOddsRatio`'s is `StrandOddsRatio.formattedValue(ratio)`, `String.format("%.3f", sor)`.
pub fn annotate_direct(
    annotation: AsStrandBias,
    vc: &VariantContext,
    likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
) -> Vec<(String, AnnotationValue)> {
    if !vc.is_variant() {
        return Vec::new();
    }
    // `g.hasAnyAttribute(GATKVCFConstants.STRAND_BIAS_BY_SAMPLE_KEY)` on any genotype. `.any(..)`
    // stops at the first genotype the closure accepts, as the Java's loop returns at the first.
    if vc.genotypes.iter().any(|g| {
        g.extended
            .iter()
            .any(|(key, _)| key == STRAND_BIAS_BY_SAMPLE_KEY)
    }) {
        // The overridden `calculateAnnotationFromGTfield`, which is an empty map.
        return Vec::new();
    }
    let Some(likelihoods) = likelihoods else {
        return Vec::new();
    };
    // `StrandBiasTest.getContingencyTable(likelihoods, vc, MIN_COUNT)`: the 2x2 table pooled over
    // every alternate, with this family's threshold of 2. `java_max` is `Math.max`, kept as a
    // named function so its `NaN` and signed-zero rules are Java's.
    let table = contingency_table(likelihoods, vc, MIN_COUNT);
    let value = match annotation {
        AsStrandBias::Fisher => format_decimals(
            phred_scale_error_rate(gatk_engine::math_utils::java_max(
                p_value_for_contingency_table(table),
                MIN_PVALUE,
            )),
            3,
        ),
        AsStrandBias::OddsRatio => format_decimals(calculate_sor(table), 3),
    };
    vec![(
        annotation.vcf_key().to_string(),
        AnnotationValue::Str(value),
    )]
}

/// `AS_StrandBiasTest.parseRawDataString`.
///
/// An entry may hold **any** number of comma-separated integers, or none; only the count of
/// pipe-separated entries is checked, and it must equal the allele count exactly.
///
/// Java's `parseRawDataString` fills the `AlleleSpecificAnnotationData` it is handed; the port
/// returns the per-allele lists instead, in allele order, every allele present.
pub fn parse_raw_data_string(
    alleles: &[Allele],
    raw: &str,
) -> Result<Vec<(Allele, Vec<i32>)>, AsStrandBiasError> {
    // `getAlleleLengthListOfStringFromRawData`: a **leading** bracket triggers stripping both ends
    // and removing every whitespace character, and split keeps trailing empties.
    //
    // ```java
    // rawDataString = rawDataString.substring(1, rawDataString.length() - 1).replaceAll("\\s", "");
    // ```
    //
    // `\s` in a Java regex is the six ASCII whitespace characters; `char::is_whitespace` is
    // Unicode's wider set, which differs only on non-ASCII input. The golden's `bracketed` and
    // `spaced` rows parse to the same `10,8|3,4`.
    let cleaned = if raw.starts_with('[') {
        raw[1..raw.len() - 1]
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>()
    } else {
        raw.to_string()
    };
    // `split(ALLELE_SPECIFIC_SPLIT_REGEX, -1)`: Rust's `split` keeps every empty piece, which is
    // what the `-1` limit asks of Java's.
    let entries: Vec<&str> = cleaned.split(RAW_DELIM).collect();
    if entries.len() != alleles.len() {
        return Err(AsStrandBiasError::AlleleCountMismatch {
            entries: entries.len(),
            alleles: alleles.len(),
        });
    }
    let mut out = Vec::with_capacity(alleles.len());
    // `.zip(..)` walks the alleles and the entries in step, Java's `myData.getAlleles().get(i)`
    // with `values.get(i)`.
    for (allele, entry) in alleles.iter().zip(entries) {
        let mut counts = Vec::new();
        // `values.get(i).split(",")` and `if (!s.isEmpty())`. Java's split drops trailing empty
        // pieces and Rust's keeps them; the emptiness test skips them either way.
        for piece in entry.split(',') {
            if piece.is_empty() {
                continue;
            }
            // `Integer.parseInt(s.trim())`; `.map_err(..)?` turns a parse failure into the
            // `NumberFormatException` variant and rethrows it.
            let value =
                piece
                    .trim()
                    .parse::<i32>()
                    .map_err(|_| AsStrandBiasError::MalformedCount {
                        raw: piece.to_string(),
                    })?;
            counts.push(value);
        }
        out.push((allele.clone(), counts));
    }
    Ok(out)
}

/// `AS_StrandBiasTest.combineRawData`: the tables added componentwise, through
/// `StrandBiasUtils.combineAttributeMap`, then written with `makeRawAnnotationString` over the
/// merged alleles (the golden's `sbcombine` rows: `two-sources` sums to `12,11|4,5`).
///
/// One allele list serves every source here. In the Java each source is an
/// `AlleleSpecificAnnotationData` built over that source's own `vcPair.getNewAlleles()` in
/// `ReferenceConfidenceVariantContextMerger`, so a source with fewer alleles than the merged list
/// parses against its own; here it would be refused as a count mismatch.
///
/// Where `combineAttributeMap` would throw an index exception (`List.get`) on a list shorter than
/// two, the port does not: it skips the addition, or stores an empty list.
pub fn combine_raw_data(
    alleles: &[Allele],
    raw_strings: &[String],
) -> Result<String, AsStrandBiasError> {
    let mut combined: Vec<(Allele, Option<Vec<i32>>)> = alleles
        .iter()
        .map(|allele| (allele.clone(), None))
        .collect();
    for raw in raw_strings {
        let parsed = parse_raw_data_string(alleles, raw)?;
        for (allele, counts) in &parsed {
            let Some(slot) = combined.iter_mut().find(|(a, _)| a == allele) else {
                continue;
            };
            // `Some(existing) if ..` is a match arm with a guard: it is taken only when the slot
            // is present **and** non-empty; every other state falls to `_`, the catch-all arm.
            match &mut slot.1 {
                Some(existing) if !existing.is_empty() => {
                    // `combineAttributeMap` reads index 0 and 1 unconditionally once the combined
                    // entry exists, so an empty combined entry followed by a non-empty source is
                    // an IndexOutOfBounds upstream. No writer produces that order.
                    if counts.len() >= 2 {
                        existing[FORWARD] += counts[FORWARD];
                        existing[REVERSE] += counts[REVERSE];
                    }
                }
                _ => {
                    // An empty source entry stores an **empty list**, which is present-but-empty
                    // and renders as nothing at all rather than as `0,0`.
                    slot.1 = Some(if counts.len() >= 2 {
                        vec![counts[FORWARD], counts[REVERSE]]
                    } else {
                        Vec::new()
                    });
                }
            }
        }
    }
    Ok(make_raw_annotation_string(&combined))
}

/// `QualityUtils.phredScaleErrorRate`, which is `phredScaleLog10ErrorRate` of the error rate's
/// base-ten logarithm.
///
/// ```java
/// // abs is necessary for edge base with errorRateLog10 = 0 producing -0.0 doubles
/// return Math.abs(-10.0 * Math.max(errorRateLog10, MIN_LOG10_SCALED_QUAL));
/// ```
///
/// The `Math.abs` is there because an error rate of exactly one has a logarithm of `+0.0`, and
/// `-10.0 * +0.0` is `-0.0`, which `%.3f` would print as `-0.000`. The golden's `phred 1.0` row is
/// the bits `0`, a positive zero.
///
/// `phredScaleLog10ErrorRate` also validates its argument as a log probability
/// (`Utils.validateArg(MathUtils.isValidLog10Probability(..))`); the port omits the check, which
/// every p-value passes.
pub fn phred_scale_error_rate(error_rate: f64) -> f64 {
    let log10 = jmath::math::log10(error_rate);
    // `MIN_LOG10_SCALED_QUAL = Math.log10(Double.MIN_VALUE)`, about -323.3.
    // `f64::from_bits(1)` is the double whose bit pattern is 1: the smallest subnormal,
    // `Double.MIN_VALUE`.
    let min_log10 = jmath::math::log10(f64::from_bits(1));
    (-10.0 * gatk_engine::math_utils::java_max(log10, min_log10)).abs()
}

/// `AS_FisherStrand.calculateReducedData` and `AS_StrandOddsRatio.calculateReducedData`.
///
/// `None` in the value position is the reference's null, which becomes the missing value. The two
/// members differ only in whether the reference allele is skipped; see the module note.
///
/// An empty list is `annotationMap.put(a, null)`; the golden's `empty-entry` row gives `AS_FS=.`
/// and `AS_SOR=.`. The table is `{refFwd, refRev}, {altFwd, altRev}`. A reference list shorter than
/// two is an index exception (`List.get`) in the Java and an index panic here.
pub fn calculate_reduced_data(
    annotation: AsStrandBias,
    per_allele: &[(Allele, Vec<i32>)],
    reference: &Allele,
) -> Vec<(Allele, Option<f64>)> {
    let reference_counts: Vec<i32> = per_allele
        .iter()
        .find(|(a, _)| a == reference)
        .map(|(_, counts)| counts.clone())
        .unwrap_or_default();
    // `!a.equals(combinedData.getRefAllele(),true)` ignores the reference flag; the port's `!=`
    // compares it. The two agree unless an alternate carries the reference's own bases, which a
    // well-formed variant does not.
    per_allele
        .iter()
        .filter(|(allele, _)| annotation == AsStrandBias::OddsRatio || allele != reference)
        .map(|(allele, counts)| {
            // `return` inside a closure leaves the closure only, giving this allele's pair.
            if counts.is_empty() {
                return (allele.clone(), None);
            }
            let table = [
                [reference_counts[FORWARD], reference_counts[REVERSE]],
                [counts[FORWARD], counts[REVERSE]],
            ];
            let value = match annotation {
                AsStrandBias::Fisher => phred_scale_error_rate(gatk_engine::math_utils::java_max(
                    p_value_for_contingency_table(table),
                    MIN_PVALUE,
                )),
                AsStrandBias::OddsRatio => calculate_sor(table),
            };
            (allele.clone(), Some(value))
        })
        .collect()
}

/// `AS_StrandBiasTest.makeReducedAnnotationString`.
///
/// An alternate the annotation data does not carry gets no value rather than the missing value.
/// That is the opposite of what the rank sums do with the same situation, and it is what the
/// reference's own log line calls out as an error without acting on it.
///
/// ```java
/// if (!annotationString.isEmpty()) {
///     annotationString += REDUCED_DELIM;
/// }
/// if (!perAltsStrandCounts.containsKey(a)) {
///     logger.warn(
///             "ERROR: VC allele not found in annotation alleles -- maybe there was trimming?");
/// ```
///
/// (The `warn` call is one line in `AS_StrandBiasTest.java`, wrapped here.)
///
/// **The port differs here**: the Java writes the delimiter **before** the lookup, so a missing
/// alternate after the first written value still leaves an empty slot (`1.721,,0.000` for a
/// missing middle alternate, a trailing `1.721,` for a missing last one); only missing alternates
/// before any value vanish. The port skips the delimiter too, giving `1.721,0.000` and `1.721`.
/// No golden row has a missing alternate.
pub fn make_reduced_annotation_string(
    alternates: &[Allele],
    per_allele: &[(Allele, Option<f64>)],
) -> String {
    let mut out = String::new();
    for allele in alternates {
        match per_allele.iter().find(|(a, _)| a == allele) {
            // "ERROR: VC allele not found in annotation alleles -- maybe there was trimming?",
            // logged and then no value appended. The Java has already appended the delimiter by
            // this point (see the item note); the port has not.
            None => continue,
            Some((_, value)) => {
                if !out.is_empty() {
                    out.push(REDUCED_DELIM);
                }
                match value {
                    None => out.push_str(MISSING_VALUE),
                    Some(value) => out.push_str(&format_decimals(*value, 3)),
                }
            }
        }
    }
    out
}

/// `AS_StrandBiasTest.finalizeRawData`: the reduced value and the raw table, both.
///
/// The raw string it writes back is built over the **current** variant's alleles while the data was
/// parsed over the original's, so an allele the merge dropped disappears from the table and one it
/// added is written as `0,0`.
///
/// The steps are the Java's: parse over `originalVC.getAlleles()` (refusing a count mismatch),
/// `calculateReducedData`, the reduced string over `vc.getAlternateAlleles()`, and
/// `StrandBiasUtils.makeRawAnnotationString(vc.getAlleles(), ..)`. Unlike the rank sums there is
/// no early exit for an empty result. The four strings are the reduced key and value, then the raw
/// key and value, as one tuple.
pub fn finalize_raw_data(
    annotation: AsStrandBias,
    vc_alleles: &[Allele],
    original_alleles: &[Allele],
    raw: Option<&str>,
) -> Result<Option<(String, String, String, String)>, AsStrandBiasError> {
    let Some(raw) = raw else { return Ok(None) };
    let reference = reference_allele(original_alleles)?.clone();
    let parsed = parse_raw_data_string(original_alleles, raw)?;
    let reduced_values = calculate_reduced_data(annotation, &parsed, &reference);

    let alternates: Vec<Allele> = vc_alleles
        .iter()
        .filter(|a| !a.is_reference())
        .cloned()
        .collect();
    let reduced = make_reduced_annotation_string(&alternates, &reduced_values);

    let raw_out = make_raw_annotation_string_preserving_empties(&parsed, vc_alleles);
    Ok(Some((
        annotation.vcf_key().to_string(),
        reduced,
        annotation.raw_key().to_string(),
        raw_out,
    )))
}

/// `makeRawAnnotationString` over a parsed map, where an entry that exists but is empty renders as
/// the empty string and only an **absent** entry falls back to `0,0`.
///
/// The same function as [`make_raw_annotation_string`] in the Java; here a second copy because the
/// parsed map holds plain lists (every original allele present) rather than `Option`s. The
/// golden's `sbfinal empty-entry` rows write back `10,8|`.
fn make_raw_annotation_string_preserving_empties(
    parsed: &[(Allele, Vec<i32>)],
    vc_alleles: &[Allele],
) -> String {
    let mut out = String::new();
    for allele in vc_alleles {
        if !out.is_empty() {
            out.push(RAW_DELIM);
        }
        match parsed.iter().find(|(a, _)| a == allele) {
            None => out.push_str("0,0"),
            Some((_, counts)) => {
                let text: Vec<String> = counts.iter().map(|c| c.to_string()).collect();
                out.push_str(&text.join(","));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles() -> Vec<Allele> {
        vec![
            Allele::from_str("A", true).expect("an allele"),
            Allele::from_str("C", false).expect("an allele"),
        ]
    }

    #[test]
    fn the_raw_string_has_an_entry_for_the_reference_and_no_leading_delimiter() {
        let alleles = alleles();
        let counts = vec![
            (alleles[0].clone(), Some(vec![10, 8])),
            (alleles[1].clone(), Some(vec![3, 4])),
        ];
        assert_eq!(make_raw_annotation_string(&counts), "10,8|3,4");
    }

    #[test]
    fn a_count_mismatch_is_refused_rather_than_padded() {
        assert_eq!(
            parse_raw_data_string(&alleles(), "10,8"),
            Err(AsStrandBiasError::AlleleCountMismatch {
                entries: 1,
                alleles: 2
            })
        );
    }

    #[test]
    fn the_odds_ratio_keeps_the_reference_and_fisher_drops_it() {
        let alleles = alleles();
        let parsed = vec![
            (alleles[0].clone(), vec![10, 8]),
            (alleles[1].clone(), vec![3, 4]),
        ];
        let fisher = calculate_reduced_data(AsStrandBias::Fisher, &parsed, &alleles[0]);
        let sor = calculate_reduced_data(AsStrandBias::OddsRatio, &parsed, &alleles[0]);
        assert_eq!(fisher.len(), 1);
        // The reference's own row survives, and is a symmetric table against itself.
        assert_eq!(sor.len(), 2);
    }
}
