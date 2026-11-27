//! `AssemblyComplexity`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.AssemblyComplexity` (GATK 4.7.0.0;
//! `AssemblyComplexity.java` is unchanged since 4.6.2.0): how a
//! variant's haplotypes relate to the germline ones.
//!
//! Every haplotype is credited with the reads that best support it, ties broken as the no-argument
//! `bestAllelesBreakingTies` breaks them: a near-tie with the reference haplotype goes to the
//! reference (`AlleleLikelihoods.java:1082`, `a -> a.isReference() ? 1.0 : 0`), and among
//! haplotypes of equal priority the first allele index wins (golden case `ties`, lines 82-91: the
//! read tied between haplotypes 1 and 2 goes to 1, so `HEC` is `4,0`). Then:
//!
//!  * `HEC`: the haplotypes grouped by their events away from the site (each event written as its
//!    start followed by its alternate bases, concatenated), each group's support, descending;
//!  * `HAPCOMP`: per alternate allele, the events away from the site that the most supported
//!    haplotype carrying it does not share with the closest germline haplotype, counted both ways.
//!    The germline haplotypes are the most supported one, and the second if its support is at least
//!    half the first's (integer division), or the reference alone in reference mode;
//!  * `HAPDOM`: per alternate allele, the share of the carriers' support the most supported carrier
//!    holds, or one over the number of haplotypes when no carrier has any.
//!
//! A symbolic allele, or one whose first base is `*`, gets 0 for both. An allele no haplotype
//! carries makes the reference's `findFirst().get()` throw.
//!
//! "Most supported" is a sort by support, descending, then by bases. The reference sorts a
//! `HashMap`'s entries, so two haplotypes with the same support **and** the same bases (a reference
//! and a non-reference copy, or two uniqueness values) would come out in hash order; here they keep
//! the matrix's order. The assembler never produces such a pair.
//!
//! # Golden rows
//!
//! `tests/data/assembly_complexity.txt.gz` holds one `case` per scenario, its haplotypes, its
//! likelihood rows and the three outputs. The rows that pin a rule rather than illustrate it are
//! named on the item they pin.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, ComplexityError>` is "a value or one of the exceptions the reference throws".
//!   A `?` after an expression returns its error to the caller at once, like an exception
//!   propagating out of this method;
//! - `.ok_or(e)` turns an `Option` (a value or nothing; Java's nullable) into a `Result`, with `e`
//!   as the error for nothing. It is how a Java `NullPointerException` is modelled here;
//! - `x as i32` converts between number types and silently truncates if the value does not fit.
//!   Here it narrows htsjdk-rs's 64-bit `start` to the reference's `int`;
//! - haplotypes are referred to by their index in the likelihood matrix (`usize`), and `support`,
//!   `maps` and `haplotypes` are parallel arrays over that index. This stands in for the
//!   reference's `Map<Haplotype, MutableInt>`, so no hashing of haplotypes is involved;
//! - `|&a, &b| ...` is a lambda whose parameters are unpacked from references on the way in.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::context::ReferenceContext;
use gatk_engine::event_map::{Event, EventMap};
use gatk_engine::fragment::Fragment;
use gatk_engine::haplotype::Haplotype;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

use crate::haplotype_filtering::{HaplotypeLikelihoods, JumboInfoAnnotation};
use crate::info_annotation::AnnotationValue;

/// `GATKVCFConstants.HAPLOTYPE_EQUIVALENCE_COUNTS_KEY`.
pub const HAPLOTYPE_EQUIVALENCE_COUNTS_KEY: &str = "HEC";
/// `GATKVCFConstants.HAPLOTYPE_COMPLEXITY_KEY`.
pub const HAPLOTYPE_COMPLEXITY_KEY: &str = "HAPCOMP";
/// `GATKVCFConstants.HAPLOTYPE_DOMINANCE_KEY`.
pub const HAPLOTYPE_DOMINANCE_KEY: &str = "HAPDOM";

/// What `annotate` throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ComplexityError {
    /// `NoSuchElementException("No value present")`: an allele no haplotype carries.
    NoCarrier,
    /// `NullPointerException`: a haplotype without an event map, or a read with no best allele.
    NullPointer,
    /// `IndexOutOfBoundsException`: reference mode over a matrix without a reference haplotype.
    NoReference,
}

/// The three arrays, `Triple<int[], int[], double[]>`: `getLeft()` is `HEC`, `getMiddle()` is
/// `HAPCOMP`, `getRight()` is `HAPDOM` (`AssemblyComplexity.java:50-52`).
#[derive(Debug, Clone, PartialEq)]
pub struct Complexity {
    pub equivalence_counts: Vec<i32>,
    pub edit_distances: Vec<i32>,
    pub dominance: Vec<f64>,
}

/// `AssemblyComplexity`, with its one argument.
#[derive(Debug, Clone, Default)]
pub struct AssemblyComplexity {
    /// `--assembly-complexity-reference-mode`.
    pub germline_mode: bool,
}

/// `annotate(vc, haplotypeLikelihoods, germlineMode)`, the static one.
///
/// `E` is the evidence type (read or fragment); the reference's `? extends Locatable`. The
/// event maps are read from the matrix's own haplotypes, as the reference does.
pub fn annotate<E: Clone + PartialEq>(
    vc: &VariantContext,
    likelihoods: &AlleleLikelihoods<E, Haplotype>,
    germline_mode: bool,
) -> Result<Complexity, ComplexityError> {
    let n = likelihoods.number_of_alleles();
    let haplotypes: Vec<&Haplotype> = (0..n)
        .map(|a| likelihoods.get_allele(a).expect("an allele index"))
        .collect();
    annotate_with_event_maps(vc, likelihoods, &haplotypes, germline_mode)
}

/// [`annotate`], reading each allele's event map from `haplotypes[allele]` rather than from the
/// matrix's own allele. HaplotypeCaller's matrix and its genotyping haplotypes are ONE set of
/// objects in the reference, so the event maps the genotyping built are the ones read here; a port
/// that holds them as two copies hands over the copy that carries them.
pub fn annotate_with_event_maps<E: Clone + PartialEq>(
    vc: &VariantContext,
    likelihoods: &AlleleLikelihoods<E, Haplotype>,
    haplotypes: &[&Haplotype],
    germline_mode: bool,
) -> Result<Complexity, ComplexityError> {
    let n = likelihoods.number_of_alleles();
    // Every haplotype's event map up front. `collect::<Result<_, _>>()` stops at the first `Err`,
    // so one haplotype without a map is the reference's `NullPointerException` at line 69, where
    // `hap.getEventMap().getEvents()` is reached for every haplotype.
    let maps: Vec<&EventMap> = haplotypes
        .iter()
        .map(|h| h.event_map().ok_or(ComplexityError::NullPointer))
        .collect::<Result<_, _>>()?;
    let vc_start = vc.start as i32;

    // Lines 63-64 of the reference:
    //
    //     haplotypeLikelihoods.bestAllelesBreakingTies()
    //     .forEach(bestHaplotype -> haplotypeSupportCounts.get(bestHaplotype.allele).increment());
    //
    // A read with no best allele makes `get(null)` return null, hence `NullPointer`. `None` is
    // passed for the tie-break priorities because the reference calls the no-argument overload.
    // Best-read support for each haplotype.
    let mut support = vec![0i32; n];
    for best in likelihoods.best_alleles_breaking_ties(None) {
        let allele = best.allele.as_ref().ok_or(ComplexityError::NullPointer)?;
        let index = likelihoods
            .index_of_allele(allele)
            .ok_or(ComplexityError::NullPointer)?;
        support[index] += 1;
    }

    // Haplotypes grouped by their events away from the site; the sums, descending.
    let mut groups: Vec<(String, i32)> = Vec::new();
    for (a, map) in maps.iter().enumerate() {
        let key: String = map
            .events()
            .filter(|e| e.start() != vc_start)
            .map(|e| format!("{}{}", e.start(), e.alt_allele().base_string()))
            .collect();
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some((_, sum)) => *sum += support[a],
            None => groups.push((key, support[a])),
        }
    }
    // The reference's `groupingBy` builds a map with no order guarantee. That does not matter,
    // because only the sums survive and they are sorted. `b.cmp(a)` sorts descending,
    // which is `Comparator.reverseOrder()`.
    let mut equivalence_counts: Vec<i32> = groups.into_iter().map(|(_, sum)| sum).collect();
    equivalence_counts.sort_unstable_by(|a, b| b.cmp(a));

    // Haplotypes by descending support, then by bases.
    // Lines 86-89 of the reference:
    //
    //     Comparator.comparingInt(entry -> -entry.getValue().intValue());
    //     Comparator.comparing(entry -> entry.getKey().getBaseString());
    //     .sorted(countComparator.thenComparing(alleleComparator))
    //
    // `getBaseString()` compares as `String.compareTo`; on the ASCII bases of a haplotype that is
    // the byte order `[u8]::cmp` gives. `sort_by` is stable, so the matrix order breaks what both
    // keys leave tied; the reference's order there is the `HashMap`'s (see the module note).
    let mut by_support: Vec<usize> = (0..n).collect();
    by_support.sort_by(|&a, &b| {
        support[b]
            .cmp(&support[a])
            .then_with(|| haplotypes[a].bases().cmp(&haplotypes[b].bases()))
    });

    // Reference mode takes `getAllele(indexOfReference())`; with no reference haplotype
    // `indexOfReference()` returns -1 (`AlleleList.java:127`) and `getAllele(-1)` throws, which
    // `NoReference` stands for. Otherwise line 99:
    //
    //     haplotypeSupportCounts.get(haplotypesByDescendingSupport.get(1)).intValue()
    //         >= haplotypeSupportCounts.get(haplotypesByDescendingSupport.get(0)).intValue() / 2
    //
    // (one line in the reference, wrapped here)
    // `/` on two `i32`s truncates, as Java's `int` division does. Golden
    // `assembly_complexity.txt.gz` case `germline-het-odd-half` (lines 64-73) pins it: supports
    // 5, 2 and 1, so 2 >= 5 / 2 = 2 admits the second haplotype and `HAPCOMP` is 0 (line 72);
    // a floating-point half, 2.5, would have excluded it and given 1.
    let germline: Vec<usize> = if germline_mode {
        vec![likelihoods
            .index_of_reference()
            .ok_or(ComplexityError::NoReference)?]
    } else {
        let mut germline = vec![by_support[0]];
        if by_support.len() > 1 && support[by_support[1]] >= support[by_support[0]] / 2 {
            germline.push(by_support[1]);
        }
        germline
    };

    let alternates = vc.alternate_alleles();
    // Line 105 (one line in the reference, wrapped here):
    //
    //     if (vc.getAlternateAllele(altAlleleIndex).isSymbolic()
    //         || vc.getAlternateAllele(altAlleleIndex).getBases()[0] == '*') {
    //
    // Golden lines 158-159 (`star-allele`) and 166-167 (`non-ref-allele`) show 0 for such an
    // allele in both `HAPCOMP` and `HAPDOM`.
    let skipped = |allele: &Allele| {
        allele.is_symbolic() || allele.display_string().as_bytes().first() == Some(&b'*')
    };
    let mut edit_distances = Vec::with_capacity(alternates.len());
    let mut dominance = Vec::with_capacity(alternates.len());
    for alt in alternates {
        if skipped(alt) {
            edit_distances.push(0);
            continue;
        }
        // `findFirst().get()` on the carriers in support order: none is `NoSuchElementException`
        // (golden line 173, case `allele-not-carried`).
        let carrier = *by_support
            .iter()
            .find(|&&h| contains_alt_allele(maps[h], vc, alt))
            .ok_or(ComplexityError::NoCarrier)?;
        let distance = germline
            .iter()
            .map(|&g| edit_distance(maps[g], maps[carrier], vc_start))
            .min()
            .expect("a germline haplotype");
        edit_distances.push(distance);
    }
    // The reference computes `HAPCOMP` for every alternate before any `HAPDOM`, so an allele
    // without a carrier has already returned `NoCarrier` above; `counts` is never empty here and
    // the `unwrap_or(i32::MIN)` (a default for an empty maximum) is never taken.
    for alt in alternates {
        if skipped(alt) {
            dominance.push(0.0);
            continue;
        }
        let counts: Vec<i32> = by_support
            .iter()
            .filter(|&&h| contains_alt_allele(maps[h], vc, alt))
            .map(|&h| support[h])
            .collect();
        let max = counts.iter().copied().max().unwrap_or(i32::MIN);
        let sum: i32 = counts.iter().sum();
        // Line 127: `maxCount == 0 ? (1 / (double) haplotypesByDescendingSupport.size()) :
        // maxCount / (double) MathUtils.sum(counts)`. The `max == 0` branch is golden case
        // `unsupported-carrier`: one carrier with no reads, two haplotypes, `HAPDOM` 0.5
        // (line 123).
        // `f64::from` is a lossless widening, Java's implicit `int` to `double`.
        dominance.push(if max == 0 {
            1.0 / by_support.len() as f64
        } else {
            f64::from(max) / f64::from(sum)
        });
    }
    Ok(Complexity {
        equivalence_counts,
        edit_distances,
        dominance,
    })
}

/// `containsAltAllele`: the first event overlapping the site starts there, and its alternate is the
/// context's with as many trailing bases dropped as the context's reference is longer.
///
/// The reference's own comment gives the case this serves (`AssemblyComplexity.java:141-143`): a
/// merged context can write a haplotype's `CAC -> C` as `CACAC -> CAC`, so the context's alleles
/// carry a suffix the event map's minimal alleles do not. Golden case `deletion-merged-context`
/// (`ATAGCT` to `ACT,A`, lines 142-151) has both alternates found, `HAPCOMP` `0,0`.
fn contains_alt_allele(map: &EventMap, vc: &VariantContext, alt: &Allele) -> bool {
    let overlapping = map.overlapping_events(vc.start as i32);
    // `let Some(x) = .. else { .. }` unwraps or takes the early exit: here the reference's
    // `overlapping.isEmpty()` test followed by `overlapping.get(0)`.
    let Some(event) = overlapping.first() else {
        return false;
    };
    if event.start() != vc.start as i32 {
        return false;
    }
    // `excessBases` is an `int` and can be negative; `i64` holds the same values without the
    // unsigned `usize` lengths wrapping on the subtraction.
    let excess = vc.reference().len() as i64 - event.ref_allele().len() as i64;
    equal_bases_excluding_suffix(
        event.alt_allele().display_string().as_bytes(),
        alt.display_string().as_bytes(),
        excess,
    )
}

/// `equalBasesExcludingSuffix`: the event's alternate is a prefix of the context's, exactly
/// `suffix` bases shorter.
fn equal_bases_excluding_suffix(event_bases: &[u8], vc_bases: &[u8], suffix: i64) -> bool {
    if event_bases.len() as i64 + suffix != vc_bases.len() as i64 {
        return false;
    }
    if event_bases.len() > vc_bases.len() {
        // The event map's allele is longer, though minimal, because it is an MNP.
        return false;
    }
    // The reference's byte-by-byte loop. `&vc_bases[..k]` is a view of the first `k` bytes, and
    // `==` on two slices compares length and contents, like `Arrays.equals`.
    event_bases == &vc_bases[..event_bases.len()]
}

/// `uniqueVariants`: the events of the first map, away from the excluded position, that the second
/// does not hold at the same start.
///
/// `AssemblyComplexity.java:186`, where a missing event is `null` and never equal:
///
/// ```java
/// .filter(event -> !event.equals(eventMap2.get(event.getStart())))
/// ```
///
/// `second.get(..)` is `Some(event)` or `None`, so comparing it to `Some(e)` is that same test
/// with `None` playing `null`. `.count()` returns a `usize`; the `as i32` is the reference's
/// `(int)` cast of a `long` count.
fn unique_variants(first: &EventMap, second: &EventMap, excluded: i32) -> i32 {
    first
        .events()
        .filter(|e| e.start() != excluded)
        .filter(|e| second.get(e.start()) != Some(*e as &Event))
        .count() as i32
}

/// `editDistance`: the events either haplotype has and the other lacks, away from the site.
fn edit_distance(first: &EventMap, second: &EventMap, excluded: i32) -> i32 {
    unique_variants(first, second, excluded) + unique_variants(second, first, excluded)
}

impl JumboInfoAnnotation for AssemblyComplexity {
    fn key_names(&self) -> Vec<&'static str> {
        vec![
            HAPLOTYPE_EQUIVALENCE_COUNTS_KEY,
            HAPLOTYPE_COMPLEXITY_KEY,
            HAPLOTYPE_DOMINANCE_KEY,
        ]
    }

    /// The three arrays under their keys. A refusal is the reference's exception, which the
    /// annotation engine does not catch; it is reported here as no annotation.
    ///
    /// The fragment-grouped matrix is used when present, the read-grouped one otherwise (see
    /// `HaplotypeLikelihoods`); the reference receives whichever its caller passed as
    /// `haplotypeLikelihoods`.
    fn annotate(
        &self,
        _reference: Option<&ReferenceContext>,
        vc: &VariantContext,
        _likelihoods: Option<&AlleleLikelihoods<BamRecord>>,
        _fragment_likelihoods: Option<&AlleleLikelihoods<Fragment>>,
        haplotype_likelihoods: &HaplotypeLikelihoods<'_>,
    ) -> Vec<(String, AnnotationValue)> {
        let result = match haplotype_likelihoods {
            HaplotypeLikelihoods::ByFragment(l) => annotate(vc, l, self.germline_mode),
            HaplotypeLikelihoods::ByRead(l) => annotate(vc, l, self.germline_mode),
        };
        let Ok(c) = result else {
            return Vec::new();
        };
        let ints =
            |v: &[i32]| AnnotationValue::List(v.iter().map(|&x| AnnotationValue::Int(x)).collect());
        vec![
            (
                HAPLOTYPE_EQUIVALENCE_COUNTS_KEY.to_string(),
                ints(&c.equivalence_counts),
            ),
            (
                HAPLOTYPE_COMPLEXITY_KEY.to_string(),
                ints(&c.edit_distances),
            ),
            (
                HAPLOTYPE_DOMINANCE_KEY.to_string(),
                AnnotationValue::List(
                    c.dominance
                        .iter()
                        .map(|&x| AnnotationValue::Double(x))
                        .collect(),
                ),
            ),
        ]
    }
}
