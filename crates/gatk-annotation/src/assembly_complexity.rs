//! `AssemblyComplexity`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.AssemblyComplexity` (GATK 4.6.2.0): how a
//! variant's haplotypes relate to the germline ones.
//!
//! Every haplotype is credited with the reads that best support it, ties broken by the matrix's own
//! rule (the first allele index). Then:
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

/// The three arrays, `Triple<int[], int[], double[]>`.
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

/// `annotate(vc, haplotypeLikelihoods, germlineMode)`.
pub fn annotate<E: Clone + PartialEq>(
    vc: &VariantContext,
    likelihoods: &AlleleLikelihoods<E, Haplotype>,
    germline_mode: bool,
) -> Result<Complexity, ComplexityError> {
    let n = likelihoods.number_of_alleles();
    let haplotypes: Vec<&Haplotype> = (0..n)
        .map(|a| likelihoods.get_allele(a).expect("an allele index"))
        .collect();
    let maps: Vec<&EventMap> = haplotypes
        .iter()
        .map(|h| h.event_map().ok_or(ComplexityError::NullPointer))
        .collect::<Result<_, _>>()?;
    let vc_start = vc.start as i32;

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
    let mut equivalence_counts: Vec<i32> = groups.into_iter().map(|(_, sum)| sum).collect();
    equivalence_counts.sort_unstable_by(|a, b| b.cmp(a));

    // Haplotypes by descending support, then by bases.
    let mut by_support: Vec<usize> = (0..n).collect();
    by_support.sort_by(|&a, &b| {
        support[b]
            .cmp(&support[a])
            .then_with(|| haplotypes[a].bases().cmp(&haplotypes[b].bases()))
    });

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
fn contains_alt_allele(map: &EventMap, vc: &VariantContext, alt: &Allele) -> bool {
    let overlapping = map.overlapping_events(vc.start as i32);
    let Some(event) = overlapping.first() else {
        return false;
    };
    if event.start() != vc.start as i32 {
        return false;
    }
    let excess = vc.reference().len() as i64 - event.ref_allele().len() as i64;
    equal_bases_excluding_suffix(
        event.alt_allele().display_string().as_bytes(),
        alt.display_string().as_bytes(),
        excess,
    )
}

/// `equalBasesExcludingSuffix`.
fn equal_bases_excluding_suffix(event_bases: &[u8], vc_bases: &[u8], suffix: i64) -> bool {
    if event_bases.len() as i64 + suffix != vc_bases.len() as i64 {
        return false;
    }
    if event_bases.len() > vc_bases.len() {
        // The event map's allele is longer, though minimal, because it is an MNP.
        return false;
    }
    event_bases == &vc_bases[..event_bases.len()]
}

/// `uniqueVariants`: the events of the first map, away from the excluded position, that the second
/// does not hold at the same start.
fn unique_variants(first: &EventMap, second: &EventMap, excluded: i32) -> i32 {
    first
        .events()
        .filter(|e| e.start() != excluded)
        .filter(|e| second.get(e.start()) != Some(*e as &Event))
        .count() as i32
}

/// `editDistance`.
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
