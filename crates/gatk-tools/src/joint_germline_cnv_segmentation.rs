//! `JointGermlineCNVSegmentation`: how single-sample gCNV calls become one joint call set.
//!
//! Two engines run in sequence. A defragmenter joins ONE sample's adjacent segments, padding each
//! by a fraction of its own length; a max-clique cluster then joins different samples' events into
//! one site. More than one input sample skips the defragmenter entirely, the input being assumed
//! pre-clustered, so the two are never both visible in the same run.
//!
//! This module is the golden's model of the rules. The tool itself, with the streaming engines of
//! [`crate::sv_cluster_engine`], the VCFs read and written and the genotypes squared off, is the
//! `JointGermlineCNVSegmentation` runner in `gatk-cli`.
//!
//! Reference: GATK 4.7.0.0, `tools/walkers/sv/JointGermlineCNVSegmentation.java` (unchanged since
//! 4.6.2.0) and `tools/sv/cluster/CNVLinkage.java`.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<i32>` is "an `int` or nothing", how an absent FORMAT attribute is modelled; in
//!   `alleles`, `None` is `Allele.NO_CALL`. `.unwrap_or(0)` is `getAttributeAsInt(g, key, 0)`.
//! - `Result<i32, SegmentationError>` is "the value or the exception the reference would throw".
//! - `let Some(x) = y else { ... };` unwraps `y` or takes the early exit in the braces.
//! - `matches!(value, pattern)` is a boolean test of a value's shape, and a `match` arm may join
//!   several shapes with `|`.
//! - `x as i32` from a double truncates toward zero, Java's `(int)` cast.

use crate::sv_stratify::SvType;

/// `PedigreeValidationType.STRICT`'s only refusal, and the two the genotype rules make.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SegmentationError {
    /// A sample in the VCFs and not in the pedigree.
    MissingFromPedigree { sample: String },
    /// A genotype carrying more alleles than the ploidy it was given.
    PloidyMismatch { ploidy: usize, alleles: usize },
    /// An allosomal contig this tool has no sex rule for.
    UnknownAllosomalContig { contig: String },
}

impl SegmentationError {
    /// The reference's exception message for each refusal, copied from `SampleDBBuilder.java:103`
    /// and `JointGermlineCNVSegmentation.java:586-587` and `:747-748`.
    pub fn message(&self) -> String {
        match self {
            SegmentationError::MissingFromPedigree { sample } => format!(
                "Sample {sample} found in data sources but not in pedigree files with STRICT \
                 pedigree validation"
            ),
            SegmentationError::PloidyMismatch { ploidy, alleles } => {
                format!("Encountered genotype with ploidy {ploidy} but {alleles} alleles.")
            }
            SegmentationError::UnknownAllosomalContig { contig } => format!(
                "Encountered unknown allosomal contig: {contig}. This tool only supports \
                 mammalian genomes with XX/XY sex determination."
            ),
        }
    }
}

/// The contigs whose ploidy comes from the pedigree rather than from an argument:
/// `Arrays.asList("X","Y","chrX","chrY")` (`JointGermlineCNVSegmentation.java:110`).
pub const ALLOSOMAL_CONTIGS: &[&str] = &["X", "Y", "chrX", "chrY"];

/// `Sex`, as the pedigree's fifth column spells it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    Unknown,
}

impl Sex {
    /// The pedigree column: 1 male, 2 female, anything else unknown.
    pub fn parse(code: &str) -> Sex {
        match code {
            "1" => Sex::Male,
            "2" => Sex::Female,
            _ => Sex::Unknown,
        }
    }
}

/// The pedigree, reduced to the one column the ploidy rule reads.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Pedigree {
    pub samples: Vec<(String, Sex)>,
}

impl Pedigree {
    /// `sampleDB.getSample(sampleName).getSex()`, or `None` for a sample with no entry.
    pub fn sex(&self, sample: &str) -> Option<Sex> {
        self.samples
            .iter()
            .find(|(name, _)| name == sample)
            .map(|(_, sex)| *sex)
    }

    /// `SampleDB.createSampleDBFromPedigreeAndDataSources` with STRICT validation, which refuses
    /// before a record is read.
    pub fn validate(&self, samples: &[String]) -> Result<(), SegmentationError> {
        for sample in samples {
            if self.sex(sample).is_none() {
                return Err(SegmentationError::MissingFromPedigree {
                    sample: sample.clone(),
                });
            }
        }
        Ok(())
    }
}

/// One sample's call on one segment. `None` in `alleles` is a no-call.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    pub sample: String,
    pub alleles: Vec<Option<i32>>,
    pub copy_number: Option<i32>,
    /// `QS`, the quality the entry filter reads.
    pub quality_some: Option<i32>,
    /// `ECN`, when the input already carried one.
    pub expected_copy_number: Option<i32>,
}

impl Genotype {
    /// htsjdk's `isHomRef`: at least one allele, and every one the reference.
    pub fn is_hom_ref(&self) -> bool {
        !self.alleles.is_empty() && self.alleles.iter().all(|allele| *allele == Some(0))
    }

    /// htsjdk's `isNoCall`: at least one allele, and every one a no-call. A genotype mixing
    /// called and uncalled alleles is neither this nor hom ref.
    pub fn is_no_call(&self) -> bool {
        !self.alleles.is_empty() && self.alleles.iter().all(Option::is_none)
    }

    pub fn ploidy(&self) -> usize {
        self.alleles.len()
    }

    /// `isNullCall`: a no-call whose copy number is exactly zero, which is a call on a contig the
    /// model had nothing for.
    pub fn is_null_call(&self) -> bool {
        self.copy_number == Some(0) && self.is_no_call()
    }
}

/// One gCNV segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub sv_type: SvType,
    pub genotypes: Vec<Genotype>,
}

impl Segment {
    /// The closed interval's length, `end - start + 1`.
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }

    /// The samples carrying a non-reference allele: a called allele whose index is above 0.
    ///
    /// This is the called-genotype branch of `SVCallRecord.isCarrier` only. The reference also
    /// counts a no-call genotype as a carrier when its copy number departs from its expected copy
    /// number (`SVCallRecord.java:281-290`), which this model does not.
    pub fn carriers(&self) -> Vec<String> {
        self.genotypes
            .iter()
            .filter(|genotype| {
                genotype
                    .alleles
                    .iter()
                    .any(|allele| matches!(allele, Some(index) if *index > 0))
            })
            .map(|genotype| genotype.sample.clone())
            .collect()
    }
}

/// `createDepthOnlyFromGCNVWithOriginalGenotypes`'s entry test, which applies ONLY to a
/// single-genotype record.
///
/// Four separate reasons drop one, and they are asked in this order: a hom-ref genotype, a no-call
/// with no copy number at all, a quality strictly below the threshold, and a null call. A record
/// carrying more than one genotype is never tested.
pub fn keeps_record(segment: &Segment, min_quality: i32) -> bool {
    if segment.genotypes.len() != 1 {
        return true;
    }
    let genotype = &segment.genotypes[0];
    if genotype.is_hom_ref() {
        return false;
    }
    if genotype.is_no_call() && genotype.copy_number.is_none() {
        return false;
    }
    // Absent QS reads as 0, and the comparison is strictly less than, so a call exactly at the
    // threshold survives.
    if genotype.quality_some.unwrap_or(0) < min_quality {
        return false;
    }
    !genotype.is_null_call()
}

/// `getSamplePloidy` (`JointGermlineCNVSegmentation.java:553-590`).
///
/// An `ECN` already on the genotype is read FIRST, before the contig is even looked at. That value
/// does not reach the output, though: the ploidy is derived again for every sample when the site is
/// written, and by then the genotype no longer carries the input's own.
pub fn sample_ploidy(
    reference_autosomal_copy_number: i32,
    pedigree: Option<&Pedigree>,
    sample: &str,
    contig: &str,
    genotype: Option<&Genotype>,
) -> Result<i32, SegmentationError> {
    if let Some(expected) = genotype.and_then(|g| g.expected_copy_number) {
        return Ok(expected);
    }
    if !ALLOSOMAL_CONTIGS.contains(&contig) {
        return Ok(reference_autosomal_copy_number);
    }
    let sex = pedigree.and_then(|pedigree| pedigree.sex(sample));
    let Some(sex) = sex else {
        // No pedigree entry: the genotype's own ploidy stands in.
        return match genotype {
            Some(genotype) => Ok(genotype.ploidy() as i32),
            None => Err(SegmentationError::MissingFromPedigree {
                sample: sample.to_string(),
            }),
        };
    };
    match contig {
        "X" | "chrX" => Ok(match sex {
            Sex::Female => 2,
            // An unknown sex is given the MALE answer, not the female one.
            Sex::Male | Sex::Unknown => 1,
        }),
        "Y" | "chrY" => Ok(match sex {
            Sex::Female => 0,
            Sex::Male | Sex::Unknown => 1,
        }),
        other => Err(SegmentationError::UnknownAllosomalContig {
            contig: other.to_string(),
        }),
    }
}

/// `correctGenotypePloidy` (`JointGermlineCNVSegmentation.java:740-755`): the alleles are padded
/// to the ploidy with reference alleles, and a genotype with more alleles than that is refused.
///
/// A SINGLE no-call allele is a special case: it becomes a no-call of the full ploidy rather than
/// one no-call padded with reference alleles.
pub fn correct_genotype_ploidy(
    genotype: &Genotype,
    ploidy: i32,
) -> Result<Vec<Option<i32>>, SegmentationError> {
    let ploidy = ploidy.max(0) as usize;
    if genotype.alleles.len() == 1 && genotype.alleles[0].is_none() {
        return Ok(vec![None; ploidy]);
    }
    if genotype.alleles.len() > ploidy {
        return Err(SegmentationError::PloidyMismatch {
            ploidy,
            alleles: genotype.alleles.len(),
        });
    }
    let mut alleles = genotype.alleles.clone();
    while alleles.len() < ploidy {
        alleles.push(Some(0));
    }
    Ok(alleles)
}

/// `CNVLinkage.DEFAULT_PADDING_FRACTION` and `DEFAULT_SAMPLE_OVERLAP` (`CNVLinkage.java:27-28`).
pub const DEFAULT_PADDING_FRACTION: f64 = 0.25;
pub const DEFAULT_SAMPLE_OVERLAP: f64 = 0.8;

/// `getPaddedRecordInterval`: the padding is a fraction of the record's OWN length, truncated
/// toward zero, and it is clipped to the contig (`CNVLinkage.java:115-118`):
///
/// ```java
/// return new SimpleInterval(contig, start, end)
///         .expandWithinContig((int) (paddingFraction * (end - start + 1)), dictionary);
/// ```
pub fn padded_interval(
    start: i32,
    end: i32,
    padding_fraction: f64,
    contig_length: i32,
) -> (i32, i32) {
    let padding = (padding_fraction * (end - start + 1) as f64) as i32;
    ((start - padding).max(1), (end + padding).min(contig_length))
}

/// `hasSampleOverlap`: a fraction of the SMALLER carrier set, and zero asks nothing.
///
/// Not what GATK 4.7.0.0 computes. `SVClusterLinkage.getSampleSetOverlap` divides by the LARGER
/// set and passes two empty sets (`SVClusterLinkage.java:51-57`):
///
/// ```java
/// final double denom = Math.max(a.size(), b.size());
/// if (denom == 0) {
///     return 1;
/// }
/// ```
///
/// and for a record of type `CNV`, `computeSampleOverlap` compares copy states over every sample
/// rather than carrier sets (`:71-93`).
pub fn has_sample_overlap(a: &Segment, b: &Segment, threshold: f64) -> bool {
    if threshold <= 0.0 {
        return true;
    }
    let carriers_a = a.carriers();
    let carriers_b = b.carriers();
    let shared = carriers_a
        .iter()
        .filter(|sample| carriers_b.contains(sample))
        .count();
    let smaller = carriers_a.len().min(carriers_b.len());
    if smaller == 0 {
        return false;
    }
    shared as f64 / smaller as f64 >= threshold
}

/// `CNVLinkage.areClusterable` (`CNVLinkage.java:43-95`): type, padded overlap, sample overlap,
/// then the singleton rule, each a refusal in that order.
///
/// The last test is the one the class documentation calls out: a SINGLETON is only joined to
/// another record of the same copy state, compared as the difference from the ploidy rather than
/// as the copy number itself, so a haploid CN of 0 and a diploid CN of 1 are both one lost copy
/// and still do not match.
pub fn are_clusterable(
    a: &Segment,
    b: &Segment,
    padding_fraction: f64,
    min_sample_overlap: f64,
    contig_length: i32,
) -> bool {
    if !matches!(a.sv_type, SvType::Del | SvType::Dup | SvType::Cnv)
        || !matches!(b.sv_type, SvType::Del | SvType::Dup | SvType::Cnv)
    {
        return false;
    }
    if a.contig != b.contig || a.sv_type != b.sv_type {
        return false;
    }
    let (left_a, right_a) = padded_interval(a.start, a.end, padding_fraction, contig_length);
    let (left_b, right_b) = padded_interval(b.start, b.end, padding_fraction, contig_length);
    if left_a > right_b || left_b > right_a {
        return false;
    }
    if !has_sample_overlap(a, b, min_sample_overlap) {
        return false;
    }

    // `carriersA.size() == 1 && carriersA.equals(carriersB)`.
    let carriers_a = a.carriers();
    let carriers_b = b.carriers();
    if carriers_a.len() == 1 && carriers_a == carriers_b {
        let sample = &carriers_a[0];
        let (Some(genotype_a), Some(genotype_b)) = (
            a.genotypes.iter().find(|g| g.sample == *sample),
            b.genotypes.iter().find(|g| g.sample == *sample),
        ) else {
            return true;
        };
        // Both carry `CN`: compare the copy-number deltas. Otherwise compare the allele lists,
        // sorted (`SVCallRecordUtils.sortAlleles`), which only equality reads here.
        match (genotype_a.copy_number, genotype_b.copy_number) {
            (Some(copy_a), Some(copy_b)) => {
                let delta_a = genotype_a.ploidy() as i32 - copy_a;
                let delta_b = genotype_b.ploidy() as i32 - copy_b;
                if delta_a != delta_b {
                    return false;
                }
            }
            _ => {
                let mut sorted_a = genotype_a.alleles.clone();
                let mut sorted_b = genotype_b.alleles.clone();
                sorted_a.sort();
                sorted_b.sort();
                if sorted_a != sorted_b {
                    return false;
                }
            }
        }
    }
    true
}

/// One defragmented segment: where it ended up and which inputs went into it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Defragmented {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub members: Vec<usize>,
}

/// The defragmenter over one sample's records, in input order: single linkage, with each new
/// record tested against every MEMBER of a cluster rather than against the cluster's extent.
///
/// That distinction is the whole of it. The padding is a fraction of the record's own length, so a
/// long record reaches much further than a short one, and a run of joined records reaches as far
/// as its LONGEST member does rather than as far as the joined span would. Testing against the
/// grown extent instead absorbs a record the reference leaves alone: at the default fraction the
/// span 40000-80100 would reach past 90000, but neither of the two records that made it does.
pub fn defragment(
    segments: &[Segment],
    padding_fraction: f64,
    min_sample_overlap: f64,
    contig_length: i32,
) -> Vec<Defragmented> {
    // Each cluster is the list of indices of its member segments.
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for (index, segment) in segments.iter().enumerate() {
        let linked: Vec<usize> = clusters
            .iter()
            .enumerate()
            .filter(|(_, members)| {
                members.iter().any(|member| {
                    are_clusterable(
                        &segments[*member],
                        segment,
                        padding_fraction,
                        min_sample_overlap,
                        contig_length,
                    )
                })
            })
            .map(|(at, _)| at)
            .collect();
        // `split_first` is `None` for an empty list, or `Some((first, rest))`.
        match linked.split_first() {
            None => clusters.push(vec![index]),
            Some((first, rest)) => {
                // Single linkage: a record that reaches two clusters merges them.
                let mut merged = clusters[*first].clone();
                for other in rest {
                    merged.extend(clusters[*other].iter().copied());
                }
                merged.push(index);
                merged.sort();
                for at in linked.iter().rev() {
                    clusters.remove(*at);
                }
                clusters.push(merged);
            }
        }
    }
    // A cluster's extent is its members' smallest start and largest end; the clusters come out
    // ordered by start.
    let mut out: Vec<Defragmented> = clusters
        .into_iter()
        .map(|members| Defragmented {
            contig: segments[members[0]].contig.clone(),
            start: members
                .iter()
                .map(|member| segments[*member].start)
                .min()
                .expect("a member"),
            end: members
                .iter()
                .map(|member| segments[*member].end)
                .max()
                .expect("a member"),
            members,
        })
        .collect();
    out.sort_by_key(|record| record.start);
    out
}

/// More than one sample in the inputs means the defragmenter never runs.
pub fn is_multi_sample(samples: &[String]) -> bool {
    samples.len() != 1
}
