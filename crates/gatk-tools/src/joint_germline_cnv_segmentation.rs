//! `JointGermlineCNVSegmentation`: how single-sample gCNV calls become one joint call set.
//!
//! Two engines run in sequence. A defragmenter joins ONE sample's adjacent segments, padding each
//! by a fraction of its own length; a max-clique cluster then joins different samples' events into
//! one site. More than one input sample skips the defragmenter entirely, the input being assumed
//! pre-clustered, so the two are never both visible in the same run.
//!
//! Reading and writing the VCFs and the cross-sample clustering, which is
//! [`crate::sv_cluster`]'s, are not ported. The entry filter, the ploidy rules, the genotype
//! padding and the defragmenter's own linkage are.

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

/// The contigs whose ploidy comes from the pedigree rather than from an argument.
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
    pub fn is_hom_ref(&self) -> bool {
        !self.alleles.is_empty() && self.alleles.iter().all(|allele| *allele == Some(0))
    }

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
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }

    /// The samples carrying a non-reference allele.
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

/// `getSamplePloidy`.
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

/// `correctGenotypePloidy`: the alleles are padded to the ploidy with reference alleles, and a
/// genotype with more alleles than that is refused.
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

/// `CNVLinkage.DEFAULT_PADDING_FRACTION` and `DEFAULT_SAMPLE_OVERLAP`.
pub const DEFAULT_PADDING_FRACTION: f64 = 0.25;
pub const DEFAULT_SAMPLE_OVERLAP: f64 = 0.8;

/// `getPaddedRecordInterval`: the padding is a fraction of the record's OWN length, truncated
/// toward zero, and it is clipped to the contig.
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

/// `CNVLinkage.areClusterable`.
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
