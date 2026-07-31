//! The read side of `org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods`
//! (GATK 4.6.2.0): the matrix, the evidence counts, and the best-allele search that every
//! likelihood-reading annotation goes through.
//!
//! The matrix is `[sample][allele][evidence]`, and the two axes are the ones in
//! [`crate::allele_list`]. What is ported here is what a consumer reads; the mutators
//! (`marginalize`, `groupEvidence`, `addMissingAlleles`, the downsamplers) are their own slices.
//!
//! # `secondBestIndex` starts at zero, and is not "the second allele"
//!
//! ```java
//! int bestAlleleIndex = canBeReference || referenceAlleleIndex != 0 ? 0 : 1;
//! int secondBestIndex = 0;
//! ```
//!
//! Both start at index 0, so before any comparison the "second best" is allele 0 whether or not
//! anything is second. The loop only ever moves `secondBestIndex` when a candidate beats the
//! current best or the current second best, so a matrix with one allele ends with
//! `secondBestIndex == bestAlleleIndex`, which the tail then turns into a second-best likelihood of
//! negative infinity. A port that initialised it to "missing" would answer the same thing by
//! accident here and differently as soon as the priority pass runs.
//!
//! # Ties go to the lower index, and the comparison is strict
//!
//! `candidateLikelihood > bestLikelihood` never displaces an equal one, so among equal likelihoods
//! the earliest allele wins. That is what makes the allele **order** observable in a result that
//! looks order-independent.
//!
//! # `isInformative` and the tie-breaking threshold disagree once the log base changes
//!
//! ```java
//! private double getInformativeThreshold() { return isNaturalLog ? NATURAL_LOG_INFORMATIVE_THRESHOLD : LOG_10_INFORMATIVE_THRESHOLD; }
//! public boolean isInformative() { return confidence > LOG_10_INFORMATIVE_THRESHOLD; }
//! ```
//!
//! The search converts its threshold to natural log when the matrix has been switched;
//! `BestAllele.isInformative` does not, and keeps comparing against 0.2 whatever the base. After
//! `switch_to_natural_log`, the two questions "was this tie broken by priority" and "is this call
//! informative" are asked at two different thresholds. Both are transcribed as they are.
//!
//! # `confidence` guards the infinite case rather than the equal case
//!
//! ```java
//! confidence = likelihood == secondBestLikelihood ? 0 : likelihood - secondBestLikelihood;
//! ```
//!
//! The guard reads as "equal likelihoods mean no confidence", and its load-bearing use is the pair
//! of negative infinities that a one-allele matrix produces, where the subtraction would be `NaN`.

use crate::allele_list::{AlleleList, SampleList};
use htsjdk_vcf::allele::Allele;

/// `LOG_10_INFORMATIVE_THRESHOLD`.
pub const LOG_10_INFORMATIVE_THRESHOLD: f64 = 0.2;

/// `MathUtils.LOG_10`, which is `Math.log(10)` and not a decimal literal.
pub const LOG_10: f64 = std::f64::consts::LN_10;

/// `MathUtils.log10ToLog`.
pub fn log10_to_log(log10: f64) -> f64 {
    log10 * LOG_10
}

/// `NATURAL_LOG_INFORMATIVE_THRESHOLD`, which is the constant above put through that conversion
/// rather than written out.
pub fn natural_log_informative_threshold() -> f64 {
    log10_to_log(LOG_10_INFORMATIVE_THRESHOLD)
}

/// `AlleleLikelihoods.BestAllele`.
///
/// The two allele fields are `None` where the reference's are null, which happens only when the
/// matrix has no allele to offer.
#[derive(Debug, Clone, PartialEq)]
pub struct BestAllele {
    pub sample: String,
    pub evidence_index: usize,
    pub allele: Option<Allele>,
    pub second_best_allele: Option<Allele>,
    pub likelihood: f64,
    pub second_best_likelihood: f64,
    pub confidence: f64,
}

impl BestAllele {
    /// `isInformative()`: the confidence against the **log10** threshold, whatever base the matrix
    /// is in. See the module note.
    pub fn is_informative(&self) -> bool {
        self.confidence > LOG_10_INFORMATIVE_THRESHOLD
    }
}

/// What the constructor refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LikelihoodsError {
    /// A sample in the evidence map that the sample list does not hold. The reference reaches an
    /// `IllegalArgumentException` through `Utils.nonNull` on the lookup.
    UnknownSample(String),
    /// A likelihood row whose length is not the sample's evidence count.
    WrongRowLength { sample: String, allele: usize },
    /// `switchToNaturalLog` called twice.
    AlreadyNaturalLog,
}

/// `AlleleLikelihoods`, over an evidence type the caller names.
///
/// The reference's evidence is a `GATKRead`; here it is whatever identifies a piece of evidence to
/// the caller, because nothing in this slice reads the evidence itself. The annotations that do
/// (`CountNs` looks at a base, `MappingQualityZero` at a mapping quality) carry their own type.
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleLikelihoods<E: Clone + PartialEq> {
    samples: SampleList,
    alleles: AlleleList,
    /// `evidenceBySampleIndex`.
    evidence_by_sample: Vec<Vec<E>>,
    /// `valuesBySampleIndex`, indexed `[sample][allele][evidence]`.
    values: Vec<Vec<Vec<f64>>>,
    /// `referenceAlleleIndex`, `None` where the reference holds `MISSING_INDEX`.
    reference_allele_index: Option<usize>,
    is_natural_log: bool,
}

impl<E: Clone + PartialEq> AlleleLikelihoods<E> {
    /// The public constructor, with the likelihood values supplied rather than left at zero.
    ///
    /// `evidence_by_sample` and `values` are given per sample **index**, not per name: the
    /// reference's map is immediately turned into index order by `setupIndexes`, and a caller that
    /// passed a name the sample list does not hold gets a null there.
    pub fn new(
        samples: SampleList,
        alleles: AlleleList,
        evidence_by_sample: Vec<Vec<E>>,
        values: Vec<Vec<Vec<f64>>>,
    ) -> Result<Self, LikelihoodsError> {
        let sample_count = samples.number_of_samples();
        let allele_count = alleles.number_of_alleles();
        if evidence_by_sample.len() != sample_count || values.len() != sample_count {
            return Err(LikelihoodsError::UnknownSample(String::new()));
        }
        for (index, sample_values) in values.iter().enumerate() {
            let evidence_count = evidence_by_sample[index].len();
            if sample_values.len() != allele_count {
                return Err(LikelihoodsError::WrongRowLength {
                    sample: samples
                        .get_sample(index)
                        .cloned()
                        .unwrap_or_else(String::new),
                    allele: sample_values.len(),
                });
            }
            for (allele, row) in sample_values.iter().enumerate() {
                if row.len() != evidence_count {
                    return Err(LikelihoodsError::WrongRowLength {
                        sample: samples
                            .get_sample(index)
                            .cloned()
                            .unwrap_or_else(String::new),
                        allele,
                    });
                }
            }
        }

        let reference_allele_index = alleles.index_of_reference();
        Ok(Self {
            samples,
            alleles,
            evidence_by_sample,
            values,
            reference_allele_index,
            is_natural_log: false,
        })
    }

    pub fn number_of_samples(&self) -> usize {
        self.samples.number_of_samples()
    }

    pub fn number_of_alleles(&self) -> usize {
        self.alleles.number_of_alleles()
    }

    pub fn index_of_sample(&self, sample: &str) -> Option<usize> {
        self.samples.index_of_sample(sample)
    }

    pub fn index_of_allele(&self, allele: &Allele) -> Option<usize> {
        self.alleles.index_of_allele(allele)
    }

    pub fn get_allele(&self, index: usize) -> Option<&Allele> {
        self.alleles.get_allele(index)
    }

    pub fn get_sample(&self, index: usize) -> Option<&String> {
        self.samples.get_sample(index)
    }

    /// `sampleEvidence(sampleIndex)`.
    pub fn sample_evidence(&self, sample_index: usize) -> Option<&[E]> {
        self.evidence_by_sample
            .get(sample_index)
            .map(|list| list.as_slice())
    }

    /// `sampleEvidenceCount(sampleIndex)`.
    pub fn sample_evidence_count(&self, sample_index: usize) -> usize {
        self.evidence_by_sample
            .get(sample_index)
            .map_or(0, |list| list.len())
    }

    /// `evidenceCount()`: the total across samples, which is what `Coverage` writes as `DP`.
    pub fn evidence_count(&self) -> usize {
        self.evidence_by_sample.iter().map(Vec::len).sum()
    }

    pub fn is_natural_log(&self) -> bool {
        self.is_natural_log
    }

    /// `getInformativeThreshold()`: converted with the base, unlike `BestAllele.isInformative`.
    pub fn informative_threshold(&self) -> f64 {
        if self.is_natural_log {
            natural_log_informative_threshold()
        } else {
            LOG_10_INFORMATIVE_THRESHOLD
        }
    }

    /// One likelihood.
    pub fn value(&self, sample_index: usize, allele_index: usize, evidence_index: usize) -> f64 {
        self.values[sample_index][allele_index][evidence_index]
    }

    /// `switchToNaturalLog()`, which refuses to run twice.
    pub fn switch_to_natural_log(&mut self) -> Result<(), LikelihoodsError> {
        if self.is_natural_log {
            return Err(LikelihoodsError::AlreadyNaturalLog);
        }
        for sample_values in &mut self.values {
            for row in sample_values.iter_mut() {
                for value in row.iter_mut() {
                    *value = log10_to_log(*value);
                }
            }
        }
        self.is_natural_log = true;
        Ok(())
    }

    /// `searchBestAllele(sampleIndex, evidenceIndex, canBeReference, priorities)`.
    ///
    /// `priorities` is indexed by allele and is consulted only when the best and second best are
    /// within the informative threshold of each other, which is the tie-breaking pass
    /// `bestAllelesBreakingTies` exists for.
    pub fn search_best_allele(
        &self,
        sample_index: usize,
        evidence_index: usize,
        can_be_reference: bool,
        priorities: Option<&[f64]>,
    ) -> BestAllele {
        let allele_count = self.alleles.number_of_alleles();
        let reference_is_first = self.reference_allele_index == Some(0);

        if allele_count == 0 || (allele_count == 1 && reference_is_first && !can_be_reference) {
            return self.best_allele(
                sample_index,
                evidence_index,
                None,
                f64::NEG_INFINITY,
                None,
                f64::NEG_INFINITY,
            );
        }

        let sample_values = &self.values[sample_index];
        let mut best_allele_index = if can_be_reference || !reference_is_first {
            0
        } else {
            1
        };
        // Both indices start at zero. Nothing here means "no second best".
        let mut second_best_index = 0usize;
        let mut best_likelihood = sample_values[best_allele_index][evidence_index];
        let mut second_best_likelihood = f64::NEG_INFINITY;

        // Indexed rather than iterated: the loop compares against `reference_allele_index` and
        // writes the index it is on into two variables, so the index is the subject and not an
        // artefact of the traversal.
        #[allow(clippy::needless_range_loop)]
        for a in (best_allele_index + 1)..allele_count {
            if !can_be_reference && self.reference_allele_index == Some(a) {
                continue;
            }
            let candidate = sample_values[a][evidence_index];
            // Strictly greater: an equal likelihood never displaces the earlier allele.
            if candidate > best_likelihood {
                second_best_index = best_allele_index;
                best_allele_index = a;
                second_best_likelihood = best_likelihood;
                best_likelihood = candidate;
            } else if candidate > second_best_likelihood {
                second_best_index = a;
                second_best_likelihood = candidate;
            }
        }

        if let Some(priorities) = priorities {
            if best_likelihood - second_best_likelihood < self.informative_threshold() {
                let mut best_priority = priorities[best_allele_index];
                let mut second_best_priority = priorities[second_best_index];
                for a in 0..allele_count {
                    let candidate = sample_values[a][evidence_index];
                    if a == best_allele_index
                        || (!can_be_reference && self.reference_allele_index == Some(a))
                        || best_likelihood - candidate > self.informative_threshold()
                    {
                        continue;
                    }
                    let candidate_priority = priorities[a];
                    if candidate_priority > best_priority {
                        second_best_index = best_allele_index;
                        best_allele_index = a;
                        second_best_priority = best_priority;
                        best_priority = candidate_priority;
                    } else if candidate_priority > second_best_priority {
                        second_best_index = a;
                        second_best_priority = candidate_priority;
                    }
                }
                // `bestPriority` and `secondBestPriority` are dead after the loop in the reference
                // too: the likelihoods are re-read below from the indices.
                let _ = (best_priority, second_best_priority);
            }
        }

        let best = sample_values[best_allele_index][evidence_index];
        let second_best = if second_best_index != best_allele_index {
            sample_values[second_best_index][evidence_index]
        } else {
            f64::NEG_INFINITY
        };

        self.best_allele(
            sample_index,
            evidence_index,
            Some(best_allele_index),
            best,
            Some(second_best_index),
            second_best,
        )
    }

    fn best_allele(
        &self,
        sample_index: usize,
        evidence_index: usize,
        best_allele_index: Option<usize>,
        likelihood: f64,
        second_best_allele_index: Option<usize>,
        second_best_likelihood: f64,
    ) -> BestAllele {
        BestAllele {
            sample: self
                .samples
                .get_sample(sample_index)
                .cloned()
                .unwrap_or_default(),
            evidence_index,
            allele: best_allele_index.and_then(|i| self.alleles.get_allele(i).cloned()),
            second_best_allele: second_best_allele_index
                .and_then(|i| self.alleles.get_allele(i).cloned()),
            likelihood,
            second_best_likelihood,
            // The guard's real job is the pair of negative infinities, whose difference is NaN.
            confidence: if likelihood == second_best_likelihood {
                0.0
            } else {
                likelihood - second_best_likelihood
            },
        }
    }

    /// `bestAllelesBreakingTies(tieBreakingPriority)` over every sample, in sample then evidence
    /// order.
    pub fn best_alleles_breaking_ties(&self, priorities: Option<&[f64]>) -> Vec<BestAllele> {
        (0..self.number_of_samples())
            .flat_map(|sample| self.best_alleles_breaking_ties_for_sample(sample, priorities))
            .collect()
    }

    /// `marginalize(newToOldAlleleMap)`: fewer alleles, each taking the **maximum** likelihood of
    /// the old ones it stands for.
    ///
    /// ```java
    /// newSampleValues[newAllele][r] = oldAlleleSet.stream()
    ///         .mapToDouble(oldA -> oldSampleValues[oldA][r]).max().orElse(NEGATIVE_INFINITY);
    /// ```
    ///
    /// A maximum, not a sum: marginalising two alleles a read supports equally leaves that
    /// likelihood unchanged rather than doubling it.
    ///
    /// The **order** of the new alleles is the caller's, and in the reference the caller is a
    /// `Collectors.toMap` whose key set iterates in `HashMap` order. That order is observable,
    /// because `searchBestAllele` breaks a tie by taking the first index, so which allele a tied
    /// read is attributed to depends on `Allele.hashCode`. [`crate::java_hash::hash_map_order`]
    /// reproduces it, and the caller passes the result here.
    ///
    /// An old allele the map does not mention contributes to nothing, which the reference supports
    /// and calls "typically not the case".
    pub fn marginalize(
        &self,
        new_to_old: &[(Allele, Vec<Allele>)],
    ) -> Result<AlleleLikelihoods<E>, LikelihoodsError> {
        let mut new_values: Vec<Vec<Vec<f64>>> = Vec::with_capacity(self.number_of_samples());
        for sample in 0..self.number_of_samples() {
            let evidence_count = self.sample_evidence_count(sample);
            let mut sample_values: Vec<Vec<f64>> = Vec::with_capacity(new_to_old.len());
            for (_, old_alleles) in new_to_old {
                let mut row = vec![f64::NEG_INFINITY; evidence_count];
                for old in old_alleles {
                    let Some(old_index) = self.index_of_allele(old) else {
                        // `IllegalArgumentException("missing old allele ...")`.
                        return Err(LikelihoodsError::UnknownSample(String::new()));
                    };
                    for (evidence, slot) in row.iter_mut().enumerate() {
                        let value = self.value(sample, old_index, evidence);
                        if value > *slot {
                            *slot = value;
                        }
                    }
                }
                sample_values.push(row);
            }
            new_values.push(sample_values);
        }

        let new_alleles: Vec<Allele> = new_to_old
            .iter()
            .map(|(allele, _)| allele.clone())
            .collect();
        let mut result = AlleleLikelihoods::new(
            self.samples.clone(),
            AlleleList::new(&new_alleles),
            self.evidence_by_sample.clone(),
            new_values,
        )?;
        // `result.isNaturalLog = isNaturalLog`, which the informative threshold depends on.
        result.is_natural_log = self.is_natural_log;
        Ok(result)
    }

    /// `bestAllelesBreakingTies(sample, tieBreakingPriority)`.
    pub fn best_alleles_breaking_ties_for_sample(
        &self,
        sample_index: usize,
        priorities: Option<&[f64]>,
    ) -> Vec<BestAllele> {
        (0..self.sample_evidence_count(sample_index))
            .map(|evidence| self.search_best_allele(sample_index, evidence, true, priorities))
            .collect()
    }
}

impl AlleleLikelihoods<htsjdk_bam::record::BamRecord> {
    /// `groupEvidence(GATKRead::getName, Fragment::createAndAvoidFailure)`, the only instantiation
    /// the annotations use.
    ///
    /// The log likelihoods of a group are **summed**, not averaged, and the resulting evidence order
    /// is a `HashMap`'s over the read names. Both are recorded in [`crate::fragment`].
    pub fn group_by_fragment(
        &self,
    ) -> Result<AlleleLikelihoods<crate::fragment::Fragment>, LikelihoodsError> {
        let allele_count = self.number_of_alleles();
        let mut evidence_by_sample: Vec<Vec<crate::fragment::Fragment>> = Vec::new();
        let mut values: Vec<Vec<Vec<f64>>> = Vec::new();

        for sample in 0..self.number_of_samples() {
            let reads = self.sample_evidence(sample).unwrap_or(&[]);
            let groups = crate::fragment::group_by_read_name(reads);
            let mut fragments = Vec::with_capacity(groups.len());
            let mut sample_values: Vec<Vec<f64>> = vec![vec![0.0; groups.len()]; allele_count];
            for (new_index, group) in groups.iter().enumerate() {
                for (allele, row) in sample_values.iter_mut().enumerate() {
                    for old_index in group {
                        row[new_index] += self.value(sample, allele, *old_index);
                    }
                }
                let group_reads: Vec<htsjdk_bam::record::BamRecord> =
                    group.iter().map(|index| reads[*index].clone()).collect();
                // `createAndAvoidFailure` never fails on a non-empty group, and a group produced by
                // the grouping is never empty.
                match crate::fragment::Fragment::create_and_avoid_failure(&group_reads) {
                    Ok(fragment) => fragments.push(fragment),
                    Err(_) => return Err(LikelihoodsError::UnknownSample(String::new())),
                }
            }
            evidence_by_sample.push(fragments);
            values.push(sample_values);
        }

        let mut grouped = AlleleLikelihoods::new(
            self.samples.clone(),
            self.alleles.clone(),
            evidence_by_sample,
            values,
        )?;
        grouped.is_natural_log = self.is_natural_log;
        Ok(grouped)
    }
}
