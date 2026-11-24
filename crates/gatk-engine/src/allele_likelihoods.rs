//! The read side of `org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods`
//! (GATK 4.7.0.0; `AlleleLikelihoods.java` is unchanged since 4.6.2.0): the matrix, the evidence
//! counts, and the best-allele search that every likelihood-reading annotation goes through.
//!
//! The matrix is `[sample][allele][evidence]`, and the two axes are the ones in
//! [`crate::allele_list`]. Besides what a consumer reads, the file holds the mutators the callers
//! reach: `marginalize`, `groupEvidence` (for fragments only), `normalizeLikelihoods`,
//! `filterPoorlyModeledEvidence`, `contaminationDownsampling`, `retainEvidence`, `addEvidence`
//! and `addNonReferenceAllele`. `addMissingAlleles` is not ported in general, only the one call
//! `addNonReferenceAllele` makes of it.
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
//! negative infinity. The second-best **allele** is still read from index 0, so it is the best
//! allele itself, not null: `allele_likelihoods.txt.gz` rows `best one-allele`,
//! `best all-infinite` and `best one-infinite` all report `A*` as both. A port that initialised
//! the index to "missing" would get the likelihood right here and the allele wrong.
//!
//! # Ties go to the lower index, and the comparison is strict
//!
//! `candidateLikelihood > bestLikelihood` never displaces an equal one, so among equal likelihoods
//! the earliest allele wins: row `best tie-between-alts` gives `C` over `G`, both at -1. That is
//! what makes the allele **order** observable in a result that looks order-independent. Two limits
//! on how far it reaches:
//!
//! - with the default priorities of `bestAllelesBreakingTies` a reference less than the
//!   threshold below the best wins wherever it sits, so the index decides only among alleles of
//!   equal priority;
//! - an exact tie has a confidence of 0 and is never `isInformative`, so a count that keeps only
//!   informative evidence skips the tied read rather than crediting it to the first allele. In
//!   `crates/gatk-annotation/tests/data/depth_per_allele.txt.gz` the `triallelic` case has four
//!   reads, the fourth tied between `C` and `G` (the `marginal` rows), and `ad triallelic` is
//!   `1,1,1` with `dp` 3.
//!
//! A NaN is the exception to "ties only": it loses every comparison, so it stays best only at
//! index 0 (row `best with-nan`, a NaN confidence). With the same two values the other way round
//! the finite allele is best with an infinite confidence, which is informative.
//!
//! # `isInformative` and the tie-breaking threshold disagree once the log base changes
//!
//! ```java
//! private double getInformativeThreshold() {
//!     return isNaturalLog ? NATURAL_LOG_INFORMATIVE_THRESHOLD : LOG_10_INFORMATIVE_THRESHOLD;
//! }
//! ```
//!
//! ```java
//! public boolean isInformative() {
//!     return confidence > LOG_10_INFORMATIVE_THRESHOLD;
//! }
//! ```
//!
//! The search converts its threshold to natural log when the matrix has been switched;
//! `BestAllele.isInformative` does not, and keeps comparing against 0.2 whatever the base. After
//! `switch_to_natural_log`, the two questions "was this tie broken by priority" and "is this call
//! informative" are asked at two different thresholds. Both are transcribed as they are. Row
//! `natural switch` shows it: a read 0.1 apart in log10 is not informative, and the same read
//! after the switch (0.2303 apart) is.
//!
//! # `confidence` guards the infinite case rather than the equal case
//!
//! ```java
//! confidence = likelihood == secondBestLikelihood ? 0 : likelihood - secondBestLikelihood;
//! ```
//!
//! The guard reads as "equal likelihoods mean no confidence", and its load-bearing use is a pair
//! of negative infinities, where the subtraction would be `NaN`. That pair comes from a best
//! likelihood that is itself negative infinity (row `best all-infinite`, confidence 0) or from a
//! matrix with no allele to offer (row `best no-alleles`). A one-allele matrix with a finite
//! likelihood does not produce it: row `best one-allele` has a confidence of positive infinity.
//!
//! The confidence is signed and can be negative. The priority pass can make the best an allele
//! with a lower likelihood than the second best (a reference 0.1 below an alternate), and the
//! difference is then below zero.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<usize>` is "an index or nothing"; `None` stands where the reference writes
//!   `MISSING_INDEX` (-1) or returns null, and the compiler makes every reader handle it;
//! - `Result<T, LikelihoodsError>` is "a `T` or one of the refusals below", how this codebase
//!   models a throw. `?` after a call returns its error to the caller at once, like an uncaught
//!   exception at that one call site;
//! - `AlleleLikelihoods<E, A: AlleleType = Allele>` is the generic class: `E` is the evidence type,
//!   `A` the allele type, and `= Allele` is a default used when a caller names only `E`;
//! - an `impl AlleleLikelihoods<BamRecord, A>` block (below) adds methods that exist only for that
//!   instantiation. Java has no counterpart; there `groupEvidence` is generic over a function;
//! - `&self` borrows the matrix read-only and `&mut self` borrows it writably; `&[E]` is a
//!   read-only view of a list that the caller keeps;
//! - `impl Fn(&E) -> bool` is a parameter that accepts any closure of that shape, Java's
//!   `Predicate<EVIDENCE>`; `|x| ...` writes a closure, Java's `x -> ...`;
//! - `let Some(x) = y else { ... }` unwraps an `Option` or takes the exit in the braces;
//! - a body's last expression, without a semicolon, is its return value.

use crate::allele_list::{AlleleList, AlleleType, SampleList};
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
/// matrix has no allele to offer. Like the matrix it comes from, it is generic in the allele type
/// and defaults to `Allele`, which is what every caller here instantiates it at.
#[derive(Debug, Clone, PartialEq)]
pub struct BestAllele<A: AlleleType = Allele> {
    pub sample: String,
    /// Where the reference keeps the evidence object itself (`evidence`), this keeps its index in
    /// the sample's evidence list; the caller reads the object back with `sample_evidence`.
    pub evidence_index: usize,
    pub allele: Option<A>,
    pub second_best_allele: Option<A>,
    pub likelihood: f64,
    pub second_best_likelihood: f64,
    pub confidence: f64,
}

impl<A: AlleleType> BestAllele<A> {
    /// `isInformative()`: the confidence against the **log10** threshold, whatever base the matrix
    /// is in. See the module note. The comparison is on the signed confidence, so a negative one
    /// (a priority win against a better likelihood) is never informative.
    pub fn is_informative(&self) -> bool {
        self.confidence > LOG_10_INFORMATIVE_THRESHOLD
    }
}

/// What the constructor and the mutators refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LikelihoodsError {
    /// The catch-all for a structural refusal that has no variant of its own, raised in this file
    /// only and always with an empty name: per-sample lists whose count is not the sample count
    /// (`new`), an old allele the matrix does not hold (`marginalize`, the reference's
    /// `IllegalArgumentException("missing old allele ...")`), and a fragment that cannot be built
    /// (`group_by_fragment_with`).
    ///
    /// The first has no counterpart in the reference constructor, which looks each listed sample
    /// up in its map and never throws on the result: a key the sample list lacks is never read,
    /// and a listed sample the map lacks gets an empty list, from `setupIndexes`:
    /// `sampleEvidences == null ? new ArrayList<>() : new ArrayList<>(sampleEvidences)`.
    UnknownSample(String),
    /// A sample whose likelihoods do not have the matrix's shape. `allele` is the allele index of a
    /// row whose length is not the sample's evidence count, or, when the sample has the wrong
    /// number of rows, that number of rows.
    WrongRowLength { sample: String, allele: usize },
    /// `switchToNaturalLog` called twice.
    AlreadyNaturalLog,
    /// `normalizeLikelihoods` with a cap that is not negative.
    NonNegativeCap,
    /// `filterPoorlyModeledEvidence` on a matrix with no alleles.
    NoAlleles,
}

/// `AlleleLikelihoods<EVIDENCE, A>`, over an evidence type and an allele type the caller names.
///
/// The reference's evidence is a `GATKRead`; here it is whatever identifies a piece of evidence to
/// the caller, because nothing in this slice reads the evidence itself. The annotations that do
/// (`CountNs` looks at a base, `MappingQualityZero` at a mapping quality) carry their own type.
///
/// The allele axis defaults to `Allele`, which is what the annotations read. `HaplotypeCaller`
/// builds the same class over [`crate::haplotype::Haplotype`], and
/// `HaplotypeFilteringAnnotation` is handed that instantiation alongside the ordinary one; both
/// are this type.
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleLikelihoods<E: Clone + PartialEq, A: AlleleType = Allele> {
    samples: SampleList,
    alleles: AlleleList<A>,
    /// `evidenceBySampleIndex`.
    evidence_by_sample: Vec<Vec<E>>,
    /// `valuesBySampleIndex`, indexed `[sample][allele][evidence]`.
    values: Vec<Vec<Vec<f64>>>,
    /// `referenceAlleleIndex`, `None` where the reference holds `MISSING_INDEX`.
    reference_allele_index: Option<usize>,
    is_natural_log: bool,
    /// `filteredHaplotypeCount`, a plain `int` field of the class whatever it is instantiated over.
    ///
    /// It is written by `AlleleFiltering.subsetHaplotypesByAlleles`, which belongs to the assembler
    /// and therefore to Milestone G3, and read by `HaplotypeFilteringAnnotation`. Between the two
    /// it is carried, not computed, so it is a field here for the same reason it is one there: a
    /// matrix nobody filtered reports zero, and that zero is a real answer rather than a missing
    /// one.
    filtered_haplotype_count: i32,
    /// `filteredEvidenceBySampleIndex`: what `filterPoorlyModeledEvidence` removed, kept for
    /// genotyping.
    filtered_evidence: Vec<Vec<E>>,
}

impl<E: Clone + PartialEq, A: AlleleType> AlleleLikelihoods<E, A> {
    /// The public constructor, with the likelihood values supplied rather than left at zero.
    ///
    /// `evidence_by_sample` and `values` are given per sample **index**, not per name: the
    /// reference's map is immediately turned into index order by `setupIndexes`, which reads only
    /// the names the sample list holds and gives a listed sample missing from the map an empty
    /// list. A caller here passes that empty list itself; a list of the wrong length is refused
    /// (see [`LikelihoodsError::UnknownSample`]).
    ///
    /// The argument order is the public constructor's, samples before alleles; the reference's
    /// private constructor and `createAlleleLikelihoods` take alleles first. `values` is
    /// `[sample][allele][evidence]`, the reference's `valuesBySampleIndex`, while
    /// `evidence_by_sample` has no allele axis. The new matrix is flagged as not natural log, with
    /// no filtered evidence and a filtered haplotype count of zero, as the reference's starts.
    pub fn new(
        samples: SampleList,
        alleles: AlleleList<A>,
        evidence_by_sample: Vec<Vec<E>>,
        values: Vec<Vec<Vec<f64>>>,
    ) -> Result<Self, LikelihoodsError> {
        let sample_count = samples.number_of_samples();
        let allele_count = alleles.number_of_alleles();
        if evidence_by_sample.len() != sample_count || values.len() != sample_count {
            return Err(LikelihoodsError::UnknownSample(String::new()));
        }
        // `.iter().enumerate()` walks a list as `(index, element)` pairs, Java's indexed `for`
        // with the element fetched for you.
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
            filtered_haplotype_count: 0,
            filtered_evidence: vec![Vec::new(); sample_count],
        })
    }

    /// `getFilteredHaplotypeCount()`.
    pub fn filtered_haplotype_count(&self) -> i32 {
        self.filtered_haplotype_count
    }

    /// `setFilteredHaplotypeCount(int)`.
    pub fn set_filtered_haplotype_count(&mut self, count: i32) {
        self.filtered_haplotype_count = count;
    }

    /// `retainEvidence(predicate)`: every sample keeps, in order, the evidence the predicate
    /// accepts, with its likelihoods; the rest is dropped, not filtered.
    ///
    /// Unlike the reference, this does **not** also narrow the filtered evidence:
    ///
    /// ```java
    /// final List<EVIDENCE> sampleFiltered = filteredEvidenceBySampleIndex.get(s).stream()
    ///         .filter(predicate).collect(Collectors.toList());
    /// ```
    ///
    /// In the reference that list is shared by every matrix `marginalize` makes from one source,
    /// so narrowing it for one event narrows it for the next. `hc_genotyping` in `gatk-tools`
    /// keeps its own copy of the list and applies each event's predicate to it, which is where
    /// that sharing is reproduced.
    pub fn retain_evidence(&mut self, keep: impl Fn(&E) -> bool) {
        // `.iter_mut().zip(..)` walks the evidence lists and the likelihood tables side by side,
        // one writable pair per sample.
        for (evidence, sample_values) in self
            .evidence_by_sample
            .iter_mut()
            .zip(self.values.iter_mut())
        {
            // The predicate is asked once per piece of evidence, and the answers are replayed in
            // order to `retain` on the evidence and on every allele row, so all of them drop the
            // same positions. `retain` visits elements in order, which is what makes the replay
            // line up.
            let kept: Vec<bool> = evidence.iter().map(&keep).collect();
            let mut flags = kept.iter();
            evidence.retain(|_| *flags.next().expect("one flag per evidence"));
            for row in sample_values.iter_mut() {
                let mut flags = kept.iter();
                row.retain(|_| *flags.next().expect("one flag per value"));
            }
        }
    }

    /// `changeEvidence`, for one piece of evidence: the likelihoods stay, the evidence is replaced.
    pub fn replace_evidence(&mut self, sample_index: usize, evidence_index: usize, replacement: E) {
        self.evidence_by_sample[sample_index][evidence_index] = replacement;
    }

    /// `addEvidence(evidenceBySample, initialLikelihood)` for one sample: each new piece of
    /// evidence not already present is appended with the initial likelihood for every allele.
    ///
    /// "Already present" includes an earlier copy in the same call, as in the reference's
    /// `appendEvidence`, which reverts its index map when `put` finds a previous entry. The test
    /// here is `==` on the evidence, where the reference's is its hash map's `equals`.
    pub fn add_evidence(&mut self, sample_index: usize, evidence: &[E], initial_likelihood: f64) {
        for item in evidence {
            if self.evidence_by_sample[sample_index].contains(item) {
                continue;
            }
            self.evidence_by_sample[sample_index].push(item.clone());
            for row in self.values[sample_index].iter_mut() {
                row.push(initial_likelihood);
            }
        }
    }

    /// `indexOfReference()`: the first reference allele, `-1` (here `None`) without one.
    ///
    /// The value is the one cached at construction. The reference caches it with
    /// `findReferenceAllele`, whose sequential stream ends in `findAny()`; this takes the first
    /// match, as `AlleleList.indexOfReference` does by its loop.
    pub fn index_of_reference(&self) -> Option<usize> {
        self.reference_allele_index
    }

    /// `numberOfSamples()`.
    pub fn number_of_samples(&self) -> usize {
        self.samples.number_of_samples()
    }

    /// `numberOfAlleles()`.
    pub fn number_of_alleles(&self) -> usize {
        self.alleles.number_of_alleles()
    }

    /// `indexOfSample(sample)`, `None` for the reference's -1.
    pub fn index_of_sample(&self, sample: &str) -> Option<usize> {
        self.samples.index_of_sample(sample)
    }

    /// `indexOfAllele(allele)`, `None` for the reference's -1.
    pub fn index_of_allele(&self, allele: &A) -> Option<usize> {
        self.alleles.index_of_allele(allele)
    }

    /// `getAllele(alleleIndex)`, `None` past the end where the reference throws.
    pub fn get_allele(&self, index: usize) -> Option<&A> {
        self.alleles.get_allele(index)
    }

    /// `getSample(sampleIndex)`, `None` past the end where the reference throws.
    pub fn get_sample(&self, index: usize) -> Option<&String> {
        self.samples.get_sample(index)
    }

    /// `sampleEvidence(sampleIndex)`.
    pub fn sample_evidence(&self, sample_index: usize) -> Option<&[E]> {
        self.evidence_by_sample
            .get(sample_index)
            .map(|list| list.as_slice())
    }

    /// `sampleEvidenceCount(sampleIndex)`. A bad index answers 0 here, where the reference's
    /// `Utils.validIndex` throws.
    pub fn sample_evidence_count(&self, sample_index: usize) -> usize {
        self.evidence_by_sample
            .get(sample_index)
            .map_or(0, |list| list.len())
    }

    /// `evidenceCount()`: the total across samples, which is what `Coverage` writes as `DP`.
    pub fn evidence_count(&self) -> usize {
        self.evidence_by_sample.iter().map(Vec::len).sum()
    }

    /// `isNaturalLog()`.
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

    /// One likelihood, `valuesBySampleIndex[s][a][e]`: the reference reaches it as
    /// `sampleMatrix(s).get(a, e)`, so allele before evidence here too.
    pub fn value(&self, sample_index: usize, allele_index: usize, evidence_index: usize) -> f64 {
        self.values[sample_index][allele_index][evidence_index]
    }

    /// `filteredSampleEvidence(sampleIndex)`: the evidence `filterPoorlyModeledEvidence` removed
    /// from a sample, in removal order. Only the matrix that ran the filter holds it here; see
    /// [`Self::retain_evidence`] and [`Self::marginalize`] for where the reference shares it.
    pub fn filtered_evidence(&self, sample_index: usize) -> &[E] {
        &self.filtered_evidence[sample_index]
    }

    /// `normalizeLikelihoods(maximumLikelihoodDifferenceCap,
    /// symmetricallyNormalizeAllelesToReference)`: no likelihood of a piece of evidence may fall
    /// more than the cap below its best allele's. The best allele may be the reference only when
    /// `symmetric`; negative infinity disables the cap, and a matrix of fewer than two alleles is
    /// left alone.
    ///
    /// The cap is applied to every allele, the reference included, even when the reference was
    /// barred from being the best: `normalizeLikelihoodsPerEvidence` loops over all of them.
    pub fn normalize_likelihoods(
        &mut self,
        cap: f64,
        symmetric: bool,
    ) -> Result<(), LikelihoodsError> {
        if cap >= 0.0 || cap.is_nan() {
            return Err(LikelihoodsError::NonNegativeCap);
        }
        if cap == f64::NEG_INFINITY || self.alleles.number_of_alleles() < 2 {
            return Ok(());
        }
        for sample in 0..self.values.len() {
            for evidence in 0..self.evidence_by_sample[sample].len() {
                let best = self.search_best_allele(sample, evidence, symmetric, None);
                let worst = best.likelihood + cap;
                for row in self.values[sample].iter_mut() {
                    if row[evidence] < worst {
                        row[evidence] = worst;
                    }
                }
            }
        }
        Ok(())
    }

    /// `maximumLikelihoodOverAllAlleles`: a running maximum from negative infinity with a strict
    /// `>`, as in the reference, so a NaN never wins and an all-NaN column answers negative
    /// infinity rather than NaN.
    pub fn maximum_likelihood_over_all_alleles(
        &self,
        sample_index: usize,
        evidence_index: usize,
    ) -> f64 {
        let mut result = f64::NEG_INFINITY;
        for row in &self.values[sample_index] {
            if row[evidence_index] > result {
                result = row[evidence_index];
            }
        }
        result
    }

    /// `filterPoorlyModeledEvidence(log10MinTrueLikelihood)`: remove, sample by sample, the
    /// evidence whose best likelihood is under its threshold, and keep it aside. The threshold is
    /// asked of `(sample, evidence index, evidence)`, since what it reads (the read's HMM
    /// qualities) is a transient attribute in the reference. The reference's is a
    /// `ToDoubleFunction<EVIDENCE>` of the evidence alone; `impl FnMut(..)` accepts a closure that
    /// may update its own state between calls.
    ///
    /// The comparison is a strict `<`, so evidence exactly at its threshold stays.
    pub fn filter_poorly_modeled_evidence(
        &mut self,
        log10_min_true_likelihood: impl FnMut(usize, usize, &E) -> f64,
    ) -> Result<(), LikelihoodsError> {
        self.filter_poorly_modeled_evidence_named(log10_min_true_likelihood, None)
    }

    /// `filterPoorlyModeledEvidence` with `--debug-genotyper-output`, which prints each removed
    /// read by `name` and asks the threshold for it a second time.
    ///
    /// The second asking is the reference's: its debug line calls
    /// `log10MinTrueLikelihood.applyAsDouble(sampleEvidence.get(idx))` again rather than reusing
    /// the first answer, so a closure with state sees two calls per removed read when the debugger
    /// is on. The doubles go through `java_double_to_string` because the reference concatenates
    /// them with `+`. `&dyn Fn(&E) -> String` is a borrowed function object, Java's
    /// `Function<EVIDENCE, String>`.
    pub fn filter_poorly_modeled_evidence_named(
        &mut self,
        mut log10_min_true_likelihood: impl FnMut(usize, usize, &E) -> f64,
        name: Option<&dyn Fn(&E) -> String>,
    ) -> Result<(), LikelihoodsError> {
        if self.alleles.number_of_alleles() == 0 {
            return Err(LikelihoodsError::NoAlleles);
        }
        for sample in 0..self.samples.number_of_samples() {
            let remove: Vec<usize> = (0..self.evidence_by_sample[sample].len())
                .filter(|&i| {
                    self.maximum_likelihood_over_all_alleles(sample, i)
                        < log10_min_true_likelihood(sample, i, &self.evidence_by_sample[sample][i])
                })
                .collect();
            for &i in &remove {
                // `Option::filter` keeps the `Some` only when the test holds, so this branch runs
                // when a name function was given and the debugger is on.
                if let Some(name) = name.filter(|_| crate::genotyper_debug::is_enabled()) {
                    let evidence = &self.evidence_by_sample[sample][i];
                    let maximum = self.maximum_likelihood_over_all_alleles(sample, i);
                    let threshold = log10_min_true_likelihood(sample, i, evidence);
                    crate::genotyper_debug::println(&format!(
                        "disqualified read: {i} {} with max likelihood {} and threshold {}",
                        name(evidence),
                        crate::tsv_table::java_double_to_string(maximum),
                        crate::tsv_table::java_double_to_string(threshold)
                    ));
                }
                let evidence = self.evidence_by_sample[sample][i].clone();
                self.filtered_evidence[sample].push(evidence);
            }
            self.remove_evidence_by_index(sample, &remove);
        }
        Ok(())
    }

    /// `contaminationDownsampling(perSampleDownsamplingFraction)`: for each sample whose fraction
    /// is above zero (and not NaN), all its evidence at one or more, otherwise what
    /// `selectAlleleBiasedEvidence` picks from its informative evidence grouped by best allele
    /// (reference allowed), in allele order.
    ///
    /// `fraction` stands for the map lookup `getOrDefault(samples.getSample(s), 0.0)`, so a caller
    /// answers 0 for a sample it has no fraction for. The grouping searches without priorities,
    /// as `evidenceByBestAlleleMap` does, so the reference allele gets no preference, and keeps
    /// only informative evidence; its groups are in allele order because the reference fills a
    /// `LinkedHashMap` in that order.
    pub fn contamination_downsampling(
        &mut self,
        fraction: impl Fn(&str) -> f64,
        random: &mut crate::well19937c::Well19937c,
    ) -> Result<(), crate::permutation::PermutationError> {
        for sample in 0..self.number_of_samples() {
            let f = fraction(self.get_sample(sample).map_or("", String::as_str));
            if f.is_nan() || f <= 0.0 {
                continue;
            }
            let count = self.evidence_by_sample[sample].len();
            let mut remove: Vec<usize> = if f >= 1.0 {
                (0..count).collect()
            } else {
                // `evidenceByBestAlleleMap(s)`: evidence indices by best allele.
                let mut by_allele: Vec<Vec<usize>> = vec![Vec::new(); self.number_of_alleles()];
                for r in 0..count {
                    let best = self.search_best_allele(sample, r, true, None);
                    if !best.is_informative() {
                        continue;
                    }
                    if let Some(index) = best.allele.as_ref().and_then(|a| self.index_of_allele(a))
                    {
                        by_allele[index].push(r);
                    }
                }
                crate::allele_biased_downsampling::select_allele_biased_evidence(
                    &by_allele, f, random,
                )?
            };
            // `removeEvidence` maps the chosen evidence to indices and takes them
            // `.sorted().distinct()`, which is what `remove_evidence_by_index` assumes.
            remove.sort_unstable();
            remove.dedup();
            self.remove_evidence_by_index(sample, &remove);
        }
        Ok(())
    }

    /// `removeEvidenceByIndex`, for sorted distinct indices. The reference leaves the freed tail of
    /// each likelihood row as NaN; here the rows are shortened, which no reader can tell apart.
    fn remove_evidence_by_index(&mut self, sample_index: usize, remove: &[usize]) {
        if remove.is_empty() {
            return;
        }
        // `binary_search` is only correct on a sorted list, the same assumption the reference
        // states for its `evidencesToRemove`.
        let keep = |i: &usize| remove.binary_search(i).is_err();
        // `std::mem::take` moves the list out and leaves an empty one in its place, so the old
        // list can be consumed while the field is rebuilt from it.
        let evidence = std::mem::take(&mut self.evidence_by_sample[sample_index]);
        self.evidence_by_sample[sample_index] = evidence
            .into_iter()
            .enumerate()
            .filter(|(i, _)| keep(i))
            .map(|(_, e)| e)
            .collect();
        for row in self.values[sample_index].iter_mut() {
            *row = row
                .iter()
                .enumerate()
                .filter(|(i, _)| keep(i))
                .map(|(_, v)| *v)
                .collect();
        }
    }

    /// `switchToNaturalLog()`, which refuses to run twice: each value is multiplied by `LOG_10` in
    /// place, and the flag that `informative_threshold` reads is set.
    pub fn switch_to_natural_log(&mut self) -> Result<(), LikelihoodsError> {
        if self.is_natural_log {
            return Err(LikelihoodsError::AlreadyNaturalLog);
        }
        for sample_values in &mut self.values {
            for row in sample_values.iter_mut() {
                // `iter_mut` hands out writable references; `*value = ...` writes through one.
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
    /// `priorities` is indexed by allele (`None` is the reference's null array) and is consulted
    /// only when the best and second best are within the informative threshold of each other,
    /// which is the tie-breaking pass `bestAllelesBreakingTies` exists for.
    ///
    /// The two threshold tests are not mirror images. The pass starts on
    /// `bestLikelihood - secondBestLikelihood < getInformativeThreshold()` and skips a candidate on
    /// `bestLikelihood - candidateLikelihood > getInformativeThreshold()`. So a pair exactly the
    /// threshold apart does not start the pass, while inside a started pass a candidate exactly
    /// the threshold below the best is still eligible. Both are kept.
    ///
    /// With `can_be_reference` false and the reference at index 0, the first pass starts at index
    /// 1 while the second-best index stays 0, the barred reference. If no later allele displaces
    /// it (two alleles, say), the tail reads the reference's likelihood as the second best. The
    /// only caller passing false is `normalize_likelihoods`, which reads only `likelihood`.
    pub fn search_best_allele(
        &self,
        sample_index: usize,
        evidence_index: usize,
        can_be_reference: bool,
        priorities: Option<&[f64]>,
    ) -> BestAllele<A> {
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
                // The window stays anchored on the first pass's best: `best_likelihood` is not
                // updated when a priority moves `best_allele_index`, in the reference either.
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
                // too: the likelihoods are re-read below from the indices. `let _ = ...` reads and
                // discards them, which keeps the compiler from warning about a dead assignment.
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

    /// The private `BestAllele` constructor: an index of `None` (the reference's `MISSING_INDEX`)
    /// gives a `None` allele, where the reference stores null, and the confidence is computed
    /// here, from the two likelihoods it is handed.
    fn best_allele(
        &self,
        sample_index: usize,
        evidence_index: usize,
        best_allele_index: Option<usize>,
        likelihood: f64,
        second_best_allele_index: Option<usize>,
        second_best_likelihood: f64,
    ) -> BestAllele<A> {
        BestAllele {
            sample: self
                .samples
                .get_sample(sample_index)
                .cloned()
                .unwrap_or_default(),
            evidence_index,
            // `.and_then(..)` maps a `Some(index)` to the allele there and leaves `None` alone.
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
    ///
    /// The reference takes a `Function<A, Double>` and evaluates it over the alleles in order; a
    /// caller here passes those values already evaluated, one per allele, or `None` for the
    /// no-argument form. `flat_map` concatenates the per-sample lists, as the reference's
    /// `flatMap` does.
    pub fn best_alleles_breaking_ties(&self, priorities: Option<&[f64]>) -> Vec<BestAllele<A>> {
        (0..self.number_of_samples())
            .flat_map(|sample| self.best_alleles_breaking_ties_for_sample(sample, priorities))
            .collect()
    }

    /// `marginalize(newToOldAlleleMap)`: fewer alleles, each taking the **maximum** likelihood of
    /// the old ones it stands for.
    ///
    /// ```java
    /// newSampleValues[newAllele][finalR] = oldAlleleSet.stream()
    ///         .mapToDouble(oldA -> oldSampleValues[oldA][finalR])
    ///         .max()
    ///         .orElse(Double.NEGATIVE_INFINITY);
    /// ```
    ///
    /// A maximum, not a sum: marginalising two alleles a read supports equally leaves that
    /// likelihood unchanged rather than doubling it. The running maximum here starts at negative
    /// infinity, which is also the reference's answer for a new allele with no old one.
    ///
    /// The **order** of the new alleles is the caller's. The reference takes a `Map` and uses its
    /// `keySet()` order, so the order is whatever map the caller built: the annotations'
    /// `Collectors.toMap` (`DepthPerAlleleBySample`) iterates in `HashMap` order over
    /// `Allele.hashCode`, which [`crate::java_hash::hash_map_order`] reproduces, while
    /// `createAlleleMapper` builds a `LinkedHashMap` in insertion order. Here `new_to_old` is an
    /// ordered list of `(new allele, old alleles)` pairs standing in for that map, so the caller
    /// states the order rather than inheriting it.
    ///
    /// That order decides which allele a tied read's `BestAllele` names, since `searchBestAllele`
    /// keeps the first of equal likelihoods. For an exact tie it does not move an informative
    /// count, since a tie has a confidence of 0; a NaN likelihood is different (see the module
    /// note on ties).
    ///
    /// An old allele the map does not mention contributes to nothing, which the reference supports
    /// and calls "typically not the case".
    ///
    /// The result carries the same evidence and the same log base as `self` (the reference's
    /// `result.isNaturalLog = isNaturalLog`). It does **not** carry the filtered evidence: the
    /// reference hands the new matrix its own `filteredEvidenceBySampleIndex` list, shared, and
    /// here the new matrix starts with none. See [`Self::retain_evidence`] for the caller that
    /// keeps the shared list itself.
    ///
    /// The new alleles may be of another type than the old ones: `HaplotypeCaller` marginalizes a
    /// matrix over haplotypes into one over the alleles each site genotypes. `<B: AlleleType>` is
    /// that second type parameter, chosen by the caller.
    pub fn marginalize<B: AlleleType>(
        &self,
        new_to_old: &[(B, Vec<A>)],
    ) -> Result<AlleleLikelihoods<E, B>, LikelihoodsError> {
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

        let new_alleles: Vec<B> = new_to_old
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
    ///
    /// `None` is the no-argument form, `a -> a.isReference() ? 1.0 : 0`: a tie goes to the
    /// reference allele wherever it sits in the matrix. The reference never searches with null
    /// priorities here, and taking the first maximum instead agrees only while the reference is
    /// the first allele, which a haplotype matrix or a hash-ordered one need not be.
    pub fn best_alleles_breaking_ties_for_sample(
        &self,
        sample_index: usize,
        priorities: Option<&[f64]>,
    ) -> Vec<BestAllele<A>> {
        // Declared before the `match` and filled inside it, so the list lives to the end of the
        // function and the borrow `&reference_first` taken in the `None` arm stays valid. A list
        // created inside the arm would be freed when the arm ends.
        let reference_first: Vec<f64>;
        let priorities = match priorities {
            Some(given) => given,
            None => {
                reference_first = self
                    .alleles
                    .as_slice()
                    .iter()
                    .map(|a| if a.is_reference() { 1.0 } else { 0.0 })
                    .collect();
                &reference_first
            }
        };
        (0..self.sample_evidence_count(sample_index))
            .map(|evidence| self.search_best_allele(sample_index, evidence, true, Some(priorities)))
            .collect()
    }
}

impl<A: AlleleType> AlleleLikelihoods<htsjdk_bam::record::BamRecord, A> {
    /// `groupEvidence(GATKRead::getName, Fragment::createAndAvoidFailure)`, the only instantiation
    /// the annotations use.
    ///
    /// The log likelihoods of a group are **summed**, not averaged, and the resulting evidence
    /// order is a `HashMap`'s over the read names. Both are recorded in [`crate::fragment`].
    ///
    /// ```java
    /// newLikelihoodValues[s][a][newEvidenceIndex] += oldSampleValues[a][oldEvidenceIndex];
    /// ```
    ///
    /// Two consequences a reader would not guess:
    ///
    /// - the sum covers **every** read of the group, while the fragment built from it holds at
    ///   most two (`createAndAvoidFailure` drops supplementary, secondary and duplicate records,
    ///   then keeps the first two). In `crates/gatk-annotation/tests/data/fragment_counts.txt.gz`
    ///   the `three-reads-one-name` case
    ///   has three reads named `p` (the third supplementary, flag 0x841) at -1 for `A` and -10 for
    ///   `C`: its `frag` row shows fragment `p` holding `p/65,p/129`, and its `grouped` row gives
    ///   that fragment -3.0 and -30.0, the three-read sums;
    /// - the per-read differences add up. With two alleles, mates that each favour one allele by
    ///   0.15 give the fragment 0.3, which clears the 0.2 informative threshold that neither read
    ///   clears alone; mates favouring opposite alleles cancel instead.
    pub fn group_by_fragment(
        &self,
    ) -> Result<AlleleLikelihoods<crate::fragment::Fragment, A>, LikelihoodsError> {
        self.group_by_fragment_with(false)
    }

    /// `groupEvidence(MTAC.independentMates ? read -> read : GATKRead::getName,
    /// Fragment::createAndAvoidFailure)`, as `SomaticGenotypingEngine` calls it. Here, with
    /// independent mates, every read is a fragment of its own, **in evidence order**.
    ///
    /// The reference does not keep that order. Its grouping is
    /// `.collect(Collectors.groupingBy(groupingFunction)).values());`, a `HashMap`, and with
    /// `read -> read` the key is the read itself, whose hash is `Objects.hashCode(samRecord)`
    /// (`SAMRecordToGATKReadAdapter`), a hash over the record's fields. So the reference's
    /// fragments come out in hash order over those records, and two records equal field for
    /// field would share one group. This port does neither.
    ///
    /// Within a group the sum runs over the reads in sample order, from 0.0, as the reference's
    /// `+=` into a fresh `double[]` does, so the floating-point sum is formed in the same order.
    /// The new matrix keeps the log base and starts with no filtered evidence, as the reference's
    /// does (it passes `null` for that list).
    pub fn group_by_fragment_with(
        &self,
        independent_mates: bool,
    ) -> Result<AlleleLikelihoods<crate::fragment::Fragment, A>, LikelihoodsError> {
        let allele_count = self.number_of_alleles();
        let mut evidence_by_sample: Vec<Vec<crate::fragment::Fragment>> = Vec::new();
        let mut values: Vec<Vec<Vec<f64>>> = Vec::new();

        for sample in 0..self.number_of_samples() {
            let reads = self.sample_evidence(sample).unwrap_or(&[]);
            let groups = if independent_mates {
                (0..reads.len()).map(|index| vec![index]).collect()
            } else {
                crate::fragment::group_by_read_name(reads)
            };
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
                // the grouping is never empty. The `Err` arm exists because `match` must cover
                // every case the type allows, not because it is reachable.
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

impl<E: Clone + PartialEq> AlleleLikelihoods<E, Allele> {
    /// `addNonReferenceAllele(Allele.NON_REF_ALLELE)`: `<NON_REF>` appended as the last allele
    /// when it is missing, its row filled by [`Self::update_non_ref_allele_likelihoods`] over
    /// every allele. A matrix that already holds it is left alone.
    ///
    /// The reference gets there through `addMissingAlleles(Collections.singleton(nonRefAllele),
    /// Double.NEGATIVE_INFINITY)` followed by `updateNonRefAlleleLikelihoods()`, which considers
    /// `alleles`, the list that now includes `<NON_REF>`. This is that one call written out: the
    /// new row starts at negative infinity, and the reference index cannot move because the
    /// appended allele is not a reference.
    pub fn add_non_reference_allele(&mut self) {
        // `.expect(..)` unwraps a success or stops the program with the message, an unchecked
        // assertion; parsing this fixed literal cannot fail.
        let non_ref = Allele::from_str("<NON_REF>", false).expect("a symbolic allele");
        if self.alleles.index_of_allele(&non_ref).is_some() {
            return;
        }
        let mut alleles: Vec<Allele> = self.alleles.as_slice().to_vec();
        alleles.push(non_ref);
        self.alleles = AlleleList::new(&alleles);
        for sample in &mut self.values {
            // The evidence count read off the first allele row, or 0 when there is none.
            let evidence = sample.first().map_or(0, Vec::len);
            sample.push(vec![f64::NEG_INFINITY; evidence]);
        }
        let all = alleles.clone();
        self.update_non_ref_allele_likelihoods(&all);
    }

    /// `updateNonRefAlleleLikelihoods(allelesToConsider)`: for each piece of evidence, the
    /// `<NON_REF>` likelihood is the median (Commons Math's, so the legacy percentile) of the
    /// concrete alleles' likelihoods strictly below the best one, among the alleles to consider.
    /// With none below the best, it is the best itself, or NaN when there are fewer than two
    /// concrete alleles.
    ///
    /// ```java
    /// sampleValues[nonRefAlleleIndex][r] = !Double.isNaN(nonRefLikelihood) ? nonRefLikelihood
    ///         : nonSymbolicAlleleCount <= 1 ? Double.NaN : bestAllele.likelihood;
    /// ```
    ///
    /// The best is searched with the `<NON_REF>` row in place, as in the reference, and "concrete"
    /// is every allele but `<NON_REF>` (`alleleCount - 1`), whatever else is symbolic.
    pub fn update_non_ref_allele_likelihoods(&mut self, to_consider: &[Allele]) {
        let non_ref = Allele::from_str("<NON_REF>", false).expect("a symbolic allele");
        let Some(non_ref_index) = self.alleles.index_of_allele(&non_ref) else {
            return;
        };
        let allele_count = self.alleles.number_of_alleles();
        let concrete = allele_count - 1;
        for s in 0..self.values.len() {
            let evidence = self.evidence_by_sample[s].len();
            for r in 0..evidence {
                let best = self.search_best_allele(s, r, true, None).likelihood;
                let mut qualified: Vec<f64> = Vec::with_capacity(concrete);
                for a in 0..allele_count {
                    let likelihood = self.values[s][a][r];
                    let allele = self.alleles.get_allele(a).expect("an allele");
                    if a != non_ref_index
                        && likelihood < best
                        && !likelihood.is_nan()
                        && to_consider.contains(allele)
                    {
                        qualified.push(likelihood);
                    }
                }
                let median = if qualified.is_empty() {
                    f64::NAN
                } else {
                    jmath::percentile::evaluate(
                        &qualified,
                        50.0,
                        jmath::percentile::EstimationType::Legacy,
                    )
                };
                self.values[s][non_ref_index][r] = if !median.is_nan() {
                    median
                } else if concrete <= 1 {
                    f64::NAN
                } else {
                    best
                };
            }
        }
    }
}
