//! `FamilyLikelihoods`, ported from GATK 4.6.2.0.
//!
//! The family half of `CalculateGenotypePosteriors`: a trio's three genotypes recomputed together,
//! every one of the twenty-seven combinations weighted by how many Mendelian violations it
//! implies. The population half is [`crate::calculate_genotype_posteriors`]; reading the pedigree
//! and the VCF is not ported.
//!
//! # The non-violation coefficient is not one
//!
//! ```java
//! mvCoeff = mvCount>0 ? Math.pow(deNovoPrior,mvCount) : (1.0-10*deNovoPrior-deNovoPrior*deNovoPrior);
//! ```
//!
//! A combination that violates nothing is scaled DOWN by ten times the de novo prior, plus a
//! square term. Nothing in the tool's documentation says so, and at the default prior of 1e-6 it
//! costs every consistent configuration about 4.3e-5 of a log10.
//!
//! # A violation is overturned rather than reported
//!
//! At that default, two hom-ref parents and a hom-var child all come out HET: the de novo weight
//! of 1e-6 is a heavier penalty than the likelihood difference. Raising `--de-novo-prior` to
//! 0.001 leaves the violation standing.
//!
//! # An uncalled parent is a uniform third, and an uncalled child stops the trio
//!
//! ```java
//! if (!hasCalledGT(child.getType()) || (!hasCalledGT(mother.getType()) && !hasCalledGT(father.getType()))) { return; }
//! ```
//!
//! A missing parent's likelihoods become `log10(1/3)` each and the pair is still processed; a
//! missing child returns before the matrix is filled and the whole trio is left alone.
//!
//! # The joint likelihood is taken at the POSTERIOR'S argmax
//!
//! ```java
//! jointTrioLikelihood = motherLikelihoods[maxElementIndex(motherPosteriors)] * ...
//! ```
//!
//! So `JL` reports the likelihood of the configuration the prior chose, not of the one the data
//! preferred. Both joint tags are computed only when all three members are called, and are -1
//! otherwise.
//!
//! # And `log10(1/3)` is a rounded constant
//!
//! ```java
//! private static final double LOG10_OF_ONE_THIRD = -0.4771213;
//! ```
//!
//! Seven digits, not the double nearest to log10(1/3). It is the value an uncalled member's three
//! likelihoods are set to, so it reaches the output.

/// `FamilyLikelihoods.LOG10_OF_ONE_THIRD`, as written: seven digits rather than the exact double.
pub const LOG10_OF_ONE_THIRD: f64 = -0.4771213;
/// `NUM_CALLED_GENOTYPETYPES`.
pub const NUM_CALLED_GENOTYPE_TYPES: usize = 3;
/// `FamilyLikelihoods.NO_JOINT_VALUE`.
pub const NO_JOINT_VALUE: f64 = -1.0;
/// The tool's `--de-novo-prior` default.
pub const DEFAULT_DE_NOVO_PRIOR: f64 = 1e-6;

/// `GenotypeType`, over the three the family engine works with plus the ones it refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenotypeType {
    HomRef,
    Het,
    HomVar,
    /// No call, unavailable or mixed, which `genotypeTypeToValue` all answer -1 for.
    Uncalled,
}

impl GenotypeType {
    /// `genotypeTypeToValue`.
    pub fn value(self) -> Option<usize> {
        match self {
            GenotypeType::HomRef => Some(0),
            GenotypeType::Het => Some(1),
            GenotypeType::HomVar => Some(2),
            GenotypeType::Uncalled => None,
        }
    }

    /// `hasCalledGT`, which excludes MIXED where htsjdk's own `isCalled` would not.
    pub fn is_called(self) -> bool {
        self.value().is_some()
    }

    pub const CALLED: [GenotypeType; 3] = [
        GenotypeType::HomRef,
        GenotypeType::Het,
        GenotypeType::HomVar,
    ];
}

/// One member of a trio, reduced to what the engine reads.
#[derive(Debug, Clone, PartialEq)]
pub struct Member {
    pub sample: String,
    pub genotype_type: GenotypeType,
    /// The `PL`s, absent when the member has none.
    pub likelihoods: Option<Vec<i32>>,
    /// The `PP`s, which are preferred over the `PL`s when present.
    pub posteriors: Option<Vec<i32>>,
}

impl Member {
    fn has_likelihoods(&self) -> bool {
        self.likelihoods.is_some()
    }
}

/// `getCombinationMVCount`: how many Mendelian violations one combination implies.
///
/// A missing parent is simply left out of the count, so a parent/child pair is judged against one
/// parent rather than two.
pub fn combination_mv_count(
    mother: GenotypeType,
    father: GenotypeType,
    child: GenotypeType,
) -> i32 {
    if !child.is_called() {
        return 0;
    }
    let mut parents = Vec::new();
    if mother.is_called() {
        parents.push(mother);
    }
    if father.is_called() {
        parents.push(father);
    }
    if parents.is_empty() {
        return 0;
    }

    let mut reference_alleles = 0;
    let mut alternate_alleles = 0;
    for parent in &parents {
        match parent {
            GenotypeType::HomRef => reference_alleles += 1,
            GenotypeType::Het => {
                reference_alleles += 1;
                alternate_alleles += 1;
            }
            GenotypeType::HomVar => alternate_alleles += 1,
            GenotypeType::Uncalled => {}
        }
    }
    let parent_count = parents.len() as i32;

    if child == GenotypeType::HomRef {
        return if reference_alleles == parent_count {
            0
        } else {
            parent_count - reference_alleles
        };
    }
    if child == GenotypeType::HomVar {
        return if alternate_alleles == parent_count {
            0
        } else {
            parent_count - alternate_alleles
        };
    }
    // A het child needs one of each, unless there is only one parent to ask.
    if child == GenotypeType::Het
        && ((reference_alleles > 0 && alternate_alleles > 0) || parent_count < 2)
    {
        return 0;
    }
    1
}

/// `getLikelihoodMatrixIndex`, over the 3×3×3 matrix laid out flat.
pub fn likelihood_matrix_index(
    mother: GenotypeType,
    father: GenotypeType,
    child: GenotypeType,
) -> Option<usize> {
    let mother = mother.value()?;
    let father = father.value()?;
    let child = child.value()?;
    Some(
        mother * NUM_CALLED_GENOTYPE_TYPES * NUM_CALLED_GENOTYPE_TYPES
            + father * NUM_CALLED_GENOTYPE_TYPES
            + child,
    )
}

/// `GeneralUtils.normalizeFromLog10(values, takeLog10OfOutput = true, keepInLogSpace = true)`.
fn normalize_log10(values: &[f64]) -> Vec<f64> {
    let maximum = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = values.iter().map(|value| 10f64.powf(value - maximum)).sum();
    values
        .iter()
        .map(|value| value - maximum - total.log10())
        .collect()
}

/// `getLikelihoodsAsMapSafeNull`: the `PP` if there is one, then the `PL`, and `log10(1/3)` for an
/// uncalled member.
pub fn likelihoods_as_map(member: Option<&Member>) -> [f64; NUM_CALLED_GENOTYPE_TYPES] {
    let Some(member) = member else {
        return [LOG10_OF_ONE_THIRD; NUM_CALLED_GENOTYPE_TYPES];
    };
    if member.genotype_type.is_called() {
        if let Some(posteriors) = &member.posteriors {
            let vector: Vec<f64> = posteriors
                .iter()
                .map(|value| f64::from(*value) / -10.0)
                .collect();
            let normalized = normalize_log10(&vector);
            return [normalized[0], normalized[1], normalized[2]];
        }
    }
    if !member.genotype_type.is_called() || member.likelihoods.is_none() {
        return [LOG10_OF_ONE_THIRD; NUM_CALLED_GENOTYPE_TYPES];
    }
    let vector: Vec<f64> = member
        .likelihoods
        .as_ref()
        .expect("likelihoods")
        .iter()
        .map(|value| f64::from(*value) / -10.0)
        .collect();
    let normalized = normalize_log10(&vector);
    [normalized[0], normalized[1], normalized[2]]
}

/// Which member of the trio a marginal is being taken for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FamilyMember {
    Mother,
    Father,
    Child,
}

/// `MathUtils.log10sumLog10`.
fn log10_sum_log10(values: &[f64]) -> f64 {
    let maximum = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if maximum.is_infinite() {
        return maximum;
    }
    let total: f64 = values.iter().map(|value| 10f64.powf(value - maximum)).sum();
    maximum + total.log10()
}

/// `MathUtils.scaleLogSpaceArrayForNumericalStability`: the maximum becomes zero.
fn scale_log_space(values: &[f64]) -> Vec<f64> {
    let maximum = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    values.iter().map(|value| value - maximum).collect()
}

/// `getPosteriors`: marginalise the configuration matrix over the other two members.
pub fn posteriors(matrix: &[f64], member: FamilyMember) -> Vec<f64> {
    let mut marginals = [[0.0f64; 9]; NUM_CALLED_GENOTYPE_TYPES];
    let mut counter = 0;
    for outer in GenotypeType::CALLED {
        for inner in GenotypeType::CALLED {
            for (slot, changed) in GenotypeType::CALLED.iter().enumerate() {
                let index = match member {
                    FamilyMember::Mother => likelihood_matrix_index(*changed, outer, inner),
                    FamilyMember::Father => likelihood_matrix_index(outer, *changed, inner),
                    FamilyMember::Child => likelihood_matrix_index(outer, inner, *changed),
                }
                .expect("a called combination");
                marginals[slot][counter] = matrix[index];
            }
            counter += 1;
        }
    }
    let summed: Vec<f64> = marginals.iter().map(|row| log10_sum_log10(row)).collect();
    scale_log_space(&summed)
}

/// The matrix `updateFamilyGenotypes` fills, and whether it filled it at all.
#[derive(Debug, Clone, PartialEq)]
pub struct Configuration {
    pub matrix: Vec<f64>,
}

/// `updateFamilyGenotypes`'s first half: the twenty-seven weighted combinations, or nothing when
/// the child is uncalled or neither parent is.
pub fn configuration_likelihoods(
    mother: Option<&Member>,
    father: Option<&Member>,
    child: Option<&Member>,
    de_novo_prior: f64,
) -> Option<Configuration> {
    let child_type = child.map_or(GenotypeType::Uncalled, |member| member.genotype_type);
    let mother_type = mother.map_or(GenotypeType::Uncalled, |member| member.genotype_type);
    let father_type = father.map_or(GenotypeType::Uncalled, |member| member.genotype_type);
    if !child_type.is_called() || (!mother_type.is_called() && !father_type.is_called()) {
        return None;
    }

    let mother_likelihoods = likelihoods_as_map(mother);
    let father_likelihoods = likelihoods_as_map(father);
    let child_likelihoods = likelihoods_as_map(child);

    let mut matrix = vec![0.0; 27];
    for (child_index, child_genotype) in GenotypeType::CALLED.iter().enumerate() {
        for (mother_index, mother_genotype) in GenotypeType::CALLED.iter().enumerate() {
            for (father_index, father_genotype) in GenotypeType::CALLED.iter().enumerate() {
                let violations =
                    combination_mv_count(*mother_genotype, *father_genotype, *child_genotype);
                let joint = mother_likelihoods[mother_index]
                    + father_likelihoods[father_index]
                    + child_likelihoods[child_index];
                // The consistent case is not 1.0.
                let coefficient = if violations > 0 {
                    de_novo_prior.powi(violations)
                } else {
                    1.0 - 10.0 * de_novo_prior - de_novo_prior * de_novo_prior
                };
                let index =
                    likelihood_matrix_index(*mother_genotype, *father_genotype, *child_genotype)
                        .expect("a called combination");
                matrix[index] = coefficient.log10() + joint;
            }
        }
    }
    Some(Configuration { matrix })
}

/// What one member comes out as.
#[derive(Debug, Clone, PartialEq)]
pub struct Updated {
    pub sample: String,
    /// The `PP`, absent when the member was not updated.
    pub posteriors: Option<Vec<i32>>,
    /// `JL` and `JP`, both -1 when the three were not all called.
    pub joint_likelihood: Option<i32>,
    pub joint_posterior: Option<i32>,
    /// The log10 posteriors, which the caller recalls the genotype from.
    pub log10_posteriors: Option<Vec<f64>>,
}

/// `getUpdatedGenotypes`, for a trio whose matrix has been filled.
pub fn updated_genotypes(
    mother: Option<&Member>,
    father: Option<&Member>,
    child: Option<&Member>,
    configuration: &Configuration,
) -> Vec<Updated> {
    let called = |member: Option<&Member>| {
        member.is_some_and(|member| member.genotype_type.is_called() && member.has_likelihoods())
    };
    let uninformative = [1.0 / 3.0; NUM_CALLED_GENOTYPE_TYPES];
    // Note the asymmetry: these are LINEAR probabilities, where the matrix was built from log10
    // ones, and an uncalled member's third is 1/3 rather than log10(1/3).
    let linear = |member: Option<&Member>| -> [f64; 3] {
        if !called(member) {
            return uninformative;
        }
        let vector: Vec<f64> = member
            .expect("a called member")
            .likelihoods
            .as_ref()
            .expect("likelihoods")
            .iter()
            .map(|value| f64::from(*value) / -10.0)
            .collect();
        let normalized = normalize_log10(&vector);
        [
            10f64.powf(normalized[0]),
            10f64.powf(normalized[1]),
            10f64.powf(normalized[2]),
        ]
    };

    let mother_linear = linear(mother);
    let father_linear = linear(father);
    let child_linear = linear(child);

    let mother_log10 = posteriors(&configuration.matrix, FamilyMember::Mother);
    let father_log10 = posteriors(&configuration.matrix, FamilyMember::Father);
    let child_log10 = posteriors(&configuration.matrix, FamilyMember::Child);

    let linearise = |values: &[f64]| -> Vec<f64> {
        normalize_log10(values)
            .iter()
            .map(|value| 10f64.powf(*value))
            .collect()
    };
    let mother_posteriors = linearise(&mother_log10);
    let father_posteriors = linearise(&father_log10);
    let child_posteriors = linearise(&child_log10);

    let mut joint_likelihood = NO_JOINT_VALUE;
    let mut joint_posterior = NO_JOINT_VALUE;
    if called(child) && called(mother) && called(father) {
        // The likelihood is read at the POSTERIOR'S argmax, not at its own.
        joint_likelihood = mother_linear[argmax(&mother_posteriors)]
            * father_linear[argmax(&father_posteriors)]
            * child_linear[argmax(&child_posteriors)];
        joint_posterior =
            maximum(&mother_posteriors) * maximum(&father_posteriors) * maximum(&child_posteriors);
    }

    vec![
        update(mother, joint_likelihood, joint_posterior, &mother_log10),
        update(father, joint_likelihood, joint_posterior, &father_log10),
        update(child, joint_likelihood, joint_posterior, &child_log10),
    ]
}

fn argmax(values: &[f64]) -> usize {
    let mut best = 0;
    for (index, value) in values.iter().enumerate() {
        if *value > values[best] {
            best = index;
        }
    }
    best
}

fn maximum(values: &[f64]) -> f64 {
    values.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
}

/// `getUpdatedGenotype`, which leaves an uncalled member exactly as it was.
fn update(
    member: Option<&Member>,
    joint_likelihood: f64,
    joint_posterior: f64,
    log10_posteriors: &[f64],
) -> Updated {
    let Some(member) = member else {
        return Updated {
            sample: String::new(),
            posteriors: None,
            joint_likelihood: None,
            joint_posterior: None,
            log10_posteriors: None,
        };
    };
    if !member.genotype_type.is_called() {
        return Updated {
            sample: member.sample.clone(),
            posteriors: None,
            joint_likelihood: None,
            joint_posterior: None,
            log10_posteriors: None,
        };
    }
    Updated {
        sample: member.sample.clone(),
        posteriors: Some(crate::calculate_genotype_posteriors::gls_to_pls(
            log10_posteriors,
        )),
        joint_likelihood: Some(phred_scale_joint(joint_likelihood)),
        joint_posterior: Some(phred_scale_joint(joint_posterior)),
        log10_posteriors: Some(log10_posteriors.to_vec()),
    }
}

/// `QualityUtils.phredScaleLog10ErrorRate(log10(1 - value))`, capped at `Byte.MAX_VALUE` and left
/// at -1 when the joint value was never computed.
fn phred_scale_joint(value: f64) -> i32 {
    if value == NO_JOINT_VALUE {
        return -1;
    }
    let phred = -10.0 * (1.0 - value).log10();
    if phred < f64::from(i8::MAX) {
        phred as i32
    } else {
        i32::from(i8::MAX)
    }
}
