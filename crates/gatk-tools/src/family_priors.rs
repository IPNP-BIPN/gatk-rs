//! `FamilyLikelihoods`, ported from GATK 4.6.2.0.
//!
//! The family half of `CalculateGenotypePosteriors`: a trio's three genotypes recomputed together,
//! every one of the twenty-seven combinations weighted by how many Mendelian violations it
//! implies. The population half is [`crate::calculate_genotype_posteriors`]; the record-level
//! step that runs before it is [`calculate_posterior_gls`].
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

// ---------------------------------------------------------------------------------------------
// The record-level step `CalculateGenotypePosteriors.apply` runs before the population priors
// ---------------------------------------------------------------------------------------------

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotype_type::GenotypeType as HtsjdkType;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

/// A trio `setTrios` keeps: the mother's, the father's and the child's sample names.
#[derive(Debug, Clone, PartialEq)]
pub struct TrioIds {
    pub mother: String,
    pub father: String,
    pub child: String,
}

/// `FamilyLikelihoods.setTrios`: per family of `getFamilies(vcfSamples)`, which already holds only
/// the samples in the VCF, a family of exactly three whose member has both parents in the
/// database and both of them in the family.
///
/// So a family of four (two children) gives no trio at all, and a child whose parents are in the
/// PED but not in the VCF gives none either, where `getTrios()` would count both.
pub fn set_trios(db: &gatk_engine::sample_db::SampleDb, vcf_samples: &[String]) -> Vec<TrioIds> {
    let mut trios = Vec::new();
    for (_, family) in db.families(Some(vcf_samples)) {
        if family.len() != 3 {
            continue;
        }
        for member in &family {
            let parents = db.parents(member);
            if parents.len() == 2 && parents.iter().all(|parent| family.contains(parent)) {
                trios.push(TrioIds {
                    mother: member.maternal_id.clone().expect("a mother"),
                    father: member.paternal_id.clone().expect("a father"),
                    child: member.id.clone(),
                });
            }
        }
    }
    trios
}

/// `hasCalledGT`'s three types, everything else (no call, unavailable, MIXED) uncalled.
fn genotype_type(genotype: &Genotype) -> GenotypeType {
    match genotype.genotype_type() {
        HtsjdkType::HomRef => GenotypeType::HomRef,
        HtsjdkType::Het => GenotypeType::Het,
        HtsjdkType::HomVar => GenotypeType::HomVar,
        _ => GenotypeType::Uncalled,
    }
}

/// The `PP` as integers, which `getLikelihoodsAsMapSafeNull` splits on commas.
fn posteriors_of(genotype: &Genotype) -> Option<Vec<i32>> {
    let text = genotype.get("PP")?.format()?;
    Some(
        text.split(',')
            .filter_map(|piece| piece.trim().parse().ok())
            .collect(),
    )
}

fn member_of(genotype: &Genotype) -> Member {
    Member {
        sample: genotype.sample_name.clone(),
        genotype_type: genotype_type(genotype),
        likelihoods: genotype.pl.clone(),
        posteriors: posteriors_of(genotype),
    }
}

/// `getLikelihoodsAsMapSafeNull`'s refusal of anything but three values, which it reaches for each
/// member before the trio is judged at all, so an uncalled child does not save a haploid parent.
fn ploidy_refusal(genotype: &Genotype) -> Option<String> {
    let called = genotype_type(genotype).is_called();
    let count = match (called, posteriors_of(genotype), &genotype.pl) {
        (true, Some(posteriors), _) => posteriors.len(),
        (true, None, Some(likelihoods)) => likelihoods.len(),
        _ => NUM_CALLED_GENOTYPE_TYPES,
    };
    if count == NUM_CALLED_GENOTYPE_TYPES {
        return None;
    }
    let key = if genotype.get("PP").is_some() {
        "PP"
    } else {
        "PL"
    };
    Some(format!(
        "{} has {count} {key} values, should be {NUM_CALLED_GENOTYPE_TYPES} since only the diploid \
         case is supported when applying family priors.",
        crate::variant_annotator_engine::java_genotype_string(genotype)
    ))
}

/// `calculatePosteriorGLs`: every trio read from the record as it came in and written into the
/// copy, so a sample in two trios (a parent who is also a child) keeps the LAST trio's answer.
/// The error is the `UserException` message a member that is not diploid raises.
pub fn calculate_posterior_gls(
    vc: &VariantContext,
    trios: &[TrioIds],
    de_novo_prior: f64,
) -> Result<VariantContext, String> {
    let mut out = vc.clone();
    let find = |id: &str| vc.genotypes.iter().find(|g| g.sample_name == id);
    let ploidy = vc.genotypes.iter().map(Genotype::ploidy).max().unwrap_or(2);
    for trio in trios {
        let (mother, father, child) = (find(&trio.mother), find(&trio.father), find(&trio.child));
        if mother.is_none() && father.is_none() || child.is_none() {
            continue;
        }
        for genotype in [mother, father, child].into_iter().flatten() {
            if let Some(message) = ploidy_refusal(genotype) {
                return Err(message);
            }
        }
        let members = [
            mother.map(member_of),
            father.map(member_of),
            child.map(member_of),
        ];
        let Some(configuration) = configuration_likelihoods(
            members[0].as_ref(),
            members[1].as_ref(),
            members[2].as_ref(),
            de_novo_prior,
        ) else {
            continue;
        };
        let updated = updated_genotypes(
            members[0].as_ref(),
            members[1].as_ref(),
            members[2].as_ref(),
            &configuration,
        );
        for (original, update) in [mother, father, child].into_iter().zip(updated) {
            let (Some(original), Some(log10)) = (original, update.log10_posteriors.as_ref()) else {
                continue;
            };
            let rebuilt = updated_genotype(vc, original, ploidy, log10, &update);
            if let Some(slot) = out
                .genotypes
                .iter_mut()
                .find(|g| g.sample_name == rebuilt.sample_name)
            {
                *slot = rebuilt;
            }
        }
    }
    Ok(out)
}

/// `put` on the genotype's attribute map: a key already there keeps its place.
fn put(extended: &mut Vec<(String, Value)>, key: &str, value: Value) {
    match extended.iter_mut().find(|(k, _)| k == key) {
        Some(slot) => slot.1 = value,
        None => extended.push((key.to_string(), value)),
    }
}

/// `getUpdatedGenotype` for a called member: `makeGenotypeCall(USE_PLS_TO_ASSIGN)` over the
/// posteriors, then the `PP`, then the original attributes with `JL` and `JP`.
///
/// That last step is `builder.attributes(original + JL + JP)` AFTER `attribute(PP)`, so an input
/// that already carried a `PP` keeps ITS value and the one just computed is lost.
fn updated_genotype(
    vc: &VariantContext,
    original: &Genotype,
    ploidy: usize,
    log10: &[f64],
    update: &Updated,
) -> Genotype {
    let mut genotype = original.clone();
    let no_call = vec![Allele::no_call(); ploidy];
    let no_data = original.is_hom_ref() && original.gq == Some(0) && original.dp == Some(0);
    if no_data {
        // The "no data" case: a plain no-call with everything but the attributes dropped.
        genotype.alleles = no_call;
        genotype.pl = None;
        genotype.dp = None;
        genotype.ad = None;
        genotype.gq = None;
        genotype.extended.clear();
    } else {
        let best = argmax(log10);
        // The biallelic diploid-or-not order: `best` alternate copies.
        let alleles: Vec<Allele> = (0..ploidy)
            .map(|copy| vc.alleles[usize::from(copy >= ploidy - best)].clone())
            .collect();
        let gq = crate::calculate_genotype_posteriors::gq_log10_from_likelihoods(best, log10);
        if alleles
            .iter()
            .any(|allele| allele.display_string() == "<NON_REF>")
        {
            genotype.alleles = vec![vc.alleles[0].clone(); ploidy];
            genotype.pl = Some(vec![0; log10.len()]);
        } else if best == 0 && gq > SUM_GL_THRESH_NOCALL {
            genotype.alleles = no_call;
        } else {
            genotype.alleles = alleles;
        }
        if vc.alleles.len() > 1 {
            // `log10PError`: `Math.round(gq * -10)`, with no cap.
            genotype.gq = Some((gq * -10.0 + 0.5).floor() as i32);
        }
    }
    let pp = Value::List(
        update
            .posteriors
            .as_ref()
            .expect("a called member's posteriors")
            .iter()
            .map(|value| Value::Int(i64::from(*value)))
            .collect(),
    );
    put(&mut genotype.extended, "PP", pp);
    for (key, value) in &original.extended {
        put(&mut genotype.extended, key, value.clone());
    }
    let joint = |value: Option<i32>| Value::Int(i64::from(value.unwrap_or(-1)));
    put(&mut genotype.extended, "JL", joint(update.joint_likelihood));
    put(&mut genotype.extended, "JP", joint(update.joint_posterior));
    genotype
}

/// `GATKVariantContextUtils.SUM_GL_THRESH_NOCALL`.
const SUM_GL_THRESH_NOCALL: f64 = -0.1;
