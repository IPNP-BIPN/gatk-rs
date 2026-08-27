//! `LearnReadOrientationModel`: which of twelve states a site is in.
//!
//! Eight of the states are orientation artefacts, one per alternate base per read orientation, and
//! four are real: hom ref, germline het, somatic het and hom var. The tool fits the prior over
//! those states by EM over a counts file.
//!
//! The EM loop and the counts file are not ported here; the counts file is
//! [`crate::collect_f1r2_counts`]'s. The two functions the fit is built from are: the flat prior it
//! starts from, and the responsibilities one site takes given a prior.

use gatk_engine::beta_binomial::BetaBinomialDistribution;
use gatk_engine::natural_log_utils::normalize_from_log_to_linear_space;

/// `ArtifactState`, in the order every prior array is indexed in: the eight artefacts first, the
/// four real states after.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State {
    F1R2A,
    F1R2C,
    F1R2G,
    F1R2T,
    F2R1A,
    F2R1C,
    F2R1G,
    F2R1T,
    HomRef,
    GermlineHet,
    SomaticHet,
    HomVar,
}

/// `F1R2FilterConstants.NUM_STATES`.
pub const NUM_STATES: usize = 12;

/// The four bases, as the states name them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Base {
    A,
    C,
    G,
    T,
}

impl Base {
    pub fn name(self) -> &'static str {
        match self {
            Base::A => "A",
            Base::C => "C",
            Base::G => "G",
            Base::T => "T",
        }
    }

    pub fn parse(text: &str) -> Option<Base> {
        match text {
            "A" => Some(Base::A),
            "C" => Some(Base::C),
            "G" => Some(Base::G),
            "T" => Some(Base::T),
            _ => None,
        }
    }
}

impl State {
    pub fn all() -> &'static [State; NUM_STATES] {
        &[
            State::F1R2A,
            State::F1R2C,
            State::F1R2G,
            State::F1R2T,
            State::F2R1A,
            State::F2R1C,
            State::F2R1G,
            State::F2R1T,
            State::HomRef,
            State::GermlineHet,
            State::SomaticHet,
            State::HomVar,
        ]
    }

    pub fn index(self) -> usize {
        State::all()
            .iter()
            .position(|state| *state == self)
            .expect("a known state")
    }

    pub fn name(self) -> &'static str {
        match self {
            State::F1R2A => "F1R2_A",
            State::F1R2C => "F1R2_C",
            State::F1R2G => "F1R2_G",
            State::F1R2T => "F1R2_T",
            State::F2R1A => "F2R1_A",
            State::F2R1C => "F2R1_C",
            State::F2R1G => "F2R1_G",
            State::F2R1T => "F2R1_T",
            State::HomRef => "HOM_REF",
            State::GermlineHet => "GERMLINE_HET",
            State::SomaticHet => "SOMATIC_HET",
            State::HomVar => "HOM_VAR",
        }
    }

    /// `getAltAlleleOfArtifact`, which is absent for the four real states.
    pub fn artifact_base(self) -> Option<Base> {
        match self {
            State::F1R2A | State::F2R1A => Some(Base::A),
            State::F1R2C | State::F2R1C => Some(Base::C),
            State::F1R2G | State::F2R1G => Some(Base::G),
            State::F1R2T | State::F2R1T => Some(Base::T),
            _ => None,
        }
    }

    pub fn is_artifact(self) -> bool {
        self.artifact_base().is_some()
    }
}

/// `getRefToRefArtifacts`: the two artefact states whose alternate base IS the reference base, and
/// which are therefore just hom ref wearing another name.
pub fn ref_to_ref_artifacts(reference: Base) -> [State; 2] {
    match reference {
        Base::A => [State::F1R2A, State::F2R1A],
        Base::C => [State::F1R2C, State::F2R1C],
        Base::G => [State::F1R2G, State::F2R1G],
        Base::T => [State::F1R2T, State::F2R1T],
    }
}

/// `getFlatPrior`, which is not flat over twelve states: the two ref-to-ref artefacts take zero and
/// the remaining TEN share the mass, so no state ever gets a twelfth and the prior depends on the
/// reference base.
pub fn flat_prior(reference: Base) -> [f64; NUM_STATES] {
    let skipped = ref_to_ref_artifacts(reference);
    let mut prior = [1.0 / (NUM_STATES - skipped.len()) as f64; NUM_STATES];
    for state in skipped {
        prior[state.index()] = 0.0;
    }
    prior
}

/// The two beta shapes each state carries: one over the alternate fraction, one over the alternate
/// F1R2 fraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BetaShape {
    pub alpha: f64,
    pub beta: f64,
}

// The hyperparameters, named as the reference names them.
const ALT_PSEUDOCOUNT: f64 = 1.0;
const REF_PSEUDOCOUNT: f64 = 9.0;
const PSEUDOCOUNT_OF_HOM_LIKELY: f64 = 10000.0;
const PSEUDOCOUNT_OF_HOM_UNLIKELY: f64 = 3.0;
const BALANCED_HET_PSEUDOCOUNT: f64 = 5.0;
const BALANCED_F1R2_PRIOR: f64 = 10.0;
const PSEUDOCOUNT_OF_SOMATIC_ALT: f64 = 2.0;
const PSEUDOCOUNT_OF_SOMATIC_REF: f64 = 5.0;
const PSEUDOCOUNT_OF_LIKELY_OUTCOME: f64 = 100.0;
const PSEUDOCOUNT_OF_RARE_OUTCOME: f64 = 1.0;

/// `getPseudoCountsForAlleleFraction`. Both orientations of one base share a shape, because the
/// alternate fraction knows nothing about read orientation.
pub fn allele_fraction_shape(state: State) -> BetaShape {
    match state {
        State::HomRef => BetaShape {
            alpha: PSEUDOCOUNT_OF_HOM_UNLIKELY,
            beta: PSEUDOCOUNT_OF_HOM_LIKELY,
        },
        State::GermlineHet => BetaShape {
            alpha: BALANCED_HET_PSEUDOCOUNT,
            beta: BALANCED_HET_PSEUDOCOUNT,
        },
        State::SomaticHet => BetaShape {
            alpha: PSEUDOCOUNT_OF_SOMATIC_ALT,
            beta: PSEUDOCOUNT_OF_SOMATIC_REF,
        },
        State::HomVar => BetaShape {
            alpha: PSEUDOCOUNT_OF_HOM_LIKELY,
            beta: PSEUDOCOUNT_OF_HOM_UNLIKELY,
        },
        _ => BetaShape {
            alpha: ALT_PSEUDOCOUNT,
            beta: REF_PSEUDOCOUNT,
        },
    }
}

/// `getPseudoCountsForAltF1R2Fraction`. This is where the two orientations differ: an F1R2 artefact
/// expects its alternate reads to be F1R2 and an F2R1 one expects the opposite, in the same
/// hundred-to-one proportion.
pub fn alt_f1r2_fraction_shape(state: State) -> BetaShape {
    match state {
        State::F1R2A | State::F1R2C | State::F1R2G | State::F1R2T => BetaShape {
            alpha: PSEUDOCOUNT_OF_LIKELY_OUTCOME,
            beta: PSEUDOCOUNT_OF_RARE_OUTCOME,
        },
        State::F2R1A | State::F2R1C | State::F2R1G | State::F2R1T => BetaShape {
            alpha: PSEUDOCOUNT_OF_RARE_OUTCOME,
            beta: PSEUDOCOUNT_OF_LIKELY_OUTCOME,
        },
        _ => BetaShape {
            alpha: BALANCED_F1R2_PRIOR,
            beta: BALANCED_F1R2_PRIOR,
        },
    }
}

/// `computeLogPosterior`: the state's own prior, then the alternate depth given the depth, then the
/// alternate F1R2 count given the alternate depth.
pub fn log_posterior(
    alt_depth: i32,
    f1r2_alt_count: i32,
    depth: i32,
    state_prior: f64,
    allele_fraction: BetaShape,
    alt_f1r2_fraction: BetaShape,
) -> f64 {
    let over_depth =
        BetaBinomialDistribution::new(allele_fraction.alpha, allele_fraction.beta, depth)
            .expect("a valid shape")
            .log_probability(alt_depth)
            .expect("a valid count");
    let over_alt =
        BetaBinomialDistribution::new(alt_f1r2_fraction.alpha, alt_f1r2_fraction.beta, alt_depth)
            .expect("a valid shape")
            .log_probability(f1r2_alt_count)
            .expect("a valid count");
    state_prior.ln() + over_depth + over_alt
}

/// `computeResponsibilities`.
///
/// Two whole classes of state are ruled out before any arithmetic. A ref-to-ref artefact is just
/// hom ref, and an artefact state whose base is not the OBSERVED alternate has an indicator of
/// zero, so at most TWO of the eight artefact states can ever be non-zero for one site: F1R2 and
/// F2R1 of the observed base. Which of the two takes the mass is decided by the F1R2 count alone.
///
/// `given_not_hom_ref` zeroes hom ref AFTER the posteriors are computed, so what is left is
/// renormalised rather than rescaled.
pub fn compute_responsibilities(
    reference: Base,
    alternate: Base,
    alt_depth: i32,
    f1r2_alt_count: i32,
    depth: i32,
    artifact_prior: &[f64; NUM_STATES],
    given_not_hom_ref: bool,
) -> [f64; NUM_STATES] {
    let ref_to_ref = ref_to_ref_artifacts(reference);
    let mut log_unnormalized = [f64::NEG_INFINITY; NUM_STATES];
    for state in State::all() {
        let index = state.index();
        if ref_to_ref.contains(state) {
            continue;
        }
        if state.is_artifact() && state.artifact_base() != Some(alternate) {
            continue;
        }
        log_unnormalized[index] = log_posterior(
            alt_depth,
            f1r2_alt_count,
            depth,
            artifact_prior[index],
            allele_fraction_shape(*state),
            alt_f1r2_fraction_shape(*state),
        );
    }
    if given_not_hom_ref {
        log_unnormalized[State::HomRef.index()] = f64::NEG_INFINITY;
    }
    let normalized = normalize_from_log_to_linear_space(&log_unnormalized)
        .expect("a finite sum over the twelve states");
    let mut out = [0.0; NUM_STATES];
    out.copy_from_slice(&normalized);
    out
}

/// What the engine's constructor refuses, and what it does not refuse and then crashes on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelError {
    ContextLength {
        context: String,
    },
    NonCanonicalKmer {
        context: String,
    },
    /// Not the constructor's own: an empty alt design matrix passes both validations and then asks
    /// for a matrix with zero rows, which Apache Commons refuses.
    EmptyDesignMatrix,
}

impl ModelError {
    pub fn message(&self) -> String {
        match self {
            ModelError::ContextLength { context } => format!(
                "reference context must have length {} but got {context}",
                REFERENCE_CONTEXT_SIZE
            ),
            ModelError::NonCanonicalKmer { context } => {
                format!("{context} is not in the set of canonical kmers")
            }
            ModelError::EmptyDesignMatrix => {
                "0 is smaller than, or equal to, the minimum (0)".to_string()
            }
        }
    }
}

/// `F1R2FilterConstants.REFERENCE_CONTEXT_SIZE`.
pub const REFERENCE_CONTEXT_SIZE: usize = 3;

/// A kmer is canonical when its MIDDLE base is A or C: the other half of the alphabet is reached
/// by reverse complement, so only one of each pair is kept.
pub fn is_canonical_kmer(context: &str) -> bool {
    let bases: Vec<char> = context.chars().collect();
    bases.len() == REFERENCE_CONTEXT_SIZE
        && bases.iter().all(|base| "ACGT".contains(*base))
        && matches!(bases[REFERENCE_CONTEXT_SIZE / 2], 'A' | 'C')
}

/// The constructor's validations, in its own order, and the crash that follows them.
pub fn validate_context(context: &str, alt_design_matrix_size: usize) -> Result<(), ModelError> {
    if context.len() != REFERENCE_CONTEXT_SIZE {
        return Err(ModelError::ContextLength {
            context: context.to_string(),
        });
    }
    if !is_canonical_kmer(context) {
        return Err(ModelError::NonCanonicalKmer {
            context: context.to_string(),
        });
    }
    if alt_design_matrix_size == 0 {
        return Err(ModelError::EmptyDesignMatrix);
    }
    Ok(())
}
