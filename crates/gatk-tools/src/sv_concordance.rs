//! `SVConcordance`: which truth variant an evaluation variant is judged to be, and what that
//! judgement writes on it.
//!
//! Not the clustering of [`crate::sv_cluster`]. Each eval record takes ONE truth record, the
//! closest by total breakend distance, under a linkage whose single override is asymmetric.
//!
//! Reading and writing the VCFs and the output sort are the runner's. Which truth record each eval
//! record takes, in which group, the order the records finish in, and every value the annotator
//! writes on them, are here.
//!
//! Since 4.7.0.0 the records are stratified: each eval record is matched once per stratum it falls
//! in and once in the default group, each under its group's linkage, and the highest group that
//! found a match is the one written. Which records two groups finish when is what the output's
//! sort by position keeps between two records that start together, so the streaming flush is
//! ported as the reference has it.
//!
//! The per-sample contingency string and the per-record metrics come from two DIFFERENT schemes:
//! `SVConcordanceAnnotator` uses its own cut-down table for the string, and
//! `GenotypeConcordanceSummaryMetrics` builds a fresh `GA4GHSchemeWithMissingAsHomRef` for the
//! numbers. The two agree on every cell this tool can reach, so the split is not observable here;
//! both are ported because only one of them is the one being read at a time.

use crate::sv_cluster::{CallRecord, Linkage, LinkageResult};
use crate::sv_stratify::SvType;

/// `GenotypeConcordanceStates.TruthState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TruthState {
    Missing,
    HomRef,
    HetRefVar1,
    HetVar1Var2,
    HomVar1,
    NoCall,
    LowGq,
    LowDp,
    VcFiltered,
    GtFiltered,
    IsMixed,
}

/// `GenotypeConcordanceStates.CallState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallState {
    Missing,
    HomRef,
    HetRefVar1,
    HetVar1Var2,
    HomVar1,
    HetRefVar2,
    HetRefVar3,
    HetVar1Var3,
    HetVar3Var4,
    HomVar2,
    HomVar3,
    NoCall,
    LowGq,
    LowDp,
    VcFiltered,
    GtFiltered,
    IsMixed,
}

/// The shared code both enums draw from. Two states are "on the diagonal" when their codes are
/// equal, and every call state that has no truth counterpart shares the one INCOMPARABLE code.
const INCOMPARABLE: i32 = -1;

impl TruthState {
    pub fn code(self) -> i32 {
        match self {
            TruthState::Missing => 0,
            TruthState::HomRef => 1,
            TruthState::HetRefVar1 => 2,
            TruthState::HetVar1Var2 => 3,
            TruthState::HomVar1 => 4,
            TruthState::NoCall => 5,
            TruthState::LowGq => 6,
            TruthState::LowDp => 7,
            TruthState::VcFiltered => 8,
            TruthState::GtFiltered => 9,
            TruthState::IsMixed => 10,
        }
    }

    pub fn all() -> &'static [TruthState] {
        &[
            TruthState::Missing,
            TruthState::HomRef,
            TruthState::HetRefVar1,
            TruthState::HetVar1Var2,
            TruthState::HomVar1,
            TruthState::NoCall,
            TruthState::LowGq,
            TruthState::LowDp,
            TruthState::VcFiltered,
            TruthState::GtFiltered,
            TruthState::IsMixed,
        ]
    }

    /// The column this state indexes in a scheme row.
    fn column(self) -> usize {
        match self {
            TruthState::Missing => 0,
            TruthState::HomRef => 1,
            TruthState::HetRefVar1 => 2,
            TruthState::HetVar1Var2 => 3,
            TruthState::HomVar1 => 4,
            TruthState::NoCall => 5,
            TruthState::LowGq => 6,
            TruthState::LowDp => 7,
            TruthState::VcFiltered => 8,
            TruthState::GtFiltered => 9,
            TruthState::IsMixed => 10,
        }
    }
}

impl CallState {
    pub fn code(self) -> i32 {
        match self {
            CallState::Missing => 0,
            CallState::HomRef => 1,
            CallState::HetRefVar1 => 2,
            CallState::HetVar1Var2 => 3,
            CallState::HomVar1 => 4,
            CallState::NoCall => 5,
            CallState::LowGq => 6,
            CallState::LowDp => 7,
            CallState::VcFiltered => 8,
            CallState::GtFiltered => 9,
            CallState::IsMixed => 10,
            _ => INCOMPARABLE,
        }
    }

    pub fn all() -> &'static [CallState] {
        &[
            CallState::Missing,
            CallState::HomRef,
            CallState::HetRefVar1,
            CallState::HetVar1Var2,
            CallState::HomVar1,
            CallState::HetRefVar2,
            CallState::HetRefVar3,
            CallState::HetVar1Var3,
            CallState::HetVar3Var4,
            CallState::HomVar2,
            CallState::HomVar3,
            CallState::NoCall,
            CallState::LowGq,
            CallState::LowDp,
            CallState::VcFiltered,
            CallState::GtFiltered,
            CallState::IsMixed,
        ]
    }
}

/// `GenotypeConcordanceStates.ContingencyState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Contingency {
    Tp,
    Fp,
    Tn,
    Fn,
    Empty,
    /// Not applicable: a cell the scheme declares impossible.
    Na,
}

impl Contingency {
    fn name(self) -> &'static str {
        match self {
            Contingency::Tp => "TP",
            Contingency::Fp => "FP",
            Contingency::Tn => "TN",
            Contingency::Fn => "FN",
            Contingency::Empty => "EMPTY",
            Contingency::Na => "NA",
        }
    }
}

// The cell values the two schemes are written from.
const TP_ONLY: &[Contingency] = &[Contingency::Tp];
const TP_FN: &[Contingency] = &[Contingency::Tp, Contingency::Fn];
const TP_TN: &[Contingency] = &[Contingency::Tp, Contingency::Tn];
const TP_FP: &[Contingency] = &[Contingency::Tp, Contingency::Fp];
const TP_FP_FN: &[Contingency] = &[Contingency::Tp, Contingency::Fp, Contingency::Fn];
const FP_ONLY: &[Contingency] = &[Contingency::Fp];
const FP_TN: &[Contingency] = &[Contingency::Fp, Contingency::Tn];
const FP_FN: &[Contingency] = &[Contingency::Fp, Contingency::Fn];
const FP_TN_FN: &[Contingency] = &[Contingency::Fp, Contingency::Tn, Contingency::Fn];
const TN_ONLY: &[Contingency] = &[Contingency::Tn];
const TN_FN: &[Contingency] = &[Contingency::Tn, Contingency::Fn];
const FN_ONLY: &[Contingency] = &[Contingency::Fn];
const EMPTY: &[Contingency] = &[Contingency::Empty];
const NA: &[Contingency] = &[Contingency::Na];

/// One row of a scheme: a call state and its eleven truth columns.
type Row = (CallState, [&'static [Contingency]; 11]);

/// A `GenotypeConcordanceScheme`: the rows it declares, with every undeclared cell reading `NA`.
#[derive(Debug, Clone, Copy)]
pub struct Scheme {
    rows: &'static [Row],
}

impl Scheme {
    /// `getConcordanceStateArray`.
    pub fn states(&self, truth: TruthState, call: CallState) -> &'static [Contingency] {
        self.rows
            .iter()
            .find(|(state, _)| *state == call)
            .map(|(_, columns)| columns[truth.column()])
            .unwrap_or(NA)
    }

    /// `getContingencyStateString`, which joins on a comma and never sees an empty array.
    pub fn contingency_string(&self, truth: TruthState, call: CallState) -> String {
        self.states(truth, call)
            .iter()
            .map(|state| state.name())
            .collect::<Vec<&str>>()
            .join(",")
    }
}

/// `SVConcordanceAnnotator.SVGenotypeConcordanceScheme`, which is the GA4GH one with the rows it
/// cannot reach removed, and with `NA` where that one has `EMPTY`.
pub const SV_SCHEME: Scheme = Scheme {
    rows: &[
        (
            CallState::Missing,
            [
                TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, NA, NA, NA, NA, NA,
            ],
        ),
        (
            CallState::HomRef,
            [
                TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, NA, NA, NA, NA, NA,
            ],
        ),
        (
            CallState::HetRefVar1,
            [FP_TN, FP_TN, TP_TN, TP_FN, TP_FN, EMPTY, NA, NA, NA, NA, NA],
        ),
        (
            CallState::HetRefVar2,
            [
                FP_TN, FP_TN, FP_TN_FN, TP_FP_FN, FP_FN, EMPTY, NA, NA, NA, NA, NA,
            ],
        ),
        (
            CallState::HetVar1Var2,
            [
                FP_ONLY, FP_ONLY, TP_FP, TP_ONLY, TP_FP_FN, EMPTY, NA, NA, NA, NA, NA,
            ],
        ),
        (
            CallState::HomVar1,
            [
                FP_ONLY, FP_ONLY, TP_FP, TP_FN, TP_ONLY, EMPTY, NA, NA, NA, NA, NA,
            ],
        ),
        (
            CallState::HomVar2,
            [
                FP_ONLY, FP_ONLY, FP_FN, TP_FN, FP_FN, EMPTY, NA, NA, NA, NA, NA,
            ],
        ),
        (
            CallState::NoCall,
            [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, NA, NA, NA, NA, NA],
        ),
        (
            CallState::IsMixed,
            [EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, NA, NA, NA, NA, NA],
        ),
    ],
};

/// `GA4GHSchemeWithMissingAsHomRef`, which is what the summary metrics build for themselves
/// whatever scheme the caller was using for the strings.
pub const GA4GH_SCHEME: Scheme = Scheme {
    rows: &[
        (
            CallState::Missing,
            [
                TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::HomRef,
            [
                TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::HetRefVar1,
            [
                FP_TN, FP_TN, TP_TN, TP_FN, TP_FN, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::HetRefVar2,
            [NA, NA, FP_TN_FN, NA, FP_FN, NA, NA, NA, NA, NA, NA],
        ),
        (
            CallState::HetRefVar3,
            [NA, NA, NA, FP_FN, NA, NA, NA, NA, NA, NA, NA],
        ),
        (
            CallState::HetVar1Var2,
            [
                FP_ONLY, FP_ONLY, TP_FP, TP_ONLY, TP_FP_FN, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
                EMPTY,
            ],
        ),
        (
            CallState::HetVar1Var3,
            [NA, NA, NA, TP_FP_FN, NA, NA, NA, NA, NA, NA, NA],
        ),
        (
            CallState::HetVar3Var4,
            [
                FP_ONLY, FP_ONLY, FP_FN, FP_FN, FP_FN, NA, NA, NA, NA, NA, NA,
            ],
        ),
        (
            CallState::HomVar1,
            [
                FP_ONLY, FP_ONLY, TP_FP, TP_FN, TP_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::HomVar2,
            [NA, NA, FP_FN, TP_FN, FP_FN, NA, NA, NA, NA, NA, NA],
        ),
        (
            CallState::HomVar3,
            [NA, NA, NA, FP_FN, NA, NA, NA, NA, NA, NA, NA],
        ),
        (
            CallState::NoCall,
            [
                EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::VcFiltered,
            [
                TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::GtFiltered,
            [
                TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::LowGq,
            [
                TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::LowDp,
            [
                TN_ONLY, TN_ONLY, TN_FN, FN_ONLY, FN_ONLY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
        (
            CallState::IsMixed,
            [
                EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY,
            ],
        ),
    ],
};

/// The truth states the metrics group over.
pub const HET_TRUTH_STATES: &[TruthState] = &[TruthState::HetRefVar1, TruthState::HetVar1Var2];
pub const HOM_VAR_TRUTH_STATES: &[TruthState] = &[TruthState::HomVar1];
/// Note that this one includes HOM_REF and MISSING, which is not what its name suggests.
pub const VAR_TRUTH_STATES: &[TruthState] = &[
    TruthState::HetRefVar1,
    TruthState::HetVar1Var2,
    TruthState::HomVar1,
    TruthState::HomRef,
    TruthState::Missing,
];
pub const HET_CALL_STATES: &[CallState] = &[
    CallState::HetRefVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var2,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
];
pub const HOM_VAR_CALL_STATES: &[CallState] =
    &[CallState::HomVar1, CallState::HomVar2, CallState::HomVar3];
pub const VAR_CALL_STATES: &[CallState] = &[
    CallState::HetRefVar1,
    CallState::HetRefVar2,
    CallState::HetRefVar3,
    CallState::HetVar1Var2,
    CallState::HetVar1Var3,
    CallState::HetVar3Var4,
    CallState::HomVar1,
    CallState::HomVar2,
    CallState::HomVar3,
];

/// `GenotypeConcordanceCounts`: how many samples fell in each (truth, call) cell.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Counts {
    cells: Vec<((TruthState, CallState), u64)>,
}

impl Counts {
    pub fn increment(&mut self, truth: TruthState, call: CallState) {
        match self
            .cells
            .iter_mut()
            .find(|((t, c), _)| *t == truth && *c == call)
        {
            Some(cell) => cell.1 += 1,
            None => self.cells.push(((truth, call), 1)),
        }
    }

    pub fn get(&self, truth: TruthState, call: CallState) -> u64 {
        self.cells
            .iter()
            .find(|((t, c), _)| *t == truth && *c == call)
            .map(|(_, count)| *count)
            .unwrap_or(0)
    }

    /// `calculateGenotypeConcordanceUtil`. `missing_sites` is always true here, so no cell is
    /// skipped for being missing.
    fn concordance_util(&self, missing_sites: bool, include_hom_ref: bool) -> f64 {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for truth in TruthState::all() {
            for call in CallState::all() {
                if !missing_sites && is_missing(*truth, *call) {
                    continue;
                }
                if !(include_hom_ref || is_var(*truth, *call)) {
                    continue;
                }
                let count = self.get(*truth, *call) as f64;
                if truth.code() == call.code() {
                    numerator += count;
                }
                denominator += count;
            }
        }
        if denominator > 0.0 {
            numerator / denominator
        } else {
            f64::NAN
        }
    }

    pub fn genotype_concordance(&self, missing_sites: bool) -> f64 {
        self.concordance_util(missing_sites, true)
    }

    pub fn non_ref_genotype_concordance(&self, missing_sites: bool) -> f64 {
        self.concordance_util(missing_sites, false)
    }

    /// `getSensitivity`: TP over TP plus FN, with a cell counted once per contingency it carries.
    pub fn sensitivity(&self, scheme: &Scheme, truth_states: &[TruthState]) -> f64 {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for truth in truth_states {
            for call in CallState::all() {
                let count = self.get(*truth, *call) as f64;
                for state in scheme.states(*truth, *call) {
                    match state {
                        Contingency::Tp => {
                            numerator += count;
                            denominator += count;
                        }
                        Contingency::Fn => denominator += count,
                        _ => {}
                    }
                }
            }
        }
        numerator / denominator
    }

    /// `Ppv`: TP over TP plus FP.
    pub fn ppv(&self, scheme: &Scheme, call_states: &[CallState]) -> f64 {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for call in call_states {
            for truth in TruthState::all() {
                let count = self.get(*truth, *call) as f64;
                for state in scheme.states(*truth, *call) {
                    match state {
                        Contingency::Tp => {
                            numerator += count;
                            denominator += count;
                        }
                        Contingency::Fp => denominator += count,
                        _ => {}
                    }
                }
            }
        }
        numerator / denominator
    }

    /// `getSpecificity`: TN over TN plus FP.
    pub fn specificity(&self, scheme: &Scheme, truth_states: &[TruthState]) -> f64 {
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for truth in truth_states {
            for call in CallState::all() {
                let count = self.get(*truth, *call) as f64;
                for state in scheme.states(*truth, *call) {
                    match state {
                        Contingency::Tn => {
                            numerator += count;
                            denominator += count;
                        }
                        Contingency::Fp => denominator += count,
                        _ => {}
                    }
                }
            }
        }
        numerator / denominator
    }
}

fn is_missing(truth: TruthState, call: CallState) -> bool {
    truth == TruthState::Missing || call == CallState::Missing
}

/// `isVar`, which is an OR over two explicit lists rather than a test that neither side is
/// hom-ref: a hom-ref TRUTH against a variant CALL counts, and so does a variant truth against a
/// hom-ref call. Only the pair where both are non-variant is dropped.
fn is_var(truth: TruthState, call: CallState) -> bool {
    matches!(
        truth,
        TruthState::HomVar1 | TruthState::HetRefVar1 | TruthState::HetVar1Var2
    ) || VAR_CALL_STATES.contains(&call)
}

/// `GenotypeConcordanceSummaryMetrics`, which builds its OWN GA4GH scheme rather than using the
/// one the caller annotated the samples with.
#[derive(Debug, Clone, PartialEq)]
pub struct SummaryMetrics {
    pub genotype_concordance: f64,
    pub non_ref_genotype_concordance: f64,
    pub het_sensitivity: f64,
    pub het_ppv: f64,
    pub homvar_sensitivity: f64,
    pub homvar_ppv: f64,
    pub var_sensitivity: f64,
    pub var_ppv: f64,
    pub var_specificity: f64,
}

impl SummaryMetrics {
    pub fn new(counts: &Counts) -> SummaryMetrics {
        let scheme = &GA4GH_SCHEME;
        SummaryMetrics {
            genotype_concordance: counts.genotype_concordance(true),
            non_ref_genotype_concordance: counts.non_ref_genotype_concordance(true),
            het_sensitivity: counts.sensitivity(scheme, HET_TRUTH_STATES),
            het_ppv: counts.ppv(scheme, HET_CALL_STATES),
            homvar_sensitivity: counts.sensitivity(scheme, HOM_VAR_TRUTH_STATES),
            homvar_ppv: counts.ppv(scheme, HOM_VAR_CALL_STATES),
            var_sensitivity: counts.sensitivity(scheme, VAR_TRUTH_STATES),
            var_ppv: counts.ppv(scheme, VAR_CALL_STATES),
            var_specificity: counts.specificity(scheme, VAR_TRUTH_STATES),
        }
    }
}

/// One sample's call. `None` in `alleles` is a no-call.
#[derive(Debug, Clone, PartialEq)]
pub struct Genotype {
    pub sample: String,
    pub alleles: Vec<Option<i32>>,
    pub copy_number: Option<i32>,
}

impl Genotype {
    pub fn is_no_call(&self) -> bool {
        !self.alleles.is_empty() && self.alleles.iter().all(Option::is_none)
    }

    pub fn is_mixed(&self) -> bool {
        self.alleles.iter().any(Option::is_none) && self.alleles.iter().any(Option::is_some)
    }

    pub fn is_hom_ref(&self) -> bool {
        !self.alleles.is_empty() && self.alleles.iter().all(|allele| *allele == Some(0))
    }

    pub fn is_hom_var(&self) -> bool {
        !self.alleles.is_empty()
            && self
                .alleles
                .iter()
                .all(|allele| matches!(allele, Some(index) if *index > 0))
    }

    /// Two called alleles that are not the same.
    pub fn is_het(&self) -> bool {
        self.alleles.iter().all(Option::is_some)
            && self.alleles.windows(2).any(|pair| pair[0] != pair[1])
    }

    /// `Genotype.sameGenotype`, which compares the allele multiset.
    pub fn same_genotype(&self, other: &Genotype) -> bool {
        let mut mine: Vec<Option<i32>> = self.alleles.clone();
        let mut theirs: Vec<Option<i32>> = other.alleles.clone();
        mine.sort();
        theirs.sort();
        mine == theirs
    }
}

/// `getTruthState`, whose order matters: a genotype that is both hom-ref and fully called is
/// answered before the no-call test is ever reached, and an ABSENT genotype is NO_CALL rather
/// than MISSING.
pub fn truth_state(genotype: Option<&Genotype>) -> TruthState {
    let Some(genotype) = genotype else {
        return TruthState::NoCall;
    };
    if genotype.is_hom_ref() {
        TruthState::HomRef
    } else if genotype.is_het() {
        TruthState::HetRefVar1
    } else if genotype.is_hom_var() {
        TruthState::HomVar1
    } else if genotype.is_no_call() || genotype.alleles.is_empty() {
        TruthState::NoCall
    } else {
        TruthState::IsMixed
    }
}

/// `getEvalState`, the same ladder without the absent case: an eval genotype is never null.
pub fn eval_state(genotype: &Genotype) -> CallState {
    if genotype.is_hom_ref() {
        CallState::HomRef
    } else if genotype.is_het() {
        CallState::HetRefVar1
    } else if genotype.is_hom_var() {
        CallState::HomVar1
    } else if genotype.is_no_call() || genotype.alleles.is_empty() {
        CallState::NoCall
    } else {
        CallState::IsMixed
    }
}

/// AC, AF and AN, as `SVAlleleCounter` computes them or as a VCF already carried them.
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleCounts {
    pub count: Vec<i32>,
    pub frequency: Vec<f64>,
    pub number: i32,
    /// The text the input VCF carried, when these were read rather than counted. A copied value is
    /// written back verbatim, which is how a copy is told from a recount: `0.5` against `0.500`.
    pub verbatim: Option<(String, String, String)>,
}

/// `SVAlleleCounter`: AN counts every CALLED allele, reference ones included, and AF is NaN when
/// there are none.
pub fn count_alleles(alternate_count: usize, genotypes: &[Genotype]) -> AlleleCounts {
    let mut number = 0;
    let mut count = vec![0; alternate_count];
    for genotype in genotypes {
        for allele in genotype.alleles.iter().flatten() {
            number += 1;
            // Allele 0 is the reference, and alternate `n` is index `n - 1` of the count array.
            if *allele > 0 && (*allele as usize) <= alternate_count {
                count[*allele as usize - 1] += 1;
            }
        }
    }
    let frequency = if number == 0 {
        vec![f64::NAN; alternate_count]
    } else {
        count
            .iter()
            .map(|value| *value as f64 / number as f64)
            .collect()
    };
    AlleleCounts {
        count,
        frequency,
        number,
        verbatim: None,
    }
}

/// One record: what the linkage reads plus what the annotator reads.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub call: CallRecord,
    pub genotypes: Vec<Genotype>,
    /// The counts the VCF already carried, if it carried all three.
    pub allele_counts: Option<AlleleCounts>,
}

impl Record {
    pub fn genotype(&self, sample: &str) -> Option<&Genotype> {
        self.genotypes
            .iter()
            .find(|genotype| genotype.sample == sample)
    }

    pub fn is_cnv(&self) -> bool {
        self.call.sv_type == SvType::Cnv
    }
}

/// `SVConcordanceLinkage.areClusterable`.
///
/// Its comment says CNV/DEL and CNV/DUP matching is not allowed. The condition it uses,
/// `(aType == CNV || bType != CNV) && aType != bType`, only refuses the pair when the CNV is the
/// EVAL record: with a non-CNV eval and a CNV truth both halves of the disjunction are false, and
/// the pair falls through to the base linkage. The asymmetry is the behaviour, not the comment.
pub fn are_clusterable(linkage: &Linkage, eval: &CallRecord, truth: &CallRecord) -> LinkageResult {
    if (eval.sv_type == SvType::Cnv || truth.sv_type != SvType::Cnv)
        && eval.sv_type != truth.sv_type
    {
        return LinkageResult::default();
    }
    linkage.linkage(eval, truth)
}

/// `totalDistance`: the sum of both ends.
pub fn total_distance(a: &CallRecord, b: &CallRecord) -> i32 {
    (a.position_a - b.position_a).abs() + (a.position_b - b.position_b).abs()
}

/// `minDistance`: the closer of the two ends, which is the first tiebreaker.
pub fn min_distance(a: &CallRecord, b: &CallRecord) -> i32 {
    (a.position_a - b.position_a)
        .abs()
        .min((a.position_b - b.position_b).abs())
}

/// `genotypeDistance`: the NEGATIVE number of matching genotypes, so more matches sorts first.
pub fn genotype_distance(a: &Record, b: &Record) -> i32 {
    let mut matches = 0;
    for genotype in &a.genotypes {
        if let Some(other) = b.genotype(&genotype.sample) {
            if genotype.same_genotype(other) {
                matches += 1;
            }
        }
    }
    -matches
}

/// The comparator chain `getClosestItem` and `ActiveClosestPair.update` share, in order. The last
/// two exist only so that a tie resolves the same way twice: a candidate whose id equals the eval
/// record's wins, and after that the smaller id does.
fn closeness<'a>(eval: &Record, truth: &'a Record) -> (i32, i32, i32, bool, &'a str) {
    (
        total_distance(&eval.call, &truth.call),
        min_distance(&eval.call, &truth.call),
        genotype_distance(eval, truth),
        truth.call.id != eval.call.id,
        &truth.call.id,
    )
}

/// `getClosestItem`: the closest candidate the linkage allows, keyed by whatever the caller keys
/// them by, the first of two equal ones in the order given.
pub fn closest<'a, K>(
    linkage: &Linkage,
    eval: &Record,
    candidates: impl IntoIterator<Item = (K, &'a Record)>,
) -> Option<(K, LinkageResult)> {
    candidates
        .into_iter()
        .map(|(key, truth)| {
            (
                key,
                truth,
                are_clusterable(linkage, &eval.call, &truth.call),
            )
        })
        .filter(|(_, _, result)| result.result)
        .min_by(|left, right| closeness(eval, left.1).cmp(&closeness(eval, right.1)))
        .map(|(key, _, result)| (key, result))
}

/// `ConcordanceState`, as the two abbreviations that reach the output.
pub const TRUE_POSITIVE: &str = "TP";
pub const FALSE_POSITIVE: &str = "FP";

/// What one annotated eval record carries.
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    pub id: String,
    pub status: &'static str,
    pub truth_variant_id: Option<String>,
    /// The per-sample contingency string, `None` for a sample outside the common set.
    pub contingency: Vec<(String, Option<String>)>,
    /// The per-sample copy-state answer, for a multiallelic CNV only.
    pub truth_copy_number_equal: Vec<(String, Option<bool>)>,
    pub copy_number_concordance: Option<f64>,
    pub metrics: Option<SummaryMetrics>,
    pub allele_counts: Option<AlleleCounts>,
    pub truth_allele_counts: Option<AlleleCounts>,
    /// What the linkage measured against the truth record, written as the four `TRUTH_*`
    /// measurements; each is absent when it was not measured.
    pub truth_linkage: Option<LinkageResult>,
}

/// `copyNumbersMatch`: null unless BOTH sides carry a copy number.
pub fn copy_numbers_match(sample: &str, eval: &Record, truth: Option<&Record>) -> Option<bool> {
    let truth = truth?;
    let evaluated = eval.genotype(sample)?;
    let expected = truth.genotype(sample)?;
    Some(evaluated.copy_number? == expected.copy_number?)
}

/// `SVConcordanceAnnotator.annotate`.
///
/// `common` is the intersection of the two VCFs' sample sets: a sample outside it is annotated on
/// neither side, and is not counted either.
pub fn annotate(
    eval: &Record,
    matched: Option<(&Record, LinkageResult)>,
    common: &[String],
) -> Annotation {
    let truth = matched.map(|(record, _)| record);
    let is_cnv = eval.is_cnv();
    let mut counts = Counts::default();
    let mut contingency = Vec::new();
    let mut copy_number_equal = Vec::new();
    let mut cnv_matches = 0;
    let mut cnv_comparisons = 0;

    for genotype in &eval.genotypes {
        let counted = common.contains(&genotype.sample);
        if !counted {
            contingency.push((genotype.sample.clone(), None));
            copy_number_equal.push((genotype.sample.clone(), None));
            continue;
        }
        if is_cnv {
            let result = copy_numbers_match(&genotype.sample, eval, truth);
            if let Some(matched) = result {
                cnv_comparisons += 1;
                if matched {
                    cnv_matches += 1;
                }
            }
            copy_number_equal.push((genotype.sample.clone(), result));
            contingency.push((genotype.sample.clone(), None));
        } else {
            let truth_state =
                truth_state(truth.and_then(|record| record.genotype(&genotype.sample)));
            let call_state = eval_state(genotype);
            counts.increment(truth_state, call_state);
            contingency.push((
                genotype.sample.clone(),
                Some(SV_SCHEME.contingency_string(truth_state, call_state)),
            ));
            copy_number_equal.push((genotype.sample.clone(), None));
        }
    }

    let copy_number_concordance = if is_cnv {
        if cnv_comparisons == 0 {
            None
        } else {
            Some(cnv_matches as f64 / cnv_comparisons as f64)
        }
    } else {
        None
    };
    let metrics = if !is_cnv && truth.is_some() {
        Some(SummaryMetrics::new(&counts))
    } else {
        None
    };

    // A multiallelic CNV gets no allele counts of its own and no truth ones either.
    let (allele_counts, truth_allele_counts) = if is_cnv {
        (None, None)
    } else {
        let alternate_count = 1;
        let own = match &eval.allele_counts {
            Some(existing) => Some(existing.clone()),
            None if eval.genotypes.is_empty() => None,
            None => Some(count_alleles(alternate_count, &eval.genotypes)),
        };
        let theirs = truth.map(|record| match &record.allele_counts {
            // Carried on the truth record: copied across, text and all.
            Some(existing) => existing.clone(),
            // Absent: recounted over the TRUTH genotypes but against the EVAL record's alternate
            // alleles, which is the eval record's arity rather than the truth record's.
            None => count_alleles(alternate_count, &record.genotypes),
        });
        (own, theirs)
    };

    Annotation {
        id: eval.call.id.clone(),
        status: if truth.is_some() {
            TRUE_POSITIVE
        } else {
            FALSE_POSITIVE
        },
        truth_variant_id: truth.map(|record| record.call.id.clone()),
        contingency,
        truth_copy_number_equal: copy_number_equal,
        copy_number_concordance,
        metrics,
        allele_counts,
        truth_allele_counts,
        truth_linkage: matched.map(|(_, linkage)| linkage),
    }
}

/// What the engine refuses, each reaching the user as the reference's own exception, unwrapped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConcordanceError {
    /// `ClosestSVFinder.add`: a record starting before the last one its group was given.
    Unsorted,
    /// `flushEngineToBuffer`: an eval record finished by its group before the engine registered
    /// it, which only one reaching no further than its own start can be.
    Unregistered { id: u64 },
}

impl ConcordanceError {
    pub fn class(&self) -> &'static str {
        match self {
            ConcordanceError::Unsorted => "java.lang.IllegalArgumentException",
            ConcordanceError::Unregistered { .. } => "java.lang.IllegalStateException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            ConcordanceError::Unsorted => {
                "Items must be added in dictionary-sorted order".to_string()
            }
            ConcordanceError::Unregistered { id } => format!("Unregistered variant id: {id}"),
        }
    }
}

/// A `java.util.HashMap<Long, _>`, kept for the order it iterates in: the buckets in turn, each in
/// insertion order. A bucket is `(h ^ h >>> 16) & (n - 1)` for a table of `n`, which starts at
/// sixteen, doubles when the size passes three quarters of it or when a bucket reaches nine
/// entries below sixty-four, and never shrinks, `clear` included.
///
/// shortcut: a bucket of nine or more on a table of sixty-four or more is a tree whose iteration
/// order this does not follow; upgrade if a group ever holds that many colliding records at once.
#[derive(Debug, Clone)]
struct LongMap<T> {
    capacity: usize,
    entries: Vec<(u64, T)>,
}

impl<T> LongMap<T> {
    fn new() -> Self {
        LongMap {
            capacity: 16,
            entries: Vec::new(),
        }
    }

    fn bucket(&self, key: u64) -> usize {
        let hash = (key ^ (key >> 32)) as u32;
        (hash ^ (hash >> 16)) as usize & (self.capacity - 1)
    }

    fn insert(&mut self, key: u64, value: T) {
        self.entries.push((key, value));
        let bucket = self.bucket(key);
        // `treeifyBin` resizes a table under sixty-four instead of building the tree.
        if self.capacity < 64
            && self
                .entries
                .iter()
                .filter(|(other, _)| self.bucket(*other) == bucket)
                .count()
                > 8
        {
            self.capacity *= 2;
        }
        if self.entries.len() > self.capacity * 3 / 4 {
            self.capacity *= 2;
        }
    }

    /// The positions in `entries`, in the map's iteration order.
    fn order(&self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.entries.len()).collect();
        order.sort_by_key(|&index| self.bucket(self.entries[index].0));
        order
    }
}

/// One eval record waiting in a group: how far it reaches, and the closest truth record so far.
#[derive(Debug, Clone, Copy)]
struct Active {
    max_start: i32,
    closest: Option<(u64, LinkageResult)>,
}

/// `ClosestSVFinder`, unsorted, as every group of `StratifiedConcordanceEngine` builds it.
#[derive(Debug, Clone)]
struct Finder {
    linkage: Linkage,
    /// `truthIdToItemMap`, each truth record with how far it reaches.
    truth: LongMap<i32>,
    /// `idToClusterMap`.
    active: LongMap<Active>,
    last_contig: Option<String>,
    last_start: Option<i32>,
}

impl Finder {
    fn new(linkage: Linkage) -> Finder {
        Finder {
            linkage,
            truth: LongMap::new(),
            active: LongMap::new(),
            last_contig: None,
            last_start: None,
        }
    }

    /// `add`: a truth record updates every eval record waiting, an eval record takes the closest
    /// truth record already held.
    fn add(
        &mut self,
        id: u64,
        is_truth: bool,
        items: &[Record],
        sequences: &[(String, i32)],
    ) -> Result<(), ConcordanceError> {
        let item = &items[id as usize];
        if self
            .last_start
            .is_some_and(|start| start > item.call.position_a)
        {
            return Err(ConcordanceError::Unsorted);
        }
        self.last_contig = Some(item.call.contig_a.clone());
        self.last_start = Some(item.call.position_a);
        let length = sequences
            .iter()
            .find(|(name, _)| *name == item.call.contig_a)
            .map_or(0, |(_, length)| *length);
        let reach = self.linkage.max_clusterable_start(&item.call, length);
        if is_truth {
            self.truth.insert(id, reach);
            for (eval_id, active) in &mut self.active.entries {
                let eval = &items[*eval_id as usize];
                let result = are_clusterable(&self.linkage, &eval.call, &item.call);
                if result.result {
                    // `update`: the new one only when it is strictly closer.
                    let replaces = active.closest.is_none_or(|(current, _)| {
                        closeness(eval, item) < closeness(eval, &items[current as usize])
                    });
                    if replaces {
                        active.closest = Some((id, result));
                    }
                }
            }
        } else {
            let order = self.truth.order();
            let candidates = order.iter().map(|&index| {
                let truth_id = self.truth.entries[index].0;
                (truth_id, &items[truth_id as usize])
            });
            let closest = closest(&self.linkage, item, candidates);
            self.active.insert(
                id,
                Active {
                    max_start: reach,
                    closest,
                },
            );
        }
        Ok(())
    }

    /// `flush` through `flushClusters`: every eval record no later record can reach, or all of
    /// them, in the map's order.
    fn flush(&mut self, force: bool) -> Vec<(u64, Option<(u64, LinkageResult)>)> {
        let order = self.active.order();
        let take: Vec<usize> = if force {
            order
        } else {
            let Some(last) = self.last_start else {
                return Vec::new();
            };
            self.truth.entries.retain(|(_, reach)| *reach >= last);
            order
                .into_iter()
                .filter(|&index| self.active.entries[index].1.max_start < last)
                .collect()
        };
        let out = take
            .iter()
            .map(|&index| {
                let (id, active) = self.active.entries[index];
                (id, active.closest)
            })
            .collect();
        if force {
            self.truth.entries.clear();
            self.active.entries.clear();
            self.last_contig = None;
            self.last_start = None;
        } else {
            let mut index = 0;
            self.active.entries.retain(|_| {
                index += 1;
                !take.contains(&(index - 1))
            });
        }
        out
    }
}

/// `StratifiedConcordanceEngine.ItemTracker`: the groups one eval record went to, and the answer
/// kept so far.
#[derive(Debug, Clone)]
struct Tracker {
    /// In priority order: the strata as the configuration lists them, then the default.
    groups: Vec<String>,
    pending: Vec<String>,
    out: Option<(usize, Option<(u64, LinkageResult)>)>,
}

impl Tracker {
    /// `eject`. The first answer is kept whatever it is, and a match replaces it when it comes from
    /// a higher group or the answer kept is no match; so the highest group that matched wins, and
    /// an eval record nothing matched is the same false positive from every group.
    fn eject(&mut self, name: &str, answer: Option<(u64, LinkageResult)>) {
        self.pending.retain(|group| group != name);
        let priority = self
            .groups
            .iter()
            .position(|group| group == name)
            .expect("a registered group");
        let replaces = match &self.out {
            None => true,
            Some((kept, kept_answer)) => {
                answer.is_some()
                    && (priority < *kept || (priority > *kept && kept_answer.is_none()))
            }
        };
        if replaces {
            self.out = Some((priority, answer));
        }
    }
}

/// One eval record finished by every group it went to.
#[derive(Debug, Clone, PartialEq)]
pub struct Finished {
    /// The eval record, by the order it was added among every record.
    pub id: u64,
    /// The truth record the winning group matched, by the same numbering, with its measurements.
    pub truth: Option<(u64, LinkageResult)>,
    /// `STRAT`: every group it went to, sorted, `default` included.
    pub strata: Vec<String>,
}

/// `StratifiedConcordanceEngine`: one finder per stratum and one for the default group, a truth
/// record given to all of them and an eval record to the strata it matched and the default.
///
/// The groups are walked in the clustering configuration's `HashMap` order and the default last,
/// and every add flushes the group it went to. The order records finish in is therefore the
/// groups' and the maps' own, and it is what the output's stable sort by position keeps between
/// two records that start together.
#[derive(Debug, Clone)]
pub struct Engine {
    /// Every record added, numbered by `nextItemId`.
    items: Vec<Record>,
    groups: Vec<(String, Finder)>,
    default: Finder,
    trackers: std::collections::HashMap<u64, Tracker>,
    finished: Vec<Finished>,
    sequences: Vec<(String, i32)>,
}

impl Engine {
    /// `groups` are the clustering configuration's distinct names in the order they were first
    /// read, each with the linkage of its last row, as `HashMap.put` leaves them.
    pub fn new(
        groups: &[(String, Linkage)],
        default: Linkage,
        sequences: &[(String, i32)],
    ) -> Self {
        let names: Vec<String> = groups.iter().map(|(name, _)| name.clone()).collect();
        let groups = crate::local_assembler_engine::java_hash_map_order(&names)
            .into_iter()
            .map(|index| (groups[index].0.clone(), Finder::new(groups[index].1)))
            .collect();
        Engine {
            items: Vec::new(),
            groups,
            default: Finder::new(default),
            trackers: std::collections::HashMap::new(),
            finished: Vec::new(),
            sequences: sequences.to_vec(),
        }
    }

    /// The record numbered `id`.
    pub fn item(&self, id: u64) -> &Record {
        &self.items[id as usize]
    }

    /// `addTruthVariant`: matched in every group.
    pub fn add_truth(&mut self, record: Record) -> Result<(), ConcordanceError> {
        let id = self.items.len() as u64;
        self.items.push(record);
        for index in 0..self.groups.len() {
            self.add_to(Some(index), id, true)?;
        }
        self.add_to(None, id, true)
    }

    /// `addEvalVariant`, with the names of the strata it matched in the stratification's order.
    pub fn add_eval(&mut self, record: Record, strata: &[String]) -> Result<(), ConcordanceError> {
        let id = self.items.len() as u64;
        self.items.push(record);
        for stratum in strata {
            let index = self
                .groups
                .iter()
                .position(|(name, _)| name == stratum)
                .expect("a validated stratum");
            self.add_to(Some(index), id, false)?;
        }
        self.add_to(None, id, false)?;
        let mut groups = strata.to_vec();
        groups.push(crate::sv_stratify::DEFAULT_STRATUM.to_string());
        self.trackers.insert(
            id,
            Tracker {
                pending: groups.clone(),
                groups,
                out: None,
            },
        );
        Ok(())
    }

    /// `flush`: every group in turn, then the default, and what finished since the last call.
    pub fn flush(&mut self, force: bool) -> Result<Vec<Finished>, ConcordanceError> {
        for index in 0..self.groups.len() {
            self.flush_group(Some(index), force)?;
        }
        self.flush_group(None, force)?;
        Ok(std::mem::take(&mut self.finished))
    }

    /// `addToEngine`: a record on a new contig first flushes the group of everything it holds.
    fn add_to(
        &mut self,
        group: Option<usize>,
        id: u64,
        is_truth: bool,
    ) -> Result<(), ConcordanceError> {
        let contig = &self.items[id as usize].call.contig_a;
        let finder = match group {
            Some(index) => &self.groups[index].1,
            None => &self.default,
        };
        if finder
            .last_contig
            .as_ref()
            .is_some_and(|last| last != contig)
        {
            self.flush_group(group, true)?;
        }
        let finder = match group {
            Some(index) => &mut self.groups[index].1,
            None => &mut self.default,
        };
        finder.add(id, is_truth, &self.items, &self.sequences)?;
        self.flush_group(group, false)
    }

    /// `flushEngineToBuffer`.
    fn flush_group(&mut self, group: Option<usize>, force: bool) -> Result<(), ConcordanceError> {
        let (name, finder) = match group {
            Some(index) => {
                let (name, finder) = &mut self.groups[index];
                (name.clone(), finder)
            }
            None => (
                crate::sv_stratify::DEFAULT_STRATUM.to_string(),
                &mut self.default,
            ),
        };
        for (id, answer) in finder.flush(force) {
            let tracker = self
                .trackers
                .get_mut(&id)
                .ok_or(ConcordanceError::Unregistered { id })?;
            tracker.eject(&name, answer);
            if tracker.pending.is_empty() {
                let tracker = self.trackers.remove(&id).expect("just found");
                let mut strata = tracker.groups;
                strata.sort();
                self.finished.push(Finished {
                    id,
                    truth: tracker.out.and_then(|(_, answer)| answer),
                    strata,
                });
            }
        }
        Ok(())
    }
}
