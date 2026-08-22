//! `MergeMutect2CallsWithMC3`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.validation.MergeMutect2CallsWithMC3`
//! (GATK 4.6.2.0).
//!
//! Mutect2 calls merged into the MC3 pan-cancer call set. The walk is
//! [`gatk_engine::concordance_walker`]; what is here is what each of the five states writes.
//!
//! # The two sides are not treated alike
//!
//! A true positive and a filtered false negative keep the MC3 record whole and only add to it, so
//! MC3's annotations are authoritative. A false positive is rebuilt from scratch:
//!
//! ```java
//! new VariantContextBuilder(m2.getSource(), m2.getContig(), m2.getStart(), m2.getEnd(), m2.getAlleles())
//!     .attribute(CENTERS_KEY, M2_CENTER_NAME).genotypes(genotypes).make()
//! ```
//!
//! so the M2 record's ID, QUAL, FILTER column and every INFO field are dropped, and only
//! `CENTERS=M2` survives. And a filtered true negative is skipped outright, so an eval-only
//! filtered call disappears without a trace.
//!
//! A false negative never learns that M2 looked at it: emitted unchanged except for its genotype,
//! with no `CENTERS` added. Three of the five states add the field; two do not.
//!
//! # The genotype's ploidy is the number of alleles at the site
//!
//! ```java
//! new GenotypeBuilder(tumorSample, truthVersusEval.getTruthIfPresentElseEval().getAlleles())
//! ```
//!
//! Every allele rather than a called pair, so a multiallelic false positive comes out `0/1/2`. And
//! since the depths are taken straight from the M2 genotype, a biallelic record can end up
//! carrying an AD of three entries.
//!
//! An M2 genotype without `AD` leaves the output genotype without one: `getAD()` answers null and
//! `GenotypeBuilder.AD(null)` sets nothing, rather than throwing or writing zeroes. The `NREF` and
//! `NALT` fallback does write zeroes, and only a truth-only record reaches it.

use gatk_engine::concordance_walker::{ConcordanceRecord, ConcordanceState, TruthVersusEval};

/// `CENTERS_KEY`.
pub const CENTERS_KEY: &str = "CENTERS";
/// `M2_CENTER_NAME`.
pub const M2_CENTER_NAME: &str = "M2";
/// `M2_FILTERS_KEY`.
pub const M2_FILTERS_KEY: &str = "M2_FILTERS";
/// `MC3_REF_COUNT_KEY`.
pub const MC3_REF_COUNT_KEY: &str = "NREF";
/// `MC3_ALT_COUNT_KEY`.
pub const MC3_ALT_COUNT_KEY: &str = "NALT";

/// One input record, reduced to what this tool reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub contig: String,
    pub start: i32,
    pub id: String,
    pub quality: String,
    pub filters: Vec<String>,
    /// Reference first, then the alternates, as the VCF declares them.
    pub alleles: Vec<String>,
    /// The INFO fields, in the order the file gives them.
    pub info: Vec<(String, String)>,
    /// The first genotype's `AD`, absent when the genotype has none.
    pub allele_depths: Option<Vec<i32>>,
}

/// The walk orders by contig index then start, and asks each record whether it is filtered.
///
/// `PASS` is not a filter: htsjdk's `isFiltered` asks whether the filter SET is non-empty, and
/// `PASS` leaves it empty.
impl ConcordanceRecord for Variant {
    fn contig(&self) -> &str {
        &self.contig
    }

    fn start(&self) -> i32 {
        self.start
    }

    fn is_filtered(&self) -> bool {
        !self.filters.is_empty() && self.filters != ["PASS"]
    }
}

impl Variant {
    fn attribute(&self, key: &str) -> Option<&str> {
        self.info
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }

    /// `getAttributeAsInt(key, 0)`.
    fn attribute_as_int(&self, key: &str) -> i32 {
        self.attribute(key)
            .and_then(|value| value.parse().ok())
            .unwrap_or(0)
    }
}

/// What the tool writes for one step, or nothing at all.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Merged {
    pub contig: String,
    pub start: i32,
    pub id: String,
    pub quality: String,
    pub filters: Vec<String>,
    pub alleles: Vec<String>,
    /// The INFO fields as the builder leaves them, before the writer sorts them by key.
    pub info: Vec<(String, String)>,
    /// The genotype's called alleles, one index per allele at the site.
    pub genotype: Vec<usize>,
    pub allele_depths: Option<Vec<i32>>,
}

/// `areVariantsAtSameLocusConcordant`: the same reference allele, and the eval containing the
/// truth's FIRST alternate.
pub fn concordant(truth: &Variant, eval: &Variant) -> bool {
    let same_reference = truth.alleles.first() == eval.alleles.first();
    let truth_alternate = match truth.alleles.get(1) {
        Some(allele) => allele,
        None => return false,
    };
    same_reference && eval.alleles[1..].contains(truth_alternate)
}

/// `makeVariantContextBuilderWithM2Center`: appended to, not replaced.
fn with_m2_center(mc3: &Variant) -> Vec<(String, String)> {
    let mut info = mc3.info.clone();
    match info.iter_mut().find(|(name, _)| name == CENTERS_KEY) {
        Some((_, value)) => {
            value.push(',');
            value.push_str(M2_CENTER_NAME);
        }
        None => info.push((CENTERS_KEY.to_string(), M2_CENTER_NAME.to_string())),
    }
    info
}

/// The genotype the tool builds: every allele at the site, and the depths from whichever side has
/// them.
fn genotype_for(truth: Option<&Variant>, eval: Option<&Variant>) -> (Vec<usize>, Option<Vec<i32>>) {
    // `getAD()` when M2 is there, the MC3 counts when it is not.
    let depths = match eval {
        Some(eval) => eval.allele_depths.clone(),
        None => truth.map(|truth| {
            vec![
                truth.attribute_as_int(MC3_REF_COUNT_KEY),
                truth.attribute_as_int(MC3_ALT_COUNT_KEY),
            ]
        }),
    };
    // `getTruthIfPresentElseEval().getAlleles()`, which is every allele at the site.
    let site = truth.or(eval).expect("a step with at least one record");
    ((0..site.alleles.len()).collect(), depths)
}

/// `apply` for one step. `None` is a filtered true negative, which writes nothing.
pub fn apply(step: &TruthVersusEval, truth: &[Variant], eval: &[Variant]) -> Option<Merged> {
    let truth_record = step.truth.map(|index| &truth[index]);
    let eval_record = step.eval.map(|index| &eval[index]);
    let (genotype, allele_depths) = genotype_for(truth_record, eval_record);

    match step.state {
        ConcordanceState::TruePositive => {
            let mc3 = truth_record.expect("a true positive has truth");
            Some(Merged {
                contig: mc3.contig.clone(),
                start: mc3.start,
                id: mc3.id.clone(),
                quality: mc3.quality.clone(),
                filters: mc3.filters.clone(),
                alleles: mc3.alleles.clone(),
                info: with_m2_center(mc3),
                genotype,
                allele_depths,
            })
        }
        ConcordanceState::FalsePositive => {
            // Rebuilt from the site and the alleles alone: no ID, no QUAL, no filters, no INFO.
            let m2 = eval_record.expect("a false positive has eval");
            Some(Merged {
                contig: m2.contig.clone(),
                start: m2.start,
                id: ".".to_string(),
                quality: ".".to_string(),
                filters: Vec::new(),
                alleles: m2.alleles.clone(),
                info: vec![(CENTERS_KEY.to_string(), M2_CENTER_NAME.to_string())],
                genotype,
                allele_depths,
            })
        }
        ConcordanceState::FalseNegative => {
            // Unchanged except for the genotype: no CENTERS, nothing else added.
            let mc3 = truth_record.expect("a false negative has truth");
            Some(Merged {
                contig: mc3.contig.clone(),
                start: mc3.start,
                id: mc3.id.clone(),
                quality: mc3.quality.clone(),
                filters: mc3.filters.clone(),
                alleles: mc3.alleles.clone(),
                info: mc3.info.clone(),
                genotype,
                allele_depths,
            })
        }
        ConcordanceState::FilteredTrueNegative => None,
        ConcordanceState::FilteredFalseNegative => {
            let mc3 = truth_record.expect("a filtered false negative has truth");
            let m2 = eval_record.expect("a filtered false negative has eval");
            let mut info = with_m2_center(mc3);
            info.push((M2_FILTERS_KEY.to_string(), m2.filters.join(",")));
            Some(Merged {
                contig: mc3.contig.clone(),
                start: mc3.start,
                id: mc3.id.clone(),
                quality: mc3.quality.clone(),
                filters: mc3.filters.clone(),
                alleles: mc3.alleles.clone(),
                info,
                genotype,
                allele_depths,
            })
        }
    }
}

/// Every step of a whole run, in the order the writer receives them.
pub fn merge(steps: &[TruthVersusEval], truth: &[Variant], eval: &[Variant]) -> Vec<Merged> {
    steps
        .iter()
        .filter_map(|step| apply(step, truth, eval))
        .collect()
}
