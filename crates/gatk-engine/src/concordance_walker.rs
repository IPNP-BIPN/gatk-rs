//! `AbstractConcordanceWalker.ConcordanceIterator`, ported from
//! `org.broadinstitute.hellbender.engine.AbstractConcordanceWalker` (GATK 4.6.2.0), with
//! `ConcordanceState` from `org.broadinstitute.hellbender.tools.walkers.validation`.
//!
//! Two VCFs walked in lockstep, each step labelled with one of five states. The labels are the
//! visible part; which of the two iterators a step consumes is what decides everything after it.
//!
//! # A same-locus disagreement advances truth alone
//!
//! ```java
//! } else {
//!     // advance truth in case of same-locus discordance -- we could equally well advance eval
//!     return TruthVersusEval.falseNegative(truthIterator.next());
//! }
//! ```
//!
//! The eval record stays where it is and is compared against the **next** truth record. In the
//! golden that is three steps from one disagreement: truth `chr1:200 AT/A` and eval `chr1:200 A/C`
//! disagree, so truth 200 is a false negative; eval 200 then sits before truth 210 and comes out a
//! false positive; and truth 210 is a false negative of its own. A port pairing locally would emit
//! one step here and disagree with the reference on everything downstream.
//!
//! # A filtered eval record is labelled by what truth has
//!
//! At a truth locus it is a filtered **false negative** and consumes both iterators without ever
//! asking whether the two agree. Alone it is a filtered **true negative**. The same record, the
//! same filter, two different labels.
//!
//! # Two of the five states depend on the walker's filters
//!
//! The base class drops filtered **truth** records and keeps every eval record, which is what makes
//! the two filtered states reachable. A walker that drops filtered records on both sides, as
//! `EvaluateInfoFieldConcordance` does, can never produce them: the golden runs the same two files
//! through both filter sets and the second one has no filtered state in it at all.
//!
//! # The order is the dictionary's
//!
//! `VariantContextComparator` compares the contig's index in the sequence dictionary and then the
//! start. Nothing looks at the end of a record, so a spanning deletion sorts by where it begins.

/// `ConcordanceState`, with the abbreviations the summary tables print.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConcordanceState {
    TruePositive,
    FalsePositive,
    FalseNegative,
    FilteredTrueNegative,
    FilteredFalseNegative,
}

impl ConcordanceState {
    /// `getAbbreviation()`.
    pub fn abbreviation(&self) -> &'static str {
        match self {
            ConcordanceState::TruePositive => "TP",
            ConcordanceState::FalsePositive => "FP",
            ConcordanceState::FalseNegative => "FN",
            ConcordanceState::FilteredTrueNegative => "FTN",
            ConcordanceState::FilteredFalseNegative => "FFN",
        }
    }

    /// The enum constant's own name, which is what an unexpected-state message quotes.
    pub fn name(&self) -> &'static str {
        match self {
            ConcordanceState::TruePositive => "TRUE_POSITIVE",
            ConcordanceState::FalsePositive => "FALSE_POSITIVE",
            ConcordanceState::FalseNegative => "FALSE_NEGATIVE",
            ConcordanceState::FilteredTrueNegative => "FILTERED_TRUE_NEGATIVE",
            ConcordanceState::FilteredFalseNegative => "FILTERED_FALSE_NEGATIVE",
        }
    }
}

/// One step of the traversal: `TruthVersusEval`, as indices into the two filtered inputs.
///
/// Indices rather than records, so that a caller can keep whatever record type it already has and
/// so that "the same record twice" is visible as the same index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TruthVersusEval {
    pub truth: Option<usize>,
    pub eval: Option<usize>,
    pub state: ConcordanceState,
}

/// As much of a record as the iterator itself looks at.
pub trait ConcordanceRecord {
    fn contig(&self) -> &str;
    fn start(&self) -> i32;
    fn is_filtered(&self) -> bool;
}

/// `VariantContextComparator`: the contig's index in the dictionary, then the start.
fn compare<T: ConcordanceRecord>(left: &T, right: &T, dictionary: &[String]) -> std::cmp::Ordering {
    let index = |record: &T| {
        dictionary
            .iter()
            .position(|contig| contig == record.contig())
            // A contig the dictionary does not name is a refusal in the reference, thrown when the
            // comparator is built rather than here; ordering it last keeps this function total.
            .unwrap_or(usize::MAX)
    };
    index(left)
        .cmp(&index(right))
        .then(left.start().cmp(&right.start()))
}

/// The whole traversal: every step, in order.
///
/// `truth` and `eval` are the records that survived the walker's own filters, in file order.
/// `concordant` is `areVariantsAtSameLocusConcordant`, which each tool defines for itself.
pub fn concordance<T: ConcordanceRecord>(
    truth: &[T],
    eval: &[T],
    dictionary: &[String],
    concordant: impl Fn(&T, &T) -> bool,
) -> Vec<TruthVersusEval> {
    let mut steps = Vec::new();
    let mut t = 0usize;
    let mut e = 0usize;

    while t < truth.len() || e < eval.len() {
        if t >= truth.len() {
            // Eval alone: filtered or not decides the label.
            steps.push(eval_only(e, eval[e].is_filtered()));
            e += 1;
            continue;
        }
        if e >= eval.len() {
            steps.push(truth_only(t));
            t += 1;
            continue;
        }

        match compare(&truth[t], &eval[e], dictionary) {
            std::cmp::Ordering::Greater => {
                steps.push(eval_only(e, eval[e].is_filtered()));
                e += 1;
            }
            std::cmp::Ordering::Less => {
                steps.push(truth_only(t));
                t += 1;
            }
            std::cmp::Ordering::Equal => {
                if eval[e].is_filtered() {
                    // Both consumed, and the agreement is never tested.
                    steps.push(TruthVersusEval {
                        truth: Some(t),
                        eval: Some(e),
                        state: ConcordanceState::FilteredFalseNegative,
                    });
                    t += 1;
                    e += 1;
                } else if concordant(&truth[t], &eval[e]) {
                    steps.push(TruthVersusEval {
                        truth: Some(t),
                        eval: Some(e),
                        state: ConcordanceState::TruePositive,
                    });
                    t += 1;
                    e += 1;
                } else {
                    // The choice the reference's own comment calls arbitrary: truth alone.
                    steps.push(truth_only(t));
                    t += 1;
                }
            }
        }
    }
    steps
}

fn eval_only(index: usize, filtered: bool) -> TruthVersusEval {
    TruthVersusEval {
        truth: None,
        eval: Some(index),
        state: if filtered {
            ConcordanceState::FilteredTrueNegative
        } else {
            ConcordanceState::FalsePositive
        },
    }
}

fn truth_only(index: usize) -> TruthVersusEval {
    TruthVersusEval {
        truth: Some(index),
        eval: None,
        state: ConcordanceState::FalseNegative,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Record {
        contig: &'static str,
        start: i32,
        reference: &'static str,
        alternates: Vec<&'static str>,
        filtered: bool,
    }

    impl ConcordanceRecord for Record {
        fn contig(&self) -> &str {
            self.contig
        }
        fn start(&self) -> i32 {
            self.start
        }
        fn is_filtered(&self) -> bool {
            self.filtered
        }
    }

    fn record(
        contig: &'static str,
        start: i32,
        reference: &'static str,
        alternate: &'static str,
        filtered: bool,
    ) -> Record {
        Record {
            contig,
            start,
            reference,
            alternates: vec![alternate],
            filtered,
        }
    }

    /// The rule the golden's probe uses: the same reference allele, and truth's first alternate
    /// among eval's.
    fn concordant(truth: &Record, eval: &Record) -> bool {
        truth.reference == eval.reference && eval.alternates.contains(&truth.alternates[0])
    }

    fn dictionary() -> Vec<String> {
        vec!["chr1".to_string(), "chr2".to_string()]
    }

    #[test]
    fn one_disagreement_puts_the_two_sides_out_of_step() {
        // The golden's records at 200 and 210.
        let truth = [
            record("chr1", 200, "AT", "A", false),
            record("chr1", 210, "A", "G", false),
        ];
        let eval = [record("chr1", 200, "A", "C", false)];
        let steps = concordance(&truth, &eval, &dictionary(), concordant);
        assert_eq!(
            steps.iter().map(|step| step.state).collect::<Vec<_>>(),
            vec![
                ConcordanceState::FalseNegative,
                ConcordanceState::FalsePositive,
                ConcordanceState::FalseNegative
            ]
        );
        // The eval record survives the first step and is spent on the second.
        assert_eq!(steps[0].eval, None);
        assert_eq!(steps[1].eval, Some(0));
    }

    #[test]
    fn a_filtered_eval_record_is_labelled_by_what_truth_has() {
        let truth = [record("chr1", 140, "A", "C", false)];
        let eval = [record("chr1", 140, "A", "C", true)];
        let both = concordance(&truth, &eval, &dictionary(), concordant);
        assert_eq!(both[0].state, ConcordanceState::FilteredFalseNegative);
        assert_eq!((both[0].truth, both[0].eval), (Some(0), Some(0)));

        let alone = concordance(&[], &eval, &dictionary(), concordant);
        assert_eq!(alone[0].state, ConcordanceState::FilteredTrueNegative);
    }

    #[test]
    fn dropping_filtered_records_on_both_sides_removes_two_states() {
        let truth = [record("chr1", 140, "A", "C", false)];
        let eval = [record("chr1", 140, "A", "C", true)];
        // What EvaluateInfoFieldConcordance's filters leave the iterator.
        let kept: Vec<&Record> = eval.iter().filter(|record| !record.is_filtered()).collect();
        assert!(kept.is_empty());

        let steps = concordance(&truth, &[], &dictionary(), concordant);
        assert_eq!(steps[0].state, ConcordanceState::FalseNegative);
    }

    #[test]
    fn the_order_is_the_dictionarys() {
        let truth = [
            record("chr1", 500, "A", "C", false),
            record("chr2", 100, "A", "C", false),
        ];
        let eval = [record("chr2", 100, "A", "C", false)];
        let steps = concordance(&truth, &eval, &dictionary(), concordant);
        assert_eq!(
            steps.iter().map(|step| step.state).collect::<Vec<_>>(),
            vec![
                ConcordanceState::FalseNegative,
                ConcordanceState::TruePositive
            ]
        );
    }

    #[test]
    fn the_abbreviations_are_what_the_tables_print() {
        assert_eq!(ConcordanceState::TruePositive.abbreviation(), "TP");
        assert_eq!(ConcordanceState::FilteredTrueNegative.abbreviation(), "FTN");
        assert_eq!(
            ConcordanceState::FilteredFalseNegative.name(),
            "FILTERED_FALSE_NEGATIVE"
        );
    }
}
