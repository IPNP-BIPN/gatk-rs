//! `AbstractConcordanceWalker.ConcordanceIterator`, ported from
//! `org.broadinstitute.hellbender.engine.AbstractConcordanceWalker` (GATK 4.7.0.0;
//! `AbstractConcordanceWalker.java` and `ConcordanceState.java` are unchanged since 4.6.2.0), with
//! `ConcordanceState` from `org.broadinstitute.hellbender.tools.walkers.validation`.
//!
//! The golden is `tests/data/concordance_walker.txt.gz`, replayed by `tests/concordance_walker.rs`;
//! the step numbers below are its `state` rows.
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
//! golden (`default-filters`, steps 5 to 7) that is three steps from one disagreement: truth
//! `chr1:200 AT/A` and eval `chr1:200 A/C` disagree, so truth 200 is a false negative; eval 200
//! then sits before truth 210 and comes out a false positive; and truth 210 is a false negative of
//! its own. A port pairing locally would emit one step here and disagree with the reference on
//! everything downstream.
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
//! through both filter sets and the second one (`filtered-dropped-both-sides`) has no filtered
//! state in it at all. The base class's two filters are:
//!
//! ```java
//! protected Predicate<VariantContext> makeTruthVariantFilter() {
//!     return vc -> !vc.isFiltered();
//! }
//! ```
//!
//! and `protected Predicate<VariantContext> makeEvalVariantFilter() { return vc -> true; }`.
//!
//! # The order is the dictionary's
//!
//! `VariantContextComparator` compares the contig's index in the sequence dictionary and then the
//! start. Nothing looks at the end of a record, so a spanning deletion sorts by where it begins.
//! The golden's `across-contigs` case is the dictionary order at work: truth `chr1:500` comes out
//! before the `chr2:100` pair although 500 is the larger start.
//!
//! # Reading this file without reading Rust
//!
//! - `trait ConcordanceRecord` is a Java interface; `<T: ConcordanceRecord>` makes a function
//!   generic over any record type that implements it, as `<T extends ConcordanceRecord>` would.
//! - `impl Fn(&T, &T) -> bool` is a function value, Java's `BiPredicate<T, T>`; here it stands in
//!   for the abstract `areVariantsAtSameLocusConcordant` a subclass overrides.
//! - `Option<usize>` is "an index or nothing", the reference's `Optional<VariantContext>`.
//! - `std::cmp::Ordering` is the sign of a `compare`: `Less`, `Equal` or `Greater`, and
//!   `a.cmp(&b)` is `Integer.compare(a, b)`.

/// `ConcordanceState`, with the abbreviations the summary tables print.
///
/// ```java
/// TRUE_POSITIVE("TP"), FALSE_POSITIVE("FP"), FALSE_NEGATIVE("FN"),
///         FILTERED_TRUE_NEGATIVE("FTN"), FILTERED_FALSE_NEGATIVE("FFN");
/// ```
///
/// (One line in the reference, wrapped here.)
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

    /// The enum constant's own name, which is what an unexpected-state message quotes: Java's
    /// `toString()` of an enum, as in `TruthVersusEval.getTruth`'s
    /// `"This is a " + concordanceState.toString() + " and has no truth VariantContext."`.
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
    /// `truth`: the index into the filtered truth records, `None` for an eval-only step.
    pub truth: Option<usize>,
    /// `eval`: the index into the filtered eval records, `None` for a truth-only step.
    pub eval: Option<usize>,
    /// `concordanceState`.
    pub state: ConcordanceState,
}

/// As much of a record as the iterator itself looks at.
pub trait ConcordanceRecord {
    /// `VariantContext.getContig()`.
    fn contig(&self) -> &str;
    /// `VariantContext.getStart()`.
    fn start(&self) -> i32;
    /// `VariantContext.isFiltered()`, which is htsjdk's
    /// `filters != null && !filters.isEmpty()`.
    fn is_filtered(&self) -> bool;
}

/// `VariantContextComparator`: the contig's index in the dictionary, then the start.
///
/// The reference subtracts (htsjdk 5.0.0, wrapped here):
///
/// ```java
/// final int contigCompare = this.contigIndexLookup.get(firstVariantContext.getContig())
///         - this.contigIndexLookup.get(secondVariantContext.getContig());
/// return contigCompare != 0 ? contigCompare
///         : firstVariantContext.getStart() - secondVariantContext.getStart();
/// ```
///
/// Only the sign is used, and for positions the subtraction cannot overflow, so `cmp` gives the
/// same answer.
fn compare<T: ConcordanceRecord>(left: &T, right: &T, dictionary: &[String]) -> std::cmp::Ordering {
    let index = |record: &T| {
        dictionary
            .iter()
            .position(|contig| contig == record.contig())
            // **Divergence.** A contig the dictionary does not name is not refused when the
            // reference's comparator is built (that checks only for an empty list and duplicates):
            // its `contigIndexLookup.get(..)` answers `null` and unboxing it throws
            // `NullPointerException` inside `compare`, as htsjdk's own comment there says. Here it
            // is ordered last instead, which keeps this function total and does not throw.
            .unwrap_or(usize::MAX)
    };
    index(left)
        .cmp(&index(right))
        .then(left.start().cmp(&right.start()))
}

/// The whole traversal: every step, in order. `ConcordanceIterator` drained.
///
/// `truth` and `eval` are the records that survived the walker's own filters, in file order.
/// `concordant` is `areVariantsAtSameLocusConcordant`, which each tool defines for itself.
///
/// `t` and `e` are the two iterators' positions; `truth[t]` is `truthIterator.peek()`, and
/// `t += 1` is the `next()` that consumes it. The branches are `next()`'s, in its order:
///
/// ```java
/// } else if(evalIterator.peek().isFiltered()) {
///     return TruthVersusEval.filteredFalseNegative(truthIterator.next(), evalIterator.next());
/// } else if (areVariantsAtSameLocusConcordant(truthIterator.peek(), evalIterator.peek())) {
///     return TruthVersusEval.truePositive(truthIterator.next(), evalIterator.next());
/// ```
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

        // `positionCompare > 0`: eval is behind, so it goes alone.
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

/// `nextEvalOnlyVariant`: "it is either a false positive or a filtered true negative".
///
/// ```java
/// return evalVariant.isFiltered() ? TruthVersusEval.filteredTrueNegative(evalVariant)
///         : TruthVersusEval.falsePositive(evalVariant);
/// ```
///
/// (One line in the reference, wrapped here.)
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

/// `nextTruthOnlyVariant`: always a false negative, filtered or not. Under the base class's truth
/// filter no filtered truth record reaches the iterator at all.
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
