//! Ported from `org.broadinstitute.hellbender.utils.read.CigarBuilder` (GATK 4.7.0.0;
//! `CigarBuilder.java` is unchanged since 4.6.2.0).
//!
//! Not a list of elements. Every cigar GATK produces by clipping is assembled through this
//! builder, and the builder rewrites what it is given:
//!
//!  * consecutive elements with the same operator merge, so `3M` then `4M` is `7M`;
//!  * a zero-length element is dropped;
//!  * a deletion that lands at either end is **removed**, and the bases it covered are counted;
//!  * a deletion arriving after an insertion is **moved before it**, because the order of the two
//!    is arbitrary in SAM and the reference standardises it;
//!  * the section order (hard clip, soft clip, middle, soft clip, hard clip) is validated, and a
//!    cigar that violates it, or that is entirely soft-clipped, is an error rather than a value.
//!
//! A port that concatenated elements would produce a cigar that is a valid string, describes the
//! same alignment, and differs byte for byte from the reference's in the output BAM.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<(), CigarError>` is a `void` method that may throw. `Ok(())` is the normal return,
//!   `Err(..)` the throw, and `?` after a call passes that call's throw on to the caller.
//! - `Option<Op>` is the reference's nullable `lastOperator`: `None` before the first element.
//! - `matches!(x, A | B)` is `x == A || x == B`.
//! - A `match` arm written `Some(last) if condition =>` is taken only when the value is present
//!   **and** the condition holds; otherwise the next arm is tried, top to bottom, which is how the
//!   reference's `if / else if` chain is reproduced here.
//! - `x.is_some_and(|e| ..)` is `x != null && test(x)`.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

/// What the reference throws rather than returning. Every one is `Utils.validate`, so an
/// `IllegalStateException`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CigarError {
    /// `advanceSectionAndValidateCigarOrder`: a soft clip after the right hard clip ("cigar has
    /// already reached its right hard clip", `CigarBuilder.java:155`), or any other non-hard-clip
    /// element after a right clip ("cigar has already reached right clip", line 162).
    SectionOutOfOrder,
    /// `make`: the cigar consists of a single soft clip ("cigar is completely soft-clipped",
    /// `CigarBuilder.java:129`).
    CompletelySoftClipped,
    /// `make(false)`: nothing left after removing flanking deletions ("No cigar elements left
    /// after removing leading and trailing deletions.", `CigarBuilder.java:138`).
    Empty,
}

/// `CigarBuilder.Section`: the order a cigar's parts must appear in. The builder only ever
/// moves forward through it, so a section once left cannot be re-entered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Section {
    LeftHardClip,
    LeftSoftClip,
    Middle,
    RightSoftClip,
    RightHardClip,
}

/// `CigarOperator.isClipping`: `S` or `H`.
fn is_clipping(op: Op) -> bool {
    matches!(op, Op::S | Op::H)
}

/// `CigarBuilder`: elements are given one at a time to [`CigarBuilder::add`], which rewrites as it
/// goes, and [`CigarBuilder::make`] finishes the cigar.
pub struct CigarBuilder {
    /// `cigarElements`. A `Vec` is Java's `ArrayList`.
    elements: Vec<CigarElement>,
    /// `lastOperator`, which is not always the last element's operator. When a right clip replaces
    /// a trailing `D, I` pair the list ends `I, S` and this stays `I`, as the reference's does
    /// (`CigarBuilder.java:89-92` does not assign it). So `5M 3D 2I 10S 5S` builds `5M2I10S5S`: the
    /// second soft clip is appended, not merged.
    last_operator: Option<Op>,
    section: Section,
    /// `true` for the default builder. `false` keeps leading and trailing deletions.
    remove_deletions_at_ends: bool,
    /// Deletion bases dropped from the front by `add`.
    leading_deletion_bases_removed: u32,
    /// Deletion bases dropped from the back by `add`, when a right clip arrived after them.
    trailing_deletion_bases_removed: u32,
    /// Deletion bases dropped from the back by the last `make`, reset by each call.
    trailing_deletion_bases_removed_in_make: u32,
}

/// `new CigarBuilder()`, which is `this(true)`: deletions at the ends are removed.
impl Default for CigarBuilder {
    fn default() -> Self {
        CigarBuilder::new(true)
    }
}

impl CigarBuilder {
    /// `new CigarBuilder(removeDeletionsAtEnds)`.
    pub fn new(remove_deletions_at_ends: bool) -> CigarBuilder {
        CigarBuilder {
            elements: Vec::new(),
            last_operator: None,
            section: Section::LeftHardClip,
            remove_deletions_at_ends,
            leading_deletion_bases_removed: 0,
            trailing_deletion_bases_removed: 0,
            trailing_deletion_bases_removed_in_make: 0,
        }
    }

    /// `CigarBuilder.add`: append one element, applying the rewrites in the module note.
    ///
    /// The order of the checks is the reference's and decides the result: a leading deletion is
    /// dropped **before** the section is advanced, so it never moves the builder out of the left
    /// clips; then merging, then the right-clip and deletion-after-insertion rewrites.
    pub fn add(&mut self, element: CigarElement) -> Result<(), CigarError> {
        if element.length == 0 {
            return Ok(());
        }
        let operator = element.op;

        // A deletion at the start of the read is dropped, including the edge case of a deletion
        // that follows a *leading insertion*: `10S 2I 5D 5M` keeps no deletion at all.
        if self.remove_deletions_at_ends && operator == Op::D {
            let leading = match self.last_operator {
                None => true,
                Some(last) if is_clipping(last) => true,
                // An insertion that is itself leading: the only element, or right after a clip.
                // `wrapping_sub(2)` only keeps the index arithmetic from panicking on a short list;
                // the `len() == 1` test before it means it is never asked about one.
                Some(Op::I) => {
                    self.elements.len() == 1
                        || self
                            .elements
                            .get(self.elements.len().wrapping_sub(2))
                            .is_some_and(|e| is_clipping(e.op))
                }
                _ => false,
            };
            if leading {
                self.leading_deletion_bases_removed += element.length;
                return Ok(());
            }
        }

        self.advance_section_and_validate_cigar_order(operator)?;

        // `operator == lastOperator`: merge into the last element. After a deletion was moved
        // behind an insertion the last element is that insertion, and so is the last operator.
        if Some(operator) == self.last_operator {
            let last = self.elements.len() - 1;
            self.elements[last].length += element.length;
            return Ok(());
        }

        match self.last_operator {
            None => {
                self.elements.push(element);
                self.last_operator = Some(operator);
            }
            Some(last) if is_clipping(operator) => {
                // Clipping has just started on the right: a deletion sitting immediately before it
                // is meaningless and is replaced, and so is a deletion hiding behind an insertion.
                if self.remove_deletions_at_ends
                    && !last.consumes_read_bases()
                    && !is_clipping(last)
                {
                    let n = self.elements.len() - 1;
                    self.trailing_deletion_bases_removed += self.elements[n].length;
                    self.elements[n] = element;
                    self.last_operator = Some(operator);
                } else if self.remove_deletions_at_ends
                    && self.last_two_elements_were_deletion_and_insertion()
                {
                    let n = self.elements.len();
                    self.trailing_deletion_bases_removed += self.elements[n - 2].length;
                    self.elements[n - 2] = self.elements[n - 1];
                    self.elements[n - 1] = element;
                } else {
                    self.elements.push(element);
                    self.last_operator = Some(operator);
                }
            }
            Some(Op::I) if operator == Op::D => {
                // Deletions shift left past an insertion, and merge into a deletion already there.
                // Note the last operator stays an insertion: the deletion went *behind* it.
                let size = self.elements.len();
                if size > 1 && self.elements[size - 2].op == Op::D {
                    self.elements[size - 2].length += element.length;
                } else {
                    self.elements.insert(size - 1, element);
                }
            }
            Some(_) => {
                self.elements.push(element);
                self.last_operator = Some(operator);
            }
        }
        Ok(())
    }

    /// `lastTwoElementsWereDeletionAndInsertion`: the list ends `..., D, I`, which is what a
    /// deletion added after an insertion leaves behind.
    fn last_two_elements_were_deletion_and_insertion(&self) -> bool {
        self.last_operator == Some(Op::I)
            && self.elements.len() > 1
            && self.elements[self.elements.len() - 2].op == Op::D
    }

    /// `CigarBuilder.make(allowEmpty)`: drop a trailing deletion (or one hidden behind a trailing
    /// insertion), record its length, and hand back the cigar.
    ///
    /// The removal mutates the builder, as the reference's does. The returned cigar is a copy
    /// (`clone`), so the builder can still be asked for its counters afterwards.
    pub fn make(&mut self, allow_empty: bool) -> Result<Cigar, CigarError> {
        if self.section == Section::LeftSoftClip
            && self.elements.first().is_some_and(|e| e.op == Op::S)
        {
            return Err(CigarError::CompletelySoftClipped);
        }
        self.trailing_deletion_bases_removed_in_make = 0;
        if self.remove_deletions_at_ends && self.last_operator == Some(Op::D) {
            let last = self.elements.len() - 1;
            self.trailing_deletion_bases_removed_in_make = self.elements[last].length;
            self.elements.remove(last);
        } else if self.remove_deletions_at_ends
            && self.last_two_elements_were_deletion_and_insertion()
        {
            let n = self.elements.len() - 2;
            self.trailing_deletion_bases_removed_in_make = self.elements[n].length;
            self.elements.remove(n);
        }
        if !allow_empty && self.elements.is_empty() {
            return Err(CigarError::Empty);
        }
        Ok(Cigar::new(self.elements.clone()))
    }

    /// `CigarBuilder.getLeadingDeletionBasesRemoved`.
    pub fn leading_deletion_bases_removed(&self) -> u32 {
        self.leading_deletion_bases_removed
    }

    /// `CigarBuilder.getTrailingDeletionBasesRemoved`: both the ones removed while adding and the
    /// one removed by `make`.
    pub fn trailing_deletion_bases_removed(&self) -> u32 {
        self.trailing_deletion_bases_removed + self.trailing_deletion_bases_removed_in_make
    }

    /// `advanceSectionAndValidateCigarOrder`: move to the section `operator` belongs in, or refuse.
    ///
    /// A hard clip either leaves the builder in the left hard clip or moves it to the right one,
    /// so `5H5H` and `5H5S` are both accepted while `5M5H5S` is not.
    fn advance_section_and_validate_cigar_order(&mut self, operator: Op) -> Result<(), CigarError> {
        match operator {
            Op::H => {
                if matches!(
                    self.section,
                    Section::LeftSoftClip | Section::Middle | Section::RightSoftClip
                ) {
                    self.section = Section::RightHardClip;
                }
            }
            Op::S => {
                if self.section == Section::RightHardClip {
                    return Err(CigarError::SectionOutOfOrder);
                }
                if self.section == Section::LeftHardClip {
                    self.section = Section::LeftSoftClip;
                } else if self.section == Section::Middle {
                    self.section = Section::RightSoftClip;
                }
            }
            _ => {
                if matches!(
                    self.section,
                    Section::RightSoftClip | Section::RightHardClip
                ) {
                    return Err(CigarError::SectionOutOfOrder);
                }
                if matches!(self.section, Section::LeftHardClip | Section::LeftSoftClip) {
                    self.section = Section::Middle;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build(elements: &[(u32, Op)]) -> Result<String, CigarError> {
        let mut builder = CigarBuilder::default();
        for (length, op) in elements {
            builder.add(CigarElement {
                length: *length,
                op: *op,
            })?;
        }
        Ok(builder.make(false)?.to_text())
    }

    #[test]
    fn consecutive_identical_operators_merge() {
        assert_eq!(build(&[(3, Op::M), (4, Op::M)]).unwrap(), "7M");
        // A zero-length element between them does not stop the merge.
        assert_eq!(build(&[(3, Op::M), (0, Op::M), (4, Op::M)]).unwrap(), "7M");
    }

    #[test]
    fn a_deletion_after_an_insertion_moves_before_it() {
        assert_eq!(
            build(&[(5, Op::M), (2, Op::I), (3, Op::D), (5, Op::M)]).unwrap(),
            "5M3D2I5M"
        );
        // And merges into a deletion that is already there.
        assert_eq!(
            build(&[(5, Op::M), (3, Op::D), (2, Op::I), (2, Op::D), (5, Op::M)]).unwrap(),
            "5M5D2I5M"
        );
    }

    #[test]
    fn deletions_at_either_end_are_removed_and_counted() {
        let mut builder = CigarBuilder::default();
        for (length, op) in [(10, Op::S), (5, Op::D), (5, Op::M)] {
            builder.add(CigarElement { length, op }).unwrap();
        }
        assert_eq!(builder.make(false).unwrap().to_text(), "10S5M");
        assert_eq!(builder.leading_deletion_bases_removed(), 5);

        let mut builder = CigarBuilder::default();
        for (length, op) in [(5, Op::M), (5, Op::D)] {
            builder.add(CigarElement { length, op }).unwrap();
        }
        assert_eq!(builder.make(false).unwrap().to_text(), "5M");
        assert_eq!(builder.trailing_deletion_bases_removed(), 5);
    }

    #[test]
    fn a_completely_soft_clipped_cigar_is_refused() {
        assert_eq!(
            build(&[(10, Op::S)]),
            Err(CigarError::CompletelySoftClipped)
        );
    }

    #[test]
    fn a_soft_clip_in_the_middle_is_refused() {
        assert_eq!(
            build(&[(5, Op::M), (3, Op::S), (5, Op::M)]),
            Err(CigarError::SectionOutOfOrder)
        );
    }
}
