//! Ported from `org.broadinstitute.hellbender.utils.read.CigarBuilder` (GATK 4.6.2.0).
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

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};

/// What the reference throws rather than returning.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CigarError {
    /// `advanceSectionAndValidateCigarOrder`: a soft clip after the right hard clip, or an aligned
    /// element after a right clip.
    SectionOutOfOrder,
    /// `make`: the cigar consists of a single soft clip.
    CompletelySoftClipped,
    /// `make(false)`: nothing left after removing flanking deletions.
    Empty,
}

/// `CigarBuilder.Section`: the order a cigar's parts must appear in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Section {
    LeftHardClip,
    LeftSoftClip,
    Middle,
    RightSoftClip,
    RightHardClip,
}

fn is_clipping(op: Op) -> bool {
    matches!(op, Op::S | Op::H)
}

/// `CigarBuilder`.
pub struct CigarBuilder {
    elements: Vec<CigarElement>,
    last_operator: Option<Op>,
    section: Section,
    remove_deletions_at_ends: bool,
    leading_deletion_bases_removed: u32,
    trailing_deletion_bases_removed: u32,
    trailing_deletion_bases_removed_in_make: u32,
}

impl Default for CigarBuilder {
    fn default() -> Self {
        CigarBuilder::new(true)
    }
}

impl CigarBuilder {
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

    /// `CigarBuilder.add`.
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

    fn last_two_elements_were_deletion_and_insertion(&self) -> bool {
        self.last_operator == Some(Op::I)
            && self.elements.len() > 1
            && self.elements[self.elements.len() - 2].op == Op::D
    }

    /// `CigarBuilder.make(allowEmpty)`.
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
