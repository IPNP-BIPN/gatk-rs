//! `FlowAnnotatorBase` and its eight members, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.flow` (GATK 4.7.0.0; every file in that
//! package is unchanged since 4.6.2.0), with `FlowBasedKeyCodec.baseArrayToKey`.
//!
//! `VARIANT_TYPE`, `INDEL_CLASSIFY`, `INDEL_LENGTH`, `HMER_INDEL_LENGTH`, `HMER_INDEL_NUC`,
//! `LEFT_MOTIF`/`RIGHT_MOTIF`, `GC_CONTENT` and `CYCLESKIP_STATUS`: the annotations for flow-based
//! sequencing, which describe a variant in **flow space** rather than in base space.
//!
//! # A flow key is a run-length encoding against a repeating flow order
//!
//! ```java
//! final byte flowBase = flowOrderBytes[flowNumber%period];
//! if ((bases[loc]!=flowBase) && ( bases[loc]!= BaseUtils.Base.N.base)) {
//!     result.add(0);
//!     if ( ++periodGuard > period )
//! ```
//!
//! The order cycles (`TGCA` by default), and each flow emits the number of bases it read, zero
//! included. `N` matches **every** flow base, so an ambiguous base is absorbed into whatever run it
//! sits in rather than breaking it.
//!
//! And the guard: more consecutive zeroes than the period means a base the flow order does not
//! contain, which throws rather than looping. That is the only way this function fails.
//!
//! # An hmer indel is one whose flow keys differ in exactly one position
//!
//! The reference and the alternate are each padded with one base before and a homopolymer run
//! plus five bases after, converted to flow keys, and compared. Different lengths, or more than
//! one differing flow, means it is not an hmer indel. The length reported is the **larger** of the
//! two flows, and the nucleotide is the flow order's base at that index.
//!
//! # `VARIANT_TYPE` asks three questions in order
//!
//! `snp` if every alternate is the same length as the reference; otherwise `h-indel` if **every**
//! alternate has a non-zero hmer length; otherwise `non-h-indel`. A site with one hmer indel and
//! one that is not is `non-h-indel`, and a spanning deletion or `<NON_REF>` is skipped in both
//! loops rather than deciding either.
//!
//! # The motifs are five bases, and the left one shifts for an indel
//!
//! ```java
//! if ( a.length() != refLength ) {
//!     motif = motif.substring(1) + vc.getReference().getBaseString().substring(0, 1);
//! }
//! ```
//!
//! For an indel the left motif drops its first base and appends the reference's first, so the two
//! motifs of a mixed site are taken from different windows. A motif that would run off the contig
//! is not truncated: the annotation is **dropped entirely**, and the reference logs it once.
//!
//! # Where the port differs
//!
//! Each of these is noted again at the line it concerns:
//!
//! - `isHmerIndel` throws where the port skips the allele: a `before` base outside the window
//!   (`Utils.validIndex` in `getReferenceNucleotide`) and a flow key whose period guard trips.
//! - When `getReferenceHmerPlus` answers empty, `isHmerIndel` returns before it puts
//!   `HMER_INDEL_LENGTH` and `HMER_INDEL_NUC` and before it sets the local context's lists. The
//!   port returns the lists built so far.
//! - The reference tests `Allele.length()`, which is zero for a symbolic allele, in the left motif,
//!   the right motif and `cycleSkip`. The port tests the length of the allele's text there, which
//!   differs only for a symbolic allele that is not `*` or `<NON_REF>`, or for `<NON_REF>` beside
//!   a nine-base reference in the left motif.
//!
//! # Reading this file without reading Rust
//!
//! - `&[u8]` is a borrowed, read-only view of a `byte[]`; `&str` is the same for a `String`.
//!   `.as_bytes()` views a string's bytes without copying, as `getBytes()` would with a copy.
//! - `Option<T>` is "a `T` or nothing". Here it stands for three Java things: a `null` list entry
//!   (`Option<String>` inside a `Vec`), an attribute the reference never puts (an `Option<Vec<..>>`
//!   field that is `None`), and an exception (`base_array_to_key`'s `None`).
//! - `let Some(x) = expr else { continue; };` unwraps an `Option` or takes the exit in the block.
//! - `usize` is an unsigned index type, so a position that may be negative is kept in an `i64` and
//!   tested before `as usize` converts it, because a negative value would wrap to a huge index.
//! - `Window<'a>` holds a borrowed slice; `'a` names how long the borrow lasts and has no Java
//!   counterpart: the window cannot outlive the bases it points into.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

/// `FlowBasedRead.DEFAULT_FLOW_ORDER`.
pub const DEFAULT_FLOW_ORDER: &str = "TGCA";

/// The constants `FlowAnnotatorBase` names, with the reference's names and values. They are the
/// text written into the VCF.
pub const C_INSERT: &str = "ins";
pub const C_DELETE: &str = "del";
pub const C_NA: &str = "NA";
pub const C_CSS_CS: &str = "cycle-skip";
pub const C_CSS_PCS: &str = "possible-cycle-skip";
pub const C_CSS_NS: &str = "non-skip";
pub const C_SNP: &str = "snp";
pub const C_NON_H_MER: &str = "non-h-indel";
pub const C_H_MER: &str = "h-indel";

/// `FlowAnnotatorBase.MOTIF_SIZE`: the length of each motif, and the bases after an hmer.
const MOTIF_SIZE: usize = 5;
/// `FlowAnnotatorBase.GC_CONTENT_SIZE`: the window `GC_CONTENT` is measured over.
const GC_CONTENT_SIZE: usize = 10;

/// `FlowBasedKeyCodec.baseArrayToKey`.
///
/// **What**: the bases as a flow key, one count per flow of the cycling `flow_order`.
///
/// **How**: walk the bases once. At each flow, a base that is neither the flow's base nor `N`
/// emits a zero and moves to the next flow; otherwise the whole run of matching bases (and `N`s)
/// is consumed and its length emitted.
///
/// `None` is the reference's `GATKException`: the period guard tripped, which means a base the
/// flow order does not contain. The reference never returns `null`. The golden's `key X TGCA` row
/// (`flow.txt.gz`) is that error.
pub fn base_array_to_key(bases: &[u8], flow_order: &str) -> Option<Vec<i32>> {
    let flow = flow_order.as_bytes();
    let period = flow.len();
    let mut result: Vec<i32> = Vec::new();
    let mut loc = 0usize;
    let mut flow_number = 0usize;
    let mut period_guard = 0usize;
    while loc < bases.len() {
        let flow_base = flow[flow_number % period];
        // `N` matches every flow base, so it is absorbed rather than breaking a run.
        if bases[loc] != flow_base && bases[loc] != b'N' {
            result.push(0);
            // `++periodGuard > period`: increment first, then compare, so a full period of zeroes
            // is allowed and one more is not.
            period_guard += 1;
            if period_guard > period {
                return None;
            }
        } else {
            let mut count = 0;
            while loc < bases.len() && (bases[loc] == flow_base || bases[loc] == b'N') {
                loc += 1;
                count += 1;
            }
            result.push(count);
            period_guard = 0;
        }
        flow_number += 1;
    }
    Some(result)
}

/// `isSpecial`: the spanning deletion and `<NON_REF>`, which every loop here skips.
///
/// ```java
/// return a.equals(Allele.SPAN_DEL) || a.equals(Allele.NON_REF_ALLELE);
/// ```
///
/// Compared by text here. `Allele.equals` also compares the reference flag, but neither constant
/// is a reference allele and neither text can be one.
fn is_special(allele: &Allele) -> bool {
    let bases = allele.display_string();
    bases == "*" || bases == "<NON_REF>"
}

/// `Allele.length()`, which is **zero** for a symbolic allele rather than the length of its text.
///
/// So `<NON_REF>` is shorter than any reference and classifies as a deletion, even though it is
/// skipped everywhere the classification is then read. The golden's `non-ref` row is that.
fn allele_length(allele: &Allele) -> usize {
    let bases = allele.display_string();
    if bases.starts_with('<') && bases.ends_with('>') {
        0
    } else {
        bases.len()
    }
}

/// The reference window a flow annotation reads from: `ref.getBases()` and
/// `ref.getWindow().getStart()`, the 1-based position of `bases[0]`.
pub struct Window<'a> {
    pub start: i64,
    pub bases: &'a [u8],
}

impl Window<'_> {
    /// `getReferenceNucleotide`, which does **not** catch: an index outside the window is a
    /// programming error to the reference and a `None` here.
    ///
    /// ```java
    /// Utils.validIndex(index, bases.length);
    /// ```
    ///
    /// The Java comment on that line begins `do not catch`. `Utils.validIndex` throws
    /// `IllegalArgumentException`, which the reference lets escape.
    fn nucleotide(&self, position: i64) -> Option<u8> {
        let index = position - self.start;
        if index < 0 || index as usize >= self.bases.len() {
            return None;
        }
        Some(self.bases[index as usize])
    }

    /// `getReferenceHmerPlus`: the homopolymer run starting at `position`, plus `additional` more
    /// bases, truncated at the window's end rather than refused.
    ///
    /// Empty only when `position` itself is outside the window, which is the reference's caught
    /// `validIndex`: `return new byte[0];`.
    fn hmer_plus(&self, position: i64, additional: usize) -> Vec<u8> {
        let mut index = position - self.start;
        if index < 0 || index as usize >= self.bases.len() {
            return Vec::new();
        }
        let mut out = Vec::new();
        let base0 = self.bases[index as usize];
        out.push(base0);
        index += 1;
        while (index as usize) < self.bases.len() && self.bases[index as usize] == base0 {
            out.push(self.bases[index as usize]);
            index += 1;
        }
        for _ in 0..additional {
            if index as usize >= self.bases.len() {
                break;
            }
            out.push(self.bases[index as usize]);
            index += 1;
        }
        out
    }

    /// `getRefMotif`, which answers the empty string rather than truncating: both `startIndex`
    /// and `endIndex-1` go through `Utils.validIndex`, and either failure is caught as `""`.
    fn motif(&self, start: i64, length: usize) -> String {
        let start_index = start - self.start;
        if start_index < 0 {
            return String::new();
        }
        let end_index = start_index as usize + length;
        if start_index as usize >= self.bases.len() || end_index > self.bases.len() {
            return String::new();
        }
        String::from_utf8_lossy(&self.bases[start_index as usize..end_index]).into_owned()
    }
}

/// Everything the eight annotations compute, in one pass, exactly as the base class does.
///
/// The reference computes these per annotation, each member calling the base class steps it
/// needs into a fresh `LocalContext`. The steps are deterministic and read only the variant, the
/// window and the flow order, so computing all of them once gives each member the values it would
/// compute itself; the caller picks the fields a member writes.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowAnnotations {
    /// `FLOW_INDEL_CLASSIFY`: `ins`, `del` or `NA`, one per alternate.
    pub indel: Vec<String>,
    /// `FLOW_INDEL_LENGTH`, `None` where the reference puts a null.
    pub indel_length: Vec<Option<i32>>,
    /// `FLOW_HMER_INDEL_LENGTH`.
    pub hmer_indel_length: Vec<i32>,
    /// `FLOW_HMER_INDEL_NUC`.
    pub hmer_indel_nuc: Vec<Option<String>>,
    /// `FLOW_LEFT_MOTIF`, absent when the window could not supply five bases.
    pub left_motif: Option<Vec<String>>,
    /// `FLOW_RIGHT_MOTIF`.
    pub right_motif: Option<Vec<Option<String>>>,
    /// `FLOW_GC_CONTENT`, a **float** rather than a double.
    pub gc_content: Option<f32>,
    /// `FLOW_CYCLESKIP_STATUS`, or the refusal a missing motif produces.
    pub cycle_skip: Result<Vec<String>, CycleSkipError>,
    /// `FLOW_VARIANT_TYPE`.
    pub variant_type: Option<String>,
}

/// `indelClassify`: `ins`, `del` or `NA` per alternate, and the length difference.
///
/// ```java
/// indelClassify.add(refLength == a.length() ? C_NA
///         : (refLength < a.length() ? C_INSERT : C_DELETE));
/// ```
///
/// (One line in `FlowAnnotatorBase.java`, wrapped here.)
///
/// Lengths are `Allele.length()` ([`allele_length`]), so `<NON_REF>` is a `del`. The length list
/// holds a null for a special allele and for a same-length one.
fn indel_classify(vc: &VariantContext) -> (Vec<String>, Vec<Option<i32>>) {
    let ref_length = vc
        .alleles
        .iter()
        .find(|a| a.is_reference())
        .map(allele_length)
        .unwrap_or(0);
    let mut classify = Vec::new();
    let mut length = Vec::new();
    for allele in vc.alleles.iter().filter(|a| !a.is_reference()) {
        let allele_length = allele_length(allele);
        classify.push(
            if ref_length == allele_length {
                C_NA
            } else if ref_length < allele_length {
                C_INSERT
            } else {
                C_DELETE
            }
            .to_string(),
        );
        if !is_special(allele) && allele_length != ref_length {
            // `Math.abs(refLength - a.length())`, in `int` arithmetic.
            length.push(Some((ref_length as i32 - allele_length as i32).abs()));
        } else {
            // A null in the list, which the encoder writes as a missing value rather than as 0.
            length.push(None);
        }
    }
    (classify, length)
}

/// `isHmerIndel`, which also produces the right motif for an indel.
///
/// **What**: per alternate, the hmer length (0 if not an hmer indel), its nucleotide (`None` if
/// not), and the right motif for an hmer indel whose alleles differ in length (`None` otherwise,
/// filled in later by [`right_motif`]).
///
/// **How**: pad both alleles with the base before the variant and the hmer plus five bases after
/// its end, convert both to flow keys, and accept only keys of equal length that differ in exactly
/// one flow, where that flow is not zero on both sides.
///
/// The returned tuple is the reference's three `LinkedList`s, `hmerIndelLength`, `hmerIndelNuc`
/// and `rightMotif`, which it fills by appending a default and then replacing the last element.
#[allow(clippy::type_complexity)]
fn is_hmer_indel(
    vc: &VariantContext,
    window: &Window,
    flow_order: &str,
) -> (Vec<i32>, Vec<Option<String>>, Vec<Option<String>>) {
    let mut hmer_length = Vec::new();
    let mut hmer_nuc = Vec::new();
    let mut right_motif = Vec::new();
    let Some(reference) = vc.alleles.iter().find(|a| a.is_reference()) else {
        return (hmer_length, hmer_nuc, right_motif);
    };

    for allele in vc.alleles.iter().filter(|a| !a.is_reference()) {
        hmer_length.push(0);
        hmer_nuc.push(None);
        right_motif.push(None);
        if is_special(allele) {
            continue;
        }

        // The reference throws here (`Utils.validIndex`); the port skips the allele, leaving its
        // defaults.
        let Some(before) = window.nucleotide(vc.start - 1) else {
            continue;
        };
        let after = window.hmer_plus(vc.stop + 1, MOTIF_SIZE);
        if after.is_empty() {
            // "probably because the variant is very close to the end of the chromosome". The
            // reference's `return;` comes before it puts `FLOW_HMER_INDEL_LENGTH` and
            // `FLOW_HMER_INDEL_NUC` and before it sets `localContext.rightMotif`, so it writes
            // neither attribute and `cycleSkip` later finds a null list. The port returns the lists
            // built so far, which the caller writes. The golden has no row for this case.
            return (hmer_length, hmer_nuc, right_motif);
        }

        // `buildHaplotype(before, bases, after)` as a closure: `|bases: &[u8]| { .. }` is a Java
        // lambda that also reads `before` and `after` from the enclosing scope.
        let build = |bases: &[u8]| {
            let mut hap = Vec::with_capacity(1 + bases.len() + after.len());
            hap.push(before);
            hap.extend_from_slice(bases);
            hap.extend_from_slice(&after);
            hap
        };
        let ref_hap = build(reference.display_string().as_bytes());
        let alt_hap = build(allele.display_string().as_bytes());

        // A tripped period guard is a `GATKException` thrown out of `baseArrayToKey` in the
        // reference (the `refKey == null` test after it can never fire); the port skips the
        // allele instead.
        let (Some(ref_key), Some(alt_key)) = (
            base_array_to_key(&ref_hap, flow_order),
            base_array_to_key(&alt_hap, flow_order),
        ) else {
            continue;
        };
        if ref_key.len() != alt_key.len() {
            continue;
        }

        // `diffIndex` stays a signed `int` with -1 for "none yet", as in the reference, rather than
        // an `Option`, so the loop reads like the Java. `refBasesCountUpInclHmer` sums the
        // reference key up to and including the differing flow: the offset of the right motif.
        let mut diff_index: i32 = -1;
        let mut ref_bases_up_incl_hmer = 0i32;
        for n in 0..ref_key.len() {
            if diff_index < 0 {
                ref_bases_up_incl_hmer += ref_key[n];
            }
            if ref_key[n] != alt_key[n] {
                if diff_index >= 0 {
                    // A second differing flow: not an hmer indel at all.
                    diff_index = -1;
                    break;
                }
                diff_index = n as i32;
            }
        }
        if diff_index < 0 {
            continue;
        }
        let index = diff_index as usize;
        // `Math.max(refKey[diffIndex], altKey[diffIndex]) == 0`; `.max` on two integers is the
        // same function.
        if ref_key[index].max(alt_key[index]) == 0 {
            continue;
        }

        let length = ref_key[index].max(alt_key[index]);
        let nuc = flow_order.as_bytes()[index % flow_order.len()];
        // "replace last element of list": the slot this allele appended above.
        let last = hmer_length.len() - 1;
        hmer_length[last] = length;
        hmer_nuc[last] = Some((nuc as char).to_string());

        // `a.length() != ref.length()`, measured here on the alleles' text. Special alleles were
        // skipped above, and any other symbolic allele has already failed at its flow key.
        if allele.display_string().len() != reference.display_string().len() {
            let from = ref_bases_up_incl_hmer as usize;
            let to = (from + MOTIF_SIZE).min(ref_hap.len());
            right_motif[last] = Some(String::from_utf8_lossy(&ref_hap[from..to]).into_owned());
        }
    }
    (hmer_length, hmer_nuc, right_motif)
}

/// `variantType`, which asks three questions in order.
///
/// `indel` and `hmer_length` are indexed by the alternate's position, as the reference indexes
/// `localContext.indel.get(i)` over `vc.getAlternateAlleles()`. The reference also treats a null
/// `hmerIndelLength` list as non-hmer; an empty slice, or a short one, has the same effect here
/// because a missing entry fails the match below.
fn variant_type(vc: &VariantContext, indel: &[String], hmer_length: &[i32]) -> String {
    let alternates: Vec<&Allele> = vc.alleles.iter().filter(|a| !a.is_reference()).collect();

    let mut is_snp = true;
    for (index, allele) in alternates.iter().enumerate() {
        if is_special(allele) {
            continue;
        }
        if indel.get(index).map(|c| c != C_NA).unwrap_or(false) {
            is_snp = false;
        }
    }
    if is_snp {
        return C_SNP.to_string();
    }

    let mut is_hmer = true;
    for (index, allele) in alternates.iter().enumerate() {
        if is_special(allele) {
            continue;
        }
        match hmer_length.get(index) {
            Some(length) if *length != 0 => {}
            // A null or a zero makes the whole site non-hmer, so one ordinary indel beside an
            // hmer one decides for both.
            _ => is_hmer = false,
        }
    }
    if is_hmer {
        C_H_MER.to_string()
    } else {
        C_NON_H_MER.to_string()
    }
}

/// `getLeftMotif`, including the shift an indel applies.
///
/// `None` is the reference's early `return;` before `FLOW_LEFT_MOTIF` is put, which also leaves
/// `localContext.leftMotif` null for `cycleSkip`. The `?` after `find(..)` returns `None` at once
/// when there is no reference allele.
fn left_motif(vc: &VariantContext, window: &Window) -> Option<Vec<String>> {
    let reference = vc.alleles.iter().find(|a| a.is_reference())?;
    let ref_bases = reference.display_string();
    let ref_length = ref_bases.len();
    let mut motifs = Vec::new();
    for allele in vc.alleles.iter().filter(|a| !a.is_reference()) {
        let mut motif = window.motif(vc.start - MOTIF_SIZE as i64, MOTIF_SIZE);
        if motif.len() != MOTIF_SIZE {
            // Dropped entirely rather than truncated.
            return None;
        }
        // The reference tests `a.length()`, zero for a symbolic allele; the text length differs
        // from it only for a symbolic allele. `&motif[1..]` is `substring(1)` and `&ref_bases[..1]`
        // `substring(0, 1)`; `format!` concatenates.
        if allele.display_string().len() != ref_length {
            motif = format!("{}{}", &motif[1..], &ref_bases[..1]);
        }
        motifs.push(motif);
    }
    Some(motifs)
}

/// `getRightMotif`: the window's motif fills every slot the hmer pass left empty.
///
/// The motif starts right after the reference allele, at `vc.getStart() + refLength`, and is the
/// same for every allele; an hmer indel keeps the motif `isHmerIndel` cut from its own padded
/// haplotype.
fn right_motif(
    vc: &VariantContext,
    window: &Window,
    from_hmer: &[Option<String>],
) -> Option<Vec<Option<String>>> {
    let reference = vc.alleles.iter().find(|a| a.is_reference())?;
    let ref_length = reference.display_string().len() as i64;
    let motif = window.motif(vc.start + ref_length, MOTIF_SIZE);
    if motif.len() != MOTIF_SIZE {
        // The reference returns before putting RIGHT_MOTIF into the attribute map, but the list
        // itself was already built by `isHmerIndel` and stays in the local context. So the
        // annotation is absent while `cycleSkip` still sees a list, with nulls in it.
        return None;
    }
    Some(
        from_hmer
            .iter()
            // `if (rightMotif.get(i) == null) rightMotif.set(i, motif)`: `or_else` keeps a present
            // value and supplies the motif for a missing one.
            .map(|existing| existing.clone().or_else(|| Some(motif.clone())))
            .collect(),
    )
}

/// `gcContent`, which is a `float` and is dropped near a contig edge.
///
/// The ten bases start at `vc.getStart() - 5 + 1`, so the variant's own base is the fifth.
fn gc_content(vc: &VariantContext, window: &Window) -> Option<f32> {
    let begin = vc.start - (GC_CONTENT_SIZE as i64 / 2);
    let seq = window.motif(begin + 1, GC_CONTENT_SIZE);
    if seq.len() != GC_CONTENT_SIZE {
        return None;
    }
    let gc = seq.bytes().filter(|b| *b == b'G' || *b == b'C').count();
    // `(float) gcCount / seq.length()`: a float division, so the value carries single precision
    // into the VCF.
    Some(gc as f32 / seq.len() as f32)
}

/// What `cycleSkip` raises rather than answering, and they are two different exceptions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CycleSkipError {
    /// The motif **list** is null, which the reference dereferences: `NullPointerException`.
    MissingMotifList,
    /// A motif **entry** is null, which Java concatenates as the four characters `null`. Those are
    /// not bases, so the flow key's period guard trips: `GATKException`.
    KeyFromNullMotif,
}

/// `cycleSkip`.
///
/// **What**: per alternate, `NA` for a special allele or an indel; otherwise whether replacing the
/// reference with the alternate between the two motifs changes the flow key's length
/// (`cycle-skip`), empties or fills a flow (`possible-cycle-skip`), or neither (`non-skip`).
///
/// The motif lists are dereferenced without a guard, so a variant whose left motif ran off the
/// window makes this throw a `NullPointerException` **before** the annotation gets the chance to
/// decline for want of a flow order. The golden's `CycleSkipStatus near-window-start` row
/// (`flow.txt.gz`) is that. Its `near-window-end` row is the other failure: a right motif that ran
/// off the window leaves a null entry, which becomes the text `null`, and the flow key throws a
/// `GATKException`.
fn cycle_skip(
    vc: &VariantContext,
    flow_order: &str,
    left: &Option<Vec<String>>,
    right: &Option<Vec<Option<String>>>,
) -> Result<Vec<String>, CycleSkipError> {
    let mut css = Vec::new();
    let Some(reference) = vc.alleles.iter().find(|a| a.is_reference()) else {
        return Ok(css);
    };
    let ref_bases = reference.display_string();
    let ref_length = ref_bases.len();

    for allele in vc.alleles.iter().filter(|a| !a.is_reference()) {
        if is_special(allele) || allele.display_string().len() != ref_length {
            css.push(C_NA.to_string());
            continue;
        }
        // `final int i = css.size();   // always working on the last`.
        let index = css.len();
        // `localContext.leftMotif.get(i)` with a null list: the reference throws here.
        let (Some(left), Some(right)) = (left.as_ref(), right.as_ref()) else {
            return Err(CycleSkipError::MissingMotifList);
        };
        // An index past the end of either list is an `IndexOutOfBoundsException` in the
        // reference; it is reported as the same refusal as a null list.
        let (Some(left_motif), Some(right_entry)) = (left.get(index), right.get(index)) else {
            return Err(CycleSkipError::MissingMotifList);
        };
        // `left + ref + right` with a null right motif is the string "null" in Java, not an empty
        // string, and those four characters are not bases.
        let right_motif = match right_entry {
            Some(motif) => motif.clone(),
            None => "null".to_string(),
        };
        let ref_key = base_array_to_key(
            format!("{left_motif}{ref_bases}{right_motif}").as_bytes(),
            flow_order,
        );
        let alt_key = base_array_to_key(
            format!("{left_motif}{}{right_motif}", allele.display_string()).as_bytes(),
            flow_order,
        );
        let (Some(ref_key), Some(alt_key)) = (ref_key, alt_key) else {
            return Err(CycleSkipError::KeyFromNullMotif);
        };
        // The reference compares `cssValue == C_CSS_NS` by reference identity, which holds because
        // both are the same constant; comparing text here gives the same answer.
        let mut value = if ref_key.len() != alt_key.len() {
            C_CSS_CS
        } else {
            C_CSS_NS
        };
        if value == C_CSS_NS {
            for n in 0..ref_key.len() {
                // An exclusive or on "this flow is empty": one side reads nothing where the other
                // reads something, which is a *possible* cycle skip rather than a certain one.
                if (ref_key[n] == 0) ^ (alt_key[n] == 0) {
                    value = C_CSS_PCS;
                    break;
                }
            }
        }
        css.push(value.to_string());
    }
    Ok(css)
}

/// Everything the eight annotations report, computed in the order the base class computes it.
///
/// `flow_order` is the one `establishFlowOrder` settles on; the caller passes the default `TGCA`
/// when none was given, and decides whether an annotation that requires an actual flow order
/// (`CycleSkipStatus`) is written at all.
pub fn annotate(vc: &VariantContext, window: &Window, flow_order: &str) -> FlowAnnotations {
    let (indel, indel_length) = indel_classify(vc);
    let (hmer_indel_length, hmer_indel_nuc, hmer_right) = is_hmer_indel(vc, window, flow_order);
    let left = left_motif(vc, window);
    let right = right_motif(vc, window, &hmer_right);
    // What `cycleSkip` sees is the local context's list, which survives even when the attribute
    // was not written. That holds when `getRightMotif` returned early; when `isHmerIndel` did,
    // the reference's list is null and `cycleSkip` would throw `NullPointerException`.
    let right_for_cycle_skip = right.clone().or(Some(hmer_right.clone()));
    let gc = gc_content(vc, window);
    let css = cycle_skip(vc, flow_order, &left, &right_for_cycle_skip);
    let variant = variant_type(vc, &indel, &hmer_indel_length);

    FlowAnnotations {
        indel,
        indel_length,
        hmer_indel_length,
        hmer_indel_nuc,
        left_motif: left,
        right_motif: right,
        gc_content: gc,
        cycle_skip: css,
        variant_type: Some(variant),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_flow_key_counts_runs_against_a_cycling_order() {
        // TGCA: T reads one, G reads none, C reads none, A reads two.
        assert_eq!(base_array_to_key(b"TAA", "TGCA"), Some(vec![1, 0, 0, 2]));
        // N matches whatever flow it lands in.
        assert_eq!(base_array_to_key(b"TNA", "TGCA"), Some(vec![2, 0, 0, 1]));
    }

    #[test]
    fn a_base_outside_the_flow_order_trips_the_period_guard() {
        assert_eq!(base_array_to_key(b"X", "TGCA"), None);
    }
}
