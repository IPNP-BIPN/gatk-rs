//! `FlowAnnotatorBase` and its eight members, ported from
//! `org.broadinstitute.hellbender.tools.walkers.annotator.flow` (GATK 4.6.2.0).
//!
//! `VARIANT_TYPE`, `INDEL_CLASSIFY`, `INDEL_LENGTH`, `HMER_INDEL_LENGTH`, `HMER_INDEL_NUC`,
//! `LEFT_MOTIF`/`RIGHT_MOTIF`, `GC_CONTENT` and `CYCLESKIP_STATUS`: the annotations for flow-based
//! sequencing, which describe a variant in **flow space** rather than in base space.
//!
//! # A flow key is a run-length encoding against a repeating flow order
//!
//! ```java
//! final byte flowBase = flowOrderBytes[flowNumber % period];
//! if ((bases[loc] != flowBase) && (bases[loc] != 'N')) { result.add(0); ... }
//! else { count consecutive matching bases; result.add(count); }
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
//! if (a.length() != refLength) { motif = motif.substring(1) + ref.getBaseString().substring(0, 1); }
//! ```
//!
//! For an indel the left motif drops its first base and appends the reference's first, so the two
//! motifs of a mixed site are taken from different windows. A motif that would run off the contig
//! is not truncated: the annotation is **dropped entirely**, and the reference logs it once.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;

/// `FlowBasedRead.DEFAULT_FLOW_ORDER`.
pub const DEFAULT_FLOW_ORDER: &str = "TGCA";

/// The constants `FlowAnnotatorBase` names.
pub const C_INSERT: &str = "ins";
pub const C_DELETE: &str = "del";
pub const C_NA: &str = "NA";
pub const C_CSS_CS: &str = "cycle-skip";
pub const C_CSS_PCS: &str = "possible-cycle-skip";
pub const C_CSS_NS: &str = "non-skip";
pub const C_SNP: &str = "snp";
pub const C_NON_H_MER: &str = "non-h-indel";
pub const C_H_MER: &str = "h-indel";

const MOTIF_SIZE: usize = 5;
const GC_CONTENT_SIZE: usize = 10;

/// `FlowBasedKeyCodec.baseArrayToKey`.
///
/// `None` is the reference's `GATKException`: the period guard tripped, which means a base the
/// flow order does not contain.
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
fn is_special(allele: &Allele) -> bool {
    let bases = allele.display_string();
    bases == "*" || bases == "<NON_REF>"
}

/// The reference window a flow annotation reads from.
pub struct Window<'a> {
    pub start: i64,
    pub bases: &'a [u8],
}

impl Window<'_> {
    /// `getReferenceNucleotide`, which does **not** catch: an index outside the window is a
    /// programming error to the reference and a `None` here.
    fn nucleotide(&self, position: i64) -> Option<u8> {
        let index = position - self.start;
        if index < 0 || index as usize >= self.bases.len() {
            return None;
        }
        Some(self.bases[index as usize])
    }

    /// `getReferenceHmerPlus`: the homopolymer run starting at `position`, plus `additional` more
    /// bases, truncated at the window's end rather than refused.
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

    /// `getRefMotif`, which answers the empty string rather than truncating.
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
#[derive(Debug, Clone, Default, PartialEq)]
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
    /// `FLOW_CYCLESKIP_STATUS`.
    pub cycle_skip: Vec<String>,
    /// `FLOW_VARIANT_TYPE`.
    pub variant_type: Option<String>,
}

/// `indelClassify`.
fn indel_classify(vc: &VariantContext) -> (Vec<String>, Vec<Option<i32>>) {
    let ref_length = vc
        .alleles
        .iter()
        .find(|a| a.is_reference())
        .map(|a| a.display_string().len())
        .unwrap_or(0);
    let mut classify = Vec::new();
    let mut length = Vec::new();
    for allele in vc.alleles.iter().filter(|a| !a.is_reference()) {
        let allele_length = allele.display_string().len();
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
            length.push(Some((ref_length as i32 - allele_length as i32).abs()));
        } else {
            // A null in the list, which the encoder writes as a missing value rather than as 0.
            length.push(None);
        }
    }
    (classify, length)
}

/// `isHmerIndel`, which also produces the right motif for an indel.
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

        let Some(before) = window.nucleotide(vc.start - 1) else {
            continue;
        };
        let after = window.hmer_plus(vc.stop + 1, MOTIF_SIZE);
        if after.is_empty() {
            // "probably because the variant is very close to the end of the chromosome": the
            // whole annotation stops here, so the alleles after this one get nothing either.
            return (hmer_length, hmer_nuc, right_motif);
        }

        let build = |bases: &[u8]| {
            let mut hap = Vec::with_capacity(1 + bases.len() + after.len());
            hap.push(before);
            hap.extend_from_slice(bases);
            hap.extend_from_slice(&after);
            hap
        };
        let ref_hap = build(reference.display_string().as_bytes());
        let alt_hap = build(allele.display_string().as_bytes());

        let (Some(ref_key), Some(alt_key)) = (
            base_array_to_key(&ref_hap, flow_order),
            base_array_to_key(&alt_hap, flow_order),
        ) else {
            continue;
        };
        if ref_key.len() != alt_key.len() {
            continue;
        }

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
        if ref_key[index].max(alt_key[index]) == 0 {
            continue;
        }

        let length = ref_key[index].max(alt_key[index]);
        let nuc = flow_order.as_bytes()[index % flow_order.len()];
        let last = hmer_length.len() - 1;
        hmer_length[last] = length;
        hmer_nuc[last] = Some((nuc as char).to_string());

        if allele.display_string().len() != reference.display_string().len() {
            let from = ref_bases_up_incl_hmer as usize;
            let to = (from + MOTIF_SIZE).min(ref_hap.len());
            right_motif[last] = Some(String::from_utf8_lossy(&ref_hap[from..to]).into_owned());
        }
    }
    (hmer_length, hmer_nuc, right_motif)
}

/// `variantType`, which asks three questions in order.
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
        if allele.display_string().len() != ref_length {
            motif = format!("{}{}", &motif[1..], &ref_bases[..1]);
        }
        motifs.push(motif);
    }
    Some(motifs)
}

/// `getRightMotif`: the window's motif fills every slot the hmer pass left empty.
fn right_motif(
    vc: &VariantContext,
    window: &Window,
    from_hmer: &[Option<String>],
) -> Option<Vec<Option<String>>> {
    let reference = vc.alleles.iter().find(|a| a.is_reference())?;
    let ref_length = reference.display_string().len() as i64;
    let motif = window.motif(vc.start + ref_length, MOTIF_SIZE);
    if motif.len() != MOTIF_SIZE {
        return None;
    }
    Some(
        from_hmer
            .iter()
            .map(|existing| existing.clone().or_else(|| Some(motif.clone())))
            .collect(),
    )
}

/// `gcContent`, which is a `float` and is dropped near a contig edge.
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

/// `cycleSkip`.
fn cycle_skip(
    vc: &VariantContext,
    flow_order: &str,
    left: &Option<Vec<String>>,
    right: &Option<Vec<Option<String>>>,
) -> Vec<String> {
    let mut css = Vec::new();
    let Some(reference) = vc.alleles.iter().find(|a| a.is_reference()) else {
        return css;
    };
    let ref_bases = reference.display_string();
    let ref_length = ref_bases.len();

    for allele in vc.alleles.iter().filter(|a| !a.is_reference()) {
        if is_special(allele) || allele.display_string().len() != ref_length {
            css.push(C_NA.to_string());
            continue;
        }
        let index = css.len();
        let (Some(left), Some(right)) = (left.as_ref(), right.as_ref()) else {
            // The reference dereferences the motif lists here, so a missing motif is a
            // NullPointerException there. Here the status is simply not produced.
            css.push(C_NA.to_string());
            continue;
        };
        let (Some(left_motif), Some(Some(right_motif))) = (left.get(index), right.get(index))
        else {
            css.push(C_NA.to_string());
            continue;
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
            css.push(C_NA.to_string());
            continue;
        };
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
    css
}

/// Everything the eight annotations report, computed in the order the base class computes it.
pub fn annotate(vc: &VariantContext, window: &Window, flow_order: &str) -> FlowAnnotations {
    let (indel, indel_length) = indel_classify(vc);
    let (hmer_indel_length, hmer_indel_nuc, hmer_right) = is_hmer_indel(vc, window, flow_order);
    let left = left_motif(vc, window);
    let right = right_motif(vc, window, &hmer_right);
    let gc = gc_content(vc, window);
    let css = cycle_skip(vc, flow_order, &left, &right);
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
