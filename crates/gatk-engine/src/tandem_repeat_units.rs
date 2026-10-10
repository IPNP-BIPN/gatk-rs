//! The tandem-repeat arithmetic of `GATKVariantContextUtils` (GATK 4.6.2.0): the repeat unit of an
//! indel and how many copies of it each allele carries in the reference context that follows.
//!
//! `TandemRepeat` (the `STR`, `RU` and `RPA` annotation) and `AssemblyRegionTrimmer` (the padding
//! of an indel in a short tandem repeat) both rest on it.
//!
//! # The counts are differences, so a repeat unit that already appears in the allele cancels
//!
//! ```java
//! int repetitionsInRef = findNumberOfRepetitions(repeatUnit, refBases, true);
//! repetitionCount[0] = findNumberOfRepetitions(repeatUnit, refBases + remainingRefContext, true) - repetitionsInRef;
//! repetitionCount[1] = findNumberOfRepetitions(repeatUnit, altBases + remainingRefContext, true) - repetitionsInRef;
//! ```
//!
//! Both counts subtract the repeats found in the **reference allele alone**, so `RPA` measures
//! copies in the surrounding context rather than in the allele. A count of zero on either side
//! means the allele is not a tandem expansion of its context, and the whole annotation is then
//! dropped for **every** allele, not only that one.
//!
//! # `findRepeatedSubstring` cannot see a partial trailing repeat, by construction
//!
//! ```java
//! final byte[] basePiece = Arrays.copyOfRange(bases,start,start+candidateRepeatUnit.length);
//! ```
//!
//! `copyOfRange` past the end of the array pads with zero bytes rather than throwing, and a zero
//! byte is not a base, so a candidate unit that divides the string only partially always fails the
//! comparison. `ACTACTAC` is therefore reported as a repeat of itself, length 8, not as `ACT`.
//!
//! And on an empty input the loop body never runs, so it returns 1 and the repeat unit is a single
//! **zero byte**. That is reachable: a multiallelic site where one alternate has the reference's
//! length leaves both base arrays empty after the padding base is stripped. The zero byte then
//! matches nothing, the count comes out zero, and the annotation is dropped, which is the right
//! answer arrived at by an accident of `Arrays.copyOf`.

use htsjdk_vcf::allele::Allele;

/// `GATKVariantContextUtils.findRepeatedSubstring`: the length of the shortest unit the string is a
/// whole number of copies of, or the string's own length when there is none.
pub fn find_repeated_substring(bases: &[u8]) -> usize {
    let mut rep_length = 1usize;
    while rep_length <= bases.len() {
        let candidate = &bases[..rep_length];
        let mut all_match = true;
        let mut start = rep_length;
        while start < bases.len() {
            // `Arrays.copyOfRange` pads with zero bytes past the end, so a short tail can never
            // compare equal to a unit made of bases.
            let matches = (0..rep_length).all(|i| {
                let byte = bases.get(start + i).copied().unwrap_or(0);
                byte == candidate[i]
            });
            if !matches {
                all_match = false;
                break;
            }
            start += rep_length;
        }
        if all_match {
            return rep_length;
        }
        rep_length += 1;
    }
    rep_length
}

/// `GATKVariantContextUtils.findNumberOfRepetitions`.
///
/// `None` is the reference's `IllegalArgumentException` on an empty repeat unit. An empty test
/// string is zero repetitions and not an error.
pub fn find_number_of_repetitions(
    repeat_unit: &[u8],
    test_string: &[u8],
    leading_repeats: bool,
) -> Option<i32> {
    if repeat_unit.is_empty() {
        return None;
    }
    if test_string.is_empty() {
        return Some(0);
    }
    let unit = repeat_unit.len() as i64;
    let length_difference = test_string.len() as i64 - unit;
    let mut repeats = 0i32;
    if leading_repeats {
        let mut start = 0i64;
        while start <= length_difference {
            if test_string[start as usize..(start + unit) as usize] == *repeat_unit {
                repeats += 1;
            } else {
                return Some(repeats);
            }
            start += unit;
        }
    } else {
        // Backwards, so `GATAT` has two trailing repeats of `AT` and no leading ones.
        let mut start = length_difference;
        while start >= 0 {
            if test_string[start as usize..(start + unit) as usize] == *repeat_unit {
                repeats += 1;
            } else {
                return Some(repeats);
            }
            start -= unit;
        }
    }
    Some(repeats)
}

/// `GATKVariantContextUtils.getNumTandemRepeatUnits(byte[], byte[], byte[])`.
///
/// Returns the `(reference count, alternate count)` pair and the repeat unit, both measured against
/// the context that follows the variant.
pub fn num_tandem_repeat_units_for_bases(
    ref_bases: &[u8],
    alt_bases: &[u8],
    remaining_ref_context: &[u8],
) -> Option<([i32; 2], Vec<u8>)> {
    // The longer of the two decides the unit, so an insertion of `ATAT` into `AT` is described as
    // `(AT)` rather than as `(ATAT)`.
    let long_b: &[u8] = if alt_bases.len() > ref_bases.len() {
        alt_bases
    } else {
        ref_bases
    };
    let repeat_unit_length = find_repeated_substring(long_b);
    // `Arrays.copyOf` pads, so an empty input gives a unit of one zero byte. See the module note.
    let repeat_unit: Vec<u8> = (0..repeat_unit_length)
        .map(|i| long_b.get(i).copied().unwrap_or(0))
        .collect();

    let repetitions_in_ref = find_number_of_repetitions(&repeat_unit, ref_bases, true)?;
    let mut ref_with_context = ref_bases.to_vec();
    ref_with_context.extend_from_slice(remaining_ref_context);
    let mut alt_with_context = alt_bases.to_vec();
    alt_with_context.extend_from_slice(remaining_ref_context);
    let counts = [
        find_number_of_repetitions(&repeat_unit, &ref_with_context, true)? - repetitions_in_ref,
        find_number_of_repetitions(&repeat_unit, &alt_with_context, true)? - repetitions_in_ref,
    ];
    Some((counts, repeat_unit))
}

/// `GATKVariantContextUtils.getNumTandemRepeatUnits(Allele, List<Allele>, byte[])`.
///
/// `None` covers three different reference outcomes: no alternate differs in length from the
/// reference (an explicit `null`), an allele that is not a tandem expansion of its context (also
/// `null`), and a symbolic alternate, whose zero length makes `Arrays.copyOfRange(bases, 1, 0)`
/// throw `IllegalArgumentException`. Only the last is an error rather than an answer, and it is
/// unreachable from `TandemRepeat`, whose `isIndel` guard rejects a symbolic site first.
pub fn num_tandem_repeat_units(
    reference: &Allele,
    alternates: &[Allele],
    ref_bases_starting_at_vc_without_pad: &[u8],
) -> Option<(Vec<i32>, Vec<u8>)> {
    if alternates
        .iter()
        .all(|allele| allele.len() == reference.len())
    {
        return None;
    }
    let ref_allele_bases = strip_padding_base(reference)?;

    let mut repeat_unit: Vec<u8> = Vec::new();
    let mut lengths: Vec<i32> = Vec::new();
    for allele in alternates {
        let alt_bases = strip_padding_base(allele)?;
        let (counts, unit) = num_tandem_repeat_units_for_bases(
            &ref_allele_bases,
            &alt_bases,
            ref_bases_starting_at_vc_without_pad,
        )?;
        // Zero on either side drops the annotation for the whole site.
        if counts[0] == 0 || counts[1] == 0 {
            return None;
        }
        if lengths.is_empty() {
            lengths.push(counts[0]);
        }
        lengths.push(counts[1]);
        // The last alternate's unit wins, and nothing checks that the alternates agreed on it.
        repeat_unit = unit;
    }
    Some((lengths, repeat_unit))
}

/// `Arrays.copyOfRange(allele.getBases(), 1, allele.length())`, which throws when the length is
/// zero, that is, when the allele is symbolic or a no-call.
fn strip_padding_base(allele: &Allele) -> Option<Vec<u8>> {
    let bases = allele.display_string().into_bytes();
    if allele.is_symbolic() || bases.is_empty() {
        return None;
    }
    Some(bases[1..].to_vec())
}
