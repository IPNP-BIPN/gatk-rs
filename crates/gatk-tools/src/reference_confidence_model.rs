//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel`
//! (`calculateRefConfidence` and everything under it) and `AssemblyBasedCallerUtils.getPileupsOverReference`
//! (GATK 4.6.2.0): one reference-confidence record per base of a region, as GVCF mode writes it.
//!
//! At each base of the region's span the reads are piled up again (re-sorted by coordinate, since
//! realignment moved them) and the base gets either the call that starts there or a hom-ref
//! record over `<NON_REF>`. That record's likelihoods are the less confident of two: the
//! ref-vs-any likelihoods of the pileup (bases of quality 6 or below ignored, a deletion counted
//! at the reference-model deletion quality, a base alternate when it differs from the reference
//! or is a deletion), capped by the hom-ref likelihood; and the indel likelihoods of the number
//! of reads with no plausible indel of up to ten bases there, capped at forty reads.
//!
//! What a reader would not guess: whether a read has a plausible indel is computed **once per
//! read and call**, at the first base it is asked about, and the answer for every later base is
//! read from that bitset; and the read's bases are laid one to one against the reference by an
//! array that counts soft clips but never fills them, so a read with both a soft clip and an
//! indel ends in zero bytes.

use std::collections::HashMap;

use gatk_engine::genotype_likelihood_calculator::gls_to_pls;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::locus_iterator::{self, LocusIteratorOptions};
use gatk_engine::pair_hmm::approximate_log10_sum_log10;
use gatk_engine::pileup::PileupElement;
use gatk_engine::qual_quantizer::qual_to_prob_log10;
use gatk_engine::read_states::{DownsamplingInfo, ReadStateManager};
use gatk_engine::read_utils;
use htsjdk_bam::cigar::Op;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::variant::{Genotype, VariantContext};

use crate::calculate_genotype_posteriors::gq_log10_from_likelihoods;
use crate::genotyping_engine::EngineError;
use crate::gvcf_blocks::calculate_gq_from_pls;

/// `ReferenceConfidenceModel.BASE_QUAL_THRESHOLD`.
const BASE_QUAL_THRESHOLD: u8 = 6;
/// `MAX_N_INDEL_INFORMATIVE_READS`.
const MAX_N_INDEL_INFORMATIVE_READS: usize = 40;
/// `INDEL_QUAL`, `Math.round(-4.5 * -10.0)`.
const INDEL_QUAL: u8 = 45;
/// `AlignmentUtils.GAP_CHARACTER`.
const GAP_CHARACTER: u8 = b'-';

/// The model's switches.
#[derive(Debug, Clone)]
pub struct RefConfidenceArguments {
    /// `--indel-size-to-eliminate-in-ref-model`.
    pub indel_informative_depth_indel_size: usize,
    /// `--reference-model-deletion-quality`.
    pub ref_model_deletion_quality: u8,
    /// `!--dont-use-soft-clipped-bases`.
    pub use_soft_clipped_bases: bool,
}

impl Default for RefConfidenceArguments {
    fn default() -> Self {
        RefConfidenceArguments {
            indel_informative_depth_indel_size: 10,
            ref_model_deletion_quality: 30,
            use_soft_clipped_bases: true,
        }
    }
}

/// `calculateRefConfidence(refHaplotype, ..., activeRegion, readLikelihoods, ploidyModel,
/// variantCalls)` for the one sample: `ref_bases` are the reference haplotype's, which span the
/// region's padded span starting at `padded_start`; `evidence` is the likelihoods' evidence.
#[allow(clippy::too_many_arguments)]
pub fn calculate_ref_confidence(
    ref_bases: &[u8],
    span: &SimpleInterval,
    padded_start: i32,
    evidence: &[BamRecord],
    sample: &str,
    ploidy: usize,
    variant_calls: &[VariantContext],
    header: &SamHeader,
    arguments: &RefConfidenceArguments,
) -> Result<Vec<VariantContext>, EngineError> {
    let mut reads: Vec<BamRecord> = evidence.to_vec();
    reads.sort_by(read_utils::compare_read_coordinate);
    let partition = vec![Some(sample.to_string())];
    let states = ReadStateManager::new(partition.clone(), DownsamplingInfo::NONE)
        .map_err(|e| runtime(format!("{e:?}")))?;
    let contexts = locus_iterator::contexts(
        &reads,
        partition,
        header,
        LocusIteratorOptions {
            include_deletions: true,
            include_ns: false,
        },
        states,
    )
    .map_err(|e| runtime(format!("{e:?}")))?;
    let by_position: HashMap<i32, usize> = contexts
        .iter()
        .enumerate()
        .filter(|(_, context)| context.contig == span.contig)
        .map(|(index, context)| (context.position, index))
        .collect();

    let global_ref_offset = span.start - padded_start;
    let mut cache: HashMap<usize, Vec<bool>> = HashMap::new();
    let mut results = Vec::with_capacity(span.end as usize - span.start as usize + 1);
    let empty: Vec<PileupElement<'_>> = Vec::new();
    for position in span.start..=span.end {
        let offset = position - span.start;
        let overlapping = overlapping_variant_context(&span.contig, position, variant_calls);
        if let Some(site) = overlapping.filter(|site| site.start == i64::from(position)) {
            results.push(site.clone());
            continue;
        }
        let pileup = by_position
            .get(&position)
            .map(|&index| &contexts[index].pileup.elements)
            .unwrap_or(&empty);
        results.push(reference_confidence_context(
            ploidy,
            ref_bases,
            sample,
            global_ref_offset,
            pileup,
            &span.contig,
            position,
            offset,
            arguments,
            &mut cache,
        )?);
    }
    Ok(results)
}

fn runtime(message: String) -> EngineError {
    EngineError::Runtime {
        class: "IllegalStateException".to_string(),
        message,
    }
}

/// `GATKVariantContextUtils.getOverlappingVariantContext`: of the calls overlapping the base, the
/// one starting last.
fn overlapping_variant_context<'v>(
    contig: &str,
    position: i32,
    calls: &'v [VariantContext],
) -> Option<&'v VariantContext> {
    let mut overlaps: Option<&VariantContext> = None;
    for vc in calls {
        let position = i64::from(position);
        if vc.contig == contig
            && vc.start <= position
            && position <= vc.stop
            && overlaps.is_none_or(|o| vc.start > o.start)
        {
            overlaps = Some(vc);
        }
    }
    overlaps
}

/// `makeReferenceConfidenceVariantContext`, without priors.
#[allow(clippy::too_many_arguments)]
fn reference_confidence_context(
    ploidy: usize,
    ref_bases: &[u8],
    sample: &str,
    global_ref_offset: i32,
    pileup: &[PileupElement<'_>],
    contig: &str,
    position: i32,
    offset: i32,
    arguments: &RefConfidenceArguments,
    cache: &mut HashMap<usize, Vec<bool>>,
) -> Result<VariantContext, EngineError> {
    let ref_offset = (offset + global_ref_offset) as usize;
    let ref_base = ref_bases[ref_offset];
    let result = ref_vs_any(ploidy, pileup, position, ref_base, arguments)?;

    let ref_allele = Allele::create(&[ref_base], true).map_err(|e| runtime(e.to_string()))?;
    let non_ref = Allele::from_str("<NON_REF>", false).expect("a symbolic allele");
    let mut vc = VariantContext::new(
        contig,
        i64::from(position),
        vec![ref_allele.clone(), non_ref],
    );
    vc.stop = i64::from(position);
    let mut genotype = Genotype::new(sample, vec![ref_allele; ploidy]);
    genotype.ad = Some(vec![result.ref_depth, result.non_ref_depth]);
    genotype.dp = Some(result.ref_depth + result.non_ref_depth);

    // `doIndelRefConfCalc`: the SNP likelihoods capped by the hom-ref one, against the indel
    // likelihoods of the informative reads, and the less confident of the two.
    let capped: Vec<f64> = result
        .likelihoods
        .iter()
        .map(|&l| l.min(result.likelihoods[0]))
        .collect();
    let informative = reads_with_no_plausible_indels(
        pileup,
        ref_offset,
        ref_bases,
        arguments.indel_informative_depth_indel_size,
        cache,
    );
    let indel = indel_likelihoods(ploidy, informative.min(MAX_N_INDEL_INFORMATIVE_READS));
    let worst = if gq_log10_from_likelihoods(0, &indel) > gq_log10_from_likelihoods(0, &capped) {
        indel
    } else {
        capped
    };
    let pls = gls_to_pls(&worst);
    genotype.gq = Some(calculate_gq_from_pls(&pls)?);
    genotype.pl = Some(pls);
    vc.genotypes = GenotypesContext::new(vec![genotype]);
    Ok(vc)
}

/// `RefVsAnyResult`: the likelihoods and the two depths.
struct RefVsAny {
    likelihoods: Vec<f64>,
    ref_depth: i32,
    non_ref_depth: i32,
}

/// `calcGenotypeLikelihoodsOfRefVsAny(ploidy, pileup, refBase, BASE_QUAL_THRESHOLD, null, true)`.
fn ref_vs_any(
    ploidy: usize,
    pileup: &[PileupElement<'_>],
    position: i32,
    ref_base: u8,
    arguments: &RefConfidenceArguments,
) -> Result<RefVsAny, EngineError> {
    let count = ploidy + 1;
    let log10_ploidy = jmath::math::log10(ploidy as f64);
    let log10_one_third = -jmath::math::log10(3.0);
    let mut result = RefVsAny {
        likelihoods: vec![0.0; count],
        ref_depth: 0,
        non_ref_depth: 0,
    };
    let mut read_count = 0usize;
    for element in pileup {
        let qual = if element.is_deletion() {
            arguments.ref_model_deletion_quality
        } else {
            element.qual()
        };
        if (qual as i8) <= (BASE_QUAL_THRESHOLD as i8) && !element.is_deletion() {
            continue;
        }
        if !arguments.use_soft_clipped_bases {
            let start = original_soft_clip(element.read, b"os")?;
            let end = original_soft_clip(element.read, b"oe")?;
            if start > position || end < position {
                continue;
            }
        }
        read_count += 1;
        // `isAltAfterAssembly`.
        let is_alt = element.base() != ref_base || element.is_deletion();
        let qual_to_error_log10 = f64::from(qual) * -0.1;
        let (reference, non_ref) = if is_alt {
            result.non_ref_depth += 1;
            (
                qual_to_error_log10 + log10_one_third,
                qual_to_prob_log10(qual),
            )
        } else {
            result.ref_depth += 1;
            (
                qual_to_prob_log10(qual),
                qual_to_error_log10 + log10_one_third,
            )
        };
        result.likelihoods[0] += reference + log10_ploidy;
        result.likelihoods[count - 1] += non_ref + log10_ploidy;
        let mut j = count as i64 - 2;
        for (i, slot) in result
            .likelihoods
            .iter_mut()
            .enumerate()
            .take(count - 1)
            .skip(1)
        {
            *slot += approximate_log10_sum_log10(
                reference + jmath::math::log10(j as f64),
                non_ref + jmath::math::log10(i as f64),
            );
            j -= 1;
        }
    }
    let denominator = read_count as f64 * log10_ploidy;
    for value in &mut result.likelihoods {
        *value -= denominator;
    }
    Ok(result)
}

/// `getOriginalSoftStart` and `getOriginalSoftEnd`: the `os` and `oe` tags finalization wrote.
fn original_soft_clip(read: &BamRecord, tag: &[u8; 2]) -> Result<i32, EngineError> {
    match read.tags.get(Tag::new(tag)) {
        Some(TagValue::Int(value)) => Ok(*value as i32),
        _ => Err(EngineError::Runtime {
            class: "GATKException".to_string(),
            message: format!(
                "Attempt to read soft clip {} that was not saved",
                if tag == b"os" { "start" } else { "end" }
            ),
        }),
    }
}

/// `getIndelPLs(ploidy, n)`, as log10 likelihoods: none informative is all zeros, otherwise each
/// read contributes the no-indel likelihood to hom-ref and a mixture to every other genotype.
fn indel_likelihoods(ploidy: usize, informative: usize) -> Vec<f64> {
    let mut likelihoods = vec![0.0; ploidy + 1];
    if informative == 0 {
        return likelihoods;
    }
    let no_indel = qual_to_prob_log10(INDEL_QUAL);
    let indel = f64::from(INDEL_QUAL) / -10.0;
    let denominator = -jmath::math::log10(ploidy as f64);
    let n = informative as f64;
    likelihoods[0] = n * no_indel;
    for (alt_count, slot) in likelihoods.iter_mut().enumerate().skip(1) {
        let reference = no_indel + jmath::math::log10((ploidy - alt_count) as f64);
        let alternate = indel + jmath::math::log10(alt_count as f64);
        *slot = n * (approximate_log10_sum_log10(reference, alternate) + denominator);
    }
    likelihoods
}

/// `calcNReadsWithNoPlausibleIndelsReads`: the reads at the base with no plausible indel of up to
/// `max_indel_size` bases, stopping at one past the cap.
fn reads_with_no_plausible_indels(
    pileup: &[PileupElement<'_>],
    ref_offset: usize,
    ref_bases: &[u8],
    max_indel_size: usize,
    cache: &mut HashMap<usize, Vec<bool>>,
) -> usize {
    let mut informative = 0usize;
    for element in pileup {
        if element.is_before_deletion_start()
            || element.is_before_insertion()
            || element.is_deletion()
        {
            continue;
        }
        let offset = cigar_modified_offset(element);
        let key = element.read as *const BamRecord as usize;
        let bits = cache.entry(key).or_insert_with(|| {
            informative_bases(element.read, offset, ref_bases, ref_offset, max_indel_size)
        });
        if bits.get(offset).copied().unwrap_or(false) {
            informative += 1;
            if informative > MAX_N_INDEL_INFORMATIVE_READS {
                return MAX_N_INDEL_INFORMATIVE_READS;
            }
        }
    }
    informative
}

/// `getCigarModifiedOffset`: the offset counting the bases of every element before the current
/// one that consumes the reference or is a soft clip, insertions left out.
fn cigar_modified_offset(element: &PileupElement<'_>) -> usize {
    let counts = |op: Op| op.consumes_reference_bases() || op == Op::S;
    let mut offset = if counts(element.current_cigar_element.op) {
        element.offset_in_current_cigar as usize
    } else {
        0
    };
    for e in element
        .read
        .cigar
        .elements
        .iter()
        .take(element.current_cigar_offset as usize)
    {
        if counts(e.op) {
            offset += e.length as usize;
        }
    }
    offset
}

/// `Nucleotide`'s one-bit-per-base mask: the four bases, the IUPAC codes as their unions, and
/// zero for anything else.
fn nucleotide_mask(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => 0,
    }
}

/// `isMismatchAndNotAnAlignmentGap`.
fn is_mismatch(read_base: u8, ref_base: u8) -> bool {
    nucleotide_mask(read_base) & nucleotide_mask(ref_base) == 0 && read_base != GAP_CHARACTER
}

/// `AlignmentUtils.getBasesAndBaseQualitiesAlignedOneToOne`: the read untouched when it has no
/// indel; otherwise an array as long as its reference bases and soft clips, deletions as gaps of
/// quality zero, insertions and soft clips skipped.
fn aligned_one_to_one(read: &BamRecord) -> (Vec<u8>, Vec<u8>) {
    let elements = &read.cigar.elements;
    if !elements.iter().any(|e| matches!(e.op, Op::I | Op::D)) {
        return (read.read_bases.clone(), read.base_qualities.clone());
    }
    let length: usize = elements
        .iter()
        .filter(|e| e.op.consumes_reference_bases() || e.op == Op::S)
        .map(|e| e.length as usize)
        .sum();
    let mut bases = vec![0u8; length];
    let mut quals = vec![0u8; length];
    let mut literal = 0usize;
    let mut padded = 0usize;
    for e in elements {
        let n = e.length as usize;
        if e.op.consumes_read_bases() {
            if e.op.consumes_reference_bases() {
                bases[padded..padded + n].copy_from_slice(&read.read_bases[literal..literal + n]);
                quals[padded..padded + n]
                    .copy_from_slice(&read.base_qualities[literal..literal + n]);
                padded += n;
            }
            literal += n;
        } else if e.op.consumes_reference_bases() {
            for _ in 0..n {
                bases[padded] = GAP_CHARACTER;
                quals[padded] = 0;
                padded += 1;
            }
        }
    }
    (bases, quals)
}

/// `readHasNoPlausibleIdealsOfSize`'s bitset, computed at `read_start` against the reference from
/// `ref_start`: true where no indel of up to `max_indel_size` bases explains the read better than
/// its alignment.
fn informative_bases(
    read: &BamRecord,
    read_start: usize,
    ref_bases: &[u8],
    ref_start: usize,
    max_indel_size: usize,
) -> Vec<bool> {
    let read_length = read.read_bases.len();
    let mut bits = vec![false; read_length.max(1)];
    let mut set = |bits: &mut Vec<bool>, index: usize, value: bool| {
        if index >= bits.len() {
            bits.resize(index + 1, false);
        }
        bits[index] = value;
    };
    if read_length as i64 - (read_start as i64) < max_indel_size as i64
        || ref_bases.len() as i64 - (ref_start as i64) < max_indel_size as i64
    {
        return bits;
    }
    let secondary_break = read_length as i64 - max_indel_size as i64;
    let (bases, quals) = aligned_one_to_one(read);
    if bases.len() as i64 - read_start as i64 <= max_indel_size as i64 {
        return bits;
    }
    let (last, reference_shorter) = if bases.len() < ref_bases.len() - ref_start + read_start + 1 {
        (bases.len() as i64 - max_indel_size as i64, false)
    } else {
        (
            ref_bases.len() as i64 - ref_start as i64 + read_start as i64 - max_indel_size as i64
                + 1,
            true,
        )
    };
    let baseline = baseline_mismatch_qualities(&bases, &quals, read_start, ref_bases, ref_start);
    for indel_size in 1..=max_indel_size {
        for insertion in [false, true] {
            traverse_for_indel_mismatches(
                &mut bits,
                &mut set,
                read_start,
                &bases,
                &quals,
                last,
                secondary_break,
                ref_start,
                ref_bases,
                &baseline,
                indel_size,
                insertion,
            );
        }
    }
    let flip_to = if last <= secondary_break {
        last
    } else {
        secondary_break + 1
    };
    for index in 0..flip_to.max(0) as usize {
        let current = bits.get(index).copied().unwrap_or(false);
        set(&mut bits, index, !current);
    }
    if last <= secondary_break && reference_shorter {
        set(&mut bits, (last - 1) as usize, false);
    }
    bits
}

/// `calculateBaselineMMQualities`: from each position on, the summed quality of the mismatches
/// against the reference.
fn baseline_mismatch_qualities(
    bases: &[u8],
    quals: &[u8],
    read_start: usize,
    ref_bases: &[u8],
    ref_start: usize,
) -> Vec<i32> {
    let n = (bases.len() - read_start).min(ref_bases.len() - ref_start);
    let mut results = vec![0; n];
    let mut sum = 0i32;
    for i in (0..n).rev() {
        if is_mismatch(bases[read_start + i], ref_bases[ref_start + i]) {
            sum += i32::from(quals[read_start + i] as i8);
        }
        results[i] = sum;
    }
    results
}

/// `traverseEndOfReadForIndelMismatches`: walking back from the end of the read with an indel of
/// `indel_size` assumed, mark each base where the indel costs no more than the alignment.
#[allow(clippy::too_many_arguments)]
fn traverse_for_indel_mismatches(
    bits: &mut Vec<bool>,
    set: &mut dyn FnMut(&mut Vec<bool>, usize, bool),
    read_start: usize,
    bases: &[u8],
    quals: &[u8],
    last: i64,
    secondary_break: i64,
    ref_start: usize,
    ref_bases: &[u8],
    baseline: &[i32],
    indel_size: usize,
    insertion: bool,
) {
    let global = baseline[0];
    let mut quality_sum = 0i32;
    let insertion_length = if insertion { indel_size as i64 } else { 0 };
    let deletion_length = if insertion { 0 } else { indel_size as i64 };
    let direct = (bases.len() as i64 - read_start as i64 - insertion_length)
        .min(ref_bases.len() as i64 - ref_start as i64 - deletion_length);
    let mut read_offset = direct + insertion_length - 1;
    let mut ref_offset = direct + deletion_length - 1;
    while read_offset >= 0 && ref_offset >= 0 {
        let read_base = bases[read_start + read_offset as usize];
        let ref_base = ref_bases[ref_start + ref_offset as usize];
        if is_mismatch(read_base, ref_base) {
            quality_sum += i32::from(quals[read_start + read_offset as usize] as i8);
            if quality_sum > global {
                break;
            }
        }
        let site = read_offset.min(ref_offset);
        let at = read_start as i64 + site;
        if bases[at as usize] != GAP_CHARACTER
            && at < last
            && at <= secondary_break
            && baseline[site as usize] >= quality_sum
        {
            set(bits, at as usize, true);
        }
        read_offset -= 1;
        ref_offset -= 1;
    }
}
