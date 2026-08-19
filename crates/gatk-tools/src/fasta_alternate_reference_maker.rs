//! Ported from `org.broadinstitute.hellbender.tools.walkers.fasta.FastaAlternateReferenceMaker`
//! (GATK 4.6.2.0).
//!
//! [`crate::fasta_reference_maker`] with a VCF applied: every locus of the traversal contributes
//! zero, one or several bases depending on the record at it and on state carried from the loci
//! before it.
//!
//! # The deletion counter is the only state that crosses an apply
//!
//! `deletionBasesRemaining` is set when a simple deletion is seen and decremented at every
//! following locus, which is what makes the bases removed the ones AFTER the record's position
//! rather than at it. A traversal that started inside a deletion never sees the record, so nothing
//! is dropped: the counter is not recovered from the interval.
//!
//! # The mask and the call meet at the same locus
//!
//! With `--snp-mask-priority` the mask is consulted first and wins; without it the called variants
//! are consulted first and a mask at the same site is never reached. So the flag is not "prefer the
//! mask", it is "look at the mask first", and the difference shows only where both exist.
//!
//! # Where the reference crashes
//!
//! `--use-iupac-sample` on a sample homozygous for a spanning deletion returns `EMPTY_BASE`, the
//! string `" "`, which the FASTA writer then refuses. The port reproduces the refusal rather than
//! inventing a base: `fasta-alternate-reference-maker`'s golden carries the exception, so this is
//! measured behaviour and not a guess about intent.

use htsjdk_bam::fasta_writer::{FastaOutputs, FastaReferenceWriter};
use htsjdk_vcf::variant::VariantContext;

use gatk_engine::interval::SimpleInterval;
use gatk_engine::interval_args::IntervalArguments;
use gatk_engine::reference::ReferenceFileSource;

use crate::fasta_reference_maker::MakerError;
use crate::reference_walker::{self};

/// `FastaAlternateReferenceMaker.EMPTY_BASE`, which is a space and not an empty string.
pub const EMPTY_BASE: u8 = b' ';

/// What the argument checks refuse, before the traversal runs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ArgumentError {
    /// `--snp-mask-priority` with no `--snp-mask`.
    PriorityWithoutMask,
    /// `--use-iupac-sample` naming a sample the VCF header does not carry.
    UnknownIupacSample,
}

impl ArgumentError {
    /// The exception class the reference throws.
    pub fn java_class(&self) -> &'static str {
        match self {
            ArgumentError::PriorityWithoutMask => {
                "org.broadinstitute.barclay.argparser.CommandLineException"
            }
            ArgumentError::UnknownIupacSample => {
                "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
            }
        }
    }

    /// The message, including the double space the reference's concatenation leaves in the first.
    pub fn message(&self) -> String {
        match self {
            ArgumentError::PriorityWithoutMask => {
                "Cannot specify --snp-mask-priority without  --snp-mask".to_string()
            }
            ArgumentError::UnknownIupacSample => {
                "Bad input: the IUPAC sample specified is not present in the provided VCF file"
                    .to_string()
            }
        }
    }
}

/// What stopped a run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlternateError {
    /// An argument check, which fires before anything is read.
    Argument(ArgumentError),
    /// The traversal or the writer, which is where the space lands.
    Maker(MakerError),
}

/// `BaseUtils.basesToIUPAC`.
///
/// The two bases are ordered first, so the table is only half written; anything that is not a plain
/// `ACGT` gives `N`, and two equal bases give themselves rather than a code.
pub fn bases_to_iupac(first: u8, second: u8) -> u8 {
    if second < first {
        return bases_to_iupac(second, first);
    }
    if !is_regular_base(first) || !is_regular_base(second) {
        return b'N';
    }
    // `basesAreEqual`, which is case-insensitive and answers the FIRST base rather than a code.
    if first.eq_ignore_ascii_case(&second) {
        return first;
    }
    match first.to_ascii_uppercase() {
        b'A' => match second.to_ascii_uppercase() {
            b'C' => b'M',
            b'G' => b'R',
            _ => b'W',
        },
        b'C' => match second.to_ascii_uppercase() {
            b'G' => b'S',
            _ => b'Y',
        },
        // The only pair left is G with T.
        _ => b'K',
    }
}

/// `BaseUtils.isRegularBase`: one of `ACGT`, in either case.
fn is_regular_base(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// `VariantContext.isSimpleIndel`.
fn is_simple_indel(variant: &VariantContext) -> bool {
    if variant.alleles.len() != 2 {
        return false;
    }
    let reference = variant.alleles[0].base_string();
    let alternate = variant.alleles[1].base_string();
    if reference.len() == alternate.len() {
        return false;
    }
    !reference.is_empty()
        && !alternate.is_empty()
        && reference.as_bytes()[0] == alternate.as_bytes()[0]
        && (reference.len() == 1 || alternate.len() == 1)
}

/// `VariantContext.isSimpleDeletion`.
fn is_simple_deletion(variant: &VariantContext) -> bool {
    is_simple_indel(variant) && variant.alleles[1].base_string().len() == 1
}

/// `VariantContext.isSimpleInsertion`.
fn is_simple_insertion(variant: &VariantContext) -> bool {
    is_simple_indel(variant) && variant.alleles[0].base_string().len() == 1
}

/// `VariantContext.isSNP`, which is the type being SNP: every alternate one base against a
/// one-base reference.
fn is_snp(variant: &VariantContext) -> bool {
    crate::remove_nearby_indels::variant_type(variant)
        == crate::remove_nearby_indels::VariantType::Snp
}

/// `getFirstConcreteAltAllele`: the first alternate that is neither symbolic nor the spanning
/// deletion.
fn first_concrete_alternate(variant: &VariantContext) -> Option<String> {
    variant
        .alleles
        .iter()
        .skip(1)
        .find(|allele| !allele.is_symbolic() && allele.base_string() != "*")
        .map(|allele| allele.base_string())
}

/// `getIUPACBase(genotype)`.
///
/// A spanning deletion in the genotype is handled before anything else, and a HOM VAR spanning
/// deletion answers `EMPTY_BASE` -- a space, which the writer refuses. That is the reference's
/// behaviour and the golden carries the crash.
fn iupac_base(alleles: &[String]) -> Vec<u8> {
    let spanning = alleles.iter().any(|allele| allele == "*");
    if spanning {
        let hom_var = alleles.iter().all(|allele| allele == "*");
        if hom_var {
            return vec![EMPTY_BASE];
        }
        let kept = if alleles[0] == "*" {
            alleles[1].clone()
        } else {
            alleles[0].clone()
        };
        return kept.into_bytes();
    }
    // `isHet` is false for a hom ref and a hom var alike, and both take the first allele whole.
    let het = alleles.len() == 2 && alleles[0] != alleles[1];
    if !het {
        return alleles[0].clone().into_bytes();
    }
    vec![bases_to_iupac(
        alleles[0].as_bytes()[0],
        alleles[1].as_bytes()[0],
    )]
}

/// The genotype of one sample, as this tool reads it.
fn genotype_alleles(variant: &VariantContext, sample: &str) -> Option<Vec<String>> {
    variant
        .genotypes
        .iter()
        .find(|genotype| genotype.sample_name == sample)
        .map(|genotype| {
            genotype
                .alleles
                .iter()
                .map(|allele| allele.base_string())
                .collect()
        })
}

/// The arguments this tool adds to the maker's.
#[derive(Debug, Clone, Default)]
pub struct AlternateArguments<'a> {
    /// `--snp-mask`, whose records write `N`.
    pub mask: Option<&'a [VariantContext]>,
    /// `--snp-mask-priority`.
    pub mask_priority: bool,
    /// `--use-iupac-sample`.
    pub iupac_sample: Option<String>,
}

/// `onTraversalStart`'s two checks, which run before the reference is opened.
///
/// `samples` is the VCF header's sample list, which is where the second check looks: a sample that
/// appears in no record but is declared in the header passes.
pub fn check_arguments(
    arguments: &AlternateArguments,
    samples: &[String],
) -> Result<(), ArgumentError> {
    if arguments.mask_priority && arguments.mask.is_none() {
        return Err(ArgumentError::PriorityWithoutMask);
    }
    if let Some(sample) = &arguments.iupac_sample {
        if !samples.iter().any(|name| name == sample) {
            return Err(ArgumentError::UnknownIupacSample);
        }
    }
    Ok(())
}

/// `doWork`: the FASTA, its index and its dictionary.
pub fn run(
    reference: &mut ReferenceFileSource,
    intervals: &IntervalArguments,
    bases_per_line: usize,
    variants: &[VariantContext],
    arguments: &AlternateArguments,
    samples: &[String],
) -> Result<FastaOutputs, AlternateError> {
    check_arguments(arguments, samples).map_err(AlternateError::Argument)?;

    let mut writer = FastaReferenceWriter::new(bases_per_line, true)
        .map_err(|error| AlternateError::Maker(MakerError::Writer(error)))?;

    let applied =
        reference_walker::traverse(reference, intervals, |locus: &SimpleInterval| locus.clone())
            .map_err(|error| AlternateError::Maker(MakerError::Traversal(error)))?;

    let mut deletion_bases_remaining = 0i32;
    let mut count = 0usize;
    let mut last: Option<SimpleInterval> = None;
    let mut start_position = 0;
    let mut sequence: Vec<u8> = Vec::new();

    for call in &applied {
        let interval = &call.window;
        let bases = handle_position(
            interval,
            call.bases[0],
            variants,
            arguments,
            &mut deletion_bases_remaining,
        );

        // `apply`: an empty answer advances the position without adding a base, which keeps the
        // sequence unbroken across a deletion rather than starting a new one after it.
        let new_sequence = match &last {
            None => true,
            Some(previous) => !within_distance_of(previous, interval),
        };
        if new_sequence {
            if last.is_some() {
                finalize(
                    &mut writer,
                    count,
                    &last,
                    start_position,
                    &sequence,
                    bases_per_line,
                )?;
            }
            count += 1;
            start_position = interval.start;
            sequence.clear();
        }
        last = Some(interval.clone());
        sequence.extend_from_slice(&bases);
    }

    if last.is_some() {
        finalize(
            &mut writer,
            count,
            &last,
            start_position,
            &sequence,
            bases_per_line,
        )?;
    }

    writer
        .close()
        .map_err(|error| AlternateError::Maker(MakerError::Writer(error)))
}

/// `handlePosition`: what this locus contributes.
fn handle_position(
    interval: &SimpleInterval,
    base: u8,
    variants: &[VariantContext],
    arguments: &AlternateArguments,
    deletion_bases_remaining: &mut i32,
) -> Vec<u8> {
    if *deletion_bases_remaining > 0 {
        *deletion_bases_remaining -= 1;
        return Vec::new();
    }

    if arguments.mask_priority && is_masked(interval, arguments.mask) {
        return vec![b'N'];
    }

    for variant in variants {
        if variant.contig != interval.contig || variant.start as i32 != interval.start {
            continue;
        }
        if variant.is_filtered() {
            continue;
        }
        if is_simple_deletion(variant) {
            // The next n bases go, not this one.
            *deletion_bases_remaining = variant.alleles[0].base_string().len() as i32 - 1;
            return vec![base];
        }
        if is_simple_insertion(variant) || is_snp(variant) {
            let allele = first_concrete_alternate(variant).unwrap_or_else(|| {
                String::from_utf8(vec![EMPTY_BASE]).expect("a space is a string")
            });
            if is_simple_insertion(variant) {
                return allele.into_bytes();
            }
            return match &arguments.iupac_sample {
                Some(sample) => match genotype_alleles(variant, sample) {
                    Some(alleles) => iupac_base(&alleles),
                    // `getGenotype` returning null is a `Utils.nonNull` failure upstream; nothing
                    // in a decoded VCF reaches it, because every record carries every sample.
                    None => allele.into_bytes(),
                },
                None => allele.into_bytes(),
            };
        }
    }

    if !arguments.mask_priority && is_masked(interval, arguments.mask) {
        return vec![b'N'];
    }

    vec![base]
}

/// `isMasked`: any mask record starting at this locus, filtered or not.
fn is_masked(interval: &SimpleInterval, mask: Option<&[VariantContext]>) -> bool {
    let Some(mask) = mask else {
        return false;
    };
    mask.iter()
        .any(|variant| variant.contig == interval.contig && variant.start as i32 == interval.start)
}

/// `finalizeSequence`, which is the maker's.
fn finalize(
    writer: &mut FastaReferenceWriter,
    count: usize,
    last: &Option<SimpleInterval>,
    start_position: i32,
    sequence: &[u8],
    bases_per_line: usize,
) -> Result<(), AlternateError> {
    let last = last.as_ref().expect("a sequence is open");
    let description = format!("{}:{}-{}", last.contig, start_position, last.end);
    writer
        .start_sequence_with(&count.to_string(), &description, bases_per_line)
        .map_err(|error| AlternateError::Maker(MakerError::Writer(error)))?;
    writer
        .append_bases(sequence)
        .map_err(|error| AlternateError::Maker(MakerError::Writer(error)))
}

/// `Locatable.withinDistanceOf(interval, 1)`.
fn within_distance_of(left: &SimpleInterval, right: &SimpleInterval) -> bool {
    left.contig == right.contig && left.start <= right.end + 1 && right.start <= left.end + 1
}
