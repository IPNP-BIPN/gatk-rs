//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyRegionTrimmer`
//! (GATK 4.7.0.0; `AssemblyRegionTrimmer.java` is unchanged since 4.6.2.0): the part of an
//! assembled region that is genotyped.
//!
//! The default trim takes the events overlapping the region's active span, spans them, and pads
//! each one: a SNP by `snpPaddingForGenotyping`, an indel by `indelPaddingForGenotyping`, or, when
//! the indel is an expansion or contraction of a tandem repeat in the reference that follows it, by
//! `strPaddingForGenotyping` plus the longest run of the repeat either allele carries. The padded
//! span is cut to the region's padded span; the variant span is cut to the active span.
//!
//! The legacy trim (`enableLegacyAssemblyRegionTrimming`) pads the variant span once, by the indel
//! padding if any overlapping event is not a SNP, and caps it at the active span extended by
//! `maxExtensionIntoRegionPadding`.
//!
//! The reference context the trimmer reads is, in every caller, the region's padded span: it is
//! taken here as those bases and the position of the first.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<SimpleInterval>` is "an interval or nothing", the reference's nullable
//!   `SimpleInterval` field, checked by the compiler;
//! - `Result<T, TrimmerError>` is "a value or an error", how the reference's exceptions are
//!   modelled; a trailing `?` returns the error to the caller at once, `.ok_or(e)` and
//!   `.ok_or_else(|| e)` turn a `None` into that error, and `.map_err(f)` wraps another module's
//!   error in this one's;
//! - `let (Some(a), Some(b)) = (x, y) else { ... };` unwraps both or takes the early exit;
//! - `.iter().filter(..).collect()` is Java's stream `filter` then `collect(toList())`, and
//!   `|e| ...` is a lambda;
//! - `x as i32` and `x as usize` convert between integer widths, as a Java cast does.
//!
//! Long Java lines quoted in this file are broken at token boundaries to fit its width, and a
//! trailing comment on a quoted line is left out.

use htsjdk_bam::header::SamHeader;

use crate::assembly_region::{AssemblyRegion, RegionError};
use crate::event_map::Event;
use crate::interval::SimpleInterval;
use crate::tandem_repeat_units::num_tandem_repeat_units;

/// The `AssemblyRegionArgumentCollection` fields the trimmer reads (the class is
/// `org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection`, the one the
/// trimmer imports).
#[derive(Debug, Clone)]
pub struct TrimmerArguments {
    /// `assemblyRegionPadding`: the padding of a flank region (see
    /// [`AssemblyRegionTrimmer::non_variant_left_flank_region`]).
    pub assembly_region_padding: i32,
    /// `indelPaddingForGenotyping`.
    pub indel_padding_for_genotyping: i32,
    /// `snpPaddingForGenotyping`.
    pub snp_padding_for_genotyping: i32,
    /// `strPaddingForGenotyping`, added to the repeat length for an indel in a tandem repeat.
    pub str_padding_for_genotyping: i32,
    /// `maxExtensionIntoRegionPadding`, read only by the legacy trim.
    pub max_extension_into_region_padding: i32,
    /// `enableLegacyAssemblyRegionTrimming`.
    pub enable_legacy_assembly_region_trimming: bool,
}

impl Default for TrimmerArguments {
    /// The argument collection's defaults: `DEFAULT_ASSEMBLY_REGION_PADDING` (100) and the field
    /// initialisers `indelPaddingForGenotyping = 75`, `snpPaddingForGenotyping = 20`,
    /// `strPaddingForGenotyping = 75`, `maxExtensionIntoRegionPadding = 25` and
    /// `enableLegacyAssemblyRegionTrimming = false`.
    fn default() -> Self {
        TrimmerArguments {
            assembly_region_padding: 100,
            indel_padding_for_genotyping: 75,
            snp_padding_for_genotyping: 20,
            str_padding_for_genotyping: 75,
            max_extension_into_region_padding: 25,
            enable_legacy_assembly_region_trimming: false,
        }
    }
}

/// What the trimmer refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TrimmerError {
    /// `CommandLineException.BadArgumentValue` from `validate()`.
    BadArgument(String),
    /// `IllegalArgumentException`: intervals that do not overlap, a padded span that does not
    /// contain the variant span, or a reference context the event starts past (the last is
    /// `Arrays.copyOfRange` with a start past the end).
    IllegalArgument(String),
    /// `SimpleInterval.mergeWithContiguous`'s
    /// `throw new GATKException("The two intervals need to be contiguous: " + this + " " + that);`
    /// (the variant does not carry the two intervals).
    NotContiguous,
    /// `IllegalStateException("There is no variation present.")`.
    NoVariation,
    /// `NullPointerException`: the right flank of a result with no variation.
    NullPointer,
    /// `ArrayIndexOutOfBoundsException`: an event before the reference context.
    IndexOutOfBounds,
    /// A refusal from [`AssemblyRegion`]'s own `trim`.
    Region(RegionError),
}

/// `AssemblyRegionTrimmer`.
#[derive(Debug, Clone)]
pub struct AssemblyRegionTrimmer {
    arguments: TrimmerArguments,
    /// The length of the region's contig, for `expandWithinContig`. The reference holds the
    /// whole `SAMSequenceDictionary` and looks the contig up on each call; the trimmer only
    /// ever expands spans on the region's contig, so its length is what is kept.
    contig_length: i32,
}

/// `AssemblyRegionTrimmer.Result`: the variant span and its padded span, both null without
/// variation.
///
/// The reference's `Result` is an inner class holding the original region and reaching the
/// trimmer's arguments through its enclosing instance; here the region and the trimmer are
/// passed to the methods that need them instead.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrimResult {
    /// `variantSpan`; `isVariationPresent()` is whether this is `Some`.
    pub variant_span: Option<SimpleInterval>,
    /// `paddedSpan`.
    pub padded_span: Option<SimpleInterval>,
}

impl AssemblyRegionTrimmer {
    /// The constructor, which runs `assemblyRegionArgs.validate()` on the fields it reads. The
    /// region size, read count and probability checks are the walker's and are made there.
    ///
    /// The three checks below are `validate()`'s, in its order:
    ///
    /// ```java
    /// if ( assemblyRegionPadding < 0 ) {
    ///     throw new CommandLineException.BadArgumentValue("assemblyRegionPadding must be >= 0");
    /// }
    /// ```
    ///
    /// and the two padding checks, which use the two-argument `BadArgumentValue(name, value)`:
    /// `"" + snpPaddingForGenotyping + "< 0"`, with no space before `<`.
    pub fn new(
        arguments: TrimmerArguments,
        contig_length: i32,
    ) -> Result<AssemblyRegionTrimmer, TrimmerError> {
        if arguments.assembly_region_padding < 0 {
            return Err(TrimmerError::BadArgument(
                "assemblyRegionPadding must be >= 0".to_string(),
            ));
        }
        if arguments.snp_padding_for_genotyping < 0 {
            return Err(TrimmerError::BadArgument(format!(
                "Argument paddingAroundSNPs has a bad value: {}< 0",
                arguments.snp_padding_for_genotyping
            )));
        }
        if arguments.indel_padding_for_genotyping < 0 {
            return Err(TrimmerError::BadArgument(format!(
                "Argument paddingAroundIndels has a bad value: {}< 0",
                arguments.indel_padding_for_genotyping
            )));
        }
        Ok(AssemblyRegionTrimmer {
            arguments,
            contig_length,
        })
    }

    /// `trim(region, events, referenceContext)`, the reference context being `reference_bases`
    /// from `reference_start` on.
    ///
    /// `events` is the reference's `SortedSet<Event>`, given in its iteration order. The variant
    /// span is the events overlapping the active span, cut to it; the padded span is each of
    /// those events padded as the module note says, spanned, and cut to the region's padded span.
    pub fn trim(
        &self,
        region: &AssemblyRegion,
        events: &[Event],
        reference_bases: &[u8],
        reference_start: i32,
    ) -> Result<TrimResult, TrimmerError> {
        if self.arguments.enable_legacy_assembly_region_trimming {
            return self.trim_legacy(region, events);
        }
        let span = region.span();
        let in_region: Vec<&Event> = events
            .iter()
            .filter(|e| span.overlaps(e.contig(), e.start(), e.end()))
            .collect();
        if in_region.is_empty() {
            return Ok(no_variation());
        }
        // `.min()` answers an `Option` (nothing for an empty list); the list was just tested
        // non-empty, so `expect` cannot fire. Java's `.getAsInt()`.
        let mut min_start = in_region.iter().map(|e| e.start()).min().expect("an event");
        let mut max_end = in_region.iter().map(|e| e.end()).max().expect("an event");
        let variant_span = intersect(&interval(&span.contig, min_start, max_end)?, span)?;
        for event in &in_region {
            let mut padding = self.arguments.snp_padding_for_genotyping;
            if event.is_indel() {
                padding = self.arguments.indel_padding_for_genotyping;
                // `numRepeatsAndUnit != null && numRepeatsAndUnit.getRight() != null`: both are the
                // `None` here.
                if let Some((counts, unit)) =
                    tandem_repeat_units(event, reference_bases, reference_start)?
                {
                    let most_repeats = counts.iter().copied().max().unwrap_or(0);
                    padding = self.arguments.str_padding_for_genotyping
                        + most_repeats * unit.len() as i32;
                }
            }
            // `Math.min(minStart, Math.max(event.getStart() - padding,1))`: the start is held at 1,
            // the end is not held at the contig's end; the cut to the padded span bounds it.
            min_start = min_start.min((event.start() - padding).max(1));
            max_end = max_end.max(event.end() + padding);
        }
        let padded = intersect(
            &interval(&span.contig, min_start, max_end)?,
            region.padded_span(),
        )?;
        // `HaplotypeCallerGenotypingDebugger.println`, which writes only when the debugger is
        // enabled.
        crate::genotyper_debug::println(&format!(
            "Padded and trimmed the region to this span: {}:{}-{}",
            padded.contig, padded.start, padded.end
        ));
        result(Some(variant_span), Some(padded))
    }

    /// `trimLegacy(originalRegion, allEventsInExtendedRegion)`.
    ///
    /// ```java
    /// final SimpleInterval maximumSpan = originalRegionRange.expandWithinContig(
    ///         assemblyRegionArgs.maxExtensionIntoRegionPadding, sequenceDictionary);
    /// final SimpleInterval idealSpan =
    ///         variantSpan.expandWithinContig(padding, sequenceDictionary);
    /// final SimpleInterval finalSpan =
    ///         maximumSpan.intersect(idealSpan).mergeWithContiguous(variantSpan);
    /// ```
    ///
    /// One padding for the whole span, the indel padding if any overlapping event is not a SNP.
    /// The variant span is not cut to the active span here, unlike the default trim, so it may
    /// reach past it by an event that starts inside and ends outside.
    pub fn trim_legacy(
        &self,
        region: &AssemblyRegion,
        events: &[Event],
    ) -> Result<TrimResult, TrimmerError> {
        if events.is_empty() {
            return Ok(no_variation());
        }
        let range = region.span();
        let mut found_non_snp = false;
        let mut variant_span: Option<SimpleInterval> = None;
        let mut within: Vec<&Event> = Vec::new();
        for event in events {
            let location = interval(event.contig(), event.start(), event.end())?;
            if range.overlaps_interval(&location) {
                within.push(event);
                found_non_snp = found_non_snp || !event.is_snp();
                variant_span = Some(match variant_span {
                    None => location,
                    Some(span) => span_with(&span, &location)?,
                });
            }
        }
        let padding = if found_non_snp {
            self.arguments.indel_padding_for_genotyping
        } else {
            self.arguments.snp_padding_for_genotyping
        };
        let Some(variant_span) = variant_span else {
            return Ok(no_variation());
        };
        let maximum =
            self.expand_within_contig(range, self.arguments.max_extension_into_region_padding)?;
        let ideal = self.expand_within_contig(&variant_span, padding)?;
        let final_span = merge_with_contiguous(&intersect(&maximum, &ideal)?, &variant_span)?;
        // The reference's debug lines, each printing the object's `toString()`.
        if crate::genotyper_debug::is_enabled() {
            let show = |i: &SimpleInterval| format!("{}:{}-{}", i.contig, i.start, i.end);
            let events: Vec<String> = within.iter().map(|e| e.java_string()).collect();
            crate::genotyper_debug::println(&format!("events       : [{}]", events.join(", ")));
            crate::genotyper_debug::println(&format!(
                "region       : AssemblyRegion {} active?={} nReads={}",
                show(range),
                region.is_active(),
                region.reads().len()
            ));
            crate::genotyper_debug::println(&format!("callableSpan : {}", show(&variant_span)));
            crate::genotyper_debug::println(&format!("padding      : {padding}"));
            crate::genotyper_debug::println(&format!("idealSpan    : {}", show(&ideal)));
            crate::genotyper_debug::println(&format!("maximumSpan  : {}", show(&maximum)));
            crate::genotyper_debug::println(&format!("finalSpan    : {}", show(&final_span)));
        }
        result(Some(variant_span), Some(final_span))
    }

    /// `SimpleInterval.expandWithinContig(padding, dictionary)`, which pads both ends and clamps
    /// them to the contig.
    ///
    /// The `None` it can answer is `trimIntervalToContig`'s `null`, which happens only when the
    /// span starts past the contig's end; the reference would carry that `null` on to a
    /// `NullPointerException` at the next `intersect`, where this reports an
    /// `IllegalArgumentException` instead. No span on the region's contig reaches it.
    fn expand_within_contig(
        &self,
        span: &SimpleInterval,
        padding: i32,
    ) -> Result<SimpleInterval, TrimmerError> {
        span.expand_within_contig(padding, self.contig_length)
            .ok_or_else(|| TrimmerError::IllegalArgument("expandWithinContig".to_string()))
    }

    /// `Result.getVariantRegion()`: the original region trimmed to the two spans.
    ///
    /// ```java
    /// Utils.validate(isVariationPresent(), "There is no variation present.");
    /// return originalRegion.trim(variantSpan, paddedSpan);
    /// ```
    pub fn variant_region(
        &self,
        result: &TrimResult,
        region: &AssemblyRegion,
        header: &SamHeader,
    ) -> Result<AssemblyRegion, TrimmerError> {
        let (Some(variant), Some(padded)) = (&result.variant_span, &result.padded_span) else {
            return Err(TrimmerError::NoVariation);
        };
        region
            .trim(variant, padded, header)
            .map_err(TrimmerError::Region)
    }

    /// `Result.nonVariantLeftFlankRegion()`: the whole region without variation, the part before
    /// the variant span padded by `assemblyRegionPadding`, or nothing.
    ///
    /// The outer `Option` is Java's `Optional<AssemblyRegion>`. Without variation the reference
    /// returns `Optional.of(originalRegion)`, the same object; here it is a copy.
    pub fn non_variant_left_flank_region(
        &self,
        result: &TrimResult,
        region: &AssemblyRegion,
        header: &SamHeader,
    ) -> Result<Option<AssemblyRegion>, TrimmerError> {
        let Some(variant) = &result.variant_span else {
            return Ok(Some(region.clone()));
        };
        let span = region.span();
        if span.start < variant.start {
            let flank = interval(&span.contig, span.start, variant.start - 1)?;
            region
                .trim_with_padding(&flank, self.arguments.assembly_region_padding, header)
                // `.map(Some)` wraps a successful region in `Some`, Java's `Optional.of`.
                .map(Some)
                .map_err(TrimmerError::Region)
        } else {
            Ok(None)
        }
    }

    /// `Result.nonVariantRightFlankRegion()`, which dereferences the variant span unchecked.
    ///
    /// ```java
    /// if (variantSpan.getEnd() < originalRegion.getEnd()) {
    /// ```
    ///
    /// Unlike the left flank there is no `isVariationPresent()` test first, so a result without
    /// variation is a `NullPointerException`, not the whole region.
    pub fn non_variant_right_flank_region(
        &self,
        result: &TrimResult,
        region: &AssemblyRegion,
        header: &SamHeader,
    ) -> Result<Option<AssemblyRegion>, TrimmerError> {
        let variant = result
            .variant_span
            .as_ref()
            .ok_or(TrimmerError::NullPointer)?;
        let span = region.span();
        if variant.end < span.end {
            let flank = interval(&span.contig, variant.end + 1, span.end)?;
            region
                .trim_with_padding(&flank, self.arguments.assembly_region_padding, header)
                .map(Some)
                .map_err(TrimmerError::Region)
        } else {
            Ok(None)
        }
    }
}

/// `noVariation(region)`: `new Result(targetRegion, null, null)`.
fn no_variation() -> TrimResult {
    TrimResult {
        variant_span: None,
        padded_span: None,
    }
}

/// The `Result` constructor's check:
///
/// ```java
/// Utils.validateArg(paddedSpan == null || variantSpan == null
///         || paddedSpan.contains(variantSpan), "the padded span must include the variant span");
/// ```
fn result(
    variant_span: Option<SimpleInterval>,
    padded_span: Option<SimpleInterval>,
) -> Result<TrimResult, TrimmerError> {
    if let (Some(v), Some(p)) = (&variant_span, &padded_span) {
        if !p.contains(v) {
            return Err(TrimmerError::IllegalArgument(
                "the padded span must include the variant span".to_string(),
            ));
        }
    }
    Ok(TrimResult {
        variant_span,
        padded_span,
    })
}

/// The repeat count of each allele, and the repeat unit.
type RepeatUnits = (Vec<i32>, Vec<u8>);

/// `TandemRepeat.getNumTandemRepeatUnits(referenceContext, event)`: the context after the event's
/// padding base, from the reference context's bases.
///
/// ```java
/// final int startIndex = event.getStart() + 1 - ref.getWindow().getStart();
/// return GATKVariantContextUtils.getNumTandemRepeatUnits(event.refAllele(),
///         Collections.singletonList(event.altAllele()),
///         Arrays.copyOfRange(refBases, startIndex, refBases.length));
/// ```
///
/// The two refusals are `Arrays.copyOfRange`'s: a negative start is out of bounds, and a start
/// past the end is `from > to`, an `IllegalArgumentException`. A start exactly at the end is an
/// empty context, not a refusal. `std::slice::from_ref` makes a one-element list of the
/// alternate allele without copying it, `Collections.singletonList`.
fn tandem_repeat_units(
    event: &Event,
    reference_bases: &[u8],
    reference_start: i32,
) -> Result<Option<RepeatUnits>, TrimmerError> {
    let start_index = event.start() + 1 - reference_start;
    if start_index < 0 {
        return Err(TrimmerError::IndexOutOfBounds);
    }
    let start_index = start_index as usize;
    if start_index > reference_bases.len() {
        return Err(TrimmerError::IllegalArgument(format!(
            "{start_index} > {}",
            reference_bases.len()
        )));
    }
    Ok(num_tandem_repeat_units(
        event.ref_allele(),
        std::slice::from_ref(event.alt_allele()),
        &reference_bases[start_index..],
    ))
}

/// `new SimpleInterval(contig, start, end)`, which refuses an end before the start.
fn interval(contig: &str, start: i32, end: i32) -> Result<SimpleInterval, TrimmerError> {
    SimpleInterval::new(contig, start, end).ok_or_else(|| {
        TrimmerError::IllegalArgument(format!(
            "Invalid interval. Contig:{contig} start:{start} end:{end}"
        ))
    })
}

/// `SimpleInterval.intersect`, which refuses intervals that do not overlap.
///
/// The reference's message goes on to name both intervals,
/// `"SimpleInterval::intersect(): The two intervals need to overlap " + this + " " + that`;
/// the one carried here stops before them.
fn intersect(a: &SimpleInterval, b: &SimpleInterval) -> Result<SimpleInterval, TrimmerError> {
    a.intersect(b).ok_or_else(|| {
        TrimmerError::IllegalArgument(
            "SimpleInterval::intersect(): The two intervals need to overlap".to_string(),
        )
    })
}

/// `SimpleInterval.spanWith`.
fn span_with(a: &SimpleInterval, b: &SimpleInterval) -> Result<SimpleInterval, TrimmerError> {
    if a.contig != b.contig {
        return Err(TrimmerError::IllegalArgument(
            "Cannot get span for intervals on different contigs".to_string(),
        ));
    }
    interval(&a.contig, a.start.min(b.start), a.end.max(b.end))
}

/// `SimpleInterval.mergeWithContiguous`: the span of two intervals that overlap or touch, and a
/// `GATKException` for two that do not.
fn merge_with_contiguous(
    a: &SimpleInterval,
    b: &SimpleInterval,
) -> Result<SimpleInterval, TrimmerError> {
    let contiguous = a.contig == b.contig && a.start <= b.end + 1 && b.start <= a.end + 1;
    if !contiguous {
        return Err(TrimmerError::NotContiguous);
    }
    span_with(a, b)
}
