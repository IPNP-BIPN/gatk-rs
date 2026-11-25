//! `GVCFWriter`'s block combining, and the reblocking variant `ReblockGVCF` writes through.
//!
//! Ported from `org.broadinstitute.hellbender.utils.variant.writers` (GATK 4.7.0.0; every file
//! named here is unchanged since 4.6.2.0): `GVCFBlockCombiner`, `HomRefBlock`, `GVCFBlock`,
//! `ReblockingGVCFBlockCombiner` and the two writers around them, and Mutect2's
//! `SomaticGVCFBlockCombiner` and `TLODBlock`. The combiner takes records one at a time and hands
//! back the records to write, in order: a hom-ref record joins the open block when its GQ falls in
//! the block's band, and anything else closes the block and is written as it is.
//!
//! # The reblocking buffer
//!
//! `ReblockingGVCFBlockCombiner` holds the reference records back in a buffer before they reach
//! the combining, so that a variant arriving later can cut them: a block overlapping a variant is
//! trimmed to end before it, and one that runs past it is split, its tail moved to start after it
//! with the reference base of its new first position. `vcfOutputEnd` is the last position written,
//! and a reference record starting at or before it is moved past it or dropped.
//!
//! # Posteriors
//!
//! A genotype carrying `PP` switches `HomRefBlock` to phred-scaled posteriors for its GQ: the
//! block keeps the element-wise minimum of the `PP`s as it does of the `PL`s, its `GQ` is taken
//! from that minimum instead, and the block record carries it as `PP`.
//!
//! # How a band is chosen
//!
//! `parse_partitions` turns the thresholds `[A, B, C]` into the bands `[0, A)`,
//! `[A, B)`, `[B, C)` and, when `C <= 99`, `[C, 100)`. A block's band is fixed by its FIRST
//! record: that record's `GQ` field, capped at 99, or 0 when it has none. The GQ the block later
//! WRITES is something else (see below), so a block's `GQ` need not be the one its band was
//! chosen by.
//!
//! # Where a block breaks
//!
//! The open block is closed and written, and a new one started, when the next record:
//!
//! - is not within one base of it (`vc.withinDistanceOf(block, 1)`: a gap, or another contig);
//! - has a `GQ` field (capped at 99) outside the band. A MISSING `GQ` reads as -1, which is in no
//!   band, so in the plain combiner every record without a `GQ` starts a block of its own (in
//!   band `[0, A)`);
//! - has another ploidy, or a `PL` of another length when both have one;
//! - is not a hom-ref (or a no-call whose first PL is 0) with exactly one ALT, `<NON_REF>`. Such
//!   a record closes the block and is written as it is.
//!
//! `withinDistanceOf` admits an OVERLAPPING record. If it would otherwise merge (same band,
//! ploidy and `PL` length), `HomRefBlock.add` then throws "overlaps previous end" on it; if not,
//! it simply starts a new block.
//! And after a variant, hom-ref records on its contig that start and end at or before the
//! variant's end are dropped, not blocked, until the first one that does not.
//!
//! # What a block record carries
//!
//! The first record's site (contig, start, ID, alleles, QUAL, FILTER) with its INFO cleared to
//! one `END`, the last merged record's `END` (or its start when it had none). One genotype,
//! rebuilt from nothing: ploidy copies of the reference allele, then
//!
//! - `PL`: the element-wise minimum of every merged `PL`;
//! - `GQ`: recomputed from that minimum `PL` (second smallest less smallest; from the minimum
//!   `PP` when there is one), or the smallest `GQ` field when no record had a `PL`. The genotype
//!   may hold more than 99; htsjdk caps it only as it writes the line (`GQAccessor`:
//!   `Math.min(g.getGQ(), VCFConstants.MAX_GENOTYPE_QUAL)`);
//! - `DP`: the median of the records' `DP`s, rounded half up. One value per RECORD, not per
//!   position, so a merged 1000-base record counts once;
//! - `MIN_DP`: the smallest of the same `DP`s. An input record's own `MIN_DP` is not read.
//!
//! With `--floor-blocks` the genotype has no `PL` and no `MIN_DP`, and its `GQ` is the band's lower
//! bound. The order the fields are set in does not matter for the written line: htsjdk's
//! `VariantContext.calcVCFGenotypeKeys` puts `GT` first and sorts the rest, so a block reads
//! `GT:DP:GQ:MIN_DP:PL` (`:PP` after it when there are posteriors).
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is "a `T` or nothing", Java's nullable field checked by the compiler. `Some(x)`
//!   is a value, `None` the `null`. `.as_ref()` and `.as_mut()` look inside without taking it,
//!   `.take()` moves the value out and leaves `None`, `is_some_and(|x| ..)` is
//!   `x != null && ..`.
//! - `Result<T, EngineError>` is "a `T` or the exception the reference throws"; `?` after a call
//!   passes the error up, Java's `throws`.
//! - `match (a, b) { (Some(x), _) => .., (None, Some(y)) => .., .. }` tests several values at
//!   once, top arm first, and must cover every case; `_` matches anything.
//! - `.clone()` is a deep copy, where the Java builds a new object from another (`new
//!   VariantContextBuilder(vc)`). The Java's buffer holds mutable builders and trims them through
//!   references; here the buffer owns its records and changes them in place.
//! - `&mut dyn ReferenceBase` is a writable reference to "any type implementing this trait",
//!   Java's parameter of an interface type.
//! - `x as i32` and `x as usize` are casts. From a float they truncate toward zero and saturate
//!   at the type's bounds, with NaN going to 0, which is what Java's `(int)` cast does too.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotypes_context::GenotypesContext;
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

use crate::genotyping_engine::EngineError;

/// `VCFConstants.MAX_GENOTYPE_QUAL`.
pub const MAX_GENOTYPE_QUAL: i32 = 99;
/// `GVCFWriter.GVCF_BLOCK`.
pub const GVCF_BLOCK: &str = "GVCFBlock";
/// `GATKVCFConstants.MIN_DP_FORMAT_KEY`.
pub const MIN_DP_FORMAT_KEY: &str = "MIN_DP";
/// `GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY`.
const PHRED_SCALED_POSTERIORS_KEY: &str = "PP";

/// The reference's exception classes, as `EngineError::Runtime` carries them: the fully qualified
/// class name and the message, so the caller can print what GATK prints.
fn illegal_argument(message: String) -> EngineError {
    EngineError::Runtime {
        class: "java.lang.IllegalArgumentException".to_string(),
        message,
    }
}

fn should_never_reach_here(message: String) -> EngineError {
    EngineError::Runtime {
        class:
            "org.broadinstitute.hellbender.exceptions.GATKException$ShouldNeverReachHereException"
                .to_string(),
        message,
    }
}

fn gatk_exception(message: String) -> EngineError {
    EngineError::Runtime {
        class: "org.broadinstitute.hellbender.exceptions.GATKException".to_string(),
        message,
    }
}

/// `allele.equals(Allele.NON_REF_ALLELE)`: the symbolic `<NON_REF>`.
fn is_non_ref(allele: &Allele) -> bool {
    allele.is_symbolic() && allele.display_string() == "<NON_REF>"
}

/// `vc.getAttributeAsInt(END, default)`, the default left to the caller (`None` is "absent").
///
/// htsjdk's `CommonInfo.getAttributeAsInt` (5.0.0):
///
/// ```java
/// if (x == null || x == VCFConstants.MISSING_VALUE_v4) return defaultValue;
/// if (x instanceof Integer) return (Integer) x;
/// return Integer.parseInt((String) x); // throws an exception if this isn't a string
/// ```
///
/// This is more forgiving than that: a `Double`, a one-element list, surrounding spaces and an
/// unparsable string are read or treated as absent here where the Java would throw. An integer
/// or a plain integer string reads the same in both.
pub fn end_attribute(vc: &VariantContext) -> Option<i64> {
    vc.attributes
        .iter()
        .find(|(key, _)| key == "END")
        .and_then(|(_, value)| match value {
            Value::Int(end) => Some(*end),
            Value::Double(end) => Some(*end as i64),
            Value::Str(text) => text.trim().parse().ok(),
            Value::List(values) if values.len() == 1 => match &values[0] {
                Value::Int(end) => Some(*end),
                Value::Str(text) => text.trim().parse().ok(),
                _ => None,
            },
            _ => None,
        })
}

/// `builder.attribute(END, value)`: replaced in place, or added at the end.
///
/// Where it lands does not matter for the written line: htsjdk's `VCFEncoder` collects the INFO
/// fields in a `TreeMap` and writes them sorted. `iter_mut().find(..)` gives a writable reference
/// to the matching pair, and `*slot = value` writes through it.
pub fn put_attribute(vc: &mut VariantContext, key: &str, value: Value) {
    match vc.attributes.iter_mut().find(|(name, _)| name == key) {
        Some((_, slot)) => *slot = value,
        None => vc.attributes.push((key.to_string(), value)),
    }
}

/// `GATKVariantContextUtils.calculateGQFromPLs`: the second smallest PL less the smallest.
///
/// Kept as the Java's single pass (seed with the first two, then `candidate >= second` skips and
/// `candidate <= first` demotes `first`) rather than a sort; the result is the same either way,
/// and a tie for smallest gives 0. The refusal is the Java's
/// `Utils.validateArg(plValues.length >= 2, ..)`.
pub fn calculate_gq_from_pls(pls: &[i32]) -> Result<i32, EngineError> {
    if pls.len() < 2 {
        return Err(illegal_argument(
            "Array of PL values must contain at least two elements.".to_string(),
        ));
    }
    let (mut first, mut second) = (pls[0], pls[1]);
    if first > second {
        second = first;
        first = pls[1];
    }
    for &candidate in &pls[2..] {
        if candidate >= second {
            continue;
        }
        if candidate <= first {
            second = first;
            first = candidate;
        } else {
            second = candidate;
        }
    }
    Ok(second - first)
}

/// `MathUtils.median(Collection)`: commons-math `Median`, the legacy percentile estimate, so an
/// even count is the mean of its two middle values.
///
/// `GVCFBlock.getMedianDP` is the caller, and its Javadoc states the even case ("the median is
/// the mean of the two middle values"):
///
/// ```java
/// return (int) Math.round(MathUtils.median(DPs));
/// ```
///
/// **How**: sort, take the 1-based position `(n + 1) / 2`, and interpolate between the values
/// either side of it. For whole depths that is the middle value, or the mean of the two middle
/// values (`20, 22` gives 21, and `20, 21` gives 20.5, which the caller rounds to 21).
fn median(values: &[i32]) -> f64 {
    // `f64::from` widens each `int` to `double` exactly. A float has no total order in Rust
    // (NaN), so the sort takes a comparison that may fail; depths are whole numbers and never do.
    let mut sorted: Vec<f64> = values.iter().map(|value| f64::from(*value)).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite depths"));
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    // `floor as usize - 1` is the 0-based index of the value at or below the position.
    let position = (n as f64 + 1.0) * 0.5;
    let floor = position.floor();
    let fraction = position - floor;
    let lower = sorted[floor as usize - 1];
    if floor as usize >= n {
        return sorted[n - 1];
    }
    let upper = sorted[floor as usize];
    lower + fraction * (upper - lower)
}

/// `GVCFBlockCombiner.parsePartitions`: the bands, each closed below and open above, with a last
/// one up to `MAX_GENOTYPE_QUAL + 1` when the list stops short of it.
///
/// The Java keys a Guava `RangeMap` by the bands so that `createNewBlock` can look a GQ up; here
/// the bands are a list of `(lower, upper)` pairs and the lookup is a linear search.
///
/// Two edges follow from the threshold starting at 0:
///
/// ```java
/// int lastThreshold = 0;
/// ...
/// } else if (value == lastThreshold) {
///     throw new IllegalArgumentException(String.format(
///             "The value %d appears more than once in the list of GQ partitions.", value));
/// ```
///
/// so a band list starting with `0` is refused as a repeat, though it appears once; and `100` is
/// accepted, and then `if (lastThreshold <= MAX_GENOTYPE_QUAL)` adds no tail band after it.
pub fn parse_partitions(bands: &[i32]) -> Result<Vec<(i32, i32)>, String> {
    if bands.is_empty() {
        return Err("The collection is empty: collection must not be null or empty.".to_string());
    }
    let mut result = Vec::new();
    let mut last = 0;
    for &value in bands {
        if value < 0 {
            return Err("The list of GQ partitions contains a non-positive integer.".to_string());
        } else if value > MAX_GENOTYPE_QUAL + 1 {
            return Err(format!(
                "The value {value} in the list of GQ partitions is greater than \
                 VCFConstants.MAX_GENOTYPE_QUAL + 1 = {}.",
                MAX_GENOTYPE_QUAL + 1
            ));
        } else if value < last {
            return Err(format!(
                "The list of GQ partitions is out of order. Previous value is {last} but the next \
                 is {value}."
            ));
        } else if value == last {
            return Err(format!(
                "The value {value} appears more than once in the list of GQ partitions."
            ));
        }
        result.push((last, value));
        last = value;
    }
    if last <= MAX_GENOTYPE_QUAL {
        result.push((last, MAX_GENOTYPE_QUAL + 1));
    }
    Ok(result)
}

/// `GVCFBlockCombiner.addRangesToHeader`: `END`, `MIN_DP` and one `GVCFBlock` line per band.
///
/// Each band line is `rangeToVCFHeaderLine`:
///
/// ```java
/// final String key = String.format(GVCF_BLOCK+"%d-%d", genotypeQualityBand.lowerEndpoint(),
///         genotypeQualityBand.upperEndpoint());
/// return new VCFHeaderLine(key, "minGQ=" + genotypeQualityBand.lowerEndpoint()
///         + "(inclusive),maxGQ=" + genotypeQualityBand.upperEndpoint() + "(exclusive)");
/// ```
///
/// so `##GVCFBlock20-60=minGQ=20(inclusive),maxGQ=60(exclusive)`. The band is in the KEY because
/// htsjdk keeps only one line per key.
///
/// htsjdk's `VCFHeader.addMetaDataLine` keeps the FIRST line for an INFO or FORMAT ID (or for
/// any other line's key) and drops a later one. The test here is whole-line equality, which
/// agrees whenever an existing `END` or `MIN_DP` line is the standard one; a header that already
/// declares either with other text would keep both here.
pub fn add_ranges_to_header(header: &mut VcfHeader, partitions: &[(i32, i32)]) {
    // A closure that captures `header` mutably, so the lines below each go through one check.
    let mut add = |line: HeaderLine| {
        if !header.lines.contains(&line) {
            header.lines.push(line);
        }
    };
    add(htsjdk_vcf::standard_header_lines::standard_info_line("END").expect("END is standard"));
    add(HeaderLine::Compound {
        key: "FORMAT".to_string(),
        id: MIN_DP_FORMAT_KEY.to_string(),
        number: htsjdk_vcf::header::Cardinality::Fixed(1),
        line_type: htsjdk_vcf::header::LineType::Integer,
        description: "Minimum DP observed within the GVCF block".to_string(),
        extra: Vec::new(),
    });
    for (lower, upper) in partitions {
        add(HeaderLine::Unstructured {
            key: format!("{GVCF_BLOCK}{lower}-{upper}"),
            value: format!("minGQ={lower}(inclusive),maxGQ={upper}(exclusive)"),
        });
    }
}

/// `HomRefBlock`, with the fields of its parent `GVCFBlock` folded in.
#[derive(Debug, Clone)]
struct HomRefBlock {
    /// `startingVC`: the first record, whose site the block record is written at.
    starting: VariantContext,
    /// The band, `[min_gq_bound, max_gq_bound)`. `GVCFBlock` calls these `minGQ` and `maxGQ`, a
    /// different thing from `HomRefBlock.minGQ` below.
    min_gq_bound: i32,
    max_gq_bound: i32,
    /// `startingVC.getReference()`, repeated ploidy times as the block's genotype.
    reference: Allele,
    /// `DPs`: one `max(DP, 0)` per merged record that had a `DP`.
    dps: Vec<i32>,
    /// The last position covered, the `END` written.
    end: i64,
    ploidy: usize,
    /// The element-wise minimum of the `PL`s; `None` until a record has one.
    min_pls: Option<Vec<i32>>,
    /// The element-wise minimum of the `PP`s, which decides `GQ` where there is one.
    min_pps: Option<Vec<i32>>,
    /// `HomRefBlock.minGQ`, the `GQ` written: from `min_pps`, else `min_pls`, else the smallest
    /// `GQ` field seen. `-1` when no genotype has had a PL or a GQ, which writes no `GQ`.
    min_gq: i32,
}

/// `PosteriorProbabilitiesUtils.parsePosteriorsIntoPhredSpace` for a genotype that carries `PP`:
/// its values, or its `PL` when they start with a `.`.
///
/// ```java
/// //samples not in trios will have PP tag like ".,.,." if family priors are applied
/// return PPstring.charAt(0)=='.' ? genotype.getPL() :
///         Arrays.stream(PPstring.split(",")).mapToInt(i -> Integer.parseInt(i)).toArray();
/// ```
///
/// `None` is the genotype without `PP`, which the Java callers test with
/// `hasExtendedAttribute` before calling. One difference: an entry that is not an integer reads
/// as 0 here, where `Integer.parseInt` throws.
fn posteriors(genotype: &Genotype) -> Option<Vec<i32>> {
    // `?` on an `Option` returns `None` from the function when the key is absent.
    let (_, value) = genotype
        .extended
        .iter()
        .find(|(key, _)| key == PHRED_SCALED_POSTERIORS_KEY)?;
    let text = value.format().unwrap_or_default();
    if text.starts_with('.') {
        return genotype.pl.clone();
    }
    Some(
        text.split(',')
            .map(|part| part.trim().parse().unwrap_or(0))
            .collect(),
    )
}

impl HomRefBlock {
    /// `new HomRefBlock(startingVC, lowerGQBound, upperGQBound, defaultPloidy)`.
    ///
    /// The first record seeds the minima, and its GQ is computed the way `add` recomputes it:
    ///
    /// ```java
    /// if (minPPs != null) {
    ///     minGQ = GATKVariantContextUtils.calculateGQFromPLs(minPPs);
    /// } else if (minPLs != null) {
    ///     minGQ = GATKVariantContextUtils.calculateGQFromPLs(minPLs);
    /// } else if (g.hasGQ()) {
    ///     minGQ = g.getGQ();
    /// }
    /// ```
    ///
    /// The constructor's own checks (`upperGQBound <= 100`, `lower <= upper`) are not repeated:
    /// the callers build the bands with `parse_partitions`, which cannot produce one that fails
    /// them.
    fn new(
        starting: &VariantContext,
        lower: i32,
        upper: i32,
        default_ploidy: usize,
    ) -> Result<HomRefBlock, EngineError> {
        let g = &starting.genotypes[0];
        // `startingVC.getMaxPloidy(defaultPloidy)`.
        let max = starting
            .genotypes
            .iter()
            .map(Genotype::ploidy)
            .max()
            .unwrap_or(0);
        // `getMaxPloidy` falls back to the default only when every genotype has ploidy 0.
        let ploidy = if max == 0 { default_ploidy } else { max };
        let min_pls = g.pl.clone();
        let min_pps = posteriors(g);
        let min_gq = match (&min_pps, &min_pls, g.gq) {
            (Some(pps), _, _) => calculate_gq_from_pls(pps)?,
            (None, Some(pls), _) => calculate_gq_from_pls(pls)?,
            (None, None, Some(gq)) => gq,
            (None, None, None) => -1,
        };
        Ok(HomRefBlock {
            starting: starting.clone(),
            min_gq_bound: lower,
            max_gq_bound: upper,
            reference: starting.alleles[0].clone(),
            // `GVCFBlock`'s constructor: `if (g.hasDP()) DPs.add(Math.max(g.getDP(), 0))`. An
            // `Option` iterates as zero or one items, so this is a list of zero or one depth.
            dps: g.dp.map(|dp| dp.max(0)).into_iter().collect(),
            end: end_attribute(starting).unwrap_or(starting.start),
            ploidy,
            min_pls,
            min_pps,
            min_gq,
        })
    }

    /// `GVCFBlock.withinBounds`: `GQ >= minGQ && GQ < maxGQ`.
    fn within_bounds(&self, gq: i32) -> bool {
        gq >= self.min_gq_bound && gq < self.max_gq_bound
    }

    /// `HomRefBlock.add(pos, newEnd, genotype)`: extend the block to `new_end` and fold the
    /// genotype's `PL`, `PP`, `GQ` and `DP` in.
    ///
    /// The checks run in the Java's order, so the same record raises the same message. The band
    /// test caps the GQ first:
    ///
    /// ```java
    /// // Make sure the GQ is within the bounds of this band. Treat GQs > 99 as 99.
    /// if ( !withinBounds(Math.min(genotype.getGQ(), VCFConstants.MAX_GENOTYPE_QUAL))) {
    /// ```
    fn add(&mut self, position: i64, new_end: i64, genotype: &Genotype) -> Result<(), EngineError> {
        if position > self.end + 1 {
            return Err(illegal_argument(format!(
                "adding genotype at pos {position} isn't contiguous with previous end {}",
                self.end
            )));
        }
        if position < self.end + 1 {
            return Err(illegal_argument(format!(
                "adding genotype at pos {position} overlaps previous end {}",
                self.end
            )));
        }
        if genotype.ploidy() != self.ploidy {
            return Err(illegal_argument(format!(
                "cannot add a genotype with a different ploidy: {} != {}",
                genotype.ploidy(),
                self.ploidy
            )));
        }
        // `getGQ()` is -1 when the genotype has none.
        let gq = genotype.gq.unwrap_or(-1);
        if !self.within_bounds(gq.min(MAX_GENOTYPE_QUAL)) {
            return Err(illegal_argument(format!(
                "cannot add a genotype with GQ={gq} because it's not within bounds [{},{})",
                self.min_gq_bound, self.max_gq_bound
            )));
        }
        // `if (minPLs == null) minPLs = genotype.getPL()`, which may still be null; otherwise the
        // element-wise minimum, a record without a PL leaving it alone. `zip` walks the two
        // lists side by side, and `*slot = ..` writes into the block's own list.
        match (&mut self.min_pls, &genotype.pl) {
            (None, pls) => self.min_pls = pls.clone(),
            (Some(min), Some(pls)) => {
                if pls.len() != min.len() {
                    return Err(gatk_exception(format!(
                        "trying to merge different PL array sizes: {} != {}",
                        pls.len(),
                        min.len()
                    )));
                }
                for (slot, pl) in min.iter_mut().zip(pls) {
                    *slot = (*slot).min(*pl);
                }
            }
            (Some(_), None) => {}
        }
        if let Some(pps) = posteriors(genotype) {
            match &mut self.min_pps {
                None => self.min_pps = Some(pps),
                Some(min) => {
                    if pps.len() != min.len() {
                        return Err(gatk_exception(format!(
                            "trying to merge different PP array sizes: {} != {}",
                            pps.len(),
                            min.len()
                        )));
                    }
                    for (slot, pp) in min.iter_mut().zip(&pps) {
                        *slot = (*slot).min(*pp);
                    }
                }
            }
        }
        // Recomputed from the minima after every record, so the block's GQ is the GQ of the
        // minimum PL vector, not the minimum of the records' GQ fields. Nothing caps it here; the
        // VCF encoder writes at most 99.
        self.min_gq = match (&self.min_pps, &self.min_pls) {
            (Some(pps), _) => calculate_gq_from_pls(pps)?,
            (None, Some(pls)) => calculate_gq_from_pls(pls)?,
            (None, None) if self.min_gq == -1 => gq,
            (None, None) => self.min_gq.min(gq),
        };
        self.end = new_end;
        if let Some(dp) = genotype.dp {
            self.dps.push(dp.max(0));
        }
        Ok(())
    }

    /// `GVCFBlock.toVariantContext(sampleName, floorBlocks)` with
    /// `HomRefBlock.createHomRefGenotype`.
    ///
    /// ```java
    /// final VariantContextBuilder vcb = new VariantContextBuilder(getStartingVC());
    /// vcb.attributes(new LinkedHashMap<>(2)); // clear the attributes
    /// vcb.stop(getEnd());
    /// vcb.attribute(VCFConstants.END_KEY, getEnd());
    /// ```
    ///
    /// and the genotype, built from nothing (`gb.noAD().noPL().noAttributes()`):
    ///
    /// ```java
    /// if (!floorBlocks) {
    ///     if (minPLs != null) {
    ///         gb.PL(minPLs);
    ///     }
    ///     gb.GQ(getMinGQ());
    ///     if (!DPs.isEmpty()) {
    ///         gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, getMinDP());
    ///     }
    /// }
    /// else {
    ///     gb.GQ(getGQLowerBound());
    /// }
    /// if (!DPs.isEmpty()) {
    ///     gb.DP(getMedianDP());
    /// }
    /// ```
    ///
    /// `gb.GQ(-1)` is a genotype without `GQ` (htsjdk's `hasGQ()` is `getGQ() != -1`), which is
    /// the `then_some` below: `Some(min_gq)` when the condition holds, `None` otherwise.
    fn to_variant_context(&self, sample_name: &str, floor_blocks: bool) -> VariantContext {
        let mut vc = self.starting.clone();
        vc.attributes = vec![("END".to_string(), Value::Int(self.end))];
        vc.stop = self.end;
        // `vec![x; n]` is `n` copies of `x`, `Collections.nCopies(getPloidy(), getRef())`.
        let mut genotype = Genotype::new(sample_name, vec![self.reference.clone(); self.ploidy]);
        if !floor_blocks {
            genotype.pl = self.min_pls.clone();
            genotype.gq = (self.min_gq != -1).then_some(self.min_gq);
            if let Some(min) = self.dps.iter().min() {
                genotype
                    .extended
                    .push((MIN_DP_FORMAT_KEY.to_string(), Value::Int(i64::from(*min))));
            }
        } else {
            genotype.gq = Some(self.min_gq_bound);
        }
        if !self.dps.is_empty() {
            genotype.dp =
                Some(crate::reference_confidence_merger::java_round(median(&self.dps)) as i32);
        }
        // `Utils.listFromPrimitives(minPPs)`, a `List<Integer>` that the encoder writes joined by
        // commas; the joined string here is the same text.
        if let Some(pps) = &self.min_pps {
            let text: Vec<String> = pps.iter().map(i32::to_string).collect();
            genotype.extended.push((
                PHRED_SCALED_POSTERIORS_KEY.to_string(),
                Value::Str(text.join(",")),
            ));
        }
        vc.genotypes = GenotypesContext::new(vec![genotype]);
        vc
    }

    /// `isContiguous(vc)`: `vc.withinDistanceOf(block, 1)`.
    ///
    /// htsjdk's `Locatable.withinDistanceOf` is `contigsMatch(other) && CoordMath.overlaps(..,
    /// other.getStart() - distance, other.getEnd() + distance)`, which for intervals with start at
    /// or before end is the two comparisons below: the record starts no later than the base after
    /// the block's end, and ends no earlier than the base before the block's start.
    fn is_contiguous(&self, vc: &VariantContext) -> bool {
        vc.contig == self.starting.contig
            && vc.start <= self.end + 1
            && self.starting.start - 1 <= vc.stop
    }
}

/// `ReblockingOptions`: `--drop-low-quals`, `--allow-missing-hom-ref-data` and
/// `--rgq-threshold-to-no-call` as `ReblockGVCF` passes them to the combiner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReblockingOptions {
    pub drop_low_quals: bool,
    pub allow_missing_hom_ref_data: bool,
    pub rgq_threshold: f64,
}

/// `GVCFBlockCombiner`, with `ReblockingGVCFBlockCombiner.addHomRefSite` when `reblocking` is set.
///
/// The Java has a subclass override `addHomRefSite`; here one type carries an `Option` of the
/// reblocking settings and `add_hom_ref_site` dispatches on it. The rest of the reblocking
/// combiner (its buffer) is [`ReblockingWriter`], which feeds this one.
pub struct BlockCombiner {
    partitions: Vec<(i32, i32)>,
    floor_blocks: bool,
    to_output: Vec<VariantContext>,
    /// `nextAvailableStart`: the end of the last variant written, `-1` once a later record has
    /// cleared it. Hom-ref records ending at or before it on its contig are dropped.
    next_available_start: i64,
    contig_of_next_available_start: Option<String>,
    /// The first record's sample, which every block record is written for.
    sample_name: Option<String>,
    current_block: Option<HomRefBlock>,
    reblocking: Option<ReblockingOptions>,
    /// `logger.warn` lines, in order.
    pub warnings: Vec<String>,
}

impl BlockCombiner {
    pub fn new(
        partitions: Vec<(i32, i32)>,
        floor_blocks: bool,
        reblocking: Option<ReblockingOptions>,
    ) -> BlockCombiner {
        BlockCombiner {
            partitions,
            floor_blocks,
            to_output: Vec::new(),
            next_available_start: -1,
            contig_of_next_available_start: None,
            sample_name: None,
            current_block: None,
            reblocking,
            warnings: Vec::new(),
        }
    }

    /// `genotypeCanBeMergedInCurrentBlock`:
    ///
    /// ```java
    /// return currentHomRefBlock != null
    ///         && currentHomRefBlock.withinBounds(Math.min(g.getGQ(), MAX_GENOTYPE_QUAL))
    ///         && currentHomRefBlock.getPloidy() == g.getPloidy()
    ///         && (currentHomRefBlock.getMinPLs() == null || !g.hasPL()
    ///             || (currentHomRefBlock.getMinPLs().length == g.getPL().length));
    /// ```
    ///
    /// Position is not tested here: `submit` has already closed a block the record does not
    /// touch, and `HomRefBlock::add` refuses one that overlaps.
    fn can_merge(&self, g: &Genotype) -> bool {
        self.current_block.as_ref().is_some_and(|block| {
            block.within_bounds(g.gq.unwrap_or(-1).min(MAX_GENOTYPE_QUAL))
                && block.ploidy == g.ploidy()
                && (block.min_pls.is_none()
                    || g.pl.is_none()
                    || block.min_pls.as_ref().map(Vec::len) == g.pl.as_ref().map(Vec::len))
        })
    }

    /// `createNewBlock`: the band is looked up by the record's GQ, capped at 99, or 0 when it has
    /// none (`g.hasGQ() ? Math.min(g.getGQ(), MAX_GENOTYPE_QUAL) : 0`).
    ///
    /// The refusal cannot be reached from `parse_partitions`'s bands, which cover 0 to 99 without
    /// a gap. Its message is not the Java's, which concatenates the whole genotype and record
    /// (`"GQ " + g + " from " + vc + ..`).
    fn create_new_block(
        &self,
        vc: &VariantContext,
        g: &Genotype,
    ) -> Result<HomRefBlock, EngineError> {
        let gq = g.gq.map_or(0, |gq| gq.min(MAX_GENOTYPE_QUAL));
        // `let Some(..) = .. else { .. };` binds the band when `find` found one and otherwise
        // runs the `else` block, which must leave the function. `.copied()` turns the borrowed
        // pairs into values.
        let Some((lower, upper)) = self
            .partitions
            .iter()
            .copied()
            .find(|(lower, upper)| *lower <= gq && gq < *upper)
        else {
            return Err(gatk_exception(format!(
                "GQ {gq} from {}:{} didn't fit into any partition",
                vc.contig, vc.start
            )));
        };
        HomRefBlock::new(vc, lower, upper, g.ploidy())
    }

    /// `GVCFBlockCombiner.addHomRefSite`: merge into the open block, or close it (returned, to be
    /// written) and open a new one; `Ok(None)` when nothing is ready to write.
    ///
    /// ```java
    /// if (nextAvailableStart != -1) {
    ///     ...
    ///     if ( vc.getStart() <= nextAvailableStart
    ///             && vc.getContig().equals(contigOfNextAvailableStart) ) {
    ///         if (vc.getEnd() <= nextAvailableStart) {
    ///             return null;
    ///         }
    ///     }
    ///     // otherwise, reset to non-relevant
    ///     nextAvailableStart = -1;
    ///     contigOfNextAvailableStart = null;
    /// }
    /// ```
    ///
    /// The elided line is the Java's reason: "there's a use case here related to ReblockGVCFs for
    /// overlapping deletions on different haplotypes". A dropped record leaves the marker set, so
    /// every following record inside the variant is dropped too, and the first one past it
    /// clears it.
    fn add_hom_ref_site_base(
        &mut self,
        vc: &VariantContext,
        g: &Genotype,
    ) -> Result<Option<VariantContext>, EngineError> {
        if self.next_available_start != -1 {
            if vc.start <= self.next_available_start
                && self.contig_of_next_available_start.as_deref() == Some(vc.contig.as_str())
                && vc.stop <= self.next_available_start
            {
                return Ok(None);
            }
            self.next_available_start = -1;
            self.contig_of_next_available_start = None;
        }
        if self.can_merge(g) {
            let end = end_attribute(vc).unwrap_or(vc.start);
            self.current_block
                .as_mut()
                .expect("a block to merge into")
                .add(vc.start, end, g)?;
            Ok(None)
        } else {
            // `currentBlock != null ? currentBlock.toVariantContext(..) : null`, then the new
            // block replaces it. `unwrap_or_default` is the empty string for a missing name.
            let sample = self.sample_name.clone().unwrap_or_default();
            let result = self
                .current_block
                .as_ref()
                .map(|block| block.to_variant_context(&sample, self.floor_blocks));
            self.current_block = Some(self.create_new_block(vc, g)?);
            Ok(result)
        }
    }

    /// `ReblockingGVCFBlockCombiner.addHomRefSite`, or the plain one.
    ///
    /// The reblocking version first drops low-quality reference sites, then decides what counts
    /// as hom-ref by its own rule ([`is_hom_ref_reblocking`]), and mends a hom-ref without `PL`:
    ///
    /// - with a `GQ`: warns and blocks it by that `GQ`;
    /// - with neither: refuses it as `UserException.BadInput`, or with
    ///   `--allow-missing-hom-ref-data` warns and blocks it as `GQ` 0 with `PL` `[0,0,0]`. That
    ///   is three PLs whatever the ploidy (`gBuilder.GQ(0).PL(new int[]{0,0,0})`).
    ///
    /// Anything else is dropped. The Java has one more branch between the two, for a no-call
    /// with a first PL of 0 ("some external data has no-called genotypes with good
    /// likelihoods"), but `isHomRef` is already true for any genotype whose first PL is 0, so it
    /// is never taken; it is kept below as the Java has it.
    fn add_hom_ref_site(
        &mut self,
        vc: &VariantContext,
        g: &Genotype,
    ) -> Result<Option<VariantContext>, EngineError> {
        let Some(options) = self.reblocking else {
            return self.add_hom_ref_site_base(vc, g);
        };
        let genotype = &vc.genotypes[0];
        // `!genotype.hasGQ() || genotype.getGQ() < rgqThreshold || genotype.getGQ() == 0`. The
        // `int` is compared to the `double` threshold after widening, as Java's `<` does.
        if options.drop_low_quals
            && genotype
                .gq
                .is_none_or(|gq| f64::from(gq) < options.rgq_threshold || gq == 0)
        {
            return Ok(None);
        }
        if is_hom_ref_reblocking(g) {
            if genotype.pl.is_none() {
                if genotype.gq.is_some() {
                    self.warnings.push(format!(
                        "PL is missing for hom ref genotype at at least one position for sample \
                         {}: {}:{}.  Using GQ to determine quality.",
                        genotype.sample_name, vc.contig, vc.start
                    ));
                    return self.add_hom_ref_site_base(vc, genotype);
                }
                let message = format!(
                    "Homozygous reference genotypes must contain GQ or PL. Both are missing for \
                     hom ref genotype at {}:{} for sample {}.",
                    vc.contig, vc.start, genotype.sample_name
                );
                if !options.allow_missing_hom_ref_data {
                    return Err(EngineError::Runtime {
                        class: "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
                            .to_string(),
                        message,
                    });
                }
                self.warnings.push(message);
                let mut filled = genotype.clone();
                filled.gq = Some(0);
                filled.pl = Some(vec![0, 0, 0]);
                let mut rebuilt = vc.clone();
                rebuilt.genotypes = GenotypesContext::new(vec![filled.clone()]);
                return self.add_hom_ref_site_base(&rebuilt, &filled);
            }
            return self.add_hom_ref_site_base(vc, genotype);
        }
        if !genotype.is_called()
            && genotype
                .pl
                .as_ref()
                .is_some_and(|pl| pl.first() == Some(&0))
        {
            return self.add_hom_ref_site_base(vc, genotype);
        }
        Ok(None)
    }

    /// `emitCurrentBlock`: write the open block, if any, and leave none open.
    fn emit_current_block(&mut self) {
        if let Some(block) = self.current_block.take() {
            let sample = self.sample_name.clone().unwrap_or_default();
            self.to_output
                .push(block.to_variant_context(&sample, self.floor_blocks));
        }
    }

    /// `GVCFBlockCombiner.submit`: one single-sample record in, zero or more records queued out.
    ///
    /// Which records are blocked:
    ///
    /// ```java
    /// if ((g.isHomRef()
    ///         || (g.isNoCall() && g.hasPL() && g.getPL()[0] == 0))
    ///         && vc.hasAlternateAllele(Allele.NON_REF_ALLELE) && vc.isBiallelic()) {
    ///         //genotypes with PL=0,0,0 get treated as zero confidence hom refs
    /// ```
    ///
    /// Anything else closes the open block, is queued as it is, and sets the marker that drops
    /// the hom-ref records it covers.
    pub fn submit(&mut self, vc: VariantContext) -> Result<(), EngineError> {
        if vc.genotypes.is_empty() {
            return Err(illegal_argument(
                "GVCF assumes that the VariantContext has genotypes".to_string(),
            ));
        }
        if vc.genotypes.len() != 1 {
            return Err(illegal_argument(format!(
                "GVCF assumes that the VariantContext has exactly one genotype but saw {}",
                vc.genotypes.len()
            )));
        }
        if self.sample_name.is_none() {
            self.sample_name = Some(vc.genotypes[0].sample_name.clone());
        }
        // "we've made a non-contiguous step (across interval, onto another chr), so finalize".
        if self
            .current_block
            .as_ref()
            .is_some_and(|block| !block.is_contiguous(&vc))
        {
            self.emit_current_block();
        }
        let g = vc.genotypes[0].clone();
        let hom_ref_like = g.is_hom_ref()
            || (g.is_no_call() && g.pl.as_ref().is_some_and(|pl| pl.first() == Some(&0)));
        if hom_ref_like && vc.alleles.len() == 2 && vc.alleles.iter().any(is_non_ref) {
            if let Some(band) = self.add_hom_ref_site(&vc, &g)? {
                self.to_output.push(band);
            }
        } else {
            self.emit_current_block();
            self.next_available_start = vc.stop;
            self.contig_of_next_available_start = Some(vc.contig.clone());
            self.to_output.push(vc);
        }
        Ok(())
    }

    /// `consumeFinalizedItems`. `std::mem::take` moves the queue out and leaves an empty one,
    /// the Java's copy-then-`clear()`.
    pub fn drain(&mut self) -> Vec<VariantContext> {
        std::mem::take(&mut self.to_output)
    }

    /// `signalEndOfInput`.
    pub fn end_of_input(&mut self) {
        self.emit_current_block();
    }
}

/// `ReblockingGVCFBlockCombiner.isHomRef`: "consider ./. with no PL a GQ0 hom ref".
///
/// ```java
/// //funky logic for DRAGEN GVCFs where call may not match PLs
/// private boolean isHomRef(final Genotype g) {
///     //consider ./. with no PL a GQ0 hom ref, as if for [0,0,0] PLs
///     return (g.hasPL() && g.getPL()[0] == 0)
///             || (!g.hasPL() && (g.isHomRef() || g.isNoCall()));
/// }
/// ```
///
/// So with a `PL` the call is ignored: this predicate calls a `0/1` whose first PL is 0 hom-ref.
pub fn is_hom_ref_reblocking(g: &Genotype) -> bool {
    match &g.pl {
        Some(pl) => pl.first() == Some(&0),
        None => g.is_hom_ref() || g.is_no_call(),
    }
}

/// The reference, one base at a time, as `ReferenceUtils.getRefBaseAtPosition` reads it
/// (`reference.getSubsequenceAt(contig, start, start).getBases()`).
///
/// A `trait` is Java's interface. The caller passes whatever reads its reference; `&mut self`
/// lets an implementation change its own state while it reads.
pub trait ReferenceBase {
    fn base(&mut self, contig: &str, position: i64) -> u8;
}

/// `ReblockingGVCFBlockCombiner.moveBuilderStart`: a new start, with the reference base there as
/// the reference allele of the record and of every genotype.
///
/// ```java
/// final byte[] newRef = ReferenceUtils.getRefBaseAtPosition(referenceReader, builder.getContig(),
///         newStart);
/// final Allele newRefAllele = Allele.create(newRef, true);
/// ```
///
/// Only the start moves; the stop and the `END` attribute are the caller's to set.
pub fn move_start(vc: &mut VariantContext, new_start: i64, reference: &mut dyn ReferenceBase) {
    let base = reference.base(&vc.contig, new_start);
    let new_reference = Allele::create(&[base], true).expect("a reference base");
    // `iter_mut` walks the list with a writable reference to each element, so the assignment
    // below replaces the allele in place (the Java builds new lists and new genotypes).
    for genotype in vc.genotypes.iter_mut() {
        for allele in genotype.alleles.iter_mut() {
            if allele.is_reference() {
                *allele = new_reference.clone();
            }
        }
    }
    for allele in vc.alleles.iter_mut() {
        if allele.is_reference() {
            *allele = new_reference.clone();
        }
    }
    vc.start = new_start;
}

/// `ReblockingGVCFWriter` over a `ReblockingGVCFBlockCombiner`.
///
/// The Java's combiner is one subclass holding both the reference buffer and the blocking; here
/// the buffer half is this type and the blocking half is the [`BlockCombiner`] it submits to,
/// whose output is collected in `written` instead of going to a `VariantContextWriter`.
pub struct ReblockingWriter {
    pub combiner: BlockCombiner,
    /// `homRefBlockBuffer`: reference records held back, kept sorted by start.
    buffer: Vec<VariantContext>,
    /// `vcfOutputEnd`: the last position handed to the combiner on this contig, 0 for none.
    vcf_output_end: i64,
    /// `bufferEnd`: the furthest stop of a record added to the buffer since the last flush.
    buffer_end: i64,
    current_contig: Option<String>,
    /// What reached the underlying writer, in order.
    pub written: Vec<VariantContext>,
}

impl ReblockingWriter {
    pub fn new(
        partitions: Vec<(i32, i32)>,
        floor_blocks: bool,
        options: ReblockingOptions,
    ) -> Self {
        ReblockingWriter {
            combiner: BlockCombiner::new(partitions, floor_blocks, Some(options)),
            buffer: Vec::new(),
            vcf_output_end: 0,
            buffer_end: 0,
            current_contig: None,
            written: Vec::new(),
        }
    }

    /// `getVcfOutputEnd()`: `None` before anything on this contig was written.
    ///
    /// ```java
    /// if (contig == null) {
    ///     return null;
    /// }
    /// final int position = ((ReblockingGVCFBlockCombiner)gvcfBlockCombiner).getVcfOutputEnd();
    /// if (position == 0) {
    ///     return null;
    /// }
    /// ```
    pub fn vcf_output_end(&self) -> Option<(String, i64)> {
        let contig = self.current_contig.clone()?;
        (self.vcf_output_end != 0).then_some((contig, self.vcf_output_end))
    }

    /// `siteOverlapsBuffer(vc)`: the record starts between the buffer's first start and its
    /// furthest stop. The lower test is there, per the Java's comment, "in case there's a tail
    /// block from a trimmed deletion, but the deletion still overlaps".
    pub fn site_overlaps_buffer(&self, vc: &VariantContext) -> bool {
        let Some(first) = self.buffer.first() else {
            return false;
        };
        self.current_contig.as_deref() == Some(vc.contig.as_str())
            && vc.start <= self.buffer_end
            && vc.start >= first.start
    }

    /// `super.submit(vc)` followed by the writer's `output()`: block it, and collect whatever the
    /// combiner finished.
    fn submit_to_combiner(&mut self, vc: VariantContext) -> Result<(), EngineError> {
        self.combiner.submit(vc)?;
        let drained = self.combiner.drain();
        self.written.extend(drained);
        Ok(())
    }

    fn sort_buffer(&mut self) {
        // `Comparator.comparingLong(getStart)` over a stable `List.sort`. `sort_by_key` is
        // stable too, so records with one start keep their order.
        self.buffer.sort_by_key(|builder| builder.start);
    }

    /// `ReblockingGVCFWriter.add`, which is `ReblockingGVCFBlockCombiner.submit`.
    ///
    /// Every record, reference or not, cuts the buffered blocks it overlaps (`trimBlockToVariant`,
    /// over the record's own span, not the moved start below). A reference record (by
    /// [`is_hom_ref_reblocking`]) then goes into the buffer, moved past `vcf_output_end` first,
    /// or is dropped if it ends at or before it. Only a variant hands the blocks that start
    /// before it to the combiner ("only flush ref blocks if we're outputting a variant"), and is
    /// then handed over itself: at once when a buffered block starts after its end, otherwise
    /// after the walk when the buffer is empty or `bufferEnd` reaches its end.
    ///
    /// The negative-length message names `currentContig` BEFORE this record updates it, so the
    /// first record of a run prints `null`, as the Java's does:
    ///
    /// ```java
    /// Utils.validate(variantContextToOutput.getStart() <= variantContextToOutput.getEnd(),
    ///         () -> "Input variant context at position " + currentContig + ":"
    ///                 + variantContextToOutput.getStart() + " has negative length: start="
    ///                 + variantContextToOutput.getStart() + " end="
    ///                 + variantContextToOutput.getEnd());
    /// ```
    pub fn add(
        &mut self,
        vc: VariantContext,
        reference: &mut dyn ReferenceBase,
    ) -> Result<(), EngineError> {
        if vc.start > vc.stop {
            return Err(EngineError::Runtime {
                class: "java.lang.IllegalStateException".to_string(),
                message: format!(
                    "Input variant context at position {}:{} has negative length: start={} end={}",
                    self.current_contig
                        .clone()
                        .unwrap_or_else(|| "null".to_string()),
                    vc.start,
                    vc.start,
                    vc.stop
                ),
            });
        }
        match &self.current_contig {
            None => self.current_contig = Some(vc.contig.clone()),
            Some(contig) if *contig != vc.contig => {
                self.flush_buffer()?;
                self.current_contig = Some(vc.contig.clone());
                self.vcf_output_end = 0;
            }
            Some(_) => {}
        }
        // `newHomRefBlock`, a copy of the record that only a reference record goes on to use.
        let mut new_block = vc.clone();
        let g = vc.genotypes[0].clone();
        let hom_ref = is_hom_ref_reblocking(&g);
        if hom_ref && vc.start <= self.vcf_output_end {
            if vc.stop <= self.vcf_output_end {
                return Ok(());
            }
            move_start(&mut new_block, self.vcf_output_end + 1, reference);
        }

        // `completedBlocks` and `tailBuffer`. The Java collects the builders themselves and
        // removes them with `removeAll`; here the completed blocks are their positions in the
        // buffer, and the loop is over an index because `trim_block` changes the buffer's
        // records while it runs.
        let mut completed: Vec<usize> = Vec::new();
        let mut tails: Vec<VariantContext> = Vec::new();
        let variant_start = vc.start;
        let variant_end = vc.stop;
        let mut index = 0;
        while index < self.buffer.len() {
            let block_start = self.buffer[index].start;
            // The first block wholly after the record ends the walk. The Java returns from inside
            // its for-each here; a reference record is added WITHOUT raising `bufferEnd`, which
            // only the path after the loop does.
            if block_start > variant_end {
                if !hom_ref {
                    self.submit_to_combiner(vc.clone())?;
                    self.vcf_output_end = self.vcf_output_end.max(variant_end);
                } else {
                    self.buffer.push(new_block);
                }
                self.remove_completed(&completed);
                self.buffer.extend(tails);
                self.sort_buffer();
                return Ok(());
            }
            // "then trim out overlap".
            let mut block_end = self.buffer[index].stop;
            if block_end >= variant_start {
                block_end = self.trim_block(&vc, index, &mut completed, &mut tails, reference)?;
            }
            // "only flush ref blocks if we're outputting a variant"
            if block_start < variant_start && !hom_ref {
                let block = self.buffer[index].clone();
                self.submit_to_combiner(block)?;
                self.vcf_output_end = block_end;
                completed.push(index);
            }
            index += 1;
        }
        self.remove_completed(&completed);
        if hom_ref {
            // `ReblockGVCF.isHomRefBlock(vc)`: one ALT and it is `<NON_REF>`.
            let is_block = vc.alleles.len() == 2 && is_non_ref(&vc.alleles[1]);
            if is_block && new_block.start < self.vcf_output_end {
                return Err(EngineError::Runtime {
                    class: "java.lang.IllegalStateException".to_string(),
                    message: format!(
                        "Reference positions added to buffer should not overlap positions already \
                         output to VCF. {} overlaps position {}:{} already emitted.",
                        vc.start,
                        self.current_contig.clone().unwrap_or_default(),
                        self.vcf_output_end
                    ),
                });
            }
            self.buffer_end = self.buffer_end.max(new_block.stop);
            self.buffer.push(new_block);
        } else if self.buffer_end >= vc.stop || self.buffer.is_empty() {
            // The Java's `else if` has no `else`: a variant meeting neither test is not handed
            // over on this path.
            self.submit_to_combiner(vc.clone())?;
            self.vcf_output_end = self.vcf_output_end.max(vc.stop);
        }
        self.buffer.extend(tails);
        self.sort_buffer();
        Ok(())
    }

    /// `homRefBlockBuffer.removeAll(completedBlocks)`, by position in the buffer.
    ///
    /// `retain` keeps the elements its closure answers `true` for, in order; the closure counts
    /// the positions itself.
    fn remove_completed(&mut self, completed: &[usize]) {
        let mut index = 0;
        self.buffer.retain(|_| {
            let keep = !completed.contains(&index);
            index += 1;
            keep
        });
    }

    /// `trimBlockToVariant`: remove the overlap between the buffered block at `index` and the
    /// record `vc`, returning the block's new end.
    ///
    /// Three shapes:
    ///
    /// - the block starts before the record and ends after it: SPLIT. A copy moved to start the
    ///   base after the record (with the reference base there) goes to `tails`, and the block
    ///   itself is cut to end the base before the record;
    /// - the block starts before the record: RIGHT TRIM, to end the base before it;
    /// - otherwise: dropped when the record covers it, else LEFT TRIM, moved to start the base
    ///   after the record.
    ///
    /// The right trim sets `END` and the stop together, as `builder.attribute(END_KEY, ..)` and
    /// `builder.stop(..)` do, so the written `END` follows.
    fn trim_block(
        &mut self,
        vc: &VariantContext,
        index: usize,
        completed: &mut Vec<usize>,
        tails: &mut Vec<VariantContext>,
        reference: &mut dyn ReferenceBase,
    ) -> Result<i64, EngineError> {
        let block_start = self.buffer[index].start;
        let variant_end = vc.stop;
        let mut block_end = self.buffer[index].stop;
        let variant_start = vc.start;
        if block_end > variant_end && block_start < variant_start {
            let mut tail = self.buffer[index].clone();
            move_start(&mut tail, variant_end + 1, reference);
            tails.push(tail);
            let builder = &mut self.buffer[index];
            builder.stop = variant_start - 1;
            put_attribute(builder, "END", Value::Int(variant_start - 1));
            block_end = variant_start - 1;
        }
        if block_start < variant_start {
            let builder = &mut self.buffer[index];
            put_attribute(builder, "END", Value::Int(variant_start - 1));
            builder.stop = variant_start - 1;
        } else if vc.contig == self.buffer[index].contig
            && vc.start <= block_start
            && block_end <= vc.stop
        {
            // `variantContextToOutput.contains(new SimpleInterval(currentContig, blockStart,
            // blockEnd))`: a buffered block is always on the current contig.
            // "if block is entirely overlapped by VC to output, then remove it from the buffer"
            completed.push(index);
        } else {
            if block_end < variant_end + 1 {
                return Err(should_never_reach_here(format!(
                    "ref block end overlaps variant end; current builder: {} to {}",
                    self.buffer[index].start, self.buffer[index].stop
                )));
            }
            move_start(&mut self.buffer[index], variant_end + 1, reference);
        }
        let builder = &self.buffer[index];
        if builder.start > builder.stop {
            return Err(should_never_reach_here(format!(
                "builder start follows stop; current builder: {} to {}",
                builder.start, builder.stop
            )));
        }
        Ok(block_end)
    }

    /// `flushRefBlockBuffer`: every buffered block to the combiner, in order.
    ///
    /// ```java
    /// try {
    ///     super.submit(builder.make());
    /// } catch (Exception e) {
    ///     throw new GATKException("builder threw an exception at " + builder.getContig() + ":"
    ///             + builder.getStart(), e);
    /// }
    /// vcfOutputEnd = (int)builder.getStop();
    /// ```
    ///
    /// The Java wraps any exception and keeps it as the cause; this puts the cause's class and
    /// message in the wrapped message's parentheses. A limitation of the port is not a reference
    /// exception, so it is passed through unwrapped.
    fn flush_buffer(&mut self) -> Result<(), EngineError> {
        // Take the whole buffer out (leaving it empty, the Java's `clear()` after the loop), so
        // the loop owns each record and can hand it over.
        let buffer = std::mem::take(&mut self.buffer);
        for builder in buffer {
            let (contig, start, stop) = (builder.contig.clone(), builder.start, builder.stop);
            // `map_err` rewrites the error of a failed call and leaves a success alone.
            self.submit_to_combiner(builder).map_err(|error| {
                let cause = match &error {
                    EngineError::Runtime { class, message } => format!("{class}: {message}"),
                    EngineError::Limitation(what) => what.clone(),
                };
                match error {
                    EngineError::Limitation(_) => error,
                    EngineError::Runtime { .. } => gatk_exception(format!(
                        "builder threw an exception at {contig}:{start} ({cause})"
                    )),
                }
            })?;
            self.vcf_output_end = stop;
        }
        self.buffer_end = 0;
        Ok(())
    }

    /// `GVCFWriter.close`: `signalEndOfInput`, which for the reblocking combiner flushes the
    /// buffer and then closes the open block, then `output()`.
    pub fn close(&mut self) -> Result<(), EngineError> {
        self.flush_buffer()?;
        self.combiner.end_of_input();
        let drained = self.combiner.drain();
        self.written.extend(drained);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partitions_close_below_and_open_above_with_a_last_band_to_one_hundred() {
        assert_eq!(
            parse_partitions(&[20, 100]).unwrap(),
            vec![(0, 20), (20, 100)]
        );
        assert_eq!(parse_partitions(&[60]).unwrap(), vec![(0, 60), (60, 100)]);
        assert!(parse_partitions(&[20, 20]).is_err());
        assert!(parse_partitions(&[101]).is_err());
    }

    #[test]
    fn the_gq_is_the_gap_between_the_two_smallest_pls() {
        assert_eq!(calculate_gq_from_pls(&[0, 25, 250]).unwrap(), 25);
        assert_eq!(calculate_gq_from_pls(&[40, 0, 10]).unwrap(), 10);
    }

    #[test]
    fn an_even_count_takes_the_mean_of_its_middle_depths() {
        assert_eq!(median(&[20, 22]), 21.0);
        assert_eq!(median(&[20, 30, 22]), 22.0);
    }
}

/// `SomaticGVCFBlockCombiner.calculatePartitionPrecision`: the decimal places the smallest gap
/// between two bands needs, `ceil(-log10(gap))`, after refusing a repeated or out-of-order band.
///
/// The first gap is measured from `Double.NEGATIVE_INFINITY` and is infinite, so a ONE-band list
/// has no finite gap: `-log10` of infinity is minus infinity, its `ceil` too, and the `(int)`
/// cast (like `as i32` here) saturates it to `Integer.MIN_VALUE`. The default bands, -2.5 to 1
/// in steps of 0.5 (`M2ArgumentCollection`), give 1.
///
/// The Java's two checks are `value != lastLOD` and `value > lastLOD`; the second is written
/// `value <= last` here, which is not its negation for NaN: a NaN band is refused by the Java
/// as out of order and passes here. The messages are the Java's `String.format("%f")`, written
/// `{:.6}` here.
fn lod_partition_precision(bands: &[f64]) -> Result<i32, String> {
    let mut smallest = f64::INFINITY;
    let mut last = f64::NEG_INFINITY;
    for &value in bands {
        if value == last {
            return Err(format!(
                "The value {value:.6} appears more than once in the list of TLOD partitions."
            ));
        }
        if value <= last {
            return Err(format!(
                "The list of TLOD partitions is out of order. Previous value is {last:.6} but the \
                 next is {value:.6}."
            ));
        }
        let delta = value - last;
        if delta < smallest {
            smallest = delta;
        }
        last = value;
    }
    Ok((-jmath::math::log10(smallest)).ceil() as i32)
}

/// `Math.pow(10, p)` for the integral exponents a band precision takes: a product of tens,
/// exact up to `10^22`, and for a negative exponent one correctly rounded division, which is the
/// value `pow` returns there.
///
/// The exponent is small for any list of two or more bands, but a one-band list gives
/// `i32::MIN` (see [`lod_partition_precision`]). The loop then runs 2^31 times: the product
/// overflows to infinity after 309 of them and stays there, so the result is `1 / inf`, 0, which
/// is also what `Math.pow(10, Integer.MIN_VALUE)` underflows to. Every LOD then bins to 0.
fn ten_to(precision: i32) -> f64 {
    let mut power = 1.0f64;
    for _ in 0..precision.unsigned_abs() {
        power *= 10.0;
    }
    if precision < 0 {
        1.0 / power
    } else {
        power
    }
}

/// `convertLODtoInt`: the LOD binned to the partitions' precision.
///
/// ```java
/// return (int)Math.floor(LOD * Math.pow(10, partitionPrecision));
/// ```
///
/// The cast saturates: a LOD whose binned value reaches `Integer.MAX_VALUE` is in no band,
/// because the last band is open at `Integer.MAX_VALUE`, and `add_hom_ref_site` refuses it.
fn lod_to_int(lod: f64, precision: i32) -> i32 {
    (lod * ten_to(precision)).floor() as i32
}

/// `SomaticGVCFBlockCombiner.parsePartitions`: the bands as binned integers, each closed below and
/// open above, from `Integer.MIN_VALUE` up to the first and from the last up to
/// `Integer.MAX_VALUE`, with the precision they were binned at.
///
/// The precision is worked out before the empty check, as the Java's `parsePartitions` sets
/// `partitionPrecision = calculatePartitionPrecision(gqPartitions)` before `Utils.nonEmpty`,
/// so an out-of-order list is refused for its order whatever else is wrong with it.
pub fn parse_lod_partitions(bands: &[f64]) -> Result<(Vec<(i32, i32)>, i32), String> {
    let precision = lod_partition_precision(bands)?;
    if bands.is_empty() {
        return Err("The collection is empty: collection must not be null or empty.".to_string());
    }
    let mut result = Vec::new();
    let mut last = i32::MIN;
    for &value in bands {
        let threshold = lod_to_int(value, precision);
        result.push((last, threshold));
        last = threshold;
    }
    result.push((last, i32::MAX));
    Ok((result, precision))
}

/// The `TLOD` a reference-confidence genotype carries, `Double.parseDouble(attribute.toString())`.
///
/// Absent, the Java's `.toString()` on `null` is a `NullPointerException`, which is the class
/// reported here. An unparsable string throws `NumberFormatException` in the Java, a subclass of
/// `IllegalArgumentException`; this reports the parent class with the same message.
fn genotype_tlod(g: &Genotype) -> Result<f64, EngineError> {
    match g.extended.iter().find(|(key, _)| key == "TLOD") {
        Some((_, Value::Double(lod))) => Ok(*lod),
        Some((_, Value::Str(text))) => text
            .parse::<f64>()
            .map_err(|_| illegal_argument(format!("For input string: \"{text}\""))),
        _ => Err(EngineError::Runtime {
            class: "java.lang.NullPointerException".to_string(),
            message: "a genotype without TLOD in a somatic GVCF".to_string(),
        }),
    }
}

/// `TLODBlock`: a run of reference sites whose binned LOD falls in one band, written with the
/// smallest LOD it saw and the median and minimum of its depths.
///
/// The Java stores the band as binned integers "to facilitate GVCF class polymorphism and to
/// avoid double precision issues", reusing `GVCFBlock`'s GQ bounds for them.
#[derive(Debug, Clone)]
struct TlodBlock {
    starting: VariantContext,
    /// The band, `[lower, upper)`, in binned units.
    lower: i32,
    upper: i32,
    precision: i32,
    reference: Allele,
    dps: Vec<i32>,
    end: i64,
    /// `minBlockLOD`: `Double.POSITIVE_INFINITY` until a genotype with a `TLOD` is seen.
    min_block_lod: f64,
}

impl TlodBlock {
    fn new(
        starting: &VariantContext,
        lower: i32,
        upper: i32,
        precision: i32,
    ) -> Result<TlodBlock, EngineError> {
        let g = &starting.genotypes[0];
        let mut block = TlodBlock {
            starting: starting.clone(),
            lower,
            upper,
            precision,
            reference: starting.alleles[0].clone(),
            dps: g.dp.map(|dp| dp.max(0)).into_iter().collect(),
            end: end_attribute(starting).unwrap_or(starting.start),
            min_block_lod: f64::INFINITY,
        };
        // `if (g.hasExtendedAttribute(TUMOR_LOG_10_ODDS_KEY)) setMinBlockLod(g)`.
        if g.extended.iter().any(|(key, _)| key == "TLOD") {
            block.set_min_block_lod(g)?;
        }
        Ok(block)
    }

    fn within_bounds(&self, binned: i32) -> bool {
        binned >= self.lower && binned < self.upper
    }

    /// `setMinBlockLod`: refuse a LOD outside the band, then keep the smaller.
    ///
    /// The refusal is unreachable from this file's callers, which pick the band by the same
    /// binned LOD. Its message is not the Java's, which concatenates doubles with
    /// `Double.toString` (`LOD=1.0E-5`) where `{}` writes `0.00001`.
    fn set_min_block_lod(&mut self, g: &Genotype) -> Result<(), EngineError> {
        let lod = genotype_tlod(g)?;
        if !self.within_bounds(lod_to_int(lod, self.precision)) {
            return Err(illegal_argument(format!(
                "cannot add a genotype with LOD={lod} because it's not within bounds [{},{})",
                f64::from(self.lower) / ten_to(self.precision),
                f64::from(self.upper) / ten_to(self.precision),
            )));
        }
        if self.min_block_lod == f64::INFINITY || lod < self.min_block_lod {
            self.min_block_lod = lod;
        }
        Ok(())
    }

    /// `TLODBlock.add`: unlike `HomRefBlock.add` it tests neither ploidy nor `PL`, and one
    /// message covers both a gap and an overlap
    /// (`Utils.validateArg(pos == end + 1 , "adding genotype at pos " + ..)`).
    fn add(&mut self, position: i64, new_end: i64, g: &Genotype) -> Result<(), EngineError> {
        if position != self.end + 1 {
            return Err(illegal_argument(format!(
                "adding genotype at pos {position} isn't contiguous with previous end {}",
                self.end
            )));
        }
        self.set_min_block_lod(g)?;
        self.end = new_end;
        if let Some(dp) = g.dp {
            self.dps.push(dp.max(0));
        }
        Ok(())
    }

    /// `toVariantContext` with `createHomRefGenotype`: diploid whatever the sample, its only
    /// attribute the block's smallest TLOD, then `DP` and `MIN_DP` when any depth was seen.
    ///
    /// ```java
    /// final GenotypeBuilder gb = new GenotypeBuilder(sampleName,
    ///         Collections.nCopies(2, getRef()));
    /// gb.noAD().noPL().noAttributes(); // clear all attributes
    ///
    /// gb.attribute(GATKVCFConstants.TUMOR_LOG_10_ODDS_KEY, minBlockLOD);
    /// if (!DPs.isEmpty()) {
    ///     gb.DP(getMedianDP());
    ///     gb.attribute(GATKVCFConstants.MIN_DP_FORMAT_KEY, getMinDP());
    /// }
    /// ```
    ///
    /// The first line carries the Java's note "for somatic stuff we output the genotype as
    /// diploid because that's familiar for human". No `GQ` and no `PL`, so a somatic block reads
    /// `GT:DP:MIN_DP:TLOD`. The `TLOD` is a
    /// `Double` and is formatted by the VCF encoder, not here.
    fn to_variant_context(&self, sample_name: &str) -> VariantContext {
        let mut vc = self.starting.clone();
        vc.attributes = vec![("END".to_string(), Value::Int(self.end))];
        vc.stop = self.end;
        let mut genotype = Genotype::new(sample_name, vec![self.reference.clone(); 2]);
        genotype
            .extended
            .push(("TLOD".to_string(), Value::Double(self.min_block_lod)));
        if !self.dps.is_empty() {
            genotype.dp =
                Some(crate::reference_confidence_merger::java_round(median(&self.dps)) as i32);
            let min = *self.dps.iter().min().expect("a depth");
            genotype
                .extended
                .push((MIN_DP_FORMAT_KEY.to_string(), Value::Int(i64::from(min))));
        }
        vc.genotypes = GenotypesContext::new(vec![genotype]);
        vc
    }

    /// `GVCFBlock.isContiguous`, as for [`HomRefBlock`].
    fn is_contiguous(&self, vc: &VariantContext) -> bool {
        vc.contig == self.starting.contig
            && vc.start <= self.end + 1
            && self.starting.start - 1 <= vc.stop
    }
}

/// `SomaticGVCFBlockCombiner`: `GVCFBlockCombiner` with the blocks banded by `TLOD` instead of GQ.
///
/// The Java subclass overrides `genotypeCanBeMergedInCurrentBlock`, `createNewBlock` and
/// `parsePartitions` and inherits the rest; Rust has no inheritance, so the inherited parts
/// (`submit`, the variant marker, `emitCurrentBlock`) are repeated here over `TlodBlock`.
/// `Mutect2` builds it through `new SomaticGVCFWriter(vcfWriter, ..)`, the constructor that does
/// not override the computed precision.
pub struct SomaticBlockCombiner {
    partitions: Vec<(i32, i32)>,
    precision: i32,
    to_output: Vec<VariantContext>,
    next_available_start: i64,
    contig_of_next_available_start: Option<String>,
    sample_name: Option<String>,
    current_block: Option<TlodBlock>,
}

impl SomaticBlockCombiner {
    pub fn new(partitions: Vec<(i32, i32)>, precision: i32) -> SomaticBlockCombiner {
        SomaticBlockCombiner {
            partitions,
            precision,
            to_output: Vec::new(),
            next_available_start: -1,
            contig_of_next_available_start: None,
            sample_name: None,
            current_block: None,
        }
    }

    fn emit_current_block(&mut self) {
        if let Some(block) = self.current_block.take() {
            let sample = self.sample_name.clone().unwrap_or_default();
            self.to_output.push(block.to_variant_context(&sample));
        }
    }

    /// `addHomRefSite`, with `genotypeCanBeMergedInCurrentBlock` and `createNewBlock` reading the
    /// genotype's `TLOD`.
    ///
    /// The `TLOD` is parsed before the block test (`genotypeCanBeMergedInCurrentBlock` parses it
    /// before its `currentTLODBlock != null`), so every record that reaches here must carry one,
    /// the first included. The refusal's message is not the Java's, which concatenates the whole
    /// genotype and record (`"LOD for genotype " + g + " from " + vc + ..`).
    fn add_hom_ref_site(
        &mut self,
        vc: &VariantContext,
        g: &Genotype,
    ) -> Result<Option<VariantContext>, EngineError> {
        if self.next_available_start != -1 {
            if vc.start <= self.next_available_start
                && self.contig_of_next_available_start.as_deref() == Some(vc.contig.as_str())
                && vc.stop <= self.next_available_start
            {
                return Ok(None);
            }
            self.next_available_start = -1;
            self.contig_of_next_available_start = None;
        }
        let binned = lod_to_int(genotype_tlod(g)?, self.precision);
        if self
            .current_block
            .as_ref()
            .is_some_and(|block| block.within_bounds(binned))
        {
            let end = end_attribute(vc).unwrap_or(vc.start);
            self.current_block
                .as_mut()
                .expect("a block to merge into")
                .add(vc.start, end, g)?;
            return Ok(None);
        }
        let sample = self.sample_name.clone().unwrap_or_default();
        let result = self
            .current_block
            .as_ref()
            .map(|block| block.to_variant_context(&sample));
        let Some((lower, upper)) = self
            .partitions
            .iter()
            .copied()
            .find(|(lower, upper)| *lower <= binned && binned < *upper)
        else {
            return Err(gatk_exception(format!(
                "LOD for genotype {} from {}:{} didn't fit into any partition",
                g.sample_name, vc.contig, vc.start
            )));
        };
        self.current_block = Some(TlodBlock::new(vc, lower, upper, self.precision)?);
        Ok(result)
    }

    /// `GVCFBlockCombiner.submit`, inherited unchanged: the same as [`BlockCombiner::submit`].
    pub fn submit(&mut self, vc: VariantContext) -> Result<(), EngineError> {
        if vc.genotypes.is_empty() {
            return Err(illegal_argument(
                "GVCF assumes that the VariantContext has genotypes".to_string(),
            ));
        }
        if vc.genotypes.len() != 1 {
            return Err(illegal_argument(format!(
                "GVCF assumes that the VariantContext has exactly one genotype but saw {}",
                vc.genotypes.len()
            )));
        }
        if self.sample_name.is_none() {
            self.sample_name = Some(vc.genotypes[0].sample_name.clone());
        }
        if self
            .current_block
            .as_ref()
            .is_some_and(|block| !block.is_contiguous(&vc))
        {
            self.emit_current_block();
        }
        let g = vc.genotypes[0].clone();
        let hom_ref_like = g.is_hom_ref()
            || (g.is_no_call() && g.pl.as_ref().is_some_and(|pl| pl.first() == Some(&0)));
        if hom_ref_like && vc.alleles.len() == 2 && vc.alleles.iter().any(is_non_ref) {
            if let Some(band) = self.add_hom_ref_site(&vc, &g)? {
                self.to_output.push(band);
            }
        } else {
            self.emit_current_block();
            self.next_available_start = vc.stop;
            self.contig_of_next_available_start = Some(vc.contig.clone());
            self.to_output.push(vc);
        }
        Ok(())
    }

    /// `consumeFinalizedItems`.
    pub fn drain(&mut self) -> Vec<VariantContext> {
        std::mem::take(&mut self.to_output)
    }

    /// `signalEndOfInput`.
    pub fn end_of_input(&mut self) {
        self.emit_current_block();
    }
}
