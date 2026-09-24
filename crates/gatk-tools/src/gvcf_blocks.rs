//! `GVCFWriter`'s block combining, and the reblocking variant `ReblockGVCF` writes through.
//!
//! Ported from `org.broadinstitute.hellbender.utils.variant.writers` (GATK 4.6.2.0):
//! `GVCFBlockCombiner`, `HomRefBlock`, `GVCFBlock`, `ReblockingGVCFBlockCombiner` and the two
//! writers around them. The combiner takes records one at a time and hands back the records to
//! write, in order: a hom-ref record joins the open block when its GQ falls in the block's band,
//! and anything else closes the block and is written as it is.
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
//! A genotype carrying `PP` switches `HomRefBlock` to phred-scaled posteriors for its GQ. No input
//! this port reads carries them, so a genotype with `PP` is refused as a limitation.

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

fn is_non_ref(allele: &Allele) -> bool {
    allele.is_symbolic() && allele.display_string() == "<NON_REF>"
}

/// `vc.getAttributeAsInt(END, default)`.
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
pub fn put_attribute(vc: &mut VariantContext, key: &str, value: Value) {
    match vc.attributes.iter_mut().find(|(name, _)| name == key) {
        Some((_, slot)) => *slot = value,
        None => vc.attributes.push((key.to_string(), value)),
    }
}

/// `GATKVariantContextUtils.calculateGQFromPLs`: the second smallest PL less the smallest.
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
fn median(values: &[i32]) -> f64 {
    let mut sorted: Vec<f64> = values.iter().map(|value| f64::from(*value)).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite depths"));
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
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
pub fn add_ranges_to_header(header: &mut VcfHeader, partitions: &[(i32, i32)]) {
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

/// `HomRefBlock`.
#[derive(Debug, Clone)]
struct HomRefBlock {
    starting: VariantContext,
    min_gq_bound: i32,
    max_gq_bound: i32,
    reference: Allele,
    dps: Vec<i32>,
    end: i64,
    ploidy: usize,
    min_pls: Option<Vec<i32>>,
    /// `-1` when no genotype has had a PL or a GQ.
    min_gq: i32,
}

fn refuse_posteriors(genotype: &Genotype) -> Result<(), EngineError> {
    if genotype
        .extended
        .iter()
        .any(|(key, _)| key == PHRED_SCALED_POSTERIORS_KEY)
    {
        return Err(EngineError::Limitation(
            "a reference block whose genotype carries PP (phred-scaled posteriors) is not ported."
                .to_string(),
        ));
    }
    Ok(())
}

impl HomRefBlock {
    fn new(
        starting: &VariantContext,
        lower: i32,
        upper: i32,
        default_ploidy: usize,
    ) -> Result<HomRefBlock, EngineError> {
        let g = &starting.genotypes[0];
        refuse_posteriors(g)?;
        // `startingVC.getMaxPloidy(defaultPloidy)`.
        let max = starting
            .genotypes
            .iter()
            .map(Genotype::ploidy)
            .max()
            .unwrap_or(0);
        let ploidy = if max == 0 { default_ploidy } else { max };
        let min_pls = g.pl.clone();
        let min_gq = match (&min_pls, g.gq) {
            (Some(pls), _) => calculate_gq_from_pls(pls)?,
            (None, Some(gq)) => gq,
            (None, None) => -1,
        };
        Ok(HomRefBlock {
            starting: starting.clone(),
            min_gq_bound: lower,
            max_gq_bound: upper,
            reference: starting.alleles[0].clone(),
            dps: g.dp.map(|dp| dp.max(0)).into_iter().collect(),
            end: end_attribute(starting).unwrap_or(starting.start),
            ploidy,
            min_pls,
            min_gq,
        })
    }

    fn within_bounds(&self, gq: i32) -> bool {
        gq >= self.min_gq_bound && gq < self.max_gq_bound
    }

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
        let gq = genotype.gq.unwrap_or(-1);
        if !self.within_bounds(gq.min(MAX_GENOTYPE_QUAL)) {
            return Err(illegal_argument(format!(
                "cannot add a genotype with GQ={gq} because it's not within bounds [{},{})",
                self.min_gq_bound, self.max_gq_bound
            )));
        }
        refuse_posteriors(genotype)?;
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
        self.min_gq = match &self.min_pls {
            Some(pls) => calculate_gq_from_pls(pls)?,
            None if self.min_gq == -1 => gq,
            None => self.min_gq.min(gq),
        };
        self.end = new_end;
        if let Some(dp) = genotype.dp {
            self.dps.push(dp.max(0));
        }
        Ok(())
    }

    /// `GVCFBlock.toVariantContext(sampleName, floorBlocks)`.
    fn to_variant_context(&self, sample_name: &str, floor_blocks: bool) -> VariantContext {
        let mut vc = self.starting.clone();
        vc.attributes = vec![("END".to_string(), Value::Int(self.end))];
        vc.stop = self.end;
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
        vc.genotypes = GenotypesContext::new(vec![genotype]);
        vc
    }

    /// `isContiguous(vc)`: `vc.withinDistanceOf(block, 1)`.
    fn is_contiguous(&self, vc: &VariantContext) -> bool {
        vc.contig == self.starting.contig
            && vc.start <= self.end + 1
            && self.starting.start - 1 <= vc.stop
    }
}

/// `ReblockingOptions`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReblockingOptions {
    pub drop_low_quals: bool,
    pub allow_missing_hom_ref_data: bool,
    pub rgq_threshold: f64,
}

/// `GVCFBlockCombiner`, with `ReblockingGVCFBlockCombiner.addHomRefSite` when `reblocking` is set.
pub struct BlockCombiner {
    partitions: Vec<(i32, i32)>,
    floor_blocks: bool,
    to_output: Vec<VariantContext>,
    next_available_start: i64,
    contig_of_next_available_start: Option<String>,
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

    fn can_merge(&self, g: &Genotype) -> bool {
        self.current_block.as_ref().is_some_and(|block| {
            block.within_bounds(g.gq.unwrap_or(-1).min(MAX_GENOTYPE_QUAL))
                && block.ploidy == g.ploidy()
                && (block.min_pls.is_none()
                    || g.pl.is_none()
                    || block.min_pls.as_ref().map(Vec::len) == g.pl.as_ref().map(Vec::len))
        })
    }

    fn create_new_block(
        &self,
        vc: &VariantContext,
        g: &Genotype,
    ) -> Result<HomRefBlock, EngineError> {
        let gq = g.gq.map_or(0, |gq| gq.min(MAX_GENOTYPE_QUAL));
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

    /// `GVCFBlockCombiner.addHomRefSite`.
    fn add_hom_ref_site_base(
        &mut self,
        vc: &VariantContext,
        g: &Genotype,
    ) -> Result<Option<VariantContext>, EngineError> {
        if self.next_available_start != -1 {
            // "overlapping deletions on different haplotypes"
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
    fn add_hom_ref_site(
        &mut self,
        vc: &VariantContext,
        g: &Genotype,
    ) -> Result<Option<VariantContext>, EngineError> {
        let Some(options) = self.reblocking else {
            return self.add_hom_ref_site_base(vc, g);
        };
        let genotype = &vc.genotypes[0];
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

    fn emit_current_block(&mut self) {
        if let Some(block) = self.current_block.take() {
            let sample = self.sample_name.clone().unwrap_or_default();
            self.to_output
                .push(block.to_variant_context(&sample, self.floor_blocks));
        }
    }

    /// `GVCFBlockCombiner.submit`.
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

/// `ReblockingGVCFBlockCombiner.isHomRef`: "consider ./. with no PL a GQ0 hom ref".
pub fn is_hom_ref_reblocking(g: &Genotype) -> bool {
    match &g.pl {
        Some(pl) => pl.first() == Some(&0),
        None => g.is_hom_ref() || g.is_no_call(),
    }
}

/// The reference, one base at a time, as `ReferenceUtils.getRefBaseAtPosition` reads it.
pub trait ReferenceBase {
    fn base(&mut self, contig: &str, position: i64) -> u8;
}

/// `ReblockingGVCFBlockCombiner.moveBuilderStart`: a new start, with the reference base there as
/// the reference allele of the record and of every genotype.
pub fn move_start(vc: &mut VariantContext, new_start: i64, reference: &mut dyn ReferenceBase) {
    let base = reference.base(&vc.contig, new_start);
    let new_reference = Allele::create(&[base], true).expect("a reference base");
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
pub struct ReblockingWriter {
    pub combiner: BlockCombiner,
    buffer: Vec<VariantContext>,
    vcf_output_end: i64,
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
    pub fn vcf_output_end(&self) -> Option<(String, i64)> {
        let contig = self.current_contig.clone()?;
        (self.vcf_output_end != 0).then_some((contig, self.vcf_output_end))
    }

    /// `siteOverlapsBuffer(vc)`.
    pub fn site_overlaps_buffer(&self, vc: &VariantContext) -> bool {
        let Some(first) = self.buffer.first() else {
            return false;
        };
        self.current_contig.as_deref() == Some(vc.contig.as_str())
            && vc.start <= self.buffer_end
            && vc.start >= first.start
    }

    fn submit_to_combiner(&mut self, vc: VariantContext) -> Result<(), EngineError> {
        self.combiner.submit(vc)?;
        let drained = self.combiner.drain();
        self.written.extend(drained);
        Ok(())
    }

    fn sort_buffer(&mut self) {
        // `Comparator.comparingLong(getStart)` over a stable `List.sort`.
        self.buffer.sort_by_key(|builder| builder.start);
    }

    /// `ReblockingGVCFWriter.add`, which is `ReblockingGVCFBlockCombiner.submit`.
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
        let mut new_block = vc.clone();
        let g = vc.genotypes[0].clone();
        let hom_ref = is_hom_ref_reblocking(&g);
        if hom_ref && vc.start <= self.vcf_output_end {
            if vc.stop <= self.vcf_output_end {
                return Ok(());
            }
            move_start(&mut new_block, self.vcf_output_end + 1, reference);
        }

        let mut completed: Vec<usize> = Vec::new();
        let mut tails: Vec<VariantContext> = Vec::new();
        let variant_start = vc.start;
        let variant_end = vc.stop;
        let mut index = 0;
        while index < self.buffer.len() {
            let block_start = self.buffer[index].start;
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
            self.submit_to_combiner(vc.clone())?;
            self.vcf_output_end = self.vcf_output_end.max(vc.stop);
        }
        self.buffer.extend(tails);
        self.sort_buffer();
        Ok(())
    }

    /// `homRefBlockBuffer.removeAll(completedBlocks)`, by position in the buffer.
    fn remove_completed(&mut self, completed: &[usize]) {
        let mut index = 0;
        self.buffer.retain(|_| {
            let keep = !completed.contains(&index);
            index += 1;
            keep
        });
    }

    /// `trimBlockToVariant`.
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

    /// `flushRefBlockBuffer`.
    fn flush_buffer(&mut self) -> Result<(), EngineError> {
        let buffer = std::mem::take(&mut self.buffer);
        for builder in buffer {
            let (contig, start, stop) = (builder.contig.clone(), builder.start, builder.stop);
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

    /// `GVCFWriter.close`: the buffer, then the open block.
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
