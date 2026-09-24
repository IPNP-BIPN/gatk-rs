//! `CombineGVCFs`: how several single-sample GVCFs become one multi-sample one.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.CombineGVCFs` and the
//! `MultiVariantWalkerGroupedOnStart` it extends (GATK 4.6.2.0). The merge of the records that meet
//! at a variant is [`crate::reference_confidence_merger`]; what is here is the walk that decides
//! where the output's records start and stop, and the cheap merge it uses where only reference
//! blocks meet.
//!
//! # Every sample's edges cut every other sample's blocks
//!
//! The walker keeps the records that overlap the position it has reached. Each new group of records
//! first closes, at the base before it, whatever the old ones had not yet written; each record's
//! end is a stop, and a record carrying a real alternate stops at EVERY base it covers. So the
//! output's records are the union of every input's edges, and a variant in one sample splits the
//! blocks of all the others around it.
//!
//! # The reference base a closed block takes is read through a window that moves
//!
//! A block that starts after the last written position takes its reference allele from the base
//! after that position, and that base was read when the position was written: from the group's
//! reference window when a new group closed it, and from a window stretched over the whole span
//! being closed when the walker closed it itself. `--ref-padding` moves the first of those windows,
//! and with it which base is read, since the window's SECOND base is the one taken.
//!
//! # The last block written picks its reference from the last record overlapping
//!
//! The stopped records are collected walking the overlap list BACKWARDS, and the block merge takes
//! the first of them. So when two samples' blocks close at one position, the block's reference allele
//! and its start come from the one added last.
//!
//! # `--convert-to-base-pair-resolution` is `--break-bands-at-multiples-of 1`
//!
//! ```java
//! if ( multipleAtWhichToBreakBands == 1 || useBpResolution) {
//!     useBpResolution = true;
//!     multipleAtWhichToBreakBands = 1;
//! }
//! ```
//!
//! Both arguments together are therefore base-pair resolution whatever the grid says, and a grid of
//! one also drops every `END`, which a grid of two does not.

use crate::reference_confidence_merger::{
    is_unmixed_mnp_ignoring_non_ref, make_genotype_call, non_ref, MergeError, MergeNotes, Merger,
    Method,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

/// What the tool's own arguments and its walker's say.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Arguments {
    /// `--convert-to-base-pair-resolution`.
    pub base_pair_resolution: bool,
    /// `--break-bands-at-multiples-of`.
    pub break_bands_at_multiples_of: i64,
    /// `--call-genotypes`.
    pub call_genotypes: bool,
    /// `--ignore-variants-starting-outside-interval`.
    pub ignore_variants_starting_outside_interval: bool,
    /// `--combine-variants-distance`.
    pub combine_variants_distance: i64,
    /// `--max-distance`.
    pub max_distance: i64,
    /// `--ref-padding`.
    pub ref_padding: i64,
}

impl Default for Arguments {
    fn default() -> Self {
        Arguments {
            base_pair_resolution: false,
            break_bands_at_multiples_of: 0,
            call_genotypes: false,
            ignore_variants_starting_outside_interval: false,
            combine_variants_distance: 0,
            max_distance: i32::MAX as i64,
            ref_padding: 1,
        }
    }
}

/// The reference the walker reads, one-based and inclusive.
pub trait Reference {
    /// The contig's length, which clamps every window.
    fn length(&self, contig: &str) -> i64;
    /// The bases from `start` to `end`.
    fn bases(&mut self, contig: &str, start: i64, end: i64) -> Vec<u8>;
}

/// A `ReferenceContext`: an interval, and the window around it that `getBases` reads.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Context {
    contig: String,
    interval: (i64, i64),
    window: (i64, i64),
}

impl Context {
    /// `setWindow(left, right)`, clamped to the contig.
    fn set_window(&mut self, left: i64, right: i64, length: i64) {
        self.window = (
            (self.interval.0 - left).max(1),
            (self.interval.1 + right).min(length),
        );
    }

    fn bases(&self, reference: &mut dyn Reference) -> Vec<u8> {
        if self.window.1 < self.window.0 {
            return Vec::new();
        }
        reference.bases(&self.contig, self.window.0, self.window.1)
    }
}

/// `Arrays.copyOfRange`: zeros past the end, and a failure when the range starts past it.
fn copy_of_range(bases: &[u8], from: i64, to: i64) -> Result<Vec<u8>, Failure> {
    if from < 0 || from > bases.len() as i64 {
        return Err(Failure::Runtime {
            class: "java.lang.ArrayIndexOutOfBoundsException".to_string(),
            message: format!("Array index out of range: {from}"),
        });
    }
    Ok((from..to)
        .map(|index| bases.get(index as usize).copied().unwrap_or(0))
        .collect())
}

/// Why the walk stopped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Failure {
    /// A `UserException`, with the whole message it prints: a `BadInput` carries its own
    /// `Bad input: ` prefix here.
    User(String),
    /// Anything else the reference throws, by class.
    Runtime { class: String, message: String },
}

impl From<MergeError> for Failure {
    fn from(error: MergeError) -> Self {
        match error {
            MergeError::MissingNonRef { .. } => Failure::User(error.message()),
            other => Failure::Runtime {
                class: other.java_class().to_string(),
                message: other.message(),
            },
        }
    }
}

/// One record and the name of the input it came from, which is `VariantContext.getSource()`.
#[derive(Debug, Clone)]
pub struct Sourced {
    pub record: VariantContext,
    pub source: String,
}

/// The sample names a record carries, `getSampleNames()`.
fn samples_of(vc: &VariantContext) -> Vec<String> {
    vc.genotypes.iter().map(|g| g.sample_name.clone()).collect()
}

/// The walk, fed one record at a time in the merged order.
pub struct Walker<'a, 'r> {
    merger: Merger<'a>,
    arguments: Arguments,
    reference: &'r mut dyn Reference,
    /// The dictionary's contig order, which `IntervalUtils.isAfter` compares by.
    contigs: Vec<String>,
    /// The user's intervals, when there were any: `isWithinInterval`'s overlap detector.
    intervals: Option<Vec<(String, i64, i64)>>,

    // `MultiVariantWalkerGroupedOnStart`.
    current: Vec<Sourced>,
    first_current_start: i64,
    last_current_start: i64,

    // `CombineGVCFs`.
    overlapping: Vec<Sourced>,
    samples: Vec<String>,
    previous: Option<(String, i64)>,
    reference_after_previous: u8,
    stored: Option<Context>,

    /// What has been written, in order.
    pub written: Vec<VariantContext>,
    pub notes: MergeNotes,
    /// Whether `onTraversalSuccess` found state left over and warned about it.
    pub cut_blocks: bool,
}

impl<'a, 'r> Walker<'a, 'r> {
    pub fn new(
        merger: Merger<'a>,
        mut arguments: Arguments,
        reference: &'r mut dyn Reference,
        contigs: Vec<String>,
        intervals: Option<Vec<(String, i64, i64)>>,
    ) -> Self {
        if arguments.break_bands_at_multiples_of == 1 || arguments.base_pair_resolution {
            arguments.base_pair_resolution = true;
            arguments.break_bands_at_multiples_of = 1;
        }
        Walker {
            merger,
            arguments,
            reference,
            contigs,
            intervals,
            current: Vec::new(),
            first_current_start: 0,
            last_current_start: 0,
            overlapping: Vec::new(),
            samples: Vec::new(),
            previous: None,
            reference_after_previous: 0,
            stored: None,
            written: Vec::new(),
            notes: MergeNotes::default(),
            cut_blocks: false,
        }
    }

    /// `isWithinInterval`: true when there are no intervals, or when one overlaps the locus.
    fn is_within_interval(&self, contig: &str, start: i64, end: i64) -> bool {
        match &self.intervals {
            None => true,
            Some(intervals) => intervals
                .iter()
                .any(|(c, s, e)| c == contig && *s <= end && start <= *e),
        }
    }

    /// `MultiVariantWalkerGroupedOnStart.apply(variant, ...)`: collect the records that start
    /// together, and hand the previous group over when this one starts elsewhere.
    pub fn accept(&mut self, variant: Sourced) -> Result<(), Failure> {
        let start = variant.record.start;
        if self.arguments.ignore_variants_starting_outside_interval
            && !self.is_within_interval(&variant.record.contig, start, start)
        {
            return Ok(());
        }
        if self.current.is_empty() {
            self.first_current_start = start;
        } else if self.current[0].record.contig != variant.record.contig
            || self.last_current_start < start - self.arguments.combine_variants_distance
            || self.first_current_start < start - self.arguments.max_distance
        {
            let group = std::mem::take(&mut self.current);
            self.apply_group(group)?;
            self.first_current_start = start;
        }
        self.current.push(variant);
        self.last_current_start = start;
        Ok(())
    }

    /// `afterTraverse` and then `onTraversalSuccess`.
    pub fn finish(&mut self) -> Result<(), Failure> {
        if !self.current.is_empty() {
            let group = std::mem::take(&mut self.current);
            self.apply_group(group)?;
        }
        if self.stored.is_none() {
            return Ok(());
        }
        if !self.overlapping.is_empty() {
            let first = &self.overlapping[0].record;
            let contig = first.contig.clone();
            let start = first.start;
            let end = self
                .overlapping
                .iter()
                .map(|s| s.record.stop)
                .max()
                .expect("an overlapping record");
            self.create_intermediate_variants(&contig, start, end)?;
            if !self.overlapping.is_empty() {
                self.cut_blocks = true;
            }
        }
        Ok(())
    }

    /// `apply(List<VariantContext>, List<ReadsContext>)`, which builds the spanning reference
    /// context and hands both to the tool.
    fn apply_group(&mut self, group: Vec<Sourced>) -> Result<(), Failure> {
        let contig = group[0].record.contig.clone();
        let min_start = group
            .iter()
            .map(|s| s.record.start)
            .min()
            .expect("a record");
        let max_end = group.iter().map(|s| s.record.stop).max().expect("a record");
        // `ReferenceContext.setWindow` refuses a negative offset, the start one first.
        if self.arguments.ref_padding < 0 {
            return Err(Failure::Runtime {
                class: "org.broadinstitute.hellbender.exceptions.GATKException".to_string(),
                message: "Reference window starts after the current interval".to_string(),
            });
        }
        let length = self.reference.length(&contig);
        let mut context = Context {
            contig: contig.clone(),
            interval: (min_start, max_end),
            window: (min_start, max_end),
        };
        context.set_window(
            self.arguments.ref_padding,
            self.arguments.ref_padding,
            length,
        );
        self.apply(group, context)
    }

    /// `CombineGVCFs.apply`.
    fn apply(&mut self, group: Vec<Sourced>, context: Context) -> Result<(), Failure> {
        for sourced in &group {
            if is_unmixed_mnp_ignoring_non_ref(&sourced.record) {
                return Err(Failure::User(format!(
                    "Bad input: Combining gVCFs containing MNPs is not supported. {} contained a MNP at {}:{}",
                    sourced.source, sourced.record.contig, sourced.record.start
                )));
            }
        }

        if !self.overlapping.is_empty() {
            let first = &self.overlapping[0].record;
            let (last_contig, last_start) = match &self.previous {
                Some((contig, start)) if *contig == first.contig => (contig.clone(), *start),
                _ => (first.contig.clone(), first.start),
            };
            let end = if last_contig == context.contig {
                context.interval.0 - 1
            } else {
                self.overlapping
                    .iter()
                    .map(|s| s.record.stop)
                    .max()
                    .expect("an overlapping record")
            };
            if end < last_start {
                return Err(Failure::Runtime {
                    class: "java.lang.IllegalArgumentException".to_string(),
                    message: format!(
                        "Invalid interval. Contig:{last_contig} start:{last_start} end:{end}"
                    ),
                });
            }
            self.create_intermediate_variants(&last_contig, last_start, end)?;
        }

        self.merge_with_new(group, &context)?;

        let replace = match &self.stored {
            None => true,
            Some(stored) => stored.contig != context.contig || stored.window.1 < context.window.1,
        };
        if replace {
            self.stored = Some(context);
        }
        Ok(())
    }

    /// `createIntermediateVariants`: every stop inside the span, each closed with the base the
    /// stretched stored window holds there.
    fn create_intermediate_variants(
        &mut self,
        contig: &str,
        start: i64,
        end: i64,
    ) -> Result<(), Failure> {
        // `resizeReferenceIfNeeded`.
        let mut stored = self.stored.clone().expect("a stored reference context");
        let left = stored.interval.0 - start;
        let right = end - stored.interval.1;
        let length = self.reference.length(&stored.contig);
        stored.set_window(left.max(1), right.max(1), length);
        self.stored = Some(stored.clone());

        let mut stops =
            intermediate_stop_sites(start, end, self.arguments.break_bands_at_multiples_of);
        for sourced in &self.overlapping {
            let vc = &sourced.record;
            if vc.alleles.len() > 2 {
                for position in vc.start..=vc.stop {
                    stops.insert(position);
                }
            } else if vc.stop <= end {
                stops.insert(vc.stop);
            }
        }

        let bases = stored.bases(self.reference);
        for stop in stops {
            if stop <= end && stop >= start && self.is_within_interval(contig, stop, stop) {
                let offset = stop - stored.window.0;
                let reference_bases = copy_of_range(&bases, offset, offset + 2)?;
                self.end_previous_states(contig, stop, &reference_bases, &[], true)?;
            }
        }
        Ok(())
    }

    /// `mergeWithNewVCs`.
    fn merge_with_new(&mut self, group: Vec<Sourced>, context: &Context) -> Result<(), Failure> {
        if group.is_empty() {
            return Ok(());
        }
        if !self.okay_to_skip(&group, context) {
            let start = context.interval.0;
            if start - 1 > 0 {
                let bases = context.bases(self.reference);
                let reference_bases = copy_of_range(&bases, 1, bases.len() as i64)?;
                if reference_bases.is_empty() {
                    return Err(Failure::Runtime {
                        class: "java.lang.ArrayIndexOutOfBoundsException".to_string(),
                        message: "Index 0 out of bounds for length 0".to_string(),
                    });
                }
                self.end_previous_states(
                    &context.contig,
                    start - 1,
                    &reference_bases,
                    &group,
                    false,
                )?;
            }
        }
        self.overlapping.extend(group);
        for sourced in &self.overlapping {
            for sample in samples_of(&sourced.record) {
                if !self.samples.contains(&sample) {
                    self.samples.push(sample);
                }
            }
        }
        Ok(())
    }

    /// `okayToSkipThisSite`: the group starts right after the last write and shares no sample with
    /// the records still open. The contigs are not compared.
    fn okay_to_skip(&self, group: &[Sourced], context: &Context) -> bool {
        let shares = group
            .iter()
            .flat_map(|s| samples_of(&s.record))
            .any(|sample| self.samples.contains(&sample));
        matches!(&self.previous, Some((_, start)) if context.interval.0 == start + 1) && !shares
    }

    /// `IntervalUtils.isAfter(pos, prevPos, dictionary)`.
    fn is_after(
        &self,
        contig: &str,
        start: i64,
        previous: &(String, i64),
    ) -> Result<bool, Failure> {
        let index = |name: &str| self.contigs.iter().position(|c| c == name);
        let (Some(first), Some(second)) = (index(contig), index(&previous.0)) else {
            return Err(Failure::Runtime {
                class: "java.lang.IllegalArgumentException".to_string(),
                message: "Can't do comparison because Locatables' contigs not found in sequence \
                          dictionary"
                    .to_string(),
            });
        };
        Ok(first > second || (first == second && start > previous.1))
    }

    /// `endPreviousStates`: stop every open record at `position`, write their merge, and keep the
    /// ones that go on.
    fn end_previous_states(
        &mut self,
        contig: &str,
        position: i64,
        reference_bases: &[u8],
        group: &[Sourced],
        force: bool,
    ) -> Result<(), Failure> {
        let new_samples: Vec<String> = group.iter().flat_map(|s| samples_of(&s.record)).collect();
        let reference_base = reference_bases[0];
        let next_base = if force {
            reference_bases.get(1).copied().unwrap_or(b'N')
        } else {
            reference_base
        };

        let mut stopped: Vec<Sourced> = Vec::new();
        let mut index = self.overlapping.len();
        while index > 0 {
            index -= 1;
            let vc = &self.overlapping[index].record;
            if vc.start <= position || vc.contig != contig {
                stopped.push(self.overlapping[index].clone());
                let names = samples_of(vc);
                let covered = !group.is_empty()
                    && !force
                    && names.iter().all(|name| new_samples.contains(name));
                if vc.stop == position || covered {
                    self.samples.retain(|sample| !names.contains(sample));
                    self.overlapping.remove(index);
                }
            }
        }

        let after = match &self.previous {
            None => true,
            Some(previous) => self.is_after(contig, position, previous)?,
        };
        if !stopped.is_empty() && after {
            let closing_contig = stopped[0].record.contig.clone();
            let records: Vec<VariantContext> = stopped.iter().map(|s| s.record.clone()).collect();
            let merged = if records.iter().any(|vc| vc.alleles.len() > 2) {
                self.merger
                    .merge(
                        &records,
                        &closing_contig,
                        position,
                        Some(reference_base),
                        false,
                        &mut self.notes,
                    )?
                    .expect("a reference base was given")
            } else {
                self.reference_block_merge(&records, position)?
            };
            self.written.push(merged);
            self.previous = Some((closing_contig, position));
            self.reference_after_previous = next_base;
        }
        Ok(())
    }

    /// `referenceBlockMerge`: one block over every stopped record, each genotype kept and called
    /// again against the reference and `<NON_REF>`.
    fn reference_block_merge(
        &self,
        records: &[VariantContext],
        end: i64,
    ) -> Result<VariantContext, Failure> {
        let first = &records[0];
        let (start, reference) = match &self.previous {
            Some((contig, previous)) if *contig == first.contig && first.start < previous + 1 => (
                previous + 1,
                Allele::create(&[self.reference_after_previous], true).map_err(|error| {
                    Failure::Runtime {
                        class: "java.lang.IllegalArgumentException".to_string(),
                        message: error.to_string(),
                    }
                })?,
            ),
            _ => (first.start, first.reference().clone()),
        };
        let alleles = vec![reference, non_ref()];

        let method = if self.merger.call_genotypes {
            Method::PreferPls
        } else {
            Method::SetToNoCall
        };
        let mut genotypes: Vec<Genotype> = Vec::new();
        for vc in records {
            for original in vc.genotypes.iter() {
                let mut genotype = original.clone();
                make_genotype_call(
                    original.ploidy(),
                    &mut genotype,
                    method,
                    original.pl.as_deref(),
                    &alleles,
                    original,
                );
                genotypes.push(genotype);
            }
        }

        let mut merged = VariantContext::new(&first.contig, start, alleles);
        merged.stop = end;
        if !self.arguments.base_pair_resolution && end != start {
            merged.attributes = vec![("END".to_string(), Value::Str(end.to_string()))];
        }
        merged.genotypes = genotypes.into();
        Ok(merged)
    }
}

/// `getIntermediateStopSites`: the base before every multiple of the grid inside the span, and the
/// first one no earlier than the grid itself (or two, for a grid of one).
pub fn intermediate_stop_sites(
    start: i64,
    end: i64,
    multiple: i64,
) -> std::collections::BTreeSet<i64> {
    let mut sites = std::collections::BTreeSet::new();
    if multiple > 0 {
        let mut block_end = if start < multiple + 1 {
            multiple.max(2)
        } else {
            (start / multiple) * multiple
        };
        while block_end <= end {
            sites.insert(block_end - 1);
            block_end += multiple;
        }
    }
    sites
}

/// htsjdk's `MergingIterator` over one record list per input, compared by
/// `VariantContextComparator`: the contig's index in the merged dictionary, then the start.
///
/// It is a `java.util.PriorityQueue` of iterators, each keyed by the record it would hand out next,
/// and `next()` polls the head, takes its record, and offers the SAME iterator back keyed by its
/// following record. Two records that compare equal come out in whatever order the heap's shape
/// puts them, which is neither input order nor stable: [`JavaQueue`] is `siftUp` and `siftDown`
/// transcribed so that shape is the reference's.
///
/// The answer is `(input, index)` for every record in the order handed out. A record that compares
/// before the one handed out just before it is the reference's `IllegalStateException`, thrown by
/// the iterator itself and so never wrapped by the walker.
pub fn merging_order(keys: &[Vec<(i64, i64)>]) -> Result<Vec<(usize, usize)>, Failure> {
    let mut cursors = vec![0usize; keys.len()];
    let mut queue = JavaQueue::default();
    for (input, list) in keys.iter().enumerate() {
        if let Some(key) = list.first() {
            queue.offer((*key, input));
        }
    }
    let mut order = Vec::new();
    let mut last: Option<(i64, i64)> = None;
    while let Some((key, input)) = queue.poll() {
        let index = cursors[input];
        cursors[input] += 1;
        if last.is_some_and(|previous| previous > key) {
            return Err(Failure::Runtime {
                class: "java.lang.IllegalStateException".to_string(),
                message: "The elements of the input Iterators are not sorted according to the \
                          comparator htsjdk.variant.variantcontext.VariantContextComparator"
                    .to_string(),
            });
        }
        if let Some(next) = keys[input].get(cursors[input]) {
            queue.offer((*next, input));
        }
        last = Some(key);
        order.push((input, index));
    }
    Ok(order)
}

/// `java.util.PriorityQueue`, ordered by the key alone: an input whose next record compares equal
/// to another's is not preferred for being named first.
#[derive(Default)]
struct JavaQueue {
    queue: Vec<((i64, i64), usize)>,
}

impl JavaQueue {
    /// `offer`: append, then `siftUp`, which stops at an equal parent.
    fn offer(&mut self, element: ((i64, i64), usize)) {
        let mut k = self.queue.len();
        self.queue.push(element);
        while k > 0 {
            let parent = (k - 1) >> 1;
            if element.0 >= self.queue[parent].0 {
                break;
            }
            self.queue[k] = self.queue[parent];
            k = parent;
        }
        self.queue[k] = element;
    }

    /// `poll`: take the head, move the last element to the root, then `siftDown`, which prefers the
    /// left child on a tie and stops at an equal child.
    fn poll(&mut self) -> Option<((i64, i64), usize)> {
        if self.queue.is_empty() {
            return None;
        }
        let result = self.queue[0];
        let last = self.queue.pop().expect("a non-empty queue");
        let size = self.queue.len();
        if size > 0 {
            let mut k = 0;
            let half = size >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let right = child + 1;
                if right < size && self.queue[child].0 > self.queue[right].0 {
                    child = right;
                }
                if last.0 <= self.queue[child].0 {
                    break;
                }
                self.queue[k] = self.queue[child];
                k = child;
            }
            self.queue[k] = last;
        }
        Some(result)
    }
}
