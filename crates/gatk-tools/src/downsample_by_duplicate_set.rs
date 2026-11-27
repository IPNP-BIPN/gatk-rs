//! `DownsampleByDuplicateSet`, ported from
//! `org.broadinstitute.hellbender.tools.walkers.consensus.DownsampleByDuplicateSet` and
//! `org.broadinstitute.hellbender.engine.DuplicateSetWalker` (GATK 4.6.2.0).
//!
//! Whole molecules dropped rather than reads, so that a mixture keeps its family-size
//! distribution. The tool is a hundred and five lines over a hundred and fifty of walker, and the
//! walker is where the surprises are.
//!
//! # The last duplicate set of the file escapes every rejection rule
//!
//! ```java
//! private void processLastReadSet(){
//!     if (currentReadsWithSameUMI.getReads().size() > 0){
//!         apply(currentReadsWithSameUMI, ...);
//!     }
//! }
//! ```
//!
//! No `rejectSet` call. A trailing set that is too small, or has an odd number of reads, or is
//! short on one strand, is offered to the tool anyway: at `--min-reads 4` over ten two-read
//! molecules, exactly one molecule is written, and it is the last.
//!
//! # A rejected set does not consume a random draw
//!
//! The rejection happens before `apply`, so a rejectable molecule at the FRONT of a file leaves
//! every later molecule with the decision it had before. Adding one anywhere the walker accepts it
//! would have shifted the whole sequence instead.
//!
//! # A set with an odd number of reads is rejected at the defaults
//!
//! `size() % 2 == 1` is one of the three rules, described in the source as checking that the set is
//! paired, so a three-read molecule is dropped whatever the minimums say.

use gatk_engine::java_random::JavaRandom;

/// `RANDOM_SEED`.
pub const RANDOM_SEED: i64 = 142;

/// The arguments, with the walker's defaults beside the tool's own.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Arguments {
    pub fraction_to_keep: f64,
    /// `DEFAULT_MINIMUM_READS_PER_SET`.
    pub minimum_reads: usize,
    /// `DEFAULT_MINIMUM_READS_PER_STRAND`.
    pub minimum_reads_per_strand: usize,
}

impl Arguments {
    /// The defaults, for a given fraction. `--fraction-to-keep` is required, so it has none.
    pub fn keeping(fraction_to_keep: f64) -> Self {
        Arguments {
            fraction_to_keep,
            minimum_reads: 1,
            minimum_reads_per_strand: 0,
        }
    }
}

/// One read, reduced to what the walker reads from it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Read {
    pub name: String,
    /// The molecule number, the part of `MI:Z:<number>/<strand>` before the slash.
    pub molecule: i32,
    /// The strand suffix, which only feeds the per-strand minimum.
    pub strand: String,
}

/// What the walk refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DownsampleError {
    /// A molecule number went backwards.
    NotSortedByMoleculeId,
    /// No read reached the walker, so the last set was never started.
    NoReads,
    /// A JVM exception thrown while a read's `MI` is parsed or added to its set.
    Jvm {
        class: &'static str,
        message: String,
    },
}

impl DownsampleError {
    pub fn java_class(&self) -> &'static str {
        match self {
            DownsampleError::NotSortedByMoleculeId => {
                "org.broadinstitute.hellbender.exceptions.UserException"
            }
            DownsampleError::NoReads => "java.lang.NullPointerException",
            DownsampleError::Jvm { class, .. } => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            DownsampleError::NotSortedByMoleculeId => {
                "The input bam must be sorted by the molecule ID (MI) tag.".to_string()
            }
            // `processLastReadSet` asks the set it never started for its reads.
            DownsampleError::NoReads => "Cannot invoke \"org.broadinstitute.hellbender.tools.walkers.consensus.ReadsWithSameUMI.getReads()\" because \"this.currentReadsWithSameUMI\" is null".to_string(),
            DownsampleError::Jvm { message, .. } => message.clone(),
        }
    }
}

/// `rejectSet`: the three rules, in the order the reference applies them.
pub fn reject_set(set: &[Read], arguments: &Arguments) -> bool {
    // `MoleculeID.countStrands`, which counts the two suffixes separately.
    let a = set.iter().filter(|read| read.strand == "A").count();
    let b = set.iter().filter(|read| read.strand == "B").count();
    if a.min(b) < arguments.minimum_reads_per_strand {
        return true;
    }
    // "Check that the read set is paired", which an odd count fails whatever the minimums say.
    if set.len() % 2 == 1 {
        return true;
    }
    set.len() < arguments.minimum_reads
}

/// The whole traversal: reads grouped by molecule number, each group offered or rejected, and the
/// last group offered without ever being asked.
///
/// Returns the reads written, in order.
pub fn run(reads: &[Read], arguments: &Arguments) -> Result<Vec<Read>, DownsampleError> {
    let tagged: Vec<Tagged> = reads
        .iter()
        .map(|read| Tagged {
            mi: Some(format!("{}/{}", read.molecule, read.strand)),
            contig: String::new(),
        })
        .collect();
    match walk(&tagged, arguments) {
        // An empty walk writes nothing; the refusal is the tool's, not the grouping's.
        (written, None | Some(DownsampleError::NoReads)) => {
            Ok(written.iter().map(|&index| reads[index].clone()).collect())
        }
        (_, Some(error)) => Err(error),
    }
}

/// A written BAM and its index, with the walk's refusal if it stopped. `closeTool` closes the
/// writer whatever happened, so a walk that throws still leaves the sets applied before it.
/// `None` is an `MI` the reference renders as an identity hash, which no port can reproduce.
pub type Downsampled = Option<((Vec<u8>, Option<Vec<u8>>), Option<DownsampleError>)>;

/// One read as the walker meets it: its `MI` as `getAttributeAsString` renders it, and its contig.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tagged {
    pub mi: Option<String>,
    pub contig: String,
}

/// `MITag.split("/")`: no slash leaves the whole string, and trailing empty fields are dropped.
fn split(tag: &str) -> Vec<&str> {
    if !tag.contains('/') {
        return vec![tag];
    }
    let mut fields: Vec<&str> = tag.split('/').collect();
    while fields.last() == Some(&"") {
        fields.pop();
    }
    fields
}

fn field(mi: &Option<String>, index: usize) -> Result<&str, DownsampleError> {
    let tag = mi.as_deref().ok_or_else(|| DownsampleError::Jvm {
        class: "java.lang.NullPointerException",
        message: "Cannot invoke \"String.split(String)\" because \"MITag\" is null".to_string(),
    })?;
    let fields = split(tag);
    fields
        .get(index)
        .copied()
        .ok_or_else(|| DownsampleError::Jvm {
            class: "java.lang.ArrayIndexOutOfBoundsException",
            message: format!("Index {index} out of bounds for length {}", fields.len()),
        })
}

/// `MoleculeID.getMoleculeNumberOfRead`.
fn number(read: &Tagged) -> Result<i32, DownsampleError> {
    let text = field(&read.mi, 0)?;
    // shortcut: Rust's i32 parse is Integer.parseInt for ASCII digits only; parseInt also takes
    // other Unicode digits, which an MI tag has never been seen to carry.
    text.parse().map_err(|_| DownsampleError::Jvm {
        class: "java.lang.NumberFormatException",
        message: format!("For input string: \"{text}\""),
    })
}

/// `MoleculeID.getStrandOfRead`.
fn strand(read: &Tagged) -> Result<String, DownsampleError> {
    field(&read.mi, 1).map(str::to_string)
}

/// `DuplicateSetWalker.apply` over every read, then `processLastReadSet`, parsing `MI` where and
/// when the reference does: a read's number as it arrives, its strand when it starts a set or
/// when `rejectSet` counts the strands of the set it joined, and never for a later read of the
/// last set. Returns the indices of the reads written and the exception that stopped the walk.
pub fn walk(reads: &[Tagged], arguments: &Arguments) -> (Vec<usize>, Option<DownsampleError>) {
    let mut rng = JavaRandom::new(RANDOM_SEED);
    let mut written = Vec::new();
    let mut current: Vec<usize> = Vec::new();
    let mut held = 0;
    let mut offer = |set: &[usize], written: &mut Vec<usize>| {
        if rng.next_double() < arguments.fraction_to_keep {
            written.extend_from_slice(set);
        }
    };
    // `new ReadsWithSameUMI(read)`, whose `MoleculeID` reads both halves.
    let start = |index: usize| -> Result<i32, DownsampleError> {
        let molecule = number(&reads[index])?;
        strand(&reads[index])?;
        Ok(molecule)
    };
    for (index, read) in reads.iter().enumerate() {
        if current.is_empty() {
            match start(index) {
                Ok(molecule) => held = molecule,
                Err(error) => return (written, Some(error)),
            }
            current.push(index);
            continue;
        }
        let molecule = match number(read) {
            Ok(molecule) => molecule,
            Err(error) => return (written, Some(error)),
        };
        if held > molecule {
            return (written, Some(DownsampleError::NotSortedByMoleculeId));
        }
        if held < molecule {
            // `MoleculeID.countStrands`, which reads every strand of the closing set.
            let mut set = Vec::with_capacity(current.len());
            for &member in &current {
                match strand(&reads[member]) {
                    Ok(strand) => set.push(Read {
                        name: member.to_string(),
                        molecule: held,
                        strand,
                    }),
                    Err(error) => return (written, Some(error)),
                }
            }
            if !reject_set(&set, arguments) {
                offer(&current, &mut written);
            }
            match start(index) {
                Ok(molecule) => held = molecule,
                Err(error) => return (written, Some(error)),
            }
            current = vec![index];
            continue;
        }
        // `addRead`, whose second check is `interval.contigsMatch(read)`.
        let contig = &reads[current[0]].contig;
        if *contig != read.contig {
            return (
                written,
                Some(DownsampleError::Jvm {
                    class: "java.lang.IllegalStateException",
                    message: format!(
                        "Adding a read from another contig: set contig = {contig}, read contig = {}",
                        read.contig
                    ),
                }),
            );
        }
        current.push(index);
    }
    if current.is_empty() {
        return (written, Some(DownsampleError::NoReads));
    }
    offer(&current, &mut written);
    (written, None)
}

/// The tool end to end: the reads the traversal hands over, walked, and what was written.
pub fn downsample_with(
    source: &gatk_engine::reads::ReadsDataSource,
    options: &crate::sam_output::Options,
    filter: &dyn Fn(&htsjdk_bam::record::BamRecord) -> bool,
    level: u32,
    deflater: htsjdk_bgzf::Deflater,
    arguments: &Arguments,
) -> Result<Downsampled, gatk_engine::reads::ReadsError> {
    let records = crate::read_walker::traverse(source, &options.intervals, filter)?;
    let sequences = &source.header().sequences;
    let mut reads = Vec::with_capacity(records.len());
    for record in &records {
        let mi = crate::transfer_read_tags::attribute_as_string(record, "MI");
        if mi.is_none() && record.tags.get(htsjdk_bam::tag::Tag::new(b"MI")).is_some() {
            return Ok(None);
        }
        let contig = usize::try_from(record.reference_index)
            .ok()
            .and_then(|index| sequences.get(index))
            .map(|sequence| sequence.name.clone())
            .unwrap_or_default();
        reads.push(Tagged { mi, contig });
    }
    let (written, error) = walk(&reads, arguments);
    let kept: Vec<htsjdk_bam::record::BamRecord> = written
        .iter()
        .map(|&index| records[index].clone())
        .collect();
    let header = crate::sam_output::header_for_sam_writer(
        source.header(),
        "GATK DownsampleByDuplicateSet",
        options,
    );
    let files = crate::sam_output::write_records_with(
        &header,
        &kept,
        options.create_output_bam_index,
        level,
        deflater,
    )?;
    Ok(Some((files, error)))
}
