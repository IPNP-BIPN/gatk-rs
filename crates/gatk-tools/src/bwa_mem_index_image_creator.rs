//! `BwaMemIndexImageCreator`: where the image lands, and why its bytes are not a claim.
//!
//! Building the index is BWA's, through JNI, and is not ported. What is ported is the naming and
//! the refusal, which is the whole of the tool's own code.
//!
//! Ported from `org.broadinstitute.hellbender.tools.BwaMemIndexImageCreator` and
//! `org.broadinstitute.hellbender.utils.bwa.BwaMemIndex` in GATK 4.6.2.0.

/// The extension `doWork` appends when `--output` is not given.
pub const IMAGE_EXTENSION: &str = ".img";

/// `doWork`'s first line: the default output is the input's WHOLE name plus `.img`.
///
/// `reference.fasta` becomes `reference.fasta.img`, not `reference.img`, which is the same rule
/// [`crate::create_hadoop_bam_splitting_index`] follows and the opposite of `BuildBamIndex`'s.
pub fn default_output(reference: &str) -> String {
    format!("{reference}{IMAGE_EXTENSION}")
}

/// `BwaMemIndex.createIndexImageFromFastaFile`, on a reference it cannot read.
///
/// The message is the native side's and names the file and the reason it gave.
pub fn cannot_read_reference(path: &str, reason: &str) -> String {
    format!("cannot read the reference file '{path}': {reason}")
}

/// Whether two images of one reference may be compared byte for byte.
///
/// They may not. The file carries in-process pointers, so two builds of one reference in one
/// process differ in a handful of bytes and two runs differ again under another address layout.
/// The constant is here so a caller asking for a byte comparison finds the answer rather than the
/// silence that would let it write one.
pub const IMAGE_BYTES_ARE_REPRODUCIBLE: bool = false;
