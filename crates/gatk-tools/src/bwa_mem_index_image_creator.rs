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

/// `BwaMemIndex.FASTA_FILE_EXTENSIONS`, in its own order.
pub const FASTA_EXTENSIONS: [&str; 2] = [".fasta", ".fa"];

/// `BwaMemIndex.resolveFastaFileExtension`: the name must END in one of the two extensions, so a
/// gzipped FASTA is refused by name before anything is read.
pub fn check_extension(reference: &str) -> Result<(), String> {
    if FASTA_EXTENSIONS.iter().any(|ext| reference.ends_with(ext)) {
        Ok(())
    } else {
        Err(format!(
            "the fasta file provided '{reference}' does not have any of the standard fasta extensions: {}",
            FASTA_EXTENSIONS.join(", ")
        ))
    }
}

/// `CouldNotReadReferenceException`'s class, which `handleNonUserException` prints.
pub const COULD_NOT_READ_REFERENCE: &str =
    "org.broadinstitute.hellbender.utils.bwa.CouldNotReadReferenceException";
/// `InvalidFileFormatException`'s class.
pub const INVALID_FILE_FORMAT: &str =
    "org.broadinstitute.hellbender.utils.bwa.InvalidFileFormatException";
/// `CouldNotCreateIndexImageException`'s class.
pub const COULD_NOT_CREATE_INDEX_IMAGE: &str =
    "org.broadinstitute.hellbender.utils.bwa.CouldNotCreateIndexImageException";

/// `nonEmptyReadableFile` failing: a missing file, a directory and an empty file are one refusal.
pub fn unreachable_reference(reference: &str) -> String {
    cannot_read_reference(reference, "input file unreachable or not a file")
}

/// `InvalidFileFormatException(file, detail)`, whose message names the file and not the detail.
pub fn invalid_format(reference: &str) -> String {
    format!("file {reference}: invalid format")
}

/// `CouldNotCreateIndexImageException(image, reason)`.
pub fn cannot_create_image(image: &str, reason: &str) -> String {
    format!("could not create index image '{image}': {reason}")
}

/// What `assertCanCreateOrOverwriteImageFile` says of an image path that exists and is a directory
/// or cannot be written.
pub const EXISTS_UNWRITABLE: &str = "already exists as a non-regular or unwritable file";

/// `Character.isSpaceChar`: the Unicode space, line and paragraph separators, which is NOT white
/// space: a newline or a tab is not one.
pub fn is_java_space_char(c: char) -> bool {
    matches!(
        c,
        '\u{20}' | '\u{a0}' | '\u{1680}' | '\u{2000}'
            ..='\u{200a}' | '\u{2028}' | '\u{2029}' | '\u{202f}' | '\u{205f}' | '\u{3000}'
    )
}

/// `assertLooksLikeFastaFile`'s scan: among the first 4092 characters, the first that is not a
/// space character must be `>`. A file of spaces alone passes.
pub fn looks_like_fasta(text: &str) -> bool {
    for c in text.chars().take(4092) {
        if is_java_space_char(c) {
            continue;
        }
        return c == '>';
    }
    true
}
