//! `CheckTerminatorBlock`, ported from `picard.sam.CheckTerminatorBlock` (Picard 3.4.0).
//!
//! Four lines over `BlockCompressedInputStream.checkTermination`, which
//! [`htsjdk_bgzf::termination::check_termination`] already carries. What is here is the tool's own
//! decision: which of the three answers is a failure.
//!
//! ```java
//! final FileTermination term = BlockCompressedInputStream.checkTermination(INPUT);
//! System.err.println(term.name());
//! if (term == FileTermination.DEFECTIVE) { return 100; } else { return 0; }
//! ```
//!
//! # A file with no terminator passes
//!
//! Only `DEFECTIVE` is a failure. A file whose terminator was cut off still has a healthy last
//! block, and this tool exits zero on it: the question it answers is whether the file was
//! TRUNCATED MID-BLOCK, not whether it was closed properly. A file that was never gzip at all also
//! answers `DEFECTIVE`, because the backwards search simply finds no preamble.
//!
//! # The check never decompresses
//!
//! A complete file with a flipped payload byte answers `HAS_TERMINATOR_BLOCK` and exits zero. The
//! `corrupt-payload` row of the golden is what says so, and a port that verified the CRC would
//! disagree with the reference on a real file.

use htsjdk_bgzf::termination::{check_termination, FileTermination};

/// The tool's return code for a defective file.
pub const DEFECTIVE_RETURN_CODE: i32 = 100;

/// The enum constant's name, which the tool prints on standard error.
pub fn termination_name(termination: FileTermination) -> &'static str {
    match termination {
        FileTermination::HasTerminatorBlock => "HAS_TERMINATOR_BLOCK",
        FileTermination::HasHealthyLastBlock => "HAS_HEALTHY_LAST_BLOCK",
        FileTermination::Defective => "DEFECTIVE",
    }
}

/// `doWork()`: the termination, and the return code that follows from it.
pub fn check(data: &[u8]) -> (FileTermination, i32) {
    let termination = check_termination(data);
    let code = if termination == FileTermination::Defective {
        DEFECTIVE_RETURN_CODE
    } else {
        0
    };
    (termination, code)
}
