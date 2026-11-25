//! `SAMPileupCodec`, `SAMPileupFeature` and `SAMPileupElement`, ported from
//! `org.broadinstitute.hellbender.utils.codecs.sampileup` (GATK 4.7.0.0; all three Java files are
//! unchanged since 4.6.2.0).
//!
//! The samtools mpileup format, which `CheckPileup` reads as its truth. One line is a locus:
//! contig, position, reference base, coverage, the read bases and their qualities.
//!
//! # The bases column is a little language
//!
//! `.` and `,` are the reference base, `*` is a deletion, `$` ends a read and consumes **no**
//! quality, `^` eats the **next** character as a mapping quality, and `+`/`-` introduce an indel
//! whose length is the number that follows and whose bases are skipped whole. So the quality
//! string is consumed at a different rate from the bases string: one quality per emitted element
//! and none for the markers, which is why the two indices are independent and why leftover
//! qualities are their own error.
//!
//! What an element's base comes out as is also not what the column says. `.` and `,` become the
//! reference column's base (strand is dropped), `*` becomes `D` (`BaseUtils.Base.D`), `N`/`n`
//! become `N`, and every other base is routed through `BaseUtils.simpleBaseToBaseIndex` and back,
//! which upper-cases it. The same routing applies to the **reference** column, where it has a
//! surprise of its own: `baseIndexMap['*'] = Base.A.ordinal();` (`BaseUtils.java:49`), so a
//! reference base of `*` decodes as `A`.
//!
//! # The field count is not the format's
//!
//! ```java
//! if(tokens.length < MINIMUM_FIELDS || tokens.length > MAXIMUM_FIELDS) {
//! ```
//!
//! `MINIMUM_FIELDS` is 4 and `MAXIMUM_FIELDS` is 6, where a samtools mpileup line has six columns
//! and may have more. The reference's own reason for the 4 is a comment above it
//! (`SAMPileupCodec.java:33`):
//! `if the coverage is 0 sometimes bases and qualities does not appear`.
//!
//! A **seven**-column line is therefore refused. A **five**-column line is worse: it passes this
//! check and then reads `tokens[5]`, which is an `ArrayIndexOutOfBoundsException` rather than the
//! codec's own exception (a four-column line with a non-zero coverage does the same at
//! `tokens[4]`). Both are reproduced, because a caller that catches the codec's error type would
//! not catch the second.
//!
//! Two of the messages carry typos, "THe SAM pileup line" and "this codes is only valid". They are
//! kept: a port that fixes the spelling of a message changes what a caller matching on it sees.
//!
//! # Reading this file without reading Rust
//!
//! - `Result<T, PileupCodecError>` is "a `T` or the exception the reference would throw";
//!   `Err(..)` is the throw. A `?` after a call rethrows that call's `Err` to this function's
//!   caller, which is Java's exception propagation written at the one call site.
//! - `Option<T>` is "a value or nothing". `slice.get(i)` returns one instead of throwing on a bad
//!   index, and `.ok_or(e)` / `.ok_or_else(|| e)` turn the nothing into the `Err(e)` the reference
//!   would have thrown. `|| ...` is a closure, Java's lambda.
//! - `&str` is a borrowed string; `.as_bytes()` views it as the bytes `u8` (Java's `byte`, but
//!   unsigned). `b'$'` is the byte value of the character `$`, Java's `(byte) '$'`.
//! - `match` is Java's `switch` without fall-through, and the compiler requires it to cover every
//!   value, which is what the final `other =>` arm (Java's `default:`) is for.
//! - A function's last expression is its return value; `return` is only written for early exits.

/// `MINIMUM_FIELDS`, which is 4 rather than the format's 6: see the module note.
pub const MINIMUM_FIELDS: usize = 4;
/// `MAXIMUM_FIELDS`.
pub const MAXIMUM_FIELDS: usize = 6;

/// `SAMPileupElement`: one base of one read at this locus, with its quality.
///
/// `#[derive(...)]` asks the compiler to write `equals`, copying and debug printing field by
/// field, so two elements are equal exactly when base and quality are. `qual` is already the
/// Phred value (the file's character minus 33), not the character.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SamPileupElement {
    pub base: u8,
    pub qual: u8,
}

/// `SAMPileupFeature`: everything one line describes.
///
/// The coverage column is not kept: `decode` checks it against the element count and then drops
/// it: `return new SAMPileupFeature(tokens[0], pos, ref, pileupElements);` takes no coverage.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SamPileupFeature {
    pub contig: String,
    /// 1-based, as the file writes it.
    pub position: i32,
    pub reference_base: u8,
    pub elements: Vec<SamPileupElement>,
}

impl SamPileupFeature {
    /// `size()`, which is the number of elements and not the coverage column.
    pub fn size(&self) -> usize {
        self.elements.len()
    }

    /// `getBasesString()`: the decoded bases, so `.` and `,` read as the reference base and `*`
    /// as `D`, not the bases column as the file wrote it.
    pub fn bases_string(&self) -> String {
        // `.iter().map(..).collect()` is Java's `stream().map(..).collect(..)`; `as char` reads
        // each byte as the character with that code, which is what `StringUtil.bytesToString`
        // does for these ASCII bases.
        self.elements.iter().map(|e| e.base as char).collect()
    }

    /// `getBaseQuals()`.
    pub fn base_quals(&self) -> Vec<u8> {
        self.elements.iter().map(|e| e.qual).collect()
    }
}

/// What the codec refuses, and how.
///
/// An `enum` here is a closed set of variants that may each carry data, closer to a sealed
/// class hierarchy than to a Java `enum`. Two variants because the reference throws two different
/// classes, and the class is part of what the golden `sam_pileup.txt.gz` records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PileupCodecError {
    /// `CodecLineParsingException` with the codec's own message.
    Parsing(String),
    /// The raw `ArrayIndexOutOfBoundsException` a five-column line produces, with Java's own text.
    IndexOutOfBounds { index: usize, length: usize },
}

impl PileupCodecError {
    /// The message the exception carries. For the index case that is the JDK's own text, which
    /// the golden row `five-columns` holds as `Index 5 out of bounds for length 5`.
    pub fn message(&self) -> String {
        match self {
            PileupCodecError::Parsing(text) => text.clone(),
            PileupCodecError::IndexOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
        }
    }

    /// The Java class the reference throws, which is not the same for the two.
    pub fn java_class(&self) -> &'static str {
        match self {
            PileupCodecError::Parsing(_) => "htsjdk.tribble.exception.CodecLineParsingException",
            PileupCodecError::IndexOutOfBounds { .. } => "java.lang.ArrayIndexOutOfBoundsException",
        }
    }
}

/// `SAMPileupCodec.decode(line)`.
///
/// **What**: one mpileup line as a feature, or the exception the reference throws for it.
///
/// **How**: the Java's steps in the Java's order: split, count the columns, parse position,
/// reference base and coverage, stop at coverage zero, then walk bases and qualities together and
/// check the element count against the coverage. The order is what decides which of two faults in
/// one line is reported.
pub fn decode(line: &str) -> Result<SamPileupFeature, PileupCodecError> {
    // `SPLIT_PATTERN.split(line.trim(), -1)`, keeping trailing empty fields as `-1` does.
    // The Java pattern is `Pattern.compile("\\t|( +)")`, a tab OR a run of spaces, and this
    // splits on tabs only: a space-separated line the reference decodes is refused here as having
    // too few columns. No input line in the golden is space-delimited.
    let tokens: Vec<&str> = line.trim().split('\t').collect();
    if tokens.len() < MINIMUM_FIELDS || tokens.len() > MAXIMUM_FIELDS {
        return Err(PileupCodecError::Parsing(format!(
            "The SAM pileup line didn't have the expected number of columns ({MINIMUM_FIELDS}-{MAXIMUM_FIELDS}): {line}. \
Note that this codes is only valid for single-sample pileups",
        )));
    }

    let contig = tokens[0].to_string();
    let position = parse_integer(tokens[1], "position")?;
    let reference_base = parse_base_token(tokens[2], "reference")?;
    let coverage = parse_integer(tokens[3], "coverage")?;

    // Coverage zero returns before the two columns after it are looked at, so a line whose bases
    // and qualities are nonsense parses as long as its coverage says zero.
    if coverage == 0 {
        return Ok(SamPileupFeature {
            contig,
            position,
            reference_base,
            elements: Vec::new(),
        });
    }

    // Five columns reach here and index past the end, which the reference does not guard.
    // `tokens[4]` and `tokens[5]` in the Java, evaluated in that order, so four columns fail on
    // index 4 and five on index 5.
    let bases = tokens.get(4).ok_or(PileupCodecError::IndexOutOfBounds {
        index: 4,
        length: tokens.len(),
    })?;
    let qualities = tokens.get(5).ok_or(PileupCodecError::IndexOutOfBounds {
        index: 5,
        length: tokens.len(),
    })?;

    let elements = parse_bases_and_quals(bases, qualities, reference_base)?;
    // `cov != pileupElements.size()`. `as usize` reinterprets the `i32` as an unsigned index, so
    // a negative coverage becomes a huge number; it still differs from any real count, and the
    // line is refused with the same message as in Java.
    if coverage as usize != elements.len() {
        return Err(PileupCodecError::Parsing(format!(
            "THe SAM pileup line didn't have the same number of elements as the expected coverage: {coverage}"
        )));
    }
    Ok(SamPileupFeature {
        contig,
        position,
        reference_base,
        elements,
    })
}

/// `parseBasesAndQuals(bases, qualities, ref)`.
///
/// **What**: the elements of one locus, one per base that consumes a quality.
///
/// **How**: two cursors, `i` over the bases column and `j` over the qualities, advanced at their
/// own rates as the module note describes. The Java wraps the whole walk in
/// `catch(IndexOutOfBoundsException e)`, so running out of qualities is reported with the codec's
/// own message; here each quality read is `qualities.get(j)`, which yields nothing past the end,
/// and that nothing becomes the same message.
pub fn parse_bases_and_quals(
    bases: &str,
    qualities: &str,
    reference_base: u8,
) -> Result<Vec<SamPileupElement>, PileupCodecError> {
    let bases = bases.as_bytes();
    let qualities = qualities.as_bytes();
    let mut elements = Vec::with_capacity(qualities.len());
    let mut i = 0usize;
    let mut j = 0usize;

    // The one error the reference reports by catching IndexOutOfBounds, so it has a message of its
    // own rather than Java's.
    let missing = || {
        PileupCodecError::Parsing(
            "Malformed SAM pileup: Different number of bases and qualities found.".to_string(),
        )
    };

    while i < bases.len() {
        match bases[i] {
            // End of read: no element and no quality.
            b'$' => {}
            // Start of read: the next character is a mapping quality, not a base.
            b'^' => i += 1,
            b'.' | b',' => {
                // `qualities.charAt(j++)`. `*` reads the byte out of the reference `get` returns.
                let qual = *qualities.get(j).ok_or_else(missing)?;
                j += 1;
                elements.push(SamPileupElement {
                    base: reference_base,
                    qual: fastq_to_phred(qual),
                });
            }
            b'*' => {
                let qual = *qualities.get(j).ok_or_else(missing)?;
                j += 1;
                // `BaseUtils.Base.D.base`, which is not a base the reads column would accept.
                elements.push(SamPileupElement {
                    base: b'D',
                    qual: fastq_to_phred(qual),
                });
            }
            b'+' | b'-' => {
                // The Java matches `Pattern.compile("([0-9]+).*")` against everything after the
                // marker and takes group 1: the leading run of ASCII digits, which is what the
                // `take_while` counts. `&bases[i + 1..]` is a view of the tail, not a copy.
                let rest = &bases[i + 1..];
                let digits = rest.iter().take_while(|b| b.is_ascii_digit()).count();
                if digits == 0 {
                    return Err(PileupCodecError::Parsing(
                        "The SAM pileup line has an indel marker (+/-) without length".to_string(),
                    ));
                }
                let length_text = std::str::from_utf8(&rest[..digits]).expect("ASCII digits");
                let indel_length = parse_integer(length_text, "indel-length")?;
                // The number's own characters plus that many bases, and the `+`/`-` comes from the
                // increment at the end of the loop. Nothing checks that the indel's bases are
                // there: a length past the end of the column leaves `i` past the end too.
                i += indel_length as usize + digits;
            }
            other => {
                let base = parse_base(other, "reads String")?;
                let qual = *qualities.get(j).ok_or_else(missing)?;
                j += 1;
                elements.push(SamPileupElement {
                    base,
                    qual: fastq_to_phred(qual),
                });
            }
        }
        i += 1;
    }

    // The Java tests `i != bases.length() || j != qualities.length()`, and only the second half
    // is ported. The first half CAN fail: a `^` as the last character, or an indel length that
    // runs past the end, leaves `i` beyond `bases.length()` when the loop stops. For
    // `chr1\t10\tA\t1\t.^\tI` the reference throws "Not all bases/qualities have been parsed
    // because of a malformed line" and this returns one element. No golden row covers it.
    if j != qualities.len() {
        return Err(PileupCodecError::Parsing(
            "Not all bases/qualities have been parsed because of a malformed line".to_string(),
        ));
    }
    Ok(elements)
}

/// `SAMPileupCodec.canDecode(path)`: the extension and nothing else.
///
/// One block-compressed extension is stripped first, and the comparison is case-insensitive, so
/// `x.mpileup.gz` and `x.PILEUP` decode while `x.pileup.txt` and a bare `pileup` do not.
///
/// ```java
/// if (IOUtil.hasBlockCompressedExtension(path)) {
///     noBlockCompressedPath = FilenameUtils.removeExtension(path).toLowerCase();
/// ```
///
/// `IOUtil.hasBlockCompressedExtension` tests `FileExtensions.BLOCK_COMPRESSED`, which in htsjdk
/// 5.0.0 is `".gz", ".gzip", ".bgz", ".bgzf"`. The list below lacks `.gzip`, so
/// `x.pileup.gzip` is `true` in the reference and `false` here. That helper also strips an
/// `http://` or `https://` query string before testing; this does not. No golden row has either.
pub fn can_decode(path: &str) -> bool {
    // Lower-casing first and then cutting at the last `.` gives the same stem as the Java's
    // cut-then-lower-case: the cut point is the `.` of the compression extension either way.
    let lower = path.to_lowercase();
    let stripped = if lower.ends_with(".gz") || lower.ends_with(".bgz") || lower.ends_with(".bgzf")
    {
        match lower.rfind('.') {
            Some(index) => lower[..index].to_string(),
            None => lower.clone(),
        }
    } else {
        lower.clone()
    };
    stripped.ends_with(".pileup") || stripped.ends_with(".mpileup")
}

/// `SAMUtils.fastqToPhred`: the quality character minus 33.
///
/// The Java refuses a character outside the printable range before subtracting:
///
/// ```java
/// if (ch < 33 || ch > 126) {
///     throw new IllegalArgumentException("Invalid fastq character: " + ch);
/// ```
///
/// That check is not ported, and the exception is not an `IndexOutOfBoundsException`, so the
/// codec's `catch` would not have turned it into a codec message. A quality of DEL (127) is
/// thrown on by the reference and decodes as 94 here; a byte below 33 underflows the `u8`
/// subtraction, which panics in a debug build and wraps in a release one (the workspace sets no
/// `overflow-checks`). No golden row has such a quality.
fn fastq_to_phred(character: u8) -> u8 {
    character - 33
}

/// `parseInteger(token, parsedValue)`: `Integer.parseInt`, with its `NumberFormatException`
/// rethrown as the codec's message. `parse::<i32>()` has the same range and, like `parseInt`,
/// accepts one leading `+` or `-`; `.map_err` swaps the error for the codec's.
fn parse_integer(token: &str, parsed_value: &str) -> Result<i32, PileupCodecError> {
    token.parse::<i32>().map_err(|_| {
        PileupCodecError::Parsing(format!(
            "The SAM pileup line had unexpected {parsed_value}: {token}"
        ))
    })
}

/// `parseBase(String, parsedValue)`: the reference column, which must be exactly one character.
///
/// `token.len()` counts bytes where Java's `length()` counts UTF-16 units; the two agree for
/// ASCII, which is every input the golden has.
fn parse_base_token(token: &str, parsed_value: &str) -> Result<u8, PileupCodecError> {
    if token.len() != 1 {
        return Err(PileupCodecError::Parsing(format!(
            "The SAM pileup line had unexpected base at {parsed_value}: {token}"
        )));
    }
    parse_base(token.as_bytes()[0], parsed_value)
}

/// `parseBase(byte, parsedValue)`: an N folds to N, everything else goes through the base index
/// and comes back upper case.
///
/// The index map knows `ACGT`, `acgt` and `*` (as A), and answers `-1` for anything else, which
/// is the refusal. That is why the golden row `bad-reference` refuses `Z`, and why a reference
/// column of `*` is accepted as `A`.
fn parse_base(base: u8, parsed_value: &str) -> Result<u8, PileupCodecError> {
    if base == b'N' || base == b'n' {
        return Ok(b'N');
    }
    let index = crate::base_utils::simple_base_to_base_index(base);
    if index == -1 {
        return Err(PileupCodecError::Parsing(format!(
            "The SAM pileup line had wrong base at {parsed_value}: {}",
            base as char
        )));
    }
    Ok(crate::base_utils::base_index_to_simple_base(index))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_dot_is_the_reference_base_and_a_star_is_a_d() {
        let feature = decode("chr1\t10\tA\t3\t.*.\tIII").expect("it parses");
        assert_eq!(feature.bases_string(), "ADA");
        assert_eq!(feature.base_quals(), vec![40, 40, 40]);
    }

    #[test]
    fn the_markers_consume_no_quality() {
        // A read start eats the character after it; a read end eats nothing.
        assert_eq!(
            decode("chr1\t10\tA\t2\t^I.,\tII").expect("parses").size(),
            2
        );
        assert_eq!(decode("chr1\t10\tA\t2\t.$,\tII").expect("parses").size(), 2);
        // An indel's own bases are skipped whole, however long the number is.
        assert_eq!(
            decode("chr1\t10\tA\t2\t.-10ACGTACGTAC,\tII")
                .expect("parses")
                .size(),
            2
        );
    }

    #[test]
    fn seven_columns_are_refused_and_five_index_out_of_bounds() {
        let seven = decode("chr1\t10\tA\t2\t.,\tII\t~~").unwrap_err();
        assert_eq!(
            seven.java_class(),
            "htsjdk.tribble.exception.CodecLineParsingException"
        );
        let five = decode("chr1\t10\tA\t2\t.,").unwrap_err();
        assert_eq!(
            five.java_class(),
            "java.lang.ArrayIndexOutOfBoundsException"
        );
        assert_eq!(five.message(), "Index 5 out of bounds for length 5");
    }

    #[test]
    fn a_coverage_of_zero_never_looks_at_the_columns_after_it() {
        let feature = decode("chr1\t10\tA\t0\tZZZZ\t!!!!").expect("it parses");
        assert_eq!(feature.size(), 0);
        assert_eq!(feature.bases_string(), "");
    }

    #[test]
    fn the_extension_is_all_can_decode_looks_at() {
        assert!(can_decode("x.pileup"));
        assert!(can_decode("x.mpileup.gz"));
        assert!(can_decode("x.PILEUP"));
        assert!(!can_decode("x.pileup.txt"));
        assert!(!can_decode("pileup"));
    }
}
