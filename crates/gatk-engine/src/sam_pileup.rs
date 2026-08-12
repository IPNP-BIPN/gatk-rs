//! `SAMPileupCodec`, `SAMPileupFeature` and `SAMPileupElement`, ported from
//! `org.broadinstitute.hellbender.utils.codecs.sampileup` (GATK 4.6.2.0).
//!
//! The samtools mpileup format, which `CheckPileup` reads as its truth. One line is a locus: contig,
//! position, reference base, coverage, the read bases and their qualities.
//!
//! # The bases column is a little language
//!
//! `.` and `,` are the reference base, `*` is a deletion, `$` ends a read and consumes **no**
//! quality, `^` eats the **next** character as a mapping quality, and `+`/`-` introduce an indel
//! whose length is the number that follows and whose bases are skipped whole. So the quality string
//! is consumed at a different rate from the bases string: one quality per emitted element and none
//! for the markers, which is why the two indices are independent and why leftover qualities are
//! their own error.
//!
//! # The field count is not the format's
//!
//! ```java
//! if (tokens.length < MINIMUM_FIELDS || tokens.length > MAXIMUM_FIELDS) { ... }
//! ```
//!
//! `MINIMUM_FIELDS` is 4 and `MAXIMUM_FIELDS` is 6, where a samtools mpileup line has six columns
//! and may have more. A **seven**-column line is therefore refused. A **five**-column line is worse:
//! it passes this check and then reads `tokens[5]`, which is an `ArrayIndexOutOfBoundsException`
//! rather than the codec's own exception. Both are reproduced, because a caller that catches the
//! codec's error type would not catch the second.
//!
//! Two of the messages carry typos, "THe SAM pileup line" and "this codes is only valid". They are
//! kept: a port that fixes the spelling of a message changes what a caller matching on it sees.

/// `MINIMUM_FIELDS`.
pub const MINIMUM_FIELDS: usize = 4;
/// `MAXIMUM_FIELDS`.
pub const MAXIMUM_FIELDS: usize = 6;

/// `SAMPileupElement`: one base of one read at this locus, with its quality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SamPileupElement {
    pub base: u8,
    pub qual: u8,
}

/// `SAMPileupFeature`: everything one line describes.
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

    /// `getBasesString()`.
    pub fn bases_string(&self) -> String {
        self.elements.iter().map(|e| e.base as char).collect()
    }

    /// `getBaseQuals()`.
    pub fn base_quals(&self) -> Vec<u8> {
        self.elements.iter().map(|e| e.qual).collect()
    }
}

/// What the codec refuses, and how.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PileupCodecError {
    /// `CodecLineParsingException` with the codec's own message.
    Parsing(String),
    /// The raw `ArrayIndexOutOfBoundsException` a five-column line produces, with Java's own text.
    IndexOutOfBounds { index: usize, length: usize },
}

impl PileupCodecError {
    /// The message the exception carries.
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
pub fn decode(line: &str) -> Result<SamPileupFeature, PileupCodecError> {
    // `SPLIT_PATTERN.split(line.trim(), -1)`: tabs, keeping trailing empty fields.
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
    let bases = tokens.get(4).ok_or(PileupCodecError::IndexOutOfBounds {
        index: 4,
        length: tokens.len(),
    })?;
    let qualities = tokens.get(5).ok_or(PileupCodecError::IndexOutOfBounds {
        index: 5,
        length: tokens.len(),
    })?;

    let elements = parse_bases_and_quals(bases, qualities, reference_base)?;
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
                // increment at the end of the loop.
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

    // `i != bases.length()` cannot fail here, because the loop runs to the end; only leftover
    // qualities can.
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
pub fn can_decode(path: &str) -> bool {
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

/// `SAMUtils.fastqToPhred`.
fn fastq_to_phred(character: u8) -> u8 {
    character - 33
}

/// `parseInteger(token, parsedValue)`.
fn parse_integer(token: &str, parsed_value: &str) -> Result<i32, PileupCodecError> {
    token.parse::<i32>().map_err(|_| {
        PileupCodecError::Parsing(format!(
            "The SAM pileup line had unexpected {parsed_value}: {token}"
        ))
    })
}

/// `parseBase(String, parsedValue)`.
fn parse_base_token(token: &str, parsed_value: &str) -> Result<u8, PileupCodecError> {
    if token.len() != 1 {
        return Err(PileupCodecError::Parsing(format!(
            "The SAM pileup line had unexpected base at {parsed_value}: {token}"
        )));
    }
    parse_base(token.as_bytes()[0], parsed_value)
}

/// `parseBase(byte, parsedValue)`: an N folds to N, everything else goes through the base index and
/// comes back upper case.
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
