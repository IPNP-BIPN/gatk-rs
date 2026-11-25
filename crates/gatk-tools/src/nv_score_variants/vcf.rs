//! What pysam 0.22 and its htslib 1.18 do to a VCF that goes through them unchanged.
//!
//! `create_output_vcf.py` does not copy lines. It parses each record into BCF and formats it back,
//! so a float in INFO or FORMAT comes out with six significant digits, an integer loses its leading
//! zeros, a sample short of fields gets them as dots, and a tag the header does not define ends the
//! run. This is the part of `vcf.c` and `kstring.c` that a text VCF round trip goes through.
//!
//! The reference is GATK 4.7.0.0's Python. `create_output_vcf.py`, `models/wrapper.py`,
//! `encoders.py` and `nvscorevariants.py` are unchanged since 4.6.2.0, the image the port was
//! measured against (`docs/nvscorevariants.md`). The 4.7.0.0 tree pins `bioconda::pysam=0.22.0`
//! (`scripts/gatkcondaenv.yml.template`). No htslib source is in this checkout, so the htslib
//! functions each item models are named and not quoted. `NVScoreVariants.java` checks the Python
//! environment, extracts the model files, builds the script's command line and returns the
//! script's exit value; it opens neither the input VCF nor the output.
//!
//! # What is read and what is written
//!
//! The VCF half of the reference is `create_output_vcf.py`. Its loop, lines 23 to 37:
//!
//! ```python
//!     with open(scores_file) as scoredVariants:
//!         sv = next(scoredVariants)
//!         for variant in variant_file:
//!             scoredVariant = sv.split('\t')
//!             if variant.contig == scoredVariant[CONTIG_INDEX] and \
//!                variant.pos == int(scoredVariant[POS_INDEX]) and \
//!                variant.ref == scoredVariant[REF_INDEX] and \
//!                ', '.join(variant.alts or []) == re.sub('[\[\]]', '', scoredVariant[ALT_INDEX]):
//!
//!                     if len(scoredVariant) > KEY_INDEX:
//!                         variant.info.update({label: float(scoredVariant[KEY_INDEX])})
//!
//!                     vcfWriter.write(variant)
//!
//!                     sv = next(scoredVariants, None)
//! ```
//!
//! - **Read**: CHROM, POS, REF and ALT for the match above and for the encoders; INFO values by
//!   key, through `encoders.py:161`, `encoding[i] = variant.info.get(key, 0.0)` (which keys is
//!   `encode.rs`'s business); and `variant.stop` (`encoders.py:541`), which is [`Record::stop`]
//!   and so depends on an INFO `END`. FORMAT is never read by the tool.
//! - **Written by the tool**: one INFO value, `CNN_1D` or `CNN_2D` ([`Record::set_float`]), and one
//!   header line for it ([`Header::add_info`]).
//! - **Written by htslib on its own**: `END`, added, rewritten or dropped by `bcf_sync_end`
//!   (`Record::sync_end`), and every other field, which is parsed into BCF and printed back, not
//!   copied. That is where the bytes move, as the next four items say.
//! - Every Float, in INFO, FORMAT and QUAL, is held as a 32-bit float and printed by [`kputd`]
//!   with at most six significant digits.
//! - Every Integer is printed from its value, so `007` and `+7` come out as `7`, and a value out
//!   of range comes out as `.` (`parse_integers`).
//! - A sample with fewer fields than FORMAT names gets `.` for the rest (`render_sample`).
//! - A contig, FILTER, INFO or FORMAT tag the header does not define, or a symbolic allele that
//!   needs an `END` the header does not define, ends the run with `OSError` after the lines
//!   before it were written ([`Record::format_line`]). The docs record `Tool returned 1` and a
//!   truncated VCF on each of these.
//!
//! # The score's bytes
//!
//! The score goes through text twice. `models/wrapper.py` line 43 writes it to the temporary file:
//!
//! ```python
//!                 self.tmp_file.write('%s\t%d\t%s\t[%s]\t%.3f\n' % (batch['chrom'][i],
//! ```
//!
//! and line 33 above reads that text back with `float(...)`. pysam stores the Python `float` in
//! the record as a BCF Float, which is 32 bits, and htslib prints that with [`kputd`]. So what
//! the VCF shows is not the `%.3f` text: trailing zeros go (`1.500` is `1.5`, `0.000` is `0`,
//! `-0.000` is `-0`), and a score of 1000 or more keeps six significant digits, so `1234.568` in
//! the score file is `1234.57` in the VCF. `nan`, `inf` and `-inf` come out as themselves.
//!
//! # Header lines
//!
//! See [`Header::parse`] and [`Header::add_info`]: a `##fileformat=VCFv4.2` line is supplied when
//! the first line is not a `##fileformat`; `##FILTER=<ID=PASS,...>` with htslib's own description
//! is put after the input's first line, and a PASS line later in the input is dropped; any line
//! whose key and ID (or whole value) repeat an earlier one is dropped; the score's `##INFO` line
//! goes last, unless the input already defines that ID. The `#CHROM` line is rebuilt, not copied.
//!
//! The `HashMap`s in [`Header`] are lookups only. Output order comes from `Header::lines` and
//! from each record's `info` list, which keep the input's order and append what is new.
//!
//! # Reading this file without reading Rust
//!
//! - `Option<T>` is a value or nothing (`Some(x)` or `None`); `Result<T, String>` is a value or an
//!   error message (`Ok(x)` or `Err(why)`). `?` after an expression returns its `None` or `Err`
//!   from the enclosing function at once, Java's `throws` at one call site.
//! - `&str` is a borrowed, read-only view of text; `text[a..b]` is the bytes `a` to `b - 1` of it,
//!   and `&[T]` is the same for an array. `bytes.get(i)` is `None` past the end where `bytes[i]`
//!   would panic, which is how the parsers below read "the next character, if any".
//! - `match` must cover every case, and the compiler checks that it does. `matches!(x, A | B)` is
//!   a one-line `match` that answers true or false.
//! - `let Some(x) = y else { ... };` unwraps `y` or runs the block, which must leave the
//!   enclosing function or closure.
//! - `x as f32` from an `f64` rounds to nearest; `x as u64` from an `f64` truncates toward zero.
//! - `.iter().map(..).collect()` is Java's `stream().map(..).collect(..)`; `.any(..)` and
//!   `.find(..)` are `anyMatch` and `findFirst` with the predicate inline.
//! - A function's last expression, with no `;`, is its return value.

use std::collections::HashMap;

/// A header line's `Type`.
///
/// `String` and `Character`, and any `Type` this does not know, are all [`Kind::Text`]: the value
/// is kept as the text it was and printed back unchanged. (`#[derive(...)]` asks the compiler for
/// copying, printing and `==`, which Java would write by hand or get from a record.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Flag,
    Integer,
    Float,
    Text,
}

/// A header line's `Number`: a count, `A`, `R`, `G`, or `.`. A missing or unreadable `Number` is
/// [`Number::Variable`]. Only `encode.rs` reads it; nothing in this file's output depends on it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Number {
    Fixed(i64),
    PerAlternate,
    PerAllele,
    PerGenotype,
    Variable,
}

/// A tag's definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Definition {
    pub kind: Kind,
    pub number: Number,
}

/// `bcf_hdr_t`, as far as a text round trip uses it.
///
/// `Vec<T>` is Java's `ArrayList<T>` and `HashMap` is Java's `HashMap`, iteration order undefined
/// in both languages.
#[derive(Clone, Debug, Default)]
pub struct Header {
    /// The `##` lines in the order `bcf_hdr_format` prints them.
    pub lines: Vec<String>,
    /// The `#CHROM` line's columns after the ninth (`FORMAT`), in order.
    pub samples: Vec<String>,
    /// The `##INFO` IDs and their definitions, the first line for each ID winning.
    pub info: HashMap<String, Definition>,
    /// The `##FORMAT` IDs, as [`Header::info`].
    pub format: HashMap<String, Definition>,
    /// The `##FILTER` IDs, `PASS` always among them.
    pub filters: Vec<String>,
    /// The `##contig` IDs. A record on any other contig cannot be written.
    pub contigs: Vec<String>,
}

/// The value of `key` in a structured line `##X=<A=1,B="two, three",...>`, with its quotes taken
/// off, or `None` when the line has no `=<` or no such key.
///
/// A quoted value runs to its closing quote, so a comma inside `Description="..."` does not end
/// it; a backslash skips the character after it. `'a` says the text returned is a slice of `line`
/// and lives as long as it does, which is how Rust returns part of a string without copying it.
fn structured_field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let open = line.find("=<")?;
    let body = line[open + 2..]
        .strip_suffix('>')
        .unwrap_or(&line[open + 2..]);
    let bytes = body.as_bytes();
    let mut at = 0;
    while at < bytes.len() {
        let name_end = at + body[at..].find('=')?;
        let name = body[at..name_end].trim();
        let mut value_end = name_end + 1;
        let value_start = value_end;
        if bytes.get(value_end) == Some(&b'"') {
            value_end += 1;
            while value_end < bytes.len() && bytes[value_end] != b'"' {
                if bytes[value_end] == b'\\' {
                    value_end += 1;
                }
                value_end += 1;
            }
            value_end = (value_end + 1).min(bytes.len());
        } else {
            while value_end < bytes.len() && bytes[value_end] != b',' {
                value_end += 1;
            }
        }
        if name == key {
            return Some(body[value_start..value_end].trim_matches('"'));
        }
        at = value_end + 1;
    }
    None
}

// `impl Header { ... }` holds the methods of `Header`; Java would write them inside the class.
impl Header {
    /// `bcf_hdr_parse`: the `##` lines in file order, and the `#CHROM` line.
    ///
    /// **What comes out that did not go in**: a `##fileformat=VCFv4.2` line when the first line is
    /// not a `##fileformat`, and htslib's own PASS line,
    /// `##FILTER=<ID=PASS,Description="All filters passed">`, straight after the first line.
    /// Because it is added before the rest of the input and `Header::add_line` keeps the first of
    /// two lines with one ID, a PASS line anywhere but first in the input is dropped, whatever its
    /// description says. A file with no `##` lines at all gets the same two lines.
    pub fn parse(meta: &[&str], columns: &str) -> Header {
        // `..Header::default()` fills every field not named with its empty value.
        let mut header = Header {
            samples: columns.split('\t').skip(9).map(str::to_string).collect(),
            ..Header::default()
        };
        // The (key, ID) pairs already taken, for the duplicate test in `add_line`.
        let mut seen: Vec<(String, String)> = Vec::new();
        let mut first = true;
        for line in meta {
            if first {
                first = false;
                if line.starts_with("##fileformat=") {
                    header.lines.push((*line).to_string());
                } else {
                    header.lines.push("##fileformat=VCFv4.2".to_string());
                    header.add_line(line, &mut seen);
                }
                // The filter every header has, whatever it says about it later.
                header.add_line(
                    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
                    &mut seen,
                );
                continue;
            }
            header.add_line(line, &mut seen);
        }
        // Still `first`: there were no `##` lines at all.
        if first {
            header.lines.push("##fileformat=VCFv4.2".to_string());
            header.add_line(
                "##FILTER=<ID=PASS,Description=\"All filters passed\">",
                &mut seen,
            );
        }
        header
    }

    /// `bcf_hdr_add_hrec`, which ignores what it already has.
    ///
    /// **What**: append one `##key=value` line, and record what it defines.
    ///
    /// **How**: a structured `FILTER`, `INFO`, `FORMAT` or `contig` line is identified by its
    /// `ID`; any other line by its whole value. A line whose key and identity were seen before is
    /// dropped, so the **first** definition of an ID wins, both for what is printed and for the
    /// `Type` the records are parsed with. A line with no `=` is dropped.
    fn add_line(&mut self, line: &str, seen: &mut Vec<(String, String)>) {
        let Some(equals) = line.find('=') else { return };
        // `line[2..equals]`: the text between the leading `##` and the first `=`.
        let key = line[2..equals].to_string();
        let structured = line[equals + 1..].starts_with('<');
        // A `match` arm with `if` (a guard) only matches when the condition also holds.
        let identity = match key.as_str() {
            "FILTER" | "INFO" | "FORMAT" | "contig" if structured => {
                structured_field(line, "ID").unwrap_or_default().to_string()
            }
            _ => line[equals + 1..].to_string(),
        };
        // `|(k, i)| ...` is a closure, Java's `(k, i) -> ...`; `*k` reads through the reference.
        if seen.iter().any(|(k, i)| *k == key && *i == identity) {
            return;
        }
        match key.as_str() {
            "FILTER" if structured => self.filters.push(identity.clone()),
            "contig" if structured => self.contigs.push(identity.clone()),
            "INFO" | "FORMAT" if structured => {
                let kind = match structured_field(line, "Type") {
                    Some("Integer") => Kind::Integer,
                    Some("Float") => Kind::Float,
                    Some("Flag") => Kind::Flag,
                    _ => Kind::Text,
                };
                let number = match structured_field(line, "Number") {
                    Some("A") => Number::PerAlternate,
                    Some("R") => Number::PerAllele,
                    Some("G") => Number::PerGenotype,
                    Some(".") | None => Number::Variable,
                    // A count that is not a number is `.`; `map_or(a, f)` is `a` on failure.
                    Some(text) => text.parse().map_or(Number::Variable, Number::Fixed),
                };
                let definition = Definition { kind, number };
                if key == "INFO" {
                    self.info.insert(identity.clone(), definition);
                } else {
                    self.format.insert(identity.clone(), definition);
                }
            }
            _ => {}
        }
        seen.push((key, identity));
        self.lines.push(line.to_string());
    }

    /// `header.info.add`, which leaves a line it already has where it is.
    ///
    /// `create_output_vcf.py` lines 18 and 19 (line 18 is one line in the reference, wrapped here
    /// before `description=`; line 19 is verbatim):
    ///
    /// ```python
    ///     variant_file.header.info.add(id=label, number=1, type='Float',
    ///                                  description='Log odds of being a true variant versus \
    ///                     being false under the trained Convolutional Neural Network')
    /// ```
    ///
    /// The backslash continues the string literal, so the space before it and the twenty that open
    /// line 19 are all in the description (`SCORE_DESCRIPTION` in `mod.rs` carries the 21). The
    /// line is appended after every input line, so it is the last `##` line of the output. It is
    /// added to the reader's header before `header.copy()` (line 20) makes the writer's, which is
    /// why the writer knows it.
    ///
    /// When the input already defines the ID, its line and its `Type` stay as they were, and
    /// [`Record::set_float`] still writes a Float.
    pub fn add_info(&mut self, id: &str, description: &str) {
        if self.info.contains_key(id) {
            return;
        }
        self.info.insert(
            id.to_string(),
            Definition {
                kind: Kind::Float,
                number: Number::Fixed(1),
            },
        );
        self.lines.push(format!(
            "##INFO=<ID={id},Number=1,Type=Float,Description=\"{description}\">"
        ));
    }

    /// `bcf_hdr_format`: the `##` lines, then a `#CHROM` line built from the sample names.
    ///
    /// The `#CHROM` line is not the input's. `FORMAT` is printed only when there are samples, so a
    /// sites-only input whose `#CHROM` line ends in `FORMAT` loses that column.
    pub fn format_text(&self) -> String {
        let mut text = String::new();
        for line in &self.lines {
            text.push_str(line);
            text.push('\n');
        }
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            text.push_str("\tFORMAT");
            for sample in &self.samples {
                text.push('\t');
                text.push_str(sample);
            }
        }
        text.push('\n');
        text
    }
}

/// `hts_str2dbl`: short decimals by integer arithmetic, and `strtod` for what does not fit.
///
/// Returns the value and how much of the text it read, or `None` where nothing could be read.
/// (`(f64, usize)` is a pair; Java would return a small object.) This is how htslib reads every
/// Float in INFO and FORMAT.
///
/// **How**: skip spaces and one sign; skip leading zeros; read the digits, then a point and more
/// digits, into one integer `n`; divide by ten to the number of digits after the point. A budget
/// of 15 is spent on **every check**, the one that finds a non-digit included, so the integer
/// path takes at most 13 digits without a point and 12 with one (leading zeros before the point
/// do not count, zeros after it do). Text that spends the budget, has an exponent, is hex
/// (`0x...`) or does not start with a digit (`.5`, `inf`, `nan`) goes to `strtod` instead.
///
/// **Why both paths give one answer**: `n` is below 2^53 and every divisor used is a power of ten
/// no larger than `1e12`, so both are exact doubles and the one division rounds once, to the
/// double nearest the decimal, which is what a correct `strtod` returns too. Where the split
/// matters is hex, which goes to `strtod`, and which that function here does not model.
pub fn hts_str2dbl(text: &str) -> Option<(f64, usize)> {
    let bytes = text.as_bytes();
    let mut v = 0;
    while bytes.get(v).is_some_and(u8::is_ascii_whitespace) {
        v += 1;
    }
    let negative = bytes.get(v) == Some(&b'-');
    if negative || bytes.get(v) == Some(&b'+') {
        v += 1;
    }
    match bytes.get(v) {
        Some(b'1'..=b'9') => {}
        Some(b'0') if !matches!(bytes.get(v + 1), Some(b'x' | b'X')) => {}
        _ => return strtod(text),
    }
    while bytes.get(v) == Some(&b'0') {
        v += 1;
    }
    // `start` is the first significant digit, or the point, after the leading zeros.
    let start = v;
    // `15i32`: the literal 15 as a 32-bit `int`.
    let mut budget = 15i32;
    let mut n = 0u64;
    let mut point: Option<usize> = None;
    loop {
        budget -= 1;
        if budget == 0 || !bytes.get(v).is_some_and(u8::is_ascii_digit) {
            break;
        }
        n = n * 10 + u64::from(bytes[v] - b'0');
        v += 1;
    }
    if budget != 0 && bytes.get(v) == Some(&b'.') {
        point = Some(v - start);
        v += 1;
        loop {
            budget -= 1;
            if budget == 0 || !bytes.get(v).is_some_and(u8::is_ascii_digit) {
                break;
            }
            n = n * 10 + u64::from(bytes[v] - b'0');
            v += 1;
        }
    }
    // No point: as if it were just past the last digit, so the divisor below is one.
    let point = point.unwrap_or(v - start);
    if budget == 0 || matches!(bytes.get(v), Some(b'e' | b'E')) {
        return strtod(text);
    }
    const POWERS: [f64; 22] = [
        1.0, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15,
        1e16, 1e17, 1e18, 1e19, 1e20,
    ];
    // `v - start - point` counts the point itself as well as the digits after it, so the table
    // is shifted by one: slot 0 is the no-point case, slot 1 a point with no digit after it
    // (`5.`), and slot k is 10^(k - 1), which is why the table starts with two ones.
    // `n as f64` is exact here (see the item comment).
    let value = n as f64 / POWERS[v - start - point];
    Some((if negative { -value } else { value }, v))
}

/// glibc's `strtod` on the decimal forms, `inf` and `nan`, which it rounds correctly as `parse`
/// does.
///
/// **How**: find the longest prefix `strtod` would take (a sign, digits with at most one point, and
/// an exponent only when a digit follows the `e`), and hand that prefix to Rust's own correctly
/// rounded `parse`. `inf`, `infinity` and `nan` are matched without regard to case.
///
/// Not modelled: glibc's hexadecimal floats. `0x10` reads as `0`, with one character taken, where
/// glibc reads 16. Every caller here drops the length, so only the value matters.
fn strtod(text: &str) -> Option<(f64, usize)> {
    let bytes = text.as_bytes();
    let mut at = 0;
    while bytes.get(at).is_some_and(u8::is_ascii_whitespace) {
        at += 1;
    }
    let number_start = at;
    if matches!(bytes.get(at), Some(b'+' | b'-')) {
        at += 1;
    }
    let lower = text[at..].to_ascii_lowercase();
    for (word, value) in [
        ("infinity", f64::INFINITY),
        ("inf", f64::INFINITY),
        ("nan", f64::NAN),
    ] {
        if lower.starts_with(word) {
            let signed = if bytes[number_start] == b'-' {
                -value
            } else {
                value
            };
            return Some((signed, at + word.len()));
        }
    }
    // `is_some_and(f)`: there is a next byte and `f` holds for it.
    let digits_start = at;
    while bytes.get(at).is_some_and(u8::is_ascii_digit) {
        at += 1;
    }
    let mut digits = at - digits_start;
    if bytes.get(at) == Some(&b'.') {
        at += 1;
        let fraction = at;
        while bytes.get(at).is_some_and(u8::is_ascii_digit) {
            at += 1;
        }
        digits += at - fraction;
    }
    if digits == 0 {
        return None;
    }
    if matches!(bytes.get(at), Some(b'e' | b'E')) {
        let mut exponent = at + 1;
        if matches!(bytes.get(exponent), Some(b'+' | b'-')) {
            exponent += 1;
        }
        if bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
            while bytes.get(exponent).is_some_and(u8::is_ascii_digit) {
                exponent += 1;
            }
            at = exponent;
        }
    }
    // `.ok()` turns a parse error into `None`; Rust's `parse` accepts `1.` and `.5` as `strtod`
    // does, and gives infinity on overflow as `strtod` does.
    text[number_start..at].parse().ok().map(|value| (value, at))
}

/// `hts_str2int` over a C `long`, and how much of the text it read.
///
/// Returns the value, whether it overflowed a 64-bit `long`, and the length read; `None` when
/// there is no digit after the optional sign. Unlike [`hts_str2dbl`] it skips no spaces. The
/// digits are gathered unsigned and the sign applied last, so `-9223372036854775808` (Java's
/// `Long.MIN_VALUE`) is in range and one more is not.
fn hts_str2int(text: &str) -> Option<(i64, bool, usize)> {
    let bytes = text.as_bytes();
    let mut at = 0;
    let negative = bytes.first() == Some(&b'-');
    if negative || bytes.first() == Some(&b'+') {
        at = 1;
    }
    let digits = at;
    let mut value: u64 = 0;
    let mut overflow = false;
    // `while let Some(x) = ...` loops for as long as the expression is `Some`. `checked_mul` and
    // `checked_add` give `None` on overflow where Java's `long` arithmetic would wrap silently.
    while let Some(digit) = bytes.get(at).filter(|b| b.is_ascii_digit()) {
        match value
            .checked_mul(10)
            .and_then(|v| v.checked_add(u64::from(digit - b'0')))
        {
            Some(next) => value = next,
            None => overflow = true,
        }
        at += 1;
    }
    if at == digits {
        return None;
    }
    let limit = if negative {
        1u64 << 63
    } else {
        (1u64 << 63) - 1
    };
    if value > limit {
        overflow = true;
    }
    // Two's complement: negating `2^63` as unsigned and reading the bits as `i64` gives
    // `i64::MIN`, which a signed negation could not reach. `as` between integer types keeps the
    // bits, never checks the range.
    let signed = if negative {
        value.wrapping_neg() as i64
    } else {
        value as i64
    };
    Some((signed, overflow, at))
}

/// C's `printf("%g", d)` with the default precision of six, for the values [`kputd`] hands it.
///
/// **How**: C's own rule. Write the value in scientific form with six significant digits; if its
/// exponent is below -4 or at least 6, print that form, otherwise print fixed-point with
/// `5 - exponent` decimals; then strip trailing zeros and a trailing point. Rust's `{:.5e}` and
/// `{:.N}` are both correctly rounded from the exact binary value, as glibc's `printf` is.
///
/// The exponent is rewritten because Rust prints `1.23457e6` where C prints `1.23457e+06`: a sign
/// always, and at least two digits. The sign of a NaN is not printed; [`kputd`] only sends this
/// non-negative numbers, infinity and NaN.
fn percent_g(d: f64) -> String {
    if d.is_nan() {
        return "nan".to_string();
    }
    if d.is_infinite() {
        return "inf".to_string();
    }
    let scientific = format!("{d:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    if !(-4..6).contains(&exponent) {
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        let sign = if exponent < 0 { '-' } else { '+' };
        return format!("{mantissa}e{sign}{:02}", exponent.abs());
    }
    // `{:.*}` takes the number of decimals from the argument before the value.
    let fixed = format!("{d:.*}", (5 - exponent) as usize);
    if fixed.contains('.') {
        fixed
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    } else {
        fixed
    }
}

/// `kputd`, which is not `%g` between 0.0001 and 999999: it scales by 1e10, in a double, and rounds
/// by adding half of the last digit it keeps.
///
/// **What**: how htslib prints every Float of a record (QUAL, INFO, FORMAT) and so the score:
/// six significant digits, trailing zeros and a trailing point removed. The callers pass a
/// 32-bit float widened to `f64` (`f64::from`, which is exact), so `12.345` arrives as
/// `12.345000267028809` and the six digits are what hide the 32-bit error.
///
/// **How**, between 0.0001 and 999999: `i = (d * 1e10) as u64`, so `i` holds the value with ten
/// decimals, truncated. Adding 5 at the position just past the sixth significant digit, then
/// keeping six digits, rounds half up on that truncated double. Outside the range, `percent_g`.
///
/// **Why not `%g` everywhere**: the two disagree on ties. `590384.5` (exact in 32 bits too) prints
/// `590385` here and `590384` under `%g`, which rounds half to even; and the product `d * 1e10`
/// is itself rounded, so a value just below a tie can land on it. Zero is `0` and negative zero
/// is `-0`.
pub fn kputd(value: f64) -> String {
    if value == 0.0 {
        return if value.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    let mut out = String::new();
    let mut d = value;
    if d < 0.0 {
        out.push('-');
        d = -d;
    }
    if !(0.0001..=999999.0).contains(&d) {
        out.push_str(&percent_g(d));
        return out;
    }
    // Truncated toward zero. `d` is at most 999999 here, so `d * 1e10` is below 2^64 and the
    // cast never saturates.
    let mut i = (d * 1e10) as u64;
    i += if d < 0.001 {
        5
    } else if d < 0.01 {
        50
    } else if d < 0.1 {
        500
    } else if d < 1.0 {
        5000
    } else if d < 10.0 {
        50000
    } else if d < 100.0 {
        500000
    } else if d < 1000.0 {
        5000000
    } else if d < 10000.0 {
        50000000
    } else if d < 100000.0 {
        500000000
    } else {
        5000000000
    };
    // The decimal digits of `i`, as bytes. Ten of them are below the point.
    let digits = i.to_string().into_bytes();
    let count = digits.len();
    let mut text: Vec<u8> = Vec::new();
    if count <= 10 {
        // Below one: six digits from the first that is not a zero, once the zeros go in front.
        text.extend_from_slice(b"0.");
        text.extend(std::iter::repeat_n(b'0', 10 - count));
        let kept = 6.min(count);
        text.extend_from_slice(&digits[..kept]);
    } else {
        // Seven characters, the point among them, and a point left last goes.
        let whole = count - 10;
        text.extend_from_slice(&digits[..whole]);
        text.push(b'.');
        text.extend_from_slice(&digits[whole..]);
        text.truncate(7);
        if text.last() == Some(&b'.') {
            text.pop();
        }
    }
    // Trailing zeros only after a point: `100` keeps its zeros, `1.50` loses one.
    if text.contains(&b'.') {
        while text.last() == Some(&b'0') {
            text.pop();
        }
        if text.last() == Some(&b'.') {
            text.pop();
        }
    }
    out.push_str(&String::from_utf8_lossy(&text));
    out
}

/// A value of an INFO key, as encoded in BCF.
///
/// A Rust `enum` variant can carry data, so this is one type with four shapes, which Java would
/// write as a sealed interface with four records. `None` in a list is a `.` in the text.
#[derive(Clone, Debug, PartialEq)]
pub enum InfoValue {
    /// A key with no value, which is how a Flag looks, and a text value that is empty.
    Nothing,
    /// A String, Character or Flag value, or the value of a key the header does not define: kept
    /// as the text it was.
    Text(String),
    /// An Integer list. BCF holds 32-bit integers, hence `i32`.
    Integers(Vec<Option<i32>>),
    /// A Float list. BCF holds 32-bit floats, hence `f32`: Java's `float`, not `double`.
    Floats(Vec<Option<f32>>),
}

/// One INFO entry, and whether the header defines it.
///
/// An undefined key is kept, as its text, until the record is written: [`Record::format_line`]
/// is where it stops the run.
#[derive(Clone, Debug, PartialEq)]
pub struct InfoEntry {
    pub key: String,
    pub value: InfoValue,
    pub defined: bool,
}

/// `bcf1_t` after `vcf_parse`: pysam's `VariantRecord`.
///
/// The fields are what [`Record::format_line`] prints back. FORMAT values are rendered once, at
/// parse time (`render_sample`), and nothing in the tool reads them.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub contig: String,
    /// 1-based, as `variant.pos`.
    pub pos: i64,
    pub id: String,
    pub reference: String,
    pub alternates: Vec<String>,
    /// `None` is a `.`; a number is held in 32 bits, as BCF holds it.
    pub quality: Option<f32>,
    pub filters: Vec<String>,
    /// In the order of the line, followed by what [`Record::set_float`] and
    /// `Record::sync_end` add.
    pub info: Vec<InfoEntry>,
    /// The FORMAT keys; empty for a sites-only file or a FORMAT of `.`.
    pub format: Vec<String>,
    /// For each sample, its fields as rendered.
    pub samples: Vec<Vec<String>>,
    /// `rlen`: the length of REF unless an INFO/END says otherwise.
    pub reference_length: i64,
    /// A contig, filter or FORMAT tag the header does not define.
    pub undefined: Option<String>,
}

/// The lowest integer an INFO or FORMAT Integer may hold here; anything below is read as `.`.
/// It is the port's value: no htslib source of the reference's version is in this checkout to
/// quote it from.
const INT32_MISSING_BELOW: i64 = -2147483392;
/// Java's `Integer.MAX_VALUE`.
const INT32_MAX: i64 = 2147483647;

/// The comma-separated values of an Integer field, as BCF would hold them.
///
/// A value that is `.`, has no digits, overflows a `long`, or falls outside
/// `INT32_MISSING_BELOW..=INT32_MAX` becomes `None`, so the output shows `.` where the input had
/// a number. Text after the digits is ignored: `12abc` is 12.
fn parse_integers(text: &str) -> Vec<Option<i32>> {
    text.split(',')
        .map(|one| match hts_str2int(one) {
            Some((value, false, _)) if (INT32_MISSING_BELOW..=INT32_MAX).contains(&value) => {
                Some(value as i32)
            }
            _ => None,
        })
        .collect()
}

/// The comma-separated values of a Float field, as BCF would hold them: [`hts_str2dbl`], then
/// rounded to 32 bits. `.` and other text with no number become `None`.
fn parse_floats(text: &str) -> Vec<Option<f32>> {
    text.split(',')
        .map(|one| hts_str2dbl(one).map(|(value, _)| value as f32))
        .collect()
}

/// A Float list as htslib prints it: [`kputd`] for each value, `.` for each missing one.
fn render_floats(values: &[Option<f32>]) -> String {
    values
        .iter()
        .map(|v| v.map_or(".".to_string(), |v| kputd(f64::from(v))))
        .collect::<Vec<_>>()
        .join(",")
}

/// An Integer list as htslib prints it: plain decimal, `.` for each missing one.
fn render_integers(values: &[Option<i32>]) -> String {
    values
        .iter()
        .map(|v| v.map_or(".".to_string(), |v| v.to_string()))
        .collect::<Vec<_>>()
        .join(",")
}

impl Record {
    /// `vcf_parse`. Why a record cannot be parsed is what the reference's `OSError` would be about.
    ///
    /// **What**: one tab-separated data line, read against `header`. ID, REF, ALT, FILTER and
    /// text values are kept as written; numbers are parsed into what BCF would hold, so they print
    /// back differently ([`kputd`], `parse_integers`).
    ///
    /// **Not an error here, an error later**: a contig, FILTER or FORMAT key the header does not
    /// define is noted in `undefined` (the first one found, in that order) and an undefined INFO
    /// key in its entry. The record still parses, because the reader reads it and the encoders
    /// score it; it is `VariantFile.write` that refuses it ([`Record::format_line`]).
    pub fn parse(header: &Header, line: &str) -> Result<Record, String> {
        let mut fields = line.split('\t');
        // A closure that takes the next column or names the missing one as an error. It is
        // `mut` because each call advances `fields`, which it borrows.
        let mut next = |name: &str| {
            fields
                .next()
                .ok_or_else(|| format!("no {name} in the VCF line: {line}"))
        };
        let contig = next("CHROM")?.to_string();
        let position = next("POS")?;
        let pos: i64 = position
            .parse()
            .map_err(|_| format!("Could not parse the position '{position}'"))?;
        let id = next("ID")?.to_string();
        let reference = next("REF")?.to_string();
        let alternates: Vec<String> = match next("ALT")? {
            "." => Vec::new(),
            text => text.split(',').map(str::to_string).collect(),
        };
        let quality = match next("QUAL")? {
            "." => None,
            // `atof`, into the record's float.
            text => Some(strtod(text).map_or(0.0, |(value, _)| value) as f32),
        };
        let filters: Vec<String> = match next("FILTER")? {
            "." => Vec::new(),
            text => text
                .split(';')
                .filter(|f| !f.is_empty())
                .map(str::to_string)
                .collect(),
        };
        let mut undefined = None;
        if !header.contigs.contains(&contig) {
            undefined = Some(format!("contig {contig}"));
        }
        for filter in &filters {
            if !header.filters.contains(filter) && undefined.is_none() {
                undefined = Some(format!("FILTER {filter}"));
            }
        }
        // `rlen`. `.len()` of a `String` is its length in bytes.
        let mut reference_length = reference.len() as i64;
        let mut info = Vec::new();
        for entry in next("INFO")?.split(';') {
            if entry.is_empty() || entry == "." {
                continue;
            }
            let (key, text) = match entry.split_once('=') {
                Some((key, text)) => (key, Some(text)),
                None => (entry, None),
            };
            // How a value is held follows the header's `Type`, not the text: `AF=0.50` under a
            // Float definition prints back as `AF=0.5`. A key the header does not define is kept
            // as text. The `match` on a pair tries the arms top to bottom.
            let definition = header.info.get(key);
            let value = match (text, definition.map(|d| d.kind)) {
                (None, _) => InfoValue::Nothing,
                (Some(""), None | Some(Kind::Flag | Kind::Text)) => InfoValue::Nothing,
                (Some(text), None | Some(Kind::Flag | Kind::Text)) => {
                    InfoValue::Text(text.to_string())
                }
                (Some(text), Some(Kind::Integer)) => InfoValue::Integers(parse_integers(text)),
                (Some(text), Some(Kind::Float)) => InfoValue::Floats(parse_floats(text)),
            };
            // An `END` changes `rlen`, and so `variant.stop` and the reads fetched for the 2D
            // model, only when the header makes it an Integer, it is one value, and it is at or
            // after POS. `if let [Some(end)] = values[..]` matches a list of exactly one present
            // value.
            if key == "END" {
                if let InfoValue::Integers(values) = &value {
                    if let [Some(end)] = values[..] {
                        if i64::from(end) > pos - 1 {
                            reference_length = i64::from(end) - (pos - 1);
                        }
                    }
                }
            }
            info.push(InfoEntry {
                key: key.to_string(),
                value,
                defined: definition.is_some(),
            });
        }
        // FORMAT and the samples. The header's sample list decides how many sample columns
        // there are: extra columns on the line are ignored, a missing one reads as `.`, and with
        // no samples in the header the FORMAT column is ignored altogether.
        let mut format = Vec::new();
        let mut samples = Vec::new();
        if let Ok(layout) = next("FORMAT") {
            if !header.samples.is_empty() && layout != "." {
                format = layout.split(':').map(str::to_string).collect();
                for key in &format {
                    if !header.format.contains_key(key) && undefined.is_none() {
                        undefined = Some(format!("FORMAT {key}"));
                    }
                }
            }
            for _ in 0..header.samples.len() {
                let text = next("sample").unwrap_or(".");
                samples.push(render_sample(header, &format, text));
            }
        } else if !header.samples.is_empty() {
            // No FORMAT column at all: each sample is empty, and prints as `.`.
            samples = vec![Vec::new(); header.samples.len()];
        }
        Ok(Record {
            contig,
            pos,
            id,
            reference,
            alternates,
            quality,
            filters,
            info,
            format,
            samples,
            reference_length,
            undefined,
        })
    }

    /// `variant.stop`: the 0-based end, exclusive, which is the 1-based end, inclusive.
    /// `encoders.py:541`, `stop = variant.stop`, bounds the reads the 2D model fetches.
    pub fn stop(&self) -> i64 {
        self.pos - 1 + self.reference_length
    }

    /// `has_symbolic_allele`: some ALT is `<...>`, such as `<DEL>` or `<NON_REF>`. A breakend is
    /// not symbolic by this test.
    fn has_symbolic_allele(&self) -> bool {
        self.alternates
            .iter()
            .any(|a| a.starts_with('<') && a.ends_with('>'))
    }

    /// `variant.info.update({key: value})` for a Float of `Number=1`.
    ///
    /// `create_output_vcf.py` line 33:
    /// `variant.info.update({label: float(scoredVariant[KEY_INDEX])})`. An INFO key already on
    /// the record keeps its place and takes the new value; a new one goes last, so the score is
    /// normally the last INFO entry (before an `END` that `Record::sync_end` adds). The value is
    /// held in 32 bits, which is why the VCF shows six digits of it and not the `%.3f` text
    /// (module header, "The score's bytes"). `defined` is set because [`Header::add_info`] has
    /// put the ID in the header.
    pub fn set_float(&mut self, key: &str, value: f32) {
        let replacement = InfoValue::Floats(vec![Some(value)]);
        if let Some(entry) = self.info.iter_mut().find(|e| e.key == key) {
            entry.value = replacement;
            entry.defined = true;
        } else {
            self.info.push(InfoEntry {
                key: key.to_string(),
                value: replacement,
                defined: true,
            });
        }
    }

    /// `bcf_sync_end`, which `VariantFile.write` does first. `Err` is the tag its `END` would need,
    /// added to the record's header and not to the writer's copy of it.
    ///
    /// **What it does to the bytes**: a record with no symbolic ALT and a span equal to its REF
    /// loses its `END`, even one that was correct, when the header defines `END`. Any other
    /// record gets `END` set to [`Record::stop`], appended if it was not there, or the run ends
    /// when the header does not define `END`.
    fn sync_end(&mut self, header: &Header) -> Result<(), String> {
        let defined = header.info.contains_key("END");
        if !self.has_symbolic_allele()
            && (self.reference.is_empty() || self.reference_length == self.reference.len() as i64)
        {
            if defined {
                self.info.retain(|entry| entry.key != "END");
            }
            return Ok(());
        }
        if !defined {
            return Err("END".to_string());
        }
        // `as i32` from `i64` keeps the low 32 bits; a stop past `Integer.MAX_VALUE` would wrap.
        let end = InfoValue::Integers(vec![Some((self.stop()) as i32)]);
        if let Some(entry) = self.info.iter_mut().find(|e| e.key == "END") {
            entry.value = end;
        } else {
            self.info.push(InfoEntry {
                key: "END".to_string(),
                value: end,
                defined: true,
            });
        }
        Ok(())
    }

    /// `VariantFile.write`: the line, or which tag ended the run. A tag the reader's header took in
    /// while reading this record was never in the copy the writer has.
    ///
    /// `create_output_vcf.py` line 35, `vcfWriter.write(variant)`, with the writer opened on
    /// line 21 from a copy of the header taken on line 20:
    ///
    /// ```python
    ///     header = variant_file.header.copy()
    ///     vcfWriter = VariantFile(vcf_out, 'w', header=header)
    /// ```
    ///
    /// **The checks, in the order they fire**: the `END` that `Record::sync_end` needs; then the
    /// first undefined contig, FILTER or FORMAT key from [`Record::parse`]; then the first INFO key
    /// the header does not define. Each `Err` becomes the reference's `OSError`, and the caller
    /// stops there, after the lines already written.
    ///
    /// **The line**: the columns in VCF order; `.` for no ALT, no QUAL, no FILTER and no INFO;
    /// QUAL and every Float through [`kputd`]; a Flag, and a text key given an empty value
    /// (`X=`), as the bare key; the FORMAT column and the samples only when the header has
    /// samples, and `.` for each sample when there are no FORMAT keys.
    pub fn format_line(&mut self, header: &Header) -> Result<String, String> {
        // `&mut self`: this changes the record, because `bcf_sync_end` does. `map_err` rewrites
        // the error into a message before `?` returns it.
        self.sync_end(header).map_err(|tag| {
            format!(
                "Invalid BCF, the INFO tag {tag} is too large at {}:{}",
                self.contig, self.pos
            )
        })?;
        if let Some(what) = &self.undefined {
            return Err(format!("Invalid BCF, {what} not present in the header"));
        }
        if let Some(entry) = self.info.iter().find(|entry| !entry.defined) {
            return Err(format!(
                "Invalid BCF, the INFO tag {} is too large at {}:{}",
                entry.key, self.contig, self.pos
            ));
        }
        let mut line = format!(
            "{}\t{}\t{}\t{}\t",
            self.contig, self.pos, self.id, self.reference
        );
        if self.alternates.is_empty() {
            line.push('.');
        } else {
            line.push_str(&self.alternates.join(","));
        }
        line.push('\t');
        match self.quality {
            None => line.push('.'),
            Some(quality) => line.push_str(&kputd(f64::from(quality))),
        }
        line.push('\t');
        if self.filters.is_empty() {
            line.push('.');
        } else {
            line.push_str(&self.filters.join(";"));
        }
        line.push('\t');
        if self.info.is_empty() {
            line.push('.');
        } else {
            for (at, entry) in self.info.iter().enumerate() {
                if at > 0 {
                    line.push(';');
                }
                line.push_str(&entry.key);
                match &entry.value {
                    InfoValue::Nothing => {}
                    InfoValue::Text(text) => {
                        line.push('=');
                        line.push_str(text);
                    }
                    InfoValue::Integers(values) => {
                        line.push('=');
                        line.push_str(&render_integers(values));
                    }
                    InfoValue::Floats(values) => {
                        line.push('=');
                        line.push_str(&render_floats(values));
                    }
                }
            }
        }
        if !self.samples.is_empty() {
            line.push('\t');
            if self.format.is_empty() {
                line.push('.');
                for _ in 0..self.samples.len() {
                    line.push_str("\t.");
                }
            } else {
                line.push_str(&self.format.join(":"));
                for sample in &self.samples {
                    line.push('\t');
                    line.push_str(&sample.join(":"));
                }
            }
        }
        line.push('\n');
        Ok(line)
    }
}

/// `vcf_parse_format` and what `bcf_fmt_array` prints of it, for one sample.
///
/// One string per FORMAT key, in FORMAT's order: `GT` as written; an Integer or Float key parsed
/// and printed back as in INFO (`render_integers`, `render_floats`); anything else as
/// written. A key with no value left in the sample is `.`, and values past the last key are
/// dropped. This runs at parse time and the strings are what [`Record::format_line`] joins.
fn render_sample(header: &Header, format: &[String], text: &str) -> Vec<String> {
    let mut values = text.split(':');
    format
        .iter()
        .map(|key| {
            // `return` inside a closure leaves the closure only: this key's string is `.`.
            let Some(value) = values.next() else {
                return ".".to_string();
            };
            // `_ if ...` matches anything when the guard holds, so `GT` is copied whatever the
            // header says its `Type` is.
            match header.format.get(key).map(|d| d.kind) {
                _ if key == "GT" => value.to_string(),
                Some(Kind::Integer) => render_integers(&parse_integers(value)),
                Some(Kind::Float) => render_floats(&parse_floats(value)),
                _ => value.to_string(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kputd_keeps_six_significant_digits_and_no_more() {
        for (value, text) in [
            (0.0, "0"),
            (-0.0, "-0"),
            (1.5, "1.5"),
            (100.0, "100"),
            (12.345000267028809, "12.345"),
            (577.4400024414062, "577.44"),
            (1234567.875, "1.23457e+06"),
            (0.00005000000237487257, "5e-05"),
            (123456.0, "123456"),
            (-16.118000030517578, "-16.118"),
            (0.03200000151991844, "0.032"),
        ] {
            assert_eq!(kputd(value), text, "{value}");
        }
    }

    #[test]
    fn the_fast_parser_and_strtod_divide_where_htslib_does() {
        assert_eq!(hts_str2dbl("12.5"), Some((12.5, 4)));
        assert_eq!(hts_str2dbl("-0.25"), Some((-0.25, 5)));
        assert_eq!(hts_str2dbl("1.5e-7"), Some((1.5e-7, 6)));
        assert_eq!(hts_str2dbl("7abc"), Some((7.0, 1)));
        assert_eq!(hts_str2dbl("abc"), None);
        assert_eq!(hts_str2dbl("."), None);
    }
}
