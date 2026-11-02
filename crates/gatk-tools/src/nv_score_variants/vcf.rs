//! What pysam 0.22 and its htslib 1.18 do to a VCF that goes through them unchanged.
//!
//! `create_output_vcf.py` does not copy lines. It parses each record into BCF and formats it back,
//! so a float in INFO or FORMAT comes out with six significant digits, an integer loses its leading
//! zeros, a sample short of fields gets them as dots, and a tag the header does not define ends the
//! run. This is the part of `vcf.c` and `kstring.c` that a text VCF round trip goes through.

use std::collections::HashMap;

/// A header line's `Type`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Flag,
    Integer,
    Float,
    Text,
}

/// A header line's `Number`.
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
#[derive(Clone, Debug, Default)]
pub struct Header {
    /// The `##` lines in the order `bcf_hdr_format` prints them.
    pub lines: Vec<String>,
    pub samples: Vec<String>,
    pub info: HashMap<String, Definition>,
    pub format: HashMap<String, Definition>,
    pub filters: Vec<String>,
    pub contigs: Vec<String>,
}

/// One `key=value` of a structured line, with the quotes taken off.
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

impl Header {
    /// `bcf_hdr_parse`.
    pub fn parse(meta: &[&str], columns: &str) -> Header {
        let mut header = Header {
            samples: columns.split('\t').skip(9).map(str::to_string).collect(),
            ..Header::default()
        };
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
    fn add_line(&mut self, line: &str, seen: &mut Vec<(String, String)>) {
        let Some(equals) = line.find('=') else { return };
        let key = line[2..equals].to_string();
        let structured = line[equals + 1..].starts_with('<');
        let identity = match key.as_str() {
            "FILTER" | "INFO" | "FORMAT" | "contig" if structured => {
                structured_field(line, "ID").unwrap_or_default().to_string()
            }
            _ => line[equals + 1..].to_string(),
        };
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

    /// `bcf_hdr_format`.
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

/// `hts_str2dbl`: fifteen digits by integer arithmetic, and `strtod` for what does not fit.
///
/// Returns the value and how much of the text it read, or `None` where nothing could be read.
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
    let start = v;
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
    let point = point.unwrap_or(v - start);
    if budget == 0 || matches!(bytes.get(v), Some(b'e' | b'E')) {
        return strtod(text);
    }
    const POWERS: [f64; 22] = [
        1.0, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15,
        1e16, 1e17, 1e18, 1e19, 1e20,
    ];
    let value = n as f64 / POWERS[v - start - point];
    Some((if negative { -value } else { value }, v))
}

/// glibc's `strtod` on the decimal forms, `inf` and `nan`, which it rounds correctly as `parse` does.
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
    text[number_start..at].parse().ok().map(|value| (value, at))
}

/// `hts_str2int` over a C `long`, and how much of the text it read.
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
    let signed = if negative {
        value.wrapping_neg() as i64
    } else {
        value as i64
    };
    Some((signed, overflow, at))
}

/// `%g` with the default precision.
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
#[derive(Clone, Debug, PartialEq)]
pub enum InfoValue {
    /// A key with no value, which is how a Flag looks, and a text value that is empty.
    Nothing,
    Text(String),
    Integers(Vec<Option<i32>>),
    Floats(Vec<Option<f32>>),
}

/// One INFO entry, and whether the header defines it.
#[derive(Clone, Debug, PartialEq)]
pub struct InfoEntry {
    pub key: String,
    pub value: InfoValue,
    pub defined: bool,
}

/// `bcf1_t` after `vcf_parse`.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub contig: String,
    /// 1-based, as `variant.pos`.
    pub pos: i64,
    pub id: String,
    pub reference: String,
    pub alternates: Vec<String>,
    pub quality: Option<f32>,
    pub filters: Vec<String>,
    pub info: Vec<InfoEntry>,
    pub format: Vec<String>,
    /// For each sample, its fields as rendered.
    pub samples: Vec<Vec<String>>,
    /// `rlen`: the length of REF unless an INFO/END says otherwise.
    pub reference_length: i64,
    /// A contig, filter or FORMAT tag the header does not define.
    pub undefined: Option<String>,
}

const INT32_MISSING_BELOW: i64 = -2147483392;
const INT32_MAX: i64 = 2147483647;

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

fn parse_floats(text: &str) -> Vec<Option<f32>> {
    text.split(',')
        .map(|one| hts_str2dbl(one).map(|(value, _)| value as f32))
        .collect()
}

fn render_floats(values: &[Option<f32>]) -> String {
    values
        .iter()
        .map(|v| v.map_or(".".to_string(), |v| kputd(f64::from(v))))
        .collect::<Vec<_>>()
        .join(",")
}

fn render_integers(values: &[Option<i32>]) -> String {
    values
        .iter()
        .map(|v| v.map_or(".".to_string(), |v| v.to_string()))
        .collect::<Vec<_>>()
        .join(",")
}

impl Record {
    /// `vcf_parse`. Why a record cannot be parsed is what the reference's `OSError` would be about.
    pub fn parse(header: &Header, line: &str) -> Result<Record, String> {
        let mut fields = line.split('\t');
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

    /// `variant.stop`.
    pub fn stop(&self) -> i64 {
        self.pos - 1 + self.reference_length
    }

    /// `has_symbolic_allele`.
    fn has_symbolic_allele(&self) -> bool {
        self.alternates
            .iter()
            .any(|a| a.starts_with('<') && a.ends_with('>'))
    }

    /// `variant.info.update({key: value})` for a Float of `Number=1`.
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
    pub fn format_line(&mut self, header: &Header) -> Result<String, String> {
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
fn render_sample(header: &Header, format: &[String], text: &str) -> Vec<String> {
    let mut values = text.split(':');
    format
        .iter()
        .map(|key| {
            let Some(value) = values.next() else {
                return ".".to_string();
            };
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
