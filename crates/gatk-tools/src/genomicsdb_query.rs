//! A GenomicsDB workspace read back as GATK 4.7.0.0 reads it: `GenomicsDBFeatureReader` over
//! `GenomicsDBQueryStream` with the VCF codec, which is GenomicsDB 1.5.5's
//! `GenomicsDBBCFGenerator` writing VCF text through htslib.
//!
//! The text is the composition of three native stages, each named after the code it stands for:
//! `VariantQueryProcessor::scan_and_operate` cuts the cells of an array into gVCF intervals,
//! `BroadCombinedGVCFOperator::operate` merges the calls of an interval into one record, and
//! htslib's `vcf_format` prints it. Where the native code leans on `std::unordered_set` or
//! `std::unordered_map`, the iteration order of libstdc++ (GCC 7) is reproduced, since it decides
//! the order of IDs and of the histogram annotations. What this port has not measured is refused
//! with a reason rather than guessed at.

use std::fs;
use std::path::{Path, PathBuf};

use serde_json::Value as Json;

use crate::gvcf_blocks::ReferenceBase;

/// What the reference does that this port does not carry, in a sentence.
pub type Reason = String;

const NULL_INT32: i32 = i32::MAX;
const NULL_FLOAT32: u32 = 0x7f7f_ffff;
const NULL_CHAR: u8 = 0x7f;
const INT32_MISSING: i32 = i32::MIN;
const INT32_VECTOR_END: i32 = i32::MIN + 1;
const INT64_MISSING: i64 = i64::MIN;
const INT64_VECTOR_END: i64 = i64::MIN + 1;
const FLOAT_MISSING: u32 = 0x7f80_0001;
const FLOAT_VECTOR_END: u32 = 0x7f80_0002;
const STR_MISSING: u8 = 0x07;
const STR_VECTOR_END: u8 = 0;

/// `GenomicsDBExportConfiguration` as `GATKGenomicsDBUtils.createExportConfiguration` fills it.
#[derive(Debug, Clone, Copy)]
pub struct ExportOptions {
    pub produce_gt: bool,
    pub max_diploid_alt_alleles: u32,
    pub max_genotype_count: u32,
}

// ---------------------------------------------------------------------------------------------
// The vid mapper
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
enum Elem {
    Int,
    Float,
    Char,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Len {
    Fixed(usize),
    Var,
    A,
    R,
    G,
    P,
    PP,
}

impl Len {
    fn allele_dependent(self) -> bool {
        matches!(self, Len::A | Len::R | Len::G)
    }

    fn ploidy(self, elements: usize) -> usize {
        if self == Len::PP {
            elements.div_ceil(2)
        } else {
            elements
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Op {
    None,
    Unknown,
    Sum,
    Mean,
    Median,
    Dp,
    Concatenate,
    ElementWiseSum,
    HistogramSum,
    MoveToFormat,
}

impl Op {
    fn is_sum(self) -> bool {
        matches!(
            self,
            Op::Sum | Op::Dp | Op::ElementWiseSum | Op::HistogramSum
        )
    }
}

/// A `FieldInfo` of the flattened table `VidMapper` builds from `vidmap.json`.
#[derive(Debug, Clone)]
struct VidField {
    name: String,
    vcf_name: String,
    info: bool,
    format: bool,
    filter: bool,
    types: Vec<Elem>,
    lengths: Vec<Len>,
    delimiters: Vec<u8>,
    op: Op,
    remap_missing_with_non_ref: bool,
    parent: Option<usize>,
    tuple_index: usize,
}

impl VidField {
    fn elem(&self) -> Elem {
        self.types[0]
    }

    fn len(&self) -> Len {
        self.lengths[0]
    }

    fn two_dimensional(&self) -> bool {
        self.lengths.len() > 1
    }

    fn multi(&self) -> bool {
        self.types.len() > 1 || self.two_dimensional()
    }

    fn delimiter(&self, dimension: usize) -> u8 {
        self.delimiters
            .get(dimension)
            .copied()
            .unwrap_or(if dimension == 0 { b'|' } else { b',' })
    }
}

/// The `KnownVariantFieldsEnum` values the operator distinguishes.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Known {
    End,
    Ref,
    Alt,
    Id,
    Qual,
    Filter,
    Dp,
    DpFormat,
    MinDp,
    Gq,
    Pl,
    Gt,
    Other,
    None,
}

fn known(name: &str) -> Known {
    match name {
        "END" => Known::End,
        "REF" => Known::Ref,
        "ALT" => Known::Alt,
        "ID" => Known::Id,
        "QUAL" => Known::Qual,
        "FILTER" => Known::Filter,
        "DP" => Known::Dp,
        "DP_FORMAT" => Known::DpFormat,
        "MIN_DP" => Known::MinDp,
        "GQ" => Known::Gq,
        "PL" => Known::Pl,
        "GT" => Known::Gt,
        "BaseQRankSum" | "ClippingRankSum" | "MQRankSum" | "ReadPosRankSum" | "MQ" | "RAW_MQ"
        | "MQ0" | "SB" | "AD" | "AF" | "AN" | "AC" | "PS" | "PGT" | "PID" | "ExcessHet" => {
            Known::Other
        }
        _ => Known::None,
    }
}

/// `KnownFieldInfo`'s combine operation for a known field without one in the JSON.
fn known_op(name: &str) -> Op {
    match name {
        "BaseQRankSum" | "ClippingRankSum" | "MQRankSum" | "ReadPosRankSum" | "MQ" | "MQ0"
        | "ExcessHet" => Op::Median,
        "RAW_MQ" => Op::Sum,
        "DP" => Op::Dp,
        _ => Op::Unknown,
    }
}

/// `KnownFieldInfo`'s length for a known field without one in the JSON.
fn known_len(name: &str) -> Len {
    match name {
        "AD" => Len::R,
        "PL" => Len::G,
        "GT" => Len::PP,
        "AF" | "AC" | "MLEAC" | "MLEAF" => Len::A,
        "SB" => Len::Fixed(4),
        "PGT" | "PID" | "REF" | "ALT" | "ID" | "FILTER" => Len::Var,
        _ => Len::Fixed(1),
    }
}

fn json_i64(value: &Json) -> Option<i64> {
    match value {
        Json::Number(number) => number.as_i64(),
        Json::String(text) => text.parse().ok(),
        _ => None,
    }
}

fn parse_elem(name: &str) -> Result<Elem, Reason> {
    match name {
        "int" | "Integer" => Ok(Elem::Int),
        "float" | "Float" => Ok(Elem::Float),
        "char" | "String" | "Flag" => Ok(Elem::Char),
        other => Err(format!(
            "the vid map has a field of type {other}, which this port's GenomicsDB reader does not read."
        )),
    }
}

fn parse_len(value: &Json) -> Result<Len, Reason> {
    let text = match value {
        Json::Object(map) => map
            .get("variable_length_descriptor")
            .and_then(Json::as_str)
            .map(str::to_string)
            .or_else(|| {
                map.get("fixed_length")
                    .and_then(json_i64)
                    .map(|n| n.to_string())
            }),
        Json::String(text) => Some(text.clone()),
        Json::Number(number) => Some(number.to_string()),
        _ => None,
    };
    match text.as_deref() {
        Some("var") => Ok(Len::Var),
        Some("A") => Ok(Len::A),
        Some("R") => Ok(Len::R),
        Some("G") => Ok(Len::G),
        Some("P") => Ok(Len::P),
        Some("PP") => Ok(Len::PP),
        Some(other) => other.parse().map(Len::Fixed).map_err(|_| {
            format!("the vid map has a length {other}, which this port's GenomicsDB reader does not read.")
        }),
        None => Err("the vid map has a length this port's GenomicsDB reader does not read.".to_string()),
    }
}

fn parse_op(name: &str) -> Result<Op, Reason> {
    match name {
        "none" => Ok(Op::None),
        "sum" => Ok(Op::Sum),
        "mean" => Ok(Op::Mean),
        "median" => Ok(Op::Median),
        "move_to_FORMAT" => Ok(Op::MoveToFormat),
        "element_wise_sum" | "elementwise_sum" => Ok(Op::ElementWiseSum),
        "concatenate" => Ok(Op::Concatenate),
        "histogram_sum" => Ok(Op::HistogramSum),
        other => Err(format!(
            "the vid map combines a field with {other}, which this port's GenomicsDB reader does not read."
        )),
    }
}

/// `VidMapper::parse_vidmap_json` with `flatten_field` and `add_mandatory_fields`.
fn vid_fields(json: &Json) -> Result<Vec<VidField>, Reason> {
    let mut fields: Vec<VidField> = Vec::new();
    let entries = json
        .get("fields")
        .and_then(Json::as_array)
        .ok_or("the workspace's vidmap.json has no fields.")?;
    for entry in entries {
        let name = entry
            .get("name")
            .and_then(Json::as_str)
            .ok_or("the workspace's vidmap.json has a field without a name.")?
            .to_string();
        let classes: Vec<&str> = entry
            .get("vcf_field_class")
            .and_then(Json::as_array)
            .map(|classes| classes.iter().filter_map(Json::as_str).collect())
            .unwrap_or_default();
        let types = match entry.get("type") {
            Some(Json::Array(types)) => types
                .iter()
                .map(|t| parse_elem(t.as_str().unwrap_or("")))
                .collect::<Result<Vec<_>, _>>()?,
            Some(Json::String(text)) => vec![parse_elem(text)?],
            _ => return Err(format!("vid field {name} has no type this port reads.")),
        };
        let lengths = match entry.get("length") {
            Some(Json::Array(items)) => {
                items.iter().map(parse_len).collect::<Result<Vec<_>, _>>()?
            }
            Some(other) => vec![parse_len(other)?],
            None => vec![known_len(&name)],
        };
        let delimiters = match entry.get("vcf_delimiter") {
            Some(Json::Array(items)) => items
                .iter()
                .filter_map(|item| item.as_str().and_then(|text| text.bytes().next()))
                .collect(),
            Some(Json::String(text)) => text.bytes().take(1).collect(),
            _ => Vec::new(),
        };
        let op = match entry
            .get("VCF_field_combine_operation")
            .and_then(Json::as_str)
        {
            Some(text) => parse_op(text)?,
            None => known_op(&name),
        };
        let vcf_name = entry
            .get("vcf_name")
            .and_then(Json::as_str)
            .map(str::to_string)
            .unwrap_or_else(|| name.clone());
        let original = VidField {
            name,
            vcf_name,
            info: classes.contains(&"INFO"),
            format: classes.contains(&"FORMAT"),
            filter: classes.contains(&"FILTER"),
            types,
            lengths,
            delimiters,
            op,
            remap_missing_with_non_ref: !entry
                .get("disable_remap_missing_with_non_ref")
                .and_then(Json::as_bool)
                .unwrap_or(false),
            parent: None,
            tuple_index: 0,
        };
        let original_index = fields.len();
        let both = original.info && original.format;
        fields.push(original.clone());
        let mut format_index = original_index;
        if both {
            fields[original_index].format = false;
            let mut copy = original.clone();
            copy.name.push_str("_FORMAT");
            copy.info = false;
            copy.op = Op::Unknown;
            format_index = fields.len();
            fields.push(copy);
        }
        if original.types.len() > 1 {
            for parent in std::iter::once(original_index).chain(both.then_some(format_index)) {
                for element in 0..original.types.len() {
                    let mut copy = fields[parent].clone();
                    copy.name = format!("{}_tuple_element_{element}", fields[parent].name);
                    copy.types = vec![original.types[element]];
                    copy.parent = Some(parent);
                    copy.tuple_index = element;
                    fields.push(copy);
                }
            }
        }
    }
    for (name, elem, len) in [
        ("END", Elem::Int, Len::Fixed(1)),
        ("REF", Elem::Char, Len::Var),
        ("ALT", Elem::Char, Len::Var),
        ("QUAL", Elem::Float, Len::Fixed(1)),
        ("FILTER", Elem::Int, Len::Var),
    ] {
        if !fields.iter().any(|field| field.name == name) {
            fields.push(VidField {
                name: name.to_string(),
                vcf_name: name.to_string(),
                info: false,
                format: false,
                filter: false,
                types: vec![elem],
                lengths: vec![len],
                delimiters: Vec::new(),
                op: Op::Unknown,
                remap_missing_with_non_ref: true,
                parent: None,
                tuple_index: 0,
            });
        }
    }
    Ok(fields)
}

// ---------------------------------------------------------------------------------------------
// The workspace and its arrays
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct WorkspaceContig {
    name: String,
    length: i64,
    offset: i64,
}

/// One array directory, named `contig$start$end` by `GenomicsDBImport`.
#[derive(Debug, Clone)]
pub struct ArrayDir {
    pub contig: String,
    pub start: i64,
    pub end: i64,
    path: PathBuf,
}

/// What `createExportConfiguration` points the native library at.
#[derive(Debug)]
pub struct Workspace {
    fields: Vec<VidField>,
    contigs: Vec<WorkspaceContig>,
    samples: Vec<String>,
    header_lines: Vec<String>,
    arrays: Vec<ArrayDir>,
}

impl Workspace {
    /// The workspace at `path`: its vid map, callsets, merged header and arrays, the arrays in
    /// `ChrArrayFolderComparator` order.
    pub fn open(path: &Path) -> Result<Workspace, Reason> {
        let read_json = |name: &str| -> Result<Json, Reason> {
            let text = fs::read_to_string(path.join(name))
                .map_err(|_| format!("the workspace has no readable {name}."))?;
            serde_json::from_str(&text).map_err(|_| format!("the workspace's {name} is not JSON."))
        };
        let vid = read_json("vidmap.json")?;
        let fields = vid_fields(&vid)?;
        let mut contigs = Vec::new();
        for contig in vid
            .get("contigs")
            .and_then(Json::as_array)
            .into_iter()
            .flatten()
        {
            let name = contig.get("name").and_then(Json::as_str);
            let length = contig.get("length").and_then(json_i64);
            let offset = contig.get("tiledb_column_offset").and_then(json_i64);
            match (name, length, offset) {
                (Some(name), Some(length), Some(offset)) => contigs.push(WorkspaceContig {
                    name: name.to_string(),
                    length,
                    offset,
                }),
                _ => {
                    return Err(
                        "the workspace's vidmap.json has a contig this port does not read.".into(),
                    )
                }
            }
        }

        let callsets = read_json("callset.json")?;
        let mut rows: Vec<(i64, String)> = Vec::new();
        for callset in callsets
            .get("callsets")
            .and_then(Json::as_array)
            .ok_or("the workspace's callset.json has no callsets.")?
        {
            let name = callset.get("sample_name").and_then(Json::as_str);
            let row = callset.get("row_idx").and_then(json_i64);
            match (name, row) {
                (Some(name), Some(row)) => rows.push((row, name.to_string())),
                _ => {
                    return Err(
                        "the workspace's callset.json has a callset this port does not read."
                            .into(),
                    )
                }
            }
        }
        rows.sort();
        if rows
            .iter()
            .enumerate()
            .any(|(index, (row, _))| *row != index as i64)
        {
            return Err(
                "the workspace's callsets do not fill rows 0 to n-1, which this port's GenomicsDB \
                 reader has not measured."
                    .into(),
            );
        }
        let samples = rows.into_iter().map(|(_, name)| name).collect();

        let header = fs::read_to_string(path.join("vcfheader.vcf"))
            .map_err(|_| "the workspace has no readable vcfheader.vcf.".to_string())?;
        let header_lines = header
            .lines()
            .filter(|line| line.starts_with("##"))
            .map(str::to_string)
            .collect();

        let mut arrays = Vec::new();
        let entries =
            fs::read_dir(path).map_err(|_| "the workspace cannot be listed.".to_string())?;
        for entry in entries.flatten() {
            let array_path = entry.path();
            if !array_path.join("__array_schema.tdb").is_file() {
                continue;
            }
            let name = entry.file_name().to_string_lossy().to_string();
            let mut parts = name.rsplitn(3, '$');
            let end = parts.next().and_then(|text| text.parse::<i64>().ok());
            let start = parts.next().and_then(|text| text.parse::<i64>().ok());
            let contig = parts.next();
            match (contig, start, end) {
                (Some(contig), Some(start), Some(end)) => arrays.push(ArrayDir {
                    contig: contig.to_string(),
                    start,
                    end,
                    path: array_path,
                }),
                _ => {
                    return Err(format!(
                        "the workspace has an array named {name}, which is not contig$start$end \
                         as this port's GenomicsDB reader expects."
                    ))
                }
            }
        }
        arrays.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.start.cmp(&b.start)));
        Ok(Workspace {
            fields,
            contigs,
            samples,
            header_lines,
            arrays,
        })
    }

    pub fn arrays(&self) -> &[ArrayDir] {
        &self.arrays
    }

    fn contig(&self, name: &str) -> Result<&WorkspaceContig, Reason> {
        self.contigs
            .iter()
            .find(|contig| contig.name == name)
            .ok_or_else(|| format!("contig {name} is not in the workspace's vid map."))
    }

    fn field(&self, name: &str) -> Option<usize> {
        self.fields.iter().position(|field| field.name == name)
    }

    /// The text of one `GenomicsDBQueryStream` over `array`, from contig position `begin` to
    /// `end`, both 1-based and inclusive: the header when `header` is set, then the records.
    pub fn stream(
        &self,
        array: &ArrayDir,
        begin: i64,
        end: i64,
        options: &ExportOptions,
        reference: &mut dyn ReferenceBase,
        out: &mut String,
    ) -> Result<(), Reason> {
        let schema = read_schema(&array.path.join("__array_schema.tdb"))?;
        let config = QueryConfig::new(self, &schema)?;
        let contig = self.contig(&array.contig)?;
        let cells = read_cells(&array.path, &schema, &config)?;
        let column_begin = contig.offset + begin - 1;
        let column_end = contig.offset + end.min(contig.length) - 1;
        let mut operator = Operator::new(self, &config, options, contig, reference);
        scan(
            &config,
            &cells,
            self.samples.len(),
            column_begin,
            column_end,
            &mut operator,
            out,
        )
    }

    /// The header `BroadCombinedGVCFOperator`'s constructor leaves in the template, as htslib
    /// prints it, `#CHROM` line included.
    pub fn header(&self) -> Result<String, Reason> {
        let array = self.arrays.first().ok_or("the workspace has no arrays.")?;
        let schema = read_schema(&array.path.join("__array_schema.tdb"))?;
        let config = QueryConfig::new(self, &schema)?;
        let mut lines: Vec<String> = Vec::new();
        let fileformat = self
            .header_lines
            .iter()
            .find(|line| line.starts_with("##fileformat="))
            .cloned()
            .unwrap_or_else(|| "##fileformat=VCFv4.2".to_string());
        lines.push(fileformat);
        lines.push("##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string());
        for line in &self.header_lines {
            if line.starts_with("##fileformat=") || header_id(line) == Some(("FILTER", "PASS")) {
                continue;
            }
            lines.push(line.clone());
        }
        for slot in &config.slots {
            let field = &self.fields[slot.vid];
            if config.info_vector(slot) && !matches!(field.op, Op::None | Op::Unknown) {
                self.add_if_missing(&mut lines, &field.vcf_name, "INFO");
            }
            if config.format_vector(slot) {
                if field.format || field.op == Op::MoveToFormat {
                    self.add_if_missing(&mut lines, &field.vcf_name, "FORMAT");
                } else {
                    self.add_if_missing(&mut lines, "DP", "INFO");
                }
            }
        }
        if config.slots.iter().any(|slot| slot.known == Known::Filter) {
            for field in self.fields.iter().filter(|field| field.filter) {
                self.add_if_missing(&mut lines, &field.vcf_name, "FILTER");
            }
        }
        for contig in &self.contigs {
            if !lines
                .iter()
                .any(|line| header_id(line) == Some(("contig", &contig.name)))
            {
                lines.push(format!(
                    "##contig=<ID={},length={}>",
                    contig.name, contig.length
                ));
            }
        }
        let mut text = String::new();
        for line in lines {
            text.push_str(&line);
            text.push('\n');
        }
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for sample in &self.samples {
            text.push('\t');
            text.push_str(sample);
        }
        text.push('\n');
        Ok(text)
    }

    /// `VCFAdapter::add_field_to_hdr_if_missing`.
    fn add_if_missing(&self, lines: &mut Vec<String>, name: &str, kind: &str) {
        let field = self.fields.iter().find(|field| {
            field.vcf_name == name
                && match kind {
                    "INFO" => field.info,
                    "FORMAT" => field.format,
                    _ => field.filter,
                }
                && field.parent.is_none()
        });
        let multi = field.is_some_and(VidField::multi);
        let existing = lines
            .iter()
            .position(|line| header_id(line) == Some((kind, name)));
        let mut description = None;
        if let Some(index) = existing {
            if !multi {
                return;
            }
            description = header_description(&lines[index]);
            lines.remove(index);
        }
        let mut line = format!("##{kind}=<ID={name}");
        if kind != "FILTER" {
            if kind == "FORMAT" && name == "GT" {
                line.push_str(",Number=1,Type=String,Description=\"Genotype\"");
            } else if let Some(field) = field {
                let flag = false;
                let number = if flag || multi {
                    "1".to_string()
                } else {
                    match field.len() {
                        Len::Fixed(count) => count.to_string(),
                        Len::Var | Len::P | Len::PP => ".".to_string(),
                        Len::A => "A".to_string(),
                        Len::R => "R".to_string(),
                        Len::G => "G".to_string(),
                    }
                };
                let kind_name = if multi {
                    "String"
                } else {
                    match field.elem() {
                        Elem::Int => "Integer",
                        Elem::Float => "Float",
                        Elem::Char => "String",
                    }
                };
                line.push_str(&format!(",Number={number},Type={kind_name}"));
            }
        }
        match description {
            Some(description) if !description.is_empty() => {
                line.push_str(",Description=");
                line.push_str(&description);
            }
            _ => line.push_str(&format!(",Description=\"{name}\"")),
        }
        line.push('>');
        lines.push(line);
    }
}

/// The key and `ID` of a structured header line.
fn header_id(line: &str) -> Option<(&str, &str)> {
    let rest = line.strip_prefix("##")?;
    let (key, value) = rest.split_once("=<ID=")?;
    let id = value.split([',', '>']).next()?;
    Some((key, id))
}

/// The `Description` value of a header line, quotes included, as htslib keeps it.
fn header_description(line: &str) -> Option<String> {
    let start = line.find(",Description=")? + ",Description=".len();
    let rest = &line[start..];
    if let Some(quoted) = rest.strip_prefix('"') {
        let mut escaped = false;
        for (index, byte) in quoted.bytes().enumerate() {
            match byte {
                b'\\' if !escaped => escaped = true,
                b'"' if !escaped => return Some(rest[..index + 2].to_string()),
                _ => escaped = false,
            }
        }
        None
    } else {
        rest.split([',', '>']).next().map(str::to_string)
    }
}

// ---------------------------------------------------------------------------------------------
// TileDB: the array schema and the one fragment
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct SchemaAttribute {
    name: String,
    cell_type: u8,
    values: i32,
}

struct ByteReader<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl ByteReader<'_> {
    fn take(&mut self, count: usize) -> Result<&[u8], Reason> {
        let slice = self
            .bytes
            .get(self.position..self.position + count)
            .ok_or("an array schema is shorter than TileDB writes it.")?;
        self.position += count;
        Ok(slice)
    }

    fn i32(&mut self) -> Result<i32, Reason> {
        Ok(i32::from_le_bytes(
            self.take(4)?.try_into().expect("four bytes"),
        ))
    }

    fn string(&mut self) -> Result<String, Reason> {
        let length = self.i32()?.max(0) as usize;
        Ok(String::from_utf8_lossy(self.take(length)?).to_string())
    }
}

/// `ArraySchema::deserialize`, keeping the attributes.
fn read_schema(path: &Path) -> Result<Vec<SchemaAttribute>, Reason> {
    let bytes = fs::read(path).map_err(|_| "an array has no readable schema.".to_string())?;
    let mut reader = ByteReader {
        bytes: &bytes,
        position: 0,
    };
    if reader.take(4)? != 0xffff_fffdu32.to_le_bytes() {
        return Err(
            "an array schema has a version this port's GenomicsDB reader does not read.".into(),
        );
    }
    reader.string()?;
    reader.string()?;
    reader.take(3 + 8)?;
    let count = reader.i32()?.max(0) as usize;
    let names = (0..count)
        .map(|_| reader.string())
        .collect::<Result<Vec<_>, _>>()?;
    let dimensions = reader.i32()?.max(0) as usize;
    for _ in 0..dimensions {
        reader.string()?;
    }
    let domain = reader.i32()?.max(0) as usize;
    reader.take(domain)?;
    let extents = reader.i32()?.max(0) as usize;
    reader.take(extents)?;
    let types = reader.take(count + 1)?.to_vec();
    let values = (0..count)
        .map(|_| reader.i32())
        .collect::<Result<Vec<_>, _>>()?;
    let compression = reader.take(count + 1)?.to_vec();
    if compression[..count].iter().any(|code| *code != 1) || compression[count] != 0x11 {
        return Err(
            "an array is compressed other than as GenomicsDBImport writes it, which this port's \
             GenomicsDB reader does not read."
                .into(),
        );
    }
    Ok(names
        .into_iter()
        .zip(types)
        .zip(values)
        .map(|((name, cell_type), values)| SchemaAttribute {
            name,
            cell_type,
            values,
        })
        .collect())
}

fn inflate_one(data: &[u8]) -> Result<(Vec<u8>, usize), Reason> {
    use libz_sys::*;
    use std::mem::MaybeUninit;

    let corrupt = || "an array tile is not a zlib stream.".to_string();
    let mut stream = MaybeUninit::<z_stream>::zeroed();
    // SAFETY: `inflateInit_` sets the allocators of a zeroed stream before it is assumed
    // initialised; every pointer handed to zlib outlives the calls, and `out` has `avail_out`
    // bytes of spare capacity at `next_out`.
    unsafe {
        if inflateInit_(
            stream.as_mut_ptr(),
            zlibVersion(),
            std::mem::size_of::<z_stream>() as i32,
        ) != Z_OK
        {
            return Err(corrupt());
        }
        let stream = stream.assume_init_mut();
        stream.next_in = data.as_ptr() as *mut Bytef;
        stream.avail_in = data.len() as uInt;
        let mut out: Vec<u8> = Vec::with_capacity(data.len() * 4 + 64);
        let result = loop {
            if out.capacity() == out.len() {
                out.reserve(out.capacity());
            }
            let spare = out.capacity() - out.len();
            stream.next_out = out.as_mut_ptr().add(out.len());
            stream.avail_out = spare as uInt;
            let status = inflate(stream, Z_NO_FLUSH);
            out.set_len(out.len() + spare - stream.avail_out as usize);
            match status {
                Z_STREAM_END => break Ok(stream.total_in as usize),
                Z_OK => {}
                Z_BUF_ERROR if stream.avail_in > 0 => {}
                _ => break Err(corrupt()),
            }
        };
        inflateEnd(stream);
        result.map(|used| (out, used))
    }
}

/// Each tile of a TileDB file, inflated on its own.
fn tiles(path: &Path) -> Result<Vec<Vec<u8>>, Reason> {
    let bytes = fs::read(path).map_err(|_| format!("{} cannot be read.", path.display()))?;
    let mut rest = &bytes[..];
    let mut out = Vec::new();
    while !rest.is_empty() {
        let (tile, used) = inflate_one(rest)?;
        out.push(tile);
        rest = &rest[used..];
    }
    Ok(out)
}

/// A tile of TileDB's delta filter undone.
fn delta_decode(tile: &[u8], stride: usize) -> Vec<i64> {
    let mut previous = vec![0i64; stride];
    tile.chunks_exact(8)
        .enumerate()
        .map(|(index, chunk)| {
            let slot = index % stride;
            let value =
                previous[slot].wrapping_add(i64::from_le_bytes(chunk.try_into().expect("8 bytes")));
            previous[slot] = value;
            value
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
enum Data {
    Int(Vec<i32>),
    /// The bits, since the missing values are NaNs that must survive.
    Float(Vec<u32>),
    Str(Vec<u8>),
    Alt(Vec<String>),
    /// A multi-D field's entries along its first dimension, each the bytes of its elements.
    Multi(Vec<Vec<u8>>),
}

impl Data {
    fn length(&self) -> usize {
        match self {
            Data::Int(values) => values.len(),
            Data::Float(values) => values.len(),
            Data::Str(values) => values.len(),
            Data::Alt(values) => values.len(),
            Data::Multi(values) => values.len(),
        }
    }
}

#[derive(Debug, Clone)]
struct Cell {
    row: usize,
    column: i64,
    end: i64,
    fields: Vec<Option<Data>>,
}

/// `GenomicsDBMultiDVectorIdx` over the bytes `[size][data][#entries][offsets]`.
fn parse_multi(bytes: &[u8]) -> Option<Vec<Vec<u8>>> {
    let read = |at: usize| -> Option<u64> {
        Some(u64::from_le_bytes(bytes.get(at..at + 8)?.try_into().ok()?))
    };
    let size = read(0)? as usize;
    let entries = read(8 + size)? as usize;
    let offsets: Vec<usize> = (0..=entries)
        .map(|index| read(16 + size + index * 8).map(|value| value as usize))
        .collect::<Option<_>>()?;
    offsets
        .windows(2)
        .map(|pair| bytes.get(8 + pair[0]..8 + pair[1]).map(<[u8]>::to_vec))
        .collect()
}

/// `VariantFieldBase::copy_data_from_tile` for one attribute of one cell: `None` when every
/// element is TileDB's null or there is none.
fn field_data(slot: &Slot, field: &VidField, bytes: &[u8]) -> Result<Option<Data>, Reason> {
    if slot.known == Known::Alt {
        let text = String::from_utf8_lossy(bytes);
        let text = text.split('\0').next().unwrap_or("");
        return Ok(Some(Data::Alt(
            text.split('|')
                .filter(|token| !token.is_empty())
                .map(str::to_string)
                .collect(),
        )));
    }
    if field.two_dimensional() || field.parent.is_some() {
        if bytes.is_empty() || bytes.iter().all(|byte| *byte == 0xff) {
            return Ok(None);
        }
        return parse_multi(bytes)
            .map(|entries| Some(Data::Multi(entries)))
            .ok_or_else(|| {
                format!(
                    "field {} holds a multi-D value this port cannot read.",
                    field.name
                )
            });
    }
    let data = match field.elem() {
        Elem::Int => {
            let values: Vec<i32> = bytes
                .chunks_exact(4)
                .map(|chunk| i32::from_le_bytes(chunk.try_into().expect("4 bytes")))
                .collect();
            if values.iter().all(|value| *value == NULL_INT32) {
                return Ok(None);
            }
            Data::Int(values)
        }
        Elem::Float => {
            let values: Vec<u32> = bytes
                .chunks_exact(4)
                .map(|chunk| u32::from_le_bytes(chunk.try_into().expect("4 bytes")))
                .collect();
            if values.iter().all(|value| *value == NULL_FLOAT32) {
                return Ok(None);
            }
            Data::Float(values)
        }
        Elem::Char => {
            if bytes.iter().all(|value| *value == NULL_CHAR) {
                return Ok(None);
            }
            Data::Str(bytes.to_vec())
        }
    };
    Ok(Some(data))
}

/// The cells of the array's one fragment, in TileDB's column-major order.
fn read_cells(
    array: &Path,
    schema: &[SchemaAttribute],
    config: &QueryConfig,
) -> Result<Vec<Cell>, Reason> {
    let fragments: Vec<PathBuf> = fs::read_dir(array)
        .map_err(|_| "an array cannot be listed.".to_string())?
        .flatten()
        .map(|entry| entry.path())
        .filter(|path| path.join("__tiledb_fragment.tdb").is_file())
        .collect();
    let fragment =
        match fragments.as_slice() {
            [] => return Ok(Vec::new()),
            [one] => one,
            _ => return Err(
                "an array holds more than one fragment, from more than one import batch, which \
                 this port's GenomicsDB reader has not measured."
                    .into(),
            ),
        };
    let mut coordinates = Vec::new();
    for tile in tiles(&fragment.join("__coords.tdb"))? {
        coordinates.extend(delta_decode(&tile, 2));
    }
    let count = coordinates.len() / 2;
    let mut columns: Vec<Vec<Vec<u8>>> = Vec::with_capacity(schema.len());
    for attribute in schema {
        let size = match attribute.cell_type {
            0 | 2 => 4,
            1 => 8,
            4 => 1,
            _ => {
                return Err(format!(
                    "attribute {} has a TileDB type this port's GenomicsDB reader does not read.",
                    attribute.name
                ))
            }
        };
        let file = fragment.join(format!("{}.tdb", attribute.name));
        let mut cells = Vec::with_capacity(count);
        if attribute.values == i32::MAX {
            let mut offsets = Vec::with_capacity(count);
            for tile in tiles(&file)? {
                offsets.extend(delta_decode(&tile, 1));
            }
            let data: Vec<u8> =
                tiles(&fragment.join(format!("{}_var.tdb", attribute.name)))?.concat();
            for (index, offset) in offsets.iter().enumerate() {
                let end = offsets.get(index + 1).copied().unwrap_or(data.len() as i64);
                cells.push(
                    data.get(*offset as usize..end as usize)
                        .ok_or("an array's variable-length offsets run past its data.")?
                        .to_vec(),
                );
            }
        } else {
            let width = size * attribute.values.max(0) as usize;
            let data = tiles(&file)?.concat();
            cells.extend(data.chunks_exact(width.max(1)).map(<[u8]>::to_vec));
        }
        if cells.len() != count {
            return Err(format!(
                "attribute {} does not have a value for every cell.",
                attribute.name
            ));
        }
        columns.push(cells);
    }
    let mut out = Vec::with_capacity(count);
    for index in 0..count {
        let mut fields = Vec::with_capacity(config.slots.len());
        let mut end = 0;
        for (slot_index, slot) in config.slots.iter().enumerate() {
            let bytes = &columns[slot_index][index];
            if slot.known == Known::End {
                end = i64::from_le_bytes(
                    bytes
                        .get(..8)
                        .ok_or("an END value is short.")?
                        .try_into()
                        .expect("8 bytes"),
                );
                fields.push(None);
            } else {
                fields.push(field_data(slot, &config.fields[slot.vid], bytes)?);
            }
        }
        let row = coordinates[2 * index];
        if row < 0 {
            return Err("a cell has a negative row.".into());
        }
        out.push(Cell {
            row: row as usize,
            column: coordinates[2 * index + 1],
            end,
            fields,
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// VariantQueryConfig
// ---------------------------------------------------------------------------------------------

/// A queried attribute, in schema order: `VariantQueryConfig`'s query index is its position.
#[derive(Debug, Clone)]
struct Slot {
    vid: usize,
    known: Known,
}

struct QueryConfig {
    fields: Vec<VidField>,
    slots: Vec<Slot>,
    ref_slot: usize,
    alt_slot: usize,
    id_slot: Option<usize>,
    gt_slot: Option<usize>,
    pl_slot: Option<usize>,
}

impl QueryConfig {
    fn new(workspace: &Workspace, schema: &[SchemaAttribute]) -> Result<QueryConfig, Reason> {
        let mut slots = Vec::with_capacity(schema.len());
        for attribute in schema {
            let vid = workspace.field(&attribute.name).ok_or_else(|| {
                format!(
                    "attribute {} is not in the workspace's vid map.",
                    attribute.name
                )
            })?;
            slots.push(Slot {
                vid,
                known: known(&attribute.name),
            });
        }
        let find = |wanted: Known| slots.iter().position(|slot| slot.known == wanted);
        if find(Known::End) != Some(0) {
            return Err("an array's first attribute is not END.".into());
        }
        let config = QueryConfig {
            fields: workspace.fields.clone(),
            ref_slot: find(Known::Ref).ok_or("an array has no REF attribute.")?,
            alt_slot: find(Known::Alt).ok_or("an array has no ALT attribute.")?,
            id_slot: find(Known::Id),
            gt_slot: find(Known::Gt),
            pl_slot: find(Known::Pl),
            slots,
        };
        for slot in &config.slots {
            let field = &config.fields[slot.vid];
            if config.info_vector(slot)
                && field.len().allele_dependent()
                && matches!(field.op, Op::Sum | Op::Mean | Op::Median | Op::Concatenate)
            {
                return Err(format!(
                    "INFO field {} is allele-dependent and combined by {:?}, which this port's \
                     GenomicsDB reader has not measured.",
                    field.name, field.op
                ));
            }
            if slot.known == Known::Qual && field.op != Op::Unknown {
                return Err("the vid map combines QUAL, which this port's GenomicsDB reader has not measured.".into());
            }
            if config.format_vector(slot) && field.multi() {
                return Err(format!(
                    "FORMAT field {} is multi-D, which this port's GenomicsDB reader has not measured.",
                    field.name
                ));
            }
        }
        Ok(config)
    }

    fn field(&self, slot: usize) -> &VidField {
        &self.fields[self.slots[slot].vid]
    }

    /// `add_to_INFO_vector` in the operator's constructor.
    fn info_vector(&self, slot: &Slot) -> bool {
        let field = &self.fields[slot.vid];
        field.info
            && slot.known != Known::End
            && !(slot.known == Known::Dp && field.op == Op::Dp)
            && field.op != Op::MoveToFormat
    }

    /// `add_to_FORMAT_vector` in the operator's constructor, for a query that is not sites-only.
    fn format_vector(&self, slot: &Slot) -> bool {
        let field = &self.fields[slot.vid];
        field.format
            || (field.info
                && ((slot.known == Known::Dp && field.op == Op::Dp)
                    || field.op == Op::MoveToFormat))
    }
}

// ---------------------------------------------------------------------------------------------
// scan_and_operate
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
struct Call {
    valid: bool,
    begin: i64,
    end: i64,
    deletion_or_mnv: bool,
    reference_block: bool,
    fields: Vec<Option<Data>>,
}

fn is_non_ref(allele: &str) -> bool {
    allele == "&" || allele == "<NON_REF>"
}

/// `VariantUtils::is_symbolic_allele`.
fn is_symbolic(allele: &str) -> bool {
    is_non_ref(allele)
        || allele == "*"
        || (allele.starts_with('<') && allele.ends_with('>') && !allele.is_empty())
        || allele.contains('[')
        || allele.contains(']')
}

fn is_deletion(reference: &str, allele: &str) -> bool {
    reference.len() > 1
        && (allele == "*" || (!is_symbolic(allele) && allele.len() < reference.len()))
}

fn is_mnv(reference: &str, allele: &str) -> bool {
    reference.len() > 1
        && !is_symbolic(allele)
        && allele.len() == reference.len()
        && reference
            .bytes()
            .zip(allele.bytes())
            .filter(|(a, b)| a != b)
            .count()
            > 1
}

fn str_of(data: &Option<Data>) -> &str {
    match data {
        Some(Data::Str(bytes)) => std::str::from_utf8(bytes).unwrap_or(""),
        _ => "",
    }
}

fn alts_of(data: &Option<Data>) -> &[String] {
    match data {
        Some(Data::Alt(alleles)) => alleles,
        _ => &[],
    }
}

/// `gt_fill_row`: the call of `cell`'s row from that cell.
fn fill_row(config: &QueryConfig, call: &mut Call, cell: &Cell) {
    call.valid = true;
    let (begin, end) = if cell.column > cell.end {
        (cell.end, cell.column)
    } else {
        (cell.column, cell.end)
    };
    call.begin = begin;
    call.end = end;
    call.fields = cell.fields.clone();
    let reference = str_of(&call.fields[config.ref_slot]);
    let alts = alts_of(&call.fields[config.alt_slot]);
    let deletion = reference.len() > 1 && alts.iter().any(|alt| is_deletion(reference, alt));
    let mnv = reference.len() > 1 && alts.iter().any(|alt| is_mnv(reference, alt));
    call.deletion_or_mnv = deletion || mnv;
    call.reference_block = reference.len() == 1 && alts.len() == 1 && is_non_ref(&alts[0]);
}

/// `handle_gvcf_ranges`.
#[allow(clippy::too_many_arguments)]
fn handle_ranges(
    calls: &mut [Call],
    current: &mut i64,
    next: i64,
    deletions: &mut usize,
    operator: &mut Operator,
    out: &mut String,
) -> Result<(), Reason> {
    while calls.iter().any(|call| call.valid) && *current < next {
        let top = calls
            .iter()
            .filter(|call| call.valid)
            .map(|call| call.end)
            .min()
            .expect("a valid call");
        let mut end = if top < next - 1 { top } else { next - 1 };
        if *deletions > 0 {
            end = *current;
        }
        operator.operate(calls, *current, end, out)?;
        for call in calls
            .iter_mut()
            .filter(|call| call.valid && call.end == end)
        {
            if call.deletion_or_mnv {
                *deletions -= 1;
            }
            call.valid = false;
        }
        *current = end + 1;
    }
    Ok(())
}

/// `scan_and_operate` over one queried column interval, with spanning deletions handled.
fn scan(
    config: &QueryConfig,
    cells: &[Cell],
    rows: usize,
    begin: i64,
    end: i64,
    operator: &mut Operator,
    out: &mut String,
) -> Result<(), Reason> {
    let mut calls = vec![Call::default(); rows];
    for cell in cells {
        if cell.row >= rows {
            return Err("a cell's row is not one of the workspace's callsets.".into());
        }
    }
    // gt_get_column: the first cell at or after `begin` in each row.
    let first = cells.partition_point(|cell| cell.column < begin);
    let mut seen = vec![false; rows];
    for cell in &cells[first..] {
        if seen[cell.row] {
            continue;
        }
        seen[cell.row] = true;
        if cell.column > cell.end || cell.column == begin {
            fill_row(config, &mut calls[cell.row], cell);
        }
    }
    let mut deletions = calls
        .iter()
        .filter(|call| call.valid && call.deletion_or_mnv)
        .count();
    let mut current = if calls.iter().any(|call| call.valid) {
        begin
    } else {
        -1
    };
    let forward = cells.partition_point(|cell| cell.column < begin + 1);
    if current < 0 {
        if let Some(cell) = cells.get(forward) {
            current = cell.column;
        }
    }
    for cell in &cells[forward..] {
        if cell.column > cell.end {
            continue;
        }
        if cell.column > end {
            break;
        }
        if cell.column != current {
            handle_ranges(
                &mut calls,
                &mut current,
                cell.column,
                &mut deletions,
                operator,
                out,
            )?;
            current = cell.column;
        }
        let call = &mut calls[cell.row];
        if call.valid && call.end >= cell.column {
            if !call.deletion_or_mnv && !call.reference_block {
                return Err(format!(
                    "the workspace has overlapping variants in one sample around column {}, which \
                     the native query stops on.",
                    cell.column
                ));
            }
            if call.deletion_or_mnv {
                deletions -= 1;
            }
        }
        *call = Call::default();
        fill_row(config, call, cell);
        if call.deletion_or_mnv {
            deletions += 1;
        }
    }
    let next = if end == i64::MAX { end } else { end + 1 };
    handle_ranges(
        &mut calls,
        &mut current,
        next,
        &mut deletions,
        operator,
        out,
    )
}

// ---------------------------------------------------------------------------------------------
// Allele remapping
// ---------------------------------------------------------------------------------------------

/// `CombineAllelesLUT` for one call: input allele to merged allele and back.
#[derive(Debug, Clone, Default)]
struct Lut {
    to_merged: Vec<Option<usize>>,
    to_input: Vec<Option<usize>>,
}

impl Lut {
    fn add(&mut self, input: usize, merged: usize) {
        if self.to_merged.len() <= input {
            self.to_merged.resize(input + 1, None);
        }
        if self.to_input.len() <= merged {
            self.to_input.resize(merged + 1, None);
        }
        self.to_merged[input] = Some(merged);
        self.to_input[merged] = Some(input);
    }

    fn merged(&self, input: usize) -> Option<usize> {
        self.to_merged.get(input).copied().flatten()
    }

    fn input(&self, merged: usize) -> Option<usize> {
        self.to_input.get(merged).copied().flatten()
    }
}

/// `KnownFieldInfo::get_number_of_genotypes`.
fn genotype_count(num_alt: usize, ploidy: usize) -> u128 {
    match ploidy {
        1 => num_alt as u128 + 1,
        2 => (num_alt as u128 + 1) * (num_alt as u128 + 2) / 2,
        _ => binomial(ploidy + num_alt, num_alt),
    }
}

fn binomial(n: usize, k: usize) -> u128 {
    let k = k.min(n.saturating_sub(k));
    let mut result: u128 = 1;
    for i in 0..k {
        result = result.saturating_mul((n - i) as u128) / (i as u128 + 1);
    }
    result
}

/// `VariantOperations::get_genotype_index` of sorted allele indices.
fn genotype_index(alleles: &[usize]) -> u128 {
    alleles
        .iter()
        .enumerate()
        .map(|(i, allele)| {
            if *allele == 0 {
                0
            } else {
                binomial(i + allele, i + 1)
            }
        })
        .sum()
}

/// Every genotype of `ploidy` over `alleles` alleles, in VCF order, as sorted allele lists.
fn genotypes(alleles: usize, ploidy: usize) -> Vec<Vec<usize>> {
    if ploidy == 0 {
        return vec![Vec::new()];
    }
    let mut out = Vec::new();
    for last in 0..alleles {
        for mut prefix in genotypes(last + 1, ploidy - 1) {
            prefix.push(last);
            out.push(prefix);
        }
    }
    out
}

/// `VariantFieldHandler::remap_vector_data`: an `A`, `R` or `G` vector in the merged allele
/// order of `lut`, a merged allele the input lacks read from its `<NON_REF>` when `non_ref`.
fn remap_vector(
    data: &Data,
    lut: &Lut,
    merged_alleles: usize,
    non_ref: bool,
    ploidy: usize,
    len: Len,
) -> Data {
    let input_non_ref = if non_ref {
        lut.input(merged_alleles - 1)
    } else {
        None
    };
    let map = |merged: usize| lut.input(merged).or(input_non_ref);
    let pick = |source: Option<usize>| -> Option<usize> { source };
    let indices: Vec<Option<usize>> = match len {
        Len::G => {
            if ploidy == 0 {
                // The remap writes nothing: the resized copy of the input stays.
                return resize(data, 1);
            }
            genotypes(merged_alleles, ploidy)
                .into_iter()
                .map(|genotype| {
                    let mut input = Vec::with_capacity(ploidy);
                    for allele in genotype {
                        input.push(map(allele)?);
                    }
                    input.sort_unstable();
                    usize::try_from(genotype_index(&input)).ok()
                })
                .collect()
        }
        Len::A => (1..merged_alleles)
            .map(|merged| pick(map(merged)).map(|input| input - 1))
            .collect(),
        _ => (0..merged_alleles).map(map).collect(),
    };
    match data {
        Data::Int(values) => Data::Int(
            indices
                .iter()
                .map(|index| {
                    index
                        .and_then(|i| values.get(i).copied())
                        .unwrap_or(INT32_MISSING)
                })
                .collect(),
        ),
        Data::Float(values) => Data::Float(
            indices
                .iter()
                .map(|index| {
                    index
                        .and_then(|i| values.get(i).copied())
                        .unwrap_or(FLOAT_MISSING)
                })
                .collect(),
        ),
        other => other.clone(),
    }
}

fn resize(data: &Data, length: usize) -> Data {
    match data {
        Data::Int(values) => {
            let mut values = values.clone();
            values.resize(length, 0);
            Data::Int(values)
        }
        Data::Float(values) => {
            let mut values = values.clone();
            values.resize(length, 0);
            Data::Float(values)
        }
        other => other.clone(),
    }
}

/// `remap_allele_specific_annotations`: a multi-D field's entries in the merged allele order.
fn remap_multi(
    entries: &[Vec<u8>],
    lut: &Lut,
    merged_alleles: usize,
    non_ref: bool,
    len: Len,
) -> Vec<Vec<u8>> {
    let input_non_ref = if non_ref {
        lut.input(merged_alleles - 1)
    } else {
        None
    };
    let alt_only = len == Len::A;
    let alleles: Vec<usize> = if alt_only {
        (1..merged_alleles).collect()
    } else {
        (0..merged_alleles).collect()
    };
    alleles
        .into_iter()
        .map(|merged| {
            lut.input(merged)
                .or(input_non_ref)
                .map(|input| {
                    if alt_only {
                        input.wrapping_sub(1)
                    } else {
                        input
                    }
                })
                .and_then(|input| entries.get(input).cloned())
                .unwrap_or_default()
        })
        .collect()
}

/// `VariantOperations::remap_GT_field`.
fn remap_gt(
    input: &[i32],
    lut: &Lut,
    merged_alleles: usize,
    non_ref: bool,
    phased: bool,
) -> Vec<i32> {
    let step = if phased { 2 } else { 1 };
    let mut out = input.to_vec();
    for index in (0..input.len()).step_by(step) {
        let value = input[index];
        if value == NULL_INT32 || value == -1 || value == INT32_MISSING || value < 0 {
            continue;
        }
        out[index] = match lut.merged(value as usize) {
            Some(merged) => merged as i32,
            None if non_ref => merged_alleles as i32 - 1,
            None => -1,
        };
    }
    out
}

// ---------------------------------------------------------------------------------------------
// BroadCombinedGVCFOperator
// ---------------------------------------------------------------------------------------------

fn int_valid(value: i32) -> bool {
    value != INT32_MISSING && value != INT32_VECTOR_END
}

fn float_valid(bits: u32) -> bool {
    bits != FLOAT_MISSING && bits != FLOAT_VECTOR_END
}

/// An INFO value as `bcf_update_info` stores it.
enum InfoValue {
    Ints(Vec<i64>),
    Floats(Vec<u32>),
    Text(String),
}

/// A FORMAT column as `bcf_update_format` stores it: values per sample.
enum FormatValue {
    Gt(Vec<i32>),
    Ints(Vec<i32>),
    Floats(Vec<u32>),
    Text(Vec<u8>),
}

struct Operator<'a> {
    config: &'a QueryConfig,
    options: ExportOptions,
    contig: &'a WorkspaceContig,
    reference: &'a mut dyn ReferenceBase,
    info_slots: Vec<usize>,
    /// (bin slot, count slot, parent field) in `std::unordered_map` iteration order.
    histograms: Vec<(usize, usize, usize)>,
    format_slots: Vec<usize>,
    dp_info_slot: Option<usize>,
    allele_slots: Vec<usize>,
}

impl<'a> Operator<'a> {
    fn new(
        workspace: &Workspace,
        config: &'a QueryConfig,
        options: &ExportOptions,
        contig: &'a WorkspaceContig,
        reference: &'a mut dyn ReferenceBase,
    ) -> Operator<'a> {
        let _ = workspace;
        let mut info_slots = Vec::new();
        let mut format_slots = Vec::new();
        let mut dp_info_slot = None;
        let mut parents: Vec<(usize, [usize; 2])> = Vec::new();
        for (index, slot) in config.slots.iter().enumerate() {
            let field = &config.fields[slot.vid];
            if config.info_vector(slot) {
                match field.op {
                    Op::None | Op::Unknown => {}
                    Op::HistogramSum => {
                        info_slots.push(index);
                        let parent = field.parent.unwrap_or(slot.vid);
                        let at = match parents.iter().position(|(key, _)| *key == parent) {
                            Some(at) => at,
                            None => {
                                parents.push((parent, [0, 0]));
                                parents.len() - 1
                            }
                        };
                        parents[at].1[field.tuple_index.min(1)] = index;
                    }
                    _ => info_slots.push(index),
                }
            }
            if config.format_vector(slot) {
                if field.format || field.op == Op::MoveToFormat {
                    format_slots.push(index);
                } else {
                    dp_info_slot = Some(index);
                }
            }
        }
        let order = unordered_order(
            &parents
                .iter()
                .map(|(key, _)| *key as u64)
                .collect::<Vec<_>>(),
        )
        .unwrap_or_else(|_| (0..parents.len()).collect());
        let histograms = order
            .into_iter()
            .map(|at| (parents[at].1[0], parents[at].1[1], parents[at].0))
            .collect();
        let allele_slots = (0..config.slots.len())
            .filter(|index| config.field(*index).len().allele_dependent())
            .collect();
        Operator {
            config,
            options: *options,
            contig,
            reference,
            info_slots,
            histograms,
            format_slots,
            dp_info_slot,
            allele_slots,
        }
    }

    fn gt_phased(&self) -> bool {
        self.config
            .gt_slot
            .is_some_and(|slot| self.config.field(slot).len() == Len::PP)
    }

    fn ploidy_of(&self, call: &Call) -> usize {
        match self
            .config
            .gt_slot
            .and_then(|slot| call.fields[slot].as_ref())
        {
            Some(gt) => self
                .config
                .field(self.config.gt_slot.expect("GT"))
                .len()
                .ploidy(gt.length()),
            None => 0,
        }
    }

    /// `BroadCombinedGVCFOperator::handle_deletions`: a deletion begun before this position
    /// becomes `*` (and `<NON_REF>`), its INFO fields dropped.
    fn handle_deletions(&self, calls: &mut [Call], begin: i64) {
        let config = self.config;
        for call in calls.iter_mut() {
            if !(call.valid && call.deletion_or_mnv && begin > call.begin) {
                continue;
            }
            for slot in &self.info_slots {
                call.fields[*slot] = None;
            }
            let alts = alts_of(&call.fields[config.alt_slot]).to_vec();
            if alts.first().map(String::as_str) == Some("*")
                && (alts.len() == 1 || (alts.len() == 2 && is_non_ref(&alts[1])))
            {
                continue;
            }
            let reference = str_of(&call.fields[config.ref_slot]).to_string();
            let ploidy = self.ploidy_of(call);
            let mut lut = Lut::default();
            lut.add(0, 0);
            let pl: Option<Vec<i32>> =
                match config.pl_slot.and_then(|slot| call.fields[slot].as_ref()) {
                    Some(Data::Int(values)) => Some(values.clone()),
                    _ => None,
                };
            let mut lowest: Option<usize> = None;
            let mut lowest_pl = i32::MAX;
            let mut has_non_ref = false;
            for (i, alt) in alts.iter().enumerate() {
                let allele = i + 1;
                if is_deletion(&reference, alt) || is_mnv(&reference, alt) {
                    if lowest.is_none() {
                        lowest = Some(allele);
                    }
                    if let Some(pl) = &pl {
                        let index = genotype_index(&vec![allele; ploidy]);
                        if let Some(value) = usize::try_from(index).ok().and_then(|i| pl.get(i)) {
                            if *value < lowest_pl {
                                lowest_pl = *value;
                                lowest = Some(allele);
                            }
                        }
                    }
                } else if is_non_ref(alt) {
                    lut.add(allele, 2);
                    has_non_ref = true;
                }
            }
            let Some(lowest) = lowest else { continue };
            lut.add(lowest, 1);
            let reduced_alts: Vec<String> = if has_non_ref {
                vec!["*".to_string(), "&".to_string()]
            } else {
                vec!["*".to_string()]
            };
            let reduced = reduced_alts.len() + 1;
            call.fields[config.ref_slot] = Some(Data::Str(b"N".to_vec()));
            call.fields[config.alt_slot] = Some(Data::Alt(reduced_alts));
            for slot in &self.allele_slots {
                let field = config.field(*slot);
                let Some(data) = call.fields[*slot].take() else {
                    continue;
                };
                call.fields[*slot] = Some(match &data {
                    Data::Multi(entries) => Data::Multi(remap_multi(
                        entries,
                        &lut,
                        reduced,
                        has_non_ref,
                        field.len(),
                    )),
                    _ => remap_vector(&data, &lut, reduced, has_non_ref, ploidy, field.len()),
                });
            }
            if let Some(slot) = config.gt_slot {
                if let Some(Data::Int(gt)) = &call.fields[slot] {
                    let remapped = remap_gt(gt, &lut, reduced, has_non_ref, self.gt_phased());
                    call.fields[slot] = Some(Data::Int(remapped));
                }
            }
        }
    }

    fn operate(
        &mut self,
        calls: &mut [Call],
        begin: i64,
        end: i64,
        out: &mut String,
    ) -> Result<(), Reason> {
        let config = self.config;
        self.handle_deletions(calls, begin);

        // SingleVariantOperatorBase::operate: the merged REF and ALT alleles.
        let mut merged_ref = String::from("N");
        for call in calls
            .iter()
            .filter(|call| call.valid && call.begin >= begin)
        {
            let current = str_of(&call.fields[config.ref_slot]);
            if current.len() > merged_ref.len() {
                if merged_ref.starts_with('N') {
                    merged_ref = current.to_string();
                } else {
                    merged_ref.push_str(&current[merged_ref.len()..]);
                }
            } else if merged_ref.starts_with('N') && !current.starts_with('N') {
                merged_ref = current.to_string();
            }
        }
        let mut merged_alts: Vec<String> = Vec::new();
        let mut luts = vec![Lut::default(); calls.len()];
        let mut input_non_ref: Vec<Option<usize>> = vec![None; calls.len()];
        let mut non_ref_exists = false;
        for (index, call) in calls.iter().enumerate().filter(|(_, call)| call.valid) {
            let reference = str_of(&call.fields[config.ref_slot]);
            luts[index].add(0, 0);
            for (i, allele) in alts_of(&call.fields[config.alt_slot]).iter().enumerate() {
                let input = i + 1;
                if is_non_ref(allele) {
                    input_non_ref[index] = Some(input);
                    non_ref_exists = true;
                    continue;
                }
                let mut allele = allele.clone();
                if reference.len() < merged_ref.len() && !is_symbolic(&allele) {
                    allele.push_str(&merged_ref[reference.len()..]);
                }
                match merged_alts.iter().position(|seen| *seen == allele) {
                    Some(position) => luts[index].add(input, position + 1),
                    None => {
                        merged_alts.push(allele);
                        luts[index].add(input, merged_alts.len());
                    }
                }
            }
        }
        if non_ref_exists {
            merged_alts.push("<NON_REF>".to_string());
            for (index, input) in input_non_ref.iter().enumerate() {
                if let Some(input) = input {
                    luts[index].add(*input, merged_alts.len());
                }
            }
        }
        let reference_block_only =
            merged_ref.len() == 1 && merged_alts.len() == 1 && merged_alts[0] == "<NON_REF>";
        let remapping = !reference_block_only;
        let merged_alleles = merged_alts.len() + 1;
        let too_many_alts = merged_alts.len() > self.options.max_diploid_alt_alleles as usize;

        // GA4GHOperator::operate: GT and the allele-dependent fields in the merged order.
        let mut remapped: Vec<Vec<Option<Data>>> =
            calls.iter().map(|call| call.fields.clone()).collect();
        let mut ploidy = vec![0usize; calls.len()];
        if remapping {
            if let Some(slot) = config.gt_slot {
                let non_ref = non_ref_exists && config.field(slot).remap_missing_with_non_ref;
                for (index, call) in calls.iter().enumerate().filter(|(_, call)| call.valid) {
                    if let Some(Data::Int(gt)) = &call.fields[slot] {
                        remapped[index][slot] = Some(Data::Int(remap_gt(
                            gt,
                            &luts[index],
                            merged_alleles,
                            non_ref,
                            self.gt_phased(),
                        )));
                        ploidy[index] = config.field(slot).len().ploidy(gt.len());
                    }
                }
            }
            for slot in &self.allele_slots {
                let field = config.field(*slot);
                if field.info && field.op.is_sum() {
                    continue;
                }
                if field.len() == Len::G && too_many_alts {
                    continue;
                }
                for (index, call) in calls.iter().enumerate().filter(|(_, call)| call.valid) {
                    remapped[index][*slot] = self.remap_field(
                        call.fields[*slot].as_ref(),
                        field,
                        &luts[index],
                        merged_alleles,
                        non_ref_exists,
                        ploidy[index],
                    );
                }
            }
        }
        let source = |slot: usize| -> bool {
            remapping
                && (config.field(slot).len().allele_dependent() || Some(slot) == config.gt_slot)
        };

        // The record.
        let line_start = out.len();
        out.push_str(&self.contig.name);
        out.push('\t');
        let position = begin - self.contig.offset + 1;
        out.push_str(&position.to_string());
        out.push('\t');
        let mut ids = Vec::new();
        if let Some(slot) = config.id_slot {
            for call in calls.iter().filter(|call| call.valid) {
                if let Some(Data::Str(bytes)) = &call.fields[slot] {
                    let text = String::from_utf8_lossy(bytes).to_string();
                    for id in text.split(';').filter(|id| !id.is_empty()) {
                        if !ids.iter().any(|seen| seen == id) {
                            ids.push(id.to_string());
                        }
                    }
                }
            }
        }
        if ids.is_empty() {
            out.push('.');
        } else {
            let codes: Vec<u64> = ids.iter().map(|id| hash_bytes(id.as_bytes())).collect();
            let order = unordered_order(&codes)?;
            let ordered: Vec<&str> = order.iter().map(|at| ids[*at].as_str()).collect();
            out.push_str(&ordered.join(";"));
        }
        out.push('\t');
        if merged_ref == "N" {
            let base = self.reference.base(&self.contig.name, position);
            let base = if b"ACGTacgt".contains(&base) {
                base.to_ascii_uppercase()
            } else {
                b'N'
            };
            merged_ref = (base as char).to_string();
        }
        out.push_str(&merged_ref);
        out.push('\t');
        if merged_alts.is_empty() {
            out.push('.');
        } else {
            out.push_str(&merged_alts.join(","));
        }
        out.push_str("\t.\t.\t");

        // INFO.
        let mut info: Vec<(String, InfoValue)> = Vec::new();
        if end > begin {
            info.push((
                "END".to_string(),
                InfoValue::Ints(vec![end - self.contig.offset + 1]),
            ));
        }
        for slot in &self.info_slots {
            let field = config.field(*slot);
            if field.len() == Len::G && too_many_alts {
                continue;
            }
            let rows: Vec<Option<&Data>> = calls
                .iter()
                .enumerate()
                .filter(|(_, call)| call.valid)
                .map(|(index, call)| {
                    if source(*slot) {
                        remapped[index][*slot].as_ref()
                    } else {
                        call.fields[*slot].as_ref()
                    }
                })
                .collect();
            let value = match field.op {
                Op::Sum | Op::Mean => sum_or_mean(&rows, field),
                Op::Median => median(&rows),
                Op::Concatenate => concatenate(&rows),
                Op::ElementWiseSum => {
                    let remap = remapping && field.len().allele_dependent();
                    if remap && field.len() == Len::G && too_many_alts {
                        None
                    } else {
                        let inputs: Vec<Option<Data>> = calls
                            .iter()
                            .enumerate()
                            .filter(|(_, call)| call.valid)
                            .map(|(index, call)| {
                                if remap {
                                    self.remap_field(
                                        call.fields[*slot].as_ref(),
                                        field,
                                        &luts[index],
                                        merged_alleles,
                                        non_ref_exists,
                                        ploidy[index],
                                    )
                                } else {
                                    call.fields[*slot].clone()
                                }
                            })
                            .collect();
                        if field.two_dimensional() {
                            element_wise_sum_2d(&inputs, field)
                        } else {
                            element_wise_sum(&inputs)
                        }
                    }
                }
                _ => None,
            };
            if let Some(value) = value {
                info.push((field.vcf_name.clone(), value));
            }
        }
        for (bin_slot, count_slot, parent) in &self.histograms {
            let remap = remapping && config.field(*bin_slot).len().allele_dependent();
            let mut sums: Vec<Vec<(f64, HistogramCount)>> = Vec::new();
            let mut found = false;
            for (index, call) in calls.iter().enumerate().filter(|(_, call)| call.valid) {
                let get = |slot: usize| {
                    if remap {
                        self.remap_field(
                            call.fields[slot].as_ref(),
                            config.field(slot),
                            &luts[index],
                            merged_alleles,
                            non_ref_exists,
                            ploidy[index],
                        )
                    } else {
                        call.fields[slot].clone()
                    }
                };
                let (Some(Data::Multi(bins)), Some(Data::Multi(counts))) =
                    (get(*bin_slot), get(*count_slot))
                else {
                    continue;
                };
                found = true;
                if bins.len() > sums.len() {
                    sums.resize(bins.len(), Vec::new());
                }
                let bin_elem = config.field(*bin_slot).elem();
                let count_elem = config.field(*count_slot).elem();
                for (entry, (bin_bytes, count_bytes)) in bins.iter().zip(counts.iter()).enumerate()
                {
                    for (bin, count) in bin_bytes.chunks_exact(4).zip(count_bytes.chunks_exact(4)) {
                        let bin = u32::from_le_bytes(bin.try_into().expect("4 bytes"));
                        let count = u32::from_le_bytes(count.try_into().expect("4 bytes"));
                        let (Some(bin), Some(count)) = (
                            element_value(bin, bin_elem),
                            element_value(count, count_elem),
                        ) else {
                            continue;
                        };
                        let count = match count_elem {
                            Elem::Float => HistogramCount::Float(count as f32),
                            _ => HistogramCount::Int(count as i64),
                        };
                        let map = &mut sums[entry];
                        match map.binary_search_by(|(key, _)| {
                            key.partial_cmp(&bin).unwrap_or(std::cmp::Ordering::Less)
                        }) {
                            Ok(at) => map[at].1 = map[at].1.add(count),
                            Err(at) => map.insert(at, (bin, count)),
                        }
                    }
                }
            }
            if !found {
                continue;
            }
            let parent_field = &config.fields[*parent];
            let bin_float = config.field(*bin_slot).elem() == Elem::Float;
            let mut text = String::new();
            for (entry, map) in sums.iter().enumerate() {
                if entry > 0 {
                    text.push(parent_field.delimiter(0) as char);
                }
                for (i, (bin, count)) in map.iter().enumerate() {
                    if i > 0 {
                        text.push(parent_field.delimiter(1) as char);
                    }
                    if bin_float {
                        text.push_str(&format!("{:.3}", *bin as f32 as f64));
                    } else {
                        text.push_str(&format!("{}", *bin as i64));
                    }
                    text.push(parent_field.delimiter(1) as char);
                    text.push_str(&count.text());
                }
            }
            info.push((parent_field.vcf_name.clone(), InfoValue::Text(text)));
        }

        // FORMAT.
        let mut format: Vec<(String, FormatValue)> = Vec::new();
        let samples = calls.len();
        let mut gq_zero: Vec<usize> = Vec::new();
        let mut pl_zero: Vec<usize> = Vec::new();
        let mut min_dp: Option<Vec<i32>> = None;
        let mut dp_format: Option<Vec<i32>> = None;
        let mut dp_info: Option<Vec<i32>> = None;
        let slots: Vec<usize> = self
            .format_slots
            .iter()
            .copied()
            .chain(self.dp_info_slot)
            .collect();
        for slot in slots {
            let field = config.field(slot);
            if field.len() == Len::G && too_many_alts {
                continue;
            }
            let known = config.slots[slot].known;
            let values: Vec<Option<&Data>> = calls
                .iter()
                .enumerate()
                .map(|(index, call)| {
                    if !call.valid {
                        None
                    } else if source(slot) {
                        remapped[index][slot].as_ref()
                    } else {
                        call.fields[slot].as_ref()
                    }
                })
                .collect();
            let Some(column) = collect_and_extend(&values, known == Known::Gt) else {
                continue;
            };
            let per_sample = column_length(&column) / samples.max(1);
            let mut insert = true;
            let column = match (known, column) {
                (Known::Gt, FormatValue::Ints(values)) => {
                    let ploidy = field.len().ploidy(per_sample);
                    let phased = field.len() == Len::PP;
                    let mut encoded = Vec::with_capacity(ploidy * samples);
                    for sample in values.chunks(per_sample.max(1)) {
                        encode_gt(sample, phased, self.options.produce_gt, &mut encoded);
                    }
                    FormatValue::Gt(encoded)
                }
                (Known::Gq, FormatValue::Ints(values)) => {
                    gq_zero = (0..samples)
                        .filter(|i| values.get(*i) == Some(&0))
                        .collect();
                    FormatValue::Ints(values)
                }
                (Known::MinDp, FormatValue::Ints(values)) => {
                    min_dp = Some(values[..samples.min(values.len())].to_vec());
                    FormatValue::Ints(values)
                }
                (Known::DpFormat, FormatValue::Ints(values)) => {
                    dp_format = Some(values[..samples.min(values.len())].to_vec());
                    insert = false;
                    FormatValue::Ints(values)
                }
                (Known::Dp, FormatValue::Ints(values)) => {
                    dp_info = Some(values[..samples.min(values.len())].to_vec());
                    insert = false;
                    FormatValue::Ints(values)
                }
                (Known::Pl, FormatValue::Ints(values)) => {
                    pl_zero = (0..samples)
                        .filter(|i| values.get(i * per_sample) == Some(&0))
                        .collect();
                    FormatValue::Ints(values)
                }
                (_, column) => column,
            };
            if insert {
                let key = field.vcf_name.clone();
                if known == Known::Gt {
                    format.insert(0, (key, column));
                } else {
                    format.push((key, column));
                }
            }
        }
        if self.options.produce_gt && !gq_zero.is_empty() && !pl_zero.is_empty() {
            if let Some((_, FormatValue::Gt(gt))) = format.first_mut() {
                let ploidy = gt.len() / samples.max(1);
                for sample in gq_zero.iter().filter(|sample| pl_zero.contains(sample)) {
                    for value in gt[sample * ploidy..(sample + 1) * ploidy].iter_mut() {
                        if *value == INT32_VECTOR_END {
                            break;
                        }
                        *value = 0;
                    }
                }
            }
        }
        if dp_info.is_some() || dp_format.is_some() {
            let mut sum: i64 = 0;
            let mut any_format = false;
            let mut format_values = vec![INT32_MISSING; samples];
            for (j, slot) in format_values.iter_mut().enumerate() {
                let mut info_value = dp_info.as_ref().map_or(INT32_MISSING, |values| values[j]);
                let format_value = dp_format.as_ref().map_or(INT32_MISSING, |values| values[j]);
                if !int_valid(info_value) {
                    info_value = match min_dp.as_ref().map(|values| values[j]) {
                        Some(value) if int_valid(value) => value,
                        _ => format_value,
                    };
                }
                *slot = format_value;
                any_format |= int_valid(format_value);
                if int_valid(info_value) {
                    sum += info_value as i64;
                }
            }
            if any_format {
                format.push(("DP".to_string(), FormatValue::Ints(format_values)));
            }
            if sum > 0 && !reference_block_only {
                info.push(("DP".to_string(), InfoValue::Ints(vec![sum])));
            }
        }

        if info.is_empty() {
            out.push('.');
        }
        for (index, (key, value)) in info.iter().enumerate() {
            if index > 0 {
                out.push(';');
            }
            out.push_str(key);
            out.push('=');
            match value {
                InfoValue::Ints(values) => push_ints64(values, out),
                InfoValue::Floats(values) => push_floats(values, out),
                InfoValue::Text(text) => out.push_str(text),
            }
        }
        if !format.is_empty() {
            out.push('\t');
            out.push_str(
                &format
                    .iter()
                    .map(|(key, _)| key.as_str())
                    .collect::<Vec<_>>()
                    .join(":"),
            );
            for sample in 0..samples {
                out.push('\t');
                for (index, (_, column)) in format.iter().enumerate() {
                    if index > 0 {
                        out.push(':');
                    }
                    let per_sample = column_length(column) / samples.max(1);
                    let range = sample * per_sample..(sample + 1) * per_sample;
                    match column {
                        FormatValue::Gt(values) => push_gt(&values[range], out),
                        FormatValue::Ints(values) => {
                            push_ints(&values[range], out);
                        }
                        FormatValue::Floats(values) => push_floats(&values[range], out),
                        FormatValue::Text(values) => push_text(&values[range], out),
                    }
                }
            }
        } else if samples > 0 {
            out.push_str("\t.");
            for _ in 0..samples {
                out.push_str("\t.");
            }
        }
        out.push('\n');
        let _ = line_start;
        Ok(())
    }

    /// `GA4GHOperator::remap_if_needed`.
    fn remap_field(
        &self,
        data: Option<&Data>,
        field: &VidField,
        lut: &Lut,
        merged_alleles: usize,
        non_ref_exists: bool,
        ploidy: usize,
    ) -> Option<Data> {
        let data = data?;
        let num_alt = merged_alleles - 1;
        if field.len() == Len::G
            && genotype_count(num_alt, ploidy) > self.options.max_genotype_count as u128
        {
            return None;
        }
        let non_ref = non_ref_exists && field.remap_missing_with_non_ref;
        Some(match data {
            Data::Multi(entries) => Data::Multi(remap_multi(
                entries,
                lut,
                merged_alleles,
                non_ref,
                field.len(),
            )),
            _ => remap_vector(data, lut, merged_alleles, non_ref, ploidy, field.len()),
        })
    }
}

#[derive(Debug, Clone, Copy)]
enum HistogramCount {
    Int(i64),
    Float(f32),
}

impl HistogramCount {
    fn add(self, other: HistogramCount) -> HistogramCount {
        match (self, other) {
            (HistogramCount::Int(a), HistogramCount::Int(b)) => HistogramCount::Int(a + b),
            (HistogramCount::Float(a), HistogramCount::Float(b)) => HistogramCount::Float(a + b),
            (a, _) => a,
        }
    }

    fn text(self) -> String {
        match self {
            HistogramCount::Int(value) => value.to_string(),
            HistogramCount::Float(value) => format!("{:.3}", value as f64),
        }
    }
}

/// A valid element of a multi-D field's bytes as a number.
fn element_value(bits: u32, elem: Elem) -> Option<f64> {
    match elem {
        Elem::Float => float_valid(bits).then(|| f32::from_bits(bits) as f64),
        _ => {
            let value = bits as i32;
            int_valid(value).then_some(value as f64)
        }
    }
}

/// `get_valid_sum` and `get_valid_mean` over the first element of each call's field.
fn sum_or_mean(rows: &[Option<&Data>], field: &VidField) -> Option<InfoValue> {
    let mut count = 0usize;
    match field.elem() {
        Elem::Int => {
            let mut sum: i64 = 0;
            for data in rows.iter().flatten() {
                if let Data::Int(values) = data {
                    if let Some(value) = values.first().filter(|value| int_valid(**value)) {
                        sum += *value as i64;
                        count += 1;
                    }
                }
            }
            if count == 0 {
                return None;
            }
            Some(InfoValue::Ints(vec![if field.op == Op::Mean {
                (sum / count as i64) as i32 as i64
            } else {
                sum
            }]))
        }
        Elem::Float => {
            let mut sum: f32 = 0.0;
            for data in rows.iter().flatten() {
                if let Data::Float(values) = data {
                    if let Some(value) = values.first().filter(|value| float_valid(**value)) {
                        sum += f32::from_bits(*value);
                        count += 1;
                    }
                }
            }
            if count == 0 {
                return None;
            }
            let value = if field.op == Op::Mean {
                sum / count as f32
            } else {
                sum
            };
            Some(InfoValue::Floats(vec![value.to_bits()]))
        }
        Elem::Char => None,
    }
}

/// `get_valid_median`: `std::nth_element` at the middle of the valid first elements, which is
/// the order statistic of that rank.
fn median(rows: &[Option<&Data>]) -> Option<InfoValue> {
    let mut ints: Vec<i32> = Vec::new();
    let mut floats = Vec::new();
    for data in rows.iter().flatten() {
        match data {
            Data::Int(values) => ints.extend(values.first().filter(|value| int_valid(**value))),
            Data::Float(values) => floats.extend(
                values
                    .first()
                    .filter(|value| float_valid(**value))
                    .map(|bits| f32::from_bits(*bits)),
            ),
            _ => {}
        }
    }
    if !ints.is_empty() {
        ints.sort_unstable();
        return Some(InfoValue::Ints(vec![ints[ints.len() / 2] as i64]));
    }
    if !floats.is_empty() {
        floats.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return Some(InfoValue::Floats(vec![floats[floats.len() / 2].to_bits()]));
    }
    None
}

fn concatenate(rows: &[Option<&Data>]) -> Option<InfoValue> {
    let mut ints = Vec::new();
    let mut floats = Vec::new();
    let mut text = Vec::new();
    for data in rows.iter().flatten() {
        match data {
            Data::Int(values) => ints.extend(values.iter().map(|value| *value as i64)),
            Data::Float(values) => floats.extend(values),
            Data::Str(bytes) => text.extend(bytes),
            _ => {}
        }
    }
    if !ints.is_empty() {
        Some(InfoValue::Ints(ints))
    } else if !floats.is_empty() {
        Some(InfoValue::Floats(floats))
    } else if !text.is_empty() {
        Some(InfoValue::Text(String::from_utf8_lossy(&text).to_string()))
    } else {
        None
    }
}

/// `compute_valid_element_wise_sum` over the calls of a 1-D field.
fn element_wise_sum(inputs: &[Option<Data>]) -> Option<InfoValue> {
    let mut ints: Vec<i64> = Vec::new();
    let mut floats: Vec<u32> = Vec::new();
    let mut is_float = false;
    for data in inputs.iter().flatten() {
        match data {
            Data::Int(values) => {
                let previous = ints.len();
                if values.len() > ints.len() {
                    ints.resize(values.len(), INT64_MISSING);
                }
                for (i, value) in values.iter().enumerate() {
                    let old = ints[i];
                    let new_valid = int_valid(*value);
                    ints[i] = if i < previous {
                        let old_valid = old != INT64_MISSING && old != INT64_VECTOR_END;
                        match (old_valid, new_valid) {
                            (true, true) => old + *value as i64,
                            (_, true) => *value as i64,
                            _ => old,
                        }
                    } else if *value == INT32_MISSING {
                        INT64_MISSING
                    } else if *value == INT32_VECTOR_END {
                        INT64_VECTOR_END
                    } else {
                        *value as i64
                    };
                }
            }
            Data::Float(values) => {
                is_float = true;
                let previous = floats.len();
                if values.len() > floats.len() {
                    floats.resize(values.len(), FLOAT_MISSING);
                }
                for (i, value) in values.iter().enumerate() {
                    let old = floats[i];
                    floats[i] = if i < previous {
                        match (float_valid(old), float_valid(*value)) {
                            (true, true) => {
                                (f32::from_bits(old) + f32::from_bits(*value)).to_bits()
                            }
                            (_, true) => *value,
                            _ => old,
                        }
                    } else {
                        *value
                    };
                }
            }
            _ => {}
        }
    }
    if is_float {
        (!floats.is_empty()).then_some(InfoValue::Floats(floats))
    } else {
        (!ints.is_empty()).then_some(InfoValue::Ints(ints))
    }
}

/// `compute_valid_element_wise_sum_2D_vector` then `stringify_2D_vector`.
fn element_wise_sum_2d(inputs: &[Option<Data>], field: &VidField) -> Option<InfoValue> {
    let float = field.elem() == Elem::Float;
    let mut sums: Vec<Vec<Option<f64>>> = Vec::new();
    let mut valid = false;
    for data in inputs.iter().flatten() {
        let Data::Multi(entries) = data else { continue };
        if entries.len() > sums.len() {
            sums.resize(entries.len(), Vec::new());
        }
        for (entry, bytes) in entries.iter().enumerate() {
            let values: Vec<u32> = bytes
                .chunks_exact(4)
                .map(|chunk| u32::from_le_bytes(chunk.try_into().expect("4 bytes")))
                .collect();
            if values.len() > sums[entry].len() {
                sums[entry].resize(values.len(), None);
            }
            for (i, bits) in values.iter().enumerate() {
                let Some(new) = element_value(*bits, field.elem()) else {
                    continue;
                };
                valid = true;
                sums[entry][i] = Some(match sums[entry][i] {
                    Some(old) if float => (old as f32 + new as f32) as f64,
                    Some(old) => old + new,
                    None => new,
                });
            }
        }
    }
    if !valid {
        return None;
    }
    let mut text = String::new();
    for (entry, values) in sums.iter().enumerate() {
        if entry > 0 {
            text.push(field.delimiter(0) as char);
        }
        for (i, value) in values.iter().enumerate() {
            if i > 0 {
                text.push(field.delimiter(1) as char);
            }
            if let Some(value) = value {
                if float {
                    text.push_str(&format!("{:.3}", *value as f32 as f64));
                } else {
                    text.push_str(&format!("{}", *value as i64));
                }
            }
        }
    }
    Some(InfoValue::Text(text))
}

fn column_length(column: &FormatValue) -> usize {
    match column {
        FormatValue::Gt(values) | FormatValue::Ints(values) => values.len(),
        FormatValue::Floats(values) => values.len(),
        FormatValue::Text(values) => values.len(),
    }
}

/// `collect_and_extend_fields` with VCF output's vector-end padding: `None` when no call has
/// the field.
fn collect_and_extend(values: &[Option<&Data>], gt: bool) -> Option<FormatValue> {
    let width = values.iter().flatten().map(|data| data.length()).max()?;
    let kind = values.iter().flatten().next()?;
    Some(match kind {
        Data::Float(_) => {
            let mut out = Vec::with_capacity(width * values.len());
            for value in values {
                let start = out.len();
                if let Some(Data::Float(items)) = value {
                    out.extend_from_slice(items);
                }
                if out.len() == start {
                    out.push(FLOAT_MISSING);
                }
                out.resize(start + width, FLOAT_VECTOR_END);
            }
            FormatValue::Floats(out)
        }
        Data::Str(_) => {
            let mut out = Vec::with_capacity(width * values.len());
            for value in values {
                let start = out.len();
                if let Some(Data::Str(items)) = value {
                    out.extend_from_slice(items);
                }
                if out.len() == start {
                    out.push(STR_MISSING);
                }
                out.resize(start + width, STR_VECTOR_END);
            }
            FormatValue::Text(out)
        }
        _ => {
            let mut out = Vec::with_capacity(width * values.len());
            for value in values {
                let start = out.len();
                if let Some(Data::Int(items)) = value {
                    out.extend_from_slice(items);
                }
                if out.len() == start {
                    out.push(if gt { INT32_VECTOR_END } else { INT32_MISSING });
                }
                out.resize(start + width, INT32_VECTOR_END);
            }
            FormatValue::Ints(out)
        }
    })
}

/// `encode_GT_vector` for one sample: htslib's genotype encoding, phase stored between alleles
/// when the length is `PP`.
fn encode_gt(sample: &[i32], phased_ploidy: bool, produce: bool, out: &mut Vec<i32>) {
    let element = |value: i32, is_phased: bool, with_phase: bool| -> i32 {
        if !int_valid(value) {
            return value;
        }
        match (produce, with_phase && is_phased) {
            (true, true) => ((value + 1) << 1) | 1,
            (true, false) => (value + 1) << 1,
            (false, true) => 1,
            (false, false) => 0,
        }
    };
    if sample.is_empty() {
        return;
    }
    out.push(element(sample[0], false, false));
    if phased_ploidy {
        let mut k = 2;
        while k < sample.len() {
            out.push(element(sample[k], sample[k - 1] > 0, true));
            k += 2;
        }
    } else {
        for value in &sample[1..] {
            out.push(element(*value, false, false));
        }
    }
}

// ---------------------------------------------------------------------------------------------
// htslib's vcf_format
// ---------------------------------------------------------------------------------------------

fn push_ints(values: &[i32], out: &mut String) {
    if values.is_empty() {
        out.push('.');
        return;
    }
    for (j, value) in values.iter().enumerate() {
        if *value == INT32_VECTOR_END {
            break;
        }
        if j > 0 {
            out.push(',');
        }
        if *value == INT32_MISSING {
            out.push('.');
        } else {
            out.push_str(&value.to_string());
        }
    }
}

fn push_ints64(values: &[i64], out: &mut String) {
    for (j, value) in values.iter().enumerate() {
        if *value == INT64_VECTOR_END {
            break;
        }
        if j > 0 {
            out.push(',');
        }
        if *value == INT64_MISSING {
            out.push('.');
        } else {
            out.push_str(&value.to_string());
        }
    }
}

fn push_floats(values: &[u32], out: &mut String) {
    if values.is_empty() {
        out.push('.');
        return;
    }
    for (j, bits) in values.iter().enumerate() {
        if *bits == FLOAT_VECTOR_END {
            break;
        }
        if j > 0 {
            out.push(',');
        }
        if *bits == FLOAT_MISSING {
            out.push('.');
        } else {
            kputd(f32::from_bits(*bits) as f64, out);
        }
    }
}

fn push_text(values: &[u8], out: &mut String) {
    if values.is_empty() {
        out.push('.');
        return;
    }
    for byte in values {
        if *byte == STR_VECTOR_END {
            break;
        }
        out.push(if *byte == STR_MISSING {
            '.'
        } else {
            *byte as char
        });
    }
}

/// `bcf_format_gt`.
fn push_gt(values: &[i32], out: &mut String) {
    let mut written = 0;
    for value in values {
        if *value == INT32_VECTOR_END {
            break;
        }
        if written > 0 {
            out.push(if value & 1 == 1 { '|' } else { '/' });
        }
        if value >> 1 == 0 {
            out.push('.');
        } else {
            out.push_str(&((value >> 1) - 1).to_string());
        }
        written += 1;
    }
    if written == 0 {
        out.push('.');
    }
}

/// htslib's `kputd`: six significant digits, trailing zeros culled, `%#g` outside
/// [0.0001, 999999].
fn kputd(d: f64, out: &mut String) {
    if d == 0.0 {
        out.push_str(if d.is_sign_negative() { "-0" } else { "0" });
        return;
    }
    let mut d = d;
    if d < 0.0 {
        out.push('-');
        d = -d;
    }
    if !(0.0001..=999_999.0).contains(&d) {
        out.push_str(&format_sharp_g(d));
        return;
    }
    let mut i = (d * 10_000_000_000.0) as u64;
    i += if d < 0.001 {
        5
    } else if d < 0.01 {
        50
    } else if d < 0.1 {
        500
    } else if d < 1.0 {
        5000
    } else if d < 10.0 {
        50_000
    } else if d < 100.0 {
        500_000
    } else if d < 1000.0 {
        5_000_000
    } else if d < 10_000.0 {
        50_000_000
    } else if d < 100_000.0 {
        500_000_000
    } else {
        5_000_000_000
    };
    // buf[0..21], cp walking down from 20.
    let mut buf = [0u8; 21];
    let mut cp = 20usize;
    loop {
        cp -= 1;
        buf[cp] = b'0' + (i % 10) as u8;
        i /= 10;
        if i < 1 {
            break;
        }
    }
    buf[20] = 0;
    let mut p = 20 - cp;
    let mut ep;
    if p <= 10 {
        buf[cp + 6] = 0;
        ep = cp + 5;
        while p < 10 {
            cp -= 1;
            buf[cp] = b'0';
            p += 1;
        }
        cp -= 1;
        buf[cp] = b'.';
        cp -= 1;
        buf[cp] = b'0';
    } else {
        cp -= 1;
        let mut xp = cp;
        while p > 10 {
            buf[xp] = buf[xp + 1];
            p -= 1;
            xp += 1;
        }
        buf[xp] = b'.';
        buf[cp + 7] = 0;
        ep = cp + 6;
    }
    while buf[ep] == b'0' && ep > cp {
        ep -= 1;
    }
    let z = ep + 1;
    while ep > cp {
        if buf[ep] == b'.' {
            if buf[z - 1] == b'.' {
                buf[z] = b'0';
                buf[z + 1] = 0;
            } else {
                buf[z] = 0;
            }
            break;
        }
        ep -= 1;
    }
    let length = buf[cp..]
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(21 - cp);
    out.push_str(std::str::from_utf8(&buf[cp..cp + length]).unwrap_or(""));
}

/// C's `%#g` of a positive finite double.
fn format_sharp_g(d: f64) -> String {
    if d.is_nan() {
        return "nan".to_string();
    }
    if d.is_infinite() {
        return "inf".to_string();
    }
    let scientific = format!("{d:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').expect("an exponent");
    let exponent: i32 = exponent.parse().expect("a decimal exponent");
    if (-4..6).contains(&exponent) {
        let decimals = (5 - exponent) as usize;
        format!("{d:.decimals$}")
    } else {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exponent.abs())
    }
}

// ---------------------------------------------------------------------------------------------
// libstdc++'s unordered containers
// ---------------------------------------------------------------------------------------------

/// `std::_Hash_bytes` (64-bit MurmurHash2) with `std::hash<std::string>`'s seed.
fn hash_bytes(bytes: &[u8]) -> u64 {
    const MUL: u64 = (0xc6a4_a793u64 << 32) + 0x5bd1_e995;
    let shift_mix = |v: u64| v ^ (v >> 47);
    let mut hash = 0xc70f_6907u64 ^ (bytes.len() as u64).wrapping_mul(MUL);
    let aligned = bytes.len() & !7;
    for chunk in bytes[..aligned].chunks_exact(8) {
        let data =
            shift_mix(u64::from_le_bytes(chunk.try_into().expect("8 bytes")).wrapping_mul(MUL))
                .wrapping_mul(MUL);
        hash ^= data;
        hash = hash.wrapping_mul(MUL);
    }
    if bytes.len() & 7 != 0 {
        let mut data = 0u64;
        for byte in bytes[aligned..].iter().rev() {
            data = (data << 8) + *byte as u64;
        }
        hash ^= data;
        hash = hash.wrapping_mul(MUL);
    }
    hash = shift_mix(hash).wrapping_mul(MUL);
    shift_mix(hash)
}

/// `_Prime_rehash_policy::_M_next_bkt` (GCC 7): the bucket count and the next resize.
fn next_bucket_count(n: u64) -> Result<(u64, u64), Reason> {
    const FAST: [u64; 14] = [2, 2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11, 13, 13];
    // shortcut: libstdc++'s __prime_list only as far as 199, which holds 199 keys; a set larger
    // than that would need the rest of the list.
    const PRIMES: [u64; 31] = [
        17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 103, 109, 113,
        127, 137, 139, 149, 157, 167, 179, 193, 199,
    ];
    if n < FAST.len() as u64 {
        return Ok((FAST[n as usize], FAST[n as usize]));
    }
    PRIMES
        .iter()
        .find(|prime| **prime >= n)
        .map(|prime| (*prime, *prime))
        .ok_or_else(|| {
            "a record merges more IDs or histograms than this port orders as libstdc++ does."
                .to_string()
        })
}

/// The iteration order of a libstdc++ `unordered_set` (or map) after inserting keys with these
/// distinct hash codes in this order: indices into `codes`.
fn unordered_order(codes: &[u64]) -> Result<Vec<usize>, Reason> {
    const NONE: usize = usize::MAX;
    const HEAD: usize = usize::MAX - 1;
    let mut next = vec![NONE; codes.len()];
    let mut head = NONE;
    let mut buckets: Vec<usize> = vec![NONE];
    let mut bucket_count = 1u64;
    let mut next_resize = 0u64;
    let link_after = |next: &mut Vec<usize>, head: &mut usize, before: usize, node: usize| {
        if before == HEAD {
            next[node] = *head;
            *head = node;
        } else {
            next[node] = next[before];
            next[before] = node;
        }
    };
    for (node, code) in codes.iter().enumerate() {
        let elements = node as u64;
        if elements + 1 > next_resize {
            let minimum = (elements + 1) as f64;
            if minimum >= bucket_count as f64 {
                let (count, resize) = next_bucket_count((elements + 2).max(bucket_count * 2))?;
                next_resize = resize;
                // _M_rehash_aux for unique keys.
                let mut new_buckets = vec![NONE; count as usize];
                let mut p = head;
                head = NONE;
                let mut begin_bucket = 0usize;
                while p != NONE {
                    let following = next[p];
                    let bucket = (codes[p] % count) as usize;
                    if new_buckets[bucket] == NONE {
                        next[p] = head;
                        head = p;
                        new_buckets[bucket] = HEAD;
                        if next[p] != NONE {
                            new_buckets[begin_bucket] = p;
                        }
                        begin_bucket = bucket;
                    } else {
                        let before = new_buckets[bucket];
                        link_after(&mut next, &mut head, before, p);
                    }
                    p = following;
                }
                buckets = new_buckets;
                bucket_count = count;
            } else {
                next_resize = bucket_count;
            }
        }
        // _M_insert_bucket_begin.
        let bucket = (code % bucket_count) as usize;
        if buckets[bucket] != NONE {
            let before = buckets[bucket];
            link_after(&mut next, &mut head, before, node);
        } else {
            next[node] = head;
            head = node;
            if next[node] != NONE {
                buckets[(codes[next[node]] % bucket_count) as usize] = node;
            }
            buckets[bucket] = HEAD;
        }
    }
    let mut order = Vec::with_capacity(codes.len());
    let mut p = head;
    while p != NONE {
        order.push(p);
        p = next[p];
    }
    Ok(order)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kputd_text(d: f64) -> String {
        let mut out = String::new();
        kputd(d, &mut out);
        out
    }

    #[test]
    fn kputd_matches_htslib() {
        assert_eq!(kputd_text(3.0103f32 as f64), "3.0103");
        assert_eq!(kputd_text(-0.524f32 as f64), "-0.524");
        assert_eq!(kputd_text(0.0), "0");
        assert_eq!(kputd_text(1.253f32 as f64), "1.253");
        assert_eq!(kputd_text(123456.7), "123457.0");
        assert_eq!(kputd_text(0.00001), "1.00000e-05");
        assert_eq!(kputd_text(1234567.0), "1.23457e+06");
    }

    #[test]
    fn unordered_iteration_is_libstdcxx() {
        // The measured order of the histogram annotations, keyed by vid index.
        assert_eq!(unordered_order(&[15, 19, 22]).unwrap(), vec![2, 0, 1]);
        // `rs100` then `rs102` inserted; the reference prints `rs102;rs100`.
        let codes = [hash_bytes(b"rs100"), hash_bytes(b"rs102")];
        assert_eq!(unordered_order(&codes).unwrap(), vec![1, 0]);
    }

    #[test]
    fn the_bucket_count_follows_the_prime_list() {
        assert_eq!(next_bucket_count(13).unwrap(), (13, 13));
        assert_eq!(next_bucket_count(130).unwrap(), (137, 137));
        assert_eq!(next_bucket_count(194).unwrap(), (199, 199));
        assert!(next_bucket_count(200).is_err());
        // 199 keys fit, the 200th does not.
        assert!(unordered_order(&(0..199).collect::<Vec<u64>>()).is_ok());
        assert!(unordered_order(&(0..200).collect::<Vec<u64>>()).is_err());
    }

    #[test]
    fn genotypes_are_in_vcf_order() {
        let order = genotypes(3, 2);
        for (index, genotype) in order.iter().enumerate() {
            assert_eq!(genotype_index(genotype), index as u128);
        }
        assert_eq!(order[3], vec![0, 2]);
    }
}
