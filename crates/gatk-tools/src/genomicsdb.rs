//! The files `GenomicsDBImport` writes into a workspace, as `org.genomicsdb` 1.5 writes them.
//!
//! The reference hands every record to a native library as a VCF stream: htsjdk's `VCFEncoder`,
//! with QUAL at two decimals and the sample column as the raw text it read,
//! writes the `VariantContext`, a fork of htslib parses it, `VCF2Binary` turns it into a cell,
//! `LoaderArrayWriter` puts the cells in column-major order and TileDB writes them as compressed
//! tiles. Nothing in those files depends on the run but the names of the fragments and of the
//! metadata files, which hold a random UUID, a thread id and the wall clock. Two runs of the
//! reference differ in those names and in the `Date` of the command line in `vcfheader.vcf`, and
//! in no other byte.
//!
//! Each conversion below is therefore the composition of those stages, and each is named after
//! the stage it stands for. What the stream cannot carry, or what makes the native code throw, is
//! returned as a reason rather than guessed at, and the runner refuses with it.

use std::collections::{BTreeMap, HashMap};
use std::rc::Rc;

use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

use crate::variant_recalibrator_model::java_parse_double;

/// What the reference does that this port does not carry, in a sentence.
pub type Reason = String;

/// `TILEDB_EMPTY_INT32`, the cell value of a fixed-length field the record does not have.
const NULL_INT32: i32 = i32::MAX;
/// `TILEDB_EMPTY_FLOAT32`.
const NULL_FLOAT32: u32 = 0x7f7f_ffff;
/// `TILEDB_EMPTY_CHAR`.
const NULL_CHAR: u8 = 0x7f;
/// `bcf_int32_missing`, what a `.` inside a vector decodes to.
const BCF_INT32_MISSING: i32 = i32::MIN;
/// `bcf_float_missing`.
const BCF_FLOAT_MISSING: u32 = 0x7f80_0001;
/// `TILEDB_VAR_NUM`, the cell value count of a variable-length attribute.
const VAR_NUM: i32 = i32::MAX;
/// `GenomicsDBImportConfig`'s default, which `GenomicsDBImport` does not change.
const CELLS_PER_TILE: usize = 1000;

// ---------------------------------------------------------------------------------------------
// vidmap.json and callset.json
// ---------------------------------------------------------------------------------------------

/// `GenomicsDBVidMapProto.GenomicsDBFieldInfo`, holding only what a reference run can set.
#[derive(Debug, Clone, PartialEq)]
pub struct Field {
    pub name: String,
    pub types: Vec<String>,
    pub classes: Vec<String>,
    pub lengths: Vec<String>,
    pub delimiters: Vec<String>,
    pub combine: Option<String>,
    pub disable_remap: bool,
}

/// `GenomicsDBVidMapProto.Chromosome`.
#[derive(Debug, Clone, PartialEq)]
pub struct Contig {
    pub name: String,
    pub length: i64,
    pub offset: i64,
}

/// `GenomicsDBVidMapProto.VidMappingPB`.
#[derive(Debug, Clone, PartialEq)]
pub struct VidMap {
    pub fields: Vec<Field>,
    pub contigs: Vec<Contig>,
}

const HISTOGRAM_FIELDS: [&str; 3] = [
    "AS_RAW_BaseQRankSum",
    "AS_RAW_MQRankSum",
    "AS_RAW_ReadPosRankSum",
];
const TWO_DIM_FLOAT_FIELDS: [&str; 1] = ["AS_RAW_MQ"];
const TWO_DIM_INT_FIELDS: [&str; 3] = ["AS_SB_TABLE", "AS_QUALapprox", "AS_VarDP"];

/// A `##FILTER` line in either of the two shapes it has: built in code, or read from a file.
fn filter_id(line: &HeaderLine) -> Option<&str> {
    match line {
        HeaderLine::Filter { id, .. } => Some(id),
        HeaderLine::Structured { key, fields } if key == "FILTER" => fields
            .iter()
            .find(|(field, _)| field == "ID")
            .map(|(_, value)| value.as_str()),
        _ => None,
    }
}

fn line_type_name(line_type: LineType) -> &'static str {
    match line_type {
        LineType::Integer => "Integer",
        LineType::Float => "Float",
        LineType::String => "String",
        LineType::Character => "Character",
        LineType::Flag => "Flag",
    }
}

/// `VidMapExtensions.getVcfHeaderLength`.
fn header_length(number: &Cardinality, line_type: LineType) -> String {
    match number {
        Cardinality::Unbounded => "VAR".to_string(),
        Cardinality::A => "A".to_string(),
        Cardinality::R => "R".to_string(),
        Cardinality::G => "G".to_string(),
        // Weird Flag fields - Number=0 in the VCF header.
        Cardinality::Fixed(0) if line_type == LineType::Flag => "1".to_string(),
        Cardinality::Fixed(count) => count.to_string(),
    }
}

/// `VidMapExtensions.generateVidMapFromMergedHeader`, then
/// `GATKGenomicsDBUtils.updateImportProtobufVidMapping`, which is what the file ends up holding.
pub fn vid_map(lines: &[HeaderLine]) -> Result<VidMap, Reason> {
    let mut fields = vec![Field {
        name: "ID".to_string(),
        types: vec!["char".to_string()],
        classes: Vec::new(),
        lengths: vec!["var".to_string()],
        delimiters: Vec::new(),
        combine: None,
        disable_remap: false,
    }];
    // `checkForDuplicateVcfFieldNames`: per name, per class, the type, the length and the entry.
    let mut seen: HashMap<String, Vec<(String, String, String, usize)>> = HashMap::new();
    let mut contigs = Vec::new();
    let mut offset = 0i64;

    for line in lines {
        match line {
            HeaderLine::Compound {
                key,
                id,
                number,
                line_type,
                ..
            } if key == "FORMAT" => {
                if *line_type == LineType::Character {
                    return Err(format!(
                        "FORMAT field {id} has type Character, which GenomicsDB's vid map does \
                         not name."
                    ));
                }
                let is_gt = id == "GT";
                let genomicsdb_type = if is_gt {
                    "int".to_string()
                } else {
                    line_type_name(*line_type).to_string()
                };
                let length = if is_gt {
                    "PP".to_string()
                } else if *line_type == LineType::String {
                    "var".to_string()
                } else {
                    header_length(number, *line_type)
                };
                if !duplicate_check(
                    id,
                    "FORMAT",
                    &genomicsdb_type,
                    &length,
                    &mut fields,
                    &mut seen,
                )? {
                    fields.push(Field {
                        name: id.clone(),
                        types: vec![genomicsdb_type],
                        classes: vec!["FORMAT".to_string()],
                        lengths: vec![length],
                        delimiters: Vec::new(),
                        combine: None,
                        disable_remap: false,
                    });
                }
            }
            HeaderLine::Compound {
                key,
                id,
                number,
                line_type,
                ..
            } if key == "INFO" => {
                if *line_type == LineType::Character {
                    return Err(format!(
                        "INFO field {id} has type Character, which GenomicsDB's vid map does not \
                         name."
                    ));
                }
                let length = if *line_type == LineType::String {
                    "var".to_string()
                } else {
                    header_length(number, *line_type)
                };
                let type_name = line_type_name(*line_type);
                if duplicate_check(id, "INFO", type_name, &length, &mut fields, &mut seen)? {
                    continue;
                }
                let name = id.as_str();
                let mut field = Field {
                    name: name.to_string(),
                    types: Vec::new(),
                    classes: vec!["INFO".to_string()],
                    lengths: Vec::new(),
                    delimiters: Vec::new(),
                    combine: None,
                    disable_remap: false,
                };
                if HISTOGRAM_FIELDS.contains(&name) {
                    set_allele_specific(&mut field, &["float", "int"], "histogram_sum");
                } else if TWO_DIM_FLOAT_FIELDS.contains(&name) {
                    set_allele_specific(&mut field, &["float"], "element_wise_sum");
                } else if TWO_DIM_INT_FIELDS.contains(&name) {
                    set_allele_specific(&mut field, &["int"], "element_wise_sum");
                } else {
                    field.types.push(type_name.to_string());
                    field.lengths.push(length);
                }
                fields.push(field);
            }
            line if filter_id(line).is_some() => {
                let id = filter_id(line).expect("just tested");
                if !duplicate_check(id, "FILTER", "Integer", "1", &mut fields, &mut seen)? {
                    fields.push(Field {
                        name: id.to_string(),
                        types: vec!["Integer".to_string()],
                        classes: vec!["FILTER".to_string()],
                        lengths: vec!["1".to_string()],
                        delimiters: Vec::new(),
                        combine: None,
                        disable_remap: false,
                    });
                }
            }
            HeaderLine::Contig { fields: pairs, .. } => {
                let name = pairs
                    .iter()
                    .find(|(key, _)| key == "ID")
                    .map(|(_, value)| value.clone())
                    .unwrap_or_default();
                let length = pairs
                    .iter()
                    .find(|(key, _)| key == "length")
                    .and_then(|(_, value)| value.parse::<i64>().ok())
                    .ok_or_else(|| {
                        format!(
                            "##contig line {name} has no length, which GenomicsDB would lay out \
                             from a length this port has not measured."
                        )
                    })?;
                contigs.push(Contig {
                    name,
                    length,
                    offset,
                });
                offset += length;
            }
            _ => {}
        }
    }

    // `updateImportProtobufVidMapping`.
    let combine = |fields: &mut Vec<Field>, name: &str, operation: &str| {
        if let Some(field) = fields.iter_mut().find(|field| field.name == name) {
            field.combine = Some(operation.to_string());
        }
    };
    combine(&mut fields, "RAW_MQandDP", "element_wise_sum");
    if let Some(field) = fields.iter_mut().find(|field| field.name == "AS_RAW_MQ") {
        set_allele_specific(field, &["float"], "element_wise_sum");
    }
    combine(&mut fields, "MQ_DP", "sum");
    combine(&mut fields, "QUALapprox", "sum");
    combine(&mut fields, "VarDP", "sum");
    combine(&mut fields, "RAW_GT_COUNT", "element_wise_sum");
    if let Some(field) = fields
        .iter_mut()
        .find(|field| field.name == "AS_QUALapprox")
    {
        set_allele_specific(field, &["int"], "element_wise_sum");
    }
    for name in ["AS_RAW_MQ", "AS_SB_TABLE"] {
        if let Some(field) = fields.iter_mut().find(|field| field.name == name) {
            field.disable_remap = true;
        }
    }
    combine(&mut fields, "TREE_SCORE", "sum");
    combine(&mut fields, "HAPCOMP", "element_wise_sum");
    combine(&mut fields, "HAPDOM", "element_wise_sum");
    combine(&mut fields, "ASSEMBLED_HAPS", "sum");
    combine(&mut fields, "FILTERED_HAPS", "sum");

    Ok(VidMap { fields, contigs })
}

fn set_allele_specific(field: &mut Field, types: &[&str], operation: &str) {
    field.lengths = vec!["R".to_string(), "var".to_string()];
    field.delimiters = vec!["|".to_string(), ",".to_string()];
    field.types = types.iter().map(|name| name.to_string()).collect();
    field.combine = Some(operation.to_string());
}

/// `checkForDuplicateVcfFieldNames`. True when no new entry is to be made.
fn duplicate_check(
    name: &str,
    class: &str,
    genomicsdb_type: &str,
    length: &str,
    fields: &mut [Field],
    seen: &mut HashMap<String, Vec<(String, String, String, usize)>>,
) -> Result<bool, Reason> {
    let lines = seen.entry(name.to_string()).or_default();
    if lines.iter().any(|(seen_class, ..)| seen_class == class) {
        return Ok(true);
    }
    if let Some((_, other_type, other_length, entry)) = lines.first().cloned() {
        if other_type == genomicsdb_type && other_length == length {
            fields[entry].classes.push(class.to_string());
            lines.push((
                class.to_string(),
                genomicsdb_type.to_string(),
                length.to_string(),
                entry,
            ));
            return Ok(true);
        }
        return Err(format!(
            "header lines of different classes share the name {name} with different types or \
             lengths, which GenomicsDB renames with a vcf_name this port does not carry."
        ));
    }
    lines.push((
        class.to_string(),
        genomicsdb_type.to_string(),
        length.to_string(),
        fields.len(),
    ));
    Ok(false)
}

/// Gson's `toJson` of a string, with the HTML-safe escapes on as `JsonFormat` leaves them.
fn json_string(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 2);
    out.push('"');
    for character in text.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\u{c}' => out.push_str("\\f"),
            '<' | '>' | '&' | '=' | '\'' | '\u{2028}' | '\u{2029}' => {
                out.push_str(&format!("\\u{:04x}", character as u32))
            }
            control if (control as u32) < 0x20 => {
                out.push_str(&format!("\\u{:04x}", control as u32))
            }
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

fn json_string_list(items: &[String]) -> String {
    let rendered: Vec<String> = items.iter().map(|item| json_string(item)).collect();
    format!("[{}]", rendered.join(", "))
}

impl VidMap {
    /// `JsonFormat.printer().print`, in the field order of the message's descriptor.
    pub fn to_json(&self) -> String {
        let mut entries = Vec::with_capacity(self.fields.len());
        for field in &self.fields {
            let mut members = vec![format!("    \"name\": {}", json_string(&field.name))];
            members.push(format!("    \"type\": {}", json_string_list(&field.types)));
            if !field.classes.is_empty() {
                members.push(format!(
                    "    \"vcf_field_class\": {}",
                    json_string_list(&field.classes)
                ));
            }
            let lengths: Vec<String> = field
                .lengths
                .iter()
                .map(|length| {
                    format!(
                        "{{\n      \"variable_length_descriptor\": {}\n    }}",
                        json_string(length)
                    )
                })
                .collect();
            members.push(format!("    \"length\": [{}]", lengths.join(", ")));
            if !field.delimiters.is_empty() {
                members.push(format!(
                    "    \"vcf_delimiter\": {}",
                    json_string_list(&field.delimiters)
                ));
            }
            if let Some(combine) = &field.combine {
                members.push(format!(
                    "    \"VCF_field_combine_operation\": {}",
                    json_string(combine)
                ));
            }
            if field.disable_remap {
                members.push("    \"disable_remap_missing_with_non_ref\": true".to_string());
            }
            entries.push(format!("{{\n{}\n  }}", members.join(",\n")));
        }
        let contigs: Vec<String> = self
            .contigs
            .iter()
            .map(|contig| {
                format!(
                    "{{\n    \"name\": {},\n    \"length\": \"{}\",\n    \
                     \"tiledb_column_offset\": \"{}\"\n  }}",
                    json_string(&contig.name),
                    contig.length,
                    contig.offset
                )
            })
            .collect();
        format!(
            "{{\n  \"fields\": [{}],\n  \"contigs\": [{}]\n}}",
            entries.join(", "),
            contigs.join(", ")
        )
    }
}

/// `GenomicsDBCallsetsMapProto`, one row per sample in the sorted order of the sample map.
pub fn callset_json(samples: &[String]) -> String {
    let entries: Vec<String> = samples
        .iter()
        .enumerate()
        .map(|(row, sample)| {
            format!(
                "{{\n    \"sample_name\": {},\n    \"row_idx\": \"{row}\",\n    \
                 \"idx_in_file\": \"0\",\n    \"stream_name\": {}\n  }}",
                json_string(sample),
                json_string(&format!("{sample}_stream"))
            )
        })
        .collect();
    format!("{{\n  \"callsets\": [{}]\n}}", entries.join(", "))
}

/// `genomicsdb_column_bounds.json`, as RapidJSON's pretty writer leaves it.
pub fn column_bounds_json(begin: i64, end: i64) -> String {
    format!("{{\n    \"min_column\": {begin},\n    \"max_column\": {end}\n}}")
}

/// `genomicsdb_meta_<uuid>.json`.
pub fn row_bounds_json(lower: i64, upper: i64) -> String {
    format!("{{\n    \"lb_row_idx\": {lower},\n    \"max_valid_row_idx_in_array\": {upper}\n}}")
}

// ---------------------------------------------------------------------------------------------
// __array_schema.tdb
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
enum CellType {
    Int32,
    Int64,
    Float32,
    Char,
}

impl CellType {
    fn code(self) -> u8 {
        match self {
            CellType::Int32 => 0,
            CellType::Int64 => 1,
            CellType::Float32 => 2,
            CellType::Char => 4,
        }
    }
}

#[derive(Debug, Clone)]
enum Source {
    End,
    Ref,
    Alt,
    Id,
    Qual,
    Filter,
    Info {
        vcf_name: String,
        /// The multi-D layout's delimiters are all `|` then `,`; the element types and which one
        /// this attribute is.
        multi_d: Option<(Vec<CellType>, usize)>,
    },
    Format {
        vcf_name: String,
    },
}

#[derive(Debug, Clone)]
struct Attribute {
    name: String,
    cell_type: CellType,
    /// `None` is `TILEDB_VAR_NUM`.
    values: Option<i32>,
    source: Source,
}

impl Attribute {
    fn is_var(&self) -> bool {
        self.values.is_none()
    }
}

/// `VidMapper`'s field table as it stands after `flatten_field`, reduced to what the schema and
/// the cells read: the attributes in `build_tiledb_array_schema`'s order, and the global field
/// index the `FILTER` attribute stores.
#[derive(Debug, Clone)]
pub struct Schema {
    attributes: Vec<Attribute>,
    filter_index: HashMap<String, i32>,
    info: HashMap<String, (LineType, Cardinality)>,
    format: HashMap<String, (LineType, Cardinality)>,
    filters: Vec<String>,
}

fn genomicsdb_cell_type(name: &str) -> Result<CellType, Reason> {
    match name {
        "int" | "Integer" => Ok(CellType::Int32),
        "float" | "Float" => Ok(CellType::Float32),
        "char" | "String" | "Flag" => Ok(CellType::Char),
        other => Err(format!(
            "vid field type {other} is not one this port lays out."
        )),
    }
}

fn fixed_length(lengths: &[String]) -> Option<i32> {
    if lengths.len() == 1 {
        lengths[0].parse::<i32>().ok()
    } else {
        None
    }
}

impl Schema {
    /// `VidMapper::parse_infofields_from_vidmap` then `build_tiledb_array_schema`.
    pub fn new(vid: &VidMap, lines: &[HeaderLine]) -> Result<Self, Reason> {
        let mut info = HashMap::new();
        let mut format = HashMap::new();
        let mut filters = Vec::new();
        for line in lines {
            match line {
                HeaderLine::Compound {
                    key,
                    id,
                    number,
                    line_type,
                    ..
                } => {
                    let table = if key == "INFO" {
                        &mut info
                    } else if key == "FORMAT" {
                        &mut format
                    } else {
                        continue;
                    };
                    table.entry(id.clone()).or_insert((*line_type, *number));
                }
                line => filters.extend(filter_id(line).map(str::to_string)),
            }
        }

        // The global index counts what `flatten_field` inserts after each field.
        let mut filter_index = HashMap::new();
        let mut next_index = 0i32;
        let mut info_attributes = Vec::new();
        let mut format_attributes = Vec::new();
        for field in &vid.fields {
            let index = next_index;
            next_index += 1;
            let is_info = field.classes.iter().any(|class| class == "INFO");
            let is_format = field.classes.iter().any(|class| class == "FORMAT");
            let is_filter = field.classes.iter().any(|class| class == "FILTER");
            let both = is_info && is_format;
            if both {
                next_index += 1;
            }
            if is_filter {
                filter_index.insert(field.name.clone(), index);
            }
            let tuple = field.types.len();
            if tuple > 1 {
                next_index += (tuple as i32) * if both { 2 } else { 1 };
            }
            if field.name == "END" {
                continue;
            }
            let two_dimensional = field.lengths.len() > 1;
            if two_dimensional {
                let declared = if is_info { info.get(&field.name) } else { None };
                if !matches!(declared, Some((LineType::String, _))) {
                    return Err(format!(
                        "the multi-D field {} is not declared as an INFO String, which the native \
                         conversion reads differently.",
                        field.name
                    ));
                }
                if is_format {
                    return Err(format!(
                        "the multi-D field {} is also a FORMAT field, which this port has not \
                         measured.",
                        field.name
                    ));
                }
            }
            let element_types: Vec<CellType> = field
                .types
                .iter()
                .map(|name| genomicsdb_cell_type(name))
                .collect::<Result<_, _>>()?;
            let (cell_type, values) = if two_dimensional {
                (CellType::Char, None)
            } else {
                (element_types[0], fixed_length(&field.lengths))
            };
            if values == Some(0) {
                return Err(format!(
                    "field {} declares Number=0 without being a Flag, which is a zero-length \
                     attribute this port has not measured.",
                    field.name
                ));
            }
            if is_info {
                if tuple > 1 {
                    for element in 0..tuple {
                        info_attributes.push(Attribute {
                            name: format!("{}_tuple_element_{element}", field.name),
                            cell_type,
                            values,
                            source: Source::Info {
                                vcf_name: field.name.clone(),
                                multi_d: Some((element_types.clone(), element)),
                            },
                        });
                    }
                } else {
                    info_attributes.push(Attribute {
                        name: field.name.clone(),
                        cell_type,
                        values,
                        source: Source::Info {
                            vcf_name: field.name.clone(),
                            multi_d: two_dimensional.then(|| (element_types.clone(), 0)),
                        },
                    });
                }
            }
            if is_format {
                format_attributes.push(Attribute {
                    name: if is_info {
                        format!("{}_FORMAT", field.name)
                    } else {
                        field.name.clone()
                    },
                    cell_type,
                    values,
                    source: Source::Format {
                        vcf_name: field.name.clone(),
                    },
                });
            }
        }

        let mandatory = |name: &str, cell_type, values, source| Attribute {
            name: name.to_string(),
            cell_type,
            values,
            source,
        };
        let mut attributes = vec![
            mandatory("END", CellType::Int64, Some(1), Source::End),
            mandatory("REF", CellType::Char, None, Source::Ref),
            mandatory("ALT", CellType::Char, None, Source::Alt),
            mandatory("ID", CellType::Char, None, Source::Id),
            mandatory("QUAL", CellType::Float32, Some(1), Source::Qual),
            mandatory("FILTER", CellType::Int32, None, Source::Filter),
        ];
        attributes.extend(info_attributes);
        attributes.extend(format_attributes);
        Ok(Schema {
            attributes,
            filter_index,
            info,
            format,
            filters,
        })
    }

    /// `ArraySchema::serialize`. `array_path` is the array's absolute directory.
    pub fn serialize(&self, array_path: &str) -> Vec<u8> {
        let mut out = Vec::new();
        // The version tag, the empty array_workspace_, then the name.
        out.extend_from_slice(&0xffff_fffdu32.to_le_bytes());
        out.extend_from_slice(&0i32.to_le_bytes());
        out.extend_from_slice(&(array_path.len() as i32).to_le_bytes());
        out.extend_from_slice(array_path.as_bytes());
        // Sparse, row-major tiles, column-major cells.
        out.extend_from_slice(&[0, 0, 1]);
        out.extend_from_slice(&(CELLS_PER_TILE as i64).to_le_bytes());
        out.extend_from_slice(&(self.attributes.len() as i32).to_le_bytes());
        for attribute in &self.attributes {
            out.extend_from_slice(&(attribute.name.len() as i32).to_le_bytes());
            out.extend_from_slice(attribute.name.as_bytes());
        }
        out.extend_from_slice(&2i32.to_le_bytes());
        for dimension in ["samples", "position"] {
            out.extend_from_slice(&(dimension.len() as i32).to_le_bytes());
            out.extend_from_slice(dimension.as_bytes());
        }
        out.extend_from_slice(&32i32.to_le_bytes());
        for _ in 0..2 {
            out.extend_from_slice(&0i64.to_le_bytes());
            out.extend_from_slice(&(i64::MAX - 1).to_le_bytes());
        }
        // No tile extents: the array is sparse.
        out.extend_from_slice(&0i32.to_le_bytes());
        for attribute in &self.attributes {
            out.push(attribute.cell_type.code());
        }
        out.push(CellType::Int64.code());
        for attribute in &self.attributes {
            out.extend_from_slice(&attribute.values.unwrap_or(VAR_NUM).to_le_bytes());
        }
        // Gzip for every attribute, gzip over a delta filter for the coordinates; the offsets of
        // a variable-length attribute likewise.
        out.extend(std::iter::repeat_n(1u8, self.attributes.len()));
        out.push(0x11);
        out.extend(std::iter::repeat_n(0xffu8, self.attributes.len() + 1));
        for attribute in &self.attributes {
            out.push(if attribute.is_var() { 0x11 } else { 0 });
        }
        for attribute in &self.attributes {
            out.push(if attribute.is_var() { 0xff } else { 0 });
        }
        out
    }
}

// ---------------------------------------------------------------------------------------------
// VariantContext -> VCF stream -> VCF2Binary
// ---------------------------------------------------------------------------------------------

/// What the VCF stream holds for one field, as `bcf_get_*_values` hands it back.
#[derive(Debug, Clone, PartialEq)]
enum Encoded {
    /// `encodeType(0, ...)`, or a key the record does not carry.
    Absent,
    Ints(Vec<i32>),
    /// The bits, since a missing float is a signalling NaN that must survive.
    Floats(Vec<u32>),
    Text(String),
    Flag,
}

/// `VariantContext.decodeOne`, with the file's own spelling of the value.
enum Decoded {
    Int(i32),
    Float(f64),
    Text(String),
    True,
}

/// `VCFUtils.parseVcfDouble`.
fn parse_vcf_double(text: &str) -> Option<f64> {
    java_parse_double(text)
}

fn decode_one(field: &str, text: &str, line_type: LineType) -> Result<Option<Decoded>, Reason> {
    if text == "." {
        return Ok(None);
    }
    match line_type {
        LineType::Integer => text
            .parse::<i32>()
            .map(|value| Some(Decoded::Int(value)))
            .map_err(|_| format!("field {field} holds \"{text}\" where the header says Integer.")),
        LineType::Float => parse_vcf_double(text)
            .map(|value| Some(Decoded::Float(value)))
            .ok_or_else(|| format!("field {field} holds \"{text}\" where the header says Float.")),
        LineType::Flag => {
            if text.eq_ignore_ascii_case("true") || text == "1" {
                Ok(Some(Decoded::True))
            } else {
                Err(format!("Flag field {field} holds \"{text}\"."))
            }
        }
        LineType::String | LineType::Character => Ok(Some(Decoded::Text(text.to_string()))),
    }
}

/// What `VCF2Binary` holds a numeric field's length to: `G` and `.` are not checked, and the
/// stream carries as many values as the record does, so a mismatch reaches the native check.
fn check_count(
    field: &str,
    held: usize,
    number: &Cardinality,
    vc: &VariantContext,
) -> Result<(), Reason> {
    let expected = match number {
        Cardinality::Fixed(count) => (*count).max(0) as usize,
        Cardinality::A => vc.alleles.len().saturating_sub(1),
        Cardinality::R => vc.alleles.len(),
        Cardinality::G | Cardinality::Unbounded => return Ok(()),
    };
    if held != expected {
        return Err(format!(
            "field {field} holds {held} values where its header line asks for {expected}, which \
             the native loader fails on."
        ));
    }
    Ok(())
}

/// `VCFEncoder`'s text as htslib reads it: as many values as the record holds, each in its place,
/// and a list of strings joined back with its commas, a missing one as the `.` it was.
fn encode_decoded(
    field: &str,
    values: Vec<Option<Decoded>>,
    line_type: LineType,
    number: &Cardinality,
    vc: &VariantContext,
) -> Result<Encoded, Reason> {
    match line_type {
        LineType::String | LineType::Character => {
            let strings: Vec<String> = values
                .into_iter()
                .map(|value| match value {
                    Some(Decoded::Text(text)) => text,
                    _ => ".".to_string(),
                })
                .collect();
            let text = strings.join(",");
            Ok(if text.is_empty() {
                Encoded::Absent
            } else {
                Encoded::Text(text)
            })
        }
        LineType::Flag => Ok(if values.is_empty() {
            Encoded::Absent
        } else {
            Encoded::Flag
        }),
        LineType::Integer | LineType::Float => {
            check_count(field, values.len(), number, vc)?;
            Ok(if line_type == LineType::Integer {
                let ints: Vec<i32> = values
                    .iter()
                    .map(|value| match value {
                        Some(Decoded::Int(value)) => *value,
                        _ => BCF_INT32_MISSING,
                    })
                    .collect();
                if ints.is_empty() {
                    Encoded::Absent
                } else {
                    Encoded::Ints(ints)
                }
            } else {
                let floats: Vec<u32> = values
                    .iter()
                    .map(|value| match value {
                        Some(Decoded::Float(value)) => (*value as f32).to_bits(),
                        _ => BCF_FLOAT_MISSING,
                    })
                    .collect();
                if floats.is_empty() {
                    Encoded::Absent
                } else {
                    Encoded::Floats(floats)
                }
            })
        }
    }
}

/// `VariantContext.decodeInfo` and the `GenericSiteWriter`.
fn encode_site(
    field: &str,
    raw: Option<&Value>,
    line_type: LineType,
    number: &Cardinality,
    vc: &VariantContext,
) -> Result<Encoded, Reason> {
    let Some(raw) = raw else {
        return Ok(Encoded::Absent);
    };
    let is_string = matches!(line_type, LineType::String | LineType::Character);
    let (values, was_list) = match raw {
        // A string keeps the `.`, and an empty value is read as one.
        Value::Missing if is_string => return Ok(Encoded::Text(".".to_string())),
        Value::Str(text) if is_string && text == "." => return Ok(Encoded::Text(".".to_string())),
        Value::Missing => return Ok(Encoded::Absent),
        Value::Bool(true) => (vec![Some(Decoded::True)], false),
        Value::Str(text) => (vec![decode_one(field, text, line_type)?], false),
        Value::List(items) => {
            let mut decoded = Vec::with_capacity(items.len());
            for item in items {
                match item {
                    Value::Str(text) => decoded.push(decode_one(field, text, line_type)?),
                    Value::Missing => decoded.push(None),
                    other => {
                        return Err(format!(
                            "field {field} holds the unexpected value {other:?}."
                        ))
                    }
                }
            }
            (decoded, true)
        }
        other => {
            return Err(format!(
                "field {field} holds the unexpected value {other:?}."
            ))
        }
    };
    if !was_list && values[0].is_none() {
        return Ok(Encoded::Absent);
    }
    encode_decoded(field, values, line_type, number, vc)
}

/// The `decodeValue` of a genotype attribute, whose comma split is `String.split`'s.
fn encode_extended(
    field: &str,
    text: &str,
    line_type: LineType,
    number: &Cardinality,
    vc: &VariantContext,
) -> Result<Encoded, Reason> {
    if text.contains(',') {
        let mut parts: Vec<&str> = text.split(',').collect();
        while parts.last().is_some_and(|part| part.is_empty()) {
            parts.pop();
        }
        let values = parts
            .iter()
            .map(|part| decode_one(field, part, line_type))
            .collect::<Result<Vec<_>, _>>()?;
        encode_decoded(field, values, line_type, number, vc)
    } else {
        let value = decode_one(field, text, line_type)?;
        if value.is_none() {
            return Ok(Encoded::Absent);
        }
        encode_decoded(field, vec![value], line_type, number, vc)
    }
}

/// `IntGenotypeFieldAccessors`: the array as it is.
fn encode_int_accessor(
    field: &str,
    values: &[i32],
    number: &Cardinality,
    vc: &VariantContext,
) -> Result<Encoded, Reason> {
    check_count(field, values.len(), number, vc)?;
    Ok(if values.is_empty() {
        Encoded::Absent
    } else {
        Encoded::Ints(values.to_vec())
    })
}

/// `strtoll(text, &endptr, 0)` narrowed to an int, for what a well-formed annotation can hold.
fn strtol_base_zero(text: &str) -> i32 {
    let bytes = text.as_bytes();
    let mut position = 0;
    while position < bytes.len() && bytes[position].is_ascii_whitespace() {
        position += 1;
    }
    let negative = position < bytes.len() && bytes[position] == b'-';
    if position < bytes.len() && (bytes[position] == b'-' || bytes[position] == b'+') {
        position += 1;
    }
    let mut radix = 10i64;
    if position + 1 < bytes.len()
        && bytes[position] == b'0'
        && (bytes[position + 1] == b'x' || bytes[position + 1] == b'X')
        && position + 2 < bytes.len()
        && bytes[position + 2].is_ascii_hexdigit()
    {
        radix = 16;
        position += 2;
    } else if position < bytes.len() && bytes[position] == b'0' {
        radix = 8;
    }
    let mut value = 0i64;
    while position < bytes.len() {
        let digit = match (bytes[position] as char).to_digit(radix as u32) {
            Some(digit) => digit as i64,
            None => break,
        };
        value = value.saturating_mul(radix).saturating_add(digit);
        position += 1;
    }
    (if negative { -value } else { value }) as i32
}

/// `strtof`, for the decimal forms an annotation holds. `None` for what is not plainly one, which
/// glibc reads by rules this port has not carried.
fn strtof_plain(text: &str) -> Option<f32> {
    let is_plain = !text.is_empty()
        && text
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'.' | b'-' | b'+' | b'e' | b'E'));
    if is_plain {
        text.parse::<f32>().ok()
    } else {
        None
    }
}

/// `GenomicsDBMultiDVectorField::parse_and_store_numeric` over `|` then `,`, in the default
/// operator since a stream here holds one sample. One buffer per tuple element.
fn parse_multi_d(field: &str, text: &str, types: &[CellType]) -> Result<Vec<Vec<u8>>, Reason> {
    let tuple = types.len();
    let mut data: Vec<Vec<u8>> = vec![Vec::new(); tuple];
    let mut offsets: Vec<Vec<u64>> = vec![vec![0]; tuple];
    let mut read_in_row = vec![0usize; tuple];
    let mut most_in_row = 0usize;
    let mut element_in_tuple = 0usize;

    let bytes = text.as_bytes();
    let mut begin = 0usize;
    let mut position = 0usize;
    while position <= bytes.len() {
        let separator = if position == bytes.len() {
            Some(b'\0')
        } else if bytes[position] == b'|' || bytes[position] == b',' {
            Some(bytes[position])
        } else {
            None
        };
        if let Some(separator) = separator {
            let element = &text[begin..position];
            let is_missing = element.is_empty() || element.eq_ignore_ascii_case("nan");
            match types[element_in_tuple] {
                CellType::Int32 => {
                    let value = if is_missing {
                        BCF_INT32_MISSING
                    } else {
                        strtol_base_zero(element)
                    };
                    data[element_in_tuple].extend_from_slice(&value.to_le_bytes());
                }
                CellType::Float32 => {
                    let bits = if is_missing {
                        BCF_FLOAT_MISSING
                    } else {
                        strtof_plain(element)
                            .ok_or_else(|| {
                                format!(
                                    "the multi-D field {field} holds \"{element}\", which strtof \
                                     reads by rules this port does not carry."
                                )
                            })?
                            .to_bits()
                    };
                    data[element_in_tuple].extend_from_slice(&bits.to_le_bytes());
                }
                _ => {
                    return Err(format!(
                        "multi-D field {field} has an unexpected element type."
                    ))
                }
            }
            read_in_row[element_in_tuple] += 1;
            most_in_row = most_in_row.max(read_in_row[element_in_tuple]);
            if separator == b',' {
                element_in_tuple = (element_in_tuple + 1) % tuple;
            } else {
                for (index, cell_type) in types.iter().enumerate() {
                    for _ in read_in_row[index]..most_in_row {
                        match cell_type {
                            CellType::Int32 => {
                                data[index].extend_from_slice(&BCF_INT32_MISSING.to_le_bytes())
                            }
                            _ => data[index].extend_from_slice(&BCF_FLOAT_MISSING.to_le_bytes()),
                        }
                    }
                    let last = *offsets[index].last().expect("starts with the zero offset");
                    offsets[index].push(last + (most_in_row * 4) as u64);
                }
                element_in_tuple = 0;
                read_in_row = vec![0; tuple];
                most_in_row = 0;
            }
            begin = position + 1;
        }
        position += 1;
    }

    Ok((0..tuple)
        .map(|index| {
            let mut out = Vec::with_capacity(data[index].len() + 16 + offsets[index].len() * 8);
            out.extend_from_slice(&(data[index].len() as u64).to_le_bytes());
            out.extend_from_slice(&data[index]);
            out.extend_from_slice(&((offsets[index].len() - 1) as u64).to_le_bytes());
            for offset in &offsets[index] {
                out.extend_from_slice(&offset.to_le_bytes());
            }
            out
        })
        .collect())
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum VariantType {
    Indel,
    SpanningDeletion,
    NonRef,
    Other,
}

/// htslib's `bcf_set_variant_type`, in the fork GenomicsDB builds, for the three types read here.
fn variant_type(reference: &[u8], alternate: &[u8]) -> VariantType {
    if alternate == b"*" {
        return VariantType::SpanningDeletion;
    }
    if reference.len() == 1 && alternate.len() == 1 {
        return VariantType::Other;
    }
    if alternate.first() == Some(&b'<') {
        if alternate.starts_with(b"<NON_REF>") {
            return VariantType::NonRef;
        }
        return VariantType::Other;
    }
    let mut common = 0;
    while common < reference.len()
        && common < alternate.len()
        && reference[common].eq_ignore_ascii_case(&alternate[common])
    {
        common += 1;
    }
    let reference_left = reference.len() - common;
    let alternate_left = alternate.len() - common;
    if alternate_left > 0 && reference_left == 0 {
        if alternate[common] == b']' || alternate[common] == b'[' {
            return VariantType::Other;
        }
        return VariantType::Indel;
    }
    if reference_left > 0 && alternate_left == 0 {
        return VariantType::Indel;
    }
    if reference_left == 0 && alternate_left == 0 {
        return VariantType::Other;
    }
    let mut reference_end = reference.len() - 1;
    let mut alternate_end = alternate.len() - 1;
    while reference_end > common
        && alternate_end > common
        && reference[reference_end].eq_ignore_ascii_case(&alternate[alternate_end])
    {
        reference_end -= 1;
        alternate_end -= 1;
    }
    let bases_match = reference[reference_end].eq_ignore_ascii_case(&alternate[alternate_end]);
    if alternate_end == common {
        if reference_end == common {
            return VariantType::Other;
        }
        return if bases_match {
            VariantType::Indel
        } else {
            VariantType::Other
        };
    }
    if reference_end == common {
        return if bases_match {
            VariantType::Indel
        } else {
            VariantType::Other
        };
    }
    VariantType::Other
}

/// One cell in `VCF2Binary`'s layout. `attributes[0]`, the `END`, is written from `end` when the
/// cell goes to disk, since the copy at the end of an interval swaps the two.
#[derive(Debug, Clone)]
pub struct Cell {
    row: i64,
    column: i64,
    end: i64,
    attributes: Rc<Vec<Vec<u8>>>,
}

impl Schema {
    /// `VCF2Binary::convert_VCF_to_binary_for_callset` for the single sample of a stream.
    pub fn cell(&self, vc: &VariantContext, row: i64, contig_offset: i64) -> Result<Cell, Reason> {
        // What `VCFEncoder` refuses before anything is written.
        for (key, _) in &vc.attributes {
            if !self.info.contains_key(key) {
                return Err(format!(
                    "the record at {}:{} carries INFO field {key}, which the merged header does \
                     not declare and which the reference's VCF encoder fails on.",
                    vc.contig, vc.start
                ));
            }
        }
        let genotype = vc.genotypes.iter().next();
        if let Some(genotype) = genotype {
            self.check_format_keys(vc, genotype)?;
        }

        let column = contig_offset + vc.start - 1;
        let reference = vc
            .alleles
            .first()
            .map(|allele| allele.display_string())
            .unwrap_or_default();
        let alternates: Vec<String> = vc.alleles[1.min(vc.alleles.len())..]
            .iter()
            .map(|allele| allele.display_string())
            .collect();

        let end_attribute = vc
            .attributes
            .iter()
            .rev()
            .find(|(key, _)| key == "END")
            .map(|(_, value)| value);
        let end = match end_attribute {
            Some(Value::Str(text)) => {
                let end: i64 = text.parse().map_err(|_| {
                    format!(
                        "the record at {}:{} has END={text}, which the reference fails on.",
                        vc.contig, vc.start
                    )
                })?;
                contig_offset + end - 1
            }
            Some(Value::Missing) | None => {
                // `m_treat_deletions_as_intervals`.
                let mut end = column;
                for alternate in &alternates {
                    let kind = variant_type(reference.as_bytes(), alternate.as_bytes());
                    if matches!(kind, VariantType::Indel | VariantType::SpanningDeletion)
                        && reference.len() > alternate.len()
                    {
                        end = column + reference.len() as i64 - 1;
                        break;
                    }
                }
                end
            }
            Some(other) => {
                return Err(format!(
                    "the record at {}:{} has END={other:?}, which the reference fails on.",
                    vc.contig, vc.start
                ))
            }
        };
        if end < column {
            return Err(format!(
                "the record at {}:{} has an END before its start, which breaks an assertion of \
                 the native loader.",
                vc.contig, vc.start
            ));
        }

        let mut attributes = Vec::with_capacity(self.attributes.len());
        for attribute in &self.attributes {
            attributes.push(self.attribute_bytes(
                attribute,
                vc,
                genotype,
                &reference,
                &alternates,
            )?);
        }
        Ok(Cell {
            row,
            column,
            end,
            attributes: Rc::new(attributes),
        })
    }

    fn check_format_keys(&self, vc: &VariantContext, genotype: &Genotype) -> Result<(), Reason> {
        let undeclared = |key: &str| {
            format!(
                "the record at {}:{} carries FORMAT field {key}, which the merged header does \
                 not declare and which the reference's VCF encoder fails on.",
                vc.contig, vc.start
            )
        };
        let mut keys: Vec<&str> = genotype
            .extended
            .iter()
            .map(|(key, _)| key.as_str())
            .collect();
        if !genotype.alleles.is_empty() {
            keys.push("GT");
        }
        if genotype.gq.is_some() {
            keys.push("GQ");
        }
        if genotype.dp.is_some() {
            keys.push("DP");
        }
        if genotype.ad.is_some() {
            keys.push("AD");
        }
        if genotype.pl.is_some() {
            keys.push("PL");
        }
        if genotype.filters.is_some() {
            keys.push("FT");
        }
        for key in keys {
            if !self.format.contains_key(key) {
                return Err(undeclared(key));
            }
        }
        for (key, line_type) in [
            ("GQ", "Integer"),
            ("DP", "Integer"),
            ("AD", "Integer"),
            ("PL", "Integer"),
        ] {
            if let Some((declared, _)) = self.format.get(key) {
                if line_type_name(*declared) != line_type {
                    return Err(format!(
                        "FORMAT field {key} is declared {} where the reference's VCF encoder writes \
                         integers.",
                        line_type_name(*declared)
                    ));
                }
            }
        }
        if let Some((declared, _)) = self.format.get("GT") {
            if *declared != LineType::String {
                return Err(
                    "FORMAT field GT is not declared a String, which htslib refuses to read back."
                        .to_string(),
                );
            }
        }
        Ok(())
    }

    fn attribute_bytes(
        &self,
        attribute: &Attribute,
        vc: &VariantContext,
        genotype: Option<&Genotype>,
        reference: &str,
        alternates: &[String],
    ) -> Result<Vec<u8>, Reason> {
        match &attribute.source {
            // Filled in as the cell is written.
            Source::End => Ok(Vec::new()),
            Source::Ref => Ok(reference.as_bytes().to_vec()),
            Source::Alt => {
                let mut out = Vec::new();
                for (index, alternate) in alternates.iter().enumerate() {
                    if index > 0 {
                        out.push(b'|');
                    }
                    if variant_type(reference.as_bytes(), alternate.as_bytes())
                        == VariantType::NonRef
                    {
                        out.push(b'&');
                    } else {
                        out.extend_from_slice(alternate.as_bytes());
                    }
                }
                Ok(out)
            }
            Source::Id => Ok(if vc.id == "." {
                Vec::new()
            } else {
                vc.id.as_bytes().to_vec()
            }),
            Source::Qual => Ok(if vc.has_log10_p_error() {
                // `VCFEncoder.formatQualValue`, two decimals, then what htslib reads back.
                let text = htsjdk_vcf::encoder::format_qual_value(vc.phred_scaled_qual());
                let qual = text
                    .parse::<f64>()
                    .map_err(|_| format!("the QUAL {text} is not one htslib reads as a number."))?;
                (qual as f32).to_bits().to_le_bytes().to_vec()
            } else {
                NULL_FLOAT32.to_le_bytes().to_vec()
            }),
            Source::Filter => self.filter_bytes(vc),
            Source::Info { vcf_name, multi_d } => {
                let (line_type, number) = self
                    .info
                    .get(vcf_name)
                    .ok_or_else(|| format!("INFO field {vcf_name} has no header line."))?;
                let raw = vc
                    .attributes
                    .iter()
                    .rev()
                    .find(|(key, _)| key == vcf_name)
                    .map(|(_, value)| value);
                let encoded = encode_site(vcf_name, raw, *line_type, number, vc)?;
                match multi_d {
                    Some((types, element)) => match encoded {
                        Encoded::Text(text) => {
                            Ok(parse_multi_d(vcf_name, &text, types)?.swap_remove(*element))
                        }
                        _ => Ok(Vec::new()),
                    },
                    None => cell_bytes(vcf_name, encoded, *line_type, number, attribute),
                }
            }
            Source::Format { vcf_name } => {
                let (line_type, number) = self
                    .format
                    .get(vcf_name)
                    .ok_or_else(|| format!("FORMAT field {vcf_name} has no header line."))?;
                let Some(genotype) = genotype else {
                    return cell_bytes(vcf_name, Encoded::Absent, *line_type, number, attribute);
                };
                if vcf_name == "GT" {
                    return gt_bytes(vc, genotype);
                }
                let encoded = match vcf_name.as_str() {
                    "GQ" => match genotype.gq {
                        Some(value) => encode_int_accessor(vcf_name, &[value], number, vc)?,
                        None => Encoded::Absent,
                    },
                    "DP" => match genotype.dp {
                        Some(value) => encode_int_accessor(vcf_name, &[value], number, vc)?,
                        None => Encoded::Absent,
                    },
                    "AD" => match &genotype.ad {
                        Some(values) => encode_int_accessor(vcf_name, values, number, vc)?,
                        None => Encoded::Absent,
                    },
                    "PL" => match &genotype.pl {
                        Some(values) => encode_int_accessor(vcf_name, values, number, vc)?,
                        None => Encoded::Absent,
                    },
                    "FT" => match &genotype.filters {
                        // The FT as written, which the runner put back: the stream is raw text.
                        Some(text) => Encoded::Text(text.clone()),
                        None => Encoded::Absent,
                    },
                    _ => match genotype
                        .extended
                        .iter()
                        .rev()
                        .find(|(key, _)| key == vcf_name)
                        .map(|(_, value)| value)
                    {
                        Some(Value::Str(text)) => {
                            encode_extended(vcf_name, text, *line_type, number, vc)?
                        }
                        Some(Value::Missing) | None => Encoded::Absent,
                        Some(other) => {
                            return Err(format!(
                                "FORMAT field {vcf_name} holds the unexpected value {other:?}."
                            ))
                        }
                    },
                };
                cell_bytes(vcf_name, encoded, *line_type, number, attribute)
            }
        }
    }

    /// The filters of the text, and the local-to-global index the cell stores.
    fn filter_bytes(&self, vc: &VariantContext) -> Result<Vec<u8>, Reason> {
        let names: Vec<String> = match &vc.filters {
            None => Vec::new(),
            Some(filters) if filters.is_empty() => vec!["PASS".to_string()],
            Some(filters) => gatk_engine::java_hash::hash_set_order(filters).map_err(|_| {
                format!(
                    "the record at {}:{} has more filters than this port lays out as a Java set.",
                    vc.contig, vc.start
                )
            })?,
        };
        let mut out = Vec::with_capacity(names.len() * 4);
        for name in names {
            if !self.filters.contains(&name) && name != "PASS" {
                return Err(format!(
                    "the record at {}:{} carries filter {name}, which the merged header does not \
                     declare and which the reference's VCF encoder fails on.",
                    vc.contig, vc.start
                ));
            }
            // `PASS` with no header line of its own has no global index.
            let index = self.filter_index.get(&name).copied().unwrap_or(-1);
            out.extend_from_slice(&index.to_le_bytes());
        }
        Ok(out)
    }
}

/// The GT of the text and `is_GT_field` in `convert_field_to_tiledb`: allele, then phase
/// and allele for each further one.
fn gt_bytes(vc: &VariantContext, genotype: &Genotype) -> Result<Vec<u8>, Reason> {
    if vc.alleles.len() > 127 {
        return Err(format!(
            "the record at {}:{} has more alleles than the reference's stream carries \
             genotypes for.",
            vc.contig, vc.start
        ));
    }
    let mut out = Vec::with_capacity(genotype.alleles.len() * 8);
    for (index, allele) in genotype.alleles.iter().enumerate() {
        if index > 0 {
            let phase: i32 = if genotype.phased { 1 } else { 0 };
            out.extend_from_slice(&phase.to_le_bytes());
        }
        let value: i32 = if allele.is_no_call() {
            -1
        } else {
            vc.alleles
                .iter()
                .position(|candidate| candidate == allele)
                .ok_or_else(|| {
                    format!(
                        "the record at {}:{} has a genotype allele it does not list.",
                        vc.contig, vc.start
                    )
                })? as i32
        };
        out.extend_from_slice(&value.to_le_bytes());
    }
    Ok(out)
}

/// The rest of `convert_field_to_tiledb`: a missing fixed field is that many nulls, a missing
/// variable one is empty, a lone missing value is missing, and strings stop at a NUL.
fn cell_bytes(
    field: &str,
    encoded: Encoded,
    line_type: LineType,
    number: &Cardinality,
    attribute: &Attribute,
) -> Result<Vec<u8>, Reason> {
    let is_string = matches!(line_type, LineType::String | LineType::Character);
    let header_fixed = match number {
        Cardinality::Fixed(count) if !is_string => Some(if line_type == LineType::Flag {
            1
        } else {
            *count
        }),
        _ => None,
    };
    let is_missing = match &encoded {
        Encoded::Absent => true,
        Encoded::Ints(values) => values.len() == 1 && values[0] == BCF_INT32_MISSING,
        Encoded::Floats(values) => values.len() == 1 && values[0] == BCF_FLOAT_MISSING,
        Encoded::Text(_) | Encoded::Flag => false,
    };
    let out = if is_missing {
        match header_fixed {
            Some(count) => {
                let mut out = Vec::new();
                for _ in 0..count {
                    match line_type {
                        LineType::Integer => out.extend_from_slice(&NULL_INT32.to_le_bytes()),
                        LineType::Float => out.extend_from_slice(&NULL_FLOAT32.to_le_bytes()),
                        _ => out.push(NULL_CHAR),
                    }
                }
                out
            }
            None => Vec::new(),
        }
    } else {
        match encoded {
            Encoded::Ints(values) => values
                .iter()
                .flat_map(|value| value.to_le_bytes())
                .collect(),
            Encoded::Floats(values) => values
                .iter()
                .flat_map(|value| value.to_le_bytes())
                .collect(),
            Encoded::Text(text) => {
                let bytes = text.as_bytes();
                let length = bytes
                    .iter()
                    .position(|byte| *byte == 0)
                    .unwrap_or(bytes.len());
                bytes[..length].to_vec()
            }
            Encoded::Flag => vec![1],
            Encoded::Absent => Vec::new(),
        }
    };
    if let Some(values) = attribute.values {
        let size = match attribute.cell_type {
            CellType::Int32 | CellType::Float32 => 4,
            CellType::Int64 => 8,
            CellType::Char => 1,
        };
        if out.len() != values as usize * size {
            return Err(format!(
                "field {field} does not fill its fixed-length attribute, which the native \
                 conversion refuses with a length mismatch."
            ));
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// LoaderArrayWriter
// ---------------------------------------------------------------------------------------------

/// `VCF2TileDBLoader::produce_cells_in_column_major_order` over the rows of a batch, then
/// `LoaderArrayWriter::operate` and `finish`: the cells in the order TileDB is handed them.
///
/// A row's cells come in stream order. `rows` is indexed by nothing: each cell knows its row.
pub fn load(
    rows: Vec<Vec<Cell>>,
    row_count: usize,
    partition_begin: i64,
) -> Result<Vec<Cell>, Reason> {
    let mut loader = Loader {
        partition_begin,
        crossed: false,
        last_end: vec![-1; row_count],
        copies: vec![None; row_count],
        pending: BTreeMap::new(),
        written: Vec::new(),
    };
    let mut next: Vec<usize> = vec![0; rows.len()];
    let mut previous: Option<(i64, i64)> = None;
    loop {
        let mut best: Option<(i64, i64, usize)> = None;
        for (index, row) in rows.iter().enumerate() {
            if let Some(cell) = row.get(next[index]) {
                let key = (cell.column, cell.row, index);
                if best.is_none_or(|candidate| key < candidate) {
                    best = Some(key);
                }
            }
        }
        let Some((column, row, index)) = best else {
            break;
        };
        if let Some(last) = previous {
            if (column, row) <= last {
                return Err(format!(
                    "two records of a sample share a position or are out of order at column \
                     {column}, which the native loader fails on with \"Incorrect cell order \
                     found\"."
                ));
            }
        }
        previous = Some((column, row));
        let cell = rows[index][next[index]].clone();
        next[index] += 1;
        loader.operate(cell)?;
    }
    loader.finish()?;
    Ok(loader.written)
}

struct Pending {
    /// The swapped value once this is the copy at an interval's end.
    other_end: i64,
    cell: Cell,
}

struct Loader {
    partition_begin: i64,
    crossed: bool,
    last_end: Vec<i64>,
    copies: Vec<Option<Cell>>,
    /// `m_cell_wrapper_pq`, by column then row.
    pending: BTreeMap<(i64, i64), Pending>,
    written: Vec<Cell>,
}

impl Loader {
    /// `handle_intervals_spanning_partition_begin`, whose copies stand until the first cell past
    /// the partition's begin.
    fn handle_spanning(&mut self, cell: Option<&Cell>) -> Result<(), Reason> {
        let (begin, end) = match cell {
            Some(cell) => (cell.column, cell.end),
            None => (i64::MAX - 1, i64::MAX - 1),
        };
        if begin > self.partition_begin {
            if !self.crossed {
                self.crossed = true;
                let mut copies = Vec::new();
                for row in 0..self.last_end.len() {
                    if self.last_end[row] >= 0 {
                        if let Some(copy) = self.copies[row].take() {
                            copies.push(copy);
                        }
                    }
                    self.last_end[row] = -1;
                }
                copies.sort_by_key(|copy| (copy.column, copy.row));
                for copy in copies {
                    self.operate(copy)?;
                }
            }
            return Ok(());
        }
        let cell = cell.expect("the sentinel is always past the partition's begin");
        let row = cell.row as usize;
        if end >= self.partition_begin {
            self.last_end[row] = end;
            self.copies[row] = Some(cell.clone());
        } else {
            self.last_end[row] = -1;
        }
        Ok(())
    }

    fn operate(&mut self, cell: Cell) -> Result<(), Reason> {
        if !self.crossed {
            self.handle_spanning(Some(&cell))?;
            if !self.crossed {
                return Ok(());
            }
        }
        let row = cell.row;
        let begin = cell.column;
        // Everything less than (row+1, begin-1) goes to disk.
        while let Some((&key, _)) = self.pending.first_key_value() {
            if key < (begin - 1, row + 1) {
                self.write_top();
            } else {
                break;
            }
        }
        if self.last_end[row as usize] >= begin {
            let key = self
                .pending
                .keys()
                .find(|(_, pending_row)| *pending_row == row)
                .copied()
                .ok_or_else(|| format!("no pending cell for row {row} at column {begin}."))?;
            let pending = self.pending.remove(&key).expect("the key was just found");
            // Always the copy at an interval's end: its begin cell already went to disk.
            if pending.other_end >= key.0 {
                return Err(format!(
                    "two cells in incorrect order at row {row}, column {begin}, which the native \
                     loader fails on."
                ));
            }
            let truncated = begin - 1;
            if truncated != pending.other_end {
                let mut copy = pending.cell;
                copy.column = truncated;
                copy.end = pending.other_end;
                self.written.push(copy);
            }
        }
        let end = cell.end;
        if self
            .pending
            .insert(
                (begin, row),
                Pending {
                    other_end: end,
                    cell,
                },
            )
            .is_some()
        {
            return Err(format!(
                "two cells at row {row}, column {begin}, which the native loader orders by a \
                 heap this port does not reproduce."
            ));
        }
        self.last_end[row as usize] = end;
        Ok(())
    }

    /// `write_top_element_to_disk`.
    fn write_top(&mut self) {
        let ((column, row), pending) = self
            .pending
            .pop_first()
            .expect("only called with something pending");
        let mut cell = pending.cell;
        cell.column = column;
        cell.end = pending.other_end;
        self.written.push(cell.clone());
        if pending.other_end > column {
            self.pending.insert(
                (pending.other_end, row),
                Pending {
                    other_end: column,
                    cell,
                },
            );
        }
    }

    fn finish(&mut self) -> Result<(), Reason> {
        if !self.crossed {
            self.handle_spanning(None)?;
        }
        while !self.pending.is_empty() {
            self.write_top();
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------------------------
// TileDB's WriteState and BookKeeping
// ---------------------------------------------------------------------------------------------

/// TileDB's `CodecGzip`: `deflateInit` at the default level, one `deflate(Z_FINISH)`.
///
/// This is the C zlib on purpose: the workspace's flate2 is built on a different deflate whose
/// output differs from the library the reference links.
fn zlib_deflate(data: &[u8], window_bits: i32) -> Vec<u8> {
    use libz_sys::*;
    use std::mem::MaybeUninit;

    let mut stream = MaybeUninit::<z_stream>::zeroed();
    // SAFETY: a zeroed z_stream with null allocators is what `deflateInit2_` documents as the
    // start state, every pointer passed outlives the calls, and the bound leaves room for the
    // whole of the output so `deflate` finishes in one call.
    unsafe {
        let status = deflateInit2_(
            stream.as_mut_ptr(),
            Z_DEFAULT_COMPRESSION,
            Z_DEFLATED,
            window_bits,
            8,
            Z_DEFAULT_STRATEGY,
            zlibVersion(),
            std::mem::size_of::<z_stream>() as i32,
        );
        assert_eq!(status, Z_OK, "deflateInit2_ takes these constant arguments");
        let stream = stream.assume_init_mut();
        let bound = deflateBound(stream, data.len() as uLong) as usize + 64;
        let mut out = vec![0u8; bound];
        stream.next_in = data.as_ptr() as *mut Bytef;
        stream.avail_in = data.len() as uInt;
        stream.next_out = out.as_mut_ptr();
        stream.avail_out = bound as uInt;
        let status = deflate(stream, Z_FINISH);
        assert_eq!(status, Z_STREAM_END, "the output buffer holds the bound");
        out.truncate(stream.total_out as usize);
        deflateEnd(stream);
        out
    }
}

/// The delta filter TileDB applies before gzip to a tile of offsets (`stride` 1) or of coordinates
/// (`stride` 2), each stream starting from zero.
fn delta_encode(values: &[i64], stride: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 8);
    let mut previous = vec![0i64; stride];
    for (index, value) in values.iter().enumerate() {
        let slot = index % stride;
        out.extend_from_slice(&value.wrapping_sub(previous[slot]).to_le_bytes());
        previous[slot] = *value;
    }
    out
}

/// `WriteState` and `BookKeeping` for one fragment: each file in it, by name. An empty list when
/// no cell reached it, since TileDB then leaves no fragment behind.
pub fn fragment_files(schema: &Schema, cells: &[Cell]) -> Vec<(String, Vec<u8>)> {
    if cells.is_empty() {
        return Vec::new();
    }
    let mut files = Vec::new();
    let tiles: Vec<&[Cell]> = cells.chunks(CELLS_PER_TILE).collect();

    let mut book_keeping = Vec::new();
    let put_i64 = |out: &mut Vec<u8>, value: i64| out.extend_from_slice(&value.to_le_bytes());
    // The domain.
    book_keeping.extend_from_slice(&32u64.to_le_bytes());
    for _ in 0..2 {
        put_i64(&mut book_keeping, 0);
        put_i64(&mut book_keeping, i64::MAX - 1);
    }
    // The MBRs, then the bounding coordinates.
    put_i64(&mut book_keeping, tiles.len() as i64);
    for tile in &tiles {
        let rows = tile.iter().map(|cell| cell.row);
        let columns = tile.iter().map(|cell| cell.column);
        put_i64(&mut book_keeping, rows.clone().min().unwrap_or(0));
        put_i64(&mut book_keeping, rows.max().unwrap_or(0));
        put_i64(&mut book_keeping, columns.clone().min().unwrap_or(0));
        put_i64(&mut book_keeping, columns.max().unwrap_or(0));
    }
    put_i64(&mut book_keeping, tiles.len() as i64);
    for tile in &tiles {
        let first = &tile[0];
        let last = &tile[tile.len() - 1];
        put_i64(&mut book_keeping, first.row);
        put_i64(&mut book_keeping, first.column);
        put_i64(&mut book_keeping, last.row);
        put_i64(&mut book_keeping, last.column);
    }

    let mut tile_offsets_section = Vec::new();
    let mut var_offsets_section = Vec::new();
    let mut var_sizes_section = Vec::new();

    for (index, attribute) in schema.attributes.iter().enumerate() {
        let mut file = Vec::new();
        let mut var_file = Vec::new();
        let mut tile_offsets = Vec::with_capacity(tiles.len());
        let mut var_offsets = Vec::with_capacity(tiles.len());
        let mut var_sizes = Vec::with_capacity(tiles.len());
        let mut position = 0i64;
        for tile in &tiles {
            tile_offsets.push(file.len() as i64);
            if attribute.is_var() {
                let mut offsets = Vec::with_capacity(tile.len());
                let mut data = Vec::new();
                for cell in tile.iter() {
                    offsets.push(position + data.len() as i64);
                    data.extend_from_slice(&cell.attributes[index]);
                }
                position += data.len() as i64;
                file.extend(zlib_deflate(&delta_encode(&offsets, 1), 15));
                var_offsets.push(var_file.len() as i64);
                var_sizes.push(data.len() as i64);
                if !data.is_empty() {
                    var_file.extend(zlib_deflate(&data, 15));
                }
            } else {
                let mut data = Vec::new();
                for cell in tile.iter() {
                    if matches!(attribute.source, Source::End) {
                        data.extend_from_slice(&cell.end.to_le_bytes());
                    } else {
                        data.extend_from_slice(&cell.attributes[index]);
                    }
                }
                file.extend(zlib_deflate(&data, 15));
            }
        }
        put_i64(&mut tile_offsets_section, tile_offsets.len() as i64);
        for offset in &tile_offsets {
            put_i64(&mut tile_offsets_section, *offset);
        }
        files.push((format!("{}.tdb", attribute.name), file));
        // A fixed-length attribute has an empty list in both sections, not no list.
        put_i64(&mut var_offsets_section, var_offsets.len() as i64);
        for offset in &var_offsets {
            put_i64(&mut var_offsets_section, *offset);
        }
        put_i64(&mut var_sizes_section, var_sizes.len() as i64);
        for size in &var_sizes {
            put_i64(&mut var_sizes_section, *size);
        }
        if attribute.is_var() {
            files.push((format!("{}_var.tdb", attribute.name), var_file));
        }
    }

    let mut coordinates = Vec::new();
    let mut coordinate_offsets = Vec::with_capacity(tiles.len());
    for tile in &tiles {
        coordinate_offsets.push(coordinates.len() as i64);
        let mut values = Vec::with_capacity(tile.len() * 2);
        for cell in tile.iter() {
            values.push(cell.row);
            values.push(cell.column);
        }
        coordinates.extend(zlib_deflate(&delta_encode(&values, 2), 15));
    }
    put_i64(&mut tile_offsets_section, coordinate_offsets.len() as i64);
    for offset in &coordinate_offsets {
        put_i64(&mut tile_offsets_section, *offset);
    }
    files.push(("__coords.tdb".to_string(), coordinates));

    book_keeping.extend_from_slice(&tile_offsets_section);
    book_keeping.extend_from_slice(&var_offsets_section);
    book_keeping.extend_from_slice(&var_sizes_section);
    // `last_tile_cell_num_`, with a full last tile recorded as itself.
    put_i64(&mut book_keeping, tiles[tiles.len() - 1].len() as i64);

    // windowBits 15 + 16 is the gzip wrapper with a zero time and the OS byte of the build.
    let mut gzipped = zlib_deflate(&book_keeping, 15 + 16);
    gzipped[9] = 3;
    files.push(("__book_keeping.tdb.gz".to_string(), gzipped));
    files.push(("__tiledb_fragment.tdb".to_string(), Vec::new()));
    files
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_tuple_field_keeps_each_row_the_same_length() {
        let parsed = parse_multi_d(
            "AS_RAW_BaseQRankSum",
            "|0.5,1,1.5|",
            &[CellType::Float32, CellType::Int32],
        )
        .unwrap();
        // Empty, then a row of two pairs' worth, then empty, each empty row holding one missing.
        let floats = &parsed[0];
        assert_eq!(u64::from_le_bytes(floats[0..8].try_into().unwrap()), 16);
        assert_eq!(
            u32::from_le_bytes(floats[8..12].try_into().unwrap()),
            BCF_FLOAT_MISSING
        );
    }

    #[test]
    fn a_comma_list_is_joined_back_with_its_commas() {
        let vc = VariantContext::new(
            "chr1",
            1,
            vec![htsjdk_vcf::allele::Allele::from_str("A", true).unwrap()],
        );
        let encoded = encode_site(
            "AS_SB_TABLE",
            Some(&Value::List(vec![
                Value::Str("1".into()),
                Value::Str("2|3".into()),
                Value::Str("4".into()),
            ])),
            LineType::String,
            &Cardinality::Fixed(1),
            &vc,
        )
        .unwrap();
        assert_eq!(encoded, Encoded::Text("1,2|3,4".to_string()));
    }

    #[test]
    fn the_variant_types_are_htslibs() {
        assert_eq!(variant_type(b"AT", b"A"), VariantType::Indel);
        assert_eq!(variant_type(b"A", b"*"), VariantType::SpanningDeletion);
        assert_eq!(variant_type(b"A", b"<NON_REF>"), VariantType::NonRef);
        assert_eq!(variant_type(b"AT", b"CG"), VariantType::Other);
    }
}
