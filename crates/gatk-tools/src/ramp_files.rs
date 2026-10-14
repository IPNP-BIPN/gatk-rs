//! The state files `RampedHaplotypeCaller`'s off ramps write.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.OffRampBase`,
//! `PreFilterOffRamp` and `AssemblerOffRamp`, and from
//! `org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter` with its
//! `HaplotypeBAMDestination` and `SAMFileDestination`, at GATK 4.6.2.0.
//!
//! An off ramp is a zip. Each region it is handed adds an object to the `regions` array of an
//! `info.json` written last, and, for the pre-filter ramp, CSV entries under a directory named for
//! the region; the assembler ramps add the region's haplotypes to one `haplotypes.bam`, sorted by
//! coordinate, and its index. The run itself goes on, every region answering no calls.
//!
//! `info.json` is the text the JSON library GATK bundles writes with an indent of two, and that
//! library's licence is not one this repository can port from. The layout here is what its output
//! shows, measured against the reference's files rather than taken from its source: an object's
//! keys come out in the iteration order of a `java.util.HashMap` of strings (see
//! [`gatk_engine::java_hash`] for what that order rests on), an object or array with exactly one
//! member is written on one line, and anything longer puts each member on its own line, indented.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_hash::{byte_array_hash_code, JavaHashMap};
use gatk_engine::tsv_table::java_double_to_string;
use htsjdk_bam::header::{ProgramRecord, ReadGroup, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// A value in `info.json`.
#[derive(Debug, Clone)]
pub enum Json {
    Int(i64),
    Str(String),
    Object(JsonObject),
    Array(Vec<Json>),
}

/// A JSON object whose keys iterate as a `java.util.HashMap`'s would.
#[derive(Debug, Clone, Default)]
pub struct JsonObject(JavaHashMap<String, Json>);

impl JsonObject {
    pub fn new() -> JsonObject {
        JsonObject(JavaHashMap::new())
    }

    /// `put`: a new key takes its bucket's last place, an existing one keeps its own.
    pub fn put(&mut self, key: &str, value: Json) {
        self.0.insert(key.to_string(), value);
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// `OffRampBase.addInfo`: a dotted name walks down nested objects, making the missing ones.
    pub fn put_path(&mut self, name: &str, value: Json) {
        let mut parts: Vec<&str> = name.split('.').collect();
        let last = parts.pop().expect("a name has at least one part");
        let mut object = self;
        for part in parts {
            if !object.0.contains_key(&part.to_string()) {
                object.put(part, Json::Object(JsonObject::new()));
            }
            object = match object.0.get_mut(&part.to_string()) {
                Some(Json::Object(inner)) => inner,
                _ => unreachable!("a path step is an object"),
            };
        }
        object.put(last, value);
    }

    /// The text, at an indent of two.
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        write_object(self, 2, 0, &mut out);
        out
    }
}

fn quote(text: &str, out: &mut String) {
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

fn write_value(value: &Json, factor: usize, indent: usize, out: &mut String) {
    match value {
        Json::Int(i) => out.push_str(&i.to_string()),
        Json::Str(s) => quote(s, out),
        Json::Object(o) => write_object(o, factor, indent, out),
        Json::Array(a) => write_array(a, factor, indent, out),
    }
}

fn write_object(object: &JsonObject, factor: usize, indent: usize, out: &mut String) {
    out.push('{');
    let entries: Vec<(&String, &Json)> = object.0.iter().collect();
    if entries.len() == 1 {
        quote(entries[0].0, out);
        out.push_str(": ");
        write_value(entries[0].1, factor, indent, out);
    } else if !entries.is_empty() {
        let inner = indent + factor;
        for (position, (key, value)) in entries.iter().enumerate() {
            if position > 0 {
                out.push(',');
            }
            out.push('\n');
            out.push_str(&" ".repeat(inner));
            quote(key, out);
            out.push_str(": ");
            write_value(value, factor, inner, out);
        }
        out.push('\n');
        out.push_str(&" ".repeat(indent));
    }
    out.push('}');
}

fn write_array(array: &[Json], factor: usize, indent: usize, out: &mut String) {
    out.push('[');
    if array.len() == 1 {
        write_value(&array[0], factor, indent, out);
    } else if !array.is_empty() {
        let inner = indent + factor;
        for (position, value) in array.iter().enumerate() {
            if position > 0 {
                out.push(',');
            }
            out.push('\n');
            out.push_str(&" ".repeat(inner));
            write_value(value, factor, inner, out);
        }
        out.push('\n');
        out.push_str(&" ".repeat(indent));
    }
    out.push(']');
}

/// `AssemblyRegion.toString()`: its span, whether it is active, and how many reads it holds.
pub fn region_description(region: &AssemblyRegion) -> String {
    format!(
        "AssemblyRegion {} active?={} nReads={}",
        interval_text(region.span()),
        region.is_active(),
        region.reads().len()
    )
}

/// `SimpleInterval.toString()`.
fn interval_text(interval: &SimpleInterval) -> String {
    format!("{}:{}-{}", interval.contig, interval.start, interval.end)
}

/// `RampBase.getLocFilenameSuffix`.
fn location_suffix(location: &SimpleInterval) -> String {
    format!("{}-{}-{}", location.contig, location.start, location.end)
}

/// Which off ramp a file is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OffRampType {
    PreFilter,
    PreAssembler,
    PostAssembler,
}

impl OffRampType {
    /// `OffRampTypeEnum`'s names; `NONE` is no ramp.
    pub fn parse(name: &str) -> Option<Option<OffRampType>> {
        match name {
            "NONE" => Some(None),
            "PRE_FILTER_OFF" => Some(Some(OffRampType::PreFilter)),
            "PRE_ASSEMBLER_OFF" => Some(Some(OffRampType::PreAssembler)),
            "POST_ASSEMBLER_OFF" => Some(Some(OffRampType::PostAssembler)),
            _ => None,
        }
    }
}

/// An off ramp being filled: the regions' objects, the entries, and the haplotype reads.
pub struct OffRamp {
    pub kind: OffRampType,
    regions: Vec<Json>,
    entries: Vec<(String, Vec<u8>)>,
    haplotype_header: Option<SamHeader>,
    haplotype_records: Vec<BamRecord>,
    /// `HaplotypeBAMWriter.uniqueNameCounter`, which starts at one.
    name_counter: i64,
    /// The BGZF level and deflater the haplotype BAM is written with: GATK's, as for any BAM
    /// a run writes, and the index's offsets follow from them.
    compression: (u32, htsjdk_bgzf::Deflater),
}

/// `HaplotypeBAMWriter.DEFAULT_HAPLOTYPE_READ_GROUP_ID`.
const HAPLOTYPE_READ_GROUP: &str = "ArtificialHaplotypeRG";
/// `HaplotypeBAMDestination.haplotypeSampleTag`.
const HAPLOTYPE_SAMPLE_TAG: &str = "HC";

impl OffRamp {
    pub fn new(kind: OffRampType, compression: (u32, htsjdk_bgzf::Deflater)) -> OffRamp {
        OffRamp {
            compression,
            kind,
            regions: Vec::new(),
            entries: Vec::new(),
            haplotype_header: None,
            haplotype_records: Vec::new(),
            name_counter: 1,
        }
    }

    fn region_object(location: &SimpleInterval) -> JsonObject {
        let mut region = JsonObject::new();
        region.put("contig", Json::Str(location.contig.clone()));
        region.put("start", Json::Int(location.start as i64));
        region.put("end", Json::Int(location.end as i64));
        region
    }

    fn add_entry(&mut self, location: Option<&SimpleInterval>, name: &str, bytes: Vec<u8>) {
        let prefix = location.map_or(String::new(), |l| format!("{}/", location_suffix(l)));
        self.entries.push((format!("{prefix}{name}"), bytes));
    }

    /// `OffRampBase.addHaplotypes`: the haplotypes as CSV, one row each.
    fn add_haplotypes(&mut self, location: &SimpleInterval, name: &str, haplotypes: &[Haplotype]) {
        let mut text =
            String::from("contig,start,end,ref,cigar,bases,score,alignmentStartHapwrtRef\n");
        for haplotype in haplotypes {
            let place = haplotype
                .genome_location()
                .expect("an assembled haplotype has a location");
            text.push_str(&format!(
                "{},{},{},{},{},{},{},{}\n",
                place.contig,
                place.start,
                place.end,
                if haplotype.is_reference() { 1 } else { 0 },
                haplotype.cigar().map_or(String::new(), |c| c.to_text()),
                String::from_utf8_lossy(&haplotype.bases()),
                java_double_to_string(haplotype.score()),
                haplotype.alignment_start_hap_wrt_ref(),
            ));
        }
        self.add_entry(Some(location), name, text.into_bytes());
    }

    /// `PreFilterOffRamp.add`: the likelihoods' haplotypes, and for each sample its reads and its
    /// matrix, under the region; the counts and the genotyping spans in `info.json`.
    pub fn add_likelihoods(
        &mut self,
        location: &SimpleInterval,
        name: &str,
        likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>,
        region_for_genotyping: &AssemblyRegion,
    ) {
        let mut region = OffRamp::region_object(location);
        region.put_path(
            &format!("{name}.haplotypeCount"),
            Json::Int(likelihoods.number_of_alleles() as i64),
        );
        region.put_path(
            &format!("{name}.readCount"),
            Json::Int(likelihoods.evidence_count() as i64),
        );
        region.put_path(
            &format!("{name}.r4g"),
            Json::Str(interval_text(region_for_genotyping.span())),
        );
        region.put_path(
            &format!("{name}.r4gp"),
            Json::Str(interval_text(region_for_genotyping.padded_span())),
        );
        let haplotypes: Vec<Haplotype> = (0..likelihoods.number_of_alleles())
            .filter_map(|index| likelihoods.get_allele(index).cloned())
            .collect();
        self.add_haplotypes(location, &format!("{name}.haplotypes"), &haplotypes);
        for sample_index in 0..likelihoods.number_of_samples() {
            let sample = likelihoods
                .get_sample(sample_index)
                .cloned()
                .unwrap_or_default();
            let base = format!("{name}.samples.{sample}");
            let reads = likelihoods.sample_evidence(sample_index).unwrap_or(&[]);
            let mut counts = JsonObject::new();
            counts.put("readCount", Json::Int(reads.len() as i64));
            region.put_path(&base, Json::Object(counts));

            let supplementary = |read: &BamRecord| if read.flags & 0x800 != 0 { 1 } else { 0 };
            let mut names = String::from("name,supp\n");
            for read in reads {
                names.push_str(&format!("{},{}\n", read.read_name, supplementary(read)));
            }
            self.add_entry(Some(location), &format!("{base}.reads"), names.into_bytes());

            let mut matrix = String::from("read,supp");
            for allele in 0..likelihoods.number_of_alleles() {
                matrix.push_str(&format!(",h{allele}"));
            }
            matrix.push('\n');
            for (read_index, read) in reads.iter().enumerate() {
                matrix.push_str(&format!("{},{}", read.read_name, supplementary(read)));
                for allele in 0..likelihoods.number_of_alleles() {
                    matrix.push(',');
                    matrix.push_str(&java_double_to_string(likelihoods.value(
                        sample_index,
                        allele,
                        read_index,
                    )));
                }
                matrix.push('\n');
            }
            self.add_entry(
                Some(location),
                &format!("{base}.matrix"),
                matrix.into_bytes(),
            );
        }
        self.regions.push(Json::Object(region));
    }

    /// `AssemblerOffRamp.add`: the assembly's spans in `info.json` and its haplotypes in the BAM;
    /// before the assembler there is no assembly, and the region's object holds its place alone.
    #[allow(clippy::too_many_arguments)]
    pub fn add_assembly(
        &mut self,
        location: &SimpleInterval,
        description: &str,
        name: &str,
        assembly: Option<(&SimpleInterval, &AssemblyRegion, &[Haplotype])>,
        left_flank: Option<&AssemblyRegion>,
        right_flank: Option<&AssemblyRegion>,
        header: &SamHeader,
    ) {
        let mut region = OffRamp::region_object(location);
        let mut assembler = JsonObject::new();
        let spans = |r: &AssemblyRegion| {
            let mut object = JsonObject::new();
            object.put("activeSpan", Json::Str(interval_text(r.span())));
            object.put("paddedSpan", Json::Str(interval_text(r.padded_span())));
            Json::Object(object)
        };
        if let Some((padded, genotyping, _)) = assembly {
            assembler.put("paddedReferenceLoc", Json::Str(interval_text(padded)));
            assembler.put("regionForGenotyping", spans(genotyping));
        }
        if let Some(flank) = left_flank {
            assembler.put("nonVariantLeftFlankRegion", spans(flank));
        }
        if let Some(flank) = right_flank {
            assembler.put("nonVariantRightFlankRegion", spans(flank));
        }
        if let Some((padded, _, haplotypes)) = assembly {
            self.write_haplotypes(location, description, padded, haplotypes, header);
        }
        if !assembler.is_empty() {
            region.put(name, Json::Object(assembler));
        }
        self.regions.push(Json::Object(region));
    }

    /// `HaplotypeBAMWriter.writeHaplotypesAsReads` for `ALL_POSSIBLE_HAPLOTYPES` with the special
    /// tag on and no best haplotypes: every haplotype as a reverse-strand read at mapping quality
    /// zero, named and tagged as the writer names and tags them.
    fn write_haplotypes(
        &mut self,
        location: &SimpleInterval,
        description: &str,
        padded: &SimpleInterval,
        haplotypes: &[Haplotype],
        header: &SamHeader,
    ) {
        let output_header = self
            .haplotype_header
            .get_or_insert_with(|| haplotype_bam_header(header))
            .clone();
        let reference_index = output_header
            .sequences
            .iter()
            .position(|s| s.name == padded.contig)
            .map_or(-1, |i| i as i32);
        // The callable region is the assembly region itself: `CR` is its `toString()`, `XS` its
        // coordinates.
        let callable = interval_text(location);
        for (index, haplotype) in haplotypes.iter().enumerate() {
            let bases = haplotype.bases();
            let mut record = BamRecord {
                read_name: format!("{HAPLOTYPE_SAMPLE_TAG}_{}_{index}", self.name_counter),
                flags: 0x10,
                reference_index,
                alignment_start: padded.start + haplotype.alignment_start_hap_wrt_ref(),
                mapping_quality: 0,
                cigar: haplotype.cigar().cloned().unwrap_or_default(),
                read_bases: bases.to_vec(),
                // `Utils.dupBytes((byte) '!', length)` as the qualities: the byte 33, a phred of 33.
                base_qualities: vec![b'!'; bases.len()],
                ..BamRecord::default()
            };
            self.name_counter += 1;
            record.tags.insert(
                Tag::new(b"HC"),
                TagValue::Int(byte_array_hash_code(&bases) as i64),
            );
            record.tags.insert(
                Tag::new(b"RG"),
                TagValue::Str(HAPLOTYPE_READ_GROUP.to_string()),
            );
            record
                .tags
                .insert(Tag::new(b"CR"), TagValue::Str(description.to_string()));
            record.tags.insert(
                Tag::new(b"XS"),
                TagValue::Str(format!(
                    "{},{},{},{}",
                    if haplotype.is_reference() { "1" } else { "0" },
                    java_double_to_string(haplotype.score()),
                    haplotype.alignment_start_hap_wrt_ref(),
                    callable
                )),
            );
            self.haplotype_records.push(record);
        }
    }

    /// The zip's entries in the order the reference writes them: the regions' entries as they were
    /// added, then the haplotype BAM and its index, then `info.json`.
    pub fn entries(&self) -> Result<Vec<(String, Vec<u8>)>, String> {
        let mut entries = self.entries.clone();
        if let Some(header) = &self.haplotype_header {
            let mut records = self.haplotype_records.clone();
            records.sort_by(htsjdk_bam::coordinate::compare);
            // The writer indexes as it writes (`createBamOutIndex`), and an index built that way
            // ends a chunk inside the block its last record closed, where one built from the
            // finished file points at the start of the next.
            let mut writer = htsjdk_bam::writer::BamWriter::with_compression(
                Vec::new(),
                header,
                self.compression.0,
                self.compression.1,
            )
            .map_err(|e| format!("{e:?}"))?
            .with_index();
            for record in &records {
                writer.write(record).map_err(|e| format!("{e:?}"))?;
            }
            let (bam, bai) = writer.finish_with_index().map_err(|e| format!("{e:?}"))?;
            entries.push(("haplotypes.bam".to_string(), bam));
            entries.push(("haplotypes.bai".to_string(), bai));
        }
        let mut info = JsonObject::new();
        info.put("regions", Json::Array(self.regions.clone()));
        entries.push(("info.json".to_string(), info.to_text().into_bytes()));
        Ok(entries)
    }

    /// `close()`: the zip, deflated, at `path`.
    pub fn write(&self, path: &std::path::Path) -> Result<(), String> {
        use std::io::Write;
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
        let mut zip = zip::ZipWriter::new(file);
        let options = zip::write::SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Deflated);
        for (name, bytes) in self.entries()? {
            zip.start_file(name, options).map_err(|e| e.to_string())?;
            zip.write_all(&bytes).map_err(|e| e.to_string())?;
        }
        zip.finish().map_err(|e| e.to_string())?;
        Ok(())
    }
}

/// `HaplotypeBAMDestination`'s header: the source's dictionary, coordinate order, its read groups
/// and the haplotypes' own, and its program records and the writer's.
fn haplotype_bam_header(source: &SamHeader) -> SamHeader {
    let mut header = SamHeader::new();
    header.set_sort_order("coordinate");
    header.sequences = source.sequences.clone();
    header.read_groups = source.read_groups.clone();
    let mut group = ReadGroup::new(HAPLOTYPE_READ_GROUP);
    group.attributes.set("SM", HAPLOTYPE_SAMPLE_TAG);
    group.attributes.set("CN", "BI");
    header.read_groups.push(group);
    header.programs = source.programs.clone();
    header
        .programs
        .push(ProgramRecord::new("HaplotypeBAMWriter"));
    header
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_region_object_lists_its_keys_in_hash_order_and_a_lone_member_stays_inline() {
        let mut region = JsonObject::new();
        region.put("contig", Json::Str("chr1".to_string()));
        region.put("start", Json::Int(963));
        region.put("end", Json::Int(1037));
        let mut samples = JsonObject::new();
        samples.put_path("tumor.readCount", Json::Int(30));
        region.put("samples", Json::Object(samples));
        let mut info = JsonObject::new();
        info.put("regions", Json::Array(vec![Json::Object(region)]));
        assert_eq!(
            info.to_text(),
            "{\"regions\": [{\n  \"start\": 963,\n  \"end\": 1037,\n  \"samples\": {\"tumor\": {\"readCount\": 30}},\n  \"contig\": \"chr1\"\n}]}"
        );
    }
}

/// Which on ramp a file is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnRampType {
    PostFilter,
    PostAssembler,
}

impl OnRampType {
    /// `OnRampTypeEnum`'s names; `NONE` is no ramp.
    pub fn parse(name: &str) -> Option<Option<OnRampType>> {
        match name {
            "NONE" => Some(None),
            "POST_FILTER_ON" => Some(Some(OnRampType::PostFilter)),
            "POST_ASSEMBLER_ON" => Some(Some(OnRampType::PostAssembler)),
            _ => None,
        }
    }
}

/// A value read back from `info.json`: what the off ramps write, objects in file order.
#[derive(Debug, Clone, PartialEq)]
pub enum JsonValue {
    Int(i64),
    Str(String),
    Object(Vec<(String, JsonValue)>),
    Array(Vec<JsonValue>),
}

impl JsonValue {
    /// `JSONObject.has(key)` then `get(key)`.
    pub fn get(&self, key: &str) -> Option<&JsonValue> {
        match self {
            JsonValue::Object(members) => members.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    /// `getString(key)`.
    pub fn string(&self, key: &str) -> Result<&str, String> {
        match self.get(key) {
            Some(JsonValue::Str(text)) => Ok(text),
            _ => Err(format!("JSONObject[\"{key}\"] not a string.")),
        }
    }

    /// `getInt(key)`.
    fn int(&self, key: &str) -> Result<i64, String> {
        match self.get(key) {
            Some(JsonValue::Int(value)) => Ok(*value),
            _ => Err(format!("JSONObject[\"{key}\"] is not a int.")),
        }
    }
}

/// The JSON the off ramps write: objects, arrays, strings and integers.
fn parse_json(text: &str) -> Result<JsonValue, String> {
    struct Parser<'a> {
        bytes: &'a [u8],
        at: usize,
    }
    impl Parser<'_> {
        fn blank(&mut self) {
            while self.at < self.bytes.len() && self.bytes[self.at].is_ascii_whitespace() {
                self.at += 1;
            }
        }
        fn expect(&mut self, byte: u8) -> Result<(), String> {
            self.blank();
            if self.bytes.get(self.at) == Some(&byte) {
                self.at += 1;
                Ok(())
            } else {
                Err(format!("expected '{}' at {}", byte as char, self.at))
            }
        }
        fn string(&mut self) -> Result<String, String> {
            self.expect(b'"')?;
            let mut out = Vec::new();
            while let Some(&byte) = self.bytes.get(self.at) {
                self.at += 1;
                match byte {
                    b'"' => return String::from_utf8(out).map_err(|e| e.to_string()),
                    b'\\' => {
                        let escaped = *self.bytes.get(self.at).ok_or("an unfinished escape")?;
                        self.at += 1;
                        match escaped {
                            b'n' => out.push(b'\n'),
                            b't' => out.push(b'\t'),
                            b'r' => out.push(b'\r'),
                            b'b' => out.push(8),
                            b'f' => out.push(12),
                            b'u' => {
                                let hex = std::str::from_utf8(
                                    self.bytes.get(self.at..self.at + 4).ok_or("a short \\u")?,
                                )
                                .map_err(|e| e.to_string())?;
                                let code =
                                    u32::from_str_radix(hex, 16).map_err(|e| e.to_string())?;
                                self.at += 4;
                                let c = char::from_u32(code).unwrap_or('\u{fffd}');
                                out.extend_from_slice(c.to_string().as_bytes());
                            }
                            other => out.push(other),
                        }
                    }
                    other => out.push(other),
                }
            }
            Err("an unterminated string".to_string())
        }
        fn value(&mut self) -> Result<JsonValue, String> {
            self.blank();
            match self.bytes.get(self.at) {
                Some(b'{') => {
                    self.at += 1;
                    let mut members = Vec::new();
                    self.blank();
                    if self.bytes.get(self.at) == Some(&b'}') {
                        self.at += 1;
                        return Ok(JsonValue::Object(members));
                    }
                    loop {
                        let key = self.string()?;
                        self.expect(b':')?;
                        members.push((key, self.value()?));
                        self.blank();
                        match self.bytes.get(self.at) {
                            Some(b',') => self.at += 1,
                            Some(b'}') => {
                                self.at += 1;
                                return Ok(JsonValue::Object(members));
                            }
                            _ => return Err(format!("expected ',' or '}}' at {}", self.at)),
                        }
                    }
                }
                Some(b'[') => {
                    self.at += 1;
                    let mut items = Vec::new();
                    self.blank();
                    if self.bytes.get(self.at) == Some(&b']') {
                        self.at += 1;
                        return Ok(JsonValue::Array(items));
                    }
                    loop {
                        items.push(self.value()?);
                        self.blank();
                        match self.bytes.get(self.at) {
                            Some(b',') => self.at += 1,
                            Some(b']') => {
                                self.at += 1;
                                return Ok(JsonValue::Array(items));
                            }
                            _ => return Err(format!("expected ',' or ']' at {}", self.at)),
                        }
                    }
                }
                Some(b'"') => Ok(JsonValue::Str(self.string()?)),
                _ => {
                    let start = self.at;
                    while self.at < self.bytes.len()
                        && (self.bytes[self.at] == b'-' || self.bytes[self.at].is_ascii_digit())
                    {
                        self.at += 1;
                    }
                    std::str::from_utf8(&self.bytes[start..self.at])
                        .map_err(|e| e.to_string())?
                        .parse::<i64>()
                        .map(JsonValue::Int)
                        .map_err(|_| format!("an unexpected value at {start}"))
                }
            }
        }
    }
    let mut parser = Parser {
        bytes: text.as_bytes(),
        at: 0,
    };
    parser.value()
}

/// `new SimpleInterval(String)` for the `contig:start-end` the ramps write.
fn parse_interval(text: &str) -> Result<SimpleInterval, String> {
    let bad = || format!("not an interval: {text}");
    let (contig, range) = text.rsplit_once(':').ok_or_else(bad)?;
    let (start, end) = range.split_once('-').ok_or_else(bad)?;
    let number = |s: &str| s.replace(',', "").parse::<i32>().map_err(|_| bad());
    SimpleInterval::new(contig, number(start)?, number(end)?).ok_or_else(bad)
}

/// `TextCigarCodec.decode`.
fn parse_cigar(text: &str) -> Result<htsjdk_bam::cigar::Cigar, String> {
    use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
    if text == "*" {
        return Ok(Cigar::new(Vec::new()));
    }
    let mut elements = Vec::new();
    let mut length = 0u32;
    for byte in text.bytes() {
        if byte.is_ascii_digit() {
            length = length * 10 + u32::from(byte - b'0');
            continue;
        }
        let op = match byte {
            b'M' => Op::M,
            b'I' => Op::I,
            b'D' => Op::D,
            b'N' => Op::N,
            b'S' => Op::S,
            b'H' => Op::H,
            b'P' => Op::P,
            b'=' => Op::Eq,
            b'X' => Op::X,
            _ => return Err(format!("Malformed CIGAR string: {text}")),
        };
        elements.push(CigarElement { length, op });
        length = 0;
    }
    Ok(Cigar::new(elements))
}

/// `OnRampBase.CSVReader`: a header of column names, then rows split on commas.
struct CsvTable {
    columns: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl CsvTable {
    fn read(bytes: &[u8]) -> Result<CsvTable, String> {
        let text = String::from_utf8_lossy(bytes);
        let mut lines = text.lines();
        let columns = lines
            .next()
            .ok_or("an empty table")?
            .split(',')
            .map(str::to_string)
            .collect();
        let rows = lines
            .map(|line| line.split(',').map(str::to_string).collect())
            .collect();
        Ok(CsvTable { columns, rows })
    }

    /// `getColumn(name)`: the last column of that name, as the `LinkedHashMap` keeps it.
    fn column(&self, name: &str) -> Result<usize, String> {
        self.columns
            .iter()
            .rposition(|c| c == name)
            .ok_or_else(|| format!("java.io.IOException: no such column: {name}"))
    }
}

/// `Double.parseDouble`, for the `Double.toString` the ramps write.
fn parse_double(text: &str) -> Result<f64, String> {
    text.trim()
        .parse::<f64>()
        .map_err(|_| format!("java.lang.NumberFormatException: For input string: \"{text}\""))
}

/// `Integer.parseInt`.
fn parse_int(text: &str) -> Result<i32, String> {
    text.parse::<i32>()
        .map_err(|_| format!("java.lang.NumberFormatException: For input string: \"{text}\""))
}

/// An on ramp: the zip an off ramp wrote, read back. `OnRampBase` with `PostFilterOnRamp` and
/// `PostAssemblerOnRamp`.
pub struct OnRamp {
    pub kind: OnRampType,
    entries: Vec<(String, Vec<u8>)>,
    regions: Vec<JsonValue>,
    /// `regionIndex`: `contig:start-end` to the region's place in `regions`, the last one winning.
    region_index: std::collections::HashMap<String, usize>,
    /// The haplotype BAM's records in file order, with its contig names, for the post-assembler
    /// ramp.
    haplotype_records: Vec<BamRecord>,
    haplotype_contigs: Vec<String>,
}

impl OnRamp {
    /// `new OnRampBase(filename)`: the zip opened and `info.json` read; the post-assembler ramp
    /// also copies out its haplotype BAM and index.
    pub fn open(kind: OnRampType, path: &std::path::Path) -> Result<OnRamp, String> {
        use std::io::Read;
        let file = std::fs::File::open(path).map_err(|e| format!("java.io.IOException: {e}"))?;
        let mut archive =
            zip::ZipArchive::new(file).map_err(|e| format!("java.util.zip.ZipException: {e}"))?;
        let mut entries = Vec::new();
        for index in 0..archive.len() {
            let mut entry = archive.by_index(index).map_err(|e| e.to_string())?;
            let name = entry.name().map_err(|e| e.to_string())?.to_string();
            let mut bytes = Vec::new();
            entry.read_to_end(&mut bytes).map_err(|e| e.to_string())?;
            entries.push((name, bytes));
        }
        let mut ramp = OnRamp {
            kind,
            entries,
            regions: Vec::new(),
            region_index: std::collections::HashMap::new(),
            haplotype_records: Vec::new(),
            haplotype_contigs: Vec::new(),
        };
        let info = parse_json(&String::from_utf8_lossy(ramp.entry(None, "info.json")?))?;
        let Some(JsonValue::Array(regions)) = info.get("regions") else {
            return Err("JSONObject[\"regions\"] is not a JSONArray.".to_string());
        };
        for (index, region) in regions.iter().enumerate() {
            let key = format!(
                "{}:{}-{}",
                region.string("contig")?,
                region.int("start")?,
                region.int("end")?
            );
            ramp.region_index.insert(key, index);
        }
        ramp.regions = regions.clone();
        if kind == OnRampType::PostAssembler {
            let bam = ramp.entry(None, "haplotypes.bam")?.to_vec();
            ramp.entry(None, "haplotypes.bai")?;
            let plain = htsjdk_bgzf::read::decompress_all(&bam).map_err(|e| format!("{e:?}"))?;
            let reader =
                htsjdk_bam::reader::BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
            ramp.haplotype_contigs = reader
                .header
                .text
                .sequences
                .iter()
                .map(|r| r.name.clone())
                .collect();
            for record in reader {
                ramp.haplotype_records
                    .push(record.map_err(|e| format!("{e:?}"))?);
            }
        }
        Ok(ramp)
    }

    /// `getEntry(loc, name)`: the first entry of that name, under the region's directory.
    fn entry(&self, location: Option<&SimpleInterval>, name: &str) -> Result<&[u8], String> {
        let prefix = location.map_or(String::new(), |l| format!("{}/", location_suffix(l)));
        let full = format!("{prefix}{name}");
        self.entries
            .iter()
            .find(|(n, _)| *n == full)
            .map(|(_, bytes)| bytes.as_slice())
            .ok_or_else(|| format!("java.io.IOException: no such: {full}"))
    }

    /// `hasRegion(loc)`.
    pub fn has_region(&self, location: &SimpleInterval) -> bool {
        self.region_index.contains_key(&interval_text(location))
    }

    /// `getRegion(loc, name)`: the region's object, then down the dotted path; nothing where a
    /// level is missing.
    pub fn region(&self, location: &SimpleInterval, name: &str) -> Option<&JsonValue> {
        let index = self.region_index.get(&interval_text(location))?;
        let mut object = self.regions.get(*index)?;
        for key in name.split('.') {
            object = object.get(key)?;
        }
        Some(object)
    }

    /// `getHaplotypes(loc, name)`: the haplotype table, one haplotype per row.
    pub fn haplotypes(
        &self,
        location: &SimpleInterval,
        name: &str,
    ) -> Result<Vec<Haplotype>, String> {
        let table = CsvTable::read(self.entry(Some(location), name)?)?;
        let contig = table.column("contig")?;
        let start = table.column("start")?;
        let end = table.column("end")?;
        let reference = table.column("ref")?;
        let cigar = table.column("cigar")?;
        let bases = table.column("bases")?;
        let score = table.column("score")?;
        let offset = table.column("alignmentStartHapwrtRef")?;
        let mut haplotypes = Vec::new();
        for row in &table.rows {
            let cell = |i: usize| -> Result<&str, String> {
                row.get(i).map(String::as_str).ok_or_else(|| {
                    format!("java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}", row.len())
                })
            };
            let mut haplotype =
                Haplotype::new(cell(bases)?.as_bytes(), parse_int(cell(reference)?)? != 0)
                    .map_err(|e| format!("{e:?}"))?;
            haplotype.set_genome_location(
                SimpleInterval::new(
                    cell(contig)?,
                    parse_int(cell(start)?)?,
                    parse_int(cell(end)?)?,
                )
                .ok_or("a bad haplotype location")?,
            );
            haplotype.set_score(parse_double(cell(score)?)?);
            haplotype
                .set_cigar(&parse_cigar(cell(cigar)?)?)
                .map_err(|e| format!("{e:?}"))?;
            haplotype.set_alignment_start_hap_wrt_ref(parse_int(cell(offset)?)?);
            haplotypes.push(haplotype);
        }
        Ok(haplotypes)
    }

    /// `PostAssemblerOnRamp.readHaplotypes(region)`: the haplotype reads overlapping the region,
    /// in file order, keeping those whose special tag names this very region.
    pub fn bam_haplotypes(&self, location: &SimpleInterval) -> Result<Vec<Haplotype>, String> {
        let wanted = interval_text(location);
        let mut haplotypes = Vec::new();
        for record in &self.haplotype_records {
            let Some(contig) = usize::try_from(record.reference_index)
                .ok()
                .and_then(|i| self.haplotype_contigs.get(i))
            else {
                continue;
            };
            let end = record.alignment_end();
            if *contig != location.contig
                || end < location.start
                || record.alignment_start > location.end
            {
                continue;
            }
            let mut is_reference = false;
            let mut score = f64::NAN;
            let mut offset = 0;
            if let Some(TagValue::Str(special)) = record.tags.get(Tag::new(b"XS")) {
                let tokens: Vec<&str> = special.split(',').collect();
                is_reference = parse_int(tokens.first().copied().unwrap_or(""))? != 0;
                score = parse_double(tokens.get(1).copied().unwrap_or(""))?;
                offset = parse_int(tokens.get(2).copied().unwrap_or(""))?;
                if tokens.get(3).copied() != Some(wanted.as_str()) {
                    continue;
                }
            }
            let mut haplotype =
                Haplotype::new(&record.read_bases, is_reference).map_err(|e| format!("{e:?}"))?;
            haplotype.set_genome_location(
                SimpleInterval::new(contig, record.alignment_start, end)
                    .ok_or("a bad haplotype location")?,
            );
            haplotype
                .set_cigar(&record.cigar)
                .map_err(|e| format!("{e:?}"))?;
            haplotype.set_alignment_start_hap_wrt_ref(offset);
            if !score.is_nan() {
                haplotype.set_score(score);
            }
            haplotypes.push(haplotype);
        }
        Ok(haplotypes)
    }

    /// `PostAssemblerOnRamp.parseAssemblyRegion`: an active region of the object's two spans,
    /// with no reads.
    pub fn assembly_region(object: &JsonValue) -> Result<AssemblyRegion, String> {
        AssemblyRegion::new(
            parse_interval(object.string("activeSpan")?)?,
            parse_interval(object.string("paddedSpan")?)?,
            true,
        )
        .map_err(|e| format!("{e:?}"))
    }

    /// `getRegion(span, name).getString(key)` as an interval.
    pub fn region_interval(
        &self,
        location: &SimpleInterval,
        name: &str,
        key: &str,
    ) -> Result<SimpleInterval, String> {
        let object = self
            .region(location, name)
            .ok_or("java.lang.NullPointerException")?;
        parse_interval(object.string(key)?)
    }

    /// `PostFilterOnRamp.getAlleleLikelihoods`: the haplotypes, then per sample the reads the
    /// table names (looked up by name and supplementary flag, first among the genotyping reads
    /// then among all the region's, a read found in neither skipped) and their matrix.
    pub fn allele_likelihoods(
        &self,
        location: &SimpleInterval,
        name: &str,
        samples: &[String],
        reads: &[Vec<BamRecord>],
        all_reads: &[Vec<BamRecord>],
    ) -> Result<AlleleLikelihoods<BamRecord, Haplotype>, String> {
        let haplotypes = self.haplotypes(location, &format!("{name}.haplotypes"))?;
        let supp_name = |read_name: &str, supplementary: bool| {
            format!("{read_name},{}", if supplementary { 1 } else { 0 })
        };
        let by_name = |groups: &[Vec<BamRecord>]| {
            let mut map: std::collections::HashMap<String, BamRecord> =
                std::collections::HashMap::new();
            for read in groups.iter().flatten() {
                map.insert(
                    supp_name(&read.read_name, read.flags & 0x800 != 0),
                    read.clone(),
                );
            }
            map
        };
        let first = by_name(reads);
        let second = by_name(all_reads);
        let mut evidence = Vec::new();
        let mut values = Vec::new();
        for sample in samples {
            let base = format!("{name}.samples.{sample}");
            let table = CsvTable::read(self.entry(Some(location), &format!("{base}.reads"))?)?;
            let name_column = table.column("name")?;
            let supp_column = table.column("supp")?;
            let mut enriched = Vec::new();
            for row in &table.rows {
                let key = supp_name(
                    row.get(name_column).map_or("", String::as_str),
                    parse_int(row.get(supp_column).map_or("", String::as_str))? != 0,
                );
                if let Some(read) = first.get(&key).or_else(|| second.get(&key)) {
                    enriched.push(read.clone());
                }
            }

            let matrix = CsvTable::read(self.entry(Some(location), &format!("{base}.matrix"))?)?;
            const HAPLOTYPE_START_COLUMN: usize = 2;
            if matrix.columns.len() != haplotypes.len() + HAPLOTYPE_START_COLUMN {
                return Err(format!(
                    "java.io.IOException: column count {} does not match number of haplotypes {}",
                    matrix.columns.len(),
                    haplotypes.len()
                ));
            }
            let mut sample_values = vec![vec![0.0; enriched.len()]; haplotypes.len()];
            let mut row_index = 0;
            for row in &matrix.rows {
                let read_name = row.first().map_or("", String::as_str);
                let expected: &BamRecord = enriched.get(row_index).ok_or_else(|| {
                    format!(
                        "java.lang.IndexOutOfBoundsException: Index {row_index} out of bounds for length {}",
                        enriched.len()
                    )
                })?;
                // IGNORE_DANGLING_READS: a row for a read that was not found is passed over.
                if read_name != expected.read_name {
                    continue;
                }
                for (haplotype, values) in sample_values.iter_mut().enumerate() {
                    values[row_index] = parse_double(
                        row.get(haplotype + HAPLOTYPE_START_COLUMN)
                            .map_or("", String::as_str),
                    )?;
                }
                row_index += 1;
            }
            if row_index != enriched.len() {
                return Err(format!(
                    "java.io.IOException: {row_index} rows found, but there are {} reads",
                    enriched.len()
                ));
            }
            evidence.push(enriched);
            values.push(sample_values);
        }
        AlleleLikelihoods::new(
            gatk_engine::allele_list::SampleList::new(samples),
            gatk_engine::allele_list::AlleleList::new(&haplotypes),
            evidence,
            values,
        )
        .map_err(|e| format!("{e:?}"))
    }
}

/// `RampUtils.compareHaplotypes`: the same count, and pairwise equal once both are sorted by
/// `HaplotypeComparator` (reference flag, contig, start, end, score, bases).
pub fn compare_haplotypes(
    first: &[gatk_engine::haplotype::Haplotype],
    second: &[gatk_engine::haplotype::Haplotype],
) -> Result<(), String> {
    if first.len() != second.len() {
        return Err(format!(
            "java.lang.RuntimeException: haplotype size verification failed: {} vs {}",
            first.len(),
            second.len()
        ));
    }
    let compare = |a: &gatk_engine::haplotype::Haplotype, b: &gatk_engine::haplotype::Haplotype| {
        let location = |h: &gatk_engine::haplotype::Haplotype| {
            h.genome_location()
                .map(|l| (l.contig.clone(), l.start, l.end))
                .unwrap_or_default()
        };
        let (ca, sa, ea) = location(a);
        let (cb, sb, eb) = location(b);
        (a.is_reference() as i32 - b.is_reference() as i32)
            .cmp(&0)
            .then_with(|| gatk_engine::java_hash::compare_strings(&ca, &cb))
            .then_with(|| sa.cmp(&sb))
            .then_with(|| ea.cmp(&eb))
            .then_with(|| {
                let delta = a.score() - b.score();
                if delta < 0.0 {
                    std::cmp::Ordering::Less
                } else if delta > 0.0 {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Equal
                }
            })
            .then_with(|| a.bases().cmp(&b.bases()))
    };
    let mut l1: Vec<_> = first.iter().collect();
    let mut l2: Vec<_> = second.iter().collect();
    l1.sort_by(|a, b| compare(a, b));
    l2.sort_by(|a, b| compare(a, b));
    for (i, (a, b)) in l1.iter().zip(&l2).enumerate() {
        if compare(a, b) != std::cmp::Ordering::Equal {
            return Err(format!(
                "java.lang.RuntimeException: haplotype failed verification on index {i}"
            ));
        }
    }
    Ok(())
}

/// `RampUtils.compareReads`: the same count, and pairwise equal once both are sorted by
/// `GATKReadComparator` (strand, name and span, bases, qualities, soft, aligned and unclipped
/// ends).
pub fn compare_reads(
    first: &[htsjdk_bam::record::BamRecord],
    second: &[htsjdk_bam::record::BamRecord],
    header: &htsjdk_bam::header::SamHeader,
) -> Result<(), String> {
    use gatk_engine::read_utils as ru;
    if first.len() != second.len() {
        return Err(format!(
            "java.lang.RuntimeException: reads size verification failed: {} vs {}",
            first.len(),
            second.len()
        ));
    }
    let common = |r: &htsjdk_bam::record::BamRecord| {
        let name = r.read_name.clone();
        if r.flags & 0x4 != 0 || r.cigar.elements.is_empty() {
            format!("{name} UNMAPPED")
        } else {
            let contig = header
                .sequences
                .get(r.reference_index.max(0) as usize)
                .map(|s| s.name.clone())
                .unwrap_or_default();
            format!("{name} {contig}:{}-{}", ru::start(r), ru::end(r))
        }
    };
    let compare = |a: &htsjdk_bam::record::BamRecord, b: &htsjdk_bam::record::BamRecord| {
        let reverse = |r: &htsjdk_bam::record::BamRecord| (r.flags & 0x10 != 0) as i32;
        (reverse(a) - reverse(b))
            .cmp(&0)
            .then_with(|| gatk_engine::java_hash::compare_strings(&common(a), &common(b)))
            .then_with(|| a.read_bases.cmp(&b.read_bases))
            .then_with(|| a.base_qualities.cmp(&b.base_qualities))
            .then_with(|| ru::soft_start(a).cmp(&ru::soft_start(b)))
            .then_with(|| ru::soft_end(a).cmp(&ru::soft_end(b)))
            .then_with(|| ru::start(a).cmp(&ru::start(b)))
            .then_with(|| ru::end(a).cmp(&ru::end(b)))
            .then_with(|| ru::unclipped_start(a).cmp(&ru::unclipped_start(b)))
            .then_with(|| ru::unclipped_end(a).cmp(&ru::unclipped_end(b)))
    };
    let mut l1: Vec<_> = first.iter().collect();
    let mut l2: Vec<_> = second.iter().collect();
    l1.sort_by(|a, b| compare(a, b));
    l2.sort_by(|a, b| compare(a, b));
    for (i, (a, b)) in l1.iter().zip(&l2).enumerate() {
        if compare(a, b) != std::cmp::Ordering::Equal {
            return Err(format!(
                "java.lang.RuntimeException: reads failed verification on index {i}"
            ));
        }
    }
    Ok(())
}
