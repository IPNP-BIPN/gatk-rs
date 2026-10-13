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
