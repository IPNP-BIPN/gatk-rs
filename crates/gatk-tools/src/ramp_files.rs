//! The state files `RampedHaplotypeCaller`'s off ramps write.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ramps.OffRampBase`,
//! `PreFilterOffRamp` and `AssemblerOffRamp`, and from
//! `org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter` with its
//! `HaplotypeBAMDestination` and `SAMFileDestination`, at GATK 4.7.0.0; those files, and the on
//! ramps' `OnRampBase`, `PostFilterOnRamp` and `PostAssemblerOnRamp`, are unchanged since 4.6.2.0.
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
//!
//! The second half of the file reads such a zip back ([`OnRamp`]): `info.json` with a small parser
//! of its own, the CSV entries split on commas, and the haplotype BAM decoded in memory.
//!
//! # Reading this file without reading Rust
//!
//! - `&str` and `&[T]` are borrowed, read-only views of a string and of an array; the caller keeps
//!   ownership. `&mut self` is a method that may change its object, and only one such borrow can
//!   exist at a time. `String` and `Vec<T>` are the owned, growable versions;
//! - `Option<T>` is "a value or nothing" (Java's nullable reference, checked by the compiler) and
//!   `Result<T, String>` "a value or an error message"; this file models each Java exception as an
//!   error string that starts with the exception's class name;
//! - `?` after a call returns its error to the caller at once, Java's `throws` at one call site;
//!   `.map_err(..)` rewrites that error first, `.ok_or(..)` turns a "nothing" into an error;
//! - `let Some(x) = y else { ... };` takes the value out of `y` or runs the block, which must leave
//!   (`continue`, `return`); Java has no single construct for it;
//! - `|a| ...` is a closure, Java's lambda; `.iter().map(..).collect()` is Java's stream, map and
//!   collect, and `.then_with(..)` is `Comparator.thenComparing`, consulted only on a tie;
//! - `match` must cover every case of what it inspects, so an `enum`'s cases are listed in full;
//! - `x as i64` is a numeric cast; widening casts such as `i32` to `i64` are exact;
//! - a function's last expression is its return value; `return` is only written for early exits.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::assembly_region::AssemblyRegion;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_hash::{byte_array_hash_code, JavaHashMap};
use gatk_engine::tsv_table::java_double_to_string;
use htsjdk_bam::header::{ProgramRecord, ReadGroup, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// A value in `info.json`, as an off ramp builds it: the four kinds the ramps `put` (an `int`,
/// a `String`, a `JSONObject`, a `JSONArray`). An interval is put as its `toString()`, a string.
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

    /// `JSONObject.put`: a new key takes its bucket's last place, an existing one keeps its own.
    /// That is `java.util.HashMap.put`, which decides the order the keys are written in.
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
    ///
    /// ```java
    /// if (!parent.has(toks[tokIndex])) {
    ///     parent.put(toks[tokIndex], new JSONObject());
    /// }
    /// parent = parent.getJSONObject(toks[tokIndex]);
    /// ```
    ///
    /// Java's `name.split("\\.")` drops trailing empty parts and Rust's `split('.')` keeps them;
    /// the names the ramps pass (`<name>.readCount`, `<name>.samples.<sample>`) end in no dot.
    pub fn put_path(&mut self, name: &str, value: Json) {
        let mut parts: Vec<&str> = name.split('.').collect();
        let last = parts.pop().expect("a name has at least one part");
        // `object` is a writable borrow that is moved one level down on each step, the Java's
        // `parent = parent.getJSONObject(...)`.
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

    /// The text, at an indent of two: `JSONObject.toString(2)`, as `OffRampBase.close` writes it.
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        write_object(self, 2, 0, &mut out);
        out
    }
}

/// A JSON string literal. The escapes are those the reference's files show; the strings the ramps
/// write (contig names, intervals, sample names) are not expected to need any.
fn quote(text: &str, out: &mut String) {
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            // A guard: this arm takes any other control character, written as `\u` and four
            // lower-case hex digits.
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

/// One value at the current indent. `factor` is the indent step (two), `indent` the column the
/// enclosing member starts at.
fn write_value(value: &Json, factor: usize, indent: usize, out: &mut String) {
    match value {
        Json::Int(i) => out.push_str(&i.to_string()),
        Json::Str(s) => quote(s, out),
        Json::Object(o) => write_object(o, factor, indent, out),
        Json::Array(a) => write_array(a, factor, indent, out),
    }
}

/// An object: `{}` when empty, `{"key": value}` on one line for one member, otherwise one member
/// per line at `indent + factor` and the closing brace back at `indent`. A lone member is written
/// at the enclosing indent, not one step in, so its own children line up with the parent's.
fn write_object(object: &JsonObject, factor: usize, indent: usize, out: &mut String) {
    out.push('{');
    // The `HashMap`'s iteration order, collected once so the count is known before writing.
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

/// An array, laid out by the same rule as [`write_object`]'s.
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
///
/// ```java
/// return "AssemblyRegion "  + activeSpan.toString() + " active?=" + isActive
///     + " nReads=" + reads.size();
/// ```
///
/// (One Java line, `AssemblyRegion.java:157`, wrapped here at the `+`.)
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

/// `RampBase.getLocFilenameSuffix`: `contig-start-end`, the directory a region's entries go in.
fn location_suffix(location: &SimpleInterval) -> String {
    format!("{}-{}-{}", location.contig, location.start, location.end)
}

/// Which off ramp a file is: `PreFilter` is `PreFilterOffRamp`, the two assembler ramps are both
/// `AssemblerOffRamp`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OffRampType {
    PreFilter,
    PreAssembler,
    PostAssembler,
}

impl OffRampType {
    /// `OffRampTypeEnum`'s names; `NONE` is no ramp.
    ///
    /// The outer `Option` is whether the name is known at all (`None`: an unknown name), the inner
    /// one whether it names a ramp (`Some(None)`: `NONE`).
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
    /// `new OffRampBase(filename)`: an empty `regions` array. Nothing is written until
    /// [`OffRamp::write`], where the Java opens the zip at once and fills it as it goes.
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

    /// The region's object both ramps start from:
    ///
    /// ```java
    /// regionObj.put("contig", loc.getContig());
    /// regionObj.put("start", loc.getStart());
    /// regionObj.put("end", loc.getEnd());
    /// ```
    fn region_object(location: &SimpleInterval) -> JsonObject {
        let mut region = JsonObject::new();
        region.put("contig", Json::Str(location.contig.clone()));
        region.put("start", Json::Int(location.start as i64));
        region.put("end", Json::Int(location.end as i64));
        region
    }

    /// `OffRampBase.addEntry`: the entry is named under the region's directory, or at the root when
    /// there is no region. The bytes are kept in order until the zip is written.
    fn add_entry(&mut self, location: Option<&SimpleInterval>, name: &str, bytes: Vec<u8>) {
        let prefix = location.map_or(String::new(), |l| format!("{}/", location_suffix(l)));
        self.entries.push((format!("{prefix}{name}"), bytes));
    }

    /// `OffRampBase.addHaplotypes`: the haplotypes as CSV, one row each.
    ///
    /// ```java
    /// pw.println(String.format("%s,%d,%d,%d,%s,%s,%s,%d",
    /// ```
    ///
    /// The score is `Double.toString(haplotype.getScore())`, hence `java_double_to_string`.
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
    ///
    /// The reads table is the `WIDE_READS = false` branch, `name,supp`. Both tables are written
    /// from the sample matrix's own evidence list, in its order. Unlike the Java, which builds the
    /// `msg` text for the wide branch on every sample, nothing is built for the branch not taken.
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
        // `value.alleles()`, one clone per allele index.
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

            // `read.isSupplementaryAlignment() ? 1 : 0`; 0x800 is the SAM supplementary flag.
            let supplementary = |read: &BamRecord| if read.flags & 0x800 != 0 { 1 } else { 0 };
            let mut names = String::from("name,supp\n");
            for read in reads {
                names.push_str(&format!("{},{}\n", read.read_name, supplementary(read)));
            }
            self.add_entry(Some(location), &format!("{base}.reads"), names.into_bytes());

            // `addSampleMatrix`: a header `read,supp,h0,h1,...`, then per read its name, its flag
            // and its likelihood against each haplotype, each `Double.toString`.
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
    ///
    /// ```java
    /// if ( assemblerObj.length() != 0 ) {
    ///     regionObj.put(name, assemblerObj);
    /// }
    /// ```
    ///
    /// `assembly` is the result set's padded reference span, its region for genotyping and its
    /// haplotypes, together or not at all, as the Java's `value != null`. `description` is the
    /// assembly region's `toString()`, which becomes the haplotype reads' `CR` tag.
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
    ///
    /// ```java
    /// record.setReadName(output.getHaplotypeSampleTag() + "_" + uniqueNameCounter++ + "_"
    ///     + haplotype_count);
    /// record.setAttribute(output.getHaplotypeSampleTag(), haplotype.hashCode());
    /// ```
    ///
    /// (The first statement is one Java line, `HaplotypeBAMWriter.java:235`, wrapped at the `+`.)
    ///
    /// The counter runs across the whole file and the index restarts at each region. The `HC` tag
    /// is `Haplotype.hashCode`, which is `Arrays.hashCode(getBases())`. The Java's `XC` tag, set
    /// only for a haplotype `isCollapsed()`, is never written here.
    fn write_haplotypes(
        &mut self,
        location: &SimpleInterval,
        description: &str,
        padded: &SimpleInterval,
        haplotypes: &[Haplotype],
        header: &SamHeader,
    ) {
        // The writer is opened on the first region that has an assembly (`haplotypeBAMWriter ==
        // null`), so the header is built from that region's source header and kept.
        let output_header = self
            .haplotype_header
            .get_or_insert_with(|| haplotype_bam_header(header))
            .clone();
        // `getSequenceIndex`, -1 for a contig the dictionary does not hold.
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
                // `Utils.dupBytes((byte) '!', length)` as the qualities: the byte 33, a phred of
                // 33.
                base_qualities: vec![b'!'; bases.len()],
                ..BamRecord::default()
            };
            self.name_counter += 1;
            // `as i64` widens the `int` hash exactly; a tag holds any integer width.
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
            // `EXT_SPECIAL_TAG`: reference flag, score, offset and the callable region's
            // coordinates, so the post-assembler on ramp can restore what a BAM record lacks.
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
            // The destination's header is coordinate-sorted, so the writer sorts the records
            // before they reach the file.
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

    /// `close()`: the zip, deflated, at `path`, its parent directories made first as the
    /// constructor's `getParentFile().mkdirs()` makes them.
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
///
/// ```java
/// final SAMReadGroupRecord rgRec = new SAMReadGroupRecord(haplotypeReadGroupID.get());
/// rgRec.setSample(haplotypeSampleTag);
/// rgRec.setSequencingCenter("BI");
/// ```
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

/// Which on ramp a file is: `PostFilterOnRamp` or `PostAssemblerOnRamp`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OnRampType {
    PostFilter,
    PostAssembler,
}

impl OnRampType {
    /// `OnRampTypeEnum`'s names; `NONE` is no ramp. The nesting is [`OffRampType::parse`]'s.
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
///
/// File order is enough here because the on ramp only looks keys up; it never writes the object
/// back out.
#[derive(Debug, Clone, PartialEq)]
pub enum JsonValue {
    Int(i64),
    Str(String),
    Object(Vec<(String, JsonValue)>),
    Array(Vec<JsonValue>),
}

impl JsonValue {
    /// `JSONObject.has(key)` then `get(key)`: the first member of that name, or nothing.
    pub fn get(&self, key: &str) -> Option<&JsonValue> {
        match self {
            JsonValue::Object(members) => members.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    /// `getString(key)`. The error text is the one this port reports for a missing or non-string
    /// member.
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
///
/// Written here rather than taken from a library for the reason the module header gives. It is
/// not a general parser: there are no floats, booleans or `null`, and a number is a run of `-`
/// and digits read as an `i64`.
fn parse_json(text: &str) -> Result<JsonValue, String> {
    // A cursor over the text's bytes. `'a` is a lifetime: it says the parser borrows `bytes` and
    // cannot outlive the text it was made from.
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
                            // `\uXXXX`: four hex digits, a code that is no character (a lone
                            // surrogate half) becoming U+FFFD.
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

/// `new SimpleInterval(String)` for the `contig:start-end` the ramps write: the LAST colon
/// splits the contig off, commas in the numbers are dropped. A contig name holding a colon or a
/// bare contig with no range is beyond what the ramps write and is refused here.
fn parse_interval(text: &str) -> Result<SimpleInterval, String> {
    let bad = || format!("not an interval: {text}");
    let (contig, range) = text.rsplit_once(':').ok_or_else(bad)?;
    let (start, end) = range.split_once('-').ok_or_else(bad)?;
    let number = |s: &str| s.replace(',', "").parse::<i32>().map_err(|_| bad());
    SimpleInterval::new(contig, number(start)?, number(end)?).ok_or_else(bad)
}

/// `TextCigarCodec.decode` (htsjdk 5.0.0), for the CIGARs `addHaplotypes` writes.
///
/// It is looser than the Java on text the ramps never write: the Java requires a digit before
/// every operator and throws on trailing digits, where this reads a bare operator as length zero
/// and drops trailing digits.
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
///
/// ```java
/// row = line.split(",");
/// ```
///
/// Java's `split` drops trailing empty fields and Rust's `split` keeps them, so a row ending in a
/// comma is one field longer here. The tables the off ramps write end in no empty field.
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

    /// `getColumn(name)`: the last column of that name, as the `LinkedHashMap` keeps it
    /// (`columns.put(toks[n], n)` overwrites an earlier position).
    fn column(&self, name: &str) -> Result<usize, String> {
        self.columns
            .iter()
            .rposition(|c| c == name)
            .ok_or_else(|| format!("java.io.IOException: no such column: {name}"))
    }
}

/// `Double.parseDouble`, for the `Double.toString` the ramps write. `parseDouble` ignores
/// surrounding whitespace, hence `trim`; Rust's parse also reads `NaN` and `Infinity`, the two
/// non-finite texts `Double.toString` produces.
fn parse_double(text: &str) -> Result<f64, String> {
    text.trim()
        .parse::<f64>()
        .map_err(|_| format!("java.lang.NumberFormatException: For input string: \"{text}\""))
}

/// `Integer.parseInt`, which, unlike `parseDouble`, does not trim.
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
    /// The Java map is a `LinkedHashMap`; it is only looked up, so a plain `HashMap` serves.
    region_index: std::collections::HashMap<String, usize>,
    /// The haplotype BAM's records in file order, with its contig names, for the post-assembler
    /// ramp.
    haplotype_records: Vec<BamRecord>,
    haplotype_contigs: Vec<String>,
}

impl OnRamp {
    /// `new OnRampBase(filename)`: the zip opened and `info.json` read; the post-assembler ramp
    /// also copies out its haplotype BAM and index.
    ///
    /// Every entry is read into memory at once, where the Java reads each on demand. The index is
    /// only checked for presence: the records are read in full and filtered by span in
    /// [`OnRamp::bam_haplotypes`] rather than queried through it.
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
        // `info.getJSONArray("regions")`: anything else there is the Java's `JSONException`.
        let Some(JsonValue::Array(regions)) = info.get("regions") else {
            return Err("JSONObject[\"regions\"] is not a JSONArray.".to_string());
        };
        // `regionKey`, `String.format("%s:%d-%d", ...)` of the object's own fields.
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
            // The Java copies the BAM out to a file and opens it with its index; here it is
            // decompressed and read in memory.
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
    ///
    /// ```java
    /// throw new IOException("no such: " + name);
    /// ```
    fn entry(&self, location: Option<&SimpleInterval>, name: &str) -> Result<&[u8], String> {
        let prefix = location.map_or(String::new(), |l| format!("{}/", location_suffix(l)));
        let full = format!("{prefix}{name}");
        self.entries
            .iter()
            .find(|(n, _)| *n == full)
            .map(|(_, bytes)| bytes.as_slice())
            .ok_or_else(|| format!("java.io.IOException: no such: {full}"))
    }

    /// `hasRegion(loc)`: whether `loc.toString()`, `contig:start-end`, is a key of the index.
    pub fn has_region(&self, location: &SimpleInterval) -> bool {
        self.region_index.contains_key(&interval_text(location))
    }

    /// `getRegion(loc, name)`: the region's object, then down the dotted path; nothing where a
    /// level is missing.
    ///
    /// ```java
    /// if ( !obj.has(key) )
    ///     return null;
    /// obj = obj.getJSONObject(key);
    /// ```
    ///
    /// The Java's `getJSONObject` throws when a step, the last one included, is not an object;
    /// here a step that is not an object answers nothing and the last step may be any value.
    pub fn region(&self, location: &SimpleInterval, name: &str) -> Option<&JsonValue> {
        let index = self.region_index.get(&interval_text(location))?;
        let mut object = self.regions.get(*index)?;
        for key in name.split('.') {
            object = object.get(key)?;
        }
        Some(object)
    }

    /// `getHaplotypes(loc, name)`: the haplotype table, one haplotype per row, its fields set in
    /// the Java's order (location, score, CIGAR, offset).
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
            // `row[column]`: a short row is the Java's `ArrayIndexOutOfBoundsException`.
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
    ///
    /// ```java
    /// if ( !specialAttrToks[3].equals(locStr) )
    ///     continue;;
    /// ```
    ///
    /// The Java's `query(contig, start, end, false)` is an overlap query; the condition below is
    /// the same overlap test over the records in file order, which a coordinate-sorted file
    /// returns in the same order. A record with no `XS` tag is kept with the defaults the Java
    /// starts from: not reference, no score, offset zero.
    pub fn bam_haplotypes(&self, location: &SimpleInterval) -> Result<Vec<Haplotype>, String> {
        let wanted = interval_text(location);
        let mut haplotypes = Vec::new();
        for record in &self.haplotype_records {
            // An unplaced record (reference index -1, which `try_from` refuses) is in no region.
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
            // `tokens.get(n)` is `None` past the end, where the Java's `specialAttrToks[n]` throws;
            // the tag this file writes always has four fields.
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
    ///
    /// ```java
    /// return new AssemblyRegion(
    ///         new SimpleInterval(obj.getString("activeSpan")),
    ///         new SimpleInterval(obj.getString("paddedSpan")),
    ///         true,
    ///         header);
    /// ```
    pub fn assembly_region(object: &JsonValue) -> Result<AssemblyRegion, String> {
        AssemblyRegion::new(
            parse_interval(object.string("activeSpan")?)?,
            parse_interval(object.string("paddedSpan")?)?,
            true,
        )
        .map_err(|e| format!("{e:?}"))
    }

    /// `getRegion(span, name).getString(key)` as an interval; a missing region or path is the
    /// Java's `NullPointerException` on the `null` that `getRegion` returns.
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
    ///
    /// `IGNORE_DANGLING_READS` is `true` in the reference, so both of its branches here are the
    /// skipping ones. A matrix row whose name is not the next expected read is passed over:
    ///
    /// ```java
    /// if ( !readName.equals(reads.get(row).getName()) ) {
    ///     if (IGNORE_DANGLING_READS) {
    /// ```
    ///
    /// The Java allocates one more column than there are reads and copies the last read's values
    /// into it (`matrix[i][row] = matrix[i][row - 1]`); this allocates exactly one per read.
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
        // `collectReadNames`: name plus supplementary flag to read, a later read replacing an
        // earlier one of the same key, as `LinkedHashMap.put` does.
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
                // `readsByName` first, then `readsByName2`; `or_else` consults the second map
                // only when the first has nothing.
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
            // `vec![x; n]` is n copies of x: a haplotypes-by-reads matrix of zeros.
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
///
/// Unlike [`crate::ramped_haplotype_caller::compare_haplotypes`], this one works on the engine's
/// `Haplotype` and returns an `Ordering`, so only the comparison's sign exists. Contigs compare
/// as `String.compareTo` does, by UTF-16 code unit (`compare_strings`).
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
        // `bool as i32` is 0 or 1, the Java's `isReference() ? 1 : 0`.
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
///
/// The qualities are compared as bytes where the Java compares
/// `new String(o1.getBaseQualitiesNoCopy())`; for bytes below 128, every Phred score a SAM can
/// carry, the two orders agree.
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
    // `GATKRead.commonToString`:
    // `String.format("%s %s:%d-%d", getName(), getContig(), getStart(), getEnd())`, or
    // `"<name> UNMAPPED"` for an unmapped read or one with an empty CIGAR.
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
