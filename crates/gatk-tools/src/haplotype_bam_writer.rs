//! `--bam-output`, ported from `utils.haplotype.HaplotypeBAMWriter`, `HaplotypeBAMDestination`
//! and `SAMFileDestination`, with the read tags `AssemblyBasedCallerUtils` sets for it
//! (GATK 4.6.2.0).
//!
//! Each called region adds its haplotypes as reads (named `HC_<n>_<index>`, base qualities of
//! 33, mapping quality 60 for the best and 0 for the others, reverse strand, an `HC` tag holding
//! the haplotype's hash, the artificial read group and the callable region in `CR`) and then its
//! reads as they were scored: realigned to their best haplotype, their original soft-clip tags
//! dropped, `AR` naming the location of their best haplotype, `CR` the callable region and `XA`
//! the alleles they were informative for. The writer is not presorted, so htsjdk sorts the lot by
//! coordinate when it closes, stably.
//!
//! # Orders
//!
//! With `CALLED_HAPLOTYPES` the haplotypes are written in the order of a `HashSet` copied from the
//! genotyping engine's own, which is their `hashCode` (the bases') laid out by that copy's capacity.
//! The order reaches the read names only, since the output is sorted.

use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::java_hash::{JavaHashCode, JavaHashMap};
use htsjdk_bam::header::{ProgramRecord, ReadGroup, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// `DEFAULT_HAPLOTYPE_READ_GROUP_ID`.
pub const HAPLOTYPE_READ_GROUP_ID: &str = "ArtificialHaplotypeRG";
/// `haplotypeSampleTag`, which is the sample, the tag and the read-name prefix.
const HAPLOTYPE_SAMPLE_TAG: &str = "HC";
const BEST_HAPLOTYPE_MQ: u8 = 60;
const OTHER_MQ: u8 = 0;
/// `ReferenceConfidenceModel.ORIGINAL_SOFTCLIP_START_TAG` and `_END_TAG`.
const ORIGINAL_SOFTCLIP_TAGS: [&[u8; 2]; 2] = [b"OS", b"OE"];

/// `HaplotypeBAMWriter.WriterType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WriterType {
    AllPossibleHaplotypes,
    CalledHaplotypes,
    NoHaplotypes,
    CalledHaplotypesNoReads,
}

impl WriterType {
    pub fn parse(name: &str) -> Option<WriterType> {
        Some(match name {
            "ALL_POSSIBLE_HAPLOTYPES" => WriterType::AllPossibleHaplotypes,
            "CALLED_HAPLOTYPES" => WriterType::CalledHaplotypes,
            "NO_HAPLOTYPES" => WriterType::NoHaplotypes,
            "CALLED_HAPLOTYPES_NO_READS" => WriterType::CalledHaplotypesNoReads,
            _ => return None,
        })
    }
}

/// A haplotype as a `HashSet` key: its `hashCode` and `equals`.
#[derive(Clone, PartialEq)]
struct HaplotypeKey(Haplotype);

impl JavaHashCode for HaplotypeKey {
    fn java_hash_code(&self) -> i32 {
        self.0.java_hash_code()
    }
}

/// The writer, holding what it was given until it closes.
#[derive(Debug, Clone)]
pub struct HaplotypeBamWriter {
    kind: WriterType,
    header: SamHeader,
    records: Vec<BamRecord>,
    counter: i64,
}

impl HaplotypeBamWriter {
    /// `new HaplotypeBAMWriter(type, path, ...)`: the source's dictionary, coordinate order, its
    /// read groups and the artificial one (unless `NO_HAPLOTYPES`), its programs and one more.
    pub fn new(kind: WriterType, source: &SamHeader) -> HaplotypeBamWriter {
        let mut header = SamHeader::default();
        header.sequences = source.sequences.clone();
        header.attributes.set("SO", "coordinate");
        header.read_groups = source.read_groups.clone();
        if kind != WriterType::NoHaplotypes {
            let mut group = ReadGroup::new(HAPLOTYPE_READ_GROUP_ID);
            group.attributes.set("SM", HAPLOTYPE_SAMPLE_TAG);
            group.attributes.set("CN", "BI");
            header.read_groups.push(group);
        }
        header.programs = source.programs.clone();
        header
            .programs
            .push(ProgramRecord::new("HaplotypeBAMWriter"));
        HaplotypeBamWriter {
            kind,
            header,
            records: Vec::new(),
            counter: 1,
        }
    }

    /// `writeReadsAlignedToHaplotypes(haplotypes, paddedReferenceLoc, haplotypes,
    /// calledHaplotypes, readLikelihoods, callableRegion)`: `called` are the genotyping engine's
    /// called haplotypes in the order its set was filled, `extra_called` one added to the copy (the
    /// reference, under `--disable-optimizations`), `supported` each read's `XA`.
    pub fn write_region(
        &mut self,
        haplotypes: &[Haplotype],
        padded_reference_loc: &SimpleInterval,
        called: &[Haplotype],
        extra_called: Option<&Haplotype>,
        likelihoods: &AlleleLikelihoods<BamRecord, Haplotype>,
        callable: &SimpleInterval,
        supported: &dyn Fn(&BamRecord) -> Option<String>,
    ) -> Result<(), String> {
        match self.kind {
            WriterType::CalledHaplotypes | WriterType::CalledHaplotypesNoReads => {
                if called.is_empty() && extra_called.is_none() {
                    return Ok(());
                }
                // `new HashSet<>(calledHaplotypes.getCalledHaplotypes())`: the engine's set,
                // then its copy.
                let mut engine: JavaHashMap<HaplotypeKey, ()> = JavaHashMap::new();
                for haplotype in called {
                    engine.insert(HaplotypeKey(haplotype.clone()), ());
                }
                let mut copy: JavaHashMap<HaplotypeKey, ()> = JavaHashMap::with_capacity(
                    JavaHashMap::<HaplotypeKey, ()>::copy_capacity(engine.len()),
                );
                for key in engine.keys() {
                    copy.insert(key.clone(), ());
                }
                // `calledHaplotypeSet.add(refHaplotype)` under `--disable-optimizations`.
                if let Some(reference) = extra_called {
                    copy.insert(HaplotypeKey(reference.clone()), ());
                }
                let ordered: Vec<Haplotype> = copy.keys().map(|k| k.0.clone()).collect();
                for (index, haplotype) in ordered.iter().enumerate() {
                    self.write_haplotype(haplotype, padded_reference_loc, true, callable, index)?;
                }
            }
            WriterType::AllPossibleHaplotypes => {
                for (index, haplotype) in haplotypes.iter().enumerate() {
                    self.write_haplotype(haplotype, padded_reference_loc, true, callable, index)?;
                }
            }
            WriterType::NoHaplotypes => {}
        }
        if self.kind == WriterType::CalledHaplotypesNoReads {
            return Ok(());
        }
        // `annotateReadLikelihoodsWithRegions`: `AR` from each read's best haplotype, the
        // reference and the simplest CIGAR winning ties, then `CR`.
        let priorities: Vec<f64> = (0..likelihoods.number_of_alleles())
            .map(|a| {
                let h = likelihoods.get_allele(a).expect("an allele");
                let cigar_term = h.cigar().map_or(0, |c| 1 - c.elements.len() as i64);
                f64::from(i32::from(h.is_reference())) + cigar_term as f64
            })
            .collect();
        for sample in 0..likelihoods.number_of_samples() {
            let evidence = likelihoods.sample_evidence(sample).unwrap_or(&[]);
            let best = likelihoods.best_alleles_breaking_ties_for_sample(sample, Some(&priorities));
            for (index, read) in evidence.iter().enumerate() {
                let mut read = read.clone();
                let region = best
                    .iter()
                    .find(|b| b.evidence_index == index)
                    .and_then(|b| b.allele.as_ref())
                    .and_then(|h| h.genome_location())
                    .map(|l| format!("{}:{}-{}", l.contig, l.start, l.end));
                if let Some(region) = region {
                    set_tag(&mut read, b"AR", TagValue::Str(region));
                }
                set_tag(&mut read, b"CR", TagValue::Str(render(callable)));
                if let Some(xa) = supported(&read) {
                    set_tag(&mut read, b"XA", TagValue::Str(xa));
                }
                for tag in ORIGINAL_SOFTCLIP_TAGS {
                    read.tags.remove(Tag::new(tag));
                }
                self.records.push(read);
            }
        }
        Ok(())
    }

    /// `writeHaplotype`.
    fn write_haplotype(
        &mut self,
        haplotype: &Haplotype,
        padded_reference_loc: &SimpleInterval,
        best: bool,
        callable: &SimpleInterval,
        index: usize,
    ) -> Result<(), String> {
        let bases = haplotype.bases();
        let reference_index = self
            .header
            .sequences
            .iter()
            .position(|s| s.name == padded_reference_loc.contig)
            .map_or(-1, |i| i as i32);
        let mut record = BamRecord::default();
        record.read_name = format!("{HAPLOTYPE_SAMPLE_TAG}_{}_{index}", self.counter);
        self.counter += 1;
        record.base_qualities = vec![b'!'; bases.len()];
        record.read_bases = bases;
        record.alignment_start =
            padded_reference_loc.start + haplotype.alignment_start_hap_wrt_ref();
        record.cigar = haplotype
            .cigar()
            .cloned()
            .ok_or_else(|| "java.lang.NullPointerException".to_string())?;
        record.mapping_quality = if best { BEST_HAPLOTYPE_MQ } else { OTHER_MQ };
        record.reference_index = reference_index;
        record.mate_reference_index = -1;
        record.flags = 0x10;
        set_tag(
            &mut record,
            HAPLOTYPE_SAMPLE_TAG
                .as_bytes()
                .try_into()
                .expect("two bytes"),
            TagValue::Int(i64::from(haplotype.java_hash_code())),
        );
        set_tag(
            &mut record,
            b"RG",
            TagValue::Str(HAPLOTYPE_READ_GROUP_ID.to_string()),
        );
        set_tag(&mut record, b"CR", TagValue::Str(render(callable)));
        if haplotype.is_collapsed() {
            set_tag(&mut record, b"XC", TagValue::Str("1".to_string()));
        }
        self.records.push(record);
        Ok(())
    }

    /// `close()`: the records sorted by coordinate, stably, and written.
    pub fn finish(
        mut self,
        create_index: bool,
        level: u32,
        deflater: htsjdk_bgzf::Deflater,
    ) -> Result<(Vec<u8>, Option<Vec<u8>>), String> {
        self.records.sort_by(htsjdk_bam::coordinate::compare);
        crate::sam_output::write_records_with(
            &self.header,
            &self.records,
            create_index,
            level,
            deflater,
        )
        .map_err(|e| format!("{e:?}"))
    }
}

/// `SimpleInterval.toString()`.
fn render(interval: &SimpleInterval) -> String {
    format!("{}:{}-{}", interval.contig, interval.start, interval.end)
}

/// `setAttribute`: replaced where the tag is, appended where it is not.
fn set_tag(record: &mut BamRecord, tag: &[u8; 2], value: TagValue) {
    record.tags.insert(Tag::new(tag), value);
}
