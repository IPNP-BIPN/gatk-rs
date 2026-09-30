//! `JointGermlineCNVSegmentation`: single-sample gCNV calls combined into one joint call set.
//!
//! The tool is a `MultiVariantWalkerGroupedOnStart` around two streaming cluster engines
//! ([`gatk_tools::sv_cluster_engine`]): a defragmenter that joins ONE sample's neighbouring
//! segments, and a max-clique engine that joins different samples' events. What it adds around
//! them is its own, and is all here:
//!
//! * **four entry filters drop a single-genotype record** before anything else: a hom-ref call, a
//!   no-call without CN, a QS strictly below the threshold, and a no-call whose CN is 0;
//! * **each genotype is padded to its ploidy** with reference alleles, the ploidy being the
//!   genotype's own `ECN`, the autosomal argument, or the pedigree's sex on an allosome; a single
//!   no-call allele becomes a no-call of the full ploidy;
//! * **more than one sample in the inputs skips the defragmenter**, the calls being assumed
//!   pre-clustered;
//! * **records are written as they leave the engines at each new contig**, sorted by start, and
//!   each group of overlapping records has its genotypes squared off: every sample gets a genotype,
//!   a sample inside an earlier overlapping event takes that event's copy number, and `AC`, `AF`
//!   and `AN` are recomputed. `AF` of a record with both `<DEL>` and `<DUP>` is the COUNT, not a
//!   frequency, as the reference writes it.
//!
//! The grouping arguments of the walker (`--combine-variants-distance`, `--max-distance`,
//! `--ref-padding`) only change how records are handed over, never their order, and the tool reads
//! no reference window, so none of them reaches the output.

use super::*;

use gatk_tools::sv_cluster_engine::{
    CanonicalLinkage, ClusteringType, CnvLinkage, Engine, EngineError, Padding,
};
use gatk_tools::sv_collapser::{self as collapser, genotype_int, Member};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

const TOOL: &str = "JointGermlineCNVSegmentation";
const ALLOSOMAL_CONTIGS: &[&str] = &["X", "Y", "chrX", "chrY"];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Sex {
    Male,
    Female,
    Unknown,
}

fn del_allele() -> Allele {
    Allele::from_str("<DEL>", false).expect("a symbolic allele")
}

fn dup_allele() -> Allele {
    Allele::from_str("<DUP>", false).expect("a symbolic allele")
}

fn engine_error(error: EngineError) -> Thrown {
    Thrown::non_user(error.class, error.message)
}

/// `PedReader`: one line per sample, the sex in the fifth column.
fn read_pedigree(path: &str) -> Result<Vec<(String, Sex)>, Thrown> {
    let text = std::fs::read_to_string(path)
        .map_err(|_| Thrown::user(format!("Couldn't read file {path}")))?;
    let mut samples = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(Thrown::user(format!(
                "Bad PED line: expected 6 fields but found {}",
                fields.len()
            )));
        }
        let sex = match fields[4] {
            "1" => Sex::Male,
            "2" => Sex::Female,
            _ => Sex::Unknown,
        };
        samples.push((fields[1].to_string(), sex));
    }
    Ok(samples)
}

fn has_extended(genotype: &Genotype, key: &str) -> bool {
    genotype.extended.iter().any(|(name, _)| name == key)
}

fn set_extended(genotype: &mut Genotype, key: &str, value: Value) {
    match genotype.extended.iter_mut().find(|(name, _)| name == key) {
        Some(slot) => slot.1 = value,
        None => genotype.extended.push((key.to_string(), value)),
    }
}

fn set_attribute(variant: &mut VariantContext, key: &str, value: Value) {
    match variant.attributes.iter_mut().find(|(name, _)| name == key) {
        Some(slot) => slot.1 = value,
        None => variant.attributes.push((key.to_string(), value)),
    }
}

/// `Genotype.isHomRef`.
fn is_hom_ref(genotype: &Genotype) -> bool {
    !genotype.alleles.is_empty()
        && genotype
            .alleles
            .iter()
            .all(|allele| !allele.is_no_call() && allele.is_reference())
}

/// `Genotype.isNoCall`: every allele uncalled, and at least one allele.
fn is_no_call(genotype: &Genotype) -> bool {
    !genotype.alleles.is_empty() && genotype.alleles.iter().all(Allele::is_no_call)
}

/// What the tool reads off its arguments once.
struct Settings {
    min_quality: i32,
    ref_autosomal_copy_number: i32,
    pedigree: Vec<(String, Sex)>,
    breakpoints: collapser::BreakpointSummary,
    alternates: collapser::AltAlleleSummary,
}

impl Settings {
    /// `getSamplePloidy`: the genotype's own `ECN` first, then the contig.
    fn ploidy(
        &self,
        sample: &str,
        contig: &str,
        genotype: Option<&Genotype>,
    ) -> Result<i32, Thrown> {
        if let Some(genotype) = genotype {
            if has_extended(genotype, "ECN") {
                return Ok(genotype_int(genotype, "ECN", 0));
            }
        }
        if !ALLOSOMAL_CONTIGS.contains(&contig) {
            return Ok(self.ref_autosomal_copy_number);
        }
        let sex = self
            .pedigree
            .iter()
            .find(|(name, _)| name == sample)
            .map(|(_, sex)| *sex);
        let Some(sex) = sex else {
            return match genotype {
                Some(genotype) => Ok(genotype.ploidy() as i32),
                None => Err(Thrown::non_user(
                    "java.lang.IllegalStateException",
                    format!("Sample {sample} is missing from the pedigree"),
                )),
            };
        };
        match contig {
            "X" | "chrX" => Ok(if sex == Sex::Female { 2 } else { 1 }),
            "Y" | "chrY" => Ok(if sex == Sex::Female { 0 } else { 1 }),
            other => Err(Thrown::non_user(
                "java.lang.IllegalArgumentException",
                format!(
                    "Encountered unknown allosomal contig: {other}. This tool only supports \
                     mammalian genomes with XX/XY sex determination."
                ),
            )),
        }
    }
}

/// `createDepthOnlyFromGCNVWithOriginalGenotypes`: `None` for a record the entry filters drop.
fn depth_only_record(
    record: &VariantContext,
    settings: &Settings,
    sequences: &[(String, i32)],
    wrapped: &dyn Fn(&VariantContext) -> Thrown,
) -> Result<Option<Member>, Thrown> {
    let genotypes: Vec<Genotype> = record.genotypes.iter().cloned().collect();
    if genotypes.len() == 1 {
        let genotype = &genotypes[0];
        let null_call = has_extended(genotype, "CN")
            && genotype_int(genotype, "CN", 0) == 0
            && is_no_call(genotype);
        if is_hom_ref(genotype)
            || (is_no_call(genotype) && !has_extended(genotype, "CN"))
            || genotype_int(genotype, "QS", 0) < settings.min_quality
            || null_call
        {
            return Ok(None);
        }
    }
    let mut variant = record.clone();
    set_attribute(
        &mut variant,
        gatk_tools::sv_call_record::ALGORITHMS_ATTRIBUTE,
        Value::List(vec![Value::Str(
            gatk_tools::sv_cluster::DEPTH_ALGORITHM.to_string(),
        )]),
    );
    let reference = record.reference().clone();
    let mut prepared = Vec::with_capacity(genotypes.len());
    for genotype in &genotypes {
        let ploidy = settings.ploidy(&genotype.sample_name, &record.contig, Some(genotype))?;
        let mut rebuilt = genotype.clone();
        if genotype.alleles.len() == 1 && genotype.alleles[0].is_no_call() {
            rebuilt.alleles = vec![Allele::no_call(); ploidy.max(0) as usize];
        } else {
            if genotype.alleles.len() as i32 > ploidy {
                return Err(wrapped(record));
            }
            while (rebuilt.alleles.len() as i32) < ploidy {
                rebuilt.alleles.push(reference.clone());
            }
        }
        set_extended(&mut rebuilt, "ECN", Value::Int(i64::from(ploidy)));
        prepared.push(rebuilt);
    }
    variant.genotypes = prepared.into();
    let call =
        gatk_tools::sv_call_record::create(&variant, sequences).map_err(|_| wrapped(record))?;
    let kept: Vec<Genotype> = variant
        .genotypes
        .iter()
        .filter(|g| !(is_hom_ref(g) || (is_no_call(g) && !has_extended(g, "CN"))))
        .cloned()
        .collect();
    Ok(Some(Member {
        call,
        alleles: variant.alleles.clone(),
        genotypes: kept,
        filters: variant.filters.clone().unwrap_or_default(),
    }))
}

/// `GATKVariantContextUtils.makePloidyLengthAlleleList`.
fn ploidy_length(ploidy: i32, allele: &Allele) -> Vec<Allele> {
    if ploidy == 0 {
        return vec![Allele::no_call()];
    }
    vec![allele.clone(); ploidy.max(0) as usize]
}

/// `GATKSVVariantContextUtils.makeGenotypeAllelesFromCopyNumber`.
fn alleles_from_copy_number(
    copy_number: i32,
    ref_copy_number: i32,
    reference: &Allele,
) -> Vec<Allele> {
    let allele = if copy_number > ref_copy_number {
        dup_allele()
    } else if copy_number < ref_copy_number {
        del_allele()
    } else {
        reference.clone()
    };
    if ref_copy_number == 0 {
        return vec![Allele::no_call()];
    }
    if ref_copy_number == 1 {
        return vec![allele];
    }
    if allele == dup_allele() {
        return vec![Allele::no_call(); ref_copy_number.max(0) as usize];
    }
    if ref_copy_number == 2 {
        let first = if copy_number == 0 {
            allele.clone()
        } else {
            reference.clone()
        };
        return vec![first, allele];
    }
    let mut out = Vec::new();
    for _ in 0..copy_number {
        out.push(reference.clone());
    }
    for _ in copy_number..ref_copy_number {
        out.push(del_allele());
    }
    out
}

/// `updateGenotypes`: every sample given a genotype of its ploidy, and the counts recomputed.
fn update_genotypes(
    settings: &Settings,
    samples: &[String],
    vc: &VariantContext,
    copy_numbers: &[(String, (i32, i32))],
) -> Result<VariantContext, Thrown> {
    let del = del_allele();
    let dup = dup_allele();
    let alternates: Vec<Allele> = vc.alternate_alleles().to_vec();
    if alternates.iter().any(|a| *a != del && *a != dup) {
        let listed: Vec<String> = alternates.iter().map(Allele::display_string).collect();
        return Err(Thrown::non_user(
            "java.lang.IllegalArgumentException",
            format!(
                "At site {}:{} variant context contains alternate alleles in addition to CNV <DEL> and <DUP> alleles: [{}]",
                vc.contig,
                vc.start,
                listed.join(", ")
            ),
        ));
    }
    let reference = vc.reference().clone();
    let (mut del_count, mut dup_count, mut allele_number) = (0i64, 0i64, 0i32);
    let mut genotypes = Vec::with_capacity(samples.len());
    for sample in samples {
        let genotype = vc.genotype(sample);
        let ploidy = settings.ploidy(sample, &vc.contig, genotype)?;
        allele_number += ploidy;
        let known = copy_numbers.iter().find(|(name, _)| name == sample);
        if known.is_none() && genotype.is_none() {
            let mut fresh = Genotype::new(sample, ploidy_length(ploidy, &reference));
            set_extended(&mut fresh, "CN", Value::Int(i64::from(ploidy)));
            genotypes.push(fresh);
            continue;
        }
        let (copy_number, alleles) = match (known, genotype) {
            (Some((_, (copy_number, end))), _) if i64::from(*end) > vc.start => {
                (*copy_number, ploidy_length(ploidy, &reference))
            }
            (_, Some(genotype)) => {
                let copy_number = genotype_int(genotype, "CN", ploidy);
                let alleles = if ploidy as usize == genotype.ploidy() {
                    genotype.alleles.clone()
                } else {
                    alleles_from_copy_number(copy_number, ploidy, &reference)
                };
                (copy_number, alleles)
            }
            _ => (ploidy, alleles_from_copy_number(ploidy, ploidy, &reference)),
        };
        let mut rebuilt = match genotype {
            Some(genotype) => genotype.clone(),
            None => Genotype::new(sample, Vec::new()),
        };
        set_extended(&mut rebuilt, "CN", Value::Int(i64::from(copy_number)));
        rebuilt.alleles = alleles.clone();
        genotypes.push(rebuilt);
        if genotype.is_some() {
            if alleles.contains(&del) {
                // `Allele::isNonReference` is true of a no-call too.
                del_count += alleles.iter().filter(|a| !a.is_reference()).count() as i64;
            } else if copy_number > ploidy {
                dup_count += 1;
            }
        }
    }
    let mut out = vc.clone();
    out.genotypes = genotypes.into();
    if allele_number > 0 {
        let count_of = |allele: &Allele| if *allele == del { del_count } else { dup_count };
        if alternates.len() == 1 {
            let count = count_of(&alternates[0]);
            set_attribute(&mut out, "AC", Value::Int(count));
            set_attribute(
                &mut out,
                "AF",
                Value::Double(count as f64 / f64::from(allele_number)),
            );
            set_attribute(&mut out, "AN", Value::Int(i64::from(allele_number)));
        } else {
            let counts: Vec<i64> = out.alleles[1..].iter().map(count_of).collect();
            set_attribute(
                &mut out,
                "AC",
                Value::List(counts.iter().map(|c| Value::Int(*c)).collect()),
            );
            set_attribute(
                &mut out,
                "AF",
                Value::List(counts.iter().map(|c| Value::Double(*c as f64)).collect()),
            );
            set_attribute(&mut out, "AN", Value::Int(i64::from(allele_number)));
        }
    }
    Ok(out)
}

/// `resolveVariantContexts` over one group of overlapping records.
fn resolve(
    settings: &Settings,
    samples: &[String],
    group: &[VariantContext],
) -> Result<Vec<VariantContext>, Thrown> {
    let mut copy_numbers: Vec<(String, (i32, i32))> = Vec::new();
    let mut resolved = Vec::with_capacity(group.len());
    for current in group {
        resolved.push(update_genotypes(settings, samples, current, &copy_numbers)?);
        let end = current
            .attributes
            .iter()
            .find(|(key, _)| key == "END")
            .and_then(|(_, value)| match value {
                Value::Int(end) => Some(*end as i32),
                other => other.format().and_then(|text| text.parse().ok()),
            })
            .unwrap_or(current.start as i32);
        for genotype in current.genotypes.iter() {
            if has_extended(genotype, "CN") {
                let copy_number = genotype_int(genotype, "CN", settings.ref_autosomal_copy_number);
                match copy_numbers
                    .iter_mut()
                    .find(|(name, _)| *name == genotype.sample_name)
                {
                    Some(slot) => slot.1 = (copy_number, end),
                    None => copy_numbers.push((genotype.sample_name.clone(), (copy_number, end))),
                }
            }
        }
    }
    Ok(resolved)
}

/// `buildAndSanitizeRecord`: the record as a variant, without its members, its algorithms or ECN.
fn sanitized(member: Member) -> VariantContext {
    let mut variant = gatk_tools::sv_call_record::to_variant(
        &member.call,
        member.alleles,
        member.genotypes,
        &member.filters,
    );
    variant.attributes.retain(|(key, _)| {
        key != collapser::CLUSTER_MEMBER_IDS_KEY
            && key != gatk_tools::sv_call_record::ALGORITHMS_ATTRIBUTE
    });
    let genotypes: Vec<Genotype> = variant
        .genotypes
        .iter()
        .cloned()
        .map(|mut genotype| {
            genotype.extended.retain(|(key, _)| key != "ECN");
            genotype
        })
        .collect();
    variant.genotypes = genotypes.into();
    variant
}

/// What the two engines and the output have in common while the records stream through.
struct Traversal<'a> {
    settings: &'a Settings,
    samples: &'a [String],
    reference: gatk_engine::reference::ReferenceFileSource,
    defragmenter: Engine<CnvLinkage>,
    clusterer: Engine<CanonicalLinkage>,
    defragmented: Vec<Member>,
    clustered: Vec<Member>,
    written: Vec<VariantContext>,
}

impl Traversal<'_> {
    fn collapse(
        &mut self,
        group: Vec<Member>,
        breakpoints: collapser::BreakpointSummary,
    ) -> Result<Member, Thrown> {
        let reference = &mut self.reference;
        let collapsed = collapser::collapse(
            &group,
            breakpoints,
            self.settings.alternates,
            collapser::FlagFieldLogic::Or,
            &mut |contig, position| {
                reference
                    .query(contig, position, position)
                    .ok()
                    .and_then(|bases| bases.first().copied())
            },
        )
        .map_err(|error| Thrown::non_user(error.class(), error.message()))?;
        Ok(Member {
            call: collapsed.call,
            alleles: collapsed.alleles,
            genotypes: collapsed.genotypes,
            filters: collapsed.filters,
        })
    }

    fn defragment(&mut self, groups: Vec<Vec<Member>>) -> Result<(), Thrown> {
        for group in groups {
            let member = self.collapse(group, collapser::BreakpointSummary::MinStartMaxEnd)?;
            self.defragmented.push(member);
        }
        Ok(())
    }

    fn cluster(&mut self, groups: Vec<Vec<Member>>) -> Result<(), Thrown> {
        for group in groups {
            let member = self.collapse(group, self.settings.breakpoints)?;
            self.clustered.push(member);
        }
        Ok(())
    }

    /// `processClusters`: the defragmenter flushed into the clusterer, the clusterer flushed, and
    /// what came out written in start order.
    fn process_clusters(&mut self) -> Result<(), Thrown> {
        let flushed = self.defragmenter.flush();
        self.defragment(flushed)?;
        let mut defragmented = std::mem::take(&mut self.defragmented);
        defragmented.sort_by_key(|member| member.call.position_a);
        for member in defragmented {
            let groups = self.clusterer.add_and_flush(member).map_err(engine_error)?;
            self.cluster(groups)?;
        }
        let flushed = self.clusterer.flush();
        self.cluster(flushed)?;
        let mut clustered = std::mem::take(&mut self.clustered);
        clustered.sort_by_key(|member| member.call.position_a);
        self.write(clustered)
    }

    /// `write`: overlapping records grouped, and each group's genotypes squared off.
    fn write(&mut self, calls: Vec<Member>) -> Result<(), Thrown> {
        let records: Vec<VariantContext> = calls.into_iter().map(sanitized).collect();
        let mut group: Vec<VariantContext> = Vec::new();
        let mut cluster_end: i64 = -1;
        let mut cluster_contig: Option<String> = None;
        for current in records {
            let joins = (cluster_end == -1 || current.start < cluster_end)
                && cluster_contig
                    .as_deref()
                    .is_none_or(|c| c == current.contig);
            if joins {
                if current.stop > cluster_end {
                    cluster_end = current.stop;
                }
                if cluster_contig.is_none() {
                    cluster_contig = Some(current.contig.clone());
                }
                group.push(current);
            } else {
                let resolved = resolve(self.settings, self.samples, &group)?;
                self.written.extend(resolved);
                cluster_end = current.stop;
                cluster_contig = Some(current.contig.clone());
                group = vec![current];
            }
        }
        if !group.is_empty() {
            let resolved = resolve(self.settings, self.samples, &group)?;
            self.written.extend(resolved);
        }
        Ok(())
    }
}

/// The model's call intervals: `featureFileToIntervals` over the file, intersected with the
/// traversal intervals by `mergeListsBySetOperator`.
fn model_call_intervals(
    path: &str,
    sequences: &[(String, i32)],
    traversal: Option<&[gatk_engine::interval::SimpleInterval]>,
) -> Result<Vec<(usize, i32, i32)>, Thrown> {
    let text = std::fs::read_to_string(path)
        .map_err(|_| Thrown::user(format!("Couldn't read file {path}")))?;
    let index_of = |contig: &str| -> Result<usize, Thrown> {
        sequences
            .iter()
            .position(|(name, _)| name == contig)
            .ok_or_else(|| {
                Thrown::user(format!(
                    "Contig {contig} not found in the sequence dictionary"
                ))
            })
    };
    let mut coverage = Vec::new();
    for line in text.lines() {
        if line.starts_with('@') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let start: i32 = fields.get(1).and_then(|f| f.parse().ok()).unwrap_or(0);
        let end: i32 = fields.get(2).and_then(|f| f.parse().ok()).unwrap_or(0);
        coverage.push((index_of(fields[0])?, start, end));
    }
    let traversal: Vec<(usize, i32, i32)> = match traversal {
        Some(intervals) => intervals
            .iter()
            .map(|interval| Ok((index_of(&interval.contig)?, interval.start, interval.end)))
            .collect::<Result<_, Thrown>>()?,
        None => sequences
            .iter()
            .enumerate()
            .map(|(index, (_, length))| (index, 1, *length))
            .collect(),
    };
    if coverage.is_empty() || traversal.is_empty() {
        return Ok(if coverage.is_empty() {
            traversal
        } else {
            coverage
        });
    }
    let is_before =
        |a: &(usize, i32, i32), b: &(usize, i32, i32)| a.0 < b.0 || (a.0 == b.0 && a.2 < b.1);
    let mut out = Vec::new();
    let (mut one, mut two) = (0, 0);
    while two < traversal.len() && one < coverage.len() {
        if is_before(&traversal[two], &coverage[one]) {
            two += 1;
        } else if is_before(&coverage[one], &traversal[two]) {
            one += 1;
        } else {
            let (a, b) = (coverage[one], traversal[two]);
            out.push((a.0, a.1.max(b.1), a.2.min(b.2)));
            if a.2 < b.2 {
                one += 1;
            } else {
                two += 1;
            }
        }
    }
    if out.is_empty() {
        return Err(Thrown::user("There was an empty intersection"));
    }
    Ok(out)
}

fn parameter<T: std::str::FromStr>(parser: &Parser, name: &str, default: T) -> T {
    scalar(parser, name)
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

/// `JointGermlineCNVSegmentation`.
pub fn joint_germline_cnv_segmentation(parser: &Parser) -> Outcome {
    use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType};

    let _ = resolve_read_filters(parser, TOOL)?;
    let inputs = arguments(parser, "variant");
    if inputs.len() > 1 {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "More than one --variant is a GATK feature that this port does not carry yet. This message is the port's own and not GATK's.",
        ));
    }
    let input = inputs.into_iter().next().ok_or_else(|| {
        Thrown::command_line("Argument variant was missing: Argument 'variant' is required")
    })?;
    let VariantWalkerStart {
        input,
        text,
        intervals,
        ..
    } = variant_walker_startup_over(parser, input)?;
    let output = argument(parser, "output").ok_or_else(|| {
        Thrown::command_line("Argument output was missing: Argument 'output' is required")
    })?;
    let pedigree_path = argument(parser, "pedigree").ok_or_else(|| {
        Thrown::command_line("Argument pedigree was missing: Argument 'pedigree' is required")
    })?;
    let Some(reference_path) = argument(parser, "reference") else {
        return Err(Thrown::non_user(
            PORT_LIMITATION,
            "A walker without --reference is refused by the engine before it starts, which this port does not word yet. This message is the port's own and not GATK's.",
        ));
    };
    let Some(dictionary) = reference_dictionary(parser)? else {
        return Err(Thrown::user("Reference sequence dictionary required"));
    };
    let reference =
        gatk_engine::reference::ReferenceFileSource::open(std::path::Path::new(&reference_path))
            .map_err(|error| Thrown::user(format!("{error:?}")))?;
    let sequences: Vec<(String, i32)> = dictionary
        .sequences
        .iter()
        .map(|sequence| (sequence.name.clone(), sequence.length))
        .collect();
    let file = htsjdk_vcf::reader::read_vcf(&text)
        .map_err(|failure| Thrown::user(format!("{:?}", failure.error)))?;
    let mut samples = file.header.samples.clone();
    samples.sort();
    samples.dedup();

    // `onTraversalStart`: the pedigree, validated strictly, then the model's intervals.
    let pedigree = read_pedigree(&pedigree_path)?;
    for sample in &samples {
        if !pedigree.iter().any(|(name, _)| name == sample) {
            return Err(Thrown::user(format!(
                "Sample {sample} found in data sources but not in pedigree files with STRICT pedigree validation"
            )));
        }
    }
    let padding_fraction = parameter(parser, "defragmentation-padding-fraction", 0.25);
    let min_sample_overlap = parameter(parser, "min-sample-set-fraction-overlap", 0.0);
    let padding = match argument(parser, "model-call-intervals") {
        Some(path) => Padding::Bins(model_call_intervals(
            &path,
            &sequences,
            intervals.as_deref(),
        )?),
        None => Padding::Fraction,
    };
    let settings = Settings {
        min_quality: parameter(parser, "minimum-qs-score", 20),
        ref_autosomal_copy_number: parameter(parser, "autosomal-ref-copy-number", 2),
        pedigree,
        breakpoints: scalar(parser, "breakpoint-summary-strategy")
            .and_then(|name| collapser::BreakpointSummary::value_of(&name))
            .unwrap_or(collapser::BreakpointSummary::MedianStartMedianEnd),
        alternates: match scalar(parser, "alt-allele-summary-strategy").as_deref() {
            Some("MOST_SPECIFIC_SUBTYPE") => collapser::AltAlleleSummary::MostSpecificSubtype,
            _ => collapser::AltAlleleSummary::CommonSubtype,
        },
    };
    let defragmenter = Engine::new(
        ClusteringType::SingleLinkage,
        CnvLinkage {
            dictionary: sequences.clone(),
            padding_fraction,
            min_sample_overlap,
            padding,
        },
    );
    let clusterer = Engine::new(
        ClusteringType::MaxClique,
        CanonicalLinkage {
            dictionary: sequences.clone(),
            linkage: gatk_tools::sv_cluster::Linkage {
                depth: gatk_tools::sv_cluster::ClusteringParameters::depth(
                    parameter(parser, "clustering-interval-overlap", 0.8),
                    parameter(parser, "clustering-size-similarity", 0.0),
                    parameter(parser, "clustering-breakend-window", 10_000_000),
                    0.0,
                ),
                mixed: gatk_tools::sv_cluster::default_mixed_parameters(),
                pesr: gatk_tools::sv_cluster::default_pesr_parameters(),
                cluster_del_with_dup: true,
            },
        },
    );

    // `getVCFWriter`: the input's lines, the tool's own, and the samples sorted.
    let mut header = file.header.clone();
    header.samples = samples.clone();
    let compound =
        |id: &str, number: Cardinality, line_type: LineType, text: &str| HeaderLine::Compound {
            key: "INFO".to_string(),
            id: id.to_string(),
            number,
            line_type,
            description: text.to_string(),
            extra: Vec::new(),
        };
    let mut added = default_tool_vcf_header_lines(parser, TOOL);
    added.extend([
        compound(
            "SVLEN",
            Cardinality::Unbounded,
            LineType::Integer,
            "Difference in length between REF and ALT alleles",
        ),
        compound(
            "SVTYPE",
            Cardinality::Fixed(1),
            LineType::String,
            "Type of structural variant",
        ),
        compound(
            "AF",
            Cardinality::A,
            LineType::Float,
            "Allele Frequency, for each ALT allele, in the same order as listed",
        ),
        compound(
            "AC",
            Cardinality::A,
            LineType::Integer,
            "Allele count in genotypes, for each ALT allele, in the same order as listed",
        ),
        compound(
            "AN",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Total number of alleles in called genotypes",
        ),
    ]);
    for line in added {
        if !header.lines.iter().any(|existing| existing == &line) {
            header.lines.push(line);
        }
    }

    let multi_sample = samples.len() != 1;
    let records = variants_in_traversal(&file.records, intervals.as_deref(), &input)?;
    let only_starting_inside = flag(parser, "ignore-variants-starting-outside-interval");
    let mut traversal = Traversal {
        settings: &settings,
        samples: &samples,
        reference,
        defragmenter,
        clusterer,
        defragmented: Vec::new(),
        clustered: Vec::new(),
        written: Vec::new(),
    };
    let wrapped = |record: &VariantContext| {
        Thrown::non_user(
            "org.broadinstitute.hellbender.exceptions.GATKException",
            format!(
                "Exception thrown at {}:{} {}",
                record.contig,
                record.start,
                java_variant_context_string(record, &input)
            ),
        )
    };
    let result = (|| -> Result<(), Thrown> {
        let mut current: Option<String> = None;
        for record in records {
            if only_starting_inside {
                if let Some(intervals) = intervals.as_deref() {
                    let start = record.start as i32;
                    let inside = intervals.iter().any(|interval| {
                        interval.contig == record.contig
                            && interval.start <= start
                            && start <= interval.end
                    });
                    if !inside {
                        continue;
                    }
                }
            }
            if current.as_deref().is_some_and(|c| c != record.contig) {
                traversal.process_clusters()?;
            }
            current = Some(record.contig.clone());
            let Some(member) = depth_only_record(record, &settings, &sequences, &wrapped)? else {
                continue;
            };
            if multi_sample {
                let groups = traversal
                    .clusterer
                    .add_and_flush(member)
                    .map_err(engine_error)?;
                traversal.cluster(groups)?;
            } else {
                let groups = traversal
                    .defragmenter
                    .add_and_flush(member)
                    .map_err(engine_error)?;
                traversal.defragment(groups)?;
            }
        }
        traversal.process_clusters()
    })();

    // `closeTool`: whatever was written before a refusal stays written.
    let mut written = std::mem::take(&mut traversal.written);
    apply_sites_only(parser, &mut header, &mut written);
    let out = write_vcf_honouring_lenient(parser, &header, &written)?;
    write_variant_output(parser, &output, &out)?;
    result?;
    // `onTraversalSuccess` returns null, so `handleResult` prints nothing.
    Ok(None)
}
