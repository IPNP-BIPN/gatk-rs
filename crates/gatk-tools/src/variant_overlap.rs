//! `VariantOverlapAnnotator` with `--dbsnp` as its one input: the rsID and the `DB` flag.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.annotator.VariantOverlapAnnotator`
//! (GATK 4.6.2.0). The overlap map holds `--dbsnp` under `DB` and nothing else, so
//! `annotateOverlaps` asks the same question `annotateRsID` does and sets the flag exactly when an
//! rsID was found.
//!
//! # What "the same variant" means
//!
//! Both records are split into `Event`s with `splitVariantContextToEvents(vc, true, ...)`, and a
//! source matches when any of its events equals any of the annotated record's. An `Event` compares
//! its start, reference and alternate, and NOT its contig. The split keeps a biallelic record as
//! it stands, gives a record with no alternate no event at all, and trims every pair it cuts out
//! of a multi-allelic record at both ends. The `Event` then cuts the suffix its two alleles
//! share, so an untrimmed biallelic `GGCT>GT` still matches the `GGC>G` a split trims to.

use gatk_engine::variant_context_utils::{
    trim_alleles, Allele as EngineAllele, Variant as EngineVariant,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Value, VariantContext};

/// `Event`: the start and the two alleles, in their minimal representation.
#[derive(Debug, Clone, PartialEq)]
struct Event {
    start: i64,
    reference: Allele,
    alternate: Allele,
}

/// `new Event(contig, start, ref, alt)`, which cuts the two alleles to their minimal
/// representation: the suffix they share, unless either is one base long or their last bases
/// differ. A pair with no base left on one side is `Allele.create` of nothing, which throws.
fn event(start: i64, reference: &Allele, alternate: &Allele) -> Result<Event, String> {
    let (ours, theirs) = (reference.display_string(), alternate.display_string());
    let (ours, theirs) = (ours.as_bytes(), theirs.as_bytes());
    if reference.len() == 1 || alternate.len() == 1 || ours.last() != theirs.last() {
        return Ok(Event {
            start,
            reference: reference.clone(),
            alternate: alternate.clone(),
        });
    }
    if ours == theirs {
        return Err("ref and alt alleles are identical".to_string());
    }
    let shared = ours
        .iter()
        .rev()
        .zip(theirs.iter().rev())
        .take_while(|(a, b)| a == b)
        .count();
    let create = |bases: &[u8], is_ref: bool| {
        Allele::create(bases, is_ref).map_err(|_| "Null alleles are not supported".to_string())
    };
    Ok(Event {
        start,
        reference: create(&ours[..ours.len() - shared], true)?,
        alternate: create(&theirs[..theirs.len() - shared], false)?,
    })
}

/// `splitVariantContextToEvents(vc, true, SET_TO_NO_CALL_NO_ANNOTATIONS, true)`, as far as an
/// event reads it: genotypes and attributes never reach the comparison.
fn events(vc: &VariantContext) -> Result<Vec<Event>, String> {
    if vc.alleles.len() <= 1 {
        return Ok(Vec::new());
    }
    if vc.alleles.len() == 2 {
        return Ok(vec![event(vc.start, &vc.alleles[0], &vc.alleles[1])?]);
    }
    vc.alleles[1..]
        .iter()
        .map(|alternate| {
            let pair = EngineVariant {
                contig: vc.contig.clone(),
                start: vc.start as i32,
                stop: vc.stop as i32,
                alleles: vec![
                    EngineAllele::new(vc.alleles[0].display_string().as_bytes(), true),
                    EngineAllele::new(alternate.display_string().as_bytes(), false),
                ],
                genotypes: Vec::new(),
                attributes: Vec::new(),
            };
            let trimmed = trim_alleles(&pair, true, true).map_err(|error| format!("{error:?}"))?;
            // A symbolic allele or a spanning deletion comes back as it went in.
            let cut = |index: usize, original: &Allele| {
                if original.is_symbolic()
                    || trimmed.alleles[index].bases == pair.alleles[index].bases
                {
                    Ok(original.clone())
                } else {
                    Allele::create(&trimmed.alleles[index].bases, index == 0)
                        .map_err(|_| "Null alleles are not supported".to_string())
                }
            };
            event(
                i64::from(trimmed.start),
                &cut(0, &vc.alleles[0])?,
                &cut(1, alternate)?,
            )
        })
        .collect()
}

/// `getRsID(rsIDSourceVCs, vcToAnnotate)`: the IDs of every unfiltered source sharing an event
/// with the record, joined with `;`, in the sources' order. A source with no ID contributes `.`.
pub fn rs_id(sources: &[&VariantContext], vc: &VariantContext) -> Result<Option<String>, String> {
    let annotated = events(vc)?;
    let mut ids: Vec<String> = Vec::new();
    for source in sources {
        if source.is_filtered() {
            continue;
        }
        if source.contig != vc.contig {
            return Err(format!(
                "source rsID VariantContext {}:{} is not on same chromosome as vcToAnnotate {}:{}",
                source.contig, source.start, vc.contig, vc.start
            ));
        }
        if events(source)?
            .iter()
            .any(|event| annotated.contains(event))
        {
            ids.push(source.id.clone());
        }
    }
    Ok((!ids.is_empty()).then(|| ids.join(";")))
}

/// `annotateOverlaps(features, annotateRsID(features, vc))` with `--dbsnp` alone: the rsID into
/// the ID column, appended with `;` unless the column already contains it, and `DB` set.
///
/// `sources` are the `--dbsnp` records starting where `vc` starts, `getValues(dbSNP, start)`.
pub fn annotate(
    sources: &[&VariantContext],
    vc: &VariantContext,
) -> Result<VariantContext, String> {
    let Some(id) = rs_id(sources, vc)? else {
        return Ok(vc.clone());
    };
    let mut out = vc.clone();
    if vc.id == "." {
        out.id = id;
    } else if !vc.id.contains(&id) {
        out.id = format!("{};{id}", vc.id);
    }
    match out.attributes.iter_mut().find(|(key, _)| key == "DB") {
        Some((_, value)) => *value = Value::Bool(true),
        None => out.attributes.push(("DB".to_string(), Value::Bool(true))),
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(start: i64, alleles: &[(&str, bool)], id: &str) -> VariantContext {
        let mut vc = VariantContext::new(
            "20",
            start,
            alleles
                .iter()
                .map(|(bases, is_ref)| Allele::from_str(bases, *is_ref).unwrap())
                .collect(),
        );
        vc.id = id.to_string();
        vc
    }

    #[test]
    fn a_matching_snp_takes_the_id_and_the_flag() {
        let source = record(100, &[("A", true), ("G", false)], "rs1");
        let vc = record(100, &[("A", true), ("G", false)], ".");
        let out = annotate(&[&source], &vc).unwrap();
        assert_eq!(out.id, "rs1");
        assert!(out
            .attributes
            .contains(&("DB".to_string(), Value::Bool(true))));
    }

    #[test]
    fn a_split_pair_and_an_untrimmed_source_meet_in_the_minimal_representation() {
        // GGCT -> GT,G: the first pair trims to GGC>G, and the untrimmed biallelic GGCT>GT is
        // cut to the same by the event. The second pair keeps GGCT>G, one base long.
        let vc = record(100, &[("GGCT", true), ("GT", false), ("G", false)], ".");
        let trimmed = record(100, &[("GGC", true), ("G", false)], "rs1");
        let untrimmed = record(100, &[("GGCT", true), ("GT", false)], "rs2");
        let shorter = record(100, &[("GGC", true), ("TC", false)], "rs3");
        assert_eq!(
            rs_id(&[&trimmed, &untrimmed, &shorter], &vc)
                .unwrap()
                .as_deref(),
            Some("rs1;rs2")
        );
    }

    #[test]
    fn a_suffix_that_eats_an_allele_throws() {
        // GGC>GC shares GC, which leaves the alternate with no base: `Allele.create` refuses it.
        let source = record(100, &[("GGC", true), ("GC", false)], "rs1");
        let vc = record(100, &[("A", true), ("G", false)], ".");
        assert_eq!(
            rs_id(&[&source], &vc),
            Err("Null alleles are not supported".to_string())
        );
    }

    #[test]
    fn an_existing_id_is_appended_to_unless_it_holds_the_rsid() {
        let source = record(100, &[("A", true), ("G", false)], "rs1");
        let other = record(100, &[("A", true), ("G", false)], "rs9");
        assert_eq!(annotate(&[&source], &other).unwrap().id, "rs9;rs1");
        let holding = record(100, &[("A", true), ("G", false)], "rs1");
        assert_eq!(annotate(&[&source], &holding).unwrap().id, "rs1");
    }

    #[test]
    fn a_filtered_source_or_a_different_alternate_is_no_match() {
        let mut filtered = record(100, &[("A", true), ("G", false)], "rs1");
        filtered.filters = Some(vec!["q10".to_string()]);
        let other = record(100, &[("A", true), ("T", false)], "rs2");
        let vc = record(100, &[("A", true), ("G", false)], ".");
        let out = annotate(&[&filtered, &other], &vc).unwrap();
        assert_eq!(out.id, ".");
        assert!(out.attributes.is_empty());
    }
}
