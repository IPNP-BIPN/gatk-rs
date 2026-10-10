//! Conformance for the haplotype-to-allele mapping against the oracle.
//!
//! Golden from `tools/readfilter-conformance/AlleleMappingDump.java`. Each case gives the
//! haplotypes and the likelihood rows; this test builds the event maps, walks every event start in
//! the window as `HaplotypeCallerGenotypingEngine` does, and renders the `starts`, `event`,
//! `merged`, `mapper`, `marginal` and `error` rows, which must be the golden's in order.

use gatk_corpus as corpus;
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::allele_mapping::{
    create_allele_mapper, make_merged_variant_context, replace_span_dels,
    variants_from_active_haplotypes,
};
use gatk_engine::event_map::{build_event_maps_for_haplotypes, event_start_positions};
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;

const REF: &str = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGTAAAAAGTC";
const START: i32 = 1000;

/// `Allele.toString()`: the bases, and `*` after a reference allele.
fn text(a: &Allele) -> String {
    let mut s = a.display_string();
    if a.is_reference() {
        s.push('*');
    }
    s
}

fn run(case: &[&str], haps: &[Vec<&str>], lk: &[Vec<&str>]) -> Vec<String> {
    let label = case[1];
    let max_mnp: i32 = case[2].parse().unwrap();
    let spanning = case[3] == "true";
    let (window_start, window_end): (i32, i32) =
        (case[4].parse().unwrap(), case[5].parse().unwrap());
    let ref_loc = SimpleInterval::new("chr1", START, START + REF.len() as i32 - 1).unwrap();
    let mut haplotypes: Vec<Haplotype> = haps
        .iter()
        .map(|f| {
            let mut h = Haplotype::new(f[4].as_bytes(), f[3] == "true").unwrap();
            h.set_cigar(&htsjdk_bam::text_parse::parse_cigar(f[5]).unwrap())
                .unwrap();
            h.set_alignment_start_hap_wrt_ref(0);
            h.set_genome_location(ref_loc.clone());
            h
        })
        .collect();
    build_event_maps_for_haplotypes(&mut haplotypes, REF.as_bytes(), &ref_loc, max_mnp).unwrap();
    let values: Vec<Vec<f64>> = lk
        .iter()
        .map(|f| f[3].split(',').map(|v| v.parse().unwrap()).collect())
        .collect();
    let reads: Vec<BamRecord> = (0..values[0].len())
        .map(|r| BamRecord {
            read_name: format!("r{r}"),
            ..Default::default()
        })
        .collect();
    let likelihoods = AlleleLikelihoods::new(
        SampleList::new(&["s1".to_string()]),
        AlleleList::new(&haplotypes),
        vec![reads],
        vec![values],
    )
    .unwrap();

    let mut out = Vec::new();
    let starts = event_start_positions(&haplotypes).unwrap();
    out.push(format!(
        "starts\t{label}\t{}",
        starts
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    ));
    for &loc in &starts {
        if loc < window_start || loc > window_end {
            continue;
        }
        let events = variants_from_active_haplotypes(loc, &haplotypes, spanning).unwrap();
        let ref_base =
            Allele::create(&REF.as_bytes()[(loc - START) as usize..][..1], true).unwrap();
        let replaced = replace_span_dels(events, &ref_base, loc);
        for c in &replaced {
            out.push(format!(
                "event\t{label}\t{loc}\t{}\t{}\t{}\t{}\t{}",
                c.source,
                c.start,
                c.end,
                text(&c.reference),
                text(&c.alternate)
            ));
        }
        let merged = match make_merged_variant_context(&replaced) {
            Ok(Some(m)) => m,
            Ok(None) => {
                out.push(format!("merged\t{label}\t{loc}\tnull"));
                continue;
            }
            Err(e) => {
                out.push(format!("error\t{label}\t{loc}\t{e:?}"));
                continue;
            }
        };
        out.push(format!(
            "merged\t{label}\t{loc}\t{}\t{}\t{}\t{}",
            merged.source,
            merged.start,
            merged.end,
            merged
                .alleles
                .iter()
                .map(text)
                .collect::<Vec<_>>()
                .join(",")
        ));
        let mapper = create_allele_mapper(&merged, loc, &haplotypes, spanning).unwrap();
        for (allele, haps) in &mapper {
            out.push(format!(
                "mapper\t{label}\t{loc}\t{}\t{}",
                text(allele),
                haps.iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }
        let new_to_old: Vec<(Allele, Vec<Haplotype>)> = mapper
            .iter()
            .map(|(a, hs)| {
                (
                    a.clone(),
                    hs.iter().map(|&h| haplotypes[h].clone()).collect(),
                )
            })
            .collect();
        let marginal = likelihoods.marginalize(&new_to_old).unwrap();
        for a in 0..marginal.number_of_alleles() {
            let row: Vec<String> = (0..marginal.sample_evidence_count(0))
                .map(|r| format!("{:x}", marginal.value(0, a, r).to_bits()))
                .collect();
            out.push(format!(
                "marginal\t{label}\t{loc}\t{a}\t{}\t{}",
                text(marginal.get_allele(a).unwrap()),
                row.join(",")
            ));
        }
    }
    out
}

#[test]
fn every_site_matches_the_reference() {
    let golden = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/allele_mapping.txt.gz"),
    );
    type Case<'a> = (Vec<&'a str>, Vec<Vec<&'a str>>, Vec<Vec<&'a str>>);
    let mut cases: Vec<Case> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "case" => cases.push((f, Vec::new(), Vec::new())),
            "haplotype" => cases.last_mut().unwrap().1.push(f),
            "lk" => cases.last_mut().unwrap().2.push(f),
            _ => expected.push(line),
        }
    }
    let got: Vec<String> = cases.iter().flat_map(|(c, h, l)| run(c, h, l)).collect();
    assert_eq!(got.len(), expected.len(), "{}", got.join("\n"));
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g, e);
    }
}
