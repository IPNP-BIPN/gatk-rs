//! Conformance for `AssemblyComplexity` against the oracle.
//!
//! Golden from `tools/readfilter-conformance/AssemblyComplexityDump.java`. Each case gives the
//! site, the haplotypes (whose event maps are built here against the dump's reference) and the
//! likelihood matrix; this test runs the annotation and renders the `hec`, `hapcomp`, `hapdom` and
//! `error` rows, which must be the golden's in order.

use gatk_annotation::assembly_complexity::{annotate, ComplexityError};
use gatk_engine::allele_likelihoods::AlleleLikelihoods;
use gatk_engine::allele_list::{AlleleList, SampleList};
use gatk_engine::event_map::EventMap;
use gatk_engine::haplotype::Haplotype;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::tsv_table::java_double_to_string;
use htsjdk_bam::record::BamRecord;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::VariantContext;
use std::io::Read;

const REF: &str = "GATTACACGTTGCAAGTCCGATAGCTTAGGCATCGATCGGATCCATGCAAGTCAGT";
const START: i32 = 1000;

#[derive(Default)]
struct Case<'a> {
    fields: Vec<&'a str>,
    haplotypes: Vec<Vec<&'a str>>,
    lk: Vec<Vec<&'a str>>,
}

fn allele(text: &str) -> Allele {
    match text {
        "*" => Allele::create(b"*", false).unwrap(),
        "<NON_REF>" => Allele::create(b"<NON_REF>", false).unwrap(),
        bases => Allele::create(bases.as_bytes(), false).unwrap(),
    }
}

fn run(case: &Case) -> Vec<String> {
    let label = case.fields[1];
    let germline = case.fields[2] == "true";
    let ref_loc = SimpleInterval::new("chr1", START, START + REF.len() as i32 - 1).unwrap();
    let haplotypes: Vec<Haplotype> = case
        .haplotypes
        .iter()
        .map(|f| {
            let mut h = Haplotype::new(f[4].as_bytes(), f[3] == "true").unwrap();
            h.set_cigar(&htsjdk_bam::text_parse::parse_cigar(f[5]).unwrap())
                .unwrap();
            h.set_alignment_start_hap_wrt_ref(0);
            let map = EventMap::from_haplotype(&h, REF.as_bytes(), &ref_loc, 0).unwrap();
            h.set_event_map(map);
            h
        })
        .collect();
    let sample_count = case
        .lk
        .iter()
        .map(|f| f[2].parse::<usize>().unwrap() + 1)
        .max()
        .unwrap();
    let samples: Vec<String> = (1..=sample_count).map(|s| format!("s{s}")).collect();
    let mut values: Vec<Vec<Vec<f64>>> = vec![Vec::new(); sample_count];
    for f in &case.lk {
        let s: usize = f[2].parse().unwrap();
        values[s].push(f[4].split(',').map(|v| v.parse::<f64>().unwrap()).collect());
    }
    let evidence: Vec<Vec<BamRecord>> = values
        .iter()
        .enumerate()
        .map(|(s, v)| {
            (0..v[0].len())
                .map(|r| BamRecord {
                    read_name: format!("s{}r{r}", s + 1),
                    ..Default::default()
                })
                .collect()
        })
        .collect();
    let likelihoods = AlleleLikelihoods::new(
        SampleList::new(&samples),
        AlleleList::new(&haplotypes),
        evidence,
        values,
    )
    .unwrap();
    let mut alleles = vec![Allele::create(case.fields[4].as_bytes(), true).unwrap()];
    alleles.extend(case.fields[5].split(',').map(allele));
    let vc = VariantContext::new("chr1", case.fields[3].parse().unwrap(), alleles);
    match annotate(&vc, &likelihoods, germline) {
        Ok(c) => {
            let join = |v: &[i32]| v.iter().map(i32::to_string).collect::<Vec<_>>().join(",");
            vec![
                format!("hec\t{label}\t{}", join(&c.equivalence_counts)),
                format!("hapcomp\t{label}\t{}", join(&c.edit_distances)),
                format!(
                    "hapdom\t{label}\t{}",
                    c.dominance
                        .iter()
                        .map(|d| format!("{:x}:{}", d.to_bits(), java_double_to_string(*d)))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
            ]
        }
        Err(ComplexityError::NoCarrier) => {
            vec![format!(
                "error\t{label}\tNoSuchElementException: No value present"
            )]
        }
        Err(e) => vec![format!("error\t{label}\t{e:?}")],
    }
}

#[test]
fn every_case_matches_the_reference() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/assembly_complexity.txt.gz");
    let mut golden = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("golden"))
        .read_to_string(&mut golden)
        .expect("golden is gzip");
    let mut cases: Vec<Case> = Vec::new();
    let mut expected: Vec<&str> = Vec::new();
    for line in golden
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "case" => cases.push(Case {
                fields: f,
                ..Default::default()
            }),
            "haplotype" => cases.last_mut().unwrap().haplotypes.push(f),
            "lk" => cases.last_mut().unwrap().lk.push(f),
            _ => expected.push(line),
        }
    }
    let got: Vec<String> = cases.iter().flat_map(run).collect();
    assert_eq!(got.len(), expected.len(), "{}", got.join("\n"));
    for (g, e) in got.iter().zip(&expected) {
        assert_eq!(g, e);
    }
}
