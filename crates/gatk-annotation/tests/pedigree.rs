//! Conformance for the de novo calls and the transmitted singletons, against the oracle.
//!
//! Golden from `tools/annotation-conformance/PedigreeAnnotationDump.java`.
//!
//! ```text
//! denovo    shallow-child    hiConfDeNovo=[kid][java.util.ArrayList]
//! denovo    no-depth         loConfDeNovo=[kid][java.util.ArrayList]
//! singleton shallow-parents  transmittedSingleton=[mom][java.util.ArrayList]
//! ```
//!
//! Those three rows are the two depth facts. `shallow-child` passes the high-confidence branch
//! because the default depth threshold is zero and a depth of zero clears it; `no-depth` fails it
//! because an absent `DP` reads as minus one; and `shallow-parents` is emitted because all three of
//! `TransmittedSingleton`'s depth tests read the **child**.

use std::io::Read;

use gatk_annotation::pedigree::{
    genotype_type, is_violation, possible_de_novo, transmitted_singleton, GenotypeType, Trio,
};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::{Genotype, VariantContext};

const START: i64 = 105;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/pedigree.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn reference() -> Allele {
    Allele::from_str("A", true).expect("an allele")
}

fn alternate() -> Allele {
    Allele::from_str("C", false).expect("an allele")
}

fn alleles_for(kind: &str) -> Vec<Allele> {
    match kind {
        "hom-ref" => vec![reference(), reference()],
        "het" => vec![reference(), alternate()],
        "hom-var" => vec![alternate(), alternate()],
        "no-call" => vec![Allele::no_call(), Allele::no_call()],
        "half-called" => vec![reference(), Allele::no_call()],
        "empty" => vec![],
        other => panic!("{other} is not a genotype kind"),
    }
}

fn genotype(name: &str, kind: &str, gq: Option<i32>, dp: Option<i32>) -> Genotype {
    let mut g = Genotype::new(name, alleles_for(kind));
    g.gq = gq;
    g.dp = dp;
    g
}

fn trio_kinds(label: &str) -> [&'static str; 3] {
    match label {
        "de-novo-het" => ["hom-ref", "hom-ref", "het"],
        "de-novo-hom-var" => ["hom-ref", "hom-ref", "hom-var"],
        "inherited-het" => ["het", "hom-ref", "het"],
        "all-hom-ref" => ["hom-ref", "hom-ref", "hom-ref"],
        "mother-no-call" => ["no-call", "hom-ref", "hom-var"],
        "father-no-call" => ["hom-ref", "no-call", "hom-var"],
        "child-no-call" => ["hom-ref", "hom-ref", "no-call"],
        "mother-half-called" => ["half-called", "hom-ref", "het"],
        "both-parents-no-call" => ["no-call", "no-call", "het"],
        "hom-var-parent" => ["hom-var", "hom-ref", "hom-ref"],
        other => panic!("{other} has no trio fixture"),
    }
}

fn trios(count: usize) -> Vec<Trio> {
    (0..count)
        .map(|i| {
            let suffix = if i == 0 { String::new() } else { i.to_string() };
            Trio {
                family_id: format!("fam{i}"),
                maternal_id: format!("mom{suffix}"),
                paternal_id: format!("dad{suffix}"),
                child_id: format!("kid{suffix}"),
            }
        })
        .collect()
}

fn de_novo_context(label: &str) -> VariantContext {
    let alleles = if label == "multiallelic" {
        vec![
            reference(),
            alternate(),
            Allele::from_str("G", false).expect("an allele"),
        ]
    } else {
        vec![reference(), alternate()]
    };
    let mut vc = VariantContext::new("chr1", START, alleles);
    vc.stop = START;
    let g = &mut vc.genotypes;
    match label {
        "high-confidence" | "multiallelic" => {
            g.push(genotype("mom", "hom-ref", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
        }
        "low-confidence-gq" => {
            g.push(genotype("mom", "hom-ref", Some(5), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(5), Some(30)));
            g.push(genotype("kid", "het", Some(15), Some(30)));
        }
        "no-depth" => {
            g.push(genotype("mom", "hom-ref", Some(50), None));
            g.push(genotype("dad", "hom-ref", Some(50), None));
            g.push(genotype("kid", "het", Some(50), None));
        }
        "shallow-child" => {
            g.push(genotype("mom", "hom-ref", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(0)));
        }
        "low-parent-gq" => {
            g.push(genotype("mom", "hom-ref", Some(15), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(15), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
        }
        "inherited" => {
            g.push(genotype("mom", "het", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
        }
        "no-gq" => {
            g.push(genotype("mom", "hom-ref", None, Some(30)));
            g.push(genotype("dad", "hom-ref", None, Some(30)));
            g.push(genotype("kid", "het", None, Some(30)));
        }
        "common-allele" => {
            g.push(genotype("mom", "hom-ref", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
            for i in 0..4 {
                g.push(genotype(&format!("x{i}"), "het", Some(50), Some(30)));
            }
        }
        "two-trios" => {
            g.push(genotype("mom", "hom-ref", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
            g.push(genotype("mom1", "hom-ref", Some(50), Some(30)));
            g.push(genotype("dad1", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid1", "het", Some(15), Some(30)));
        }
        other => panic!("{other} has no de novo fixture"),
    }
    vc
}

/// The dump's singleton fixtures, with the `AC` attribute they carry.
fn singleton_context(label: &str) -> (VariantContext, i32) {
    let mut vc = VariantContext::new("chr1", START, vec![reference(), alternate()]);
    vc.stop = START;
    let mut allele_count = 2;
    let g = &mut vc.genotypes;
    match label {
        "transmitted" | "ac-three" => {
            g.push(genotype("mom", "het", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
            if label == "ac-three" {
                allele_count = 3;
            }
        }
        "non-transmitted" => {
            g.push(genotype("mom", "het", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "hom-ref", Some(50), Some(30)));
            allele_count = 1;
        }
        "shallow-parents" => {
            g.push(genotype("mom", "het", Some(50), Some(1)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(1)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
        }
        "shallow-child" => {
            g.push(genotype("mom", "het", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(1)));
        }
        "low-call-rate" => {
            g.push(genotype("mom", "het", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
            for i in 0..7 {
                g.push(genotype(&format!("x{i}"), "hom-ref", Some(5), Some(30)));
            }
        }
        "child-hom-var" => {
            g.push(genotype("mom", "het", Some(50), Some(30)));
            g.push(genotype("dad", "hom-ref", Some(50), Some(30)));
            g.push(genotype("kid", "hom-var", Some(50), Some(30)));
        }
        "both-parents-het" => {
            g.push(genotype("mom", "het", Some(50), Some(30)));
            g.push(genotype("dad", "het", Some(50), Some(30)));
            g.push(genotype("kid", "het", Some(50), Some(30)));
        }
        other => panic!("{other} has no singleton fixture"),
    }
    (vc, allele_count)
}

/// `key=[a, b][java.util.ArrayList]` joined with `;`, in the map's insertion order.
fn rendered(entries: &[(String, Vec<String>)]) -> String {
    entries
        .iter()
        .map(|(key, values)| format!("{key}=[{}][java.util.ArrayList]", values.join(", ")))
        .collect::<Vec<_>>()
        .join(";")
}

#[test]
fn every_genotype_type_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("type\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label and a type");
        let ours = match genotype_type(&genotype("s", label, Some(30), Some(20))) {
            GenotypeType::Unavailable => "UNAVAILABLE",
            GenotypeType::NoCall => "NO_CALL",
            GenotypeType::HomRef => "HOM_REF",
            GenotypeType::Het => "HET",
            GenotypeType::HomVar => "HOM_VAR",
            GenotypeType::Mixed => "MIXED",
        };
        assert_eq!(ours, expected, "type of {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no type rows");
    println!("{count} genotype types identical");
}

#[test]
fn every_mendelian_violation_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("violation\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label and a verdict");
        let kinds = trio_kinds(label);
        let mother = genotype("mom", kinds[0], Some(50), Some(30));
        let father = genotype("dad", kinds[1], Some(50), Some(30));
        let child = genotype("kid", kinds[2], Some(50), Some(30));
        assert_eq!(
            is_violation(&mother, &father, &child).to_string(),
            expected,
            "violation of {label}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no violation rows");
    println!("{count} Mendelian verdicts identical");
}

#[test]
fn every_de_novo_call_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("denovo\t") else {
            continue;
        };
        let (label, expected) = match rest.split_once('\t') {
            Some(pair) => pair,
            None => (rest, ""),
        };
        let vc = de_novo_context(label);
        let trio_count = if label == "two-trios" { 2 } else { 1 };
        // The defaults: a parent GQ threshold of twenty and a depth threshold of zero.
        let ours = rendered(&possible_de_novo(&vc, &trios(trio_count), 20, 0));
        assert_eq!(ours, expected, "de novo on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no de novo rows");
    println!("{count} de novo answers identical");
}

#[test]
fn every_singleton_call_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("singleton\t") else {
            continue;
        };
        let (label, expected) = match rest.split_once('\t') {
            Some(pair) => pair,
            None => (rest, ""),
        };
        let (vc, allele_count) = singleton_context(label);
        let ours = rendered(&transmitted_singleton(&vc, &trios(1), allele_count));
        assert_eq!(ours, expected, "singleton on {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no singleton rows");
    println!("{count} singleton answers identical");
}

#[test]
fn the_parents_depths_are_never_read() {
    // The same trio, once with shallow parents and once with a shallow child. Only the second
    // changes the answer, which is the reference's copy-paste and not an accident of this fixture.
    let (shallow_parents, ac) = singleton_context("shallow-parents");
    let (shallow_child, _) = singleton_context("shallow-child");
    assert!(!transmitted_singleton(&shallow_parents, &trios(1), ac).is_empty());
    assert!(transmitted_singleton(&shallow_child, &trios(1), ac).is_empty());
}
