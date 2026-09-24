//! Conformance for `AlleleFrequencyCalculator` against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/AlleleFrequencyCalculatorDump.java`. Every double is
//! compared as its raw bits, so a fit that converges one round later or sums in another order
//! fails here rather than in a QUAL three tools away.
//!
//! # What this suite is for
//!
//!  * **the flat first round and the pseudocounts after it**, which decide every QUAL;
//!  * **the SNP and indel pseudocounts chosen by length**;
//!  * **the biallelic short-circuit**, where the allele-absent probability is the no-variant one;
//!  * **a spanning deletion counted as no variant**;
//!  * **hom-refs with only a GQ given likelihoods, and other genotypes without them skipped**;
//!  * **the MLE counts as rounded effective counts, and the epsilon in `passesThreshold`**;
//!  * **and the three refusals.**

use gatk_corpus as corpus;
use gatk_engine::allele_frequency_calculator::{AlleleFrequencyCalculator, Priors};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::variant::Genotype;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/allele_frequency_calculator.txt.gz"),
    )
}

fn allele(bases: &str, reference: bool) -> Allele {
    Allele::create(bases.as_bytes(), reference).expect("an allele")
}

fn a() -> Allele {
    allele("A", true)
}
fn c() -> Allele {
    allele("C", false)
}
fn g() -> Allele {
    allele("G", false)
}
fn t() -> Allele {
    allele("T", false)
}
fn at() -> Allele {
    allele("AT", false)
}
fn span() -> Allele {
    allele("*", false)
}

fn pl(sample: &str, alleles: Vec<Allele>, pls: &[i32]) -> Genotype {
    let mut genotype = Genotype::new(sample, alleles);
    genotype.pl = Some(pls.to_vec());
    genotype
}

/// A hom-ref of the given ploidy carrying a GQ and no likelihoods.
fn gq(sample: &str, ploidy: usize, quality: i32) -> Genotype {
    let mut genotype = Genotype::new(sample, vec![a(); ploidy]);
    genotype.gq = Some(quality);
    genotype
}

fn no_call(sample: &str, pls: Option<&[i32]>) -> Genotype {
    let mut genotype = Genotype::missing(sample, 2);
    genotype.pl = pls.map(|p| p.to_vec());
    genotype
}

fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn row(
    calculator: &AlleleFrequencyCalculator,
    label: &str,
    alleles: Vec<Allele>,
    genotypes: Vec<Genotype>,
) -> String {
    match calculator.calculate("chr1", 100, &alleles, &genotypes) {
        Ok(result) => {
            let counts: Vec<String> = result
                .allele_counts_of_mle
                .iter()
                .map(|c| c.to_string())
                .collect();
            let absent: Vec<String> = result
                .log10_p_ref_by_allele
                .iter()
                .map(|&x| bits(x))
                .collect();
            let passes: Vec<String> = (0..alleles.len() - 1)
                .map(|alt| {
                    [10.0, 30.0, 50.0]
                        .iter()
                        .map(|&q| {
                            if result.passes_threshold(alt, q) {
                                "1"
                            } else {
                                "0"
                            }
                        })
                        .collect()
                })
                .collect();
            format!(
                "calc\t{label}\t{}\t{}\t{}\t{}\t{}",
                bits(result.log10_posterior_of_no_variant),
                bits(result.log10_prob_variant_present()),
                counts.join(","),
                absent.join(","),
                passes.join(",")
            )
        }
        Err(error) => format!("error\t{label}\t{}:{}", error.java_class(), error.message()),
    }
}

fn produced() -> Vec<String> {
    let standard = AlleleFrequencyCalculator::make_calculator(&Priors::default());
    let wide = AlleleFrequencyCalculator::make_calculator(&Priors {
        snp_heterozygosity: 0.01,
        indel_heterozygosity: 0.0005,
        heterozygosity_standard_deviation: 0.05,
        ..Priors::default()
    });
    let s = &standard;
    let mut rows = vec![
        row(
            s,
            "het",
            vec![a(), c()],
            vec![pl("s1", vec![a(), c()], &[300, 0, 500])],
        ),
        row(
            s,
            "hom-var",
            vec![a(), c()],
            vec![pl("s1", vec![c(), c()], &[900, 60, 0])],
        ),
        row(
            s,
            "reference-best",
            vec![a(), c()],
            vec![pl("s1", vec![a(), a()], &[0, 30, 400])],
        ),
        row(
            s,
            "marginal",
            vec![a(), c()],
            vec![pl("s1", vec![a(), c()], &[8, 0, 300])],
        ),
        row(
            s,
            "flat",
            vec![a(), c()],
            vec![pl("s1", vec![a(), c()], &[0, 0, 0])],
        ),
        row(
            s,
            "cohort",
            vec![a(), c()],
            vec![
                pl("s1", vec![a(), c()], &[250, 0, 480]),
                pl("s2", vec![a(), a()], &[0, 45, 700]),
                pl("s3", vec![c(), c()], &[820, 55, 0]),
                pl("s4", vec![a(), c()], &[120, 0, 350]),
            ],
        ),
    ];
    let mut rare = vec![pl("s0", vec![a(), c()], &[200, 0, 400])];
    for i in 1..20 {
        rare.push(pl(&format!("s{i}"), vec![a(), a()], &[0, 30 + i, 450]));
    }
    rows.push(row(s, "rare-het", vec![a(), c()], rare));
    rows.push(row(
        s,
        "all-hom-ref",
        vec![a(), c()],
        vec![
            pl("s1", vec![a(), a()], &[0, 60, 900]),
            pl("s2", vec![a(), a()], &[0, 42, 630]),
            pl("s3", vec![a(), a()], &[0, 21, 315]),
        ],
    ));
    rows.push(row(
        s,
        "snp-and-snp",
        vec![a(), c(), g()],
        vec![pl("s1", vec![c(), g()], &[900, 400, 500, 450, 0, 800])],
    ));
    rows.push(row(
        s,
        "snp-and-indel",
        vec![a(), c(), at()],
        vec![
            pl("s1", vec![a(), c()], &[300, 0, 500, 280, 450, 900]),
            pl("s2", vec![a(), at()], &[260, 350, 700, 0, 380, 640]),
        ],
    ));
    let del_ref = allele("ACG", true);
    let del_alt = allele("A", false);
    rows.push(row(
        s,
        "indel-only",
        vec![del_ref.clone(), del_alt.clone()],
        vec![pl("s1", vec![del_ref, del_alt], &[400, 0, 600])],
    ));
    rows.push(row(
        s,
        "five-alleles",
        vec![a(), c(), g(), t(), at()],
        vec![
            pl(
                "s1",
                vec![c(), g()],
                &[
                    900, 400, 500, 450, 0, 800, 600, 700, 650, 900, 620, 710, 660, 910, 990,
                ],
            ),
            pl(
                "s2",
                vec![a(), t()],
                &[
                    300, 350, 600, 380, 610, 640, 0, 400, 420, 500, 330, 440, 460, 520, 800,
                ],
            ),
        ],
    ));
    rows.push(row(
        s,
        "haploid",
        vec![a(), c()],
        vec![pl("s1", vec![c()], &[200, 0])],
    ));
    rows.push(row(
        s,
        "triploid",
        vec![a(), c()],
        vec![pl("s1", vec![a(), a(), c()], &[150, 0, 90, 400])],
    ));
    rows.push(row(
        s,
        "tetraploid",
        vec![a(), c()],
        vec![pl("s1", vec![a(), a(), c(), c()], &[400, 60, 0, 70, 500])],
    ));
    rows.push(row(
        s,
        "mixed-ploidy",
        vec![a(), c()],
        vec![
            pl("s1", vec![c()], &[200, 0]),
            pl("s2", vec![a(), c()], &[300, 0, 500]),
        ],
    ));
    rows.push(row(
        s,
        "span-del-only",
        vec![a(), span()],
        vec![pl("s1", vec![a(), span()], &[300, 0, 500])],
    ));
    rows.push(row(
        s,
        "span-del-and-snp",
        vec![a(), span(), c()],
        vec![
            pl("s1", vec![a(), span()], &[300, 0, 500, 320, 520, 900]),
            pl("s2", vec![a(), c()], &[280, 400, 700, 0, 350, 600]),
        ],
    ));
    rows.push(row(
        s,
        "gq-hom-ref",
        vec![a(), c()],
        vec![pl("s1", vec![a(), c()], &[300, 0, 500]), gq("s2", 2, 40)],
    ));
    rows.push(row(
        s,
        "gq-hom-ref-triallelic",
        vec![a(), c(), g()],
        vec![
            pl("s1", vec![c(), g()], &[900, 400, 500, 450, 0, 800]),
            gq("s2", 2, 25),
        ],
    ));
    rows.push(row(
        s,
        "gq-haploid-skipped",
        vec![a(), c()],
        vec![pl("s1", vec![a(), c()], &[300, 0, 500]), gq("s2", 1, 40)],
    ));
    rows.push(row(
        s,
        "no-call-with-pl",
        vec![a(), c()],
        vec![
            pl("s1", vec![a(), c()], &[300, 0, 500]),
            no_call("s2", Some(&[0, 10, 100])),
        ],
    ));
    rows.push(row(
        s,
        "no-call-without-pl",
        vec![a(), c()],
        vec![
            pl("s1", vec![a(), c()], &[300, 0, 500]),
            no_call("s2", None),
        ],
    ));
    rows.push(row(
        &wide,
        "wide-prior-het",
        vec![a(), c()],
        vec![pl("s1", vec![a(), c()], &[300, 0, 500])],
    ));
    rows.push(row(
        &wide,
        "wide-prior-snp-and-indel",
        vec![a(), c(), at()],
        vec![
            pl("s1", vec![a(), c()], &[300, 0, 500, 280, 450, 900]),
            pl("s2", vec![a(), at()], &[260, 350, 700, 0, 380, 640]),
        ],
    ));
    rows.push(row(
        s,
        "only-reference",
        vec![a()],
        vec![pl("s1", vec![a(), a()], &[0])],
    ));
    rows.push(row(
        s,
        "no-likelihoods",
        vec![a(), c()],
        vec![gq("s1", 2, 40)],
    ));
    rows.push(row(
        s,
        "inconsistent-pl",
        vec![a(), c(), g()],
        vec![pl("s1", vec![a(), c()], &[300, 0, 500])],
    ));
    rows
}

#[test]
fn every_fit_matches_the_golden() {
    let expected: Vec<String> = golden()
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(str::to_string)
        .collect();
    let produced = produced();
    assert_eq!(produced.len(), expected.len(), "the rows the dump wrote");
    for (produced, expected) in produced.iter().zip(&expected) {
        assert_eq!(produced, expected);
    }
}
