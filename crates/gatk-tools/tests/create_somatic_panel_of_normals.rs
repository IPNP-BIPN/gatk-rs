//! Conformance for `CreateSomaticPanelOfNormals` against GATK 4.6.2.0, compared as every panel
//! entry of every run.
//!
//! Golden from `tools/readfilter-conformance/CreateSomaticPanelOfNormalsDump.java`.
//!
//! # What this suite is for
//!
//!  * **a site with no alternate, or only the spanning deletion, being skipped**;
//!  * **a multiallelic site skipping the germline test entirely**;
//!  * **a germline probability of exactly zero where the resource says nothing**, which no
//!    threshold can drop;
//!  * **a high germline frequency removing a het-looking genotype and leaving a low-fraction one**;
//!  * **`--min-sample-count` being compared against the survivors**;
//!  * **FRACTION being over all samples in the header**;
//!  * **and the BETA being fitted by a Brent search over a scale.**

use gatk_corpus as corpus;
use gatk_tools::create_somatic_panel_of_normals::{
    build_panel, fit_beta, germline_probability, has_artifact, variant_genotypes, Genotype,
    PanelEntry, Site, DEFAULT_MAX_GERMLINE_PROBABILITY, DEFAULT_MIN_SAMPLE_COUNT,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/create_somatic_panel_of_normals.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn section(text: &str, kind: &str, name: &str) -> String {
    unescape(
        text.lines()
            .find_map(|line| line.strip_prefix(&format!("{kind}\t{name}=")))
            .unwrap_or_else(|| panic!("the golden carries {kind}/{name}")),
    )
}

/// One input VCF, read as the sites the walker sees.
fn sites(text: &str, which: &str) -> (Vec<Site>, usize) {
    let vcf = section(text, "vcf", which);
    let mut samples: Vec<String> = Vec::new();
    let mut out = Vec::new();
    for line in vcf.lines() {
        if line.starts_with("#CHROM") {
            samples = line.split('\t').skip(9).map(str::to_string).collect();
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        let keys: Vec<&str> = columns[8].split(':').collect();
        let depth_index = keys.iter().position(|key| *key == "AD");
        out.push(Site {
            contig: columns[0].to_string(),
            position: columns[1].parse().expect("a position"),
            reference: columns[3].to_string(),
            alternates: if columns[4] == "." {
                Vec::new()
            } else {
                columns[4].split(',').map(str::to_string).collect()
            },
            genotypes: samples
                .iter()
                .enumerate()
                .map(|(index, sample)| {
                    let parts: Vec<&str> = columns[9 + index].split(':').collect();
                    Genotype {
                        sample: sample.clone(),
                        allele_depths: depth_index.and_then(|at| parts.get(at)).map(|text| {
                            text.split(',')
                                .map(|value| value.parse().expect("a depth"))
                                .collect()
                        }),
                    }
                })
                .collect(),
        });
    }
    (out, samples.len())
}

/// The germline resource, as the frequency it gives one position.
fn germline(text: &str) -> Vec<(i32, f64)> {
    section(text, "germline", "main")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            (
                columns[1].parse().expect("a position"),
                columns[7]
                    .strip_prefix("AF=")
                    .expect("a frequency")
                    .parse()
                    .expect("a frequency"),
            )
        })
        .collect()
}

/// The entries one run wrote, as position, FRACTION and BETA in the writer's own text.
fn measured(text: &str, label: &str) -> Vec<(i32, String, String)> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            let field = |key: &str| {
                columns[7]
                    .split(';')
                    .find_map(|part| part.strip_prefix(&format!("{key}=")))
                    .unwrap_or_else(|| panic!("{label} carries {key}"))
                    .to_string()
            };
            (
                columns[1].parse().expect("a position"),
                field("FRACTION"),
                field("BETA"),
            )
        })
        .collect()
}

/// `VCFEncoder.formatVCFDouble`, which is how the writer renders both fields.
fn format(value: f64) -> String {
    htsjdk_vcf::variant::format_vcf_double(value)
}

fn rendered(entry: &PanelEntry) -> (i32, String, String) {
    (
        entry.position,
        format(entry.fraction),
        format!("{},{}", format(entry.beta.alpha), format(entry.beta.beta)),
    )
}

/// label, input, germline resource used, minimum samples, maximum germline probability.
type Run = (&'static str, &'static str, bool, usize, f64);

fn runs() -> Vec<Run> {
    vec![
        (
            "default",
            "main",
            false,
            DEFAULT_MIN_SAMPLE_COUNT,
            DEFAULT_MAX_GERMLINE_PROBABILITY,
        ),
        (
            "germline",
            "main",
            true,
            DEFAULT_MIN_SAMPLE_COUNT,
            DEFAULT_MAX_GERMLINE_PROBABILITY,
        ),
        (
            "germline-permissive",
            "main",
            true,
            DEFAULT_MIN_SAMPLE_COUNT,
            1.0,
        ),
        (
            "germline-strict",
            "main",
            true,
            DEFAULT_MIN_SAMPLE_COUNT,
            0.0001,
        ),
        (
            "min-one",
            "main",
            false,
            1,
            DEFAULT_MAX_GERMLINE_PROBABILITY,
        ),
        (
            "min-three",
            "main",
            false,
            3,
            DEFAULT_MAX_GERMLINE_PROBABILITY,
        ),
        (
            "four-samples",
            "four",
            false,
            DEFAULT_MIN_SAMPLE_COUNT,
            DEFAULT_MAX_GERMLINE_PROBABILITY,
        ),
    ]
}

#[test]
fn every_panel_entry_matches_the_golden() {
    let text = golden();
    let resource = germline(&text);
    let mut compared = 0;
    for (label, which, use_resource, minimum, maximum) in runs() {
        let (sites, sample_count) = sites(&text, which);
        let frequency = |site: &Site| -> f64 {
            if !use_resource {
                return 0.0;
            }
            resource
                .iter()
                .filter(|(position, _)| *position == site.position)
                .map(|(_, value)| *value)
                .sum()
        };
        let produced: Vec<(i32, String, String)> =
            build_panel(&sites, sample_count, frequency, minimum, maximum)
                .iter()
                .map(rendered)
                .collect();
        assert_eq!(produced, measured(&text, label), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 7, "the runs the golden carries");
}

/// No alternate at all, or only the spanning deletion.
#[test]
fn two_site_shapes_are_skipped() {
    let text = golden();
    let (sites, _) = sites(&text, "main");
    let at = |position: i32| {
        sites
            .iter()
            .find(|site| site.position == position)
            .unwrap_or_else(|| panic!("a site at {position}"))
    };
    assert!(at(1000).is_skipped(), "no alternate");
    assert!(at(2000).is_skipped(), "only the spanning deletion");
    assert!(!at(3000).is_skipped());
    // Which the golden shows by their absence from every run.
    for (label, ..) in runs() {
        let positions: Vec<i32> = measured(&text, label)
            .into_iter()
            .map(|(position, _, _)| position)
            .collect();
        assert!(!positions.contains(&1000), "{label}");
        assert!(!positions.contains(&2000), "{label}");
    }
}

/// Every genotype counts, whatever its counts say.
#[test]
fn a_multiallelic_site_skips_the_germline_test() {
    let text = golden();
    let (sites, _) = sites(&text, "main");
    let multiallelic = sites
        .iter()
        .find(|site| site.position == 7000)
        .expect("the multiallelic site");
    assert!(multiallelic.is_multiallelic());
    // All three genotypes, including the one carrying no alternate at all.
    assert_eq!(variant_genotypes(multiallelic, 0.4, 0.5).len(), 3);
    // The biallelic site with the same genotype shape keeps only the two that carry an alternate.
    let biallelic = sites
        .iter()
        .find(|site| site.position == 6000)
        .expect("the half-fraction site");
    assert_eq!(variant_genotypes(biallelic, 0.0, 0.5).len(), 2);
    // And at the same germline frequency the multiallelic one survives while it does not.
    let positions = |label: &str| -> Vec<i32> {
        measured(&text, label)
            .into_iter()
            .map(|(position, _, _)| position)
            .collect()
    };
    assert!(positions("germline").contains(&7000));
    assert!(!positions("germline").contains(&6000));
    assert!(positions("default").contains(&6000));
}

/// Exactly zero, which is below every threshold.
#[test]
fn a_site_the_resource_does_not_mention_can_never_be_dropped() {
    let text = golden();
    assert_eq!(germline_probability(0.0, 2, 20), 0.0);
    assert_eq!(germline_probability(1e-9, 2, 20), 0.0, "below negligible");
    assert_eq!(germline_probability(1.5, 2, 20), 0.0, "and above one");
    let genotype = Genotype {
        sample: "n1".to_string(),
        allele_depths: Some(vec![18, 2]),
    };
    assert!(
        has_artifact(&genotype, 0.0, 0.0001),
        "zero is below any threshold"
    );

    // Which is why the strictest run drops nothing the resource had not already removed.
    assert_eq!(
        measured(&text, "germline"),
        measured(&text, "germline-strict")
    );
    // A genotype with no alternate read never counts, whatever the frequency.
    let no_alt = Genotype {
        sample: "n3".to_string(),
        allele_depths: Some(vec![20, 0]),
    };
    assert!(!has_artifact(&no_alt, 0.0, 1.0));
}

/// The half-fraction site goes and the tenth-fraction one stays, at the same frequency.
#[test]
fn a_germline_frequency_removes_a_het_and_leaves_a_low_fraction() {
    let text = golden();
    let het = Genotype {
        sample: "n1".to_string(),
        allele_depths: Some(vec![10, 10]),
    };
    let low = Genotype {
        sample: "n1".to_string(),
        allele_depths: Some(vec![1800, 200]),
    };
    assert!(!has_artifact(&het, 0.4, DEFAULT_MAX_GERMLINE_PROBABILITY));
    assert!(has_artifact(&low, 0.4, DEFAULT_MAX_GERMLINE_PROBABILITY));
    // And raising the threshold to one brings the het back.
    assert!(has_artifact(&het, 0.4, 1.0));

    let positions = |label: &str| -> Vec<i32> {
        measured(&text, label)
            .into_iter()
            .map(|(position, _, _)| position)
            .collect()
    };
    assert!(!positions("germline").contains(&6000));
    assert!(positions("germline-permissive").contains(&6000));
    assert!(
        positions("germline").contains(&9000),
        "the deep low-fraction site"
    );
}

/// Against the survivors, and FRACTION is over the header's samples.
#[test]
fn the_minimum_counts_survivors_and_the_fraction_counts_samples() {
    let text = golden();
    let positions = |label: &str| -> Vec<i32> {
        measured(&text, label)
            .into_iter()
            .map(|(position, _, _)| position)
            .collect()
    };
    // The singleton site appears only when the minimum is one.
    assert!(!positions("default").contains(&4000));
    assert!(positions("min-one").contains(&4000));
    // And the two-sample sites go when the minimum is three.
    assert!(positions("default").contains(&3000));
    assert!(!positions("min-three").contains(&3000));

    // A fourth sample that carries nothing moves FRACTION and leaves BETA alone.
    let three = measured(&text, "default");
    let four = measured(&text, "four-samples");
    let entry = |rows: &[(i32, String, String)], position: i32| {
        rows.iter()
            .find(|(at, _, _)| *at == position)
            .expect("an entry")
            .clone()
    };
    let (_, fraction_three, beta_three) = entry(&three, 3000);
    let (_, fraction_four, beta_four) = entry(&four, 3000);
    assert_eq!(fraction_three, "0.667");
    assert_eq!(fraction_four, "0.500");
    assert_eq!(beta_three, beta_four, "the fit sees only the survivors");
}

/// The base shape is the empirical mean and the scale comes from a Brent search.
#[test]
fn the_beta_is_fitted_by_a_brent_search() {
    let text = golden();
    // The scale cancels out of the ratio, so whatever the search returns, alpha over beta is the
    // alternate total plus one over the reference total plus one. NOT the empirical mean of the
    // counts: the two pseudocounts move it.
    let shape = fit_beta(&[(2, 18), (3, 17)]);
    assert!((shape.alpha / shape.beta - 6.0 / 36.0).abs() < 1e-12);
    assert!(
        (shape.alpha / (shape.alpha + shape.beta) - 5.0 / 40.0).abs() > 0.01,
        "and that is not the fraction 5/40 the counts themselves give"
    );
    // Which is what the golden wrote for that site.
    let (_, _, beta) = measured(&text, "default")
        .into_iter()
        .find(|(position, _, _)| *position == 3000)
        .expect("the site");
    assert_eq!(
        beta,
        format!(
            "{},{}",
            htsjdk_vcf::variant::format_vcf_double(shape.alpha),
            htsjdk_vcf::variant::format_vcf_double(shape.beta)
        )
    );
    // Deeper counts at the same fraction keep the ratio rule and land on a different total, which
    // is the search rather than the base shape: the deeper site's beta is the SMALLER of the two.
    let deep = fit_beta(&[(200, 1800), (300, 1700)]);
    assert!((deep.alpha / deep.beta - 501.0 / 3501.0).abs() < 1e-12);
    assert!(deep.alpha + deep.beta < shape.alpha + shape.beta);

    // A genotype without AD counted for the site and contributes nothing to the fit: the site with
    // one such genotype is fitted from the two that have counts.
    let without = fit_beta(&[(2, 18), (4, 16)]);
    let (_, _, measured_beta) = measured(&text, "default")
        .into_iter()
        .find(|(position, _, _)| *position == 8000)
        .expect("the site");
    assert_eq!(
        measured_beta,
        format!(
            "{},{}",
            htsjdk_vcf::variant::format_vcf_double(without.alpha),
            htsjdk_vcf::variant::format_vcf_double(without.beta)
        )
    );
}
