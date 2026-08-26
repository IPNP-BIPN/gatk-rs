//! Conformance for `SVCluster` against GATK 4.6.2.0, compared as the member set of every cluster
//! of every run.
//!
//! Golden from `tools/readfilter-conformance/SVClusterDump.java`.
//!
//! Reading the VCF and collapsing a cluster into a representative record are not ported. Which
//! records belong together is, which is what `MEMBERS` holds.
//!
//! # What this suite is for
//!
//!  * **the parameter set being chosen by the pair**, not by the user;
//!  * **single linkage and max clique disagreeing on a chain** the type rule makes;
//!  * **`--enable-cnv` being observable only under max clique**;
//!  * **overlap and size similarity being separate tests**;
//!  * **sample overlap being zero by default**;
//!  * **and `--omit-members` removing the only field that shows any of it**.

use gatk_corpus as corpus;
use gatk_tools::sv_cluster::{
    cluster, has_sample_overlap, is_reciprocal_overlap, test_size_similarity, Algorithm,
    CallRecord, ClusteringParameters, Linkage,
};
use gatk_tools::sv_stratify::SvType;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sv_cluster.txt.gz"),
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

/// The measured input, read out of the VCF the golden carries.
fn records(text: &str) -> Vec<CallRecord> {
    let vcf = section(text, "vcf", "input");
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
        let info: Vec<(&str, &str)> = columns[7]
            .split(';')
            .filter_map(|part| part.split_once('='))
            .collect();
        let field = |key: &str| {
            info.iter()
                .find(|(name, _)| *name == key)
                .map(|(_, value)| *value)
        };
        let sv_type = SvType::parse(field("SVTYPE").expect("a type")).expect("a known type");
        let end: i32 = field("END").expect("an end").parse().expect("an end");
        let strands = field("STRANDS");
        // The carriers are the samples whose GT holds a non-reference allele.
        let keys: Vec<&str> = columns[8].split(':').collect();
        let genotype_index = keys.iter().position(|key| *key == "GT").expect("a GT");
        let carriers: Vec<String> = samples
            .iter()
            .enumerate()
            .filter(|(index, _)| {
                columns[9 + index]
                    .split(':')
                    .nth(genotype_index)
                    .expect("a genotype")
                    .split(['/', '|'])
                    .any(|allele| allele != "0" && allele != ".")
            })
            .map(|(_, sample)| sample.clone())
            .collect();
        out.push(CallRecord {
            id: columns[2].to_string(),
            sv_type,
            contig_a: columns[0].to_string(),
            position_a: columns[1].parse().expect("a position"),
            contig_b: field("CHR2").unwrap_or(columns[0]).to_string(),
            position_b: field("END2")
                .map(|value| value.parse().expect("a second position"))
                .unwrap_or(end),
            strand_a: strands.map(|value| value.starts_with('+')),
            strand_b: strands.map(|value| value.ends_with('+')),
            length: match sv_type {
                SvType::Ins | SvType::Bnd | SvType::Ctx => None,
                _ => Some(field("SVLEN").expect("a length").parse().expect("a length")),
            },
            algorithms: field("ALGORITHMS")
                .expect("an algorithms field")
                .split(',')
                .map(str::to_string)
                .collect(),
            carriers,
        });
    }
    out
}

/// Every cluster one run wrote, as its `MEMBERS` set, in output order.
fn measured(text: &str, label: &str) -> Vec<Vec<String>> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            columns[7]
                .split(';')
                .find_map(|part| part.strip_prefix("MEMBERS="))
                .unwrap_or_else(|| panic!("{label} carries MEMBERS"))
                .split(',')
                .map(str::to_string)
                .collect()
        })
        .collect()
}

/// Cluster sets compared without their order, which is the collapser's business rather than the
/// engine's.
fn sorted(mut clusters: Vec<Vec<String>>) -> Vec<Vec<String>> {
    clusters.sort();
    clusters
}

fn defaults() -> Linkage {
    Linkage::default()
}

/// label, linkage, algorithm.
fn runs() -> Vec<(&'static str, Linkage, Algorithm)> {
    let base = defaults();
    vec![
        ("single-linkage", base, Algorithm::SingleLinkage),
        ("max-clique", base, Algorithm::MaxClique),
        (
            "enable-cnv",
            Linkage {
                cluster_del_with_dup: true,
                ..base
            },
            Algorithm::SingleLinkage,
        ),
        (
            "max-clique-enable-cnv",
            Linkage {
                cluster_del_with_dup: true,
                ..base
            },
            Algorithm::MaxClique,
        ),
        (
            "depth-overlap-high",
            Linkage {
                depth: ClusteringParameters::depth(0.99, 0.0, 10_000_000, 0.0),
                ..base
            },
            Algorithm::SingleLinkage,
        ),
        (
            "depth-overlap-low",
            Linkage {
                depth: ClusteringParameters::depth(0.1, 0.1, 10_000_000, 0.0),
                ..base
            },
            Algorithm::SingleLinkage,
        ),
        (
            "sample-overlap",
            Linkage {
                depth: ClusteringParameters::depth(0.8, 0.0, 10_000_000, 0.5),
                ..base
            },
            Algorithm::SingleLinkage,
        ),
        (
            "wide-window",
            Linkage {
                pesr: ClusteringParameters::pesr(0.5, 0.0, 5000, 0.0),
                ..base
            },
            Algorithm::SingleLinkage,
        ),
    ]
}

#[test]
fn every_cluster_matches_the_golden() {
    let text = golden();
    let records = records(&text);
    let mut compared = 0;
    for (label, linkage, algorithm) in runs() {
        assert_eq!(
            sorted(cluster(&records, &linkage, algorithm)),
            sorted(measured(&text, label)),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 8, "the runs that wrote member ids");
}

/// The chain is made by the type rule: a deletion and a duplication do not cluster, both cluster
/// with the CNV between them, and the two algorithms answer differently.
#[test]
fn the_two_algorithms_disagree_on_a_chain_the_types_make() {
    let text = golden();
    let single = sorted(measured(&text, "single-linkage"));
    let clique = sorted(measured(&text, "max-clique"));
    assert!(single.contains(&vec!["c1".to_string(), "e1".to_string(), "u1".to_string()]));
    assert!(clique.contains(&vec!["c1".to_string(), "u1".to_string()]));
    assert!(clique.contains(&vec!["c1".to_string(), "e1".to_string()]));
    assert!(!clique.contains(&vec!["c1".to_string(), "e1".to_string(), "u1".to_string()]));

    // And the direct link the CNV was standing in for.
    let records = records(&text);
    let of = |id: &str| records.iter().find(|r| r.id == id).expect(id);
    assert!(
        !defaults().are_clusterable(of("e1"), of("u1")),
        "DEL with DUP"
    );
    assert!(
        defaults().are_clusterable(of("c1"), of("u1")),
        "CNV with DUP"
    );
    assert!(
        defaults().are_clusterable(of("c1"), of("e1")),
        "CNV with DEL"
    );
    assert!(
        Linkage {
            cluster_del_with_dup: true,
            ..defaults()
        }
        .are_clusterable(of("e1"), of("u1")),
        "DEL with DUP once enabled"
    );
}

/// Enabling it changes nothing under single linkage, because the CNV already chained them.
#[test]
fn enable_cnv_is_only_observable_under_max_clique() {
    let text = golden();
    assert_eq!(
        sorted(measured(&text, "single-linkage")),
        sorted(measured(&text, "enable-cnv"))
    );
    assert_ne!(
        sorted(measured(&text, "max-clique")),
        sorted(measured(&text, "max-clique-enable-cnv"))
    );
    assert_eq!(
        sorted(measured(&text, "max-clique-enable-cnv")),
        sorted(measured(&text, "single-linkage"))
    );
}

/// Two deletions sharing a start and differing twentyfold in length overlap and still do not
/// cluster, because size similarity is a separate test.
#[test]
fn overlap_and_size_similarity_are_separate_tests() {
    let text = golden();
    let records = records(&text);
    let of = |id: &str| records.iter().find(|r| r.id == id).expect(id);
    let (z1, z2) = (of("z1"), of("z2"));
    assert_eq!(z1.position_a, z2.position_a, "the same start");
    assert!(
        is_reciprocal_overlap(50000, 51000, 50000, 51000, 0.8),
        "a control"
    );
    assert!(!defaults().are_clusterable(z1, z2));

    // Lowering both thresholds is what joins them, and the golden shows they are apart by default.
    assert!(sorted(measured(&text, "single-linkage")).contains(&vec!["z1".to_string()]));
    assert!(test_size_similarity(1001, 20001, 0.04));
    assert!(!test_size_similarity(1001, 20001, 0.1));
}

/// Sample overlap is zero by default, so two records sharing no carrier cluster until it is
/// raised.
#[test]
fn sample_overlap_is_zero_by_default() {
    let text = golden();
    let records = records(&text);
    let of = |id: &str| records.iter().find(|r| r.id == id).expect(id);
    let (s1, s2) = (of("s1"), of("s2"));
    assert!(s1
        .carriers
        .iter()
        .all(|sample| !s2.carriers.contains(sample)));
    assert!(has_sample_overlap(s1, s2, 0.0), "zero asks nothing");
    assert!(!has_sample_overlap(s1, s2, 0.5));

    assert!(sorted(measured(&text, "single-linkage"))
        .contains(&vec!["s1".to_string(), "s2".to_string()]));
    assert!(sorted(measured(&text, "sample-overlap")).contains(&vec!["s1".to_string()]));
}

/// Widening the PESR window twentyfold changes nothing, because every factory requires overlap AND
/// proximity whatever the class documentation says.
#[test]
fn proximity_alone_never_suffices() {
    let text = golden();
    assert_eq!(
        sorted(measured(&text, "single-linkage")),
        sorted(measured(&text, "wide-window"))
    );
    for parameters in [
        ClusteringParameters::depth(0.8, 0.0, 1000, 0.0),
        ClusteringParameters::mixed(0.8, 0.0, 1000, 0.0),
        ClusteringParameters::pesr(0.5, 0.0, 500, 0.0),
    ] {
        assert!(
            parameters.requires_overlap_and_proximity,
            "every factory passes true"
        );
    }
}

/// The run that proves the others were reading something real: with the field omitted there is
/// nothing to read.
#[test]
fn omitting_the_members_removes_the_only_evidence() {
    let text = golden();
    let body = section(&text, "out", "omit-members");
    assert!(!body.contains("MEMBERS="));
    // And the same run still wrote the same number of records.
    let rows = body
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .count();
    assert_eq!(rows, measured(&text, "single-linkage").len());
}
