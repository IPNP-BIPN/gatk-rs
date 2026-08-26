//! Conformance for `GroupedSVCluster` against GATK 4.6.2.0, compared as the stratum and member set
//! of every cluster of every run, and as the message of every refusal.
//!
//! Golden from `tools/readfilter-conformance/GroupedSVClusterDump.java`.
//!
//! # What this suite is for
//!
//!  * **each stratum clustering with its own thresholds**, in one run and from one engine;
//!  * **a record matching two strata being refused**, with a different message from `SVStratify`'s;
//!  * **an unmatched record not being clustered at all**;
//!  * **the two configurations being checked twice**, once by count and once by name;
//!  * **an empty stratification configuration being refused** before either engine is built;
//!  * **and the doubled-number column message appearing a second time.**

use gatk_corpus as corpus;
use gatk_tools::grouped_sv_cluster::{
    check_columns, linkage_for, run, Cluster, Engines, GroupedError, StratumParameters,
};
use gatk_tools::sv_cluster::{is_reciprocal_overlap, Algorithm, CallRecord};
use gatk_tools::sv_stratify::{
    parse_integer_maybe_null, Engine, Stratum, SvType, Thresholds, Tracks,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/grouped_sv_cluster.txt.gz"),
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

/// `error\t<label>\t<class>:<message>`, as the class and the message.
fn refusal(text: &str, label: &str) -> (String, String) {
    let line = text
        .lines()
        .find_map(|line| line.strip_prefix(&format!("error\t{label}\t")))
        .unwrap_or_else(|| panic!("the golden carries error/{label}"));
    let (class, message) = line.split_once(':').expect("a class and a message");
    (class.to_string(), unescape(message))
}

/// A stratification table, read as the engine reads it.
fn strata(text: &str, name: &str) -> Engine {
    let table = section(text, "strata", name);
    let mut rows = table.lines();
    rows.next().expect("a header");
    let strata: Vec<Stratum> = rows
        .filter(|line| !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Stratum::new(
                columns[0],
                SvType::parse(columns[1]).expect("a known type"),
                parse_integer_maybe_null(columns[2]),
                parse_integer_maybe_null(columns[3]),
                Vec::new(),
            )
            .expect("a valid stratum")
        })
        .collect();
    Engine::new(strata, Tracks::new(&[], &[]).expect("no tracks")).expect("a valid engine")
}

/// A clustering table, read as `StratifiedClusteringTableParser` reads it.
fn clustering(text: &str, name: &str) -> Engines {
    let table = section(text, "cluster", name);
    let mut rows = table.lines();
    let header: Vec<String> = rows
        .next()
        .expect("a header")
        .split('\t')
        .map(str::to_string)
        .collect();
    check_columns(&header).expect("a well formed header");
    Engines::new(
        &rows
            .filter(|line| !line.is_empty())
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                StratumParameters {
                    name: columns[0].to_string(),
                    reciprocal_overlap: columns[1].parse().expect("an overlap"),
                    size_similarity: columns[2].parse().expect("a similarity"),
                    breakend_window: columns[3].parse().expect("a window"),
                    sample_overlap: columns[4].parse().expect("a sample overlap"),
                }
            })
            .collect::<Vec<StratumParameters>>(),
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
            contig_b: columns[0].to_string(),
            position_b: end,
            strand_a: None,
            strand_b: None,
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

/// Every record one run wrote, as its stratum and its member set.
fn measured(text: &str, label: &str) -> Vec<Cluster> {
    section(text, "out", label)
        .lines()
        .filter(|line| !line.starts_with("#CHROM") && !line.is_empty())
        .map(|line| {
            let info = line.split('\t').nth(7).expect("an INFO column");
            let field = |key: &str| {
                info.split(';')
                    .find_map(|part| part.strip_prefix(&format!("{key}=")))
                    .unwrap_or_else(|| panic!("{label} carries {key}"))
            };
            Cluster {
                stratum: field("STRAT").to_string(),
                members: field("MEMBERS").split(',').map(str::to_string).collect(),
            }
        })
        .collect()
}

/// Compared without their order, which is the writer's business: the tool disables its own output
/// index because it does not guarantee one.
fn sorted(mut clusters: Vec<Cluster>) -> Vec<Cluster> {
    clusters.sort_by(|a, b| (&a.stratum, &a.members).cmp(&(&b.stratum, &b.members)));
    clusters
}

fn thresholds() -> Thresholds {
    Thresholds {
        overlap_fraction: 1.0,
        num_breakpoint_overlaps: 0,
        num_breakpoint_overlaps_interchrom: 0,
    }
}

#[test]
fn every_cluster_matches_the_golden() {
    let text = golden();
    let records = records(&text);
    let engine = strata(&text, "main");
    let engines = clustering(&text, "main");
    let mut compared = 0;
    for (label, algorithm) in [
        ("default", Algorithm::SingleLinkage),
        ("max-clique", Algorithm::MaxClique),
    ] {
        let produced =
            run(&records, &engine, &engines, thresholds(), algorithm, false).expect(label);
        assert_eq!(sorted(produced), sorted(measured(&text, label)), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 2, "the runs that wrote records");
}

/// Two large deletions that a half-overlap threshold would join stay apart under their own
/// stratum's 0.99, in the same run in which two small ones cluster.
#[test]
fn each_stratum_clusters_with_its_own_thresholds() {
    let text = golden();
    let produced = sorted(measured(&text, "default"));
    assert!(produced.contains(&Cluster {
        stratum: "DEL_small".to_string(),
        members: vec!["small1".to_string(), "small2".to_string()],
    }));
    for id in ["large1", "large2"] {
        assert!(produced.contains(&Cluster {
            stratum: "DEL_large".to_string(),
            members: vec![id.to_string()],
        }));
    }

    // And the thresholds really are the difference: the pair that clusters is the one whose own
    // row is loose enough, and the other row would have kept it apart.
    let records = records(&text);
    let engines = clustering(&text, "main");
    let of = |id: &str| {
        records
            .iter()
            .find(|record| record.id == id)
            .expect(id)
            .clone()
    };
    let under = |name: &str, a: &CallRecord, b: &CallRecord| {
        linkage_for(engines.get(name).expect(name), false).are_clusterable(a, b)
    };
    let (small1, small2) = (of("small1"), of("small2"));
    assert!(
        under("DEL_small", &small1, &small2),
        "its own row joins them"
    );
    assert!(
        !under("DEL_large", &small1, &small2),
        "the other row would not: 0.70 overlap against 0.99"
    );

    // The large pair fails on the window rather than on the overlap, so neither row would have
    // joined it: what its stratum decides is only that it is judged at all.
    let (large1, large2) = (of("large1"), of("large2"));
    assert!(!under("DEL_large", &large1, &large2));
    assert!(!under("DEL_small", &large1, &large2));
    assert!(
        is_reciprocal_overlap(20000, 30000, 21000, 31000, 0.5),
        "0.90 overlap clears 0.5"
    );
    assert!(
        large2.position_a - large1.position_a > 500,
        "and the window is 500"
    );
}

/// It is written straight out, with its own id as its only member.
#[test]
fn an_unmatched_record_is_not_clustered_at_all() {
    let text = golden();
    assert!(sorted(measured(&text, "default")).contains(&Cluster {
        stratum: "default".to_string(),
        members: vec!["unmatched".to_string()],
    }));

    // Nothing in the configuration names it, which is why.
    let engine = strata(&text, "main");
    assert!(!engine
        .strata
        .iter()
        .any(|stratum| stratum.sv_type == SvType::Ins));
}

/// A different message from `SVStratify`'s, which offers `--allow-multiple-matches` instead, and
/// the groups are listed in `java.util.HashMap` order over their names.
#[test]
fn a_record_matching_two_strata_is_refused_outright() {
    let text = golden();
    let (class, message) = refusal(&text, "multiple-matches");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.GATKException"
    );
    let produced = run(
        &records(&text),
        &strata(&text, "overlapping"),
        // The clustering table this run was given is not reported beside the strata, because the
        // refusal happens before any record is clustered. Restated here from the run's arguments,
        // where all it has to do is name the same two groups.
        &Engines::new(
            &["DEL_a", "DEL_b"]
                .iter()
                .map(|name| StratumParameters {
                    name: (*name).to_string(),
                    reciprocal_overlap: 0.5,
                    size_similarity: 0.0,
                    breakend_window: 500,
                    sample_overlap: 0.0,
                })
                .collect::<Vec<StratumParameters>>(),
        ),
        thresholds(),
        Algorithm::SingleLinkage,
        false,
    )
    .expect_err("two strata match");
    assert_eq!(produced.message(), message);
    assert!(
        matches!(&produced, GroupedError::MultipleMatches { names, .. } if names == &["DEL_b", "DEL_a"]),
        "the HashMap order over the names"
    );
    assert!(!message.contains("--allow-multiple-matches"), "not a flag");
}

/// Once by count, once by name, with two different messages.
#[test]
fn the_two_configurations_must_name_the_same_groups() {
    let text = golden();
    let records = records(&text);
    let engine = strata(&text, "main");
    for (label, name) in [("too-few-groups", "short"), ("group-not-found", "renamed")] {
        let (class, message) = refusal(&text, label);
        assert_eq!(class, "java.lang.IllegalStateException", "{label}");
        let produced = run(
            &records,
            &engine,
            &clustering(&text, name),
            thresholds(),
            Algorithm::SingleLinkage,
            false,
        )
        .expect_err(label);
        assert_eq!(produced.message(), message, "{label}");
    }
    // The renamed table holds the right NUMBER of groups, which is what makes the second check a
    // second check rather than the same one.
    assert_eq!(clustering(&text, "renamed").entries.len(), 3);
    assert_eq!(clustering(&text, "short").entries.len(), 2);
}

/// Refused before either engine is built, so no record is ever consulted.
#[test]
fn an_empty_stratification_configuration_is_refused() {
    let text = golden();
    let (class, message) = refusal(&text, "no-strata");
    assert_eq!(class, "java.lang.IllegalStateException");
    let engine = strata(&text, "empty");
    assert!(engine.strata.is_empty());
    let produced = run(
        &[],
        &engine,
        &clustering(&text, "main"),
        thresholds(),
        Algorithm::SingleLinkage,
        false,
    )
    .expect_err("no strata");
    assert_eq!(produced.message(), message);
    // An empty table is refused for being empty, not for having the wrong count against a
    // three-group clustering table.
    assert_eq!(produced, GroupedError::NoStrata);
}

/// A second copy of the mistake the stratification parser makes: both halves read the count that
/// was found, so the numbers can never differ.
#[test]
fn the_column_count_message_prints_the_same_number_twice() {
    let text = golden();
    let (class, message) = refusal(&text, "extra-column");
    assert_eq!(
        class,
        "org.broadinstitute.hellbender.exceptions.UserException$BadInput"
    );
    let header: Vec<String> = [
        "NAME",
        "RECIPROCAL_OVERLAP",
        "SIZE_SIMILARITY",
        "BREAKEND_WINDOW",
        "SAMPLE_OVERLAP",
        "EXTRA",
    ]
    .iter()
    .map(|column| column.to_string())
    .collect();
    let produced = check_columns(&header).expect_err("an extra column");
    assert_eq!(produced, GroupedError::ColumnCount { count: 6 });
    assert!(message.ends_with(produced.message().as_str()), "{message}");
    assert_eq!(produced.message(), "Expected 6 columns but found 6");

    // The positive control: the same header without the extra column is accepted.
    assert!(check_columns(&header[..5]).is_ok());
    // And a missing column is reported by name, in the order the column set holds them.
    assert_eq!(
        check_columns(&header[1..]).expect_err("no NAME"),
        GroupedError::MissingColumn {
            column: "NAME".to_string()
        }
    );
}

/// Unlike `SVCluster`, one row feeds all three parameter sets, so which set a pair takes stops
/// mattering.
#[test]
fn one_row_becomes_all_three_parameter_sets() {
    let text = golden();
    let engines = clustering(&text, "main");
    let linkage = linkage_for(engines.get("DEL_large").expect("the large group"), false);
    assert_eq!(linkage.depth.reciprocal_overlap, 0.99);
    assert_eq!(linkage.mixed.reciprocal_overlap, 0.99);
    assert_eq!(linkage.pesr.reciprocal_overlap, 0.99);
    assert_eq!(linkage.depth.window, 100);
    assert_eq!(linkage.pesr.window, 100);
}
