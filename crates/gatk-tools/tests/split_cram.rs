//! Conformance for `SplitCRAM` against GATK 4.6.2.0, compared as the shard each container lands in
//! and the name each shard is given.
//!
//! Golden from `tools/readfilter-conformance/SplitCRAMDump.java`.
//!
//! # What this suite is for
//!
//!  * **the threshold being a minimum**, so a shard overshoots by up to one container, and being
//!    tested strictly, so a threshold exactly the size of a container does not;
//!  * **`--shard-max-output-count` limiting nothing above one**, the counter it compares being
//!    reset for every shard;
//!  * **the names**, which are `String.format` on a template whose counter starts at zero;
//!  * **and the two refusals**, which happen before the input is opened.
//!
//! The bytes of a shard are htsjdk's, not the tool's, so what the golden carries for each shard is
//! the record count of every container it holds and the read names it gives back when read. The
//! counts are what this compares against; the names are what says the shard is a whole CRAM.

use gatk_corpus as corpus;
use gatk_tools::split_cram::{accepts_template, format_name, plan, Shard, SplitError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/split_cram.txt.gz"),
    )
}

/// The record count of every container of the input, which is what the tool plans over.
fn containers(text: &str) -> Vec<i32> {
    counts(
        text.lines()
            .find_map(|line| line.strip_prefix("fixture\tinput\tcontainers="))
            .expect("the golden carries the input's containers"),
    )
}

fn counts(field: &str) -> Vec<i32> {
    if field.is_empty() {
        return Vec::new();
    }
    field
        .split(',')
        .map(|count| count.parse().expect("a record count"))
        .collect()
}

/// Every shard of one run, in the order the dump listed them, which is by name.
fn shards(text: &str, label: &str) -> Vec<Shard> {
    let prefix = format!("shard\t{label}\t");
    text.lines()
        .filter_map(|line| line.strip_prefix(prefix.as_str()))
        .map(|rest| {
            let fields: Vec<&str> = rest.split('\t').collect();
            let containers = counts(
                fields[1]
                    .strip_prefix("records=")
                    .expect("a shard names its record counts"),
            );
            let names = fields[2]
                .strip_prefix("names=")
                .expect("a shard names its reads");
            // Every record of every container of the shard came back when it was read, which is
            // what says the shard is a CRAM on its own rather than a slice of one.
            let read: i64 = names.split(',').filter(|name| !name.is_empty()).count() as i64;
            assert_eq!(
                read,
                containers.iter().map(|count| *count as i64).sum::<i64>(),
                "shard {} of {label} reads back every record it holds",
                fields[0]
            );
            Shard {
                name: fields[0].to_string(),
                containers,
            }
        })
        .collect()
}

fn refusal(text: &str, label: &str) -> String {
    let prefix = format!("error\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .expect("the golden carries the refusal")
        .to_string()
}

#[test]
fn a_threshold_of_one_gives_one_container_per_shard() {
    let text = golden();
    assert_eq!(
        plan(&containers(&text), 1, 0, "shard_%d.cram"),
        Ok(shards(&text, "one-per-shard"))
    );
}

#[test]
fn a_threshold_a_container_reaches_exactly_ends_the_shard() {
    let text = golden();
    let planned = plan(&containers(&text), 3, 0, "exact_%d.cram").expect("a plan");
    assert_eq!(planned, shards(&text, "exact-threshold"));
    assert!(planned.iter().all(|shard| shard.containers.len() == 1));
}

#[test]
fn a_threshold_one_above_it_takes_a_second_container() {
    let text = golden();
    let planned = plan(&containers(&text), 4, 0, "pair_%d.cram").expect("a plan");
    assert_eq!(planned, shards(&text, "overshoot"));
    // Six records for a threshold of four, which is the overshoot.
    assert_eq!(planned[0].records(), 6);
}

#[test]
fn a_threshold_no_shard_reaches_is_one_shard() {
    let text = golden();
    assert_eq!(
        plan(&containers(&text), 1000, 0, "all_%d.cram"),
        Ok(shards(&text, "one-shard"))
    );
}

#[test]
fn only_a_maximum_of_one_ever_stops_the_run() {
    let text = golden();
    let input = containers(&text);
    assert_eq!(
        plan(&input, 1, 1, "max1_%d.cram"),
        Ok(shards(&text, "max-one"))
    );
    assert_eq!(plan(&input, 1, 1, "max1_%d.cram").expect("a plan").len(), 1);
    // Two and three leave every shard standing, the counter being reset for each of them.
    assert_eq!(
        plan(&input, 1, 2, "max2_%d.cram"),
        Ok(shards(&text, "max-two"))
    );
    assert_eq!(
        plan(&input, 1, 3, "max3_%d.cram"),
        Ok(shards(&text, "max-three"))
    );
    assert_eq!(plan(&input, 1, 3, "max3_%d.cram").expect("a plan").len(), 5);
}

#[test]
fn a_padded_template_is_accepted_and_padded() {
    let text = golden();
    assert_eq!(
        plan(&containers(&text), 1000, 0, "padded_%04d.cram"),
        Ok(shards(&text, "padded-template"))
    );
    assert_eq!(
        format_name("padded_%04d.cram", 0).as_deref(),
        Some("padded_0000.cram")
    );
}

#[test]
fn a_width_with_a_flag_is_not_a_formatter() {
    let text = golden();
    // The tool is handed the whole path, and names it back whole.
    let template = "<dir>/flagged-template/flagged_%-4d.cram";
    let error = plan(&containers(&text), 1000, 0, template).expect_err("a refusal");
    assert_eq!(
        error,
        SplitError::TemplateMissingFormatter {
            template: template.to_string()
        }
    );
    assert!(!accepts_template(template));
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "flagged-template")
    );
}

#[test]
fn a_template_with_no_formatter_is_refused_before_anything_is_read() {
    let text = golden();
    let template = "<dir>/no-formatter/plain.cram";
    let error = plan(&containers(&text), 1000, 0, template).expect_err("a refusal");
    assert_eq!(
        format!("{}:{}", error.java_class(), error.message()),
        refusal(&text, "no-formatter")
    );
}

#[test]
fn an_empty_cram_produces_no_shard() {
    let text = golden();
    let empty = text
        .lines()
        .find_map(|line| line.strip_prefix("fixture\tempty\tcontainers="))
        .expect("the golden carries the empty input");
    assert!(counts(empty).is_empty());
    assert_eq!(plan(&counts(empty), 1, 0, "empty_%d.cram"), Ok(Vec::new()));
    assert!(shards(&text, "empty-input").is_empty());
}
