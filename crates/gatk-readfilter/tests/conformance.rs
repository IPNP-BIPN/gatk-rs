//! Conformance for the ported read filters against GATK 4.6.2.0.
//!
//! The golden is a decision matrix produced by `tools/readfilter-conformance/ReadFilterDump.java`
//! in the pinned container: one row per filter, one character per record, taken by the reference
//! through `SAMRecordToGATKReadAdapter`. The corpus travels in the same file, field by field, so
//! this test judges the records the reference judged rather than a reconstruction of them.
//!
//! Rows whose label carries parameters (`MappingQualityReadFilter(min=30,max=60)`) are instances
//! the reference built; the port rebuilds them from the label, so one list of instantiations
//! drives both sides.
//!
//! What this catches that a unit test does not: `NotProperlyPairedReadFilter` is
//! `isPaired() && !isProperlyPaired()`, not the negation of `ProperlyPairedReadFilter`. The first
//! version of the port used the negation, which keeps every unpaired read. The decision matrix
//! disagreed on five records of nineteen the first time it ran.

use gatk_readfilter::{by_name, with_header, Parameterized, PORTED};

use gatk_corpus as corpus;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_filters.txt.gz"),
    )
}

#[test]
fn every_filter_matches_the_reference_decision_for_decision() {
    let text = golden();
    let records = corpus::records(&text);
    let header = corpus::header(&text);

    let mut checked = 0;
    let mut compared = 0;
    for line in text.lines() {
        let mut parts = line.splitn(3, '\t');
        if parts.next() != Some("filter") {
            continue;
        }
        let name = parts.next().expect("a filter row has a name");
        let expected = parts.next().expect("a filter row has decisions");
        assert_eq!(
            expected.len(),
            records.len(),
            "{name}: the golden has {} decisions for {} records",
            expected.len(),
            records.len()
        );

        // A filter in the golden that the port does not implement is a failure, not a skip:
        // silently ignoring it is how a suite shrinks without anyone noticing.
        //
        // A parameterised filter carries its parameters in the label, and the port rebuilds the
        // instance from them, so the reference's own instantiation drives the comparison rather
        // than a second list on this side that could drift from it.
        let ours: String = if let Some(filter) = by_name(name) {
            records
                .iter()
                .map(|read| if filter(read) { '1' } else { '0' })
                .collect()
        } else if let Some(filter) = Parameterized::parse(name) {
            // `E` is the golden's third outcome: the reference threw rather than deciding. A
            // filter that throws stops the tool, so folding it into `0` would hide a crash behind
            // a dropped read.
            records
                .iter()
                .map(|read| match filter.decide(read) {
                    Some(true) => '1',
                    Some(false) => '0',
                    None => 'E',
                })
                .collect()
        } else {
            // The header-dependent family: the label names the filter and its arguments, and the
            // header comes from the golden, so both sides resolve against the same @RG lines.
            let (label, args) = name.split_once('(').expect("a filter label");
            let args = args.strip_suffix(')').expect("a filter label ends with )");
            let values: Vec<String> = args
                .split_once('=')
                .map(|(_, list)| list)
                .unwrap_or("")
                .split(';')
                .filter(|v| !v.is_empty())
                .map(str::to_string)
                .collect();
            records
                .iter()
                .map(|read| {
                    let kept: Option<bool> = match label {
                        "HasReadGroupWithHeader" => Some(gatk_readfilter::has_read_group(read)),
                        "AlignmentAgreesWithHeaderReadFilter" => {
                            Some(with_header::alignment_agrees_with_header(read, &header))
                        }
                        "WellformedReadFilter" => Some(with_header::wellformed(read, &header)),
                        "LibraryReadFilter" => Some(with_header::library(read, &header, &values)),
                        "SampleReadFilter" => Some(with_header::sample(read, &header, &values)),
                        "PlatformReadFilter" => Some(with_header::platform(read, &header, &values)),
                        "PlatformUnitReadFilter" => {
                            Some(with_header::platform_unit(read, &header, &values))
                        }
                        "ReadGroupBlackListReadFilter" => {
                            Some(with_header::read_group_black_list(read, &header, &values))
                        }
                        "IntervalOverlapReadFilter" => {
                            with_header::interval_overlap(read, &header, &values)
                        }
                        "ReadGroupHasFlowOrderReadFilter" => {
                            Some(with_header::read_group_has_flow_order(read, &header))
                        }
                        "WellformedFlowBasedReadFilter" => {
                            with_header::wellformed_flow_based(read, &header)
                        }
                        _ => panic!("{name} is in the golden but not ported; add it or remove it"),
                    };
                    match kept {
                        Some(true) => '1',
                        Some(false) => '0',
                        None => 'E',
                    }
                })
                .collect()
        };

        if ours != expected {
            let differing: Vec<usize> = expected
                .chars()
                .zip(ours.chars())
                .enumerate()
                .filter(|(_, (a, b))| a != b)
                .map(|(i, _)| i)
                .collect();
            let names: Vec<&str> = differing
                .iter()
                .map(|&i| records[i].read_name.as_str())
                .collect();
            panic!(
                "{name} differs on {} of {} records: {names:?}\n  reference: {expected}\n  port     : {ours}",
                differing.len(),
                records.len()
            );
        }
        checked += 1;
        compared += records.len();
    }

    assert!(checked > 0, "the golden carries no filter rows");
    println!("{checked} filters, {compared} decisions, all identical");
}

/// Every ported filter must appear in the golden, or it is untested against the reference.
#[test]
fn no_ported_filter_is_missing_from_the_golden() {
    let text = golden();
    let in_golden: Vec<&str> = text
        .lines()
        .filter_map(|line| {
            let mut parts = line.splitn(3, '\t');
            (parts.next() == Some("filter")).then(|| parts.next().unwrap())
        })
        .collect();

    let missing: Vec<&&str> = PORTED
        .iter()
        .filter(|name| !in_golden.contains(*name))
        .collect();
    assert!(
        missing.is_empty(),
        "ported but not in the golden, so never compared to the reference: {missing:?}. \
         Add them to ReadFilterDump.filters() and regenerate."
    );
}
