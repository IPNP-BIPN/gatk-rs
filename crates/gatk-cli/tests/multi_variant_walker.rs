//! The `multi-variant-walker` suite: `MultiVariantWalker`, its data source and its two grouping
//! subclasses, measured through the three probe walkers of `MultiVariantWalkerDump`.
//!
//! The probes are not tools, so each run is replayed here through the same calls a tool makes:
//! the inputs written to disk from the golden and indexed by the port's own `IndexFeatureFile`
//! where the dump indexed them, opened through `runners::multi_variant::open`, merged, and grouped
//! by `gatk_tools::multi_variant`. Every row of every run is compared, in order.

use gatk_cli::runners::multi_variant::{self, DrivingVariants};
use gatk_engine::interval::SimpleInterval;
use gatk_tools::multi_variant::{self as engine, FeatureInput, Located};
use htsjdk_vcf::variant::VariantContext;

fn golden() -> String {
    gatk_corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/multi_variant_walker.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

#[derive(Clone, Copy, PartialEq)]
enum Probe {
    Plain,
    Throwing,
    OnStart,
    Overlap,
}

/// One run of the dump, as its command line set it up.
struct Run {
    label: &'static str,
    probe: Probe,
    inputs: Vec<(Option<&'static str>, &'static str)>,
    cross_validate: bool,
    interval: Option<(&'static str, i32, i32)>,
    distance: i64,
    max_distance: i64,
    padding: i64,
    ignore_outside: bool,
    ignore_reference_blocks: bool,
}

fn run(label: &'static str, probe: Probe, inputs: &[&'static str]) -> Run {
    Run {
        label,
        probe,
        inputs: inputs.iter().map(|input| (None, *input)).collect(),
        cross_validate: true,
        interval: None,
        distance: if probe == Probe::Overlap { 10000 } else { 0 },
        max_distance: i64::from(i32::MAX),
        padding: 1,
        ignore_outside: false,
        ignore_reference_blocks: false,
    }
}

fn runs() -> Vec<Run> {
    use Probe::*;
    let grouped = ["g1", "g2"];
    let overlapping = ["o1", "o2"];
    vec![
        run("interleaved", Plain, &["a", "b"]),
        run("one-input", Plain, &["b"]),
        run("ties-two", Plain, &["t1", "t2"]),
        run("ties-three", Plain, &["t1", "t2", "t3"]),
        run("ties-three-reversed", Plain, &["t3", "t2", "t1"]),
        Run {
            interval: Some(("chr1", 150, 350)),
            ..run("intervals", Plain, &["a", "b"])
        },
        run("dictionaries-union", Plain, &["d1", "d2"]),
        run("dictionaries-union-reversed", Plain, &["d2", "d1"]),
        run("dictionaries-reordered", Plain, &["d1", "reordered"]),
        run("length-mismatch", Plain, &["a", "longer"]),
        Run {
            cross_validate: false,
            ..run("length-mismatch-not-validated", Plain, &["a", "longer"])
        },
        run("md5-mismatch", Plain, &["md5a", "md5b"]),
        run("md5-one-side", Plain, &["md5a", "a"]),
        run("no-contigs-indexed", Plain, &["no-contigs-indexed", "a"]),
        run("no-contigs-indexed-alone", Plain, &["no-contigs-indexed"]),
        run("no-contigs-beside-one", Plain, &["no-contigs", "a"]),
        run(
            "no-dictionary-anywhere",
            Plain,
            &["no-contigs", "no-contigs-too"],
        ),
        run("unknown-contig", Plain, &["a", "unknown-contig"]),
        run("unsorted", Plain, &["unsorted", "b"]),
        run("samples", Plain, &["samples-a", "samples-b"]),
        run("header-type-conflict", Plain, &["samples-b", "samples-c"]),
        run("duplicate-input", Plain, &["a", "a"]),
        Run {
            inputs: vec![(Some("one"), "a"), (Some("two"), "a")],
            ..run("same-file-two-tags", Plain, &[])
        },
        Run {
            inputs: vec![(Some("same"), "a"), (Some("same"), "b")],
            ..run("two-files-one-tag", Plain, &[])
        },
        run("throws-in-apply", Throwing, &["a", "b"]),
        run("on-start-default", OnStart, &grouped),
        Run {
            distance: 5,
            ..run("on-start-distance-5", OnStart, &grouped)
        },
        Run {
            distance: 40,
            ..run("on-start-distance-40", OnStart, &grouped)
        },
        Run {
            distance: 40,
            max_distance: 3,
            ..run("on-start-max-distance-3", OnStart, &grouped)
        },
        Run {
            padding: 3,
            ..run("on-start-ref-padding-3", OnStart, &grouped)
        },
        Run {
            interval: Some(("chr1", 101, 200)),
            ..run("on-start-intervals", OnStart, &grouped)
        },
        Run {
            interval: Some(("chr1", 101, 200)),
            ignore_outside: true,
            ..run("on-start-ignore-outside", OnStart, &grouped)
        },
        Run {
            interval: Some(("chr1", 400, 500)),
            ..run("on-start-empty", OnStart, &grouped)
        },
        run("overlap-default", Overlap, &overlapping),
        Run {
            ignore_reference_blocks: true,
            ..run("overlap-ignore-reference-blocks", Overlap, &overlapping)
        },
        Run {
            interval: Some(("chr1", 105, 300)),
            ignore_outside: true,
            ..run("overlap-ignore-outside", Overlap, &overlapping)
        },
    ]
}

/// The inputs written to a scratch directory, indexed by the port where the dump indexed them.
fn write_inputs(text: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join("gatk-cli-multivariantwalker-dump");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("a scratch directory");
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("input\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.splitn(3, '\t').collect();
        let path = dir.join(format!("{}.vcf", fields[0]));
        std::fs::write(&path, unescape(fields[2])).expect("an input");
        if fields[1] == "indexed" {
            let path = path.to_string_lossy().to_string();
            let outcome = gatk_cli::run(&["IndexFeatureFile".to_string(), "-I".to_string(), path]);
            assert_eq!(outcome.status, 0, "{}", outcome.stderr);
        }
    }
    dir
}

fn reference(text: &str) -> Vec<(String, String)> {
    let fasta = unescape(
        text.lines()
            .find_map(|line| line.strip_prefix("fasta\t"))
            .expect("the reference"),
    );
    let mut contigs: Vec<(String, String)> = Vec::new();
    for line in fasta.lines() {
        match line.strip_prefix('>') {
            Some(name) => contigs.push((name.to_string(), String::new())),
            None => contigs.last_mut().expect("a contig").1.push_str(line),
        }
    }
    contigs
}

struct Record<'a>(&'a VariantContext);

impl Located for Record<'_> {
    fn contig(&self) -> &str {
        &self.0.contig
    }
    fn start(&self) -> i64 {
        self.0.start
    }
    fn end(&self) -> i64 {
        self.0.stop
    }
}

fn mine(setup: &Run, dir: &std::path::Path, bases: &[(String, String)]) -> Vec<String> {
    let label = setup.label;
    let masked = |message: &str| message.replace(&dir.to_string_lossy().to_string(), "<dir>");
    let mut rows = Vec::new();
    let error = |rows: &mut Vec<String>, thrown: gatk_tools::main_entry::Thrown| {
        rows.push(format!(
            "error\t{label}\t{}:{}",
            thrown.exception,
            escape(&masked(thrown.message.as_deref().unwrap_or("null")))
        ));
    };
    let inputs: Vec<FeatureInput> = setup
        .inputs
        .iter()
        .map(|(tag, file)| {
            let path = dir
                .join(format!("{file}.vcf"))
                .to_string_lossy()
                .to_string();
            FeatureInput {
                tag: tag.map(str::to_string),
                raw: path.clone(),
                absolute: path,
            }
        })
        .collect();
    let driving: DrivingVariants = match multi_variant::open(inputs, setup.cross_validate) {
        Ok(driving) => driving,
        Err(thrown) => {
            error(&mut rows, thrown);
            return rows;
        }
    };

    // `onTraversalStart`.
    rows.push(format!(
        "dict\t{label}\t{}",
        match &driving.merged.dictionary {
            None => "null".to_string(),
            Some(dictionary) => dictionary
                .iter()
                .map(|contig| format!("{}:{}", contig.name, contig.length))
                .collect::<Vec<_>>()
                .join(","),
        }
    ));
    rows.push(format!(
        "samples\t{label}\t{}",
        driving.merged.samples.join(",")
    ));
    let mut header = String::from("##fileformat=VCFv4.2\n");
    for line in engine::sorted_metadata(&driving.merged.lines) {
        header.push_str("##");
        header.push_str(&line.render());
        header.push('\n');
    }
    rows.push(format!("header\t{label}\t{}", escape(&masked(&header))));

    let spans: Option<Vec<SimpleInterval>> = setup.interval.map(|(contig, start, end)| {
        vec![SimpleInterval::new(contig, start, end).expect("an interval")]
    });
    let (order, failure) = driving
        .merged(spans.as_deref())
        .expect("every input is indexed");
    let user_intervals: Option<Vec<(String, i64, i64)>> = setup
        .interval
        .map(|(contig, start, end)| vec![(contig.to_string(), i64::from(start), i64::from(end))]);

    let mut step = 0;
    let mut group = |rows: &mut Vec<String>, members: Vec<Record>| {
        let contig = members[0].contig().to_string();
        let sequence = &bases
            .iter()
            .find(|(name, _)| *name == contig)
            .expect("a contig")
            .1;
        let (start, end) = engine::spanning_window(&members, setup.padding, sequence.len() as i64);
        let ids: Vec<&str> = members.iter().map(|member| member.0.id.as_str()).collect();
        rows.push(format!(
            "group\t{label}\t{step}\t{}\t{contig}:{start}-{end}\t{}",
            ids.join(","),
            &sequence[(start - 1) as usize..end as usize]
        ));
        step += 1;
    };
    let mut record_step = 0;
    let mut on_start = engine::GroupedOnStart::new(setup.distance, setup.max_distance);
    let mut by_overlap = engine::GroupedByOverlap::default();
    for (input, record) in order {
        match setup.probe {
            Probe::Plain | Probe::Throwing => {
                if setup.probe == Probe::Throwing && record.id == "b200" {
                    error(&mut rows, driving.apply_failure(input, record));
                    return rows;
                }
                rows.push(format!(
                    "record\t{label}\t{record_step}\t{}:{}-{}\t{}",
                    record.contig, record.start, record.stop, record.id
                ));
                record_step += 1;
            }
            Probe::OnStart | Probe::Overlap => {
                if setup.ignore_outside
                    && !engine::starts_within(
                        user_intervals.as_deref(),
                        &record.contig,
                        record.start,
                    )
                {
                    continue;
                }
                let closed = if setup.probe == Probe::OnStart {
                    on_start.accept(Record(record))
                } else {
                    let alternates = record.alternate_alleles();
                    if setup.ignore_reference_blocks
                        && alternates.len() == 1
                        && alternates[0].display_string() == "<NON_REF>"
                    {
                        continue;
                    }
                    by_overlap.accept(Record(record))
                };
                if let Some(members) = closed {
                    group(&mut rows, members);
                }
            }
        }
    }
    if let Some(thrown) = failure {
        error(&mut rows, thrown);
        return rows;
    }
    let last = match setup.probe {
        Probe::OnStart => on_start.finish(),
        Probe::Overlap => by_overlap.finish(),
        _ => None,
    };
    if let Some(members) = last {
        group(&mut rows, members);
    }
    rows
}

fn expected(text: &str, label: &str) -> Vec<String> {
    text.lines()
        .filter(|line| {
            ["dict", "samples", "header", "record", "group", "error"]
                .iter()
                .any(|kind| line.starts_with(&format!("{kind}\t{label}\t")))
        })
        .map(str::to_string)
        .collect()
}

#[test]
fn every_run_matches_the_golden_row_for_row() {
    let text = golden();
    let dir = write_inputs(&text);
    let bases = reference(&text);
    let all = runs();
    let mut labels: Vec<&str> = text
        .lines()
        .filter_map(|line| {
            let mut fields = line.splitn(3, '\t');
            let kind = fields.next()?;
            ["dict", "error"].contains(&kind).then(|| fields.next())?
        })
        .collect();
    labels.dedup();
    assert_eq!(
        labels,
        all.iter().map(|setup| setup.label).collect::<Vec<_>>(),
        "every run of the golden has a setup here"
    );
    for setup in &all {
        assert_eq!(
            mine(setup, &dir, &bases),
            expected(&text, setup.label),
            "the rows of {}",
            setup.label
        );
    }
}
