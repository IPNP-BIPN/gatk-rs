//! Conformance for `ReferenceContext` against GATK 4.6.2.0.
//!
//! One row per accessor per context, replayed against the FASTA that travels in the golden. The
//! per-accessor shape is what makes the failures legible: a composite row would collapse to `E`
//! the moment one call throws, and *which* call throws is the measurement.

use gatk_corpus as corpus;
use gatk_engine::context::ReferenceContext;
use gatk_engine::interval::SimpleInterval;
use gatk_engine::reference::ReferenceFileSource;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/reference_context.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => out.push(other),
            None => {}
        }
    }
    out
}

fn field<'a>(text: &'a str, kind: &str) -> &'a str {
    text.lines()
        .find_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .unwrap_or_else(|| panic!("the golden carries no {kind} row"))
}

/// A label is `<contig>:<start>-<end>+<leading>,<trailing>/<accessor>`.
fn parse_label(label: &str) -> Option<(String, i32, i32, i32, i32, String)> {
    let (span, accessor) = label.split_once('/')?;
    let (location, window) = span.split_once('+')?;
    let (contig, range) = location.split_once(':')?;
    let (start, end) = range.split_once('-')?;
    let (leading, trailing) = window.split_once(',')?;
    Some((
        contig.to_string(),
        start.parse().ok()?,
        end.parse().ok()?,
        leading.parse().ok()?,
        trailing.parse().ok()?,
        accessor.to_string(),
    ))
}

fn bytes(value: Vec<u8>) -> String {
    String::from_utf8(value).expect("the reference returns ASCII bases")
}

#[test]
fn every_accessor_answers_what_the_reference_answers() {
    let text = golden();

    let dir = std::env::temp_dir().join(format!("gatk-rs-refcontext-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("ref.fasta");
    std::fs::write(&fasta, unescape(field(&text, "fasta"))).unwrap();
    std::fs::write(dir.join("ref.fasta.fai"), unescape(field(&text, "fai"))).unwrap();
    let mut source = ReferenceFileSource::open(&fasta).expect("the fixture opens");

    let mut compared = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("result\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').unwrap_or((rest, ""));

        let Some((contig, start, end, leading, trailing, accessor)) = parse_label(label) else {
            // The named cases at the end of the golden are checked separately below.
            continue;
        };
        let interval = SimpleInterval::new(&contig, start, end);
        let ours = context_answer(
            &mut source,
            interval,
            leading,
            trailing,
            &accessor,
            start,
            end,
        );
        assert_eq!(ours, expected, "{label}");
        compared += 1;
    }

    // The cases whose labels are names rather than coordinates.
    let named = |label: &str| -> String {
        text.lines()
            .find_map(|l| {
                l.strip_prefix("result\t")
                    .and_then(|r| r.strip_prefix(label))
                    .and_then(|r| r.strip_prefix('\t'))
            })
            .unwrap_or_else(|| panic!("no {label} row"))
            .to_string()
    };

    // A negative window is a GATKException, not a clamp to zero.
    assert_eq!(
        ReferenceContext::new(&source, SimpleInterval::new("chr1", 5, 5), -1, 0)
            .err()
            .map_or("no throw", |_| "E"),
        named("negativeWindow")
    );
    // A contig the reference does not carry is a UserException from trimToContigLength.
    assert_eq!(
        ReferenceContext::new(&source, SimpleInterval::new("chr9", 5, 5), 3, 3)
            .err()
            .map_or("no throw", |_| "E"),
        named("unknownContig")
    );
    // The explicit window must contain the interval.
    assert_eq!(
        ReferenceContext::with_window(
            &source,
            SimpleInterval::new("chr1", 5, 10),
            SimpleInterval::new("chr1", 6, 7),
        )
        .err()
        .map_or("no throw", |_| "E"),
        named("windowInsideInterval")
    );
    // No interval at all: empty answers rather than failures.
    let mut empty = ReferenceContext::new(&source, None, 0, 0).expect("a contextless context");
    assert_eq!(
        format!(
            "bases={} lead={} trail={} backing={}",
            bytes(empty.bases(&mut source).expect("empty bases")),
            empty.num_window_leading_bases(),
            empty.num_window_trailing_bases(),
            false
        ),
        named("noInterval")
    );

    // The copy constructor carries the window *sizes*, and at a contig edge those are the cropped
    // ones: a context built at position 1 with a lead of 10 has a lead of 0, so the context
    // derived from it is asymmetric.
    let edge = ReferenceContext::new(&source, SimpleInterval::new("chr1", 1, 1), 10, 10)
        .expect("the edge context");
    let mut moved = edge
        .with_interval(&source, SimpleInterval::new("chr1", 20, 20))
        .expect("the moved context");
    let window = moved.window().expect("a window").clone();
    assert_eq!(
        format!(
            "edgeLead={} edgeTrail={} movedWindow={}:{}-{} bases={}",
            edge.num_window_leading_bases(),
            edge.num_window_trailing_bases(),
            window.contig,
            window.start,
            window.end,
            bytes(moved.bases(&mut source).expect("moved bases")),
        ),
        named("copyFromEdge")
    );

    std::fs::remove_dir_all(&dir).ok();
    assert!(compared > 0, "the golden carries no contexts");
    println!("{compared} reference-context answers, all identical");
}

/// One accessor on one freshly built context, rendered the way the harness renders it.
fn context_answer(
    source: &mut ReferenceFileSource,
    interval: Option<SimpleInterval>,
    leading: i32,
    trailing: i32,
    accessor: &str,
    start: i32,
    end: i32,
) -> String {
    let contig = interval
        .as_ref()
        .map(|i| i.contig.clone())
        .unwrap_or_default();
    let build = |source: &ReferenceFileSource| {
        ReferenceContext::new(source, interval.clone(), leading, trailing)
    };

    let threw = "E".to_string();
    match accessor {
        "window" => match build(source) {
            Err(_) => threw,
            Ok(context) => match context.window() {
                None => "null".to_string(),
                Some(window) => format!("{}:{}-{}", window.contig, window.start, window.end),
            },
        },
        "lead" => build(source).map_or(threw, |c| c.num_window_leading_bases().to_string()),
        "trail" => build(source).map_or(threw, |c| c.num_window_trailing_bases().to_string()),
        "bases" => match build(source) {
            Err(_) => threw,
            Ok(mut context) => context.bases(source).map_or(threw, bytes),
        },
        "forward" => match build(source) {
            Err(_) => threw,
            Ok(mut context) => context.forward_bases(source).map_or(threw, bytes),
        },
        "base" => match build(source) {
            Err(_) => threw,
            Ok(mut context) => context
                .base(source)
                .map_or(threw, |b| (b as char).to_string()),
        },
        "expand5" => match build(source) {
            Err(_) => threw,
            Ok(context) => context.bases_expanded(source, 5, 5).map_or(threw, bytes),
        },
        "expand0" => match build(source) {
            Err(_) => threw,
            Ok(context) => context.bases_expanded(source, 0, 0).map_or(threw, bytes),
        },
        "kmer3" | "kmer20" => {
            let each_side = if accessor == "kmer3" { 3 } else { 20 };
            match build(source) {
                Err(_) => threw,
                Ok(context) => match context.kmer_around(source, start, each_side) {
                    Err(_) => threw,
                    // The reference returns null rather than a shorter kmer at a contig edge, and
                    // String.valueOf(null) is the text "null".
                    Ok(None) => "null".to_string(),
                    Ok(Some(kmer)) => bytes(kmer),
                },
            }
        }
        // The window given as coordinates rather than as two counts. Building the argument is
        // itself where the reference throws when the padding runs off the front of the contig,
        // which is why the interval constructor is checked.
        "explicit" => {
            let window = SimpleInterval::new(&contig, start - leading, end + trailing);
            if window.is_none() {
                return threw;
            }
            match ReferenceContext::with_window(source, interval.clone(), window) {
                Err(_) => threw,
                Ok(mut context) => context.bases(source).map_or(threw, bytes),
            }
        }
        other => panic!("{other} is in the golden but not ported"),
    }
}
