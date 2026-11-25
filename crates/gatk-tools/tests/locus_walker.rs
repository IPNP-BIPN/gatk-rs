//! Conformance for the `LocusWalker` traversal against GATK 4.6.2.0.
//!
//! One row per `apply` call: the locus, its depth, its bases and the reference base under it.
//!
//! The rows that earn the suite are the pair `gap` and `gap-emptyloci`, which run the same
//! interval both ways: with the default the walker makes six calls, one per covered position, and
//! with `emitEmptyLoci` overridden it makes fifty-one, the uncovered ones carrying a real context
//! with zero depth and the reference base still present.
//!
//! Two cases are named rather than compared, and for opposite reasons:
//!
//!  * `negative-depth` is a refusal on both sides, so what is asserted is the refusal;
//!  * `depth-1` is a run the port refuses and the reference completes. Its rows are **identical**
//!    to the undownsampled run, because this fixture never exceeds one read per sample at any
//!    locus, so the run does not measure downsampling at all. That is stated here rather than
//!    counted as coverage, and a probe that actually downsamples belongs with the `java.util.Random`
//!    port that the reservoir needs.

use gatk_corpus as corpus;
use gatk_engine::interval::{self, SimpleInterval};
use gatk_engine::reads::ReadsDataSource;
use gatk_engine::reference::ReferenceFileSource;
use gatk_tools::locus_walker::{self, LocusWalkerError, Options};
use htsjdk_bam::record::BamRecord;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/locus_walker.txt.gz"),
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

/// The arguments each labelled run was given: intervals, empty loci, reference, depth cap.
fn configuration(label: &str) -> (Option<Vec<&'static str>>, bool, bool, i32) {
    match label {
        "all" => (None, false, true, 0),
        "chr1:100-130" => (Some(vec!["chr1:100-130"]), false, true, 0),
        "gap" => (Some(vec!["chr1:20-70"]), false, true, 0),
        "gap-emptyloci" => (Some(vec!["chr1:20-70"]), true, true, 0),
        "chr2" => (Some(vec!["chr2"]), false, true, 0),
        "all-noref" => (None, false, false, 0),
        "negative-depth" => (None, false, true, -1),
        "depth-1" => (None, false, true, 1),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// The run whose rows the reference produced but the port refuses, with the reason.
// `depth-1` used to be asserted as a REFUSAL: the reference completed the run and the port did
// not, because the locus iterator had no downsampler. It has one now (#1102), so the case is
// compared like every other label. It still thins nothing -- the fixture's depth is one, so a
// target of one is already met, which is why the golden's rows for it equal the undownsampled
// run's -- and a corpus where a target below the depth is observed is #1135.

#[test]
fn every_apply_call_matches_the_reference() {
    let text = golden();

    // The ReadWalker fixture travels in the ReadWalker golden; this suite shares it, so it is read
    // from there rather than duplicated.
    let read_walker = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_walker.txt.gz"),
    );
    let dir = std::env::temp_dir().join(format!("gatk-rs-locuswalker-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bam = dir.join("reads.bam");
    let bai = dir.join("reads.bai");
    let fasta = dir.join("ref.fasta");
    std::fs::write(&bam, corpus::decode_base64(field(&read_walker, "bam"))).unwrap();
    std::fs::write(&bai, corpus::decode_base64(field(&read_walker, "bai"))).unwrap();
    std::fs::write(&fasta, unescape(field(&read_walker, "fasta"))).unwrap();
    std::fs::write(
        dir.join("ref.fasta.fai"),
        unescape(field(&read_walker, "fai")),
    )
    .unwrap();

    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();
    let reads = source.iter_all().expect("the fixture reads");

    let mut labels: Vec<String> = Vec::new();
    let mut rows: std::collections::HashMap<String, Vec<String>> = Default::default();
    let mut outcomes: std::collections::HashMap<String, String> = Default::default();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("apply\t") {
            let mut parts = rest.splitn(3, '\t');
            let label = parts.next().expect("a label").to_string();
            let _index = parts.next();
            rows.entry(label)
                .or_default()
                .push(parts.next().unwrap_or("").to_string());
        } else if let Some(rest) = line.strip_prefix("summary\t") {
            let (label, outcome) = rest.split_once('\t').expect("a label and an outcome");
            labels.push(label.to_string());
            outcomes.insert(label.to_string(), outcome.to_string());
        }
    }
    assert!(!labels.is_empty(), "the golden carries no summary rows");

    let mut compared = 0;
    let mut refused = 0;
    for label in &labels {
        let (interval_strings, emit_empty_loci, has_reference, max_depth) = configuration(label);
        let intervals: Option<Vec<SimpleInterval>> = interval_strings.map(|strings| {
            strings
                .iter()
                .map(|text| interval::parse_interval(text, &header).expect("a parsable interval"))
                .collect()
        });

        let options = Options {
            emit_empty_loci,
            max_depth_per_sample: max_depth,
            ..Options::default()
        };
        let mut reference = if has_reference {
            Some(ReferenceFileSource::open(&fasta).expect("the fixture reference opens"))
        } else {
            None
        };
        let header_for_filter = header.clone();
        let filter = move |read: &BamRecord| locus_walker::default_filter(&header_for_filter)(read);

        let result = locus_walker::traverse(
            &reads,
            &header,
            reference.as_mut(),
            intervals.as_deref(),
            options,
            &filter,
        );

        let outcome = &outcomes[label];
        if outcome.starts_with("E:") {
            assert!(
                matches!(result, Err(LocusWalkerError::NegativeMaxDepth(_))),
                "{label}: the reference raised {outcome}, the port gave {:?}",
                result.as_ref().map(|a| a.len())
            );
            refused += 1;
            continue;
        }

        let applied = result.unwrap_or_else(|e| panic!("{label}: the port refused: {e:?}"));
        let expected = rows.get(label).cloned().unwrap_or_default();
        assert_eq!(applied.len(), expected.len(), "{label}: apply count");
        for (index, mut call) in applied.into_iter().enumerate() {
            let bases: String = call
                .context
                .pileup
                .bases()
                .iter()
                .map(|b| *b as char)
                .collect();
            let reference_bases = match reference.as_mut() {
                Some(source) => call.reference.bases(source).expect("bases"),
                None => Vec::new(),
            };
            let ours = format!(
                "{}:{}|{}|{}|{}",
                call.context.contig,
                call.context.position,
                call.context.pileup.size(),
                if bases.is_empty() {
                    "-".to_string()
                } else {
                    bases
                },
                if reference_bases.is_empty() {
                    "-".to_string()
                } else {
                    String::from_utf8(reference_bases).expect("ASCII bases")
                },
            );
            assert_eq!(ours, expected[index], "{label} apply {index}");
            compared += 1;
        }
    }

    std::fs::remove_dir_all(&dir).ok();
    println!(
        "{compared} apply calls over {} traversals, {refused} refused as the reference refuses or \
         as the port declares",
        labels.len()
    );
}

/// `for_each_streamed` hands `apply` the loci `traverse` collects, the same pileups in the same
/// order, however small its segments, on every configuration it accepts.
#[test]
fn streamed_segments_are_the_whole_traversal() {
    let read_walker = corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/read_walker.txt.gz"),
    );
    let dir = std::env::temp_dir().join(format!("gatk-rs-locusstream-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let (bam, bai, fasta) = (dir.join("r.bam"), dir.join("r.bai"), dir.join("ref.fasta"));
    std::fs::write(&bam, corpus::decode_base64(field(&read_walker, "bam"))).unwrap();
    std::fs::write(&bai, corpus::decode_base64(field(&read_walker, "bai"))).unwrap();
    std::fs::write(&fasta, unescape(field(&read_walker, "fasta"))).unwrap();
    std::fs::write(
        dir.join("ref.fasta.fai"),
        unescape(field(&read_walker, "fai")),
    )
    .unwrap();
    let source = ReadsDataSource::open(&bam, &bai).expect("the fixture opens");
    let header = source.header().clone();
    let reads = source.iter_all().expect("the fixture reads");
    let row = |call: &locus_walker::Applied<'_>| {
        let names: Vec<&str> = call
            .context
            .pileup
            .elements
            .iter()
            .map(|element| element.read.read_name.as_str())
            .collect();
        format!(
            "{}:{}|{:?}|{:?}",
            call.context.contig,
            call.context.position,
            call.context.pileup.bases(),
            names
        )
    };

    let mut compared = 0;
    let configurations = [None, Some("chr1:100-130"), Some("chr1:20-70"), Some("chr2")]
        .into_iter()
        .flat_map(|bounds| [(bounds, false), (bounds, true)]);
    for (bounds, emit_empty_loci) in configurations {
        let intervals: Option<Vec<SimpleInterval>> =
            bounds.map(|text| vec![interval::parse_interval(text, &header).unwrap()]);
        let filter = |read: &BamRecord| locus_walker::default_filter(&header)(read);
        let options = Options {
            emit_empty_loci,
            ..Options::default()
        };
        let whole: Vec<String> = locus_walker::traverse(
            &reads,
            &header,
            None,
            intervals.as_deref(),
            options,
            &filter,
        )
        .unwrap()
        .iter()
        .map(row)
        .collect();
        assert!(!whole.is_empty());
        for segment in [1, 2, 3, 5, 1000] {
            let mut streamed = Vec::new();
            locus_walker::for_each_streamed(
                reads.iter().cloned().map(Ok::<_, ()>),
                &header,
                None,
                intervals.as_deref(),
                options,
                &filter,
                segment,
                |call| {
                    streamed.push(row(&call));
                    true
                },
            )
            .unwrap();
            assert_eq!(
                streamed, whole,
                "{bounds:?}, empty loci {emit_empty_loci}, segments of {segment}"
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 40);

    // A read out of order stops it, for the caller to answer from `traverse`.
    let mut unsorted = reads.clone();
    unsorted.swap(0, reads.len() - 1);
    let stopped = locus_walker::for_each_streamed(
        unsorted.into_iter().map(Ok::<_, ()>),
        &header,
        None,
        None,
        Options::default(),
        &|read: &BamRecord| locus_walker::default_filter(&header)(read),
        2,
        |_| true,
    );
    assert!(
        matches!(stopped, Err(locus_walker::Stopped::Unsorted)),
        "{stopped:?}"
    );
}

/// The same on reads that overlap across every cut, from two samples, through deletions,
/// insertions, clips and skips, with reads sharing a start.
#[test]
fn streamed_segments_carry_the_reads_that_reach_the_cut() {
    let mut sam = String::from(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n\
         @RG\tID:a\tSM:s1\n@RG\tID:b\tSM:s2\n",
    );
    let cigars = ["10M", "4M2D6M", "3S7M", "5M3I2M", "4M20N6M", "2M1D8M"];
    for index in 0..40 {
        let contig = if index < 30 { "chr1" } else { "chr2" };
        let start = 100 + (index / 2) * 3;
        let cigar = cigars[index % cigars.len()];
        // The bases the cigar reads: its M, I and S lengths.
        let (mut length, mut number) = (0usize, 0usize);
        for c in cigar.chars() {
            match c.to_digit(10) {
                Some(digit) => number = number * 10 + digit as usize,
                None => {
                    if "MIS".contains(c) {
                        length += number;
                    }
                    number = 0;
                }
            }
        }
        let group = if index % 3 == 0 { "b" } else { "a" };
        sam.push_str(&format!(
            "r{index}\t0\t{contig}\t{start}\t60\t{cigar}\t*\t0\t0\t{}\t{}\tRG:Z:{group}\n",
            "ACGT".repeat(length).get(..length).unwrap(),
            "I".repeat(length),
        ));
    }
    let (header, reads) = htsjdk_bam::sam_file::read_sam(&sam).expect("the SAM parses");
    let row = |call: &locus_walker::Applied<'_>| {
        let names: Vec<&str> = call
            .context
            .pileup
            .elements
            .iter()
            .map(|element| element.read.read_name.as_str())
            .collect();
        format!(
            "{}:{}|{names:?}",
            call.context.contig, call.context.position
        )
    };
    let filter = |read: &BamRecord| locus_walker::default_filter(&header)(read);
    let whole: Vec<String> =
        locus_walker::traverse(&reads, &header, None, None, Options::default(), &filter)
            .unwrap()
            .iter()
            .map(row)
            .collect();
    assert!(
        whole.iter().any(|row| row.matches('r').count() > 3),
        "reads overlap"
    );
    for segment in 1..=6 {
        let mut streamed = Vec::new();
        locus_walker::for_each_streamed(
            reads.iter().cloned().map(Ok::<_, ()>),
            &header,
            None,
            None,
            Options::default(),
            &filter,
            segment,
            |call| {
                streamed.push(row(&call));
                true
            },
        )
        .unwrap();
        assert_eq!(streamed, whole, "segments of {segment}");
    }

    // Empty loci over a request that starts before the reads, spans a gap and ends past them.
    let requested = vec![
        interval::parse_interval("chr1:90-125", &header).unwrap(),
        interval::parse_interval("chr1:140-160", &header).unwrap(),
        interval::parse_interval("chr2:140-170", &header).unwrap(),
    ];
    let options = Options {
        emit_empty_loci: true,
        ..Options::default()
    };
    let whole: Vec<String> =
        locus_walker::traverse(&reads, &header, None, Some(&requested), options, &filter)
            .unwrap()
            .iter()
            .map(row)
            .collect();
    assert_eq!(whole.len(), 36 + 21 + 31);
    for segment in 1..=6 {
        let mut streamed = Vec::new();
        locus_walker::for_each_streamed(
            reads.iter().cloned().map(Ok::<_, ()>),
            &header,
            None,
            Some(&requested),
            options,
            &filter,
            segment,
            |call| {
                streamed.push(row(&call));
                true
            },
        )
        .unwrap();
        assert_eq!(streamed, whole, "empty loci, segments of {segment}");
    }
}
