//! Conformance for `DumpTabixIndex` against GATK 4.6.2.0, compared as its printed text.
//!
//! Golden from `tools/readfilter-conformance/DumpTabixIndexDump.java`. The three `.tbi` files travel
//! in full, base64, so the port reads the same bytes the reference read.
//!
//! # What this suite is for
//!
//!  * **the bin ladder**, whose unit changes from M to K at 585;
//!  * **the linear index in 16K steps**;
//!  * **the format number is the file's**, 2 for a VCF and 65536 for a BED;
//!  * **the magic's fourth byte is a number**, so a file ending it with the letter `1` is refused;
//!  * **and a file that is not gzipped fails before the magic is checked**, which is why the gzip
//!    layer is the caller's.

use gatk_corpus as corpus;
use gatk_tools::dump_tabix_index::{dump_tabix_index, TabixError};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/dump_tabix_index.txt.gz"),
    )
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .filter_map(|line| line.strip_prefix(kind).and_then(|r| r.strip_prefix('\t')))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

/// The decompressed index of one label, which is what the tool's own routine reads.
fn index_bytes(text: &str, label: &str) -> Vec<u8> {
    let encoded = rows(text, "tbi")
        .into_iter()
        .find(|row| row[0] == label)
        .unwrap_or_else(|| panic!("no index {label}"))[1]
        .to_string();
    let compressed = corpus::decode_base64(&encoded);
    // A tabix index is BGZF, which is gzip with the block extra field: the same decompressor the
    // BAM reader uses reads it.
    htsjdk_bgzf::read::decompress_all(&compressed).expect("the index is BGZF")
}

#[test]
fn every_dump_is_the_reference() {
    let text = golden();
    let dumps = rows(&text, "dump");
    assert_eq!(dumps.len(), 3, "three indexes are dumped");

    for row in dumps {
        let label = row[0];
        let ours = dump_tabix_index(&index_bytes(&text, label))
            .unwrap_or_else(|error| panic!("{label} was refused: {}", error.message()));
        assert_eq!(ours, unescape(row[1]), "dump/{label}");
    }
}

/// The VCF and the BED do not declare the same format, and the columns follow.
#[test]
fn the_format_number_is_the_files_and_not_the_tools() {
    let text = golden();
    let header_of = |label: &str| -> String {
        let dump = rows(&text, "dump")
            .into_iter()
            .find(|row| row[0] == label)
            .expect("the run")[1]
            .to_string();
        unescape(&dump)
            .lines()
            .nth(1)
            .expect("the values row")
            .to_string()
    };
    assert_eq!(header_of("small"), "1\t2\t1\t2\t0\t#\t0");
    assert_eq!(header_of("regions"), "1\t65536\t1\t2\t3\t#\t0");
}

/// The magic's fourth byte is compared against the number, which the golden records as a refusal.
#[test]
fn the_wrong_magic_is_refused_with_the_references_message() {
    let text = golden();
    let expected = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "wrong-magic")
        .expect("the refusal")[1]
        .to_string();

    let letter = [b'T', b'B', b'I', b'1', 0, 0, 0, 0];
    let error = dump_tabix_index(&letter).expect_err("this magic is refused");
    assert_eq!(error, TabixError::WrongMagic);
    assert_eq!(
        format!(
            "org.broadinstitute.hellbender.exceptions.UserException:{}",
            error.message()
        ),
        expected
    );
}

/// The gzip layer belongs to the caller, which is what the second refusal is about.
#[test]
fn a_file_that_is_not_gzipped_fails_before_the_magic() {
    let text = golden();
    let expected = rows(&text, "error")
        .into_iter()
        .find(|row| row[0] == "not-gzipped")
        .expect("the refusal")[1]
        .to_string();
    // The reference never reaches its own check: the failure is java.util.zip's.
    assert_eq!(expected, "java.util.zip.ZipException:Not in GZIP format");

    // The same bytes, handed over decompressed, pass the magic and then run out of index: the
    // failure moves from the gzip layer to the reader's own end of stream.
    let plain = [b'T', b'B', b'I', 1, 0, 0, 0, 0];
    assert_eq!(dump_tabix_index(&plain).unwrap_err(), TabixError::PastEof);
}
