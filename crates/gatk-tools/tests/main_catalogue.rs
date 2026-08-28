//! Conformance for `Main`'s tool catalogue against GATK 4.6.2.0.
//!
//! Golden from `tools/readfilter-conformance/MainCatalogueDump.java`, which performs the scan the
//! way `Main.extractCommandLineProgram` does and then asks the suggestion search of the whole set.
//!
//! # What this suite is for
//!
//!  * **the catalogue being what the port carries, name for name**;
//!  * **it being bigger than the documented tool list**;
//!  * **the two classes excluded by name being absent**;
//!  * **a deprecated tool being absent, which is why the registry has to answer first**;
//!  * **the Spark tools and the Picard tools being present**;
//!  * **and the search of #947 answering the same over three hundred tools as over five.**

use gatk_corpus as corpus;
use gatk_tools::main_catalogue::{resolves, CATALOGUE};
use gatk_tools::main_dispatch::{
    distance, suggested_alternate_command, tool_deprecation_info, unknown_command_message,
    HELP_SIMILARITY_FLOOR,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/main_catalogue.txt.gz"),
    )
}

fn unescape(text: &str) -> String {
    text.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// The catalogue the golden printed, reassembled from its hundred-name lines.
fn catalogue(text: &str) -> Vec<String> {
    let mut names = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("catalogue\t") {
            let payload = rest.split_once('\t').expect("a chunk").1;
            names.extend(payload.split(',').map(str::to_string));
        }
    }
    names
}

/// The names one query suggested, sorted, as the dump emitted them.
fn suggests(text: &str, query: &str) -> Option<Vec<String>> {
    field(text, "suggests", query).map(|names| names.split(',').map(str::to_string).collect())
}

/// The port's catalogue is the reference's, name for name.
#[test]
fn the_catalogue_is_the_references() {
    let text = golden();
    let theirs = catalogue(&text);
    assert_eq!(theirs.len(), 331);
    assert_eq!(field(&text, "count", "catalogue").as_deref(), Some("331"));
    assert_eq!(CATALOGUE.to_vec(), theirs);
    // Sorted, which is what makes the lookup a binary search.
    let mut sorted = theirs.clone();
    sorted.sort();
    assert_eq!(theirs, sorted);
}

/// It is bigger than the documented tool list: the CLI reports 311 tools and the scan finds 331.
#[test]
fn the_catalogue_is_bigger_than_the_documentation() {
    assert_eq!(CATALOGUE.len(), 331);
    assert!(CATALOGUE.len() > 311);
}

/// The names the golden asked about are present or absent as it says.
#[test]
fn the_catalogue_holds_what_the_golden_says_it_holds() {
    let text = golden();
    let holds = field(&text, "holds", "names").expect("the holds row");
    let mut checked = 0;
    for part in holds.split(',') {
        let (name, present) = part.split_once('=').expect("a name and a verdict");
        let present: bool = present.parse().expect("a boolean");
        assert_eq!(resolves(name), present, "{name}");
        checked += 1;
    }
    assert_eq!(checked, 9);
    // The two the reference excludes by name, whatever their annotation says.
    assert!(!resolves("PicardCommandLineProgramExecutor"));
    assert!(!resolves("CommandLineArgumentValidator"));
    // The Spark tools and the Picard tools are both in it.
    assert!(resolves("PrintReadsSpark"));
    assert!(resolves("SortVcf"));
    assert!(resolves("MarkDuplicates"));
}

/// A deprecated tool is not in the catalogue at all, which is the only reason its notice is seen:
/// the name falls through to the search, and the registry answers before the search does.
#[test]
fn a_deprecated_tool_is_not_in_the_catalogue() {
    let text = golden();
    for name in ["CNNScoreVariants", "IndelRealigner"] {
        assert!(!resolves(name), "{name}");
        assert!(tool_deprecation_info(name).is_some(), "{name}");
    }
    let message = field(&text, "message", "IndelRealigner").expect("its notice");
    assert_eq!(
        message,
        tool_deprecation_info("IndelRealigner").expect("a notice")
    );
    // The registry is the only thing that has anything to say about it: over this catalogue the
    // search finds no neighbour under the floor either, so without the registry the message would
    // have been the bare `is not a valid command` line.
    let names: Vec<String> = CATALOGUE.iter().map(|name| name.to_string()).collect();
    assert_eq!(
        suggested_alternate_command(&names, "IndelRealigner").expect("no refusal"),
        "'IndelRealigner' is not a valid command.\n"
    );
}

/// A name that resolves is refused rather than searched, over the whole catalogue as over five.
#[test]
fn a_name_that_resolves_is_refused() {
    let text = golden();
    let names: Vec<String> = CATALOGUE.iter().map(|name| name.to_string()).collect();
    assert!(resolves("SortVcf"));
    assert_eq!(
        unknown_command_message(&names, "SortVcf"),
        Err("Command matches: SortVcf".to_string())
    );
    assert_eq!(
        field(&text, "error", "SortVcf").as_deref(),
        Some("java.lang.RuntimeException:Command matches: SortVcf")
    );
}

/// The search answers the same over three hundred tools as it did over five: every query finds
/// exactly the tools the golden says it finds.
#[test]
fn the_search_finds_the_same_tools() {
    let text = golden();
    let names: Vec<String> = CATALOGUE.iter().map(|name| name.to_string()).collect();
    let mut compared = 0;
    for query in [
        "PrintRead",
        "PrintReadz",
        "HaplotypeCallr",
        "MarkDuplicate",
        "CollectQuality",
        "PathSeq",
        "Fingerprint",
    ] {
        let expected = suggests(&text, query).unwrap_or_else(|| panic!("{query}"));
        let message = suggested_alternate_command(&names, query).expect("no refusal");
        let mut ours: Vec<String> = names
            .iter()
            .filter(|name| {
                distance(query, name).expect("a distance")
                    == names
                        .iter()
                        .map(|other| distance(query, other).expect("a distance"))
                        .min()
                        .expect("a best")
            })
            .cloned()
            .collect();
        ours.sort();
        ours.dedup();
        assert_eq!(ours, expected, "{query}");
        assert!(message.contains("Did you mean"), "{query}");
        compared += 1;
    }
    assert_eq!(compared, 7);
}

/// A name far from everything finds nothing, even against three hundred tools.
#[test]
fn a_name_far_from_everything_finds_nothing() {
    let text = golden();
    let names: Vec<String> = CATALOGUE.iter().map(|name| name.to_string()).collect();
    let query = "Zzzzzzzzzzzzzzzzzzzz";
    assert_eq!(
        field(&text, "message", query).as_deref(),
        Some("'Zzzzzzzzzzzzzzzzzzzz' is not a valid command.\n")
    );
    let message = suggested_alternate_command(&names, query).expect("no refusal");
    assert_eq!(message, "'Zzzzzzzzzzzzzzzzzzzz' is not a valid command.\n");
    // Nothing is under the floor.
    let best = names
        .iter()
        .map(|name| distance(query, name).expect("a distance"))
        .min()
        .expect("a best");
    assert!(best >= HELP_SIMILARITY_FLOOR, "{best}");
}
