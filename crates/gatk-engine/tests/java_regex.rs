//! Conformance for [`gatk_engine::java_regex`] against `java.util.regex` in GATK 4.6.2.0, compared
//! as whether every pattern compiles, what its refusal says, and what `find()` answers for every
//! input.
//!
//! Golden from `tools/readfilter-conformance/JavaRegexDump.java`.
//!
//! # What this suite is for
//!
//!  * **`find()` is a search**, which is why `-se s1` selects `xs10`;
//!  * **the anchors and `.` know four line terminators**, not one;
//!  * **the predefined classes are ASCII**, so an Arabic-Indic digit is not a digit;
//!  * **and a refusal's index has a rule per construct**, measured rather than inferred.
//!
//! # The boundary is asserted, not assumed
//!
//! The port is a subset: it refuses what it cannot represent instead of matching it differently.
//! Every pattern the reference accepts and the port refuses is listed in [`BOUNDARY`], so adding
//! one is a deliberate edit and losing one is a failure. Everything else must agree exactly.

use gatk_corpus as corpus;
use gatk_engine::java_regex::Pattern;

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/java_regex.txt.gz"),
    )
}

/// The reverse of the dump's `escape`, including its UTF-16 escapes.
fn unescape(text: &str) -> String {
    let mut units: Vec<u16> = Vec::new();
    let mut characters = text.chars().peekable();
    while let Some(character) = characters.next() {
        if character != '\\' {
            let mut buffer = [0u16; 2];
            units.extend_from_slice(character.encode_utf16(&mut buffer));
            continue;
        }
        match characters.next() {
            Some('t') => units.push(u16::from(b'\t')),
            Some('n') => units.push(u16::from(b'\n')),
            Some('r') => units.push(u16::from(b'\r')),
            Some('\\') => units.push(u16::from(b'\\')),
            Some('u') => {
                let digits: String = (0..4).filter_map(|_| characters.next()).collect();
                units.push(u16::from_str_radix(&digits, 16).expect("four hex digits"));
            }
            Some(other) => {
                units.push(u16::from(b'\\'));
                let mut buffer = [0u16; 2];
                units.extend_from_slice(other.encode_utf16(&mut buffer));
            }
            None => units.push(u16::from(b'\\')),
        }
    }
    String::from_utf16(&units).expect("valid UTF-16")
}

/// What the reference did with one pattern.
struct Row {
    pattern: String,
    /// `None` where it compiled; the description, index and message where it did not.
    refusal: Option<(String, i32, String)>,
    /// `(input, whether find() matched)`, empty where the pattern did not compile.
    finds: Vec<(String, bool)>,
}

fn rows(text: &str) -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::new();
    for line in text.lines().filter(|line| !line.starts_with('#')) {
        let field: Vec<&str> = line.split('\t').collect();
        match field[0] {
            "compile" => {
                let refusal = if field[2] == "error" {
                    Some((
                        unescape(field[3]),
                        field[4].parse().expect("an index"),
                        unescape(field[5]),
                    ))
                } else {
                    None
                };
                rows.push(Row {
                    pattern: unescape(field[1]),
                    refusal,
                    finds: Vec::new(),
                });
            }
            "find" => {
                let last = rows.last_mut().expect("a compile row first");
                assert_eq!(
                    last.pattern,
                    unescape(field[1]),
                    "a find without its compile"
                );
                last.finds.push((unescape(field[2]), field[3] == "true"));
            }
            other => panic!("unknown row {other}"),
        }
    }
    rows
}

/// The constructs the reference has and this port deliberately does not: refused at compile rather
/// than matched differently. Each was measured to compile in the reference.
const BOUNDARY: [&str; 17] = [
    "\\bNA",
    "(?i)na",
    "(?:ab)+",
    "(?=NA)N",
    "(?!NA)s",
    "(\\w)\\1",
    "(?>a*)b",
    "[a-z&&[^aeiou]]",
    "\\p{Alpha}+",
    "\\p{Lower}",
    "\\Qa.c\\E",
    "\\x41",
    "\\070",
    "\\cA",
    "^\\p{IsLatin}+$",
    "\\ANA",
    "NA\\z",
];

/// The rest of the boundary, kept separate only because the array above is already long.
const BOUNDARY_MORE: [&str; 7] = [
    "NA\\Z",
    "(?s).",
    "(?m)^x$",
    "\\h",
    "\\R",
    "(?<name>NA)\\k<name>",
    "[[a-c]&&[b-d]]",
];

fn is_boundary(pattern: &str) -> bool {
    BOUNDARY.contains(&pattern) || BOUNDARY_MORE.contains(&pattern)
}

#[test]
fn every_pattern_compiles_or_refuses_as_the_reference_does() {
    let text = golden();
    for row in rows(&text) {
        let ours = Pattern::compile(&row.pattern);
        match (&row.refusal, &ours) {
            (Some((description, index, message)), Err(error)) => {
                assert_eq!(
                    &error.description, description,
                    "description/{:?}",
                    row.pattern
                );
                assert_eq!(error.index, *index, "index/{:?}", row.pattern);
                assert_eq!(&error.message(), message, "message/{:?}", row.pattern);
            }
            (Some(_), Ok(_)) => panic!("{:?} compiles here and is refused there", row.pattern),
            (None, Err(error)) => assert!(
                is_boundary(&row.pattern),
                "{:?} compiles there and is refused here as {}, which is not a declared boundary",
                row.pattern,
                error.description
            ),
            (None, Ok(_)) => {}
        }
    }
}

#[test]
fn every_search_answers_what_the_reference_answered() {
    let text = golden();
    let mut compared = 0;
    for row in rows(&text) {
        let Ok(pattern) = Pattern::compile(&row.pattern) else {
            continue;
        };
        for (input, expected) in &row.finds {
            assert_eq!(
                pattern.find(input),
                *expected,
                "find({:?}, {:?})",
                row.pattern,
                input
            );
            compared += 1;
        }
    }
    // The boundary patterns contribute nothing, so this asserts the suite is still doing work.
    assert!(compared > 1200, "only {compared} searches compared");
}

/// The four behaviours a crate cannot be configured into, each read off the golden and then
/// asserted of the port.
#[test]
fn the_four_reasons_a_general_engine_cannot_stand_in() {
    let text = golden();
    let answer = |pattern: &str, input: &str| {
        rows(&text)
            .into_iter()
            .find(|row| row.pattern == pattern)
            .unwrap_or_else(|| panic!("no row for {pattern}"))
            .finds
            .into_iter()
            .find(|(candidate, _)| candidate == input)
            .unwrap_or_else(|| panic!("no input {input:?}"))
            .1
    };

    // `$` stands before a final line terminator, and takes CRLF as one.
    assert!(answer("^s1$", "s1\n"));
    assert!(answer("^s1$", "s1\r\n"));
    assert!(!answer("^s1$", "s1\nx"));
    assert!(Pattern::compile("^s1$").expect("compiled").find("s1\r\n"));
    assert!(!Pattern::compile("^s1$").expect("compiled").find("s1\nx"));

    // `.` refuses a carriage return as well as a newline.
    assert!(!answer("a.c", "a\rc"));
    assert!(!Pattern::compile("a.c").expect("compiled").find("a\rc"));

    // The predefined classes are ASCII.
    assert!(!answer("\\d", "\u{663}"));
    assert!(!answer("\\w", "\u{e9}"));
    assert!(!Pattern::compile("\\d").expect("compiled").find("\u{663}"));
    assert!(!Pattern::compile("\\w").expect("compiled").find("\u{e9}"));

    // A possessive quantifier gives nothing back.
    assert!(answer("^.*1$", "s1"));
    assert!(!answer("^.*+1$", "s1"));
    assert!(!Pattern::compile("^.*+1$").expect("compiled").find("s1"));
}

/// The index rule differs per construct, and two of them omit the caret line entirely.
#[test]
fn a_refusals_index_is_the_reference_s_index() {
    let text = golden();
    let refusal = |pattern: &str| {
        rows(&text)
            .into_iter()
            .find(|row| row.pattern == pattern)
            .unwrap_or_else(|| panic!("no row for {pattern}"))
            .refusal
            .unwrap_or_else(|| panic!("{pattern} was not refused"))
    };

    // A bare `)` has nothing to point at: index -1, no "near index" clause, no caret line.
    let (description, index, message) = refusal(")");
    assert_eq!(index, -1);
    assert_eq!(message, format!("{description}\n)"));
    let ours = Pattern::compile(")").expect_err("refused");
    assert_eq!(ours.index, -1);
    assert_eq!(ours.message(), message);

    // With a prefix it points at the character BEFORE the parenthesis.
    assert_eq!(refusal("abc)").1, 2);
    assert_eq!(Pattern::compile("abc)").expect_err("refused").index, 2);

    // An unclosed group points at the character after the parenthesis.
    assert_eq!(refusal("abc(").1, 4);
    assert_eq!(Pattern::compile("abc(").expect_err("refused").index, 4);

    // An unclosed class points at the bracket itself.
    assert_eq!(refusal("abc[").1, 3);
    assert_eq!(Pattern::compile("abc[").expect_err("refused").index, 3);

    // A repetition range whose low is above its high points at the closing brace.
    assert_eq!(refusal("ab{3,2}").1, 6);
    let ours = Pattern::compile("ab{3,2}").expect_err("refused");
    assert_eq!(ours.description, "Illegal repetition range");
    assert_eq!(ours.index, 6);

    // A trailing backslash points at the end, where there is no character to carry a caret.
    let (_, index, message) = refusal("a\\");
    assert_eq!(index, 2);
    assert!(!message.ends_with('^'));
    assert_eq!(
        Pattern::compile("a\\").expect_err("refused").message(),
        message
    );
}
