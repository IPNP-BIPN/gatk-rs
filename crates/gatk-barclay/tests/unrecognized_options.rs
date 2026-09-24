//! The token an unrecognized option is reported with, and what jopt-simple counts as an option.
//!
//! Measured on 2026-09-24 against `gatk-rs-oracle:4.6.2.0`, `CountVariants -V va.vcf <shape>`:
//!
//! ```text
//! -Q foo     -Q is not a recognized option
//! -Qfoo      -Qfoo is not a recognized option
//! -foo       -foo is not a recognized option
//! -5         -5 is not a recognized option
//! -1.5       -1.5 is not a recognized option
//! --foo      foo is not a recognized option
//! --Q foo    Q is not a recognized option
//! --         (accepted: the option terminator)
//! -- -5      Positional arguments were provided ',-5}' ...
//! ```
//!
//! A long option goes through jopt-simple's `handleLongOptionToken`, which reports the key without
//! its dashes. A short one goes through `BarclayOptionParser.handleShortOptionToken`, which throws
//! `new UnrecognizedOptionException(candidate)`: the whole token, dash included. A negative number
//! is a short option token to `ParserRules`, and `--` ends the options.

use gatk_barclay::{Annotation, Definition, Parser, Value, ValueClass};

/// `-V/--variant`, the one argument the measured command line gives.
fn parser() -> Parser {
    Parser::new(vec![Definition::new(
        Annotation {
            full_name: "variant",
            short_name: "V",
            doc: "the input",
            optional: true,
            ..Annotation::default()
        },
        "variant",
        ValueClass::Text,
        false,
        false,
        Value::Null,
    )])
}

fn refusal(argv: &[&str]) -> String {
    parser()
        .parse_arguments(argv)
        .expect_err("the command line is refused")
        .message
}

#[test]
fn a_short_option_is_reported_with_its_dash() {
    for (shape, token) in [
        (&["-V", "va.vcf", "-Q", "foo"][..], "-Q"),
        (&["-X"][..], "-X"),
        (&["-Qfoo"][..], "-Qfoo"),
        (&["-foo"][..], "-foo"),
    ] {
        assert_eq!(
            refusal(shape),
            format!("{token} is not a recognized option"),
            "{shape:?}"
        );
    }
}

#[test]
fn a_long_option_is_reported_without_its_dashes() {
    assert_eq!(refusal(&["--foo"]), "foo is not a recognized option");
    assert_eq!(refusal(&["--foo", "bar"]), "foo is not a recognized option");
    assert_eq!(refusal(&["--Q", "foo"]), "Q is not a recognized option");
}

#[test]
fn a_negative_number_is_a_short_option() {
    assert_eq!(refusal(&["-5"]), "-5 is not a recognized option");
    assert_eq!(refusal(&["-1.5"]), "-1.5 is not a recognized option");
}

#[test]
fn a_double_dash_ends_the_options() {
    parser()
        .parse_arguments(&["-V", "va.vcf", "--"])
        .expect("a trailing terminator is accepted");
    let positional = refusal(&["--", "-5"]);
    assert!(
        positional.contains("Positional arguments were provided ',-5}'"),
        "{positional}"
    );
}
