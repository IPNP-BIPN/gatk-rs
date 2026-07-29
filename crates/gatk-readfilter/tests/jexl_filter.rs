//! Conformance for `JexlExpressionReadTagValueFilter` and the commons-jexl 2.1.1 arithmetic under
//! it, against GATK 4.6.2.0.
//!
//! Two rows per (expression, read): what the engine evaluated to, and what the filter made of it.
//! They are separate because the filter's `!v.equals(Boolean.TRUE)` turns a null into a
//! `NullPointerException` and a non-boolean into a quiet false, and a port could get the second
//! right while getting the first wrong.
//!
//! The golden corrected the port on one point. An absent tag does **not** evaluate to null: the
//! interpreter throws `JexlException.Variable`, because `setLenient(false)` makes it strict, so
//! even `ZZ == null` throws rather than answering true.

use gatk_corpus as corpus;
use gatk_engine::jexl::{self, Context, JexlError, Value};
use gatk_readfilter::jexl_filter::{self, Decision};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/jexl_filter.txt.gz"),
    )
}

/// The context each labelled read presents, rebuilt from the `attr` rows of the golden.
///
/// Rebuilt rather than reconstructed from the tag values: `getAttributeAsString` is the adapter's
/// rendering, and comparing against the port's own rendering of the same tag would compare two
/// reimplementations instead of the reference.
fn contexts(text: &str) -> Vec<(String, Context)> {
    let mut order: Vec<String> = Vec::new();
    let mut contexts: std::collections::HashMap<String, Context> = Default::default();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("attr\t") else {
            continue;
        };
        let (label, pair) = rest.split_once('\t').expect("a label and a pair");
        let (tag, value) = pair.split_once('=').expect("tag=value");
        if !contexts.contains_key(label) {
            order.push(label.to_string());
        }
        let context = contexts.entry(label.to_string()).or_default();
        // `null` is the row for a tag the read does not carry, and an absent key is exactly what
        // makes the identifier throw.
        if value != "null" {
            context.insert(tag.to_string(), value.to_string());
        }
    }
    order
        .into_iter()
        .map(|label| {
            let context = contexts.remove(&label).expect("a context");
            (label, context)
        })
        .collect()
}

/// How the reference printed an evaluated value: `ok:<simple class>:<toString>`, or `ok:null`.
fn render(value: &Value) -> String {
    match value {
        Value::Null => "ok:null".to_string(),
        Value::Bool(b) => format!("ok:Boolean:{b}"),
        Value::Int(i) => format!("ok:Integer:{i}"),
        Value::Long(l) => format!("ok:Long:{l}"),
        Value::Float(f) => format!("ok:Float:{}", render_float(*f)),
        Value::Double(d) => format!("ok:Double:{}", render_double(*d)),
        Value::Str(s) => format!("ok:String:{s}"),
    }
}

fn render_float(value: f32) -> String {
    if value == value.trunc() && value.is_finite() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

fn render_double(value: f64) -> String {
    if value == value.trunc() && value.is_finite() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

/// The exception class the reference raised, for the refusal the port produced.
fn class_of(error: &JexlError) -> String {
    match error {
        JexlError::UndefinedVariable(_) => {
            "E:org.apache.commons.jexl2.JexlException$Variable".to_string()
        }
        JexlError::NumberFormat(_) => "E:java.lang.NumberFormatException".to_string(),
        JexlError::Arithmetic(_) => "E:java.lang.ArithmeticException".to_string(),
        JexlError::Parse(message) => panic!("the port failed to parse: {message}"),
        JexlError::Unsupported(message) => panic!("the port refuses: {message}"),
    }
}

#[test]
fn every_expression_answers_what_the_reference_answers() {
    let text = golden();
    let contexts = contexts(&text);
    assert!(!contexts.is_empty(), "the golden carries no attr rows");

    let mut evaluated = 0;
    let mut filtered = 0;
    for line in text.lines() {
        let (kind, rest) = match line.split_once('\t') {
            Some(("eval", rest)) => ("eval", rest),
            Some(("filter", rest)) => ("filter", rest),
            _ => continue,
        };
        let mut parts = rest.split('\t');
        let expression = parts.next().expect("an expression");
        let label = parts.next().expect("a read label");
        let expected = parts.next().expect("an outcome");

        let context = contexts
            .iter()
            .find(|(l, _)| l == label)
            .map(|(_, c)| c)
            .unwrap_or_else(|| panic!("{label} has no attr rows"));

        let compiled = jexl::create_expression(expression)
            .unwrap_or_else(|e| panic!("{expression}: the port did not parse it: {e:?}"));

        if kind == "eval" {
            let ours = match compiled.evaluate(context) {
                Ok(value) => render(&value),
                Err(error) => class_of(&error),
            };
            assert_eq!(ours, expected, "eval {expression} on {label}");
            evaluated += 1;
            continue;
        }

        // The filter's own answer. `Decision::NullResult` is the NullPointerException that
        // `!v.equals(Boolean.TRUE)` raises on a null, which no absent tag can reach.
        let names: Vec<String> = context.keys().cloned().collect();
        let ours = match jexl_filter::test_context(context, &[compiled], &names) {
            Decision::Keep => "true".to_string(),
            Decision::Drop => "false".to_string(),
            Decision::NullResult => "E:java.lang.NullPointerException".to_string(),
            Decision::Failed(error) => class_of(&error),
        };
        assert_eq!(ours, expected, "filter {expression} on {label}");
        filtered += 1;
    }

    println!(
        "{evaluated} expression evaluations and {filtered} filter decisions over {} reads",
        contexts.len()
    );
}
