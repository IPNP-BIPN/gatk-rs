//! Conformance for the per-tool argument declarations against GATK 4.6.2.0.
//!
//! Golden from `tools/argument-conformance/ToolArgumentDeclarationDump.java`, which asks each tool
//! for its own parser and prints every named argument it reports. The module beside this test is
//! generated from that golden by `tools/declarations/generate.py`, which checks it against a
//! second, independent reading of the same reference on the way through: the inventory, taken
//! from the tool's usage text rather than from its reflection. This file holds the facts a
//! command line depends on.
//!
//! # What this suite is for
//!
//!  * **the counts being the reference's, and the two readings staying side by side**;
//!  * **a read walker and a variant walker not mirroring each other**;
//!  * **the parse decisions the golden holds following from the declarations**;
//!  * **the type being the UNDERLYING field's, which for a collection is its element class**;
//!  * **and the four other declarations a parser and a usage text are built from: the visibility
//!    flags, the bounds, the controlling plugin and the documentation.**

use gatk_corpus as corpus;
use gatk_tools::tool_declarations::{
    declarations, Declaration, APPLYBQSR, COUNTREADS, COUNTVARIANTS, GATHERVCFSCLOUD,
    INDEXFEATUREFILE, PRINTREADS, SELECTVARIANTS,
};

fn golden() -> String {
    corpus::read_golden(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/tool_declarations.txt.gz"),
    )
}

fn field(text: &str, kind: &str, tool: &str) -> Vec<String> {
    let prefix = format!("{kind}\t{tool}\t");
    text.lines()
        .filter(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
        .collect()
}

fn find<'a>(list: &'a [Declaration], name: &str) -> Option<&'a Declaration> {
    list.iter().find(|d| d.long_name == name)
}

/// The seven tools the golden and the module both carry.
const TOOLS: [(&str, &[Declaration]); 7] = [
    ("CountReads", COUNTREADS),
    ("CountVariants", COUNTVARIANTS),
    ("PrintReads", PRINTREADS),
    ("ApplyBQSR", APPLYBQSR),
    ("SelectVariants", SELECTVARIANTS),
    ("IndexFeatureFile", INDEXFEATUREFILE),
    ("GatherVcfsCloud", GATHERVCFSCLOUD),
];

/// The counts the golden holds are the ones the module carries.
#[test]
fn the_counts_are_the_reference_ones() {
    let text = golden();
    for (tool, list) in TOOLS {
        let counts = field(&text, "count", tool);
        let row = counts.first().unwrap_or_else(|| panic!("{tool}"));
        let tool_count: usize = row
            .split_once("tool=")
            .expect("a tool count")
            .1
            .parse()
            .expect("a number");
        assert_eq!(list.len(), tool_count, "{tool}");
        // And the instance-built parser's shorter list is named in the same row, so the two
        // readings stay side by side rather than one replacing the other.
        let instance_count: usize = row
            .split_once("instance=")
            .expect("an instance count")
            .1
            .split(' ')
            .next()
            .expect("a number")
            .parse()
            .expect("a number");
        // A walker's instance-built list is shorter; a non-walker's is the same list.
        assert!(instance_count <= tool_count, "{tool}");
        assert_eq!(
            field(&text, "only-on-the-tool", tool).len(),
            tool_count - instance_count
        );
    }
    assert!(declarations("CountReads").is_some());
    assert!(declarations("NoSuchTool").is_none());
    // The gap is the WALKER surface: every walker gains exactly thirty-two, and the two tools that
    // are no walkers gain nothing at all.
    for (tool, list) in TOOLS {
        let row = field(&text, "count", tool).remove(0);
        let instance: usize = row
            .split_once("instance=")
            .expect("a count")
            .1
            .split(' ')
            .next()
            .expect("a number")
            .parse()
            .expect("a number");
        let gap = list.len() - instance;
        let expected = if matches!(tool, "IndexFeatureFile" | "GatherVcfsCloud") {
            0
        } else {
            32
        };
        assert_eq!(gap, expected, "{tool}");
    }
}

/// A read walker and a variant walker do not mirror each other, and neither is a non-walker.
#[test]
fn the_two_archetypes_do_not_mirror_each_other() {
    let text = golden();
    // The read walker requires its input and has no `--variant` at all.
    assert!(find(COUNTREADS, "input").expect("input").required);
    assert!(find(COUNTREADS, "variant").is_none());
    // The variant walker requires `--variant` and takes an optional input.
    assert!(find(COUNTVARIANTS, "variant").expect("variant").required);
    assert!(!find(COUNTVARIANTS, "input").expect("input").required);
    // Which is exactly what the parser did with those command lines.
    let refusal = field(&text, "parse", "CountReads")
        .into_iter()
        .find(|row| row.starts_with("a-variant-argument\t"))
        .expect("the refusal");
    assert!(
        refusal.contains("-V is not a recognized option"),
        "{refusal}"
    );
    let accepted = field(&text, "parse", "CountVariants")
        .into_iter()
        .find(|row| row.starts_with("an-input\t"))
        .expect("the acceptance");
    assert!(accepted.ends_with("ok"), "{accepted}");
    // A tool that is no walker has neither: fourteen arguments, an input among them, and no
    // interval argument at all.
    assert_eq!(INDEXFEATUREFILE.len(), 14);
    assert!(find(INDEXFEATUREFILE, "input").expect("input").required);
    assert!(find(INDEXFEATUREFILE, "intervals").is_none());
    let refused = field(&text, "parse", "IndexFeatureFile")
        .into_iter()
        .find(|row| row.starts_with("an-interval\t"))
        .expect("the refusal");
    assert!(refused.contains("not a recognized option"), "{refused}");
}

/// The type is the underlying field's, which for a collection is its element class.
#[test]
fn the_type_is_the_underlying_fields() {
    let text = golden();
    // Every declaration in the module says what the golden's own line says.
    for (tool, list) in TOOLS {
        for row in field(&text, "def", tool) {
            let (index, body) = row.split_once('\t').expect("an index and a body");
            let fields: Vec<&str> = body.split('|').collect();
            let index: usize = index.parse().expect("a number");
            let declaration = &list[index];
            assert_eq!(declaration.long_name, fields[0], "{tool}");
            assert_eq!(declaration.type_name, fields[5], "{tool}");
            assert_eq!(declaration.primitive, fields[6] == "primitive", "{tool}");
            assert_eq!(declaration.flag, fields[7] == "flag", "{tool}");
            assert_eq!(declaration.hidden, fields[8] == "hidden", "{tool}");
            assert_eq!(declaration.advanced, fields[9] == "advanced", "{tool}");
            assert_eq!(declaration.common, fields[10] == "common", "{tool}");
            assert_eq!(
                declaration.controlled_by.is_none(),
                fields[18] == "none",
                "{tool}"
            );
        }
    }
    // `--input` is a `List<GATKPath>`, so it is a collection whose ELEMENT is the path: the
    // conversion a value goes through is the element's, and the collection is only how many.
    let input = find(COUNTREADS, "input").expect("input");
    assert!(input.collection);
    assert_eq!(input.type_name, "GATKPath");
    // A flag is a boolean that takes no value, and it is the one place primitiveness shows.
    let flag = find(COUNTREADS, "create-output-bam-index").expect("the index flag");
    assert!(flag.flag);
    assert_eq!(flag.type_name, "Boolean");
    // An enum argument reports the enum, which is what the parser converts the string with.
    assert_eq!(
        find(COUNTREADS, "interval-set-rule")
            .expect("the set rule")
            .type_name,
        "IntervalSetRule"
    );
}

/// The rest of what a parser and a usage text are built from.
#[test]
fn the_declarations_carry_what_a_parser_and_a_usage_need() {
    let text = golden();
    // The documentation is the annotation's own string, and it is what the usage wraps.
    let doc = |tool: &str, index: usize| {
        field(&text, "doc", tool)
            .into_iter()
            .find(|row| row.starts_with(&format!("{index}\t")))
            .map(|row| row.split_once('\t').expect("an index").1.to_string())
            .unwrap_or_else(|| panic!("doc/{tool}/{index}"))
    };
    for (tool, list) in TOOLS {
        for (index, declaration) in list.iter().enumerate() {
            assert_eq!(declaration.doc, doc(tool, index), "{tool}");
        }
    }
    // A bound is four doubles and not two: an undeclared range is the infinities, and the
    // recommended one is declared BESIDE the hard one.
    let plain = find(COUNTREADS, "output").expect("output");
    assert_eq!(plain.min_value, f64::NEG_INFINITY);
    assert_eq!(plain.max_recommended_value, f64::INFINITY);
    // The read-filter SELECTOR is not itself controlled by the descriptor: it is the common
    // argument that names filters, and it belongs to the tool's own namespace.
    let selector = find(COUNTREADS, "read-filter").expect("the read filter argument");
    assert!(selector.controlled_by.is_none());
    assert!(selector.common);
    assert!(selector.collection);
    // What the descriptor controls is the arguments the individual FILTERS declare, and they are
    // in the parser whether or not anybody selected the filter that owns them. Two of them read
    // as required, which is the descriptor's business and not a command line's: nothing on a
    // plain `CountReads` line supplies them and the parse succeeds anyway.
    let controlled = find(COUNTREADS, "read-filter-expression").expect("a filter's argument");
    assert_eq!(
        controlled.controlled_by,
        Some("GATKReadFilterPluginDescriptor")
    );
    assert!(controlled.required);
    assert!(!controlled.common);
    // And a tool that is no walker has no plugin-controlled argument at all.
    assert!(INDEXFEATUREFILE
        .iter()
        .all(|declaration| declaration.controlled_by.is_none()));
}

/// The parse decisions follow from the declarations: a collection may be repeated, a scalar not.
#[test]
fn the_parse_decisions_follow_from_the_declarations() {
    let text = golden();
    assert!(find(COUNTREADS, "input").expect("input").collection);
    assert!(!find(PRINTREADS, "output").expect("output").collection);
    let outcome = |tool: &str, case: &str| {
        field(&text, "parse", tool)
            .into_iter()
            .find(|row| row.starts_with(&format!("{case}\t")))
            .unwrap_or_else(|| panic!("{tool}/{case}"))
    };
    assert!(outcome("CountReads", "input-twice").ends_with("ok"));
    assert!(outcome("PrintReads", "output-twice").contains("BadArgumentValue"));
    // A required argument left out is refused by name, before the tool runs.
    assert!(outcome("CountReads", "no-arguments").contains("Argument input was missing"));
    assert!(outcome("PrintReads", "input-only").contains("Argument output was missing"));
    assert!(find(PRINTREADS, "output").expect("output").required);
    // And the intervals argument is a collection, which is why two of them parse.
    assert!(find(COUNTREADS, "intervals").expect("intervals").collection);
    assert!(outcome("CountReads", "two-intervals").ends_with("ok"));
}
