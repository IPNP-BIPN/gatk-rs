//! `htsjdk.variant.variantcontext.filter.JavascriptVariantFilter` (htsjdk 5.0.0), which runs a
//! script in Nashorn with `variant` and `header` bound, here run in Boa, a JavaScript engine
//! written in Rust, with the two bound to stand-ins written in JavaScript
//! (`javascript_variant_filter.js`) over what this module works out with htsjdk-rs.
//!
//! # What `AbstractJavascriptFilter` does, and so what this does
//!
//! The script is compiled once, when the filter is built, and a compile error is a
//! `RuntimeScriptException` "Script error in input". It is then evaluated once per record with the
//! same `SimpleBindings`, which is where Nashorn keeps the global it creates on the first
//! evaluation: a top-level `var` survives from one record to the next, and does so here because
//! one Boa context serves the whole run. The last expression's value accepts the record when it is
//! `true`, or a number whose `intValue()` is 1, so `1.5` accepts and `"true"` does not.
//!
//! # What the stand-ins do not do
//!
//! Measured against Nashorn on the coverage corpus (`docs/filtervcf-javascript.md`): a method of
//! the Java objects the stand-ins do not carry, a Java class reached through `Java.type`, a string
//! read as `s.length()` (a JavaScript string's `length` is a number, which Nashorn answers both
//! ways), or any other run-time script error stops the run here as a port limitation, where
//! Nashorn's message would name the error in its own words. A refusal is the outcome where this
//! engine and Nashorn could disagree, not a different answer.
use boa_engine::object::{FunctionObjectBuilder, JsObject};
use boa_engine::{js_string, Context, JsResult, JsValue, NativeFunction, Script, Source};
use htsjdk_vcf::attributes::{as_boolean, as_double, as_int, as_string};
use htsjdk_vcf::genotype_type::{determine_type, sample_names_ordered_by_name};
use htsjdk_vcf::header::{HeaderLine, VcfHeader};
use htsjdk_vcf::variant::{Value, VariantContext};
use serde_json::{json, Value as Json};

use crate::local_assembler_engine::java_hash_map_order;
use crate::remove_nearby_indels::{variant_type, VariantType};

const PRELUDE: &str = include_str!("javascript_variant_filter.js");

/// A script that would not compile, which `AbstractJavascriptFilter` throws as
/// `new RuntimeScriptException("Script error in input", err)`.
pub const SCRIPT_ERROR: &str = "Script error in input";

pub struct JavascriptVariantFilter {
    context: Context,
    script: Script,
    wrap_variant: JsObject,
    wrap_header: JsObject,
}

impl JavascriptVariantFilter {
    /// The constructor: the script compiled, and nothing evaluated yet. `Err` is a compile error.
    pub fn new(source: &str) -> Result<Self, String> {
        let mut context = Context::default();
        let refuse =
            FunctionObjectBuilder::new(context.realm(), NativeFunction::from_fn_ptr(refuse))
                .build();
        let wrappers = context
            .eval(Source::from_bytes(PRELUDE))
            .and_then(|prelude| prelude.to_object(&mut context))
            .and_then(|prelude| prelude.call(&JsValue::undefined(), &[refuse.into()], &mut context))
            .and_then(|value| value.to_object(&mut context))
            .map_err(|error| error.to_string())?;
        let wrap_variant = function(&wrappers, "variant", &mut context)?;
        let wrap_header = function(&wrappers, "header", &mut context)?;
        // Nashorn's `print`, which writes to the engine's standard output.
        context
            .register_global_builtin_callable(
                js_string!("print"),
                0,
                NativeFunction::from_fn_ptr(print),
            )
            .map_err(|error| error.to_string())?;
        let script = Script::parse(Source::from_bytes(source), None, &mut context)
            .and_then(|script| script.codeblock(&mut context).map(|_| script))
            .map_err(|error| error.to_string())?;
        Ok(Self {
            context,
            script,
            wrap_variant,
            wrap_header,
        })
    }

    /// `bindings.put("header", header)`: the header as the run will have it by the first record,
    /// with the tool's own FILTER and FORMAT lines added.
    pub fn bind_header(&mut self, header: &VcfHeader) -> Result<(), String> {
        let ids = |wanted: fn(&HeaderLine) -> Option<&String>| {
            header.lines.iter().filter_map(wanted).collect::<Vec<_>>()
        };
        let data = json!({
            "samples": header.samples,
            "info": ids(|line| compound(line, "INFO")),
            "format": ids(|line| compound(line, "FORMAT")),
            "filter": ids(|line| match line {
                HeaderLine::Filter { id, .. } => Some(id),
                _ => None,
            }),
        });
        let wrap = self.wrap_header.clone();
        self.bind("header", &wrap, &data)
            .map_err(|error| error.to_string())
    }

    /// `accept(record)`. `Err` is a run-time script error, with the engine's message, or what a
    /// stand-in refused, even where the script caught it.
    pub fn accept(&mut self, record: &VariantContext) -> Result<bool, String> {
        REFUSED.take();
        let run = |filter: &mut Self| -> JsResult<JsValue> {
            let wrap = filter.wrap_variant.clone();
            filter.bind("variant", &wrap, &record_json(record))?;
            filter.script.evaluate(&mut filter.context)
        };
        let result = run(self);
        if let Some(refusal) = REFUSED.take() {
            return Err(refusal);
        }
        // The first line: Boa's message goes on with the script's stack.
        let result = result.map_err(|error| {
            error
                .to_string()
                .lines()
                .next()
                .unwrap_or_default()
                .to_string()
        })?;
        Ok(match result.as_boolean() {
            Some(accepted) => accepted,
            // `((Number) result).intValue() == 1`, a truncation, so [1, 2) and nothing else.
            None => result.as_number().is_some_and(|n| (1.0..2.0).contains(&n)),
        })
    }

    fn bind(&mut self, name: &str, wrap: &JsObject, data: &Json) -> JsResult<()> {
        let data = JsValue::from_json(data, &mut self.context)?;
        let object = wrap.call(&JsValue::undefined(), &[data], &mut self.context)?;
        self.context.global_object().set(
            boa_engine::JsString::from(name),
            object,
            false,
            &mut self.context,
        )?;
        Ok(())
    }
}

fn function(object: &JsObject, name: &str, context: &mut Context) -> Result<JsObject, String> {
    object
        .get(boa_engine::JsString::from(name), context)
        .and_then(|value| value.to_object(context))
        .map_err(|error| error.to_string())
}

thread_local! {
    /// What a stand-in refused during the current record.
    static REFUSED: std::cell::Cell<Option<String>> = const { std::cell::Cell::new(None) };
}

fn refuse(_: &JsValue, args: &[JsValue], context: &mut Context) -> JsResult<JsValue> {
    let what = args
        .first()
        .cloned()
        .unwrap_or_default()
        .to_string(context)?;
    REFUSED.set(Some(what.to_std_string_escaped()));
    Ok(JsValue::undefined())
}

fn print(_: &JsValue, args: &[JsValue], context: &mut Context) -> JsResult<JsValue> {
    let parts = args
        .iter()
        .map(|arg| arg.to_string(context).map(|s| s.to_std_string_escaped()))
        .collect::<JsResult<Vec<_>>>()?;
    crate::tool_stdout::println(&parts.join(" "));
    Ok(JsValue::undefined())
}

fn compound<'a>(line: &'a HeaderLine, wanted: &str) -> Option<&'a String> {
    match line {
        HeaderLine::Compound { key, id, .. } if key == wanted => Some(id),
        _ => None,
    }
}

/// `VariantContext.Type`'s constant names.
fn type_name(kind: VariantType) -> &'static str {
    match kind {
        VariantType::NoVariation => "NO_VARIATION",
        VariantType::Snp => "SNP",
        VariantType::Mnp => "MNP",
        VariantType::Indel => "INDEL",
        VariantType::Symbolic => "SYMBOLIC",
        VariantType::Mixed => "MIXED",
    }
}

/// A double that may be infinite or NaN, which JSON cannot carry, as `Number()` reads it.
fn number_text(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else if value.is_infinite() {
        if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string()
    } else {
        format!("{value:?}")
    }
}

/// What `getAttribute(key)` hands a script: a String, a List of them, or a Boolean.
fn value_json(value: &Value) -> Json {
    match value {
        Value::Missing => json!("."),
        Value::Int(n) => json!(n),
        Value::Double(d) => json!(number_text(*d)),
        Value::Bool(b) => json!(b),
        Value::Str(s) => json!(s),
        Value::List(items) => Json::Array(items.iter().map(value_json).collect()),
    }
}

/// One attribute, and what each typed accessor answers for it: a value, `"default"` where the
/// accessor returns the caller's default, or null where it throws.
fn attribute_json(key: &str, value: &Value) -> Json {
    let int = match value {
        Value::Missing => json!("default"),
        _ => as_int(Some(value), 0).map_or(Json::Null, |n| json!(n)),
    };
    json!({
        "key": key,
        "value": value_json(value),
        "int": int,
        "double": as_double(Some(value), 0.0).map_or(Json::Null, |d| json!(number_text(d))),
        "string": as_string(Some(value), ""),
        "boolean": as_boolean(Some(value), false).map_or(Json::Null, |b| json!(b)),
    })
}

fn record_json(record: &VariantContext) -> Json {
    let alleles: Vec<Json> = record
        .alleles
        .iter()
        .map(|allele| {
            json!({
                "display": allele.display_string(),
                "ref": allele.is_reference(),
                "symbolic": allele.is_symbolic(),
                "noCall": allele.is_no_call(),
            })
        })
        .collect();
    let genotypes: Vec<Json> = record
        .genotypes
        .iter()
        .map(|genotype| {
            // A decoded genotype holds the record's own alleles, so `==` between them is true.
            let indices: Vec<i64> = genotype
                .alleles
                .iter()
                .map(|allele| {
                    record
                        .alleles
                        .iter()
                        .position(|a| a == allele)
                        .map_or(-1, |i| i as i64)
                })
                .collect();
            json!({
                "name": genotype.sample_name,
                "alleles": indices,
                "phased": genotype.phased,
                "type": determine_type(genotype).name(),
                "gq": genotype.gq.unwrap_or(-1),
                "dp": genotype.dp.unwrap_or(-1),
                "ad": genotype.ad,
                "pl": genotype.pl,
                "ft": genotype.filters,
                "extended": genotype
                    .extended
                    .iter()
                    .map(|(key, value)| json!({"key": key, "value": value_json(value)}))
                    .collect::<Vec<_>>(),
            })
        })
        .collect();
    // `builder.filters(new HashSet<>(filters))`, copied into a LinkedHashSet in that order.
    // shortcut: a HashSet built from twelve or more filters starts at a larger table than one
    // filled from empty; nothing in a record carries that many.
    let filters = record.filters.as_ref().map(|filters| {
        let mut distinct: Vec<String> = Vec::new();
        for filter in filters {
            if !distinct.contains(filter) {
                distinct.push(filter.clone());
            }
        }
        java_hash_map_order(&distinct)
            .into_iter()
            .map(|i| distinct[i].clone())
            .collect::<Vec<_>>()
    });
    json!({
        "contig": record.contig,
        "start": record.start,
        "end": record.stop,
        "id": record.id,
        "alleles": alleles,
        "type": type_name(variant_type(record)),
        "qual": number_text(record.phred_scaled_qual()),
        "hasQual": record.has_log10_p_error(),
        "filters": filters,
        "attributes": record
            .attributes
            .iter()
            .map(|(key, value)| attribute_json(key, value))
            .collect::<Vec<_>>(),
        "genotypes": genotypes,
        "orderedByName": sample_names_ordered_by_name(record),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::reader::read_vcf;

    const VCF: &str = "##fileformat=VCFv4.2\n\
        ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
        chr1\t100\t.\tA\tC\t50\tq10;a1\tDP=7\tGT:DP\t0/1:120\t1/1:150\n\
        chr1\t1001\trs1\tAT\tA\t.\tPASS\tDP=.\tGT:DP\t0/0:3\t./.:.\n";

    fn verdicts(script: &str) -> Result<Vec<bool>, String> {
        let file = read_vcf(VCF).map_err(|failure| failure.error.message())?;
        let mut filter = JavascriptVariantFilter::new(script)?;
        filter.bind_header(&file.header)?;
        file.records
            .iter()
            .map(|record| filter.accept(record))
            .collect()
    }

    #[test]
    fn the_last_value_decides_as_abstract_javascript_filter_reads_it() {
        assert_eq!(
            verdicts("variant.isSNP() && variant.getStart() > 10;"),
            Ok(vec![true, false])
        );
        assert_eq!(
            verdicts("variant.start % 2 == 0 ? 1.5 : 'true'"),
            Ok(vec![true, false])
        );
        assert_eq!(verdicts("2"), Ok(vec![false, false]));
        // A top-level var lives in the one global the run keeps.
        assert_eq!(
            verdicts("if (typeof n == 'undefined') var n = 0; n++; n == 2"),
            Ok(vec![false, true])
        );
    }

    #[test]
    fn the_stand_ins_answer_as_the_java_objects_do() {
        let script = "var g = variant.getGenotypes();\n\
            variant.getReference() + '' == (variant.isSNP() ? 'A*' : 'AT*')\n\
            && variant.getAlleles().toString() == (variant.isSNP() ? '[A*, C]' : '[AT*, A]')\n\
            && g.get(0).getAllele(0) == variant.getReference()\n\
            && header.genotypeSamples.get(1) == 's2' && header.hasFilterLine('q10') == false\n\
            && variant.getFilters().toString() == (variant.isSNP() ? '[a1, q10]' : '[]')\n\
            && variant.getAttributeAsInt('DP', -1) == (variant.isSNP() ? 7 : -1) == variant.isSNP()\n\
            && variant.hasID() == !variant.isSNP() && g.get('s2').getDP() == (variant.isSNP() ? 150 : -1)\n\
            && (variant.getHetCount() == 1) == variant.isSNP() && !(variant instanceof Object)";
        // `DP=.` is the text ".", which getAttributeAsInt parses and throws on.
        assert!(verdicts(script)
            .unwrap_err()
            .contains("getAttributeAsInt throws on DP"));
        let script = script.replace(
            "variant.getAttributeAsInt('DP', -1) == (variant.isSNP() ? 7 : -1) == variant.isSNP()",
            "true",
        );
        assert_eq!(verdicts(&script), Ok(vec![true, true]));
    }

    #[test]
    fn filter_vcf_sets_the_script_name_on_what_it_rejects_and_declares_it() {
        use crate::filter_vcf::{filter_to, ScriptFilter, Thresholds};
        let mut script = ScriptFilter {
            name: "keep_snps".to_string(),
            path: "/x/keep_snps.js".to_string(),
            script: JavascriptVariantFilter::new("variant.isSNP()").unwrap(),
        };
        let out = filter_to(VCF, &Thresholds::default(), false, Some(&mut script)).unwrap();
        assert!(out.contains(
            "##FILTER=<ID=keep_snps,Description=\"Variant Filtered by Javascript file /x/keep_snps.js\">"
        ));
        // No GQ is -1, below MIN_GQ=0, so the first site's variant genotypes are all filtered.
        let filters: Vec<&str> = out
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').nth(6).unwrap())
            .collect();
        assert_eq!(filters, ["AllGtsFiltered", "keep_snps"]);
    }

    #[test]
    fn a_script_that_does_not_compile_is_refused_when_the_filter_is_built() {
        assert!(JavascriptVariantFilter::new("record.getReadName( ;").is_err());
        assert!(verdicts("variant.getCalledChrCount() > 0").is_err());
        // A member the stand-ins lack is refused even where the script catches the error.
        assert!(
            verdicts("try { variant.getCalledChrCount(); } catch (e) {} true")
                .unwrap_err()
                .contains("getCalledChrCount is not mirrored")
        );
    }
}
