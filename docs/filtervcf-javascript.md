# FilterVcf's JAVASCRIPT_FILE: Nashorn's semantics in a Rust engine

Picard's `FilterVcf --JAVASCRIPT_FILE` hands the script to htsjdk's `JavascriptVariantFilter`,
which compiles it in Nashorn (15.4 in the GATK 4.7.0.0 jar) and evaluates it once per record with
`variant` (a `VariantContext`) and `header` (a `VCFHeader`) bound. The port runs the script in
[Boa](https://github.com/boa-dev/boa) (`boa_engine` 0.22, no default features, MIT), with the two
bound to stand-ins written in JavaScript over what htsjdk-rs works out for the record
(`crates/gatk-tools/src/javascript_variant_filter.{rs,js}`).

## What the reference does, and the port with it

- **Compile once, when the filter is built**: after the header is read and before the dictionary
  check, so a script that does not compile is `RuntimeScriptException: Script error in input`
  even on an input the dictionary check would refuse. A missing file is
  `PicardException: javascript-related error`.
- **One global for the run**: Nashorn keeps the global it creates in the `SimpleBindings` the
  filter reuses, so a top-level `var`, or an undeclared assignment, carries from record to record.
- **The verdict**: `true`, or a number whose `intValue()` is 1 (so `1.5` accepts), accepts; any
  other value, the string `"true"` included, rejects, and the record is filtered under the
  script's base name, declared as `Variant Filtered by Javascript file <path>`.
- **`print`** writes to standard output, joined by spaces.

## The stand-ins

`VariantContext`, `Allele`, `Genotype`, `GenotypesContext`, `VCFHeader`, `List`, `Set`, `int[]` and
the two `Type` enums carry the methods the scripts in Picard's and htsjdk's tests reach and the
common accessors around them (positions, alleles, type, QUAL, FILTER, INFO with the typed
accessors, genotypes with GQ/DP/AD/PL/FT and extended attributes, genotype counts, the header's
samples and its INFO/FORMAT/FILTER ids), plus Nashorn's bean properties (`variant.start`,
`header.genotypeSamples`). They follow Nashorn where it differs from plain JavaScript, each point
measured: a Java object is not `instanceof Object` and has no Array methods; a `List` indexes,
has `length` and enumerates its indices, a `Set` has `length` but neither; `==` between an object
and a string compares `toString()`; alleles of a genotype are `==` the record's; a filter set
iterates in `HashSet` order; `getAttributeAsInt` throws on `DP=.` and returns the default on a bare
`DP`. Java strings' `equals`, `equalsIgnoreCase`, `isEmpty` and `contains` are on
`String.prototype`, since Nashorn's strings are `java.lang.String`.

## What is excluded, and why it is a refusal

A script that reaches past the stand-ins stops the run as `gatk_rs::PortLimitation`, never with a
different answer. The stand-ins report what they refuse to the Rust side as well as throwing it,
so a script that catches the error is refused all the same:

- **a member the stand-ins do not carry**, looked up on any of them (the end of their prototype
  chain answers every other name with a refusal). Nashorn would call the Java method, or answer
  `undefined` for a name that is not one; `list[i]` past the end is refused for the same reason.
- **an exception a Java method throws** (`List.get` out of range, `getAttributeAsInt` on `DP=.`),
  whose class and wording are Java's.
- **`toString` of a VariantContext, Genotype, GenotypesContext, VCFHeader or `int[]`**: htsjdk's
  formatting, or an identity hash, is not reproduced.
- **an uncaught run-time error of the script's own**. Nashorn wraps it as `RuntimeException:
  ScriptException: <Nashorn's message>`; the wording is Nashorn's and Boa's differs.
- **`s.length()` on a string**. Nashorn answers both `s.length` and `s.length()`; a JavaScript
  string's `length` is a number and cannot also be a method, so the call is a `TypeError` here.
- **Java interop** (`Java.type`, `java.lang.*`, `JavaImporter`), which no engine without a JVM has.

Four differences are silent, and are excluded here rather than refused: a script's own error that
the script catches and inspects (the message is each engine's own), reflection over a host
object (`Object.keys(variant)` is a `TypeError` in Nashorn and `[]` here; `JSON.stringify(variant)`
is `undefined` there and `"{}"` here), `typeof Java` (`object` there), and syntax: Nashorn is
ECMAScript 5.1 with its own extensions (`for each`, expression closures), Boa is ES2023, so a
script using `let`, arrow functions or a Nashorn extension compiles in one and not the other. No
script in Picard's or htsjdk's tests, nor in the coverage array, does any of these.

## Measured

`tools/coverage/compare_javascript_filter.py` runs every script on every VCF in both, in the pinned
image, and compares exit status, first error line, standard output and the output VCF.

On eight VCFs (the array's five, Picard's `testFiltering.vcf` and `testFilteringSitesOnly.vcf`, and
htsjdk's 1,182-sample `ILLUMINA.wex...genotypes.1000.vcf`), on 2026-11-19:

| script | agree | refused | differ |
|---|---|---|---|
| the array's `filter_variants.js` and `filter_broken.js` | 16 | 0 | 0 |
| Picard's `TestFilterVcf.testJavaScript` (`variant.getStart()%5 != 0`) | 8 | 0 | 0 |
| htsjdk's `variantFilter01.js` and `variantFilter02.js` | 16 | 0 | 0 |
| a probe printing about 110 expressions per record over every stand-in | 6 | 2 | 0 |
| two probes written to be refused (`s.length()`, a caught missing method) | 0 | 16 | 0 |

The probe's two refusals are the sites-only VCFs, where it catches the `IndexOutOfBoundsException`
of `getGenotypes().get(0)`. Picard's other scripts (`filterOddStarts.js`,
`filterReadsWithout5primeSoftClip.js`) and htsjdk's `samFilter0*.js` are `FilterSamReads`' and
bind a `SAMRecord`, so they are not this tool's. The coverage array's 20 rows match, 11 of them
`--JAVASCRIPT_FILE`.
