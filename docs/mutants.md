# Mutants: where the goldens do not look

The goldens and the covering arrays say the port matches the reference where they look.
[cargo-mutants](https://github.com/sourcefrog/cargo-mutants) says where they do not look: it puts
one bug into the port, runs the tests, and reports each bug no test caught (#1392). It never ships
in the binary, so it puts no byte at risk.

```bash
cargo install --locked cargo-mutants
cargo mutants -p gatk-readfilter -j 4
cargo mutants -p gatk-barclay -j 4
```

`.cargo/mutants.toml` holds the exclusions, each with its reason. A surviving mutant has two
possible fates: a real gap gets a record or a case from the oracle, and a mutant no input can
tell apart from the port gets an exclusion. Nothing else goes in that file.

## Survivors per crate

| crate | mutants | caught | missed | unviable | timeouts | wall time |
|---|---:|---:|---:|---:|---:|---:|
| gatk-readfilter, before | 518 | 438 | **63** | 15 | 2 | 7 min |
| gatk-readfilter, now | 484 | 469 | **0** | 13 | 2 | 5 min |
| gatk-barclay, before | 382 | 226 | **132** | 19 | 5 | 3 min |
| gatk-barclay, now | 367 | 342 | **0** | 18 | 7 | 3 min |

Measured on cargo-mutants 27.1.0, Homebrew cargo 1.99.1, `-j 4` on a shared 16-core Apple-silicon
machine. A timeout is a mutant that turns a loop infinite (`offset += length` into `*=` from 0):
the test hangs instead of passing, so it counts as caught.

The other crates have not been run. Their sizes say how: `gatk-corpus` 53 mutants,
`gatk-annotation` 2,259, `gatk-fermilite` 3,236, `gatk-cli` 5,740, `gatk-engine` 15,001,
`gatk-tools` 44,064.

## How the 63 of gatk-readfilter were sorted

**35 were real gaps**, and the oracle now closes each one. `ReadFilterDump` gained 18 records,
three filter instances and, in `CountingFilterDump`, one composition:

- Cigars (18 mutants of `is_valid`, `is_real`, `is_padding`, `starts_or_ends_with_deletion`):
  the corpus had no P at all, no plain internal deletion (`5M1D5M`), no `2S2S6M`, no same indel
  twice with an M or a P between, and no zero-length element. A port that lost the whole
  padding rule passed.
- `MateOnSameContigOrNoMappedMate`, `NonZeroReferenceLengthAlignment`, `SeqIsStored` (4): no
  unmapped mate placed on another contig, no `0M10I`, no read without a sequence.
- `ReadTagValueFilter` (5): no instance used `GREATER` or `LESS_OR_EQUAL`, so deleting either
  parse arm changed nothing. `AmbiguousBaseReadFilter` (1): no limit equal to a read's count.
- Flow filters (5): no asymmetry in the last hmer, no trailing hard clip, no `tp` on either
  bound of `[0, maxHmer]`, no `tp` stored as text.
- `PlatformUnitReadFilter` (1): no read carrying its own PU.
- The counting summary (1): no nested tree whose left child rejects nothing.

The corpus is shared, so the records also reach `ReadCoordinateDump` (1,074 lines to 1,401),
`ReadClipperDump` (3,211 to 4,189), `AssemblyRegionDump` (4,058 to 6,650) and `CovariatesDump`
(2,903 to 3,632). The port agreed with the oracle on every new line the first time, in all six
goldens: these were gaps in the evidence, not bugs. The one test that needed a change was the
covariates one, which read the dump's `empty` marker for a read with no bases as an exception.

**16 were dead code**, removed rather than tested. `with_header::header_declares_read_group` had
no caller. The JEXL filter's read path (`attribute_as_string`, `format_java_float`,
`context_for`, `compile`, `test`) had none either: no runner builds that filter, and the
conformance test calls `test_context` with the context the reference built. That path also held
a known divergence on array tags, recorded in `transfer_read_tags.rs`, now gone with it.

**10 are equivalent** and excluded: `ValidAlignmentStartReadFilter` and
`ValidAlignmentEndReadFilter` cannot answer false on a read the adapter presents (5), "is
padding" answering true takes no branch the others do not (1), and four `is_valid` mutants
differ only on cigars with two indels in a row, which `is_good` rejects before it matters (4).
Four of those are pinned by line, so moving the code brings them back as survivors.

**2 are caught downstream**: `Counting::filtered_count` returning 0 or 1 makes
`gatk-tools/tests/multi_pass.rs` fail. A run scoped to the crate does not build that test, so
they are excluded with that reason.

## How the 132 of gatk-barclay were sorted

**51 were caught by a test in the wrong crate.** `value_classes.rs` exercises nothing but
`gatk_barclay` (`Float.valueOf` and `Long.valueOf`'s grammars against the
`tool-argument-value-classes` golden) and sat in `gatk-cli`, which a run scoped to the crate does
not build. It moved to `gatk-barclay/tests`, unchanged.

**54 were real gaps**, closed with oracle rows in four dumps. The port agreed with the oracle on
every new row the first time, except where it says otherwise below.

- The order values are set in (22, `propagation_order` and the `HashMap` it transcribes): no
  case had two fields failing at once, so a port that set them in any order passed. Now every
  pair of the value model's failing fields is a case, and two more targets sit on a resize
  threshold: eleven names, which with `[arguments]` make twelve and size `toJavaUtilMap()`'s copy
  at thirty-two buckets rather than sixteen, and twelve specs with `CallableLoci`'s
  `{mlmq, max-low-mapq}` pair, which `arrangeOptions` sorts long name first.
- Positional arguments (10): no tool in the crate declared any, and no case gave two positional
  values to one that did not (`,a{b}`). A lone `-` as a positional value or as a tagged value (2).
- A bounded `Long` (1), which `--splitting-index-granularity` is on every Spark tool.
- `Long.valueOf` below zero, at and past `Long.MIN_VALUE`, and forty digits (5); `Float.valueOf`
  on `+Infinity` and on hexadecimal with several digits or a fraction (10).
- Plugin selection (4): the `allowed` rows of `plugin-argument-ownership` are replayed through
  `plugin_is_selected`, and three annotation rows were added, a tool-default group and a group
  named on the command line.

**1 was a divergence.** `is_number` left out `Byte`, so barclay never range-checked
`DepthOfCoverage`'s `--min-base-quality`, and the runner did it itself after parsing, which
reports a missing required argument before an out-of-range quality where the reference reports
the quality. A bounded `byte` is now in the value model, and barclay refuses `-1` with
`allowed range [0.0, 127.0].` as the reference does (bounds as doubles: a `Byte` is not an
`Integer`). Once `Byte` was in, the "only numeric types" guard decided nothing, since only numeric
types can have bounds, and it went too.

**14 are port-side and got unit tests**: `Math.rint`'s half-to-even (7), `hasRecommendedRange`
on one bound (4), `Filesystem::read` (2), `override_value` (1).

**11 were dead code**, removed: `ValueClass::simple_name` and `long_names` had no caller, every
caller dropped the list `set_mode_values` returned, GATK never sets `APPEND_TO_COLLECTIONS`, a
float suffix guard tested for `infinity` on a string that ends in a suffix letter, and the flag's
"does not look like an option" test can never decide, since no boolean spelling starts with a
dash.

**1 is equivalent** and excluded: `java_long` stopping its digit loop at `u64::MAX` rather than
past it refuses the same strings. Pinned by line.

## What it costs, and where it runs

Scoped to its own crate, a mutant of `gatk-readfilter` costs about a second of build and a
second of test after a 15-second baseline. Adding the downstream tests (`gatk-tools`, `gatk-cli`,
debug profile) costs about ten minutes a mutant, so a whole crate cannot be run that way; only
named survivors are checked downstream, by hand.

Decision for now: not in CI on every push. It runs by hand for a crate whose code or corpus
changes, and for a validation report (#87). `gatk-barclay` costs three minutes scoped to itself,
and only `gatk-cli` depends on it, so a named survivor can be checked downstream in seconds
(`cargo test -p gatk-cli --test plugin_ownership ...`). The next candidate is `gatk-corpus`, then `gatk-engine` split with `--shard k/n` across a scheduled
job, once a first shard says what its tests cost.
