# Roadmap to 100% byte-identical reproduction

Tracking document for the goal: a byte-for-byte reproduction of the **entire** GATK 4.6.2.0,
Picard 3.4.0 and htsjdk 4.2.0 tool set in Rust. Deliberately honest about scale (the plan sizes
the whole program at 40 to 100 person-years), and meant to be followed and ticked off.

Status legend: `[ ]` not started, `[~]` in progress or partial, `[x]` done **and** oracle-backed.

A box is only `[x]` when a golden produced by the pinned container on a real x86-64 runner is
committed, a Rust test compares against it, and CI re-derives it. Code that works but has no
golden stays `[~]`.

---

## Where we are (measured 2026-07-30)

| Repo | Scope | State |
|---|---|---|
| **htsjdk-rs** | the I/O and math foundation | substantially built; CRAM, GKL-exact deflate, full VCF and the jmath conformance corpus remain |
| **picard-rs** | 109 tools | ~50 tools have a first slice, ~43 with an oracle-backed conformance suite; many are partial (default paths only). The harness is generated from a manifest, the fuzzer and the determinism gate run in CI, and argument coverage is measured for 2 tools |
| **gatk-rs** | 202 tools | 6 crates, **58 conformance suites, all oracle-backed**; 3 tools byte-identical, and the annotation archetype opened with 53 of 54 annotations measured. **No performance number exists yet for any of it** — see Milestone S |

Totals from the generated inventory (`tools/inventory`): **311 tools** (202 GATK-origin,
109 Picard-origin), **39 Spark**, ~13,130 arguments. Non-Spark: 163 GATK + 109 Picard.

### What "100% repro" means, per tool

A tool is **bit-identical** (the target) when, under the declared canonicalization, it is
byte-equal to the oracle for:

1. every applicable input under the pinned clones' `src/test/resources/**`;
2. a t-wise covering array over its arguments (t=2 everywhere, t=3 on the critical path);
3. the coverage-guided differential fuzzer reaching its branch-coverage threshold with no
   divergence;

with **zero quarantined fields**. Any quarantined field downgrades it to **bio-identical** with
the quarantine list attached. Everything ported so far reproduces the paths its suite covers,
which is not yet this bar for any tool.

---

## Milestone G1: the GATK engine

The single biggest unlock: 163 non-Spark GATK tools stand on it.

**Closed.** Every box below is `[x]` except four, and none of the four is G1 work:

| open item | where it lives now | why |
|---|---|---|
| the Tribble index (G1.3) | Milestone H | it is an htsjdk capability GATK consumes |
| the multi-input walkers (G1.6) | Milestone H | they merge with `VCFUtils.smartMergeHeaders` and `VariantContextComparator` |
| `AllelePseudoDepth` (G1.7) | **done, G1.9** | the licence was never what blocked it. Both values leave through a `DecimalFormat` at two and four decimals, and the fixed point turned out bit-identical anyway. Oracle-backed, 46 calls, no divergence |
| `AssemblyComplexity` (G1.7) | Milestone G3 | it needs `Haplotype.getEventMap()`, the assembly event model |

57 conformance suites carry it, all oracle-backed.

One of those four came back and is now closed. **G1.9** below is `AllelePseudoDepth`, reopened not
because more effort was found for it but because the reason it was refused did not survive reading
the annotation — and then measured all the way through, ending oracle-backed with no divergence.
Three of the four remain open, none of them G1 work.

### G1.1 Read filters

- [x] 55 of the 56 `ReadFilterLibrary` filters, oracle-backed (79 instances, 59 records, 4,661
      decisions)
- [x] `CountingReadFilter` summary text, compared as bytes (the only place a conjunction's order
      is observable)
- [x] `JexlExpressionReadTagValueFilter`, the 56th, with commons-jexl 2.1.1's arithmetic ported
      under it (186 evaluations and 186 filter decisions). **All 56 read filters are now
      oracle-backed.**

### G1.2 Coordinates, cigars and clipping

- [x] `SimpleInterval` and the interval query parser, including the ambiguity rules
- [x] `SAMRecordToGATKReadAdapter` semantics (`isUnmapped` is three criteria, `isFirstOfPair`
      requires pairing)
- [x] `ReadUtils` reference-to-read coordinate mapping (872 probed positions)
- [x] `CigarBuilder` and `CigarUtils.clipCigar` (604 clips)
- [x] `ReadClipper`, all 14 entry points (3,068 clipped reads)
- [x] `ReadUtils` indel qualities (`BI`/`BD` tags, and the flat Q45 fallback)

### G1.3 Data sources and contexts

- [x] `ReferenceDataSource` through `CachingIndexedFastaSequenceFile` (45 queries, upper-casing
      and IUPAC flattening included)
- [x] `ReadsDataSource` interval queries against a fixture BAM and its `.bai` (29 queries,
      htsjdk's stateful single-pass filter included)
- [x] `ReferenceContext` and `ReadsContext` (352 window answers)
- [~] `FeatureDataSource` and `FeatureContext`: the lookahead cache, the trim that preserves file
      order, and the window arithmetic are ported and oracle-backed (20 queries at two lookahead
      settings), and the **Tribble index it needed is now in** (htsjdk-rs #83). The suite pins what a tool sees and does **not** distinguish the cache from a
      fresh query per call, which the manifest states. All three codecs `-L` reaches are ported
      and oracle-backed in htsjdk-rs (**BED**, **IntervalList**, and **VCF**'s `canDecode`, the one
      that decides by content). Still missing: the **Tribble index**, which is what turns a
      Feature file into a random-access source rather than a linear read, and belongs in
      Milestone H

### G1.4 Interval arguments

- [x] `-L` / `-XL`, padding, `--interval-set-rule`, `--interval-merging-rule`, subtraction
      (measured through `IntervalWalker`, 24 combinations)
- [x] interval **files**: `.list` and `.intervals` (lower-cased extension test, blank-line
      skipping, the empty-file refusal, and the order of the four tests in
      `parseIntervalArguments`), 13 arguments compared
- [x] the Feature-file path, all three codecs. `.bed`, `.interval_list` and VCF landed in
      htsjdk-rs with their own oracle-backed suites (90, 133 and 24 rows) and are registered in the
      seam, so all **25** measured `-L` arguments resolve identically and none is pending. The
      rows that cost the most: the two dictionaries an interval list is validated against, and the
      fact that VCF is the only codec deciding by **content**, so a `.list` holding a VCF body is a
      Feature file where one holding a BED body is not
- [x] `-L unmapped` end to end, measured through `ReadWalker` (5 runs): the tail comes after
      every interval, `-L unmapped` alone is a bounded traversal of nothing else, and an unmapped
      read carrying its mate's position is not in the tail at all

### G1.5 The pileup floor

- [x] `AlignmentStateMachine` (244 stops over 26 cigars)
- [x] `PileupElement` (217 elements, plus 231 `createPileupForReadAndOffset` calls including the
      offsets it refuses)
- [x] `ReadPileup`: the per-locus collection, its sorting, its sample split and the samtools
      overlap fix (3 pileups and 24 quality pairs)
- [x] `ReadStateManager`, `PerSampleReadStateManager` and `SamplePartitioner`, downsampling
      excepted (56 traversal steps over 5 runs)
- [x] `LIBSDownsamplingInfo` and the downsampling itself (`--max-depth-per-sample`), over both
      of GATK's static generators. `java.util.Random` from its published contract (98 sequences)
      carries `ReservoirDownsampler` (7 reservoirs, with the shared stream's position compared
      after each); commons-math3's `Well19937c` (255 sequences, 150 `nextBytes` lengths, 15
      interleaved streams) carries `RandomDataGenerator.nextPermutation` and `LevelingDownsampler`
      on top of it (17 permutations and 24 leveling cases, each with the stream position after).
      `nextGaussian` is refused rather than approximated: it runs through commons-math's own
      `FastMath`, so it belongs with jmath
- [x] `LocusIteratorByState`: one pileup per covered locus, both exclusions and the per-base
      adaptor test (148 pileups over 12 runs)
- [x] `IntervalAlignmentContextIterator`, `IntervalLocusIterator`, `IntervalOverlappingIterator`
      and `AlignmentContextIteratorBuilder`: 7 routing decisions and 99 contexts across every
      route, plus the sample-set iteration order the pileup element order depends on (7 orders,
      29 hashes), which had to be reproduced as a measured observable rather than ported

### G1.6 Walkers

- [x] `ReadWalker` (49 `apply` calls over 9 traversals)
- [x] `IntervalWalker` (25 `apply` calls over 24 argument combinations)
- [x] `LocusWalker`: 217 `apply` calls over 8 traversals, including the same interval run with
      and without `emitEmptyLoci` (6 calls against 51)
- [x] `VariantWalker` (12 traversals, 30 `apply` calls). The VCF reader it stands on landed in
      htsjdk-rs first: the header frame, the typed header lines, the site columns and the genotype
      columns, four oracle-backed suites there
- [x] `AssemblyRegionWalker` (the base of G3), and the four layers under it, each oracle-backed on
      its own: the activity profile that decides where a region starts and stops (20 kernels
      compared as raw bits, 8 profiles); `AssemblyRegion` itself, its two spans, its reads and its
      trimming (3481 comparator pairs, and the javadoc's own worked example shown to disagree with
      the code it documents); the locus iterators that manufacture an empty pileup on a coverage
      gap, without which a region ends at the last covered base; the traversal that turns loci into
      regions; and the walker, run through the real command line, where `--force-active` is shown to
      rewrite ten flags without moving one boundary
- [x] the multi-pass walkers (`MultiplePassVariantWalker`, `TwoPassVariantWalker`,
      `MultiplePassReadWalker`), 144 rows over ten runs. Three things happen between passes and
      none is guessable: the variant walker builds one counting filter before the loop and reuses
      it, so its counts accumulate (1, 2 and 3 drops for one, two and three passes over one
      filtered record) where the read walker builds a new one per pass and reports 3, 3, 3;
      `afterNthPass` runs after the **last** pass too, which its own javadoc denies; and
      `TwoPassVariantWalker`'s `afterNthPass` guard is `n == 0` then `n > 1`, so the call after the
      second pass matches neither branch and does nothing at all. Zero passes is legal throughout
- [ ] the multi-input walkers (`MultiVariantWalker` and `MultiVariantDataSource`), which merge
      several VCFs into one position-ordered stream. Not effort-blocked but layer-blocked: the
      merged header comes from htsjdk's `VCFUtils.smartMergeHeaders` and the merged order from
      `VariantContextComparator` over the merged dictionary, so it lands once those are in
      htsjdk-rs. `MultiVariantWalkerGroupedOnStart` and `MultiVariantWalkerGroupedByOverlap` sit
      on top of it

### G1.7 Annotations

- [~] the 54-annotation library. **Fifty-three** are ported and oracle-backed: the counting family
      (`ChromosomeCounts`, `SampleList`, `RawGtCount`, `Coverage`, `MappingQualityZero`,
      `CountNs`, `OriginalAlignment`) and the median family (`BaseQuality`, `MappingQuality`,
      `ReadPosition`, `FragmentLength`, i.e. MBQ/MMQ/MPOS/MFRL) and the rank-sum family
      (`BaseQRankSum`, `MQRankSum`, `ReadPosRankSum`, `ClippingRankSum`) and the strand-bias
      family (`FS`, `SOR`, `SB`) and three that read the matrix or the reference window alone
      (`AS_UNIQ_ALT_READ_COUNT`, `BQHIST`, `REF_BASES`) the depth family (`AD`, `AF`, `DP`) and the
      eight flow-based ones (`VARIANT_TYPE`, `INDEL_CLASSIFY`, `INDEL_LENGTH`, `HMER_INDEL_LENGTH`,
      `HMER_INDEL_NUC`, the motifs, `GC_CONTENT`, `CYCLESKIP_STATUS`) and three site statistics
      (`QD`, the genotype summaries, `LikelihoodRankSum`) and the cohort-heterozygosity pair
      (`ExcessHet`, `InbreedingCoeff`) with `MQ` and the tandem repeats (`STR`/`RU`/`RPA`) and
      the three allele-specific rank sums (`AS_BaseQRankSum`, `AS_MQRankSum`,
      `AS_ReadPosRankSum`, with the `Histogram` and `CompressedDataList` they travel as) and the
      allele-specific strand-bias pair (`AS_FS`, `AS_SOR`, over the `AS_SB_TABLE` they share) and
      the last three allele-specific ones (`AS_QD`, `AS_MQ`, `AS_InbreedingCoeff`, with the
      heterozygosity calculator under the third) and the pedigree pair (`hiConfDeNovo`/
      `loConfDeNovo`, `transmittedSingleton`/`nonTransmittedSingleton`, over a ported
      `MendelianViolation`) and the fragment pair (`F1R2`/`F2R1`, `FAD`, over a ported `Fragment`
      and `groupEvidence`) and `HaplotypeFilteringAnnotation` (`ASSEMBLED_HAPS`/`FILTERED_HAPS`,
      over a `Haplotype` and a likelihood matrix whose allele axis is one), together with the
      `InfoFieldAnnotation` interface and the machinery underneath it (`AlleleList`/`SampleList`
      and their permutation, the `AlleleLikelihoods` matrix and its best-allele search,
      `VariantContextGetters`, **`MannWhitneyU`** and **`FisherExactTest`**, with commons-math3's
      `FastMath.exp`, `FastMath.log`, `Gamma`, `Erf`, `NormalDistribution` and the saddle-point
      expansion under them, all oracle-backed in htsjdk-rs)
- [ ] the remaining 2, **neither of which is blocked on effort**. **The claim that
      most of the 54 wait on jmath did not survive a grep**: 10 of the 57 files in
      `tools/walkers/annotator` mention `MathUtils` at all, and what the annotators reach through
      `java.lang.Math` is `log`, `log10`, `sqrt` and `round`, all four already exact. What the
      first fifty-two actually waited on was engine machinery, and all of it is now ported. The
      third of the three came in with the haplotype-typed matrix: `HaplotypeFilteringAnnotation`
      is oracle-backed above, and `Haplotype` is ported as far as an allele axis needs it, which
      stops short of `EventMap`. The two that are left cannot be closed by working harder:
  - `AllelePseudoDepth` is **refused on the licence boundary**, not deferred. It ends in
    `SomaticLikelihoodsEngine.alleleFractionsPosterior`, whose fixed point runs through
    `NaturalLogUtils.normalizeFromLogToLinearSpace` and `logSumExp`, and both of those call
    `java.lang.Math.exp`. htsjdk-rs decision 0014 withdrew `Math.exp`: its only faithful port was
    an operation-by-operation transcription of HotSpot's x86 intrinsic, whose source is GPL2 with
    no Classpath Exception and therefore cannot be published under that crate's MIT licence. The
    `Gamma.digamma` this annotation also needs **is** ported (htsjdk-rs #62); the exponential is
    what stops it, and no amount of work on this side changes that
  - `AssemblyComplexity` needs `Haplotype.getEventMap()`, which is the assembly event model, so it
    belongs to **G3** with `HaplotypeCaller` rather than to G1. Listing it here would be counting
    the assembler as an annotation

### G1.9 `AllelePseudoDepth`, and the premise that refused it

G1.7 lists this annotation as **refused on the licence boundary**: it ends in
`SomaticLikelihoodsEngine.alleleFractionsPosterior`, which reaches `Math.exp` through
`NaturalLogUtils`, and the only exact port of `Math.exp` is a transcription of GPL2-only HotSpot
source (htsjdk-rs decision 0014).

**That framing does not survive reading the annotation.** Both values that depend on `exp` leave
through a formatter:

```java
private static DecimalFormat DEPTH_FORMAT    = new DecimalFormat("#.##");
private static DecimalFormat FRACTION_FORMAT = new DecimalFormat("#.####");
```

Two and four decimal places. htsjdk-rs decision 0025 measured the worst divergence between a
permissively-licensed `exp` and `Math.exp` at **1 ulp** — a relative difference near 2⁻⁵², about
twelve orders of magnitude below what rounding to two decimals discards. Byte-identity of the
*output* does not require byte-identity of every *intermediate*, and this annotation was blocked by
an assumption rather than by a licence.

It is an assumption on this side too until the corpus says otherwise, which is what these boxes are
for.

- [x] `NaturalLogUtils` on `jmath::strict_math::exp` (#97), 42 rows compared as raw bit patterns.
      **9 of the 55 values are exact by construction**: `logSumExp`'s accumulator starts at **1.0**
      because the maximum's own term is folded in as that 1, and `sum != 1.0` then skips the `log`
      entirely, so a one-element array, a maximum with everything else at `-Infinity`, and a
      difference large enough that `1 + exp(diff)` rounds back to `1` all return `maxValue`
      untouched. The suite computes that property from the inputs rather than listing labels. Every
      other value lands within **1 ulp**, and the test prints the worst it saw rather than only
      passing. The refusal is on the accumulator, not the inputs: it fires after the loop, on
      `sum`, so a `NaN` reaches it and a `-Infinity` does not
- [x] `alleleFractionsPosterior` and the `Dirichlet` under it (#98), 29 rows over twelve fixed-point
      runs. This was where the risk was, and it was two risks: amplification across iterations, and
      a **different iteration count**, since convergence is a threshold test on `distance1/sum` that
      a difference too small to see in the values can land on the other side of. The suite asserts
      the count exactly *before* comparing any value, because otherwise every value comparison
      measures the wrong thing, and the count is the engine's own — the harness replays the loop and
      checks its replay against the engine on every case. **Neither risk materialised: every
      iteration count matches and the worst value divergence is zero ulp.** The fixed point is
      bit-identical, not merely bounded. The `weights` rows are the control, `digamma` with no `exp`
      in it, asserted bit-identical for that reason
- [x] `AllelePseudoDepth` itself (#99), with the suite comparing the **formatted strings**.
      Comparing the doubles would re-measure the `exp` gap, which is already measured; comparing
      the strings measures whether it reaches the output. Two suites, both `golden-pending`.

      **`DecimalFormat` was not the ten-line job this line used to claim.** "HALF_EVEN to two
      places with trailing zeros dropped" is what the Javadoc says and it is not what the class
      does. Three facts had to be measured, over 5,699,818 formatted values: it rounds the
      **shortest decimal form** and not the value, so `0.1` at forty places is `0.1`; an apparent
      tie is settled by which side of the halfway point the double sits on, so `0.155` goes down
      and `0.165` goes up; and the two patterns **disagree with each other** where the rounding
      position falls before the first digit, because that is where the class's internal fast path
      stops applying. `0.005` and `5e-5` are the same shape one decade apart and go opposite ways.
      None of it needed transcribing, all of it needed measuring, and the port is exact for values
      below 2^53 with at most fifteen significant digits — decades away from anything this
      annotation emits

      Two findings in the annotation are not arithmetic at all, and neither is in any test the
      reference ships. `composePriorPseudoCounts` memoises one array per allele count and hands out
      **that array**; on the empty-evidence branch the posteriors are it, so the closing
      `posteriors[i] -= prior[i]` **zeroes the memo** and the next genotype with the same allele
      count gets a prior of zeros. The reference's own second answer is `NaN,NaN`. And the log10
      branch's visitor looks the mapping quality up at `evidence().get(row)` with `row` the
      **allele** index, so a site with more alleles than reads throws `IndexOutOfBoundsException`
      and one with fewer floors each allele's row using an unrelated read's quality

**No divergence.** All 46 calls match, all 69 rows, and the library goes to **53 of 54** — leaving
only `AssemblyComplexity` and its G3 dependency. The rounding argument this section opened with
turned out not to be needed: G1.9.2 measured the fixed point bit-identical, so nothing had to
survive the formatter. What the formatter did instead was cost a day, because reproducing it needed
three undocumented facts and 5.7M measured values.

The goldens are byte-identical to the same container on Apple Silicon, so nothing in this chain
depends on the host — which is a stronger statement than the suite needed to make and is worth
having on the record for a chain built on an unported `exp` and an unported `pow`.

**Below 2^53 the formatter reproduces the reference exactly**, on 903,121 measured values. It did
not at first: two of the four divergences were this port's fault rather than Java's, an equidistant
pair of shortest forms that the specification resolves toward the even digit and Rust's formatter
does not. That is fixed, in both repos (htsjdk-rs #74).

What is left is above the line, and it is not a rule. On a sweep of 493 such values Java gave the
shortest form 472 times, the double's exact value 9 times, and neither 12 times — `2^60` comes out
as its exact value rounded to eighteen significant digits. Those are branches inside Java 17's
pre-Schubfach `FloatingDecimal`, so closing them would mean transcribing GPL2 source or fitting an
implementation to measurements, and both are refused.

Two ways round were then examined and both came back **negative, by measurement**.

*A permissively-licensed implementation of the same algorithm.* htsjdk-rs decision 0013 listed this
and noted it had never been searched. It has now, and the answer follows from
[JDK-4511638](https://bugs.openjdk.org/browse/JDK-4511638): Java 17's behaviour **violates its own
Javadoc**, so any clean-room implementation implements the specification and therefore agrees with
this port rather than with the oracle. The concrete instance is Android's Apache-2.0 `RealToString`,
which produced `100000000000000000000000` for `1e23` where the reference produces
`99999999999999990000000` — the same two answers this suite records.

*A newer oracle JDK.* Decision 0013 said cause A disappears against a JDK 19+ oracle. That is right
for `String.format` and **wrong for `DecimalFormat`**, which is the path these divergences are on:

| JDK | `Double.toString(1e23)` | `new DecimalFormat("#.##").format(1e23)` |
|---|---|---|
| 17, the pinned oracle | `9.999999999999999E22` | `99999999999999990000000` |
| 21 LTS | **`1.0E23`** | `99999999999999990000000` |
| 22, 24, 25 | `1.0E23` | `99999999999999990000000` |
| **26** | `1.0E23` | **`100000000000000000000000`** |

The Schubfach rewrite replaced `Double.toString` and `Float.toString` and left `DecimalFormat` on
the old converter until **JDK 26**, released this year. Pinning the oracle nine major versions past
what GATK 4.6.2.0 is shipped against would make the goldens represent a runtime nobody runs these
tools on.

The cost of a bump was measured while the image existed, and it is not the obstacle: all **57
oracle-backed suites, 62 cases, 32,604 compared values** were replayed against a JDK 21 oracle and
**not one moved**. The obstacle is that the version which would help is too new to be the
reference. (The image's own probe refused JDK 21 outright — `java major is '21', expected '17'` —
so the pin cannot change by accident; it had to be relaxed in a copy to take the measurement.)

So the line stays where it is, and it is worth saying which side of it is right: **the port
implements the documented behaviour and the pinned oracle runs a version whose defect Oracle has
since fixed.** JDK-4511638 was opened in 2001 and closed in 19 and 26. Nothing this programme
produces comes within eleven orders of magnitude of 2^53 anyway — pseudo-depths are read counts and
fractions are in `[0, 1]`.

**Not in scope:** reproducing `Math.exp` bit-for-bit. It has no specification to implement
against — beyond "within 1 ulp and semi-monotone", its only definition is its own code — so
recovering its bits from black-box measurement would be reverse-engineering toward a functional
copy, which is a worse position than reading the source rather than a better one.

### Explaining the code as it is written

- [~] `docs/COMMENTING.md`: every item that is not self-evident answers what it computes, how, and
      **why it is written this way rather than the obvious way**, and Rust idioms are explained
      where they are not guessable from Java. The third question is the point: in a byte-identity
      port the obvious way is usually wrong, and a reader who can check the Java should not have to
      learn Rust to check the port. `tools/audit/comment_density.py --check` is a ratchet, not a
      floor: it fails only when a file already on the list loses its explanations
- [ ] the tranches, in the order a reader needs them: the engine's numeric primitives, then the
      likelihood matrix, then the annotations, then the readers and writers. **3 of 122 files** in
      this repository, 332 across the three

### G1.8 The argument layer

- [~] Barclay's argument model and validation at library level, so covering-array vectors are
      interpreted as upstream interprets them. The unified CLI dispatcher stays out of scope
- [x] the **value model**: `NamedArgumentDefinition` and the grammar under it, 368 rows over 41
      command lines. Six of its rules are not what the annotation names suggest. `optional()` does
      not decide optionality (`isOptional` is the annotation **or** a default that renders as
      something other than `"null"`, and an empty collection renders as `"null"`, so an
      initialised-but-empty `List` is required and an initialised-and-non-empty one is not);
      `"null"` is a value with three outcomes (clear a collection, throw on a non-optional
      argument, throw a *different* exception on a scalar whose **raw** field is primitive);
      `isValueOutOfRange` begins with `value == null ||`, so a null on a bounded numeric argument
      is out of range, and the message formats the bounds by the **value's** type, so one argument
      reports `allowed range [1, 10].` for a rejected `0` and `allowed range [1.0, 10.0].` for a
      rejected `null`; the **recommended** range is checked with that same method against the
      **hard** bounds, so its warning is unreachable for any non-null value; a scalar refuses a
      second occurrence rather than taking the last; and a collection is cleared before the first
      value unless the parser is in `APPEND_TO_COLLECTIONS` mode. The grammar under all of it is
      **jopt-simple's**, not Barclay's: `--name=value` is refused outright, a flag consumes a
      following token only when `StrictBooleanConverter` accepts it, and the positional-argument
      message is built by `Collectors.joining` given a delimiter where a prefix belongs
- [x] **tagged arguments and collection-file expansion**, 150 rows over 30 command lines. A tag is
      a rewrite that happens **before** the grammar runs: the option name is peeled off and the
      pair of tag and value is stored under a surrogate key built from the option string **and**
      the value, so the same option with the same tag and the same value twice is "duplicated on
      the command line" while the same tag with two different values is two values. A tag on a
      field whose type does not implement `TaggedArgument` is refused when the value is *set*, and
      the message is `getShortName() + "/" + getFullName()` with no guard, so an argument with no
      short name reports itself as `/plain-scalar`. Expansion is **collection-only**: the same
      `.list` path becomes three values on a collection, stays a path on a scalar, and stays a path
      on a collection declaring `suppressFileExpansion`; a tag is written onto every value the file
      produced
- [x] **`@ArgumentCollection` flattening**, 42 rows. This is how `-L`, `-XL` and the read-filter
      arguments reach a tool: not declared on it, but on collection objects it holds, flattened
      into one namespace where nothing records which object an argument came from. Two orderings
      fall out and neither is stated anywhere: `getAllFields` adds a class's own fields and *then*
      climbs to its superclass, so a subclass's required argument is reported missing before its
      base class's; and the recursion is depth-first **at the position of the field**, so a nested
      collection splices its arguments between the two it sits between. Three refusals happen
      before any command line exists, of which a duplicate alias is the one worth naming: it is a
      construction failure rather than a shadowing rule
- [x] **`--arguments_file`**, 76 rows over 14 command lines. The only argument that changes the
      command line rather than a field. The file's arguments come **first**, wherever
      `--arguments_file` sat, because the original command line is appended to the expansion
      rather than the other way round: a collection reads `[a, b, cli]` either way round, and a
      scalar given in both a file and the command line is a duplicate rather than an override. The
      recursion is bounded by a **set of file names**, not a depth, and every file *named* enters
      it including ones skipped for already being there, so a self-including file and a pair that
      include each other are each read once. The argument is not removed between passes: the
      recursion ends because the second expansion is empty
- [x] **plugin descriptors**, 48 rows, measured through GATK's own `GATKReadFilterPluginDescriptor`
      over GATK's own filters rather than a stand-in. Every read filter's arguments are in the
      parser on every run; `validatePluginArgumentValues` runs **before** the required check and
      **removes** each controlled argument that nobody set and whose filter nobody named, so a
      required argument of an unselected filter would not fire at all. Set without its filter, the
      same argument is an error naming the filter **class** rather than `--read-filter`, built
      from the short name, a slash and the long name with no guard. The refusal of an unknown name
      is the descriptor's own `validateAndResolvePlugins`, not the argument layer. **Discovery is
      not ported**: which definitions exist is a property of the filter library and belongs with
      the tool that has plugins
- [x] **scope closed.** The usage text is deliberately **not** in this box: it is help output, not
      the argument model, and no covering-array vector reads it. It belongs with the unified CLI
      dispatcher, which G1.8 excludes by its own wording, and it is listed under G2 with the tools
      that print it

---

## Milestone G2: the 163 non-Spark GATK tools, by archetype

- [x] `PrintReads`, byte-identical: six output BAMs and five `.bai` indexes, under the JDK
      deflater
- [~] the rest of `record-transform` (56 tools, the largest archetype). **The calibration gate is
      answered**: `UnmarkDuplicates` and `RevertBaseQualityScores` are ported, with one suite
      covering both. What the second and third members cost is not the transform — both `apply`
      bodies are two lines — it is what the archetype hides:

      * both **replace the default read filters** with `ALLOW_ALL_READS`, where `PrintReads` takes
        `GATKTool`'s default of `WellformedReadFilter`. Three tools in one archetype, and their
        default traversal is not the same set of reads;
      * `RevertBaseQualityScores` **aborts the whole run** on a read with no `OQ` — not skips, not
        passes through. A port that passed it through would emit a larger and healthier-looking
        file than the reference, which is the worst shape a divergence can take;
      * an **empty `OQ` is the same as an absent one**, because `getOriginalBaseQualities` returns
        null for both even though `fastqToPhred("")` returns an empty array happily. Measured, not
        inferred.

      The measured marginal cost: `PrintReads` needed 152 lines of harness and its own header
      logic; the second and third needed **95 and 235 lines of Rust between them and no new
      harness**, because extracting `sam_output` left the `@PG` handling, the ID suffixing and the
      writer shared, and `print_reads.rs` fell from 152 lines to 52 by delegating. So the
      archetype's 54 remaining members are bounded by their `apply` and their filter overrides
      rather than by the engine — which is what the gate existed to find out.

      Oracle-backed: **10 output BAMs and their indexes byte for byte, and both refusals
      reproduced** with the class and message the reference throws. The golden is byte-identical to
      the same container on Apple Silicon
- [ ] `reporting-walker` (56 tools)
- [ ] reference, interval, coverage, CNV/SV and genotyping-array utilities
- [ ] full parameter coverage per tool (covering arrays plus fuzzing), not the default path
- [ ] the **usage text** (`CommandLineArgumentParser.usage`), moved here from G1.8: it is help
      output rather than the argument model, no covering-array vector reads it, and it is printed
      by the tools rather than by the layer underneath them

---

## Milestone G3: the variant callers (individually multi-month)

- [ ] `HaplotypeCaller`, `Mutect2` and relatives: assembly graph, genotyping
- [ ] **PairHMM** targeting the pure-Java `LoglessPairHMM` semantics. The implementation must be
      pinned in the oracle contract, since `FASTEST_AVAILABLE` resolves per host

---

## Milestone H: finish the htsjdk-rs foundation

Everything downstream inherits these, so they are front-loaded.

- [x] BAM/SAM/tags/index read and write. **Measured before working, and the entry was stale.** The
      write-side BAI is oracle-backed over nine shapes (`empty`, `unmapped` with both a placed and
      an unplaced read, `window_boundary`, `all_levels`, and the rest); the read-side index has its
      own six, with the dump recording why the two forms differ — a chunk ending on a BGZF block
      boundary is `(nextBlockAddress, 0)` to the reader and `(blockAddress, blockLength)` to the
      writer. Tags cover every writable type, `B` arrays and the empty one included. **CRAI moved
      to CRAM**, where it belongs: it is the CRAM index and cannot precede CRAM
- [x] VCF and Tribble. Allele, variant, header, encoder, the Tribble index in both directions and
      both layouts, the whole-file read, the whole-file write with its index, and the typed
      attribute accessors. **Both named consumers are done**, so nothing in G1 waits on this entry,
      and the field types that closed it turned out to be a measurement rather than a port.

      The **Tribble index** landed with htsjdk-rs #83: **both** layouts are read byte for byte, and
      what the port does not reproduce is narrower than this entry once said — the comparator the
      interval-tree *query* sorts with, which htsjdk itself calls "a little cryptic" and which is
      not a consistent order at all: `compare(a, a)` is `-1` and blocks one byte apart compare
      equal. For any pair whose starts differ by two or more it is an ordinary ascending sort, and
      no index seen here contains such a pair. The dump was written before the
      port and earned that order three times. The type identifiers cannot be read out of the Java
      at all — `LinearIndex.INDEX_TYPE` reads a field of the `IndexType` enum whose own constructor
      is handed `LinearIndex.INDEX_TYPE` — so `1` and `2` are measured rather than cited. The bin
      width is **per contig**, 16000 and 8000 in one file, so a reader assuming the creator's
      default would answer every query wrongly and never fail. And the header carries the source
      file's modification time, so the raw bytes differ on every run: the golden masks those eight
      bytes and reports the offset rather than being quietly unstable.

      **Writing one landed with htsjdk-rs #99**, which is what GATK does beside every VCF it emits,
      and it is a different problem: reading is a layout, writing is a set of decisions the layout
      only records the outcome of. The per-contig bin width above now has its **cause**.
      `LinearIndex.optimize` doubles the width per contig, merging blocks pairwise, until the most
      dense block is estimated to hold more than a hundred features, or one block is left, or the
      width goes bad — and it keeps the *last* width still under the threshold rather than the
      first one over it. The estimate is the largest block's size in **bytes** over the mean bytes
      per feature, never compared to the feature count the same object carries, so two files with
      identical feature counts and different line lengths index differently.

      **The index type is chosen from the data**, and measured, the choice flips both ways on the
      same two files: sparse data gets a linear index under `FOR_SEEK_TIME` and an interval tree
      under `FOR_SIZE`, and dense data gets the opposite, with nothing on a command line to say
      which a run produced. And the header's `FEATURE_LENGTH_MEAN` is **not the mean feature
      length**: the statistics are pushed the running maximum at each step, so features of lengths
      10, 10, 600, 10, 10 write 364.0, the mean of 10, 10, 600, 600, 600. A port computing the
      honest mean writes a different file.

      **Both layouts are written** (htsjdk-rs #100). The first pass wrote only the linear one and
      justified that as symmetry with the reader; the justification was **false**, because the
      reader parses and queries both, so the interval tree was simply not done. It is now, and it
      was the harder half: `IntervalTreeIndex.ChrIndex.write` writes `tree.getIntervals()`, and
      that is a **pre-order walk of a red-black tree** rather than a sorted list, so the byte order
      is the order the rotations left the nodes in. Agreeing on those bytes means reproducing the
      CLRS insert, both rotations and the `min`/`max` update that walks to the root after each one,
      and the insert comparator sending **equal starts left** is part of the file.

      **The two halves are joined** (htsjdk-rs #101): the VCF and its `.idx` are written in one
      pass, which is what every GATK tool emitting a VCF does. Indexing is **on by default**, so a
      caller who asks for nothing gets an index and one who supplies no dictionary gets an exception
      rather than a file. The position recorded for a record is the one **before** it and is
      absolute in the stream, so **the header is counted**: the first record sits at 242 because
      that is the header's length, and an index built from a feature list that forgot the header is
      uniformly off by it. The sequence dictionary becomes `DICT:` **properties** rather than the
      flag it used to be, before the four statistics, and `flags` stays zero.

      And **the layout is not the caller's choice**: the writer always uses the dynamic creator with
      `FOR_SEEK_TIME`, so three of the six indexed files measured are interval trees — one record
      gets a linear index, two thousand records get a tree, and so does a header-only file, whose
      density is a division of zero by zero. That is why writing the tree was not optional

      **The other is done.** `VariantContextComparator` and `VCFUtils.smartMergeHeaders` are ported
      and oracle-backed (htsjdk-rs #81 and #82, 67 and 16 rows), so `MultiVariantDataSource` has
      what it needs and **the multi-input walkers G1.6 handed over are no longer blocked**.

      Between them the two dumps corrected the port **four times before a single golden row was
      compared**, and all four are the same shape — behaviour that reading the source carefully
      does not reveal. The comparator's two constructors word the *empty* case differently, and
      nothing in the class explains why. Every merge output carries a `fileformat` line **no source
      wrote**. The merge's version comes from a **field** set at parse time rather than from a
      `##fileformat` line, so a header assembled in memory never reaches the version policy at all.
      And both of the merge's Integer/Float promotion arms are **no-ops**: the Java says "promote
      key to Float" twice and neither arm does it, because the `put` writes back what the map
      already holds.

      **The whole-file read is in** (htsjdk-rs #98), and it is the loop the three existing slices
      were missing rather than a fourth slice beside them. What it measured is that the codec is
      stateful and the state is not visible in any single line. htsjdk-rs decision 0035:

      **the header a reader hands back is not the header the file contains.** The `VCFHeader`
      constructor deletes the file's own `fileformat` line and `getMetaDataInInputOrder` prepends a
      synthesized one saying `VCFv4.2` for everything below 4.3, so a v4.0 file reads back claiming
      to be v4.2. And `doOnTheFlyModifications` defaults to **true**, so for the eighteen IDs htsjdk
      holds a standard for, a count or type that disagrees replaces the whole line, description
      included; a description that disagrees on its own is kept, which is what makes the rewrite
      hard to notice.

      That rebuild re-attaches the version **only from 4.3 up**, so below 4.3 the codec knows the
      version and the header does not — and both consequences are load-bearing. `VCFWriter` refuses
      a 4.3 header by testing that field, so a v4.3 file **can be read and cannot be written back**;
      and `smartMergeHeaders` reads the same field, so the version policy above is unreachable
      through the read path, which is the other half of the note two paragraphs up.

      The **line counter** is shared and incremented in two places, so the same malformed line
      reports two different numbers: `Line 12` from the column check, which runs before
      `parseVCFLine`'s increment, and `line number 13` from `generateException`, which runs after
      it. A `#` line in the body increments nothing and decodes to null, making it a **silently
      dropped record** rather than a refusal. One upstream message is wrong and reproduced as it is:
      a sites-only file is checked against 8 columns and told it was expecting 9.

      **And what closed the entry was a measurement, not a port** (htsjdk-rs #102). The last open
      item read "full field-type coverage", which sounds like work per type; measured, **the
      declared Type converts nothing**. Every `Type` in the format stores a String, or a list of
      Strings when `Number` is not 1, with the single exception of `Flag`. `Integer`, `Float`,
      `Character` and `String` are indistinguishable in a decoded record. The Type only tells a
      caller which accessor to reach for, and the conversion happens there, once per call.

      Those accessors are where the surprises live, and they are the layer G1's annotations do
      arithmetic on. `getAttributeAsInt` tests missing with `==` on a String, so it is true only for
      the constant the codec assigns to a bare key and to `KEY=`; a value written `KEY=.` arrives as
      a substring and throws instead. Three spellings of missing, identical in every rendering, and
      the two outcomes are a number and an exception. `getAttributeAsDouble` has no such test at all
      while `getAttributeAsDoubleList` does, so two accessors over the same conversion disagree. A
      scalar accessor reaching a list or a flag fails the **cast** rather than the parse. And
      `parseVcfDouble` accepts `1f`, `0x1p3`, `" 1"`, `inf` and `nan`, so a VCF may carry numbers no
      reading of the specification predicts
- [~] **CRAM** (container model, all encodings, codec negotiation, reference-based compression),
      **and CRAI with it**: a sub-project on its own, built floor upward. Four suites are
      oracle-backed and the scope is smaller than 169 Java files suggested.

      **The floor is the integers** (htsjdk-rs #103). A container header is a run of ITF8s, so
      nothing above them can be checked until they are, and they lie in two directions. The
      five-byte ITF8 **stores four bits twice** and the reader takes byte four whole while masking
      byte five to its low nibble, so `f000000112` and `f0000001f2` both read 18 and a stream whose
      two copies disagree resolves silently. And a **truncated stream is a number, not a refusal**:
      `InputStream.read()` returns -1 past the end and nothing checks it after the first byte, so
      `80` reads -1. Adding the refusal a port would naturally reach for would be the divergence.

      **The file definition and container header** (htsjdk-rs #107). The file id is padded to
      exactly 20 bytes and **truncated in silence**. The checksum covers the header rather than the
      container, is little-endian, which is the opposite of the CRC in a BGZF block one crate away,
      and is **absent below version 3**, so the same container is four bytes shorter in a 2.1 file.
      Every file has at least two containers, and the last one is recognised by a magic number in a
      **coordinate field**: its `alignmentStart` is 4542278, which is `0x454F46`, which is `EOF`.

      **The block** (htsjdk-rs #108). The CRC covers the header **and** the content together, so a
      block cannot be verified without re-reading its own header and the four checksum bytes sit
      **outside** the `compressedSize` the header declares. A port that walks a container by adding
      header plus compressed size lands four bytes short and misreads every block after the first.

      **The scope is bounded by a measurement, not by the file count.** htsjdk 4.2.0 **ships the
      CRAM 3.1 codecs and refuses 3.1 files**: `isSupportedVersion` answers false, and opening a 3.1
      definition throws while a 3.0 one gets past the version gate. So rANS **4x8** is required and
      rANS Nx16, the range coder, fqzcomp and the name tokeniser are not, because no file htsjdk
      will open can contain them. That is 21 of the 169 files out of scope, and they are the hardest
      21. htsjdk-rs decision 0038.

      Measured on five ordinary files, a four-read CRAM uses **RAW, GZIP and rANS** over 29 blocks,
      which is what makes rANS 4x8 required rather than optional. Two external compressors are
      **unreachable from the oracle as it stands**: Commons Compress is not on its classpath, so
      bzip2 and LZMA blocks cannot be produced to compare against at all.

      **rANS 4x8 order 0 is ported and byte-identical** (htsjdk-rs #114), over 1782 rows and
      eighteen inputs. It is arithmetic rather than layout, and four of its properties are not in
      the specification. The **requested order is not always the written order**: below four bytes
      `compress` uses order 0 whatever the parameters say and the order byte records what it used.
      The four final states are written big-endian and the whole blob is then **reversed**, two
      reversals that cancel, so a port doing only one of them puts sixteen plausible bytes in the
      wrong place. The normalisation is fixed point, and **one symbol absorbs the whole rounding
      residue**: on an input holding symbol `i` exactly `i` times, symbol 254 normalises to 31 and
      symbol 255, whose fair share is the same 31, is written as **152**. And the frequency table's
      run marker is **inferred, never signalled**: the decoder peeks at whether the next symbol byte
      is the current symbol plus one.

      The suite compares the encoding symbols field by field, not only the output bytes, because a
      stream comparison says a port is wrong and a symbol comparison says which multiplication is.
      **Order 1 is the next slice**
- [x] **GKL-exact deflate**, and all nine levels reproduce GKL, by two routes with the boundary
      stated. Levels 3 to 9 are a pure Rust port: htsjdk-rs `crates/gkl-deflate` reproduces GKL
      byte for byte there, **28 of 28** (fixture, level) pairs against the column the real
      library produced in the pinned container.

      The scoping was measured rather than assumed, and it moved twice. Decision 0028 read the
      level table and inferred "levels 1 to 6 are igzip"; decision 0029 read the branch out of the
      library and found the split elsewhere. **Levels 1 and 2 are igzip; 3 to 9 are a zlib 1.2.13
      carrying Intel's `deflate_medium` patch**, which the JDK's zlib 1.3.2 disagrees with below
      level 7. The BGZF default is 5, so the default path is that patched zlib, not igzip. Two
      pieces reproduce it: `deflate_medium` itself, and the CRC-32C positional hash Intel
      substitutes for zlib's multiplicative rolling one.

      **Levels 1 and 2 are done too, by linking ISA-L rather than porting it**, and the reason is
      a measurement. Decision 0034: ISA-L ships *two* implementations of its own compressor,
      readable C and hand-written SIMD kernels, and they disagree — 19749 bytes where the assembly
      gives 19044, and GKL ships the assembly. A translation of the readable version would have
      been a confident wrong answer. Linking is also the trade this programme already makes for
      the JDK deflater, where htsjdk-rs decision 0001 pins `flate2` to a *vendored C zlib*.

      What makes that safe rather than convenient is a canary. ISA-L falls back to that same C when
      built without an assembler, on a CPU without SSE4.2, and on any architecture with no kernels;
      in all three states it returns **valid deflate that decompresses correctly**, so a round-trip
      test passes and a length check passes. The crate carries 2048 bytes GKL was given and the 694
      it returned, compresses the first and compares the second before answering, and **refuses**
      if they disagree. A green CI that skipped the comparison would have looked like one that ran
      it, so on x86-64 with SSE4.2 an unavailable igzip fails outright.

      The pure-Rust route for these two levels stays open and unchosen: it means reproducing about
      2,400 lines of SIMD assembly, since the readable version is known to produce different bytes.

      **Nothing here is licence-blocked**: ISA-L is BSD-3-Clause, GKL is MIT, Intel's zlib fork is
      under the zlib licence. The first byte-deciding component in the programme whose reference
      implementation is *permissively* licensed, unlike `Math.exp` where the obstacle is law
      rather than effort.

      **"Prove it on real x86-64" came first rather than last, and the answer is one bit.**
      Decision 0033: both backends cut at **SSE4.2** and nowhere else. Above it there is one
      behaviour, and SSE4.2, AVX and AVX2 hosts all produce it; below it a second, which zlib
      reaches by using the multiplicative hash and igzip by falling back to pure C. So a byte claim
      over BGZF is a claim about a CPU class, and the class every oracle run has been in is
      "reports SSE4.2". What this retires: 0028 treated "igzip's AVX2 kernels might differ from its
      SSE ones" as the thing that could put this entry beside `Math.pow`. Measured, they do not.
      AVX512 remains untested, because QEMU's TCG drops the feature bits and no available host has
      it.
- [~] **jmath**. The target is **not** "the corpus reaches 100%": its columns are `java.lang.Math`,
      whose remaining divergent functions can only be made exact by transcribing GPL2 source, so
      that is unreachable by construction. htsjdk-rs decision 0023 replaced it with "every function
      a ported call site reaches is exact, and every one that cannot be is named at the call site".

      **The rule was written down and the list was not.** It is now, in
      `docs/numeric-functions-a-ported-call-site-reaches.md`, built by walking the call sites rather
      than the library. Eight functions are reached through jmath and every one is exact or
      bounded; **six call sites reach the host libm directly**, and those are the whole of the gap:
      four `powf`, one `exp` in the activity profile that should be `strict_math::exp`, and
      `log1mexp`'s `ln_1p`/`expm1` pair, which nothing in G1 reaches. `StrictMath.exp` and
      `StrictMath.pow` are exact and `Math.exp`/`Math.pow` are each bounded at **1 ulp**
      (decisions 0025 and 0027).

      `BinomialDistribution` and SVD, the other two names this entry used to carry, are **not
      reached by anything ported**: both are Mutect2-family and therefore Milestone G3. The honest
      entry for them is "waits for G3" rather than "remaining" 
- [x] the BGZF surface that exists today is cross-checked on real x86-64. The `bgzf` and
      `bgzf-termination` suites run on `ubuntu-latest`, which is real hardware, so every BGZF
      golden is re-derived off emulation on every push; decision 0007's addendum measured it
      directly on an AMD EPYC host as well. The **igzip** half folded into the entry above, where
      its subject lives

---

## Milestone P: finish Picard (109 tools)

- [~] **MergeBamAlignment** (in flight): transfer, `PG` linkage, `NM`/`MD`/`UQ` and the
      whole-file coordinate-sorted producer are done and oracle-backed. Remaining:
  - [ ] merged-header construction (`@SQ` from the reference `.dict` with a `UR` absolute-path
        canonicalization rule, `@RG` from the unmapped BAM, `@PG` from the aligned BAM via the
        ported `SamFileHeaderMerger`)
  - [ ] paired mate-info, proper-pair and `ClippedPairFixer`
  - [ ] off-end-of-reference cigar clipping, `UNMAP_CONTAMINANT_READS`, adapter and overlap
        clipping
  - [ ] multi-hit selection (`MultiHitAlignedReadIterator` and the primary-alignment strategies)
- [ ] the remaining ~59 Picard tools, by archetype
- [ ] full parameter coverage on every ported Picard tool

---

## Milestone X: the four hard problems

- [ ] **Spark** (39 tools): establish bit-identity against the `--spark-master local[1]`
      single-partition oracle first, then prove any parallel path byte-equal to it in CI
- [ ] **ML inference** (CNNScoreVariants, NVScoreVariants, relatives): reproduce the exact
      kernels and accumulation order, validated layer by layer against dumped intermediate
      tensors, weights pinned by hash, quarantined and reported as bio-identical if a model
      proves intractable
- [ ] **CRAM**, if not closed under H
- [ ] **GKL-exact deflate**, if not closed under H

---

## Milestone GPU: accelerate without losing the byte

There is nothing upstream to port here. GATK 4.6.2.0 ships no GPU kernel: `PairHMM` resolves to
`AVX_LOGLESS_CACHING`, `AVX_LOGLESS_CACHING_OMP` or the pure-Java `LOGLESS_CACHING` through Intel
GKL, and the only CUDA-adjacent tool is `NVScoreVariants`, which is ML inference and already sits
in Milestone X.

A GPU path is therefore a **second implementation of a path this programme has already made
byte-identical**, and it earns its place only by producing the same bytes as the first one. The
gate is not "close", not "concordant", and not "within tolerance": same golden, same bytes, or the
kernel does not ship. That is what separates this from every existing GPU reimplementation of
GATK, which targets concordance and says so.

### GPU.1 What makes a kernel eligible

- [ ] a written eligibility rule, and a guard that enforces it per kernel:
  - integer and byte work with a fixed traversal order is eligible as it stands;
  - floating point is eligible only when the operation order, the rounding and the transcendental
    implementations are pinned to the CPU path's. In practice: no fast-math, no FMA contraction
    unless the CPU path contracts identically, no TF32 or tensor-core substitution, and a fixed
    reduction tree rather than an atomics-order-dependent one;
  - `StrictMath`-equivalent transcendentals have to be computed in software on the device, because
    a vendor libm is not the one jmath reproduces;
  - anything that cannot meet the above is quarantined and reported as bio-identical, exactly as a
    CPU quarantine is
- [ ] a determinism gate across **two different GPU architectures**, not one: a reduction that is
      deterministic at one warp size or occupancy can stop being so at another, and a single-device
      green run does not establish that

### GPU.2 The targets, in order of tractability

- [ ] **BGZF deflate and inflate**: integer work, block-independent, and already the hottest path
      in every tool that writes a BAM. Byte-identity is decidable here because the output is
      defined bit for bit
- [ ] **the read filters and the pileup floor**: integer and byte predicates over many reads, no
      floating point, and a fixed order. The cheapest place to prove the harness works end to end
- [ ] **`PairHMM`**: the reason anyone wants a GPU here. Byte-identity against `LoglessPairHMM`
      requires reproducing the pure-Java accumulation order, which is the same constraint the CPU
      port already carries, so the kernel is a transliteration rather than a redesign
- [ ] **the assembly graph and genotyping** (G3): larger, and floating point throughout
- [ ] **CRAM codecs** (Milestone H), whose arithmetic is integral and whose output is defined

### GPU.3 How it is verified

- [ ] every conformance suite runs twice in CI, CPU and GPU, and the goldens are compared with the
      same comparator: a GPU run is not a separate claim, it is the same claim on other hardware
- [ ] the oracle contract records the device, the driver and the toolkit version alongside the
      container digest, because a kernel's result is a property of all three
- [ ] a measured speedup per kernel, published beside its byte-equality, so an accelerated path
      that is not faster gets deleted rather than kept

**Scope note.** This milestone is optional and runs in parallel: nothing in G1, G2, G3, H or P
depends on it, and no tool's byte-identity claim may rest on a GPU path alone.

---

## Milestone S: speed, once the bytes are settled

**Every box is gated on a conformance suite existing for the path it touches.** Speed work on a
path with no golden is not optimisation, it is unmeasured change. Tracking issue #107.

Fifty-seven conformance suites say what the port *produces*. Not one says what it **costs**. That
is the right order — a fast wrong answer is worthless — but it has been followed far enough that
the second question is now worth asking, and asking it late is an advantage: it gets asked on paths
whose correctness is already pinned.

### The constraint that makes this unlike ordinary optimisation

A byte-identity port cannot buy speed the usual way, because the arithmetic is not ours to
rearrange:

- **no reassociation.** Floating-point addition is not associative, and the port transcribes the
  reference's summation order precisely because reordering changes the double. `MathUtils.sum`
  accumulates in index order; `sumArrayFunction` accumulates over reads in index order;
- **no FMA contraction**, unless the reference contracts identically. `a * b + c` fused is not
  `a * b + c` rounded twice;
- **no fast-math, and no faster transcendental.** `jmath` exists because the host libm is not the
  JVM's;
- **no reordering of collection traversal**, which reaches output bytes in several places already
  on the record.

So the wins available are not arithmetic. They are allocation, copying, memory layout, I/O, and
work that need not happen at all. That is a narrower space than a normal port has, and this section
says so rather than promising a number.

### What is actually in dispute

Folklore says a Rust port is faster than the JVM. The answerable version is narrower:

- **JVM startup** is paid once per invocation and dominates short runs, and GATK is very often
  invoked once per file. The largest and least interesting win, measured separately so it does not
  flatter everything else;
- **steady state**, after the JIT warms up, is where the claim is genuinely uncertain. The
  reference has had two decades of tuning on these paths;
- **memory** is where a byte-identical port can lose, because it clones where Java aliases.

### The boxes

- [ ] **S.1** (#108) a harness measuring both sides on the same inputs — wall clock, CPU time and
      peak RSS, cold start reported separately from steady state, ratios rather than absolute
      seconds, on the conformance fixtures so a benchmark cannot drift onto an easier case than the
      correctness claim covers
- [ ] **S.2** (#109) the byte-neutrality gate, enforced rather than intended. A perf PR is not
      exempt from its suites, the workspace forbids the flags that would let arithmetic be
      reordered, and the "slow on purpose" list is written down so nobody attacks the contract by
      mistake
- [ ] **S.3** (#110) a published baseline before any optimisation, with the predictions written
      **first** so the measurement can contradict them
- [ ] **S.4** (#111) the costs already on the record — the clone in `AllelePseudoDepth`, the exact
      decimal expansion in `decimal_format`, the fixed reduction order, the software transcendentals
      — each marked removable or not
- [ ] **S.5** (#112) the first targets, chosen from the baseline rather than from reading the code.
      Blocked on S.3 by construction

Milestone GPU is a different thing: a second implementation on other hardware. This one is about
the CPU path that already exists, and no bit-identity claim may weaken to buy a speedup.

## Milestone V: program-level validation and reproducibility

- [x] one conformance manifest per repository, with the CI generated from it and a guard that
      fails if the generated YAML is edited by hand
- [x] the per-tool status dashboard, generated from the 311-tool inventory so it cannot drift
- [x] the golden audit: every committed golden is declared, every declared golden exists, and a
      keyed suite whose golden collapses is refused
- [x] the provenance guard (decision 0014), running in all three repositories
- [~] covering arrays generated **and run** per tool. Two tools are measured against both the
      reference and the port on every CI run (`AddOATag` 0/9, `CollectAlignmentSummaryMetrics`
      0/16 at t=2), with the corpus recording each side of a mismatching row and the dashboard
      printing the fraction. Zero is the honest figure for binaries written for the throughput
      benchmark; what remains is the other 42 tools with a suite, and a port binary to run each
      array against
- [~] the coverage-guided differential fuzzer (running in picard-rs CI, seeded from the arrays;
      not yet in gatk-rs or htsjdk-rs)
- [~] determinism gates: same input twice, different `-Xmx` / `TMPDIR` / `LC_ALL`, two clean
      builds (running in picard-rs; not yet in the other two)
- [ ] per-run provenance (input hashes, reference tags, toolchain, container digest, GKL provider
      state)
- [ ] per-tool validation reports (inputs, records, byte-equal file count, argument coverage, t
      level, fuzzer branch coverage, quarantined fields)

---

## The critical path, in order

1. **Close G1.5**, the pileup floor, then `LocusWalker`. It is what the largest GATK archetype
   needs and the last big engine piece before tools can fan out.
2. **`FeatureDataSource`** (G1.3), which unblocks `VariantWalker` and every VCF-reading tool.
3. **Fan out G2 by archetype**, measuring the marginal cost of the second and third member of
   each archetype (the calibration gate) and re-sizing from real numbers.
4. Close the **htsjdk-rs** gaps that block downstream tools: full VCF, write-side BAI, then the
   jmath functions ported call sites reach. Not "the corpus to 100%": its columns are
   `java.lang.Math`, whose remaining divergent functions can only be made exact by transcribing
   GPL2 source, so that target is unreachable by construction (htsjdk-rs decision 0023).
5. Then the **callers** (G3), then the **hard problems** (X).
6. **CRAM** and **GKL-exact deflate** proceed in parallel, being self-contained in htsjdk-rs.
7. **GPU** (Milestone GPU) is off the critical path by construction: it accelerates paths that are
   already byte-identical, and a kernel that cannot match the CPU bytes is not merged.

The honest bottom line: the tool-by-tool fan-out is tractable and amortizes well; the schedule is
dominated by the engine, the callers, and the four hard problems. Full parameter coverage on
every tool, not the count of tools with a first slice, is what the 100% claim rests on.

---

## How this is tracked

- This file is the human-readable milestone list. Each `[ ]` is a checkable unit, ticked in the
  same commit that makes it true.
- The per-tool status dashboard is generated from `tools/inventory`, so it cannot drift from the
  311-tool ground truth. This file is the plan; `docs/STATUS.md` is the measurement.
- Progress lands as one commit per slice (one ported symbol or feature), CI-gated on the
  digest-pinned oracle, merged only when green.
