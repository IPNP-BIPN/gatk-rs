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
| **picard-rs** | 121 tools | ~50 tools have a first slice, ~43 with an oracle-backed conformance suite; many are partial (default paths only). The harness is generated from a manifest, the fuzzer and the determinism gate run in CI, and argument coverage is measured for 2 tools |
| **gatk-rs** | 190 tools | 6 crates, **58 conformance suites, all oracle-backed**; 3 tools byte-identical, and the annotation archetype opened with 53 of 54 annotations measured. **No performance number exists yet for any of it** — see Milestone S |

Totals from the generated inventory (`tools/inventory`): **311 tools** (190 GATK-origin,
121 Picard-origin), **39 Spark**, ~13,130 arguments. Non-Spark: 151 GATK + 121 Picard.

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

The single biggest unlock: 151 non-Spark GATK tools stand on it.

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

## Milestone G2: the 151 non-Spark GATK tools, by archetype

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
- [x] **CRAM** (container model, all encodings, codec negotiation, reference-based compression),
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

      **Order 1 is ported too** (htsjdk-rs #117), 876 rows over sixteen inputs, and it is **not
      order 0 with more tables**. Its frequency table **counts three bytes that are not bigrams**:
      each of the four lanes starts with context 0 while only one of them follows nothing, so
      `calcFrequenciesOrder1` counts the byte at each quarter boundary as if it did. Measured on the
      four bytes `ACGT`, context 0 holds A, C, G and T at 1024 apiece, so the table says C is as
      likely as A after nothing and the input never shows it. Its normalisation is **floating
      point** where order 0's is fixed point, so one class carries two arithmetics. And **a
      frequency byte of zero means 4096 on the way in**, a reader-only rule that no input can
      produce, measured on a five-byte table built by hand.

      Both orders were byte-identical on their first run against the reference's own output.
      Measured aside: the x86-64 runner's dump is byte-identical to the same dump under emulation
      on Apple Silicon, so the rANS arithmetic does not depend on the silicon. That does not change
      where a golden may come from.

      **The compression header is open, and two of its three maps are done.** It is the RAW block
      cram-block measured behind the GZIP header block in every file, and it is three
      length-prefixed maps in a row.

      **The preservation map** (htsjdk-rs #118). Its size field is a **hardcoded 5**, not a count:
      `internalWrite` writes the literal and then exactly RN, AP, RR, SM, TD in that order whatever
      the header holds. The order is htsjdk's rather than the specification's and the reader accepts
      any order, so it is **invisible to a round trip** and visible only against the reference's own
      bytes. A boolean is `== 1` and not `!= 0`, so a 2 reads as false and nothing complains. An
      unknown key is a plain `RuntimeException` while a missing mandatory key is a `CRAMException`,
      and the latter names both keys whichever one is absent. And the tag dictionary's first group is
      always empty, present even where every record carries tags. Measured aside: the whole
      compression header is the same 160 bytes across four files differing in record count and read
      length, and only tags move it.

      **The substitution matrix** (htsjdk-rs #119), which is the five bytes the preservation map
      carries under SM. They pack four two-bit codes per reference base, and the codes are **ranks
      by observed frequency**, so the commonest substitution gets the shortest ITF8. The ranking has
      an overflow: the comparator is `(int) (o2.freq - o1.freq)`, a **long difference narrowed to an
      int**, so two frequencies whose difference is a non-zero multiple of 2^32 compare equal.
      Measured, with `C` substituted 4294967296 times and nothing else substituted at all, reference
      base `G` packs 27 in which **C ranks second behind a substitution that never happened**; one
      more occurrence packs 75 in which it ranks first. The sort also runs **twice**, the second
      pass over zeroed frequencies, so it is a sort by alphabet wearing the frequency comparator's
      clothes. TimSort's small-array path is ported rather than delegated to Rust's sort, because
      that comparator is not a total order and two sorts may legitimately disagree on one.

      **The data series encoding map** (htsjdk-rs #120), the second map, and the first thing in
      CRAM that describes a record rather than a container. **Its size is a real count where the
      preservation map's is the literal 5**: three maps in one header and two counting conventions,
      so a port that hardcodes both or computes both is wrong on one. **The write order is the
      enum's ordinal order**, not the alphabetical order the constructor populates in, and only the
      first is in the bytes. htsjdk **writes 26 of the 32 series**, and **reads `TC` and `TN` only
      to drop them**, so the reader's map can hold fewer entries than the file declared. The content
      ids are htsjdk's rather than the specification's. And an unknown encoding id is an **array
      index**, not a CRAM error, read **signed**, so an id byte of 255 arrives as index -1.

      **The tag encoding map** (htsjdk-rs #121), the third, which closes the header. **The key is
      the tag itself**, two name bytes and the type packed into twenty-four bits, so the type is
      part of the key and one name at two types is two entries. **The write order is that key's**:
      measured, two files whose records introduce the same three tags in opposite orders produce
      **byte-identical** maps, so the order the data arrived in is not in the file. Its collision
      guard covers ids 1 to 32 while the smallest printable tag packs to 2105376, which makes it
      real code no input can reach. And one finding that is not about this map but is only visible
      through it: **htsjdk narrows an integer attribute to the smallest type that holds it**, so
      `NM` at 1 to 4 is written `NMc` and at 100000 is written `NMi` from the same Java `Integer`,
      which changes the key and therefore the map.

      **The slice header** (htsjdk-rs #122), the last frame before CRAM becomes reads. **Its block
      count does not count the header block**: `1 + numberOfExternalBlocks`, which measured equals
      exactly the blocks that *follow* it, so a reader counting the header among them stops one
      short. **Six tags ride in the header and four of them digest nothing**: `B1` and `S1` a SHA-1,
      `B5` and `S5` a SHA-512, and on an unmapped slice all four are the digest of the **empty
      string**, identical in every file, which is **168 bytes of constant per slice**. Only `BD` and
      `SD` move with the reads, and two files differing only in whether their records carry tags
      have equal `BD`, so record tags do not enter the slice digest. **The tag section carries no
      length** and is read to the end of the block, so a header with no tags is indistinguishable
      from one whose tags are zero bytes long; it is ported as opaque bytes, which is what
      byte-identity needs. The MD5 is sixteen zeroes rather than absent, and an absent embedded
      reference is -1, whose ITF8 is the five-byte form: the commonest value of that field is its
      longest encoding.

      That finding came from correcting a misreading. The tag bytes were decoded by eye as one tag
      and are six; the dump was changed to decode them with the reader's own codec.

      **Ten CRAM suites are oracle-backed.** The last five are not layout but implementation
      decisions the specification does not prescribe, and **none of the five is visible to a round
      trip**: a hardcoded count beside a real one in the same header, a write order that comes from
      a `TreeMap` rather than from the data, a comparator that loses 32 bits, a tag type derived
      from the magnitude of a value, and 168 bytes of digests of nothing.

      **The read side of the tag codec** (htsjdk-rs #123), which is what that opaque section waits
      on. The write side was pinned five slices ago; reading is **not its inverse**. **Every narrow
      integer widens to one type**: `c`, `C`, `s`, `S` and `i` all come back as a Java `Integer`,
      so the width the file chose stops existing on the way in and the type a rewrite picks comes
      from the value alone. **Exactly two forms break the round trip**, and both are ones htsjdk
      never writes, so only a foreign file carries them: an `I` holding 5 is rewritten as `c`, four
      bytes shorter, and an `H` is rewritten as a `B` array a byte longer. **There is no in-memory
      `H` at all**: it decodes into the same `byte[]` a signed `B` array does, so every branch that
      would write one back is dead and `TextTagCodec`'s says so in a comment. `'A'` is a signed
      byte cast to a char, so `0xE9` becomes `U+FFE9` and the character in memory is not the one in
      the file, though the bytes survive because the write truncates it back: the same cast the
      substitution matrix goes through. The `I` range check **cannot fire**, the value being masked
      to 32 bits before it is compared against the 32-bit range. The unsigned flag of an array is
      **only the case of its type letter**, the elements staying signed, so a `C` array holding
      `0xFF` comes back as `-1`. And a repeated tag replaces rather than duplicates, the last one
      winning.

      Two things in already-merged code were wrong and are fixed by that measurement. The text
      encoder wrote `XX:H:48656C` back as an `H` line the reference never emits; it emits
      `XX:B:c,72,101,108`. And the record codec carried a **second decoder** for the same bytes,
      raising messages of its own invention, which is how the `H` handling came to be wrong in
      three places at once; ninety-nine lines of it are gone.

      Also closed: the candidate-golden artefact was named by the job's **position in the CI
      matrix**, so adding a suite renumbered it and one slice's golden took the name the previous
      slice's had. It is named after its suites now.

      **The read features** (htsjdk-rs #124), the first half of the record model and the point
      where CRAM stops being frames. A slice's records are not bases and a cigar: they are an
      alignment start, a read length, and a list of one-letter features, and **everything that
      matches the reference is stored as nothing at all**. That is where the compression comes
      from, and it makes every rule below a rule about what counts as a match.

      **The positions are one-based, and the interface says they are not**: every construction site
      passes `zeroBasedPositionInRead + 1` while `ReadFeature.getPosition`'s javadoc says
      "zero-based position in the read", so a port that believes the documentation is off by one on
      every feature it writes. **An insertion of n bases becomes n features and a soft clip of n
      becomes one**, decided five lines apart in the same loop: htsjdk's own comment says the
      insertion should use a `Bases` feature and does not, because that would need a
      `ByteArrayLenEncoding` and therefore a frequency distribution over lengths, so two of the
      twelve features are read and never written. **A mismatch splits on the alphabet, not on the
      cigar**: ACGTN against ACGTN is a substitution, anything else is a `ReadBase` carrying the
      quality score a second time. **The cigar's own claim is not consulted at all**, so an `X` over
      bases that match emits nothing and an `=` over bases that differ emits substitutions; the
      operator only says how far to walk. **Past the end of the reference every base is compared
      against `N`**, so an `N` out there matches and is stored as nothing: four `N`s at the end of
      the reference produce three features, not four. **`SEQ="*"` manufactures `N`s** that then
      mismatch like any other base, one substitution per position. And **the missing-quality test
      is an identity test**: `baseQualities.equals(NULL_QUALS)` is `Object.equals` on a `byte[]`, so
      an equal but distinct empty array takes the other branch and indexes it, which is measured as
      `Index 3 out of bounds for length 0` on a record htsjdk itself will hold.

      `htsjdk-cram` took its first dependency here, on `htsjdk-bam`. CRAM is defined in terms of SAM
      records, and the cigar walk was already ported.

      **The cigar rebuilt from those features** (htsjdk-rs #126), which is the way back. The
      interesting part is that **the cigar is not stored anywhere**: it is rebuilt from the feature
      positions and the read length, and `gap = position - (lastOpPos + lastOpLen)` is the only
      source of `M` in the output. The matches thrown away on the way in **come back as the gaps
      between what was kept**.

      **A substitution and a `ReadBase` are both `M`**, so the rebuilt cigar never emits `X` or `=`:
      a record written with `8X` comes back as `8M`, and over thirteen round trips that and `8=` are
      the only two that change. **A feature that consumes no read bases winds the read cursor back**,
      the bookkeeping being in read space, so a deletion at the first position leaves the whole read
      after it and `D@1 len=2` over a read of eight rebuilds as `2D8M`. **The switch silently ignores
      what it does not name**, dropping `BaseQualityScore`, `Scores` and `Bases`, and `Bases` carries
      read bases: a list holding one rebuilds as though it were not there. **The read length says
      where the read ends and it wins**, absorbing a feature positioned past it, and a read length of
      0 takes the accumulated length instead. Two guards in the reference cannot fire: a null feature
      list reaches the same single `M` through the empty-list check at the end, and the operator it
      compares against null never is.

      **The bases restored** (htsjdk-rs #128), which closes the reverse direction: a record can now
      be read whole. The bases come from **three sources at once**, the features saying what
      differs, the reference supplying everything else, and the substitution matrix turning a code
      back into a base.

      **It is two passes, and the second one overwrites.** `ReadBase` and `Bases` are skipped in the
      main loop under a comment saying to defer them, then applied straight into the array, so they
      win over the reference fill **and** over what the features before them wrote: an insertion of
      `GGG` at position 1 followed by a `Bases` of `TTT` at position 1 restores `TTT`. **The trailing
      fill stops at the end of the reference** and leaves the array's own zeros, so nothing ever
      writes an `N` past the end: a zero becomes one on the way through the lookup.
      **`toBamReadBasesInPlace` indexes a 127-byte table with a signed byte**, so a base of `0xE9` is
      index -23 and a base of 127 is one past the end, and because the table is built by adding 32 to
      every BAM read base, **`]` is an `=`** rather than an `N`. A substitution resolves against a
      **normalized** reference base, so an IUPAC code there resolves as though it were `N`. And
      unknown bases or a read length of zero return the empty sequence with the features not looked
      at at all.

      The comment above that loop says read features are 0-based. They are one-based, which the
      forward direction measured from the other side, and this is the function that consumes them.

      One thing about the corpus rather than the code: its reference is **aperiodic on purpose**.
      Against `ACGTACGT...` a reference cursor off by four reads the same bases as a correct one,
      and every case in the suite would pass anyway.

      **Fifty-three of fifty-four suites are oracle-backed**, thirteen of them CRAM. The record model
      is now pinned in both directions; **what remains for H.3 is the encodings that carry it**:
      the data series codecs the encoding map names, and the slice blocks they are read from

      **The codecs that carry it** (htsjdk-rs #129, #132, #133, #135, #137, #139), ten slices that
      close that gap. Every encoding identifier the map can name now has a codec behind it, and
      every one was measured before it was written.

      **The bit stream** under the three core codecs: bits go in most significant first, and the
      flush pads the partial byte with zeros **on the right**, so the padding is data as far as
      anything downstream can tell. A stream of one `true` bit and a stream of `0x80` in eight bits
      are the same byte; only the count of values expected says where the stream ends.

      **Beta, Gamma and Subexponential**, of which only Beta has an upper bound. Gamma and
      Subexponential derive a bit length from `Math.log(v) / Math.log(2)`, so the bytes written
      depend on a double division landing on the right side of an integer; the corpus walks every
      power of two to `2^31 - 1` and the runner and an Apple Silicon laptop agree on all of them.

      **Canonical Huffman**, where a file carries no tree at all: an alphabet and one code word
      length per symbol, from which both sides rebuild the same codes. **Byte symbols sort signed**,
      so `0x80` takes the first code word. **A one-symbol alphabet has length zero** and writes no
      bits, leaving a core block from which the number of symbols written cannot be recovered. The
      overflow check counts **set bits** rather than width, so three symbols at one bit are accepted
      and the third given code word 2. And an unmatched code word runs off the end of a table sized
      to the largest code word, which makes the codec's own "unable to map" message reachable only
      with an empty alphabet.

      **The external codecs**, where two codecs naming the same content id share one block and
      interleave in it. **External byte cannot see the end of its block**: past the end it returns
      -1, which a byte of `0xFF` that is really there also produces. **Byte array stop trusts the
      data**: written `01 00 02` with a stop byte of zero, the block reads back as one array of
      `01`, and nothing reports it. And **ByteArrayLen cannot wrap ByteArrayStop for reading**
      though the format allows it: the length is read, then `read(length)` throws.

      **Golomb, Golomb-Rice and Golomb-Long**, experimental in the reference and reachable from the
      factory, so a port that skips them cannot claim to read every legal file. **Golomb does not
      round-trip a value whose offset sum is negative and does not say so**: with `m` 4, `-1` is
      written `60` and read back `3`, because the quotient is written by counting up to it.
      **Golomb-Rice's parameter is not `m`**: the encoding calls it `m` and hands it to the codec as
      `log2m`, so one built with 8 divides by 256. And the divisor Golomb refuses is one Golomb-Rice
      takes without a word.

      **The encoding factory**, forty lines of Java with **one missing `break`**. Only the `BYTE`
      arm ends in one, so an `INT` that matches nothing falls into the `LONG` arm and then into the
      `BYTE_ARRAY` arm: an `INT` data series named with `BYTE_ARRAY_LEN` gets a byte array encoding
      rather than the refusal the method's last line promises. The suite is exhaustive by
      construction, four types by ten identifiers, so a new identifier upstream fails it on the row
      count.

      **The compression header, the record's flag words, a slice's blocks and the record reader**
      (htsjdk-rs #141, #143, #145, #147), which take those codecs up to a whole record.

      A compression header **read and written again is byte-identical** in all six measured cases,
      and a version 3 block carries a four-byte CRC-32 that a 2.1 block does not, **outside** the
      compressed size. A record carries **three flag words**, and two bits live in two of them at
      different positions: mate unmapped is `0x2` in the mate flags and `0x8` in the BAM flags.
      Restoring mate info **walks a ring**, and the template length is computed once and negated
      once, so the middle record of a triple keeps the zero it was built with. A slice's blocks are
      written **by content id and not by insertion**: added 3, 2, 1 the reference writes 1, 2, 3.

      **Reading a record** is the first suite where the port reads bytes **the reference wrote from
      records** rather than bytes a dump built by hand. Every field comes from its own data series
      and the series share streams, so the read order is not an implementation detail: read two in
      the wrong order and both come back wrong with nothing to say so. The alignment start is a
      delta and **the delta may be negative**; an unmapped record does not read zero read features,
      it does not consult the series at all; and a multi-reference slice reads a reference index per
      record where a single-reference slice takes the slice's own.

      That suite earned its keep on its first run: it found the port had left the **mapping quality
      series** out of the table it resolves encodings from. Nine other series would have masked it;
      the tenth did not.

      **Writing a record back** (htsjdk-rs #152), the other half of that round trip, and the
      correction it forced (htsjdk-rs #150). Porting the writer showed what the reader was missing:
      **an unmapped record's bases**, read one at a time, and the quality scores a record keeps as
      an array. The suite could not see either, because the corpus recorded neither, so a reader
      that skipped both passed it.

      Two rows fix that, and the pairing is the point. **An unmapped record followed by a mapped one
      is the only arrangement where the gap shows**: the unmapped record's bases come out of the
      same series the next record's read feature reads from, so skipping them makes the following
      record read the wrong byte and nothing says so. The other row is a record with
      `CF_QS_PRESERVED_AS_ARRAY` set, the only one that touches the quality score series as an array
      at all, through **the reference's oddest reader**: it takes the QS series' own encoding
      descriptor and hands it a data series type of `BYTE_ARRAY` instead of `BYTE`, which is the
      same external block read a different way, and htsjdk builds it by hand with a comment
      wondering why.

      The writer itself keeps three things as measured. **Two features can be read and not
      written**: `Bases` and `Scores` fall to its default arm, which throws, while the reader has
      branches for both. **A substitution whose code is negative** is resolved against the
      compression header's substitution matrix on the way out. And **an unmapped record writes its
      bases one at a time** with its quality scores nested inside that branch, the same asymmetry
      the read side has.

      Both directions share one corpus of eight slices, and each record row carries the record in
      full, feature payloads included, so the port rebuilds the input rather than taking it from a
      label. The blocks have to come out byte for byte what the reference produced.

      **The CRAM index** (htsjdk-rs #154), which closes what H.1 handed over. A `.crai` is not a
      structure, it is a sorted text file: one line per slice, six tab-separated integers, gzipped.
      So what is worth pinning is not a layout but four decisions. **Unmapped-unplaced sorts last**
      whatever its alignment start says, and its start is not consulted at all. **An unmapped entry
      never intersects, not even with itself**, which is a special case in the code rather than
      something the arithmetic produces. **The overlap test is a midpoint comparison**,
      `|a0 + b0 - a1 - b1| < span0 + span1`, which is not the expression `a0 < b1 && a1 < b0` and
      **does not agree with it on a zero span**: two identical entries of span zero do not
      intersect. And **a query with a start or a span below one matches the whole sequence**, so 0
      and -1 are a wildcard rather than an empty range, which is why an unmapped query finds
      unmapped entries that would refuse to intersect.

      **Codec negotiation** (htsjdk-rs #156), the writer's side of the choice. The reader takes what
      a file names; the writer chooses, and the choice is what makes one CRAM of a set of records
      rather than another. **The data series table is fixed** rather than derived, and **six of the
      thirty-two series are not in it at all**, so a reader that expects every series to be named
      finds nothing for them. **The compressor is chosen by running all three**, GZIP and rANS at
      both orders, and **the tie-break is the order of the comparisons** rather than of the
      compressions: a thousand identical bytes compress to 29 under GZIP and 29 under rANS 0, and
      rANS wins. **A tag of one value size gets a zero-bit Huffman length**, which writes no bits at
      all: the size lives in the encoding rather than in the data. **A `Z` of several sizes gets a
      stop byte of TAB**, chosen rather than searched for, so a `Z` whose text contains a tab is
      split by its own encoding, while a `B` over a hundred bytes searches for a byte its data never
      uses. And **two records whose tags differ only in order share one dictionary entry**.

      That measurement corrected the port a second time: `name3BytesToInt` packs a tag id **high
      byte first**, and the record reader had it the other way round. No corpus reached it, because
      no record measured until then carried a tag.

      The gzip length is the one thing here the port cannot produce: htsjdk compresses with the
      JDK's `Deflater`, whose output length is its zlib's business. So the rule takes the three
      lengths rather than the data, and every row carries all three. Measured, the runner and the
      laptop agree on those lengths too.

      **A whole file** (htsjdk-rs #158), which is the piece that shows the rest fit. The definition,
      the SAM header container, the container header, the compression header, the slice header and
      the slice's blocks, in that order and at those offsets, over `ce#5.2.1.cram` from htsjdk's own
      test resources with every block inflated. **The first container is not like the others**: it
      holds the SAM header in a `FILE_HEADER` block rather than a compression header, and a reader
      that treats it as an ordinary container is refused with the compression header's own message,
      which is the first thing the walk hit. **The EOF container is a container**, parsing as one
      whose record count is zero rather than as a byte pattern. And **the sort order htsjdk reports
      is not in the file**: it says `unsorted` for a header carrying no `SO` tag at all.

      `htsjdk-cram` takes `flate2` there, for inflate only: deflate output depends on the zlib
      behind it, and a block's GZIP content has exactly one correct expansion. bzip2 and LZMA are
      legal in a CRAM and are not ported, so a file using either is refused rather than read wrong.

      **H.3 is closed. Sixty-seven of sixty-eight suites are oracle-backed, twenty-seven of them
      CRAM.** The container model, every encoding, codec negotiation, reference-based compression
      and the index are all measured and ported, and a real file walks end to end.
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

## Milestone P: finish Picard (121 tools)

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
- [ ] the reference version moves to GATK 4.7.0.0, htsjdk 5.0.0 and Picard 3.5.0, once G2 and P
      are closed (see below)

---

## The reference has moved on, and this one deliberately has not

GATK [4.7.0.0](https://github.com/broadinstitute/gatk/releases/tag/4.7.0.0) was released on
2026-08-18. It pins **htsjdk 5.0.0**, **Picard 3.5.0**, GKL 0.9.1 and commons-beanutils 1.11.0.

**The target stays at 4.6.2.0 until G2 and P are closed.** 216 suites are oracle-backed against
4.6.2.0, and moving the target does not adjust them, it re-opens them: every golden has to be
re-measured on a new oracle and either survives the bump or is replaced with the difference
explained. The three ports are also coherent only because they name one set of pins, and htsjdk
4.2.0 to 5.0.0 is a major bump, so htsjdk-rs would have to land first. Finishing against a frozen
target and moving once is cheaper than chasing a moving one.

What will be waiting, in short. The full delta, per pull request, is in the tracking issue.

| Where | What changes |
|---|---|
| GATK 4.7.0.0 | one new tool, `ConvertCountsToDepthFile`, so the inventory's 311 becomes 312; behaviour changes in `SVConcordance`, `SVStratify`, `SVAnnotate`, `PrintSVEvidence`, `CollectSVEvidence`, `HaplotypeCaller`, `GenotypeGVCFs` and `Funcotator`, all still open; `--output-cram-version`, defaulting to 3.1 |
| htsjdk 5.0.0 | `jlibdeflate` becomes the default DEFLATE engine, so every BGZF byte depends on a pinning that must be re-verified; CRAM 3.1 write support with a trial-compression codec choice; `SAMRecord.toString()` returns the full SAM line; three fixed bugs that this port reproduced faithfully and would now have to stop reproducing |
| Picard 3.5.0 | `MarkDuplicates` physical location moves from `short` to `int`, changing optical-duplicate counts; `FilterVcf` plumbs `CREATE_INDEX`; fixes in `RevertSam`, `CollectRnaSeqMetrics`, `MergeBamAlignment` and `CrosscheckFingerprints` |

The deflater line is the one to watch. The dumps already pin the factory and print
`deflater\t<class>` into every golden, which is the guard put there after a Picard call silently
replaced the factory and the goldens after it were GKL bytes. It is what makes this bump
survivable at all, and on the move it must be re-verified against the new default rather than
assumed to still hold.

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
