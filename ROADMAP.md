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
| **gatk-rs** | 202 tools | 5 crates, **33 conformance suites, all oracle-backed**; 1 tool byte-identical, and the annotation archetype opened with 11 of 54 annotations measured |

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

The single biggest unlock: 163 non-Spark GATK tools stand on it. This is the active work.

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
      settings). The suite pins what a tool sees and does **not** distinguish the cache from a
      fresh query per call, which the manifest states. The **BED** and **IntervalList** codecs are
      ported and oracle-backed in htsjdk-rs. Still missing: the VCF codec and the Tribble index,
      which are htsjdk's and belong in Milestone H

### G1.4 Interval arguments

- [x] `-L` / `-XL`, padding, `--interval-set-rule`, `--interval-merging-rule`, subtraction
      (measured through `IntervalWalker`, 24 combinations)
- [x] interval **files**: `.list` and `.intervals` (lower-cased extension test, blank-line
      skipping, the empty-file refusal, and the order of the four tests in
      `parseIntervalArguments`), 13 arguments compared
- [~] the Feature-file path: `.bed` and `.interval_list` are closed. Both codecs landed in
      htsjdk-rs with their own oracle-backed suites (90 and 133 rows) and are registered in the
      seam, so all 19 measured `-L` arguments now resolve identically, four of them new rows for
      the two dictionaries an interval list is validated against. VCF arrives with the Tribble index
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
- [ ] the multi-pass and multi-input walker variants

### G1.7 Annotations

- [~] the 54-annotation library. **Eleven** are ported and oracle-backed: the counting family
      (`ChromosomeCounts`, `SampleList`, `RawGtCount`, `Coverage`, `MappingQualityZero`,
      `CountNs`, `OriginalAlignment`) and the median family (`BaseQuality`, `MappingQuality`,
      `ReadPosition`, `FragmentLength`, i.e. MBQ/MMQ/MPOS/MFRL), together with the
      `InfoFieldAnnotation` interface and the machinery underneath it (`AlleleList`/`SampleList`
      and their permutation, the `AlleleLikelihoods` matrix and its best-allele search, and
      `VariantContextGetters`)
- [ ] the remaining 43. **The claim that most of them wait on jmath does not survive a grep**: 10
      of the 57 files in `tools/walkers/annotator` mention `MathUtils` at all, and what the
      annotators reach through `java.lang.Math` is `log`, `log10`, `sqrt` and `round`, all four
      already exact. What they actually wait on is engine machinery: the rank-sum and
      strand-bias statistics, the pileup, and the genotype likelihoods. The median family landed
      once commons-math3 `Percentile` was ported, which htsjdk-rs decision 0023 made possible by
      separating an Apache 2.0 source from a GPL2 one

### G1.8 The argument layer

- [ ] Barclay's argument model and validation at library level, so covering-array vectors are
      interpreted as upstream interprets them. The unified CLI dispatcher stays out of scope

---

## Milestone G2: the 163 non-Spark GATK tools, by archetype

- [x] `PrintReads`, byte-identical: six output BAMs and five `.bai` indexes, under the JDK
      deflater
- [ ] the rest of `record-transform` (56 tools, the largest archetype)
- [ ] `reporting-walker` (56 tools)
- [ ] reference, interval, coverage, CNV/SV and genotyping-array utilities
- [ ] full parameter coverage per tool (covering arrays plus fuzzing), not the default path

---

## Milestone G3: the variant callers (individually multi-month)

- [ ] `HaplotypeCaller`, `Mutect2` and relatives: assembly graph, genotyping
- [ ] **PairHMM** targeting the pure-Java `LoglessPairHMM` semantics. The implementation must be
      pinned in the oracle contract, since `FASTEST_AVAILABLE` resolves per host

---

## Milestone H: finish the htsjdk-rs foundation

Everything downstream inherits these, so they are front-loaded.

- [~] BAM/SAM/tags/index read and write (core done; write-side BAI, CRAI and corner cases remain)
- [~] VCF and Tribble (allele, variant, header, encoder exist; full read, write, index and all
      field types remain)
- [ ] **CRAM** (container model, all encodings, codec negotiation, reference-based compression):
      a sub-project on its own
- [ ] **GKL-exact deflate** (ISA-L / igzip byte-exact), the default non-JDK path. Until it
      exists, every byte claim over BGZF must name the deflater it is a claim about
- [~] **jmath**: bit-exact Java `Math` / `StrictMath` / commons-math3 `FastMath`, plus `Gamma`,
      `BinomialDistribution`, SVD. `Percentile` and `Median` are ported and oracle-backed. The
      target is **not** "the corpus reaches 100%": its columns are `java.lang.Math`, whose
      remaining divergent functions can only be made exact by transcribing GPL2 source, so that
      is unreachable by construction. htsjdk-rs decision 0023 replaced it with "every function a
      ported call site reaches is exact, and every one that cannot be is named at the call site" 
- [ ] BGZF GKL/igzip surface cross-checked on real x86-64 hardware, not under emulation

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
