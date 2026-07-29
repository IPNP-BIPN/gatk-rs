# Roadmap to 100% byte-identical reproduction

Tracking document for the goal: a byte-for-byte reproduction of the **entire** GATK 4.6.2.0,
Picard 3.4.0 and htsjdk 4.2.0 tool set in Rust. Deliberately honest about scale (the plan sizes
the whole program at 40 to 100 person-years), and meant to be followed and ticked off.

Status legend: `[ ]` not started, `[~]` in progress or partial, `[x]` done **and** oracle-backed.

A box is only `[x]` when a golden produced by the pinned container on a real x86-64 runner is
committed, a Rust test compares against it, and CI re-derives it. Code that works but has no
golden stays `[~]`.

---

## Where we are (measured 2026-07-29)

| Repo | Scope | State |
|---|---|---|
| **htsjdk-rs** | the I/O and math foundation | substantially built; CRAM, GKL-exact deflate, full VCF and the jmath conformance corpus remain |
| **picard-rs** | 109 tools | ~50 tools have a first slice, ~43 with an oracle-backed conformance suite; many are partial (default paths only) |
| **gatk-rs** | 202 tools | 4 crates, **8 conformance suites, 13 goldens, all oracle-backed**; 1 tool byte-identical |

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
      fresh query per call, which the manifest states. Still missing: the codecs (VCF, BED,
      interval_list) and the Tribble index, which are htsjdk's and belong in Milestone H

### G1.4 Interval arguments

- [x] `-L` / `-XL`, padding, `--interval-set-rule`, `--interval-merging-rule`, subtraction
      (measured through `IntervalWalker`, 24 combinations)
- [x] interval **files**: `.list` and `.intervals` (lower-cased extension test, blank-line
      skipping, the empty-file refusal, and the order of the four tests in
      `parseIntervalArguments`), 13 arguments compared
- [ ] the Feature-file path (`.interval_list`, `.bed`, VCF): the seam exists and is named
      (`NoFeatureSources`), the codecs arrive with G1.3
- [x] `-L unmapped` end to end, measured through `ReadWalker` (5 runs): the tail comes after
      every interval, `-L unmapped` alone is a bounded traversal of nothing else, and an unmapped
      read carrying its mate's position is not in the tail at all

### G1.5 The pileup floor

- [x] `AlignmentStateMachine` (244 stops over 26 cigars)
- [x] `PileupElement` (217 elements, plus 231 `createPileupForReadAndOffset` calls including the
      offsets it refuses)
- [x] `ReadPileup`: the per-locus collection, its sorting, its sample split and the samtools
      overlap fix (3 pileups and 24 quality pairs)
- [ ] `ReadStateManager` and `PerSampleReadStateManager` (partition by sample)
- [ ] `SamplePartitioner`
- [ ] `LIBSDownsamplingInfo` and the downsampling itself (`--max-depth-per-sample`)
- [ ] `LocusIteratorByState` (merging the per-read machines into one pileup per locus)
- [ ] `IntervalAlignmentContextIterator` and `AlignmentContextIteratorBuilder`, including
      `emitEmptyLoci`, `includeDeletions` and `includeNs`

### G1.6 Walkers

- [x] `ReadWalker` (49 `apply` calls over 9 traversals)
- [x] `IntervalWalker` (25 `apply` calls over 24 argument combinations)
- [ ] `LocusWalker` (needs all of G1.5)
- [ ] `VariantWalker` (needs G1.3's feature sources)
- [ ] `AssemblyRegionWalker` (the base of G3)
- [ ] the multi-pass and multi-input walker variants

### G1.7 Annotations

- [ ] the 54-annotation library. Many go through jmath, so any annotation touching an unfinished
      jmath function waits rather than being approximated

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
      `BinomialDistribution`, `Percentile`, `Median`, SVD. The corpus must reach 100%
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

## Milestone V: program-level validation and reproducibility

- [x] one conformance manifest per repository, with the CI generated from it and a guard that
      fails if the generated YAML is edited by hand
- [x] the per-tool status dashboard, generated from the 311-tool inventory so it cannot drift
- [x] the golden audit: every committed golden is declared, every declared golden exists, and a
      keyed suite whose golden collapses is refused
- [x] the provenance guard (decision 0014), running in all three repositories
- [~] covering arrays generated and verified per tool (generated; not yet driving the oracle runs)
- [ ] the coverage-guided differential fuzzer
- [ ] determinism gates: same input twice, different `-Xmx` / `TMPDIR` / `LC_ALL`, two clean
      builds
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
   jmath corpus, which blocks every floating-point tool and therefore most annotations.
5. Then the **callers** (G3), then the **hard problems** (X).
6. **CRAM** and **GKL-exact deflate** proceed in parallel, being self-contained in htsjdk-rs.

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
