# Roadmap to 100% byte-identical reproduction

This is the tracking document for the goal stated in [the plan](docs) : a byte-for-byte
reproduction of the **entire** GATK 4.6.2.0, Picard 3.4.0 and htsjdk 4.2.0 tool set in Rust.
It is deliberately honest about scale (the plan sizes the whole program at 40 to 100 person-years)
and is meant to be followed and checked off, not admired.

Status legend: `[ ]` not started, `[~]` in progress / partial, `[x]` done and oracle-backed.

---

## Where we are (measured 2026-07-23)

| Repo | Scope | State |
|---|---|---|
| **htsjdk-rs** | the I/O + math foundation | substantially built; CRAM, GKL-exact deflate, full VCF and the jmath conformance corpus remain |
| **picard-rs** | 109 tools | ~40 tools have a first slice, ~37 with an oracle-backed conformance suite; many are partial (default paths only) |
| **gatk-rs** | 202 tools | inventory + oracle harness only; **0 tools ported** |

Totals from the generated inventory (`gatk-rs/tools/inventory`): **311 tools** (202 GATK-origin,
109 Picard-origin), **39 Spark**, ~13,130 arguments. Non-Spark: 163 GATK + 109 Picard.

### What "100% repro" means, per tool

A tool is **bit-identical** (the target) when, under the declared canonicalization, it is byte-equal to
the oracle for:

1. every applicable input under the pinned clones' `src/test/resources/**`;
2. a t-wise covering array over its arguments (t=2 everywhere, t=3 on the critical path);
3. the coverage-guided differential fuzzer reaching its branch-coverage threshold with no divergence;

with **zero quarantined fields**. Any quarantined field downgrades it to **bio-identical** with the
quarantine list attached. Most currently-ported tools reproduce only the default/common paths, so they
are not yet at this bar.

---

## Milestone H — finish the htsjdk-rs foundation

Everything downstream inherits these, so they are front-loaded. Gaps today:

- [~] BAM/SAM/tags/index read + write (core done; write-side BAI, CRAI, and full corner cases remain)
- [~] VCF / tribble (allele, variant, header, encoder exist; full read + write + index + all field types remain)
- [ ] **CRAM** (container model, all encodings, codec negotiation, reference-based compression) — a sub-project on its own
- [ ] **GKL-exact deflate** (ISA-L / igzip byte-exact) — the default non-JDK deflater path
- [~] **jmath**: bit-exact Java `Math`/`StrictMath`/commons-math3 `FastMath`, plus `Gamma`, `BinomialDistribution`, `Percentile`, `Median`, SVD — the conformance corpus must reach 100%
- [ ] BGZF GKL/igzip surface cross-checked on real x86-64 hardware (not just emulation)

---

## Milestone P — finish Picard (109 tools)

Ordered by archetype for amortization; floating-point-free tools first.

- [~] **MergeBamAlignment** (in flight): transfer, `PG` linkage and `NM`/`MD`/`UQ` and the whole-file
      coordinate-sorted producer are done and oracle-backed (PRs #85 to #87). Remaining slices:
  - [ ] merged-header construction (`@SQ` from the reference `.dict` with a `UR` absolute-path
        canonicalization rule, `@RG` from the unmapped BAM, `@PG` from the aligned BAM via the ported
        `SamFileHeaderMerger`) — brings the first full-file MergeBamAlignment oracle
  - [ ] paired mate-info + proper-pair + `ClippedPairFixer`
  - [ ] off-end-of-reference cigar clipping, `UNMAP_CONTAMINANT_READS`, adapter/overlap clipping
  - [ ] multi-hit selection (`MultiHitAlignedReadIterator` + the primary-alignment strategies)
- [ ] the remaining ~68 Picard tools, by archetype (metrics collectors, read transforms, interval and
      VCF utilities)
- [ ] **full parameter coverage** on every ported Picard tool (t-wise covering arrays + fuzzing), not
      just the default path — this is the bulk of the real work and is currently unmet almost everywhere

---

## Milestone G1 — the GATK engine (shared by everything above the tools)

Nothing in gatk-rs is ported yet; the engine must come first because the walkers build on it.

- [ ] the **read-filter** library (55 filters)
- [ ] the **annotation** library (54 annotations)
- [ ] the traversal **engine**: locus walker, variant walker, read walker, assembly-region walker
- [ ] interval parsing, reference/feature data sources, the Barclay argument layer (`--INPUT` syntax)

## Milestone G2 — the 163 non-Spark GATK tools, by archetype

- [ ] metrics collectors and reporting walkers (largest archetypes)
- [ ] read/variant transforms and locus/variant walkers
- [ ] reference, interval, coverage, CNV/SV and genotyping-array utilities
- [ ] full parameter coverage per tool

## Milestone G3 — the variant callers (individually multi-month)

- [ ] `HaplotypeCaller`, `Mutect2` and relatives: assembly graph, genotyping
- [ ] **PairHMM** targeting the pure-Java `LoglessPairHMM` semantics (the implementation must be pinned
      in the oracle contract, since `FASTEST_AVAILABLE` resolves per host)

---

## Milestone X — the four hard problems

- [ ] **Spark** (39 tools): establish bit-identity against the `--spark-master local[1]` single-partition
      oracle first; prove any parallel path byte-equal to it in CI
- [ ] **ML inference** (CNNScoreVariants, NVScoreVariants, relatives): reproduce the exact kernels and
      accumulation order, validated layer-by-layer against dumped intermediate tensors; pin model weights
      by hash; quarantine and report as bio-identical if a model proves intractable
- [ ] **CRAM** (if not completed under Milestone H)
- [ ] **GKL-exact deflate** (if not completed under Milestone H)

---

## Milestone V — program-level validation and reproducibility

- [ ] per-tool validation reports (inputs, records, byte-equal file count, argument coverage, t level,
      fuzzer branch coverage, quarantined fields)
- [ ] the **dashboard** generated from the inventory: per tool, ported / byte-equal / coverage / quarantine
- [ ] determinism gates: same input twice, different `-Xmx`/`TMPDIR`/`LC_ALL`, two clean builds
- [ ] per-run provenance (input hashes, reference tags, toolchain, container digest, GKL provider state)

---

## The critical path, in order

1. Finish **MergeBamAlignment** (Milestone P, in flight).
2. Close the **htsjdk-rs** foundation gaps that block downstream tools: full VCF, write-side BAI, then
   the jmath conformance corpus (blocks every floating-point tool).
3. Stand up the **GATK engine** (Milestone G1): read filters + annotations + the walker engine. This is
   the single biggest unlock, because it is the prerequisite for 163 GATK tools.
4. Fan out **G2 by archetype**, measuring the marginal cost of the second and third member of each
   archetype (the plan's calibration gate) and re-sizing from real numbers.
5. Then the **callers** (G3), then the **hard problems** (X).
6. **CRAM** and **GKL-exact deflate** can proceed in parallel with the tool fan-out, since they are
   self-contained within htsjdk-rs.

The honest bottom line: the tool-by-tool fan-out is tractable and amortizes well; the schedule is
dominated by (a) the GATK engine, (b) the variant callers, and (c) the four hard problems. Full
parameter coverage on every tool, not the count of tools with a first slice, is what the 100% claim
ultimately rests on.

---

## How this is tracked

- This file is the human-readable milestone list. Each `[ ]` is a checkable unit.
- GitHub milestones mirror H / P / G1 / G2 / G3 / X / V, with one tracking issue per milestone.
- The per-tool status dashboard (Milestone V) is generated from `tools/inventory` so it cannot drift
  from the 311-tool ground truth.
- Progress lands as one PR per slice (one ported symbol or feature), CI-gated on the digest-pinned
  oracle, merged only when green.
