# The move to GATK 4.7.0.0, measured before anything moves

V.10 (#810) moves the reference to GATK 4.7.0.0, htsjdk 5.0.0 and Picard 3.5.0. Moving it does not
adjust the oracle-backed suites, it re-opens them, so the first step is to find out how many it
re-opens and why. `tools/conformance/survey.py` regenerates every oracle-backed suite under a
4.7.0.0 oracle and compares it with the 4.6.2.0 golden it already has;
`.github/workflows/survey.yml` runs it in 36 shards on any `survey/**` branch. The fixtures are
still built by the 4.6.2.0 image, so the inputs stay fixed and only the reference moves.

Measured on runs 42827131133 and 42833576338, which agree suite for suite.

## 362 suites

| | suites | what it means for the port |
|---|---:|---|
| identical under 4.7.0.0 | 316 | nothing |
| the stamp only | 27 | the version, and one new global argument in every command line |
| the dumps' own SAM text | 7 | nothing: fixed in the dumps here |
| the ground truth grows | 6 | one new tool, and arguments the port has to declare |
| the reference behaves differently | 5 | port work, listed below |
| not surveyed | 1 | `r-plots` runs in the Broad's own 4.6.2.0 image, which the override does not reach |

**The stamp (27).** 4.7.0.0 adds `--output-cram-version`, defaulting to `3.1`, to every tool that
writes reads or variants, so `@PG CL:` and `##GATKCommandLine` grow ` --output-cram-version 3.1`,
and the BAMs the dumps carry change with their header. `survey.py` now takes the stamp out of both
sides (`NEW_ARGUMENTS`, the BAM header decompressed and its length recomputed, an index set aside
because it follows the header) and the 27 match. For the port the stamp is not free: the argument
has to be declared, and printed, by 155 tools.

**The dumps' SAM text (7).** `SAMRecord.getSAMString()` ends its line with `\n` up to htsjdk 4.2.0
and without it from 5.0.0. `PrintReadsDump.samLine` already hid that, but seven dumps still called
`getSAMString()` directly and concatenated their records into one line; they now ask `samLine`,
which gives the same bytes under 4.6.2.0. A property of the instrument, not of the tools.

**The ground truth grows (6).** `main-catalogue`, `main-usage` and `main-entry` count
`ConvertCountsToDepthFile` (#811), 331 classes becoming 332. `usage-text`,
`tool-argument-declarations` and `tool-argument-enums` carry the new arguments:

| argument | tools |
|---|---:|
| `--output-cram-version` | 155 |
| `--dragen-pd-pair-hmm-implementation` and the four `--native-dragen-pd-pair-*` (the PDHMM through GKL 0.9.1) | 9 |
| `--activeregion-alt-multiplier` | 5 (HaplotypeCaller and its kin) |
| `--somatic-quality-threshold` | GenotypeGVCFs |
| nine `--stratify-*`, `--track-*`, `--clustering-config` and `--max-records-in-ram`; `--do-not-sort` removed | SVConcordance |

**The reference behaves differently (5).**

- `sv-concordance`: every record gains `STRAT`, `TRUTH_DISTANCE_START`, `TRUTH_DISTANCE_END`,
  `TRUTH_RECIPROCAL_OVERLAP` and `TRUTH_SIZE_SIMILARITY`, and `--do-not-sort` is now refused as an
  unknown option (broadinstitute/gatk#9125).
- `sv-stratify`: the case the reference refused now writes an output (intersecting groups,
  broadinstitute/gatk#9185).
- `grouped-sv-cluster`: "matched multiple groups" names them as `DEL_a, DEL_b` where it said
  `DEL_b, DEL_a`; the order of a set reaches the message.
- `spark-walkers`: the `filter-no-paired` case throws `NullPointerException` where it threw
  `IllegalArgumentException`.
- `cram-issue-8768-detector`: the CRAM it reads is written as 3.1, the new default.

## What this says about the move

The engine, the readers and writers, the annotations and the callers' outputs survive it: no
golden of theirs moved except by the stamp. The port's work is one global argument across 155
tools, one tool, the PDHMM's arguments (not its kernel: no suite reaches a DRAGEN PDHMM path), the
SVConcordance and SVStratify changes, one message order and one exception class, and CRAM 3.1
where a tool writes CRAM. The covering arrays have not been surveyed yet; they measure the port's
side as well, so they come after the port has the new arguments.
