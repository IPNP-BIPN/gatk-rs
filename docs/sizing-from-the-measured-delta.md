# Sizing the program from the measured delta, not from the archetype story

The plan's cost estimate — 40 to 100 person-years — rested on an assumption about archetype
amortisation that has now been measured and found to work by a different mechanism. This
records what the measurement implies for sizing, so the next revision of the plan starts from a
number rather than from a story.

## What was measured

`picard-rs/docs/decisions/0002` reports the first within-stratum delta. Two stratum-mates,
ported over identical inputs so shared machinery is visible:

| | Java ported | Rust written (non-test) | ratio |
|---|---:|---:|---:|
| shared convention | — | 67 | — |
| `MeanQualityByCycle` | 385 | 151 | **0.39** |
| `CollectBaseDistributionByCycle` | 291 | 169 | **0.58** |

The second member cost *more* per line than the first. The amortisation is real but it is
**~18% of the second member's cost**, not the near-zero the plan's model implied.

## The corpus to be ported

Main source only, excluding tests, at the pinned tags:

| repository | files | lines of Java |
|---|---:|---:|
| htsjdk 4.2.0 | 792 | 132,150 |
| Picard 3.4.0 | 499 | 91,221 |
| GATK 4.6.2.0 | 1,589 | 292,297 |
| **total** | **2,880** | **515,668** |

## What that implies

At the measured **0.39 to 0.58 Rust lines per Java line**, ported main source alone comes to
roughly **200,000 to 300,000 lines of Rust**.

The tests are not a rounding error on top of that. Across the four collectors ported so far,
test code has run between **50% and 95%** of the port's own size, and that is *before* the
t-wise covering arrays and coverage-guided fuzzing the plan commits to. Taking the low end,
total output is on the order of **300,000 to 500,000 lines**.

Two things that number does **not** include, and both are large:

- **The four hard problems.** CRAM, Spark's 39 tools, bit-identical ML inference, and
  GKL-exact deflate are each their own sub-programme and none of them is line-count-shaped.
- **The findings.** This is the important omission. The session that produced the delta also
  produced thirteen decision records, three of which correct earlier decisions of mine. Several
  cost hours and changed almost no lines: the `EnumMap` iteration order is one line; the FPU's
  NaN sign is zero lines. Line counts measure typing, and typing is not where this work lives.

## The honest position

The 40-to-100 person-year range is **not** refuted by this measurement, and it is not confirmed
either. What the measurement does is remove one specific way the estimate could have been
optimistic: it is now known that members 2..n of an archetype are not nearly free, so any
sizing that leaned on that is wrong.

It replaces a guess with one data point, at the small end of one stratum, in one of the two
Picard-shaped archetypes. Decision 0002 lists what would strengthen it: a second delta at the
large end of the same stratum, where `CollectRnaSeqMetrics` (879 lines) and
`CollectAlignmentSummaryMetrics` (978) sit, and where the input plumbing the stratification
flagged as "environment cost" actually comes due.
