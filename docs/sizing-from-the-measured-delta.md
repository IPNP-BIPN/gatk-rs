# Sizing the program from the measured delta, not from the archetype story

The plan's cost estimate — 40 to 100 person-years — rested on an assumption about archetype
amortisation that has now been measured and found to work by a different mechanism. This
records what the measurement implies for sizing, so the next revision of the plan starts from a
number rather than from a story.

## What was measured

Two within-stratum pairs, one at each end of the size distribution, each ported over identical
inputs so that shared machinery is visible.

**The small end** (`picard-rs/docs/decisions/0002`):

| | Java ported | Rust written (non-test) | ratio |
|---|---:|---:|---:|
| shared convention | — | 67 | — |
| `MeanQualityByCycle` | 385 | 151 | **0.39** |
| `CollectBaseDistributionByCycle` | 291 | 169 | **0.58** |

**The large end** (`picard-rs/docs/decisions/0003`), the measurement this document previously
asked for:

| | Java ported | Rust written (non-test) | ratio |
|---|---:|---:|---:|
| `CollectInsertSizeMetrics` | 526 | 399 | **0.76** |
| `CollectAlignmentSummaryMetrics` | 1217 | 1030 | **0.85** |

Counting the new shared machinery the second member forced into `htsjdk-rs` — alignment blocks,
IUPAC base comparison, the FASTA reader, a histogram key-union fix, about 190 further lines —
the second member cost **1.00 Rust lines per Java line**.

In both pairs the second member cost **more** per Java line than the first. The archetype delta
is negative at both ends of the one stratum where the hypothesis was most likely to hold.

## The corpus to be ported

Main source only, excluding tests, at the pinned tags:

| repository | files | lines of Java |
|---|---:|---:|
| htsjdk 4.2.0 | 792 | 132,150 |
| Picard 3.4.0 | 499 | 91,221 |
| GATK 4.6.2.0 | 1,589 | 292,297 |
| **total** | **2,880** | **515,668** |

## What that implies, revised

The earlier revision of this document used the small-end pair's **0.39 to 0.58** and arrived at
200,000 to 300,000 lines of Rust for the ported main source. That was too low, and the reason is
now visible: the small-end pair's tools are small *because* they do little, and a tool that does
little has a low ratio. The large-end pair, which is where most of the corpus's lines actually
live, runs at **0.76 to 1.00**.

Weighting by where the lines are rather than by where the first measurement happened to land,
the ported main source comes to roughly **390,000 to 515,000 lines of Rust**, not 200,000 to
300,000. **The line estimate roughly doubles.**

The tests are not a rounding error on top of that. Across the five collectors ported so far,
test code has run between **50% and 95%** of the port's own size, and that is *before* the
t-wise covering arrays and coverage-guided fuzzing the plan commits to. Taking the low end,
total output is on the order of **600,000 to 900,000 lines**.

Two things that number does **not** include, and both are large:

- **The four hard problems.** CRAM, Spark's 39 tools, bit-identical ML inference, and
  GKL-exact deflate are each their own sub-programme and none of them is line-count-shaped.
- **The findings.** This is the important omission. The session that produced the delta also
  produced thirteen decision records, three of which correct earlier decisions of mine. Several
  cost hours and changed almost no lines: the `EnumMap` iteration order is one line; the FPU's
  NaN sign is zero lines. Line counts measure typing, and typing is not where this work lives.

## What replaces the archetype model

The plan's model was "port the shape once, then pay a small delta per member". Two measurements
now say something narrower and more useful:

> **The per-member cost tracks the member's own Java footprint, at 0.76 to 1.00 Rust lines per
> Java line. The stratum predicts almost nothing beyond that.**

The stratum signals turn out to describe the plumbing — how records arrive (`single_pass`),
whether the output has a histogram section (`histogram`), which base class the bean extends
(`multi_level`) — and all of that was paid for before any of the five tools was written. What
the two large-end tools actually share is 20 lines.

That makes `tools/stratify/stratify.py`'s **footprint** the sizing input, and the stratum a
scheduling convenience rather than a cost model. It is a worse answer for the program than the
archetype story, and it is the one the measurements give.

## The honest position

The 40-to-100 person-year range is **not** refuted, and it is not confirmed. What the two
measurements do is remove a specific way it could have been optimistic, and then move the line
estimate up by roughly a factor of two inside it.

Both pairs are in the metrics stratum, which is the most homogeneous of the twenty-four. The
variant callers, the Spark tools and the ML tools have had no measurement at all, and none of
them is line-count-shaped. So the direction is now checked at both ends of the easiest case and
unchecked everywhere else.

Two things the number still does not include:

- **The four hard problems.** CRAM, Spark's 39 tools, bit-identical ML inference, and GKL-exact
  deflate are each their own sub-programme.
- **The findings**, which remain the largest omission. The work since the last revision produced
  decisions 0017, 0018 and 0019 in `htsjdk-rs` and 0003 in `picard-rs`. The most expensive of
  them changed almost nothing: Picard's `BAD_CYCLES` is binned by the offset within an alignment
  block rather than by the read cycle, which took a probe in the oracle to establish and one
  comment plus one variable name to reproduce. Decision 0019 cost a probe and produced *no*
  code at all, because the divergence it looked for was not there.

  A sizing model built on line counts systematically under-prices exactly the work that makes
  the port bit-identical rather than merely correct.
