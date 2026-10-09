# The speed baseline

Milestone S, S.3 (#110). The first numbers the programme has on what the port costs, measured
before any optimisation, and the predictions the issue wrote down in advance, marked against them.

**Source.** `tools/speed/baseline.json` is the `speed` artefact of the `Speed` workflow, run
37917687555 on commit `26432cc8`, eight shards on GitHub's x86-64 runners, not a laptop. The
per-tool rows are in [STATUS.md](../STATUS.md), in the cost column beside each tool's byte-identity
claim. The harness and what each number means are in `tools/speed/README.md`.

**Scope.** Every gatk-rs tool whose covering array is at share 1.000: 137 tools, 135 measured. Two
have no row the reference accepts (`AlleleFrequencyQC`, which ends every accepted row on an R
error, and `FuncotatorDataSourceDownloader`, whose rows the reference refuses on this corpus), so they have no number,
rather than a number on a row the claim does not cover. The unit is the tool, because that is
what a user runs; the read filters, the pileup floor, the annotation library and the walkers are
measured through the tools that run them, on the rows those tools are proven on.

## The distribution

Ratios are `port / reference`; below 1 the port is cheaper.

| | min | 25% | median | 75% | 90% | max |
|---|---:|---:|---:|---:|---:|---:|
| cold wall clock | 0.0018 | 0.0021 | 0.0025 | 0.0031 | 0.0047 | 4.83 |
| steady wall clock | 0.028 | 0.049 | 0.055 | 0.070 | 0.115 | 18.4 |
| cold CPU time | 0.0006 | 0.0008 | 0.0009 | 0.0012 | 0.0019 | 2.01 |
| peak RSS | 0.032 | 0.052 | 0.052 | 0.053 | 0.055 | 0.21 |

In absolute terms, on the same rows: the reference takes 1.1 to 2.7 s cold (median 1.75 s, of
which the first in-JVM iteration is almost all), 38 to 526 ms once warm (median 75 ms), and
244 to 442 MiB; the port takes 2.8 ms to 9.7 s (median 4.0 ms) and 14 to 58 MiB.

## Where the port is not cheaper

Four tools, the only ones whose steady ratio reaches 1. They are S.5's (#112) candidates, and
nothing else is close: the next is `LearnReadOrientationModel` at 0.64.

| Tool | port | reference cold | reference steady | cold | steady |
|---|---:|---:|---:|---:|---:|
| `ModelSegments` | 9.67 s | 1.99 s | 0.53 s | **4.83** | **18.4** |
| `VariantRecalibrator` | 1.00 s | 2.74 s | 0.22 s | 0.37 | **4.63** |
| `HaplotypeBasedVariantRecaller` | 158 ms | 2.29 s | 121 ms | 0.07 | 1.30 |
| `FlowPairHMMAlignReadsToHaplotypes` | 120 ms | 1.34 s | 114 ms | 0.09 | 1.05 |

`ModelSegments` is slower than the reference even with the JVM's start-up charged to the
reference: nearly ten seconds of single-threaded CPU on a row the reference finishes in two.
None of the four uses more than one core (port CPU time equals its wall clock to the millisecond),
while the reference's cold CPU time is two to three times its wall clock: the JVM compiles and
collects on other cores.

## The predictions, marked

| Prediction (#110, written before the numbers) | Verdict |
|---|---|
| Short runs favour the port heavily, because JVM start-up is paid per invocation. | **Right**, and by more than anyone would quote: a median cold ratio of 0.0025, i.e. the reference spends about 1.7 s starting and the port about 4 ms working. It is also the least interesting number here, as predicted. |
| Steady-state throughput is genuinely uncertain. | **Mostly wrong for these rows**: against a warmed-up JVM the port is still cheaper on 131 of 135 tools, by a median factor of 18 (0.055) and never less than 1.5 (`LearnReadOrientationModel`, 0.64). It is right exactly where the work is numerical and long (the four above). |
| The port may be slower in places, and memory is where it is most exposed. | **Slower in places: right** (four tools). **Memory: wrong.** Peak RSS is 3 to 21% of the reference's on every tool; the largest is `BaseRecalibrator` at 0.21 (58 MiB against 272 MiB). The clones the issue names do not show at this scale. |
| I/O may dominate everything else. | **Not decided by this corpus.** The fixtures are the conformance corpus, a few kilobytes per input, and the port's median run of 4 ms is mostly process start. A row large enough for I/O to show is outside what the arrays prove, so the question needs a bigger input that is still covered by a golden. |

## How much to trust a single number

Within one run, the harness states a noise floor (90th percentile of per-tool spread): 24% cold
and 63% steady, the steady bound being conservative because a loop and separate processes cannot
be paired. Across runs it is worse than the within-run floor suggests: the first, failed run of
the same commit (its shards printed their ratios before the teardown bug) and this one landed on
different runner models (five CPU models across the eight shards), and per-tool ratios moved by a
median of 10% and at most 49%, cold and steady alike.

So a difference under about half is not a finding on one tool from one run; the four rows above
are, and they held (`ModelSegments` 20.0 then 18.4 steady, `VariantRecalibrator` 4.5 then 4.6).
Most of the noise is the port's own scale: a 4 ms run is at the mercy of the scheduler. A
regression check should compare ratios across the whole table with `bench.py --compare`, not
read one cell.
