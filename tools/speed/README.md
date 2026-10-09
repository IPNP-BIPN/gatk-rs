# tools/speed

What the port costs against the reference, measured on the rows it is already proven
byte-identical on. Milestone S (#107); this directory is S.1 (#108).

```sh
python3 tools/speed/bench.py --port target/release/gatk-rs --out speed.json
python3 tools/speed/bench.py --port target/release/gatk-rs --tool CountReads --tool FlagStat
python3 tools/speed/bench.py --merge shard-*.json --out speed.json
python3 tools/speed/bench.py --compare first.json second.json
```

The port binary must be linux/amd64, because it runs inside the pinned container: CI builds one,
and a Mac builds one in `rust:1.97.1-bookworm` with `--target-dir target-linux`. On Apple
silicon both sides run under Rosetta, so local ratios are a smoke test; the numbers that are
published come from the `Speed` workflow on real x86-64.

| File | Role |
|---|---|
| `bench.py` | picks the rows, runs the timer per tool, summarises, merges shards, compares runs |
| `timer.py` | runs inside the container: alternates the two sides, reads CPU and peak RSS from `wait4` |
| `Steady.java` | the reference's steady state: one JVM, the same command line `n` times |

## What a number means

- **Rows.** A tool is measured only if its covering array is at share 1.000 in
  `tools/coverage/measured.json`, and only on the first row that, in the same run, the reference
  accepts and the port answers identically under `run_array.py`'s comparison. A path with no
  golden gets no row; a tool with no such row is listed with the reason.
- **Cold.** One process per run, both sides, after one warm-up repetition that fills the page
  cache. This is what GATK costs invoked once per file, JVM start and JIT warm-up included.
- **Steady.** The reference's tool run repeatedly in one JVM through `Main.instanceMain`, median
  of the second half of the iterations. The port has no JIT and nothing to warm, so its steady
  state is its cold run; the steady ratio is `port cold / reference steady`. A tool that throws on
  a second in-process run has no steady number, and says so.
- **Ratios.** `port / reference` for wall clock, CPU time (user + system) and peak RSS, medians
  over the repetitions. Below 1 the port is cheaper.
- **Noise.** Each repetition pairs one run of each side; the spread (max - min over median) of
  those paired ratios is the tool's cold noise, and the run's cold noise floor is its 90th
  percentile across tools. The steady ratio is not paired (a loop on one side, processes on the
  other), so its spread is bounded by the sum of the two sides' spreads, and its floor is stated
  separately; it is the wider of the two. `--compare` says how far each ratio moved between two
  runs, which is the check that the floors are honest.
- **Container.** Timing happens inside the container, so `docker run` itself is in neither number.
  The harness measures it on the host and prints it beside them.
