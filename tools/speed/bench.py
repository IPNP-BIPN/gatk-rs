#!/usr/bin/env python3
"""Measure what the port costs against the reference, on rows it is already proven identical on.

    python3 tools/speed/bench.py --port target/release/gatk-rs --out speed.json
    python3 tools/speed/bench.py --port target/release/gatk-rs --tool CountReads --tool FlagStat
    python3 tools/speed/bench.py --merge shard-*.json --out speed.json
    python3 tools/speed/bench.py --compare first.json second.json

Milestone S, S.1 (#108). The conformance harness says what the port PRODUCES; this says what it
COSTS, and it holds itself to the same rules:

* **No row without a golden.** A tool is measured only if its covering array is at share 1.000 in
  `tools/coverage/measured.json`, and only on a row that, in THIS run, the reference accepts and
  the port answers byte for byte (through `run_array.py`'s own comparison). A tool with no such row
  gets no number, and the output says why. A benchmark cannot drift onto an easier case than the
  correctness claim covers, because it runs the correctness claim's own inputs: the corpus is
  `MakeFixtures.java`'s.
* **Both sides in the same container.** The port runs inside the pinned image at the same paths,
  as it does for the arrays, so the comparison is linux/amd64 against linux/amd64 over the same
  mounted files. The timing happens inside the container (`timer.py`), so the container's own start
  is outside both numbers; it is measured once on the host and stated beside them, not subtracted.
* **Three numbers, not one.** Wall clock, CPU time and peak RSS per run, because wall clock alone
  hides a port that wins by using four cores where the reference used one.
* **Cold separately from steady state.** A cold run is one process per invocation, which is how
  GATK is most often run and where the JVM pays class loading and JIT warm-up. The steady state is
  the reference's tool run repeatedly in one JVM (`Steady.java`), median of the second half. The
  port has no JIT and nothing to warm: its steady state is its cold run, and the output says so
  rather than inventing a second number.
* **A ratio, not a time.** `port / oracle`, per row, on one host in one session. Seconds on one
  machine say nothing about another; the ratio travels.
* **Its own noise floor.** Each repetition pairs one run of each side, so each tool has as many
  ratios as repetitions; their spread is the tool's noise, and the run's noise floor is the
  90th percentile of those spreads. `--compare` reads two runs and says how far each ratio moved,
  which is the check that a regression can be told from noise.
"""

import argparse
import json
import math
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SPEED = REPO / "tools" / "speed"
sys.path.insert(0, str(REPO / "tools" / "coverage"))
import run_array  # noqa: E402  (the array runner is the one place rows become command lines)

MEASURED = REPO / "tools" / "coverage" / "measured.json"


def proven_tools():
    """Every tool whose array the port matches on every row, in name order."""
    tools = json.loads(MEASURED.read_text())["tools"]
    return sorted(
        name
        for name, entry in tools.items()
        if entry.get("share") == 1.0 and (run_array.ARRAYS / f"{name}.t{entry['t']}.json").exists()
    )


def median(values):
    return statistics.median(values) if values else None


def spread(values):
    """(max - min) / median: how far apart the repetitions landed, as a fraction."""
    if len(values) < 2:
        return None
    middle = statistics.median(values)
    return (max(values) - min(values)) / middle if middle else None


def ratio(a, b):
    return a / b if a is not None and b else None


def compile_steady(workdir):
    """`Steady.java` against the reference's classpath, once per run, inside the container."""
    classes = workdir / "steady"
    classes.mkdir(exist_ok=True)
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", run_array.PLATFORM,
            "-v", f"{SPEED}:/harness:ro",
            "-v", f"{classes}:/out",
            "-w", "/work", run_array.IMAGE,
            'javac -cp "$ORACLE_CP" -d /out /harness/Steady.java',
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(result.stderr[-2000:])
        raise SystemExit("could not compile Steady.java")
    return classes


def container_overhead(samples):
    """What `docker run --rm <image> true` costs on this host, which neither side's number holds."""
    times = []
    for _ in range(samples):
        start = time.perf_counter()
        subprocess.run(
            ["docker", "run", "--rm", "--platform", run_array.PLATFORM, run_array.IMAGE, "true"],
            capture_output=True,
        )
        times.append(time.perf_counter() - start)
    return {"samples": samples, "median_wall": median(times), "spread": spread(times)}


def proven_row(tool, array, workdir, port, max_probe):
    """The first row the reference accepts and the port answers identically, or why there is none.

    `run_array.py`'s comparison, unchanged, decides "identically": the same canonicalisation the
    share was measured with, so a row timed here is a row the claim covers.
    """
    held = array.get("excluded", [])
    positional = run_array.positional_values(tool)
    tagged = run_array.tagged_arguments(tool)
    lists = run_array.list_arguments(tool)
    keep_output = run_array.output_on_failure(tool)
    probed = 0
    for row in array["array"]:
        args = run_array.row_arguments(row, held)
        code, text, error = run_array.run_oracle(tool, args, workdir, positional, tagged, lists)
        if code != 0:
            continue
        probed += 1
        reference = run_array.outcome(code, text, error, keep_output)
        port_code, port_text, port_error = run_array.run_port(
            port, tool, args, workdir, positional, tagged, lists
        )
        if run_array.outcome(port_code, port_text, port_error, keep_output) == reference:
            cli = " ".join(run_array.as_cli(args, positional, tagged, lists))
            return {"row": row["row"], "arguments": args, "cli": cli}, None
        if probed >= max_probe:
            return None, f"none of the first {probed} accepted rows matched in this run"
    if probed == 0:
        return None, "the reference accepts no row of the array"
    return None, f"none of the {probed} accepted rows matched in this run"


def time_tool(tool, row, workdir, port, classes, reps, iterations, timeout):
    """Run `timer.py` once in the container over the chosen row, and summarise what it measured."""
    spec_dir = workdir / "spec"
    spec_dir.mkdir(exist_ok=True)
    binary = Path(port).resolve()
    java = f"java {run_array.JAVA_OPENS}"
    spec = {
        "reps": reps,
        "timeout": timeout,
        "iterations": iterations,
        "oracle": f'{java} -cp "$ORACLE_CP" org.broadinstitute.hellbender.Main {tool} {row["cli"]}',
        "port": f"/work/port-binary/{binary.name} {tool} {row['cli']}",
        "steady": (
            f'{java} -cp "/work/steady:$ORACLE_CP" Steady {iterations} @REPORT@ {tool} {row["cli"]}'
            if iterations
            else None
        ),
    }
    (spec_dir / "spec.json").write_text(json.dumps(spec))
    out_dir = workdir / "out"
    out_dir.mkdir(exist_ok=True)
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", run_array.PLATFORM,
            "-v", f"{run_array.fixtures_of(workdir)}:/work/fixtures:ro",
            "-v", f"{out_dir}:/work/out",
            "-v", f"{binary.parent}:/work/port-binary:ro",
            "-v", f"{classes}:/work/steady:ro",
            "-v", f"{spec_dir}:/work/spec:ro",
            "-v", f"{SPEED}:/harness:ro",
            "-w", "/work", run_array.IMAGE,
            "python3 /harness/timer.py /work/spec/spec.json",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None, f"timer failed: {(result.stderr or result.stdout).strip()[-300:]}"
    raw = json.loads(result.stdout.strip().splitlines()[-1])
    return summarise(raw), None


def summarise(raw):
    """Medians over repetitions one onwards, the paired ratios, and the reference's steady state."""
    cold = raw["cold"]
    failed = [
        side for side in ("oracle", "port") if any(run["exit"] != 0 for run in cold[side])
    ]
    oracle, port = cold["oracle"][1:], cold["port"][1:]

    def side(runs):
        return {
            "wall": median([r["wall"] for r in runs]),
            "cpu": median([r["cpu"] for r in runs]),
            "rss_kib": median([r["rss_kib"] for r in runs]),
            "wall_spread": spread([r["wall"] for r in runs]),
        }

    paired = [p["wall"] / o["wall"] for o, p in zip(oracle, port) if o["wall"]]
    summary = {
        "oracle_cold": side(oracle),
        "port": side(port),
        "ratio_cold_wall": median(paired),
        "ratio_cold_cpu": ratio(side(port)["cpu"], side(oracle)["cpu"]),
        "ratio_rss": ratio(side(port)["rss_kib"], side(oracle)["rss_kib"]),
        "ratio_spread": spread(paired),
        "failed_runs": failed,
        "host": raw["host"],
        "raw": raw,
    }
    steady = raw.get("steady")
    if steady:
        loop = steady["iterations"]
        complete = len(loop) == steady["requested"] and all(i["status"] == "ok" for i in loop)
        if complete:
            warm = loop[len(loop) // 2:]
            summary["oracle_steady"] = {
                "wall": median([i["wall"] for i in warm]),
                "cpu": median([i["cpu"] for i in warm]),
                "first_wall": loop[0]["wall"],
                "iterations": len(loop),
                "wall_spread": spread([i["wall"] for i in warm]),
            }
            summary["ratio_steady_wall"] = ratio(
                summary["port"]["wall"], summary["oracle_steady"]["wall"]
            )
            # The two sides are not paired here, one being a loop and the other separate
            # processes, so the ratio's spread is bounded by the sum of the two spreads.
            if None not in (summary["oracle_steady"]["wall_spread"], summary["port"]["wall_spread"]):
                summary["ratio_steady_spread"] = (
                    summary["oracle_steady"]["wall_spread"] + summary["port"]["wall_spread"]
                )
        else:
            stopped = next((i["status"] for i in loop if i["status"] != "ok"), "no report")
            summary["oracle_steady"] = None
            summary["steady_unavailable"] = (
                f"the in-process loop stopped after {len(loop)} of {steady['requested']}: {stopped}"
            )
    return summary


def percentile(values, fraction):
    values = sorted(v for v in values if v is not None)
    if not values:
        return None
    return values[min(len(values) - 1, math.ceil(fraction * len(values)) - 1)]


def finish(document):
    """The run-level numbers that only exist once every tool is in."""
    measured = [e for e in document["tools"].values() if e.get("ratio_cold_wall") is not None]
    document["noise_floor"] = percentile([e["ratio_spread"] for e in measured], 0.9)
    document["noise_floor_steady"] = percentile(
        [e.get("ratio_steady_spread") for e in measured], 0.9
    )
    document["measured"] = len(measured)
    return document


def table(document):
    print(f"\n{'tool':44} {'cold':>7} {'steady':>7} {'cpu':>7} {'rss':>7} {'noise':>6}")
    for tool, entry in sorted(document["tools"].items()):
        if entry.get("ratio_cold_wall") is None:
            print(f"{tool:44} not measured: {entry.get('reason')}")
            continue

        def fmt(value):
            return f"{value:7.3f}" if value is not None else "      -"

        print(
            f"{tool:44} {fmt(entry['ratio_cold_wall'])} {fmt(entry.get('ratio_steady_wall'))} "
            f"{fmt(entry['ratio_cold_cpu'])} {fmt(entry['ratio_rss'])} "
            f"{fmt(entry['ratio_spread'])[1:]}"
        )
    print(
        f"\n{document['measured']} tools measured; ratios are port / reference; "
        f"noise floor (90th percentile of per-tool spread): cold {fmt_share(document['noise_floor'])}, "
        f"steady {fmt_share(document['noise_floor_steady'])}"
    )
    overhead = document.get("container_overhead") or {}
    if overhead.get("median_wall") is not None:
        print(f"container start, in neither number: {overhead['median_wall']:.3f}s median")


def fmt_share(value):
    return f"{value:.1%}" if value is not None else "-"


def merge(paths, out):
    document = {"tools": {}, "shards": []}
    for path in paths:
        part = json.loads(Path(path).read_text())
        document["tools"].update(part["tools"])
        document["shards"].append(
            {key: part.get(key) for key in ("host", "container_overhead", "reps", "iterations")}
        )
        for key in ("reps", "iterations", "commit"):
            document.setdefault(key, part.get(key))
    finish(document)
    Path(out).write_text(json.dumps(document, indent=1, sort_keys=True) + "\n")
    table(document)
    return 0


def compare(first, second):
    """How far each tool's ratio moved between two runs, against the noise each run claims."""
    a = json.loads(Path(first).read_text())["tools"]
    b = json.loads(Path(second).read_text())["tools"]
    for key in ("ratio_cold_wall", "ratio_steady_wall"):
        drifts = []
        print(f"\n{key:44} {'first':>7} {'second':>7} {'moved':>7}")
        for tool in sorted(set(a) & set(b)):
            ra, rb = a[tool].get(key), b[tool].get(key)
            if ra is None or rb is None:
                continue
            moved = abs(rb - ra) / ra
            drifts.append(moved)
            print(f"{tool:44} {ra:7.3f} {rb:7.3f} {moved:7.1%}")
        if drifts:
            print(
                f"{len(drifts)} tools: median move {statistics.median(drifts):.1%}, "
                f"90th percentile {percentile(drifts, 0.9):.1%}, largest {max(drifts):.1%}"
            )
    return 0


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--port", help="the linux/amd64 port binary (gatk-rs)")
    parser.add_argument("--tool", action="append", help="measure this tool (default: every proven one)")
    parser.add_argument("--shard", help="i/n: measure every n-th proven tool from the i-th")
    parser.add_argument("--reps", type=int, default=5, help="cold repetitions per side, after a warm-up")
    parser.add_argument("--iterations", type=int, default=20, help="in-JVM iterations; 0 skips the steady state")
    parser.add_argument("--max-probe", type=int, default=6, help="accepted rows to try before giving up")
    parser.add_argument("--timeout", type=int, default=600, help="seconds per run")
    parser.add_argument("--fixtures-dir", help="a corpus MakeFixtures.java already built")
    parser.add_argument("--out", help="write the measurement here")
    parser.add_argument("--merge", nargs="+", help="combine shard outputs into one --out")
    parser.add_argument("--compare", nargs=2, help="two measurements: how far each ratio moved")
    options = parser.parse_args(argv)

    if options.merge:
        return merge(options.merge, options.out or "speed.json")
    if options.compare:
        return compare(*options.compare)
    if not options.port:
        parser.error("--port is required to measure")

    proven = proven_tools()
    if options.tool:
        refused = [t for t in options.tool if t not in proven]
        if refused:
            raise SystemExit(
                f"no byte-identity claim to measure against for {', '.join(refused)}: "
                "a tool is measured only at share 1.000 in tools/coverage/measured.json"
            )
        tools = options.tool
    else:
        tools = proven
    if options.shard:
        index, count = (int(x) for x in options.shard.split("/"))
        tools = tools[index::count]

    measured = json.loads(MEASURED.read_text())["tools"]
    document = {
        "reps": options.reps,
        "iterations": options.iterations,
        "commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=REPO
        ).stdout.strip(),
        "tools": {},
    }
    with tempfile.TemporaryDirectory(prefix="gatk-speed-", ignore_cleanup_errors=True) as tmp:
        workdir = Path(tmp)
        try:
            if options.fixtures_dir:
                run_array.FIXTURES_DIR = options.fixtures_dir
            else:
                run_array.build_fixtures(workdir)
            (workdir / "tmp").mkdir(exist_ok=True)
            classes = compile_steady(workdir) if options.iterations else workdir
            document["container_overhead"] = container_overhead(5)
            for tool in tools:
                array = json.loads(
                    (run_array.ARRAYS / f"{tool}.t{measured[tool]['t']}.json").read_text()
                )
                row, reason = proven_row(tool, array, workdir, options.port, options.max_probe)
                if row is None:
                    document["tools"][tool] = {"reason": reason}
                    print(f"{tool}: not measured, {reason}", flush=True)
                    continue
                summary, reason = time_tool(
                    tool, row, workdir, options.port, classes,
                    options.reps, options.iterations, options.timeout,
                )
                if summary is None:
                    document["tools"][tool] = {"row": row, "reason": reason}
                    print(f"{tool}: not measured, {reason}", flush=True)
                    continue
                if summary["failed_runs"]:
                    reason = f"a timed run failed on the {' and '.join(summary['failed_runs'])} side"
                    document["tools"][tool] = {"row": row, "reason": reason, "raw": summary["raw"]}
                    print(f"{tool}: not measured, {reason}", flush=True)
                    continue
                document["tools"][tool] = {"row": row, **summary}
                document.setdefault("host", summary["host"])
                print(
                    f"{tool}: row {row['row']} cold {summary['ratio_cold_wall']:.3f} "
                    f"steady {summary.get('ratio_steady_wall') or float('nan'):.3f}",
                    flush=True,
                )
            # Written before the teardown, so a teardown that fails cannot take the
            # measurement with it.
            finish(document)
            if options.out:
                Path(options.out).write_text(json.dumps(document, indent=1, sort_keys=True) + "\n")
        finally:
            empty_as_root(workdir)
    table(document)
    return 0


def empty_as_root(workdir):
    """Remove everything the container wrote under `workdir`, as root, which is who owns it.

    The corpus is built by the container, so its subdirectories (`gtrb`, the trio's CSVs) belong
    to root and the host cannot unlink what is inside them: the temporary directory's teardown
    raised on every shard after the last tool had been measured.
    """
    subprocess.run(
        [
            "docker", "run", "--rm", "--platform", run_array.PLATFORM,
            "-v", f"{workdir}:/scratch", run_array.IMAGE, "rm -rf /scratch/* /scratch/.[!.]*",
        ],
        capture_output=True,
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
