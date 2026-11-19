#!/usr/bin/env python3
"""Measure what the port costs against the reference, on rows it is already proven identical on.

    python3 tools/speed/bench.py --port target/release/gatk --out speed.json
    python3 tools/speed/bench.py --port target/release/gatk --tool CountReads --tool FlagStat
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
import shutil
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
        if "steps" in part:
            document.setdefault("steps", {}).update(part["steps"])
            for key in ("scale", "commit", "host"):
                document.setdefault(key, part.get(key))
            continue
        document["tools"].update(part["tools"])
        document["shards"].append(
            {key: part.get(key) for key in ("host", "container_overhead", "reps", "iterations")}
        )
        for key in ("reps", "iterations", "commit"):
            document.setdefault(key, part.get(key))
    if "steps" in document:
        del document["tools"], document["shards"]
        Path(out).write_text(json.dumps(document, indent=1, sort_keys=True) + "\n")
        scale_table(document)
        return 0
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


SCALE = SPEED / "scale.json"
# The masks the comparison applies to a whole output directory, applied here one file at a time.
SCALE_TEXT = (".vcf", ".vcf.gz", ".tsv", ".table", ".metrics", ".interval_list", ".seg", ".tranches",
              ".recal", ".param", ".stats", ".hdf5")


def scale_steps(names):
    spec = json.loads(SCALE.read_text())
    steps = spec["steps"]
    unknown = set(names or ()) - {s["name"] for s in steps}
    if unknown:
        raise SystemExit(f"no such step in scale.json: {', '.join(sorted(unknown))}")
    return [s for s in steps if not names or s["name"] in names]


def scale_command(step, side, port, threads):
    """One side's command line, with the step's placeholders resolved inside the container."""
    args = step["args"].format(data="/work/data", out="/work/out", threads=threads, **{"in": "/work/in"})
    if side == "oracle":
        java = f"java {run_array.JAVA_OPENS} {step.get('java_options', '')}".strip()
        return f'{java} -cp "$ORACLE_CP" org.broadinstitute.hellbender.Main {step["name"]} {args}'
    return f"/work/port-binary/{Path(port).name} {step['name']} {args}"


def time_step(step, options):
    """Both sides of one step in the pinned container, under the memory limit, kept on disk."""
    results = Path(options.results).resolve()
    missing = [n for n in step["needs"] if not any((results / n / "oracle").glob("*"))]
    if missing:
        return None, f"no reference output from {', '.join(missing)} to run on"
    here = results / step["name"]
    for side in ("oracle", "port"):
        (here / side).mkdir(parents=True, exist_ok=True)
    binary = Path(options.port).resolve()
    spec = {
        "scale": True,
        "reps": 0,
        "timeout": options.timeout or step.get("timeout", 10800),
        "out": {side: f"/work/result/{side}" for side in ("oracle", "port")},
        "warm": [f"/work/in/{n}/{f.name}" for n in step["needs"] for f in (results / n / "oracle").iterdir()
                 if f.is_file()],
        "oracle": scale_command(step, "oracle", binary, options.threads),
        "port": scale_command(step, "port", binary, options.threads),
    }
    (here / "spec.json").write_text(json.dumps(spec, indent=1))
    mounts = []
    for need in step["needs"]:
        mounts += ["-v", f"{results / need / 'oracle'}:/work/in/{need}:ro"]
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", run_array.PLATFORM,
            f"--memory={options.memory}", f"--memory-swap={options.memory}",
            f"--cpus={options.threads}",
            "-e", f"RAYON_NUM_THREADS={options.threads}",
            "-v", f"{Path(options.data).resolve()}:/work/data:ro",
            *mounts,
            "-v", f"{here}:/work/result",
            "-v", f"{binary.parent}:/work/port-binary:ro",
            "-v", f"{SPEED}:/harness:ro",
            "-w", "/work", run_array.IMAGE,
            "python3 /harness/timer.py /work/result/spec.json; chmod -R a+rwX /work/result",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return None, f"timer failed: {(result.stderr or result.stdout).strip()[-300:]}"
    raw = json.loads(result.stdout.strip().splitlines()[-1])
    return raw, None


def one_file(path):
    """`run_array.read_output` on one file alone: the conformance masks, file by file."""
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / path.name).symlink_to(path.resolve())
        return run_array.read_output(Path(tmp))


def first_difference(a, b):
    """Where two outputs first part: a line for text (gzipped or not), a byte offset otherwise."""
    import gzip

    def opener(path):
        with open(path, "rb") as handle:
            gz = handle.read(2) == b"\x1f\x8b"
        return gzip.open(path, "rb") if gz else open(path, "rb")

    if a.name.endswith(SCALE_TEXT) or a.name.endswith(".gz") and a.name[:-3].endswith(SCALE_TEXT):
        with opener(a) as x, opener(b) as y:
            number = 0
            for number, (left, right) in enumerate(zip(x, y), 1):
                if left != right:
                    return {"line": number, "reference": left[:300].decode(errors="replace"),
                            "port": right[:300].decode(errors="replace")}
            return {"line": number + 1, "note": "one file ends first"} if (x.read(1) or y.read(1)) else None
    with opener(a) as x, opener(b) as y:
        offset = 0
        while True:
            left, right = x.read(1 << 20), y.read(1 << 20)
            if left != right:
                at = next((i for i, (p, q) in enumerate(zip(left, right)) if p != q), min(len(left), len(right)))
                return {"decompressed_offset": offset + at}
            if not left:
                return None
            offset += len(left)


def compare_step(here):
    """Per file: both sides' sha256, whether they are equal under the masks, where they part."""
    oracle, port = here / "oracle", here / "port"
    names = sorted({p.relative_to(oracle).as_posix() for p in oracle.rglob("*") if p.is_file()}
                   | {p.relative_to(port).as_posix() for p in port.rglob("*") if p.is_file()})
    files = {}
    for name in names:
        a, b = oracle / name, port / name
        entry = {side: sha256_of(p) for side, p in (("oracle", a), ("port", b)) if p.exists()}
        if not (a.exists() and b.exists()):
            entry["equal"] = False
            entry["missing_on"] = "port" if a.exists() else "reference"
        else:
            entry["equal"] = entry["oracle"] == entry["port"] or one_file(a) == one_file(b)
            if not entry["equal"]:
                entry["first_difference"] = first_difference(a, b)
        files[name] = entry
    return files


def sha256_of(path):
    import hashlib
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def log_tail(path, lines=5):
    return path.read_text(errors="replace").strip().splitlines()[-lines:] if path.exists() else []


def scale(options):
    """Milestone D, D2: the best-practices chain on chr20, step by step, both sides at scale."""
    results = Path(options.results)
    results.mkdir(parents=True, exist_ok=True)
    document = {
        "scale": {"data": json.loads(SCALE.read_text())["sources"]["reads"]["what"],
                  "threads": options.threads, "memory": options.memory, "timeout": options.timeout},
        "commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                                 cwd=REPO).stdout.strip(),
        "steps": {},
    }
    for step in scale_steps(options.step):
        name = step["name"]
        print(f"{name}: running", flush=True)
        raw, reason = time_step(step, options)
        if raw is None:
            document["steps"][name] = {"reason": reason}
            print(f"{name}: not measured, {reason}", flush=True)
            continue
        here = results / name
        entry = {"oracle": raw["cold"]["oracle"][0], "port": raw["cold"]["port"][0], "host": raw["host"]}
        for side, run in (("oracle", entry["oracle"]), ("port", entry["port"])):
            if run["exit"] != 0:
                run["log_tail"] = log_tail(here / f"{side}.log")
        if entry["oracle"]["exit"] != 0:
            # Nothing below may run on a reference output that did not finish.
            shutil.rmtree(here / "oracle", ignore_errors=True)
        if entry["oracle"]["exit"] == 0 and entry["port"]["exit"] == 0:
            entry["files"] = compare_step(here)
            entry["bytes_equal"] = all(f["equal"] for f in entry["files"].values())
        document["steps"][name] = entry
        document.setdefault("host", raw["host"])
        print(f"{name}: {scale_row(name, entry)}", flush=True)
        if options.out:
            Path(options.out).write_text(json.dumps(document, indent=1, sort_keys=True) + "\n")
        if not options.keep_port_output:
            shutil.rmtree(here / "port", ignore_errors=True)
    if options.out:
        Path(options.out).write_text(json.dumps(document, indent=1, sort_keys=True) + "\n")
    scale_table(document)
    return 0


def side_state(run):
    if run["exit"] == 0:
        return None
    if run.get("oom_killed"):
        return "killed: out of memory"
    if run.get("timed_out"):
        return "killed: timeout"
    return f"exit {run['exit']}"


def scale_row(name, entry):
    if "reason" in entry:
        return f"| {name} | not run: {entry['reason']} |"
    o, p = entry["oracle"], entry["port"]

    def cell(run, key):
        state = side_state(run)
        if key == "rss_kib":
            value = f"{run[key] / 1048576:.2f} GiB"
        else:
            value = f"{run[key]:.0f} s"
        return f"{value}{' (' + state + ')' if state and key == 'wall' else ''}"

    ratio_wall = p["wall"] / o["wall"] if o["wall"] and not side_state(o) and not side_state(p) else None
    equal = entry.get("bytes_equal")
    verdict = "-" if equal is None else ("yes" if equal else "**no**")
    return (f"| {name} | {cell(o, 'wall')} | {cell(p, 'wall')} | "
            f"{f'{ratio_wall:.2f}' if ratio_wall else '-'} | {cell(o, 'cpu')} | {cell(p, 'cpu')} | "
            f"{cell(o, 'rss_kib')} | {cell(p, 'rss_kib')} | {verdict} |")


def scale_table(document):
    print("\n| step | reference wall | port wall | port / ref | reference CPU | port CPU | "
          "reference peak RSS | port peak RSS | bytes equal |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|:-:|")
    order = [s["name"] for s in json.loads(SCALE.read_text())["steps"]]
    for name in sorted(document["steps"], key=lambda n: order.index(n) if n in order else len(order)):
        print(scale_row(name, document["steps"][name]))


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--port", help="the linux/amd64 port binary (gatk-rs)")
    parser.add_argument("--tool", action="append", help="measure this tool (default: every proven one)")
    parser.add_argument("--shard", help="i/n: measure every n-th proven tool from the i-th")
    parser.add_argument("--reps", type=int, default=5, help="cold repetitions per side, after a warm-up")
    parser.add_argument("--iterations", type=int, default=20, help="in-JVM iterations; 0 skips the steady state")
    parser.add_argument("--max-probe", type=int, default=6, help="accepted rows to try before giving up")
    parser.add_argument("--timeout", type=int, help="seconds per run (default 600; with --scale, the step's own)")
    parser.add_argument("--fixtures-dir", help="a corpus MakeFixtures.java already built")
    parser.add_argument("--out", help="write the measurement here")
    parser.add_argument("--merge", nargs="+", help="combine shard outputs into one --out")
    parser.add_argument("--compare", nargs=2, help="two measurements: how far each ratio moved")
    parser.add_argument("--scale", action="store_true", help="the chr20 chain in scale.json, not the fixtures")
    parser.add_argument("--step", action="append", help="with --scale: run this step (default: all)")
    parser.add_argument("--data", help="with --scale: the directory scale_data.py derived")
    parser.add_argument("--results", default="scale-results", help="with --scale: each step's outputs")
    parser.add_argument("--threads", type=int, default=4, help="with --scale: threads, both sides")
    parser.add_argument("--memory", default="14g", help="with --scale: the container's memory limit")
    parser.add_argument("--keep-port-output", action="store_true", help="with --scale: keep the port's files")
    options = parser.parse_args(argv)

    if options.merge:
        return merge(options.merge, options.out or "speed.json")
    if options.compare:
        return compare(*options.compare)
    if not options.port:
        parser.error("--port is required to measure")
    if options.scale:
        if not options.data:
            parser.error("--scale needs --data")
        return scale(options)
    options.timeout = options.timeout or 600

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
