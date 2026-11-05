#!/usr/bin/env python3
"""Coverage-guided differential fuzzing: seed from the covering array, fuzz the reference, replay.

1. **Seed** from the tool's covering array: its rows already cover every pair of values, so the
   search starts where combinatorial coverage stops. Each argument's DOMAIN is every value the
   array gives it, which is what a mutation draws from.
2. **Fuzz** inside the pinned container (FuzzDriver.java): one warm JVM, JaCoCo counting the
   reference's probes, a mutant kept when it reaches a probe or an outcome nothing reached before.
3. **Replay** every kept case on both sides out of process, through run_array.py, so a case is
   compared with the same canonicalization as a covering-array row. A divergence is **minimized**:
   arguments the array varies are dropped one at a time while both sides, re-run, still disagree.

    python3 tools/fuzz/run_fuzz.py --tool CountReads --iterations 60
    python3 tools/fuzz/run_fuzz.py --tool CountReads --iterations 60 --port target/release/gatk
    python3 tools/fuzz/run_fuzz.py --merge tools/fuzz/findings   # -> tools/fuzz/measured.json

Findings land in tools/fuzz/findings/ and are published by CI rather than committed; what is
committed is `measured.json`, frozen from a CI run's artefact like every other measurement.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "coverage"))
import run_array  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
FINDINGS = REPO / "tools" / "fuzz" / "findings"
MEASURED = REPO / "tools" / "fuzz" / "measured.json"
# Only the reference's own code is counted; htsjdk is measured by its own repository.
INCLUDES = "org.broadinstitute.hellbender.*:picard.*"
RNG_SEED = 20261105


def seed_file(tool, array, path):
    """The array's rows, plus what the driver needs to write them as `run_array.as_cli` would."""
    held = array.get("excluded", [])
    rows = [run_array.row_arguments(row, held) for row in array["array"]]
    domains = {}
    for row in array["array"]:
        for pair in row["arguments"]:
            name, _, value = pair.partition("=")
            domains.setdefault(name, [])
            if value not in domains[name]:
                domains[name].append(value)
    lines = [
        "\t".join(["positional", *run_array.positional_values(tool)]),
        "\t".join(["tagged", *run_array.tagged_arguments(tool)]),
        "\t".join(["lists", *run_array.list_arguments(tool)]),
    ]
    lines += ["\t".join(["domain", name, *values]) for name, values in domains.items()]
    lines += ["\t".join(["row", *args]) for args in rows]
    path.write_text("\n".join(lines) + "\n")
    return len(rows), set(domains)


def fuzz(tool, workdir, iterations, timeout):
    """Run the driver in the container; returns its SUMMARY fields and the kept cases."""
    harness = workdir / "fuzz"
    command = (
        "mkdir -p /work/tmp /work/tmp2 /work/out /tmp/fz && cp /harness/FuzzDriver.java /tmp/fz/ "
        '&& javac -cp "$ORACLE_CP:$JACOCO_AGENT:$JACOCO_ANALYSIS" -d /tmp/fz /tmp/fz/FuzzDriver.java '
        f"&& java {run_array.JAVA_OPENS} -javaagent:$JACOCO_AGENT=output=none,includes={INCLUDES} "
        '-cp "/tmp/fz:$ORACLE_CP:$JACOCO_AGENT:$JACOCO_ANALYSIS" '
        f"FuzzDriver {tool} /work/fuzz/seed.txt {iterations} /work/fuzz/corpus.txt {RNG_SEED} {timeout}"
    )
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", run_array.PLATFORM,
            "-v", f"{REPO}/tools/fuzz:/harness:ro",
            "-v", f"{run_array.fixtures_of(workdir)}:/work/fixtures:ro",
            "-v", f"{workdir / 'out'}:/work/out",
            "-v", f"{harness}:/work/fuzz",
            "-w", "/work", run_array.IMAGE, command,
        ],
        capture_output=True,
        text=True,
    )
    print(result.stdout[-3000:])
    summary = [line for line in result.stdout.splitlines() if line.startswith("SUMMARY ")]
    if result.returncode != 0 or not summary:
        print(result.stderr[-3000:])
        raise SystemExit(f"{tool}: the fuzz driver failed")
    fields = dict(word.split("=", 1) for word in summary[-1].split()[1:])
    kept = []
    for line in (harness / "corpus.txt").read_text().splitlines():
        args, _, rest = line.partition("\t#\t")
        reached, fresh = rest.rsplit("\t", 1)
        kept.append({"arguments": args.split("\t"), "reached": reached, "new_probes": int(fresh)})
    return fields, kept


def both_sides(tool, args, port, workdir, config):
    positional, tagged, lists, keep = config
    code, text, error = run_array.run_oracle(tool, args, workdir, positional, tagged, lists)
    reference = run_array.outcome(code, text, error, keep)
    code, text, error = run_array.run_port(port, tool, args, workdir, positional, tagged, lists)
    return reference, run_array.outcome(code, text, error, keep)


def minimize(tool, args, varied, port, workdir, config):
    """Drop varied arguments while the divergence survives, re-running BOTH sides each step:
    dropping an argument changes what the reference does too (picard-rs learned this the hard way).
    """
    current = list(args)
    changed = True
    while changed:
        changed = False
        for arg in list(current):
            if arg.split("=", 1)[0] not in varied:
                continue
            candidate = [a for a in current if a != arg]
            reference, ours = both_sides(tool, candidate, port, workdir, config)
            if reference != ours:
                current, changed = candidate, True
                break
    return current, *both_sides(tool, current, port, workdir, config)


def run(options):
    array_path = run_array.ARRAYS / f"{options.tool}.t{options.t}.json"
    if not array_path.exists():
        raise SystemExit(f"no array at {array_path}; generate it with covering.py")
    array = json.loads(array_path.read_text())
    run_array.FIXTURES_DIR = options.fixtures_dir
    config = (
        run_array.positional_values(options.tool),
        run_array.tagged_arguments(options.tool),
        run_array.list_arguments(options.tool),
        run_array.output_on_failure(options.tool),
    )

    with tempfile.TemporaryDirectory(prefix="gatk-fuzz-", ignore_cleanup_errors=True) as tmp:
        workdir = Path(tmp)
        for name in ("out", "fuzz", "tmp"):
            (workdir / name).mkdir()
        try:
            if run_array.FIXTURES_DIR is None:
                run_array.build_fixtures(workdir)
            seeds, varied = seed_file(options.tool, array, workdir / "fuzz" / "seed.txt")
            print(f"{options.tool}: {seeds} seed rows from the t={options.t} array")
            fields, kept = fuzz(options.tool, workdir, options.iterations, options.timeout)

            divergences = []
            if options.port:
                for case in kept:
                    reference, ours = both_sides(options.tool, case["arguments"], options.port, workdir, config)
                    case["matched"] = reference == ours
                    if not case["matched"]:
                        minimal, ref_min, ours_min = minimize(
                            options.tool, case["arguments"], varied, options.port, workdir, config)
                        divergences.append({"arguments": minimal, "reference": ref_min, "port": ours_min})
        finally:
            run_array.empty_output_in_container(workdir)

    summary = {
        "tool": options.tool,
        "t": options.t,
        "iterations": options.iterations,
        "seeds": seeds,
        "probes_from_seeds": int(fields["probes_from_seeds"]),
        "probes_total": int(fields["probes_total"]),
        "kept": len(kept),
        "timed_out": fields["timed_out"] == "true",
    }
    if options.port:
        summary["replayed"] = len(kept)
        summary["diverged"] = len(divergences)

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / f"{options.tool}.json").write_text(json.dumps(
        {"summary": summary, "kept": kept, "divergences": divergences}, indent=1) + "\n")
    print(json.dumps(summary))
    for d in divergences[:5]:
        print(f"  diverges: {' '.join(d['arguments'])}")
        print(f"    reference={d['reference'][:200]!r}")
        print(f"    port     ={d['port'][:200]!r}")
    return 0


def merge(directory):
    """Fold per-tool findings into the committed measurement the reports read."""
    tools = {}
    for path in sorted(Path(directory).glob("*.json")):
        summary = json.loads(path.read_text())["summary"]
        if "failed" in summary:
            print(f"not merged, the driver failed: {summary['tool']}: {summary['failed'][:120]}")
            continue
        tools[summary.pop("tool")] = summary
    MEASURED.write_text(json.dumps({
        "$comment": [
            "Produced by tools/fuzz/run_fuzz.py --merge from a CI run's fuzz-findings artefact.",
            "`probes_*` count the reference's JaCoCo probes reached (a proxy for branches, not a",
            "percentage); `diverged` counts kept cases on which the port, re-run out of process,",
            "did not answer exactly as the reference did, after minimization.",
        ],
        "tools": tools,
    }, indent=1, sort_keys=True) + "\n")
    print(f"wrote {MEASURED} with {len(tools)} tools")
    return 0


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tool")
    parser.add_argument("--t", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=60)
    parser.add_argument("--timeout", type=int, default=120, help="seconds per in-process case")
    parser.add_argument("--port", help="the port binary (linux/amd64); omit to fuzz the reference alone")
    parser.add_argument("--fixtures-dir", help="a corpus MakeFixtures.java already built")
    parser.add_argument("--merge", metavar="DIR", help="merge per-tool findings into measured.json")
    parser.add_argument("--shard", metavar="K/N",
                        help="fuzz the K-th of N slices of every tool whose array the port is measured on")
    options = parser.parse_args(argv)
    if options.merge:
        return merge(options.merge)
    if options.shard:
        return shard(options)
    if not options.tool:
        parser.error("--tool or --shard is required")
    return run(options)


def shard(options):
    """One CI slice: every tool `measured.json` holds an argument-coverage number for, so every
    tool fuzzed is one whose array already runs against both sides. A tool the driver cannot run is
    recorded as failed and the slice goes on; the report then reads `not measured` for it."""
    k, n = (int(x) for x in options.shard.split("/"))
    measured = json.loads((REPO / "tools" / "coverage" / "measured.json").read_text())["tools"]
    tools = sorted(t for t in measured if (run_array.ARRAYS / f"{t}.t{options.t}.json").exists())
    mine = tools[k::n]
    failed = []
    with tempfile.TemporaryDirectory(prefix="gatk-fuzz-fixtures-", ignore_cleanup_errors=True) as tmp:
        if options.fixtures_dir is None:
            options.fixtures_dir = str(run_array.build_fixtures(Path(tmp)))
        for tool in mine:
            options.tool = tool
            try:
                run(options)
            except (SystemExit, Exception) as error:  # noqa: BLE001 -- one tool must not end the slice
                failed.append(tool)
                FINDINGS.mkdir(parents=True, exist_ok=True)
                (FINDINGS / f"{tool}.json").write_text(json.dumps(
                    {"summary": {"tool": tool, "failed": str(error)[:500]}}, indent=1) + "\n")
                print(f"{tool}: FAILED {error}")
    print(f"shard {k}/{n}: {len(mine) - len(failed)} of {len(mine)} fuzzed; failed: {' '.join(failed) or 'none'}")
    return 1 if mine and len(failed) == len(mine) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
