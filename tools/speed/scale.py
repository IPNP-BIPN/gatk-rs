#!/usr/bin/env python3
"""Port against GATK at multi-GB scale, equal threads. Run by `speed-scale.yml`.

    python3 tools/speed/scale.py --gatk /opt/gatk/gatk --port target/release/gatk \
        --bam /work/in.bam --threads 4 --out scale.json
    python3 tools/speed/scale.py --compare x86_64.json aarch64.json

`bench.py` measures kilobyte covering-array rows; this measures one synthetic BAM of gigabytes
(`scale_bam.py`) through the tools whose memory the streaming work bounds. Each run goes through
`timer.run`, so CPU and peak RSS come from `wait4` as they do there, and the two sides alternate
per repetition the same way. Outputs land in /work/out, which `timer.prepare` empties before each
run, so each output is hashed as soon as its run ends.

Threads are equal: GATK's Spark tools get `--spark-master local[N]`, the port `RAYON_NUM_THREADS=N`.
CountReads and PrintReads are GATK's walkers, which use one thread whatever N is.
"""

import argparse
import hashlib
import json
import os
import platform
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import timer  # noqa: E402

OUT_BAM = "/work/out/out.bam"
SORT = "SortSamSpark --input {bam} --output " + OUT_BAM + " --sort-order {order} --num-reducers 1"
TOOLS = [
    ("SortSamSpark queryname", SORT.replace("{order}", "queryname"), True),
    ("SortSamSpark coordinate", SORT.replace("{order}", "coordinate"), True),
    ("CountReads", "CountReads --input {bam}", False),
    ("PrintReads", "PrintReads --input {bam} --output " + OUT_BAM, False),
]


def md5(path):
    if not os.path.exists(path):
        return None
    digest = hashlib.md5()
    with open(path, "rb") as f:
        while block := f.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def host():
    info = Path("/proc/cpuinfo").read_text()
    fields = {}
    for line in info.splitlines():
        key, _, value = line.partition(":")
        fields.setdefault(key.strip(), value.strip())
    flags = set((fields.get("flags") or fields.get("Features") or "").split())
    wanted = ("avx2", "avx512f", "avx512bw", "avx512vl", "bmi2", "pclmulqdq", "asimd", "sve", "sve2", "pmull", "crc32")
    return {
        "arch": platform.machine(),
        "cpu": fields.get("model name") or f"implementer {fields.get('CPU implementer')} part {fields.get('CPU part')}",
        "cores": os.cpu_count(),
        "flags": [f for f in wanted if f in flags],
    }


def measure(args):
    n = args.threads
    gatk_extra = f"--java-options '{args.java_options}' " if args.java_options else ""
    rows = []
    for name, template, spark in TOOLS:
        line = template.format(bam=args.bam)
        sides = {
            "gatk": f"{args.gatk} {gatk_extra}{line} --tmp-dir /work/tmp"
            + (f" --spark-master 'local[{n}]' --conf spark.local.dir=/work/tmp" if spark else ""),
            "port": f"env RAYON_NUM_THREADS={n} {args.port} {line}" + (" --tmp-dir /work/tmp" if spark else ""),
        }
        runs = {"gatk": [], "port": []}
        for rep in range(args.reps):
            for side in ("gatk", "port") if rep % 2 == 0 else ("port", "gatk"):
                shutil.rmtree("/work/tmp", ignore_errors=True)
                result = timer.run(sides[side], args.timeout)
                result["md5"] = md5(OUT_BAM) if "--output" in line else None
                runs[side].append(result)
                print(name, side, json.dumps(result), flush=True)
        rows.append({"tool": name, "commands": sides, "runs": runs})
    return {"host": host(), "threads": n, "bam_bytes": os.path.getsize(args.bam), "rows": rows}


def summary(side_runs):
    ok = [r for r in side_runs if r["exit"] == 0]
    if not ok:
        return None
    return {
        "wall": min(r["wall"] for r in ok),
        "cpu": min(r["cpu"] for r in ok),
        "rss_gib": max(r["rss_kib"] for r in ok) / (1 << 20),
        "md5s": sorted({r["md5"] for r in ok if r["md5"]}),
    }


def table(result):
    h = result["host"]
    lines = [
        f"**{h['arch']}**: {h['cpu']}, {h['cores']} cores, flags {' '.join(h['flags']) or 'none'}; "
        f"{result['threads']} threads, input {result['bam_bytes'] / (1 << 30):.2f} GiB",
        "",
        "| tool | GATK wall s | GATK CPU s | GATK RSS GiB | port wall s | port CPU s | port RSS GiB | CPU port/GATK | RSS port/GATK | port md5 |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in result["rows"]:
        g, p = summary(row["runs"]["gatk"]), summary(row["runs"]["port"])
        cell = lambda s, k: f"{s[k]:.1f}" if s else "failed"  # noqa: E731
        rss = lambda s: f"{s['rss_gib']:.2f}" if s else "failed"  # noqa: E731
        ratio = lambda k: f"{p[k] / g[k]:.2f}" if g and p else "-"  # noqa: E731
        md5s = ", ".join(p["md5s"]) if p and p["md5s"] else "-"
        lines.append(
            f"| {row['tool']} | {cell(g, 'wall')} | {cell(g, 'cpu')} | {rss(g)} | {cell(p, 'wall')} "
            f"| {cell(p, 'cpu')} | {rss(p)} | {ratio('cpu')} | {ratio('rss_gib')} | {md5s} |"
        )
    return "\n".join(lines)


def compare(paths):
    results = [json.loads(Path(p).read_text()) for p in paths]
    out = [table(r) + "\n" for r in results]
    out.append("| tool | " + " | ".join(r["host"]["arch"] for r in results) + " | port bytes identical |")
    out.append("|---|" + "---|" * len(results) + "---|")
    same = True
    for i, row in enumerate(results[0]["rows"]):
        sums = [summary(r["rows"][i]["runs"]["port"]) for r in results]
        md5s = [set(s["md5s"]) if s else set() for s in sums]
        if row["commands"]["port"].find("--output") < 0:
            continue
        equal = all(len(m) == 1 for m in md5s) and len(set.union(*md5s)) == 1
        same &= equal
        out.append(f"| {row['tool']} | " + " | ".join(", ".join(sorted(m)) or "-" for m in md5s) + f" | {'yes' if equal else 'NO'} |")
    print("\n".join(out))
    return 0 if same else 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gatk")
    parser.add_argument("--port")
    parser.add_argument("--bam")
    parser.add_argument("--threads", type=int, default=os.cpu_count())
    parser.add_argument("--reps", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=7200)
    parser.add_argument("--java-options", default="")
    parser.add_argument("--out")
    parser.add_argument("--compare", nargs="+")
    args = parser.parse_args()
    if args.compare:
        sys.exit(compare(args.compare))
    result = measure(args)
    Path(args.out).write_text(json.dumps(result, indent=1))
    print(table(result))


if __name__ == "__main__":
    main()
