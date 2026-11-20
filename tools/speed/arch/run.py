#!/usr/bin/env python3
"""Time CountReads on SAM text and CRAM, GATK against the port, through tools/speed/timer.py.
Both sides run under `taskset -c 0`. Prints a markdown table and the outputs' sha256.
Usage: run.py <data-dir> <gatk-jar> <port> <reps>"""
import hashlib, json, subprocess, sys
from pathlib import Path

data, jar, port, reps = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
here = Path(__file__).resolve().parent
sys.path.insert(0, str(here.parent.parent / "coverage"))
import run_array  # noqa: E402

java = f"java {run_array.JAVA_OPENS} -jar {jar}"
cases = {
    "SAM text": f"--input {data}/reads.sam",
    "CRAM + FASTA": f"--input {data}/reads.cram -R {data}/ref.fasta",
}
median = lambda xs: sorted(xs)[len(xs) // 2]
rows = []
for name, args in cases.items():
    spec = {"reps": reps, "timeout": 3600, "iterations": 0, "steady": None,
            "oracle": f"taskset -c 0 {java} CountReads {args}",
            "port": f"taskset -c 0 {port} CountReads {args}"}
    Path("/work/spec.json").write_text(json.dumps(spec))
    raw = json.loads(subprocess.run([sys.executable, str(here.parent / "timer.py"), "/work/spec.json"],
                                    capture_output=True, text=True, check=True).stdout.splitlines()[-1])
    cold = raw["cold"]
    for side in ("oracle", "port"):
        bad = [r for r in cold[side] if r["exit"] != 0]
        if bad:
            raise SystemExit(f"{name}: {side} exited {bad[0]['exit']}")
    o, p = cold["oracle"][1:], cold["port"][1:]
    sha = {}
    for side, cmd in (("java", f"{java} CountReads {args}"), ("port", f"{port} CountReads {args}")):
        target = Path(f"/work/{side}.txt")
        target.unlink(missing_ok=True)
        stdout = subprocess.run(f"{cmd} -O {target}", shell=True, capture_output=True, text=True, check=True).stdout
        returned = stdout.split("Tool returned:\n", 1)[1].splitlines()[0]
        sha[side] = hashlib.sha256(target.read_bytes() + b"|" + returned.encode()).hexdigest()[:16]
    rows.append((name,
                 median([r["cpu"] for r in o]), median([r["cpu"] for r in p]),
                 median([r["rss_kib"] for r in o]) / 1024, median([r["rss_kib"] for r in p]) / 1024,
                 median([r["wall"] for r in o]), median([r["wall"] for r in p]),
                 sha["java"], sha["port"]))
print(f"host: {raw['host']}")
print("| Input | GATK CPU | port CPU | GATK RSS | port RSS | GATK wall | port wall | sha256 GATK | sha256 port |")
print("|---|---:|---:|---:|---:|---:|---:|---|---|")
for r in rows:
    print(f"| {r[0]} | {r[1]:.2f} s | {r[2]:.2f} s | {r[3]:.0f} MB | {r[4]:.0f} MB | {r[5]:.2f} s | {r[6]:.2f} s | `{r[7]}` | `{r[8]}` |")
ok = all(r[2] < r[1] and r[4] < r[3] and r[7] == r[8] for r in rows)
print("port cheaper on CPU and RSS, outputs equal:", ok)
sys.exit(0 if ok else 1)
