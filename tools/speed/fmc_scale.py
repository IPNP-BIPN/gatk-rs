#!/usr/bin/env python3
"""FilterMutectCalls at scale: the port against the reference on 67,000 real Mutect2 calls.

    python3 tools/speed/fmc_scale.py PORT_BINARY GATK_CLASSPATH [REPS]

Run it where both sides run natively, inside a container with `/work` writable (the Speed
workflow runs it in eclipse-temurin:17-jdk on x86-64 and arm64 runners). `data/fmc-chr20.vcf.gz`
is Mutect2 over NA12878, chr20:10,000,000-10,300,000; it is tiled 200 times along chr20, and the
reference is all N, since FilterMutectCalls reads no base without orientation priors. Both sides
are timed by `timer.py`; this fails unless the outputs' sha256 are equal and the port is cheaper
on CPU time and on peak RSS.
"""

import hashlib
import json
import platform
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "coverage"))
from run_array import JAVA_OPENS  # noqa: E402

port, classpath = sys.argv[1], sys.argv[2]
reps = int(sys.argv[3]) if len(sys.argv) > 3 else 3
work = Path("/work/in")
work.mkdir(parents=True, exist_ok=True)
vcf, fasta = work / "fmc.vcf", work / "chr20.fasta"
subprocess.run(
    [sys.executable, HERE / "tile_vcf.py", HERE / "data/fmc-chr20.vcf.gz",
     "10000000", "300000", "200", vcf],
    check=True,
)
length = 64444167  # chr20's, as the input's ##contig line has it
with open(fasta, "w") as out:
    out.write(">chr20\n")
    for start in range(0, length, 60):
        out.write("N" * min(60, length - start) + "\n")
(work / "chr20.fasta.fai").write_text(f"chr20\t{length}\t7\t60\t61\n")
(work / "chr20.dict").write_text(f"@HD\tVN:1.6\n@SQ\tSN:chr20\tLN:{length}\n")

arguments = f"FilterMutectCalls -V {vcf} -R {fasta} -O /work/out/f.vcf"
sides = {
    "oracle": f"java {JAVA_OPENS} -cp {classpath} org.broadinstitute.hellbender.Main {arguments}",
    "port": f"{port} {arguments}",
}
spec = work / "spec.json"
spec.write_text(json.dumps({"reps": reps, "timeout": 1800, **sides}))
timing = json.loads(
    subprocess.run([sys.executable, HERE / "timer.py", spec], check=True,
                   capture_output=True, text=True).stdout
)


def digests(command):
    subprocess.run(["bash", "-c", f"exec {command}"], check=True, capture_output=True)
    lines = Path("/work/out/f.vcf").read_bytes().splitlines(keepends=True)
    body = b"".join(line for line in lines if not line.startswith(b"##GATKCommandLine"))
    stats = Path("/work/out/f.vcf.filteringStats.tsv").read_bytes()
    return hashlib.sha256(body).hexdigest()[:16], hashlib.sha256(stats).hexdigest()[:16]


hashes = {side: digests(command) for side, command in sides.items()}
summary = {}
for side, runs in timing["cold"].items():
    kept = runs[1:]  # repetition zero warms the page cache
    if any(run["exit"] != 0 for run in kept):
        sys.exit(f"{side} exited {[run['exit'] for run in kept]}")
    summary[side] = {key: statistics.median(run[key] for run in kept)
                     for key in ("wall", "cpu", "rss_kib")}
print(f"{timing['host']['cpu'] or platform.machine()}, {timing['host']['cores']} cores, "
      f"{reps} reps, medians")
print("| | wall | CPU | peak RSS | VCF sha256 | stats sha256 |\n|---|---:|---:|---:|---|---|")
for side, name in (("oracle", "GATK 4.7.0.0"), ("port", "port")):
    m = summary[side]
    print(f"| {name} | {m['wall']:.2f} s | {m['cpu']:.2f} s | {m['rss_kib'] / 1024:.1f} MiB "
          f"| `{hashes[side][0]}` | `{hashes[side][1]}` |")
failures = []
if hashes["port"] != hashes["oracle"]:
    failures.append("outputs differ")
for key in ("cpu", "rss_kib"):
    if summary["port"][key] >= summary["oracle"][key]:
        failures.append(f"port not cheaper on {key}")
sys.exit("; ".join(failures) or None)
