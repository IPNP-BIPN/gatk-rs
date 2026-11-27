#!/usr/bin/env python3
"""Inside broadinstitute/gatk:4.7.0.0: time the reference and the port on Train and Score.

The pinned oracle has no python, so bench.py cannot reach these two; this is their bench.
    python3 tools/python-tools/make_vets_fixture.py FX 200000  # then widen the ##contig lengths
    docker run --platform linux/amd64 -v PORT_DIR:/port -v W:/w broadinstitute/gatk:4.7.0.0 \
        python3 /w/vets.py /port/gatk /w/fx /w/out 3

Each command runs pinned to one core (taskset), alternating sides, N repetitions. CPU and the peak
RSS of the largest single process come from wait4; the peak of the whole tree's summed RSS comes
from polling /proc every 20 ms. Usage: measure.py <port binary> <fixture dir> <work dir> <reps>
"""
import hashlib, json, os, subprocess, sys, time

port, fx, work, reps = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
JAVA = "java -Xmx2g -XX:ActiveProcessorCount=1 -jar /gatk/gatk.jar"
ANN = "-A QD -A MQ -A FS -A SOR"
RES = f"--resource:train,training=true,calibration=true {fx}/resource.vcf"
VCF = "--add-output-vcf-command-line false --create-output-variant-index false --do-not-gzip-vcf-output"


def tree_rss(pid):
    kids, total = [pid], 0
    seen = set()
    while kids:
        p = kids.pop()
        if p in seen:
            continue
        seen.add(p)
        try:
            with open(f"/proc/{p}/statm") as f:
                total += int(f.read().split()[1]) * 4
            with open(f"/proc/{p}/task/{p}/children") as f:
                kids += [int(c) for c in f.read().split()]
        except OSError:
            pass
    return total


def run(cmd, cwd):
    os.makedirs(cwd, exist_ok=True)
    start = time.monotonic()
    child = subprocess.Popen(["taskset", "-c", "0", "bash", "-c", cmd], cwd=cwd,
                             stdout=subprocess.DEVNULL, stderr=open(os.path.join(cwd, "stderr"), "w"))
    peak = 0
    while True:
        pid, status, usage = os.wait4(child.pid, os.WNOHANG)
        if pid:
            break
        peak = max(peak, tree_rss(child.pid))
        time.sleep(0.02)
    return {"exit": os.waitstatus_to_exitcode(status), "wall": time.monotonic() - start,
            "cpu": usage.ru_utime + usage.ru_stime, "rss_kib": usage.ru_maxrss, "tree_kib": peak}


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest() if os.path.exists(path) else None


ref = f"{work}/ref-prep"
r = run(f"{JAVA} IndexFeatureFile -I {fx}/calls.vcf && {JAVA} IndexFeatureFile -I {fx}/resource.vcf && {JAVA} ExtractVariantAnnotations -V {fx}/calls.vcf {ANN} {RES} "
        f"--maximum-number-of-unlabeled-variants 400 -O extract {VCF} && "
        f"{JAVA} TrainVariantAnnotationsModel --annotations-hdf5 extract.annot.hdf5 -O train", ref)
assert r["exit"] == 0, open(f"{ref}/stderr").read()[-3000:]

cases = {
    "TrainVariantAnnotationsModel": (
        lambda gatk, d: f"cp {ref}/extract.annot.hdf5 . && {gatk} TrainVariantAnnotationsModel --annotations-hdf5 extract.annot.hdf5 -O train",
        ["train.snp.trainingScores.hdf5", "train.snp.calibrationScores.hdf5",
         "train.indel.trainingScores.hdf5", "train.indel.calibrationScores.hdf5"]),
    "ScoreVariantAnnotations": (
        lambda gatk, d: f"{gatk} ScoreVariantAnnotations -V {fx}/calls.vcf {ANN} {RES} --model-prefix {ref}/train -O score {VCF}",
        ["score.vcf"]),
}
out = {}
for tool, (command, files) in cases.items():
    sides = {"reference": [], "port": []}
    for i in range(reps):
        for side, gatk in (("reference", JAVA), ("port", port)):
            d = f"{work}/{tool}/{side}{i}"
            m = run(command(gatk, d), d)
            assert m["exit"] == 0, (tool, side, open(f"{d}/stderr").read()[-3000:])
            m["sha256"] = {f: sha(f"{d}/{f}") for f in files}
            sides[side].append(m)
    out[tool] = sides
    for side, ms in sides.items():
        med = lambda k: sorted(m[k] for m in ms)[len(ms) // 2]
        print(f"{tool:30} {side:9} wall {med('wall'):7.2f}s cpu {med('cpu'):7.2f}s "
              f"rss {med('rss_kib')/1024:7.1f} MiB tree {med('tree_kib')/1024:7.1f} MiB")
    same = all(m["sha256"] == sides["reference"][0]["sha256"] for s in sides.values() for m in s)
    print(f"{tool:30} sha256 equal across all runs: {same}")
json.dump(out, open(f"{work}/measure.json", "w"), indent=1)
