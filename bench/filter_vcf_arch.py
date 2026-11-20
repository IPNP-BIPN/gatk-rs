#!/usr/bin/env python3
"""FilterVcf, GATK 4.7.0.0's jar against the port, on one host: CPU, peak RSS, output sha256.

    python3 filter_vcf_arch.py --jar gatk.jar --port target/release/gatk --work DIR [--reps 5]

The inputs are generated here from a fixed seed, so every host times the same bytes: a sites-only
VCF with FS and QD, and the same sites with 12 samples (GT:AD:DP:GQ). Each repetition runs both
sides once, Java first on even repetitions and the port first on odd ones; repetition 0 warms the
page cache and is not summarised. CPU is user + sys of the child from wait4, RSS its ru_maxrss.
"""
import argparse, hashlib, json, os, platform, random, shutil, statistics, subprocess, sys, time
from pathlib import Path

PROPS = ["-Duser.language=en", "-Duser.country=US", "-Dsamjdk.use_libdeflate=false"]
ARGS = ["--MIN_AB", "0.2", "--MAX_FS", "30", "--MIN_QD", "2", "--MIN_GQ", "20", "--MIN_DP", "10"]
INFO = [("AC", "Integer", "A"), ("AF", "Float", "A"), ("AN", "Integer", "1"),
        ("BaseQRankSum", "Float", "1"), ("DP", "Integer", "1"), ("FS", "Float", "1"),
        ("MQ", "Float", "1"), ("MQRankSum", "Float", "1"), ("QD", "Float", "1"),
        ("ReadPosRankSum", "Float", "1"), ("SOR", "Float", "1"), ("VQSLOD", "Float", "1")]


def generate(path, records, samples, seed=1488):
    rng = random.Random(seed)
    names = [f"S{i:02d}" for i in range(1, samples + 1)]
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n")
        for key, kind, number in INFO:
            out.write(f'##INFO=<ID={key},Number={number},Type={kind},Description="{key}">\n')
        if samples:
            out.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
                      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
                      '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                  + ("\tFORMAT\t" + "\t".join(names) if samples else "") + "\n")
        pos = 60000
        for n in range(records):
            pos += rng.randint(1, 180)
            ref = rng.choice("ACGT")
            alts = [b for b in "ACGT" if b != ref]
            nalt = 2 if rng.random() < 0.05 else 1
            alt = rng.sample(alts, nalt)
            af = [round(rng.random() * 0.5, 5) for _ in alt]
            info = (f"AC={','.join(str(rng.randint(1, 200)) for _ in alt)};AF={','.join(map(str, af))};"
                    f"AN={rng.randint(200, 5000)};BaseQRankSum={rng.uniform(-5, 5):.3f};DP={rng.randint(100, 9000)};"
                    f"FS={rng.expovariate(1 / 8):.3f};MQ={rng.uniform(40, 70):.2f};MQRankSum={rng.uniform(-3, 3):.3f};"
                    f"QD={rng.uniform(0, 30):.2f};ReadPosRankSum={rng.uniform(-3, 3):.3f};SOR={rng.uniform(0, 4):.3f};"
                    f"VQSLOD={rng.uniform(-10, 20):.4f}")
            line = f"chr20\t{pos}\t{'rs%d' % n if rng.random() < 0.4 else '.'}\t{ref}\t{','.join(alt)}\t{rng.uniform(20, 3000):.2f}\tPASS\t{info}"
            if samples:
                gts = []
                for _ in names:
                    r = rng.random()
                    if r < 0.05:
                        gts.append("./."); continue
                    a, b = (0, 0) if r < 0.6 else (0, rng.randint(1, nalt)) if r < 0.9 else (rng.randint(1, nalt),) * 2
                    dp = rng.randint(0, 60)
                    ad = [0] * (nalt + 1)
                    if a == b: ad[a] = dp
                    else:
                        k = rng.randint(0, dp); ad[a] = k; ad[b] = dp - k
                    gts.append(f"{a}/{b}:{','.join(map(str, ad))}:{dp}:{rng.randint(0, 99)}")
                line += "\tGT:AD:DP:GQ\t" + "\t".join(gts)
            out.write(line + "\n")


def run(command, out_dir):
    shutil.rmtree(out_dir, ignore_errors=True); out_dir.mkdir(parents=True)
    start = time.perf_counter()
    child = subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    _, status, usage = os.wait4(child.pid, 0)
    child.returncode = os.waitstatus_to_exitcode(status)
    assert child.returncode == 0, (command, child.returncode)
    rss = usage.ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024)  # MiB
    return {"wall": time.perf_counter() - start, "cpu": usage.ru_utime + usage.ru_stime, "rss": rss}


def masked_index(raw):
    """The Tribble index with the indexed file's last-modified time (8 bytes after the path and size) zeroed."""
    end = raw.find(b"\0", 12)
    return raw[:end + 9] + b"\0" * 8 + raw[end + 17:]


def sha256(path):
    # Streamed: a file read whole here stays resident in this process, and a child forked from it
    # then reports that memory as its own ru_maxrss.
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def digests(out_dir):
    return {f.name: hashlib.sha256(masked_index(f.read_bytes())).hexdigest() if f.suffix == ".idx" else sha256(f)
            for f in sorted(out_dir.iterdir())}


def cpu_name():
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith(("model name", "Model name")):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or platform.machine()


def cpu_flags():
    try:
        text = Path("/proc/cpuinfo").read_text()
    except OSError:
        return []
    flags = set()
    for line in text.splitlines():
        if line.startswith(("flags", "Features")):
            flags |= set(line.split(":", 1)[1].split())
    return sorted(f for f in flags if f in {"avx", "avx2", "avx512f", "avx512bw", "sse4_2", "asimd", "sve", "sve2"})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--jar", required=True); p.add_argument("--port", required=True)
    p.add_argument("--work", required=True); p.add_argument("--reps", type=int, default=5)
    p.add_argument("--records", type=int, default=600000)
    o = p.parse_args()
    work = Path(o.work).resolve(); work.mkdir(parents=True, exist_ok=True)
    out = work / "out"
    rows, result = [], {"cpu": cpu_name(), "flags": cpu_flags(), "machine": platform.machine(), "cores": os.cpu_count(), "inputs": {}}
    for label, samples in (("sites", 0), ("12 samples", 12)):
        vcf = work / f"{label.replace(' ', '_')}.vcf"
        if not vcf.exists():
            generate(vcf, o.records, samples)
        result["inputs"][label] = {"bytes": vcf.stat().st_size, "sha256": sha256(vcf)}
        common = ["FilterVcf", "--INPUT", str(vcf), "--OUTPUT", str(out / "out.vcf"), *ARGS]
        sides = {"java": ["java", *PROPS, "-jar", o.jar, *common], "port": [o.port, *common]}
        runs = {"java": [], "port": []}
        for rep in range(o.reps + 1):
            for side in (("java", "port") if rep % 2 == 0 else ("port", "java")):
                runs[side].append(run(sides[side], out))
        sha = {}
        for side in ("java", "port"):
            run(sides[side], out); sha[side] = digests(out)
        med = {s: {k: statistics.median(r[k] for r in runs[s][1:]) for k in ("wall", "cpu", "rss")} for s in runs}
        equal = sha["java"] == sha["port"]
        result["inputs"][label].update(median=med, sha256=sha, equal=equal)
        rows.append(f"| {label} ({result['inputs'][label]['bytes'] / 1e6:.0f} MB) | {med['java']['wall']:.2f} s | {med['java']['cpu']:.2f} s | {med['java']['rss']:.1f} MiB "
                    f"| {med['port']['wall']:.2f} s | {med['port']['cpu']:.2f} s | {med['port']['rss']:.1f} MiB "
                    f"| {med['port']['cpu'] / med['java']['cpu']:.2f} | {med['port']['rss'] / med['java']['rss']:.3f} | {'equal' if equal else 'DIFFER'} |")
    table = (f"**{result['cpu']}** ({result['machine']}, {result['cores']} cores; {' '.join(result['flags'])}), median of {o.reps}\n\n"
             "| input | Java wall | Java CPU | Java RSS | port wall | port CPU | port RSS | CPU ratio | RSS ratio | outputs |\n"
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n" + "\n".join(rows) + "\n")
    print(table); print(json.dumps(result))
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            f.write(table)
    (work / "result.json").write_text(json.dumps(result, indent=1))
    sys.exit(0 if all(v["equal"] for v in result["inputs"].values()) else 1)


if __name__ == "__main__":
    main()
