"""Run FilterVcf's JAVASCRIPT_FILE in Picard's Nashorn and in the port, on the same scripts and VCFs.

Every script runs on every VCF, both sides in the pinned image at the same paths, and a pair agrees
when the exit status, the first error line, what the run printed and the output VCF are the same.
A port refusal (`gatk_rs::PortLimitation`) is counted apart from a disagreement: it is the port
declining a script its stand-ins do not cover, not a different answer.

    python3 tools/coverage/compare_javascript_filter.py --port target-linux/release/gatk \\
        --vcf <a.vcf> ... --script <a.js> ...

docs/filtervcf-javascript.md records which scripts were run and what came back.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_array  # noqa: E402


def run_side(work, command):
    """Every case in one container, so a JVM start is the only cost per case."""
    cases = sorted(p.name for p in (work / "cases").iterdir())
    out = work / ("port" if "port-binary" in command else "reference")
    out.mkdir()
    (out / "run.sh").write_text("".join(
        f"mkdir -p /work/out/{case} && cd /work/out/{case} && {command} "
        f"$(cat /work/cases/{case}) >stdout 2>stderr; echo $? >status\n" for case in cases
    ) + "chmod -R a+rwX /work/out\n")
    port = ["-v", f"{work / 'bin'}:/work/port-binary:ro"] if (work / "bin").exists() else []
    subprocess.run(
        ["docker", "run", "--rm", "--platform", run_array.PLATFORM,
         "-v", f"{work / 'in'}:/work/in:ro", "-v", f"{work / 'cases'}:/work/cases:ro",
         "-v", f"{out}:/work/out", *port, "-w", "/work", run_array.IMAGE, "sh /work/out/run.sh"],
        check=True, capture_output=True,
    )
    return out


def answer(directory):
    status = (directory / "status").read_text().strip()
    error = run_array.first_error((directory / "stderr").read_text())
    vcf = directory / "output.vcf"
    return status, error if status != "0" else "", (directory / "stdout").read_text(), (
        vcf.read_text() if vcf.exists() else "")


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", required=True, help="a linux/amd64 port binary")
    parser.add_argument("--vcf", nargs="+", required=True)
    parser.add_argument("--script", nargs="+", required=True)
    parser.add_argument("--keep", help="an empty directory to work in and leave the runs in")
    options = parser.parse_args(argv)
    with tempfile.TemporaryDirectory(prefix="gatk-js-") as tmp:
        work = Path(options.keep or tmp)
        (work / "in").mkdir()
        (work / "cases").mkdir()
        for path in options.vcf + options.script:
            shutil.copy(path, work / "in" / Path(path).name)
        for vcf in options.vcf:
            for script in options.script:
                case = f"{Path(script).stem}__{Path(vcf).stem}"
                (work / "cases" / case).write_text(
                    f"FilterVcf --INPUT /work/in/{Path(vcf).name} --OUTPUT output.vcf "
                    f"--CREATE_INDEX false --JAVASCRIPT_FILE /work/in/{Path(script).name}")
        java = f'java {run_array.JAVA_OPENS} -cp "$ORACLE_CP" org.broadinstitute.hellbender.Main'
        reference = run_side(work, java)
        (work / "bin").mkdir()
        shutil.copy(options.port, work / "bin" / "gatk")
        port = run_side(work, "/work/port-binary/gatk")
        tally = {"agree": 0, "refused": 0, "differ": 0}
        for case in sorted(p.name for p in (work / "cases").iterdir()):
            ours, theirs = answer(port / case), answer(reference / case)
            verdict = ("agree" if ours == theirs else
                       "refused" if "gatk_rs::PortLimitation" in ours[1] else "differ")
            tally[verdict] += 1
            print(f"{verdict:8} {case}  reference={theirs[0]} {theirs[1][:90]}  port={ours[0]} {ours[1][:90]}")
        print(" ".join(f"{k}={v}" for k, v in tally.items()))
        return 1 if tally["differ"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
