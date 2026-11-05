#!/usr/bin/env python3
"""The reference's NVScoreVariants over the fixture, as text, inside `broadinstitute/gatk:4.6.2.0`.

For each VCF of `make_fixture.py` and each model: the exit code of the process, what the tool
returned, the score file and the output VCF, or that there is none. Each line says what it belongs
to, so `tools/conformance/compare.py` can tell a score from a header. The errors go to stderr, for
whoever is reading the log: their wording is Python's and no part of the comparison.

    NVScoreVariantsDump.py [fixture directory]

With no directory the fixture is made in a temporary one. `crates/gatk-tools/tests/
nv_score_variants_output.rs` writes the same text for the port, and reads the cases from here.
"""

import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

CASES = [
    "variants",
    "symbolic",
    "empty",
    "gatk_header",
    "breakend",
    "end_tag",
    "undefined_info",
    "undefined_filter",
    "undefined_contig",
    "undefined_format",
    "not_a_number",
    "window_too_long",
    "past_the_end",
    "short_samples",
    "float_annotation_tuple",
    "multiple_values",
]


def lines_of(path):
    if not os.path.exists(path):
        return None
    with open(path) as text:
        return text.read().split("\n")[:-1]


def show(kind, name, tag, lines):
    if lines is None:
        print(f"{name}\t{tag}\t{kind}\tabsent")
        return
    for line in lines:
        print(f"{name}\t{tag}\t{kind}\t{line}")


def run(fixture, work, name, read_tensor):
    tag = "2d" if read_tensor else "1d"
    output = f"{work}/{name}.{tag}.vcf"
    scores = f"{work}/{name}.{tag}.tmp"
    command = [
        "gatk", "NVScoreVariants",
        "-V", f"{fixture}/{name}.vcf",
        "-R", f"{fixture}/reference.fasta",
        "-O", output,
        "--tmp-file", scores,
    ]
    if read_tensor:
        command += ["--tensor-type", "read_tensor", "-I", f"{fixture}/reads.bam"]
    done = subprocess.run(command, capture_output=True, text=True, timeout=1800)
    returned = re.search(r"Tool returned:\s*\n(\S+)", done.stdout + done.stderr)
    print(f"{name}\t{tag}\texit\t{done.returncode}")
    print(f"{name}\t{tag}\treturned\t{returned.group(1) if returned else '-'}")
    show("scores", name, tag, lines_of(scores))
    show("vcf", name, tag, lines_of(output))
    sys.stdout.flush()
    for line in done.stderr.splitlines():
        if re.search(r"Error|\[E::|Exception", line):
            print(f"{name} {tag}: {line}", file=sys.stderr)


def main():
    with tempfile.TemporaryDirectory() as work:
        fixture = sys.argv[1] if len(sys.argv) > 1 else f"{work}/fixture"
        if len(sys.argv) <= 1:
            subprocess.run([sys.executable, f"{HERE}/make_fixture.py", fixture], check=True)
        for name in CASES:
            run(fixture, work, name, False)
            run(fixture, work, name, True)


if __name__ == "__main__":
    main()
