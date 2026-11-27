#!/usr/bin/env bash
# The release binary reads a BAM, writes one through the BGZF deflater, and reads it back.
# `RUN` prefixes every call (e.g. "qemu-x86_64 -cpu qemu64"), so the same check proves the binary
# needs nothing a baseline x86-64 lacks: an unsupported instruction is SIGILL, not a wrong count.
set -euo pipefail
gatk=$1
here=$(cd "$(dirname "$0")" && pwd)
out=$(mktemp -d)
run=${RUN:-}

$run "$gatk" --version
count() { $run "$gatk" CountReads -I "$1" 2>/dev/null | tail -n 1; }
test "$(count "$here/smoke.bam")" = 2
$run "$gatk" PrintReads -I "$here/smoke.bam" -O "$out/out.bam" >/dev/null 2>&1
test "$(count "$out/out.bam")" = 2
echo "smoke: ok${run:+ under $run}"
