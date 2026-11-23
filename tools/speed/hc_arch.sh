#!/usr/bin/env bash
# HaplotypeCaller, GATK 4.7 against the port, on this machine: wall, user, sys, peak RSS, and the
# md5 of each output without its command line. Usage: hc_arch.sh DATA JAR PORT OUT
set -uo pipefail
DATA=$1 JAR=$2 PORT=$3 OUT=$4
mkdir -p "$OUT"
NPROC=$(nproc)
{
  echo "## $(uname -m), $NPROC CPUs"
  echo
  echo '```'
  lscpu | grep -E 'Model name|Vendor|^CPU\(s\)|Thread|Core|Socket' || true
  echo "flags: $(grep -o -w -E 'avx|avx2|fma|avx512f|avx512bw|avx512vl|sse4_2|asimd|sve|sve2' /proc/cpuinfo | sort -u | tr '\n' ' ')"
  free -g | head -2
  java -version 2>&1 | head -1
  echo '```'
  echo
  echo "| interval | mode | who | threads | wall s | user s | sys s | peak RSS MB | md5 (no command line) |"
  echo "|---|---|---|---|---|---|---|---|---|"
} > "$OUT/table.md"

row() { # interval mode who threads log out
  local wall user sys rss sum
  wall=$(grep 'Elapsed (wall clock)' "$5" | awk '{print $NF}' | awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2 }')
  user=$(grep 'User time' "$5" | awk '{print $NF}')
  sys=$(grep 'System time' "$5" | awk '{print $NF}')
  rss=$(grep 'Maximum resident' "$5" | awk '{printf "%.0f", $NF/1024}')
  if [[ -s "$6" ]]; then sum=$(grep -v '^##GATKCommandLine' "$6" | md5sum | cut -c1-12); else sum="(failed: $(tail -1 "$5" | cut -c1-60))"; fi
  echo "| $1 | $2 | $3 | $4 | $wall | $user | $sys | $rss | $sum |" | tee -a "$OUT/table.md"
}

for L in chr20:10000000-11000000 chr20:10000000-20000000; do
  for MODE in NONE GVCF; do
    tag="${L#chr20:}.$MODE"
    common=(-I "$DATA/NA12878.chr20.bam" -R "$DATA/chr20.fasta" -L "$L" -ERC "$MODE" --add-output-vcf-command-line false)
    /usr/bin/time -v java -jar "$JAR" HaplotypeCaller "${common[@]}" -O "$OUT/java.$tag.vcf" > "$OUT/java.$tag.log" 2>&1
    row "$L" "$MODE" java default "$OUT/java.$tag.log" "$OUT/java.$tag.vcf"
    for T in $(printf "%s\n" 1 2 4 "$NPROC" | sort -nu); do
      RAYON_NUM_THREADS=$T /usr/bin/time -v "$PORT" HaplotypeCaller "${common[@]}" -O "$OUT/port.$tag.$T.vcf" > "$OUT/port.$tag.$T.log" 2>&1
      row "$L" "$MODE" port "$T" "$OUT/port.$tag.$T.log" "$OUT/port.$tag.$T.vcf"
    done
  done
done
