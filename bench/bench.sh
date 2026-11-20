#!/bin/bash
# bench.sh <old-binary> <new-binary> <cpus|all> <out.json>: GATK, the old port and the new one over
# three cases against one htsget-rs server; timer.py inside the oracle image, sides alternating.
set -euo pipefail
OLD=$1; NEW=$2; CPUS=$3; OUT=$4; W=$PWD/bench-work
LIMIT=(); [ "$CPUS" != all ] && LIMIT=(--cpus "$CPUS")
P="-Duser.language=en -Duser.country=US -Dsamjdk.use_libdeflate=false"
TRUST="-Djavax.net.ssl.trustStore=/work/htsget/trust.p12 -Djavax.net.ssl.trustStorePassword=changeit -Djavax.net.ssl.trustStoreType=PKCS12"
mkdir -p $W/bin && cp "$OLD" $W/bin/old && cp "$NEW" $W/bin/new
echo "[" > "$OUT"; first=1
for c in whole interval interval-t4; do
  a="--url https://htsget:8080/reads --id big"
  [ $c != whole ] && a="$a -L chrS:20000000-30000000"
  [ $c = interval-t4 ] && a="$a --reader-threads 4"
  for pair in "java:old" "java:new" "old:new"; do
    one=${pair%%:*}; two=${pair##*:}
    cmd() { case $1 in java) echo "java $P $TRUST -cp \"\$ORACLE_CP\" org.broadinstitute.hellbender.Main HtsgetReader $a -O /work/out/o.bam";;
            *) echo "env SSL_CERT_FILE=/work/htsget/cert.pem /work/bin/$1 HtsgetReader $a -O /work/out/o.bam";; esac; }
    python3 -c 'import json,sys; json.dump({"reps":5,"timeout":1800,"oracle":sys.argv[1],"port":sys.argv[2]},open(sys.argv[3],"w"))' "$(cmd $one)" "$(cmd $two)" $W/spec.json
    docker run --rm "${LIMIT[@]}" --network htsget-bench -v $W/htsget:/work/htsget:ro -v $W/bin:/work/bin:ro -v $W/spec.json:/work/spec.json:ro \
      -v $PWD/tools/speed/timer.py:/work/timer.py:ro gatk-rs-oracle:4.7.0.0 "python3 /work/timer.py /work/spec.json" > $W/r.json
    [ $first = 1 ] || echo "," >> "$OUT"; first=0
    python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); r.update(case=sys.argv[2],a=sys.argv[3],b=sys.argv[4]); print(json.dumps(r))' $W/r.json $c $one $two >> "$OUT"
  done
  docker run --rm "${LIMIT[@]}" --network htsget-bench -v $W/htsget:/work/htsget:ro -v $W/bin:/work/bin:ro gatk-rs-oracle:4.7.0.0 \
    "cd /tmp && java $P $TRUST -cp \"\$ORACLE_CP\" org.broadinstitute.hellbender.Main HtsgetReader $a -O java.bam >/dev/null 2>&1; for b in old new; do SSL_CERT_FILE=/work/htsget/cert.pem /work/bin/\$b HtsgetReader $a -O \$b.bam >/dev/null 2>&1; done; sha256sum java.bam old.bam new.bam" | sed "s/^/$c $CPUS /" | tee -a $W/sha256.txt
done
echo "]" >> "$OUT"
