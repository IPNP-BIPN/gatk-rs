# Oracle image for gatk-rs golden generation.
#
# Same contract as picard-rs and htsjdk-rs: the base image pinned by digest, the release pinned by
# sha256, the platform linux/amd64, and a probe that runs during the build so a degraded
# environment cannot produce an image, let alone a golden.
#
# The GATK release zip is taken rather than a hand-assembled classpath. It carries the local jar
# with Picard 3.4.0, htsjdk 4.2.0, Barclay 5.0.0 and Intel GKL 0.8.11 exactly as build.gradle at
# tag 4.6.2.0 pins them, and it carries the `gatk` wrapper script, which is the entry point the
# bit-identity claim is defined against: `gatk <Tool>` fixes the argument parser to Barclay
# (`--INPUT`) rather than Picard's legacy syntax.
#
# The version is an argument so a second oracle can be built beside the first rather than in place
# of it. The move to 4.7.0.0 is not a matter of editing this file: it is a matter of building both
# images and re-running every suite against the new one, so the two have to coexist.
#
#     docker build --platform linux/amd64 -t gatk-rs-oracle:4.6.2.0 tools/oracle
#     docker build --platform linux/amd64 -t gatk-rs-oracle:4.7.0.0 \
#       --build-arg GATK_VERSION=4.7.0.0 \
#       --build-arg GATK_SHA256=d093d2693b1626361a413ca59d6d4a0bf968717f280a8fd9ce060b25eb2ed1db \
#       tools/oracle
#
# The checksum is an argument too, defaulting to 4.6.2.0's, so naming a new version without also
# naming its checksum fails at `sha256sum -c` rather than producing an unverified image.

FROM eclipse-temurin@sha256:068a8f9ae4b74d9a20de3ecca771ba1a6437f1d7a8a8ad6deaf9dbdd2274397a

ARG GATK_VERSION=4.6.2.0
ARG GATK_SHA256=32d2f90bf13fcb3a8ac765bb2cb8ec1fc9a6cc447055d0156bd1db2092d4e3e8

# python3 is not optional: the `gatk` wrapper is a Python script, so without it the entry point the
# claim is defined against cannot run at all.
#
# R is deliberately *not* installed yet. picard-rs needs it because ten of its metrics tools
# declare a chart argument as REQUIRED and shell out to Rscript; the GATK tools that plot are
# fewer, and none is on the critical path. It is a build-time cost of several minutes, so it goes
# in when the first tool needs it rather than on the chance that one might.
RUN for attempt in 1 2 3; do \
      apt-get update \
      && apt-get install -y --no-install-recommends curl ca-certificates unzip python3 \
      && break; \
      echo "apt attempt $attempt failed; retrying"; sleep 10; \
    done \
 && command -v python3 >/dev/null \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/oracle

RUN set -eux; \
    curl -fsSL -o /tmp/gatk.zip \
      https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip; \
    echo "${GATK_SHA256}  /tmp/gatk.zip" \
      | sha256sum -c -; \
    unzip -q /tmp/gatk.zip -d /opt; \
    rm /tmp/gatk.zip; \
    ln -s /opt/gatk-${GATK_VERSION} /opt/gatk

ENV GATK_HOME=/opt/gatk
ENV GATK_JAR=/opt/gatk/gatk-package-${GATK_VERSION}-local.jar
ENV ORACLE_CP=/opt/gatk/gatk-package-${GATK_VERSION}-local.jar:/opt/oracle/classes
ENV PATH=/opt/gatk:$PATH

# What the probe asserts GATK reports about itself. Passed in rather than compiled in, so one
# probe serves both images and neither can be built against the wrong jar.
ENV ORACLE_GATK_VERSION=${GATK_VERSION}

# JaCoCo for the coverage-guided fuzzer, pinned by sha256 like everything else: the fuzzer's
# stopping condition is a coverage threshold on the reference, so a silently different JaCoCo would
# silently change what "covered" means.
RUN set -eux; \
    mkdir -p /opt/jacoco/lib; cd /opt/jacoco/lib; \
    for spec in \
      "org/jacoco/org.jacoco.core/0.8.12/org.jacoco.core-0.8.12.jar fca26db37c0c5fbd5dc4985237eb82866df9799d5082af899475a73f91f5b035" \
      "org/ow2/asm/asm/9.7/asm-9.7.jar adf46d5e34940bdf148ecdd26a9ee8eea94496a72034ff7141066b3eea5c4e9d" \
      "org/ow2/asm/asm-commons/9.7/asm-commons-9.7.jar 389bc247958e049fc9a0408d398c92c6d370c18035120395d4cba1d9d9304b7a" \
      "org/ow2/asm/asm-tree/9.7/asm-tree-9.7.jar 62f4b3bc436045c1acb5c3ba2d8ec556ec3369093d7f5d06c747eb04b56d52b1" \
    ; do \
      set -- $spec; \
      curl -fsSL -O "https://repo1.maven.org/maven2/$1"; \
      echo "$2  $(basename $1)" | sha256sum -c -; \
    done; \
    curl -fsSL -o /tmp/jacoco.zip \
      https://repo1.maven.org/maven2/org/jacoco/jacoco/0.8.12/jacoco-0.8.12.zip; \
    echo "35cf2a2b8937659ecbc8d0d3902dfa7f55b7ed2250e82424036a3e1d12462cd7  /tmp/jacoco.zip" \
      | sha256sum -c -; \
    cd /opt/jacoco && jar xf /tmp/jacoco.zip lib/jacocoagent.jar lib/jacococli.jar && rm /tmp/jacoco.zip

ENV JACOCO_AGENT=/opt/jacoco/lib/jacocoagent.jar
ENV JACOCO_ANALYSIS=/opt/jacoco/lib/org.jacoco.core-0.8.12.jar:/opt/jacoco/lib/asm-9.7.jar:/opt/jacoco/lib/asm-commons-9.7.jar:/opt/jacoco/lib/asm-tree-9.7.jar

WORKDIR /opt/oracle
COPY OracleProbe.java .
RUN mkdir -p classes && javac -cp "$ORACLE_CP" -d classes OracleProbe.java

# Fail the build itself if the environment does not satisfy the contract.
#
# The deflater property is passed here exactly as the manifest passes it to every suite, because
# the probe's job is to assert the configuration the goldens are produced under, not some other
# one. It is also the only way to find out whether the property still works: under htsjdk 5.0.0
# this run is what proves `-Dsamjdk.use_libdeflate=false` actually reaches
# `Defaults.USE_LIBDEFLATE`, and under 4.2.0 it proves the property is harmlessly unknown. The
# first attempt at building a 4.7.0.0 image failed right here, which is the check working.
RUN java -Dsamjdk.use_libdeflate=false -cp "$ORACLE_CP:classes" OracleProbe

WORKDIR /work
ENTRYPOINT ["/bin/sh", "-c"]
CMD ["java -Dsamjdk.use_libdeflate=false -cp \"$ORACLE_CP:/opt/oracle/classes\" OracleProbe"]
