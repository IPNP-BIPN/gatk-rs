# fermilite-conformance

Compares `crates/gatk-fermilite` with the fermi-lite r53 that GATK 4.6.2.0 bundles. Nothing here is a
golden: the corpus is regenerated from a seed, and the reference is the JNI in the oracle image.

    python3 gen_corpus.py main corpus.bin
    python3 gen_corpus.py edge edge.bin

    FERMILITE_CORPUS=corpus.bin FERMILITE_DUMP=port.txt \
        cargo run --release -p gatk-fermilite --example fermilite_dump

    # in gatk-rs-oracle:4.6.2.0
    javac -cp "$GATK_JAR" -d . FermiLiteDump.java
    java -cp ".:$GATK_JAR" FermiLiteDump corpus.bin jni.txt        # add a thread count to call setNThreads

    python3 compare.py jni.txt port.txt

`FermiLiteDump` restarts a worker JVM when an `assert` kills one and writes an `X` line for that
assembly; the worker's messages are in `jni.txt.stderr`. `J` is the `BufferOverflowException` of
`makeReadData`. `U` is a refusal by the port, and `compare.py` counts those without comparing them.

`native_driver.c` is the fast loop, a clang build of the r53 sources with the JNI's read loop and a
fork per assembly. It disagrees with the JNI by construction in two places: it has no Java buffer
(`J` records of the main corpus), and `.data` is laid out differently (`edge` corpus). The JNI wins.
