#!/usr/bin/env python3
"""Check that every ported symbol comes from a licence-compatible source.

Decision 0014: `crates/jmath/src/exp.rs` was a line-by-line transcription of
`macroAssembler_x86_exp.cpp`, which is GPL2-only, published in an MIT crate. It was found by
reading the module headers eleven commits after it was merged. This makes that reading
mechanical, so the next one fails CI instead.

The rule: a symbol is portable into this program only if it comes from the pinned htsjdk,
Picard or GATK clones, all of which are MIT or Apache 2.0. Anything reached through the JDK is
GPL2 and is not portable, whatever the Classpath Exception says about linking.
"""

import re
import sys
from pathlib import Path

# Package prefixes that are safe to port from, with the licence that makes them safe.
ALLOWED = {
    # Fully-qualified class names, and the bare repository names used in module headers that
    # cite a whole file rather than one symbol.
    "htsjdk": "htsjdk, MIT",
    "picard": "Picard, MIT",
    "org.broadinstitute.hellbender": "GATK, Apache 2.0",
}

# Anything matching these is GPL2 and must not be transcribed. Listed explicitly rather than
# as "not allowed", so a new source has to be classified deliberately.
FORBIDDEN = {
    "openjdk": "OpenJDK, GPL2",
    "jdk/src": "OpenJDK, GPL2",
    "hotspot": "HotSpot, GPL2 only, no Classpath Exception",
    "java.lang.": "JDK class library, GPL2",
    "java.util.": "JDK class library, GPL2",
    "java.text.": "JDK class library, GPL2",
    "java.math.": "JDK class library, GPL2",
}

CLAIM = re.compile(r"Ported from\s+`?([^\s`,]+)")


def main(roots) -> int:
    violations, checked = [], 0
    for root in roots:
        for path in sorted(Path(root).rglob("*.rs")):
            text = path.read_text(errors="replace")
            for line_no, line in enumerate(text.splitlines(), 1):
                m = CLAIM.search(line)
                if not m:
                    continue
                source = m.group(1)
                checked += 1
                lowered = source.lower()
                bad = next((w for w in FORBIDDEN if w in lowered), None)
                if bad:
                    violations.append(
                        (path, line_no, source, FORBIDDEN[bad])
                    )
                    continue
                # Case-insensitive: module headers cite "Picard 3.4.0" in prose and
                # "picard.analysis.X" as a symbol, and both are the same permissive source.
                if not any(lowered.startswith(a.lower()) for a in ALLOWED):
                    violations.append(
                        (path, line_no, source, "unclassified source; add it to ALLOWED "
                         "or FORBIDDEN in this script after checking its licence")
                    )

    print(f"checked {checked} `Ported from` claims")
    for path, line_no, source, why in violations:
        print(f"  {path}:{line_no}: {source}\n      -> {why}")
    if violations:
        print(f"\n{len(violations)} provenance violations. See decision 0014.")
        return 1
    print("all claims resolve to a licence-compatible source")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["crates"]))
