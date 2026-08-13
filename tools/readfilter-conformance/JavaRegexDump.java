/*
 * java.util.regex, as GATK asks it: compile a user string and call find().
 *
 * `Utils.filterCollectionByExpressions` is what `-se` and `-xl-se` go through and what
 * `VariantEval` uses for its sample expressions; `ClipReads` compiles its `-X` sequences the same
 * way. All of them SEARCH rather than match, all of them let a compile failure reach the user
 * unwrapped, and none of them anchors anything.
 *
 * Eight behaviours this is built to catch, four of which are why a crate cannot stand in.
 *
 *   - find() IS A SEARCH, so `s1` selects `xs10`. Anchoring is the caller's business and no caller
 *     does it;
 *   - `$` MATCHES BEFORE A FINAL LINE TERMINATOR, not only at the end, and that includes a
 *     trailing "\r\n" taken as one;
 *   - `.` REFUSES FOUR LINE TERMINATORS, not one: \n, \r, ,   and  ;
 *   - THE PREDEFINED CLASSES ARE ASCII. Without UNICODE_CHARACTER_CLASS, \d is [0-9] and \w is
 *     [a-zA-Z0-9_], so an Arabic-Indic digit is not a digit and an accented letter is not a word
 *     character;
 *   - POSSESSIVE QUANTIFIERS TAKE AND NEVER GIVE BACK, so `^.*+1$` fails where `^.*1$` succeeds;
 *   - MATCHING IS BACKTRACKING AND ALTERNATION PREFERS THE FIRST BRANCH, which is observable
 *     through the leftmost result;
 *   - A CLASS WHOSE FIRST CHARACTER IS `]` CONTAINS IT rather than being empty;
 *   - AND A COMPILE FAILURE CARRIES A DESCRIPTION, AN INDEX AND A CARET. The index is not always
 *     where the reader would put it: `a)` reports index 0 and a backslash at the end reports 1, so
 *     the invalid patterns here vary the prefix to pin the rule rather than one example of it.
 *
 * Output:
 *
 *     compile\t<pattern>\tok
 *     compile\t<pattern>\terror\t<description>\t<index>\t<getMessage>
 *     find\t<pattern>\t<input>\t<true|false>
 *
 * Patterns and inputs are escaped: a backslash doubles, a tab, newline and carriage return become
 * \t, \n and \r, and anything outside printable ASCII becomes a backslash-u escape per UTF-16
 * code unit.
 *
 * Usage: JavaRegexDump
 */

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

public class JavaRegexDump {

    /** What a sample expression or a clipped sequence can be, plus what it must not silently be. */
    static final String[] PATTERNS = {
        // The shapes a sample expression actually takes.
        "s1", "^s1$", "^NA", "NA12878", "^(NA|s)1$", "[a-c]", "[^a-c]", "[]]", "[\\d_]",
        // The metacharacters, one at a time.
        ".", "a.c", "\\d", "\\D", "\\w", "\\W", "\\s", "\\S", "$", "^", "a|", "",
        // Quantifiers, in all three appetites.
        "^a*b$", "^a+b$", "^ab?$", "^a{2,3}$", "^a{2}$", "^a{2,}$", "^.*1$", "^.*+1$", "^.*?1$",
        "a??b", "(ab)+c", "^(a*)*$",
        // Constructs Java has and a subset need not: each is measured so the boundary is a fact.
        "\\bNA", "(?i)na", "(?:ab)+", "(?=NA)N", "(?!NA)s", "(\\w)\\1", "(?>a*)b",
        "[a-z&&[^aeiou]]", "\\p{Alpha}+", "\\p{Lower}", "\\Qa.c\\E", "\\x41", "\\070", "\\cA",
        "^\\p{IsLatin}+$", "\\ANA", "NA\\z", "NA\\Z", "(?s).", "(?m)^x$", "\\h", "\\R",
        "(?<name>NA)\\k<name>", "[^\\W\\d]", "[[a-c]&&[b-d]]",
        // The invalid ones, with the prefix varied so the index rule is visible and not guessed.
        "[", "a[", "abc[", "(", "a(", "abc(", ")", "a)", "abc)", "(a))",
        "*a", "a**", "{2}", "a{2,1}", "ab{3,2}", "\\", "a\\", "[z-a]", "a[z-a]", "[a-\\", "(?",
        // A brace with nothing before it, and the shapes a brace can take that are not a count.
        "{2,3}", "x{2}y", "{a}", "{2", "a{2", "{}", "a{,2}", "|{2}",
        // Ranges whose ends are not both plain characters.
        "[a-]", "[-a]", "[a-\\d]", "[\\d-a]", "[a-b-c]",
    };

    /** Sample names, a few that differ only above the BMP, and the line-terminator cases. */
    static final String[] INPUTS = {
        "s1", "xs10", "NA12878", "tumor", "a.c", "abc", "aaab", "ab", "b", "aa", "aaa", "aaaa",
        "]", "_", "A", "٣", "é", "café", "NANA", "nA", "a*c",
        "s1\n", "s1\r\n", "s1\nx", "s1\r", "a\rc", "ac", "ac", "a c", "", "x",
    };

    public static void main(final String[] args) {
        System.out.println("# JavaRegexDump: compile and find(), from the reference");
        for (final String pattern : PATTERNS) {
            final Pattern compiled;
            try {
                compiled = Pattern.compile(pattern);
            } catch (final PatternSyntaxException e) {
                System.out.printf("compile\t%s\terror\t%s\t%d\t%s%n", escape(pattern),
                        escape(e.getDescription()), e.getIndex(), escape(e.getMessage()));
                continue;
            }
            System.out.printf("compile\t%s\tok%n", escape(pattern));
            for (final String input : INPUTS) {
                final Matcher matcher = compiled.matcher(input);
                System.out.printf("find\t%s\t%s\t%s%n", escape(pattern), escape(input),
                        matcher.find() ? "true" : "false");
            }
        }
    }

    static String escape(final String text) {
        final StringBuilder out = new StringBuilder();
        for (final char c : text.toCharArray()) {
            switch (c) {
                case '\\' -> out.append("\\\\");
                case '\t' -> out.append("\\t");
                case '\n' -> out.append("\\n");
                case '\r' -> out.append("\\r");
                default -> {
                    if (c < 0x20 || c > 0x7e) {
                        out.append(String.format("\\u%04x", (int) c));
                    } else {
                        out.append(c);
                    }
                }
            }
        }
        return out.toString();
    }
}
