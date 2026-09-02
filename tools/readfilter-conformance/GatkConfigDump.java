/*
 * `GATKConfig`'s system properties, which `Main` installs before any tool runs.
 *
 * The reason this is measured at all: a covering-array row over `IndexFeatureFile` disagreed with
 * the port on a `.tbi` whose INDEX was already byte for byte, and the difference was the
 * compression. `samjdk.compression_level` is declared here with a default of TWO, where htsjdk's
 * own default is five, and two is the one level pair Intel's GKL routes through igzip rather than
 * zlib. So every BGZF byte a real `gatk` invocation writes depends on this table (gatk-rs#1032).
 *
 * Four behaviours this is built to catch.
 *
 *   - THE DEFAULTS ARE ANNOTATIONS, not a file: `@DefaultValue` on each method of the interface,
 *     and a config file only overrides them;
 *   - WHICH KEYS REACH System.getProperties is `@SystemProperty`'s answer, and here it is ALL of
 *     them: twelve keys, twelve properties. The distinction is measured rather than assumed,
 *     because a key added without the annotation would be read by GATK and never seen by htsjdk;
 *   - THE PROPERTY NAME IS `@Key`, not the method name, and the two differ by more than case;
 *   - AND A PROPERTY ALREADY SET IS NOT OVERWRITTEN, so a `-D` on the command line wins over the
 *     config's default.
 *
 * Output:
 *
 *     key\t<property>\t<default>\t<system|internal>
 *     injected\t<property>\t<the value System.getProperties holds after injection>
 *     precedence\t<property>\t<the value after injection when -D set it first>
 *
 * Usage: GatkConfigDump
 */

import org.aeonbits.owner.Config;
import org.broadinstitute.hellbender.utils.config.ConfigFactory;
import org.broadinstitute.hellbender.utils.config.GATKConfig;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

public class GatkConfigDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t').append(payload).append('\n');
    }

    /** The interface's own declaration order, which is the order the table is read in. */
    static List<Method> declared() {
        final List<Method> methods = new ArrayList<>();
        for (final Method method : GATKConfig.class.getDeclaredMethods()) {
            if (method.getAnnotation(Config.Key.class) != null) {
                methods.add(method);
            }
        }
        // getDeclaredMethods has no defined order, so the table is sorted by the property name it
        // will be printed under. A golden cannot depend on a JVM's reflection order.
        methods.sort((a, b) -> a.getAnnotation(Config.Key.class).value()
                .compareTo(b.getAnnotation(Config.Key.class).value()));
        return methods;
    }

    public static void main(final String[] args) {
        final List<Method> methods = declared();

        // The declaration: the property's name, its default, and whether it is injected at all.
        for (final Method method : methods) {
            final String key = method.getAnnotation(Config.Key.class).value();
            final Config.DefaultValue value = method.getAnnotation(Config.DefaultValue.class);
            final boolean system =
                    method.getAnnotation(Config.Sources.class) == null
                            && hasSystemProperty(method);
            emit("key", key, (value == null ? "<none>" : value.value()) + "\t"
                    + (system ? "system" : "internal"));
        }

        // What `Main` leaves behind: the config is built from no command line at all, which is the
        // default path, and then injected.
        ConfigFactory.getInstance().initializeConfigurationsFromCommandLineArgs(
                new String[0], "--gatk-config-file");
        final GATKConfig config = ConfigFactory.getInstance().getGATKConfig();
        ConfigFactory.getInstance().injectSystemPropertiesFromConfig(config);
        for (final Method method : methods) {
            final String key = method.getAnnotation(Config.Key.class).value();
            if (!hasSystemProperty(method)) {
                continue;
            }
            final String value = System.getProperty(key);
            emit("injected", key, value == null ? "<unset>" : value);
        }

        // And whether the injection overwrites a property that is already set, which decides
        // whether `-Dsamjdk.compression_level=5` on the command line means anything.
        System.setProperty("samjdk.compression_level", "9");
        ConfigFactory.getInstance().injectSystemPropertiesFromConfig(config);
        emit("precedence", "samjdk.compression_level",
                System.getProperty("samjdk.compression_level"));

        System.out.print(buf);
    }

    /** `@SystemProperty` lives in GATK's own package, so it is looked up by name. */
    static boolean hasSystemProperty(final Method method) {
        for (final java.lang.annotation.Annotation annotation : method.getAnnotations()) {
            if (annotation.annotationType().getSimpleName().equals("SystemProperty")) {
                return true;
            }
        }
        return false;
    }
}
