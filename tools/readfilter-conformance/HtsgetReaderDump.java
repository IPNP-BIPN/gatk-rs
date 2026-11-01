/*
 * HtsgetReader through Main, against two servers this dump starts on fixed ports: TLS on 18443,
 * because the tool asks for https whatever --url says, and plain on 18080 for a block's own URL.
 *
 * The routes are printed as `route` rows so the port's test can serve the same answers. A request
 * no server can answer goes to port 1, where the refusal names the URI the builder made: that is
 * how the rows over --id, --url, --tag and the rest measure the URI without a server.
 */

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpServer;
import com.sun.net.httpserver.HttpsConfigurator;
import com.sun.net.httpserver.HttpsServer;
import org.broadinstitute.hellbender.Main;

import javax.net.ssl.KeyManagerFactory;
import javax.net.ssl.SSLContext;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.KeyStore;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

public class HtsgetReaderDump {
    static final Path DIR = Path.of("htsget-dump");
    static final String TLS = "https://localhost:18443";
    static final String PLAIN = "http://localhost:18080";
    static final String UNREACHABLE = "https://localhost:1";

    /** Statuses in turn, the last one repeating; where a redirect goes; the body. */
    record Route(String statuses, String location, String body) {
        Route(final String body) {
            this("200", "-", body);
        }
    }

    static final Map<String, Route> ROUTES = new LinkedHashMap<>();
    static final Map<String, AtomicInteger> SEEN = new java.util.concurrent.ConcurrentHashMap<>();

    static String dataUri(final String text) {
        return "data:application/octet-stream;base64,"
                + Base64.getEncoder().encodeToString(text.getBytes(StandardCharsets.UTF_8));
    }

    static String answer(final String urls) {
        return "{\"htsget\":{\"format\":\"BAM\",\"urls\":[" + urls + "]}}";
    }

    static String block(final String url) {
        return "{\"url\":\"" + url + "\"}";
    }

    static void routes() {
        final String one = block(dataUri("one"));
        ROUTES.put("/data", new Route(answer(one)));
        ROUTES.put("/two", new Route(answer(one + "," + block(dataUri("two")))));
        ROUTES.put("/block", new Route("the plain block\n"));
        ROUTES.put("/secure", new Route("the secure block\n"));
        ROUTES.put("/gone", new Route("404", "-", "the body of a block that was not there\n"));
        ROUTES.put("/echo", new Route("echoes the header"));
        final String headers = ",\"headers\":{}";
        ROUTES.put("/mixed", new Route(answer(one + ",{\"url\":\"" + PLAIN + "/block\"" + headers + ",\"class\":\"body\"},"
                + "{\"url\":\"" + TLS + "/secure\"" + headers + ",\"class\":\"header\"}")));
        ROUTES.put("/block-404", new Route(answer("{\"url\":\"" + PLAIN + "/gone\"" + headers + "}")));
        ROUTES.put("/block-headers", new Route(answer(
                "{\"url\":\"" + PLAIN + "/echo\",\"headers\":{\"X-Test\":\"sent\",\"X-Number\":7}}")));
        ROUTES.put("/block-no-headers", new Route(answer("{\"url\":\"" + PLAIN + "/block\"}")));
        ROUTES.put("/block-refused", new Route(answer("{\"url\":\"http://localhost:1/x\"" + headers + "}")));
        ROUTES.put("/empty", new Route(answer("")));
        ROUTES.put("/no-urls", new Route("{\"htsget\":{\"format\":\"BAM\"}}"));
        ROUTES.put("/null-urls", new Route("{\"htsget\":{\"urls\":null}}"));
        ROUTES.put("/no-url", new Route(answer("{\"class\":\"body\"}")));
        ROUTES.put("/no-scheme", new Route(answer(block("relative/path"))));
        ROUTES.put("/ftp", new Route(answer(block("ftp://localhost/x"))));
        ROUTES.put("/not-base64", new Route(answer(block("data:text/plain,abc"))));
        ROUTES.put("/bad-base64", new Route(answer(block("data:;base64,Q!=="))));
        ROUTES.put("/bad-ending", new Route(answer(block("data:;base64,QQ=A"))));
        ROUTES.put("/short-base64", new Route(answer(block("data:;base64,Q"))));
        ROUTES.put("/two-markers", new Route(answer(block("data:;base64,x;base64,QUJD"))));
        ROUTES.put("/md5-good", new Route("{\"htsget\":{\"urls\":[" + one + "],\"md5\":\"f97c5d29941bfb1b2fdab0874906ab82\"}}"));
        ROUTES.put("/md5-bad", new Route("{\"htsget\":{\"urls\":[" + one + "],\"md5\":\"F97C5D29941BFB1B2FDAB0874906AB82\"}}"));
        ROUTES.put("/md5-number", new Route("{\"htsget\":{\"urls\":[" + one + "],\"md5\":12}}"));
        ROUTES.put("/upper", new Route("{\"htsget\":{\"URLS\":[{\"URL\":\"" + dataUri("upper") + "\"}],\"FORMAT\":\"VCF\"}}"));
        ROUTES.put("/not-json", new Route("this is not json"));
        ROUTES.put("/empty-body", new Route(""));
        ROUTES.put("/array", new Route("[]"));
        ROUTES.put("/other-root", new Route("{\"other\":{\"urls\":[]}}"));
        ROUTES.put("/extra-root", new Route("{\"htsget\":{\"urls\":[]},\"more\":1}"));
        ROUTES.put("/first-not-root", new Route("{\"more\":1,\"htsget\":{\"urls\":[]}}"));
        ROUTES.put("/upper-root", new Route("{\"HTSGET\":{\"urls\":[]}}"));
        ROUTES.put("/trailing", new Route("{\"htsget\":{\"urls\":[" + one + "]}} and more"));
        ROUTES.put("/unknown-property", new Route("{\"htsget\":{\"urls\":[],\"extra\":1}}"));
        ROUTES.put("/unknown-block-property", new Route(answer("{\"url\":\"data:;base64,QQ==\",\"extra\":1}")));
        ROUTES.put("/bad-format", new Route("{\"htsget\":{\"format\":\"bam\",\"urls\":[]}}"));
        ROUTES.put("/bad-class", new Route(answer("{\"url\":\"data:;base64,QQ==\",\"class\":\"BODY\"}")));
        ROUTES.put("/bad-url", new Route(answer(block("not a uri"))));
        ROUTES.put("/number-url", new Route(answer("{\"url\":7}")));
        ROUTES.put("/urls-object", new Route("{\"htsget\":{\"urls\":{}}}"));
        ROUTES.put("/md5-object", new Route("{\"htsget\":{\"urls\":[],\"md5\":{}}}"));
        ROUTES.put("/error-400", new Route("400", "-", "{\"htsget\":{\"error\":\"InvalidInput\",\"message\":\"bad\"}}"));
        ROUTES.put("/error-404", new Route("404", "-", "{\"htsget\":{\"error\":\"NotFound\"}}"));
        ROUTES.put("/error-499", new Route("499", "-", "{\"htsget\":{}}"));
        ROUTES.put("/error-text", new Route("403", "-", "Forbidden"));
        ROUTES.put("/error-unknown", new Route("400", "-", "{\"htsget\":{\"error\":\"x\",\"code\":1}}"));
        ROUTES.put("/no-content", new Route("204", "-", ""));
        ROUTES.put("/accepted", new Route("202", "-", answer(one)));
        ROUTES.put("/moved", new Route("302", TLS + "/data", ""));
        ROUTES.put("/flaky", new Route("500,503,200", "-", answer(one)));
        ROUTES.put("/broken", new Route("500", "-", "{}"));
    }

    static void serve(final HttpExchange exchange) throws java.io.IOException {
        final String path = exchange.getRequestURI().getPath();
        final Route route = ROUTES.get(path);
        byte[] body = (route == null ? "no route" : route.body()).getBytes(StandardCharsets.UTF_8);
        int status = 404;
        if (route != null) {
            final String[] statuses = route.statuses().split(",");
            final int seen = SEEN.computeIfAbsent(path, key -> new AtomicInteger()).getAndIncrement();
            status = Integer.parseInt(statuses[Math.min(seen, statuses.length - 1)]);
            if (!"-".equals(route.location())) {
                exchange.getResponseHeaders().add("Location", route.location());
            }
            if ("/echo".equals(path)) {
                body = (exchange.getRequestHeaders().getFirst("X-Test") + "/"
                        + exchange.getRequestHeaders().getFirst("X-Number") + "\n").getBytes(StandardCharsets.UTF_8);
            }
        }
        exchange.sendResponseHeaders(status, body.length == 0 ? -1 : body.length);
        if (body.length > 0) {
            exchange.getResponseBody().write(body);
        }
        exchange.close();
    }

    static SSLContext tls(final Path store) throws Exception {
        final char[] password = "changeit".toCharArray();
        final KeyStore keys = KeyStore.getInstance("PKCS12");
        try (var in = Files.newInputStream(store)) {
            keys.load(in, password);
        }
        final KeyManagerFactory factory = KeyManagerFactory.getInstance(KeyManagerFactory.getDefaultAlgorithm());
        factory.init(keys, password);
        final SSLContext context = SSLContext.getInstance("TLS");
        context.init(factory.getKeyManagers(), null, null);
        return context;
    }

    static void run(final String label, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>(List.of("HtsgetReader"));
        argv.addAll(Arrays.asList(rest));
        System.out.printf("argv\t%s\t%s%n", label, ReferenceQueryDump.escape(String.join("\t", argv)));
        final ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
        final PrintStream realOut = System.out;
        String failure = null;
        try {
            System.setOut(new PrintStream(outBytes, true, StandardCharsets.UTF_8));
            try {
                new Main().instanceMain(argv.toArray(new String[0]));
            } catch (final Exception e) {
                failure = e.getClass().getName() + ":" + e.getMessage();
            }
        } finally {
            System.out.flush();
            System.setOut(realOut);
        }
        System.out.printf("stdout\t%s\t%s%n", label, ReferenceQueryDump.escape(outBytes.toString(StandardCharsets.UTF_8)));
        if (failure != null) {
            System.out.printf("error\t%s\t%s%n", label, ReferenceQueryDump.escape(failure));
        }
        // The file is listed either way: a block that fails leaves what came before it.
        final Path written = DIR.resolve("out").resolve(label);
        if (Files.exists(written)) {
            System.out.printf("file\t%s\t%s%n", label, Base64.getEncoder().encodeToString(Files.readAllBytes(written)));
        }
    }

    /** One run named for the route it asks for. */
    static void fetch(final String label, final String route, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>(List.of("--url", TLS, "--id", route.substring(1),
                "--O", DIR.resolve("out").resolve(label).toString()));
        argv.addAll(Arrays.asList(rest));
        run(label, argv.toArray(new String[0]));
    }

    /** One run against nothing, which answers with the URI it was asked for. */
    static void uri(final String label, final String url, final String id, final String... rest) throws Exception {
        final List<String> argv = new ArrayList<>(List.of("--url", url, "--id", id,
                "--O", DIR.resolve("out").resolve(label).toString()));
        argv.addAll(Arrays.asList(rest));
        run(label, argv.toArray(new String[0]));
    }

    public static void main(final String[] args) throws Exception {
        Files.createDirectories(DIR.resolve("out"));
        final Path store = DIR.resolve("server.p12");
        final Process keytool = new ProcessBuilder(Path.of(System.getProperty("java.home"), "bin", "keytool").toString(),
                "-genkeypair", "-alias", "server", "-keyalg", "RSA", "-keysize", "2048", "-validity", "2",
                "-dname", "CN=localhost", "-ext", "SAN=dns:localhost", "-storetype", "PKCS12",
                "-keystore", store.toString(), "-storepass", "changeit", "-keypass", "changeit").inheritIO().start();
        if (keytool.waitFor() != 0) {
            throw new IllegalStateException("keytool failed");
        }
        // `SSLContexts.createDefault` reads the trust store the JVM was told of.
        System.setProperty("javax.net.ssl.trustStore", store.toString());
        System.setProperty("javax.net.ssl.trustStorePassword", "changeit");
        System.setProperty("javax.net.ssl.trustStoreType", "PKCS12");

        routes();
        ROUTES.forEach((path, route) -> System.out.printf("route\t%s\t%s\t%s\t%s%n", path, route.statuses(), route.location(),
                Base64.getEncoder().encodeToString(route.body().getBytes(StandardCharsets.UTF_8))));
        final HttpsServer secure = HttpsServer.create(new InetSocketAddress("localhost", 18443), 16);
        secure.setHttpsConfigurator(new HttpsConfigurator(tls(store)));
        secure.createContext("/", HtsgetReaderDump::serve);
        secure.setExecutor(java.util.concurrent.Executors.newCachedThreadPool());
        secure.start();
        final HttpServer plain = HttpServer.create(new InetSocketAddress("localhost", 18080), 16);
        plain.createContext("/", HtsgetReaderDump::serve);
        plain.setExecutor(java.util.concurrent.Executors.newCachedThreadPool());
        plain.start();

        // What the server answers.
        for (final String route : ROUTES.keySet()) {
            if (!route.equals("/flaky") && !route.equals("/broken") && !route.equals("/gone")) {
                fetch("fetch" + route.replace('/', '-'), route);
            }
        }
        fetch("fetch-flaky", "/flaky");
        fetch("fetch-broken", "/broken");
        fetch("fetch-unrouted", "/nothing-here");
        fetch("md5-checked-good", "/md5-good", "--check-md5", "true");
        fetch("md5-checked-case", "/md5-bad", "--check-md5", "true");
        fetch("md5-unchecked", "/md5-bad", "--check-md5", "false");
        fetch("md5-checked-none", "/data", "--check-md5", "true");
        fetch("md5-checked-number", "/md5-number", "--check-md5", "true");
        fetch("threads-mixed", "/mixed", "--reader-threads", "3");
        fetch("threads-two", "/two", "--reader-threads", "2");
        fetch("threads-more-than-blocks", "/data", "--reader-threads", "8");
        fetch("threads-empty", "/empty", "--reader-threads", "2");
        fetch("threads-ftp", "/ftp", "--reader-threads", "2");
        fetch("threads-no-url", "/no-url", "--reader-threads", "2");
        fetch("threads-bad-base64", "/bad-base64", "--reader-threads", "2");
        fetch("threads-refused", "/block-refused", "--reader-threads", "2");
        fetch("threads-no-urls", "/no-urls", "--reader-threads", "2");
        fetch("threads-zero", "/data", "--reader-threads", "0");
        fetch("threads-text", "/data", "--reader-threads", "x");
        run("output-nowhere", "--url", TLS, "--id", "data", "--O", DIR.resolve("nowhere").resolve("x").toString());
        run("output-nowhere-threads", "--url", TLS, "--id", "data", "--O", DIR.resolve("nowhere").resolve("x").toString(),
                "--reader-threads", "2");
        run("output-long-name", "--url", TLS, "--id", "data", "-output", DIR.resolve("out").resolve("output-long-name").toString());
        run("output-missing", "--url", TLS, "--id", "data");
        run("url-missing", "--id", "data", "--O", DIR.resolve("out").resolve("url-missing").toString());
        run("id-missing", "--url", TLS, "--O", DIR.resolve("out").resolve("id-missing").toString());

        // What the request refuses before it is sent.
        uri("header-interval", UNREACHABLE, "a", "--class", "header", "-L", "chr1");
        uri("header-field", UNREACHABLE, "a", "--class", "header", "--field", "SEQ");
        uri("header-tag", UNREACHABLE, "a", "--class", "header", "--tag", "NM");
        uri("header-notag", UNREACHABLE, "a", "--class", "header", "--notag", "NM");
        uri("header-format", UNREACHABLE, "a", "--class", "header", "--format", "BAM");
        uri("body-interval", UNREACHABLE, "a", "--class", "body", "-L", "chr1");
        for (final String format : new String[] {"BAM", "CRAM", "VCF", "BCF"}) {
            for (final String id : new String[] {"a.bam", "a.sam", "a.cram", "a.vcf", "a.vcf.gz", "a.vcf.bgz", "a.bcf", "a"}) {
                uri("format-" + format + "-" + id, UNREACHABLE, id, "--format", format);
            }
        }
        uri("format-lower", UNREACHABLE, "a", "--format", "bam");
        uri("class-upper", UNREACHABLE, "a", "--class", "BODY");
        uri("field-lower", UNREACHABLE, "a", "--field", "seq");
        uri("tags-overlap", UNREACHABLE, "a", "--tag", "NM", "--notag", "NM");
        uri("tags-overlap-two", UNREACHABLE, "a", "--tag", "XX", "--tag", "AS", "--tag", "NM", "--notag", "NM", "--notag", "AS", "--notag", "XX");

        // The URI the builder makes.
        final String[] urls = {"https://localhost:1", "https://localhost:1/", "https://localhost:1/base", "https://localhost:1/base/",
                "http://localhost:1/plain", "ftp://localhost:1/other", "localhost:1", "//localhost:1/path", "just-a-host.invalid",
                "https://localhost:1?x=1", "https://localhost:1/base?x=1&y=a%20b", "https://localhost:1/base#part",
                "https://user@localhost:1/base", "https://localhost:1/a%20b", "/only/a/path", "", "https://", "a b", "1a:b", "mailto:a@b.c",
                "https://[::1]:1/v6", "https://localhost:1/%zz"};
        int at = 0;
        for (final String url : urls) {
            uri("url-" + at, url, "id");
            uri("url-query-" + at, url, "id", "--format", "BCF");
            at++;
        }
        final String[] ids = {"plain", "a/b", "/leading", "trailing/", "a b", "a%20b", "a%zz", "a%", "a+b", "a&b", "a=b", "a?b", "a#b", "a;b", "a:b", "a@b",
                "a,b", "a'b", "a(b)", "a!b", "a$b", "a*b", "a~b", "a_b", "a-b", "a.b", "a|b", "a\"b", "a<b>", "a[b]", "a\\b", "a^b", "a`b", "é", "日本", "..", "%2F", "{x}", "a{b", "a}b", ""};
        at = 0;
        for (final String id : ids) {
            uri("id-" + at++, UNREACHABLE, id);
        }
        final String[] values = {"NM", "a b", "a%20b", "a%zz", "a+b", "a&b", "a=b", "a?b", "a#b", "a;b", "a:b", "a@b", "a,b", "a/b", "a'b", "a(b)", "a!b", "a$b", "a*b",
                "a~b", "a|b", "a\"b", "a<b>", "a[b]", "a\\b", "a^b", "a`b", "é", "{x}", ""};
        at = 0;
        for (final String value : values) {
            uri("tag-" + at, UNREACHABLE, "id", "--tag", value);
            uri("notag-" + at, UNREACHABLE, "id", "--notag", value);
            at++;
        }
        final String[] intervals = {"chr1", "chr1:5", "chr1:5+", "chr1:5-10", "chr1:1,000-2,000", "chr 1:5-10", "c&d:5-10", "HLA:A:5-10", "chr1:x-5", "chr1:5-x", "chr1:0-5", "chr1:10-5", "chr1:", ""};
        at = 0;
        for (final String interval : intervals) {
            uri("interval-" + at++, UNREACHABLE, "id", "-L", interval);
        }
        uri("interval-two", UNREACHABLE, "id", "-L", "chr1", "-L", "chr2");
        uri("fields-order", UNREACHABLE, "id", "--field", "QUAL", "--field", "QNAME", "--field", "SEQ", "--field", "QNAME");
        uri("fields-all", UNREACHABLE, "id", "--field", "QNAME", "--field", "FLAG", "--field", "RNAME", "--field", "POS", "--field", "MAPQ",
                "--field", "CIGAR", "--field", "RNEXT", "--field", "PNEXT", "--field", "TLEN", "--field", "SEQ", "--field", "QUAL");
        uri("tags-order", UNREACHABLE, "id", "--tag", "ZZ", "--tag", "AA", "--tag", "NM", "--tag", "MD", "--tag", "RG", "--tag", "AS", "--tag", "XS",
                "--tag", "BC", "--tag", "OQ", "--tag", "SA", "--tag", "CB", "--tag", "UB", "--tag", "AA", "--tag", "XA");
        uri("notags-order", UNREACHABLE, "id", "--notag", "ZZ", "--notag", "AA", "--notag", "NM", "--notag", "MD");
        uri("everything", UNREACHABLE, "id", "--format", "CRAM", "--class", "body", "-L", "chr2:7-9", "--field", "POS", "--field", "FLAG",
                "--tag", "NM", "--notag", "AS");
        secure.stop(0);
        plain.stop(0);
        // Apache's pool holds its threads, and nothing else ends the JVM.
        System.out.flush();
        System.exit(0);
    }
}
