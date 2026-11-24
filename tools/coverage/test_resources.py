#!/usr/bin/env python3
"""Replay the reference's own integration tests on both sides: criterion (1) of a bit-identical tool.

The ROADMAP asks for every applicable input under the pinned clone's `src/test/resources`. The
inputs that apply to a tool are the ones its `*IntegrationTest.java` hands it, with the arguments
it hands them with, so this reads those invocations out of the test sources and runs each one
against the pinned oracle and against the port, comparing what they write the way
`run_array.py` compares a covering array's rows (the same runners, the same masks).

The invocations are read STATICALLY: the test classes are parsed (tree-sitter) and their test
methods partially evaluated -- constants, `ArgumentsBuilder`, `IntegrationTestSpec`, String
arrays, helper methods, `@DataProvider` rows, loops and branches whose condition is known. Running
the Java tests instead would need their whole classpath and would run the reference's assertions
rather than ours. A test whose arguments cannot be evaluated is listed under `unextracted` with
the reason, so what the harness does not replay is counted rather than missed.

    # once per pinned release, on a host with the sources; the result is committed
    python3 tools/coverage/test_resources.py extract --gatk ../gatk-4.7.0.0
    # what CI runs (stdlib only), one shard of the replay
    python3 tools/coverage/test_resources.py run --gatk gatk-src --port target/release/gatk --shard 0/8
    python3 tools/coverage/test_resources.py merge shard-*.json

`extract` needs `pip install tree-sitter==0.26.1 tree-sitter-java==0.23.5`; nothing else does.
"""

import argparse
import concurrent.futures
import json
import re
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_array  # noqa: E402

DATA = HERE / "test_resources"
INVOCATIONS = DATA / "gatk.json"
EXCLUSIONS = DATA / "exclusions.json"
MEASURED = DATA / "measured.json"
MANIFEST = HERE.parent / "conformance" / "manifest.json"

# Where the reference's checkout is mounted on both sides; GATK's tests run with it as `user.dir`.
ROOT = "/work/fixtures"
OUT = "/work/out"
# TestNG groups whose tests read from a bucket or need credentials, which neither side has.
CLOUD_GROUPS = {"cloud", "bucket", "gcs", "spark_cloud", "cloud_requester_pays"}
RELATIVE_RESOURCE = re.compile(r"(?<![\w/.])src/(test|main)/resources/")


# ---------------------------------------------------------------------------------------------
# Values of the partial evaluator.


class Unresolved(Exception):
    """An expression whose value the sources alone do not determine."""


class Return(Exception):
    def __init__(self, value):
        self.value = value


class Thrown(Exception):
    """A `throw` reached: the test method ends there."""


class Poison:
    """A variable whose value could not be computed; reading it is what fails, not binding it."""

    def __init__(self, reason):
        self.reason = reason


class ClassRef:
    def __init__(self, name, info=None):
        self.name, self.info = name, info


class EnumValue:
    def __init__(self, name):
        self.name = name


class Builder:
    """`ArgumentsBuilder`: a list of command-line tokens."""

    def __init__(self, tokens=()):
        self.tokens = list(tokens)

    def raw(self, value):
        if isinstance(value, str):
            self.tokens += value.split()
        else:
            self.tokens.append(jstr(value))
        return self

    def add(self, name, value):
        # Only `add(name, File)` makes the path absolute, and a String can be `%s`; a relative
        # resource path is made absolute where the invocation is recorded instead.
        self.raw("--" + jstr(name))
        return self.raw(jstr(value))


class Spec:
    """`IntegrationTestSpec`: an argument string whose `%s` are the outputs it creates."""

    def __init__(self, args, outputs):
        self.args, self.outputs, self.extension = args, outputs, ".tmp"


class Lambda:
    def __init__(self, node, env, this):
        self.node, self.env, self.this = node, env, this


def jdouble(value):
    """`Double.toString`, for the doubles a test concatenates into an argument."""
    if value != value or value in (float("inf"), float("-inf")):
        return {True: "NaN"}.get(value != value, "Infinity" if value > 0 else "-Infinity")
    if value == 0:
        return "0.0"
    if 1e-3 <= abs(value) < 1e7:
        text = repr(value)
        return text if "." in text or "e" in text else text + ".0"
    mantissa, exponent = f"{value:.17e}".split("e")
    shortest = repr(value)
    digits = shortest.replace("-", "").split("e")[0].replace(".", "").lstrip("0").rstrip("0") or "0"
    sign = "-" if value < 0 else ""
    return f"{sign}{digits[0]}.{digits[1:] or '0'}E{int(exponent)}"


def jstr(value):
    """`String.valueOf`."""
    if isinstance(value, Poison):
        raise Unresolved(value.reason)
    if value is None:
        return "null"
    if value is True:
        return "true"
    if value is False:
        return "false"
    if isinstance(value, float):
        return jdouble(value)
    if isinstance(value, (int, str)):
        return str(value)
    if isinstance(value, EnumValue):
        return value.name
    if isinstance(value, ClassRef):
        return f"class {value.name}"
    if isinstance(value, list):
        return "[" + ", ".join(jstr(v) for v in value) + "]"
    if isinstance(value, Builder):
        return " ".join(value.tokens)
    raise Unresolved(f"no string form for {type(value).__name__}")


def absolute(value):
    """`File.getAbsolutePath()`: GATK's tests run from the checkout's root."""
    text = jstr(value)
    if text.startswith(("/", "gs://", "gendb", "hdfs:", "file:", "http")) or not text:
        return text
    return f"{ROOT}/{text}"


def normalized(path):
    """What `new File(path)` and `Paths.get(path)` keep: no repeated or trailing separator."""
    if "://" in path:
        return path
    path = re.sub(r"/+", "/", path)
    return path.rstrip("/") if len(path) > 1 else path


def join_path(parent, child):
    parent, child = jstr(parent), jstr(child)
    if child.startswith("/") or not parent:
        return normalized(child)
    return normalized(parent + "/" + child)


def jformat(fmt, args):
    """`String.format` for the conversions tests use."""
    args = list(args)
    position = 0

    def one(match):
        nonlocal position
        index, flags, width, precision, conversion = match.groups()
        if conversion == "%":
            return "%"
        if conversion == "n":
            return "\n"
        if index:
            value = args[int(index[:-1]) - 1]
        else:
            if position >= len(args):
                raise Unresolved("String.format with too few arguments")
            value = args[position]
            position += 1
        if conversion in "sS":
            text = jstr(value)
        elif conversion == "d":
            text = str(int(value))
        elif conversion == "f":
            text = f"{float(value):.{int(precision[1:]) if precision else 6}f}"
        elif conversion == "b":
            text = jstr(bool(value))
        else:
            raise Unresolved(f"String.format conversion %{conversion}")
        if width:
            text = text.ljust(int(width)) if "-" in (flags or "") else text.rjust(int(width))
        return text.upper() if conversion == "S" else text

    return re.sub(r"%(\d+\$)?([-#+ 0,(<]*)?(\d+)?(\.\d+)?([a-zA-Z%])", one, fmt)


def escape_expressions(args):
    """`Utils.escapeExpressions`: split on spaces, keeping a quoted expression whole."""
    for delimiter in ("'", '"'):
        if delimiter in args:
            out, split = [], args.split(delimiter)
            for i in range(0, len(split) - 1, 2):
                if split[i].strip():
                    out += split[i].strip().split()
                out.append(split[i + 1])
            if len(split) % 2 == 1 and split[-1].strip():
                out += split[-1].strip().split()
            return out
    return args.split()


def java_string(node):
    """A string literal's value, text blocks included."""
    raw = node.text.decode()
    if raw.startswith('"""'):
        body = raw[3:-3].split("\n", 1)[1] if "\n" in raw[3:-3] else raw[3:-3]
        lines = body.split("\n")
        indent = min((len(l) - len(l.lstrip()) for l in lines if l.strip()), default=0)
        body = "\n".join(l[indent:].rstrip() for l in lines)
    else:
        body = raw[1:-1]
    return re.sub(
        r"\\(u+[0-9a-fA-F]{4}|[0-7]{1,3}|.)",
        lambda m: (
            chr(int(m.group(1).lstrip("u"), 16)) if m.group(1)[0] == "u"
            else chr(int(m.group(1), 8)) if m.group(1)[0] in "01234567"
            else {"n": "\n", "t": "\t", "r": "\r", "b": "\b", "f": "\f", "s": " ", "\n": ""}.get(m.group(1), m.group(1))
        ),
        body,
    )


# ---------------------------------------------------------------------------------------------
# The sources: every class of the checkout's main and test trees, by simple name.


class ClassInfo:
    def __init__(self, node, package, path, outer=None, static_imports=()):
        self.node, self.package, self.path, self.outer = node, package, path, outer
        self.static_imports = static_imports
        self.name = node.child_by_field_name("name").text.decode()
        superclass = node.child_by_field_name("superclass")
        self.superclass = None
        if superclass is not None:
            type_node = kids(superclass)[0]
            if type_node.type == "generic_type":
                type_node = kids(type_node)[0]
            self.superclass = type_node.text.decode().split(".")[-1]
        self.fields, self.methods, self.enum_constants = {}, {}, []
        self.annotations = annotations(node)
        self.abstract = has_modifier(node, "abstract")
        body = node.child_by_field_name("body")
        for child in kids(body) if body is not None else ():
            if child.type == "field_declaration":
                for declarator in child.children_by_field_name("declarator"):
                    self.fields[declarator.child_by_field_name("name").text.decode()] = (child, declarator)
            elif child.type == "method_declaration":
                self.methods.setdefault(child.child_by_field_name("name").text.decode(), []).append(child)
            elif child.type == "enum_body_declarations":
                for inner in kids(child):
                    if inner.type == "method_declaration":
                        self.methods.setdefault(inner.child_by_field_name("name").text.decode(), []).append(inner)
            elif child.type == "enum_constant":
                self.enum_constants.append(child.child_by_field_name("name").text.decode())


def kids(node):
    """Named children, comments left out: a comment inside an array initializer is a child too."""
    return [c for c in node.named_children if c.type not in ("line_comment", "block_comment")]


def has_modifier(node, word):
    for child in node.children:
        if child.type == "modifiers":
            return any(c.type == word for c in child.children)
    return False


def annotations(node):
    """{name: {key: value node}} of a declaration's annotations."""
    found = {}
    for child in node.children:
        if child.type != "modifiers":
            continue
        for annotation in kids(child):
            if annotation.type not in ("annotation", "marker_annotation"):
                continue
            name = annotation.child_by_field_name("name").text.decode().split(".")[-1]
            values = {}
            arguments = annotation.child_by_field_name("arguments")
            for pair in kids(arguments) if arguments is not None else ():
                if pair.type == "element_value_pair":
                    values[pair.child_by_field_name("key").text.decode()] = pair.child_by_field_name("value")
                else:
                    values["value"] = pair
            found[name] = values
    return found


class Sources:
    def __init__(self, roots):
        import tree_sitter
        import tree_sitter_java

        self.parser = tree_sitter.Parser(tree_sitter.Language(tree_sitter_java.language()))
        self.classes = {}
        for root in roots:
            for path in sorted(Path(root).rglob("*.java")):
                self.load(path)

    def load(self, path):
        tree = self.parser.parse(path.read_bytes())
        package, static_imports = "", []
        for child in kids(tree.root_node):
            if child.type == "package_declaration":
                package = kids(child)[0].text.decode()
            if child.type == "import_declaration" and any(c.type == "static" for c in child.children):
                parts = kids(child)[0].text.decode().split(".")
                star = any(c.type == "asterisk" for c in child.children)
                static_imports.append((parts[-1], "*") if star else (parts[-2], parts[-1]))

        def walk(node, outer):
            for child in kids(node):
                if child.type in ("class_declaration", "enum_declaration", "interface_declaration"):
                    info = ClassInfo(child, package, path, outer, static_imports)
                    self.classes.setdefault(info.name, info)
                    body = child.child_by_field_name("body")
                    if body is not None:
                        walk(body, info)
                        for inner in kids(body):
                            if inner.type == "enum_body_declarations":
                                walk(inner, info)

        walk(tree.root_node, None)

    def chain(self, info):
        """The class and its superclasses that the sources hold, nearest first."""
        seen = []
        while info is not None and info not in seen:
            seen.append(info)
            info = self.classes.get(info.superclass) if info.superclass else None
        return seen


# ---------------------------------------------------------------------------------------------
# The evaluator.


class Run:
    """One test method's execution: the temporary files it made and the commands it ran."""

    def __init__(self, test, sources):
        self.test, self.sources = test, sources
        self.counter, self.touch, self.mkdir, self.invocations = 0, [], [], []
        self.steps = 0
        self.instance = {}
        # The methods being evaluated: an overload of the same arity (`f(String)` calling
        # `f(File)`) cannot be told apart without types, and must not be taken for itself.
        self.stack = []

    def temp(self, prefix, suffix, directory=False):
        self.counter += 1
        name = re.sub(r"[^\w.+-]", "_", f"{jstr(prefix)}{self.counter}{'' if suffix is None else jstr(suffix)}")
        path = f"{OUT}/{name}"
        (self.mkdir if directory else self.touch).append(path)
        return path


def gatk_extension(extension):
    """`IOUtils.createTempFile`/`createTempPath` prepend the dot an extension lacks, unlike
    `File.createTempFile`: `createTempFile("test", "fasta")` is `test<n>.fasta`."""
    extension = jstr(extension)
    return extension if extension.startswith(".") else f".{extension}"


class Interpreter:
    def __init__(self, run, this):
        self.run, self.this, self.sources = run, this, run.sources

    # --- names --------------------------------------------------------------------------------

    def tool_name(self):
        for info in self.sources.chain(self.this):
            if info.name in FRAMEWORK:
                break
            for method in info.methods.get("getTestedToolName", []):
                return jstr(self.call(method, [], info))
        return self.tested_class_name()

    def overridden(self, name, args):
        """A test class's own definition of a framework helper, which wins as it does in Java."""
        for owner in self.sources.chain(self.this):
            if owner.name in FRAMEWORK:
                return NOT_NATIVE
            for method in owner.methods.get(name, []):
                if len(kids(method.child_by_field_name("parameters"))) == len(args) and method not in self.run.stack:
                    return self.call(method, args, owner)
        return NOT_NATIVE

    def tested_class_name(self):
        value = self.overridden("getTestedClassName", [])
        if value is not NOT_NATIVE:
            return jstr(value)
        name = self.this.name
        for suffix in ("IntegrationTest", "UnitTest", "Test"):
            if name.endswith(suffix):
                return name[: -len(suffix)]
        return name

    def field(self, info, name):
        """A field's value, looked up through the class chain and the enclosing classes."""
        for owner in self.sources.chain(info):
            key = (owner.name, name)
            if key in self.run.instance:
                return self.run.instance[key]
            if name in owner.enum_constants:
                return EnumValue(name)
            if name in owner.fields:
                declaration, declarator = owner.fields[name]
                value_node = declarator.child_by_field_name("value")
                if value_node is None:
                    raise Unresolved(f"field {owner.name}.{name} has no initializer")
                try:
                    value = Interpreter(self.run, owner).eval(value_node, {"scopes": []})
                except (Unresolved, Return, Thrown, ArithmeticError, LookupError, TypeError, ValueError, AttributeError) as error:
                    value = Poison(str(error))
                self.run.instance[key] = value
                return value
        if info is not None and info.outer is not None:
            return self.field(info.outer, name)
        for owner, member in info.static_imports if info is not None else ():
            if member in (name, "*") and owner in self.sources.classes and self.sources.classes[owner] is not info:
                try:
                    return self.field(self.sources.classes[owner], name)
                except Unresolved:
                    continue
        raise Unresolved(f"unknown name {name}")

    def lookup(self, name, env):
        for scope in reversed(env["scopes"]):
            if name in scope:
                return scope[name]
        try:
            return self.field(self.this, name)
        except Unresolved:
            if name in self.sources.classes:
                return ClassRef(name, self.sources.classes[name])
            if name[:1].isupper() and not name.isupper():
                return ClassRef(name)
            raise

    # --- statements ---------------------------------------------------------------------------

    def call(self, method, args, owner, scopes=None):
        params = method.child_by_field_name("parameters")
        scope = {}
        names = []
        for param in kids(params) if params is not None else ():
            if param.type == "formal_parameter":
                names.append((param.child_by_field_name("name").text.decode(), False))
            elif param.type == "spread_parameter":
                names.append((kids(param)[-1].child_by_field_name("name").text.decode(), True))
        for i, (name, spread) in enumerate(names):
            if spread:
                rest = args[i:]
                scope[name] = rest[0] if len(rest) == 1 and isinstance(rest[0], list) else rest
            elif i < len(args):
                scope[name] = args[i]
        body = method.child_by_field_name("body")
        if body is None:
            raise Unresolved(f"abstract method {owner.name}.{method.child_by_field_name('name').text.decode()}")
        inner = Interpreter(self.run, owner if has_modifier(method, "static") else self.this)
        env = {"scopes": (scopes or []) + [scope]}
        self.run.stack.append(method)
        try:
            inner.block(body, env)
        except Return as done:
            return done.value
        finally:
            self.run.stack.pop()
        return None

    def block(self, node, env):
        env["scopes"].append({})
        try:
            for statement in kids(node):
                self.statement(statement, env)
        finally:
            env["scopes"].pop()

    def bind(self, env, name, value):
        env["scopes"][-1][name] = value

    def assign(self, env, name, value):
        for scope in reversed(env["scopes"]):
            if name in scope:
                scope[name] = value
                return
        for owner in self.sources.chain(self.this):
            if name in owner.fields:
                self.run.instance[(owner.name, name)] = value
                return
        self.bind(env, name, value)

    def poison(self, node, env, reason):
        """Every variable a statement we cannot evaluate might change."""
        for child in walk(node):
            target = None
            if child.type == "assignment_expression":
                target = child.child_by_field_name("left")
            elif child.type == "method_invocation":
                target = child.child_by_field_name("object")
            elif child.type == "variable_declarator":
                target = child.child_by_field_name("name")
            if target is not None and target.type == "identifier":
                self.assign(env, target.text.decode(), Poison(reason))
        if any(c.type == "method_invocation" and c.child_by_field_name("name").text.decode() in RUNNERS
               for c in walk(node)):
            raise Unresolved(reason)

    def statement(self, node, env):
        self.run.steps += 1
        if self.run.steps > 200_000:
            raise Unresolved("evaluation budget exhausted")
        kind = node.type
        if kind in ("line_comment", "block_comment", "empty_statement", "assert_statement",
                    "local_class_declaration", "class_declaration", "labeled_statement_"):
            return
        if kind == "block":
            return self.block(node, env)
        if kind == "local_variable_declaration":
            for declarator in node.children_by_field_name("declarator"):
                name = declarator.child_by_field_name("name").text.decode()
                value = declarator.child_by_field_name("value")
                try:
                    self.bind(env, name, None if value is None else self.eval(value, env))
                except Unresolved as error:
                    self.bind(env, name, Poison(str(error)))
                    if value is not None and any(
                        c.type == "method_invocation" and c.child_by_field_name("name").text.decode() in RUNNERS
                        for c in walk(value)
                    ):
                        raise
            return
        if kind == "expression_statement":
            try:
                self.eval(kids(node)[0], env)
            except Unresolved as error:
                self.poison(node, env, str(error))
            return
        if kind == "if_statement":
            try:
                condition = truthy(self.eval(node.child_by_field_name("condition"), env))
            except Unresolved as error:
                return self.poison(node, env, f"branch on {error}")
            branch = node.child_by_field_name("consequence" if condition else "alternative")
            if branch is not None:
                self.statement(branch, env)
            return
        if kind == "enhanced_for_statement":
            try:
                items = self.eval(node.child_by_field_name("value"), env)
                if isinstance(items, Poison):
                    raise Unresolved(items.reason)
                items = list(items)
            except (Unresolved, TypeError) as error:
                return self.poison(node, env, f"loop over {error}")
            name = node.child_by_field_name("name").text.decode()
            for item in items:
                env["scopes"].append({name: item})
                try:
                    self.statement(node.child_by_field_name("body"), env)
                finally:
                    env["scopes"].pop()
            return
        if kind == "for_statement":
            env["scopes"].append({})
            try:
                for init in node.children_by_field_name("init"):
                    self.statement(init, env) if init.type.endswith("declaration") else self.eval(init, env)
                for _ in range(10_000):
                    condition = node.child_by_field_name("condition")
                    if condition is not None and not truthy(self.eval(condition, env)):
                        break
                    self.statement(node.child_by_field_name("body"), env)
                    for update in node.children_by_field_name("update"):
                        self.eval(update, env)
            except Unresolved as error:
                self.poison(node, env, f"loop on {error}")
            finally:
                env["scopes"].pop()
            return
        if kind == "try_statement":
            try:
                self.statement(node.child_by_field_name("body"), env)
            finally:
                for child in kids(node):
                    if child.type == "finally_clause":
                        self.statement(kids(child)[-1], env)
            return
        if kind == "try_with_resources_statement":
            env["scopes"].append({})
            try:
                resources = node.child_by_field_name("resources")
                for resource in kids(resources):
                    if resource.type == "resource" and resource.child_by_field_name("name") is not None:
                        try:
                            value = self.eval(resource.child_by_field_name("value"), env)
                        except Unresolved as error:
                            value = Poison(str(error))
                        self.bind(env, resource.child_by_field_name("name").text.decode(), value)
                self.statement(node.child_by_field_name("body"), env)
            finally:
                env["scopes"].pop()
            return
        if kind == "return_statement":
            raise Return(self.eval(kids(node)[0], env) if kids(node) else None)
        if kind == "throw_statement":
            raise Thrown()
        if kind in ("while_statement", "do_statement", "switch_expression", "switch_statement",
                    "synchronized_statement", "yield_statement", "break_statement", "continue_statement"):
            return self.poison(node, env, f"{kind} is not evaluated")
        raise Unresolved(f"statement {kind}")

    # --- expressions --------------------------------------------------------------------------

    def eval(self, node, env):
        kind = node.type
        if kind == "string_literal":
            return java_string(node)
        if kind == "character_literal":
            return java_string(node)
        if kind in ("decimal_integer_literal", "hex_integer_literal", "octal_integer_literal", "binary_integer_literal"):
            return int(node.text.decode().rstrip("lL").replace("_", ""), 0 if kind != "decimal_integer_literal" else 10)
        if kind == "decimal_floating_point_literal":
            return float(node.text.decode().rstrip("fFdD").replace("_", ""))
        if kind == "true":
            return True
        if kind == "false":
            return False
        if kind == "null_literal":
            return None
        if kind == "this":
            return ClassRef("this", self.this)
        if kind == "identifier":
            value = self.lookup(node.text.decode(), env)
            if isinstance(value, Poison):
                raise Unresolved(value.reason)
            return value
        if kind in ("parenthesized_expression", "cast_expression"):
            return self.eval(kids(node)[-1] if kind == "cast_expression" else kids(node)[0], env)
        if kind == "class_literal":
            name = kids(node)[0].text.decode().split(".")[-1].split("<")[0]
            return ClassRef(name, self.sources.classes.get(name))
        if kind == "field_access":
            return self.field_access(node, env)
        if kind == "array_access":
            return self.eval(node.child_by_field_name("array"), env)[int(self.eval(node.child_by_field_name("index"), env))]
        if kind == "binary_expression":
            return self.binary(node, env)
        if kind == "unary_expression":
            operand = self.eval(node.child_by_field_name("operand"), env)
            operator = node.child_by_field_name("operator").text.decode()
            return {"!": lambda v: not v, "-": lambda v: -v, "+": lambda v: v}[operator](operand)
        if kind == "ternary_expression":
            branch = "consequence" if truthy(self.eval(node.child_by_field_name("condition"), env)) else "alternative"
            return self.eval(node.child_by_field_name(branch), env)
        if kind == "assignment_expression":
            left, right = node.child_by_field_name("left"), node.child_by_field_name("right")
            operator = node.child_by_field_name("operator").text.decode()
            value = self.eval(right, env)
            if operator == "+=":
                current = self.eval(left, env)
                value = current + value if isinstance(current, (int, float)) and isinstance(value, (int, float)) else jstr(current) + jstr(value)
            elif operator != "=":
                raise Unresolved(f"assignment {operator}")
            if left.type == "identifier":
                self.assign(env, left.text.decode(), value)
            elif left.type == "field_access" and left.child_by_field_name("object").type == "this":
                self.assign(env, left.child_by_field_name("field").text.decode(), value)
            elif left.type == "array_access":
                self.eval(left.child_by_field_name("array"), env)[int(self.eval(left.child_by_field_name("index"), env))] = value
            else:
                raise Unresolved("assignment to an expression")
            return value
        if kind == "update_expression":
            target = next(c for c in kids(node))
            delta = 1 if "++" in node.text.decode() else -1
            value = self.eval(target, env)
            self.assign(env, target.text.decode(), value + delta)
            return value
        if kind == "array_initializer":
            return [self.eval(child, env) for child in kids(node)]
        if kind == "array_creation_expression":
            value = node.child_by_field_name("value")
            if value is not None:
                return self.eval(value, env)
            dimensions = node.child_by_field_name("dimensions")
            sizes = [c for c in kids(node) if c.type == "dimensions_expr"]
            if sizes:
                return [None] * int(self.eval(sizes[0].named_children[0], env))
            return []
        if kind == "object_creation_expression":
            return self.construct(node, env)
        if kind == "method_invocation":
            return self.invoke(node, env)
        if kind == "lambda_expression":
            return Lambda(node, env, self.this)
        if kind == "instanceof_expression":
            raise Unresolved("instanceof")
        raise Unresolved(f"expression {kind}")

    def binary(self, node, env):
        operator = node.child_by_field_name("operator").text.decode()
        left = self.eval(node.child_by_field_name("left"), env)
        if operator == "&&":
            return truthy(left) and truthy(self.eval(node.child_by_field_name("right"), env))
        if operator == "||":
            return truthy(left) or truthy(self.eval(node.child_by_field_name("right"), env))
        right = self.eval(node.child_by_field_name("right"), env)
        if operator == "+":
            if isinstance(left, (int, float)) and isinstance(right, (int, float)) and not isinstance(left, bool):
                return left + right
            return jstr(left) + jstr(right)
        if operator in ("==", "!="):
            same = left is right if not isinstance(left, (int, float, str, type(None), bool)) else left == right
            if isinstance(left, EnumValue) and isinstance(right, EnumValue):
                same = left.name == right.name
            return same if operator == "==" else not same
        numbers = {"-": lambda a, b: a - b, "*": lambda a, b: a * b, "<": lambda a, b: a < b,
                   ">": lambda a, b: a > b, "<=": lambda a, b: a <= b, ">=": lambda a, b: a >= b,
                   "%": lambda a, b: a % b,
                   "/": lambda a, b: a // b if isinstance(a, int) and isinstance(b, int) else a / b}
        if operator in numbers:
            return numbers[operator](left, right)
        raise Unresolved(f"operator {operator}")

    def field_access(self, node, env):
        target, name = node.child_by_field_name("object"), node.child_by_field_name("field").text.decode()
        if target.type in ("this", "super"):
            return self.field(self.this, name)
        try:
            owner = self.eval(target, env)
        except Unresolved:
            # A qualified class name: `org.x.Tool.CONSTANT`.
            owner = ClassRef(target.text.decode().split(".")[-1], self.sources.classes.get(target.text.decode().split(".")[-1]))
        if isinstance(owner, ClassRef):
            if (owner.name, name) in (("File", "separator"), ("File", "pathSeparator")):
                return "/" if name == "separator" else ":"
            if owner.info is None:
                raise Unresolved(f"constant {owner.name}.{name} is outside the sources")
            if name in self.sources.classes and self.sources.classes[name].outer is owner.info:
                return ClassRef(name, self.sources.classes[name])
            return self.field(owner.info, name)
        if name == "length" and isinstance(owner, (list, str)):
            return len(owner)
        raise Unresolved(f"field {name} of a {type(owner).__name__}")

    def arguments(self, node, env):
        arguments = node.child_by_field_name("arguments")
        return [self.eval(child, env) for child in kids(arguments)] if arguments is not None else []

    def construct(self, node, env):
        type_node = node.child_by_field_name("type")
        name = type_node.text.decode().split("<")[0].split(".")[-1]
        if kids(node)[-1].type == "class_body":
            if name in ("ArgumentsBuilder",):
                raise Unresolved("anonymous ArgumentsBuilder")
        args = self.arguments(node, env)
        if name == "File":
            return normalized(jstr(args[0])) if len(args) == 1 else join_path(args[0], args[1])
        if name in ("GATKPath", "PicardHtsPath", "HtsPath", "URI"):
            return jstr(args[0])
        if name == "ArgumentsBuilder":
            builder = Builder()
            for arg in (args[0] if args else []):
                builder.raw(arg)
            return builder
        if name == "IntegrationTestSpec":
            outputs = len(args[1]) if isinstance(args[1], list) else int(args[1])
            return Spec(jstr(args[0]), outputs)
        if name in ("ArrayList", "LinkedList", "HashSet", "LinkedHashSet", "TreeSet"):
            items = list(args[0]) if args and isinstance(args[0], (list, tuple)) else []
            return sorted(items) if name == "TreeSet" else items
        if name in ("String", "StringBuilder"):
            return jstr(args[0]) if args else ""
        if name == "SimpleInterval":
            if len(args) == 1:
                return jstr(args[0])
            return f"{jstr(args[0])}:{jstr(args[1])}-{jstr(args[2])}"
        if name in ("Main",):
            return ClassRef("Main")
        if name in self.sources.classes and name.endswith("IntegrationTest"):
            return ClassRef("this", self.sources.classes[name])
        raise Unresolved(f"new {name}")

    def invoke(self, node, env):
        name = node.child_by_field_name("name").text.decode()
        target_node = node.child_by_field_name("object")
        arguments = node.child_by_field_name("arguments")
        arg_nodes = kids(arguments) if arguments is not None else []
        if target_node is None or target_node.type in ("this", "super"):
            return self.own_method(name, arg_nodes, env, super_call=target_node is not None and target_node.type == "super")
        try:
            target = self.eval(target_node, env)
        except Unresolved:
            if target_node.type in ("identifier", "field_access", "scoped_identifier") and target_node.text.decode()[:1].isupper():
                target = ClassRef(target_node.text.decode().split(".")[-1], self.sources.classes.get(target_node.text.decode().split(".")[-1]))
            else:
                self.lambdas(arg_nodes, env)
                raise
        if isinstance(target, ClassRef) and target.name == "this":
            return Interpreter(self.run, target.info).own_method(name, arg_nodes, env, interpreter=self)
        args = [self.eval(child, env) for child in arg_nodes] if name not in LAZY else None
        try:
            return self.member(target, name, args, arg_nodes, env)
        except Unresolved:
            self.lambdas(arg_nodes, env)
            raise

    def lambdas(self, arg_nodes, env):
        """Run the zero-argument lambdas of a call we cannot evaluate: `assertThrows(X.class, () -> runCommandLine(args))`."""
        for child in arg_nodes:
            if child.type == "lambda_expression":
                self.run_lambda(Lambda(child, env, self.this))

    def run_lambda(self, fn, args=()):
        body = fn.node.child_by_field_name("body")
        params = fn.node.child_by_field_name("parameters")
        names = [n.text.decode() for n in walk(params) if n.type == "identifier"] if params is not None else []
        names = names[-len(args):] if args else []
        env = {"scopes": fn.env["scopes"] + [dict(zip(names, args))]}
        inner = Interpreter(self.run, fn.this)
        try:
            if body.type == "block":
                inner.block(body, env)
                return None
            return inner.eval(body, env)
        except Return as done:
            return done.value
        except Thrown:
            return None

    def own_method(self, name, arg_nodes, env, super_call=False, interpreter=None):
        caller = interpreter or self
        args = [caller.eval(child, env) for child in arg_nodes]
        if name in RUNNERS and name != "executeTest":
            return self.runner(name, args)
        native = self.native(name, args)
        if native is not NOT_NATIVE:
            return native
        chain = self.sources.chain(self.this)
        if super_call:
            chain = chain[1:]
        for owner in chain + ([self.this.outer] if self.this.outer else []):
            for method in owner.methods.get(name, []):
                params = method.child_by_field_name("parameters")
                count = len([p for p in kids(params) if p.type in ("formal_parameter", "spread_parameter")])
                spread = any(p.type == "spread_parameter" for p in kids(params))
                if count == len(args) or (spread and len(args) >= count - 1):
                    if method.child_by_field_name("body") is None or method in self.run.stack:
                        continue
                    return Interpreter(self.run, self.this).call(method, args, owner)
        raise Unresolved(f"method {name}/{len(args)}")

    def native(self, name, args):
        """The test framework's own helpers, unless the test class overrides them."""
        if name not in ("getTestedClassName", "getTestedToolName", "getClass"):
            value = self.overridden(name, args)
            if value is not NOT_NATIVE:
                return value
        if name in ("createTempFile", "createTempPath"):
            return self.run.temp(args[0] if args else "default", gatk_extension(args[1]) if len(args) > 1 else ".tmp")
        if name == "createTempDir":
            return self.run.temp(args[0] if args else "tmp", "", directory=True)
        if name in ("getSafeNonExistentFile", "getSafeNonExistentPath", "getSafeNonExistentGATKPath"):
            directory = self.run.temp("nonExistentFileHolder", "", directory=True)
            return f"{directory}/{jstr(args[0])}"
        if name == "getToolTestDataDir":
            return f"src/test/resources/{self.this.package.replace('.', '/')}/{self.tested_class_name()}/"
        if name == "getTestFile":
            return join_path(self.native("getToolTestDataDir", []), args[0])
        if name == "getTestFileGATKPath":
            return f"{self.native('getToolTestDataDir', [])}/{jstr(args[0])}"
        if name == "getTestedClassName":
            return self.tested_class_name()
        if name == "getTestedToolName":
            return self.tool_name()
        if name == "getClass":
            return ClassRef(self.this.name, self.this)
        return NOT_NATIVE

    def runner(self, name, args):
        """`runCommandLine` and its relatives: the invocation this harness exists to read."""
        if name in ("runCommandLine", "runCommandLineWithoutVerbosity") or name == "runToolInNewJVM":
            if name == "runToolInNewJVM":
                tool, tokens = jstr(args[0]), args[1]
            else:
                tokens = args[0]
                tool = jstr(args[1]) if len(args) > 1 else self.tool_name()
            tokens = tokens.tokens if isinstance(tokens, Builder) else [jstr(t) for t in tokens]
            if name == "runCommandLine":
                tokens = inject_verbosity(tokens)
            self.record(tool, tokens)
            return None
        if name == "instanceMain":
            tokens = [jstr(t) for t in args[0]]
            self.record(tokens[0], tokens[1:])
            return None
        raise Unresolved(f"runner {name}")

    def record(self, tool, tokens):
        tokens = [RELATIVE_RESOURCE.sub(lambda m: f"{ROOT}/src/{m.group(1)}/resources/", jstr(t)) for t in tokens]
        self.run.invocations.append({
            "tool": tool, "args": tokens,
            "touch": list(self.run.touch), "mkdir": list(self.run.mkdir),
        })

    def member(self, target, name, args, arg_nodes, env):
        if isinstance(target, Poison):
            raise Unresolved(target.reason)
        if isinstance(target, Builder):
            return builder_method(target, name, args)
        if isinstance(target, Spec):
            if name == "executeTest":
                tester = args[1] if len(args) > 1 else None
                interpreter = Interpreter(self.run, tester.info) if isinstance(tester, ClassRef) and tester.info else self
                outputs = [self.run.temp(f"walktest.tmp_param.{i}", target.extension) for i in range(target.outputs)]
                interpreter.record_spec(target, outputs)
                return None
            if name == "setOutputFileExtension":
                target.extension = jstr(args[0])
                return None
            if name.startswith("set"):
                return None
            if name == "getArgs":
                return target.args
            raise Unresolved(f"IntegrationTestSpec.{name}")
        if isinstance(target, Lambda):
            return self.run_lambda(target, args)
        if isinstance(target, ClassRef):
            return self.static(target, name, args)
        if isinstance(target, EnumValue):
            if name in ("name", "toString"):
                return target.name
            raise Unresolved(f"enum method {name}")
        if isinstance(target, str):
            return string_method(target, name, args)
        if isinstance(target, list) and name == "forEach" and isinstance(args[0], Lambda):
            for item in list(target):
                self.run_lambda(args[0], [item])
            return None
        if isinstance(target, list):
            return list_method(target, name, args)
        if isinstance(target, (int, float)) and name in ("toString", "intValue", "doubleValue", "longValue"):
            return jstr(target) if name == "toString" else target
        raise Unresolved(f"method {name} of a {type(target).__name__}")

    def record_spec(self, spec, outputs):
        tokens = escape_expressions(jformat(spec.args, outputs))
        self.record(self.tool_name(), inject_verbosity(tokens))

    def static(self, target, name, args):
        cls = target.name
        if cls == "String" and name == "format":
            return jformat(jstr(args[0]), args[1:])
        if cls == "String" and name == "join":
            items = args[1] if len(args) == 2 and isinstance(args[1], list) else args[1:]
            return jstr(args[0]).join(jstr(v) for v in items)
        if cls == "String" and name == "valueOf" or cls in ("Integer", "Long", "Double", "Boolean") and name == "toString":
            return jstr(args[0])
        if cls in ("Integer", "Long") and name in ("parseInt", "parseLong", "valueOf"):
            return int(jstr(args[0]))
        if cls == "Double" and name in ("parseDouble", "valueOf"):
            return float(jstr(args[0]))
        if (cls, name) in (("Arrays", "asList"), ("List", "of"), ("ImmutableSet", "of"), ("Stream", "concat"), ("Lists", "newArrayList"), ("Set", "of"),
                           ("Collections", "singletonList"), ("Collections", "singleton"), ("ImmutableList", "of"),
                           ("Stream", "of"), ("Sets", "newHashSet")):
            if name == "concat":
                return list(args[0]) + list(args[1])
            return list(args[0]) if len(args) == 1 and isinstance(args[0], list) and cls in ("Arrays", "Lists", "Stream") else list(args)
        if (cls, name) in (("Collections", "emptyList"), ("Collections", "emptySet")):
            return []
        if (cls, name) in (("Paths", "get"), ("Path", "of"), ("IOUtils", "getPath")):
            path = normalized(jstr(args[0]))
            for part in args[1:]:
                path = join_path(path, part)
            return path
        if cls == "IOUtils" and name in ("createTempFile", "createTempPath"):
            return self.run.temp(args[0], gatk_extension(args[1]))
        if cls in ("IOUtils", "Files") and name in ("createTempDir", "createTempDirectory"):
            return self.run.temp(args[0] if args else "tmp", "", directory=True)
        if cls == "File" and name in ("createTempFile",):
            return self.run.temp(args[0], args[1])
        if cls in ("ArrayUtils", "ObjectArrays") and name in ("addAll", "concat"):
            return [v for part in args if isinstance(part, list) for v in part]
        if name == "values" and target.info is not None and target.info.enum_constants:
            return [EnumValue(c) for c in target.info.enum_constants]
        if cls == "Assert" or name in ("assertThrows", "expectThrows"):
            for arg in args:
                if isinstance(arg, Lambda):
                    self.run_lambda(arg)
            return None
        if cls == "ArgumentsBuilder" and name == "create":
            return Builder()
        if cls == "System" and name == "getProperty":
            if jstr(args[0]) == "user.dir":
                return ROOT
            if len(args) > 1:
                return args[1]
            raise Unresolved(f"system property {jstr(args[0])}")
        if cls == "StringUtils" and name == "join":
            return jstr(args[1] if len(args) > 1 else " ").join(jstr(v) for v in args[0])
        if cls == "Utils" and name == "escapeExpressions":
            return escape_expressions(jstr(args[0]))
        if name == "getSimpleName":
            return cls
        if name == "getName" and target.info is not None:
            return f"{target.info.package}.{cls}"
        if cls in ("CommandLineProgramTest", "GATKBaseTest", "BaseTest") or (target.info and name in target.info.methods):
            native = self.native(name, args)
            if native is not NOT_NATIVE:
                return native
            if target.info is not None:
                for method in target.info.methods.get(name, []):
                    params = method.child_by_field_name("parameters")
                    if len([p for p in kids(params) if p.type == "formal_parameter"]) == len(args) and method not in self.run.stack:
                        return Interpreter(self.run, target.info).call(method, args, target.info)
        raise Unresolved(f"{cls}.{name}")


NOT_NATIVE = object()
FRAMEWORK = {"BaseTest", "GATKBaseTest", "CommandLineProgramTest", "CommandLineProgramTester"}
RUNNERS = {"runCommandLine", "runCommandLineWithoutVerbosity", "instanceMain", "runToolInNewJVM", "executeTest"}
LAZY = set()


def inject_verbosity(tokens):
    """`CommandLineProgramTester.injectDefaultVerbosity`."""
    if any(t.lower() in ("--verbosity", "-verbosity") for t in tokens):
        return tokens
    return tokens + ["--verbosity", "ERROR"]


def truthy(value):
    if isinstance(value, Poison):
        raise Unresolved(value.reason)
    if not isinstance(value, bool):
        raise Unresolved(f"condition on a {type(value).__name__}")
    return value


def walk(node):
    yield node
    for child in node.children:
        yield from walk(child)


def builder_method(builder, name, args):
    flags = {"addInput": "input", "addOutput": "output", "addReference": "reference", "addVCF": "variant",
             "addInterval": "intervals", "addIntervals": "intervals", "addMask": "exclude-intervals"}
    if name in flags:
        return builder.add(flags[name], args[0])
    if name == "add":
        return builder.add(args[0], args[1])
    if name == "addRaw":
        return builder.raw(args[0])
    if name == "addFlag":
        return builder.raw("--" + jstr(args[0]))
    if name == "addBooleanArgument":
        return builder.add(args[0], args[1])
    if name == "addArgument":
        return builder.add(args[0], args[1])
    if name in ("getArgsList", "getArgsArray"):
        return list(builder.tokens)
    if name in ("getString", "toString"):
        return " ".join(builder.tokens)
    if name == "concat":
        builder.tokens += args[0].tokens
        return builder
    if name == "copy":
        return Builder(builder.tokens)
    raise Unresolved(f"ArgumentsBuilder.{name}")


def string_method(text, name, args):
    simple = {
        "toString": lambda: text, "trim": lambda: text.strip(), "getAbsolutePath": lambda: absolute(text),
        "getPath": lambda: text, "toPath": lambda: text, "toFile": lambda: text, "toAbsolutePath": lambda: absolute(text),
        "getCanonicalPath": lambda: absolute(text), "getAbsoluteFile": lambda: absolute(text),
        "getURIString": lambda: text, "getRawInputString": lambda: text, "toUri": lambda: absolute(text),
        "toURI": lambda: absolute(text), "getName": lambda: text.rstrip("/").rsplit("/", 1)[-1],
        "getFileName": lambda: text.rstrip("/").rsplit("/", 1)[-1], "isEmpty": lambda: text == "",
        "getParent": lambda: text.rstrip("/").rsplit("/", 1)[0], "getParentFile": lambda: text.rstrip("/").rsplit("/", 1)[0],
        "toLowerCase": lambda: text.lower(), "toUpperCase": lambda: text.upper(), "length": lambda: len(text),
        "intern": lambda: text, "delete": lambda: True, "deleteOnExit": lambda: None, "mkdirs": lambda: True,
        "mkdir": lambda: True,
    }
    if name in simple and not args:
        return simple[name]()
    if name in ("resolve", "resolveSibling"):
        base = text if name == "resolve" else text.rstrip("/").rsplit("/", 1)[0]
        return join_path(base, args[0])
    if name == "concat":
        return text + jstr(args[0])
    if name == "replace":
        return text.replace(jstr(args[0]), jstr(args[1]))
    if name == "replaceAll":
        return re.sub(jstr(args[0]), jstr(args[1]).replace("$", "\\"), text)
    if name in ("equals", "contentEquals"):
        return text == (args[0] if isinstance(args[0], str) else None)
    if name == "equalsIgnoreCase":
        return isinstance(args[0], str) and text.lower() == args[0].lower()
    if name in ("startsWith", "endsWith", "contains"):
        return {"startsWith": text.startswith, "endsWith": text.endswith, "contains": text.__contains__}[name](jstr(args[0]))
    if name == "split":
        return [p for p in re.split(jstr(args[0]), text)]
    if name == "substring":
        return text[int(args[0]): int(args[1]) if len(args) > 1 else None]
    if name == "formatted":
        return jformat(text, args)
    raise Unresolved(f"String.{name}")


def list_method(items, name, args):
    if name in ("add", "addAll") and len(args) == 1:
        if name == "addAll":
            items.extend(args[0].tokens if isinstance(args[0], Builder) else args[0])
        else:
            items.append(args[0])
        return True
    if name == "add" and len(args) == 2:
        items.insert(int(args[0]), args[1])
        return None
    if name in ("toArray", "iterator", "stream", "toList", "collect", "clone", "asList"):
        return list(items)
    if name == "size":
        return len(items)
    if name == "get":
        return items[int(args[0])]
    if name == "isEmpty":
        return not items
    if name == "contains":
        return any(jstr(v) == jstr(args[0]) for v in items)
    if name == "toString":
        return jstr(items)
    if name == "subList":
        return items[int(args[0]):int(args[1])]
    raise Unresolved(f"List.{name}")


# ---------------------------------------------------------------------------------------------
# Extraction.


def element_strings(node):
    if node is None:
        return []
    if node.type == "string_literal":
        return [java_string(node)]
    return [java_string(n) for n in walk(node) if n.type == "string_literal"]


def test_methods(sources, info):
    """(owner, method, @Test values) of every TestNG test the class runs, inherited ones included."""
    seen, found = set(), []
    class_test = None
    for owner in sources.chain(info):
        if class_test is None and "Test" in owner.annotations:
            class_test = owner.annotations["Test"]
        for name, methods in owner.methods.items():
            for method in methods:
                key = (name, len(kids(method.child_by_field_name("parameters"))))
                if key in seen:
                    continue
                seen.add(key)
                marks = annotations(method)
                if "Test" in marks:
                    found.append((owner, method, marks["Test"]))
                elif class_test is not None and has_modifier(method, "public") and not has_modifier(method, "static") \
                        and not any(m.startswith(("Before", "After", "DataProvider")) for m in marks) \
                        and method.child_by_field_name("type").type == "void_type":
                    found.append((owner, method, {}))
    return found, class_test or {}


def provider_rows(sources, info, test, values, interpreter):
    """The rows a `@DataProvider` returns, or [None] for a test that takes none."""
    if "dataProvider" not in values:
        return [None]
    provider = element_strings(values["dataProvider"])[0]
    owner_info = info
    if "dataProviderClass" in values:
        name = values["dataProviderClass"].text.decode().removesuffix(".class").split(".")[-1]
        owner_info = sources.classes.get(name)
        if owner_info is None:
            raise Unresolved(f"data provider class {name} is outside the sources")
    for owner in sources.chain(owner_info):
        for methods in owner.methods.values():
            for method in methods:
                marks = annotations(method)
                if "DataProvider" not in marks:
                    continue
                named = element_strings(marks["DataProvider"].get("name"))
                if (named[0] if named else method.child_by_field_name("name").text.decode()) == provider:
                    rows = Interpreter(interpreter.run, owner_info).call(method, [], owner)
                    if isinstance(rows, Poison):
                        raise Unresolved(rows.reason)
                    return [row if isinstance(row, list) else [row] for row in rows]
    raise Unresolved(f"data provider {provider} not found")


def setup_methods(sources, info):
    """`@BeforeClass`/`@BeforeMethod`/... methods, superclasses first."""
    found = []
    for owner in reversed(sources.chain(info)):
        for methods in owner.methods.values():
            for method in methods:
                if any(m in annotations(method) for m in ("BeforeClass", "BeforeMethod", "BeforeTest", "BeforeSuite")):
                    found.append((owner, method))
    return found


def extract_class(sources, info, relative):
    invocations, unextracted = [], []
    tests, class_test = test_methods(sources, info)
    for owner, method, values in tests:
        name = method.child_by_field_name("name").text.decode()
        base = f"{info.name}#{name}"
        merged = {**class_test, **values}
        groups = set(element_strings(merged.get("groups")))
        if merged.get("enabled") is not None and merged["enabled"].type == "false":
            unextracted.append({"id": base, "reason": "disabled (enabled = false)"})
            continue
        if groups & CLOUD_GROUPS:
            unextracted.append({"id": base, "reason": f"reads a bucket (groups {sorted(groups & CLOUD_GROUPS)})"})
            continue
        try:
            probe = Run(base, sources)
            rows = provider_rows(sources, info, method, merged, Interpreter(probe, info))
        except (Unresolved, Return, Thrown, RecursionError, ArithmeticError, LookupError, TypeError, ValueError, AttributeError) as error:
            unextracted.append({"id": base, "reason": f"data provider: {error}"})
            continue
        for index, row in enumerate(rows):
            test_id = base if row is None else f"{base}[{index}]"
            run = Run(test_id, sources)
            interpreter = Interpreter(run, info)
            try:
                for setup_owner, setup in setup_methods(sources, info):
                    try:
                        interpreter.call(setup, [], setup_owner)
                    except (Unresolved, Thrown):
                        pass
                interpreter.call(method, row or [], owner)
            except Thrown:
                pass
            except (Unresolved, RecursionError, ArithmeticError, LookupError, TypeError, ValueError, AttributeError) as error:
                unextracted.append({"id": test_id, "reason": str(error) or type(error).__name__})
                continue
            for k, invocation in enumerate(run.invocations):
                invocation_id = test_id if len(run.invocations) == 1 else f"{test_id}.{k}"
                invocations.append({"id": invocation_id, "source": relative, **invocation})
    return invocations, unextracted


def extract(options):
    gatk = Path(options.gatk).resolve()
    roots = [gatk / "src" / "main" / "java", gatk / "src" / "test" / "java", gatk / "src" / "testUtils" / "java"]
    sources = Sources([r for r in roots if r.is_dir()] + [Path(c) for c in options.constants])
    invocations, unextracted = [], []
    files = sorted((gatk / "src" / "test" / "java").rglob("*IntegrationTest.java"))
    for path in files:
        info = sources.classes.get(path.stem)
        if info is None or info.abstract or info.path != path:
            continue
        found, missed = extract_class(sources, info, path.relative_to(gatk).as_posix())
        invocations += found
        unextracted += missed
    head = sha_of(gatk)
    payload = {
        "$comment": "Generated by tools/coverage/test_resources.py extract from the pinned GATK sources; do not edit.",
        "source": {"repository": "broadinstitute/gatk", "tag": options.tag, "sha": head, "classes": len(files)},
        "invocations": invocations,
        "unextracted": unextracted,
    }
    DATA.mkdir(exist_ok=True)
    INVOCATIONS.write_text(json.dumps(payload, indent=1) + "\n")
    tools = {i["tool"] for i in invocations}
    print(f"classes={len(files)} invocations={len(invocations)} tools={len(tools)} unextracted={len(unextracted)}")
    return 0


def sha_of(checkout):
    head = (checkout / ".git" / "HEAD").read_text().strip() if (checkout / ".git" / "HEAD").is_file() else ""
    if head.startswith("ref: "):
        ref = checkout / ".git" / head[5:]
        return ref.read_text().strip() if ref.is_file() else head
    return head


# ---------------------------------------------------------------------------------------------
# Replay.


def ported_tools():
    """The tools whose port lives here: the ones a gatk-rs conformance suite names."""
    manifest = json.loads(MANIFEST.read_text())
    return {tool for suite in manifest["suites"] for tool in suite.get("tools", [])}


def exclusions():
    return json.loads(EXCLUSIONS.read_text()).get("invocations", {}) if EXCLUSIONS.is_file() else {}


def setup_command(invocation):
    words = [f"mkdir -p {d}" for d in invocation["mkdir"]] + [f"touch {f}" for f in invocation["touch"]]
    return "".join(f"{w} && " for w in words)


def replay(invocation, port, gatk):
    with tempfile.TemporaryDirectory(prefix="gatk-test-resources-", ignore_cleanup_errors=True) as tmp:
        workdir = Path(tmp)
        (workdir / "tmp").mkdir()
        try:
            setup = setup_command(invocation)
            code, text, error = run_array.run_oracle(invocation["tool"], [], workdir, invocation["args"], setup=setup)
            reference = run_array.outcome(code, text, error)
            port_code, port_text, port_error = run_array.run_port(port, invocation["tool"], [], workdir, invocation["args"], setup=setup)
            ours = run_array.outcome(port_code, port_text, port_error)
        finally:
            run_array.empty_output_in_container(workdir)
    result = {"id": invocation["id"], "tool": invocation["tool"], "reference_exit": code,
              "matched": ours == reference, "reference": reference[:4000], "port": ours[:4000]}
    if ours != reference:
        result["first_difference"] = first_difference(reference, ours)
    return result


def first_difference(reference, port):
    """Where two answers part, since a head kept to 4000 characters may not reach it."""
    left, right = reference.split("\n"), port.split("\n")
    for line, (a, b) in enumerate(zip(left, right)):
        if a != b:
            return {"line": line, "reference": a[:2000], "port": b[:2000]}
    return {"line": min(len(left), len(right)), "reference": f"{len(left)} lines", "port": f"{len(right)} lines"}


def chosen(options):
    """The invocations this shard replays: a ported tool's, not excluded, every `count`-th."""
    data = json.loads(INVOCATIONS.read_text())
    ported, excluded = ported_tools(), exclusions()
    picked = [i for i in data["invocations"] if i["tool"] in ported and i["id"] not in excluded
              and (not options.tool or i["tool"] in options.tool)]
    index, count = (int(x) for x in options.shard.split("/"))
    return picked[index::count]


def run(options):
    run_array.FIXTURES_DIR = str(Path(options.gatk).resolve())
    picked = chosen(options)
    # The R-plot tools' port runs in the image `run_array` builds for them, as their arrays do.
    if any(run_array.past_the_interpreter(i["tool"]) == "Rscript" for i in picked):
        run_array.build_port_with_r()
    results = []
    with concurrent.futures.ThreadPoolExecutor(options.jobs) as pool:
        for result in pool.map(lambda i: replay(i, options.port, options.gatk), picked):
            print(f"{result['id']}: {'match' if result['matched'] else 'DIFFERS'} (reference exit {result['reference_exit']})", flush=True)
            results.append(result)
    Path(options.out).write_text(json.dumps({"shard": options.shard, "results": results}, indent=1) + "\n")
    print(f"replayed={len(results)} matched={sum(r['matched'] for r in results)}")
    return 0


def merge(options):
    """Fold the shards' results into `measured.json`, the file VALIDATION.md reads."""
    data = json.loads(INVOCATIONS.read_text())
    ported, excluded = ported_tools(), exclusions()
    results = [r for path in options.shards for r in json.loads(Path(path).read_text())["results"]]
    tools = {}

    def entry(tool):
        return tools.setdefault(tool, {"extracted": 0, "replayed": 0, "matched": 0, "excluded": 0,
                                       "unextracted": 0, "diverged": []})

    by_class = {}
    for invocation in data["invocations"]:
        if invocation["tool"] in ported:
            entry(invocation["tool"])["extracted"] += 1
            entry(invocation["tool"])["excluded"] += invocation["id"] in excluded
            by_class.setdefault(invocation["id"].split("#")[0], invocation["tool"])
    for missed in data["unextracted"]:
        tool = by_class.get(missed["id"].split("#")[0])
        if tool:
            entry(tool)["unextracted"] += 1
    for result in results:
        found = entry(result["tool"])
        found["replayed"] += 1
        found["matched"] += result["matched"]
        if not result["matched"]:
            found["diverged"].append(result["id"])
    for found in tools.values():
        found["diverged"].sort()
    MEASURED.write_text(json.dumps({"source": data["source"], "tools": dict(sorted(tools.items()))}, indent=1) + "\n")
    print(f"tools={len(tools)} replayed={len(results)} matched={sum(r['matched'] for r in results)}")
    return 0


def lfs_include(options):
    """The LFS paths the replayed invocations name, with the files beside them that share a stem
    (an index, a dictionary), as one `git lfs pull --include` list."""
    patterns = set()
    for invocation in chosen(options):
        for token in invocation["args"]:
            for path in re.findall(rf"{ROOT}/(src/test/resources/large/[^\s,:;]+)", token):
                stem = path.rstrip("/").split(".")[0] if "/" in path else path
                patterns.update({path, f"{path}/**", f"{path}.*", f"{stem}.*"})
    print(",".join(sorted(patterns)))
    return 0


def self_check():
    assert jdouble(0.5) == "0.5" and jdouble(1e-5) == "1.0E-5" and jdouble(3.0) == "3.0" and jdouble(1.5e7) == "1.5E7"
    assert jformat("-O %s -x %d%%", ["/a", 3]) == "-O /a -x 3%"
    assert escape_expressions("-a 1 -select 'QD > 2.0' -b") == ["-a", "1", "-select", "QD > 2.0", "-b"]
    assert first_difference("a\nb", "a\nc") == {"line": 1, "reference": "b", "port": "c"}
    assert inject_verbosity(["-I", "x"]) == ["-I", "x", "--verbosity", "ERROR"]
    assert join_path("/a/b/", "c/") == "/a/b/c" and normalized("x//y/") == "x/y"
    assert Builder().add("O", "%s").tokens == ["--O", "%s"]
    assert gatk_extension("fasta") == ".fasta" and gatk_extension(".vcf") == ".vcf"


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    one = commands.add_parser("extract")
    one.add_argument("--gatk", required=True, help="the pinned GATK checkout")
    one.add_argument("--constants", action="append", default=[],
                     help="another source tree to read constants from (htsjdk's src/main/java)")
    one.add_argument("--tag", default="4.7.0.0", help="the release the checkout is at, which CI clones")
    two = commands.add_parser("run")
    two.add_argument("--gatk", required=True)
    two.add_argument("--port", required=True)
    two.add_argument("--shard", default="0/1")
    two.add_argument("--jobs", type=int, default=2)
    two.add_argument("--tool", action="append")
    two.add_argument("--out", default="test-resources-shard.json")
    three = commands.add_parser("merge")
    three.add_argument("shards", nargs="+")
    four = commands.add_parser("lfs-include")
    four.add_argument("--shard", default="0/1")
    four.add_argument("--tool", action="append")
    options = parser.parse_args(argv)
    self_check()
    return {"extract": extract, "run": run, "merge": merge, "lfs-include": lfs_include}[options.command](options)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
