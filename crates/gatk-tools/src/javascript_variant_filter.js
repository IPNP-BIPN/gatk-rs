// The Java objects htsjdk's JavascriptVariantFilter binds, as Nashorn shows them to a script.
//
// Evaluated once per run. It returns the two functions that turn the JSON the Rust side writes
// (javascript_variant_filter.rs) into `variant` and `header`; nothing here reaches the global
// object except the String methods at the end. Where a stand-in cannot answer as the Java object
// would (a method it does not carry, the toString of a VariantContext, an exception whose wording
// is Java's), it calls `refuse`, which the Rust side checks after the script whatever the script
// did with the error, so a try/catch cannot turn a refusal into a different answer.
(function (refuse) {
  var DATA = Symbol("data");
  function hidden(object, key, value) {
    Object.defineProperty(object, key, { value: value, enumerable: false });
  }
  function fail(what) {
    refuse(what);
    throw new Error(what);
  }
  function unmirrored(what) {
    return function () {
      fail(what + " is not mirrored by this port");
    };
  }
  // The end of every stand-in's prototype chain: a name none of them carries is a Java member
  // this port does not mirror. ToPrimitive's and JSON's own lookups pass through.
  var MISSING = new Proxy(Object.create(null), {
    get: function (target, key) {
      if (typeof key === "symbol" || key === "valueOf" || key === "toJSON") return undefined;
      fail("the Java member " + key + " is not mirrored by this port");
    },
  });
  // Nashorn's bean properties: `variant.start` is getStart(), `variant.SNP` isSNP(), with
  // java.beans' decapitalisation, which leaves `ID` and `NAlleles` as they are.
  function beans(proto) {
    Object.getOwnPropertyNames(proto).forEach(function (name) {
      var match = /^(get|is)([A-Z].*)$/.exec(name);
      if (!match || proto[name].length !== 0) return;
      var rest = match[2];
      var property = /^.[A-Z]/.test(rest) ? rest : rest[0].toLowerCase() + rest.slice(1);
      Object.defineProperty(proto, property, {
        get: function () {
          return this[name]();
        },
      });
    });
    return proto;
  }
  // Nashorn's Java objects are not JavaScript objects (`instanceof Object` is false, and they
  // have no Object.prototype methods), so neither are these.
  function define(methods) {
    var proto = Object.create(MISSING);
    Object.keys(methods).forEach(function (name) {
      hidden(proto, name, methods[name]);
    });
    return proto;
  }
  function javaEquals(left, right) {
    return left !== null && typeof left === "object" && "equals" in left
      ? left.equals(right)
      : left === right;
  }
  function index(i, length) {
    if (!(i >= 0 && i < length)) fail("an IndexOutOfBoundsException, whose wording is Java's");
    return i;
  }

  // java.util.List: indexable, `length`, for-in over the indices, and no Array methods.
  var LIST = beans(define({
    size: function () { return this.length; },
    get: function (i) { return this[index(i, this.length)]; },
    isEmpty: function () { return this.length === 0; },
    contains: function (o) { return this.indexOf(o) >= 0; },
    indexOf: function (o) {
      for (var i = 0; i < this.length; i++) if (javaEquals(this[i], o)) return i;
      return -1;
    },
    forEach: function (action) { for (var i = 0; i < this.length; i++) action(this[i]); },
    // AbstractCollection.toString.
    toString: function () {
      var parts = [];
      for (var i = 0; i < this.length; i++) parts.push(String(this[i]));
      return "[" + parts.join(", ") + "]";
    },
  }));
  function list(values, proto) {
    var object = Object.create(proto || LIST);
    values.forEach(function (value, i) { object[i] = value; });
    hidden(object, "length", values.length);
    return object;
  }
  // java.util.Set: `length` and the Set methods, but neither indices nor for-in keys. The order is
  // the HashSet's, which the Rust side worked out.
  var SET = define({
    size: function () { return this.length; },
    isEmpty: function () { return this.length === 0; },
    contains: function (o) { return LIST.indexOf.call(d(this), o) >= 0; },
    forEach: function (action) { d(this).forEach(function (v) { action(v); }); },
    toString: function () { return LIST.toString.call(d(this)); },
  });
  function set(values) {
    var object = host(SET, list(values));
    hidden(object, "length", values.length);
    return object;
  }
  // A Java int[]: indexable and `length`, and a toString that is an identity hash.
  var ARRAY = define({ toString: unmirrored("The toString of a Java array") });
  function array(values) {
    return values === null ? null : list(values, ARRAY);
  }

  function enumeration(names) {
    var constants = {};
    names.forEach(function (name, ordinal) {
      constants[name] = Object.create(define({
        name: function () { return name; },
        ordinal: function () { return ordinal; },
        toString: function () { return name; },
        equals: function (o) { return this === o; },
      }));
    });
    return constants;
  }
  var TYPE = enumeration(["NO_VARIATION", "SNP", "MNP", "INDEL", "SYMBOLIC", "MIXED"]);
  var GENOTYPE_TYPE = enumeration(["NO_CALL", "HOM_REF", "HET", "HOM_VAR", "UNAVAILABLE", "MIXED"]);

  function host(proto, data) {
    var object = Object.create(proto);
    hidden(object, DATA, data);
    return object;
  }
  function d(object) {
    return object[DATA];
  }

  var ALLELE = beans(define({
    getBaseString: function () { var a = d(this); return a.noCall ? "." : a.symbolic ? "" : a.display; },
    getDisplayString: function () { return d(this).display; },
    isReference: function () { return d(this).ref; },
    isNonReference: function () { return !d(this).ref; },
    isSymbolic: function () { return d(this).symbolic; },
    isNoCall: function () { return d(this).noCall; },
    isCalled: function () { return !d(this).noCall; },
    length: function () { return d(this).symbolic ? 0 : d(this).display.length; },
    equals: function (o, ignoreRefState) {
      if (this === o) return true;
      if (o === null || typeof o !== "object" || Object.getPrototypeOf(o) !== ALLELE) return false;
      var a = d(this), b = d(o);
      return (ignoreRefState === true || a.ref === b.ref) && a.noCall === b.noCall && a.display === b.display;
    },
    toString: function () { var a = d(this); return (a.noCall ? "." : a.display) + (a.ref ? "*" : ""); },
  }));
  var NO_CALL = host(ALLELE, { display: "", ref: false, symbolic: false, noCall: true });

  // An attribute's answers were worked out on the Rust side, with htsjdk-rs's own accessors:
  // `int` is a number, "default", or null where the Java accessor throws, and so on.
  function attribute(entries, key) {
    for (var i = 0; i < entries.length; i++) if (entries[i].key === key) return entries[i];
    return null;
  }
  function value(v) {
    return Array.isArray(v) ? list(v) : v;
  }
  function typed(entry, answer, fallback, accessor) {
    if (entry === null) return fallback;
    var result = entry[answer];
    if (result === "default") return fallback;
    if (result === null) fail(accessor + " throws on " + entry.key + ", in Java's words");
    return answer === "double" ? Number(result) : result;
  }

  var GENOTYPE = beans(define({
    getSampleName: function () { return d(this).name; },
    getAlleles: function () { return d(this).alleleList; },
    getAllele: function (i) { return d(this).alleleList.get(i); },
    getPloidy: function () { return d(this).alleleList.length; },
    isPhased: function () { return d(this).phased; },
    getType: function () { return GENOTYPE_TYPE[d(this).type]; },
    isHom: function () { return this.isHomRef() || this.isHomVar(); },
    isHomRef: function () { return d(this).type === "HOM_REF"; },
    isHomVar: function () { return d(this).type === "HOM_VAR"; },
    isHet: function () { return d(this).type === "HET"; },
    isHetNonRef: function () {
      return this.isHet() && this.getAllele(0).isNonReference() && this.getAllele(1).isNonReference();
    },
    isNoCall: function () { return d(this).type === "NO_CALL"; },
    isCalled: function () { var t = d(this).type; return t !== "NO_CALL" && t !== "UNAVAILABLE"; },
    isMixed: function () { return d(this).type === "MIXED"; },
    isAvailable: function () { return d(this).type !== "UNAVAILABLE"; },
    getDP: function () { return d(this).dp; },
    hasDP: function () { return d(this).dp !== -1; },
    getGQ: function () { return d(this).gq; },
    hasGQ: function () { return d(this).gq !== -1; },
    getAD: function () { return d(this).adArray; },
    hasAD: function () { return d(this).adArray !== null; },
    getPL: function () { return d(this).plArray; },
    hasPL: function () { return d(this).plArray !== null; },
    getFilters: function () { return d(this).ft; },
    isFiltered: function () { return d(this).ft !== null; },
    hasExtendedAttribute: function (key) { return attribute(d(this).extended, key) !== null; },
    getExtendedAttribute: function (key, fallback) {
      var entry = attribute(d(this).extended, key);
      return entry === null ? (fallback === undefined ? null : fallback) : value(entry.value);
    },
    toString: unmirrored("Genotype.toString"),
  }));

  var GENOTYPES = beans(define({
    get: function (which) {
      // `sampleNameToOffset`, where a name given twice maps to its last genotype.
      if (typeof which === "string") {
        var genotype = d(this).byName[which];
        return genotype === undefined ? null : genotype;
      }
      return LIST.get.call(this, which);
    },
    containsSample: function (name) { return this.get(name) !== null; },
    getSampleNamesOrderedByName: function () { return list(d(this).orderedByName); },
    toString: unmirrored("GenotypesContext.toString"),
  }));
  Object.setPrototypeOf(GENOTYPES, LIST);

  var VARIANT = beans(define({
    getContig: function () { return d(this).contig; },
    getStart: function () { return d(this).start; },
    getEnd: function () { return d(this).end; },
    getID: function () { return d(this).id; },
    hasID: function () { return d(this).id !== "."; },
    getReference: function () { return d(this).alleleList.get(0); },
    getAlleles: function () { return d(this).alleleList; },
    getAlternateAlleles: function () {
      return list(Array.prototype.slice.call(d(this).alleleList, 1));
    },
    getAlternateAllele: function (i) { return d(this).alleleList.get(i + 1); },
    getNAlleles: function () { return d(this).alleleList.length; },
    isBiallelic: function () { return d(this).alleleList.length === 2; },
    getType: function () { return TYPE[d(this).type]; },
    isSNP: function () { return d(this).type === "SNP"; },
    isMNP: function () { return d(this).type === "MNP"; },
    isIndel: function () { return d(this).type === "INDEL"; },
    isSymbolic: function () { return d(this).type === "SYMBOLIC"; },
    isMixed: function () { return d(this).type === "MIXED"; },
    isVariant: function () { return d(this).type !== "NO_VARIATION"; },
    getPhredScaledQual: function () { return Number(d(this).qual); },
    hasLog10PError: function () { return d(this).hasQual; },
    getFilters: function () { return set(d(this).filters || []); },
    isFiltered: function () { var f = d(this).filters; return f !== null && f.length > 0; },
    isNotFiltered: function () { return !this.isFiltered(); },
    filtersWereApplied: function () { return d(this).filters !== null; },
    hasAttribute: function (key) { return attribute(d(this).attributes, key) !== null; },
    getAttribute: function (key, fallback) {
      var entry = attribute(d(this).attributes, key);
      return entry === null ? (fallback === undefined ? null : fallback) : value(entry.value);
    },
    getAttributeAsInt: function (key, fallback) {
      return typed(attribute(d(this).attributes, key), "int", fallback, "getAttributeAsInt");
    },
    getAttributeAsDouble: function (key, fallback) {
      return typed(attribute(d(this).attributes, key), "double", fallback, "getAttributeAsDouble");
    },
    getAttributeAsString: function (key, fallback) {
      return typed(attribute(d(this).attributes, key), "string", fallback, "getAttributeAsString");
    },
    getAttributeAsBoolean: function (key, fallback) {
      return typed(attribute(d(this).attributes, key), "boolean", fallback, "getAttributeAsBoolean");
    },
    getGenotypes: function () { return d(this).genotypeList; },
    getGenotype: function (which) {
      var genotypes = d(this).genotypeList;
      return typeof which === "number" && which >= genotypes.length ? null : genotypes.get(which);
    },
    hasGenotype: function (name) { return d(this).genotypeList.containsSample(name); },
    hasGenotypes: function () { return d(this).genotypeList.length > 0; },
    getNSamples: function () { return d(this).genotypeList.length; },
    getSampleNamesOrderedByName: function () { return d(this).genotypeList.getSampleNamesOrderedByName(); },
    getNoCallCount: function () { return count(this, "NO_CALL"); },
    getHomRefCount: function () { return count(this, "HOM_REF"); },
    getHetCount: function () { return count(this, "HET"); },
    getHomVarCount: function () { return count(this, "HOM_VAR"); },
    getMixedCount: function () { return count(this, "MIXED"); },
    toString: unmirrored("VariantContext.toString"),
  }));
  function count(variant, type) {
    return d(variant).genotypes.filter(function (g) { return g.type === type; }).length;
  }

  var HEADER = beans(define({
    getGenotypeSamples: function () { return d(this).sampleList; },
    getNGenotypeSamples: function () { return d(this).samples.length; },
    hasGenotypingData: function () { return d(this).samples.length > 0; },
    hasInfoLine: function (id) { return d(this).info.indexOf(id) >= 0; },
    hasFormatLine: function (id) { return d(this).format.indexOf(id) >= 0; },
    hasFilterLine: function (id) { return d(this).filter.indexOf(id) >= 0; },
    toString: unmirrored("VCFHeader.toString"),
  }));

  // Nashorn hands a script java.lang.String itself, so String's own methods answer on it.
  hidden(String.prototype, "equals", function (o) { return typeof o === "string" && o === String(this); });
  hidden(String.prototype, "equalsIgnoreCase", function (o) {
    return typeof o === "string" && o.toUpperCase() === String(this).toUpperCase();
  });
  hidden(String.prototype, "isEmpty", function () { return this.length === 0; });
  hidden(String.prototype, "contains", function (o) { return String(this).indexOf(o) >= 0; });

  return {
    variant: function (data) {
      var alleles = data.alleles.map(function (a) { return host(ALLELE, a); });
      data.alleleList = list(alleles);
      var genotypes = data.genotypes.map(function (g) {
        g.alleleList = list(g.alleles.map(function (i) { return i < 0 ? NO_CALL : alleles[i]; }));
        g.adArray = array(g.ad);
        g.plArray = array(g.pl);
        return host(GENOTYPE, g);
      });
      data.genotypeList = list(genotypes, GENOTYPES);
      var byName = Object.create(null);
      genotypes.forEach(function (g) { byName[g.getSampleName()] = g; });
      hidden(data.genotypeList, DATA, { orderedByName: data.orderedByName, byName: byName });
      return host(VARIANT, data);
    },
    header: function (data) {
      data.sampleList = list(data.samples);
      return host(HEADER, data);
    },
  };
})
