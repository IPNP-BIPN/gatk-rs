//! `AlleleFrequencyQC`: the chi-squared statistic over the allele-frequency bins, and the p-value
//! it gives.
//!
//! The traversal is `VariantEval`'s with every knob preset, and is not ported here. What is ported
//! is what the tool does with the report that traversal writes: the grouping, the statistic, the
//! distribution it is read against, and the two metrics it writes.
//!
//! Ported from `org.broadinstitute.hellbender.tools.walkers.varianteval.AlleleFrequencyQC` and
//! `org.apache.commons.math3.distribution.ChiSquaredDistribution`.

/// `threshold`, under which the tool logs an error and changes nothing else.
pub const DEFAULT_THRESHOLD: f64 = 0.05;
/// `allowedVariance`, which stands in for the expected count Pearson's statistic would divide by.
pub const DEFAULT_ALLOWED_VARIANCE: f64 = 0.01;

/// The `METRIC_TYPE` column, which is a constant.
pub const METRIC_TYPE: &str = "Allele Frequency";

/// The `Filter` value the rows are cut down to before anything is grouped.
pub const CALLED: &str = "called";

/// `calculateChiSquaredStatistic`.
///
/// Not Pearson's: the expected count in the denominator is replaced by a constant variance, and
/// that variance is SQUARED, so a variance ten times larger divides the statistic by a hundred.
/// A bin holding fewer than two entries contributes nothing; on the reference's own path every bin
/// holds exactly two, one per eval track, so that guard never fires.
pub fn chi_squared_statistic(bins: &[Vec<f64>], variance: f64) -> f64 {
    let sum: f64 = bins
        .iter()
        .map(|entries| {
            if entries.len() >= 2 {
                (entries[0] - entries[1]).powi(2)
            } else {
                0.0
            }
        })
        .sum();
    sum / variance.powi(2)
}

/// The degrees of freedom: the bin count less one, counted over the bins the report holds and not
/// over the ones the data reached.
///
/// The allele-frequency stratifier emits a fixed ladder, so this is the same number whatever the
/// file holds, and a bin no variant reached contributes a term of nought to the statistic while
/// still being counted here.
pub fn degrees_of_freedom(bins: usize) -> f64 {
    bins as f64 - 1.0
}

/// `1 - ChiSquaredDistribution(df).cumulativeProbability(x)`: the upper tail.
///
/// commons-math evaluates the gamma distribution's CDF as the regularized lower incomplete gamma
/// function, which htsjdk-rs already carries, so the p-value here is that function's and not an
/// approximation of it.
pub fn p_value(statistic: f64, degrees_of_freedom: f64) -> f64 {
    if statistic <= 0.0 {
        return 1.0;
    }
    let cumulative = jmath::gamma::regularized_gamma_p(
        degrees_of_freedom / 2.0,
        statistic / 2.0,
        1e-14,
        i32::MAX,
    )
    .unwrap_or(f64::NAN);
    1.0 - cumulative
}

/// The two numbers the metrics file carries, from the bins the report gave.
pub fn metrics(bins: &[Vec<f64>], variance: f64) -> (f64, f64) {
    let statistic = chi_squared_statistic(bins, variance);
    (
        statistic,
        p_value(statistic, degrees_of_freedom(bins.len())),
    )
}

/// Whether the run complains, which is the whole of what the threshold decides.
pub fn complains(p_value: f64, threshold: f64) -> bool {
    p_value < threshold
}

/// The message that complaint carries.
pub fn complaint(p_value: f64) -> String {
    format!(
        "Allele frequencies between your array VCF and the expected VCF do not match with a \
         significant pvalue of {p_value}"
    )
}

/// The module whose table the tool reads back, which is the only module it runs.
pub const MODULE: &str = "VariantAFEvaluator";

/// The header line the sample name is taken from, and not a genotype column.
pub const SAMPLE_ALIAS_KEY: &str = "sampleAlias";

/// The `MetricBase` subclass the metrics file names.
pub const METRIC_CLASS: &str =
    "org.broadinstitute.hellbender.metrics.analysis.AlleleFrequencyQCMetric";

/// What `getOtherHeaderLine("sampleAlias").getValue()` raises when the header has no such line,
/// as HotSpot's helpful message words it.
pub const NO_ALIAS_MESSAGE: &str = "Cannot invoke \"htsjdk.variant.vcf.VCFHeaderLine.getValue()\" \
     because the return value of \"htsjdk.variant.vcf.VCFHeader.getOtherHeaderLine(String)\" is null";

/// What reading the report back refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReportError {
    /// The report holds no table of that name: `GATKReport.getTable` refuses by name.
    NoTable(String),
    /// A column the grouping reads is missing, so `indexOf` answers -1 and the row lookup fails.
    NoColumn(String),
    /// A value of the averaged column does not parse as a double.
    NotANumber(String),
}

/// The bins `onTraversalSuccess` groups, read out of the report text the traversal wrote.
///
/// `new GATKReport(outFile).getTable(module)` parses each column by the type its format declares,
/// so `avgVarAF` (`%.8f`) comes back as the `Double` its eight decimals spell and not as the
/// double the evaluator computed, and `AlleleFrequency` (`%s`) as a string. The rows are cut to
/// `Filter == called`, then grouped by that string with `Collectors.groupingBy`, whose map is a
/// `HashMap`: the bins are returned in that map's iteration order, each holding its averages in
/// row order, which is the order the statistic sums them in.
pub fn bins_from_report(report: &str, module: &str) -> Result<Vec<Vec<f64>>, ReportError> {
    let mut lines = report.lines();
    let name_line = format!("#:GATKTable:{module}:");
    let exact = format!("#:GATKTable:{module}");
    lines
        .by_ref()
        .find(|line| line.starts_with(&name_line) || *line == exact)
        .ok_or_else(|| ReportError::NoTable(module.to_string()))?;
    let header: Vec<&str> = lines.next().unwrap_or("").split_whitespace().collect();
    let index = |name: &str| {
        header
            .iter()
            .position(|column| *column == name)
            .ok_or_else(|| ReportError::NoColumn(name.to_string()))
    };
    let (frequency, filter, average) = (
        index("AlleleFrequency")?,
        index("Filter")?,
        index("avgVarAF")?,
    );
    let mut keys: Vec<String> = Vec::new();
    let mut grouped: Vec<Vec<f64>> = Vec::new();
    for line in lines.take_while(|line| !line.trim().is_empty()) {
        let columns: Vec<&str> = line.split_whitespace().collect();
        if columns.get(filter) != Some(&CALLED) {
            continue;
        }
        let key = columns.get(frequency).copied().unwrap_or("").to_string();
        let text = columns.get(average).copied().unwrap_or("");
        let value: f64 = text
            .parse()
            .map_err(|_| ReportError::NotANumber(text.to_string()))?;
        match keys.iter().position(|other| *other == key) {
            Some(slot) => grouped[slot].push(value),
            None => {
                keys.push(key);
                grouped.push(vec![value]);
            }
        }
    }
    let order = gatk_engine::java_hash::hash_set_order(&keys)
        .map_err(|error| ReportError::NoColumn(format!("{error:?}")))?;
    Ok(order
        .iter()
        .map(|key| {
            let slot = keys.iter().position(|other| other == key).unwrap_or(0);
            grouped[slot].clone()
        })
        .collect())
}

/// `calculateChiSquaredStatistic` as the JVM evaluates it: the squares summed by
/// `DoubleStream.sum`, which is compensated, over the bins in the map's order, then divided by the
/// variance squared. `Math.pow(x, 2.)` is `x * x` exactly.
pub fn java_chi_squared_statistic(bins: &[Vec<f64>], variance: f64) -> f64 {
    let squares: Vec<f64> = bins
        .iter()
        .map(|entries| {
            if entries.len() >= 2 {
                let difference = entries[0] - entries[1];
                difference * difference
            } else {
                0.0
            }
        })
        .collect();
    gatk_engine::allele_fraction_cluster::double_stream_sum(&squares) / (variance * variance)
}

/// The metric `onTraversalSuccess` writes.
pub struct AlleleFrequencyQcMetric {
    pub sample: String,
    pub metric_value: f64,
    pub chi_sq_value: f64,
}

impl htsjdk_metrics::file::MetricBean for AlleleFrequencyQcMetric {
    fn class_name(&self) -> &str {
        METRIC_CLASS
    }

    fn columns(&self) -> &[&'static str] {
        &["SAMPLE", "METRIC_TYPE", "METRIC_VALUE", "CHI_SQ_VALUE"]
    }

    fn values(&self) -> Vec<htsjdk_metrics::file::Value> {
        use htsjdk_metrics::file::Value;
        vec![
            Value::Str(self.sample.clone()),
            Value::Str(METRIC_TYPE.to_string()),
            Value::Double(self.metric_value),
            Value::Double(self.chi_sq_value),
        ]
    }
}

/// The whole metrics file: `MetricsUtils.saveMetrics` over a `MetricsFile` with no header and one
/// metric, and no histogram.
pub fn metrics_file(sample: &str, bins: &[Vec<f64>], variance: f64) -> String {
    let statistic = java_chi_squared_statistic(bins, variance);
    let metric = AlleleFrequencyQcMetric {
        sample: sample.to_string(),
        metric_value: p_value(statistic, degrees_of_freedom(bins.len())),
        chi_sq_value: statistic,
    };
    let mut file = htsjdk_metrics::file::MetricsFile::new();
    file.add_metric(&metric);
    file.write()
}
