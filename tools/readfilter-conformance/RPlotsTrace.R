# Run ahead of the reference's own, unmodified plotting scripts by RPlotsDump.java, inside the image
# that has R. It records what they hand to the graphics calls and then lets the call happen.
#
# Those scripts render a PNG and drop everything else, so this is the only way to see their
# arithmetic: every number on these rows is the one the reference's base graphics or ggplot2 was
# given, after the log2-to-linear transform, the midpoint rounding, the genomic offsetting and the
# running index. RPlotsReadBack.R writes the same rows from the port's exported .data.rds.
#
# The masks work because `source` evaluates in the global environment, which is searched ahead of
# the attached packages. Numbers are written at fifteen significant digits, which is what
# `write.table` would keep and which is the most a libm's last bit cannot move.

TraceFile = Sys.getenv("R_PLOTS_TRACE_FILE")
TraceNumber = function(x) ifelse(is.na(x) & !is.nan(x), "NA", ifelse(is.nan(x), "NaN", sprintf("%.15g", x)))
TraceRows = function(kind, ...) {
    cat(paste(kind, ..., sep="\t"), sep="\n", file=TraceFile, append=TRUE)
}

plot = function(x, ..., xlim, ylim, main) {
    TraceRows("plot", main, TraceNumber(xlim[2]), TraceNumber(ylim[2]))
    graphics::plot(x, ..., xlim=xlim, ylim=ylim, main=main)
}
rect = function(xleft, ybottom, xright, ytop, col=NA, ...) {
    TraceRows("rect", TraceNumber(xleft), TraceNumber(ybottom), TraceNumber(xright), TraceNumber(ytop), col)
    graphics::rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col=col, ...)
}
points = function(x, y, col, ...) {
    if (length(x) > 0) {
        TraceRows("points", TraceNumber(x), TraceNumber(y), col)
    }
    graphics::points(x=x, y=y, col=col, ...)
}
segments = function(x0, y0, x1, y1, ...) {
    TraceRows("segments", TraceNumber(x0), TraceNumber(y0), TraceNumber(x1), TraceNumber(y1))
    graphics::segments(x0=x0, y0=y0, x1=x1, y1=y1, ...)
}

# plotAlleleFrequencyQC.R is already ggplot2, and hands the same merged frame to both of its plots.
TraceFrameSeen = FALSE
ggplot = function(data, ...) {
    if (!TraceFrameSeen) {
        TraceFrameSeen <<- TRUE
        if (nrow(data) > 0) {
            TraceRows("frame", TraceNumber(data[["AF_bin"]]), TraceNumber(data[["avgVarAF.array"]]),
                      TraceNumber(data[["avgVarAF.thousand_g"]]), data[["sample"]])
        }
    }
    ggplot2::ggplot(data, ...)
}
