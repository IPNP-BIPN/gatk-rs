# Reads back a .data.rds the port exported and writes it, in the order the reference made its calls,
# as the rows RPlotsTrace.R records from the reference: the exported frame as a TSV, which is what
# makes the arithmetic of these tools comparable at all.
#
#     Rscript RPlotsReadBack.R <kind: denoised, modeled or qc> <.data.rds> <output>

args = commandArgs(TRUE)
kind = args[1]
frame = readRDS(args[2])
TraceFile = args[3]
invisible(file.create(TraceFile))

TraceNumber = function(x) ifelse(is.na(x) & !is.nan(x), "NA", ifelse(is.nan(x), "NaN", sprintf("%.15g", x)))
TraceRows = function(kind, ...) {
    cat(paste(kind, ..., sep="\t"), sep="\n", file=TraceFile, append=TRUE)
}

if (kind == "qc") {
    if (nrow(frame) > 0) {
        TraceRows("frame", TraceNumber(frame[["AF_bin"]]), TraceNumber(frame[["avgVarAF.array"]]),
                  TraceNumber(frame[["avgVarAF.thousand_g"]]), frame[["sample"]])
    }
    quit(save="no", status=0)
}

panels = attr(frame, "panels")
contigs = attr(frame, "contigs")
for (p in seq_len(nrow(panels))) {
    panel = panels[p, ]
    TraceRows("plot", panel[["MAIN"]], TraceNumber(panel[["X_MAX"]]), TraceNumber(panel[["Y_MAX"]]))
    for (i in seq_len(nrow(contigs))) {
        TraceRows("rect", TraceNumber(contigs[i, "START"]), TraceNumber(0), TraceNumber(contigs[i, "END"]),
                  TraceNumber(panel[["Y_MAX"]]), contigs[i, "FILL"])
    }
    rows = frame[frame[["PANEL"]] == panel[["PANEL"]], ]
    if (kind == "denoised") {
        if (nrow(rows) > 0) {
            TraceRows("points", TraceNumber(rows[["GENOMIC_COORDINATE"]]), TraceNumber(rows[["COPY_RATIO"]]), rows[["COLOR"]])
        }
        next
    }
    for (s in unique(rows[["SEGMENT"]])) {
        piece = rows[rows[["SEGMENT"]] == s, ]
        TraceRows("points", TraceNumber(piece[["GENOMIC_COORDINATE"]]), TraceNumber(piece[["VALUE"]]), piece[["COLOR"]])
        first = piece[1, ]
        x0 = TraceNumber(first[["SEGMENT_START"]])
        x1 = TraceNumber(first[["SEGMENT_END"]])
        TraceRows("segments", x0, TraceNumber(first[["POSTERIOR_50"]]), x1, TraceNumber(first[["POSTERIOR_50"]]))
        TraceRows("rect", x0, TraceNumber(first[["POSTERIOR_10"]]), x1, TraceNumber(first[["POSTERIOR_90"]]), NA)
        if (panel[["PANEL"]] == "allele_fraction") {
            TraceRows("segments", x0, TraceNumber(1 - first[["POSTERIOR_50"]]), x1, TraceNumber(1 - first[["POSTERIOR_50"]]))
            TraceRows("rect", x0, TraceNumber(1 - first[["POSTERIOR_90"]]), x1, TraceNumber(1 - first[["POSTERIOR_10"]]), NA)
        }
    }
}
