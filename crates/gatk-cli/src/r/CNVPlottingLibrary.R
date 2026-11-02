#The reference's CNVPlottingLibrary.R on ggplot2. What it draws is computed here as data frames, so
#the figure and the exported .data.rds are the same numbers. The arithmetic stays the reference's,
#including what looks like a mistake, and each such place says so.

#The reference shells out to grep for this; dropping the '@' lines here gives fread the same text.
ReadTSV = function(tsv_file) {
    lines = readLines(tsv_file)
    lines = lines[!startsWith(lines, "@")]
    return(suppressWarnings(fread(text=paste(lines, collapse="\n"), sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE)))
}

#The contigs along the genome, with the alternating background each one gets.
ContigTable = function(contig_names, contig_starts, contig_ends) {
    num_contigs = length(contig_names)
    return(data.frame(CONTIG=contig_names,
                      START=contig_starts,
                      END=contig_ends,
                      CENTER=(contig_starts + contig_ends) / 2,
                      FILL=ifelse(c(1:num_contigs) %% 2 == 1, "grey90", "white"),
                      stringsAsFactors=FALSE))
}

#A row of the table every exported frame carries in its "panels" attribute.
PanelRow = function(panel, sample_name, y.lab, x.lab, contigs, y.max) {
    return(data.frame(PANEL=panel, MAIN=sample_name, Y_LABEL=y.lab, X_LABEL=x.lab,
                      X_MAX=contigs[["END"]][nrow(contigs)], Y_MAX=y.max, stringsAsFactors=FALSE))
}

#Contig background for data. A limit that is not finite is the error base graphics raises from
#plot(ylim=), and the scripts begin the page ahead of this as plot() does, so that a run stopping
#here leaves the same file behind.
ContigPanel = function(panel, contigs, do_label_contigs) {
    if (!is.finite(panel[["Y_MAX"]])) {
        stop("need finite 'ylim' values")
    }
    #The labels alternate between two lines, so that neighbouring contigs do not run together.
    odd = c(1:nrow(contigs)) %% 2 == 1
    labels = ifelse(odd, contigs[["CONTIG"]], paste("\n", contigs[["CONTIG"]], sep=""))
    x_scale = if (do_label_contigs) {
        scale_x_continuous(breaks=contigs[["CENTER"]], labels=labels)
    } else {
        scale_x_continuous(breaks=NULL)
    }
    return(ggplot() +
        geom_rect(data=contigs, aes(xmin=START, xmax=END, fill=FILL), ymin=0, ymax=panel[["Y_MAX"]], inherit.aes=FALSE) +
        scale_fill_identity() +
        scale_colour_identity() +
        x_scale +
        coord_cartesian(xlim=c(0, panel[["X_MAX"]]), ylim=c(0, panel[["Y_MAX"]]), expand=FALSE) +
        labs(title=panel[["MAIN"]], x=panel[["X_LABEL"]], y=panel[["Y_LABEL"]]) +
        GatkPlotTheme() +
        theme(axis.ticks.x=element_blank()))
}

#A point size in the reference's cex, for geom_point(shape=".").
PointSize = function(point_size) {
    return(point_size * 12 / .pt)
}

#One point per row, at the genomic coordinate; a point on a contig that is not plotted has none.
CopyRatiosAsPlotted = function(copy_ratios_df, panel, color, contig_names, contig_starts) {
    copy_ratios_df[["PANEL"]] = rep(panel, nrow(copy_ratios_df))
    copy_ratios_df[["GENOMIC_COORDINATE"]] = contig_starts[match(copy_ratios_df[["CONTIG"]], contig_names)] + copy_ratios_df[["MIDDLE"]]
    copy_ratios_df[["COLOR"]] = rep(color, nrow(copy_ratios_df))
    return(copy_ratios_df)
}

PointsLayer = function(frame, y_column, point_size) {
    return(geom_point(data=frame, aes(x=GENOMIC_COORDINATE, y=.data[[y_column]], colour=COLOR),
                      shape=".", size=PointSize(point_size), na.rm=TRUE, inherit.aes=FALSE))
}

#Which points each segment draws, resolved. The assignment comes from the row order alone and never
#from the coordinates, and it keeps two things that are behaviour rather than styling:
#  - the last index is first + num_points, one too many, so consecutive segments share a point and
#    the last segment runs off the end of the table into a row of NAs;
#  - a segment with no points is skipped and does not advance the index.
#The coordinate is the segment's contig offset plus the point's own position, so a point shared
#across a contig boundary lands where the later contig puts it. The segment colours alternate on
#the segment's row and not on the segments that were drawn. 1:nrow is the reference's too: with no
#segments at all it walks rows 1 and 0, and the NA count stops the run.
SegmentPoints = function(points_df, modeled_segments_df, count_column, position_column, value, contig_names, contig_starts) {
    colors = c("coral", "dodgerblue")
    counts = modeled_segments_df[[count_column]]
    points_start_index = 1
    pieces = list()
    for (s in 1:nrow(modeled_segments_df)) {
        num_points = counts[s]
        if (num_points == 0) {
            next
        }
        indices = points_start_index:(points_start_index + num_points)
        offset = contig_starts[match(modeled_segments_df[s, "CONTIG"], contig_names)]
        pieces[[length(pieces) + 1]] = data.frame(
            SEGMENT=s,
            CONTIG=modeled_segments_df[s, "CONTIG"],
            SEGMENT_START=offset + modeled_segments_df[s, "START"],
            SEGMENT_END=offset + modeled_segments_df[s, "END"],
            POINT_INDEX=indices,
            GENOMIC_COORDINATE=offset + points_df[indices, position_column],
            VALUE=value[indices],
            COLOR=colors[s %% 2 + 1],
            stringsAsFactors=FALSE)
        points_start_index = points_start_index + num_points
    }
    if (length(pieces) == 0) {
        return(data.frame(SEGMENT=integer(0), CONTIG=character(0), SEGMENT_START=numeric(0), SEGMENT_END=numeric(0),
                          POINT_INDEX=integer(0), GENOMIC_COORDINATE=numeric(0), VALUE=numeric(0), COLOR=character(0),
                          stringsAsFactors=FALSE))
    }
    return(do.call(rbind, pieces))
}

#Each segment's posterior deciles beside its points, and the panel's name on every row.
WithPosteriors = function(frame, panel, modeled_segments_df, posterior_10, posterior_50, posterior_90) {
    frame = cbind(data.frame(PANEL=rep(panel, nrow(frame)), stringsAsFactors=FALSE), frame)
    frame[["POSTERIOR_10"]] = posterior_10[frame[["SEGMENT"]]]
    frame[["POSTERIOR_50"]] = posterior_50[frame[["SEGMENT"]]]
    frame[["POSTERIOR_90"]] = posterior_90[frame[["SEGMENT"]]]
    return(frame)
}

#The segments a frame holds, one row each, with the posterior median drawn at MEDIAN and the band
#between LOWER and UPPER.
SegmentBands = function(frame, median, lower, upper) {
    bands = frame[!duplicated(frame[["SEGMENT"]]), c("SEGMENT", "SEGMENT_START", "SEGMENT_END")]
    first = match(bands[["SEGMENT"]], frame[["SEGMENT"]])
    bands[["MEDIAN"]] = median[first]
    bands[["LOWER"]] = lower[first]
    bands[["UPPER"]] = upper[first]
    return(bands)
}

BandLayers = function(bands) {
    return(list(
        geom_segment(data=bands, aes(x=SEGMENT_START, xend=SEGMENT_END, y=MEDIAN, yend=MEDIAN),
                     colour="black", linewidth=0.5, na.rm=TRUE, inherit.aes=FALSE),
        geom_rect(data=bands, aes(xmin=SEGMENT_START, xmax=SEGMENT_END, ymin=LOWER, ymax=UPPER),
                  colour="black", fill=NA, linewidth=0.25, na.rm=TRUE, inherit.aes=FALSE)))
}
