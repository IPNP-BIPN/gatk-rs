#NOTE: the port first sources RPlotBootstrap.R and CNVPlottingLibrary.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

EnsurePackages(c("optparse", "data.table", "ggplot2", "cowplot"))

library(optparse)
library(data.table)
library(ggplot2)

option_list = list(
    make_option(c("--sample_name"), dest="sample_name", action="store"),
    make_option(c("--denoised_copy_ratios_file"), dest="denoised_copy_ratios_file", action="store"),
    make_option(c("--allelic_counts_file"), dest="allelic_counts_file", action="store"),
    make_option(c("--modeled_segments_file"), dest="modeled_segments_file", action="store"),
    make_option(c("--contig_names"), dest="contig_names", action="store"),      #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--contig_lengths"), dest="contig_lengths", action="store"),  #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--maximum_copy_ratio"), dest="maximum_copy_ratio", action="store", type="double"),
    make_option(c("--point_size_copy_ratio"), dest="point_size_copy_ratio", action="store", type="double"),
    make_option(c("--point_size_allele_fraction"), dest="point_size_allele_fraction", action="store", type="double"),
    make_option(c("--output_file"), dest="output_file", action="store"))

opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["sample_name"]]
denoised_copy_ratios_file = opt[["denoised_copy_ratios_file"]]
allelic_counts_file = opt[["allelic_counts_file"]]
modeled_segments_file = opt[["modeled_segments_file"]]
contig_names_string = opt[["contig_names"]]
contig_lengths_string = opt[["contig_lengths"]]
maximum_copy_ratio = opt[["maximum_copy_ratio"]]
point_size_copy_ratio = opt[["point_size_copy_ratio"]]
point_size_allele_fraction = opt[["point_size_allele_fraction"]]
output_file = opt[["output_file"]]

#check that required input files exist; if not, quit with error code that GATK will pick up
if (!file.exists(modeled_segments_file)) {
    quit(save="no", status=1, runLast=FALSE)
}

contig_names = strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]]
contig_lengths = as.numeric(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))

#plotting is extracted to a function for debugging purposes
WriteModeledSegmentsPlot = function(sample_name, allelic_counts_file, denoised_copy_ratios_file, modeled_segments_file, contig_names, output_file) {
    modeled_segments_df = ReadTSV(modeled_segments_file)
    contigs = ContigTable(contig_names, contig_starts, contig_ends)

    num_plots = ifelse(all(file.exists(c(denoised_copy_ratios_file, allelic_counts_file))), 2, 1)
    png(output_file, 12, 3.5 * num_plots, units="in", type="cairo", res=300, bg="white")

    panels = list()
    panel_rows = list()
    frames = list()

    if (file.exists(denoised_copy_ratios_file) && denoised_copy_ratios_file != "null") {
        grid::grid.newpage()
        denoised_copy_ratios_df = ReadTSV(denoised_copy_ratios_file)

        #transform to linear copy ratio
        denoised_copy_ratios_df[["COPY_RATIO"]] = 2^denoised_copy_ratios_df[["LOG2_COPY_RATIO"]]

        #determine copy-ratio midpoints
        denoised_copy_ratios_df[["MIDDLE"]] = round((denoised_copy_ratios_df[["START"]] + denoised_copy_ratios_df[["END"]]) / 2)

        #plot up to maximum_copy_ratio (or full range, if maximum_copy_ratio = Infinity)
        maximum_denoised_copy_ratio = if(is.finite(maximum_copy_ratio)) maximum_copy_ratio else 1.05 * max(denoised_copy_ratios_df[["COPY_RATIO"]])
        panel_row = PanelRow("copy_ratio", sample_name, "denoised copy ratio", "contig", contigs, maximum_denoised_copy_ratio)

        frame = SegmentPoints(denoised_copy_ratios_df, modeled_segments_df, "NUM_POINTS_COPY_RATIO", "MIDDLE", denoised_copy_ratios_df[["COPY_RATIO"]], contig_names, contig_starts)
        frame = WithPosteriors(frame, "copy_ratio", modeled_segments_df,
                               2^modeled_segments_df[["LOG2_COPY_RATIO_POSTERIOR_10"]],
                               2^modeled_segments_df[["LOG2_COPY_RATIO_POSTERIOR_50"]],
                               2^modeled_segments_df[["LOG2_COPY_RATIO_POSTERIOR_90"]])

        panels[[length(panels) + 1]] = ContigPanel(panel_row, contigs, TRUE) +
            PointsLayer(frame, "VALUE", point_size_copy_ratio) +
            BandLayers(SegmentBands(frame, frame[["POSTERIOR_50"]], frame[["POSTERIOR_10"]], frame[["POSTERIOR_90"]]))
        panel_rows[[length(panel_rows) + 1]] = panel_row
        frames[[length(frames) + 1]] = frame
    }

    if (file.exists(allelic_counts_file) && allelic_counts_file != "null") {
        if (length(panels) == 0) {
            grid::grid.newpage()
        }
        allelic_counts_df = ReadTSV(allelic_counts_file)

        panel_row = PanelRow("allele_fraction", sample_name, "alternate-allele fraction", "contig", contigs, 1.0)

        alternate_allele_fractions = allelic_counts_df[["ALT_COUNT"]] / (allelic_counts_df[["ALT_COUNT"]] + allelic_counts_df[["REF_COUNT"]])
        frame = SegmentPoints(allelic_counts_df, modeled_segments_df, "NUM_POINTS_ALLELE_FRACTION", "POSITION", alternate_allele_fractions, contig_names, contig_starts)
        frame = WithPosteriors(frame, "allele_fraction", modeled_segments_df,
                               modeled_segments_df[["MINOR_ALLELE_FRACTION_POSTERIOR_10"]],
                               modeled_segments_df[["MINOR_ALLELE_FRACTION_POSTERIOR_50"]],
                               modeled_segments_df[["MINOR_ALLELE_FRACTION_POSTERIOR_90"]])

        #the bands twice: the minor allele fraction, and mirrored at one minus it
        panels[[length(panels) + 1]] = ContigPanel(panel_row, contigs, TRUE) +
            PointsLayer(frame, "VALUE", point_size_allele_fraction) +
            BandLayers(SegmentBands(frame, frame[["POSTERIOR_50"]], frame[["POSTERIOR_10"]], frame[["POSTERIOR_90"]])) +
            BandLayers(SegmentBands(frame, 1 - frame[["POSTERIOR_50"]], 1 - frame[["POSTERIOR_90"]], 1 - frame[["POSTERIOR_10"]]))
        panel_rows[[length(panel_rows) + 1]] = panel_row
        frames[[length(frames) + 1]] = frame
    }

    #no panel draws no page, and the png device then leaves no file for the check below to find
    figure = NULL
    if (length(panels) > 0) {
        figure = cowplot::plot_grid(plotlist=panels, ncol=1, align="v")
        print(figure, newpage=FALSE)
    }
    dev.off()

    #check for created file and quit with error code if not found
    if (!file.exists(output_file)) {
        quit(save="no", status=1, runLast=FALSE)
    }

    #the segments joined to their points as plotted, and the figure with the panels it was made from
    output_base = sub("\\.png$", "", output_file)
    frame = do.call(rbind, frames)
    rownames(frame) = NULL
    attr(frame, "panels") = do.call(rbind, panel_rows)
    attr(frame, "contigs") = contigs
    saveRDS(frame, paste(output_base, ".data.rds", sep=""))
    names(panels) = attr(frame, "panels")[["PANEL"]]
    attr(figure, "panels") = panels
    saveRDS(figure, paste(output_base, ".ggplot.rds", sep=""))
}

WriteModeledSegmentsPlot(sample_name, allelic_counts_file, denoised_copy_ratios_file, modeled_segments_file, contig_names, output_file)
