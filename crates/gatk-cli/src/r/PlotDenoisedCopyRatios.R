#NOTE: the port first sources RPlotBootstrap.R and CNVPlottingLibrary.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

EnsurePackages(c("optparse", "data.table", "ggplot2", "cowplot"))

library(optparse)
library(data.table)
library(ggplot2)

option_list = list(
    make_option(c("--sample_name"), dest="sample_name", action="store"),
    make_option(c("--standardized_copy_ratios_file"), dest="standardized_copy_ratios_file", action="store"),
    make_option(c("--denoised_copy_ratios_file"), dest="denoised_copy_ratios_file", action="store"),
    make_option(c("--contig_names"), dest="contig_names", action="store"),      #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--contig_lengths"), dest="contig_lengths", action="store"),  #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--maximum_copy_ratio"), dest="maximum_copy_ratio", action="store", type="double"),
    make_option(c("--point_size_copy_ratio"), dest="point_size_copy_ratio", action="store", type="double"),
    make_option(c("--output_dir"), dest="output_dir", action="store"),
    make_option(c("--output_prefix"), dest="output_prefix", action="store"))

opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["sample_name"]]
standardized_copy_ratios_file = opt[["standardized_copy_ratios_file"]]
denoised_copy_ratios_file = opt[["denoised_copy_ratios_file"]]
contig_names_string = opt[["contig_names"]]
contig_lengths_string = opt[["contig_lengths"]]
maximum_copy_ratio = opt[["maximum_copy_ratio"]]
point_size_copy_ratio = opt[["point_size_copy_ratio"]]
output_dir = opt[["output_dir"]]
output_prefix = opt[["output_prefix"]]

#check that input files exist; if not, quit with error code that GATK will pick up
if (!all(file.exists(c(standardized_copy_ratios_file, denoised_copy_ratios_file)))) {
    quit(save="no", status=1, runLast=FALSE)
}

contig_names = strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]]
contig_lengths = as.numeric(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))

#the median of the absolute differences between neighbours, which is neither R's mad() nor a deviation from any median
CalculateMedianAbsoluteDeviation = function(dat) {
    return(median(abs(diff(dat))))
}

#plotting is extracted to a function for debugging purposes
WriteDenoisingPlots = function(sample_name, standardized_copy_ratios_file, denoised_copy_ratios_file, contig_names, output_dir, output_prefix) {
    standardized_copy_ratios_df = ReadTSV(standardized_copy_ratios_file)
    denoised_copy_ratios_df = ReadTSV(denoised_copy_ratios_file)

    #transform to linear copy ratio
    standardized_copy_ratios_df[["COPY_RATIO"]] = 2^standardized_copy_ratios_df[["LOG2_COPY_RATIO"]]
    denoised_copy_ratios_df[["COPY_RATIO"]] = 2^denoised_copy_ratios_df[["LOG2_COPY_RATIO"]]

    #determine copy-ratio midpoints
    standardized_copy_ratios_df[["MIDDLE"]] = round((standardized_copy_ratios_df[["START"]] + standardized_copy_ratios_df[["END"]]) / 2)
    denoised_copy_ratios_df[["MIDDLE"]] = round((denoised_copy_ratios_df[["START"]] + denoised_copy_ratios_df[["END"]]) / 2)

    #write the MAD files
    standardizedMAD = CalculateMedianAbsoluteDeviation(standardized_copy_ratios_df[["COPY_RATIO"]])
    denoisedMAD = CalculateMedianAbsoluteDeviation(denoised_copy_ratios_df[["COPY_RATIO"]])
    write.table(round(standardizedMAD, 3), file.path(output_dir, paste(output_prefix, ".standardizedMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(denoisedMAD, 3), file.path(output_dir, paste(output_prefix, ".denoisedMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(standardizedMAD - denoisedMAD, 3), file.path(output_dir, paste(output_prefix, ".deltaMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round((standardizedMAD - denoisedMAD) / standardizedMAD, 3), file.path(output_dir, paste(output_prefix, ".scaledDeltaMAD.txt", sep="")), col.names=FALSE, row.names=FALSE)

    #plot standardized and denoised copy ratio on top of each other
    pre_color_blue = "#3B5DFF"
    post_color_green = "#4FC601"

    contigs = ContigTable(contig_names, contig_starts, contig_ends)
    standardized_frame = CopyRatiosAsPlotted(standardized_copy_ratios_df, "standardized", pre_color_blue, contig_names, contig_starts)
    denoised_frame = CopyRatiosAsPlotted(denoised_copy_ratios_df, "denoised", post_color_green, contig_names, contig_starts)

    #plot up to maximum_copy_ratio (or full range, if maximum_copy_ratio = Infinity); the two panels then get different limits
    denoising_plot_file = file.path(output_dir, paste(output_prefix, ".denoised.png", sep=""))
    png(denoising_plot_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    grid::grid.newpage()
    maximum_standardized_copy_ratio = if(is.finite(maximum_copy_ratio)) maximum_copy_ratio else 1.05 * max(standardized_copy_ratios_df[["COPY_RATIO"]])
    standardized_row = PanelRow("standardized", sample_name, "standardized copy ratio", paste("median absolute deviation = ", round(standardizedMAD, 3), sep=""), contigs, maximum_standardized_copy_ratio)
    standardized_panel = ContigPanel(standardized_row, contigs, FALSE) +
        PointsLayer(standardized_frame, "COPY_RATIO", point_size_copy_ratio)
    maximum_denoised_copy_ratio = if(is.finite(maximum_copy_ratio)) maximum_copy_ratio else 1.05 * max(denoised_copy_ratios_df[["COPY_RATIO"]])
    denoised_row = PanelRow("denoised", sample_name, "denoised copy ratio", paste("median absolute deviation = ", round(denoisedMAD, 3), sep=""), contigs, maximum_denoised_copy_ratio)
    denoised_panel = ContigPanel(denoised_row, contigs, TRUE) +
        PointsLayer(denoised_frame, "COPY_RATIO", point_size_copy_ratio)
    figure = cowplot::plot_grid(standardized_panel, denoised_panel, ncol=1, align="v")
    print(figure, newpage=FALSE)
    dev.off()

    #check for created files and quit with error code if not found
    if (!all(file.exists(c(denoising_plot_file)))) {
        quit(save="no", status=1, runLast=FALSE)
    }

    #the frames as plotted, and the figure with the panels it was made from
    frame = rbind(standardized_frame, denoised_frame)
    rownames(frame) = NULL
    attr(frame, "panels") = rbind(standardized_row, denoised_row)
    attr(frame, "contigs") = contigs
    saveRDS(frame, file.path(output_dir, paste(output_prefix, ".denoised.data.rds", sep="")))
    attr(figure, "panels") = list(standardized=standardized_panel, denoised=denoised_panel)
    saveRDS(figure, file.path(output_dir, paste(output_prefix, ".denoised.ggplot.rds", sep="")))
}

WriteDenoisingPlots(sample_name, standardized_copy_ratios_file, denoised_copy_ratios_file, contig_names, output_dir, output_prefix)
