#NOTE: the port first sources RPlotBootstrap.R
EnsurePackages(c("dplyr", "ggplot2", "cowplot"))

library(dplyr)
library(ggplot2)

inverseLogitFn = function(score) {
    v = 10.0^(-score/10.0)
    v / (1.0 + v)
}

args <- commandArgs(TRUE)
verbose = TRUE

variant_eval_filename <- args[1]
output_basename <- tools::file_path_sans_ext(args[2])
sample <- args[3]


make_sample_df = function(filename, sample) {
    df = read.table(filename, sep="", stringsAsFactors = F, header = T)
    df = df %>% dplyr::filter(Filter == "called" ) # only called sites
    df = df %>%
	    transmute(AF_bin = inverseLogitFn(AlleleFrequency),
                      EvaluationType = EvalFeatureInput,
                      avgVarAF= avgVarAF)

    array_bins = subset(df, EvaluationType == "eval", select = -EvaluationType)
    thousand_g_bins = subset(df, EvaluationType == "thousand_genomes", select = -EvaluationType)
    merged = merge(array_bins, thousand_g_bins,
          by = c("AF_bin"), suffixes = c(".array", ".thousand_g")) %>%
      mutate(sample = sample)

    # The merge is an inner join, so a bin only one side has disappears; this is the only word of it.
    message(sprintf("%d bins from eval and %d from thousand_genomes, %d in both.",
                    nrow(array_bins), nrow(thousand_g_bins), nrow(merged)))
    merged
}

write_png = function(suffix, figure, width = 6) {
    png(paste(output_basename, suffix, ".png", sep = ""), width = width, height = 6, units = 'in', res = 300)
    print(figure)
    dev.off()
}

run_single_sample = function(filename, output_basename, sample) {
    single_samp = make_sample_df(filename, sample)

    af_differences = ggplot(single_samp, aes(x = avgVarAF.thousand_g, y = avgVarAF.array)) + geom_point() +
        geom_abline(slope = 1, intercept = 0, color = "red", alpha = 0.7) +
        labs(x = "AF from Thousand Genomes", y = "AF from Array") +
        GatkPlotTheme()
    write_png(".af_differences", af_differences)

    long_form = rbind(
        data.frame(AF_bin = single_samp[["AF_bin"]], AF = single_samp[["avgVarAF.array"]], source = rep("array", nrow(single_samp))),
        data.frame(AF_bin = single_samp[["AF_bin"]], AF = single_samp[["avgVarAF.thousand_g"]], source = rep("thousand genomes", nrow(single_samp))))
    af = ggplot(long_form, aes(x = AF_bin, y = AF, color = source)) + geom_point() +
        geom_abline(slope = 1, intercept = 0, color = "black", alpha = 0.7) +
        scale_color_manual(values = c("array" = "#F8766D", "thousand genomes" = "#00BFC4"), name = NULL) +
        labs(x = "AF Bin", y = "AF in VCFs") +
        GatkPlotTheme()
    write_png(".af", af)

    # A third figure beside the reference's two. The grid is built with the device open, because
    # outside one its layout would open R's default device and leave an Rplots.pdf behind.
    png(paste(output_basename, ".af_combined", ".png", sep = ""), width = 12, height = 6, units = 'in', res = 300)
    print(cowplot::plot_grid(af_differences, af, ncol = 2, align = "h"))
    dev.off()

    saveRDS(single_samp, paste(output_basename, ".data.rds", sep = ""))
    saveRDS(af, paste(output_basename, ".af.ggplot.rds", sep = ""))
    saveRDS(af_differences, paste(output_basename, ".af_differences.ggplot.rds", sep = ""))
}


run_single_sample(variant_eval_filename, output_basename, sample)
