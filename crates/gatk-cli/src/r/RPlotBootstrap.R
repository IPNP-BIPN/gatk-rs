#NOTE: the port sources this file first, ahead of whatever library and script the tool runs.
#It is the gatk-rs rewrite's own: the reference has no such step, and a missing package is an
#R error there.

#Install what is missing into a user library made on demand. GATK_R_REPOSITORY overrides the CRAN
#mirror and GATK_R_LIBRARY the library. The exit status is the only channel the caller reads, so a
#package that cannot be had is named on stderr and the run ends with status 1.
EnsurePackages = function(packages) {
    missing_packages = packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
    if (length(missing_packages) == 0) {
        return(invisible(NULL))
    }
    repository = Sys.getenv("GATK_R_REPOSITORY", unset="https://cloud.r-project.org")
    user_library = strsplit(Sys.getenv("GATK_R_LIBRARY", unset=Sys.getenv("R_LIBS_USER")), .Platform$path.sep)[[1]]
    user_library = if (length(user_library) > 0) path.expand(user_library[1]) else path.expand("~/R/gatk-rs-library")
    dir.create(user_library, recursive=TRUE, showWarnings=FALSE)
    .libPaths(c(user_library, .libPaths()))
    for (package in missing_packages) {
        #install.packages only warns when it fails, so the answer is whether the namespace loads.
        if (!requireNamespace(package, quietly=TRUE)) {
            tryCatch(suppressWarnings(install.packages(package, lib=user_library, repos=repository, quiet=TRUE)),
                     error=function(e) message(conditionMessage(e)))
        }
        if (!requireNamespace(package, quietly=TRUE)) {
            message(sprintf("Could not install the R package '%s' into %s from %s.", package, user_library, repository))
            quit(save="no", status=1, runLast=FALSE)
        }
    }
    invisible(NULL)
}

#The one theme of every figure these tools draw.
GatkPlotTheme = function() {
    ggplot2::theme_classic(base_size=10) +
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, face="bold"),
                       legend.position="bottom")
}
