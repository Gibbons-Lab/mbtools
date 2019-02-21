# Copyright 2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Build a configuration for raw read preprocessing.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the DADA2 workflow.
#' @export
#' @examples
#'  config <- config_preprocess(truncLen = c(240, 250))
config_preprocess <- function(...) {
    config <- list(
        threads = TRUE,
        out_dir = "preprocessed",
        trimLeft = 10,
        truncLen = 0,
        maxEE = 2
    )
    args <- list(...)
    for (arg in names(args)) {
        config[[arg]] <- args[[arg]]
    }
    return(config)
}

#' Runs preprocessing of sequencing reads.
#'
#' @param object An experiment data table as returned by
#'  \code{\link{find_read_files}} or a worflow object.
#' @param config A configuration file as returned by
#'  \code{\link{config_preprocess}}.
#' @return A list containing the workflow results:
#' \describe{
#'   \item{passed_reads}{How many reads were kept in each step. Rows are
#'     samples and columns are workflow steps.}
#'   \item{files}{Preprocessed sequencing files list.}
#' }
#' @export
#'
#' @importFrom dada2 filterAndTrim
preprocess <- function(object, config) {
    files <- get_files(object)
    files <- copy(files)
    if (!"run" %in% names(files)) {
        files[, "run" := "all"]
    }
    paired <- "reverse" %in% names(files)

    flog.info("Preprocessing reads for %d %s-end samples...",
              nrow(files), ifelse(paired, "paired", "single"))
    passed_files <- copy(files)
    passed_files$forward <- file.path(config$out_dir, "preprocessed",
                                      basename(files$forward))
    if (paired) {
        passed_files$reverse <- file.path(config$out_dir, "preprocessed",
                                          basename(files$reverse))
        passed_stats <- filterAndTrim(
            fwd = files$forward, filt = passed_files$forward,
            rev = files$reverse, filt.rev = passed_files$reverse,
            trimLeft = config$trimLeft, truncLen = config$truncLen,
            maxEE = config$maxEE, multithread = config$threads
        )
    } else {
        passed_stats <- filterAndTrim(
            fwd = files$forward, filt = passed_files$forward,
            trimLeft = config$trimLeft, truncLen = config$truncLen,
            maxEE = config$maxEE, multithread = config$threads
        )
    }
    passed_stats <- as.data.table(passed_stats) %>%
                    setNames(c("raw", "preprocessed"))
    passed_stats[, "id" := files$id]
    flog.info("%.3g/%.3g (%.2f%%) reads passed preprocessing.",
              passed_stats[, sum(preprocessed)],
              passed_stats[, sum(raw)],
              passed_stats[, 100 * mean(preprocessed / raw)])
    artifact <- list(
        files = passed_files,
        passed = passed_stats,
        steps = c(object[["steps"]], "preprocess")
    )
    return(artifact)
}
