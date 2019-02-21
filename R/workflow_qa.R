# Copyright 2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Workflow for basic quality assessment of the read files.
#'
#' @param object An experiment data table as returned by
#'  \code{\link{find_read_files}} or a worflow object.
#' @param min_score The smallest quality score still considered okayish.
#' @param n Largest number of reads for each file. Will be sampled if more
#'  are found.
#' @return A list with the following elements:
#' \describe{
#'   \item{files}{The input files to produce the report}
#'   \item{qualities}{data table of qualities scores per cycle and sample}
#'   \item{bases}{data table of base calls per cycle and sample}
#'   \item{quality_plot}{sample quality profiles}
#'   \item{length_plot}{distribution of cycles that pass the quality cutoff}
#'   \item{entropy_plot}{sample base entropy}
#' }
quality_control <- function(object, min_score = 10, n = 1e4) {
    files <- get_files(object)
    qp <- quality_profile(files, n = n) %>% suppressPackageStartupMessages
    artifact <- list(
        files = files,
        qualities = qp$qualities,
        bases = qp$bases,
        quality_plot = plot_qualities(qp, min_score = min_score),
        length_plot = plot_lengths(qp, min_score = min_score),
        entropy_plot = plot_entropy(qp),
        steps = c(object[["steps"]], "quality_control")
    )
    return(artifact)
}
