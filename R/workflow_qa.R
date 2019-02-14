# Copyright 2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.


#' Workflow for basic quality assessment of the read files.
#'
#' @param min_score The smallest quality score still considered okayish.
#' @return A list with the following elements:
#' \describe{
#'   \item{qualities}{data table of qualities scores per cycle and sample}
#'   \item{bases}{data table of base calls per cycle and sample}
#'   \item{quality_plot}{sample quality profiles}
#'   \item{length_plot}{distribution of cycles that pass the quality cutoff}
#'   \item{entropy_plot}{sample base entropy}
#' }
workflow_qa <- function(files, min_score = 10, n = 1e4) {
    qp <- quality_profile(files, n = 1e4)
    artifact <- list(
        qualities = qp$qualities,
        bases = qp$bases,
        quality_plot = plot_qualities(qp),
        length_plot = plot_lengths(qp, min_score = min_score),
        entropy_plot = plot_entropy(qp)
    )
    return(artifact)
}
