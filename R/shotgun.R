# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Aligns metagenomic shotgun reads against a reference.
#'
#' This method uses bowtie2 and should not be too sensitive to the used pre-
#' processing.
#'
#' @param object An artifact or list of files.
#' @param ... Configuration arguments or a configuration object as obtained by
#'  \code{\link{config_align}}.
#' @return A list with the generated alignments and some general diagnostics.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom pbapply pbapply pbsapply pblapply
align_short_reads <- function(object, ...) {
    config <- config_parser(list(...), config_align)
    config$preset <- "sr"
    artifact <- align(object, config)
    artifact$steps = c(object[["steps"]], "align_short_reads")
    return(artifact)
}
