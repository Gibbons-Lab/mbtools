# Copyright 2018 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Align long reads (for instance from nanopore sequencing) to a reference
#' database.
#'
#' @param object An artifact or list of files.
#' @param ... Configuration arguments or a configuration object as obtained by
#' \code{\link{config_align}}.
#' @return A list with the generated alignments and some general diagnostics.
#'
#' @export
align_long_reads <- function(object, ...) {
    config <- config_parser(list(...), config_align)
    config$preset <- "map-ont"
    artifact <- align(object, config)
    artifact$steps = c(object[["steps"]], "align_long_reads")
    return(artifact)
}
