# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' @import ShortRead ggplot2 phyloseq
#' @importFrom utils packageVersion download.file read.table
#' @importFrom stringr str_match
#' @importFrom Biostrings BStringSet
#' @importFrom R.utils gunzip
NULL

pkgs <- c("ggplot2", "dada2", "phyloseq", "ShortRead",
          "data.table", "yaml", "magrittr")

silent_lib <- function(...) suppressPackageStartupMessages(library(...))

.onAttach <- function(...) {
    is_loaded <- paste0("package:", pkgs) %in% search()
    needed <- sort(pkgs[!is_loaded])

    if (length(needed) == 0) return()

    needed <- sort(needed)
    vs <- sapply(needed, function(x) as.character(packageVersion(x)))
    packageStartupMessage(paste("Also loading:", needed, vs, collapse = "\n"))
    lapply(needed, silent_lib, character.only = TRUE, warn.conflicts = FALSE)
}
