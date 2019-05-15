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

tools <- c("minimap2", "slimm", "samtools", "bowtie2")

silent_lib <- function(...) suppressPackageStartupMessages(library(...))

#' @importFrom stringr str_split_fixed
tool_version <- function(command) {
    out <- tryCatch(
        system2(command, "--version", stdout = TRUE, stderr = FALSE)[1],
        warning = function(e) return(NULL),
        error = function(e) return(NULL)
        )
    if (!is.null(out)) {
        out <- str_split_fixed(out, " ", n = Inf)
        out <- out[length(out)]
    }
    return(out)
}

.onAttach <- function(...) {
    is_loaded <- paste0("package:", pkgs) %in% search()
    needed <- sort(pkgs[!is_loaded])

    if (length(needed) == 0) return()

    needed <- sort(needed)
    vs <- sapply(needed, function(x) as.character(packageVersion(x)))
    packageStartupMessage("Also loading:")
    packageStartupMessage(paste0("  - ", needed, "=", vs, collapse = "\n"))
    lapply(needed, silent_lib, character.only = TRUE, warn.conflicts = FALSE)
    tools_vs <- sapply(tools, tool_version)
    found <- !sapply(tools_vs, is.null)
    if (sum(found) > 0) {
        packageStartupMessage("Found tools:")
        paste0("  - ", tools[found], "=", tools_vs[found], collapse = "\n") %>%
        packageStartupMessage
    }

}
