# Copyright 2018 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Build a configuration for the long read alignment workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the long read alignment
#'  workflow.
#' @export
#' @examples
#'  config <- config_align_long(index = "refs/mouse")
config_align_long <- function(...) {
    config <- list(
        reference = NULL,
        threads = 1,
        alignment_dir = "alignments",
        max_hits = 100,
        progress = TRUE,
        use_existing = TRUE
    )
    args <- list(...)
    for (arg in names(args)) {
        config[[arg]] <- args[[arg]]
    }
    return(config)
}

#' Align long reads (for instance from nanopore sequencing) to a reference
#' database.
#'
#' @param object An artifact or list of files.
#' @param config A configuration as obtained by \code{\link{config_align_short}}.
#' @return A list with the generated alignments and some general diagnostics.
#'
#' @export
align_long_reads <- function(object, config) {
    files <- get_files(object)
    if (is.null(config$reference)) {
        stop("must specify a reference genome in configuration :/")
    }
    if (!dir.exists(config$alignment_dir)) {
        flog.info("Creating output directory %s.", config$alignment_dir)
        dir.create(config$alignment_dir, recursive = TRUE)
    }
    paired <- "reverse" %in% names(files)
    flog.info(paste("Aligning %d samples on %d threads.",
                    "Keeping up to %d secondary alignments."),
                    nrow(files), config$threads, config$max_hits)
    alns <- apply(files, 1, function(file) {
        file <- as.list(file)
        flog.info("Aligning %s...", file$id)
        reads <- file$forward
        if (paired) {
            reads <- c(reads, file$reverse)
        }
        out_path <- file.path(config$alignment_dir, paste0(file$id, ".bam"))
        log_file <- file.path(config$alignment_dir, paste0(file$id, ".log"))

        if (config$use_existing && file.exists()) {
            return(data.table(id = file$id, alignment = out_path, success = 0))
        }

        args <- c("-acx", "map-ont", "-t", config$threads, "-N",
                  config$max_hits, config$reference, reads)
        args <- append(args, c(paste0("2>", log_file), "|", "samtools",
                            "view", "-bS", "-", ">", out_path))
        success <- system2("minimap2", args = args)
        return(data.table(id = file$id, alignment = out_path,
                          success = success == 0))
    })
    logs <- lapply(files$id, function(id) {
        log_file <- file.path(config$alignment_dir, paste0(id, ".log"))
        if (!file.exists(log_file)) {
            return(NA)
        }
        content <- readChar(log_file, min(file.info(log_file)$size, 1e8))
        file.remove(log_file)
        return(content)
    })
    alns <- rbindlist(alns)
    if (alns[, any(!success)]) {
        flog.error("%d alignments failed!", alns[, sum(!success)])
    }
    artifact <- list(
        alignments = alns,
        logs = logs,
        disk_size = sum(sapply(alns$alignment,
                               function(f) file.info(f)$size)),
        steps = c(object[["steps"]], "align_long_reads")
    )
    return(artifact)
}
