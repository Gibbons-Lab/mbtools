# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' @importFrom stringr str_match
alignment_rate <- function(log_file) {
    if (file.exists(log_file)) {
        content <- readChar(log_file, file.info(log_file)$size)
        total <- single <- multi <- NA
        match <- str_match(content, "(\\d+) reads; of")
        if (nrow(match) > 0) {
            total <- as.numeric(match[, 2])
        }
        match <- str_match(content, "(\\d+) \\(\\d+.\\d+%\\) aligned >1 times")
        if (nrow(match) > 0) {
            multi <- as.numeric(match[, 2])
        }
        match <- str_match(content,
                           "(\\d+) \\(\\d+.\\d+%\\) aligned exactly 1 time")
        if (nrow(match) > 0) {
            single <- as.numeric(match[, 2])
        }
    }
    return(c(total = total, aligned = single + multi))
}

#' Build a configuration for the short read alignment workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the short read alignment
#'  workflow.
#' @export
#' @examples
#'  config <- config_count(index = "refs/mouse")
config_align_short <- function(...) {
    config <- list(
        reference = NULL,
        threads = 1,
        alignment_dir = "alignments",
        bowtie2_path = NULL,
        samtools_path = NULL,
        max_hits = 100,
        progress = TRUE
    )
    args <- list(...)
    for (arg in names(args)) {
        config[[arg]] <- args[[arg]]
    }
    return(config)
}

#' Aligns metagenomic shotgun sample against a reference
#'
#' This method uses bowtie2 and should not be too sensitive to the used pre-
#' processing.
#'
#' @param object An artifact or list of files.
#' @param config A configuration as obtained by \code{\link{config_align_short}}.
#' @return A list with the generated alignments and some general diagnostics.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom pbapply pbapply pbsapply pblapply
align_short_reads <- function(object, config) {
    files <- get_files(object)
    if (is.null(config$reference)) {
        stop("must specify index in configuration :/")
    }
    if (grepl("(\\.fa\\.gz$)|(\\.fna\\.gz$)|(\\.fasta$)",
              config$reference)) {
        flog.info("Reference is a fasta file. Building index in %s.", index)
        index <- file.path(tempdir(), "bref")
        args <- c("--threads", config$threads, config$reference, index)
        out <- system2("bowtie2-build", args = args, env = env, stdout = FALSE,
                       stderr = FALSE)
        if (out != 0) {
            stop("failed building index :(")
        }
    } else {
        index <- config$reference
    }
    env <- character()
    if (!is.null(config$bowtie2_path)) {
        env <- paste0("PATH=", config$bowtie2_path)
    }
    if (!is.null(config$samtools_path)) {
        env <- paste0(env, ":", config$samtools_path)
    }
    paired <- "reverse" %in% names(files)
    if (!dir.exists(config$alignment_dir)) {
        flog.info("Creating output directory %s.", config$alignment_dir)
        dir.create(config$alignment_dir, recursive = TRUE)
    }

    if (config$progress && interactive()) {
        apply_fun <- pbapply
    } else apply_fun <- apply

    flog.info("Aligning %d samples...", nrow(files))
    alignments <- apply_fun(files, 1, function(read) {
        read <- as.list(read)
        log_file <- file.path(config$alignment_dir, paste0(read$id, ".log"))

        if (file.exists(log_file)) {
            rate <- alignment_rate(log_file)

            if (config$bam) {
                out_path <- file.path(config$alignment_dir,
                                      paste0(read$id, ".bam"))
            } else {
                out_path <- file.path(config$alignment_dir,
                                      paste0(read$id, ".sam"))
            }
            if (!is.null(rate) && file.exists(out_path)) {
                return(data.table(id = read$id, success = TRUE,
                       log = log_file, alignment = out_path, reads = rate[1],
                       aligned = rate[2], rate = rate[2] / rate[1]))
            }
        }

        args <- c("-x", index)
        if (paired) {
            args <- append(args, c("-1", read$forward, "-2", read$reverse))
        } else {
            args <- append(args, c("-U", read$forward))
        }
        args <- append(args, c("-q", "--no-unal", "--mm", "-p", config$threads,
                               "-k", config$max_hits))


        out_path <- file.path(config$alignment_dir,
                                  paste0(read$id, ".sam"))
        args <- append(args, c("-S", out_path, "2>", log_file))
        success <- system2("bowtie2", args = args, env = env)
        rate <- c(NA, NA)
        if (success == 0) {
            rate <- alignment_rate(log_file)
        }

        return(data.table(id = read$id, success = success == 0, log = log_file,
                          alignment = out_path, reads = rate[1],
                          aligned = rate[2], rate = rate[2] / rate[1]))
    }) %>% rbindlist()
    artifact <- list(
        alignments = alignments,
        disk_size = sum(sapply(alignments$alignment,
                               function(f) file.info(f)$size)),
        steps = c(object[["steps"]], "align_short_reads")
    )
    return(artifact)
}
