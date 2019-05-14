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
    files <- get_files(object)
    if (is.null(config$reference)) {
        stop("must specify a reference in configuration :/")
    }
    if (grepl("(\\.fa\\.gz$)|(\\.fna\\.gz$)|(\\.fasta$)",
              config$reference)) {
        index <- file.path(alignment_dir, "index", "bref")
        flog.info("Reference is a fasta file. Building index in %s.", index)
        args <- c("--threads", config$threads, config$reference, index)
        out <- system2("bowtie2-build", args = args, env = env, stdout = FALSE,
                       stderr = FALSE)
        if (out != 0) {
            stop("failed building index :(")
        }
    } else {
        index <- config$reference
    }
    if (!dir.exists(config$alignment_dir)) {
        flog.info("Creating output directory %s.", config$alignment_dir)
        dir.create(config$alignment_dir, recursive = TRUE)
    }
    paired <- "reverse" %in% names(files)
    threads <- parse_threads(config$threads, FALSE)
    flog.info(paste("Aligning %d samples on %d threads.",
                    "Keeping up to %d secondary alignments."),
                    nrow(files), threads, config$max_hits)
    alns <- apply(files, 1, function(file) {
        file <- as.list(file)
        reads <- file$forward
        if (paired) {
            reads <- c(reads, file$reverse)
        }
        out_path <- file.path(config$alignment_dir, paste0(file$id, ".bam"))
        log_file <- file.path(config$alignment_dir, paste0(file$id, ".log"))

        if (config$use_existing && file.exists(out_path) &
            file.exists(log_file)) {
            rate <- alignment_rate(log_file)
            if (!is.null(rate)) {
                flog.info("Found existing alignment for %s. Will use that one.",
                          file$id)
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
                               "-k", config$max_hits, "-S", out_path,
                               "2>", log_file, "|", "samtools", "view",
                               "-bS", "-", ">", out_path))
        success <- system2("bowtie2", args = args)
        rate <- c(NA, NA)
        if (success == 0) {
            flog.info("Finished aligning %s.", file$id)
            rate <- alignment_rate(log_file)
        } else {
            flog.error("Failed aligning %s.", file$id)
        }
        return(data.table(id = read$id, success = success == 0,
                          alignment = out_path, reads = rate[1],
                          aligned = rate[2], rate = rate[2] / rate[1]))
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
    artifact <- list(
        alignments = alns,
        logs = logs,
        disk_size = sum(sapply(alns$alignment,
                               function(f) file.info(f)$size)),
        steps = c(object[["steps"]], "align_short_reads")
    )
    return(artifact)
}
