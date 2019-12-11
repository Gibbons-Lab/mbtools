# General puprose alignment code.

#' Build a configuration for the alignment workflows.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the long read alignment
#'  workflow.
#' @export
#' @examples
#'  config <- config_align(reference = "refs/mouse.fna.gz")
config_align <- config_builder(list(
    reference = NULL,
    build_index = FALSE,
    threads = getOption("mc.cores", 1),
    alignment_dir = "alignments",
    max_hits = 100,
    use_existing = TRUE,
    limited_memory = FALSE
))

#' Heuristic to get the median read length from a BAM file.access
#'
#' @param bam_file Path to a bam file.
#' @param n How many alignments to use.
#' @return The median read length in the file.
#' @importFrom Rsamtools BamFile yieldSize scanBam
read_length <- function(bam_file, n = 100) {
    bam <- BamFile(bam_file)
    yieldSize(bam) <- n
    open(bam)
    alns <- scanBam(bam)
    close(bam)
    l <- alns[[1]]$qwidth
    return(median(l, na.rm = TRUE))
}


align <- function(object, config) {
    files <- get_files(object)
    if (is.null(config$reference)) {
        stop("must specify a reference in configuration :/")
    }
    if (!dir.exists(config$alignment_dir)) {
        flog.info("Creating output directory %s.", config$alignment_dir)
        dir.create(config$alignment_dir, recursive = TRUE)
    }
    paired <- "reverse" %in% names(files)
    threads <- parse_threads(config$threads, FALSE)

    if (config$build_index) {
        index <- file.path(config$alignment_dir, "index",
                           paste0(basename(config$reference), ".mmi"))
        args <- c("-x", config$preset, "--secondary=yes", "-I", "100G",
                  "-t", threads, "-d", index, config$reference)
        flog.info("Building index in %s.", index)
        ec <- system2("minimap2", args)
    } else {
        index <- config$reference
    }

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

        if (config$use_existing && file.exists(out_path)) {
            flog.info("Found existing alignment for %s. Will use that one.",
                      file$id)
            return(data.table(id = file$id, alignment = out_path,
                              success = TRUE))
        }

        args <- c("-acx", config$preset, "-t", threads,
                  "--secondary=yes", "-N", config$max_hits)
        if (config$limited_memory) {
            args <- append(args, c("--split-prefix",
                                   file.path(config$alignment_dir, "prefix")))
        } else {
            args <- append(args, c("-I", "500G"))
        }
        args <- append(args, c(index, reads))
        args <- append(args, c(paste0("2>", log_file), "|", "samtools",
                            "sort", "--threads", threads - 1, "-o", out_path))
        success <- system2("minimap2", args = args)
        if (success == 0) {
            flog.info("Finished aligning %s.", file$id)
        } else {
            flog.error("Failed aligning %s.", file$id)
        }
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
    if (config$limited_memory) {
        unlink(file.path(config$alignment_dir, "prefix"), recursive = TRUE)
    }
    alns <- rbindlist(alns)
    disk_size = sum(sapply(alns$alignment,
                           function(f) file.info(f)$size))
    class(disk_size) <- "object_size"
    artifact <- list(
        alignments = alns,
        logs = logs,
        disk_size = disk_size
    )
    return(artifact)
}
