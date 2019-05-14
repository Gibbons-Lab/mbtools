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
    threads = TRUE,
    alignment_dir = "alignments",
    max_hits = 100,
    use_existing = TRUE
))

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
                  "-I", "64G",, "--secondary", "yes",
                  "-N", config$max_hits, config$reference, reads)
        args <- append(args, c(paste0("2>", log_file), "|", "samtools",
                            "view", "-bS", "-", ">", out_path))
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
    alns <- rbindlist(alns)
    artifact <- list(
        alignments = alns,
        logs = logs,
        disk_size = sum(sapply(alns$alignment,
                               function(f) file.info(f)$size))
    )
    return(artifact)
}
