# Copyright 2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Build a configuration for chimera removal.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the DADA2 workflow.
#' @export
#' @examples
#'  config <- config_preprocess(truncLen = c(240, 250))
config_chimera <- config_builder(list(
    threads = getOption("mc.cores", 1),
    out_dir = "non_chimeric",
    preset = "ava-ont",
    seed_distance = 500
))

#' @importFrom tools file_path_sans_ext
filter_chimeras <- function(files, output, config) {
    stats <- lapply(1:nrow(files), function(i) {
        infile <- files[i]
        spl <- split_ext(infile)
        outfile <- output[i]
        filter_file <- paste0(spl[1], "_filtered.", spl[2])
        yacfile <- file.path(
                config$out_dir,
                paste0(spl[1], ".yacrd"))
        args <- c("-x", config$preset, "-g", config$seed_distance,
                 "-t", config$threads, infile, infile, "2>",
                 "/dev/null", "|",
                 "yacrd", "chimeric", "-f", infile, ">", yacfile)
        ret <- system2("minimap2", args = args)
        if (ret != 0) {
            stop(sprintf(
                "Chimera detection failed for file %s :(", infile))
        }
        file.copy(filter_file, outfile)
        file.remove(filter_file)
        nseq <- fastq.geometry(infile)[1]
        bad <- fread(yacfile, sep = "\t") %>% nrow()
        flog.info("Finished looking for chimeras in %s. " %p%
                  "Removing %d/%d as chimeric.",
                  infile, bad, nseq)
        return(data.table(id = files$id[i], file = infile,
                          before = nseq, after = nseq - bad,
                          chimeric = bad))
    }) %>% rbindlist()

    return(stats)
}


#' Removes chimeric reads from amplicon sequencing data.
#'
#' It is recommended to run this step after preprocessing which will be a bit
#' more efficient.
#'
#' @param object An experiment data table as returned by
#'  \code{\link{find_read_files}} or a worflow object.
#' @param ... A configuration as returned by
#'  \code{\link{config_chimera}}.
#' @return A list containing the workflow results:
#' \describe{
#'   \item{passed_reads}{How many reads were kept in each step. Rows are
#'     samples and columns are workflow steps.}
#'   \item{files}{Preprocessed sequencing files list.}
#' }
#' @export
#'
#' @importFrom Biostrings fastq.geometry
remove_chimeras <- function(object, ...) {
    files <- get_files(object)
    files <- copy(files)
    config <- config_parser(list(...), config_chimera)
    if (!"run" %in% names(files)) {
        files[, "run" := "all"]
    }
    paired <- "reverse" %in% names(files)

    if (!dir.exists(config$out_dir)) {
        dir.create(config$out_dir, recursive = TRUE)
    }

    flog.info("Removing chimeras in %d %s-end samples...",
              nrow(files), ifelse(paired, "paired", "single"))
    passed_files <- copy(files)
    passed_files$forward <- file.path(config$out_dir,
                                      basename(files$forward))
    stats <- filter_chimeras(files$forward, passed_files$forward, config)
    if (paired) {
        passed_files$reverse <- file.path(config$out_dir,
                                          basename(files$reverse))
        s <- filter_chimeras(files$reverse, passed_files$reverse, config)
        stats <- rbind(stats, s)
    }
    flog.info("%.3g/%.3g (%.2f%%) reads passed preprocessing.",
              stats[, sum(after)],
              stats[, sum(before)],
              stats[, 100 * mean(after / before)])
    artifact <- list(
        files = passed_files,
        passed = stats,
        steps = c(object[["steps"]], "removes_chimeras")
    )
    return(artifact)
}
