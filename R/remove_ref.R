# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Removes sequence that map to a given reference.
#'
#' This can be used for a variety of applications the most common ones are:
#' \itemize{
#'   \item removing sequences from the host
#'   \item removing ribosomal sequences
#'   \item removing contaminants
#' }
#' This function uses minimap2 to align and identify hits and does not require
#' a prebuilt index.
#'
#' @param reads A character vector containing the read files in fastq format.
#'  Can be generated using \code{\link{find_read_files}}.
#' @param out A folder to which to save the filtered fastq files.
#' @param reference Path to a fasta file (can be gzipped) that contains the
#'  sequences to filter. Can be a genome or transcripts.
#' @param alignments Whether to keep the alignment. If not NA should be a
#'  string indicating the path to the output bam file.
#' @param threads How many threads to use for mapping.
#' @return A numeric vector with two entries. The number of sequences after
#'  filtering (non-mapped), and the number of removed sequences (mapped).
#' @examples
#'  NULL
#'
#' @importFrom data.table tstrsplit
#' @importFrom digest digest
remove_reference <- function(reads, out, reference, alignments=NA,
                             threads=3) {
    paired <- length(reads) == 2 & !any(is.na(reads))
    if (is.na(alignments)) {
        name <- digest(reads, "md5")
        alignment_file <- file.path(out, paste0(name, ".bam"))
    } else {
        alignment_file <- alignments
    }

    flog.info("[%s] Aligning reads to %s, saving in %s.",
              basename(reads[1]), reference, alignment_file)
    if (!paired) {
        system2("minimap2", c("-t", threads, "-ax", "sr", reference, reads[1],
                              "2> /dev/null | samtools view -bS - >",
                              alignment_file))
    } else {
        system2("minimap2", c("-t", threads, "-ax", "sr", reference, reads[1],
                              reads[2], "2> /dev/null | samtools view -bS - >",
                              alignment_file))
    }

    hits <- as.data.table(read_bam(alignment_file))
    ref_ids <- unique(hits$qname)
    if (is.na(alignments)) {
        unlink(alignment_file)
    }
    new_files <- file.path(out, basename(reads))
    dir.create(out, showWarnings = FALSE)

    streams <- reads
    counts <- sapply(1:length(streams), function(i) {
        reads <- readFastq(streams[i])
        n <- length(reads)
        ids <- sub("/\\d+$", "", as.character(id(reads)))
        ids <- tstrsplit(ids, " ", fixed = TRUE)[[1]]
        rem <- !(ids %in% ref_ids)
        writeFastq(reads[rem], new_files[i])
        c(n, length(reads[rem]))
    })[, 1]
    flog.info("[%s] %d/%d reads passed filtering (%.2f%%).",
              basename(reads[1]),
              counts[2], counts[1], 100 * counts[2] / counts[1])

    return(list(reads = counts[1],
                removed = counts[1] - counts[2]))
}


#' Build a configuration for the reference removal workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the DADA2 workflow.
#' @export
#' @examples
#'  config <- config_reference(reference = "refs/mouse.fna.gz")
config_reference <- config_builder(list(
        threads = getOption("mc.cores", 1),
        out_dir = "reference_removed",
        alignment_dir = NA,
        reference = NA
))


#' Filter a set of reference sequences from the data set.seed
#'
#' This will also return a data.table containing the counts of reference
#' sequences for each sample.
#'
#' @param object An experiment data table as returned by
#'  \code{\link{find_read_files}} or a worflow object.
#' @param ... A configuration as returned by
#'  \code{\link{config_reference}}.
#' @return A list with the processed files and removal counts for each sample.
#' @export
filter_reference <- function(object, ...) {
    files <- get_files(object)
    config <- config_parser(list(...), config_reference)
    if (is.na(config$reference)) {
        stop("Must specify a reference to remove in configuration :/")
    }
    paired <- "reverse" %in% names(files)
    dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
    # we will leave 3 threads for minimap2
    threads <- parse_threads(config$threads, FALSE)
    threads <- ceiling(threads / 3)
    flog.info("Actually using %d threads to filter and count.", threads * 3)
    counts <- mclapply(1:nrow(files), function(i) {
        row <- files[i]
        if ("lane" %in% names(row)) {
            lane <- as.numeric(row$lane)
        } else {
            lane <- NA
        }
        flog.info("Processing %s on lane %d.", row[, id], lane)
        r <- if (paired) row[, .(forward, reverse)] else row[, forward]
        if (is.na(config$alignment_dir)) {
            aln <- NA
        } else {
            name <- strsplit(basename(row[, id]), ".", fixed = TRUE)[[1]][1]
            aln <- file.path(config$alignment_dir,
                             paste0(basename(name), ".bam"))
        }
        res <- remove_reference(as.character(r), config$out_dir,
                                config$reference, alignments = aln,
                                threads = 3)
        res$id <- row[, id]
        res$lane <- lane
        return(res)
    }, mc.cores = threads)
    flog.info("Merging hit tables...")
    filtered <- copy(files)
    filtered[, forward := file.path(config$out_dir, basename(forward))]
    if (paired) {
        filtered[, reverse := file.path(config$out_dir, basename(reverse))]
    }
    artifact <- list(
        counts = rbindlist(counts),
        files = filtered,
        steps = c(object[["steps"]], "filter_reference")
    )
    return(artifact)
}
