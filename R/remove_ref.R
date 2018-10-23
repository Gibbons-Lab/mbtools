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
#'  Can be generated using \link{\code{find_illumina}}.
#' @param out A folder to which to save the filtered fastq files.
#' @param index Additional barcode file that should be filtered as well. Can
#'  be used to filter multiplexed samples.
#' @param reference Path to a fasta file (can be gzipped) that contains the
#'  sequences to filter. Can be a genome or transcripts.
#' @param keep_bam Whether to keep the alignment. If not FALSE should be a
#'  string indicating the path to the output bam file.
#' @return A numeric vector with two entries. The number of sequences after
#'  filtering (non-mapped), and the number of removed sequences (mapped).
#' @examples
#'  NULL
#'
#' @export
#' @importFrom GenomicAlignments readGAlignments seqnames
remove_reference <- function(reads, out, reference, index=NA, keep_bam=FALSE) {
    paired <- length(reads) == 2 & !any(is.na(reads))
    if (keep_bam == FALSE) {
        alignment_file <- file.path(out, "filtered.bam")
    } else {
        flog.info("Saving alignment in %s.", keep_bam)
        alignment_file <- keep_bam
    }

    flog.info("Aligning reads to %s...", index)
    if (!paired) {
        system2("minimap2", c("-ax", "sr", reference, reads[1], "2> /dev/null",
                              "|", "samtools", "view", "-bS", "-", ">",
                              alignment_file))
    } else {
        system2("minimap2", c("-ax", "sr", reference, reads[1], reads[2],
                              "2> /dev/null",
                              "|", "samtools", "view", "-bS", "-", ">",
                              alignment_file))
    }

    flog.info("Getting hits and saving filtered reads to %s.", out)
    hits <- readGAlignments(alignment_file, use.names = TRUE)
    ref_ids <- unique(names(hits))
    if (keep_bam == FALSE) {
        unlink(alignment_file)
    }
    new_files <- file.path(out, paste0(basename(reads), ".gz"))
    dir.create(out, showWarnings = FALSE)

    streams <- reads
    if (!is.na(index)) {
        streams <- append(streams, index)
        new_files <- append(new_files,
                            file.path(out, paste0(basename(index), ".gz")))
    }

    counts <- sapply(1:length(streams), function(i) {
        reads <- readFastq(streams[i])
        n <- length(reads)
        ids <- sub("/\\d+$", "", as.character(id(reads)))
        rem <- !(ids %in% ref_ids)
        writeFastq(reads[rem], new_files[i], mode = "a")
        c(n, length(reads[rem]))
    })[, 1]
    flog.info("%d/%d reads passed filtering.", counts[2], counts[1])

    return(c(reads = counts[1], removed = counts[1] - counts[2]))
}
