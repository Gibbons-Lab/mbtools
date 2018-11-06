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
#' @param alignments Whether to keep the alignment. If not NA should be a
#'  string indicating the path to the output bam file.
#' @param threads How many threads to use for mapping.
#' @return A numeric vector with two entries. The number of sequences after
#'  filtering (non-mapped), and the number of removed sequences (mapped).
#' @examples
#'  NULL
#'
#' @export
#' @importFrom data.table tstrsplit
#' @importFrom digest digest
remove_reference <- function(reads, out, reference, index=NA, alignments=NA,
                             threads=3) {
    paired <- length(reads) == 2 & !any(is.na(reads))
    if (is.na(alignments)) {
        name <- digest(reads, "md5")
        alignment_file <- file.path(out, paste0(name, ".bam"))
    } else {
        flog.info("Saving alignment in %s.", keep_bam)
        alignment_file <- keep_bam
    }

    flog.info("Aligning reads to %s...", reference)
    if (!paired) {
        system2("minimap2", c("-t", threads, "-ax", "sr", reference, reads[1],
                              "2> /dev/null | samtools view -bS - >",
                              alignment_file))
    } else {
        system2("minimap2", c("-t", threads, "-ax", "sr", reference, reads[1],
                              reads[2], "2> /dev/null | samtools view -bS - >",
                              alignment_file))
    }

    flog.info("Getting hits and saving filtered reads to %s.", out)
    hits <- as.data.table(read_bam(alignment_file))
    ref_ids <- unique(hits$qname)
    if (is.na(alignments)) {
        unlink(alignment_file)
    }
    new_files <- file.path(out, basename(reads))
    dir.create(out, showWarnings = FALSE)

    streams <- reads
    if (!is.na(index)) {
        streams <- append(streams, index)
        new_files <- append(new_files,
                            file.path(out, basename(index)))
    }

    counts <- sapply(1:length(streams), function(i) {
        reads <- readFastq(streams[i])
        n <- length(reads)
        ids <- sub("/\\d+$", "", as.character(id(reads)))
        ids <- tstrsplit(ids, " ", fixed = TRUE)[[1]]
        rem <- !(ids %in% ref_ids)
        writeFastq(reads[rem], new_files[i])
        c(n, length(reads[rem]))
    })[, 1]
    flog.info("%d/%d reads passed filtering (%.2f%%).",
              counts[2], counts[1], 100 * counts[2] / counts[1])

    return(list(reads = counts[1],
                removed = counts[1] - counts[2],
                counts = count_hit(hits)))
}


#' Filter a set of reference sequences from the data set.seed
#'
#' This will also return a data.table containing the counts of reference
#' sequences for each sample.
#'
#' @param reads A data frame or data table containing the read files. Can be
#'  generated with \link{\code{find_illumina}} for instance.
#' @param out The folder where to store filtered reads. Should be empty as
#'  files **will be overwritten**.
#' @param reference Fasta file (can be gzipped) containing the reference
#'  DNA sequences.
#' @param alignments Optional folder in which to store the alignments.
#' @param threads How many threads to use for mapping.
#' @export
filter_reference <- function(reads, out, reference, alignments = NA,
                             threads = 3) {
    paired <- "reverse" %in% names(reads)
    dir.create(out, showWarnings = FALSE)
    threads <- ceiling(threads / 3)
    counts <- mclapply(1:nrow(reads), function(i) {
        row <- reads[i]
        flog.info("Processing %s on lane %d.", row[, id],
                  as.numeric(row[, lane]))
        r <- if (paired) row[, .(forward, reverse)] else row[, forward]
        if (is.na(alignments)) {
            aln <- NA
        } else {
            name <- strsplit(basename(row[, id]), ".", fixed = TRUE)[[1]]
            aln <- file.path(alignments, paste0(basename(name, ".bam")))
        }
        res <- remove_reference(r, out, reference, alignments = aln,
                                threads = 3)
        res$counts[, "id" := row[, id]]
        if (!is.na(row[, lane])) {
            res$counts[, "lane" := row[, lane]]
        }
        return(res$counts)
    }, mc.cores = threads)
    flog.info("Merging hit tables.")
    return(rbindlist(counts))
}
