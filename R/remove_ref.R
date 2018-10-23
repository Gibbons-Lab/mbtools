# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

HS_INDEX <- "ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip"

#' Builds index files for human sequence identification
#'
#' This function has to be run only once to enable analysis.
#'
#' @param where Where to save index files
#' @param genome_file The human genome to be used. NA means download a recent version.
#' @return Nothing.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom utils unzip
download_index <- function(where="./bowtie2", genome_file=NA) {
    dir.create(where, showWarnings = FALSE)

    if (is.na(genome_file)) {
        genome_file <- HS_INDEX
        flog.info("No genome file given, downloading hg19...")
    }
    base <- basename(genome_file)
    loc <- file.path(where, base)
    download.file(HS_GENOME, loc)
    fileext <- strsplit(base, "\\.")[[1]]
    if (fileext[2] == "zip") {
        flog.info("Extracting %s...", base)
        unzip(loc, exdir = fileext[1])
    }
}

#' Removes sequence that map to a given reference.
#'
#' This can be used for a variety of applications the most common ones are:
#' \itemize{
#'   \item removing sequences from the host
#'   \item removing ribosomal sequences
#'   \item removing contaminants
#' }
#' This function uses bowtie2 to align and identify hits.
#' This requires a pre-built index, for instance by using \link{download_index}.
#'
#' @param reads A character vector containing the read files in fastq format.
#'  Can be generated using \link{\code{find_illumina}}.
#' @param out A folder to which to save the filtered fastq files.
#' @param index Additional barcode file that should be filtered as well.
#' @param reference The base name for the index. Defaults to human (hg19).
#' @param where Where to find the previously generated index files.
#' @param keep_bam Whether to keep the alignment. If not FALSE should be a
#'  string indicating the path to the output bam file.
#' @return A numeric vector with two entries. The number of sequences after
#'  filtering (non-mapped), and the number of removed sequences (mapped).
#' @examples
#'  NULL
#'
#' @export
#' @importFrom GenomicAlignments readGAlignments seqnames
remove_reference <- function(reads, index=NA, out, organism = "hg19",
                         where = "./bowtie2", keep_bam=FALSE) {
    paired <- length(reads) == 2 & !any(is.na(reads))
    bowtie2_index <- file.path(where, organism)
    if (keep_bam == FALSE) {
        alignment_file <- file.path(out, "filtered.bam")
    } else {
        flog.info("Saving alignment in %s.", keep_bam)
        alignment_file <- keep_bam
    }

    flog.info("Aligning reads to %s...", index)
    if (!paired) {
        system2("bowtie2", c("-x", bowtie2_index, "-U", reads[1],
                             "2>", file.path(out, "bowtie2.log"), "|",
                             "samtools", "view", "-bS", "-", ">",
                             alignment_file))
    } else {
        system2("bowtie2", c("-x", bowtie2_index, "-1", reads[1],
                             "-2", reads[2], "2>",
                             file.path(out, "bowtie2.log"), "|", "samtools",
                             "view", "-bS", "-", ">",
                             alignment_file))
    }

    flog.info("Getting hits and saving filtered reads to %s.", out)
    hits <- readGAlignments(alignment_file, use.names = TRUE)
    ref_ids <- unique(names(hits))
    if (keep_bam == FALSE) {
        unlink(alignment_file)
    }
    rates <- alignment_rate(file.path(out, "bowtie2.log"))
    unlink(file.path(out, "bowtie2.log"))
    new_files <- file.path(out, paste0(basename(reads), ".gz"))
    dir.create(out, showWarnings = FALSE)
    flog.info("%d/%d reads passed filtering.", rates[1] - rates[2], rates[1])

    streams <- list(f = FastqStreamer(reads[1]))
    if (!is.na(index)) {
        streams$i <- FastqStreamer(index)
        new_files[3] <- file.path(out, paste0(basename(index), ".gz"))
    }
    if (paired) {
        streams$r <- FastqStreamer(reads[2])
    }

    counts <- vapply(1:length(streams), function(i) {
        n <- 0
        repeat {
            reads <- yield(streams[[i]])
            if (length(reads) == 0) break
            ids <- sub("/\\d+$", "", as.character(id(reads)))
            rem <- !(ids %in% ref_ids)
            writeFastq(reads[rem], new_files[i], mode="a")
            n <- n + length(reads[rem])
        }
        n
    }, 0)
    lapply(streams, close)

    return(c(reads = unname(rates[1]), removed = unname(rates[2])))
}
