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
        cat("No genome file given, downloading hg19...\n")
    }
    base <- basename(genome_file)
    loc <- file.path(where, base)
    download.file(HS_GENOME, loc)
    fileext <- strsplit(base, "\\.")[[1]]
    if (fileext[2] == "zip") {
        cat("Extracting...\n")
        unzip(loc, exdir = fileext[1])
    }
}

#' Removes sequence contamination from another target organism (for instance
#' human sequences).
#'
#' This function uses bowtie2 to align and identify hits.
#' This requires a pre-built index, for instance by using \link{download_index}.
#'
#' @param reads A character vector containing the read files in fastq format.
#' @param out A folder to which to save the filtered fastq files.
#' @param index Additional barcode file that should be filtered as well.
#' @param organism The base name for the index. Defaults to human (hg19).
#' @param where Where to find the previously generated index files.
#' @return A numeric vector with two entries. The number of sequences after
#'  filtering (non-human), and the number of removed sequences (human).
#' @examples
#'  NULL
#'
#' @export
#' @importFrom GenomicAlignments readGAlignments seqnames
remove_organism <- function(reads, index=NA, out, organism = "hg19",
                         where = "./bowtie2") {
    paired <- length(reads) == 2 & !any(is.na(reads))
    bowtie2_index <- file.path(where, organism)

    cat("Finding human sequences...")
    if (!paired) {
        system2("bowtie2", c("-x", bowtie2_index, "-U", reads[1],
                             "2>", file.path(out, "bowtie2.log"), "|",
                             "samtools", "view", "-bS", "-", ">",
                             file.path(out, "filter.bam")))
    } else {
        system2("bowtie2", c("-x", bowtie2_index, "-1", reads[1],
                             "-2", reads[2], "2>",
                             file.path(out, "bowtie2.log"), "|", "samtools",
                             "view", "-bS", "-", ">",
                             file.path(out, "filter.bam")))
    }


    hits <- readGAlignments(file.path(out, "filter.bam"), use.names = TRUE)
    human_ids <- unique(names(hits))
    unlink(file.path(out, "filter.bam"))
    unlink(file.path(out, "bowtie2.log"))
    new_files <- file.path(out, paste0(basename(reads), ".gz"))
    dir.create(out, showWarnings = FALSE)

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
            rem <- !(ids %in% human_ids)
            writeFastq(reads[rem], new_files[i])
            n <- n + length(reads[rem])
        }
        n
    }, 0)
    lapply(streams, close)

    return(c(reads = counts[1], removed = length(human_ids)))
}
