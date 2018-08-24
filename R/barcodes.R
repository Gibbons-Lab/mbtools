# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Splits FASTQ files into individual samples.
#'
#' This function will take each of the sequences barcodes and tries to find
#' the corresponding barcode from the reference set. After that it writes new
#' fastq files containing only the reads for a single barcode.
#'
#' @param reads A character vector containing the read files in fastq format.
#' @param index The index file containing the demultiplexed barcodes for each
#'  read file.
#' @param out A folder to which to save the split fastq files.
#' @param ref A character vector or DNAStringSet containing the reference
#'  barcodes.
#' @param max_ed Maximum allowed edit distance between the sequenced and
#'  reference barcode.
#' @param n Maximum number of records to read in each iteration.
#' @return A numeric vector containing three entries, where the first defines the
#'  reads that are kept.
#'  \itemize{
#'  \item{The number of reads that could be mapped uniquely.}
#'  \item{The number of reads for which no match was found.}
#'  \item{The number of reads for which more than one reference match was found.}
#'  }
#' @examples
#'  NULL
#'
#' @export
#' @importFrom Biostrings DNAStringSet reverseComplement
split_barcodes <- function(reads, index, out, ref, n=1e5, max_ed=1) {
    if (is.null(names(ref))) snames <- paste0("S", 1:length(ref))
    else snames <- names(ref)

    ref <- DNAStringSet(ref)
    nref <- length(ref)
    if (nref < 2)
        stop("There is only one sample. Barcode filtering is not necessary!")
    ref <- c(ref, reverseComplement(ref))

    istream <- FastqStreamer(index, n = n)
    on.exit(close(istream))

    rstream <- lapply(reads, FastqStreamer, n = n)

    if (dir.exists(out)) {
        unlink(file.path(out, "*.fastq.gz"))
    } else dir.create(out)
    res <- c(unique = 0, multiple = 0, nomatch = 0)
    nseq <- 0

    repeat {
        fq <- yield(istream)
        if (length(fq) == 0) break

        ids <- sub("[/\\s].+$", "", id(fq), perl = TRUE)

        hits <- do.call(cbind, srdistance(fq, ref)) <= max_ed
        inds <- apply(hits, 1, function(x) {
            i <- unique(which(x) %% nref) + 1
            if (length(i) == 0) return(-1)
            if (length(i) > 1) return(0)
            else return(i)
        })

        for (i in 1:length(rstream)) {
            rfq <- yield(rstream[[i]])
            rids <- sub("[/\\s].+$", "", id(rfq), perl = TRUE)
            if (any(rids != ids)) {
                sapply(rstream, close)
                stop("Index file and reads do not match!")
            }
            fn <- basename(reads[i])

            for (sid in 1:nref) {
                writeFastq(rfq[inds == sid], file.path(out,
                    paste0(snames[sid], "_", fn)), "a")
            }
        }
        nseq <- nseq + length(fq)
        res <- res + c(sum(inds > 0), sum(inds == 0), sum(inds < 0))
    }
    sapply(rstream, close)

    return(res)
}
