# Copyright 2016-2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Build a configuration for the demultiplexing workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the demultiplexing workflow.
#' @export
#' @examples
#'  config <- config_demultiplex(n = 1e4)
config_demultiplex <- config_builder(list(
        out_dir = "demultiplexed",
        barcodes = NULL,
        samples = NULL,
        n = 1e5,
        max_edit = 1,
        reverse_complement = TRUE,
        threads = getOption("mc.cores", 1)
))

#' Splits FASTQ files into individual samples.
#'
#' This function will take each of the sequences barcodes and tries to find
#' the corresponding barcode from the reference set. After that it writes new
#' fastq files containing only the reads for a single barcode.
#'
#' @param object An experiment data table as returned by
#'  \code{\link{find_read_files}} or a worflow object.
#' @param ... A configuration as returned by
#'  \code{\link{config_demultiplex}}.
#' @return A list containing the split files and matching statistics.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom Biostrings DNAStringSet reverseComplement
demultiplex <- function(object, ...) {
    config <- config_parser(list(...), config_demultiplex)
    files <- get_files(object)
    if (!"index" %in% names(files)) {
        stop("must specify an index file for each sample :/")
    }
    if (is.null(config$barcodes)) {
        stop("must specify barcodes in configuration :/")
    }
    apfun <- parse_threads(config$threads)

    ref <- DNAStringSet(config$barcodes)
    nref <- length(ref)
    if (nref < 2)
        stop("There is only one sample. Barcode filtering is not necessary!")
    if (is.null(config$samples)) {
        flog.warn("No sample names specified using S1-SN...")
        snames <- paste0("S", 1:length(ref))
    } else {
        snames <- config$samples
        if (length(snames) != nref) {
            stop("Need exactly one sample name for each barcode.")
        }
    }

    if (config$reverse_complement) {
        ref <- c(ref, reverseComplement(ref))
    }

    if (dir.exists(config$out_dir)) {
        unlink(file.path(config$out_dir, "*.fastq.gz"))
    } else dir.create(config$out_dir, recursive = TRUE)
    unmatched <- 0
    multiple <- 0
    read_counts <- rep(0, length(snames))
    names(read_counts) <- snames
    nseq <- 0

    for (i in 1:nrow(files)) {
        flog.info("Processing ID %s...", files[i, id])
        istream <- FastqStreamer(files[i, index], n = config$n)
        rstream <- files[i, lapply(.(forward, reverse), FastqStreamer,
                                   n = config$n)]
        repeat {
            fq <- yield(istream)
            if (length(fq) == 0) break

            ids <- sub("[/\\s].+$", "", id(fq), perl = TRUE)

            hits <- apfun(ref, function(p) srdistance(fq, p)[[1]])
            hits <- do.call(cbind, hits)
            inds <- apply(hits, 1, function(x) {
                scores <- x[x < config$max_edit]
                i <- which(x < config$max_edit)
                # reverse complements
                i[i > nref] <- i - nref
                if (length(i) == 0) return(-1)
                if (length(i) > 1) {
                    i <- i[scores == min(scores)]
                    if (length(i) != 1) return(0)
                }
                return(i)
            })
            flog.info("Processed chunk of size %g. Found %d hits.",
                      config$n, sum(inds > 0))

            for (di in 1:length(rstream)) {
                rfq <- yield(rstream[[di]])
                rids <- sub("[/\\s].+$", "", id(rfq), perl = TRUE)
                if (any(rids != ids)) {
                    sapply(rstream, close)
                    stop("Index file and reads do not match!")
                }

                cn <- apfun(1:nref, function(sid) {
                    filename <- sprintf("%s_S%d_L%d_R%d_001.fastq.gz",
                                        snames[sid], sid, i, di)
                    matched <- rfq[inds == sid]
                    if (length(matched) > 0) {
                        writeFastq(matched,
                                   file.path(config$out_dir, filename),
                                   mode = "a", compress = TRUE)
                    }
                    return(length(matched))
                }) %>% as.numeric()
                if (di == 1) {
                    read_counts <- read_counts + cn
                    unmatched <- unmatched + sum(inds == 0)
                    multiple <- multiple + sum(inds == -1)
                }
            }
            nseq <- nseq + length(fq)
        }
        close(istream)
        sapply(rstream, close)
    }
    read_counts["unmatched"] <- unmatched
    read_counts["multiple"] <- multiple
    artifact <- list(
        files = find_read_files(config$out_dir),
        read_counts = read_counts,
        steps = c(object[["steps"]], "demultiplex")
    )
    return(artifact)
}
