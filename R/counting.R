# Helpers for counting transcript reads

#' Read alignments from a BAM file.
#'
#' @param path The file path to the BAM.
#' @param tags Additional tags to read from the BAM file.
#' @return The alignments in a GAlignments object.
#'
#' @export
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools ScanBamParam
read_bam <- function(path, tags = character(0)) {
    bam <- readGAlignments(path, param=ScanBamParam(
        what=c("qname", "mapq"), tag=tags))
    return(bam)
}

#' @useDynLib mbtools, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats quantile
count_alns <- function(alignments, txlengths, file, method="em",
                       maxit=1000, cutoff=0.01, tpm=FALSE) {
    aln <- as.data.table(alignments)
    aln[, seqnames := factor(as.character(seqnames))]
    if (is.null(txlengths)) {
        efflengths <- rep(1, aln[, length(levels(seqnames))])
    } else {
        efflengths <- effective_lengths(aln[, txlengths[levels(seqnames)]],
                                        aln[, width])
    }
    names(efflengths) <- aln[, txlengths[levels(seqnames)]]
    flog.info(paste("[%s] %d transcripts. Confidence interval for effective",
                    "lengths: [%.2f, %.2f]."),
              file, aln[, length(levels(seqnames))],
              quantile(efflengths, 0.025), quantile(efflengths, 0.975))
    if (nrow(aln) == 0) {
        return(data.table(transcript=character(), counts=integer(),
                          effective_length=integer()))
    }
    if (tpm) {
        libsize <- 1e6
    } else {
        libsize <-  aln[, uniqueN(qname)]  # return counts
    }
    if (method == "naive") {
        aln <- aln[order(-mapq), .SD[1], by="qname"]
        counts <- aln[, .(counts = .N), by="seqnames"]
        names(counts)[1] <- "transcript"
        counts[, counts := counts / efflengths[transcript]]
        counts[, counts := counts / sum(counts) * libsize]
        counts[, effective_length := efflengths[transcript]]
    } else {
        aln[, seqnames := factor(seqnames)]
        aln[, qname := factor(qname)]
        txids <- aln[, as.integer(seqnames) - 1]
        txnames <- aln[, levels(seqnames)]
        rids <- aln[, as.integer(qname) - 1]
        em_result <- em_count(cbind(txids, rids), efflengths,
                              length(txnames), max(rids) + 1, maxit,
                              cutoff, cutoff)
        flog.info(paste("[%s] Used %d EM iterations on %d equivalence classes.",
                        "Last max. abs. change was %.2g."),
                  file, em_result$iterations, em_result$num_ecs,
                  max(em_result$change))
        counts <- data.table(transcript = txnames,
                             counts = em_result$p,
                             effective_length = efflengths)
        if (tpm) {
            counts[, counts := em_result$p / (max(rids) + 1) * libsize]
            counts[counts < 1e-8, counts := 0]
        } else {
            counts[counts < cutoff, counts := 0]
        }
        counts <- counts[counts > 0]
    }
    return(counts)
}


#' Count alignment hits to a reference database.
#'
#' This will correct for effective transcript lengths as done by almost any
#' good tool those days. So the returned counts are not correlated with
#' transcript lengths. By default an expectation maximization algorithm is
#' used to resolve multiple mappings of one read to many transcripts which
#' pretty much always happens in metagenomics data sets. The optimized
#' likelihodd function is very similar to the one in kallisto
#' (https://doi.org/10.1038/nbt.3519).
#'
#' @param alignment_files Paths to BAM files.
#' @param reference Path to the reference FASTA file. Used to get transcript
#'  lengths. If NA assumes constant length transcripts and will not correct
#'  for transcript lengths.
#' @param threads Number of parallel processes.
#' @param method The counting method. Can be either "naive" for assigning each
#'  read to any transcript with the highest mapping score or "em" to resolve
#'  multimapping with an expectation maximzation algorithm.
#' @param maxit Only for method="em". Maximum iterations of EM algorithm.
#' @param cutoff Only for method="em". Stop EM if maximum relative change in
#'  transcript abundance is not at least this value. For instance a value of
#'  0.01 (default) means that the EM algorithm is stopped if transcript
#'  abundances change less than 1\% between iterations.
#' @param tpm Whether to return counts or transcripts per
#'  million. If FALSE returns counts.
#' @return A data.table with transcript names, counts, effective transcript
#'  length and sample name.
#'
#' @export
#' @importFrom data.table tstrsplit
count_transcripts <- function(alignment_files, reference=NA, threads = 1,
                              method = "em", maxit = 1000, cutoff = 0.01,
                              tpm = FALSE) {
    if (is.na(reference)) {
        flog.info(paste("No reference given so assuming constant length",
                        "transcripts. Starting counting..."))
        txlengths <- NULL
    } else {
        flog.info("Getting transcript lengths from %s...", reference)
        fasta_index <- fasta.index(reference)[, c("desc", "seqlength")]
        txlengths <- fasta_index$seqlength
        names(txlengths) <- gsub("\\s.+", "", fasta_index$desc)
        flog.info("Normalized IDs. Starting counting...")
    }
    counts <- mclapply(alignment_files, function(file) {
        bam <- read_bam(file)
        flog.info("[%s] Read %d alignments.", file, length(bam))
        cn <- count_alns(bam, txlengths, file=file)
        cn[, "sample" := strsplit(basename(file), ".bam")[[1]][1]]
        return(cn)
    }, mc.cores = threads)

    return(rbindlist(counts))
}
