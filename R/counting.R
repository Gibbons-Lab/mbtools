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
count_alns <- function(alignments, txlengths, file, method="em",
                       maxit=1000, cutoff=0.01, counts=TRUE) {
    aln <- as.data.table(alignments)
    aln[, seqnames := factor(as.character(seqnames))]
    efflengths <- effective_lengths(aln[, txlengths[levels(seqnames)]],
                                    aln[, width])
    flog.info(paste("[%s] %d transcripts. Confidence interval for effective",
                    "lengths: [%.2f, %.2f]."),
              file, aln[, length(levels(seqnames))],
              quantile(efflengths, 0.025), quantile(efflengths, 0.975))
    if (nrow(aln) == 0) {
        return(data.table(seqnames=character(), counts=integer()))
    }
    if (counts) {
        libsize <- aln[, uniqueN(qname)]
    } else {
        libsize = 1e6  # return TPM
    }
    if (method == "naive") {
        aln <- aln[order(-mapq), .SD[1], by="qname"]
        counts <- aln[, .(counts = .N), by="seqnames"]
        names(counts)[1] <- "transcript"
        counts[, counts := counts / efflengths[transcript]]
        counts[, counts := counts / sum(counts) * libsize]
    } else {
        aln[, seqnames := factor(seqnames)]
        aln[, qname := factor(qname)]
        txids <- aln[, as.integer(seqnames) - 1]
        txnames <- aln[, levels(seqnames)]
        rids <- aln[, as.integer(qname) - 1]
        em_result <- em_count(cbind(txids, rids), efflengths,
                              length(txnames), max(rids) + 1, maxit, cutoff)
        flog.info(paste("[%s] Used %d EM iterations on %d equivalence classes.",
                        "Last relative change was %.4f."),
                  file, em_result$iterations, em_result$num_ecs,
                  max(em_result$change))
        counts <- data.table(transcript = txnames,
                             counts = round(em_result$p * libsize))
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
#' @param tags Additional tags to read in the BAM file.
#' @param threads Number of parallel processes.
#' @param method The counting method. Can be either "naive" for assigning each
#'  read to any transcript with the highest mapping score or "em" to resolve
#'  multimapping with an expectation maximzation algorithm.
#' @param maxit Only for method="em". Maximum iterations of EM algorithm.
#' @param cutoff Only for method="em". Stop EM if maximum relative change in
#'  transcript abundance is not at least this value. For instance a value of
#'  0.01 (default) means that the EM algorithm is stopped if transcript
#'  abundances change less than 1% between iterations.
#' @param counts Whether to return counts. If FALSE returns transcripts per
#'  million.
#' @return A data.table with sequence names, counts and sample name.
#'
#' @export
#' @importFrom data.table tstrsplit
count_hits <- function(alignment_files, reference, threads = 1,
                       method = "em", maxit = 1000, cutoff = 0.01,
                       counts = TRUE) {
    flog.info("Getting transcript lengths from %s...", reference)
    fasta_index <- fasta.index(reference)[, c("desc", "seqlength")]
    txlengths <- fasta_index$seqlength
    names(txlengths) <- gsub("\\s.+", "", fasta_index$desc)
    flog.info("Normalized IDs. Starting counting...")
    counts <- mclapply(alignment_files, function(file) {
        bam <- read_bam(file)
        flog.info("[%s] Read %d alignments.", file, nrow(bam))
        cn <- count_alns(bam, txlengths, file=file)
        cn[, "sample" := strsplit(basename(file), ".bam")[[1]][1]]
        return(cn)
    }, mc.cores=threads)

    return(rbindlist(counts))
}
