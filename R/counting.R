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
read_bam <- function(path, tags = c("AS", "de")) {
    bam <- readGAlignments(path, param=ScanBamParam(
        what=c("qname", "mapq"), tag=tags))
    return(bam)
}

#' @useDynLib mbtools, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats quantile
count_alns <- function(alignments, reflengths, file, method = "em",
                       maxit = 1000, cutoff = 0.01, tpm = FALSE, ecs = FALSE) {
    aln <- as.data.table(alignments)
    aln[, "seqnames" := factor(as.character(seqnames))]
    if (is.null(reflengths)) {
        efflengths <- rep(1, aln[, length(levels(seqnames))])
    } else {
        efflengths <- effective_lengths(aln[, reflengths[levels(seqnames)]],
                                        aln[, width])
    }
    names(efflengths) <- aln[, reflengths[levels(seqnames)]]
    flog.info(paste("[%s] %d reference seqs. Confidence interval for effective",
                    "lengths: [%.2f, %.2f]."),
              file, aln[, length(levels(seqnames))],
              quantile(efflengths, 0.025), quantile(efflengths, 0.975))
    if (nrow(aln) == 0) {
        return(data.table(reference = character(), counts = integer(),
                          effective_length = integer()))
    }
    if (tpm) {
        libsize <- 1e6
    } else {
        libsize <-  aln[, uniqueN(qname)]  # return counts
    }
    equiv_classes <- NULL
    if (method == "naive") {
        if (aln[, !any(is.na(AS))]) {
            aln <- aln[order(-AS, -mapq), .SD[1], by = "qname"]
            flog.info(paste("[%s] Assigning reads by alignment score",
                            "and mapping quality."), file)
        } else {
            aln <- aln[order(-mapq), .SD[1], by = "qname"]
            flog.info("[%s] Assigning reads by mapping quality.", file)
        }
        counts <- aln[, .(counts = .N), by = "seqnames"]
        names(counts)[1] <- "reference"
        counts[, "counts" := counts / efflengths[reference]]
        counts[, "counts" := counts / sum(counts) * libsize]
        counts[, "effective_length" := efflengths[reference]]
    } else {
        aln[, seqnames := factor(seqnames)]
        aln[, qname := factor(qname)]
        refids <- aln[, as.integer(seqnames) - 1]
        refnames <- aln[, levels(seqnames)]
        rids <- aln[, as.integer(qname) - 1]
        em_result <- em_count(unique(cbind(refids, rids)), efflengths,
                              length(refnames), max(rids) + 1, maxit,
                              cutoff, cutoff)
        flog.info(paste("[%s] Used %d EM iterations on %d equivalence classes.",
                        "Last max. abs. change was %.2g."),
                  file, em_result$iterations, length(em_result$ecs),
                  max(em_result$change))
        equiv_classes <- em_result$ecs
        counts <- data.table(reference = refnames,
                             counts = em_result$p,
                             effective_length = efflengths)
        if (tpm) {
            counts[, counts := em_result$p]
            print(libsize == max(rids) + 1)
            counts[counts < 1e-8, "counts" := 0]
        } else {
            counts[counts < cutoff, "counts" := 0]
        }
        counts <- counts[counts > 0]
    }
    if (ecs) {
        equiv_classes <- lapply(equiv_classes, function(ec) {
            list(
                count = ec[1],
                references = refnames[ec[2:length(ec)] + 1]
            )
        })
        return(list(counts = counts, ecs = equiv_classes))
    }
    return(counts)
}

#' Build a configuration for the transcript counting workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the transcript counting
#'  workflow.
#' @export
#' @examples
#'  config <- config_count(reference = "refs/mouse.fna.gz")
config_count <- function(...) {
    config <- list(
        reference = NA,
        threads = 1,
        method = "em",
        maxit = 1000,
        cutoff = 0.01,
        tpm = FALSE
    )
    args <- list(...)
    for (arg in names(args)) {
        config[[arg]] <- args[[arg]]
    }
    return(config)
}

#' Count alignment hits to a reference database.
#'
#' This will correct for effective treference lengths as done by almost any
#' good tool those days. So the returned counts are not correlated with
#' feature lengths. By default an expectation maximization algorithm is
#' used to resolve multiple mappings of one read to many references which
#' pretty much always happens in metagenomics data sets. The optimized
#' likelihodd function is very similar to the one in kallisto
#' (https://doi.org/10.1038/nbt.3519).
#'
#' @param object An experiment data table as returned by any alignment method
#'  like \code{\link{align_short_reads}} or \code{\link{align_long_reads}} .
#' @param ... A configuration as generated by \code{\link{config_count}}.
#' @return A list containing the used alignments and the transcript counts in
#'  `counts`.
#'
#' @export
#' @importFrom data.table tstrsplit
count_references <- function(object, ...) {
    config <- config_parser(list(...), config_count)
    if (is.na(config$reference)) {
        flog.info(paste("No reference given so assuming constant length",
                        "sequences. Starting counting..."))
        reflengths <- NULL
    } else {
        flog.info("Getting reference lengths from %s...", config$reference)
        fasta_index <- fasta.index(config$reference)[, c("desc", "seqlength")]
        reflengths <- fasta_index$seqlength
        names(reflengths) <- gsub("\\s.+", "", fasta_index$desc)
        flog.info("Normalized IDs. Starting counting...")
    }
    counts <- mclapply(object$alignments$alignment, function(file) {
        bam <- read_bam(file)
        flog.info("[%s] Read %d alignments.", file, length(bam))
        cn <- count_alns(bam, reflengths, file = file, config$method)
        cn[, "sample" := strsplit(basename(file), ".bam")[[1]][1]]
        return(cn)
    }, mc.cores = config$threads)

    artifact <- list(
        alignments = object$alignments,
        counts = rbindlist(counts),
        steps = c(object[["steps"]], "count_references")
    )
    return(artifact)
}
