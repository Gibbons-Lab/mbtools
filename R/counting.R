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
                       maxit = 1000, cutoff = 0.01, tpm = FALSE,
                       ecs = FALSE, weighting = FALSE) {
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
        } else {
            aln <- aln[order(-mapq), .SD[1], by = "qname"]
            flog.info(paste("[%s] Missing or uncomplete alignments score.",
                            "Assigning reads by mapping quality only."), file)
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
        if (weighting & aln[, !any(is.na(AS))]) {
            weights <- aln[, AS - min(AS) + 1]
        } else {
            weights <- rep(1, nrow(aln))
            if (weighting) {
                flog.info(paste("[%s] Missing or incomplete alignment score.",
                                "Not using weighting."), file)
            }
        }
        em_result <- em_count(unique(cbind(refids, rids)), efflengths, weights,
                              length(refnames), max(rids) + 1, maxit,
                              cutoff, cutoff ^ 2)
        flog.info(paste("[%s] Used %d EM iterations on %d equivalence classes.",
                        "Last max. abs. change was %.2g. Database concordance",
                        "is %.2f%%."),
                  file, em_result$iterations, length(em_result$ecs),
                  max(em_result$change),
                  100 * (1 - em_result$unobserved / (max(rids) + 1)))
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
        counts <- rbind(counts, data.table(reference = NA,
                                           counts = em_result$unobserved,
                                           effective_length = NA))
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
config_count <- config_builder(list(
        reference = NA,
        threads = getOption("mc.cores", 1),
        method = "em",
        maxit = 10000,
        cutoff = 0.01,
        tpm = FALSE,
        weights = FALSE
    ))

#' Count alignment hits to a reference database.
#'
#' This will correct for effective reference lengths as done by almost any
#' good tool those days. So the returned counts are not correlated with
#' feature lengths. By default an expectation maximization algorithm is
#' used to resolve multiple mappings of one read to many references which
#' pretty much always happens in metagenomics data sets. The optimized
#' likelihood function is very similar to the one in kallisto
#' (https://doi.org/10.1038/nbt.3519).
#'
#' Note that for the EM method there will be a NA reference reported which
#' corresponds to the approximate abundance of references not contained in the
#' database.
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
    alignments <- get_alignments(object)

    apfun <- parse_threads(config$threads)
    counts <- apfun(alignments$alignment, function(file) {
        bam <- read_bam(file)
        flog.info("[%s] Read %d alignments.", file, length(bam))
        cn <- count_alns(bam, reflengths, file = file, method = config$method,
                         weighting = config$weights, maxit = config$maxit,
                         cutoff = config$cutoff)
        cn[, "sample" := strsplit(basename(file), ".bam")[[1]][1]]
        return(cn)
    })

    artifact <- list(
        alignments = object$alignments,
        counts = rbindlist(counts),
        steps = c(object[["steps"]], "count_references")
    )
    return(artifact)
}
