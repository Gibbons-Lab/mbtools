# Copyright 2018 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

types <- c(i = "as.integer", f = "as.double", A = "as.character",
           Z = "as.character")

parse_annotations <- function(vec) {
    mapped <- list()
    vec <- tstrsplit(as.character(vec), ":")
    bad <- is.na(vec[[3]])
    res <- as.list(vec[[3]][!bad])
    names(res) <- vec[[1]][!bad]
    return(res)
}

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

count_hits <- function(alignments) {
    aln <- as.data.table(alignments)
    aln <- aln[order(-mapq, -AS), .SD[1], by="qname"]
    counts <- aln[, .(counts = .N), by="seqnames"]
    return(counts)
}

#' Count nanopore hits to a 16S reference database.
#'
#' @param alignment_files Paths to BAM files.
#' @return A data.table with sequence names, counts and sample name.
#'
#' @export
count_nanopore <- function(alignment_files) {
    counts <- pblapply(alignment_files, function(file) {
        bam <- read_bam(file, tags=c("AS", "dv"))
        cn <- count_hits(bam)
        cn[, "sample" := strsplit(basename(file), ".fa")[[1]][1]]
        return(cn)
    })

    return(rbindlist(counts))
}

#' Align nanopore reads to a 16S reference database.
#'
#' @param read_files Paths to fastq files (can be gzipped).
#' @param ref Path to the reference fasta (can be gzipped).
#' @param alignments_folder Where to store the alignments.
#' @param log_file The log file used to store status messages.
#' @param threads Number of threds used for alignment.
#' @return A data.table with sequence names, counts and sample name.
#'
#' @export
align_nanopore <- function(read_files, ref, alignments_folder="./alignments",
                           log_file="minimap2.log", threads=4) {
    write("Building index...", file="")
    success <- system2("minimap2", args = c("-x", "map-ont", "-d",
                                            "index.mmi", ref, "2>", log_file))
    if (success != 0) {
        stop("Index build failed.")
    }
    if (!dir.exists(alignments_folder)) {
        dir.create(alignments_folder)
    }
    write("Aligning reads to 16S references...", file="")
    successes <- pbsapply(read_files, function(file) {
        base <- strsplit(basename(file), ".fa")[[1]][1]
        out_path <- file.path(alignments_folder,
                              paste0(base, ".bam"))
        args <- c("-acx", "map-ont", "-t", threads, "index.mmi", file)
        args <- append(args, c(paste0("2>", log_file), "|", "samtools",
                            "view", "-bS", "-", ">", out_path))
        success <- system2("minimap2", args = args)
        return(success)
    })
    if (any(successes != 0)) {
        stop("At least one alignment failed!")
    }
}
