# Copyright 2018 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.


#' Align nanopore reads to a 16S reference database.
#'
#' @param read_files Paths to fastq files (can be gzipped).
#' @param ref Path to the reference fasta (can be gzipped).
#' @param alignments_folder Where to store the alignments.
#' @param log_file The log file used to store status messages.
#' @param threads Number of threds used for alignment.
#' @param secondary Return at most that many secondary alignments.
#' @return A data.table with sample name and alignment.
#'
#' @export
align_nanopore <- function(read_files, ref, alignments_folder = "./alignments",
                           log_file="minimap2.log", threads = 4,
                           secondary = 50) {
    if (!dir.exists(alignments_folder)) {
        dir.create(alignments_folder)
    }
    flog.info("Aligning reads to 16S references")
    alns <- lapply(read_files, function(file) {
        flog.info("Aligning %s...", file)
        base <- strsplit(basename(file), ".fa")[[1]][1]
        out_path <- file.path(alignments_folder,
                              paste0(base, ".bam"))
        args <- c("-acx", "map-ont", "-t", threads, "-N", secondary, ref, file)
        args <- append(args, c(paste0("2>", log_file), "|", "samtools",
                            "view", "-bS", "-", ">", out_path))
        success <- system2("minimap2", args = args)
        return(data.table(sample = base, alignment = out_path,
               success = success == 0))
    })
    alns <- rbindlist(alns)
    if (alns[, any(!success)]) {
        stop("At least one alignment failed!")
    }
    alns[, success := NULL]
    return(alns)
}
