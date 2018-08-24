# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' @importFrom stringr str_match
alignment_rate <- function(log_file) {
    if (file.exists(log_file)) {
        content <- readChar(log_file, file.info(log_file)$size)
        match <- str_match(content, "(\\d+\\.*\\d*)% overall alignment rate")
        if (nrow(match) > 0) {
            return(as.numeric(match[ ,2])/100)
        }
    }
    return(NULL)
}


#' Aligns metagenomic shotgun sample against a reference
#'
#' This method uses bowtie2 and should not be too sensitive to the used pre-
#' processing.
#'
#' @param reads A data frame containing the paths to the forward and reverse
#'  reads in fastq format. Needs to have at least columns "id" and "forward".
#'  If the samples are paired you also need a column "reverse".
#' @param index_basename Path and basename of the index against which to align.
#' @param threads How many threads to use for bowtie2.
#' @param alignment_folder A folder to which to save log output and the
#'  generated alignment.
#' @param bam Whether to output alignments in BAM format (requires samtools).
#' @param bowtie2_path Path to the bowtie executables.
#' @param samtools_path Path to samtools.
#' @return A data frame mapping the alignments to the sample IDs. It will
#'  contain the following columns:
#'  \itemize{
#'  \item{id}{the id of the sample}
#'  \item{success}{whether alignment executed without error}
#'  \item{log}{the log output for the alignment}
#'  \item{alignment}{the file containing the alignment}
#'  }
#' @examples
#'  NULL
#'
#' @export
#' @importFrom pbapply pbapply pbsapply pblapply
align_bowtie2 <- function(reads, index_basename, threads=1,
                          alignment_folder = "alignments", bam = TRUE,
                          bowtie2_path = NULL, samtools_path = NULL) {
    env <- character()
    if (!is.null(bowtie2_path)) {
        env <- paste0("PATH=", bowtie2_path)
    }
    if (!is.null(samtools_path)) {
        env <- paste0(env, ":", samtools_path)
    }
    paired <- "reverse" %in% names(reads)

    write("Aligning reads to microbial genomes...", file="")
    alignments <- pbapply(reads, 1, function(read) {
        read <- as.list(read)
        log_file <- file.path(alignment_folder, paste0(read$id, ".log"))

        if (file.exists(log_file)) {
            rate <- alignment_rate(log_file)

            if (bam) {
                out_path <- file.path(alignment_folder,
                                      paste0(read$id, ".bam"))
            } else {
                out_path <- file.path(alignment_folder,
                                      paste0(read$id, ".sam"))
            }
            if (!is.null(rate) && file.exists(out_path)) {
                return(data.table(id = read$id, success = TRUE,
                       log = log_file, alignment = out_path, rate = rate))
            }
        }

        args <- c("-x", index_basename)
        if (paired) {
            args <- append(args, c("-1", read$forward, "-2", read$reverse))
        } else {
            args <- append(args, c("-U", read$forward))
        }
        args <- append(args, c("-q", "--no-unal", "--mm", "-p", threads, "-k",
                               "60"))

        if (bam) {
            out_path <- file.path(alignment_folder, paste0(read$id, ".bam"))
            args <- append(args, c(paste0("2>", log_file), "|", "samtools",
                           "view", "-bS", "-", ">", out_path))
        } else {
            out_path <- file.path(alignment_folder, paste0(read$id, ".sam"))
            args <- append(args, c("-S", out_path, "2>", log_file))
        }
        success <- system2("bowtie2", args = args, env = env)
        rate <- NULL
        if (success) {
            rate <- alignment_rate(log_file)
        }

        return(data.table(id = read$id, success = success == 0, log = log_file,
                          alignment = out_path, rate = rate))
    })

    return(rbindlist(alignments))
}


#' Quantifies abundances for the bacteria using SLIMM
#'
#' @param alignments A data frame as output by \code{\link{align_bowtie2}}
#' @param slimm_db Path for the SLIMM data base.
#' @param reports Path where to save the SLIMM reports. Uses a temporary
#'  directory by default.
#' @return Path to the slimm output.
#' @examples
#'  NULL
#'
#' @export
run_slimm <- function(alignments, slimm_db, reports = NULL) {
    if (!all(alignments$success)) {
        stop("some alignments were not successful!")
    }
    if (is.null(reports)) {
        reports <- tempdir()
    }

    write("Running SLIMM...", file="")
    ecodes <- pbsapply(as.character(alignments$alignment), function(al) {
        ecode <- system2("slimm", args=c("-m", slimm_db, "-o",
                         file.path(reports, ""), al),
                         stdout=file.path(reports, "slimm.log"), stderr = NULL)
        return(ecode)
    })

    if (any(ecodes != 0)) {
        stop(paste0("slimm terminated with an error, logs can be found in",
                    file.path(reports, "slimm.log")))
    }

    return(reports)
}


#' Reads SLIMM output to a data table.
#'
#' @param reports Where the slimm output is stored.
#' @return A data table counting reads and relative abundances for all found
#'  bacteria for several taxonomic ranks. It will contain the following columns:
#'  \itemize{
#'  \item{id}{id of the sample}
#'  \item{rank}{the name of the taxonomic rank, for instance "genus"}
#'  \item{name}{the value of the rank, for instance "Escherichia"}
#'  \item{reads}{the number of reads in that rank}
#'  \item{relative}{the relative abundance of the rank in [0,1]}
#'  \item{coverage}{the coverage of the rank}
#'  }
#' @examples
#'  NULL
#'
#' @export
#' @importFrom data.table fread rbindlist :=
read_slimm <- function(reports) {
    write("Summarizing results", file = "")
    tsvs <- list.files(reports, pattern = "_reported.tsv", full.names = TRUE)
    dts <- pblapply(tsvs, function(file) {
        elements <- strsplit(basename(file), "_")[[1]]
        id <- elements[1]
        rank <- elements[2]
        dt <- fread(file, drop = 1, col.names = c("name", "taxid", "reads",
                    "relative", "contributors", "coverage"))
        dt[, "rank" := rank]
        dt[, "id" := id]
        dt[, "relative" := reads / sum(reads)]
        dt
    })

    return(rbindlist(dts))
}
