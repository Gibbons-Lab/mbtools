# Copyright 2018 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.
#
# Tools for lineage/taxonomy calling of metagenomic data.


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
slimm <- function(alignments, slimm_db, reports = NULL) {
    if (!all(alignments$success)) {
        stop("some alignments were not successful!")
    }
    if (is.null(reports)) {
        reports <- tempdir()
    }

    flog.info("running SLIMM on %d alignments with database %s.",
              nrow(alignments), slimm_db)
    ecodes <- pbsapply(as.character(alignments$alignment), function(al) {
        ecode <- system2(
            "slimm",
            args = c("-m", slimm_db, "-o", file.path(reports, ""), al),
            stdout = file.path(reports, "slimm.log"),
            stderr = NULL)
        return(ecode)
    })

    if (any(ecodes != 0)) {
        paste0("slimm terminated with an error, logs can be found in ",
               file.path(reports, "slimm.log")) %>% stop
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
    flog.info("Summarizing results")
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
