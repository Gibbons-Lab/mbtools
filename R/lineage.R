# Copyright 2018 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.
#
# Tools for lineage/taxonomy calling of metagenomic data.

#' Build a configuration for the SLIMM workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the long read alignment
#'  workflow.
#' @export
#' @examples
#'  config <- config_slimm(bin_width = 2000)
config_slimm <- config_builder(list(
    reports = NULL,
    bin_width = 1000,
    relative_cutoff = 0,
    database = "refs/ABVF_SP_CMP_genomes.sldb",
    threads = TRUE,
    rank = "species"
))

#' Quantifies abundances  and coverage using SLIMM.
#'
#' @param object An alignment table or workflow artifact that includes
#'  alignments.
#' @param ... Either additional arguments specifying a config or a config
#'  created with \code{\link{config_slimm}}.
#' @return A SLIMM artifact containing the lineage calls for each sample and
#'  taxon and the coverage.
#' @examples
#'  NULL
#'
#' @export
slimm <- function(object, ...) {
    alignments <- get_alignments(object)
    conf <- config_parser(list(...), config_slimm)
    if (is.null(conf$reports)) {
        reports <- tempdir()
    } else if (!dir.exists(conf$reports)) {
        dir.create(conf$reports)
        reports <- conf$reports
    }
    apfun <- parse_threads(conf$threads)

    flog.info("running SLIMM on %d alignments with database %s.",
              nrow(alignments), conf$database)
    ecodes <- apfun(as.character(alignments$alignment), function(al) {
        id <- strsplit(basename(al), ".", fixed = TRUE)[[1]]
        flog.info("[%s] Starting SLIMM...", id)
        ecode <- system2(
            "slimm",
            args = c("-w", conf$bin_width,
                     "-r", conf$rank,
                     "-ac", conf$relative_cutoff, "-co",
                     "-o", file.path(reports, ""), conf$database, al),
            stdout = file.path(reports, "slimm.log"),
            stderr = file.path(reports, "slimm.log"))
        flog.info("[%s] Finished running SLIMM.", id)
        return(ecode)
    })

    if (any(ecodes != 0)) {
        paste0("slimm terminated with an error, logs can be found in ",
               file.path(reports, "slimm.log")) %>% stop()
    }

    flog.info("Parsing abundances on rank `%s`.", conf$rank)
    abundance <- apfun(
        file.path(reports, paste0(alignments$id, "_profile.tsv")),
        read_slimm) %>% rbindlist()

    flog.info("Parsing coverage profiles with bin width of %dbp.",
              conf$bin_width)
    coverage <- apfun(
        file.path(reports, paste0(alignments$id, "_uniq_coverage2.tsv")),
        read_slimm_coverage, conf$bin_width) %>% rbindlist()
    artifact <- list(
        alignments = alignments,
        abundance = abundance,
        coverage = coverage,
        steps = c(object[["steps"]], "slimm")
    )
    return(artifact)
}


#' Reads SLIMM output to a data table.
#'
#' @param reports Where the slimm output is stored.
#' @return A data table counting reads and relative abundances for all found
#'  bacteria for several taxonomic ranks.
#' @examples
#'  NULL
#'
#' @importFrom data.table fread rbindlist :=
read_slimm <- function(profile) {
    dt <- fread(profile, header = TRUE,
                col.names = c("rank", "taxa_id", "lineage",
                              "relative", "reads"),
                sep = "\t")
    id <- strsplit(basename(profile), "_|\\.")[[1]][1]
    dt[, "id" := id]
    dt[, "relative" := reads / sum(reads)]
    dt[, c("kingdom", "phylum", "class", "order", "family",
           "genus", "species") := tstrsplit(gsub("\\w__", "", lineage),
                                            "|", fixed = TRUE)]
    dt[, "lineage" := NULL]
    return(dt)
}

#' Reads SLIMM coverage to a data table.
#'
#' @param cofile A SLIMM coverage file.
#' @param bin_width The width of the bin in base pairs.
#' @return A data table containing the taxonomy and coverages.
#' @examples
#'  NULL
read_slimm_coverage <- function(cofile, bin_width) {
    co <- readLines(cofile) %>% strsplit(",")
    co <- lapply(co, function(l) {
        anns <- l[1:9]
        cv <- as.numeric(l[10:length(l)])
        dt <- data.table()
        dt[, "start" := bin_width * (seq_along(cv) - 1) + 1]
        dt[, "end" := bin_width * seq_along(cv)]
        dt[, "reads" := cv]
        dt[, c("genbank", "strain", "species", "genus", "family",
               "order", "class", "phylum", "kingdom") := as.list(anns)]
    }) %>% rbindlist()
    co[, "id" := strsplit(basename(cofile), "_")[[1]][1]]
    co[, "length" := .N * bin_width, by = "genbank"]
    return(co)
}
