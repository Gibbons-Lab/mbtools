# Copyright 2018 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.
#
# Tools for lineage/taxonomy calling of metagenomic data.


#' The SLIMM taxonomy databases
#' @importFrom data.table data.table
#' @export
slimm_files <- data.table(
    url = "http://ftp.mi.fu-berlin.de/pub/dadi/slimm/ABVF_SP_CMP_genomes.sldb",
    target = "refs/ABVF_SP_CMP_genomes.sldb",
    description = "SLIMM taxonomy mapping DB"
)

clean_taxa_names <- function(dt) {
    flog.info("Fixing taxa names.")
    dt[genus == "unknown_genus" & grepl("\\[.+\\]", species),
              c("genus", "species") := list(
                  tstrsplit(species, "\\[|\\]")[[2]],
                  gsub("\\[|\\]", "", species)
              )]
    dt[, "species" := gsub("Candidatus\\s", "", species)]
    dt[grepl("unknown|unclassified", genus), "genus" := NA]
    dt[grepl("unknown|unclassified", species), "species" := NA]
    return(dt)
}


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
    threads = getOption("mc.cores", 1),
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
    } else {
        reports <- conf$reports
    }
    apfun <- parse_threads(conf$threads)

    alignments[, "id" := tstrsplit(basename(alignment),
                                  ".", fixed = TRUE)[[1]]]
    flog.info("running SLIMM on %d alignments with database %s.",
              nrow(alignments), conf$database)
    ecodes <- apfun(1:nrow(alignments), function(row) {
        al <- alignments$alignment[row]
        id <- alignments$id[row]
        flog.info("[%s] Starting SLIMM...", id)
        repdir <- file.path(reports, id, "")
        if (!file.exists(repdir)) {
            dir.create(repdir)
        }
        ecode <- system2(
            "slimm",
            args = c("-w", conf$bin_width,
                     "-r", conf$rank,
                     "-ac", conf$relative_cutoff, "-co",
                     "-o", repdir, conf$database, al),
            stdout = file.path(repdir, "slimm.log"),
            stderr = file.path(repdir, "slimm.log"))
        flog.info("[%s] Finished running SLIMM.", id)
        return(ecode)
    })

    if (any(unlist(ecodes) != 0)) {
        paste0("slimm terminated with an error, logs can be found in ",
               file.path(reports, "[ID]", "slimm.log")) %>% stop()
    }

    flog.info("Parsing abundances on rank `%s`.", conf$rank)
    abundance <- apfun(
        file.path(reports, alignments$id,
                  paste0(alignments$id, "_profile.tsv")),
        read_slimm) %>% rbindlist() %>% clean_taxa_names()

    flog.info(paste("Estimating read lengths from a sample of",
                    "100 reads per alignment."))
    rlens <- apfun(alignments$alignment, read_length) %>% as.numeric()
    names(rlens) <- alignments$id
    flog.info("Estimated median read length is %d, range is [%d, %d].",
              median(rlens, na.rm = TRUE), min(rlens, na.rm = TRUE),
              max(rlens, na.rm = TRUE))

    flog.info("Parsing coverage profiles with bin width of %dbp.",
              conf$bin_width)
    coverage <- apfun(
        file.path(reports, alignments$id,
                  paste0(alignments$id, "_uniq_coverage2.tsv")),
        read_slimm_coverage, conf$bin_width) %>%
        rbindlist() %>% clean_taxa_names()

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
#' @param profile Where the slimm output is stored.
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
        dt[, "bin_width" := bin_width]
        dt[, "reads" := list(list(reads = cv))]
        dt[, c("genbank", "strain", "species", "genus", "family",
               "order", "class", "phylum", "kingdom") := as.list(anns)]
    }) %>% rbindlist()
    co[, "id" := strsplit(basename(cofile), "_")[[1]][1]]
    co[, "length" := length(reads[[1]]) * bin_width, by = "genbank"]
    return(co)
}
