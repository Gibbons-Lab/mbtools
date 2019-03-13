# Copyright 2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

classified_taxa <- function(feature_table, tax_table) {
    ranks <- colnames(tax_table)
    ranks <- ranks[ranks != "sequence"]  # if hashed
    total <- colSums(feature_table)
    classified <- lapply(ranks, function(r) {
        known <- !is.na(tax_table[, r])
        asvs <- sum(known) / nrow(tax_table)
        reads <- sum(total[known]) / sum(total)
        return(data.table(rank = r, asvs = asvs, reads = reads))
    })
    return(rbindlist(classified))
}

#' Build a configuration for the DADA2 workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the DADA2 workflow.
#' @export
#' @examples
#'  config <- config_denoise(nbases = 1e9)
config_denoise <- function(...) {
    config <- list(
        threads = TRUE,
        nbases = 2.5e8,
        pool = FALSE,
        bootstrap_confidence = 0.5,
        taxa_db = paste0("https://zenodo.org/record/1172783/files/",
                         "silva_nr_v132_train_set.fa.gz?download=1"),
        species_db = paste0("https://zenodo.org/record/1172783/files/",
                            "silva_species_assignment_v132.fa.gz?download=1"),
        hash = TRUE
    )
    args <- list(...)
    for (arg in names(args)) {
        config[[arg]] <- args[[arg]]
    }
    return(config)
}

#' @importFrom dada2 getUniques
getN <- function(x) {
    if (any(class(x) %in% c("dada", "derep", "matrix", "data.frame"))) {
        x <- list(x)
    }
    sapply(x, function(xi) sum(getUniques(xi)))
}

#' Runs a full DADA2 workflow.
#'
#' This performs the following steps: \enumerate{
#'  \item{Learning the error rate for each run separately}
#'  \item{Inferring sequence variants for each run and sample}
#'  \item{Merging of feature tables across runs}
#'  \item{Consensus chimera removal}
#'  \item{Taxa assignment with Naive Bayes}
#'  \item{Species assignment by exact alignment}
#'  \item{Diagnostic plots of the error rates}
#'  }
#' You will usually want to preprocess the read files first with
#' \code{\link{preprocess}}.
#' Depending on your config the sequences might be represented by
#' an MD5 hash. In this case the taxa table has an additional `sequence` column
#' containing the real sequences.
#'
#' @param object An experiment data table as returned by
#'  \code{\link{find_read_files}} or a worflow object.
#' @param config A configuration file as returned by
#'  \code{\link{config_denoise}}.
#' @return A list containing the workflow results:
#' \describe{
#'   \item{feature_table}{Matrix of sequence variant abundances. Samples are
#'    and sequences are columns.}
#'   \item{taxonomy}{Matrix of taxonomy assignments. Rows are sequences and
#'    columns are taxonomy ranks.}
#'   \item{errors}{The error profiles estimated by DADA2. One for each run and
#'    read direction (forward/reverse).}
#'   \item{error_plots}{Plotted error profiles. One for each run and
#'    read direction (forward/reverse).}
#'   \item{passed_reads}{How many reads were kept in each step. Rows are
#'     samples and columns are workflow steps.}
#'   \item{classified}{The proportion of sequence variants that could be
#'        where a specific taxa rank could be classified.}
#' }
#' @export
#'
#' @importFrom dada2 learnErrors dada removeBimeraDenovo
#'  mergeSequenceTables assignTaxonomy addSpecies plotErrors
#' @importFrom digest digest
denoise <- function(object, config) {
    files <- get_files(object)
    files <- copy(files)
    if (!"run" %in% names(files)) {
        files[, "run" := "all"]
    }
    paired <- "reverse" %in% names(files)

    flog.info("Running DADA2 on %d run(s) from a sample of up to %.3g bases.",
              files[, uniqueN(run)], config$nbases)
    errors <- list()
    dada_stats <- list()
    feature_table <- list()
    for (r in files[, unique(run)]) {
        fi <- files[run == r]
        dada_stats[[r]] <- data.table(id = fi$id)
        flog.info("Learning errors for run `%s` (%d samples)...", r, nrow(fi))
        errors[[r]] <- list()
        errors[[r]][["forward"]] <- learnErrors(
            fi$forward, nbases = config$nbases,
            multithread = config$threads, verbose = 0)
        if (paired) {
            errors[[r]][["reverse"]] <- learnErrors(
                fi$reverse, nbases = config$nbases,
                multithread = config$threads, verbose = 0)
        }
        flog.info("Dereplicating run `%s` (%d samples)...", r, nrow(fi))
        derep_forward <- derepFastq(fi$forward)
        dada_stats[[r]][, "derep_forward" := getN(derep_forward)]
        if (paired) {
            derep_reverse <- derepFastq(fi$reverse, )
            dada_stats[[r]][, "derep_reverse" := getN(derep_reverse)]
        }
        flog.info("Inferring sequence variants for run `%s`...", r)
        dada_forward <- dada(derep_forward, err = errors[[r]]$forward,
                             multithread = config$threads, verbose = 0,
                             pool = config$pool)
        dada_stats[[r]][, "denoised_forward" := getN(dada_forward)]
        if (paired) {
            dada_reverse <- dada(derep_reverse, err = errors[[r]]$reverse,
                                 multithread = config$threads, verbose = 0,
                                 pool = config$pool)
            dada_stats[[r]][, "denoised_reverse" := getN(dada_reverse)]
            merged <- mergePairs(dada_forward, derep_forward,
                                 dada_reverse, derep_reverse,
                                 verbose = 0)
            dada_stats[[r]][, "merged" := getN(merged)]
            feature_table[[r]] <- makeSequenceTable(merged)
        } else {
            feature_table[[r]] <- makeSequenceTable(dada_forward)
        }
        rownames(feature_table[[r]]) <- fi$id
        flog.info("Finished run `%s`.", r)
    }
    feature_table <- mergeSequenceTables(tables = feature_table)
    dada_stats <- rbindlist(dada_stats)
    if ("passed" %in% names(object)) {
        dada_stats <- object$passed[dada_stats, on = "id"]
    }

    flog.info("Merged sequence tables, now removing chimeras...")
    feature_table_nochim <- removeBimeraDenovo(feature_table,
                                               multithread = config$threads)
    flog.info("Removed %d/%d sequence variants as chimeric (%.2f%% of reads)",
              ncol(feature_table) - ncol(feature_table_nochim),
              ncol(feature_table),
              100 - 100 * sum(feature_table_nochim) / sum(feature_table))
    dada_stats[, "non_chimera" := rowSums(feature_table_nochim)[id]]

    flog.info("Assigning taxonomy to %d sequences...",
              ncol(feature_table_nochim))
    tmp <- tempdir()
    if (grepl("tp(s*)://", config$taxa_db)) {
        taxa_db <- file.path(tmp, "taxa.fna.gz")
        flog.info("Downloading taxa db to %s...", taxa_db)
        download.file(config$taxa_db, taxa_db)
    } else {
        taxa_db <- config$taxa_db
    }
    if (grepl("tp(s*)://", config$species_db)) {
        species_db <- file.path(tmp, "species.fna.gz")
        flog.info("Downloading taxa db to %s...", species_db)
        download.file(config$species_db, species_db)
    } else {
        species_db <- config$taxa_db
    }
    taxa <- assignTaxonomy(feature_table_nochim, taxa_db,
                           minBoot = config$bootstrap_confidence * 100,
                           multithread = config$threads)
    taxa <- addSpecies(taxa, species_db)
    seqs <- rownames(taxa)
    taxa <- cbind(taxa, sequence = seqs)
    classified <- classified_taxa(feature_table, taxa)
    flog.info("Reads with taxonomic classification: %s",
              classified[, paste(rank, "=", 100 * round(reads, 3),
                                 "%", collapse = ",")])
    if (config$hash) {
        flog.info("Hashing %d sequence variants.", nrow(taxa))
        seqs <- rownames(taxa)
        hashes <- sapply(seqs, digest)
        names(seqs) <- hashes
        rownames(taxa) <- hashes
        colnames(feature_table_nochim) <- hashes[
            colnames(feature_table_nochim)]
    }

    artifact <- list(
        feature_table = feature_table_nochim,
        taxonomy = taxa,
        errors = errors,
        error_plots = lapply(errors, function(x)
            lapply(x, plotErrors, nominalQ = TRUE)),
        passed_reads = dada_stats,
        classified = classified,
        steps = c(object[["steps"]], "denoise")
    )
    return(artifact)
}


#' Converts the denoise artifact to a phyloseq object.
#'
#' @param object The artifact returned by \code{\link{denoise}}.
#' @param metadata Sample metadata to add to the phyloseq object. Must have
#'  rownames that match the sample names.
#' @return A phyloseq object with additional annotations.
#' @export
as_phyloseq <- function(object, metadata = NULL) {
    if (!all(c("feature_table", "taxonomy", "passed_reads", "steps") %in%
             names(object))) {
        stop("This is not a denoise artifact :(")
    }

    ps <- phyloseq(otu_table(object$feature_table, taxa_are_rows = FALSE),
                   tax_table(object$taxonomy))
    if (!is.null(metadata)) {
        sample_data(ps) <- metadata
    }
    attr(ps, "passed_reads") <- object$passed_reads
    attr(ps, "workflow_steps") <- object$steps
    return(ps)
}

#' Quantify classification rates for each taxonomic rank.
#'
#' @param ps A phyloseq object.
#' @return A data table with 3 columns (rank, asvs, reads) specifying the
#'  fraction of ASVs or reads that were classified successfully on that
#'  taxa rank.
#' @export
#' @importFrom phyloseq otu_table tax_table
with_classification <- function(ps) {
    return(classified_taxa(otu_table(ps, taxa_are_rows = FALSE),
                           tax_table(ps)))
}
