# Copyright 2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' @export
config_dada2 <- function() {
    return(list(
        threads = TRUE,
        data_dir = "data",
        trimLeft = 10,
        truncLen = 0,
        maxEE = 2,
        nbases = 2.5e8,
        pool = "pseudo",
        bootstrap_confidence = 0.5,
        taxa_db = paste0("https://zenodo.org/record/1172783/files/",
                         "silva_nr_v132_train_set.fa.gz?download=1"),
        species_db = paste0("https://zenodo.org/record/1172783/files/",
                            "silva_species_assignment_v132.fa.gz?download=1"),
        hash = TRUE
    ))
}

#' @importFrom dada2 getUniques
getN <- function(x) sum(getUniques(x))


#' @export
#'
#' @importFrom dada2 filterAndTrim learnErrors dada removeBimeraDenovo
#'  mergeSequenceTables assignTaxonomy addSpecies
#' @importFrom digest digest
dada2 <- function(files, config) {
    files <- copy(files)
    if (!"run" %in% names(files)) {
        files[, "run" := "all"]
    }
    paired <- "reverse" %in% names(files)
    flog.info("Running DADA2 workflow for %d samples in %d runs.",
              nrow(files), files[, uniqueN(run)])

    flog.info("Preprocessing reads.")
    passed_files <- copy(files)
    passed_files$forward <- file.path(config$data_dir, "preprocessed",
                                      basename(files$forward))
    if (paired) {
        passed_files$reverse <- file.path(config$data_dir, "preprocessed",
                                          basename(files$reverse))
        passed_stats <- filterAndTrim(
            fwd = files$forward, filt = passed_files$forward,
            rev = files$reverse, filt.rev = passed_files$reverse,
            trimLeft = config$trimLeft, truncLen = config$truncLen,
            maxEE = config$maxEE, multithread = config$threads
        )
    } else {
        passed_stats <- filterAndTrim(
            fwd = files$forward, filt = passed_files$forward,
            trimLeft = config$trimLeft, truncLen = config$truncLen,
            maxEE = config$maxEE, multithread = config$threads
        )
    }
    passed_stats <- as.data.table(passed_stats)
    names(passed_stats) <- c("raw", "preprocessed")
    flog.info("%.2f%% of reads passed preprocessing.",
              passed_stats[, 100 * mean(preprocessed / raw)])

    flog.info("Running DADA2 on %d runs from a sample of %d bases.",
              files[, uniqueN(run)], config$nbases)
    errors <- list()
    dada_stats <- list()
    feature_table <- list()
    for (r in files[, unique(run)]) {
        fi <- passed_files[run == r]
        dada_stats[[r]] <- data.table(id = passed_files$id)
        flog.info("Learning errors for run `%s` (%d samples).", r, nrow(fi))
        errors[[r]] <- list()
        errors[[r]][["forward"]] <- learnErrors(
            passed_files$forward, nbases = config$nbases,
            multithread = config$threads, verbose = 0)
        if (paired) {
            errors[[r]][["reverse"]] <- learnErrors(
                passed_files$forward, nbases = config$nbases,
                multithread = config$threads, verbose = 0)
        }
        flog.info("Dereplicating run `%s` (%d samples).", r, nrow(fi))
        derep_forward <- derepFastq(passed_files$forward)
        names(derep_forward) <- passed_files$id
        if (paired) {
            derep_reverse <- derepFastq(passed_files$reverse)
            names(derep_reverse) <- passed_files$id
        }
        flog.info("Inferring sequence variants for run `%s`.", r)
        dada_forward <- dada(derep_forward, err = error[[r]]$forward,
                             multithread = config$processes, verbose = 0)
        dada_stats[[r]][, "denoised_forward" := sapply(dada_forward, getN)]
        if (paired) {
            dada_reverse <- dada(derep_reverse, err = error[[r]]$forward,
                                 multithread = config$processes, verbose = 0)
            dada_stats[[r]][, "denoised_reverse" := sapply(dada_reverse, getN)]
            merged <- mergePairs(dada_forward, derep_forward,
                                 dada_reverse, derep_reverse,
                                 verbose = 0)
            dada_stats[[r]][, "merged" := sapply(merged, getN)]
            feature_table[[r]] <- makeSequenceTable(merged)
        } else {
            feature_table[[r]] <- makeSequenceTable(dada_forward)
        }
        flog.info("Finished run `%s`.", r)
    }
    feature_table <- do.call(mergeSequenceTables, feature_table)
    dada_stats <- rbindlist(dada_stats)
    dada_stats <- dada_stats[passed_stats, on = "id"]
    flog.info("Removing chimeras.")
    feature_table_nochim <- removeBimeraDenovo(feature_table)
    flog.info("Removed %d/%d sequence variants as chimeric (%.2f%% of reads)",
              ncol(feature_table_nochim), ncol(feature_table),
              sum(feature_table_nochim) / sum(feature_table))
    dada_stats[, "non_chimera" := rowSums(feature_table_nochim)[id]]
    flog.info("Assigning taxonomy to %d sequences...",
              ncol(feature_table_nochim))
    taxa <- assignTaxonomy(feature_table_nochim, config$taxa_db,
                           minBoot = config$minBoot,
                           multithread = config$threads)
    taxa <- addSpecies(taxa, config$species_db)
    seqs <- rownames(taxa)
    taxa <- cbind(taxa, sequence = seqs)
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
        error_plots = lapply(errors, function(x) lapply(x, plot_errors)),
        passed_reads = dada_stats,
        filtered_files = passed_files
    )
    return(artifact)
}
