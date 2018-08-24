# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

# So check does not complain :()
utils::globalVariables(c("reference", "measured", "level"))

# mockrobiota base address
mb <- "https://raw.githubusercontent.com/caporaso-lab/mockrobiota/master/data/"
dl <- c("raw-data-url-forward-read", "raw-data-url-reverse-read",
        "raw-data-url-index-read")

L <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

download_reads <- function(url, folder, quiet) {
    if (is.na(url)) return(NA)
    p <- file.path(folder, basename(url))
    download.file(url, p, quiet = quiet)
    p
}

#' Downloads a complete data set from mockrobiota.
#'
#' @param name Name of the mockrobiota data set.
#' @param folder Where to save the data.
#' @param quiet Whether to show download progress.
#' @return A list with the following components
#' \describe{
#' \item{description}{A description of the data set.}
#' \item{forward}{Filepath to the forward reads.}
#' \item{reverse}{Filepath to the reverse reads.}
#' \item{index}{Filepath to the index file.}
#' \item{citation}{Reference for the data set.}
#' \item{fragment}{Which fragment(s) were sequenced.}
#' \item{equipment}{The used sequencing equipment.}
#' \item{samples}{A data frame mapping samples to barcodes.}
#' \item{ps_gg}{The reference phyloseq object using the green genes taxonomy.}
#' \item{ps_silva}{The reference phyloseq object using the silva taxonomy.}
#'}
#' @examples
#'  NULL
#'
#' @export
mockrobiota <- function(name, folder, quiet=!interactive()) {
    mock_info <- sprintf("%s/%s/dataset-metadata.tsv", mb, name)
    mock_samples <- sprintf("%s/%s/sample-metadata.tsv", mb, name)
    info <- read.table(mock_info, sep = "\t", header = TRUE)
    ivec <- as.character(info[, 2])
    names(ivec) <- as.character(info[, 1])

    dir.create(folder, showWarnings = FALSE)
    dl_list <- ivec[dl]
    dl_list <- dl_list[!is.na(dl_list)]

    downloaded <- vapply(dl_list, download_reads, "", folder = folder,
                         quiet = quiet)
    samples <- read.table(mock_samples, header = TRUE)
    rownames(samples) <- samples[, 1]

    gg <- sprintf("%s/%s/greengenes/13-8/97-otus/expected-taxonomy.tsv",
                  mb, name)
    gg <- read.table(gg, header = TRUE, sep = "\t")
    taxa <- as.character(gg[, 1])
    taxa <- do.call(rbind, strsplit(taxa, ";\\s*", perl = TRUE))
    colnames(taxa) <- L[1:ncol(taxa)]
    gg <- phyloseq(tax_table(as.matrix(taxa)),
                   otu_table(gg[, -1, drop = FALSE], taxa_are_rows = TRUE),
                   sample_data(samples))
    silva <- sprintf("%s/%s/silva/123/99-otus/expected-taxonomy.tsv", mb, name)
    silva <- read.table(silva, header = TRUE, sep = "\t")
    taxa <- as.character(silva[, 1])
    taxa <- do.call(rbind, strsplit(taxa, ";\\s*", perl = TRUE))
    colnames(taxa) <- L[1:ncol(taxa)]
    silva <- phyloseq(tax_table(as.matrix(taxa)),
                      otu_table(silva[, -1, drop = FALSE], taxa_are_rows = TRUE),
                      sample_data(samples))

    list(
        description = ivec["human-readable-description"],
        forward = downloaded[1],
        reverse = downloaded[2],
        index = ifelse(length(dl_list) == 3, downloaded[3], NA),
        citation = ivec["citation"],
        fragment = ivec["target-subfragment"],
        equipment = ivec["sequencing-instrument"],
        samples = samples,
        ps_gg = gg,
        ps_silva = silva
    )
}

taxa_str <- function(taxa, level) {
    index <- which(colnames(taxa) == level)
    bad <- apply(taxa[, 1:index, drop = FALSE], 1, function(x)
        any(is.na(x) | nchar(x) == 0))
    strs <- apply(taxa[, 1:index, drop = FALSE], 1, paste, collapse = ";")
    strs[bad] <- NA
    return(strs)
}

#' Checks whether taxa from one taxonomy table are contained in another table.
#'
#' This function takes two taxonomy tables, than looks for each sample of the first
#' table in the second one.
#'
#' @param taxa1 First taxonomy table.
#' @param taxa2 Second taxonomy table.
#' @param level At which level to compare. Must be column name in both tables.
#' @return A data frame with two columns. Has as many rows as unique values in level.
#' \describe{
#' \item{level}{The unique values found for the specified taxa level.}
#' \item{found}{A boolean indicating whether the taxa were found in the second table.}
#'}
#' @examples
#'  NULL
#'
#' @export
find_taxa <- function(taxa1, taxa2, level="Species") {
    if (!(level %in% colnames(taxa1) && level %in% colnames(taxa2)))
        stop("level must be a valid column name in both tables!")

    index <- which(colnames(taxa1) == level)
    snames <- unique(taxa_str(taxa1, level))
    snames <- snames[!is.na(snames)]
    snames_ref <- taxa_str(taxa2, level)
    found <- snames %in% snames_ref
    names(found) <- snames

    return(found)
}

#' Calculates what percentage of taxa was found in a reference set.
#'
#' Note that both arguments must of class "taxonomyTable" from the
#' phyloseq package which can be obtained with the \code{tax_table}
#' method.
#'
#' @param tax_table Measured taxonomy table.
#' @param ref Reference taxonomy table.
#' @return A data frame denoting performance metrics for taxa identification.
#' @examples
#'  NULL
#'
#' @export
taxa_metrics <- function(tax_table, ref) {
    if (any(colnames(tax_table) != colnames(ref)))
        stop("Both taxonomy tables need to have the same column names!")

    metrics <- data.frame()
    for (cn in colnames(tax_table)) {
        prec <- find_taxa(tax_table, ref, level = cn)
        precision <- sum(prec) / length(prec)
        rec <- find_taxa(ref, tax_table, level = cn)
        recall <- sum(rec) / length(rec)
        new <- data.frame(level = cn, precision = precision, recall = recall,
                          F1 = 2 * precision * recall / (precision + recall),
                          n_exp = length(prec), n_ref = length(rec))
        metrics <- rbind(metrics, new)
    }

    return(metrics)
}

#' Compares taxa quantification from a measurement to a reference ground truth.
#'
#' Note that both arguments must be phyloseq objects from the
#' phyloseq package.
#'
#' @param ps A phyloseq object describing the measurements
#' @param ref A phyloseq object describing the reference set. Can be obtained
#'  from \code{\link{mockrobiota}}.
#' @param normalize Whether to normalize taxa counts first.
#' @return A data frame with the following columns:
#'  \describe{
#'  \item{level}{The taxonomy level for the entry.}
#'  \item{name}{The taxonomy.}
#'  \item{sample}{The sample name for the quantification.}
#'  \item{measured}{The measured quantification.}
#'  \item{reference}{The reference quantification.}
#' }
#' @examples
#'  NULL
#'
#' @export
taxa_quants <- function(ps, ref, normalize = FALSE) {
    taxa1 <- tax_table(ps)
    taxa2 <- tax_table(ref)
    otu1 <- otu_table(ps)
    otu2 <- otu_table(ref)
    if (taxa_are_rows(ps)) otu1 <- t(otu1)
    if (taxa_are_rows(ref)) otu2 <- t(otu2)

    n <- ncol(taxa1)
    ns <- nrow(otu1)


    if (any(colnames(taxa1) != colnames(taxa2)))
        stop("Both taxonomy tables need to have the same column names!")

    if (any(rownames(otu1) != rownames(otu2)))
        stop("Both phyloseq objects need to have the same samples!")

    x <- data.frame()
    for (cn in colnames(taxa1)) {
        tax_m <- taxa_str(taxa1, cn)
        tax_r <- taxa_str(taxa2, cn)
        found <- find_taxa(taxa1, taxa2, level = cn)
        found <- names(found)[found]
        for (i in 1:ns) {
            measured <- as.numeric(otu1[i, ])
            reference <- as.numeric(otu2[i, ])
            if (normalize) {
                measured <- measured / sum(measured)
                reference <- reference / sum(reference)
            }
            measured <- tapply(measured, tax_m, sum, na.rm = TRUE)
            reference <- tapply(reference, tax_r, sum, na.rm = TRUE)
            if (length(found) == 0) {
                new <- NULL
            } else {
                new <- data.frame(level = cn, name = found,
                                  sample = sample_names(ps)[i],
                                  measured = measured[found],
                                  reference = reference[found])
            }
            x <- rbind(x, new)
        }
    }
    rownames(x) <- NULL

    return(x)
}

#' Creates a plot of measured taxa quantifications vs. reference quantification.
#'
#' Note that both arguments must be phyloseq objects from the
#' phyloseq package.
#'
#' @param ps A phyloseq object describing the measurements
#' @param ref A phyloseq object describing the reference set. Can be obtained
#'  from \code{\link{mockrobiota}}.
#' @return A ggplot2 plot.
#' @examples
#'  NULL
#'
#' @export
mock_plot <- function(ps, ref) {
    quants <- taxa_quants(ps, ref, normalize = TRUE)
    ggplot(quants, aes(x = reference, y = measured, col = level)) +
        geom_abline(alpha = 0.5) +
        geom_smooth(aes(group = 1), method = "lm", lty = "dashed") +
        geom_point(aes(col = level)) + facet_wrap(~ samples) + theme_bw()
}
