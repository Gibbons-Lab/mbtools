# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

#' Plots counts for several taxa across a co-variable
#'
#' @param ps A phyloseq data set.
#' @param variable The name of the co-variable.
#' @param tax_level The taxonomy level to use. Defaults to genus.
#' @param taxa A character vector denoting the taxa to be plotted. Defaults
#'  to plotting all taxa.
#' @param normalized Whether to normalize the counts using the DESeq2 size
#'  factors.
#' @param pc The pseudo count to add.
#' @param only_data Only get the raw data for the plot.
#' @return Nothing.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom ggplot2 ggplot geom_boxplot facet_wrap scale_y_log10 xlab
plot_counts <- function(ps, variable, tax_level = "genus", taxa = NULL,
                        normalized = TRUE, pc = 0.5, only_data = FALSE) {
    dts <- taxa_count(ps, lev = tax_level)
    valid_taxa <- taxa
    if (normalized) {
        dts <- normalize(dts)
    }

    if (is.null(valid_taxa)) {
        valid_taxa <- unique(dts[["taxa"]])
    }

    dts <- dts[taxa %in% valid_taxa]
    dts$variable <- variable
    dts$value <- sample_data(ps)[dts$sample, variable]
    if (only_data) return(dts)

    if (is.integer(dts$value) || is.factor(dts$value)) {
        pl <- ggplot(dts, aes(x = value, y = reads + pc, group = value)) +
              geom_boxplot() + facet_wrap(~ taxa) + scale_y_log10() +
              xlab(variable)
    } else {
        pl <- ggplot(dts, aes(x = value, y = reads + pc)) +
              geom_point(alpha = 0.5) + facet_wrap(~ taxa) + scale_y_log10() +
              xlab(variable)
    }

    return(pl)
}


shorten <- function(texts, n=40) {
    before <- sapply(texts, nchar)
    texts <- substr(texts, 1, n)
    after <- sapply(texts, nchar)
    texts[before > after] <- paste0(texts[before > after], "...")
    return(texts)
}


#' Plots relative distribution for taxa across samples.
#'
#' @param ps A phyloseq data set.
#' @param level The taxonomy level to use. Defaults to phylum.
#' @param sort Whether to sort taxa by abundance across all samples.
#' @param max_taxa Maximum number of different taxa to plot. If more than 12
#'  there is probably no color scale that can visualize them.
#' @param only_data Only get the raw data for the plot as a data table.
#' @return Nothing or a data.table containing the relative abundances.
#' @examples
#'  NULL
#'
#' @export
plot_taxa <- function(ps, level="Phylum", sort=TRUE,
                      max_taxa = 12, only_data = FALSE) {
    counts <- taxa_count(ps, lev=level)[, reads := as.double(reads)]
    counts[, reads := reads / sum(reads), by = "sample"]
    if (is.na(level)) {
        counts[, taxa := paste0(species, ": ", taxa)]
    }
    total_ord <- counts[, sum(reads, na.rm=TRUE), by = "taxa"][order(-V1), taxa]
    if (length(total_ord) > max_taxa) {
        total_ord <- total_ord[1:max_taxa]
        counts <- counts[taxa %in% total_ord]
    }
    sample_ord <- counts[taxa == total_ord[1]][order(-reads), sample]
    counts[, taxa := factor(taxa, levels=rev(total_ord))]
    counts[, sample := factor(sample, levels=sample_ord)]
    counts[, id := as.numeric(sample)]

    if (only_data) return(counts)

    pl <- ggplot(counts, aes(x=id, y=reads, fill=taxa)) +
        geom_bar(stat="identity", col=NA, width=1) +
        scale_x_continuous(expand = c(0, 1)) +
        scale_y_continuous(expand = c(0, 0.01)) +
        scale_fill_brewer(palette="Paired", direction = -1, label=shorten) +
        xlab("sample index") + ylab("% of reads") + labs(fill="") +
        theme_bw()

    return(pl)
}
