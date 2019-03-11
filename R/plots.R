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
              geom_boxplot(outlier.color = NA) +
              geom_jitter(width = 0.2, alpha = 0.5, size = 1, stroke = 0) +
              facet_wrap(~ taxa) + scale_y_log10() +
              xlab(variable)
    } else {
        pl <- ggplot(dts, aes(x = value, y = reads + pc)) +
              geom_point(alpha = 0.5, stroke = 0) +
              facet_wrap(~ taxa) + scale_y_log10() +
              xlab(variable)
    }

    return(pl)
}

#' Plots relative distribution for taxa across samples.
#'
#' @param ps A phyloseq data set.
#' @param level The taxonomy level to use. Defaults to phylum.
#' @param sort Whether to sort taxa by abundance across all samples.
#' @param show_samples Whether to show sample names.
#' @param max_taxa Maximum number of different taxa to plot. If more than 12
#'  there is probably no color scale that can visualize them.
#' @param only_data Only get the raw data for the plot as a data table.
#' @return Nothing or a data.table containing the relative abundances.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom scales percent
#' @importFrom stringr str_trunc
plot_taxa <- function(ps, level = "Phylum", show_samples = TRUE, sort = TRUE,
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
    x <- "id"
    if (show_samples) x <- "sample"

    pl <- ggplot(counts, aes_string(x=x, y="reads", fill="taxa")) +
        geom_bar(stat="identity", col=NA, width=1) +
        scale_y_continuous(expand = c(0, 0.01), labels = percent) +
        scale_fill_brewer(palette = "Paired", direction = -1,
                          label = function(x) str_trunc(x, 30)) +
        labs(x = ifelse(show_samples, "", "sample index"),
             y = "relative abundance", fill = "") +
        theme_bw()

    if (show_samples) {
        pl <- pl + theme(axis.text.x = element_text(angle = 90,
                                                    hjust = 1, vjust = 0.5))
    } else {
        pl <- pl + scale_x_continuous(limits = range(counts[[x]]),
                                      expand = c(0, 0))
    }

    return(pl)
}
