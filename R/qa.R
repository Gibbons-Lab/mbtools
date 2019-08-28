# Copyright 2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.

count_med <- function(vals, counts) {
    o <- order(vals)
    cum <- cumsum(as.numeric(counts[o]))
    idx <- which(cum > 0.5 * sum(as.numeric(counts)))[1]
    return(vals[o][idx])
}

qualities <- function(srqa) {
    cycles <- srqa[["perCycle"]]$quality %>% as.data.table
    names(cycles) <- c("cycle", "quality_symbol", "quality",
                       "count", "file")
    med_score <- cycles[, count_med(quality, count)]
    flog.info("Median per base error is %.3f%% (median score = %.2f).",
              10 ^ (-med_score / 10 + 2), med_score)
    return(cycles)
}

entropy <- function(counts) {
    p <- (counts + 1) / sum(counts + 1)
    return(-sum(p * log2(p)))
}

basecalls <- function(srqa) {
    cycles <- srqa[["perCycle"]]$baseCall %>% as.data.table
    names(cycles) <- c("cycle", "base", "count", "file")
    med_entropy <- cycles[, entropy(count), by = c("cycle", "file")][, V1]
    flog.info("Mean per cycle entropy is %.3f (in [0, 2]).",
              mean(med_entropy, na.rm = TRUE))
    return(cycles)
}

#' @importFrom stats setNames
library_size <- function(srqa) {
    sizes <- srqa[["readCounts"]]
    sizes <- data.table(file = rownames(sizes), count = sizes[, 1])
    flog.info("On average we have %.2f +- %.2f reads per file.",
              sizes[, mean(count)], sizes[, sd(count)])
    return(sizes)
}

#' Get the quality profiles and base calls for each cycle across all
#' input files.
#'
#' @param files A data table containing the files as returned by
#'  \code{\link{find_read_files}}.
#' @param n How many reads to sample from each file. If a file has less
#'  all of them will be used.
#' @return A list with two data tables giving the sampled metrics per cycle /
#'  base pair:
#'  \describe{
#'    \item{qualities}{The quality scores for each cycle and sample.}
#'    \item{bases}{The base calls for each cycle and sample.}
#'  }
#' @export
#' @importFrom ShortRead qa
quality_profile <- function(files, n = 1e4) {
    flog.info("Running quality assay for forward reads from %d files.",
              files[, .N])
    srqa <- qa(files$forward, sample = TRUE, n = n)
    cycles <- qualities(srqa)
    bases <- basecalls(srqa)
    sizes <- library_size(srqa)
    cycles[, direction := "forward"]
    bases[, direction := "forward"]
    sizes[, direction := "forward"]
    cycles <- files[, .(file = basename(forward), id)][cycles, on = "file"]
    bases <- files[, .(file = basename(forward), id)][bases, on = "file"]
    sizes <- files[, .(file = basename(forward), id)][sizes, on = "file"]
    if ("reverse" %in% names(files)) {
        flog.info("Running quality assay for reverse reads from %d files.",
              files[, .N])
        srqa <- qa(files$reverse, sample = TRUE, n = n)
        rcycles <- qualities(srqa)
        rbases <- basecalls(srqa)
        rsizes <- library_size(srqa)
        rcycles[, direction := "reverse"]
        rbases[, direction := "reverse"]
        rsizes[, direction := "reverse"]
        rcycles <- files[, .(file = basename(reverse), id)][
            rcycles, on = "file"]
        rbases <- files[, .(file = basename(reverse), id)][
            rbases, on = "file"]
        rsizes <- files[, .(file = basename(reverse), id)][
            rsizes, on = "file"]
        cycles <- rbind(cycles, rcycles)
        bases <- rbind(bases, rbases)
        sizes <- rbind(sizes, rsizes)
    }
    return(list(qualities = cycles, bases = bases, sizes = sizes))
}

#' Plots the quality profile for an entire experiment.
#'
#' @param qp A quality profile as returned by
#'  \code{\link{quality_profile}}.
#' @param min_score Smallest acceptable score. Defaults to 10 (10% per base
#'  error). A reference line will be plotted but no data will be discarded.
#' @return A ggplot2 plot mapping the read positions to mean quality for each
#'  sample. Will be facetted into fowrward and reverse if applicable.
#'
#' @export
plot_qualities <- function(qp, min_score = 10) {
    cycles <- qp$qualities
    mean_scores <- cycles[, .(quality = sum(count * quality) / sum(count),
                              direction[1]),
                          by = c("cycle", "file")]
    pl <- ggplot(mean_scores, aes(x = cycle, y = quality)) +
            geom_hline(yintercept = min_score, col = "blue") +
            geom_bin2d(bins = 50) +
            geom_smooth(col = "tomato", fill = "tomato", method = "gam",
                        formula = y ~ s(x, bs = "cs")) +
            scale_fill_viridis_c(limits = c(0, NA)) +
            labs(x = "read position / cycle [bps]",
                 y = "quality score [mean per sample]")
    n <- cycles[, uniqueN(direction)]
    if (n > 1) {
        pl <- pl + facet_wrap(~ direction)
    }
    return(pl)
}


#' Plots the distribution of sequence lengths with acceptable quality.
#'
#' @param qp A quality profile as returned by
#'  \code{\link{quality_profile}}.
#' @param min_score Smallest acceptable score. Defaults to 10 (10% per base
#'  error).
#' @return A ggplot2 plot mapping the read positions to mean quality for each
#'  sample. Will be facetted into fowrward and reverse if applicable.
#'
#' @export
plot_lengths <- function(qp, min_score = 10) {
    nice <- qp$qualities[quality > min_score,
                        .(count = sum(count), direction[1], id[1]),
                        by = c("cycle", "file")]
    pl <- ggplot(nice, aes(x = cycle, y = id, fill = count)) +
            geom_raster() +
            scale_fill_gradient(low = "white", high = "blue",
                                limits = c(0, NA)) +
            scale_x_continuous(limits = range(nice$cycle), expand = c(0, 0)) +
            labs(x = "read position / cycle [bps]", y = "")
    n <- qp$qualities[, uniqueN(direction)]
    if (n > 1) {
        pl <- pl + facet_wrap(~ direction)
    }
    return(pl)
}

#' Plots the distribution of base entropy for each cycle.
#'
#' @param qp A quality profile as returned by
#'  \code{\link{quality_profile}}.
#' @return A ggplot2 plot mapping the read positions to mean quality for each
#'  sample. Will be facetted into forward and reverse if applicable.
#'
#' @export
plot_entropy <- function(qp) {
    ent <- qp$bases[, .(entropy = entropy(count), direction[1]),
                    by = c("cycle", "file")]
    pl <- ggplot(ent, aes(x = cycle, y = entropy)) +
            geom_hline(yintercept = 2, col = "blue") +
            geom_bin2d(bins = 50) +
            scale_fill_viridis_c(limits = c(0, NA)) +
            stat_smooth(col = "tomato", fill = "tomato", method = "gam",
                        formula = y ~ s(x, bs = "cs")) +
            labs(x = "read position / cycle [bps]", y = "entropy [bits]")
    n <- qp$bases[, uniqueN(direction)]
    if (n > 1) {
        pl <- pl + facet_wrap(~ direction)
    }
    return(pl)
}
