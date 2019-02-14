# Copyright 2019 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.


qualities <- function(srqa) {
    cycles <- srqa[["perCycle"]]$quality %>% as.data.table
    names(cycles) <- c("cycle", "quality_symbol", "quality",
                       "count", "file")
    avg_score <- cycles[, sum(quality * count) / sum(count)]
    flog.info("Average per base error is %.3f%% (mean score = %.2f).",
              10 ^ (-avg_score / 10 + 2), avg_score)
    return(cycles)
}

entropy <- function(counts) {
    p <- (counts + 1) / sum(counts + 1)
    return(-sum(p * log2(p)))
}

basecalls <- function(srqa) {
    cycles <- srqa[["perCycle"]]$baseCall %>% as.data.table
    names(cycles) <- c("cycle", "base", "count", "file")
    med_entropy <- cycles[, entropy(count), by = c("cycle", "file")] %>% mean
    flog.info("Mean per cycle entropy is %.3f%% (in [0, 2]).",
              med_entropy)
    return(cycles)
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
    cycles[, direction := "forward"]
    bases[, direction := "forward"]
    cycles <- files[, .(file = basename(forward), id)][cycles, on = "file"]
    bases <- files[, .(file = basename(forward), id)][bases, on = "file"]
    if ("reverse" %in% names(files)) {
        flog.info("Running quality assay for forward reads from %d files.",
              files[, .N])
        rcycles <- qualities(files$reverse)
        rbases <- basecalls(srqa)
        rcycles[, direction := "reverse"]
        rbases[, direction := "reverse"]
        rcycles <- files[, .(file = basename(reverse), id)][
            rcycles, on = "file"]
        rbases <- files[, .(file = basename(reverse), id)][
            rbases, on = "file"]
        cycles <- rbind(cycles, rcycles)
        bases <- rbind(bases, rbases)
    }
    return(list(qualities = cycles, bases = bases))
}

#' Plots the quality profile for an entire experiment.
#'
#' @param cycles A quality profile as returned by
#'  \code{\link{quality_profile}}.
#' @return A ggplot2 plot mapping the read positions to mean quality for each
#'  sample. Will be facetted into fowrward and reverse if applicable.
#'
#' @export
plot_qualities <- function(qp) {
    cycles <- qp$qualities
    mean_scores <- cycles[, .(quality = sum(count * quality) / sum(count),
                              direction),
                          by = c("cycle", "file")]
    pl <- ggplot(mean_scores, aes(x = cycle, y = quality)) +
            geom_hline(yintercept = 10, col = "red") +
            geom_point(alpha = 0.1, stroke = 0, size = 2) +
            geom_smooth() + theme_bw() +
            labs(x = "read position / cycle [bps]", y = "quality score")
    n <- cycles[, uniqueN(direction)]
    if (n > 1) {
        pl <- pl + facet_wrap(~ direction)
    }
    return(pl)
}


#' Plots the distribution of sequence lengths with acceptable quality.
#'
#' @param cycles A quality profile as returned by
#'  \code{\link{quality_profile}}.
#' @param min_score Smallest acceptable score. Defaults to 10 (10% per base
#'  error).
#' @return A ggplot2 plot mapping the read positions to mean quality for each
#'  sample. Will be facetted into fowrward and reverse if applicable.
#'
#' @export
plot_lengths <- function(qp, min_score = 10) {
    nice <- qp$qualities[quality > min_score,
                        .(count = sum(count), direction, id),
                        by = c("cycle", "file")]
    pl <- ggplot(nice, aes(x = cycle, y = id, fill = count)) +
            geom_raster() +
            scale_fill_gradient(low = "white", high = "blue",
                                limits = c(0, NA)) +
            scale_x_continuous(limits = range(nice$cycle), expand = c(0, 0)) +
            theme_bw() +
            labs(x = "read position / cycle [bps]", y = "")
    n <- cycles[, uniqueN(direction)]
    if (n > 1) {
        pl <- pl + facet_wrap(~ direction)
    }
    return(pl)
}

#' Plots the distribution of base entropy for each cycle.
#'
#' @param cycles A quality profile as returned by
#'  \code{\link{quality_profile}}.
#' @return A ggplot2 plot mapping the read positions to mean quality for each
#'  sample. Will be facetted into forward and reverse if applicable.
#'
#' @export
plot_entropy <- function(qp) {
    ent <- qp$bases[, .(entropy = entropy(count)), by = c("cycle", "file")]
    pl <- ggplot(ent, aes(x = cycle, y = entropy)) +
            geom_hline(yintercept = 2, col = "red") +
            geom_point(alpha = 0.05, stroke = 0) +
            stat_smooth() + theme_bw() +
            labs(x = "read position / cycle [bps]", y = "entropy [bits]")
    n <- qp$bases[, uniqueN(direction)]
    if (n > 1) {
        pl <- pl + facet_wrap(~ direction)
    }
    return(pl)
}
