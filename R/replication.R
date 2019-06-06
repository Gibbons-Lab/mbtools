# Peak-to-through ratio analysis

#' Build a configuration for the rate workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the long read alignment
#'  workflow.
#' @export
#' @examples
#'  config <- config_rep(min_coverage = 0.8)
config_rep <- config_builder(list(
    remove_extremes = 0.1,
    max_median_fold = 8,
    min_coverage = 2,
    min_covered = 0.6,
    threads = getOption("mc.cores", 1)
))

#' @importFrom stats lm anova coef cor median
irep <- function(profile, conf) {
    profile <- copy(profile)
    w <- profile$bin_width
    reads <- profile$reads[[1]]
    profile[, "reads" := NULL]
    co <- reads * profile$read_length / w
    slided <- frollmean(co, 50, align = "center")
    slided[slided <= conf$min_coverage] <- NA
    m <- median(slided, na.rm = TRUE)
    slided[slided <= m / conf$max_median_fold |
           slided > m * conf$max_median_fold] <- NA
    sufficient <- ((sum(!is.na(slided)) / length(slided)) > conf$min_covered &
                   sum(!is.na(slided)) > 60)
    if (!sufficient) {
        return(list(
            rate = NULL,
            profile = NULL
        ))
    }
    so <- sort(slided[!is.na(slided)])
    k <- ceiling(conf$remove_extremes * length(slided))
    so <- so[k:(length(slided) - k)]
    pos <- seq_along(so)
    fit <- lm(log(so) ~ pos)
    coefs <- coef(fit)
    rate <- exp(coefs[2] * max(pos))
    profile[, "coverage" := list(list(coverage = slided))]
    res <- list(
        rates = data.table(intercept = coefs[1],
                           slope = coefs[2],
                           rate = rate,
                           rsquared = summary(fit)$r.squared,
                           pval = anova(fit)[["Pr(>F)"]][1]),
        profile = profile
    )
    return(res)
}

#' Estimate replication rates by the iRep method.
#'
#' Implements the method from Brown et. al.
#' (https://dx.doi.org/10.1038%2Fnbt.3704).
#'
#' @param object An artifact containing a coverage map as returned by
#' \code{link{slimm}}.
#' @param ... Additional configuration parameters or a config object as
#' returned by \code{\link{config_rep}}.
#' @return An artifact containing the peak-to-through ratios as well as the
#'  smoothed coverage profiles.
#' @export
replication_rates <- function(object, ...) {
    if (!"coverage" %in% names(object)) {
        stop("Need `coverage` in artifacts!")
    }
    co <- object$coverage
    config <- config_parser(list(...), config_rep)
    apfun <- parse_threads(config$threads)

    mc <- co[, mean(reads[[1]] * read_length / bin_width),
             by = c("id", "genbank")][, V1]
    flog.info(paste("Estimating replication rates for %d sample-genome",
                    "combinations. Median coverage is %.3g [%.3g, %.3g]."),
                    nrow(co), median(mc), min(mc), max(mc))
    rates <- apfun(1:nrow(co), function(i) {
        row <- co[i]
        res <- irep(row, config)
        if (!is.null(res$rate)) {
            res$profile[, "genbank" := row$genbank]
            res$profile[, "id" := row$id]
            res$rate[, "genbank" := row$genbank]
            res$rate[, "id" := row$id]
            taxmap <- unique(res$profile[, .(strain, species, genus, family,
                                             order, class, phylum, kingdom,
                                             id, genbank)])
            res$rate <- taxmap[res$rate, on = c("genbank", "id")]
        }
        return(res)
    })
    profiles <- lapply(rates, "[[", "profile") %>% rbindlist()
    rates <- lapply(rates, "[[", "rate") %>% rbindlist()
    flog.info(paste("Finished. %d genome-sample combinations had sufficient",
                    "coverage for obtaining replication rates."),
              nrow(rates))

    artifact <- list(
        rate = rates,
        coverage = profiles[genbank %in% rates$genbank],
        steps = c(object[["steps"]], "replication_rates")
    )
}
