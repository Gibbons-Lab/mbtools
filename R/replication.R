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
    linear_quantile = 0.8,
    min_coverage = 2,
    min_covered = 0.6,
    threads = TRUE
))

#' @importFrom stats lm anova
irep <- function(profile, conf) {
    profile <- copy(profile)
    w <- profile$bin_width
    reads <- profile$reads[[1]]
    profile[, "reads" := NULL]
    co <- reads * profile$read_length / w
    slided <- frollmean(co, 50, align = "center")
    slided[slided <= conf$min_coverage] <- NA
    sufficient <- ((sum(!is.na(slided)) / length(slided)) > conf$min_covered &
                   sum(!is.na(slided)) > 50)
    if (!sufficient) {
        return(list(
            rate = NULL,
            profile = NULL
        ))
    }
    m <- median(slided, na.rm = TRUE)
    q <- conf$linear_quantile
    d <- min(quantile(slided[slided > m] - m, q, na.rm = TRUE),
             quantile(m - slided[slided < m], q, na.rm = TRUE))
    in_range <- m + c(-1, 1) * d
    slided[!between(slided, in_range[1], in_range[2])] <- NA
    slided <- sort(slided[!is.na(slided)])
    pos <- seq_along(slided)
    fit <- lm(slided ~ pos)
    coefs <- coef(fit)
    extremes <- coefs[1] + coefs[2] * range(pos)
    rate <- extremes[2] / extremes[1]
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
    profiles <- lapply(rates, function(l) {
        pro <- l$profile
        if (is.null(pro)) return(NULL)
        co <- data.table(rank = 1:length(pro$coverage[[1]]),
                         coverage = pro$coverage[[1]])
        pro[, coverage := NULL]
        cbind(pro, co)
    }) %>% rbindlist()
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
