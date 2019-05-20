# Peak-to-through ratio analysis

#' Build a configuration for the PTR workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the long read alignment
#'  workflow.
#' @export
#' @examples
#'  config <- config_ptr(span=0.5)
config_ptr <- config_builder(list(
    max_median_fold = 8,
    min_coverage = 2,
    min_covered = 0.6,
    threads = FALSE
))

#' @importFrom mgcv gam s
ptr <- function(profile, conf, rlen) {
    profile <- copy(profile)
    w <- profile[, start[2] - start[1] + 1]
    profile <- profile[, coverage := reads * rlen / w]
    profile[coverage <= conf$min_coverage, coverage := NA]
    profile[abs(log(reads + 1) - log(median(reads + 1, na.rm = TRUE))) >
            log(conf$max_median_fold), coverage := NA]
    if (profile[, (sum(!is.na(coverage)) / .N) < conf$min_covered]) {
        return(list(
            ptr = NULL,
            profile = profile[, "smooth" := NA]
        ))
    }
    coverage <- profile[!is.na(coverage), coverage]
    pos <- profile[!is.na(coverage), start]
    data <- data.table(
        coverage = c(coverage[1:(length(coverage) - 1)],
                     coverage, coverage[2:length(coverage)]),
        start = c(-rev(pos[2:length(pos)]), pos, max(pos) + pos[2:length(pos)])
    )
    fit <- gam(coverage ~ s(start, bs = "tp", k = 18), data = data)
    profile$smooth <- predict(fit, data.frame(start = profile$start))
    m <- profile[which.max(smooth)]
    ptr <- profile[, max(smooth) / min(smooth)]

    res <- list(
        ptr = data.table(start = m$start, end = m$end, ptr = ptr,
                         rsquared = summary(fit)$r.sq,
                         pval = summary(fit)$s.pv),
        profile = profile
    )
    return(res)
}

#' Estimate replication rates by the peak-to-through ratio.
#'
#' Implements the method from Korem et. al. (DOI: 10.1126/science.aac4812).
#'
#' @param object An artifact containing a coverage map as returned by
#' \code{link{slimm}}.
#' @param ... Additional configuration parameters or a config object as
#' returned by \code{\link{config_ptr}}.
#' @return An artifact containing the peak-to-through ratios as well as the
#'  smoothed coverage profiles.
#' @export
peak_to_through <- function(object, ...) {
    if (!"coverage" %in% names(object)) {
        stop("Need `coverage` in artifacts!")
    }
    co <- object$coverage
    config <- config_parser(list(...), config_ptr)
    apfun <- parse_threads(config$threads)

    flog.info(paste("Estimating read lengths from a sample of",
                    "100 reads per alignment."))
    alns <- get_alignments(object)
    rlens <- sapply(alns$alignment, read_length)
    names(rlens) <- alns$id
    flog.info("Estimated median read length is %d, range is [%d, %d].",
              median(rlens, na.rm = TRUE), min(rlens, na.rm = TRUE),
              max(rlens, na.rm = TRUE))

    genbank_id <- unique(co[, list(genbank, id)])
    genbank_id <- lapply(1:nrow(genbank_id),
                         function(i) as.character(genbank_id[i]))
    flog.info(paste("Calculating smoothed coverage and PTR for %d genomes",
                    "and %d samples."), co[, uniqueN(genbank)],
                    co[, uniqueN(id)])
    ptrs <- apfun(genbank_id, function(row) {
        profile <- co[genbank == row[1] & id == row[2]]
        res <- ptr(profile, config, rlens[row[2]])
        res$profile[, "genbank" := row[1]]
        res$profile[, "id" := row[2]]
        if (!is.null(res$ptr)) {
            res$ptr[, "genbank" := row[1]]
            res$ptr[, "id" := row[2]]
            taxmap <- unique(res$profile[, .(strain, species, genus, family,
                                             order, class, phylum, kingdom,
                                             id, genbank)])
            res$ptr <- taxmap[res$ptr, on = c("genbank", "id")]
        }
        return(res)
    })
    profiles <- lapply(ptrs, "[[", "profile") %>% rbindlist()
    ptrs <- lapply(ptrs, "[[", "ptr") %>% rbindlist()
    flog.info(paste("Finished. %d genome-sample combinations had sufficient",
                    "coverage for obtaining PTRs."),
              ptrs[, sum(!is.na(ptr))])

    artifact <- list(
        ptr = ptrs,
        coverage = profiles[genbank %in% ptrs$genbank],
        steps = c(object[["steps"]], "peak_to_through")
    )
}
