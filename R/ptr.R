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
    max_median_dev = 4,
    min_coverage = 2,
    min_covered = 0.6,
    threads = FALSE
))

#' @importFrom mgcv gam s
ptr <- function(profile, conf, rlen) {
    profile <- copy(profile)
    w <- profile$bin_width
    reads <- profile$reads[[1]]
    profile[, "reads" := NULL]
    co <- reads * rlen / w
    co[co <= conf$min_coverage] <- NA
    iqr_range <- median(co, na.rm = TRUE) +
                 c(-1, 1) * mad(co, na.rm = TRUE) * conf$max_median_dev
    co[!between(co, iqr_range[1], iqr_range[2])] <- NA
    sufficient <- ((sum(!is.na(co)) / length(co)) > conf$min_covered &
                   sum(!is.na(co)) > 50)
    if (!sufficient) {
        return(list(
            ptr = NULL,
            profile = NULL
        ))
    }
    pos <- seq_along(co)[!is.na(co)]
    coverage <- co[!is.na(co)]
    data <- data.table(
        coverage = c(coverage[1:(length(coverage) - 1)], coverage,
                     coverage[2:length(coverage)]),
        start = c(-rev(pos[2:length(pos)]), pos, max(pos) + pos[2:length(pos)])
    )
    fit <- gam(coverage ~ s(start, bs = "gp"), data = data)
    smooth <- predict(fit, data.frame(start = seq_along(co)))
    m <- which.max(smooth)
    ptr <- profile[, max(smooth) / min(smooth)]
    profile[, "smooth" := list(list(smooth = smooth))]
    profile[, "coverage" := list(list(coverage = co))]
    res <- list(
        ptr = data.table(start = (m - 1) * w + 1, end = m * w, ptr = ptr,
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

    flog.info(paste("Calculating smoothed coverage and PTR for %d genomes",
                    "and %d samples."), co[, uniqueN(genbank)],
                    co[, uniqueN(id)])
    ptrs <- apfun(1:nrow(co), function(i) {
        row <- co[i]
        res <- ptr(row, config, rlens[row$id])
        if (!is.null(res$ptr)) {
            res$profile[, "genbank" := row$genbank]
            res$profile[, "id" := row$id]
            res$ptr[, "genbank" := row$genbank]
            res$ptr[, "id" := row$id]
            taxmap <- unique(res$profile[, .(strain, species, genus, family,
                                             order, class, phylum, kingdom,
                                             id, genbank)])
            res$ptr <- taxmap[res$ptr, on = c("genbank", "id")]
        }
        return(res)
    })
    profiles <- lapply(ptrs, function(l) {
        pro <- l$profile
        if (is.null(pro)) return(NULL)
        co <- data.table(start = (seq_along(pro$coverage[[1]]) - 1) *
                                 pro$bin_width + 1,
                         coverage = pro$coverage[[1]],
                         smooth = pro$smooth[[1]])
        pro[, smooth := NULL]
        pro[, coverage := NULL]
        cbind(pro, co)
    }) %>% rbindlist()
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
