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
    read_length = 250,
    max_median_fold = 8,
    min_coverage = 3,
    min_covered = 0.6,
    threads = FALSE
))

#' @importFrom stats loess predict.loess
#' @importFrom mgcv gam s
ptr <- function(profile, conf) {
    profile <- copy(profile)
    w <- profile[, start[2] - start[1] + 1]
    profile <- profile[, reads := reads * conf$read_length / w]
    profile[reads <= conf$min_coverage, reads := NA]
    profile[abs(log(reads + 1) - log(median(reads + 1, na.rm = TRUE))) >
            log(conf$max_median_fold), reads := NA]
    if (profile[, (sum(!is.na(reads)) / .N) < conf$min_covered]) {
        return(list(
            ptr = NULL,
            profile = profile[, "smooth" := NA]
        ))
    }
    reads <- profile[!is.na(reads), reads]
    pos <- profile[!is.na(reads), start]
    data <- data.table(
        coverage = c(reads[1:(length(reads) - 1)],
                     reads, reads[2:length(reads)]),
        start = c(-rev(pos[2:length(pos)]), pos, max(pos) + pos[2:length(pos)])
    )
    fit <- gam(coverage ~ s(start), data = data)
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

peak_to_through <- function(object, ...) {
    if (!"coverage" %in% names(object)) {
        stop("Need `coverage` in artifacts!")
    }
    co <- object$coverage
    config <- config_parser(list(...), config_ptr)
    apfun <- parse_threads(config$threads)

    genbank_id <- unique(co[, list(genbank, id)])
    genbank_id <- lapply(1:nrow(genbank_id),
                         function(i) as.character(genbank_id[i]))
    flog.info(paste("Calculating smoothed coverage and PTR for %d genomes",
                    "and %d samples."), co[, uniqueN(genbank)],
                    co[, uniqueN(id)])
    ptrs <- apfun(genbank_id, function(row) {
        profile <- co[genbank == row[1] & id == row[2]]
        res <- ptr(profile, config)
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
    flog.info("Finished. Could get PTRs for %d genomes.",
              ptrs[, sum(!is.na(ptr))])

    artifact <- list(
        ptr = ptrs,
        coverage = profiles,
        steps = c(object[["steps"]], "peak_to_through")
    )
}
