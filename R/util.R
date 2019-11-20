# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.


#' Orlitsky's diminishing attenuation estimator (q2/3).
#'
#' @param counts A vector of counts for which to approximate discrete
#'  probablities.
#' @return A named vector `p` assigning a probability to each event in data.
#'  Those do not sum up to one since there is also a remaining probability to
#'  observe a new event given as p(new) = 1 - sum(p).
#' @examples
#'  x <- sample(1:10, 100, replace=TRUE)
#'  p <- orlitsky(x)
#' @export
orlitsky <- function(counts) {
    n <- length(counts)
    phi <- c(tabulate(counts), 0) # the prevalences, denoted by phi
    cn <- ceiling((n+1)^1/3)    # a smoothing factor for the prevalences
    new <- max(cn, phi[1] + 1)
    probs <- (counts + 1) * pmax(cn, phi[counts + 1] + 1) /
             pmax(cn, phi[counts])
    if (!is.null(names(counts))) {
        names(probs) <- names(counts)
    }
    return(probs / sum(c(probs, new)))
}


#' Read a blast hit/alignment file.
#'
#' @param matches The file containing the matches like it is returned by
#'  blast or diamond.
#' @return The hits as a data table.
#' @export
read_blast <- function(matches) {
    flog.info("Reading blast file %s.", matches)
    align <- fread(matches, header = FALSE)
    names(align) <- c("query", "reference", "percent_match", "alignment_length",
                      "num_mismatches", "num_gap_open", "query_start",
                      "query_end", "ref_start", "ref_end", "evalue",
                      "bit_score")
    flog.info("Read %d hits in total (max E=%g).",
              nrow(align), align[, max(evalue)])
    return(align)
}


#' Very simple select for mbtools artifacts. Just gets a particular entry.
#'
#' @param object A mbtools artifact returned from a workflow step.
#' @param entry The name of the entry to get.
#' @return The requested entry.
#' @export
select <- function(object, entry) {
    if (!entry %in% names(object)) {
        stop(sprintf("This object has no entry named `%s` :(", entry))
    }
    return(object[[entry]])
}

#' @importFrom parallel detectCores mclapply
parse_threads <- function(th, return_fun = TRUE) {
    if (is.logical(th) & th) {
        threads <- detectCores()
    } else if (th > 1) {
        threads <- th
    } else {
        threads <- NA
    }
    if (!is.na(threads) & threads > 1) {
        apfun <- function(...) mclapply(..., mc.cores = threads)
    } else {
        apfun <- lapply
        threads <- 1
    }

    if (return_fun) {
        return(apfun)
    }
    return(threads)
}

`%p%` <- function(a, b) paste0(a, b)

split_ext <- function(path) {
    path <- basename(path)
    return(str_split_fixed(path, ".", 2) %>% as.character())
}
