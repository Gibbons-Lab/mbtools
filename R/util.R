# Copyright 2017 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0. See LICENSE for more information.


#' List samples and read files in a directory.
#'
#' @param path Where to look for files.
#' @return A data frame with columns "id", "forward" and optionally "reverse".
#' @examples
#'  NULL
#'
#' @export
#' @importFrom stringr str_match
#' @importFrom data.table data.table
sra_files <- function(path) {
    files <- list.files(path, pattern=".fastq\\.*gz*", full.names=TRUE)
    fwd <- grepl("_1.fastq", files)
    rev <- grepl("_2.fastq", files)
    if (any(fwd) && sum(fwd) != sum(rev)) {
        stop("Some paired files are missing!")
    }

    ids <- str_match(files, "([a-zA-Z\\d]+)_*\\d*\\.fastq")[, 2]

    if (any(fwd)) {
        return(data.table(id=ids[fwd], forward=files[fwd], reverse=files[rev]))
    }
    return(data.table(id=ids, forward=files))
}

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
