# Download a list of required files

#' Download a list of files.
#'
#' @param files A data frame-like object with three columns: url, target,
#'  description specifying the source file, target location (including
#'  the file name) and description.
#' @param threads Maximum number of parallel file downloads.
#' @return The list of files with indicated download success.
#' @export
download_files <- function(files, threads = getOption("mc.cores", 1)) {
    downloaded <- mclapply(1:nrow(files), function(i) {
        meta <- files[i, ]
        if (!dir.exists(dirname(meta$target))) {
            dir.create(dirname(meta$target), recursive = TRUE)
        }
        flog.info("Downloading %s to %s.", meta$description, meta$target)
        ret <- tryCatch(download.file(meta$url, meta$target,
                                     quiet = TRUE),
                       error = function(e) 1)
        meta$success <- (ret == 0)
        if (meta$success) {
            flog.info("Finished downloading %s.", meta$target)
        } else {
            flog.info("Failed downloading %s.", meta$target)
            file.remove(meta$target)
        }
        return(meta)
    }, mc.cores = threads)

    return(rbindlist(downloaded))
}

ENA <- "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"

sra_download_url <- function(id, paired = TRUE) {
    url <- paste0(ENA, "/", substr(id, 1, 6))
    if (nchar(id) > 9) {
        url <- sprintf("%s/%03d", url, as.integer(substr(id, 10, 100)))
    }
    if (paired) {
        return(paste0(url, "/", id, "/", id, "_", 1:2, ".fastq.gz"))
    } else {
        return(paste0(url, "/", id, "/", id, ".fastq.gz"))
    }
}

# nocov start
prepare_filelist <- function(runtable, files, path = "data") {
    sra <- fread(runtable, sep = "\t")
    setkey(sra, Sample_Name)
    sra <- sra[files$id]
    paired <- "reverse" %in% names(files)
    urls <- sapply(sra$Run, sra_download_url, paired = paired)
    if (paired) {
        files <- data.table(
            url = c(urls[1, ], urls[2, ]),
            target = file.path(path, c(basename(files$forward),
                                       basename(files$reverse))),
            description = c(paste("forward reads for", files$id),
                            paste("reverse reads for", files$id))
            )
    } else {
        files <- data.table(
            url = urls,
            target = file.path(path, basename(files$forward)),
            description = c(paste("forward reads for", files$id))
            )
    }
    return(files)
}
# nocov end


sra_filelist <- function(runtable, path) {
    sra <- fread(runtable, sep = "\t")
    if (sra[, uniqueN(LibraryLayout) != 1]) {
        stop("Can't handle mixed library layouts :(")
    }
    paired <- sra[, unique(LibraryLayout) %>% tolower() == "paired"]

    urls <- sapply(sra$Run, sra_download_url, paired = paired)
    if (paired) {
        files <- data.table(
            url = c(urls[1, ], urls[2, ]),
            target = file.path(path, c(paste0(sra$Run, "_1.fastq.gz"),
                                       paste0(sra$Run, "_2.fastq.gz"))),
            description = c(paste("forward reads for", sra$Run),
                            paste("reverse reads for", sra$Run))
            )
    } else {
        files <- data.table(
            url = urls,
            target = file.path(path, paste0(sra$Run, "_1.fastq.gz")),
            description = c(paste("forward reads for", sra$Run))
            )
    }
    return(files)
}


#' Download data from SRA.
#'
#' To use this function first obtain a RunInfo Table as described at
#' \url{https://www.ncbi.nlm.nih.gov/Traces/study/?go=help}.
#'
#' @param runtable Path to a runtable as selected in the SRA interface.
#' @param path The folder into which to download the FASTQ files.
#' @param threads Maximum number of parallel file downloads.
#' @return The list of files with indicated download success.
#' @export
#' @importFrom futile.logger flog.warn
download_sra <- function(runtable, path = "sra/",
                         threads = getOption("mc.cores", 1)) {
    if (!dir.exists(path)) {
        flog.info("Creating directory `%s`.", path)
        dir.create(path, recursive = TRUE)
    }
    files <- sra_filelist(runtable, path)
    dl <- download_files(files, threads)
    if (dl[, sum(success)] != nrow(dl)) {
        flog.warn(paste0("%d/%d files could not be downloaded. ",
                         "They haven been marked in the returned table ",
                         "with `success == FALSE`."),
                  dl[, sum(!success)], nrow(dl))
    }
    return(dl)
}
