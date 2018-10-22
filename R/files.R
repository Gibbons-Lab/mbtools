# Helpers to manage read files

illumina_pattern <- "([A-Za-z0-9\\-]+)_S(\\d+)_L(\\d+)_R(\\d+)_001.fastq"

illumina <- function(dir) {
    files <- list.files(dir, pattern=illumina_pattern, recursive=TRUE,
                        include.dirs=TRUE)
    annotations <- as.data.table(str_match(files, illumina_pattern))
    names(annotations) <- c("file", "id", "injection_order", "lane",
                            "direction")
    annotations$file <- files
    if (annotations[, uniqueN(direction)] == 2) {
        fwd <- annotations[direction == "001"]
        bwd <- annotations[direction == "002"]
        annotations <- bwd[fwd[, .(reverse=file, id)], on="id"]
    } else {
        names(annotations)[1] <- "forward"
        annotations[, direction := NULL]
    }
    annotations[, injection_order := as.numeric(injection_order)]
    annotations[, lane := as.numeric(lane)]
    return(annotations)
}

#' Find Illumina read files in a given directory.
#'
#' @param directory The directory in which to look.
#' @param dirs_are_runs Whether subdirctories indicate different sequencing
#'  runs.
#' @return A data table that contains the samples and their annotations.
#'
#' @export
#' @importFrom stringr str_match_all
#' @importFrom data.table setDT uniqueN
find_illumina <- function(directory, dirs_are_runs=FALSE) {
    if (dirs_are_runs) {
        files <- list()
        for (dir in list.dirs(directory, recursive=FALSE, full.names=FALSE)) {
            fi <- illumina(file.path(directory, dir))
            fi[, "run" := dir]
            fi[, forward := file.path(dir, forward)]
            if ("reverse" %in% names(fi)) {
                fi[, reverse := file.path(dir, reverse)]
            }
            files[[dir]] <- fi
        }
        files <- rbindlist(files)
    } else {
        files <- illumina(directory)
    }
    files[, forward := file.path(directory, forward)]
    if ("reverse" %in% names(files)) {
        files[, reverse := file.path(directory, reverse)]
    }
    return(files)
}
