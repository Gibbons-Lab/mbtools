# Helpers to manage read files

illumina_pattern <- "([A-Za-z0-9\\-]+)_S(\\d+)(?:_L(\\d+))*_R(\\d+)_001.fastq"
illumina_annotations <- c("id", "injection_order", "lane", "direction")

annotate_files <- function(dir, pattern, annotations) {
    if (!"id" %in% annotations) {
        stop("need at least the id for each sample")
    }
    files <- list.files(dir, pattern = pattern, recursive = TRUE,
                        include.dirs = TRUE)
    anns <- as.data.table(str_match(files, pattern))
    names(anns) <- c("file", annotations)
    if ("direction" %in% annotations) {
        anns[, direction := as.numeric(direction)]
    } else {
        anns[, direction := 1]
    }
    anns$file <- files
    names(anns)[1] <- "forward"
    if (anns[, uniqueN(direction)] == 2) {
        fwd <- anns[direction == 1]
        bwd <- anns[direction == 2]
        anns <- fwd[bwd[, .(reverse = forward, id)], on = "id"]
        other_cols <- names(anns)[!names(anns) %in% c("forward", "reverse")]
        anns <- anns[, c("forward", "reverse", other_cols), with = FALSE]
    }
    anns[, direction := NULL]
    if ("injection_order" %in% names(anns)) {
        anns[, injection_order := as.numeric(injection_order)]
    }
    if ("lane" %in% names(anns)) {
        anns[, lane := as.numeric(lane)]
    }
    return(anns)
}

get_files <- function(obj) {
    if ("data.table" %in% class(obj) &&
        all(c("forward", "id") %in% names(obj))) {
            return(obj)
    } else if ("list" %in% class(obj) && "files" %in% names(obj)) {
        return(obj[["files"]])
    }
    stop("`object` must be a file list or a workflow object :/")
}

get_alignments <- function(obj) {
    if ("data.table" %in% class(obj) &&
        all(c("alignments", "id") %in% names(obj))) {
            return(obj)
    } else if ("list" %in% class(obj) && "alignments" %in% names(obj)) {
        return(obj[["alignments"]])
    }
    stop("`object` must be an alignment list or a workflow object :/")
}

#' Find read files in a given directory.
#'
#' @param directory The directory in which to look.
#' @param dirs_are_runs Whether subdirctories indicate different sequencing
#'  runs.
#' @param pattern Regular expression pattern for the file names. Each capture
#'  group will be used as an annotation.
#' @param annotations Names for the annotations. Must contain one name for each
#'  capture group in `pattern`.
#' @return A data table that contains the samples and their annotations.
#'
#' @export
#' @importFrom stringr str_match_all
#' @importFrom data.table setDT uniqueN
find_read_files <- function(directory,
                            pattern = illumina_pattern,
                            annotations = illumina_annotations,
                            dirs_are_runs = FALSE) {
    if (dirs_are_runs) {
        files <- list()
        for (dir in list.dirs(directory, recursive = FALSE,
                              full.names = FALSE)) {
            fi <- annotate_files(file.path(directory, dir), pattern,
                                 annotations)
            fi[, "run" := dir]
            fi[, forward := file.path(dir, forward)]
            if ("reverse" %in% names(fi)) {
                fi[, reverse := file.path(dir, reverse)]
            }
            files[[dir]] <- fi
        }
        files <- rbindlist(files)
    } else {
        files <- annotate_files(directory, pattern, annotations)
    }
    files[, forward := file.path(directory, forward)]
    if ("reverse" %in% names(files)) {
        files[, reverse := file.path(directory, reverse)]
    }
    return(files)
}
