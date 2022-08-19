# Helpers to manage read files

illumina_pattern <-
    "([A-Za-z0-9\\-\\.]+)_S(\\d+)(?:_trimmed)*(?:_L(\\d+))*_R(\\d+)_001.f"
illumina_annotations <- c("id", "injection_order", "lane", "direction")

sra_pattern <- "([A-Za-z0-9\\-]+)_(\\d).fastq"
sra_annotations <- c("id", "direction")

simple_pattern <- "([A-Za-z0-9_\\-\\.]+)\\.f"
simple_annotations <- c("id")

annotate_files <- function(dir, pattern, annotations) {
    if (!"id" %in% annotations) {
        stop("need at least the id for each sample")
    }
    files <- list.files(dir, recursive = TRUE, include.dirs = TRUE)
    match <- str_match(files, pattern)
    files <- files[!is.na(match[, 1])]
    match <- match[!is.na(match[, 1]), ]
    anns <- as.data.table(match)
    names(anns) <- c("file", annotations)
    if ("direction" %in% annotations) {
        anns[, direction := as.numeric(direction)]
    } else {
        anns[, direction := 1]
    }
    matchcols <- names(anns)[!names(anns) %in% c("file", "direction")]
    anns$file <- files
    names(anns)[1] <- "forward"
    dupes <- duplicated(anns)
    if (any(dupes)) {
        flog.error(paste("Some files have duplicated metadata. Please fix",
                        "the filenames or your pattern.",
                        "Duplicated files: %s"),
                   paste(unique(anns$forward[dupes])))
        stop("Duplicated file information. Can not continue :(")
    }
    if (anns[, uniqueN(direction)] == 2) {
        fwd <- anns[direction == 1]
        bwd <- anns[direction == 2]
        names(bwd)[1] <- "reverse"
        anns <- merge(fwd, bwd, by = matchcols, all = TRUE)
        other_cols <- names(anns)[
            !names(anns) %in% c("forward", "reverse") &
            !grepl("i\\.", names(anns))]
        anns <- anns[, c("forward", "reverse", other_cols), with = FALSE]
        anns[, "direction.x" := NULL]
        anns[, "direction.y" := NULL]
    }

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
        all(c("alignment", "id") %in% names(obj))) {
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
#' @importFrom utils capture.output
find_read_files <- function(directory,
                            pattern = illumina_pattern,
                            annotations = illumina_annotations,
                            dirs_are_runs = FALSE) {
    subdirs <- list.dirs(directory, recursive = FALSE, full.names = FALSE)
    if (dirs_are_runs && (length(subdirs) > 0)) {
        files <- list()
        for (dir in subdirs) {
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
    files[, forward :=
        ifelse(is.na(forward), NA, file.path(directory, forward))]
    if ("reverse" %in% names(files)) {
        files[, reverse :=
            ifelse(is.na(reverse), NA, file.path(directory, reverse))]
    }
    tab <- table(files$id)
    if (any(tab > 1)) {
        flog.warn(paste0("There are duplicated ids, please fix them before ",
                         "advancing to other workflow steps. Duplicated: %s"),
                  paste(names(tab)[tab > 1], collapse = ", "))
    }
    missing <- files[is.na(forward) | is.na(reverse)]
    if (nrow(missing) > 0) {
        flog.warn(
            "The following samples have missing read files: %s",
            capture.output(missing)
        )
    }
    return(files)
}
