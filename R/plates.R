# Helpers to arrange samples on plates

#' Build a configuration for the plate layout.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the long read alignment
#'  workflow.
#' @export
#' @examples
#'  config <- config_align(reference = "refs/mouse.fna.gz")
config_layout <- config_builder(list(
    idcol = "id",
    blank_step = 18,
    ncol = 1,
    by = "col",
    fill = FALSE,
    ncols = 12,
    nrows = 8
))

grid <- function(ncols, nrows) {
    list(
        row = expand.grid(as.character(1:nrows), LETTERS[1:ncols])[, 2:1] %>%
            apply(1, paste, collapse = "", sep = ""),
        col = expand.grid(LETTERS[1:ncols], as.character(1:nrows)) %>%
            apply(1, paste, collapse = "", sep = "")
    )
}

#' Generate a plate layout from a sample manifest.
#'
#' @param manifest A data frame-like object containing the samples.
#' @param ... other parameters passed to \code{\link{config_layout}}.
#' @return An artifact with the annotated manifest as well as a plate map.
#'
#' @export
#' @importFrom stringr str_pad
#' @importFrom magrittr %>%
#' @importFrom data.table copy setnames
layout <- function(manifest, ...) {
    config <- config_parser(list(...), config_layout)
    blank <- data.table(layout_type = "blank")
    manifest <- as.data.table(copy(manifest))[, "layout_type" := "sample"]
    setnames(manifest, old = config$idcol, new = "id")

    g <- grid(config$ncols, config$nrows)
    n <- config$ncols * config$nrows

    if (("well" %in% names(manifest)) || ("plate" %in% names(manifest))) {
        config$blank_step <- Inf
    }

    blank[[config$idcol]] <- ""
    dt <- list()
    if (is.finite(config$blank_step)) {
        for (i in seq(1, nrow(manifest), config$blank_step)) {
            n <- min(i + config$blank_step - 1, nrow(manifest))
            if (i == 1) {
                dt <- append(dt, list(manifest[i:n]))
            } else {
                dt <- append(dt, list(blank, manifest[i:n]))
            }
        }
        manifest <- rbindlist(dt, fill = TRUE)
    }
    if (!"plate" %in% names(manifest) ) {
        manifest[, "plate" := ceiling(1:nrow(manifest) / n)]
    }
    if (!"well" %in% names(manifest)) {
        manifest[, "well" := g[[config$by]][1:.N], by = "plate"]
    }
    if (config$fill) {
        blanks <- data.table(layout_type = "blank", well = g[[config$by]])
        all_blanks = list()
        for (i in manifest[, unique(plate)]) {
            bi <- copy(blanks[!well %in% manifest$well])
            bi[, "plate" := i]
            all_blanks[[i]] <- bi
        }
        all_blanks <- rbindlist(all_blanks)
        manifest <- rbind(manifest, all_blanks, fill = TRUE)[order(plate, well)]
    }

    layout <- ggplot(
        manifest,
        aes(x = factor(substr(well, 2, 3), levels = 1:config$ncols),
            y = substr(well, 1, 1),
            fill = layout_type,
            label = id)) +
        scale_x_discrete(drop = FALSE) +
        geom_tile(color = "black", size = 0.5) +
        geom_text(hjust = 0.5, vjust = 0.5) +
        facet_wrap(~ plate, scales = "free",
                   labeller = function(x) label_both(x, sep = " "),
                   ncol = config$ncol) +
        scale_fill_manual(values = c(blank="gray", sample="white")) +
        guides(fill = guide_legend(nrow = 1)) +
        scale_y_discrete(limits = rev(LETTERS[1:config$nrows])) +
        labs(x = "", y = "", fill = "") +
        theme_minimal(base_size = 16) + theme(legend.position = "bottom")

    artifact <- list(manifest = manifest, layout = layout, steps = c("layout"))
    return(artifact)
}
