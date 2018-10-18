# Tools for dealing with 16S reference databases

#' Build the lineage from a SILVA path string.
#' importFrom stringr str_split_fixed
build_lineage <- function(silva_path, taxid, taxonomy) {
    levs <- str_split_fixed(silva_path, ";", n=Inf)[1, ]
    n <- length(levs)
    lineage <- sapply(1:(n-1), function(i) {
        p <- paste0(c(levs[1:i], ""), collapse=";")
        return(p)
    })
    lineage <- taxonomy[lineage, .(rank=rank, rank_taxid=taxid)]
    lineage[, "level" := 1:(n-1)]
    lineage[, "name" := levs[-n]]
    lineage <- rbind(lineage,
                     data.table(rank="species", rank_taxid=taxid,
                                level=n, name=levs[n]))
    return(lineage)
}

#' Build the lineages for the entire SILVA database.
#'
#' @param taxmap The SILVA taxonomy map (`taxmap_slv_*.txt`).
#' @param taxonomy The SILVA taxonomy
#' @param importFrom data.table rbind.data.table fread unique.data.table
#' @export
silva_build_taxonomy <- function(taxmap, taxonomy) {
    map <- fread(taxmap)
    map <- unique(map[, .(taxid, path, organismName)])
    tax <- fread(taxonomy,
                 col.names=c("path", "taxid", "rank", "comments", "version"),
                 key="path")
    lineage <- map[, build_lineage(paste0(path, organismName), taxid, tax),
                   by="taxid"]
    lineage[name == "", name := NA]
    return(lineage)
}
