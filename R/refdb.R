# Tools for dealing with 16S reference databases

#' Build the lineage from a SILVA path string.
#'
#' @param silva_path A string containing a taxonomy path with taxa separated
#'  by semicolons.
#' @param taxid The overall taxid for the path/taxon.
#' @param taxonomy The SILVA taxonomy as data table.
#' @param return A data table containing the individual taxa ranks with ids.
#'
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
#' @param taxonomy The SILVA taxonomy.
#' @return The SILVA lineage mapping every taxon to its upper ranks up to the
#'  root.
#'
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


ENSEMBL_CDNA <- paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/current/",
                       "fasta/%s/%s/cdna/%s.%s.cdna.all.fa.gz")

#' Helper function to download ENSEMBL transcripts (cdna).
#' @importFrom stringr str_to_title
download_ensembl_cdna <- function(out="transcripts", collection,
                                  name, assembly) {
    collection <- str_split_fixed(collection, "_core", n=2)[1]

    url <- sprintf(ENSEMBL_CDNA, collection, name,
                   str_to_title(name), assembly)
    outfile <- file.path(out, paste0(name, ".cdna.all.fa.gz"))
    if (!file.exists(outfile)) {
        ret <- download.file(url, outfile, quiet=TRUE)
    }
    if (ret != 0) {
        stop(paste("could not download", name))
    }
    return(outfile)
}

#' Download the entire ENSEMBL transcript DB for bacteria.
download_bacteria_transcripts <- function(out="transcripts", np=8,
                                          remove_redundant=TRUE) {
    dir.create(out, showWarnings=FALSE)
    db <- fread(paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/current/",
                       "species_EnsemblBacteria.txt"),
                skip=1, header=FALSE,
                col.names=c("name", "species", "division", "taxonomy_id",
                            "assembly", "assembly_accession", "genebuild",
                            "variation", "pan_compara", "peptide_compara",
                            "genome_alignments", "other_alignments",
                            "core_db", "species_id", "trash"))
    if (remove_redundant) {
        db <- db[order(-genebuild), .SD[1], by="taxonomy_id"]
        cat(sprintf("Downloading %d transcript files.", nrow(db)))
    }
    downloaded <- pbapply(db, 1, function(row) {
        download_ensembl_cdna(out, row["core_db"], row["species"],
                              row["assembly"])
    }, cl=np)
    return(downloaded)
}
