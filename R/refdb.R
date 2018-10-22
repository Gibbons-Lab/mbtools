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

ens_title <- function(name) {
    return(paste0(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name))))
}

# Try to reproduce the ENSEMBL modifications
ens_assembly <- function(assembly) {
    assembly <- gsub("[^[:alnum:]\\.\\-]", "_", assembly)
    assembly <- gsub("__", "_", assembly, fixed=TRUE)
    assembly <- gsub("Cand.", "Cand", assembly, fixed=TRUE)
    assembly <- gsub("sp.", "sp", assembly, fixed=TRUE)
    return(assembly)
}

#' Helper function to download ENSEMBL transcripts (cdna).
#' @importFrom stringr str_to_title
download_ensembl_cdna <- function(out="transcripts", collection,
                                  name, assembly) {
    collection <- str_split_fixed(collection, "_core", n=2)[1]
    assembly <- ens_assembly(assembly)
    url <- sprintf(ENSEMBL_CDNA, collection, name,
                   ens_title(name), assembly)
    outfile <- file.path(out, paste0(name, ".", assembly, ".cdna.all.fa.gz"))
    ret <- 0
    if (!file.exists(outfile)) {
        ret <- tryCatch(download.file(url, outfile, quiet=TRUE),
                        error=function(e) 1,
                        warning=function(w) 1)
    }
    if (ret != 0) {
        file.remove(outfile)
        warning(paste("could not download", url))
        return(NA)
    }
    return(name)
}

#' Download the entire ENSEMBL transcript DB for bacteria.
#'
#' @param out Where to store the downloaded transcript files.
#' @param np How many parallel download processes to use.
#' @param remove_redundant Whether to remove alternative assemblies for the
#'  same strain.
#' @return A data table containing the downloaded species and assemblies.
#' @importFrom data.table setkey
download_bacterial_transcripts <- function(out="transcripts", np=8,
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
        flog.info("Downloading %d transcript files.", nrow(db))
    }
    downloaded <- pbapply(db, 1, function(row) {
        download_ensembl_cdna(out, row["core_db"], row["species"],
                              row["assembly"])
    }, cl=np)
    setkey(db, "species")
    downloaded <- as.character(downloaded[!is.na(downloaded)])
    flog.info("Downloaded %d transcript databases.", length(downloaded))
    downloaded <- db[downloaded]
    downloaded[, "file" := paste0(species, ".", ens_assembly(assembly),
                                  ".cdna.all.fa.gz")]
    return(downloaded)
}

#' Merge transcripts into a single database
#'
#' @param transcript_files A data table as returned by
#'  \code{\link{download_bacterial_trasncripts}}.
#' @param transcript_folder In which folder to look for transcripts.
#' @param out The filename for the DB. Will be compressed and saved in the
#'  transcript folder.
#' @return The transcript counts for each file.
merge_transcripts <- function(transcripts_files, transcripts_folder,
                              out="ensembl_transcripts.fa.gz") {
    flog.info("Merging %d transcript files from %s.",
              nrow(transcripts_files), transcripts_folder)
    counts <- pbapply(transcripts_files, 1, function(row) {
        fasta <- readFasta(file.path(transcripts_folder, row["file"]))
        new_ids <- BStringSet(paste0("TAX", row["taxonomy_id"],
                                     "_", id(fasta)))
        writeFasta(ShortRead(sread(fasta), new_ids),
                   file.path(transcripts_folder, out),
                   mode="a", compress=TRUE)
        length(fasta)
    })

    return(counts)
}
