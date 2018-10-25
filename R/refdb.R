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
ENSEMBL_PROTEIN <- paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/current/",
                          "fasta/%s/%s/pep/%s.%s.pep.all.fa.gz")

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
#' @importFrom stringr str_to_title str_split_fixed
download_ensembl <- function(out="transcripts", what="cdna", collection,
                             name, assembly) {
    collection <- str_split_fixed(collection, "_core", n=2)[1]
    assembly <- ens_assembly(assembly)
    outfile <- file.path(out, paste0(name, ".", assembly))
    if (what == "cdna") {
        url <- sprintf(ENSEMBL_CDNA, collection, name,
                       ens_title(name), assembly)
        outfile <- paste0(outfile, ".cdna.all.fa.gz")
    } else {
        url <- sprintf(ENSEMBL_PROTEIN, collection, name,
                       ens_title(name), assembly)
        outfile <- paste0(outfile, ".pep.all.fa.gz")
    }

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
#' @param what Download in which format. Can be "cdna" or "protein" to download
#'  reverse transcribed DNA (cDNA) or protein sequences.
#' @param remove_redundant Whether to remove alternative assemblies for the
#'  same strain.
#' @return A data table containing the downloaded species and assemblies.
#' @importFrom data.table setkey
download_bacterial_transcripts <- function(out="transcripts", what="cdna",
                                           np=8, remove_redundant=TRUE) {
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
        flog.info("Downloading %d %s transcript files.", nrow(db), what)
    }
    downloaded <- pbapply(db, 1, function(row) {
        download_ensembl(out, what, row["core_db"], row["species"],
                         row["assembly"])
    }, cl=np)
    setkey(db, "species")
    downloaded <- as.character(downloaded[!is.na(downloaded)])
    flog.info("Downloaded %d %s transcript databases.",
              length(downloaded), what)
    downloaded <- db[downloaded]
    if (what == "cdna") {
        downloaded[, "file" := paste0(species, ".", ens_assembly(assembly),
                                      ".cdna.all.fa.gz")]
    } else {
        downloaded[, "file" := paste0(species, ".", ens_assembly(assembly),
                                      ".pep.all.fa.gz")]
    }
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
#' @importFrom Biostrings readAAStringSet readDNAStringSet writeXStringSet
merge_transcripts <- function(transcripts_files, transcripts_folder,
                              what="cdna", out="ensembl_transcripts.fa.gz") {
    flog.info("Merging %d %s transcript files from %s.",
              nrow(transcripts_files), what, transcripts_folder)
    if (what == "cdna") {
        reader <- readDNAStringSet
    } else {
        reader <- readAAStringSet
    }
    out <- file.path(transcripts_folder, out)
    if (file.exists(out)) {
        flog.info("Overwriting output %s with new data.", out)
        file.remove(out)
    }
    counts <- pbapply(transcripts_files, 1, function(row) {
        fasta <- reader(file.path(transcripts_folder, row["file"]))
        names(fasta) <- paste0("TAX", as.integer(row["taxonomy_id"]), "_",
                               names(fasta))
        writeXStringSet(fasta, out, append=TRUE, compress=TRUE)
        length(fasta)
    })

    return(counts)
}

ENSID <- paste0("TAX(\\d+)\\_(\\w+) (\\w+) supercontig:(.+) gene:(.+) ",
                "transcript:(.+) gene_biotype:(.+) transcript_biotype:(.+) ",
                "description:(.+)")

#' Parse annotations from an ENSEMBL id
#'
#' @param id The id to be parsed.
#' @param A data table containing the transcript id with annotations.
parse_ensembl_id <- function(id) {
    res <- as.data.table(str_match(id, ENSID)[, 2:10])
    names(res) <- c("taxid", "seqid", "sequence_type", "supercontig", "gene",
                    "transcript", "gene_biotype", "transcript_biotype",
                    "description")
    res[, "id" := paste0("TAX", as.integer(taxid), "_", seqid)]
    return(res)
}

#' Annotate diamond database hits with transcript info.
#'
#' @param matches Path to diamond output (*.m8).
#' @return The annotated hits as data table.
#' @importFrom Biostrings fasta.index
#' @export
annotate_contigs <- function(matches) {
    flog.info("Reading contig-reference alignments.")
    align <- fread(matches, header=FALSE)
    names(align) <- c("query", "reference", "percent_match", "alignment_length",
                      "num_mismatches", "num_gap_open", "query_start",
                      "query_end", "ref_start", "ref_end", "evalue",
                      "bit_score")
    flog.info("Getting unique reference hits...")
    ids <- align[, unique(reference)]
    flog.info("Parsing annotations for %d sequences...", length(ids))
    anns <- parse_ensembl_id(ids)
    flog.info("Merging hits with annotations...")
    merged <- anns[align, by=c(id="reference")]
    return(merged)
}
