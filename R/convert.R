# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# Apache license 2.0.

ANN_RE <- ";([\\w_.]+)\\|NN=(.+)\\|D=(.+);"

# cutoff larger than 100 means never use nearest-neighbor
hitdb_cleaner <- function(i, df, match, cutoff = 101) {
    xn <- df[i, 1]
    x <- df[i, 2]

    x <- sub("\\|.+;", ";", x)
    if (all(!is.na(match[i, ]))) {
        if (as.numeric(match[i, 4]) > cutoff) {
            x <- sub(match[i, 2], match[i, 3], x)
        } else {
            x <- sub(match[i, 2], "unclassified", x)
        }
    }

    # dada2 nows about empty fields so we do not need placeholders
    gsub("[Uu]nclassified;", "", x)
}

#' Converts taxa annotations from mothur format to dada2 format.
#'
#' @param seq_file A fasta file containing the (cluster) sequences.
#' @param taxa_file A tab-separated file with IDs on the first column and
#'  taxonomy in the second column.
#' @param out Filename for the compressed output file.
#' @return Nothing.
#' @examples
#' NULL
#'
#' @export
mothur_to_dada <- function(seq_file, taxa_file, out = "taxonomy.fa.gz") {
    taxa_df <- read.table(taxa_file, header = FALSE)
    matches <- str_match(taxa_df[, 2], ANN_RE)

    tax <- vapply(1:nrow(taxa_df), hitdb_cleaner, "",
        df = taxa_df,
        match = matches
    )
    names(tax) <- taxa_df[, 1]

    seqs <- readFasta(seq_file)
    ids <- as.character(id(seqs))
    seqs <- ShortRead(sread(seqs), BStringSet(tax[ids]))
    writeFasta(seqs, out, compress = TRUE)
}

#' Converts BRACKEN results to a phyloseq object.
#'
#' @param bracken A data.table containing the BRACKEN quantifications.
#' @param metadata A metadata in data.table format to use.
#' @param id_col Column in the metdata corresponding to the sample IDs
#' in the BRACKEN counts.
#' @return A phyloseq object for the data.
#' @examples
#' NULL
#'
#' @export
#' @importFrom data.table setkeyv
bracken_to_phyloseq <- function(
        bracken,
        metadata = NULL,
        id_col = "sample_id") {
    n_ranks <- which(names(bracken) == "reads") - 1
    lowest_rank <- names(bracken)[n_ranks]
    table <- dcast(
        bracken,
        reformulate(lowest_rank, response = "sample"),
        value.var = "reads",
        fill = 0,
        fun.aggregate = sum
    )
    samps <- table[, sample]
    table <- as.matrix(table[, !"sample"])
    rownames(table) <- samps
    taxa <- unique(bracken[, 1:n_ranks, with = FALSE])
    setkeyv(taxa, lowest_rank)
    taxa <- as.matrix(taxa[colnames(table)])
    rownames(taxa) <- taxa[, lowest_rank]
    colnames(taxa) <- names(bracken)[1:n_ranks]

    if (is.null(metadata)) {
        ps <- phyloseq(
            otu_table(table, taxa_are_rows = FALSE),
            tax_table(taxa)
        )
        return(ps)
    }

    metadata <- as.data.frame(metadata)
    rownames(metadata) <- metadata[[id_col]]
    ps <- phyloseq(
        otu_table(table, taxa_are_rows = FALSE),
        tax_table(taxa),
        sample_data(metadata)
    )
    return(ps)
}


#' Converts read count data to a phyloseq object.
#'
#' @param counts A data.table containing at least columns "reads",
#'   "sample_id", and features.
#' @param metadata A metadata in data.table format to use.
#' @param feature_col Column in the count table containing the features.
#' @param id_col Column in the metadata corresponding to the sample IDs
#' in the read counts.
#' @return A phyloseq object for the data.
#' @examples
#' NULL
#'
#' @export
reads_to_phyloseq <- function(
    counts,
    feature_col,
    metadata = NULL,
    id_col = "sample_id") {
    table <- dcast(
        counts,
        reformulate(feature_col, response = "sample_id"),
        value.var = "reads",
        fill = 0
    )
    samps <- table$sample_id
    table <- as.matrix(table[, !"sample_id"])
    rownames(table) <- samps
    features <- matrix(colnames(table), ncol = 1)
    rownames(features) <- features[, 1]
    colnames(features) <- feature_col

    if (is.null(metadata)) {
        ps <- phyloseq(
            otu_table(table, taxa_are_rows = FALSE),
            tax_table(features)
        )
        return(ps)
    }

    metadata <- as.data.frame(metadata)
    rownames(metadata) <- metadata[[id_col]]
    ps <- phyloseq(
        otu_table(table, taxa_are_rows = FALSE),
        tax_table(features),
        sample_data(metadata)
    )
    return(ps)
}
