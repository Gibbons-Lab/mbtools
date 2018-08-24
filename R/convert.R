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
        if (as.numeric(match[i, 4]) > cutoff)
            x <- sub(match[i, 2], match[i, 3], x)
        else x <- sub(match[i, 2], "unclassified", x)
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
#'  NULL
#'
#' @export
mothur_to_dada <- function(seq_file, taxa_file, out = "taxonomy.fa.gz") {
    taxa_df <- read.table(taxa_file, header = FALSE)
    matches <- str_match(taxa_df[, 2], ANN_RE)

    tax <- vapply(1:nrow(taxa_df), hitdb_cleaner, "", df = taxa_df,
                   match = matches)
    names(tax) <- taxa_df[, 1]

    seqs <- readFasta(seq_file)
    ids <- as.character(id(seqs))
    seqs <- ShortRead(sread(seqs), BStringSet(tax[ids]))
    writeFasta(seqs, out, compress = TRUE)
}
