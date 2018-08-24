context("taxonomy conversion")

REPO <- "https://github.com/microbiome/HITdb/raw/master/HITdb_v1.00/"

test_that("we can convert HITdb", {
    d <- tempdir()
    download.file(paste0(REPO, "HITdb_sequences.fna"), file.path(d, "hitdb.fa"),
        quiet = TRUE)
    download.file(paste0(REPO, "HITdb_taxonomy_mothur.txt"),
        file.path(d, "taxa.tsv"), quiet = TRUE)

    mothur_to_dada(file.path(d, "hitdb.fa"), file.path(d, "taxa.tsv"),
        file.path(d, "taxonomy.fa.gz"))

    expect_true(file.exists(file.path(d, "taxonomy.fa.gz")))
    seqs <- readFasta(file.path(d, "taxonomy.fa.gz"))
    ids <- as.character(id(seqs))
    nlev <- sapply(strsplit(ids, ";", fixed = TRUE), length)

    expect_true(all(nlev <= 6))
})
