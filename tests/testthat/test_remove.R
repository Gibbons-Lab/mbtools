context("contamination removal")

test_that("index files can be downloaded", {
    with_mock(
        download.file = function(...) print("downloading"),
        expect_output(download_index(genome_file = "blu/bla.file"),
                      "downloading")
    )
})

test_that("sequences can be removed", {
    d <- tempdir()
    seqs <- replicate(100,
        paste(sample(c("A", "C", "G", "T"), 100, replace = TRUE),
                     collapse = ""))
    quals <- replicate(100, paste(rep("@", 100), collapse = ""))
    sr <- ShortReadQ(sread = DNAStringSet(seqs),
                     quality = BStringSet(quals),
                     id = BStringSet(paste0("S", 1:100)))
    old_files <- list.files(d, "fastq", recursive = TRUE, full.names = TRUE)
    file.remove(old_files)
    writeFastq(sr, file.path(d, "f.fastq.gz"))
    writeFastq(sr, file.path(d, "r.fastq.gz"))
    writeFastq(sr, file.path(d, "i.fastq.gz"))

    reads <- file.path(d, c("f.fastq.gz", "r.fastq.gz"))
    index <- file.path(d, "i.fastq.gz")
    dir.create(file.path(d, "nh"), showWarnings = FALSE)

    index_folder <- system.file("extdata/indices", package = "mbtools")
    counts <- remove_organism(reads, index, file.path(d, "nh"),
                              organism = "lambda_virus", where = index_folder)
    expect_equal(counts, c(reads = 100, removed = 0))

})
