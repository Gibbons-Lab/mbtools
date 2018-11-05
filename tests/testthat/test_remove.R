context("contamination removal")
flog.threshold(WARN)

make_random_data <- function() {
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

    index_folder <- system.file("extdata/genomes", package = "mbtools")
    phix <- readFasta(file.path(index_folder, "phiX.fa.gz"))
    phix <- substr(as.character(sread(phix)[1]), 1, 100)
    sr <- ShortReadQ(sread = DNAStringSet(c(seqs[1:99], phix)),
                     quality = BStringSet(quals),
                     id = BStringSet(paste0("S", 1:100)))
    f2 <- file.path(d, "f2.fastq.gz")
    writeFastq(sr, f2)

    dset <- data.table(forward=file.path(d, c("f.fastq.gz", "f2.fastq.gz")),
                       reverse=file.path(d, c("r.fastq.gz", "f2.fastq.gz")),
                       index=file.path(d, c("f.fastq.gz", "f2.fastq.gz")),
                       id=c("random", "single"))
    return(dset)
}

test_that("sequences can be removed", {
    data <- make_random_data()
    d <- tempdir()
    dir.create(file.path(d, "out"), showWarnings = FALSE)
    outpath <- file.path(d, "out")

    index_folder <- system.file("extdata/genomes", package = "mbtools")
    reads <- as.character(data[1, .(forward, reverse)])
    index <- data[1, index]
    counts <- remove_reference(reads, out = outpath,
                reference = file.path(index_folder, "phiX.fa.gz"),
                index = index)
    expect_equal(counts[c("reads", "removed")], list(reads = 100, removed = 0))
    expect_equal(nrow(counts$counts), 0)

    counts <- remove_reference(reads[1], out = outpath,
                reference = file.path(index_folder, "phiX.fa.gz"))
    expect_equal(counts[c("reads", "removed")], list(reads = 100, removed = 0))
    expect_equal(nrow(counts$counts), 0)

    reads <- as.character(data[2, .(forward, reverse)])
    counts <- remove_reference(reads, out = outpath,
                reference = file.path(index_folder, "phiX.fa.gz"))
    expect_equal(counts[c("reads", "removed")], list(reads = 100, removed = 1))
    expect_equal(counts$counts[, sum(counts)], 1)
})

test_that("filtering works on full data sets", {
    data <- make_random_data()
    d <- tempdir()
    dir.create(file.path(d, "out"), showWarnings = FALSE)
    outpath <- file.path(d, "out")
    index_folder <- system.file("extdata/genomes", package = "mbtools")

    expect_output(counts <- filter_reference(data, out = outpath,
                        reference = file.path(index_folder, "phiX.fa.gz")),
                  "100% elapsed")
    expect_equal(counts[, sum(counts)], 1)
    expect_equal(counts[, uniqueN(id)], 1)
})
flog.threshold(INFO)
