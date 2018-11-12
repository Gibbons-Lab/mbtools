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
    sr_rev <- ShortReadQ(sread = reverseComplement(sread(sr)),
                         quality = reverse(quality(sr)),
                         id = id(sr))
    f2 <- file.path(d, "f2.fastq.gz")
    r2 <- file.path(d, "r2.fastq.gz")
    writeFastq(sr, f2)
    writeFastq(sr_rev, r2)

    dset <- data.table(forward=file.path(d, c("f.fastq.gz", "f2.fastq.gz")),
                       reverse=file.path(d, c("r.fastq.gz", "r2.fastq.gz")),
                       index=file.path(d, c("i.fastq.gz", "i2.fastq.gz")),
                       id=c("random", "single"))
    return(dset)
}

test_that("sequences can be removed", {
    data <- make_random_data()
    d <- tempdir()
    outpath <- file.path(d, "out")
    unlink(file.path(outpath), recursive=TRUE)
    dir.create(outpath, showWarnings = FALSE)

    index_folder <- system.file("extdata/genomes", package = "mbtools")
    reads <- as.character(data[1, .(forward, reverse)])
    index <- data[1, index]
    counts <- remove_reference(reads, out = outpath,
                reference = file.path(index_folder, "phiX.fa.gz"),
                index = index)
    expect_equal(counts, list(reads = 100, removed = 0))
    unlink(file.path(outpath, "*"), recursive=TRUE)

    counts <- remove_reference(reads[1], out = outpath,
                reference = file.path(index_folder, "phiX.fa.gz"))
    expect_equal(counts, list(reads = 100, removed = 0))
    unlink(file.path(outpath, "*"), recursive=TRUE)

    reads <- as.character(data[2, .(forward, reverse)])
    counts <- remove_reference(reads, out = outpath,
                reference = file.path(index_folder, "phiX.fa.gz"))
    expect_equal(counts, list(reads = 100, removed = 1))
    unlink(file.path(outpath, "*"), recursive=TRUE)
})

test_that("filtering works on full data sets", {
    data <- make_random_data()
    d <- tempdir()
    dir.create(file.path(d, "out"), showWarnings = FALSE)
    outpath <- file.path(d, "out")
    index_folder <- system.file("extdata/genomes", package = "mbtools")

    counts <- filter_reference(data, out = outpath,
                               reference = file.path(index_folder, "phiX.fa.gz"))
    expect_equal(counts[, sum(reads)], 200)
    expect_equal(counts[, sum(removed)], 1)
    expect_equal(counts[, uniqueN(id)], 2)
})
flog.threshold(INFO)
