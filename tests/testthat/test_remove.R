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

    dset <- data.table(forward = file.path(d, c("f.fastq.gz", "f2.fastq.gz")),
                       reverse = file.path(d, c("r.fastq.gz", "r2.fastq.gz")),
                       id = c("random", "single"))
    return(dset)
}

test_that("sequences can be removed", {
    data <- make_random_data()
    d <- tempdir()
    outpath <- file.path(d, "out")
    unlink(file.path(outpath), recursive = TRUE)
    dir.create(outpath, showWarnings = FALSE)

    index_folder <- system.file("extdata/genomes", package = "mbtools")
    reads <- data[1]
    conf <- config_reference(
        out_dir = outpath,
        reference = file.path(index_folder, "phiX.fa.gz")
    )
    counts <- filter_reference(reads, conf)
    expect_equal(counts$counts, data.table(reads = 100, removed = 0,
                                           id = "random", lane = NA))
    unlink(file.path(outpath, "*"), recursive = TRUE)

    counts <- filter_reference(reads[, .(forward, id)], conf)
    expect_equal(counts$counts, data.table(reads = 100, removed = 0,
                                           id = "random", lane = NA))
    unlink(file.path(outpath, "*"), recursive = TRUE)

    reads <- data[2]
    counts <- filter_reference(reads, conf)
    expect_equal(counts$counts, data.table(reads = 100, removed = 1,
                                           id = "single", lane = NA))
    unlink(file.path(outpath, "*"), recursive = TRUE)
})

test_that("filtering works on full data sets", {
    data <- make_random_data()
    d <- tempdir()
    dir.create(file.path(d, "out"), showWarnings = FALSE)
    outpath <- file.path(d, "out")
    index_folder <- system.file("extdata/genomes", package = "mbtools")
    conf <- config_reference(
        out_dir = outpath,
        reference = file.path(index_folder, "phiX.fa.gz")
    )
    counts <- filter_reference(data, conf)
    expect_equal(counts$counts[, sum(reads)], 200)
    expect_equal(counts$counts[, sum(removed)], 1)
    expect_equal(counts$counts[, uniqueN(id)], 2)
})
flog.threshold(INFO)
