context("contamination removal")
flog.threshold(WARN)

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

    index_folder <- system.file("extdata/genomes", package = "mbtools")
    counts <- remove_reference(reads, out = file.path(d, "filtered"),
                reference = file.path(index_folder, "phiX.fa.gz"),
                index = index)
    expect_equal(counts, c(reads = 100, removed = 0))

    counts <- remove_reference(reads[1], out = file.path(d, "filtered"),
                reference = file.path(index_folder, "phiX.fa.gz"))
    expect_equal(counts, c(reads = 100, removed = 0))

    phix <- readFasta(file.path(index_folder, "phiX.fa.gz"))
    phix <- substr(as.character(sread(phix)[1]), 1, 100)
    sr <- ShortReadQ(sread = DNAStringSet(c(seqs[1:99], phix)),
                     quality = BStringSet(quals),
                     id = BStringSet(paste0("S", 1:100)))
    writeFastq(sr, file.path(d, "f2.fastq.gz"))
    counts <- remove_reference(file.path(d, "f2.fastq.gz"), out = file.path(d, "filtered"),
                reference = file.path(index_folder, "phiX.fa.gz"))
    expect_equal(counts, c(reads = 100, removed = 1))
})

flog.threshold(INFO)
