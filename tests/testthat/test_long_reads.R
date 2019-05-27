context("long reads")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_align()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_preprocess(threads = 2)
    expect_equal(config$threads, 2)
})

path <- system.file("extdata/nanopore", package = "mbtools")
ref <- system.file("extdata/genomes/zymo_mock.fna.gz",
                   package = "mbtools")
files <- find_read_files(path)[1:2]

test_that("alignments work", {
    alns <- align_long_reads(files,
        alignment_dir = file.path(tempdir(), "aln"),
        reference = ref)

    expect_named(alns, c("alignments", "logs", "disk_size", "steps"))
    expect_s3_class(alns$alignments, "data.table")
    expect_true(alns$alignments[, all(success)])
    expect_true(alns$alignments[, all(file.exists(alignment))])
})

flog.threshold(INFO)
