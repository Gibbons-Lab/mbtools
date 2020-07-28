context("SLIMM lineage calling")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_slimm()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_preprocess(database = "hello")
    expect_equal(config$database, "hello")
})

fi <- system.file("extdata/shotgun", package = "mbtools") %>%
      find_read_files()
ref <- system.file("extdata/genomes/zymo_mock.fna.gz",
                   package = "mbtools")
alns <- align_short_reads(fi, alignment_dir = file.path(tempdir(), "aln"),
                              reference = ref)
db <- system.file("extdata/ABVF_SP_CMP_genomes.sldb", package = "mbtools")

test_that("SLIMM works", {
    sl <- slimm(alns, database = db)

    expect_named(sl, c("alignments", "abundance", "coverage", "steps"))
    expect_s3_class(sl$abundance, "data.table")
    expect_s3_class(sl$coverage, "data.table")
    expect_true(sl$coverage[, sum(coverage[[1]]), by = "contig"][, all(V1 > 0)])
    expect_equal(sl$coverage[, uniqueN(contig)], 10)
})

flog.threshold(INFO)
