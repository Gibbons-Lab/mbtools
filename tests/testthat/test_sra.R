context("SRA submissions")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_sra()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_preprocess(title = "blub")
    expect_equal(config$title, "blub")
})

out <- file.path(tempdir(), "sra")
fi <- system.file("extdata/shotgun", package = "mbtools") %>% find_read_files()

test_that("sra bundling works", {
    expect_error(sra_submission(fi))

    sra <- sra_submission(fi, out_dir = out, title = "Sequencing of X: XYZ",
                          preset = "human gut metagenome",
                          platform = "ILLUMINA",
                          instrument = "Illumina HiSeq 2000",
                          metadata = data.table(id = fi$id, lane = fi$lane,
                                                date = "2019-01-01"))

    expect_named(sra, c("files", "biosample_attributes", "sra_metadata",
                        "upload"))
    expect_true(file.exists(sra$upload))
})

flog.threshold(INFO)
