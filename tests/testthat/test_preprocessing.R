context("preprocessing")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_preprocess()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_preprocess(maxEE = 3)
    expect_equal(config$maxEE, 3)
})

path <- system.file("extdata/16S", package = "mbtools")
files <- find_read_files(path)

test_that("preprocessing works", {
    out <- file.path(tempdir(), "proc")
    pre <- preprocess(files, out_dir = out)
    expect_length(list.files(out), 10)
    expect_equal(nrow(pre$files), 5)
    expect_named(pre, c("files", "passed", "steps"))
    expect_true(pre$passed[, all(preprocessed <= raw)])
})

flog.threshold(ERROR)
