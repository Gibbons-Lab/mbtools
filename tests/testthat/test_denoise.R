context("denoise")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_denoise()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_preprocess(trimLeft = 5)
    expect_equal(config$trimLeft, 5)
})

path <- system.file("extdata/16S", package = "mbtools")
files <- find_read_files(path)[1:2] %>%
         preprocess(out_dir = file.path(tempdir(), "pre"))

test_that("denoise works", {
    den <- denoise(files[1], truncLength = c(240, 180), species_db = NULL)

    expect_named(den, c('feature_table', 'taxonomy', 'errors', 'error_plots',
                        'passed_reads', 'classified', 'steps'))
    expect_s3_class(den$passed_reads, "data.table")
    expect_true(den$passed_reads[,
        all(merged <= pmin(derep_forward, derep_reverse))])
    expect_true(den$passed_reads[, all(non_chimera <= merged)])
    expect_true(den$classified[, all(reads > 0)])
    expect_true(den$classified[rank == "Genus", reads > 0.5])
})

flog.threshold(INFO)
