context("replication rates")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_rep()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_rep(min_coverage = 3)
    expect_equal(config$min_coverage, 3)
})

sl <- ERR260132

test_that("replication rates work", {
    r <- replication_rates(sl)

    expect_named(r, c("rate", "coverage", "steps"))
    expect_gt(nrow(r$rate), 3)
    expect_s3_class(r$rate, "data.table")
    expect_true(r$rate[, all(rate > 0)])
    expect_s3_class(r$coverage, "data.table")
    expect_equal(r$coverage[, uniqueN(contig)], nrow(r$rate))
})

flog.threshold(INFO)
