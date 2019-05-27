context("quality control")

flog.threshold(WARN)

path <- system.file("extdata/16S", package = "mbtools")
files <- find_read_files(path)

test_that("we can run quality controls", {
    expect_equal(nrow(files), 5)

    quals <- files %>% quality_control(min_score = 20)
    expect_named(quals, c("files", "qualities", "bases", "quality_plot",
                          "length_plot", "entropy_plot", "steps"))
    expect_s3_class(quals$qualities, "data.table")
})

flog.threshold(INFO)
