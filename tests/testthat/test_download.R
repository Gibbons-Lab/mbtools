context("download helpers")

flog.threshold(WARN)

test_that("we can download files", {
    fi <- data.table(
        url = "https://github.com/Gibbons-Lab/mbtools/raw/master/R/000_config.R",
        target = file.path(tempdir(), "test.R"),
        description = "test file"
    )
    dld <- download_files(fi)
    expect_equal(nrow(dld), 1)
    expect_equal(ncol(dld), 4)
})

test_that("we can get SRA download links", {
    sra <- sra_download_url("ERR260132")
    expect_type(sra, "character")
    expect_equal(length(sra), 2)
})

flog.threshold(INFO)
