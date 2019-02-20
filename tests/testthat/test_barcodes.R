context("barcodes")

flog.threshold(ERROR)

test_that("barcode splitting/checking work", {
    # download an example data set
    dir <- tempdir()
    mock <- mockrobiota("mock-3", dir, quiet = TRUE)
    files <- data.table(
        forward = mock$forward,
        reverse = mock$reverse,
        index = mock$index,
        id = "mock3")
    conf <- config_demultiplex(
        barcodes = mock$samples$BarcodeSequence,
        out_dir = file.path(dir, "demultiplexed")
    )
    bc <- demultiplex(files, conf)
    expect_true(55979 == bc$matches[1])
    expect_true(0 == bc$matches[2])
    expect_true(0 == bc$matches[3])
    conf$barcodes <- "ACGAT"
    expect_error(split_barcodes(files, conf))
})

flog.threshold(INFO)
