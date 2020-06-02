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
    expect_true(55979 == sum(bc$read_counts))
    for (i in 1:4) {
        expect_gt(bc$read_counts[paste0("S", i)], 5000)
    }
    expect_true(0 == bc$read_counts["unmatched"])
    expect_true(0 == bc$read_counts["multiple"])
    conf$barcodes <- "ACGAT"
    expect_error(demultiplex(files, conf))
})

flog.threshold(INFO)
