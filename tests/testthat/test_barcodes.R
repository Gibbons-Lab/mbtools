context("barcodes")

test_that("barcode splitting/checking work", {
    # download an example data set
    dir <- tempdir()
    mock <- mockrobiota("mock-3", dir, quiet = TRUE)
    reads <- c(mock$forward, mock$reverse)

    bc <- split_barcodes(reads, mock$index, file.path(dir, "filtered"),
        as.character(mock$samples$BarcodeSequence))
    expect_true(55979 == bc[1])
    expect_true(0 == bc[2])
    expect_true(0 == bc[3])
    expect_error(split_barcodes(reads, index, file.path(dir, "filtered"),
        "ACGAT"))
})
