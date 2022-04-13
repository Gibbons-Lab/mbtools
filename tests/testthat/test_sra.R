context("SRA submissions")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_sra()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_sra(title = "blub")
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
    expect_true(file.exists(sra$upload[[1]]))
})

test_that("large sra submissions work", {
    s <- paste0("S", 1:3000)
    fi <- data.table(
        id = s,
        forward = sprintf("%s_1.fastq.gz", s),
        reverse = sprintf("%s_2.fastq.gz", s),
        lane = 1:3000
    )

    sra <- sra_submission(fi, out_dir = out, title = "Sequencing of X: XYZ",
                          preset = "human gut 16S",
                          platform = "ILLUMINA",
                          instrument = "Illumina MiSeq",
                          metadata = data.table(id = fi$id, lane = fi$lane,
                                                date = "2019-01-01"),
                          make_package = FALSE)
    expect_named(sra, c("files", "biosample_attributes", "sra_metadata",
                        "upload"))
    expect_equal(nrow(sra$biosample_attributes), 3000)
    expect_equal(nrow(sra$sra_metadata), 3000)
    expect_equal(length(sra$upload), 4)
    expect_equal(sum(sapply(sra$upload, length)), 2 * 3000)
})

flog.threshold(INFO)
