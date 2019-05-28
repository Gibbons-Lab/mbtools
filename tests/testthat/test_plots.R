context("phyloseq plots")

flog.threshold(WARN)

otus <- sapply(seq(2, 4, length.out = 20),
    function(i) rpois(6, 10 ^ i))
otus <- otus * (1:6)
rownames(otus) <- paste("sample", 1:6)
colnames(otus) <- paste("taxa", 1:20)
tax <- matrix(rep(LETTERS[1:20], 2), ncol = 2)
colnames(tax) <- c("genus", "species")
rownames(tax) <- colnames(otus)
sdata <- data.frame(val = 1:6, sex = rep(c("M", "F"), 3))
rownames(sdata) <- rownames(otus)
ps <- phyloseq(otu_table(otus, taxa_are_rows = FALSE), tax_table(tax),
               sample_data(sdata))

test_that("we can plot counts", {
    pl <- plot_counts(ps, "val", tax_level = "species")
    dt <- plot_counts(ps, "val", tax_level = "species", only_data = TRUE)
    expect_s3_class(pl, "ggplot")
    expect_equal(nrow(dt), 120)
    expect_equal(dt[, uniqueN(taxa)], 20)
})

test_that("we can plot composition", {
    pl <- plot_taxa(ps, "genus", show_samples = TRUE)
    pl <- plot_taxa(ps, "genus", show_samples = FALSE)
    expect_s3_class(pl, "ggplot")
    dt <- plot_taxa(ps, "species", only_data = TRUE)
    expect_equal(nrow(dt), 6 * 12)
    expect_equal(dt[, uniqueN(taxa)], 12)
})
