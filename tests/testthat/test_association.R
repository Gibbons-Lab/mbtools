context("association tests")

flog.threshold(WARN)

otus <- sapply(seq(2, 4, length.out = 20),
    function(i) rpois(6, 10 ^ i))
otus <- otus * (1:6)
rownames(otus) <- paste("sample", 1:6)
colnames(otus) <- paste("taxa", 1:20)
tax <- matrix(rep(LETTERS[1:20]), ncol = 1)
colnames(tax) <- "genus"
rownames(tax) <- colnames(otus)
sdata <- data.frame(val = 1:6, sex = rep(c("M", "F"), 3))
rownames(sdata) <- rownames(otus)
ps <- phyloseq(otu_table(otus, taxa_are_rows = FALSE), tax_table(tax),
               sample_data(sdata))

test_that("DESeq2 works", {
    tests <- association(ps)
    expect_equal(nrow(tests), 40)
    expect_true(tests[, all(padj > 0.01)])

    tests <- association(ps, confounders = "sex")
    expect_equal(nrow(tests), 20)
    expect_true(tests[, all(padj > 0.01)])
})

test_that("limma works", {
    tests <- association(ps, method = "voom")
    expect_equal(nrow(tests), 40)
    expect_true(tests[, all(padj > 0.01)])

    tests <- association(ps, method = "voom",
                         confounders = "sex")
    expect_equal(nrow(tests), 20)
    expect_true(tests[, all(padj > 0.01)])
})

test_that("combinatorial works", {
    sample_data(ps)$cat <- factor(rep(LETTERS[1:3], 2))
    tests <- combinatorial_association(ps, variable = "cat")
    expect_equal(nrow(tests), 60)
    expect_true(tests[, all(padj > 0.01)])
})

flog.threshold(INFO)
