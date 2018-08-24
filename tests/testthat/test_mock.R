context("mock diagnostics")

# A small mock table
taxa <- data.frame(
    kingdom = sample(letters[1:3], 100, replace = TRUE),
    genus = sample(letters[4:10], 100, replace = TRUE),
    species = sample(letters[11:26], 100, replace = TRUE)
)

otus <- replicate(4, runif(100))
colnames(otus) <- paste0("S", 1:4)

ps <- phyloseq(tax_table(as.matrix(taxa)),
               otu_table(otus, taxa_are_rows = TRUE))

test_that("taxa can be found and counted", {
    fi <- taxa_metrics(tax_table(ps), tax_table(ps))
    expect_equivalent(fi$precision, rep(1, 3))
    expect_equivalent(fi$recall, rep(1, 3))
    expect_equivalent(fi$F1, rep(1, 3))
    expect_equal(nrow(fi), 3)

    fi <- taxa_metrics(tax_table(ps), tax_table(as.matrix(taxa[-1, ])))
    expect_true(all(fi$precision <= 1))
    expect_equal(nrow(fi), 3)
})

test_that("taxa quantification can be calculated", {
    q <- taxa_quants(ps, ps)
    expect_equal(colnames(q), c("level", "name", "sample", "measured",
                                "reference"))
    expect_true(all(q$level %in% names(taxa)[-4]))
    expect_equal(q$measured, q$reference)

    p <- mock_plot(ps, ps)
    expect_true("ggplot" %in% class(p))
})
