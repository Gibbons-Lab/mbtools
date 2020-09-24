context("power analysis")

flog.threshold(WARN)

p <- rexp(10)
p <- p / sum(p)
counts <- sapply(p, function(pi) rbinom(8, 1e5, pi))
rownames(counts) <- paste0("sample_", 1:8)
colnames(counts) <- paste0("asv_", 1:10)
ps <- phyloseq(otu_table(counts, taxa_are_rows = FALSE),
               tax_table(matrix(colnames(counts), ncol = 1,
                                dimnames = list(colnames(counts), "taxa"))),
               sample_data(data.frame(x = rep(0, 8),
                                      row.names = rownames(counts))))


test_that("PERMANOVA power analysis works", {
    pa <- power_analysis(ps, n = c(4, 20), effect_size = c(0, 0.9))
    power <- pa$power
    expect_lt(power[effect != 0, max(fdr, na.rm = T)], 0.5)
    expect_gt(power[effect == max(effect), max(power)], 0.9)
    expect_gt(power[effect == max(effect) & n == 20, mean(power)], 0.5)
})

test_that("Mann-Whitney power analysis works", {
    pa <- power_analysis(ps, method = "mw", n = c(4, 20),
                         effect_size = c(0, 0.9))
    power <- pa$power
    expect_lt(power[effect != 0, max(fdr, na.rm = T)], 0.2)
    expect_gt(power[effect == max(effect), max(power)], 0.9)
    expect_gt(power[effect == max(effect) & n == 20, mean(power)], 0.5)
})

test_that("Corncob power analysis works", {
    pa <- power_analysis(ps, method = "corncob", n = c(4, 20),
                         effect_size = c(0, 0.9), n_groups = 3, n_power = 4)
    power <- pa$power
    expect_lt(power[effect != 0, max(fdr, na.rm = T)], 0.2)
    expect_gt(power[effect == max(effect), max(power)], 0.9)
    expect_gt(power[effect == max(effect) & n == 20, mean(power)], 0.5)
})

flog.threshold(INFO)
