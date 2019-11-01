context("plate layouts")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_layout()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_layout(title = "blub")
    expect_equal(config$title, "blub")
})

manifest <- data.table(
    id = paste0("S", 1:640),
    age = runif(640, 18, 97) %>% ceiling(),
    treatment = c("control", "inulin")[(runif(640) > 0.5) + 1],
    sex = c("F", "M")[(runif(640) > 0.5) + 1]
)

test_that("layout works", {
    lo <- layout(manifest, by = "row", blank_step = 14)
    expect_equal(length(lo), 3)
    expect_true(lo$manifest[, any(layout_type == "blank")])
    expect_true("well" %in% names(lo$manifest))
})
flog.threshold(INFO)
