library(testthat)
library(mbtools)

if (requireNamespace("xml2")) {
    test_check(
        "mbtools",
        reporter = MultiReporter$new(
            reporters = list(
                JunitReporter$new(file = "test-results.xml"),
                CheckReporter$new()))
    )
} else {
    test_check("mbtools")
}
