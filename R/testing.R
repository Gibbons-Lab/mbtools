# Functions for testing associations between the microbiome
# and exogenous variables

#' Run differential association tests between taxa counts and exogenous
#' factors.
#'
#' @param ps A phyloseq object containing the taxa counts.
#' @param variables Names of exogenous variables to include. Defaults to all
#'  variables.
#' @param tax The taxa level on which to run differential tests. Defaults to
#'  genus.
#' @param confounders A character vector containing the confounders that should
#'  be used.
#' @param min_count Minimum required number of average counts for a taxa.
#' @param in_samples Taxa must be present in at least this fraction of samples.
#' @param independent_weighting Whether to adjust p values by independent
#'  weighting or normal Benjamini-Hochberg.
#'  factors.
#' @param standardize Whether to standardize continuous variables to a mean
#'  of zero and a variance of 1. If True log fold changes for those variables
#'  denote are relative to a change of one standard deviation in the variable
#'  value.
#' @param shrink Whether to return shrunken log fold changes. Defaults to true.
#' @return A data.table containing the results.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom data.table set
#' @importFrom phyloseq sample_data
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink
#'  resultsNames estimateSizeFactors
association <- function(ps, variables = NULL, tax = "genus",
                        confounders = NULL, min_count = 10, in_samples = 0.1,
                        independent_weighting = TRUE, standardize = TRUE,
                        shrink = TRUE) {
    if (!is.null(confounders)) {
        missing_conf <- apply(sample_data(ps)[, confounders], 1,
                              function(row) any(is.na(row)))
        ps <- prune_samples(!missing_conf, ps)
    }
    if (!requireNamespace("IHW", quietly = TRUE)) {
        stop("independent weighting requires the IHW package!")
    }
    meta <- as(sample_data(ps), "data.frame")
    if (standardize) meta <- standardize(meta)
    if (is.null(variables)) variables <- names(meta)
    variables <- variables[!(variables %in% confounders)]
    counts <- as.matrix(taxa_count(ps, lev = tax))
    meta <- meta[rownames(counts), ]
    too_rare <- (colSums(counts >= 1) / nrow(counts)) < in_samples
    too_few <- colMeans(counts) < min_count
    counts <- counts[, !(too_rare | too_few)]

    log_file <- tempfile("mbtools", fileext = ".log")
    cat(paste("Writing logs to", log_file))
    sink(file(log_file, open = "wt"), type = "message")
    on.exit(sink(type = "message"))
    tests <- pblapply(variables, function(v) {
        good <- !is.na(meta[[v]])
        is_reg = !is.factor(meta[[v]])
        if (is.null(confounders)) {
            ref_model <- ~ 1
        } else {
            ref_model <- reformulate(confounders)
        }
        if (ncol(counts) < 50) {
            fit_type <- "mean"
        } else {
            fit_type <- "local"
        }
        dds <- DESeqDataSetFromMatrix(t(counts[good, ]), meta[good, ],
            design = reformulate(c(confounders, v)))
        dds <- estimateSizeFactors(dds, type = "poscount")
        dds <- DESeq(dds, test = "LRT", parallel = TRUE, quiet = TRUE,
                     fitType = fit_type, reduced = ref_model)
        if (is_reg) {
            totals <- colSums(counts(dds))
            mod <- glm(reformulate(c(confounders, v), "totals"),
                       data = meta[good, ])
            infl <- influence(mod)
            expected <- sum(infl$hat) / sum(good)
            outliers <- any(infl$hat > (expected * 10))
        } else {
            outliers <- FALSE
        }
        res <- results(dds)
        if (shrink) {
            res <- lfcShrink(dds, coef = length(resultsNames(dds)),
                             results = res)
        }
        res <- as.data.table(res)
        set(res, j = tax, value = colnames(counts))
        set(res, j = "variable", value = v)
        set(res, j = "robust", value = !outliers)
        if (is_reg) {
            n <- sum(good)
        } else {
            levs <- levels(meta[[v]])
            n <- min(sum(meta[[v]] == levs[1]),
                     sum(meta[[v]] == levs[length(levs)]))
        }
        set(res, j = "n_eff", value = n)
        return(res)
    })
    tests <- rbindlist(tests)

    if (length(variables) > 1) {
        if (independent_weighting) {
            weights <- IHW::ihw(pvalue ~ baseMean, tests, alpha = 0.05)
            set(tests, j = "padj", value = IHW::adj_pvalues(weights))
        } else {
            set(tests, j = "padj", value =
                p.adjust(tests$pvalues, method = "BH"))
        }
    }

    return(tests)
}


#' Run differential association tests between between all combinations of a
#' factor variable. Can be used as post-hoc test for regression.
#'
#' @param ps A phyloseq object containing the taxa counts.
#' @param variable The factor variable to permute.
#' @param tax The taxa level on which to run differential tests. Defaults to
#'  genus.
#' @param confounders A character vector containing the confounders that should
#'  be used.
#' @param min_count Minimum required number of average counts for a taxa.
#' @param in_samples Taxa must be present in at least this fraction of samples.
#' @param independent_weighting Whether to adjust p values by independent
#'  weighting or normal Benjamini-Hochberg.
#'  factors.
#' @param standardize Whether to standardize continuous variables to a mean
#'  of zero and a variance of 1. If True log fold changes for those variables
#'  denote are relative to a change of one standard deviation in the variable
#'  value.
#' @param shrink Whether to return shrunken log fold changes. Defaults to true.
#' @return A data.table containing the results.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom data.table set
#' @importFrom stats glm p.adjust influence reformulate
#' @importFrom utils combn
#' @importFrom phyloseq sample_data
combinatorial_association <- function(ps, variable, tax = "genus",
                        confounders = NULL, min_count = 10, in_samples = 0.1,
                        independent_weighting = TRUE, standardize = TRUE,
                        shrink = TRUE) {
    if (!is.null(confounders)) {
        missing_conf <- apply(sample_data(ps)[, confounders], 1,
                              function(row) any(is.na(row)))
        ps <- prune_samples(!missing_conf, ps)
    }
    if (!requireNamespace("IHW", quietly = TRUE)) {
        stop("independent weighting requires the IHW package!")
    }
    meta <- as(sample_data(ps), "data.frame")
    if (standardize) meta <- standardize(meta)
    if (!is.factor(meta[[variable]])) {
        stop("variable must be a factor.")
    }
    levs <- levels(meta[[variable]])
    if (length(levs) <= 2) {
        stop("variable must have more than 2 levels.")
    }
    counts <- as.matrix(taxa_count(ps, lev = tax))
    meta <- meta[rownames(counts), ]
    too_rare <- (colSums(counts >= 1) / nrow(counts)) < in_samples
    too_few <- colMeans(counts) < min_count
    counts <- counts[, !(too_rare | too_few)]

    good <- !is.na(meta[[variable]])
    if (is.null(confounders)) {
            ref_model <- ~ 1
        } else {
            ref_model <- reformulate(confounders)
        }
    dds <- DESeqDataSetFromMatrix(t(counts[good, ]), meta[good, ],
        design = reformulate(c(confounders, variable)))
    dds <- estimateSizeFactors(dds, type = "poscount")
    dds <- DESeq(dds, test = "LRT", parallel = TRUE, quiet = TRUE,
                 fitType = "local", reduce = ref_model)
    combinations <- combn(length(levs):1, 2)
    tests <- pbapply(combinations, 2, function(co) {
        name <- paste0(variable, levs[co[1]], "_vs_", levs[co[2]])
        res <- results(dds, contrast = c(variable, levs[co[1]], levs[co[2]]))
        if (shrink) {
            res <- lfcShrink(dds,
                             contrast = c(variable, levs[co[1]], levs[co[2]]),
                             results = res)
        }
        res <- as.data.table(res)
        set(res, j = tax, value = colnames(counts))
        set(res, j = "variable", value = name)
        n <- min(table(meta[[variable]])[co])
        set(res, j = "n_eff", value = n)
    })
    tests <- rbindlist(tests)

    if (independent_weighting) {
        if (!requireNamespace("IHW", quietly = TRUE)) {
            stop("independent weighting requires the IHW package!")
        }
        weights <- IHW::ihw(pvalue ~ baseMean, tests, alpha = 0.05)
        set(tests, j = "padj", value = IHW::adj_pvalues(weights))
    } else {
        set(tests, j = "padj", value = p.adjust(tests$pvalues, method = "BH"))
    }

    return(tests)
}

