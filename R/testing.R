# Functions for testing associations between the microbiome
# and exogenous variables

check_variables <- function(variable, counts, meta, confounders) {
    good <- !is.na(meta[[variable]])
    # Check for constant variables
    sds <- sapply(c(variable, confounders),
                 function(co) sd(as.numeric(meta[good, co])))
    if (any(is.na(sds) | sds < 1e-6)) {
        flog.info("Detected variable or confounder with zero variation.")
        return(NULL)
    }
    if (is.null(confounders)) {
        ref_model <- ~ 1
    } else {
        ref_model <- reformulate(confounders)
        num_conf <- sapply(confounders, function(co) as.numeric(meta[[co]]))
        corrs <- cor(as.numeric(meta[[variable]]), num_conf,
                     use = "pairwise.complete.obs")
        if (any(abs(corrs) > 0.9)) {
            flog.info(
                paste("Variable `%s` is correlated with the confounders.",
                       "Skipping it."),
                variable)
            return(NULL)
        }
    }
    if (sum(good) < length(confounders) + 2) {
        flog.info("Not enough degrees of freedom. Skipping variable `%s`.",
                  variable)
        return(NULL)
    }
    return(list(good = good, ref_model = ref_model))
}

iter_deseq2 <- function(variable, counts, meta, confounders, shrink, tax) {
    is_reg <- !is.factor(meta[[variable]])
    check <- check_variables(variable, counts, meta, confounders)
    if (is.null(check)) {
        return(NULL)
    } else {
        good <- check$good
        ref_model <- check$ref_model
    }
    if (ncol(counts) < 50) {
        fit_type <- "mean"
    } else {
        fit_type <- "local"
    }
    dds <- suppressMessages(
        DESeqDataSetFromMatrix(t(counts[good, ]), meta[good, ],
        design = reformulate(c(confounders, variable))))
    dds <- suppressMessages(
        DESeq(dds, test = "LRT", parallel = FALSE, quiet = TRUE,
              fitType = fit_type, reduced = ref_model,
              sfType = "poscounts"))
    res <- results(dds)
    if (shrink) {
        res <- lfcShrink(dds, coef = length(resultsNames(dds)),
                            res = res, quiet = TRUE)
    }
    res <- as.data.table(res)
    set(res, j = ifelse(is.na(tax), "variant", tax),
        value = colnames(counts))
    set(res, j = "variable", value = variable)
    set(res, j = "coef", value = resultsNames(dds)[length(resultsNames(dds))])
    if (is_reg) {
        n <- sum(good)
    } else {
        levs <- levels(meta[[variable]])
        n <- min(sum(meta[good, variable] == levs[1]),
                 sum(meta[good, variable] == levs[length(levs)]))
    }
    set(res, j = "n_eff", value = n)
    return(res)
}

iter_voom <- function(variable, counts, meta, confounders, shrink, tax) {
    is_reg <- !is.factor(meta[[variable]])
    check <- check_variables(variable, counts, meta, confounders)
    if (is.null(check)) {
        return(NULL)
    } else {
        good <- check$good
        ref_model <- check$ref_model
    }
    norm_counts <- t(suppressMessages(normalize(counts[good, ])))
    design <- model.matrix(reformulate(c(confounders, variable)), data=meta)
    model <- voom(norm_counts, design, plot=FALSE)
    fit <- lmFit(model, design)

    if (shrink) {
        fit <- eBayes(fit)
    }
    res <- topTable(fit, coef=ncol(design), sort.by="none", number=Inf)
    res <- as.data.table(res)
    names(res) <- c("log2FoldChange", "baseMean", "t", "pvalue", "padj", "B")
    res$baseMean <- 2 ^ res$baseMean
    set(res, j = ifelse(is.na(tax), "variant", tax),
        value = colnames(counts))
    set(res, j = "variable", value = variable)
    if (is_reg) {
        n <- sum(good)
    } else {
        levs <- levels(meta[[variable]])
        n <- min(sum(meta[good, variable] == levs[1]),
                 sum(meta[good, variable] == levs[length(levs)]))
    }
    set(res, j = "n_eff", value = n)
    return(res)
}

#' Build a configuration for the alignment workflows.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the long read alignment
#'  workflow.
#' @export
#' @examples
#'  config <- config_association(variable = "glucose")
config_association <- config_builder(list(
    variables = NULL,
    taxa_rank = "genus",
    confounders = NULL,
    min_count = 10,
    in_samples = 0.1,
    independent_weighting = TRUE,
    standardize = TRUE,
    shrink = TRUE,
    method = "deseq2",
    threads = getOption("mc.cores", 1)
))

#' Run differential association tests between taxa counts and exogenous
#' factors.
#'
#' @param ps A phyloseq object containing the taxa counts.
#' @param ... A configuration as returned by \code{\link{config_association}}.
#' @return An artifact containing the results.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom data.table set
#' @importFrom phyloseq sample_data
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results lfcShrink
#'  resultsNames estimateSizeFactors
#' @importFrom limma voom eBayes lmFit topTable
#' @importFrom stats model.matrix
association <- function(ps, ...) {
    config <- config_parser(list(...), config_association)
    meta <- as(sample_data(ps), "data.frame")
    if (!is.null(config$confounders)) {
        missing_conf <- apply(sample_data(ps)[, config$confounders], 1,
                              function(row) any(is.na(row)))
        ps <- prune_samples(!missing_conf, ps)
    }
    if ((!requireNamespace("IHW", quietly = TRUE)) &&
          config$independent_weighting) {
        stop("independent weighting requires the IHW package!")
    }
    if (config$standardize) meta <- standardize(meta)
    if (is.null(config$variables)) {
        variables <- names(meta)
    } else {
        variables <- config$variables
    }
    variables <- variables[!(variables %in% config$confounders)]
    counts <- as.matrix(taxa_count(ps, lev = config$taxa_rank))
    meta <- meta[rownames(counts), , drop = FALSE]
    too_rare <- (colSums(counts >= 1) / nrow(counts)) < config$in_samples
    too_few <- colMeans(counts) < config$min_count
    counts <- counts[, !(too_rare | too_few)]
    if (config$method == "deseq2") {
        iter <- iter_deseq2
    } else if (config$method == "voom") {
        iter <- iter_voom
    } else {
        stop("`%s` is not a recognized method :/", config$method)
    }

    apfun <- parse_threads(config$threads)
    report_at <- max(ceiling(length(variables) / 100),
                     ceiling(1e4 / nrow(counts)))
    tests <- apfun(1:length(variables), function(i) {
        v <- variables[i]
        test <- iter(v, counts = counts, meta = meta,
                     confounders = config$confounders,
                     shrink = config$shrink, tax = config$taxa_rank)
        if (i %% report_at == 0) {
            flog.info("Finished running %d/%d tests (%d%%).",
                      i * ncol(counts), length(variables) * ncol(counts),
                      floor(100 * i / length(variables)))
        }
        return(test)
    }) %>% rbindlist()

    if (length(variables) > 1) {
        if (config$independent_weighting) {
            weights <- IHW::ihw(pvalue ~ baseMean, tests, alpha = 0.05)
            set(tests, j = "padj", value = IHW::adj_pvalues(weights))
        } else {
            set(tests, j = "padj", value =
                p.adjust(tests$pvalue, method = "BH"))
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
    if ((!requireNamespace("IHW", quietly = TRUE)) && independent_weighting) {
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
                 fitType = "local", reduced = ref_model)
    combinations <- combn(length(levs):1, 2)
    tests <- pbapply(combinations, 2, function(co) {
        name <- paste0(variable, levs[co[1]], "_vs_", levs[co[2]])
        res <- results(dds, contrast = c(variable, levs[co[1]], levs[co[2]]))
        if (shrink) {
            res <- lfcShrink(dds,
                             contrast = c(variable, levs[co[1]], levs[co[2]]),
                             res = res, type = "normal")
        }
        res <- as.data.table(res)
        set(res, j = tax, value = colnames(counts))
        set(res, j = "variable", value = name)
        n <- min(table(meta[[variable]])[co])
        set(res, j = "n_eff", value = n)
        return(res)
    })
    tests <- rbindlist(tests)

    if (independent_weighting) {
        weights <- IHW::ihw(pvalue ~ baseMean, tests, alpha = 0.05)
        set(tests, j = "padj", value = IHW::adj_pvalues(weights))
    } else {
        set(tests, j = "padj", value = p.adjust(tests$pvalue, method = "BH"))
    }

    return(tests)
}
