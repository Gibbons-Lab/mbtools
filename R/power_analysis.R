# Runs Power Analysis for microbial abundances

config_power <- config_builder(list(
    fraction_differential = 0.5,
    method = "permanova",
    type = "categorical",
    depth = "auto",
    n = ceiling(2 ^ seq(2, 7, length.out = 16) / 2) * 2,
    effect_size = seq(0, 0.9, by = 0.1),
    threads = getOption("mc.cores", 1),
    pval = 0.05,
    n_power = 100,
    n_groups = 8
))

memory_use <- function(config, n_taxa) {
    size <- 8 * max(config$n) * config$n_power * config$n_groups *
            n_taxa * config$threads
    class(size) <- "object_size"
    return(size)
}

get_corncob_pars <- function(ps, threads) {
    if (!requireNamespace("corncob", quietly = TRUE)) {
        stop("Power Analysis requires corncob to be installed.")
    }

    clean <- corncob::clean_taxa_names(ps)
    tnames <- attr(clean, "original_names")
    names(tnames) <- taxa_names(clean)

    apfun <- parse_threads(threads)
    pars <- apfun(taxa_names(clean), function(taxon){
        r <- tryCatch(corncob::bbdml(
            formula = reformulate("1", taxon),
            phi.formula = ~ 1,
            data = clean),
            error = function(e) list(mu.resp = NA, phi.resp = NA))
        cn <- otu_table(clean)[, taxon]
        data.table(mu = r$mu.resp[1], phi = r$phi.resp[1],
                   mean_reads = mean(cn), prevalence = mean(cn > 0),
                   taxon = taxon)
    }) %>% rbindlist()
    pars[, taxon := tnames[taxon]]
    setkey(pars, "taxon")

    return(pars)
}


sample_corncob <- function(pars, sig_taxa, type, scale,
                           size = 10000, reps = 1000) {
    n <- length(scale)
    last <- pars[, unique(taxon)[uniqueN(taxon)]]
    p <- sapply(pars$taxon, function(taxon) {
        if (taxon %in% sig_taxa) {
            p <- rep(pars[taxon, mu] * scale, reps)
        } else {
            p <- rep(pars[taxon, mu], length(scale) * reps)
        }
    })
    colnames(p) <- pars$taxon
    p[, last] <- pars[, sum(mu)] - rowSums(p[, 1:(ncol(p) - 1)])
    p <- sapply(1:ncol(p), function(i){
        taxon <- colnames(p)[i]
        VGAM::rbetabinom(n = n * reps, size = size, prob = p[, i],
                         rho = pars[taxon, phi])

    })
    colnames(p) <- pars$taxon

    snames <- paste0("sample_", 1:n)
    p <- lapply(1:reps, function(i) {
        x <- p[(1 + (i - 1) * n) : (i * n), ]
        rownames(x) <- snames
        return(x)
    })

    return(p)
}

mwtest <- function(counts, sig_taxa) {
    first <- 1:(nrow(counts) / 2)
    second <- (nrow(counts) / 2) : nrow(counts)
    p <- counts / rowSums(counts)
    res <- lapply(sig_taxa, function(ta){
        ctrl <- p[first, ta] + 0.5
        case <- p[second, ta] + 0.5
        te <- suppressWarnings(
            wilcox.test(ctrl, case))
        data.table(taxa = ta, pval = te$p.value)
    }) %>% rbindlist()
    res[is.na(pval), "pval" := 1]
    return(res)
}

#' Run a power analysis for microbe abundances from a reference sample.
#'
#' @param ps A phyloseq object containing data for a reference experiment.
#'  This should assume absence of any differential effect (all variation is
#'  random).
#' @param ... A configuration as described by \code{\link{config_power}}.
#' @return An artifact with an entry `power` that quantifies the power
#'  over many n and effect size.
#' @export
power_analysis <- function(ps, ...) {
    config <- config_parser(list(...), config_power)
    apfun <- parse_threads(config$threads)
    if (config$depth == "auto") {
        config$depth <- sample_sums(ps) %>% median() %>% ceiling()
        flog.info("Using median depth from reference data (%d).",
                  config$depth)
    }

    flog.info("Estimating corncob model parameters for %d taxa...", ntaxa(ps))
    pars <- get_corncob_pars(ps, config$threads)
    ps <- prune_taxa(pars[!is.na(mu), taxon], ps)
    flog.info("Succesfully estimated parameters for %d/%d taxa.",
              ntaxa(ps), nrow(pars))
    pars <- pars[!is.na(mu)]
    sig_taxa <- sample(taxa_names(ps)[1:(ntaxa(ps) - 1)],
                       config$fraction_differential * ntaxa(ps) - 1)
    sig_taxa <- c(sig_taxa, taxa_names(ps)[ntaxa(ps)])
    comb <- expand.grid(list(n = config$n,
                             effect_size = config$effect_size))
    flog.info(paste("Estimating power for %d n/effect combinations.",
                    "%d/%d taxa are truly differential."), nrow(comb),
                    length(sig_taxa), ntaxa(ps))
    flog.info("Will need at least %s of memory.",
              memory_use(config, ntaxa(ps)) %>% format(unit = "auto"))
    if (config$method == "permanova") {
        power <- apfun(1:nrow(comb), function(i) {
            co <- as.numeric(comb[i, ])
            if (config$type == "categorical") {
                scale <- c(rep(1, co[1] / 2), rep(1 - co[2], co[1] / 2))
                v <- c(rep("control", co[1] / 2), rep("changed", co[1] / 2))
            } else {
                v <- seq(0, 1, length.out = co[1])
                scale <- 1 - v * co[2]

            }
            counts <- sample_corncob(pars, sig_taxa, config$type, scale,
                                     config$depth,
                                     config$n_power * config$n_groups)
            p <- lapply(counts, function(co) {
                res <- suppressMessages(
                    vegan::adonis2(co ~ v, data = data.frame(v = v)))
                data.table(r2 = res[1, "R2"], pval = res[1, "Pr(>F)"])
            }) %>% rbindlist()
            p[, "replicate" := rep(1:config$n_groups, config$n_power)]
            p <- p[, .(n = co[1], effect = co[2], r2 = mean(r2),
                       power = mean(pval < config$pval)), by = "replicate"]
            p <- p[, .(n = co[1], effect = co[2], r2 = mean(r2),
                       r2_sd = sd(r2),
                       power = mean(power), power_sd = sd(power))]
            flog.info(paste("Power for n=%d and effect=%.3g is %.3g ",
                            "(reported R2=%.3g)."),
                      as.integer(co[1]), co[2], p$power, p$r2)
            return(p)
        })
    } else if (config$type == "categorical") {
        power <- apfun(1:nrow(comb), function(i) {
            co <- as.numeric(comb[i, ])
            scale <- c(rep(1, co[1] / 2), rep(1 - co[2], co[1] / 2))
            counts <- sample_corncob(pars, sig_taxa, config$type, scale,
                                     config$depth,
                                     config$n_groups * config$n_power)
            p <- lapply(counts, function(co) {
                mwtest(co, sig_taxa)
            }) %>% rbindlist()
            p[, "replicate" := rep(1:config$n_groups, config$n_power),
              by = "taxa"]
            p <- p[, .(taxa, power = mean(pval < config$pval)),
                       by = c("replicate", "taxa")]
            p <- p[, .(n = co[1], effect = co[2],
                       mu = pars[taxa[1], mu], phi = pars[taxa[1], phi],
                       mean_reads = pars[taxa[1], mean_reads],
                       prevalence = pars[taxa[1], prevalence],
                       power = mean(power), power_sd = sd(power)),
                       by = "taxa"]
            flog.info("Power for n=%d and effect=%.3g is %.3g [%.3g, %.3g].",
                      as.integer(co[1]), co[2], mean(p$power),
                      quantile(p$power, 0.025), quantile(p$power, 0.975))
            return(p)
        })
    }

    power <- rbindlist(power)
    if (config$method == "permanova") {
        power[, "asym_r2" := r2[n == max(n)], by = "effect"]
        power[, "asym_r2_sd" := r2_sd[n == max(n)], by = "effect"]
    }

    artifact <- list(
        power = power,
        steps <- c("power_analysis")
    )
    return(artifact)
}
