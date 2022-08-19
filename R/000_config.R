# Helpers to manage workflow step configurations

config_builder <- function(default_parms) {
    proto <- function(...) {
        config <- default_parms
        args <- list(...)
        for (arg in names(args)) {
            config[[arg]] <- args[[arg]]
        }
        class(config) <- "config"
        return(config)
    }
    return(proto)
}

config_parser <- function(args, config_fun) {
    if (length(args) == 1 && inherits(args[[1]], "config")) {
        config <- args[[1]]
    } else {
        config <- do.call(config_fun, args)
    }
    return(config)
}
