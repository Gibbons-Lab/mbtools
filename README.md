<img src="https://github.com/Gibbons-Lab/mbtools/raw/master/inst/extdata/logo.png" width="25%" align="right">

[![Build Status](https://dev.azure.com/chdiener/cdiener/_apis/build/status/Gibbons-Lab.mbtools?branchName=master)](https://dev.azure.com/chdiener/cdiener/_build/latest?definitionId=1&branchName=master)
![Azure DevOps coverage](https://img.shields.io/azure-devops/coverage/chdiener/cdiener/1)
![Azure DevOps tests (compact)](https://img.shields.io/azure-devops/tests/chdiener/cdiener/1?compact_message)
[![Build Status](https://travis-ci.org/Gibbons-Lab/mbtools.svg?branch=master)](https://travis-ci.org/Gibbons-Lab/mbtools)
[![codecov](https://codecov.io/gh/Gibbons-Lab/mbtools/branch/master/graph/badge.svg)](https://codecov.io/gh/Gibbons-Lab/mbtools)



`mbtools` is a collection of helpers that we use to analyze microbiome
data. It makes it easier to run some common analyses and is pretty
opinionated towards our own experiences.

## General philosophy

`mbtools` sees analyses as a workflow consisting of several analysis steps
and final outputs based on those. This is pretty similar to Qiime 2 and most
of the functionality in `mbtools` is supposed to be orthogonal and not in
competition to those other excellent tools.

## Data types

`mbtools` mostly works on the following data types:

1. artifacts - compound data objects returned from an analysis step
2. phyloseq object - a [phyloseq](https://joey711.github.io/phyloseq/) object containing sequence variant abundances,
   taxonomy assignments and sample metadata
3. data tables - general data frame-like objects in [tidy format](https://r4ds.had.co.nz/tidy-data.html)

## Workflow steps

For `mbtools` a workflow step is based on input data and a configuration,
thus having the function signature `step(object, config)`.
Most steps can be chained with the pipe operator to yield workflows.
For instance, the following is a possible workflow with `mbtools`:

```r
library(mbtools)

config <- list(
    demultiplex = config_demultiplex(barcodes = c("ACGTA", "AGCTT")),
    preprocess = config_preprocess(truncLen = 200),
    denoise = config_denoise()
)

output <- find_read_files("raw") %>%
          demultiplex(config$demultiplex) %>%
          quality_control() %>%
          preprocess(config$preprocess) %>%
          denoise(config$denoise)
```

This clearly logs the used workflow and the configuration. The configuration
can also be saved and read in many formats, for instance yaml.

config.yaml:
```
preprocess:
  threads: yes
  out_dir: preprocessed
  trimLeft: 10.0
  truncLen: 200.0
  maxEE: 2.0
denoise:
  threads: yes
  nbases: 2.5e+08
  pool: no
  bootstrap_confidence: 0.5
  taxa_db: https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1
  species_db: https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1
  hash: yes
```

This can now be reused by someone else:

```r
config <- read_yaml("config.yml")

output <- find_read_files("raw") %>%
          quality_control() %>%
          preprocess(config$preprocess) %>%
          denoise(config$denoise)
```

## Other functions

All other functions are usually functions that are meant to be inside
more complex code or functions that produce plots and endpoints of
an analysis. Most of them act on phyloseq objects and some on tidy data
tables.
