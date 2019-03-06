[![Build Status](https://travis-ci.org/Gibbons-Lab/mbtools.svg?branch=master)](https://travis-ci.org/Gibbons-Lab/mbtools)

# :poop: :leaves: :earth_americas: mbtools

`mbtools` is a collection of helpers that we use to analyze microbiome
data. It makes it easier to run some common analyses and is pretty opinionated.

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

For `mbtools` a workflow step is based
