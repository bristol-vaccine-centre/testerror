# testerror: Simulations and analytical approaches to test uncertainty for prevalence calculations

<!-- badges: start -->
[![R-CMD-check](https://github.com/bristol-vaccine-centre/testerror/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bristol-vaccine-centre/testerror/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/608545107.svg)](https://zenodo.org/badge/latestdoi/608545107)
[![testerror status
badge](https://bristol-vaccine-centre.r-universe.dev/badges/testerror)](https://bristol-vaccine-centre.r-universe.dev)
<!-- badges: end -->

The goal of `testerror` is to help produce modelled prevalence estimates for 
diseases based on positive test results in the inevitable presence of test error.
This is particularly problematic when test results are based on the combination 
of many test results from a multiplex panel, in which case error in component
tests compound to cause a great deal of uncertainty and in some situations 
significant bias.

## Installation

The `testerror` package is hosted
in the [Bristol Vaccine Centre
r-universe](https://https://bristol-vaccine-centre.r-universe.dev/).
Installation from there is as follows:

``` r
options(repos = c(
  "bristol-vaccine-centre" = 'https://https://bristol-vaccine-centre.r-universe.dev/',
  CRAN = 'https://cloud.r-project.org'))

# Download and install tableone in R
install.packages("testerror")
```

You can install the development version of `testerror` from
[GitHub](https://github.com/bristol-vaccine-centre/testerror) with:

``` r
# install.packages("devtools")
devtools::install_github("bristol-vaccine-centre/testerror")
```

## Further information

For examples and usage head to
the [main documentation
website](https://bristol-vaccine-centre.github.io/testerror/). An academic paper 
is in preparation which details the methodology implemented here.