# MADPop: MHC Allele-Based Differencing Between Populations

*Martin Lysy, Peter W.J. Kim, Terin Robinson*

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/mlysy/MADPop/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mlysy/MADPop/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
---

### Description

Tools for the analysis of population differences
using the Major Histocompatibility Complex (MHC) genotypes of samples
having a variable number of alleles (1-4) recorded for each
individual.  A hierarchical Dirichlet-Multinomial model on the
genotype counts is used to pool small samples from multiple
populations for pairwise tests of equality.  Bayesian inference is
implemented via the [**rstan**](https://mc-stan.org/rstan/) package.  Bootstrapped and posterior
p-values are provided for chi-squared and likelihood ratio tests of
equal genotype probabilities.

### Installation

To install the latest R release:
```r
install.packages("MADPop")
```

To install the latest development version, first install the R package [**devtools**](https://CRAN.R-project.org/package=devtools) then run
```r
devtools::install_github("mlysy/MADPop")
```

### Usage

Please see tutorial `vignette("MADPop-quicktut")`.
