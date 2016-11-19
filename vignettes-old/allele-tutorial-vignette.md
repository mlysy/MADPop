<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Modeling of Variable-Alleled Species}
-->

Hierarchial Hardy-Weinberg Analysis of Variable Alleled Fish
========================================================
TUTORIAL
--------------------------------------------------------
There are many steps in this tutorial that are labelled STEP 1 to STEP 6. Each step has a "PURPOSE" written below for explanation of the step.
Throughout each of the steps there will be a "NOTE" as guidance.


```r
require(mapsAllele)
```

#
#
**STEP 0.**

PURPOSE: Clean up dataset. For a preloaded dataset, skip to STEP 1.

See mapsAllele/extdata/fish215-cleanup.R

#
#
**STEP 1.**
PURPOSE: Calculate some summary statistics for the data, to (1) simulate datasets which satisfy the Hardy-Weinberg equilibrium, and (2) draw MCMC samples for posterior inference. To proceed directly to inference go to step 2.


```r
# fish data
head(fish215)
```

```
##    Lake  A1  A2  A3  A4
## 1 Hogan r.3     r.1 r.6
## 2 Hogan     r.3     r.1
## 3 Hogan r.4            
## 4 Hogan r.4     r.2 r.7
## 5 Hogan     r.3     r.5
## 6 Hogan     t.1 r.8 r.1
```

```r
# Here are two lakes under consideration: Simcoe and Michipicoten
lnames <- c("Simcoe", "Michipicoten")
Xobs <- fish215[fish215$Lake %in% lnames,]
```

Pool these under the null hypothesis:

H0: the two lakes have the same Hardy-Weinberg equilibrium conditions.


```r
Xsuff <- HW.suff(X = Xobs[,-1])
# converts the fish data Xobs into statistics including:
#   A: total set of unique alleles
#   G: total set of unique full genotypes (combination of 4 alleles)
#   Xg: counts in each unique genotype from observed data
#   H: total set of possible inherited genotypes (combination of 2 alleles)
#       with repeats being 1 inherited allele
#   Y: observed data in a simplified numerical format
```

This is the end of STEP 1. In STEP 1, we have sourced all of the functions we will be using, and loaded the model, as well as the fish215. In STEP 2, we will be working with the lake data to fit models using STAN.

#
#
**STEP 2.**

PURPOSE: to make inferences for two lakes under the null hypothesis. Null hypothesis is that the two lakes are the same, and they are not treated separately (grouped into one).































