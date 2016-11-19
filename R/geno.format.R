# format a 4-column genotype matrix or data.frame into a simplified numerical format
# X is the genotype matrix or data.frame
# Y is a 4-column numeric matrix.  Each row is sorted with zeros indicating missing alleles.
# other optional output (Y.only = FALSE) is:
# A: conversion between allele names and allele numerical identifiers
# G: unique allele combinations, i.e. Y without duplicated rows
geno.format <- function(X, Y.only = TRUE, debug = FALSE) {
  # convert to numeric
  Y <- as.matrix(X)
  if(!is.numeric(Y)) {
    if(any(is.na(Y)) && any(Y[!is.na(Y)] == "")) {
      stop("Cannot supply both empty character and NA.")
    }
    Y[Y == ""] <- NA
    if(debug) browser()
    A <- factor(Y)
    Y <- matrix(as.numeric(A), ncol = 4)
    A <- levels(A)
  } else {
    A <- levels(factor(Y))
  }
  Y[is.na(Y)] <- 0
  # sort numeric
  Y <- t(apply(Y, 1, function(x) {
    x <- sort(x)
    c(x, rep(0, 4-length(x)))
  }))
  # unique allele combinations observed
  G <- unique(colSort(Y))
  if(Y.only) {
    ans <- Y
  } else {
    ans <- list(Y = Y, A = A, G = G)
  }
  ans
}
