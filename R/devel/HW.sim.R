#' Simulate data from the Hardy-Weinberg model
#'
#' Generate genotype data from a 4-allele Hardy-Weinberg model.
#'
#' @param N Number of samples to draw.
#' @param H 2-column matrix of possible chromosomes.
#' @param rho vector of the same length as \code{H} giving the probability of inheriting each chromosome in H from either parent.
#' @return An \code{N x 4} matrix giving the unique alleles for each observation.  0's indicate the presence of less than 4 alleles or duplicates.
#' @keywords internal
HW.sim <- function(N, H, rho) {
  nH <- nrow(H)
  # full genotype (i.e. which location has what and see duplicates)
  Xfull <- cbind(L1 = sample(nH, size = N, replace = TRUE, prob = rho),
                 L2 = sample(nH, size = N, replace = TRUE, prob = rho))
  # observed genotype (i.e. remove duplicates and location information)
  Xobs <- t(apply(Xfull, 1, function(x) {
    y <- unique(sort(H[x,]))
    c(rep(0, 4-length(y)), y)
  }))
  colnames(Xobs) <- paste0("A", 1:4)
  Xobs
}
