#' Equality tests for two Hardy-Weinberg multinomial samples
#'
#' Generate multinomial samples from a common Hardy-Weinberg equilibrium distribution and calculate the Chi-square and Likelihood Ratio test statistics.
#'
#' @param N1,N2 Sample sizes.
#' @param H 2-column matrix of inheritable chromosomes.
#' @param rho Vector of probabilities for inheritance of each chromosome from either parent.  Can also be a matrix, in which case each simulation randomly draws with replacement from the rows of rho.
#' @param nreps Number of replications of the simulation.
#' @return An \code{nreps x 2} matrix with the simulated chi-squared and LR values.
#' @keywords internal
HW.eqtest <- function(N1, N2, H, rho, nreps, verbose = TRUE) {
  # precomputations
  HH <- as.matrix(expand.grid(1:nrow(H), 1:nrow(H)))
  HH <- cbind(H[HH[,1],], H[HH[,2],])
  HH <- t(apply(HH, 1, function(x) {
    y <- sort(unique(x))
    c(rep(0, 4-length(y)), y)
  }))
  HH <- apply(HH, 1, paste0, collapse = ".")
  lake.id <- c(rep(1, N1), rep(2, N2))
  N <- N1+N2
  if(is.matrix(rho)) {
    Nrho <- nrow(rho)
  } else {
    Nrho <- 1
    rho.sim <- c(rho %o% rho)
  }
  # output
  boot.out <- matrix(NA, nreps, 2)
  colnames(boot.out) <- c("chi2", "LRT")
  for(ii in 1:nreps) {
    if(verbose && (ii%%5e3 == 0)) message("bootstrap sample ", ii, "/", nreps)
    if(Nrho > 1) {
      jj <- sample(Nrho, 1)
      rho.sim <- c(rho[jj,] %o% rho[jj,])
    }
    tab.sim <- table(lake.id,
                     sample(HH, N, replace = TRUE, prob = rho.sim))
    boot.out[ii,1] <- chi2.stat(tab.sim)
    boot.out[ii,2] <- LRT.stat(tab.sim)
  }
  boot.out
}
