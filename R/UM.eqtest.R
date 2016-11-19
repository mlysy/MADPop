#' @title Simulate test statistics for two samples of individuals from the same multinomial population.
#' @description Samples chi-squared and likelihood ratio test statistics two samples from an unconstrained multinomial distribution.
#' @param N1 Size of sample 1.
#' @param N2 Size of sample 2.
#' @param p0 Common probability vector from which to draw the multinomial samples.  Can also be a matrix, in which case each simulation randomly draws with replacement from the rows of p0.
#' @param nreps Number of replications of the simulation.
#' @details The chi-squared and likelihood ratio test statistics are calculated from multinomial samples \eqn{((Y_1^1, Y_2^1), \ldots, (Y_1^M, Y_2^M))}, where
#' \deqn{Y_k^m \sim \emph{Multin}(N_k, p_0^m,}
#' where \eqn{p_0^m} is the \eqn{m}th row of \code{p0}.
#' @return An \code{nreps x 2} matrix with the simulated chi-squared values in the first column and the simulated LRT values in the second.
#' @examples
#' # bootstrapped p-value calculation against equal genotype proportions
#' # in lakes Michipicoten and Simcoe
#'
#' # contingency table
#' lnames <- c("Michipicoten", "Simcoe")
#' ctable <- UM.suff(fish215[fish215$Lake %in% lnames,])$tab
#' ctable
#'
#' # MLE of probability vector
#' p.MLE <- colSums(ctable)/sum(ctable)
#' # sample sizes
#' N1 <- sum(ctable[1,]) # Michipicoten
#' N2 <- sum(ctable[2,]) # Simcoe
#'
#' # bootstrapped test statistics (chi^2 and LRT)
#' T.boot <- UM.eqtest(N1 = N1, N2 = N2, p0 = p.MLE, nreps = 1e3)
#'
#' # observed test statistics
#' T.obs <- c(chi2 = chi2.stat(ctable), LRT = LRT.stat(ctable))
#' # p-values
#' mean(T.boot[,"chi2"] > T.obs["chi2"]) # p-value of chi^2 test
#' mean(T.boot[,"LRT"] > T.obs["LRT"]) # p-value of likelihood ratio test
#' @export
UM.eqtest <- function(N1, N2, p0, nreps, verbose = TRUE, debug = FALSE) {
  N <- N1+N2
  P0 <- p0
  if(is.matrix(P0)) {
    Np0 <- nrow(P0)
  } else {
    Np0 <- 1
    p0 <- P0
  }
  # output
  boot.out <- matrix(NA, nreps, 2)
  colnames(boot.out) <- c("chi2", "LRT")
  for(ii in 1:nreps) {
    if(verbose && (ii%%5e3 == 0)) message("bootstrap sample ", ii, "/", nreps)
    if(Np0 > 1) {
      p0 <- P0[sample(Np0, 1),]
    }
    tab.sim <- rbind(c(rmultinom(1, size = N1, prob = p0)),
                     c(rmultinom(1, size = N2, prob = p0)))
    boot.out[ii,1] <- chi2.stat(tab.sim)
    boot.out[ii,2] <- LRT.stat(tab.sim)
  }
  boot.out
}
