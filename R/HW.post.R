#' @title Posterior draws from a fitted Hardy-Weinberg model
#' @param nsamples Number of posterior samples
#' @param X 4-column matrix of observations in the correct format.
#' @param full.stan.out Logical.  Whether or not to return the full stan output.  For monitoring convergence of the MCMC sampling.
#' @param ... Further arguments to be passed to Stan.
#' @return A list with elements
#' \itemize{
#'   \item \code{A}: The unique allele names.
#'   \item \code{H}: The 2-column matrix potential chromosomes.
#'   \item \code{rho}: A matrix with \code{ncol(rho) == nrow(H)}, where each row is a draw from the posterior distribution of inheritance probabilities.
#'   \item \code{sfit}: If \code{full.stan.out = TRUE}, the fitted Stan object.
#' }
HW.post <- function(nsamples, X, full.stan.out = FALSE, ...) {
  suff <- HW.suff(X)
  sfit <- sampling(object = HW.mod, data = HW.data(suff), iter = nsamples, ...)
  ans <- list(A = suff$A, H = suff$H, rho = extract(sfit, permuted = TRUE)$rho)
  if(full.stan.out) ans <- c(ans, list(sfit = sfit))
  ans
}
