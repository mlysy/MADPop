#' @title Posterior draws from a fitted hierarchical Unconstrained-Multinomial model
#' @description MCMC sampling from a Dirichlet-Multinomial model using \code{\link[rstan]{stan}}.
#' @param nsamples Number of posterior samples
#' @param X 4-column or 5-column matrix of observations in the correct format.  See \code{\link{UM.suff}}.
#' @param popId Optional vector of population identifiers.  See \code{\link{UM.suff}}.
#' @param rhoId Populations for which rho probabilities are desired.  Defaults to all populations.  Set \code{rhoId = NULL} not to compute these for any populations.
#' @param full.stan.out Logical.  Whether or not to return the full \code{stan} output.  For monitoring convergence of the MCMC sampling.
#' @param ... Further arguments to be passed to \code{stan}.
#' @return A list with elements
#' \itemize{
#'   \item \code{A}: The unique allele names.
#'   \item \code{G}: The 4-column matrix of unique genotype combinations.
#'   \item \code{rho}: A matrix with \code{ncol(rho) == nrow(G)}, where each row is a draw from the posterior distribution of inheritance probabilities.
#'   \item \code{sfit}: If \code{full.stan.out = TRUE}, the fitted \code{stan} object.
#' }
#' @details The model is given by
#' \deqn{Y_k \sim_{ind} \emph{Multin}(\rho_k, N_k), \qquad \rho_k \sim_{iid} \emph{Dir}(\alpha), \qquad (\bar \alpha, \alpha_0) \sim 1/(1+\alpha_0)^2.}
#' @examples
#' # fit hierarchical model to fish215 data
#'
#' # only output posterior samples for lake Simcoe
#' rhoId <- "Simcoe"
#' nsamples <- 500
#' hUM.fit <- hUM.post(nsamples = nsamples, X = fish215,
#'                     rhoId = rhoId,
#'                     chains = 1) # number of MCMC chains
#' rho.post <- hUM.fit$rho[,1,] # convert 3-d array to 2-d matrix
#'
#' # plot first 20 posterior probabilities in lake Simcoe
#' boxplot(rho.post[,1:20], las = 2,
#'         xlab = "Genotype", ylab = "Posterior Probability",
#'         pch = ".", col = "grey")
#' @export
hUM.post <- function(nsamples, X, popId, rhoId,
                     full.stan.out = FALSE,..., debug = FALSE) {
  suff <- UM.suff(X, popId)
  Xobs <- suff$tab
  if(missing(rhoId)) {
    rhoId <- rownames(Xobs)
  }
  if(is.null(rhoId)) {
    iLrho <- numeric(0)
  } else {
    iLrho <- sapply(rhoId, function(ri) {
      if(ri %in% rownames(Xobs)) {
        ans <- which(ri == rownames(Xobs))
      } else {
        ans <- 0
      }
      ans
    })
    if(any(iLrho == 0)) stop("rhoId must be a subset of popId.")
  }
  nLrho <- length(iLrho)
  if(debug) browser()
  sfit <- sampling(object = hUM.mod,
                   data = list(nG = ncol(Xobs), nL = nrow(Xobs),
                     nLrho = nLrho, iLrho = array(iLrho, dim = nLrho),
                     X = Xobs),
                   iter = nsamples, ...)
  #alpha <- extract(sfit, permuted = TRUE)$alpha
  #colnames(alpha) <- colnames(Xobs)
  ans <- list(A = suff$A, G = suff$G)
  if(nLrho > 0) {
    rho <- extract(sfit, permuted = TRUE)$rho
    dimnames(rho)[[2]] <- rownames(Xobs)[iLrho]
    dimnames(rho)[[3]] <- colnames(Xobs)
    ans <- c(ans, list(rho = rho))
  }
  if(full.stan.out) ans <- c(ans, list(sfit = sfit))
  ans
}
