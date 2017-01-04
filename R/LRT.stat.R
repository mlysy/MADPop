#' Likelihood ratio test statistic for contingency tables
#'
#' Calculate the likelihood ratio test statistic for general two-way contingency tables.
#'
#' @param tab A \code{K x N} matrix (contingency table) of counts. See details.
#' @return The calculated value of the LRT statistic.
#' @details The likelihood ratio test statistic is computed as
#' \deqn{2 \sum_{i=1}^K \sum_{j=1}^N O_{ij} \log(p^A_{ij}/p^0_{j}),}
#' where \eqn{O_{ij} = \code{tab[i,j]}}, \eqn{p^A_{ij} = \code{tab[i,j]/sum(tab)}}, and
#' \deqn{p^0_j = \code{sum(tab[,j])/sum(tab)}.}
#' That is, the null hypothesis is that each row of the contingency table is drawn from the same multinomial distribution.  If any column has only zeros it is removed before calculating the test statistic.
#' @examples
#' # simple contingency table
#' ctab <- rbind(pop1 = c(5, 3, 0, 3),
#'                 pop2 = c(4, 10, 2, 5))
#' colnames(ctab) <- LETTERS[1:4]
#' ctab
#' LRT.stat(ctab) # likelihood ratio statistic
#' @export
LRT.stat <- function(tab) {
  K <- nrow(tab) # number of lakes
  N <- rowSums(tab) # number of fish in each lake
  p0 <- colSums(tab)/sum(N) # MLE of common probability vector
  lp0 <- matrix(rep(log(p0), each = K), nrow = K) # log term
  L0 <- sum((tab * lp0)[lp0 > -Inf]) # loglik under H0
  LA <- sum((tab * log(tab/N))[tab > 0]) # loglik under HA
  2 * (LA - L0) # LRT (>= 0)
}
