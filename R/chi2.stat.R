#' Chi-squared test statistic for contingency tables
#'
#' Calculates the chi-squared test statistic for a two-way contingency table.
#' @param tab A \code{K x N} matrix (contingency table) of counts. See details.
#' @return The calculated value of the chi-squared statistic.
#' @details The chi-squared test statistic is computed as
#' \deqn{\sum_{i=1}^K \sum_{j=1}^N (E_{ij} - O_{ij})^2/E_{ij},}
#' where \eqn{O_{ij} = \code{tab[i,j]}} and
#' \deqn{E_{ij} = \code{sum(tab[,i])*sum(tab[,j])/sum(tab)}.}
#' That is, the null hypothesis is that each row of the contingency table is drawn from the same multinomial distribution.  If any column has only zeros it is removed before calculating the test statistic.
#' @examples
#' # simple contingency table
#' ctab <- rbind(pop1 = c(5, 3, 0, 3),
#'                 pop2 = c(4, 10, 2, 5))
#' colnames(ctab) <- LETTERS[1:4]
#' ctab
#' chi2.stat(ctab) # chi^2 test statistic
#' @export
chi2.stat <- function(tab) {
  N <- rowSums(tab) # number of fish in each lake
  p0 <- colSums(tab)/sum(N) # MLE of common probability vector
  tabE <- N %o% p0 # expected
  sum(((tab-tabE)^2/tabE)[tabE > 0]) # chi2
}
