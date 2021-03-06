#' Chi-squared test statistic for contingency tables
#'
#' Calculates the chi-squared test statistic for a two-way contingency table.
#' @param tab A \code{K x C} matrix (contingency table) of counts. See details.
#' @return The calculated value of the chi-squared statistic.
#' @details Suppose that \code{tab} consists of counts from \eqn{K} populations (rows) in \eqn{C} categories.  The chi-squared test statistic is computed as
#' \deqn{
#'   \sum_{i=1}^K \sum_{j=1}^C (E_{ij} - O_{ij})^2/E_{ij},
#' }{
#'   \sum_ij (E_ij - O_ij)^2/E_ij,
#' }
#' where \eqn{O_{ij}}{O_ij} is the observed number of counts in the \eqn{i}th row and \eqn{j}th column of \code{tab}, and \eqn{E_{ij}}{E_ij} is the expected number of counts under \eqn{H_0} that the populations have indentical proportions in each category:
#' \deqn{
#'   E_{ij} = \frac 1 N \sum_{i=1}^K O_{ij} \times \sum_{j=1}^C O_{ij}.
#' }{
#'   E_ij = \sum_i O_ij * \sum_j O_ij / N,
#' }
#' where \eqn{N} is the total number of counts in \code{tab}.
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
