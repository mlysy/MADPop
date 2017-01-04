#' Sufficient statistics for the Hardy-Weinberg model
#'
#' Converts a matrix of genotype data into the sufficient statistics required to fit a  Hardy-Weinberg model.
#' @param X Genotype adata.  Either a \code{N x 4} matrix with \code{NA}'s indicating duplicates or a \code{N x 5} column data.frame with the first column being the \code{popId}.
#' @param popId grouping variable for \code{X}. Must be supplied if \code{X} has 4 columns.
#' @return A list with elements:
#' \itemize{
#' \item \code{A}: Vector of unique alleles names.  The allele numbers in the following quantities correspond to the indices of \code{A}.
#' \item \code{G}: 4-column matrix of unique genotype combinations.  The presence of 0's indicates that there were less that a given genotype either has less than 4 distinct alleles or that some alleles are duplicated.
#' \item \code{Xg}: Vector of counts in each unique genotype combination from \code{X}.
#' \item \code{H}: 2-column matrix of all possible chromosomes that could have generated \code{X}.  Each chromosome consists of 2 alleles.
#' \item \code{Y}: Observed data in a simplified numerical format.  That is, all allele labels are given unique positive-integer identifiers and missing or duplicated alleles are represented by 0's.
#' }
#' @examples
#' # Two fish: one with 2 alleles and one with 3
#' X <- rbind(c("a", "b", NA, NA), c("c", "a", "d", NA))
#' suff <- MADPop:::HW.suff(X)
#' suff$A # representing all of the unique alleles and 4 is repeated
#' suff$G # since there are only 2 unique full genotypes
#' suff$Xg # since there is only 1 of each genotype
#' suff$H # all possible chromosomes compatible with the observed data
#' suff$Y # dataset in simplified numerical format
#' @keywords internal
HW.suff <- function(X, popId, debug = FALSE) {
  if(ncol(X) != 4) {
    stop("X must be a matrix or data.frame with 4 columns.")
  }
  Y <- as.matrix(X)
  # convert to numeric
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
  # unique observed genotype combinations
  G <- unique(colSort(Y))
  # inherited allele combinations for observed genotypes
  H <- t(matrix(apply(G, 1, combn, m = 2), 2))
  H <- H[rowSums(H) > 0,]
  H[H[,1]==0,1] <- H[H[,1]==0,2] # replace 0's by duplicates
  H <- colSort(H) # unique order
  H <- unique(H)
  # number of observations of each unique genotype per lake
  if(missing(popId)) popId <- rep(1, nrow(Y))
  popId <- tapply(1:length(popId), popId, c, simplify = FALSE)
  Xg <- sapply(popId, function(ii) Y[ii,], simplify = FALSE)
  Xg <- t(sapply(Xg, function(x) {
    apply(G, 1, function(g) sum(apply(t(x) == g, 2, all)))
  }))
  list(A = A, G = G, Xg = Xg, H = H, Y = Y)
}
