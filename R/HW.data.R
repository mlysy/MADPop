# Internal function
# convert HW sufficient statistics to STAN input
# input is the result of a cal to HW.suff
# output is a list with the following elements:
# \enumerate{
#  \item \code{nG}: number of nonzero full genotype observations
#  \item \code{nH}: number of nonzero inherited genotype observations
#  \item \code{Xg}: counts in each unique genotype from observed data
#  \item \code{Hg}: full genotype information - for each genotype, list all inherited genotypes \code{H} for each parent (unordered), at most 6 per genotype
#  \item \code{nHg}: absolute value of \code{Hg} for each observed genotype
# }
# NOTE: The potential chromosomes in the observed data can produce more genotypes than are found in the observed data, in which case an observation with 0 counts is added to the group of all these unobserved counts.
HW.data <- function(suff, debug = FALSE) {
  nG <- nrow(suff$G)
  H <- t(suff$H)
  nH <- ncol(H)
  Xg <- suff$Xg
  nL <- nrow(Xg)
  Hg <- array(0, dim = c(2, 6, nG))
  nHg <- rep(NA, nG)
  for(ii in 1:nG) {
    g <- suff$G[ii,]
    g <- g[g>0]
    singleA <- length(g) == 1
    # all potential allele combinations which could generate g
    hg <- rbind(g, g)
    if(!singleA) hg <- cbind(hg, combn(g,2))
    hg <- t(hg)
    #hg <- rbind(t(cbind(, combn(g,2))))
    if(singleA) {
      hg <- cbind(hg, hg)
    } else {
      hg <- cbind(hg[rep(1:nrow(hg), nrow(hg)),], hg[rep(1:nrow(hg), each=nrow(hg)),])
    }
    hg <- t(apply(hg, 1, function(x) {
      or <- order(x[c(1,3)], x[c(2,4)])
      if(all(or == 2:1)) x <- x[c(3:4, 1:2)]
      x
    }))
    hg <- unique(hg)
    # remove invalid combinations
    ind <- apply(hg, 1, function(y) {
      tmp <- unique(sort(y))
      length(tmp) == length(g) && all(tmp == g)
    })
    hg <- hg[ind,]
    if(singleA) hg <- t(hg)
    dimnames(hg) <- NULL
    # 2-row matrix containing the elements of Hg
    hg <- apply(hg, 1, function(y)
      c(h1 = which(colSums(y[1:2] == H) == 2),
        h2 = which(colSums(y[3:4] == H) == 2)))
    nHg[ii] <- ncol(hg)
    hg <- cbind(hg, matrix(0, 2, 6-nHg[ii]))
    Hg[,,ii] <- hg
  }
  if(nG < nH) {
    # add the "null" genotype to account for all genotypes with 0 counts implied by Hardy-Weinberg
    nG <- nG+1
    Xg <- cbind(Xg, 0)
    Hg <- array(c(Hg, rep(0, 12)), dim = c(2, 6, nG))
    nHg <- c(nHg, 0)
  }
  if(debug) browser()
  if(nL == 1) {
    ans <- list(nG = nG, nH = nH, Xg = c(Xg), Hg = Hg, nHg = nHg)
  } else {
    ans <- list(nG = nG, nH = nH, nL = nL, Xg = Xg, Hg = Hg, nHg = nHg)
  }
  ans
}
