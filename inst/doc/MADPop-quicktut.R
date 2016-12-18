## ---- echo = FALSE-------------------------------------------------------
if(FALSE) {
  #setwd("C:/Users/Peter/Dropbox/test/MADPop")
  setwd("c:/Users/Jerome/Dropbox/Shared/allele/MADPop")
  pkg.path <- getwd()
  require(rmarkdown)
  rmarkdown::render(file.path(pkg.path, "vignettes", "MADPop-quicktut.Rmd"))
  # view it by opening MADPop/vignettes/MADPop-tutorial.html
}

## ------------------------------------------------------------------------
require(MADPop)

## ---- echo = FALSE-------------------------------------------------------
nObs <- nrow(fish215)
nPop <- nlevels(fish215$Lake)
nAlleles <- length(table(c(as.matrix(fish215[,-1])))) - 1

## ------------------------------------------------------------------------
head(fish215[sample(nObs),])

## ----setup2, ref.label = "table2", echo = FALSE, results = "hide"--------
popID <- c("Michipicoten", "Simcoe")   # lakes to compare
Xsuff <- UM.suff(fish215)             # summary statistics for dataset
ctab <- Xsuff$tab[popID,]            # contingency table
ctab <- ctab[,colSums(ctab) > 0] # remove alleles with no counts
#ctab
rbind(ctab, Total = colSums(ctab))

## ----table2--------------------------------------------------------------
popID <- c("Michipicoten", "Simcoe")   # lakes to compare
Xsuff <- UM.suff(fish215)             # summary statistics for dataset
ctab <- Xsuff$tab[popID,]            # contingency table
ctab <- ctab[,colSums(ctab) > 0] # remove alleles with no counts
#ctab
rbind(ctab, Total = colSums(ctab))

## ------------------------------------------------------------------------
gtype <- colnames(ctab)[1]
gtype <- as.numeric(strsplit(gtype, "[.]")[[1]])
gtype
names(gtype) <- paste0("A", gtype)
sapply(gtype, function(ii) Xsuff$A[ii])

## ------------------------------------------------------------------------
# observed values of the test statistics
chi2.obs <- chi2.stat(ctab) # Pearson's chi^2
LRT.obs <- LRT.stat(ctab) # LR test statistic
T.obs <- c(chi2 = chi2.obs, LRT = LRT.obs)
# p-value with asymptotic calculation
C <- ncol(ctab)
pv.asy <- pchisq(q = T.obs, df = C-1, lower.tail = FALSE)
signif(pv.asy, 2)

## ---- cache = TRUE-------------------------------------------------------
N1 <- sum(ctab[1,])                     # size of first sample
N2 <- sum(ctab[2,])                     # size of second sample
rho.hat <- colSums(ctab)/(N1+N2)        # common probability vector
# bootstrap distribution of the test statistics
# set verbose = TRUE for progress output
system.time({
  T.boot <- UM.eqtest(N1 = N1, N2 = N2, p0 = rho.hat, nreps =1e4,
                      verbose = FALSE)
})
# bootstrap p-value
pv.boot <- rowMeans(t(T.boot) >= T.obs)
signif(pv.boot, 2)

## ------------------------------------------------------------------------
itab1 <- colSums(ctab) == 1             # single count genotypes
cbind(ctab[,itab1],
      Other = rowSums(ctab[,!itab1]),
      Total = rowSums(ctab))

## ---- echo = FALSE-------------------------------------------------------
c1 <- sum(itab1) # number of single-count columns
n1 <- sum(ctab[,itab1]) # number of single counts

