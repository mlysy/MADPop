#--- quick tutorial for hierarchical model --------------------------------------

require(MADPop)

# fish data
head(fish215)

# Here are two lakes under consideration: Simcoe and Michipicoten
lnames <- c("Michipicoten", "Simcoe")

#--- quick table ----------------------------------------------------------------

Xsuff <- UM.suff(fish215)
tab.obs <- Xsuff$tab[lnames,]
tab.obs <- tab.obs[,colSums(tab.obs)>0]

ftab <- function(nc = 10, rot = 90, debug = FALSE) {
  if(debug) browser()
  message("\\begin{tabular}{r", rep("l", nc), "}")
  for(II in seq(1, ncol(tab.obs), by = nc)) {
    ii <- II:min(II+nc-1,ncol(tab.obs))
    message("& ", paste("{\\small \\rotatebox{", rot, "}{", colnames(tab.obs)[ii], "}}",
                        sep = "", collapse = " & "), "\\\\")
    message("\\hline")
    message(lnames[1], " & ", paste(tab.obs[1,ii], collapse = " & "), " \\\\")
    message(lnames[2], " & ", paste(tab.obs[2,ii], collapse = " & "), " \\\\")
    message("\\\\")
  }
  message("\\end{tabular}")
}

# probably can't display this LaTeX table directly with pandoc...
# perhaps can modify?
ftab(13, rot = 90)

#--- hierarchical model ---------------------------------------------------------

# merge michipicoten and simcoe
popId0 <- paste(lnames, collapse = ".")
popId0
popId <- as.character(fish215$Lake)
popId[popId %in% lnames] <- popId0

# sample from the posterior distribution
# this has 11 lakes, with Simcoe and Michipicoten having been grouped together
# under H0: they are from the same population.
# rhoId specifies which of the lakes should the MCMC sampler return estimates for.
# in this setup we only care about Simcoe/Michipicoten, so the parameter estimates for
# the other models are neither returned nor computed (faster and less space).

nsamples <- 1e4
system.time({
  hUM.fit <- hUM.post(nsamples = nsamples, X = cbind(popId, fish215[,-1]),
                      rhoId = popId0,
                      chains = 1, warmup = min(1e4, floor(nsamples/10)))
})
rho.post <- hUM.fit$rho[,1,] # convert 3-d array to 2-d matrix

#--- boxplot of estimated probabilities -----------------------------------------

# display those with highest posterior means
# also color by count in observed lakes

# combined genotype counts in two lakes under comparison
rho.count <- colSums(Xsuff$tab[lnames,])
# order genotypes by decreasing posterior mean
rho.ord <- order(colMeans(rho.post), decreasing = TRUE)

# genotypes with highest posterior prob
nplt <- 50
clrs <- c("#999999", "#E69F00", "#56B4E9",
          "#009E73", "#F0E442", "#0072B2",
          "#D55E00", "#CC79A7")
rho.ct <- rho.count[rho.ord[1:nplt]]
rho.lgd <- unique(sort(rho.ct))
rho.box <- rho.post[,rho.ord[1:nplt]]

par(mar = c(5,4.5,.5,.5)+.1)
boxplot(x = rho.box,
        las = 2, col = clrs[rho.ct+1], pch = 16, cex = .2)
title(xlab = "Genotype", line = 4)
title(ylab = expression(p(gamma[(12)]*" | "*Y)))
legend("topright", legend = rho.lgd, fill = clrs[rho.lgd+1],
       title = "Counts")


#--- test of equality -----------------------------------------------------------

# create two-way table on desired lakes
tab.obs <- Xsuff$tab[lnames,]
tab.obs <- tab.obs[,colSums(tab.obs)>0] # remove empty columns
p0 <- colSums(tab.obs)/sum(tab.obs) # MLE of common probability vector
N <- rowSums(tab.obs) # number of fish in each lake

# observed test statistics
T.obs <- c(chi2 = chi2.stat(tab.obs), LRT = LRT.stat(tab.obs))

# sampling distribution of test statistics
nreps <- 1e5

# bootstrap
system.time({
  T.mle.boot <- UM.eqtest(N1 = N[1], N2 = N[2], p0 = p0, nreps = nreps)
})

# hierarchical model posterior
system.time({
  T.hUM.post <- UM.eqtest(N1 = N[1], N2 = N[2], p0 = rho.post, nreps = nreps)
})

# p-value table
# for bootstrap and MCMC
pv.tab <- sapply(list(Boot = T.mle.boot, Hier = T.hUM.post),
                 function(smp) rowMeans(t(smp) >= T.obs))
# for asymptotic chi^2
df0 <- length(p0)-1 # degrees of freedom
pv.tab <- cbind(Asym = pchisq(T.obs, df = df0, lower.tail=FALSE), pv.tab)
signif(pv.tab, 2)

# also plot
clrs <- c("blue", "red", "brown")
test.nm <- expression(chi, Lambda)
stat.nm <- c("chi[C-1]^2", "\"Bootstrap\"", "\"Hierarchical\"")
par(mfrow = c(1,2), mar = c(2,4,1,.2)+.1)
for(ii in 1:2) {
  # densities
  lboot <- density(T.mle.boot[,ii])
  lboot <- cbind(x = lboot$x, y = lboot$y)
  lhier <- density(T.hUM.post[,ii])
  lhier <- cbind(x = lhier$x, y = lhier$y)
  lchi2 <- seq(min(lboot[,1],lhier[,1]), max(lboot[,1],lhier[,1]),
               len = 200)
  lchi2 <- cbind(x = lchi2, y = dchisq(lchi2, df = df0))
  # plot them
  plot(lchi2[,1], lchi2[,2], type = "l", col = clrs[1],
       xlim = range(lboot[,1], lhier[,1], lchi2[,1]),
       ylim = range(lboot[,2], lhier[,2], lchi2[,2]),
       main = test.nm[ii], ylab = "Density", xlab = "")
  lines(lboot[,1], lboot[,2], col = clrs[2])
  lines(lhier[,1], lhier[,2], col = clrs[3])
  # add true test stat
  abline(v = T.obs[ii], lty = 2)
  # legend with pval
  lgd <- parse(text = paste0(stat.nm, "*\": \"*", signif(pv.tab[ii,],2)))
  legend("topright", legend = lgd, fill = clrs, title = "p-values")
}
