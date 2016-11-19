## ------------------------------------------------------------------------
require(mapsAllele)

## ------------------------------------------------------------------------
# fish data
head(fish215)

# Here are two lakes under consideration: Simcoe and Michipicoten
lnames <- c("Simcoe", "Michipicoten")
Xobs <- fish215[fish215$Lake %in% lnames,]


## ------------------------------------------------------------------------
Xsuff <- HW.suff(X = Xobs[,-1])
# converts the fish data Xobs into statistics including:
#   A: total set of unique alleles
#   G: total set of unique full genotypes (combination of 4 alleles)
#   Xg: counts in each unique genotype from observed data
#   H: total set of possible inherited genotypes (combination of 2 alleles)
#       with repeats being 1 inherited allele
#   Y: observed data in a simplified numerical format


## ------------------------------------------------------------------------
# fit the model by STAN with 100,000 iterations
# The STAN model is the HardyWeinberg.stan file
nsamples <- 1e4
system.time({
  HW.fit <- HW.post(nsamples = nsamples, X = Xobs[,-1], chains = 1,
                    warmup = min(1e4, floor(nsamples/10)))
  rho.post <- HW.fit$rho
})


## ------------------------------------------------------------------------
nHobs <- nrow(Xsuff$H)

## ---- echo=FALSE---------------------------------------------------------
if(FALSE) {
  par(mfrow = c(floor(sqrt(nHobs)), ceiling(nHobs/floor(sqrt(nHobs)))))
  par(mar = c(4.5,2,1,0)+.1)
  same.scale <- TRUE # set to false to give each plot its own x-axis
  if(same.scale) {
    rho.xlim <- range(rho.post)
  }
  for(ii in 1:nHobs) {
    nm <- paste0("rho[", paste(unique(Xsuff$H[ii,]), collapse = "*\".\"*"), "]",
                 collapse = "")
    if(!same.scale) {
      hist(rho.post[,ii], ylab = "", main = "",
           breaks = 100, freq = FALSE, xlab = parse(text = nm))
    } else {
      hist(rho.post[,ii], ylab = "", main = "", xlim = rho.xlim,
           breaks = 100, freq = FALSE, xlab = parse(text = nm))
    }
    # mean value and confidence interval (the red lines)
    abline(v = mean(rho.post[,ii]), col = "red")
    abline(v = quantile(rho.post[,ii], prob = c(.025, .975)), col = "red", lty = 2)
  }
}


## ------------------------------------------------------------------------
par(mfrow = c(1,1), mar = c(5,4,4,2)+.1)
nm <- paste0("rho[", apply(Xsuff$H, 1,
                                function(a) paste(unique(a), collapse = "*\".\"*")), "]")

## ---- echo=FALSE---------------------------------------------------------
boxplot(data.frame(rho.post), names = parse(text = nm), las = 2, ylab = "Probability")

## ------------------------------------------------------------------------
# Create a 2-way table with lakes as rows and observed genotypes as
#       columns (each observed genotype is given a unique identifier).
# tab.obs <- cbind(lake = as.character(Xobs[,1]),
#                  genotype = apply(Xsuff$Y, 1,
#                    function(x) paste0(sort(x), collapse = ".")))
# tab.obs <- table(tab.obs[,1], tab.obs[,2])
# tab.obs

tab.obs <- UM.suff(Xobs)$tab


## ------------------------------------------------------------------------
# observed statistics values
# functions used: chi2.stat and LRT.stat

T.obs <- c(chi2 = chi2.stat(tab.obs), LRT = LRT.stat(tab.obs))

# chi2 is the test statistic for chi^2
# LRT is the test statistic for likelihood ratio test
# p-values of the tests under H0: both lakes are the same

pchisq(T.obs, df = ncol(tab.obs)-1, lower.tail = FALSE)


## ------------------------------------------------------------------------
# number of repetitions will be 100,000
nreps <- 1e5

# Set p0 and N for the probability MLE and number of fish
p0 <- colSums(tab.obs)/sum(tab.obs) # MLE of common probability vector
N <- rowSums(tab.obs) # number of fish in each lake


## ------------------------------------------------------------------------
# Bootstrapping:
# Test statistic for chi2 and lrt
# 100,000 iterations
system.time({
  boot.out <- UM.eqtest(N1 = N[1], N2 = N[2], p0 = p0, nreps = nreps)
})


## ------------------------------------------------------------------------
par(mfrow = c(1,2))
df0 <- length(p0)-1 # degrees of freedom

## ---- echo=FALSE---------------------------------------------------------

# xlim is the vector with minimum and maximum in the chi2.boot and chi2.obs
# histogram of bootstrap samples - plot looks like a distribution
# red line represents the observed value - to differentiate between samples and observed value

for(stat.nm in c("chi2", "LRT")) {
  # bootstrap  values
  hist(boot.out[,stat.nm], freq = FALSE, breaks = 100,
       xlim = range(boot.out[,stat.nm], T.obs[stat.nm]),
       main = ifelse(stat.nm == "chi2", "Chi-Square Test", "Likelihood Ratio Test"),
       xlab = ifelse(stat.nm == "chi2", expression(chi^2), expression(Lambda)))
  # observed value
  abline(v = T.obs[stat.nm], col = "red")
  # theoretical chi^2 distr most ppl use to calculate p-value
  curve(dchisq(x, df = df0), col = "blue", add = TRUE)
  # (dchisq gives the density)
  # (pchisq gives the distribution function)
  # The legend shows the p-value for both cases
  legend(x = "topleft",
         legend = c(paste("theoretical:",
           signif(pchisq(T.obs[stat.nm], df = df0, lower.tail = FALSE), 2)),
           paste("bootstrap:", signif(mean(boot.out[,stat.nm] >= T.obs[stat.nm]), 2))),
         title = "p-value", bty = "n")
}

## ------------------------------------------------------------------------
# Bootstrapping:
# 100,000 iterations

system.time({
  boot.out2 <- HW.eqtest(N1 = N[1], N2 = N[2], H = Xsuff$H, rho = rho.post, nreps = nreps)
})

## ------------------------------------------------------------------------
# Same plots as above, but put into for loop for graphing purposes
par(mfrow = c(1,2))
nbreaks <- 100

# degrees of freedom
df0 <- length(p0)-1
stat.nm <- c("chi^2", "Lambda")


## ---- echo=FALSE---------------------------------------------------------
# For loop
for(ii in 1:2) {
  dens.cont <- density(boot.out[,ii]) # contingency table bootstrap
  hist.hw <- hist(boot.out2[,ii], breaks = 100, plot = FALSE)
  hist(boot.out2[,ii], freq = FALSE, breaks = 100,
       xlim = range(dens.cont$x, boot.out2[,ii], T.obs[ii]),
       ylim = range(dens.cont$y, hist.hw$density),
       main = parse(text = paste0("\"Distribution of \"*", stat.nm[ii], "*\" Statistic under \"*H[0]")),
       xlab = parse(text = stat.nm[ii])) # HW bootstrap
  lines(dens.cont$x, dens.cont$y, col = "red") # contingency table bootstrap
  # Observed value
  abline(v = T.obs[ii], col = "black", lwd = 2)
  # Theoretical chi^2 distr commonly used to calculate p-value
  curve(dchisq(x, df = df0), col = "blue", add = TRUE)
  # Legend with p-values
  pval <- signif(c(chi2 = pchisq(as.numeric(T.obs[ii]), df = df0, lower.tail = FALSE),
                   bcont = mean(boot.out[,ii] >= T.obs[ii]),
                   bHW = mean(boot.out2[,ii] >= T.obs[ii])), digits = 2)
  legend(x = "topleft",
         legend = parse(text = c("\"Test Statistic\"", "\"p-values\"",
                                 paste("\"Theoretical \"*(chi^2):", pval["chi2"]),
                                 paste("\"Bootstrap \"*(\"MLE\"):", pval["bcont"]),
                                 paste("\"Bootstrap \"*(\"Bayes-HW\"):", pval["bHW"]))),
         col = c("black", NA, "blue", "red", "black"),
         pch = c(NA, NA, NA, NA, 22),
         lty = c(1, NA, 1, 1, NA),
         lwd = c(2, NA, 1, 1, NA), cex = .8,
         bty = "o")
}


## ---- echo=FALSE---------------------------------------------------------
plot(cars)

