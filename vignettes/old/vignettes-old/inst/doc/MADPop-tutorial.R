## ------------------------------------------------------------------------
require(MADPop)

## ------------------------------------------------------------------------
head(fish215)

## ------------------------------------------------------------------------
Xsuff <- UM.suff(fish215)

## ------------------------------------------------------------------------
head(Xsuff$A)

## ------------------------------------------------------------------------
dim(Xsuff$G) # 134 unique genotypes
Xsuff$G[1,] # observed genotype is "3"
Xsuff$G[20,] # observed genotype is "1.8"

## ------------------------------------------------------------------------
Xsuff$tab[,1:10] # 11 lakes and 120 more genotypes

## ------------------------------------------------------------------------
lnames <- c("Michipicoten", "Simcoe") # two lakes under consideration
# popId vector which merges these two lakes under the null
popId0 <- paste(lnames, collapse = ".")
popId0 # common lake name
popId <- as.character(fish215$Lake)
popId[popId %in% lnames] <- popId0
table(popId) # all lake names
allele0 <- cbind(popId, fish215[,-1]) # data with merged lakes

# posterior inference in stan
# don't make too big during development or this will take forever
# eventually should use nsamples = 1e4
nsamples <- 1e3
hUM.fit <- hUM.post(nsamples = nsamples, X = allele0,
                    rhoId = popId0, # output only merged lakes to save space
                    chains = 1, # next two arguments are passed to stan
                    warmup = min(1e4, floor(nsamples/10)))
rho.post <- hUM.fit$rho[,1,] # convert 3-d array to 2-d matrix

## ---- fig.width = 7, fig.height = 3.5------------------------------------
propId <- order(colMeans(rho.post))[1:20]
par(mar = c(5,4,.1,.1))
boxplot(rho.post[,propId],
        las = 2, pch = ".", col = "grey", cex.axis = .8,
        xlab = "Genotype", ylab = "Probability")

## ------------------------------------------------------------------------
# number of observations in each lake
N1 <- sum(Xsuff$tab["Michipicoten",])
N2 <- sum(Xsuff$tab["Simcoe",])
# generate samples from p(T | data, H_0)
T.post <- UM.eqtest(N1 = N1, N2 = N2, p0 = rho.post, nreps = 1e4)
head(T.post) # computes both chi^2 and LRT
# calculate T.obs
tab.obs <- Xsuff$tab[lnames,] # observed table
# automatically removes categories where both lakes have zero observations
T.obs <- c(chi2 = chi2.stat(tab.obs),
           LRT = LRT.stat(tab.obs))
p.post <- c(chi2 = mean(T.post[,"chi2"] >= T.obs["chi2"]),
            LRT = mean(T.post[,"LRT"] >= T.obs["LRT"]))
p.post # posterior p-values

