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


