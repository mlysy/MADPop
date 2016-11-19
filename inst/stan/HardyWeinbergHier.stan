// STAN model for Hierarchical (e.g. multi-lake) Hardy-Weinberg with variable/unobserved allele counts
// {mlysy,pwjkim}@uwaterloo.ca, november 2014

data {
  int<lower=1> nL; // number of lakes
  int<lower=1> nG; // number of observed  genotype categories (total)
  int<lower=1> nH; // number of inherited genotype categories (total)
  int<lower=0> Xg[nL,nG]; // number of observations for each of the observed genotype categories per lake
  // inheritance information
  // for each observed genotype category, each (unordered) combination of inherited genotypes
  // which could make up the observation.  For two locations, there are at most 6 combinations.
  int<lower=0,upper=nH> Hg[2,6,nG];
  int<lower=0,upper=6> nHg[nG]; // number of inheritance combinations for each observed genotype category.
//  real mu0; // prior parameters on eta0
//  real<lower=0> sigma0;
  real<lower=0> eta0; // scale parameter for prior
}

transformed data {
  int ordCount[nG,6]; // how many times to count the order of each inheritance pattern (once or twice)
  for(ii in 1:(nG-1)) {
    for(jj in 1:nHg[ii]) {
      if(Hg[1,jj,ii] == Hg[2,jj,ii]) {
        ordCount[ii,jj] <- 1;
      } else {
        ordCount[ii,jj] <- 2;
      }
    }
  }
}

parameters {
  simplex[nH] rho[nL]; // probability of inheriting each of the elements in H
  simplex[nH] eta; // hyperparameters
//  real<lower=0> eta0; // scale multiple for prior
}

transformed parameters {
  simplex[nG] theta[nL]; // probability of each of the observed genotypes in each lake
  for(kk in 1:nL) {
    theta[kk,nG] <- 1;
    for(ii in 1:(nG-1)) {
      theta[kk,ii] <- 0;
      for(jj in 1:nHg[ii]) {
        theta[kk,ii] <- theta[kk,ii] + ordCount[ii,jj] * rho[kk,Hg[1,jj,ii]] * rho[kk,Hg[2,jj,ii]];
      }
      theta[kk,nG] <- theta[kk,nG] - theta[kk,ii];
    }
  }
}

model {
//eta0 ~ lognormal(mu0, sigma0);
  for(kk in 1:nL) {
    rho[kk] ~ dirichlet(eta * eta0);
    Xg[kk] ~ multinomial(theta[kk]);
  }
}
