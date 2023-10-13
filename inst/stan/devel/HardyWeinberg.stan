// STAN model for Hardy-Weinberg with variable/unobserved allele counts
// {mlysy,pwjkim}@uwaterloo.ca, november 2014

data {
  int<lower=1> nG; // number of observed  genotype categories
  int<lower=1> nH; // number of inherited genotype categories
  array[nG] int<lower=0> Xg; // number of observations for each of the observed genotype categories
  // inheritance information
  // for each observed genotype category, each (unordered) combination of inherited genotypes
  // which could make up the observation.  For two locations, there are at most 6 combinations.
  array[2,6,nG] int<lower=0,upper=nH> Hg;
  array[nG] int<lower=0,upper=6> nHg; // number of inheritance combinations for each observed genotype category.
}

transformed data {
  array[nG,6] int ordCount; // how many times to count the order of each inheritance pattern (once or twice)
  for(ii in 1:(nG-1)) {
    for(jj in 1:nHg[ii]) {
      if(Hg[1,jj,ii] == Hg[2,jj,ii]) {
        ordCount[ii,jj] = 1;
      } else {
        ordCount[ii,jj] = 2;
      }
    }
  }
}

parameters {
  simplex[nH] rho; // probability of inheriting each of the elements in H
}

transformed parameters {
  simplex[nG] theta; // probability of each of the observed genotypes
  theta[nG] = 1;
  for(ii in 1:(nG-1)) {
    theta[ii] = 0;
    for(jj in 1:nHg[ii]) {
      theta[ii] = theta[ii] + ordCount[ii,jj] * rho[Hg[1,jj,ii]] * rho[Hg[2,jj,ii]];
    }
    theta[nG] = theta[nG] - theta[ii];
  }
}

model {
  Xg ~ multinomial(theta);
}
