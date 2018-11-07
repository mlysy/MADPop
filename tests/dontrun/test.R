#--- check stan code -------------------------------------------------------

# stan code
DM.mod <- stanc(model_code = "
functions {
  // unnormalized dirichlet-multinomial distribution
  // (for likelihood only)
  real dirichlet_multinomial_lpmf(int[] x, vector eta) {
    real ans;
    ans = 0.0;
    for(ii in 1:num_elements(x)) {
      ans = ans + lgamma(x[ii] + eta[ii]) - lgamma(eta[ii]);
    }
    return ans + lgamma(sum(eta)) - lgamma(sum(x)+sum(eta));
  }
  // same thing but vectorized, i.e. accepts matrix X
  // current bug prevents overloading...
  real Dirichlet_Multinomial_lpdf(matrix X, vector eta) {
    int D[2];
    real ans;
    real seta;
    real slgeta;
    D = dims(X);
    // eta constants
    seta = sum(eta);
    slgeta = 0.0;
    for(jj in 1:D[2]) {
      slgeta = slgeta + lgamma(eta[jj]);
    }
    // all of X
    ans = 0.0;
    for(ii in 1:D[1]) {
      for(jj in 1:D[2]) {
	ans = ans + lgamma(X[ii,jj] + eta[jj]);
      }
      ans = ans - lgamma(sum(X[ii])+seta);
    }
    ans = ans + D[1] * (lgamma(seta) - slgeta);
    return ans;
  }
}
model {}
")

expose_stan_functions(DM.mod)

# sample from dirichlet distribution
rdirichlet <- function(n, alpha) {
  K <- length(alpha)
  X <- matrix(rgamma(n*K, shape = alpha), n, K, byrow = TRUE)
  sweep(X, 1, rowSums(X), "/")
}

# dirichlet-multinomial distribution
# x and alpha are vectors of the same length.
# just like dmultinom, this density is not vectorized (i.e., doesn't accept matrix inputs)
ddirmulti <- function(x, alpha, log = FALSE) {
  sx <- sum(x)
  sa <- sum(alpha)
  ans <- lgamma(sx+1) - sum(lgamma(x+1))
  ans <- ans + lgamma(sa) - lgamma(sx+sa)
  ans <- ans + sum(lgamma(x+alpha)) - sum(lgamma(alpha))
  if(!log) ans <- exp(ans)
  ans
}

# dirichlet-multinomial logliklihood
# vectorized to accept a matrix alpha and drops constants in x
# output length is nrow(alpha)
ldirmulti <- function(alpha, x) {
  sx <- sum(x)
  sa <- rowSums(alpha)
  ans <- lgamma(sa) - lgamma(sx+sa)
  ans + rowSums(lgamma(sweep(alpha, 2, x, "+")) - lgamma(alpha))
}

require(MADPop)

#expose_stan_functions(MADPop:::hUM.mod)

# generate data
nG <- 5 # number of observed genotypes
nL <- 3 # number of lakes
N <- sample(20:30, nL, replace = TRUE) # lake size
eta <- rexp(nG, .1)
Rho <- rdirichlet(nL, eta)
X <- matrix(NA, nL, nG)
for(ii in 1:nL) {
  X[ii,] <- rmultinom(1, N[ii], Rho[ii,])
}

# parameters
Eta <- replicate(10, rexp(nG, .1))

dirichlet_multinomial_lpmf(x = X[1,], eta = eta)
Dirichlet_Multinomial_lpdf(X = X[1,,drop=FALSE], eta = eta)

# loglikelihood, single lake
ll1 <- apply(Eta, 2, function(et) {
  dirichlet_multinomial_lpmf(x = X[1,], eta = et)
})
ll2 <- apply(Eta, 2, function(et) {
  ddirmulti(x = X[1,], alpha = et, log = TRUE)
})
ll3 <- apply(Eta, 2, function(et) {
  Dirichlet_Multinomial_lpmf(X = as.list(X[1,]), eta = et)
})

ll1-ll2

# loglikelihood, multiple lakes
LL1 <- apply(Eta, 2, function(et) {
  Dirichlet_Multinomial_lpdf(X = X, eta = et)
})
LL2 <- apply(Eta, 2, function(et) {
  sum(apply(X, 1, function(x) {
    ddirmulti(x = x, alpha = et, log = TRUE)
  }))
})

LL1-LL2



for(ii in 1:1000) {
  set.seed(ii)
  if(!(ii %% 5)) message("ii = ", ii)
  source("test-stanmodels.R", local = TRUE)
}
