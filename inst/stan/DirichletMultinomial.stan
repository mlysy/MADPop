/* STAN model for the Dirichlet-Multinomial Distribution

   {mlysy,pwjkim}@uwaterloo.ca, january 2015

Model is:
```
X_l | rho_l ~ Multinomial(N_l, rho_l)
rho_l ~ Dirichlet(eta * alpha),
```
where `sum(alpha = 1)`.  Thus, the prior expectation is `E[rho_l] = alpha`.

For now, let's assume that `eta` is given, though ideally it will be estimated from the data.
*/

functions {
  
  //' Unnormalized Dirichlet-Multinomial distribution.
  //'
  //' For likelihood only.
  //' No longer needed since now provided by Stan.
  // real dirichlet_multinomial_lpmf(array[] int x, vector eta) {
  //   real ans;
  //   ans = 0.0;
  //   for(ii in 1:num_elements(x)) {
  //     ans += lgamma(x[ii] + eta[ii]) - lgamma(eta[ii]);
  //   }
  //   return ans + lgamma(sum(eta)) - lgamma(sum(x)+sum(eta));
  // }
  
  //' Unnormalized Dirichlet-Multinomial distribution with vectorized inputs.
  //'
  //' @param X Matrix of DM observations, one on each row.
  //' @param eta Vector of positive DM parameters.
  //' @return The loglikelihood as a function of `eta`.  In other words, sum of the log-PDFs for each for of `X`, but dropping anything that does not involve `eta`.
  real Dirichlet_Multinomial_lpmf(array[,] int X, vector eta) {
    array[2] int D;
    real ans;
    real seta;
    real slgeta;
    D = dims(X);
    // eta constants
    seta = sum(eta);
    slgeta = 0.0;
    for(jj in 1:D[2]) {
      slgeta += lgamma(eta[jj]);
    }
    // all of X
    ans = 0.0;
    for(ii in 1:D[1]) {
      for(jj in 1:D[2]) {
	ans += lgamma(X[ii,jj] + eta[jj]);
      }
      ans -= lgamma(sum(X[ii])+seta);
    }
    ans += D[1] * (lgamma(seta) - slgeta);
    return ans;
  }
}

data {
  int<lower=1> nG; // number of observed genotype categories
  int<lower=1> nL; // number of lakes
  array[nL,nG] int<lower=0> X; // vector of counts for each lake
  //real<lower=0> eta; // hyper prior precision parameter
  int<lower=0,upper=nL> nLrho; // number of lakes for which to generate samples from rho
  array[nLrho] int<lower=1,upper=nL> iLrho; // indices of these lakes
}

parameters {
  simplex[nG] alpha; // probability parameters for dirichlet-multinomial
  real<lower=0> eta; // precision parameter
}

model {
  //vector<lower=0> kappa[nG]; // dirichlet parameters
  //kappa = eta * alpha;
  X ~ Dirichlet_Multinomial(eta * alpha);
  // avoid ridiculously small dirichlet variance
  // using p(1/sum(eta)) \propto 1
  //increment_log_prob(-2*log(1+eta));
  target += -2*log(1+eta);
}

generated quantities {
  // individual lake multinomial parameters
  array[nLrho] simplex[nG] rho;
  if(nLrho > 0) {
    for(ii in 1:nLrho) {
      rho[ii] = dirichlet_rng(to_vector(X[iLrho[ii]]) + eta * alpha);
    }
  }
}
