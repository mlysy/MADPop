
context("stanmodels")

# HMC crashes on rare occasions, so avoid this by setting the seed to a known safe result
if(identical(Sys.getenv("NOT_CRAN"), "true")) {
  set.seed(NULL)
} else {
  set.seed(1)
}
source("MADPop-testfunctions.R")
data(fish215, package = "MADPop")

# fit model
nsamples <- 100
rhoId <- "Simcoe"
nobs <- sample(100:300, 1)
X <- fish215[sample(nrow(fish215), nobs, replace = TRUE),]
invisible(capture.output({
  dm.mod <- hUM.post(nsamples, X = X, rhoId = rhoId, chains = 1,
                     full.stan.out = TRUE,
                     verbose = FALSE, show_messages = FALSE)
}))
dm.mod <- dm.mod$sfit

ntest <- nrow(as.data.frame(dm.mod))

test_that("logpost_R == logpost_Stan", {
  library(rstan)
  # stan posterior
  lp.stan <- sapply(1:ntest, function(ii) {
    upars <- unconstrain_pars(dm.mod, pars = extract.iter(dm.mod, ii))
    log_prob(dm.mod, upars = upars, adjust_transform = FALSE)
  })
  # R posterior
  Xobs <- UM.suff(X)$tab
  lp.R <- sapply(1:ntest, function(ii) {
    pars <- extract.iter(dm.mod, ii)
    ae <- pars$eta * pars$alpha
    sum(apply(Xobs, 1, function(x) {
      ddirmulti(x = x, alpha = ae, log = TRUE)
    })) - 2 * log(1+pars$eta)
  })
  max_diff <- abs(diff(lp.stan - lp.R))
  for(ii in 1:length(max_diff)) {
    expect_equal(max_diff[ii], 0)
  }
})
