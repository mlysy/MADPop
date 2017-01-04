#--- file to setup parts of the package -----------------------------------------
# won't be included in bundled version.

require(devtools)
#require(roxygen2)
require(knitr)
require(rmarkdown)
require(rstan)

# set this field to folder where the package is located
pkg.path <- getwd()
pkg.name <- "MADPop"

#--- quick install --------------------------------------------------------------

document(pkg.path) # updates documentation

# run this when you want to "compile" the vignette
# will then be visible when you install but takes longer
build_vignettes(pkg = pkg.path)

install(pkg.path) # install
# always quit and restart R after this step.

# optionally, to produce the .tar.gz that can be distributed
build(pkg.path)

# check...
check(pkg.path)

#--- vignette updating ----------------------------------------------------------
# check vignette
rmarkdown::render(file.path(pkg.path, "vignettes", "MADPop-tutorial.Rmd"))
# view it by opening MADPop/inst/doc/MADPop-tutorial.html

#--- test some of the code ------------------------------------------------------

# after quitting and restarting:
require(MADPop)

tmp <- hUM.post(nsamples = 100, X = fish215, rhoId = NULL,
                chains = 1, debug = FALSE)






#--- SCRATCH --------------------------------------------------------------------

# OK we'll need to provide the following functions/data.

# fish215: dataset

# HW model
# HW.suff: takes a dataset and extracts the HW sufficient statistics, i.e.
# everything needed for posterior inference.
# HW.post: a wrapper for the Stan model fitting.  User-friendly means more constraints.
# so I think arguments should be: suff, nsamples, stan.out, ...
# the ellipsis is for further arguments to be passed to the stan MCMC sampler.
# HW.sim: should be able to generate data from the HW model.
# arguments are sample size, H, the list of chromosomes, rho, prob vector for each of these (this is geno.sim)
# chi2.stat and LRT.stat
# boot.eq.HW, boot.eq.UM: bootstrap distributions of statistics under various nulls.

# then for DM model
# DM.post
# DM.sim
# can use boot.eq.UM for generating data from the appropriate null.

#--- documentation --------------------------------------------------------------

# run this when you modify the documentation
maps.compile <- FALSE
document(pkg = pkg.path)
# remove "undocumented" functions
undoc <- c("HW.eqtest", "HW.post", "HW.sim", "HW.suff")
if(!is.null(undoc)) {
  sapply(undoc, function(fn) {
    file.remove(file.path(pkg.path, "man", paste0(fn, ".Rd")))
  })
}


# check the pdf documentation with this command
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD", "Rd2pdf",
             shQuote(find.package(pkg.name))))


#--- vignette -------------------------------------------------------------------

# do this only once
# use_vignette("MADpop-tutorial")

# check vignette
rmarkdown::render(file.path(pkg.path, "vignettes", "MADPop-tutorial.Rmd"))

# run this when you want to "compile" the vignette
# will then be visible when you install but takes longer
build_vignettes(pkg = pkg.path)

#--- install package ------------------------------------------------------------

# run this when you want to install the package to see what the user gets.
# before running these steps, quit R, re-open, and only run the following:
#   1. require(devtools)
#   2. pkg.path <- getwd()
# then run:
file.remove(file.path(pkg.path, ".Rbuildignore"))
ignore.files <- c(".Rhistory", "Read-and-delete-me", "setup.R", "vignettes-old",
                  file.path("vignettes", "setup"))
use_build_ignore(files = ignore.files, pkg = pkg.path)
install(pkg = pkg.path)

# to then use the package, quit R, re-open, and run
require(MADPop)

#--- building the package -------------------------------------------------------

# this is what will get downloaded from CRAN
# you don't need to run this for testing/using the package
build(pkg = pkg.name)
