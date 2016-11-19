# Compile the Stan objects
if(!exists("maps.compile") || maps.compile) {
  message("--- COMPILING STAN CODE ---")
  message("HardyWeinberg.stan (1/2)")
  HW.mod <- stan_model(file = file.path("inst", "stan", "HardyWeinberg.stan"))
  message("DirichletMultinomial.stan (2/2)")
  hUM.mod <- stan_model(file = file.path("inst", "stan", "DirichletMultinomial.stan"))
  message("---------- DONE -----------")
}
