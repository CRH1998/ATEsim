source("helperfuncs.R")

# Vælg antal simulationer
n_sim <- 10

# Kør simulationer og træk estimater ud
estimates <- replicate(n_sim, runone(n = 800, correlation = T, cluster_size = 2), simplify = FALSE) # Korrelerede observationer med to observationer per cluster

# Omformater simuleringerne i et data frame
sim_ests <- flatten_sims(estimates)