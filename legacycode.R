################################################################################
#                                 LEGACY CODE                                  #
################################################################################


################################################################################
# For determining the true parameters of the simulation
################################################################################

# set.seed(03092025)
# n <- 1e5
# n_sims <- 500
# simulation_list <- vector("list", n_sims)
# for (i in 1:n_sims){
#   if (i%%10 == 0){print(i)}
#   
#   simulation_list[[i]] <- trueATE(n = n)
# }

# ## drop the sim_id column if it exists
# results_mat <- simulation_list_df[ , !(names(simulation_list_df) %in% "sim_id") ]
# 
# ## column means
# col_means <- colMeans(results_mat, na.rm = TRUE)
# 
# ## column sds
# col_sds <- apply(results_mat, 2, sd, na.rm = TRUE)
# 
# ## combine into one summary data.frame
# summary_df <- data.frame(mean = col_means,
#                          sd   = col_sds)
# 
# summary_df
# simulation_list_df <- flatten_sims(simulation_list)
# write.csv(summary_df, file = "summary_results.csv", row.names = TRUE)
