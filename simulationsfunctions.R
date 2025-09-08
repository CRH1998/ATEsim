############################################################
#                                                          #
#   The following script includes a series of              #
#   simulation functions for the ATE simulation study.     #
#                                                          #
############################################################


# Function to simulate covariates

#' Simulate covariates with user-specified distributions
#'
#' @param n integer, number of rows to simulate
#' @param specs named list describing each covariate.
#'   Each element can be ONE of:
#'   1) list(dist=<function or name>, <named args for that dist>)
#'   2) a function that takes n and returns a length-n vector
#'   3) an atomic vector to sample from (sampled with replacement unless length == n)
#'   4) a length-n atomic vector used as-is
#' @param seed optional integer for reproducibility
#' @param correlation logical, if TRUE simulates correlated outcomes
#' @param depcens numeric, if not 0 simulates correlated censoring
#' @param cluster_size integer, size of clusters if correlation = TRUE
#' @param varZ_corr numeric, variance of the correlation variable
#' @param varZ_corr_cens numeric, variance of the correlation variable for censoring
#' @return data.frame (or tibble) with n rows and length(specs) columns

simulate_covariates <- function(n, 
                                specs = list(
                                  Z1 = list(dist = rbinom, size = 1, prob = 0.5), 
                                  Z2 = list(dist = rbinom, size = 1, prob = 0.5)),
                                correlation = F, depcens = 0, 
                                cluster_size = NA, 
                                varZ_corr = 1, varZ_corr_cens = 1){
  
  
  
  # Basic checks
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("n must be a positive scalar.")

  if (correlation) {
    if (is.na(cluster_size)) stop("Please provide a cluster size if correlation = TRUE.")
    if (!is.numeric(cluster_size) || cluster_size < 1 || cluster_size != as.integer(cluster_size)) {
      stop("cluster_size must be a positive integer.")
    }
  } else {
    cluster_size <- 1L
  }
  
  if (varZ_corr <= 0) stop("varZ_corr must be > 0.")
  if (depcens != 0 && varZ_corr_cens <= 0) stop("varZ_corr_cens must be > 0.")
  
  
  
  if (correlation & is.na(cluster_size)){
    stop("Please provide a cluster size if correlation = TRUE")
  } else if (correlation) {
    
    #-------Determine number of clusters--------#
    no_clusters <- no_of_cluster(n, cluster_size)
    
    n_clusters <- no_clusters$n_clusters
    n <- no_clusters$n
    #------------------------------------------#
    
    #------Simulate correlation variable-------#
    Z_corr <- rgamma(n_clusters, shape=1/varZ_corr)*varZ_corr # Change to any other distribution if needed
    log_Z_corr <- log(Z_corr)
    #------------------------------------------#
    
    
    if (depcens != 0){
      #--Simulate censoring covariates and store in matrix--#
      Z_corr_cens <- rgamma(n_clusters, shape=1/varZ_corr_cens)*varZ_corr_cens # Change to any other distribution if needed
      log_Z_corr_cens <- log(Z_corr_cens)
      #-----------------------------------------------------#
    }
    
    
    
    #--Simulate all covariates and store in matrices - reusing z for correlation--#
    Z <- vector("list", cluster_size)
    for (i in 1:cluster_size){
      
      covars <- simulate_covariates_specs(n = n_clusters, specs = specs)
      
      Z[[i]] <- cbind(log_Z_corr, covars) 
      
      if (depcens != 0){
        Z[[i]] <- cbind(Z[[i]], log_Z_corr_cens)
      }
    }
    #----------------------------------------------------------------------------#
    
    
  } else {
    n_clusters <- n
    Z <- simulate_covariates_specs(n = n_clusters, specs = specs)
  }
  
  return(list("Z" = Z, "n" = n, "n_clusters" = n_clusters))
}






# Simulate correlated outcomes

#' Title
#'
#' @param n_clusters 
#' @param rho1 
#' @param rho2 
#' @param beta 
#' @param rc 
#' @param depcens 
#' @param rcZ 
#' @param type 
#' @param cluster_size 
#' @param rate 
#' @param Z 
#'
#' @returns
#' @export
#'
#' @examples
corr.simul.cifs <- function(n_clusters, 
                            rho1, 
                            rho2, 
                            beta, 
                            rc = 0.5, 
                            depcens = 0, 
                            rcZ = 0.5, 
                            type = c("cloglog", "logistic"), 
                            cluster_size,
                            rate = 1,
                            Z){
  
  #browser()
  # Error checks
  p <- length(beta)/2
  p_cens <- length(rcZ)
  
  if (ncol(Z[[1]]) != p + p_cens){
    stop("Number of covariates does not match number of parameters")
  }
  
  if (depcens != 0){
    
    # Construct rcZ with zeros for non-censoring covariates
    rcZ <- c(rep(0,length(beta)/2), rcZ)
    
    # Add zeros to beta for censoring covariates
    beta <- append(beta, 0, length(beta)/2)
    beta <- append(beta, 0)
  }
  
  
  dats <- vector("list", cluster_size)
  
  for (i in 1:cluster_size){
    
    Z_i <- as.matrix(Z[[i]])
    
    dats[[i]] <- simul.cifs(n = n_clusters,     # Number of clusters = n / cluster_size
                            rho1 = rho1,        # Baseline hazard parameter for cause 1
                            rho2 = rho2,        # Baseline hazard parameter for cause 2
                            beta = beta,        # Coefficients for covariates
                            rc = rc,            # Hazard for gamma distributed censoring
                            depcens = depcens,  # Independent censoring
                            rcZ = rcZ,          # Coefficients for censoring covariates
                            type = type,        # Simulate outcomes from specified model
                            rate = rate,        # Rate parameter for exponential baseline hazard
                            Z = Z_i)            # Covariates for simulation
    
    
    # Add id column for clustering
    dats[[i]]$id <- 1:n_clusters
  }
  
  
  # Combine datasets - this needs to be generalized
  dats <- do.call(rbind, dats)
  
  return(list("dats" = dats, betaO = beta, betaC = rcZ))
}