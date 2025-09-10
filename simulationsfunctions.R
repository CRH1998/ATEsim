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
                                correlation = F, depcens = 0, deptreat = 0,
                                cluster_size = NA, 
                                varZ_corr = 1, varZ_corr_cens = 1, varZ_corr_treat = 1){
  
  
  #browser()
  
  # Dependent censoring?
  if (!is.numeric(depcens) || length(depcens) != 1 || !(depcens %in% c(0,1))) {
    stop("depcens must be 0 (independent censoring) or 1 (dependent censoring).")
  }
  depcens <- depcens != 0
  deptreat <- deptreat != 0
  
  
  
  ##############################################################################
  #                               Basic checks                                 #
  ##############################################################################
  
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("n must be a positive scalar.")

  if (correlation || depcens) {
    if (is.na(cluster_size)) stop("Please provide a cluster size if correlation = TRUE.")
    if (!is.numeric(cluster_size) || cluster_size < 1 || cluster_size != as.integer(cluster_size)) {
      stop("cluster_size must be a positive integer.")
    }
  } else {
    cluster_size <- 1L
  }
  
  if (varZ_corr <= 0) stop("varZ_corr must be > 0.")
  if (depcens && varZ_corr_cens <= 0) stop("varZ_corr_cens must be > 0.")
  
  ##############################################################################
  
  
  
  ##############################################################################
  #                   Determine clusters if data is correlated                 #
  ##############################################################################
  
  
  if (correlation || depcens){
    #-------Determine number of clusters--------#
    no_clusters <- no_of_cluster(n, cluster_size)
    
    n_clusters <- no_clusters$n_clusters
    n <- no_clusters$n
    #------------------------------------------#
  } else {
    n_clusters <- n
  }
  
  ##############################################################################
  
  
  
  
  ##############################################################################
  #                  Simulate correlation variables if chosen                  #
  ##############################################################################
  
  if (correlation){
    #------Simulate correlation variable-------#
    Z_corr <- rgamma(n_clusters, shape=1/varZ_corr)*varZ_corr # Change to any other distribution if needed
    log_Z_corr <- log(Z_corr)
    #------------------------------------------#
  }
    
  
  if (depcens){
    #------Simulate censoring covariate--------#
    Z_corr_cens <- rgamma(n_clusters, shape=1/varZ_corr_cens)*varZ_corr_cens # Change to any other distribution if needed
    log_Z_corr_cens <- log(Z_corr_cens)
    #------------------------------------------#
  }
  
  
  if (deptreat){
    #------Simulate treatment covariate--------#
    Z_corr_treat <- rgamma(n_clusters, shape = 1/varZ_corr_treat)*varZ_corr_treat # Change to any other distribution if needed
  }
  
  ##############################################################################
  
  
  
  
  
  ##############################################################################
  #                 Simulate covariates based on specifications                #
  ##############################################################################
  
  
  #-----------------------------------------------------#
  if (correlation){
    #--Simulate all covariates and store in matrices - reusing log_Z_corr for correlation--#
    Z <- vector("list", cluster_size)
    for (i in 1:cluster_size){
      
      covars <- simulate_covariates_specs(n = n_clusters, specs = specs)
      
      Z[[i]] <- cbind(log_Z_corr, covars) 
      
      if (depcens){
        #--If dependent censoring is specified add censoring variable
        Z[[i]] <- cbind(Z[[i]], log_Z_corr_cens)
      }
      
      if (deptreat){
        A <- correlated_treat(n = n_clusters, Z_corr_treat = Z_corr_treat,
                              varZ_corr_treat = varZ_corr_treat, X = covars)
        
        Z[[i]] <- cbind(A, Z[[i]])
      }
    }
  }
  #-----------------------------------------------------#
  
  
  
  #-----------------------------------------------------#
  else if (depcens){
    #--Simulate censoring covariates and store in matrix--#
    Z <- vector("list", cluster_size)
    for (i in 1:cluster_size){
      
      covars <- simulate_covariates_specs(n = n_clusters, specs = specs)
      
      Z[[i]] <- cbind(covars, log_Z_corr_cens)
      
      
      if (deptreat){
        A <- correlated_treat(n = n_clusters, Z_corr_treat = Z_corr_treat,
                              varZ_corr_treat = varZ_corr_treat, X = covars)
        
        Z[[i]] <- cbind(A, Z[[i]])
      }

    }
  } 
  #-----------------------------------------------------#
  
  
  #-----------------------------------------------------#
  else if (deptreat){
    
    for (i in 1:cluster_size){
      
      covars <- simulate_covariates_specs(n = n_clusters, specs = specs)
      
      A <- correlated_treat(n = n_clusters, Z_corr_treat = Z_corr_treat,
                            varZ_corr_treat = varZ_corr_treat, X = covars)
      
      Z[[i]] <- cbind(A, Z[[i]])
    }
  }
  #-----------------------------------------------------#
  
  
  #-----------------------------------------------------#
  else {
    n_clusters <- n
    Z <- simulate_covariates_specs(n = n_clusters, specs = specs)
  }
  #-----------------------------------------------------#
  
  ##############################################################################
  
  return(list("Z" = Z, "n" = n, "n_clusters" = n_clusters))
}






# Simulate correlated outcomes

#' Title
#'
#' @param n_clusters Integer, clusters per member (rows per Z[[i]]).
#' @param rho1 Baseline hazard parameter for cause 1.
#' @param rho2 Baseline hazard parameter for cause 2.
#' @param beta Numeric vector of length 2*p (stacked cause1, cause2).
#' @param rc Hazard parameter for gamma-distributed censoring.
#' @param depcens 0 for independent censoring; nonzero to enable dependent censoring.
#' @param rcZ Numeric vector of length p_cens with coefficients for censoring covariates.
#'            Ignored when \code{depcens == 0}.
#' @param type One of "cloglog" or "logistic".
#' @param cluster_size Optional. If missing, inferred from \code{length(Z)}.
#' @param rate Rate parameter for exponential baseline hazard.
#' @param Z List of length \code{cluster_size}; each element a data.frame/matrix with
#'          \code{n_clusters} rows and \code{p + p_cens} columns (p from \code{beta}).
#'          If \code{cluster_size = 1}, \code{Z} may be a single data.frame/matrix.
#'
#' @returns list(dats = data.frame, betaO = numeric, betaC = numeric)
#' @export
#'
#' @examples
corr.simul.cifs <- function(rho1 = 0.4, 
                            rho2 = 2, 
                            beta = c(1, 0.3, -0.3, 1, -0.3, 0.3), 
                            rc = 0.5, 
                            depcens = 0, 
                            deptreat = 0,
                            rcZ = 1, 
                            trZ = 0.3,
                            type = c("cloglog", "logistic"), 
                            rate = 1,
                            Z){
  
  #browser()
  
  n_clusters <- nrow(Z[[1]])
  cluster_size <- length(Z)

  
  ##############################################################################
  #                               Basic checks                                 #
  ##############################################################################
  
  type <- match.arg(type)
  if (!is.numeric(n_clusters) || length(n_clusters) != 1L || n_clusters < 1)
    stop("n_clusters must be a positive scalar integer.")
  n_clusters <- as.integer(n_clusters)
  
  if (!is.numeric(beta) || length(beta) < 2L || length(beta) %% 2L != 0L)
    stop("beta must be a numeric vector of even length (stacked for two causes).")
  p <- length(beta) %/% 2L
  
  
  ## Z can be a single df/matrix (cluster_size = 1) or a list of them
  if (!is.list(Z)) Z <- list(Z)
  if (is.null(cluster_size)) cluster_size <- length(Z)
  if (length(Z) != cluster_size)
    stop("length(Z) must equal cluster_size.")
  
  
  ## coerce Z to matrices, check dims, and that all have same ncol
  Z_mats <- lapply(seq_len(cluster_size), function(i) {
    Zi <- Z[[i]]
    if (is.data.frame(Zi)) Zi <- as.matrix(Zi)
    if (!is.matrix(Zi)) stop(sprintf("Z[[%d]] must be a matrix/data.frame.", i))
    if (nrow(Zi) != n_clusters)
      stop(sprintf("Z[[%d]] must have %d rows (n_clusters).", i, n_clusters))
    storage.mode(Zi) <- "double"
    Zi
  })
  k_all <- vapply(Z_mats, ncol, integer(1))
  if (length(unique(k_all)) != 1L)
    stop("All elements of Z must have the same number of columns.")
  k <- k_all[1L]
  
  
  ## dependent censoring bookkeeping
  depcens <- (depcens != 0)
  p_cens <- if (depcens) length(rcZ) else 0L
  
  ## dependent treatment bookkeeping
  deptreat <- (deptreat != 0)
  p_treat <- if (deptreat) length(trZ) else 0L
  if (deptreat) beta <- c(trZ, beta[seq_len(p)], trZ, beta[p + seq_len(p)]); p <- length(beta) %/% 2
  
  ## check Z column count matches p (+ p_cens when depcens)
  expected_k <- p + p_cens
  if (k != expected_k) {
    stop(sprintf("Each Z[[i]] must have p + p_cens columns: %d + %d = %d, but got %d.",
                 p, p_cens, expected_k, k))
  }
  
  ##############################################################################
  
  

  
  
  
  ##############################################################################
  #        Constructing outcome coefficients and censoring coefficients        #
  ##############################################################################
  
  if (depcens) {
    zeros <- rep(0, p_cens)
    beta_expanded <- c(beta[seq_len(p)], zeros,
                       beta[p + seq_len(p)], zeros)
    rcZ_full <- c(rep(0, p), rcZ)  # zeros for outcome covars, then censoring covars
  } else {
    beta_expanded <- beta
    rcZ_full <- numeric(0)         # ignored downstream when depcens == 0
  }
  
  ##############################################################################
  
  
  dats <- vector("list", cluster_size)
  
  for (i in 1:cluster_size){
    
    Z_i <- as.matrix(Z[[i]])
    
    dats[[i]] <- simul.cifs(n = n_clusters,          # Number of clusters = n / cluster_size
                            rho1 = rho1,             # Baseline hazard parameter for cause 1
                            rho2 = rho2,             # Baseline hazard parameter for cause 2
                            beta = beta_expanded,    # Coefficients for covariates
                            rc = rc,                 # Hazard for gamma distributed censoring
                            depcens = depcens,       # Independent censoring
                            rcZ = rcZ_full,          # Coefficients for censoring covariates
                            type = type,             # Simulate outcomes from specified model
                            rate = rate,             # Rate parameter for exponential baseline hazard
                            Z = Z_i)                 # Covariates for simulation
    
    
    # Add id column for clustering
    dats[[i]]$id <- 1:n_clusters
  }
  
  
  # Combine datasets - this needs to be generalized
  dats <- do.call(rbind, dats)
  
  return(list("dats" = dats, "betaO" = beta_expanded, "betaC" = rcZ_full))
}








# Runone simulation

#' Title
#'
#' @param n integer, number of rows to simulate
#' @param seed integer, for reproducibility
#' @param correlation logical, if TRUE simulates correlated outcomes
#' @param cluster_size integer, size of clusters if correlation = TRUE
#' @param varZ_corr integer, variance of the correlation variable
#' @param varZ_corr_cens integer, variance of the correlation variable for censoring
#' @param rcZ vector, coefficients for censoring covariates
#' @param specs named list describing each covariate.
#' @param rho1 numeric, baseline hazard parameter for cause 1
#' @param rho2 numeric, baseline hazard parameter for cause 2
#' @param beta vector, coefficients for covariates
#' @param rc numeric, hazard for gamma distributed censoring
#' @param depcens numeric, if not 0 simulates dependent censoring
#' @param sim_type character, simulate outcomes from specified model
#' @param trZ vector, coefficient for treatment variable
#' @param deptreat numeric, if not 0 simulates dependent treatment
#' @param formula string, formula to be passed to binregATE
#'
#' @returns
#' @export
#'
#' @examples
runone_modelspec <- function(n = 2000, 
                             seed = NA, 
                             correlation = T, 
                             cluster_size = 2,
                             varZ_corr = 1, 
                             varZ_corr_cens = 1, 
                             rcZ = 1,
                             trZ = 1,
                             specs = list(
                               Z1 = list(dist = rbinom, size = 1, prob = 0.5), 
                               Z2 = list(dist = rbinom, size = 1, prob = 0.5)),
                             rho1 = 0.4, 
                             rho2 = 2, 
                             beta = c(1, 0.3, -0.3, 1, -0.3, 0.3), 
                             rc = 0.5, 
                             depcens = 1, 
                             deptreat = 1,
                             sim_type = "cloglog",
                             formula = "Event(time,status) ~ Z3 + Z4",
                             treat_model = "Z3 ~ Z4",
                             cens_model = "~ strata(Z3, Z4)"){
  
  browser()
  
  if (depcens == 0){
    rcZ <- NULL
  }
  
  #-------------------------------Simulate data---------------------------------
  
  # Set seed if provided
  if(!is.na(seed)) set.seed(seed)
  
  
  # Simulate covariates
  sim_covars <- simulate_covariates(n = n, 
                                    correlation = correlation, 
                                    depcens = depcens, 
                                    deptreat = deptreat,
                                    specs = specs,
                                    cluster_size = cluster_size, 
                                    varZ_corr = varZ_corr, 
                                    varZ_corr_cens = varZ_corr_cens)
  
  Z <- sim_covars$Z
  n <- sim_covars$n
  n_clusters <- sim_covars$n_clusters
  
  
  
  if (!correlation){
    dats <- simul.cifs(n = n,
                       rho1 = rho1,
                       rho2 = rho2,
                       beta = beta,
                       rc = rc,
                       depcens = depcens,
                       type = sim_type,
                       Z = Z)
    
    dats$id <- 1:n  # Add id column for consistency
    
  } else {
    # Simulate outcomes using corr.simul.cifs
    sims <- corr.simul.cifs(rho1 = rho1,
                            rho2 = rho2,
                            beta = beta,
                            rc = rc,
                            depcens = depcens,
                            deptreat = deptreat,
                            rcZ = rcZ,
                            trZ = trZ,
                            type = sim_type,
                            Z = Z)
    
    dats <- sims$dats
    betaO <- sims$betaO
    betaC <- sims$betaC
  }
  
  
  # Changing binary columns to factor columns
  dats <- binary_to_factor(dats)
  
  #-------------------------------ATE estimation--------------------------------
  
  # With specified formula, treat.model and cens.model
  bATEt <- binregATE(formula = as.formula(formula), data = dats, model = "logit",
                     treat.model = as.formula(treat_model),
                     cens.model = as.formula(cens_model), 
                     time = c(4), cause = 1)
  

  
  #-------------------------------Extract estimates-----------------------------
  bATEtest <- extract_binregATE_coef(bATEt)
  
  
  return(list("n" = n, "ATEests" = bATEtest))
}













