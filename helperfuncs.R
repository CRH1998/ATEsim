############################################################
#                                                          #
#   The following script includes a series of helper       #
#   functions for the ATE simulation study                 #
#                                                          #
############################################################




#---------------------------Simulation helper functions------------------------#


# Helper function to determine number of clusters given n and cluster size
no_of_cluster <- function(n, cluster_size){

  if (n %% cluster_size == 0){
    n_clusters <- n / cluster_size
  } else {
    warning("n is not a multiple of cluster_size, adjusting n to be a multiple of cluster_size")
    n <- ceiling(n / cluster_size) * cluster_size
    n_clusters <- n / cluster_size
  }
  
  return(list("n" = n, "n_clusters" = n_clusters))
}



# Helper function to allow the user to specify distribution of covariates
simulate_covariates_specs <- function(n, specs, seed = NULL) {
  if (!is.numeric(n) || length(n) != 1 || n < 1) stop("n must be a positive scalar.")
  
  # Optional reproducibility without permanently changing RNG state
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
  }
  
  out <- vector("list", length(specs))
  names(out) <- names(specs)
  
  j <- 0L
  for (nm in names(specs)) {
    j <- j + 1L
    s <- specs[[nm]]
    
    # 1) list with $dist and its args
    if (is.list(s) && !is.null(s$dist)) {
      f <- match.fun(s$dist)
      args <- s[setdiff(names(s), "dist")]
      val <- tryCatch(
        do.call(f, c(list(n), args)),                        # many RNGs take n as first arg
        error = function(e) do.call(f, c(list(n = n), args)) # others expect named 'n'
      )
      
      # 2) a function taking n and returning length n
    } else if (is.function(s)) {
      val <- s(n)
      
      # 3) atomic vector: length n (as-is), length 1 (repeat), or sample with replacement
    } else if (is.atomic(s)) {
      if (length(s) == n) {
        val <- s
      } else if (length(s) == 1L) {
        val <- rep(s, n)
      } else {
        val <- sample(s, n, replace = TRUE)
      }
      
    } else {
      stop(sprintf("Spec for '%s' is not understood.", nm))
    }
    
    if (length(val) != n) stop(sprintf("Spec '%s' returned length %d (expected %d).", nm, length(val), n))
    out[[j]] <- val
  }
  
  # data.frame
  df <- do.call(
    data.frame,
    c(out, list(check.names = TRUE, stringsAsFactors = FALSE))
  )
  as.matrix(df)
}














# Function to extract estimates and standard errors from a fitted model

extract_binregATE_coef <- function(model){
  model_summary <- summary(model)
  
  coef_est <- model_summary$coef[,1]
  coef_se <- model_summary$coef[,2]
  
  G_est <- model_summary$ateG[,1]
  G_se <- model_summary$ateG[,2]
  
  DR_est <- model_summary$ateDR[,1]
  DR_se <- model_summary$ateDR[,2]
  
  return(list(coef_est = coef_est,
              coef_se = coef_se,
              G_est = G_est,
              G_se = G_se,
              DR_est = DR_est,
              DR_se = DR_se))
}



# Flatten simulations
flatten_one <- function(sim) {
  v <- unlist(sim, use.names = TRUE)
  names(v) <- gsub("[()]", "", gsub("\\.", "_", names(v)))
  as.data.frame(as.list(v), check.names = FALSE)
}

flatten_sims <- function(sims){
  results_df <- do.call(rbind, lapply(sims, flatten_one))
  rownames(results_df) <- seq_len(nrow(results_df))
  
  return(results_df)  
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
#'
#' @returns
#' @export
#'
#' @examples
runone <- function(n = 800, 
                   seed = NA, 
                   correlation = T, 
                   cluster_size = 2,
                   varZ_corr = 1, 
                   varZ_corr_cens = 1, 
                   rcZ = 1,
                   specs = list(
                     Z1 = list(dist = rbinom, size = 1, prob = 0.5), 
                     Z2 = list(dist = rbinom, size = 1, prob = 0.5)),
                   rho1 = 0.4, 
                   rho2 = 2, 
                   beta = c(1, 0.3, -0.3, 1, -0.3, 0.3), 
                   rc = 0.5, 
                   depcens = 1, 
                   sim_type = "cloglog"){
  
  #browser()
  
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
    sims <- corr.simul.cifs(n_clusters = n_clusters,
                            rho1 = rho1,
                            rho2 = rho2,
                            beta = beta,
                            rc = rc,
                            depcens = depcens,
                            rcZ = rcZ,
                            type = sim_type,
                            cluster_size = cluster_size,
                            Z = Z)
    
    dats <- sims$dats
    betaO <- sims$betaO
    betaC <- sims$betaC
  }
  
  
  # Changing binary columns to factor columns
  dats <- binary_to_factor(dats)
  
  #-------------------------------ATE estimation--------------------------------
  
  # With specified cens.model but no cluster(id)
  bATEtBigCensNoID <- binregATE(Event(time,status) ~ Z2 + Z3, data = dats, model = "logit",
                                treat.model = Z2 ~ Z3,
                                cens.model = ~ strata(Z2, Z3), 
                                time = c(4), cause = 1)
  
  # With specified cens.model and cluster(id)
  bATEtBigCensWithID <- binregATE(Event(time,status) ~ Z2 + Z3 + cluster(id), data = dats, model = "logit",
                                  treat.model = Z2 ~ Z3,
                                  cens.model = ~ strata(Z2, Z3), 
                                  time = c(4), cause = 1)
  
  
  # Without specified cens.model and without cluster(id)
  bATEtSmallCensNoID <- binregATE(Event(time,status) ~ Z2 + Z3,data = dats, model = "logit",
                                  treat.model = Z2 ~ Z3,
                                  time = c(4), cause = 1)
  
  
  # Without specified cens.model but with cluster(id)
  bATEtSmallCensWithID <- binregATE(Event(time,status) ~ Z2 + Z3 + cluster(id), data = dats, model = "logit",
                                    treat.model = Z2 ~ Z3,
                                    time = c(4), cause = 1)
  
  
  atemodels <- list("BigCensNoID" = bATEtBigCensNoID, 
                    "BigCensWithID" = bATEtBigCensWithID, 
                    "SmallCensNoID" = bATEtSmallCensNoID, 
                    "SmallCensWithID" = bATEtSmallCensWithID)
  
  
  #-------------------------------Extract estimates-----------------------------
  atemodelests <- lapply(atemodels, extract_binregATE_coef)
  
  
  return(list("n" = n, "ATEests" = atemodelests))
}











# For estimating true ATE

trueATE <- function(n = 1e6, rho1 = 0.4, rho2 = 2, beta = c(1,0.3,-0.3,1,-0.3,0.3), sim_type = "cloglog"){
  
  # Simulate covariates
  sim_covars <- simulate_covariates(n = n, correlation = T, cluster_size = 1)
  
  Z <- sim_covars$Z
  n <- sim_covars$n
  
  if (ncol(Z[[1]]) != length(beta)/2){
    stop("Number of covariates does not match number of beta parameters")
  }
  
  # Simulate outcomes using simul.cifs
  dats <- simul.cifs(n = n,
                     rho1 = rho1,
                     rho2 = rho2,
                     beta = beta,
                     rc = 0,
                     depcens = 0,
                     type = sim_type,
                     Z = Z[[1]])
  
  dats$Z2 <- as.factor(dats$Z2)
  dats$Z3 <- as.factor(dats$Z3)
  
  # Estimate true ATE using binregATE
  bATEt <- binregATE(Event(time,status) ~ Z2 + Z3, data = dats, model = "logit",
                     treat.model = Z2 ~ Z3,
                     time = c(4), cause = 1)
  
  atemodelest <- extract_binregATE_coef(bATEt)
  
  return(atemodelest)
}



# helper: is a vector binary (exactly 2 distinct, ignoring NAs)?
is_binary <- function(x) {
  if (is.factor(x)) {
    return(nlevels(droplevels(x)) == 2L)
  }
  ux <- unique(x[!is.na(x)])
  length(ux) == 2L
}

# convert all binary columns to factor
binary_to_factor <- function(df) {
  if (!is.data.frame(df)) stop("df must be a data.frame")
  bin_cols <- vapply(df, is_binary, logical(1))
  df[bin_cols] <- lapply(df[bin_cols], function(x) {
    if (is.logical(x)) {
      # optional: label logicals nicely
      factor(x, levels = c(FALSE, TRUE), labels = c("FALSE", "TRUE"))
    } else if (is.factor(x)) {
      droplevels(x)  # ensure only the 2 levels remain
    } else {
      lv <- sort(unique(x[!is.na(x)]))
      factor(x, levels = lv)
    }
  })
  df
}

