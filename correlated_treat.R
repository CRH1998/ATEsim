



#' Title
#'
#' @param n total number of individuals
#' @param cluster_size individuals per cluster
#' @param varZ_corr_treat variance of correlation variable
#' @param X covariates
#' @param beta_coef multiplicative coefficients for covariates
#' @param alpha_coef additive coefficients for covariates
#'
#' @returns correlated treatment variable
#' @export
#'
#' @examples
correlated_treat <- function(n, Z_corr_treat, varZ_corr_treat, X,
                             beta_coef = NULL, alpha_coef = NULL){
  
  
  #browser()
  
  eta <- 1/varZ_corr_treat             # Reparameterization constant. Ensures mean 1 and variance 1/eta = varZ_corr_treat for random effect (frailty variable)

  
  
  # Define coefficients
  if (is.null(alpha_coef)){
    alpha_coef <- -0.5
  }
  if (is.null(beta_coef)){
    beta_coef <- rep(0.5, ncol(X))
  }
  
  
  
  
  ## --- Marginal target probabilities p_ik via logistic link ---
  # eta_{ik} = alpha_i + beta_i * x_{ik}
  linpred <- alpha_coef + X %*% beta_coef
  p_marg  <- plogis(linpred)  # p_{ik} in (0,1)
  
  
  # Computes the inverse laplace phi^{-1}(eta, p) = eta*(p^{-1/eta} - 1) used to get correct marginal probability
  inverse_laplace <- ilap(eta, p_marg)
  
  # Define conditional probability to get the correct marginal probability: E[exp(-Gam1 * inverse_laplace(p))] = laplace(inverse_laplace(p)) = p
  pgivenZ <- exp(-Z_corr_treat * inverse_laplace)
  
  
  # Draw binary outcomes with the correct probability. That is correct marginal probability and correlated through Z
  A <- rbinom(n, size = 1, prob = c(pgivenZ))

  return(A)
}
