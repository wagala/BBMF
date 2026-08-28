#' Update the Feature-Specific Sparsity Indicators Psi
#'
#' Updates the binary sparsity indicators \eqn{\psi_g} conditional on
#' the feature-specific probabilities \eqn{\beta_g}, the global sparsity
#' probability \eqn{\pi}, and the parameters of the two Beta mixture
#' components.
#'
#' @param beta_g A numeric vector of length \eqn{G} containing the
#'   feature-specific probabilities \eqn{\beta_g}.
#' @param b1 First shape parameter of the Beta distribution associated
#'   with \eqn{\psi_g = 1}.
#' @param b2 Second shape parameter of the Beta distribution associated
#'   with \eqn{\psi_g = 1}.
#' @param c1 First shape parameter of the Beta distribution associated
#'   with \eqn{\psi_g = 0}.
#' @param c2 Second shape parameter of the Beta distribution associated
#'   with \eqn{\psi_g = 0}.
#' @param pi_s Current value of the global sparsity probability
#'   \eqn{\pi}.
#'
#' @return A binary vector of length \eqn{G} containing the updated
#'   feature-specific sparsity indicators \eqn{\psi_g}.
#'
#' @details
#' For each feature \eqn{g}, the function computes the posterior log
#' probabilities corresponding to \eqn{\psi_g = 1} and
#' \eqn{\psi_g = 0}.
#'
#' Conditional on \eqn{\psi_g = 1}, the feature-specific probability
#' \eqn{\beta_g} is associated with a Beta distribution having shape
#' parameters \code{b1} and \code{b2}. Conditional on
#' \eqn{\psi_g = 0}, it is associated with a Beta distribution having
#' shape parameters \code{c1} and \code{c2}.
#'
#' The Beta normalizing constants are included through \code{lbeta()}.
#' The two log-probabilities are normalized using \code{sumLog()}, and
#' the updated value of \eqn{\psi_g} is sampled from a Bernoulli
#' distribution.
#'
update_psi <- function(beta_g, b1, b2, c1, c2, pi_s) {
  # Initialize a vector to store updated psi_g values
  psi <- numeric(length(beta_g))
  
  # constants from the Beta normalizing factors
  const1 <- -lbeta(b1, b2)  # -log B(b1, b2)
  const0 <- -lbeta(c1, c2)  # -log B(c1, c2)
  
  # Iterate over each beta_g to update psi_g
  for (g in 1:length(beta_g)) {
    # Compute log probabilities for psi_g = 1 and psi_g = 0
    log_prob_1 <- const1 +(b1 - 1) * log(beta_g[g]) + (b2 - 1) * log(1 - beta_g[g]) + log(pi_s)
    log_prob_0 <- const0 + (c1 - 1) * log(beta_g[g]) + (c2 - 1) * log(1 - beta_g[g]) + log(1 - pi_s)
    
    # Calculate the normalized probability for psi_g = 1
    probPsi <- log_prob_1 - sumLog(c(log_prob_0,log_prob_1))
    prob<-exp(probPsi)
    
    # Draw from Bernoulli distribution to determine psi_g
    psi[g] <- stats::rbinom(1, 1, prob)
  }
  
  return(psi)
}