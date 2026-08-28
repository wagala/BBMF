#' Update the Binary Factor Matrix H
#'
#' Performs a Gibbs update of each entry of the binary matrix \code{H}.
#' For each element \eqn{H_{rg}}, the function evaluates the
#' log-posterior probability under the two possible states,
#' \eqn{H_{rg}=0} and \eqn{H_{rg}=1}, normalizes the corresponding
#' posterior masses, and samples a new binary value.
#'
#' @param W A binary matrix containing the current factor activations.
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns.
#' @param X An observed binary data matrix.
#' @param beta A numeric vector of length \eqn{G} containing the
#'   feature-specific Bernoulli prior probabilities.
#' @param p11 Probability of observing an entry equal to 1 when the
#'   corresponding reconstructed entry is 1.
#' @param p10 Probability of observing an entry equal to 1 when the
#'   corresponding reconstructed entry is 0.
#'
#' @return An updated binary matrix \code{H} with the same dimensions
#'   as the input matrix.
#'
#' @details
#' The entries of \code{H} are updated sequentially. For each position
#' \eqn{(r,g)}, the function computes the log-posterior probability
#' after setting \eqn{H_{rg}} to 0 and then to 1.
#'
#' The two log-posterior values are normalized using \code{sumLog()}
#' to obtain the conditional probability that \eqn{H_{rg}=1}. A new
#' value is then sampled from a Bernoulli distribution.
#'
#' Because updates are sequential, each subsequent entry is updated
#' conditional on the most recently updated values of \code{H}.
#'
update_H <- function(W, H, X, beta, p11, p10) {
  G <- ncol(H)
  R <- nrow(H)
  for (r in 1:R) {
    for (g in 1:G) {
      original_h_rg <- H[r, g]
      
      # Compute the log posterior when h_rg = 0
      H[r, g] <- 0
      log_posterior_0 <- compute_log_posteriorH(W, H, X, beta, p11, p10)
      
      # Compute the log posterior when h_rg = 1
      H[r, g] <- 1
      log_posterior_1 <- compute_log_posteriorH(W, H, X, beta, p11, p10)
      
      # Reset H to its original value
      H[r, g] <- original_h_rg
      
      # Compute the normalized probabilities
      log_propH <- log_posterior_1 - sumLog(c(log_posterior_0, log_posterior_1))
      
      # Sample H
      p_H <- exp(log_propH)
      H[r, g] <- stats::rbinom(1, 1, p_H)
    }
  }
  return(H)
}