#' Update the Binary Factor-Activation Matrix W
#'
#' Performs a Gibbs update of each entry of the binary matrix \code{W}.
#' For each element \eqn{W_{kr}}, the function evaluates the
#' log-posterior probability under the two possible states,
#' \eqn{W_{kr}=0} and \eqn{W_{kr}=1}, normalizes the corresponding
#' posterior masses, and samples a new binary value.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations.
#' @param H A binary matrix containing the current factor-specific
#'   feature patterns.
#' @param X An observed binary data matrix.
#' @param alpha A scalar, vector, or matrix of Bernoulli prior
#'   probabilities compatible with \code{W}.
#' @param p11 Probability of observing an entry equal to 1 when the
#'   corresponding reconstructed entry is 1.
#' @param p10 Probability of observing an entry equal to 1 when the
#'   corresponding reconstructed entry is 0.
#'
#' @return An updated binary matrix \code{W} with the same dimensions
#'   as the input matrix.
#'
#' @details
#' The entries of \code{W} are updated sequentially. For each position
#' \eqn{(k,r)}, the function computes the log-posterior probability
#' after setting \eqn{W_{kr}} to 0 and then to 1.
#'
#' The two log-posterior values are normalized using \code{sumLog()}
#' to obtain the conditional probability that \eqn{W_{kr}=1}. A new
#' value is then sampled from a Bernoulli distribution.
#'
#' Because updates are sequential, each subsequent entry is updated
#' conditional on the most recently updated values of \code{W}.
#'
update_W <- function(W, H, X, alpha, p11, p10) {
  K <- nrow(W)
  R <- ncol(W)
  for (k in 1:K) {
    for (r in 1:R) {
      original_w_kr <- W[k, r]
      
      # Compute the log posterior when w_kr = 0
      W[k, r] <- 0
      log_posterior_0 <- compute_log_postW(W, H, X, alpha, p11, p10)
      
      # Compute the log posterior when w_kr = 1
      W[k, r] <- 1
      log_posterior_1 <- compute_log_postW(W, H, X, alpha, p11, p10)
      
      # Reset W to its original value
      W[k, r] <- original_w_kr
      
      # Compute the normalized probabilities
      log_propW <- log_posterior_1 - sumLog(c(log_posterior_1, log_posterior_0))
      p_W <- exp(log_propW) # Normalized Probability
      
      # Update W based on the normalized probabilities
      W[k, r] <- stats::rbinom(1, 1, p_W)
    }
  }
  return(W)
}