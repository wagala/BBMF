#' Compute the Log-Posterior Probability of W
#'
#' Computes the log-posterior probability of the binary matrix \code{W}
#' by combining its log-likelihood contribution with its log-prior
#' probability.
#'
#' @param W A binary matrix containing the current factor activations.
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
#' @return A numeric scalar containing the log-posterior probability
#'   of \code{W}, up to an additive normalizing constant.
#'
#' @details
#' The function computes the log-likelihood using
#' \code{compute_likelihood()} and adds the Bernoulli log-prior
#' contribution obtained from \code{compute_log_priorW()}.
#'
compute_log_postW <- function(W, H, X, alpha, p11, p10) {
  loglikeW <- compute_likelihood(W, H, X, p11, p10)
  logpriorW <- compute_log_priorW(W, alpha)
  return(loglikeW + logpriorW)
}