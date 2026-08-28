#' Compute the Log-Posterior Probability of H
#'
#' Computes the log-posterior probability of the binary matrix \code{H}
#' by combining its log-likelihood contribution with its log-prior
#' probability.
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
#' @return A numeric scalar containing the log-posterior probability
#'   of \code{H}, up to an additive normalizing constant.
#'
#' @details
#' The function computes the log-likelihood of the observed matrix
#' \code{X} using \code{compute_likelihood()} and adds the Bernoulli
#' log-prior contribution for \code{H} obtained from
#' \code{compute_log_priorH()}.
#'
compute_log_posteriorH <- function(W, H, X, beta, p11, p10) {
  loglikeH <- compute_likelihood(W, H, X, p11, p10)
  logpriorH <- compute_log_priorH(H, beta)
  return(loglikeH + logpriorH)
}