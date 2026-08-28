#' Compute the Bernoulli Log-Likelihood
#'
#' Computes the log-likelihood of an observed binary matrix \code{X}
#' given the current Boolean factorization defined by \code{W} and
#' \code{H}, and the observation probabilities \code{p11} and
#' \code{p10}.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations.
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns.
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#' @param p11 Probability of observing \eqn{X_{kg}=1} when the
#'   reconstructed entry \eqn{\tilde{X}_{kg}=1}.
#' @param p10 Probability of observing \eqn{X_{kg}=1} when the
#'   reconstructed entry \eqn{\tilde{X}_{kg}=0}.
#'
#' @return A numeric scalar containing the log-likelihood of the
#'   observed matrix \code{X}.
#'
#' @details
#' The reconstructed binary matrix is obtained using the Boolean
#' product of \code{W} and \code{H}. Conditional on the reconstructed
#' entry, the observed entry follows a Bernoulli distribution with
#' probability \code{p11} when the reconstructed entry is 1 and
#' probability \code{p10} when the reconstructed entry is 0.
#'
compute_likelihood <- function(W, H, X, p11, p10) {
  
  X_tilde <- bool_prod(W, H)
  
  loglike <- sum(
    X_tilde * X * log(p11) +
      X_tilde * (1 - X) * log(1 - p11) +
      (1 - X_tilde) * X * log(p10) +
      (1 - X_tilde) * (1 - X) * log(1 - p10)
  )
  
  return(loglike)
}