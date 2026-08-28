#' Update the Observation Probabilities
#'
#' Samples the observation probabilities \code{p11} and \code{p10}
#' from their full conditional Beta distributions given the current
#' factor matrices \code{W} and \code{H}.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations.
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the current factor-specific feature patterns.
#' @param X An observed binary matrix of dimension \eqn{K \times G}.
#'
#' @return A numeric vector of length two containing the updated values
#'   of \code{p11} and \code{p10}, in that order.
#'
#' @details
#' The function first reconstructs the binary matrix from \code{W} and
#' \code{H}. The parameter \code{p11} is the probability of observing
#' an entry equal to 1 when the corresponding reconstructed entry is 1,
#' while \code{p10} is the probability of observing an entry equal to 1
#' when the corresponding reconstructed entry is 0.
#'
#' Both parameters are sampled from their Beta full conditional
#' distributions using Beta(1, 1) priors.
#'
update_p <- function(W, H, X) {
  K <- nrow(W)
  G <- ncol(H)
  
  # Compute x_tilde using matrix multiplication and convert to binary
  X_tilde <- bool_prod(W, H)
  
  # Compute the values for p11 and p10 using vectorized operations
  m1 <- sum((1 - X_tilde) * X)
  m <- sum(1 - X_tilde)
  n1 <- sum(X_tilde * X)
  n <- sum(X_tilde)
  
  # Sample p11 and p10 from the Beta distribution
  p11 <- stats::rbeta(1, shape1 = n1 + 1, shape2 = n - n1 + 1)
  p10 <- stats::rbeta(1, shape1 = m1 + 1, shape2 = m - m1 + 1)
  return(c(p11, p10))
}