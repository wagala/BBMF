#' Compute the Log-Prior Probability of W
#'
#' Computes the log-prior probability of the binary matrix \code{W}
#' under independent Bernoulli prior probabilities specified by
#' \code{alpha}.
#'
#' @param W A binary matrix whose entries indicate factor activations.
#' @param alpha A scalar, vector, or matrix of Bernoulli prior
#'   probabilities compatible with \code{W}.
#'
#' @return A numeric scalar containing the log-prior probability of
#'   \code{W}.
#'
#' @details
#' Each entry of \code{W} is assumed to follow a Bernoulli distribution
#' with success probability given by the corresponding value of
#' \code{alpha}. The total log-prior is obtained by summing the
#' element-wise Bernoulli log-probabilities over all entries of
#' \code{W}.
#'
compute_log_priorW <- function(W, alpha) {
  logpriorW <- sum(W * log(alpha) + (1 - W) * log(1 - alpha))
  return(logpriorW)
}