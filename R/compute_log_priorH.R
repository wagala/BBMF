#' Compute the Log-Prior Probability of H
#'
#' Computes the log-prior probability of the binary matrix \code{H}
#' under independent Bernoulli prior probabilities specified by
#' \code{beta}.
#'
#' @param H A binary matrix of dimension \eqn{R \times G} containing
#'   the factor-specific feature patterns.
#' @param beta A numeric vector of length \eqn{G} containing the
#'   feature-specific Bernoulli prior probabilities.
#'
#' @return A numeric scalar containing the log-prior probability of
#'   \code{H}.
#'
#' @details
#' Each entry \eqn{H_{rg}} is assumed to follow a Bernoulli distribution
#' with success probability \eqn{\beta_g}. Thus, all entries in column
#' \eqn{g} share the same prior probability \eqn{\beta_g}.
#'
#' The log-prior is obtained by summing the Bernoulli log-probabilities
#' over all rows and columns of \code{H}.
#'
compute_log_priorH <- function(H, beta) {
  log_beta <- log(beta)
  log_one_minus_beta <- log(1 - beta)
  logpriorH <- sum(H %*% log_beta + (1 - H) %*% log_one_minus_beta)
  return(logpriorH)
}