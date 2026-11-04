#' Log-Prior Density for the Factor-Activation Matrix W
#'
#' Computes the log-density of the factor-activation matrix \code{W}
#' under its Bernoulli prior, with success probability
#' \eqn{\alpha_k \phi_r} for entry \eqn{W_{kr}}.
#'
#' @param W A binary matrix of dimension \eqn{K \times R} containing
#'   the current factor activations.
#' @param alpha A numeric vector of length \eqn{K} containing the
#'   per-sample activation probabilities \eqn{\alpha_k}.
#' @param phi A numeric vector of length \eqn{R} containing the
#'   per-factor inclusion probabilities \eqn{\phi_r}.
#'
#' @return A numeric scalar containing the log-prior density of
#'   \code{W}.
#'
#' @details
#' Each entry is assumed independent with
#' \eqn{W_{kr} \mid \alpha_k, \phi_r \sim \mathrm{Bernoulli}(\alpha_k \phi_r)}.
#' The factor-level inclusion probability \eqn{\phi_r} allows an
#' entire factor to be switched off (\eqn{\phi_r \to 0}) independently
#' of the per-sample rate \eqn{\alpha_k}.
#'
#' @export
compute_log_prior_W <- function(W, alpha, phi) {

  K <- nrow(W)
  R <- ncol(W)

  p_mat <- outer(alpha, phi)  # entry [k,r] = alpha_k * phi_r

  logprior <- sum(W * log(p_mat) + (1 - W) * log(1 - p_mat))

  return(logprior)
}
