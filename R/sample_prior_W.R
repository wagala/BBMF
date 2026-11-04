#' Sample the Factor-Activation Matrix W from Its Prior
#'
#' Draws a fresh factor-activation matrix \code{W} entrywise from its
#' Bernoulli prior, with success probability \eqn{\alpha_k \phi_r} for
#' entry \eqn{W_{kr}}.
#'
#' @param alpha A numeric vector of length \eqn{K} containing the
#'   per-sample activation probabilities \eqn{\alpha_k}.
#' @param phi A numeric vector of length \eqn{R} containing the
#'   per-factor inclusion probabilities \eqn{\phi_r}.
#'
#' @return A binary (0/1) matrix of dimension \eqn{K \times R}.
#'
#' @details
#' Entries are drawn independently as
#' \eqn{W_{kr} \mid \alpha_k, \phi_r \sim \mathrm{Bernoulli}(\alpha_k \phi_r)}.
#'
#' @export
sample_prior_W <- function(alpha, phi) {

  K <- length(alpha)
  R <- length(phi)

  p_mat <- outer(alpha, phi)
  W <- matrix(stats::rbinom(K * R, 1, p_mat), K, R)

  return(W)
}
