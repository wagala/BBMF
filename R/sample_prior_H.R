#' Sample the Factor-Pattern Matrix H from Its Prior
#'
#' Draws a fresh factor-pattern matrix \code{H} entrywise from its
#' Bernoulli prior, with success probability \eqn{\beta_g \phi_r} for
#' entry \eqn{H_{rg}}.
#'
#' @param beta A numeric vector of length \eqn{G} containing the
#'   per-feature activation probabilities \eqn{\beta_g}.
#' @param phi A numeric vector of length \eqn{R} containing the
#'   per-factor inclusion probabilities \eqn{\phi_r}.
#'
#' @return A binary (0/1) matrix of dimension \eqn{R \times G}.
#'
#' @details
#' Entries are drawn independently as
#' \eqn{H_{rg} \mid \beta_g, \phi_r \sim \mathrm{Bernoulli}(\beta_g \phi_r)}.
#'
#' @export
sample_prior_H <- function(beta, phi) {

  R <- length(phi)
  G <- length(beta)

  p_mat <- outer(phi, beta)
  H <- matrix(stats::rbinom(R * G, 1, p_mat), R, G)

  return(H)
}
